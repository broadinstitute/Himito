# Plan 2026-09-27 T7
"""Measure homopolymer length confusion from a primary-alignment BAM.

The length allele is ``denoise_indel::hp_allele``: the read must cover the
anchor and ``INDEL_FLANK`` bases past the run, an insertion inside the run
counts only when every inserted base is the run base, and a deletion inside
the run shortens it. ``obs_len = ref_len + allele``.
"""
import argparse
import bisect
import pathlib
import sys
from collections import Counter, defaultdict

import pysam

from hp_tracts import CATALOGUE, find_homopolymer_runs

INDEL_FLANK = 5

_MATCH = (0, 7, 8)  # M, =, X
_INS = 1
_DEL = 2
_REF_SKIP = 3
_SOFT = 4


def overlap(s1, e1, s2, e2):
    """Size of the overlap of half-open intervals ``[s1, e1)`` and ``[s2, e2)``."""
    return max(0, min(e1, e2) - max(s1, s2))


def match_run_len(pos, cigar, target):
    """Aligned bases of the match op covering ``target``, else 0.

    A deletion or ref-skip that covers ``target`` reports 0, matching
    ``denoise_indel::match_run_len``.
    """
    rp = max(pos, 0)
    for op, n in cigar:
        if op in _MATCH:
            if rp <= target < rp + n:
                return rp + n - target
            rp += n
        elif op in (_DEL, _REF_SKIP):
            if rp <= target < rp + n:
                return 0
            rp += n
    return 0


def hp_allele(pos, cigar, seq, refseq, a, e):
    """Net length change of the run ``[a + 1, e)`` on one alignment.

    ``pos`` is the 0-based leftmost reference coordinate. ``a`` is the 0-based
    anchor (the base before the run) and ``e`` is the 0-based exclusive end.
    ``None`` means the read is uninformative, not that it observes the
    reference length.
    """
    if e < a + 2 or a + 1 >= len(refseq):
        return None
    base = refseq[a + 1].upper()
    if max(pos, 0) > a:
        return None
    if match_run_len(pos, cigar, a) == 0:
        return None

    rp = max(pos, 0)
    qp = 0
    last_aligned = None
    k = 0
    seq_u = seq.upper()
    for op, n in cigar:
        if op in _MATCH:
            rp += n
            qp += n
            if n > 0:
                last_aligned = rp - 1
        elif op == _INS:
            inside = last_aligned is not None and a <= last_aligned < e
            if inside:
                chunk = seq_u[qp:qp + n]
                if len(chunk) != n or any(c != base for c in chunk):
                    return None
                k += n
            qp += n
        elif op == _DEL:
            k -= overlap(rp, rp + n, a + 1, e)
            rp += n
        elif op == _REF_SKIP:
            if overlap(rp, rp + n, a + 1, e) > 0:
                return None
            rp += n
        elif op == _SOFT:
            qp += n
    if rp < e + INDEL_FLANK:
        return None
    return k


def load_fasta(path):
    """``{contig: sequence}`` from a FASTA file. Sequences are uppercased."""
    seqs = {}
    name = None
    chunks = []

    def flush():
        if name is not None:
            seqs[name] = "".join(chunks).upper()

    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        flush()
    if not seqs:
        raise SystemExit(f"no sequences in {path}")
    return seqs


def parse_tract_ids(spec):
    """Tract ids to skip.

    ``None`` skips nothing. The token ``hp_tracts`` skips the catalogue.
    An existing file contributes its first column (a ``tract`` header is
    ignored). Anything else is a comma-separated id list.
    """
    if spec is None or spec == "":
        return set()
    if spec == "hp_tracts":
        return {tract.id for tract in CATALOGUE}
    path = pathlib.Path(spec)
    if path.is_file():
        ids = set()
        for line in path.read_text().splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            field = line.split("\t")[0].strip()
            if field.lower() in {"tract", "id", "tract_id"}:
                continue
            ids.add(field)
        return ids
    if "/" in spec or spec.endswith((".txt", ".tsv", ".csv")):
        raise SystemExit(f"tract file not found: {spec}")
    return {part.strip() for part in spec.split(",") if part.strip()}


class RunIndex:
    """Homopolymer runs laid out so one CIGAR walk can score every run.

    ``alleles`` returns the same ``{run_index: k}`` map ``hp_allele`` would
    for each run on this reference. Positions that no run cares about are
    absent from the maps.
    """

    def __init__(self, runs):
        self.runs = runs
        order = sorted(range(len(runs)), key=lambda i: runs[i][0])
        self._order = order
        self._anchors = [runs[i][0] for i in order]
        self._ins_at = defaultdict(list)
        self._del_at = defaultdict(list)
        for i, (anchor, exclusive, base, _length, _start) in enumerate(runs):
            for pos in range(anchor, exclusive):
                self._ins_at[pos].append(i)
            for pos in range(anchor + 1, exclusive):
                self._del_at[pos].append(i)

    def alleles(self, pos, cigar, seq):
        anchor_ok = set()
        dead = set()
        k = {}
        rp = max(pos, 0)
        qp = 0
        last_aligned = None
        seq_u = seq.upper()
        for op, n in cigar:
            if op in _MATCH:
                if n:
                    left = bisect.bisect_left(self._anchors, rp)
                    right = bisect.bisect_left(self._anchors, rp + n)
                    for j in range(left, right):
                        anchor_ok.add(self._order[j])
                    rp += n
                    last_aligned = rp - 1
                qp += n
            elif op == _INS:
                if last_aligned is not None:
                    chunk = seq_u[qp:qp + n]
                    intact = len(chunk) == n
                    for idx in self._ins_at.get(last_aligned, ()):
                        base = self.runs[idx][2]
                        if not intact or any(c != base for c in chunk):
                            dead.add(idx)
                        else:
                            k[idx] = k.get(idx, 0) + n
                qp += n
            elif op == _DEL:
                for p in range(rp, rp + n):
                    for idx in self._del_at.get(p, ()):
                        k[idx] = k.get(idx, 0) - 1
                rp += n
            elif op == _REF_SKIP:
                for p in range(rp, rp + n):
                    dead.update(self._del_at.get(p, ()))
                rp += n
            elif op == _SOFT:
                qp += n
        out = {}
        for idx in anchor_ok:
            if idx in dead:
                continue
            exclusive = self.runs[idx][1]
            if rp < exclusive + INDEL_FLANK:
                continue
            out[idx] = k.get(idx, 0)
        return out


def _runs(seq, min_len, skip_intervals):
    runs = []
    for start, end, base, length in find_homopolymer_runs(seq, min_len):
        if (start, end) in skip_intervals:
            continue
        anchor = start - 2  # 0-based base immediately before the run
        if anchor < 0:
            continue
        runs.append((anchor, end, base, length, start))
    return runs


def _primary(aln):
    return not (aln.is_unmapped or aln.is_secondary or aln.is_supplementary)


def measure(bam, ref, min_len=4, tracts=None):
    """Rows ``(base, ref_len, obs_len, strand, n)`` for reference HP runs.

    A run is omitted when it is named by ``tracts``, or when the unique modal
    observed length is not ``ref_len`` (the sample is not reference there,
    including a tie). Primary alignments only.
    """
    seqs = load_fasta(ref)
    skip_ids = parse_tract_ids(tracts)
    skip_intervals = {
        (tract.start, tract.end) for tract in CATALOGUE if tract.id in skip_ids
    }
    runs_by_name = {}
    for name, seq in seqs.items():
        runs = _runs(seq, min_len, skip_intervals)
        runs_by_name[name] = (runs, RunIndex(runs))

    # (base, ref_len, start) -> list[(obs_len, strand)]
    observed = defaultdict(list)
    with pysam.AlignmentFile(bam, "rb") as bam_fh:
        name_map = {name: name for name in bam_fh.references if name in seqs}
        if not name_map and len(seqs) == 1 and len(bam_fh.references) == 1:
            name_map[bam_fh.references[0]] = next(iter(seqs))

        def consume(aln, ref_name):
            if not _primary(aln):
                return
            runs, index = runs_by_name[ref_name]
            start = aln.reference_start
            qseq = aln.query_sequence
            cigar = aln.cigartuples
            if start is None or not qseq or not cigar:
                return
            strand = "-" if aln.is_reverse else "+"
            for idx, allele in index.alleles(start, cigar, qseq).items():
                _anchor, _exclusive, base, length, run_start = runs[idx]
                observed[(base, length, run_start)].append((length + allele, strand))

        if bam_fh.has_index():
            for bam_name, ref_name in name_map.items():
                for aln in bam_fh.fetch(bam_name):
                    consume(aln, ref_name)
        else:
            for aln in bam_fh.fetch(until_eof=True):
                ref_name = name_map.get(aln.reference_name)
                if ref_name is not None:
                    consume(aln, ref_name)

    pooled = Counter()
    for (base, ref_len, _start), items in observed.items():
        counts = Counter(obs_len for obs_len, _strand in items)
        best = max(counts.values())
        modes = [obs_len for obs_len, n in counts.items() if n == best]
        if len(modes) != 1 or modes[0] != ref_len:
            continue
        for obs_len, strand in items:
            pooled[(base, ref_len, obs_len, strand)] += 1
    return [(b, rl, ol, st, n) for (b, rl, ol, st), n in sorted(pooled.items())]


def write_confusion(rows, path):
    with open(path, "w") as fh:
        fh.write("base\tref_len\tobs_len\tstrand\tn\n")
        for base, ref_len, obs_len, strand, n in rows:
            fh.write(f"{base}\t{ref_len}\t{obs_len}\t{strand}\t{n}\n")


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Measure homopolymer length confusion from a BAM (plan T7)."
    )
    parser.add_argument("--bam", required=True)
    parser.add_argument("--ref", required=True)
    parser.add_argument("--min-len", type=int, default=4)
    parser.add_argument(
        "--tracts",
        default=None,
        help="skip these tract ids: hp_tracts, a comma list, or a file",
    )
    parser.add_argument("--output", "-o", default="hp_confusion.tsv")
    args = parser.parse_args(argv)
    rows = measure(args.bam, args.ref, min_len=args.min_len, tracts=args.tracts)
    write_confusion(rows, args.output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
