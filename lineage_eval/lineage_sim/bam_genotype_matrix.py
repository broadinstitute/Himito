#!/usr/bin/env python3
"""Reference prototype: genotype Himito-called variants per read from a BAM.

Oracle for the Rust `genotype_from_bam` path (plan 2026-10-11). For every
variant row of an existing Himito matrix (or a VCF), each read whose alignment
spans a window around the variant is re-scored against the REF haplotype and
each ALT haplotype of that site with an infix edit distance (haplotype fully
aligned, read-segment ends free). The read gets

    1   if the variant's haplotype is the unique best by >= --margin,
    0   if some other haplotype (REF or another ALT at the same site) is,
    NaN if it ties for best, the window is not fully spanned, or the best
        distance exceeds --max-ed (read too noisy here to say).

Overlapping variants (same REF span) are scored together as one multi-allelic
site. Reads with several alignments (supplementary pieces across the circular
origin, concatemers) are combined: decisive calls that disagree become NaN.

Output is a CSV in Himito's matrix format: `variant,<read1>,<read2>,...`, values
`1`/`0`, empty for NaN.

Usage:
  bam_genotype_matrix.py --bam raw.bam --ref rCRS.fasta \
      --variants sim.matrix.csv --out sim.bam_matrix.csv [--ref-column chrM]
"""
import argparse
import bisect
import csv
import re
import sys
from collections import defaultdict

import edlib
import numpy as np
import pysam

NAME_RE = re.compile(r"^m\.(\d+)([ACGTN]+)>([ACGTN]+)$")


def parse_name(name):
    """Himito names put the 1-based position of REF's first base after `m.`.

    SNP: `m.{pos0+1}{R}>{A}`; INS/DEL: `m.{VCF POS}{anchor..}>{anchor..}`, and
    VCF POS of an indel is the 1-based anchor, i.e. REF's first base too.
    """
    m = NAME_RE.match(name)
    if not m:
        return None
    return int(m.group(1)) - 1, m.group(2), m.group(3)


def read_variant_names(path):
    if path.endswith(".vcf") or path.endswith(".vcf.gz"):
        names = []
        for rec in pysam.VariantFile(path):
            for alt in rec.alts or ():
                names.append(f"m.{rec.pos}{rec.ref}>{alt}")
        return names
    with open(path) as fh:
        reader = csv.reader(fh)
        next(reader)
        return [row[0] for row in reader]


def group_sites(variants):
    """Group alleles of the same CLASS whose REF spans overlap into one site.

    SNPs group only with SNPs at the same base (multi-allelic SNV); indels group
    with overlapping indels (length alleles of one repeat). A SNP and an
    overlapping indel are deliberately scored apart: on ONT the indel is usually
    the deletion ERROR of the SNP base itself (e.g. m.8170AT>A over m.8171T>A),
    and scoring them jointly hands every read whose SNP base was deleted to the
    artifact -- a false 0 for the SNP -- where on its own the read ties REF and
    ALT and is correctly left NaN.
    """
    sites = []
    for is_indel in (False, True):
        idx = [i for i in range(len(variants))
               if (len(variants[i][1][1]) != len(variants[i][1][2])) == is_indel]
        idx.sort(key=lambda i: variants[i][1][0])
        cur, cur_end = [], -1
        for i in idx:
            start, ref, _alt = variants[i][1]
            end = start + len(ref)
            if cur and start < cur_end:
                cur.append(i)
                cur_end = max(cur_end, end)
            else:
                if cur:
                    sites.append(cur)
                cur, cur_end = [i], end
        if cur:
            sites.append(cur)
    return sites


def site_window(refseq, members, variants, flank):
    lo = min(variants[i][1][0] for i in members)
    hi = max(variants[i][1][0] + len(variants[i][1][1]) for i in members)
    left, right = lo - flank, hi + flank
    # Push each boundary out of any homopolymer it lands in, so an indel in a
    # run is never judged on a window that cuts the run.
    while 0 < left < len(refseq) and refseq[left] == refseq[left - 1]:
        left -= 1
    while 0 < right < len(refseq) and refseq[right] == refseq[right - 1]:
        right += 1
    return left - 2, right + 2


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--contig", default="chrM", help="contig to fetch from the BAM")
    ap.add_argument("--variants", required=True, help="Himito matrix.csv (row names) or VCF")
    ap.add_argument("--out", required=True)
    ap.add_argument("--flank", type=int, default=10)
    ap.add_argument("--slack", type=int, default=4)
    ap.add_argument("--margin", type=int, default=1)
    ap.add_argument("--max-ed", type=int, default=10**6)
    ap.add_argument("--ref-column", default=None,
                    help="also emit an all-0 column with this name (Himito's reference pseudo-read)")
    args = ap.parse_args()

    refseq = str(next(pysam.FastxFile(args.ref)).sequence).upper()
    names = read_variant_names(args.variants)
    variants = []
    for n in names:
        p = parse_name(n)
        if p is None:
            print(f"skip unparsable variant name {n}", file=sys.stderr)
            continue
        variants.append((n, p))

    sites = []  # (left, right, members, haps)
    for members in group_sites(variants):
        left, right = site_window(refseq, members, variants, args.flank)
        if left < 0 or right > len(refseq):
            continue  # window crosses the circular origin: left NaN (see plan)
        haps = {"REF": refseq[left:right]}
        for i in members:
            name, (s, r, a) = variants[i]
            haps[name] = refseq[left:s] + a + refseq[s + len(r):right]
        sites.append((left, right, members, haps))
    sites.sort(key=lambda x: x[0])
    site_lefts = [s[0] for s in sites]

    # One pass over the BAM: each alignment's matched (query, ref) pairs are
    # computed once and reused for every site it spans.
    per_site_read = defaultdict(lambda: defaultdict(list))  # site idx -> read -> best-sets
    bam = pysam.AlignmentFile(args.bam)
    for aln in bam.fetch(args.contig):
        if aln.is_unmapped or aln.is_secondary or aln.query_sequence is None:
            continue
        rs, re_ = aln.reference_start, aln.reference_end
        k0 = bisect.bisect_left(site_lefts, rs)
        if k0 >= len(sites) or sites[k0][0] >= re_:
            continue
        pairs = np.array(aln.get_aligned_pairs(matches_only=True), dtype=np.int64)
        qpos, rpos = pairs[:, 0], pairs[:, 1]
        seq, qlen = aln.query_sequence, aln.query_length
        for k in range(k0, len(sites)):
            left, right, members, haps = sites[k]
            if left >= re_:
                break
            if right > re_:
                continue
            a_idx = np.searchsorted(rpos, left, "left")
            b_idx = np.searchsorted(rpos, right, "left")
            if b_idx <= a_idx:
                continue
            qa = max(0, int(qpos[a_idx:b_idx].min()) - args.slack)
            qb = min(qlen, int(qpos[a_idx:b_idx].max()) + 1 + args.slack)
            seg = seq[qa:qb]
            d = {h: edlib.align(hs, seg, mode="HW")["editDistance"] for h, hs in haps.items()}
            best = min(d.values())
            if best > args.max_ed:
                per_site_read[k][aln.query_name].append(None)
                continue
            per_site_read[k][aln.query_name].append(
                frozenset(h for h, v in d.items() if v - best < args.margin))

    calls = defaultdict(dict)  # variant name -> read -> 0/1
    for k, per_read in per_site_read.items():
        members = sites[k][2]
        for read, sets in per_read.items():
            usable = [s for s in sets if s is not None]
            decisive = [s for s in usable if len(s) == 1]
            if not decisive:
                # ambiguous: still a 0 for any allele outside every best set
                if not usable:
                    continue
                for i in members:
                    name = variants[i][0]
                    if all(name not in s for s in usable):
                        calls[name][read] = 0
                continue
            winners = {next(iter(s)) for s in decisive}
            if len(winners) > 1:
                continue  # alignments of one read disagree
            win = winners.pop()
            for i in members:
                name = variants[i][0]
                calls[name][read] = 1 if win == name else 0

    reads = sorted({r for d in calls.values() for r in d})
    if args.ref_column:
        reads = sorted(set(reads) | {args.ref_column})
    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["variant"] + reads)
        for name, _ in variants:
            row = calls.get(name, {})
            out = []
            for r in reads:
                if args.ref_column and r == args.ref_column:
                    out.append("0")
                elif r in row:
                    out.append(str(row[r]))
                else:
                    out.append("")
            w.writerow([name] + out)
    n_cells = sum(len(d) for d in calls.values())
    print(f"wrote {args.out}: {len(variants)} variants x {len(reads)} reads, {n_cells} genotyped cells", file=sys.stderr)


if __name__ == "__main__":
    main()
