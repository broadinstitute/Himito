#!/usr/bin/env python3
"""Cell-level accuracy of Himito variant x read matrices on the HapMap mixture.

Each read is assigned to a cell line from the private truth alleles it carries
(>= 2 private alt bases of one line, none of another), leaving out the variant
being scored. A cell's truth is 1 if that line carries the variant. Rates are
for comparing matrices, not absolute: only ~20% of reads are assignable.

Excluded: 513GCA>G (the truth spells HG005's deletion differently), 12969 and
14861 (discordant with the data), 3106 (REF N).

2026-10-12, merge215k ont-r10 (graph vs BAM-genotyped, max-ed off / 4):
  missing 16-18% vs 1.4-2.6% / 3-4%; minor-SNV FP 0.07% vs 0.57% / 0.34%.

usage: matrix_cell_eval.py --bam X.mt.bam(.bai) --truth T.vcf.gz LABEL=matrix.csv [LABEL=matrix.csv ...]
needs: pysam
"""
import argparse, collections, csv
import pysam

ap = argparse.ArgumentParser()
ap.add_argument("--bam", required=True); ap.add_argument("--truth", required=True)
ap.add_argument("--contig", default="chrM")
ap.add_argument("matrices", nargs="+", help="LABEL=path")
a = ap.parse_args()

LIN = ["HG005", "HG02622", "HG02257", "HG002", "HG02486", "HG00438"]
SKIP = {12969, 14861, 3106, 513}
bam = pysam.AlignmentFile(a.bam)
tv = []
for r in pysam.VariantFile(a.truth):
    car = [s for s in LIN if str(r.info.get(f"{s}_GT", "")).startswith("1")]
    tv.append((r.pos, r.ref, r.alts[0], float(r.samples[0]["HL"][0]), car))
priv = collections.defaultdict(list)
for pos, ref, alt, hl, car in tv:
    if len(ref) == 1 and len(alt) == 1 and len(car) == 1 and pos not in SKIP:
        priv[car[0]].append((pos, alt))
hits = collections.defaultdict(collections.Counter); sites = collections.defaultdict(set)
for L in LIN:
    for pos, alt in priv[L]:
        for col in bam.pileup(a.contig, pos - 1, pos, truncate=True, min_base_quality=0, stepper="nofilter", max_depth=10**7):
            for p in col.pileups:
                if not (p.is_del or p.is_refskip or p.alignment.is_secondary) and p.alignment.query_sequence[p.query_position] == alt:
                    q = p.alignment.query_name; hits[q][L] += 1; sites[q].add((L, (pos, alt)))
spans = collections.defaultdict(list)
for x in bam.fetch(a.contig):
    if not (x.is_secondary or x.is_unmapped):
        spans[x.query_name].append((x.reference_start, x.reference_end))

def lineage(q, exclude):
    c = collections.Counter(hits[q])
    for L in LIN:
        if (L, exclude) in sites[q]:
            c[L] -= 1
    c = {L: v for L, v in c.items() if v > 0}
    if not c:
        return None
    L, v = max(c.items(), key=lambda kv: kv[1])
    return L if v >= 2 and sum(c.values()) == v else None

def load(path):
    with open(path) as fh:
        rd = csv.reader(fh); hdr = next(rd)[1:]
        return {q: i for i, q in enumerate(hdr)}, {r[0]: r[1:] for r in rd}

mats = {lab: load(p) for lab, p in (m.split("=", 1) for m in a.matrices)}
truth = {f"m.{pos}{ref}>{alt}": (pos, ref, alt, hl, car) for pos, ref, alt, hl, car in tv}
names = set().union(*(set(r) for _, r in mats.values()))
cells = []
for name in names:
    if name not in truth or truth[name][0] in SKIP:
        continue
    pos, ref, alt, hl, car = truth[name]
    typ = "SNV" if len(ref) == 1 and len(alt) == 1 else "indel"
    cls = ("major " if hl >= 0.5 else "minor ") + typ
    for q, sp in spans.items():
        if not any(s <= pos - 13 and e >= pos + 12 for s, e in sp):
            continue
        L = lineage(q, (pos, alt) if typ == "SNV" else None)
        if L is not None:
            cells.append((cls, name, q, 1 if L in car else 0))
print(f"cells {len(cells)} over {len({c[1] for c in cells})} truth variants; assigned reads {sum(1 for q in spans if lineage(q, None))}/{len(spans)}")
print(f"{'class':<12} {'matrix':<14} {'true1':>7} {'false1':>7} {'prec(1)':>8} {'FN':>7} {'FP':>7} {'missing':>8}")
for cls in ("major SNV", "major indel", "minor SNV", "minor indel"):
    for lab, (idx, rows) in mats.items():
        c = collections.Counter()
        for cl, name, q, t in cells:
            if cl == cls:
                v = rows[name][idx[q]] if (name in rows and q in idx) else ""
                c[(t, v or "NA")] += 1
        t1, f0, fp, t0 = c[(1, "1")], c[(1, "0")], c[(0, "1")], c[(0, "0")]
        tot = sum(c.values())
        if tot:
            print(f"{cls:<12} {lab:<14} {t1:>7} {fp:>7} {t1/max(t1+fp,1):>8.1%} {f0/max(t1+f0,1):>7.2%} {fp/max(fp+t0,1):>7.2%} {(c[(1,'NA')]+c[(0,'NA')])/tot:>8.1%}")
    print()
