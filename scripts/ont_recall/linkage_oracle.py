#!/usr/bin/env python3
"""Read-linkage rescue oracle for low-frequency SNVs (ONT-recall revision 3).

For SNV alleles at the given positions, collect the reads carrying each allele
from a raw (pre-denoise) BAM pileup. Then test every allele pair for
co-occurrence on the same reads, and report which alleles the R3 rule
accepts. This is the reference behaviour for `src/linkage.rs` (R3-T2).

Rule (defaults = plan revision 3 decisions):
  candidates: SNV alleles with min_hf <= HF < max_hf and >= min_ac reads; NO strand gate
  reads:      the top 5% by error rate (NM / aligned length) are dropped from all counts
  pair test:  partners > min_dist bp apart, co-covered by >= 50 reads.
    overall:     shared >= 5, enrichment (shared / expected) >= 10, hypergeometric p <= 1e-10
    each strand: shared >= 3, enrichment >= 10, p <= 1e-5 (forward-strand reads alone,
                 then reverse-strand reads alone)
  accept:     an allele with >= 1 partner passing both.

Measured 2026-10-12 on HG00097_ont.mt.bam (HapMap merge, 50k reads, ~5,000x),
with 400 random non-truth positions plus all truth alleles under 20% HL:
noise rescued 0 of 760; HG00438 12/12; HG02486 4/5 (11087 has 6 reads);
HG002 5/8 (309, 310 and 13376 have no linked partner). Without the per-strand
requirement, 2 forward-only artifacts 25 bp apart (2542A, 2567A) were rescued.

usage: linkage_oracle.py --bam X.mt.bam --truth T.vcf.gz [--random 400] [--seed 7]
       (needs a .bai; ~0.75 s per position at 5,000x)
needs: pysam
"""
import argparse, collections, itertools, math, random
from math import lgamma, exp
import pysam

ap = argparse.ArgumentParser()
ap.add_argument("--bam", required=True); ap.add_argument("--truth", required=True)
ap.add_argument("--contig", default="chrM")
ap.add_argument("--random", type=int, default=400, help="random non-truth positions to sample as noise controls")
ap.add_argument("--seed", type=int, default=7)
ap.add_argument("--min-hf", type=float, default=0.004); ap.add_argument("--max-hf", type=float, default=0.2)
ap.add_argument("--noise-max-hf", type=float, default=0.05)
ap.add_argument("--min-ac", type=int, default=10)
ap.add_argument("--min-dist", type=int, default=20)
ap.add_argument("--err-quantile", type=float, default=0.95)
ap.add_argument("--log10p", type=float, default=-10.0); ap.add_argument("--strand-log10p", type=float, default=-5.0)
ap.add_argument("--min-enrich", type=float, default=10.0)
ap.add_argument("--min-shared", type=int, default=5); ap.add_argument("--min-shared-strand", type=int, default=3)
a = ap.parse_args()

bam = pysam.AlignmentFile(a.bam)
samples = ["HG002", "HG00438", "HG005", "HG02257", "HG02486", "HG02622"]
truth, near = {}, set()
for r in pysam.VariantFile(a.truth):
    for d in range(-3, len(r.ref) + 3):
        near.add(r.pos + d)
    if len(r.ref) == 1 and len(r.alts[0]) == 1:
        car = [s for s in samples if str(r.info.get(f"{s}_GT", "")).startswith("1")]
        truth[(r.pos, r.alts[0])] = (float(r.samples[0]["HL"][0]), car[0] if len(car) == 1 else "multi")

err = {}
for x in bam.fetch(a.contig):
    if x.is_secondary or not x.has_tag("NM") or x.query_alignment_length == 0:
        continue
    err[x.query_name] = max(err.get(x.query_name, 0.0), x.get_tag("NM") / x.query_alignment_length)
cut = sorted(err.values())[int(a.err_quantile * (len(err) - 1))] if err else 1.0
bad = {q for q, e in err.items() if e > cut}

random.seed(a.seed)
rand_pos = set(random.sample([p for p in range(60, bam.get_reference_length(a.contig) - 60) if p not in near], a.random))
positions = sorted(rand_pos | {p for (p, _), (hl, _) in truth.items() if hl < a.max_hf})

alle = {}
for pos in positions:
    cov, alts = {}, collections.defaultdict(set)
    for col in bam.pileup(a.contig, pos - 1, pos, truncate=True, min_base_quality=0, stepper="nofilter", max_depth=10**7):
        for p in col.pileups:
            if p.is_del or p.is_refskip or p.alignment.is_secondary:
                continue
            q = p.alignment.query_name
            if q in bad:
                continue
            cov[q] = p.alignment.is_reverse
            alts[p.alignment.query_sequence[p.query_position]].add(q)
    D = len(cov)
    for base, reads in alts.items():
        hf = len(reads) / max(D, 1); key = (pos, base)
        if key in truth and truth[key][0] < a.max_hf and len(reads) >= a.min_ac and hf >= a.min_hf:
            alle[key] = dict(reads=reads, cov=cov, hf=hf, kind="truth", lab=truth[key][1])
        elif pos in rand_pos and a.min_hf <= hf <= a.noise_max_hf and len(reads) >= a.min_ac:
            alle[key] = dict(reads=reads, cov=cov, hf=hf, kind="noise", lab="noise")

def log10_hyper_sf(k, N, K, n):
    lc = lambda n, k: lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)
    v = [lc(K, x) + lc(N - K, n - x) - lc(N, n) for x in range(k, min(K, n) + 1)]
    m = max(v); return (m + math.log(sum(exp(t - m) for t in v))) / math.log(10)

def test(x, y, strand=None, min_shared=5, log10p=-10.0):
    both = [q for q in x["cov"].keys() & y["cov"].keys() if strand is None or x["cov"][q] == strand]
    N = len(both)
    if N < (50 if strand is None else 30):
        return None
    bs = set(both); A = x["reads"] & bs; B = y["reads"] & bs; k = len(A & B)
    if k < min_shared:
        return None
    enr = k / (len(A) * len(B) / N)
    lp = log10_hyper_sf(k, N, len(A), len(B))
    return (k, enr, lp) if enr >= a.min_enrich and lp <= log10p else None

partners = collections.defaultdict(list)
for k1, k2 in itertools.combinations(alle, 2):
    if abs(k1[0] - k2[0]) <= a.min_dist:
        continue
    x, y = alle[k1], alle[k2]
    if len(x["reads"] & y["reads"]) < a.min_shared:
        continue
    o = test(x, y)
    if not o:
        continue
    if not (test(x, y, False, a.min_shared_strand, a.strand_log10p) and test(x, y, True, a.min_shared_strand, a.strand_log10p)):
        continue
    partners[k1].append((k2, o)); partners[k2].append((k1, o))

noise = [k for k in alle if alle[k]["kind"] == "noise"]
print(f"BAM {a.bam}: reads {len(err)} (dropped top {100*(1-a.err_quantile):.0f}% error > {cut:.3f}); positions {len(positions)}; candidate alleles {len(alle)} ({len(noise)} noise)")
print(f"noise alleles rescued: {sum(1 for k in noise if partners.get(k))} of {len(noise)}")
for g in sorted({v["lab"] for v in alle.values() if v["kind"] == "truth"}):
    ks = sorted(k for k in alle if alle[k]["kind"] == "truth" and alle[k]["lab"] == g)
    got = [k for k in ks if partners.get(k)]
    print(f"{g:<8} rescued {len(got)}/{len(ks)}; not rescued: {[f'{p}{b}' for p, b in ks if (p, b) not in got]}")
low_truth = sorted(k for k, (hl, _) in truth.items() if hl < a.max_hf and k not in alle)
print("truth SNVs below candidate gates (min_hf/min_ac):", [f"{p}{b}" for p, b in low_truth])
for k in sorted(alle):
    if partners.get(k):
        best = min(partners[k], key=lambda t: t[1][2])
        print(f"  accept {k[0]}{k[1]:<2} hf {alle[k]['hf']:.4f} {alle[k]['kind']:<5} {alle[k]['lab']:<8} partners {len(partners[k]):>2}; best {best[0][0]}{best[0][1]} shared {best[1][0]} enrich {best[1][1]:.0f}x log10p {best[1][2]:.1f}")
