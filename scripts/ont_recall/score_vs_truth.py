#!/usr/bin/env python3
"""Exact-match scorer for Himito chrM VCFs against the SMaHT HapMap mito truth.

Counts PASS records whose (POS, REF, ALT) equals a truth record. On the two
2026-10-11 HapMap ONT runs this reproduces `rtg vcfeval --squash-ploidy`
exactly (ont-r10: 37/1/84, ont-denoised: 79/1/42), so it is a valid quick
proxy. Use vcfeval for anything you publish.

usage: score_vs_truth.py --truth T.vcf.gz VCF [VCF ...] [--exclude-pos 3106]
needs: pysam
"""
import argparse, collections, os
import pysam

ap = argparse.ArgumentParser()
ap.add_argument("--truth", required=True)
ap.add_argument("--exclude-pos", type=int, nargs="*", default=[],
                help="truth POS to drop before scoring (e.g. 3106: REF contains N, never callable)")
ap.add_argument("vcfs", nargs="+")
a = ap.parse_args()

truth = {}
for r in pysam.VariantFile(a.truth):
    if r.pos in a.exclude_pos:
        continue
    truth[(r.pos, r.ref, r.alts[0])] = round(float(r.samples[0]["HL"][0]), 3)

for vcf in a.vcfs:
    calls = {}
    for line in open(vcf):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        calls[(int(f[1]), f[3], f[4])] = (f[6], f[9])
    pas = {k for k, v in calls.items() if v[0] == "PASS"}
    tp = pas & truth.keys()
    fp = pas - tp
    fn = truth.keys() - tp
    bins = collections.defaultdict(lambda: [0, 0])
    for k, h in truth.items():
        bins[h][0 if k in tp else 1] += 1
    print(f"== {os.path.basename(vcf)}  (truth n={len(truth)})")
    print(f"   TP {len(tp)}  FP {len(fp)}  FN {len(fn)}  "
          f"recall {len(tp)/len(truth):.4f}  precision {len(tp)/max(len(pas),1):.4f}")
    print("   by truth HL (TP, FN): " + ", ".join(f"{h}:{tuple(v)}" for h, v in sorted(bins.items(), reverse=True)))
    print("   FP: " + ", ".join(f"{p}{r}>{x}" for p, r, x in sorted(fp)))
    nonpass = [(k, truth[k], calls[k][0]) for k in sorted(truth) if k in calls and k not in pas]
    print("   truth called but not PASS: " + ", ".join(f"{p}{r}>{x}(HL {h}, {flt})" for (p, r, x), h, flt in nonpass))
