#!/usr/bin/env python3
"""Per-site read evidence at low-HL truth SNVs, raw vs denoised BAM.

Shows what each FN looks like in the reads: alt count, strands and HF in the
raw (pre-denoise) BAM and in the denoised BAM, plus the Himito call. This is
how the 2026-10-11 review found that denoise erases real 2-4% alleles
(e.g. 9545A>G: 15 raw alt reads -> 2 after denoise).

usage: fn_pileup_audit.py --truth T.vcf.gz --raw X.mt.bam [--denoised X.mt.denoised.bam]
                          [--vcf X.vcf] [--max-hl 0.2] [--contig chrM]
needs: pysam (indexes temporary copies when a BAM has no index)
"""
import argparse, os, shutil, tempfile
import pysam

ap = argparse.ArgumentParser()
ap.add_argument("--truth", required=True)
ap.add_argument("--raw", required=True)
ap.add_argument("--denoised")
ap.add_argument("--vcf")
ap.add_argument("--max-hl", type=float, default=0.2)
ap.add_argument("--contig", default="chrM")
a = ap.parse_args()
tmp = tempfile.mkdtemp()

def open_indexed(path):
    try:
        b = pysam.AlignmentFile(path); b.fetch(a.contig, 0, 1); return b
    except ValueError:
        dst = os.path.join(tmp, os.path.basename(path)); shutil.copy(path, dst)
        pysam.index(dst); return pysam.AlignmentFile(dst)

def evidence(bam, pos, alt_base):
    d = f = r = 0
    for col in bam.pileup(a.contig, pos - 1, pos, truncate=True, min_base_quality=0,
                          stepper="nofilter", max_depth=10**6):
        for p in col.pileups:
            if p.is_del or p.is_refskip or p.alignment.is_secondary:
                continue
            d += 1
            if p.alignment.query_sequence[p.query_position] == alt_base:
                if p.alignment.is_reverse: r += 1
                else: f += 1
    return d, f, r

calls = {}
if a.vcf:
    for line in open(a.vcf):
        if not line.startswith("#"):
            x = line.split("\t"); calls[(int(x[1]), x[3], x[4])] = (x[6], x[9].strip().split(":")[-1])
bams = [("raw", open_indexed(a.raw))] + ([("den", open_indexed(a.denoised))] if a.denoised else [])
print("pos   ref>alt  HL     carriers            " + "  ".join(f"{n}: DP alt(f/r) hf   " for n, _ in bams) + " call")
samples = ["HG002", "HG00438", "HG005", "HG02257", "HG02486", "HG02622"]
for rec in pysam.VariantFile(a.truth):
    hl = float(rec.samples[0]["HL"][0])
    if hl > a.max_hl or len(rec.ref) != 1 or len(rec.alts[0]) != 1:
        continue
    carriers = ",".join(s for s in samples if str(rec.info.get(f"{s}_GT", "")).startswith("1"))
    row = f"{rec.pos:<5} {rec.ref}>{rec.alts[0]}      {hl:<6.3f} {carriers:<19}"
    for _, b in bams:
        d, f, r = evidence(b, rec.pos, rec.alts[0])
        row += f"  {d:>4} {f+r:>3}({f}/{r}) {(f+r)/max(d,1):.3f}"
    row += "  " + str(calls.get((rec.pos, rec.ref, rec.alts[0]), "-"))
    print(row)
shutil.rmtree(tmp, ignore_errors=True)
