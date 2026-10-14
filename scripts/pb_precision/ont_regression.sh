#!/usr/bin/env bash
# ONT non-regression gate for the PacBio precision work (plan 2026-10-13-pacbio-precision.md, T0).
#
# freeze : run every ONT case TWICE with the current binary, require the two runs to be byte-identical
#          (otherwise the gate would be flaky), then store the outputs as the golden set.
# check  : run every ONT case once with the current binary and `cmp` each output against the golden set.
#          Exit 1 on any difference and print the first differing records.
#
# usage (from the repo root or anywhere):
#   scripts/pb_precision/ont_regression.sh freeze          # once, at the starting commit, before any change
#   scripts/pb_precision/ont_regression.sh check           # after every task
#   scripts/pb_precision/ont_regression.sh check --full    # also the end-to-end quick-start runs (slow; T6)
#
# env overrides: HIMITO BENCH REF GOLDEN CASES (space-separated subset of: r10 r9 denoised r10_ds100 r10_ds1000)
#                FULL_BAM (input of the --full quick-start runs; default HapMap_ont_merge.bam)
#
# Cases (call-only, on the frozen inputs, so a difference can only come from `call`):
#   r10        -d ont-r10      on HG00097_r10.gfa + HG00097_r10.mt.bam
#   r9         -d ont-r9       on the same inputs (no r9 data; exercises the r9 preset branch)
#   denoised   -d ont-denoised on HG00097_ont.gfa + HG00097_ont.mt.bam (the raw, pre-denoise mt BAM)
#   r10_ds100  -d ont-r10      on a ~100x subset of HG00097_r10.mt.bam (subset + its GFA are frozen)
#   r10_ds1000 -d ont-r10      on a ~1000x subset
# No -m/-v is passed, so the data-type presets are what gets exercised. --strand-bias-threshold 0.05
# matches quick-start.
# --full adds quick-start (filter -> denoise -> build -> call) for ont-r10 and ont-denoised through
# scripts/ont_recall/run_hapmap_bench.sh; compares the VCF and the vcfeval summary.
set -euo pipefail
MODE=${1:?usage: ont_regression.sh freeze|check [--full]}
FULL=${2:-}
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
HIMITO=${HIMITO:-$ROOT/target/release/Himito}
BENCH=${BENCH:-$ROOT/ont_benchmark}
REF=${REF:-$ROOT/rCRS.fasta}
GOLDEN=${GOLDEN:-$BENCH/golden}
CASES=${CASES:-"r10 r9 denoised r10_ds100 r10_ds1000"}
FULL_BAM=${FULL_BAM:-HapMap_ont_merge.bam}
OUTS="vcf matrix.csv raw_matrix.csv"

[[ -x "$HIMITO" ]] || { echo "no Himito binary at $HIMITO (cargo build --release)"; exit 2; }
need_index() { [[ -f "$1.bai" || -f "${1%.bam}.bai" ]] || samtools index "$1"; }

case_inputs() {   # case -> "datatype gfa bam"
  case $1 in
    r10)        echo "ont-r10 $BENCH/HG00097_r10.gfa $BENCH/HG00097_r10.mt.bam" ;;
    r9)         echo "ont-r9 $BENCH/HG00097_r10.gfa $BENCH/HG00097_r10.mt.bam" ;;
    denoised)   echo "ont-denoised $BENCH/HG00097_ont.gfa $BENCH/HG00097_ont.mt.bam" ;;
    r10_ds100)  echo "ont-r10 $GOLDEN/inputs/r10_ds100.gfa $GOLDEN/inputs/r10_ds100.bam" ;;
    r10_ds1000) echo "ont-r10 $GOLDEN/inputs/r10_ds1000.gfa $GOLDEN/inputs/r10_ds1000.bam" ;;
    *) echo "unknown case $1" >&2; exit 2 ;;
  esac
}

make_subsets() {  # frozen once, at freeze time
  mkdir -p "$GOLDEN/inputs"
  local src=$BENCH/HG00097_r10.mt.bam; need_index "$src"
  local mean; mean=$(samtools depth -a -r chrM "$src" | awk '{s+=$3;n++} END{printf "%.1f", s/n}')
  for d in 100 1000; do
    local out=$GOLDEN/inputs/r10_ds$d.bam
    [[ -f $out ]] && continue
    local frac; frac=$(awk -v d=$d -v m="$mean" 'BEGIN{f=d/m; if(f>=1)f=0.999; printf "%.6f", f}')
    samtools view -b -s "42${frac#0}" -o "$out" "$src"
    samtools index "$out"
    "$HIMITO" build -i "$out" -r "$REF" -o "$GOLDEN/inputs/r10_ds$d.gfa" > "$GOLDEN/inputs/r10_ds$d.build.log" 2>&1
    echo "subset r10_ds$d: mean depth of source $mean, fraction $frac"
  done
}

run_case() {      # case outdir
  local c=$1 od=$2; read -r dt gfa bam <<< "$(case_inputs "$c")"
  need_index "$bam"; mkdir -p "$od"
  ( cd "$od" && "$HIMITO" call -g "$gfa" -r "$REF" -s HG00097_ont -d "$dt" -o "$c.vcf" -b "$bam" \
      --strand-bias-threshold 0.05 > "$c.log" 2>&1 ) || { echo "FAIL: call crashed for $c (see $od/$c.log)"; return 1; }
}

compare() {       # a b label
  if cmp -s "$1" "$2"; then return 0; fi
  echo "DIFF: $3"
  if [[ $1 == *.vcf ]]; then
    diff <(grep -v '^##' "$1" | cut -f1-7,10) <(grep -v '^##' "$2" | cut -f1-7,10) | head -20 || true
  else
    echo "  (sizes $(wc -c < "$1") vs $(wc -c < "$2") bytes; rows $(wc -l < "$1") vs $(wc -l < "$2"))"
  fi
  return 1
}

full_runs() {     # outdir ; end-to-end quick-start through the canonical bench script
  local od=$1; mkdir -p "$od/output"
  for m in ont-r10 ont-denoised; do
    ( cd "$BENCH" && BAM="$FULL_BAM" OUT="$od/output" HIMITO="$HIMITO" bash "$ROOT/scripts/ont_recall/run_hapmap_bench.sh" "$m" "regress_$m" > "$od/regress_$m.bench.txt" 2>&1 )
    cp "$BENCH/regress_$m.vcf" "$od/regress_$m.vcf"
    cp "$od/output/regress_$m.eval/summary.txt" "$od/regress_$m.summary.txt"
  done
}

status=0
case $MODE in
  freeze)
    [[ -e "$GOLDEN/MANIFEST" ]] && { echo "golden set exists at $GOLDEN; move it away to re-freeze"; exit 2; }
    make_subsets
    for c in $CASES; do
      run_case "$c" "$GOLDEN/run1"; run_case "$c" "$GOLDEN/run2"
      for x in $OUTS; do compare "$GOLDEN/run1/$c.$x" "$GOLDEN/run2/$c.$x" "non-deterministic: $c.$x" || status=1; done
    done
    if [[ $FULL == --full ]]; then
      full_runs "$GOLDEN/full1"; full_runs "$GOLDEN/full2"
      for m in ont-r10 ont-denoised; do
        compare "$GOLDEN/full1/regress_$m.vcf" "$GOLDEN/full2/regress_$m.vcf" "non-deterministic: quick-start $m VCF" || status=1
      done
    fi
    [[ $status == 0 ]] || { echo "freeze aborted: outputs are not reproducible run to run; fix that first"; exit 1; }
    mv "$GOLDEN/run1" "$GOLDEN/calls"; rm -rf "$GOLDEN/run2"
    [[ -d "$GOLDEN/full1" ]] && { mv "$GOLDEN/full1" "$GOLDEN/full"; rm -rf "$GOLDEN/full2"; }
    { echo "commit $(git -C "$ROOT" rev-parse --verify -q HEAD || echo unknown)$(git -C "$ROOT" diff --quiet 2>/dev/null || echo ' (+uncommitted changes)')"
      echo "frozen $(date -u +%FT%TZ)"
      ( cd "$GOLDEN" && find calls full inputs -type f \( -name '*.vcf' -o -name '*.csv' -o -name '*.txt' -o -name '*.bam' -o -name '*.gfa' \) 2>/dev/null | sort | xargs md5sum 2>/dev/null || true )
    } > "$GOLDEN/MANIFEST"
    for c in $CASES; do printf "%-11s PASS records: %s\n" "$c" "$(grep -v '^#' "$GOLDEN/calls/$c.vcf" | awk -F'\t' '$7=="PASS"' | wc -l)"; done
    echo "golden set frozen at $GOLDEN"
    ;;
  check)
    [[ -e "$GOLDEN/MANIFEST" ]] || { echo "no golden set at $GOLDEN; run freeze at the starting commit first"; exit 2; }
    tmp=$(mktemp -d "${TMPDIR:-/tmp}/ont_regress.XXXX")
    for c in $CASES; do
      run_case "$c" "$tmp" || { status=1; continue; }
      ok=1; for x in $OUTS; do compare "$GOLDEN/calls/$c.$x" "$tmp/$c.$x" "$c.$x" || { ok=0; status=1; }; done
      [[ $ok == 1 ]] && echo "OK   $c (vcf, matrix.csv, raw_matrix.csv identical)"
    done
    if [[ $FULL == --full ]]; then
      [[ -d "$GOLDEN/full" ]] || { echo "no frozen --full set; run freeze --full"; exit 2; }
      full_runs "$tmp/full"
      for m in ont-r10 ont-denoised; do
        compare "$GOLDEN/full/regress_$m.vcf" "$tmp/full/regress_$m.vcf" "quick-start $m VCF" || status=1
        compare "$GOLDEN/full/regress_$m.summary.txt" "$tmp/full/regress_$m.summary.txt" "quick-start $m vcfeval summary" || status=1
        [[ $status == 0 ]] && echo "OK   quick-start $m"
      done
    fi
    [[ $status == 0 ]] && { echo "ONT regression: PASS"; rm -rf "$tmp"; } || echo "ONT regression: FAIL (outputs kept in $tmp)"
    exit $status
    ;;
  *) echo "usage: ont_regression.sh freeze|check [--full]"; exit 2 ;;
esac
