#!/usr/bin/env python3
"""Score a reconstructed Himito mutation tree against the simulated truth.

Both the truth and the reconstruction use Himito's mutation_tree.tsv schema
(node_id, variant, parent_id, parent_variant). Nodes are keyed by their
`variant` string; ROOT is the sentinel "ROOT".
"""
import argparse
import csv
import os
import sys


def parse_mutation_tree(path: str) -> dict[str, str]:
    """variant -> parent_variant. ROOT maps to ROOT."""
    parent = {}
    with open(path) as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for row in rdr:
            parent[row["variant"]] = row["parent_variant"]
    return parent


def ancestor_sets(parent_map: dict[str, str]) -> dict[str, set[str]]:
    """variant -> set of proper variant ancestors (excludes self and ROOT).

    ROOT is deliberately left out: it is an ancestor of every variant in every
    tree, so it carries no information about a reconstruction. See `score` for
    what including it did to `ad_precision`.
    """
    anc = {}
    for v in parent_map:
        if v == "ROOT":
            continue
        acc, cur, steps = set(), parent_map.get(v, "ROOT"), 0
        while cur != "ROOT":
            acc.add(cur)
            cur = parent_map.get(cur, "ROOT")
            steps += 1
            if steps > len(parent_map):  # guard against malformed cycles
                break
        anc[v] = acc
    return anc


def _anc_pairs(anc: dict[str, set[str]], keep: set[str]) -> set[tuple[str, str]]:
    """Ordered (ancestor, descendant) pairs restricted to `keep`."""
    pairs = set()
    for desc, ancestors in anc.items():
        if desc not in keep:
            continue
        for a in ancestors:
            if a in keep:
                pairs.add((a, desc))
    return pairs


def detected_variants_with_hf_from_vcf(path: str) -> dict[str, float | None]:
    """PASS/. calls from a Himito VCF → {m.<pos><ref>><alt>: HF}.

    The FILTER column is the only gate: every record marked PASS or "." counts as a
    detected variant. HF is read from the first sample's FORMAT/HF field (first allele
    if multi-valued) and is reporting-only -- it never excludes a record, and a
    missing/unparseable HF yields None rather than dropping the call.

    No HF band and no SNV-only gate are applied here on purpose. Both used to shrink
    the ``var_precision`` denominator to the slice ``Himito lineage`` consumes, which
    made precision read 1.0 on runs whose VCF was in fact full of false positives:
    error-driven calls sat below the band's floor and near-homoplasmic artifacts above
    its ceiling, so neither reached the metric. Indels are likewise counted now.
    ``simulate_tree.py`` emits substitutions only, so every indel call is a false
    positive by construction -- expect ONT ``var_precision`` to read low whenever
    the permutation test lets indel artifacts through. ``run_himito.sh`` leaves
    ``call -p`` at the data-type default (does not pass ``-p 1``). That is a
    property of the caller output, not of lineage reconstruction; read
    ``ad_*`` for tree accuracy.
    """
    out: dict[str, float | None] = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 7:
                continue
            pos, ref, alt, filt = f[1], f[3], f[4], f[6]
            if filt not in ("PASS", "."):
                continue
            hf: float | None = None
            if len(f) >= 10:
                keys = f[8].split(":")
                vals = f[9].split(":")
                sample = dict(zip(keys, vals))
                if "HF" in sample and sample["HF"] not in (".", ""):
                    hf = float(sample["HF"].split(",")[0])
            out[f"m.{pos}{ref}>{alt}"] = hf
    return out


def parse_clone_haplotypes(path: str) -> set[frozenset[str]]:
    """Truth clones from simulate_tree.py's clones.tsv, as variant sets.

    The mutation-free `ref` clone becomes the empty set; `score_haplotypes`
    drops it. Frequencies are deliberately ignored -- see `score_haplotypes`
    for why read fractions are not a usable abundance estimate here.
    """
    haps = set()
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            haps.add(frozenset(v for v in row["variant_path"].split(",") if v))
    return haps


def parse_matrix_haplotypes(path: str) -> set[frozenset[str]]:
    """Reconstructed haplotypes from Himito's lineage matrix, as variant sets.

    The matrix is variants (rows) x reads (columns); reads sharing a genotype
    pattern are one haplotype, which is exactly the grouping Himito reports in
    <prefix>.cleaned_haplotype_map.tsv. Only '1' counts as present: the pre-clean
    matrix writes an empty cell for "read does not span this site", and reading
    that as a mutation would invent haplotypes.
    """
    with open(path) as fh:
        rows = list(csv.reader(fh))
    if not rows:
        return set()
    variants = [r[0] for r in rows[1:]]
    haps = set()
    for i in range(len(rows[0]) - 1):
        haps.add(frozenset(
            v for v, row in zip(variants, rows[1:])
            if len(row) > i + 1 and row[i + 1] == "1"
        ))
    return haps


def score_haplotypes(
    truth_haps: set[frozenset[str]],
    recon_haps: set[frozenset[str]],
) -> dict:
    """Exact variant-set agreement between truth clones and reconstructed haplotypes.

    A truth clone and a reconstructed haplotype match when their full variant sets
    are equal, so a false-positive or missed call breaks the match. That is what a
    downstream consumer of the haplotype map actually gets.

    The empty (mutation-free) haplotype is excluded from both sides. It is not a
    clone: a read lands there whenever it fails to span any variant site, so on
    seed5_mut10_depth1000 that bin holds 1226 reads of which only 590 are genuinely
    reference -- the other 636 come from mutant clones. Scoring it as a recovered
    reference clone would add a free true positive to every run.

    Haplotypes are counted, not read-weighted, and no abundance metric is reported:
    a haplotype's read count is not proportional to its clone frequency, so read
    fractions cannot be compared against the truth frequencies in clones.tsv.
    Per-site heteroplasmy is available from the VCF's HF field instead.
    """
    truth = {h for h in truth_haps if h}
    recon = {h for h in recon_haps if h}

    matched = truth & recon
    precision = len(matched) / len(recon) if recon else 0.0
    recall = len(matched) / len(truth) if truth else 0.0

    return {
        "n_truth_clones": len(truth),
        "n_recon_haps": len(recon),
        "hap_precision": precision,
        "hap_recall": recall,
        "hap_f1": _f1(precision, recall),
    }


def _f1(p: float, r: float) -> float:
    return 0.0 if (p + r) == 0 else 2 * p * r / (p + r)


def score(truth_parent, recon_parent, truth_vars, detected_vars,
          truth_haps=None, recon_haps=None) -> dict:
    # --- variant detection (tree-independent) ---
    tp = len(truth_vars & detected_vars)
    var_precision = tp / len(detected_vars) if detected_vars else 0.0
    var_recall = tp / len(truth_vars) if truth_vars else 0.0

    # --- shared variant set for tree metrics ---
    truth_tree_vars = {v for v in truth_parent if v != "ROOT"}
    recon_tree_vars = {v for v in recon_parent if v != "ROOT"}
    shared = truth_tree_vars & recon_tree_vars

    truth_anc = ancestor_sets(truth_parent)
    recon_anc = ancestor_sets(recon_parent)
    # ROOT is excluded from the pair sets. It is an ancestor of every variant in
    # both trees, so keeping it adds |shared| pairs that match by construction and
    # drags every score toward 1: on a 4-mutation chain truth, a completely flat
    # reconstruction scored ad_precision 1.00 / ad_f1 0.57 with ROOT in the sets,
    # and 0.00 with it out. ROOT was originally kept so that two identical star
    # topologies would not score 0 on empty pair sets; the both-empty branch below
    # handles that case directly, without inflating everything else.
    #
    # ad_recall is measured against EVERY truth ancestral pair, not only pairs
    # among `shared`. Restricting the denominator to `shared` meant a run that
    # detected fewer variants scored HIGHER, because the pairs it had no chance of
    # recovering simply left the denominator: over 100 runs corr(var_recall,
    # ad_f1) ran -0.40 to -0.63 within a config, so improving detection lowered the
    # tree score. The extreme case was --n-mutations 15 with denoise --vaf
    # destroying 40% of the truth set, which scored a perfect ad_f1 = 1.000 on all
    # 10 seeds; against all truth pairs it scores 0.374. With this denominator the
    # same correlations run +0.27 to +0.98.
    #
    # ad_precision stays restricted to `shared`. A pair built from a
    # false-positive variant is a *detection* error and var_precision already
    # reports it; charging it here too would double-count one mistake.
    tp_all = _anc_pairs(truth_anc, truth_tree_vars)
    tp_shared = _anc_pairs(truth_anc, shared)
    rp_pairs = _anc_pairs(recon_anc, shared)
    inter = tp_shared & rp_pairs
    if not tp_all and not rp_pairs:
        # Truth has no ancestral pairs at all and neither does the reconstruction:
        # they agree, there is simply nothing to count. Distinct from "one side is
        # flat", which scores 0 below.
        ad_precision = ad_recall = 1.0
    else:
        ad_recall = len(inter) / len(tp_all) if tp_all else 0.0
        ad_precision = len(inter) / len(rp_pairs) if rp_pairs else 0.0

    # --- haplotype (clone) recovery ---
    if truth_haps is None or recon_haps is None:
        hap = {f: "NA" for f in HAP_FIELDS}
    else:
        hap = score_haplotypes(truth_haps, recon_haps)

    return {
        "n_truth_vars": len(truth_vars),
        "n_detected_vars": len(detected_vars),
        "n_shared": len(shared),
        # Denominator of ad_recall. With recall scored against all truth pairs,
        # n_shared alone no longer says what ad_* was measured over.
        "n_truth_pairs": len(tp_all),
        "var_precision": var_precision,
        "var_recall": var_recall,
        "var_f1": _f1(var_precision, var_recall),
        "ad_precision": ad_precision,
        "ad_recall": ad_recall,
        "ad_f1": _f1(ad_precision, ad_recall),
        **hap,
    }


# Haplotype columns, kept separate so they can be filled with "NA" in one place
# when --recon-matrix/--truth-clones are not supplied.
HAP_FIELDS = ["n_truth_clones", "n_recon_haps",
              "hap_precision", "hap_recall", "hap_f1"]

# Three metric families, each with its own n_* denominators: variant detection,
# ancestor-descendant tree accuracy, and clone (haplotype) recovery. Consumers
# look columns up by header name, not by position -- see sweep_fpfn.sh.
FIELDS = ["profile", "fp", "fn", "n_truth_vars", "n_detected_vars", "n_shared",
          "n_truth_pairs",
          "var_precision", "var_recall", "var_f1",
          "ad_precision", "ad_recall", "ad_f1"] + HAP_FIELDS


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth-tree", required=True)
    ap.add_argument("--recon-tree", required=True)
    ap.add_argument("--truth-variants", required=True)
    ap.add_argument("--vcf", required=True)
    # Optional pair: Himito's per-read genotype matrix and the simulator's clone
    # table. Supply both to get the hap_* columns; omit both and they read NA.
    ap.add_argument("--recon-matrix", default="",
                    help="Himito <prefix>.cleaned_matrix.csv (reconstructed haplotypes)")
    ap.add_argument("--truth-clones", default="",
                    help="simulate_tree.py truth/clones.tsv (truth clone genomes)")
    ap.add_argument("--profile", default="NA")
    ap.add_argument("--fp", default="NA")
    ap.add_argument("--fn", default="NA")
    ap.add_argument("--metrics-out", default="")
    args = ap.parse_args()

    truth_pm = parse_mutation_tree(args.truth_tree)
    recon_pm = parse_mutation_tree(args.recon_tree)
    truth_vars = {l.strip() for l in open(args.truth_variants) if l.strip()}
    detected_hf = detected_variants_with_hf_from_vcf(args.vcf)
    detected = set(detected_hf)

    if bool(args.recon_matrix) != bool(args.truth_clones):
        ap.error("--recon-matrix and --truth-clones must be given together")
    truth_haps = recon_haps = None
    if args.recon_matrix:
        truth_haps = parse_clone_haplotypes(args.truth_clones)
        recon_haps = parse_matrix_haplotypes(args.recon_matrix)

    m = score(truth_pm, recon_pm, truth_vars, detected,
              truth_haps=truth_haps, recon_haps=recon_haps)
    row = {"profile": args.profile, "fp": args.fp, "fn": args.fn, **m}

    def fmt(x):
        return f"{x:.4f}" if isinstance(x, float) else str(x)

    line = "\t".join(fmt(row[f]) for f in FIELDS)
    if args.metrics_out:
        new = not os.path.exists(args.metrics_out) or os.path.getsize(args.metrics_out) == 0
        with open(args.metrics_out, "a") as fh:
            if new:
                fh.write("\t".join(FIELDS) + "\n")
            fh.write(line + "\n")
    else:
        # Standalone use has no file to read the result back from, so print it.
        # When --metrics-out IS given, callers (e.g. run_eval.sh) read that file
        # back themselves to display it -- printing here too would just show the
        # same single result twice in the terminal.
        print("\t".join(FIELDS))
        print(line)

    # False positives: PASS/. VCF calls not in the truth variant list.
    fps = sorted(detected - truth_vars)
    print(f"false_positives\t{len(fps)}", file=sys.stderr)
    print("variant\tHF", file=sys.stderr)
    for vid in fps:
        hf = detected_hf[vid]
        print(f"{vid}\t{hf if hf is not None else 'NA'}", file=sys.stderr)


if __name__ == "__main__":
    main()
