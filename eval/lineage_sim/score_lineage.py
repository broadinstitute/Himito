#!/usr/bin/env python3
"""Score a reconstructed Himito mutation tree against the simulated truth.

Both the truth and the reconstruction use Himito's mutation_tree.tsv schema
(node_id, variant, parent_id, parent_variant). Nodes are keyed by their
`variant` string; ROOT is the sentinel "ROOT".
"""
import argparse
import csv
import os


def parse_mutation_tree(path: str) -> dict[str, str]:
    """variant -> parent_variant. ROOT maps to ROOT."""
    parent = {}
    with open(path) as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for row in rdr:
            parent[row["variant"]] = row["parent_variant"]
    return parent


def ancestor_sets(parent_map: dict[str, str]) -> dict[str, set[str]]:
    """variant -> set of STRICT ancestors (excludes self and ROOT)."""
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


def _pc_edges(parent_map: dict[str, str], keep: set[str]) -> set[tuple[str, str]]:
    """Direct (parent, child) edges restricted to `keep` (excludes ROOT edges)."""
    edges = set()
    for child, par in parent_map.items():
        if child == "ROOT" or par == "ROOT":
            continue
        if child in keep and par in keep:
            edges.add((par, child))
    return edges


def detected_variants_from_vcf(path: str) -> set[str]:
    """PASS SNVs from a Himito VCF, as m.<pos><ref>><alt> strings."""
    out = set()
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
            out.add(f"m.{pos}{ref}>{alt}")
    return out


def _f1(p: float, r: float) -> float:
    return 0.0 if (p + r) == 0 else 2 * p * r / (p + r)


def score(truth_parent, recon_parent, truth_vars, detected_vars) -> dict:
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
    tp_pairs = _anc_pairs(truth_anc, shared)
    rp_pairs = _anc_pairs(recon_anc, shared)
    inter = tp_pairs & rp_pairs
    union_pairs = tp_pairs | rp_pairs
    ad_recall = len(inter) / len(union_pairs) if union_pairs else 0.0
    ad_precision = len(inter) / len(union_pairs) if union_pairs else 0.0

    truth_edges = _pc_edges(truth_parent, shared)
    recon_edges = _pc_edges(recon_parent, shared)
    pc_recall = (len(truth_edges & recon_edges) / len(truth_edges)) if truth_edges else 0.0

    return {
        "n_truth_vars": len(truth_vars),
        "n_detected_vars": len(detected_vars),
        "n_shared": len(shared),
        "var_precision": var_precision,
        "var_recall": var_recall,
        "var_f1": _f1(var_precision, var_recall),
        "ad_precision": ad_precision,
        "ad_recall": ad_recall,
        "ad_f1": _f1(ad_precision, ad_recall),
        "pc_recall": pc_recall,
    }


FIELDS = ["profile", "fp", "fn", "n_truth_vars", "n_detected_vars", "n_shared",
          "var_precision", "var_recall", "var_f1",
          "ad_precision", "ad_recall", "ad_f1", "pc_recall"]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth-tree", required=True)
    ap.add_argument("--recon-tree", required=True)
    ap.add_argument("--truth-variants", required=True)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--profile", default="NA")
    ap.add_argument("--fp", default="NA")
    ap.add_argument("--fn", default="NA")
    ap.add_argument("--metrics-out", default="")
    args = ap.parse_args()

    truth_pm = parse_mutation_tree(args.truth_tree)
    recon_pm = parse_mutation_tree(args.recon_tree)
    truth_vars = {l.strip() for l in open(args.truth_variants) if l.strip()}
    detected = detected_variants_from_vcf(args.vcf)

    m = score(truth_pm, recon_pm, truth_vars, detected)
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
    print("\t".join(FIELDS))
    print(line)


if __name__ == "__main__":
    main()
