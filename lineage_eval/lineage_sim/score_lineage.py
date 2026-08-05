#!/usr/bin/env python3
"""Score a reconstructed Himito mutation tree against the simulated truth.

Both the truth and the reconstruction use Himito's mutation_tree.tsv schema
(node_id, variant, parent_id, parent_variant). Nodes are keyed by their
`variant` string; ROOT is the sentinel "ROOT".
"""
import argparse
import csv
import itertools
import os
import sys
from collections import deque


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


def _tree_graph(parent_map: dict[str, str]) -> dict[str, set[str]]:
    """Undirected adjacency over all tree nodes (variants + ROOT).

    Edges are the parent<->child links; ROOT is kept as a connector node so
    distances between variants in different top-level lineages route through it.
    """
    adj: dict[str, set[str]] = {}
    for child, par in parent_map.items():
        if child == "ROOT":
            continue
        adj.setdefault(child, set()).add(par)
        adj.setdefault(par, set()).add(child)
    return adj


def _pairwise_distances(parent_map: dict[str, str], nodes: list[str]) -> dict[str, dict[str, int]]:
    """Topological (edge-count) distance from every node in `nodes` to all others.

    BFS over the full tree graph (unit edge weights), so paths correctly traverse
    non-shared intermediate mutations and ROOT.
    """
    adj = _tree_graph(parent_map)
    dist: dict[str, dict[str, int]] = {}
    for src in nodes:
        seen = {src: 0}
        queue = deque([src])
        while queue:
            u = queue.popleft()
            for w in adj.get(u, ()):
                if w not in seen:
                    seen[w] = seen[u] + 1
                    queue.append(w)
        dist[src] = seen
    return dist


def _clades(anc: dict[str, set[str]], keep: set[str]) -> set[frozenset[str]]:
    """Non-trivial rooted clades restricted to `keep`.

    Each node maps to its subtree membership (self + descendants); singletons and
    the universal set are dropped (they are shared by every tree on the same taxa).
    """
    desc: dict[str, set[str]] = {v: {v} for v in keep}
    for d, ancestors in anc.items():
        if d not in keep:
            continue
        for a in ancestors:
            if a in keep:
                desc[a].add(d)
    n = len(keep)
    return {frozenset(members) for members in desc.values() if 1 < len(members) < n}


def robinson_foulds(truth_anc, recon_anc, shared: set[str]) -> tuple[int, float]:
    """Rooted Robinson-Foulds distance on the shared variant set.

    Symmetric difference of the two trees' clade sets. Returns
    (raw_count, normalized in [0, 1]); 0 == identical topology on `shared`.
    """
    ct = _clades(truth_anc, shared)
    cr = _clades(recon_anc, shared)
    raw = len(ct ^ cr)
    denom = len(ct) + len(cr)
    return raw, (raw / denom if denom else 0.0)


def _quartet_topology(d: dict[str, dict[str, int]], q: tuple[str, str, str, str]):
    """Unrooted quartet topology of 4 taxa via the four-point condition.

    Returns the bipartition (frozenset of two frozensets) implied by the smallest
    of the three distance-sum pairings, or None when the smallest is tied (an
    unresolved / star quartet).
    """
    a, b, c, e = q
    sums = [
        (d[a][b] + d[c][e], frozenset({frozenset({a, b}), frozenset({c, e})})),
        (d[a][c] + d[b][e], frozenset({frozenset({a, c}), frozenset({b, e})})),
        (d[a][e] + d[b][c], frozenset({frozenset({a, e}), frozenset({b, c})})),
    ]
    lo = min(val for val, _ in sums)
    winners = [topo for val, topo in sums if val == lo]
    return winners[0] if len(winners) == 1 else None


def quartet_distance(truth_parent, recon_parent, shared: set[str]) -> tuple[int, float]:
    """Quartet distance on the shared variant set.

    Number of 4-taxon subsets whose induced (unrooted) topology differs between
    the trees; a resolved-vs-star mismatch counts as a difference. Returns
    (raw_count, normalized in [0, 1]); 0.0 when fewer than 4 shared variants.
    """
    verts = sorted(shared)
    if len(verts) < 4:
        return 0, 0.0
    td = _pairwise_distances(truth_parent, verts)
    rd = _pairwise_distances(recon_parent, verts)
    diff = total = 0
    for quad in itertools.combinations(verts, 4):
        total += 1
        if _quartet_topology(td, quad) != _quartet_topology(rd, quad):
            diff += 1
    return diff, (diff / total if total else 0.0)


_BASES = frozenset("ACGT")


def detected_variants_with_hf_from_vcf(
    path: str,
    min_hf: float = 0.01,
    max_hf: float = 0.95,
) -> dict[str, float | None]:
    """PASS/. SNVs from a Himito VCF → {m.<pos><ref>><alt>: HF}.

    HF is read from the first sample's FORMAT/HF field (first allele if
    multi-valued). Non-PASS filters are excluded. Keeps only variants with
    ``min_hf <= HF < max_hf`` (same band as Himito ``lineage``); records
    without a usable HF are dropped.

    Only single-base substitutions count. ``simulate_tree.py`` builds the truth set
    from substitutions alone, so an indel call can never be a true positive; letting
    indels through only deflated ``var_precision``. That mattered a lot on ONT graphs,
    where indel artifacts outnumber SNV calls by ~50x and are filtered solely by the
    permutation test -- with that test disabled (``call -p 1``) they dominated the
    denominator and made precision look like a caller regression.
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
            # SNVs only: rejects indels, multiallelic ALTs ("G,T") and non-ACGT bases.
            if ref not in _BASES or alt not in _BASES:
                continue
            hf: float | None = None
            if len(f) >= 10:
                keys = f[8].split(":")
                vals = f[9].split(":")
                sample = dict(zip(keys, vals))
                if "HF" in sample and sample["HF"] not in (".", ""):
                    hf = float(sample["HF"].split(",")[0])
            if hf is None or not (min_hf <= hf < max_hf):
                continue
            out[f"m.{pos}{ref}>{alt}"] = hf
    return out


def detected_variants_from_vcf(
    path: str,
    min_hf: float = 0.01,
    max_hf: float = 0.95,
) -> set[str]:
    """PASS/. SNVs from a Himito VCF, as m.<pos><ref>><alt> strings."""
    return set(detected_variants_with_hf_from_vcf(path, min_hf=min_hf, max_hf=max_hf))


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
    ad_recall = len(inter) / len(tp_pairs) if tp_pairs else 0.0
    ad_precision = len(inter) / len(rp_pairs) if rp_pairs else 0.0

    truth_edges = _pc_edges(truth_parent, shared)
    recon_edges = _pc_edges(recon_parent, shared)
    pc_recall = (len(truth_edges & recon_edges) / len(truth_edges)) if truth_edges else 0.0

    # --- topology distances on the shared variant set ---
    rf, rf_norm = robinson_foulds(truth_anc, recon_anc, shared)
    quartet_dist, quartet_norm = quartet_distance(truth_parent, recon_parent, shared)

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
        "rf": rf,
        "rf_norm": rf_norm,
        "quartet_dist": quartet_dist,
        "quartet_norm": quartet_norm,
    }


FIELDS = ["profile", "fp", "fn", "n_truth_vars", "n_detected_vars", "n_shared",
          "var_precision", "var_recall", "var_f1",
          "ad_precision", "ad_recall", "ad_f1", "pc_recall",
          "rf", "rf_norm", "quartet_dist", "quartet_norm"]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth-tree", required=True)
    ap.add_argument("--recon-tree", required=True)
    ap.add_argument("--truth-variants", required=True)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--profile", default="NA")
    ap.add_argument("--fp", default="NA")
    ap.add_argument("--fn", default="NA")
    ap.add_argument(
        "--min-hf", type=float, default=0.01,
        help="keep detected variants with HF >= this (default: 0.01)",
    )
    ap.add_argument(
        "--max-hf", type=float, default=0.95,
        help="keep detected variants with HF < this (default: 0.95)",
    )
    ap.add_argument("--metrics-out", default="")
    args = ap.parse_args()

    truth_pm = parse_mutation_tree(args.truth_tree)
    recon_pm = parse_mutation_tree(args.recon_tree)
    truth_vars = {l.strip() for l in open(args.truth_variants) if l.strip()}
    detected_hf = detected_variants_with_hf_from_vcf(
        args.vcf, min_hf=args.min_hf, max_hf=args.max_hf,
    )
    detected = set(detected_hf)

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

    # False positives: PASS/. VCF calls not in the truth variant list.
    fps = sorted(detected - truth_vars)
    print(f"false_positives\t{len(fps)}", file=sys.stderr)
    print("variant\tHF", file=sys.stderr)
    for vid in fps:
        hf = detected_hf[vid]
        print(f"{vid}\t{hf if hf is not None else 'NA'}", file=sys.stderr)


if __name__ == "__main__":
    main()
