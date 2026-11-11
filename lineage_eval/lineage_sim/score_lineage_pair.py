#!/usr/bin/env python3
"""Compare two Himito lineages on their original and imputed matrices."""

import argparse
from collections import Counter
import csv
import json
import sys

from score_lineage import (
    _allele_call,
    _frequency_matrix_rows,
    add_support_args,
    branch_ancestry,
    haplotype_frequency_tv,
    haplotype_support,
    matrix_haplotype_counts,
    parse_mutation_tree,
    retained_haplotype_frequencies,
    score,
    score_haplotypes,
    support_detail,
    supported_calls,
    supported_haplotype_metrics,
    tree_haplotypes,
)


def order_evidence(tree_path, raw_path, tree):
    edges = [(parent, child) for child, parent in tree.items() if child != "ROOT" and parent != "ROOT"]
    reads, rows = _frequency_matrix_rows(raw_path, {v for edge in edges for v in edge}, None)
    with open(tree_path) as handle:
        margins = {row["variant"]: row.get("order_support_ll") for row in csv.DictReader(handle, delimiter="\t")}
    result = []
    for parent, child in edges:
        counts = Counter((_allele_call(rows[parent], i), _allele_call(rows[child], i))
                         for i, read in enumerate(reads)
                         if read != "chrM" and parent in rows and child in rows)
        result.append(dict(parent=parent, child=child,
                           co_callable=sum(counts[a, b] for a in (0, 1) for b in (0, 1)),
                           parent_only=counts[1, 0], child_only=counts[0, 1], both_alt=counts[1, 1],
                           order_support_ll=margins.get(child)))
    return result


def compare_pair(reference, query, min_support, min_eligible, scope="intermediate"):
    # Each input is (mutation tree, pre-imputation matrix, cleaned matrix).
    inputs = {"reference": reference, "query": query}
    trees = {label: parse_mutation_tree(paths[0]) for label, paths in inputs.items()}
    genotypes = {label: tree_haplotypes(tree) for label, tree in trees.items()}
    cleaned = {
        label: matrix_haplotype_counts(paths[2], trees[label].keys() - {"ROOT"}, None)
        for label, paths in inputs.items()
    }
    counts = {label: result[0] for label, result in cleaned.items()}
    ref, rec = counts["reference"].keys() - {frozenset()}, counts["query"].keys() - {frozenset()}
    candidates = ref | rec
    alternatives = genotypes["reference"] | genotypes["query"]
    metrics = {k: v for k, v in score(trees["reference"], trees["query"], set(), set()).items()
               if k.startswith("ad_")}
    metrics.update(branch_ancestry(trees["reference"], trees["query"]))
    strict_haps = score_haplotypes(ref, rec)
    metrics.update(n_reference_haps=strict_haps.pop("n_truth_clones"),
                   n_query_haps=strict_haps.pop("n_recon_haps"), **strict_haps)
    evidence, retained, details, totals = {}, {}, {}, {}
    for label, paths in inputs.items():
        evidence[label], total = haplotype_support(paths[1], candidates, alternatives, None,
                                                   min_support, min_eligible, required_reads=cleaned[label][1])
        totals[label] = total
        retained[label], intermediates = supported_calls(counts[label], genotypes[label], evidence[label], scope)
        removed = (counts[label].keys() - {frozenset()}) - retained[label]
        metrics[f"{label}_filtered_read_fraction"] = (
            sum(counts[label][h] for h in removed) / total if total else None
        )
        details[label] = support_detail(evidence[label], counts[label], retained[label], intermediates)

    eligible = {h for h in candidates if all(evidence[label][h]["eligible"] >= min_eligible
                                            for label in inputs)}
    reference_calls = retained["reference"] & eligible
    query_calls = retained["query"] & eligible
    supported = supported_haplotype_metrics(reference_calls, query_calls)
    if not (reference_calls or query_calls):
        for field in ("precision", "recall", "f1"):
            supported[f"supported_hap_{field}"] = None
    metrics.update(supported, supported_hap_n_evaluable=len(eligible),
                   supported_hap_n_unresolved=len(candidates - eligible),
                   supported_hap_n_evaluable_called=len(reference_calls | query_calls))
    ref_freq, ref_missing = retained_haplotype_frequencies(counts["reference"], reference_calls, totals["reference"])
    query_freq, query_missing = retained_haplotype_frequencies(counts["query"], query_calls, totals["query"])
    metrics.update(reference_unassigned_read_fraction=ref_missing, query_unassigned_read_fraction=query_missing,
                   supported_hap_freq_tv=haplotype_frequency_tv(ref_freq, query_freq, ref_missing, query_missing))
    details["unresolved_haplotypes"] = [sorted(h) for h in sorted(candidates - eligible, key=lambda h: sorted(h))]
    details["order_evidence"] = {
        label: order_evidence(paths[0], paths[1], trees[label]) for label, paths in inputs.items()
    }
    return metrics, details


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for side in ("reference", "query"):
        parser.add_argument(f"--{side}", nargs=3, required=True,
                            metavar=("TREE_TSV", "RAW_CSV", "CLEAN_CSV"))
    add_support_args(parser)
    parser.add_argument("--metrics-out")
    args = parser.parse_args()
    if min(args.min_clone_support, args.min_clone_eligible) < 1:
        parser.error("clone thresholds must be positive")
    try:
        metrics, details = compare_pair(args.reference, args.query, args.min_clone_support,
                                        args.min_clone_eligible, args.clone_support_scope)
        if args.metrics_out:
            with open(args.metrics_out, "w", newline="") as handle:
                write_metrics(handle, metrics)
        else:
            write_metrics(sys.stdout, metrics)
        if args.support_out:
            with open(args.support_out, "w") as handle:
                json.dump(details, handle, indent=2)
                handle.write("\n")
    except (OSError, ValueError) as exc:
        parser.error(str(exc))


def write_metrics(handle, metrics):
    writer = csv.DictWriter(handle, fieldnames=list(metrics), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerow({k: "NA" if v is None else v for k, v in metrics.items()})


if __name__ == "__main__":
    main()
