#!/usr/bin/env python3
"""Score a reconstructed Himito mutation tree against the simulated truth.

Both the truth and the reconstruction use Himito's mutation_tree.tsv schema
(node_id, variant, parent_id, parent_variant). Nodes are keyed by their
`variant` string; ROOT is the sentinel "ROOT".
"""
import argparse
import csv
import gzip
import json
import math
from collections import Counter
import os
import re
import sys

from hp_tracts import canonical_name
from simulate_tree import parse_background


def parse_mutation_tree(path: str) -> dict[str, str]:
    """variant -> parent_variant. ROOT maps to ROOT."""
    with open(path) as fh:
        return {row["variant"]: row["parent_variant"]
                for row in csv.DictReader(fh, delimiter="\t")
                if row.get("node_id") not in {"ROOT_BLOCK", "UNASSIGNED", "TRACT_ANNOTATED"}}


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


def detected_variants_with_hf_from_vcf(
    path: str, *, pass_only: bool = False,
) -> dict[str, float | None]:
    """Return accepted alleles and the first sample's per-ALT HF (or AF).

    By default both PASS and unfiltered (.) records are accepted.
    Missing or invalid estimates do not remove alleles from detection scoring.
    """
    out: dict[str, float | None] = {}
    accepted_filters = {"PASS"} if pass_only else {"PASS", "."}
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 7 or f[6] not in accepted_filters:
                continue
            pos, ref = f[1], f[3]
            sample = dict(zip(f[8].split(":"), f[9].split(":"))) if len(f) >= 10 else {}
            for i, alt in enumerate(f[4].split(",")):
                if alt in ("", "."):
                    continue
                hf = None
                for field in ("HF", "AF"):
                    values = sample.get(field, "").split(",")
                    if i < len(values):
                        try:
                            value = float(values[i])
                        except ValueError:
                            continue
                        if math.isfinite(value) and 0.0 <= value <= 1.0:
                            hf = value
                            break
                out[f"m.{pos}{ref}>{alt}"] = hf
    return out


def parse_clone_haplotypes(path: str) -> set[frozenset[str]]:
    """Truth clones from simulate_tree.py's clones.tsv, as variant sets.

    The mutation-free `ref` clone becomes the empty set; `score_haplotypes`
    drops it. Frequencies are scored separately by `score_frequencies`.
    """
    haps = set()
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            haps.add(frozenset(v for v in row["variant_path"].split(",") if v))
    return haps


def _canonical_variant(name: str, background: list[str]) -> str:
    return canonical_name(name, background)


def _canonical_parent_map(parent: dict[str, str], background: list[str]) -> dict[str, str]:
    return {
        _canonical_variant(key, background): _canonical_variant(value, background)
        for key, value in parent.items()
    }


def parse_matrix_haplotypes(path: str, keep: set[str] | None = None,
                            background: list[str] | None = None) -> set[frozenset[str]]:
    """Reconstructed haplotypes from Himito's lineage matrix, as variant sets.

    The matrix is variants (rows) x reads (columns); reads sharing a genotype
    pattern are one haplotype, which is exactly the grouping Himito reports in
    <prefix>.cleaned_haplotype_map.tsv. Only '1' counts as present: the pre-clean
    matrix writes an empty cell for "read does not span this site", and reading
    that as a mutation would invent haplotypes.

    ``keep``, when given, ignores rows whose variant is not in the set.
    """
    with open(path) as fh:
        rows = list(csv.reader(fh))
    if not rows:
        return set()
    body = rows[1:]
    # Plan 2026-09-27 T6. Rename before the keep filter so an equivalent
    # spelling of a mutation-tree variant is not dropped.
    if background is not None:
        body = [
            [_canonical_variant(r[0], background), *r[1:]] if r else r
            for r in body
        ]
    if keep is not None:
        body = [r for r in body if r and r[0] in keep]
    variants = [r[0] for r in body]
    haps = set()
    for i in range(len(rows[0]) - 1):
        haps.add(frozenset(
            v for v, row in zip(variants, body)
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

    This exact-match metric is unweighted. Frequency agreement is reported
    separately, using read assignments as an estimate rather than known abundance.
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


def truth_haplotype_frequencies(path, background):
    """Terminal clone mass aggregated by canonical lineage genotype."""
    masses = Counter()
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            mass = float(row["frequency"])
            if not math.isfinite(mass) or mass < 0:
                raise ValueError(f"Invalid clone frequency in {path}: {mass}")
            hap = frozenset(_canonical_variant(v, background)
                            for v in row["variant_path"].split(",") if v)
            if mass:
                masses[hap] += mass
    total = math.fsum(masses.values())
    if total <= 0:
        raise ValueError(f"No positive clone frequencies in {path}")
    # Simulator TSVs round each clone mass independently.
    return {h: m / total for h, m in masses.items()}


def _frequency_matrix_rows(path, keep, background):
    with open(path) as fh:
        reader = csv.reader(fh)
        reads = next(reader)[1:]
        rows = {}
        for row in reader:
            name = row[0]
            if background is not None or name not in keep:
                name = _canonical_variant(name, background)
            if name in keep:
                if len(row) != len(reads) + 1:
                    raise ValueError(f"Wrong row length for {name} in {path}")
                rows.setdefault(name, []).append(row[1:])
    return reads, rows


def _allele_call(copies, i):
    """Combine canonical duplicates without counting a molecule twice."""
    values = [row[i] for row in copies]
    if any(v.isascii() and v.isdigit() and int(v) > 0 for v in values):
        return 1
    return 0 if all(v == "0" for v in values) else None


def _frequency_errors(prefix, errors):
    n = len(errors)
    return {
        f"{prefix}_n_compared": n,
        f"{prefix}_mae": math.fsum(map(abs, errors)) / n if n else None,
        f"{prefix}_rmse": math.sqrt(math.fsum(e * e for e in errors) / n) if n else None,
        f"{prefix}_bias": math.fsum(errors) / n if n else None,
    }


def score_frequencies(truth_freq, truth_vars, detected, detected_hf, cleaned_path,
                      keep, background, input_path=""):
    """Caller and matrix VAF, and attachment-frequency agreement."""
    truth_vaf = Counter()
    for hap, freq in truth_freq.items():
        for variant in hap:
            truth_vaf[variant] += freq
    reads, rows = _frequency_matrix_rows(cleaned_path, keep, background)
    if set(rows) != keep:
        raise ValueError(f"Tree variants missing from {cleaned_path}: {sorted(keep - rows.keys())}")

    clean_ids = {r for r in reads if r != "chrM"}
    if len(clean_ids) != sum(r != "chrM" for r in reads):
        raise ValueError(f"Duplicate read IDs in {cleaned_path}")
    denominator = clean_ids
    matrix_vaf = {}
    if input_path:
        raw_reads, raw_rows = _frequency_matrix_rows(input_path, truth_vars | detected, background)
        physical_indices = [i for i, r in enumerate(raw_reads) if r != "chrM"]
        input_reads = [raw_reads[i] for i in physical_indices]
        denominator = set(input_reads)
        if len(denominator) != len(input_reads) or not clean_ids <= denominator:
            raise ValueError("Input-matrix read IDs must be unique and include every cleaned read")
        for variant, copies in raw_rows.items():
            support = Counter(_allele_call(copies, i) for i in physical_indices)
            n_called = support[0] + support[1]
            matrix_vaf[variant] = (support[1] / n_called if n_called else None, n_called)

    variant_rows = []
    vcf_errors = []
    matrix_errors = []
    for variant in sorted(truth_vars | detected):
        expected = truth_vaf[variant]
        estimate = detected_hf.get(variant)
        matrix_estimate, n_called = matrix_vaf.get(variant, (None, 0))
        error = None if estimate is None else estimate - expected
        matrix_error = None if matrix_estimate is None else matrix_estimate - expected
        if variant in truth_vars:
            if error is not None:
                vcf_errors.append(error)
            if matrix_error is not None:
                matrix_errors.append(matrix_error)
        variant_rows.append({
            "variant": variant, "in_truth": variant in truth_vars,
            "called": variant in detected, "truth_vaf": expected,
            "vcf_vaf": estimate, "error": error,
            "matrix_vaf": matrix_estimate, "matrix_n_called": n_called,
            "matrix_error": matrix_error,
        })

    counts = Counter()
    for i, read in enumerate(reads):
        if read == "chrM":
            continue
        calls = {v: _allele_call(copies, i) for v, copies in rows.items()}
        if None not in calls.values():
            counts[frozenset(v for v, call in calls.items() if call == 1)] += 1
    n_reads = len(denominator)
    n_unassigned = n_reads - sum(counts.values())
    hap_rows = []
    for hap in sorted(set(truth_freq) | set(counts), key=lambda h: (len(h), sorted(h))):
        expected = truth_freq.get(hap, 0.0)
        estimate = counts[hap] / n_reads if n_reads else None
        hap_rows.append({
            "variants": sorted(hap), "truth_frequency": expected,
            "estimated_read_fraction": estimate, "n_reads": counts[hap],
            "error": None if estimate is None else estimate - expected,
        })
    unassigned = n_unassigned / n_reads if n_reads else None
    summary = {
        "vcf_vaf_n_truth": len(truth_vars),
        **_frequency_errors("vcf_vaf", vcf_errors),
        **(_frequency_errors("matrix_vaf", matrix_errors) if input_path else
           dict.fromkeys(_frequency_errors("matrix_vaf", []))),
        "hap_freq_n_reads": n_reads,
        "hap_freq_basis": "input_matrix" if input_path else "cleaned_matrix",
        "hap_freq_unassigned": unassigned,
        "hap_freq_mae": math.fsum(abs(row["error"]) for row in hap_rows) / len(hap_rows)
                        if n_reads and hap_rows else None,
        "hap_freq_tv": (math.fsum(abs(row["error"]) for row in hap_rows) + unassigned) / 2
                       if n_reads else None,
    }
    return summary, {
        "summary": summary, "variant_vaf": variant_rows, "haplotype_frequency": hap_rows,
        "n_unassigned_reads": n_unassigned,
    }


def _f1(p: float, r: float) -> float:
    return 0.0 if (p + r) == 0 else 2 * p * r / (p + r)


def branch_ancestry(reference, reconstruction):
    """Compare directed relationships after setting aside within-segment order."""
    def segments(tree):
        child_count = Counter(p for v, p in tree.items() if v != "ROOT")
        result = {}
        for variant in tree.keys() - {"ROOT"}:
            start = variant
            for _ in tree:
                parent = tree.get(start, "ROOT")
                if parent == "ROOT" or child_count[parent] != 1:
                    break
                start = parent
            else:
                raise ValueError(f"Cycle in tree at {variant}")
            result[variant] = start
        return result

    shared = (reference.keys() & reconstruction.keys()) - {"ROOT"}
    ref_pairs = _anc_pairs(ancestor_sets(reference), shared)
    recon_pairs = _anc_pairs(ancestor_sets(reconstruction), shared)
    ref_segments, recon_segments = segments(reference), segments(reconstruction)
    both_comparable = {frozenset(p) for p in ref_pairs} & {frozenset(p) for p in recon_pairs}
    skipped = {frozenset((a, b)) for a, b in both_comparable
               if ref_segments[a] == ref_segments[b] or recon_segments[a] == recon_segments[b]}
    ref_pairs = {p for p in ref_pairs if frozenset(p) not in skipped}
    recon_pairs = {p for p in recon_pairs if frozenset(p) not in skipped}
    matches = len(ref_pairs & recon_pairs)
    informative = bool(ref_pairs or recon_pairs)
    precision = matches / len(recon_pairs) if recon_pairs else 0.0
    recall = matches / len(ref_pairs) if ref_pairs else 0.0
    return dict(branch_ad_precision=precision if informative else None,
                branch_ad_recall=recall if informative else None,
                branch_ad_f1=_f1(precision, recall) if informative else None,
                branch_ad_n_ref_pairs=len(ref_pairs), branch_ad_n_query_pairs=len(recon_pairs),
                branch_ad_n_matching_pairs=matches, branch_ad_n_skipped_pairs=len(skipped),
                branch_ad_n_unshared_sites=len((reference.keys() ^ reconstruction.keys()) - {"ROOT"}),
                branch_ad_status="evaluated" if informative else "no_informative_pairs")


def tree_haplotypes(tree):
    return {frozenset(anc | {variant}) for variant, anc in ancestor_sets(tree).items()}


def matrix_haplotype_counts(path, sites, background):
    reads, rows = _frequency_matrix_rows(path, sites, background)
    if set(rows) != sites:
        raise ValueError(f"Tree variants missing from {path}: {sorted(sites - rows.keys())}")
    counts = Counter()
    for i, read in enumerate(reads):
        if read != "chrM":
            calls = {v: _allele_call(copies, i) for v, copies in rows.items()}
            if None not in calls.values():
                counts[frozenset(v for v, call in calls.items() if call == 1)] += 1
    physical_reads = [read for read in reads if read != "chrM"]
    if len(set(physical_reads)) != len(physical_reads):
        raise ValueError(f"Duplicate read IDs in {path}")
    return counts, set(physical_reads)


def haplotype_support(path, candidates, alternatives, background, min_support, min_eligible, *,
                      required_reads=()):
    """Raw evidence per candidate, using unique compatible observations only."""
    genotypes = candidates | alternatives | {frozenset()}
    panel = sorted(set().union(*genotypes))
    bits = {v: 1 << i for i, v in enumerate(panel)}
    masks = {h: sum(bits[v] for v in h) for h in genotypes}
    reads, rows = _frequency_matrix_rows(path, set(panel), background)
    indices = [i for i, read in enumerate(reads) if read != "chrM"]
    raw_ids = {reads[i] for i in indices}
    if len(raw_ids) != len(indices):
        raise ValueError(f"Duplicate read IDs in {path}")
    if not raw_ids.issuperset(required_reads):
        raise ValueError(f"{path} must include every cleaned read ID")
    patterns = Counter()
    for i in indices:
        observed = alternate = 0
        for v, copies in rows.items():
            call = _allele_call(copies, i)
            if call is not None:
                observed |= bits[v]
                if call:
                    alternate |= bits[v]
        patterns[observed, alternate] += 1

    evidence = {}
    for hap in candidates:
        expected = masks[hap]
        differences = [expected ^ mask for other, mask in masks.items() if other != hap]
        eligible = support = ambiguous = conflicts = 0
        for (observed, alternate), n in patterns.items():
            distinguishable = bool(observed) and all(observed & d for d in differences)
            eligible += n if distinguishable else 0
            if (observed & expected) != alternate:
                conflicts += n
            elif distinguishable:
                support += n
            else:
                ambiguous += n
        status = "supported"
        if eligible < min_eligible:
            status = "insufficient_coverage"
        elif support < min_support:
            status = "insufficient_support"
        evidence[hap] = dict(eligible=eligible, direct_support=support, ambiguous=ambiguous,
                             conflicting=conflicts, status=status)
    return evidence, len(indices)


def supported_calls(counts, own_tree_haps, evidence, scope="intermediate"):
    calls = counts.keys() - {frozenset()}
    intermediates = {h for h in own_tree_haps if any(h < other for other in own_tree_haps)}
    retained = {h for h in calls if evidence[h]["status"] == "supported"
                or (scope == "intermediate" and h not in intermediates)}
    return retained, intermediates


def support_detail(evidence, counts, retained, intermediates):
    return [dict(variants=sorted(hap), **evidence[hap], imputed_reads=counts[hap],
                 candidate_class=("not_reconstructed" if hap not in counts else
                                  "intermediate" if hap in intermediates else "terminal"),
                 reconstructed=hap in counts, retained=hap in retained)
            for hap in sorted(evidence, key=lambda h: sorted(h))]


BRANCH_FIELDS = list(branch_ancestry({}, {}))
SUPPORT_FIELDS = ["supported_hap_precision", "supported_hap_recall", "supported_hap_f1",
                  "supported_hap_tp", "supported_hap_fp", "supported_hap_fn",
                  "supported_hap_n_filtered", "supported_hap_filtered_read_fraction",
                  "supported_hap_freq_unassigned", "supported_hap_freq_tv"]


def supported_haplotype_metrics(reference, query):
    scores = score_haplotypes(reference, query)
    return dict(supported_hap_precision=scores["hap_precision"],
                supported_hap_recall=scores["hap_recall"], supported_hap_f1=scores["hap_f1"],
                supported_hap_tp=len(reference & query), supported_hap_fp=len(query - reference),
                supported_hap_fn=len(reference - query))


def retained_haplotype_frequencies(counts, retained, n_reads):
    if not n_reads:
        return {}, None
    frequencies = {h: count / n_reads for h, count in counts.items() if not h or h in retained}
    return frequencies, (n_reads - sum(counts[h] for h in frequencies)) / n_reads


def haplotype_frequency_tv(reference, query, ref_unassigned, query_unassigned):
    if ref_unassigned is None or query_unassigned is None:
        return None
    distance = math.fsum(abs(reference.get(h, 0.0) - query.get(h, 0.0))
                         for h in reference.keys() | query.keys())
    return (distance + abs(ref_unassigned - query_unassigned)) / 2


def add_support_args(parser):
    parser.add_argument("--min-clone-support", type=int, default=3,
                        help="minimum direct reads to retain an intermediate clone (default: 3)")
    parser.add_argument("--min-clone-eligible", type=int, default=10,
                        help="minimum eligible raw reads to evaluate a clone (default: 10)")
    parser.add_argument("--clone-support-scope", choices=("intermediate", "all"), default="intermediate",
                        help="which reconstructed clones must pass the support rule")
    parser.add_argument("--support-out", help="JSON with per-clone raw support and imputed counts")


def score(truth_parent, recon_parent, truth_vars, detected_vars,
          truth_haps=None, recon_haps=None, background=None) -> dict:
    # Plan 2026-09-27 T6. None leaves names as written. A list, including an
    # empty one, renames alleles whose REF span sits in an effective tract.
    if background is not None:
        background = list(background)
        truth_parent = _canonical_parent_map(truth_parent, background)
        recon_parent = _canonical_parent_map(recon_parent, background)
        truth_vars = {_canonical_variant(v, background) for v in truth_vars}
        detected_vars = {_canonical_variant(v, background) for v in detected_vars}
        if truth_haps is not None:
            truth_haps = {
                frozenset(_canonical_variant(v, background) for v in hap)
                for hap in truth_haps
            }
        if recon_haps is not None:
            recon_haps = {
                frozenset(_canonical_variant(v, background) for v in hap)
                for hap in recon_haps
            }

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

FREQUENCY_FIELDS = [
    "vcf_vaf_n_truth", "vcf_vaf_n_compared", "vcf_vaf_mae", "vcf_vaf_rmse", "vcf_vaf_bias",
    "matrix_vaf_n_compared", "matrix_vaf_mae", "matrix_vaf_rmse", "matrix_vaf_bias",
    "hap_freq_n_reads", "hap_freq_basis", "hap_freq_unassigned", "hap_freq_mae", "hap_freq_tv",
]

# Three metric families, each with its own n_* denominators: variant detection,
# ancestor-descendant tree accuracy, and clone (haplotype) recovery. Consumers
# look columns up by header name, not by position.
# Plan 2026-09-27 T6. Appended after n_background_calls. Do not reorder earlier columns.
HP_FIELDS = ["n_hp_unstable_truth", "n_hp_unstable_called",
             "hp_var_precision", "hp_var_recall", "hp_var_f1",
             "hp_intrusion", "hp_backbone_nodes", "hp_len_l1", "root_fraction"]

FIELDS = ["profile", "fp", "fn", "n_truth_vars", "n_detected_vars", "n_shared",
          "n_truth_pairs",
          "var_precision", "var_recall", "var_f1",
          "ad_precision", "ad_recall", "ad_f1"] + HAP_FIELDS + ["n_background_calls"] + HP_FIELDS + FREQUENCY_FIELDS + BRANCH_FIELDS + SUPPORT_FIELDS

_SNV_NAME = re.compile(r"^m\.\d+([ACGT])>([ACGT])$")


def load_truth_background(path: str) -> set[str]:
    """Names in ``truth/background.txt``. Missing or empty means exclude nothing.

    Plan 2026-09-27 T5. These are fixed homoplasmies. ``Himito lineage`` already
    drops them as ``fixed``; the VCF can still report them at HF near 1.
    """
    if not path or not os.path.isfile(path):
        return set()
    names = set()
    with open(path) as fh:
        for line in fh:
            name = line.strip()
            if name:
                names.add(name)
    return names


def _is_snv(name: str) -> bool:
    return _SNV_NAME.fullmatch(name) is not None


def load_hp_tract_rows(path: str, background: list[str]) -> list[dict]:
    """Non-zero δ rows of ``truth/hp_tracts.tsv``, names canonicalised.

    Plan 2026-09-27 T6. A missing file is an empty list. δ = 0 rows are
    not truth alleles.
    """
    if not path or not os.path.isfile(path):
        return []
    rows = []
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if int(row["delta"]) == 0:
                continue
            rows.append({
                "tract": row["tract"],
                "truth_freq": float(row["truth_freq"]),
                "canonical_name": _canonical_variant(
                    row["canonical_name"].strip(), background,
                ),
            })
    return rows


def load_hp_lineage_names(path: str, background: list[str]) -> set[str] | None:
    """Variant column of ``truth/hp_lineage.tsv``, or None if the file is missing."""
    if not path or not os.path.isfile(path):
        return None
    names = set()
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            name = row.get("variant", "").strip()
            if name:
                names.add(_canonical_variant(name, background))
    return names


def _hp_detection(hp_names: set[str] | None, detected: set[str]) -> dict:
    """Precision and recall restricted to the HP-lineage names."""
    if hp_names is None:
        return {field: "NA" for field in ("hp_var_precision", "hp_var_recall", "hp_var_f1")}
    called = detected & hp_names
    # The file is the truth list, so every call of one of its names is a true positive.
    precision = 1.0 if called else 0.0
    recall = len(called) / len(hp_names) if hp_names else 0.0
    return {
        "hp_var_precision": precision,
        "hp_var_recall": recall,
        "hp_var_f1": _f1(precision, recall),
    }


def hp_intrusion(truth_parent: dict[str, str], recon_parent: dict[str, str],
                 unstable: set[str]) -> float:
    """Fraction of truth SNVs in the recon tree whose ancestors include an unstable allele.

    Plan 2026-09-27 T6. 0.0 when no such SNV is in the reconstruction.
    """
    recon_nodes = {v for v in recon_parent if v != "ROOT"}
    snvs = [
        v for v in truth_parent
        if v != "ROOT" and _is_snv(v) and v in recon_nodes
    ]
    if not snvs:
        return 0.0
    ancestors = ancestor_sets(recon_parent)
    return sum(1 for v in snvs if ancestors.get(v, set()) & unstable) / len(snvs)


def hp_backbone_nodes(recon_parent: dict[str, str], unstable: set[str]) -> int:
    """Unstable names that are the parent of at least one recon node."""
    parents = set(recon_parent.values())
    return sum(1 for name in unstable if name in parents)


def root_fraction(path: str):
    """ROOT ``n_reads_attached`` over the per-node sum, or ``NA``."""
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames or "n_reads_attached" not in reader.fieldnames:
            return "NA"
        total = 0
        root = None
        for row in reader:
            raw = (row.get("n_reads_attached") or "").strip()
            if raw in ("", "NA"):
                continue
            n = int(raw)
            total += n
            if row.get("variant") == "ROOT":
                root = n
        if root is None or total == 0:
            return "NA"
        return root / total


def _matrix_presence(path: str, background: list[str]) -> tuple[dict[str, set[int]], int]:
    with open(path) as fh:
        rows = list(csv.reader(fh))
    if not rows:
        return {}, 0
    n_reads = max(0, len(rows[0]) - 1)
    presence: dict[str, set[int]] = {}
    for row in rows[1:]:
        if not row:
            continue
        name = _canonical_variant(row[0], background)
        hits = presence.setdefault(name, set())
        for i in range(n_reads):
            if len(row) > i + 1 and row[i + 1] == "1":
                hits.add(i)
    return presence, n_reads


def hp_length_l1(rows: list[dict], matrix_path: str, background: list[str]):
    """Mean over tracts of Σ_δ |truth_freq(δ) − recon_freq(δ)|.

    Plan 2026-09-27 T6. Reads carrying none of a tract's canonical alleles are
    δ = 0. ``NA`` when there is no non-zero δ row or no matrix to count.
    """
    if not rows or not matrix_path or not os.path.isfile(matrix_path):
        return "NA"
    presence, n_reads = _matrix_presence(matrix_path, background)
    if n_reads == 0:
        return "NA"
    by_tract: dict[str, dict[str, float]] = {}
    for row in rows:
        alleles = by_tract.setdefault(row["tract"], {})
        name = row["canonical_name"]
        alleles[name] = alleles.get(name, 0.0) + row["truth_freq"]
    total = 0.0
    for alleles in by_tract.values():
        truth0 = max(0.0, 1.0 - sum(alleles.values()))
        counts = {name: 0 for name in alleles}
        zeros = 0
        for read in range(n_reads):
            hit = [name for name in alleles if read in presence.get(name, ())]
            if not hit:
                zeros += 1
            else:
                for name in hit:
                    counts[name] += 1
        l1 = abs(truth0 - zeros / n_reads)
        for name, freq in alleles.items():
            l1 += abs(freq - counts[name] / n_reads)
        total += l1
    return total / len(by_tract)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth-tree", required=True)
    ap.add_argument("--recon-tree", required=True)
    ap.add_argument("--truth-variants", required=True)
    ap.add_argument("--vcf", required=True)
    # Optional pair: Himito's per-read genotype matrix and the simulator's clone
    # table. Supply both to get the hap_* columns; omit both and they read NA.
    ap.add_argument("--recon-matrix", default="",
                    help="Himito <prefix>.cleaned_matrix_all_variants.csv "
                         "(reconstructed haplotypes; rows outside the mutation "
                         "tree are ignored)")
    ap.add_argument("--truth-clones", default="",
                    help="simulate_tree.py truth/clones.tsv (truth clone genomes)")
    ap.add_argument("--profile", default="NA")
    ap.add_argument("--fp", default="NA")
    ap.add_argument("--fn", default="NA")
    ap.add_argument("--metrics-out", default="")
    ap.add_argument("--input-matrix", default="", help="Pre-lineage BAM genotype matrix; counts omitted reads as unassigned")
    ap.add_argument("--frequency-out", default="", help="JSON detail for truth/called VAF and haplotype read fractions")
    add_support_args(ap)
    ap.add_argument(
        "--truth-background",
        default="",
        help="fixed-homoplasmy names removed from detected before var_precision "
        "(plan 2026-09-27 T5). Missing or empty leaves detected unchanged",
    )
    ap.add_argument(
        "--background",
        default="none",
        help="fixed homoplasmies that reshape tract anchors before names are "
        "canonicalised (plan 2026-09-27 T6): none, hvs1-polyC, or "
        "comma-separated <REF><pos><ALT> tokens",
    )
    ap.add_argument(
        "--truth-hp-tracts",
        default="",
        help="truth/hp_tracts.tsv. canonical_name values, after canonicalisation, "
        "are unstable alleles excluded from the var_precision denominator "
        "(plan 2026-09-27 T6)",
    )
    ap.add_argument(
        "--truth-hp-lineage",
        default="",
        help="truth/hp_lineage.tsv. hp_var_* are computed on its variant column "
        "(plan 2026-09-27 T6)",
    )
    args = ap.parse_args()
    if min(args.min_clone_support, args.min_clone_eligible) < 1:
        ap.error("clone thresholds must be positive")
    if args.support_out and not args.input_matrix:
        ap.error("--support-out requires --input-matrix")
    try:
        background_tokens = parse_background(args.background)
    except ValueError as exc:
        ap.error(str(exc))

    # Plan 2026-09-27 T6. Both sides use the same effective tracts, so a
    # collapsed reconstruction matches an uncollapsed truth.
    truth_pm = _canonical_parent_map(parse_mutation_tree(args.truth_tree), background_tokens)
    recon_pm = _canonical_parent_map(parse_mutation_tree(args.recon_tree), background_tokens)
    truth_vars = {
        _canonical_variant(l.strip(), background_tokens)
        for l in open(args.truth_variants) if l.strip()
    }
    detected_hf: dict[str, float | None] = {}
    for name, hf in detected_variants_with_hf_from_vcf(args.vcf, pass_only=True).items():
        detected_hf.setdefault(_canonical_variant(name, background_tokens), hf)
    detected = set(detected_hf)
    vcf_names = set(detected)
    # Plan 2026-09-27 T5. Count hits before dropping them from the precision
    # denominator. An empty set leaves `detected` unchanged. Names are
    # canonicalised with the same background so a tract-resident homoplasmy
    # still matches the VCF.
    background = {
        _canonical_variant(name, background_tokens)
        for name in load_truth_background(args.truth_background)
    }
    n_background_calls = len(detected & background)
    # Plan 2026-09-27 T6. Unstable length alleles are real, so they are not
    # false positives and they are not added to truth_vars.
    tract_rows = load_hp_tract_rows(args.truth_hp_tracts, background_tokens)
    unstable = {row["canonical_name"] for row in tract_rows}
    n_hp_unstable_truth = len(unstable)
    n_hp_unstable_called = len(unstable & vcf_names)
    # A lineage allele can share its canonical name with an hp_tracts row.
    # Those calls are true positives, so they stay in the precision denominator.
    detected -= (background - truth_vars)
    detected -= (unstable - truth_vars)
    hp_det = _hp_detection(
        load_hp_lineage_names(args.truth_hp_lineage, background_tokens), vcf_names,
    )

    if bool(args.recon_matrix) != bool(args.truth_clones):
        ap.error("--recon-matrix and --truth-clones must be given together")
    truth_haps = recon_haps = None
    if args.recon_matrix:
        truth_haps = parse_clone_haplotypes(args.truth_clones)
        keep = set(recon_pm) - {"ROOT"}
        recon_haps = parse_matrix_haplotypes(
            args.recon_matrix, keep=keep, background=background_tokens,
        )

    truth_snvs = {name for name in truth_vars if _is_snv(name)}
    detected_snvs = {name for name in detected if _is_snv(name)}
    m = score(truth_pm, recon_pm, truth_snvs, detected_snvs,
              truth_haps=truth_haps, recon_haps=recon_haps,
              background=background_tokens)
    row = {"profile": args.profile, "fp": args.fp, "fn": args.fn, **m,
           "n_background_calls": n_background_calls,
           "n_hp_unstable_truth": n_hp_unstable_truth,
           "n_hp_unstable_called": n_hp_unstable_called,
           **hp_det,
           "hp_intrusion": hp_intrusion(truth_pm, recon_pm, unstable),
           "hp_backbone_nodes": hp_backbone_nodes(recon_pm, unstable),
           "hp_len_l1": hp_length_l1(tract_rows, args.recon_matrix, background_tokens),
           "root_fraction": root_fraction(args.recon_tree)}

    row.update(branch_ancestry(truth_pm, recon_pm))
    row.update(dict.fromkeys(SUPPORT_FIELDS))
    row.update(dict.fromkeys(FREQUENCY_FIELDS))
    if args.input_matrix or args.frequency_out:
        if not args.recon_matrix:
            ap.error("--input-matrix and --frequency-out require --recon-matrix and --truth-clones")
    if args.truth_clones:
        try:
            truth_freq = truth_haplotype_frequencies(args.truth_clones, background_tokens)
            summary, details = score_frequencies(
                truth_freq, truth_vars, detected, detected_hf, args.recon_matrix,
                keep, background_tokens, args.input_matrix,
            )
        except (ValueError, KeyError) as exc:
            ap.error(str(exc))
        row.update(summary)
        if args.input_matrix:
            counts, _ = matrix_haplotype_counts(args.recon_matrix, keep, background_tokens)
            candidates = counts.keys() - {frozenset()}
            tree_haps = tree_haplotypes(recon_pm)
            evidence, n_reads = haplotype_support(args.input_matrix, candidates, tree_haps,
                                                  background_tokens, args.min_clone_support,
                                                  args.min_clone_eligible)
            retained, intermediates = supported_calls(counts, tree_haps, evidence, args.clone_support_scope)
            expected = set(truth_freq) - {frozenset()}
            removed = candidates - retained
            frequencies, unassigned = retained_haplotype_frequencies(counts, retained, n_reads)
            row.update(supported_haplotype_metrics(expected, retained),
                       supported_hap_n_filtered=len(removed),
                       supported_hap_filtered_read_fraction=sum(counts[h] for h in removed) / n_reads if n_reads else None,
                       supported_hap_freq_unassigned=unassigned,
                       supported_hap_freq_tv=haplotype_frequency_tv(truth_freq, frequencies, 0.0, unassigned))
            hap_details = support_detail(evidence, counts, retained, intermediates)
            for record in hap_details:
                hap = frozenset(record["variants"])
                record["in_truth"] = hap in expected
                record["truth_frequency"] = truth_freq.get(hap, 0.0)
            details["haplotype_support"] = hap_details
            details["supported_hap_missing_truth"] = [
                dict(variants=sorted(hap), truth_frequency=truth_freq[hap],
                     reason=evidence[hap]["status"] if hap in candidates else "not_reconstructed")
                for hap in sorted(expected - retained, key=lambda h: sorted(h))
            ]
            truth_tree_haps = tree_haplotypes(truth_pm)
            truth_intermediates = {h for h in expected if any(h < other for other in truth_tree_haps)}
            intermediate_calls = (candidates & truth_intermediates) | ((candidates - expected) & intermediates)
            details["haplotype_support_by_clone_class"] = {
                "intermediate": supported_haplotype_metrics(truth_intermediates, retained & intermediate_calls),
                "terminal": supported_haplotype_metrics(expected - truth_intermediates, retained - intermediate_calls),
            }
            details["summary"].update({field: row[field] for field in SUPPORT_FIELDS})
            if args.support_out:
                with open(args.support_out, "w") as fh:
                    json.dump(hap_details, fh, indent=2)
                    fh.write("\n")
        if args.frequency_out:
            with open(args.frequency_out, "w") as fh:
                json.dump(details, fh, indent=2, allow_nan=False)
                fh.write("\n")

    def fmt(x):
        return "NA" if x is None else f"{x:.4f}" if isinstance(x, float) else str(x)

    line = "\t".join(fmt(row[f]) for f in FIELDS)
    if args.metrics_out:
        new = not os.path.exists(args.metrics_out) or os.path.getsize(args.metrics_out) == 0
        if not new:
            with open(args.metrics_out) as existing:
                if existing.readline().rstrip("\r\n").split("\t") != FIELDS:
                    ap.error("--metrics-out has a different schema; choose a new output file")
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

    # Report false positives from the same PASS SNV set used for variant precision.
    fps = sorted(detected_snvs - truth_snvs)
    print(f"false_positives\t{len(fps)}", file=sys.stderr)
    print("variant\tHF", file=sys.stderr)
    for vid in fps:
        hf = detected_hf[vid]
        print(f"{vid}\t{hf if hf is not None else 'NA'}", file=sys.stderr)


if __name__ == "__main__":
    main()
