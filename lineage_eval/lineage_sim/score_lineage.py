#!/usr/bin/env python3
"""Score a reconstructed Himito mutation tree against the simulated truth.

Both the truth and the reconstruction use Himito's mutation_tree.tsv schema
(node_id, variant, parent_id, parent_variant). Nodes are keyed by their
`variant` string; ROOT is the sentinel "ROOT".
"""
import argparse
import csv
import os
import re
import sys

from hp_tracts import canonical_name
from simulate_tree import parse_background


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


def _canonical_variant(name: str, background: list[str]) -> str:
    """Canonicalise one name, or each member of a comma-joined ROOT_BLOCK."""
    if "," not in name:
        return canonical_name(name, background)
    return ",".join(
        canonical_name(part, background) if part else part
        for part in name.split(",")
    )


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

    Plan 2026-09-26 T3. ``keep``, when given, ignores rows whose variant is not
    in the set. ``cleaned_matrix_all_variants.csv`` appends near-fixed
    homoplasmies that were not rows of the old ``cleaned_matrix.csv``; passing
    the mutation-tree variant names leaves haplotype scores unchanged.
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

# Three metric families, each with its own n_* denominators: variant detection,
# ancestor-descendant tree accuracy, and clone (haplotype) recovery. Consumers
# look columns up by header name, not by position -- see sweep_fpfn.sh.
# Plan 2026-09-27 T6. Appended after n_background_calls. Do not reorder earlier columns.
HP_FIELDS = ["n_hp_unstable_truth", "n_hp_unstable_called",
             "hp_var_precision", "hp_var_recall", "hp_var_f1",
             "hp_intrusion", "hp_backbone_nodes", "hp_len_l1", "root_fraction"]

FIELDS = ["profile", "fp", "fn", "n_truth_vars", "n_detected_vars", "n_shared",
          "n_truth_pairs",
          "var_precision", "var_recall", "var_f1",
          "ad_precision", "ad_recall", "ad_f1"] + HAP_FIELDS + ["n_background_calls"] + HP_FIELDS

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
    """ROOT ``n_reads_attached`` over the per-node sum, or ``NA``.

    Plan 2026-09-27 T6. ``write_mutation_tree`` repeats the root count on a
    ``ROOT_BLOCK`` row. That row is an annotation, not another node, so it
    is left out of the denominator. The numerator is still the ROOT row.
    """
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames or "n_reads_attached" not in reader.fieldnames:
            return "NA"
        total = 0
        root = None
        for row in reader:
            if row.get("node_id") == "ROOT_BLOCK":
                continue
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
                         "tree, including appended near-fixed homoplasmies, are ignored)")
    ap.add_argument("--truth-clones", default="",
                    help="simulate_tree.py truth/clones.tsv (truth clone genomes)")
    ap.add_argument("--profile", default="NA")
    ap.add_argument("--fp", default="NA")
    ap.add_argument("--fn", default="NA")
    ap.add_argument("--metrics-out", default="")
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
    for name, hf in detected_variants_with_hf_from_vcf(args.vcf).items():
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
        # Plan 2026-09-26 T3. These names are the rows of the old
        # cleaned_matrix.csv: every mutation-tree key except ROOT, with a
        # comma-joined ROOT_BLOCK key split into its member variants.
        keep: set[str] = set()
        for key in recon_pm:
            if key == "ROOT":
                continue
            if "," in key:
                keep.update(part for part in key.split(",") if part)
            else:
                keep.add(key)
        recon_haps = parse_matrix_haplotypes(
            args.recon_matrix, keep=keep, background=background_tokens,
        )

    m = score(truth_pm, recon_pm, truth_vars, detected,
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
