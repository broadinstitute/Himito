# eval/lineage_sim/tests/test_score_lineage.py
import pathlib
import sys

import pytest

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import score_lineage as sl

TRUTH = str(HERE / "fixtures" / "truth_tree.tsv")
RECON = str(HERE / "fixtures" / "recon_tree.tsv")


def test_parse_mutation_tree_maps_variant_to_parent():
    pm = sl.parse_mutation_tree(TRUTH)
    assert pm["m.200C>T"] == "m.100A>G"
    assert pm["m.100A>G"] == "ROOT"


def test_ancestor_sets_are_transitive():
    pm = sl.parse_mutation_tree(TRUTH)
    anc = sl.ancestor_sets(pm)
    assert anc["m.400T>C"] == {"m.200C>T", "m.100A>G"}
    assert anc["m.300G>A"] == {"m.100A>G"}
    assert anc["m.100A>G"] == set()


def test_score_ad_and_pc_against_known_trees():
    truth_pm = sl.parse_mutation_tree(TRUTH)
    recon_pm = sl.parse_mutation_tree(RECON)
    truth_vars = {"m.100A>G", "m.200C>T", "m.300G>A", "m.400T>C"}
    detected = truth_vars  # assume perfect detection for this unit test
    m = sl.score(truth_pm, recon_pm, truth_vars, detected)

    assert m["var_precision"] == 1.0
    assert m["var_recall"] == 1.0

    # Truth ancestor pairs (ancestor, descendant):
    #   (A,C),(A,G),(A,T),(C,T) -> 4 truth-anc pairs
    # Recon ancestor pairs (G misplaced under C instead of A):
    #   (A,C),(A,T),(C,T),(A,G),(C,G) -> 5 recon-anc pairs
    # Intersection: {(A,C),(A,G),(A,T),(C,T)} = 4
    # ad_recall = 4/4 = 1.0  (denominator = |truth_pairs|)
    # ad_precision = 4/5 = 0.8  (denominator = |recon_pairs|)
    assert abs(m["ad_recall"] - 1.0) < 1e-9
    assert abs(m["ad_precision"] - 0.8) < 1e-9
    # Regression guard: asymmetric definitions must not collapse to Jaccard (recall != precision here).
    assert m["ad_recall"] != m["ad_precision"]

    # Truth parent->child edges: A100->C200, A100->G300, C200->T400 (3).
    # Recon has A100->C200 and C200->T400 (G300 misplaced). Recovered = 2/3.
    assert abs(m["pc_recall"] - 2/3) < 1e-9


def test_rf_and_quartet_against_known_trees():
    truth_pm = sl.parse_mutation_tree(TRUTH)
    recon_pm = sl.parse_mutation_tree(RECON)
    truth_vars = {"m.100A>G", "m.200C>T", "m.300G>A", "m.400T>C"}
    m = sl.score(truth_pm, recon_pm, truth_vars, truth_vars)

    # Truth clades (subtree sets, excl. singletons/universal): {C,T}.
    # Recon clades (G moved under C): {C,T,G}. Symmetric difference = 2 clades.
    assert m["rf"] == 2
    assert abs(m["rf_norm"] - 1.0) < 1e-9

    # Only one quartet {A,C,G,T}. Truth resolves A,G|C,T; in recon A, G, T all
    # hang off C (a star) -> unresolved -> the two topologies differ.
    assert m["quartet_dist"] == 1
    assert abs(m["quartet_norm"] - 1.0) < 1e-9


def test_rf_and_quartet_are_zero_for_identical_trees():
    truth_pm = sl.parse_mutation_tree(TRUTH)
    truth_vars = {"m.100A>G", "m.200C>T", "m.300G>A", "m.400T>C"}
    m = sl.score(truth_pm, truth_pm, truth_vars, truth_vars)
    assert m["rf"] == 0
    assert m["rf_norm"] == 0.0
    assert m["quartet_dist"] == 0
    assert m["quartet_norm"] == 0.0
    assert m["ue_f1"] == 1.0
    assert m["clade_f1"] == 1.0
    assert m["path_undirected_f1"] == 1.0


def test_two_node_path_reversal_keeps_undirected_scores_perfect():
    # Truth: ROOT -> A -> B. Recon: ROOT -> B -> A (same edge, flipped order).
    truth = {"ROOT": "ROOT", "m.A": "ROOT", "m.B": "m.A"}
    recon = {"ROOT": "ROOT", "m.B": "ROOT", "m.A": "m.B"}
    vars_ = {"m.A", "m.B"}
    m = sl.score(truth, recon, vars_, vars_)

    assert m["ue_precision"] == 1.0
    assert m["ue_recall"] == 1.0
    assert m["ue_f1"] == 1.0
    assert m["path_undirected_precision"] == 1.0
    assert m["path_undirected_recall"] == 1.0
    assert m["path_undirected_f1"] == 1.0
    # n=2 drops the universal clade, so no non-trivial clades → F1 is 0/0 → 0.
    assert m["clade_f1"] == 0.0
    # Directed metrics must still see the order error.
    assert m["ad_f1"] < 1.0
    assert m["pc_recall"] == 0.0


def test_clade_f1_perfect_when_rooted_clades_match_despite_label_order():
    # Identical trees with non-trivial clades → clade F1 = 1.
    pm = {
        "ROOT": "ROOT",
        "m.A": "ROOT",
        "m.B": "m.A",
        "m.C": "m.A",
        "m.D": "m.B",
    }
    vars_ = {"m.A", "m.B", "m.C", "m.D"}
    m = sl.score(pm, pm, vars_, vars_)
    assert m["clade_precision"] == 1.0
    assert m["clade_recall"] == 1.0
    assert m["clade_f1"] == 1.0


def test_quartet_zero_when_fewer_than_four_shared_variants():
    # Only 3 shared variants -> no quartet can be formed.
    pm = {"ROOT": "ROOT", "m.100A>G": "ROOT", "m.200C>T": "m.100A>G", "m.300G>A": "m.100A>G"}
    shared_vars = {"m.100A>G", "m.200C>T", "m.300G>A"}
    m = sl.score(pm, pm, shared_vars, shared_vars)
    assert m["quartet_dist"] == 0
    assert m["quartet_norm"] == 0.0


def test_detected_variants_from_vcf_reads_pass_snvs(tmp_path):
    vcf = tmp_path / "x.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSIM\n"
        "chrM\t100\t.\tA\tG\t.\tPASS\t.\tGT:AD:HF\t1:3:0.3\n"
        "chrM\t250\t.\tC\tT\t.\tLowQual\t.\tGT:AD:HF\t1:2:0.2\n"
    )
    got = sl.detected_variants_from_vcf(str(vcf))
    assert got == {"m.100A>G"}


def test_detected_variants_with_hf_parses_format_hf(tmp_path):
    vcf = tmp_path / "x.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSIM\n"
        "chrM\t100\t.\tA\tG\t.\tPASS\tDP=10\tGT:AD:HF\t1:3:0.3\n"
        "chrM\t200\t.\tC\tT\t.\tPotential_Artifact\tDP=10\tGT:AD:HF\t1:2:0.2\n"
        "chrM\t300\t.\tG\tA\t.\tPASS\tDP=10\tGT:AD:HF\t1:5:0.5\n"
    )
    got = sl.detected_variants_with_hf_from_vcf(str(vcf))
    assert got == {"m.100A>G": 0.3, "m.300G>A": 0.5}


def test_detected_variants_includes_indels_and_multiallelic(tmp_path):
    # FILTER is the ONLY gate: every PASS/. record counts toward the precision
    # denominator, whatever its ALT looks like. simulate_tree.py emits substitutions
    # only, so indel/multiallelic calls are all false positives by construction --
    # that is deliberate, the metric is meant to expose them rather than exempt them.
    vcf = tmp_path / "x.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSIM\n"
        "chrM\t100\t.\tA\tG\t.\tPASS\t.\tGT:AD:HF\t1:3:0.3\n"      # SNV
        "chrM\t200\t.\tT\tTC\t.\tPASS\t.\tGT:AD:HF\t1:3:0.3\n"     # insertion
        "chrM\t300\t.\tCT\tC\t.\tPASS\t.\tGT:AD:HF\t1:3:0.3\n"     # deletion
        "chrM\t400\t.\tA\tG,T\t.\tPASS\t.\tGT:AD:HF\t1:3:0.3\n"    # multiallelic
        "chrM\t500\t.\tA\tN\t.\tPASS\t.\tGT:AD:HF\t1:3:0.3\n"      # non-ACGT
    )
    got = sl.detected_variants_with_hf_from_vcf(str(vcf))
    assert got == {
        "m.100A>G": 0.3,
        "m.200T>TC": 0.3,
        "m.300CT>C": 0.3,
        "m.400A>G,T": 0.3,
        "m.500A>N": 0.3,
    }


def test_detected_variants_applies_no_hf_band(tmp_path):
    # Regression guard for the removed [min_hf, max_hf) band: low-HF error calls and
    # near-homoplasmic calls alike must reach the precision denominator. Previously
    # both ends were silently dropped, which is what made var_precision read 1.0 on
    # runs whose VCF was full of low-HF false positives.
    vcf = tmp_path / "x.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSIM\n"
        "chrM\t100\t.\tA\tG\t.\tPASS\t.\tGT:AD:HF\t1:1:0.005\n"
        "chrM\t200\t.\tC\tT\t.\tPASS\t.\tGT:AD:HF\t1:1:0.01\n"
        "chrM\t300\t.\tG\tA\t.\tPASS\t.\tGT:AD:HF\t1:1:0.94\n"
        "chrM\t400\t.\tT\tC\t.\tPASS\t.\tGT:AD:HF\t1:1:0.99\n"
    )
    got = sl.detected_variants_with_hf_from_vcf(str(vcf))
    assert got == {
        "m.100A>G": 0.005,
        "m.200C>T": 0.01,
        "m.300G>A": 0.94,
        "m.400T>C": 0.99,
    }


def test_detected_variants_keeps_records_without_usable_hf(tmp_path):
    # A missing/"." HF used to drop the record entirely (it could not satisfy the
    # band). With no band, HF is reporting-only: the call still counts, HF is None.
    vcf = tmp_path / "x.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSIM\n"
        "chrM\t100\t.\tA\tG\t.\tPASS\t.\tGT:AD:HF\t1:3:.\n"        # HF is "."
        "chrM\t200\t.\tC\tT\t.\tPASS\t.\tGT:AD\t1:3\n"             # no HF key
        "chrM\t300\t.\tG\tA\t.\tPASS\t.\n"                          # no sample column
    )
    got = sl.detected_variants_with_hf_from_vcf(str(vcf))
    assert got == {"m.100A>G": None, "m.200C>T": None, "m.300G>A": None}


def _write_matrix(tmp_path, header, rows):
    p = tmp_path / "cleaned_matrix.csv"
    p.write_text("\n".join([",".join(header)] + [",".join(r) for r in rows]) + "\n")
    return str(p)


def _write_clones(tmp_path, paths_and_freqs):
    p = tmp_path / "clones.tsv"
    lines = ["clone_id\tnode_variant\tvariant_path\tfrequency\tn_reads"]
    for i, (vpath, freq) in enumerate(paths_and_freqs):
        node = vpath.split(",")[-1] if vpath else "ROOT"
        lines.append(f"{i}\t{node}\t{vpath}\t{freq}\t0")
    p.write_text("\n".join(lines) + "\n")
    return str(p)


def test_parse_matrix_haplotypes_groups_reads_by_pattern(tmp_path):
    path = _write_matrix(
        tmp_path,
        ["variant", "r1", "r2", "r3", "r4"],
        [["m.100A>G", "1", "1", "0", "1"],
         ["m.200C>T", "0", "1", "0", "0"]],
    )
    assert sl.parse_matrix_haplotypes(path) == {
        frozenset({"m.100A>G"}): 2,
        frozenset({"m.100A>G", "m.200C>T"}): 1,
        frozenset(): 1,
    }


def test_all_ones_row_outside_keep_does_not_change_score_haplotypes(tmp_path):
    # Plan 2026-09-26 T3. cleaned_matrix_all_variants.csv appends near-fixed
    # homoplasmies as all-1 rows. Those variants are outside the mutation-tree
    # keep set, so haplotype scores match the matrix without that row.
    header = ["variant", "r1", "r2", "r3"]
    base_rows = [
        ["m.100A>G", "1", "1", "0"],
        ["m.200C>T", "0", "1", "0"],
    ]
    extra_rows = base_rows + [["m.fixedA>G", "1", "1", "1"]]

    def dump(name, rows):
        path = tmp_path / name
        path.write_text(
            "\n".join([",".join(header)] + [",".join(r) for r in rows]) + "\n"
        )
        return str(path)

    base = dump("base.csv", base_rows)
    extra = dump("all_variants.csv", extra_rows)
    keep = {"m.100A>G", "m.200C>T"}
    truth = {frozenset({"m.100A>G"}), frozenset({"m.100A>G", "m.200C>T"})}
    without = sl.parse_matrix_haplotypes(base)
    filtered = sl.parse_matrix_haplotypes(extra, keep=keep)
    assert sl.score_haplotypes(truth, without) == sl.score_haplotypes(truth, filtered)
    assert sl.parse_matrix_haplotypes(extra) != without


def test_parse_matrix_haplotypes_treats_missing_cells_as_absent(tmp_path):
    # Himito's pre-clean matrix uses an empty cell for "read does not span this
    # site". Only '1' means present; missing must not be read as a mutation.
    path = _write_matrix(
        tmp_path,
        ["variant", "r1", "r2"],
        [["m.100A>G", "1", ""],
         ["m.200C>T", "", "0"]],
    )
    assert sl.parse_matrix_haplotypes(path) == {
        frozenset({"m.100A>G"}): 1,
        frozenset(): 1,
    }


def test_parse_clone_haplotypes_keeps_reference_clone_as_empty_set(tmp_path):
    path = _write_clones(tmp_path, [("", 0.4), ("m.100A>G", 0.2),
                                    ("m.100A>G,m.200C>T", 0.4)])
    assert sl.parse_clone_haplotypes(path) == {
        frozenset(),
        frozenset({"m.100A>G"}),
        frozenset({"m.100A>G", "m.200C>T"}),
    }


def test_haplotype_scores_are_perfect_when_every_clone_is_recovered():
    truth = {frozenset(), frozenset({"m.100A>G"}),
             frozenset({"m.100A>G", "m.200C>T"})}
    recon = {frozenset(): 5, frozenset({"m.100A>G"}): 2,
             frozenset({"m.100A>G", "m.200C>T"}): 1}
    m = sl.score_haplotypes(truth, recon, {"m.100A>G", "m.200C>T"})
    assert m["hap_precision"] == 1.0
    assert m["hap_recall"] == 1.0
    assert m["hap_f1"] == 1.0
    assert m["hap_read_frac"] == 1.0
    # The mutation-free bin is excluded from both sides, so 3 rows -> 2 clones.
    assert m["n_truth_clones"] == 2
    assert m["n_recon_haps"] == 2


def test_strict_penalizes_false_positive_variant_but_shared_does_not():
    # The reconstruction groups reads correctly; one haplotype just carries an
    # extra variant the caller invented. Strict scoring must see a miss, the
    # shared-variant restriction must not.
    truth = {frozenset({"m.100A>G"}), frozenset({"m.100A>G", "m.200C>T"})}
    recon = {frozenset({"m.100A>G", "m.999T>A"}): 3,
             frozenset({"m.100A>G", "m.200C>T"}): 1}
    m = sl.score_haplotypes(truth, recon, {"m.100A>G", "m.200C>T"})
    assert m["hap_precision"] == 0.5
    assert m["hap_recall"] == 0.5
    assert m["hap_s_precision"] == 1.0
    assert m["hap_s_recall"] == 1.0


def test_shared_restriction_merges_haplotypes_that_collapse_together():
    # Two reconstructed haplotypes differing only by a false-positive variant
    # become one key under restriction; they must merge rather than count twice
    # (which would drag precision to 0.5 for a single correct haplotype).
    truth = {frozenset({"m.100A>G"})}
    recon = {frozenset({"m.100A>G"}): 4,
             frozenset({"m.100A>G", "m.999T>A"}): 1}
    m = sl.score_haplotypes(truth, recon, {"m.100A>G"})
    assert m["hap_s_precision"] == 1.0
    assert m["hap_s_recall"] == 1.0


def test_haplotype_of_only_false_positives_is_dropped_by_shared_restriction():
    truth = {frozenset({"m.100A>G"})}
    recon = {frozenset({"m.100A>G"}): 4, frozenset({"m.999T>A"}): 1}
    m = sl.score_haplotypes(truth, recon, {"m.100A>G"})
    assert m["hap_precision"] == 0.5   # strict: the FP haplotype is a miss
    assert m["hap_s_precision"] == 1.0  # restricted: it collapses to empty, dropped


def test_empty_haplotype_is_excluded_from_matching_and_read_fraction():
    # The all-zero bin is dominated by reads that never span a variant site
    # (on seed5: 636 of its 1226 reads come from mutant clones), so scoring it
    # as a recovered reference clone would inflate every number.
    truth = {frozenset(), frozenset({"m.100A>G"})}
    recon = {frozenset(): 1000, frozenset({"m.100A>G"}): 1}
    m = sl.score_haplotypes(truth, recon, {"m.100A>G"})
    assert m["n_truth_clones"] == 1
    assert m["n_recon_haps"] == 1
    assert m["hap_precision"] == 1.0
    assert m["hap_read_frac"] == 1.0  # the 1000 unphased reads are not a denominator


def test_hap_read_frac_weights_haplotypes_by_read_support():
    # Same 1-of-2 haplotype precision, but almost all reads sit on the wrong
    # haplotype -- hap_read_frac is what distinguishes the two situations.
    truth = {frozenset({"m.100A>G"})}
    recon = {frozenset({"m.100A>G"}): 10, frozenset({"m.999T>A"}): 90}
    m = sl.score_haplotypes(truth, recon, {"m.100A>G"})
    assert m["hap_precision"] == 0.5
    assert abs(m["hap_read_frac"] - 0.1) < 1e-9


def test_hap_fields_are_na_when_haplotype_inputs_are_absent():
    truth_pm = sl.parse_mutation_tree(TRUTH)
    recon_pm = sl.parse_mutation_tree(RECON)
    truth_vars = {"m.100A>G", "m.200C>T", "m.300G>A", "m.400T>C"}
    m = sl.score(truth_pm, recon_pm, truth_vars, truth_vars)
    for f in sl.HAP_FIELDS:
        assert m[f] == "NA"


def test_score_uses_truth_intersect_detected_as_the_shared_restriction():
    # A truth variant the caller never detected must not be part of the shared
    # set -- otherwise a reconstruction is charged for phasing a variant that
    # was never available to it.
    truth_pm = sl.parse_mutation_tree(TRUTH)
    recon_pm = sl.parse_mutation_tree(RECON)
    truth_vars = {"m.100A>G", "m.200C>T"}
    detected = {"m.100A>G"}  # m.200C>T was missed entirely
    truth_haps = {frozenset({"m.100A>G", "m.200C>T"})}
    recon_haps = {frozenset({"m.100A>G"}): 5}
    m = sl.score(truth_pm, recon_pm, truth_vars, detected,
                 truth_haps=truth_haps, recon_haps=recon_haps)
    assert m["hap_precision"] == 0.0  # strict: the missed variant breaks the match
    assert m["hap_s_precision"] == 1.0


def _run_main(monkeypatch, tmp_path, extra_args):
    vcf = tmp_path / "x.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSIM\n"
        "chrM\t100\t.\tA\tG\t.\tPASS\t.\tGT:AD:HF\t1:3:0.3\n"
    )
    truth_variants = tmp_path / "truth_variants.txt"
    truth_variants.write_text("m.100A>G\n")
    argv = [
        "score_lineage.py",
        "--truth-tree", TRUTH, "--recon-tree", RECON,
        "--truth-variants", str(truth_variants), "--vcf", str(vcf),
    ] + extra_args
    monkeypatch.setattr(sys, "argv", argv)
    sl.main()


def test_main_skips_stdout_row_when_metrics_out_given(monkeypatch, tmp_path, capsys):
    # run_eval.sh invokes score_lineage.py with --metrics-out and then re-displays
    # that same file itself (`column -t "$OUTDIR/metrics.tsv"`) -- if main() also
    # prints the row to stdout unconditionally, the one result shows up twice in
    # the terminal. Only print to stdout when there is no file destination for it.
    metrics_out = tmp_path / "metrics.tsv"
    _run_main(monkeypatch, tmp_path, ["--metrics-out", str(metrics_out)])
    out, _ = capsys.readouterr()
    assert out == ""
    assert metrics_out.read_text().splitlines()[0].split("\t") == sl.FIELDS


def test_main_prints_row_to_stdout_when_no_metrics_out(monkeypatch, tmp_path, capsys):
    # Standalone use (no --metrics-out) has no other way to see the result, so
    # the stdout print must survive for that case.
    _run_main(monkeypatch, tmp_path, [])
    out, _ = capsys.readouterr()
    lines = out.splitlines()
    assert lines[0].split("\t") == sl.FIELDS
    assert len(lines) == 2


def test_main_emits_na_haplotype_columns_without_the_new_flags(monkeypatch, tmp_path, capsys):
    # Callers predating --recon-matrix/--truth-clones must keep working.
    _run_main(monkeypatch, tmp_path, [])
    out, _ = capsys.readouterr()
    row = dict(zip(*(l.split("\t") for l in out.splitlines())))
    assert [row[f] for f in sl.HAP_FIELDS] == ["NA"] * len(sl.HAP_FIELDS)


def test_main_scores_haplotypes_when_both_flags_are_given(monkeypatch, tmp_path, capsys):
    matrix = _write_matrix(
        tmp_path,
        ["variant", "r1", "r2", "r3"],
        [["m.100A>G", "1", "1", "0"]],
    )
    clones = _write_clones(tmp_path, [("", 0.5), ("m.100A>G", 0.5)])
    _run_main(monkeypatch, tmp_path, ["--recon-matrix", matrix, "--truth-clones", clones])
    out, _ = capsys.readouterr()
    row = dict(zip(*(l.split("\t") for l in out.splitlines())))
    assert row["n_truth_clones"] == "1"
    assert row["n_recon_haps"] == "1"
    assert row["hap_f1"] == "1.0000"


def test_main_rejects_one_haplotype_flag_without_the_other(monkeypatch, tmp_path):
    clones = _write_clones(tmp_path, [("m.100A>G", 1.0)])
    with pytest.raises(SystemExit):
        _run_main(monkeypatch, tmp_path, ["--truth-clones", clones])
