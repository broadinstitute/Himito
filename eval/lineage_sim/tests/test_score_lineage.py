# eval/lineage_sim/tests/test_score_lineage.py
import pathlib
import sys

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


def test_detected_variants_from_vcf_reads_pass_snvs(tmp_path):
    vcf = tmp_path / "x.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chrM\t100\t.\tA\tG\t.\tPASS\t.\n"
        "chrM\t250\t.\tC\tT\t.\tLowQual\t.\n"
    )
    got = sl.detected_variants_from_vcf(str(vcf))
    assert got == {"m.100A>G"}
