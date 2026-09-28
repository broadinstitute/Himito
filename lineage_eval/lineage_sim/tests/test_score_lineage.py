# eval/lineage_sim/tests/test_score_lineage.py
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import score_lineage as sl


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


def test_n_background_calls_is_appended_to_fields():
    """Plan 2026-09-27 T5. The new column is last; earlier columns stay put."""
    assert sl.FIELDS[-1] == "n_background_calls"
    assert sl.FIELDS[:-1] == [
        "profile", "fp", "fn", "n_truth_vars", "n_detected_vars", "n_shared",
        "n_truth_pairs",
        "var_precision", "var_recall", "var_f1",
        "ad_precision", "ad_recall", "ad_f1",
        "n_truth_clones", "n_recon_haps",
        "hap_precision", "hap_recall", "hap_f1",
    ]


def _score_fixtures(tmp_path):
    truth_tree = tmp_path / "truth_mutation_tree.tsv"
    recon_tree = tmp_path / "recon_mutation_tree.tsv"
    tree = (
        "node_id\tvariant\tparent_id\tparent_variant\n"
        "0\tROOT\t0\tROOT\n"
        "1\tm.100A>G\t0\tROOT\n"
    )
    truth_tree.write_text(tree)
    recon_tree.write_text(tree)
    variants = tmp_path / "truth_variants.txt"
    variants.write_text("m.100A>G\n")
    vcf = tmp_path / "sim.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        "chrM\t100\t.\tA\tG\t.\tPASS\t.\tGT:HF\t0/1:0.50\n"
        "chrM\t16189\t.\tT\tC\t.\tPASS\t.\tGT:HF\t1/1:0.98\n"
        "chrM\t50\t.\tA\tG\t.\tPASS\t.\tGT:HF\t0/1:0.20\n"
    )
    return truth_tree, recon_tree, variants, vcf


def _run_score(tmp_path, *extra):
    import subprocess
    truth_tree, recon_tree, variants, vcf = _score_fixtures(tmp_path)
    script = pathlib.Path(__file__).resolve().parents[1] / "score_lineage.py"
    return subprocess.run(
        [sys.executable, str(script),
         "--truth-tree", str(truth_tree),
         "--recon-tree", str(recon_tree),
         "--truth-variants", str(variants),
         "--vcf", str(vcf),
         "--profile", "ont", "--fp", "0.1", "--fn", "0.1",
         *extra],
        capture_output=True, text=True,
    )


def _metrics_row(stdout: str) -> dict[str, str]:
    lines = [ln for ln in stdout.splitlines() if ln and not ln.startswith("false")]
    header, values = lines[0].split("\t"), lines[1].split("\t")
    return dict(zip(header, values))


def test_truth_background_dropped_before_var_precision(tmp_path):
    """Plan 2026-09-27 T5. Background names leave detected before precision.

    n_background_calls counts how many of those names were in the VCF.
    A missing or empty path leaves detected unchanged and writes 0.
    """
    named = tmp_path / "named"
    named.mkdir()
    bg = named / "background.txt"
    bg.write_text("m.16189T>C\nm.16183A>C\n")
    hit = _run_score(named, "--truth-background", str(bg))
    assert hit.returncode == 0, hit.stderr
    row = _metrics_row(hit.stdout)
    assert row["n_background_calls"] == "1"
    assert row["var_precision"] == "0.5000"
    assert row["n_detected_vars"] == "2"
    assert row["var_recall"] == "1.0000"
    assert "m.16189T>C" not in hit.stderr

    plain = tmp_path / "plain"
    plain.mkdir()
    omitted = _run_score(plain)
    assert omitted.returncode == 0, omitted.stderr
    base = _metrics_row(omitted.stdout)
    assert base["n_background_calls"] == "0"
    assert base["var_precision"] == "0.3333"
    assert base["n_detected_vars"] == "3"

    empty_dir = tmp_path / "empty"
    empty_dir.mkdir()
    empty = empty_dir / "background.txt"
    empty.write_text("\n")
    empty_run = _run_score(empty_dir, "--truth-background", str(empty))
    assert empty_run.returncode == 0, empty_run.stderr
    assert _metrics_row(empty_run.stdout)["n_background_calls"] == "0"
    assert _metrics_row(empty_run.stdout)["var_precision"] == "0.3333"

    missing_dir = tmp_path / "missing"
    missing_dir.mkdir()
    missing = _run_score(
        missing_dir, "--truth-background", str(missing_dir / "nope.txt"),
    )
    assert missing.returncode == 0, missing.stderr
    assert _metrics_row(missing.stdout)["n_background_calls"] == "0"
    assert _metrics_row(missing.stdout)["var_precision"] == "0.3333"


# Plan 2026-09-27 T6.
# rCRS 16183-16193 is ACCCCCTCCCC. T16189C plus one inserted C is A + 11 C.
_HVS1_CANON = "m.16183ACCCCCTCCCC>ACCCCCCCCCCC"
# rCRS 16182-16193 is AACCCCCTCCCC. T16189C, A16183C and one inserted C
# is A followed by 12 C.
_HVS1_EXTENDED = "m.16182AACCCCCTCCCC>ACCCCCCCCCCCC"


def test_canonical_name_collapses_equivalent_hvs1_spellings():
    """Plan 2026-09-27 T6. One HVS1 allele, three spellings, two backgrounds."""
    import hp_tracts as hp

    bg = ["T16189C"]
    assert hp.canonical_name("m.16189T>TC", bg) == _HVS1_CANON
    assert hp.canonical_name("m.16183A>AC", bg) == _HVS1_CANON
    extended = ["T16189C", "A16183C"]
    assert hp.canonical_name("m.16182A>AC", extended) == _HVS1_EXTENDED
    assert hp.canonical_name("m.16183A>AC", extended) == _HVS1_EXTENDED
    assert hp.canonical_name("m.16189T>TC", extended) == _HVS1_EXTENDED


def test_canonical_name_leaves_variants_outside_tracts():
    import hp_tracts as hp

    assert hp.canonical_name("m.3144A>G", ["T16189C"]) == "m.3144A>G"


def test_canonical_name_returns_collapsed_allele_unchanged():
    import hp_tracts as hp

    assert hp.canonical_name(_HVS1_CANON, ["T16189C"]) == _HVS1_CANON


def test_score_matches_equivalent_hp_alleles_when_background_extends_tract():
    """Plan 2026-09-27 T6. Same insertion, two spellings: recall 1 with the background, 0 without."""
    truth = {"m.16183A>AC"}
    detected = {"m.16189T>TC"}
    parents = {"ROOT": "ROOT"}
    with_bg = sl.score(parents, parents, truth, detected, background=["T16189C"])
    assert with_bg["var_recall"] == 1.0
    without = sl.score(parents, parents, truth, detected)
    assert without["var_recall"] == 0.0


def test_parent_maps_and_root_block_use_canonical_names():
    """Plan 2026-09-27 T6. Keys, values, and comma-joined ROOT_BLOCK members collide."""
    bg = ["T16189C"]
    block_truth = "m.16183A>AC,m.100A>G"
    block_recon = "m.16189T>TC,m.100A>G"
    truth_parent = {
        "ROOT": "ROOT",
        block_truth: "ROOT",
        "m.200C>T": block_truth,
    }
    recon_parent = {
        "ROOT": "ROOT",
        block_recon: "ROOT",
        "m.200C>T": block_recon,
    }
    matched = sl.score(
        truth_parent, recon_parent,
        {"m.200C>T"}, {"m.200C>T"},
        background=bg,
    )
    assert matched["n_shared"] == 2
    assert matched["ad_recall"] == 1.0
    raw = sl.score(
        truth_parent, recon_parent,
        {"m.200C>T"}, {"m.200C>T"},
    )
    assert raw["n_shared"] == 1
    assert raw["ad_recall"] == 0.0


def test_parse_matrix_haplotypes_canonicalises_row_names(tmp_path):
    """Plan 2026-09-27 T6. A matrix row spelled the other way still joins the keep set."""
    import hp_tracts as hp

    bg = ["T16189C"]
    canon = hp.canonical_name("m.16189T>TC", bg)
    path = tmp_path / "matrix.csv"
    path.write_text(
        "variant,r1,r2\n"
        "m.16189T>TC,1,0\n"
        "m.100A>G,1,1\n"
    )
    haps = sl.parse_matrix_haplotypes(str(path), keep={canon, "m.100A>G"}, background=bg)
    assert frozenset({canon, "m.100A>G"}) in haps
    assert frozenset({"m.100A>G"}) in haps


def test_cli_background_canonicalises_truth_and_vcf(tmp_path):
    """Plan 2026-09-27 T6. --background makes the two HVS1 spellings one variant."""
    import subprocess

    truth_tree = tmp_path / "truth_mutation_tree.tsv"
    recon_tree = tmp_path / "recon_mutation_tree.tsv"
    truth_body = (
        "node_id\tvariant\tparent_id\tparent_variant\n"
        "0\tROOT\t0\tROOT\n"
        "1\tm.16183A>AC\t0\tROOT\n"
    )
    recon_body = (
        "node_id\tvariant\tparent_id\tparent_variant\n"
        "0\tROOT\t0\tROOT\n"
        "1\tm.16189T>TC\t0\tROOT\n"
    )
    truth_tree.write_text(truth_body)
    recon_tree.write_text(recon_body)
    variants = tmp_path / "truth_variants.txt"
    variants.write_text("m.16183A>AC\n")
    vcf = tmp_path / "sim.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        "chrM\t16189\t.\tT\tTC\t.\tPASS\t.\tGT:HF\t0/1:0.40\n"
    )
    script = pathlib.Path(__file__).resolve().parents[1] / "score_lineage.py"
    base = [
        sys.executable, str(script),
        "--truth-tree", str(truth_tree),
        "--recon-tree", str(recon_tree),
        "--truth-variants", str(variants),
        "--vcf", str(vcf),
    ]
    hit = subprocess.run([*base, "--background", "T16189C"], capture_output=True, text=True)
    assert hit.returncode == 0, hit.stderr
    assert _metrics_row(hit.stdout)["var_recall"] == "1.0000"
    assert _metrics_row(hit.stdout)["n_shared"] == "1"
    miss = subprocess.run(base, capture_output=True, text=True)
    assert miss.returncode == 0, miss.stderr
    assert _metrics_row(miss.stdout)["var_recall"] == "0.0000"
