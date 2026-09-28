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
    """Plan 2026-09-27 T5/T6. Background stays put; HP columns follow it."""
    i = sl.FIELDS.index("n_background_calls")
    assert sl.FIELDS[:i] == [
        "profile", "fp", "fn", "n_truth_vars", "n_detected_vars", "n_shared",
        "n_truth_pairs",
        "var_precision", "var_recall", "var_f1",
        "ad_precision", "ad_recall", "ad_f1",
        "n_truth_clones", "n_recon_haps",
        "hap_precision", "hap_recall", "hap_f1",
    ]
    assert sl.FIELDS[i:] == [
        "n_background_calls",
        "n_hp_unstable_truth", "n_hp_unstable_called",
        "hp_var_precision", "hp_var_recall", "hp_var_f1",
        "hp_intrusion", "hp_backbone_nodes", "hp_len_l1", "root_fraction",
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
    # 16183 and 16189 sit inside Himito's fixed 16183–16193 window, so the
    # anchor stays 16183 even when the run extends to 16182.
    fixed_on_extended = "m.16183ACCCCCTCCCC>CCCCCCCCCCCC"
    assert hp.canonical_name("m.16183A>AC", extended) == fixed_on_extended
    assert hp.canonical_name("m.16189T>TC", extended) == fixed_on_extended


def test_canonical_name_matches_himito_window_on_plain_rcrs():
    """Plan 2026-09-27 T6. Uncollapsed and --collapse-hp spellings share the fixed window."""
    import hp_tracts as hp

    collapsed = "m.16183ACCCCCTCCCC>ACCCCCCTCCCC"
    assert hp.canonical_name("m.16183A>AC", []) == collapsed
    assert hp.canonical_name(collapsed, []) == collapsed
    # REF longer than one base, and outside the bare C-run, still collapses.
    assert hp.canonical_name("m.16190CC>C", []) == "m.16183ACCCCCTCCCC>ACCCCCTCCC"


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


def test_hp_columns_default_when_flags_omitted(tmp_path):
    """Plan 2026-09-27 T6. No HP truth files: counts are 0 and rates are NA."""
    plain = tmp_path / "plain"
    plain.mkdir()
    run = _run_score(plain)
    assert run.returncode == 0, run.stderr
    row = _metrics_row(run.stdout)
    assert row["n_hp_unstable_truth"] == "0"
    assert row["n_hp_unstable_called"] == "0"
    assert row["hp_var_precision"] == "NA"
    assert row["hp_var_recall"] == "NA"
    assert row["hp_var_f1"] == "NA"
    assert row["hp_intrusion"] == "0.0000"
    assert row["hp_backbone_nodes"] == "0"
    assert row["hp_len_l1"] == "NA"
    assert row["root_fraction"] == "NA"


def _hp_score(tmp_path, *, truth_tree, recon_tree, variants, vcf, extra):
    import subprocess
    script = pathlib.Path(__file__).resolve().parents[1] / "score_lineage.py"
    return subprocess.run(
        [sys.executable, str(script),
         "--truth-tree", str(truth_tree),
         "--recon-tree", str(recon_tree),
         "--truth-variants", str(variants),
         "--vcf", str(vcf),
         *extra],
        capture_output=True, text=True,
    )


def test_unstable_calls_are_not_false_positives(tmp_path):
    """Plan 2026-09-27 T6. Unstable and background alleles leave the precision denominator."""
    truth_tree = tmp_path / "truth.tsv"
    recon_tree = tmp_path / "recon.tsv"
    header = "node_id\tvariant\tparent_id\tparent_variant\tn_reads_attached\n"
    truth_tree.write_text(
        header
        + "0\tROOT\t0\tROOT\t8\n"
        + "1\tm.100A>G\t0\tROOT\t20\n"
        + "2\tm.300G>A\t0\tROOT\t10\n"
    )
    recon_tree.write_text(
        header
        + "0\tROOT\t0\tROOT\t8\n"
        + "1\tm.16183A>AC\t0\tROOT\t2\n"
        + "2\tm.100A>G\t1\tm.16183A>AC\t20\n"
        + "3\tm.300G>A\t0\tROOT\t10\n"
        # Keyed by variant, so this row must not reuse m.100A>G.
        + "ROOT_BLOCK\tm.50A>G,m.200C>T\t0\tROOT\t8\n"
    )
    variants = tmp_path / "truth_variants.txt"
    variants.write_text("m.100A>G\nm.300G>A\n")
    vcf = tmp_path / "sim.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        "chrM\t100\t.\tA\tG\t.\tPASS\t.\tGT:HF\t0/1:0.40\n"
        "chrM\t300\t.\tG\tA\t.\tPASS\t.\tGT:HF\t0/1:0.20\n"
        "chrM\t50\t.\tA\tG\t.\tPASS\t.\tGT:HF\t0/1:0.10\n"
        "chrM\t200\t.\tC\tT\t.\tPASS\t.\tGT:HF\t1/1:0.99\n"
        "chrM\t16189\t.\tT\tTC\t.\tPASS\t.\tGT:HF\t0/1:0.30\n"
    )
    tracts = tmp_path / "hp_tracts.tsv"
    tracts.write_text(
        "tract\tanchor\tref_len\tbackground_len\tdelta\ttruth_freq\tcanonical_name\torigin\n"
        "HVS1\t16183\t5\t10\t+1\t0.400000\tm.16183A>AC\tunstable\n"
    )
    lineage = tmp_path / "hp_lineage.tsv"
    lineage.write_text(
        "node_id\tvariant\ttract\tdelta\n"
        "1\tm.300G>A\tT455\t+1\n"
        "2\tm.400A>C\tC960\t+1\n"
    )
    bg = tmp_path / "background.txt"
    bg.write_text("m.200C>T\n")
    run = _hp_score(
        tmp_path,
        truth_tree=truth_tree, recon_tree=recon_tree,
        variants=variants, vcf=vcf,
        extra=[
            "--background", "T16189C",
            "--truth-background", str(bg),
            "--truth-hp-tracts", str(tracts),
            "--truth-hp-lineage", str(lineage),
        ],
    )
    assert run.returncode == 0, run.stderr
    row = _metrics_row(run.stdout)
    assert row["n_truth_vars"] == "2"
    assert row["n_detected_vars"] == "3"
    assert row["var_precision"] == "0.6667"
    assert row["var_recall"] == "1.0000"
    assert row["n_background_calls"] == "1"
    assert row["n_hp_unstable_truth"] == "1"
    assert row["n_hp_unstable_called"] == "1"
    assert _HVS1_CANON not in run.stderr
    assert "m.200C>T" not in run.stderr
    assert "m.50A>G" in run.stderr
    assert row["hp_var_precision"] == "1.0000"
    assert row["hp_var_recall"] == "0.5000"
    assert row["hp_var_f1"] == "0.6667"
    assert row["hp_intrusion"] == "0.5000"
    assert row["hp_backbone_nodes"] == "1"
    # ROOT 8 + tract 2 + m.100 20 + m.300 10. ROOT_BLOCK repeats the root count.
    assert row["root_fraction"] == "0.2000"
    assert row["hp_len_l1"] == "NA"


def test_root_fraction_skips_root_block_annotation_row(tmp_path):
    """Plan 2026-09-27 T6. ROOT_BLOCK repeats the root's n_reads_attached."""
    path = tmp_path / "mutation_tree.tsv"
    path.write_text(
        "node_id\tvariant\tparent_id\tparent_variant\tn_reads_attached\n"
        "0\tROOT\t0\tROOT\t10\n"
        "1\tm.100A>G\t0\tROOT\t30\n"
        "ROOT_BLOCK\tm.100A>G\t0\tROOT\t10\n"
    )
    assert sl.root_fraction(str(path)) == 0.25


def test_lineage_allele_listed_in_hp_tracts_stays_in_var_recall(tmp_path):
    """Plan 2026-09-27 T6. A name in both truth and hp_tracts.tsv is still a true positive."""
    header = "node_id\tvariant\tparent_id\tparent_variant\n"
    tree = header + "0\tROOT\t0\tROOT\n1\tm.16183A>AC\t0\tROOT\n"
    truth_tree = tmp_path / "truth.tsv"
    recon_tree = tmp_path / "recon.tsv"
    truth_tree.write_text(tree)
    recon_tree.write_text(
        header + "0\tROOT\t0\tROOT\n1\tm.16189T>TC\t0\tROOT\n"
    )
    variants = tmp_path / "truth_variants.txt"
    variants.write_text("m.16183A>AC\n")
    vcf = tmp_path / "sim.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        "chrM\t16189\t.\tT\tTC\t.\tPASS\t.\tGT:HF\t0/1:0.40\n"
    )
    tracts = tmp_path / "hp_tracts.tsv"
    tracts.write_text(
        "tract\tanchor\tref_len\tbackground_len\tdelta\ttruth_freq\tcanonical_name\torigin\n"
        "HVS1\t16183\t5\t10\t+1\t0.400000\tm.16183A>AC\tunstable\n"
    )
    run = _hp_score(
        tmp_path,
        truth_tree=truth_tree, recon_tree=recon_tree,
        variants=variants, vcf=vcf,
        extra=["--background", "T16189C", "--truth-hp-tracts", str(tracts)],
    )
    assert run.returncode == 0, run.stderr
    row = _metrics_row(run.stdout)
    assert row["var_recall"] == "1.0000"
    assert row["var_precision"] == "1.0000"
    assert row["n_detected_vars"] == "1"


def test_hp_len_l1_is_mean_l1_over_tracts(tmp_path):
    """Plan 2026-09-27 T6. δ=0 is the share of reads carrying none of the tract's alleles."""
    truth_tree = tmp_path / "truth.tsv"
    recon_tree = tmp_path / "recon.tsv"
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
    )
    tracts = tmp_path / "hp_tracts.tsv"
    tracts.write_text(
        "tract\tanchor\tref_len\tbackground_len\tdelta\ttruth_freq\tcanonical_name\torigin\n"
        "HVS1\t16183\t5\t10\t+1\t0.200000\tm.16189T>TC\tunstable\n"
        "HVS1\t16183\t5\t10\t+1\t0.300000\tm.16189T>TC\tmixed\n"
        "C960\t955\t5\t5\t+1\t0.250000\tm.960A>G\tunstable\n"
        "D310\t302\t7\t7\t0\t0.100000\tm.302A>AC\tunstable\n"
    )
    matrix = tmp_path / "matrix.csv"
    matrix.write_text(
        "variant,r1,r2,r3,r4\n"
        "m.16189T>TC,1,0,0,0\n"
        "m.960A>G,0,1,0,0\n"
    )
    clones = tmp_path / "clones.tsv"
    clones.write_text("clone_id\tnode_variant\tvariant_path\tfrequency\nref\tROOT\t\t1.0\n")
    run = _hp_score(
        tmp_path,
        truth_tree=truth_tree, recon_tree=recon_tree,
        variants=variants, vcf=vcf,
        extra=[
            "--background", "T16189C",
            "--truth-hp-tracts", str(tracts),
            "--recon-matrix", str(matrix),
            "--truth-clones", str(clones),
        ],
    )
    assert run.returncode == 0, run.stderr
    row = _metrics_row(run.stdout)
    # HVS1 truth δ=+1 is 0.2+0.3=0.5, recon 1/4, L1=0.5. C960 matches, L1=0.
    # The δ=0 row is ignored. Mean over the two tracts with a non-zero δ is 0.25.
    assert row["n_hp_unstable_truth"] == "2"
    assert row["hp_len_l1"] == "0.2500"

    zero = tmp_path / "zero.tsv"
    zero.write_text(
        "tract\tanchor\tref_len\tbackground_len\tdelta\ttruth_freq\tcanonical_name\torigin\n"
        "D310\t302\t7\t7\t0\t1.000000\tm.302A>AC\tunstable\n"
    )
    zero_run = _hp_score(
        tmp_path,
        truth_tree=truth_tree, recon_tree=recon_tree,
        variants=variants, vcf=vcf,
        extra=["--truth-hp-tracts", str(zero), "--recon-matrix", str(matrix),
               "--truth-clones", str(clones)],
    )
    assert zero_run.returncode == 0, zero_run.stderr
    assert _metrics_row(zero_run.stdout)["hp_len_l1"] == "NA"


def test_sweep_rebuilds_background_tokens_from_truth_file(tmp_path):
    """Plan 2026-09-27 T6. A sweep that omits --background still canonicalises from background.txt."""
    import subprocess

    bg = tmp_path / "background.txt"
    bg.write_text("m.16189T>C\nm.16183A>C\n")
    script = pathlib.Path(__file__).resolve().parents[1] / "sweep_fpfn.sh"
    # macOS bash 3.2 drops functions defined via `source <(...)`.
    fn = tmp_path / "fn.sh"
    lines = script.read_text().splitlines()
    start = next(i for i, line in enumerate(lines) if line.startswith("background_spec_from_truth()"))
    end = next(i for i, line in enumerate(lines[start:], start) if line == "}")
    fn.write_text("\n".join(lines[start:end + 1]) + "\n")
    proc = subprocess.run(
        ["bash", "-c", f"source '{fn}'; background_spec_from_truth '{bg}'"],
        capture_output=True, text=True,
    )
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.strip() == "T16189C,A16183C"
    text = script.read_text()
    assert 'background_spec_from_truth "$OUTDIR/truth/background.txt"' in text


def test_sweep_fpfn_passes_hp_truth_files():
    """Plan 2026-09-27 T6. The sweep scorer sees the same truth files as run_eval."""
    script = (pathlib.Path(__file__).resolve().parents[1] / "sweep_fpfn.sh").read_text()
    assert '--truth-background "$OUTDIR/truth/background.txt"' in script
    assert '--truth-hp-tracts "$OUTDIR/truth/hp_tracts.tsv"' in script
    assert '--truth-hp-lineage "$OUTDIR/truth/hp_lineage.tsv"' in script
    assert 'SCORE_ARGS+=(--background "$BACKGROUND")' in script
