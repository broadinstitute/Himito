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
