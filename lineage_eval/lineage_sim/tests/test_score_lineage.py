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
