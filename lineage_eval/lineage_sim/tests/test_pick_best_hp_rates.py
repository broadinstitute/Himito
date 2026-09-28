# Plan 2026-09-27 T11. Noise-aware selection of --hp-fp-rate x --hp-fn-rate.
import csv
import math
import pathlib
import sys

import pytest

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import pick_best_hp_rates as pick


def _write(path, rows):
    cols = list(rows[0])
    with open(path, "w") as fh:
        writer = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _row(seed, hp_fp, hp_fn, intrusion, ad_f1, root, fp="0.0001", fn="0.01"):
    return {
        "seed": str(seed),
        "hp_fp_rate": str(hp_fp),
        "hp_fn_rate": str(hp_fn),
        "fp": fp,
        "fn": fn,
        "hp_intrusion": f"{intrusion:.4f}",
        "ad_f1": f"{ad_f1:.4f}",
        "root_fraction": f"{root:.4f}",
    }


def _run(tmp_path, rows, **kwargs):
    metrics = tmp_path / "metrics.tsv"
    summary = tmp_path / "hp_sweep_summary.tsv"
    best = tmp_path / "best_hp_rates.json"
    report = tmp_path / "hp_sweep_report.md"
    _write(metrics, rows)
    result = pick.pick_rates(
        str(metrics),
        summary=str(summary),
        json_out=str(best),
        report=str(report),
        **kwargs,
    )
    text = report.read_text()
    table = list(csv.DictReader(summary.open(), delimiter="\t"))
    return result, text, table


def test_clear_best_on_hp_intrusion_is_chosen(tmp_path):
    """(a) One cell is clearly best on hp_intrusion."""
    rows = []
    for seed, a, b in ((1, 0.10, 0.50), (2, 0.11, 0.51), (3, 0.09, 0.49)):
        rows.append(_row(seed, 0.0001, 0.01, a, 0.80, 0.20))
        rows.append(_row(seed, 0.0016, 0.32, b, 0.70, 0.55))
    result, text, _table = _run(tmp_path, rows)
    assert result["flat"] is False
    assert result["best"]["hp_fp_rate"] == pytest.approx(0.0001)
    assert result["best"]["hp_fn_rate"] == pytest.approx(0.01)
    assert result["hp_fp_rate"] == pytest.approx(0.0001)
    assert result["hp_fn_rate"] == pytest.approx(0.01)
    assert result["objective"] == "hp_intrusion"
    assert result["root_fraction_best"] == pytest.approx(0.20)
    assert result["root_fraction_global"] == pytest.approx(0.20)
    assert "chosen" in text.lower()


def test_primary_tie_breaks_on_ad_f1(tmp_path):
    """(b) Two cells tied on the primary objective; the better ad_f1 wins."""
    rows = []
    for seed in (1, 2, 3):
        rows.append(_row(seed, 0.0001, 0.01, 0.20, 0.50, 0.40))
        rows.append(_row(seed, 0.0004, 0.08, 0.20, 0.90, 0.25))
        rows.append(_row(seed, 0.0016, 0.32, 0.80, 0.95, 0.60))
    result, _text, _table = _run(tmp_path, rows)
    assert result["flat"] is False
    assert result["best"]["hp_fp_rate"] == pytest.approx(0.0004)
    assert result["best"]["hp_fn_rate"] == pytest.approx(0.08)


def test_flat_grid_refuses_a_best_cell(tmp_path):
    """(c) Best and worst sit inside the noise band; best is null."""
    rows = []
    for seed, a, b in ((1, 0.40, 0.45), (2, 0.50, 0.50), (3, 0.60, 0.55)):
        rows.append(_row(seed, 0.0001, 0.01, a, 0.80, 0.30))
        rows.append(_row(seed, 0.0016, 0.32, b, 0.70, 0.35))
    result, text, _table = _run(tmp_path, rows)
    assert result["best"] is None
    assert result["flat"] is True
    assert "flat" in text


def test_summary_se_is_sd_over_sqrt_n(tmp_path):
    """(d) A 3-seed cell's SE is the sample sd divided by sqrt(3)."""
    vals = (0.1, 0.2, 0.3)
    rows = [_row(seed, 0.0001, 0.01, v, 0.5, 0.2) for seed, v in enumerate(vals, start=1)]
    _result, _text, table = _run(tmp_path, rows)
    assert len(table) == 1
    mean = sum(vals) / 3
    sd = math.sqrt(sum((v - mean) ** 2 for v in vals) / 2)
    assert float(table[0]["n_seeds"]) == 3
    assert float(table[0]["hp_intrusion_se"]) == pytest.approx(sd / math.sqrt(3))
    assert float(table[0]["hp_intrusion_mean"]) == pytest.approx(mean)


def test_missing_hp_fn_rate_column_errors(tmp_path):
    """(e) A table without hp_fn_rate fails with that column named."""
    metrics = tmp_path / "metrics.tsv"
    metrics.write_text(
        "seed\thp_fp_rate\tfp\tfn\thp_intrusion\tad_f1\troot_fraction\n"
        "1\t0.0001\t0.0001\t0.01\t0.2\t0.5\t0.3\n"
    )
    with pytest.raises(ValueError, match="hp_fn_rate"):
        pick.pick_rates(
            str(metrics),
            summary=str(tmp_path / "summary.tsv"),
            json_out=str(tmp_path / "best.json"),
            report=str(tmp_path / "report.md"),
        )
