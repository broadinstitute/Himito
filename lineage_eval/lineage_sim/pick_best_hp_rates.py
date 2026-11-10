#!/usr/bin/env python3
"""Pick an HP false-positive / false-negative pair from a sweep table.

Plan 2026-09-27 T11. Selection is noise-aware: a cell within 2·SE of a better
cell is tied with it (the same band bench_standard.sh uses). A grid whose best
and worst primary means sit inside that band is flat, and no rate is declared.

Closest to the global rates means the smallest
|log2(hp_fn/fn)| + |log2(hp_fp/fp)|. On the default grid every multiplier is
>= 1, so that is the same as the smallest log2(hp_fn/fn) + log2(hp_fp/fp).
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys

# Same focus set bench_hp.sh averages. NA cells are left out of the mean.
FOCUS = [
    "var_precision", "var_recall", "var_f1",
    "ad_precision", "ad_recall", "ad_f1",
    "hap_precision", "hap_recall", "hap_f1",
    "hp_intrusion", "hp_backbone_nodes", "hp_len_l1",
    "root_fraction", "hp_var_f1",
]

# A smaller mean is the improvement. ad_f1 and hap_f1 are maximised.
MINIMIZE = {"hp_intrusion", "hp_len_l1", "hp_backbone_nodes", "root_fraction"}

SELECTABLE = {"hp_intrusion", "ad_f1", "hap_f1", "hp_len_l1"}


def _mean(vals: list[float]) -> float:
    return sum(vals) / len(vals)


def _sd(vals: list[float]) -> float:
    """Sample standard deviation. One observation has sd 0 (no scatter to estimate)."""
    n = len(vals)
    if n < 2:
        return 0.0
    m = _mean(vals)
    return math.sqrt(sum((x - m) ** 2 for x in vals) / (n - 1))


def _num(text) -> float | None:
    if text is None:
        return None
    s = str(text).strip()
    if s in ("", "NA", "nan", "NaN"):
        return None
    return float(s)


def _fmt(x, digits: str = ".6g") -> str:
    if x is None:
        return "NA"
    return format(x, digits)


def _load(path: str):
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"metrics table {path} has no header")
        fields = list(reader.fieldnames)
        missing = [c for c in ("seed", "hp_fp_rate", "hp_fn_rate") if c not in fields]
        if missing:
            raise ValueError(
                "metrics table missing column(s): " + ", ".join(missing)
            )
        rows = list(reader)
    if not rows:
        raise ValueError(f"metrics table {path} has no data rows")
    return fields, rows


def _check_objective(name: str, fields: list[str], role: str) -> None:
    allowed = SELECTABLE if role == "objective" else SELECTABLE | set(FOCUS)
    if name not in allowed:
        raise ValueError(
            f"{role} must be one of {', '.join(sorted(allowed))}, got {name!r}"
        )
    if name not in fields:
        raise ValueError(f"metrics table missing column(s): {name}")


def _aggregate(rows: list[dict], metrics: list[str]) -> list[dict]:
    groups: dict[tuple[str, str], list[dict]] = {}
    order: list[tuple[str, str]] = []
    for row in rows:
        key = (row["hp_fp_rate"].strip(), row["hp_fn_rate"].strip())
        if key not in groups:
            order.append(key)
            groups[key] = []
        groups[key].append(row)

    cells = []
    for hp_fp_s, hp_fn_s in order:
        members = groups[(hp_fp_s, hp_fn_s)]
        stats = {}
        for metric in metrics:
            vals = []
            for row in members:
                if metric not in row:
                    continue
                value = _num(row.get(metric))
                if value is not None:
                    vals.append(value)
            if not vals:
                stats[metric] = {"mean": None, "sd": None, "se": None, "n": 0}
                continue
            sd = _sd(vals)
            n = len(vals)
            stats[metric] = {
                "mean": _mean(vals),
                "sd": sd,
                "se": sd / math.sqrt(n),
                "n": n,
            }
        fp = fn = None
        for row in members:
            if fp is None and "fp" in row:
                fp = _num(row.get("fp"))
            if fn is None and "fn" in row:
                fn = _num(row.get("fn"))
        cells.append({
            "hp_fp_rate": float(hp_fp_s),
            "hp_fn_rate": float(hp_fn_s),
            "hp_fp_s": hp_fp_s,
            "hp_fn_s": hp_fn_s,
            "n_seeds": len(members),
            "stats": stats,
            "fp": fp,
            "fn": fn,
        })
    return cells


def _global_rates(cells: list[dict]) -> tuple[float, float]:
    fps = {c["fp"] for c in cells if c["fp"] is not None}
    fns = {c["fn"] for c in cells if c["fn"] is not None}
    if len(fps) != 1 or len(fns) != 1:
        raise ValueError(
            "metrics fp and fn columns must hold one global rate pair "
            f"(got fp={sorted(x for x in fps if x is not None)} "
            f"fn={sorted(x for x in fns if x is not None)}). "
            "Those are the profile rates, not the HP rates."
        )
    fp, fn = fps.pop(), fns.pop()
    if fp <= 0 or fn <= 0:
        raise ValueError(f"global fp/fn must be positive (got fp={fp} fn={fn})")
    return fp, fn


def _distance(hp_fp: float, hp_fn: float, fp: float, fn: float) -> float:
    """Absolute log2 steps away from the global rates. Zero is the 1x cell."""
    if min(hp_fp, hp_fn, fp, fn) <= 0:
        return float("inf")
    return abs(math.log2(hp_fn / fn)) + abs(math.log2(hp_fp / fp))


def _raw_log_sum(hp_fp: float, hp_fn: float, fp: float, fn: float) -> float:
    if min(hp_fp, hp_fn, fp, fn) <= 0:
        return float("inf")
    return math.log2(hp_fn / fn) + math.log2(hp_fp / fp)


def _band(a: dict, b: dict) -> float:
    """2·SE_pooled. Matches bench_standard.sh: sqrt(sd_a^2/n_a + sd_b^2/n_b) * 2."""
    if a["n"] <= 0 or b["n"] <= 0:
        return float("inf")
    se = math.sqrt((a["sd"] ** 2) / a["n"] + (b["sd"] ** 2) / b["n"])
    return 2.0 * se


def _within(a: dict, b: dict, metric: str) -> bool:
    sa, sb = a["stats"][metric], b["stats"][metric]
    if sa["mean"] is None or sb["mean"] is None:
        return False
    return abs(sa["mean"] - sb["mean"]) <= _band(sa, sb)


def _usable(cells: list[dict], metric: str) -> list[dict]:
    return [c for c in cells if c["stats"][metric]["n"] > 0]


def _extreme(cells: list[dict], metric: str, want_best: bool) -> dict:
    usable = _usable(cells, metric)
    if not usable:
        raise ValueError(f"no numeric values for {metric}")
    minimize = metric in MINIMIZE
    # best+minimize and worst+maximize both want the smaller mean.
    want_small = minimize == want_best

    def key(cell):
        mean = cell["stats"][metric]["mean"]
        return (mean if want_small else -mean, cell["hp_fp_rate"], cell["hp_fn_rate"])

    return min(usable, key=key)


def _tied_with(cells: list[dict], ref: dict, metric: str) -> list[dict]:
    return [c for c in _usable(cells, metric) if _within(c, ref, metric)]


def _closeness_key(cell: dict, fp: float, fn: float):
    return (
        _distance(cell["hp_fp_rate"], cell["hp_fn_rate"], fp, fn),
        _raw_log_sum(cell["hp_fp_rate"], cell["hp_fn_rate"], fp, fn),
        cell["hp_fn_rate"],
        cell["hp_fp_rate"],
    )


def _stat_mean(cell: dict, metric: str):
    st = cell["stats"].get(metric)
    if not st or st["n"] <= 0:
        return None
    return st["mean"]


def _sense(metric: str) -> str:
    return "minimised" if metric in MINIMIZE else "maximised"


def _summary_metrics(fields: list[str], objective: str, secondary: str) -> list[str]:
    present = [m for m in FOCUS if m in fields]
    for extra in (objective, secondary):
        if extra not in present and extra in fields:
            present.append(extra)
    return present


def _write_summary(path: str, cells: list[dict], metrics: list[str]) -> None:
    cols = ["hp_fp_rate", "hp_fn_rate", "n_seeds"]
    for metric in metrics:
        cols.append(f"{metric}_mean")
        cols.append(f"{metric}_se")
    ordered = sorted(cells, key=lambda c: (c["hp_fp_rate"], c["hp_fn_rate"]))
    with open(path, "w") as fh:
        writer = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for cell in ordered:
            row = {
                "hp_fp_rate": cell["hp_fp_s"],
                "hp_fn_rate": cell["hp_fn_s"],
                "n_seeds": str(cell["n_seeds"]),
            }
            for metric in metrics:
                st = cell["stats"][metric]
                row[f"{metric}_mean"] = "NA" if st["mean"] is None else f"{st['mean']:.10g}"
                row[f"{metric}_se"] = "NA" if st["se"] is None else f"{st['se']:.10g}"
            writer.writerow(row)


def _cell_record(cell: dict, objective: str) -> dict:
    st = cell["stats"][objective]
    return {
        "hp_fp_rate": cell["hp_fp_rate"],
        "hp_fn_rate": cell["hp_fn_rate"],
        "objective": objective,
        "value": st["mean"],
        "se": st["se"],
        "n_seeds": cell["n_seeds"],
    }


def _describe(cell: dict, metric: str) -> str:
    st = cell["stats"][metric]
    return (
        f"hp_fp={cell['hp_fp_s']} hp_fn={cell['hp_fn_s']} "
        f"mean={_fmt(st['mean'], '.4f')} se={_fmt(st['se'], '.4f')} "
        f"n_seeds={cell['n_seeds']}"
    )


def _report(cells, metrics, objective, secondary, fp, fn, best_primary, worst_primary,
            tied_primary, tied_secondary, chosen, global_cell, flat) -> str:
    lines = [
        "# HP rate sweep",
        "",
        f"Primary objective `{objective}` ({_sense(objective)}). "
        f"Secondary `{secondary}` ({_sense(secondary)}).",
        f"Global rates fp={_fmt(fp, '.6g')} fn={_fmt(fn, '.6g')}. "
        "A cell within 2·SE of the current best mean is tied with it.",
        "",
        "## Selection trace",
        "",
        "1. Primary. "
        f"Best {_describe(best_primary, objective)}. "
        f"Worst {_describe(worst_primary, objective)}. "
        f"Tied with the best ({len(tied_primary)} cell(s)): "
        + ", ".join(f"({c['hp_fp_s']}, {c['hp_fn_s']})" for c in tied_primary)
        + ".",
    ]
    if flat:
        band = _band(best_primary["stats"][objective], worst_primary["stats"][objective])
        diff = abs(
            best_primary["stats"][objective]["mean"]
            - worst_primary["stats"][objective]["mean"]
        )
        lines += [
            "",
            "4. The grid is flat. The primary objective's best and worst cells "
            f"differ by {_fmt(diff, '.4f')}, which is within 2·SE ({_fmt(band, '.4f')}). "
            "No parameter should change. Leave the default alone. "
            "Do not pick the top row of a tie.",
            "",
        ]
    else:
        lines += [
            "2. Secondary, among cells tied on the primary. "
            f"Still tied ({len(tied_secondary)} cell(s)): "
            + ", ".join(f"({c['hp_fp_s']}, {c['hp_fn_s']})" for c in tied_secondary)
            + ".",
            "3. Among cells still tied, closest to the global rates: "
            "smallest |log2(hp_fn/fn)| + |log2(hp_fp/fp)| "
            "(then the raw log2 sum, so a >=1x grid prefers the 1x corner). "
            f"Chosen cell: {_describe(chosen, objective)}.",
            "4. Not flat. Best and worst primary means are outside 2·SE, "
            "so a rate change is distinguishable from seed noise.",
            "",
        ]
    lines += [
        "## Grid",
        "",
        "| hp_fp_rate | hp_fn_rate | n_seeds | "
        + " | ".join(f"{m} mean" for m in metrics)
        + " | |",
        "|" + " --- |" * (4 + len(metrics)),
    ]
    for cell in sorted(cells, key=lambda c: (c["hp_fp_rate"], c["hp_fn_rate"])):
        mark = ""
        if chosen is not None and cell is chosen:
            mark = "chosen"
        elif (
            global_cell is not None
            and cell["hp_fp_rate"] == global_cell["hp_fp_rate"]
            and cell["hp_fn_rate"] == global_cell["hp_fn_rate"]
        ):
            mark = "global"
        means = " | ".join(_fmt(_stat_mean(cell, m), ".4f") for m in metrics)
        lines.append(
            f"| {cell['hp_fp_s']} | {cell['hp_fn_s']} | {cell['n_seeds']} | {means} | {mark} |"
        )
    lines += [
        "",
        f"root_fraction at the chosen cell: {_fmt(_stat_mean(chosen, 'root_fraction') if chosen else None, '.4f')}.",
        f"root_fraction at the global-rate cell: {_fmt(_stat_mean(global_cell, 'root_fraction') if global_cell else None, '.4f')}.",
        "",
    ]
    return "\n".join(lines)


def pick_rates(
    metrics_path: str,
    summary: str,
    json_out: str,
    report: str,
    objective: str = "hp_intrusion",
    secondary: str = "ad_f1",
) -> dict:
    """Read a combined metrics table and write summary, JSON, and the report.

    Returns the JSON object. ``best`` is null when the grid is flat.
    """
    fields, rows = _load(metrics_path)
    _check_objective(objective, fields, "objective")
    _check_objective(secondary, fields, "secondary")
    if "fp" not in fields or "fn" not in fields:
        raise ValueError(
            "metrics table missing column(s): "
            + ", ".join(c for c in ("fp", "fn") if c not in fields)
        )

    shown = _summary_metrics(fields, objective, secondary)
    # Always aggregate the objectives plus root_fraction, even if a focus
    # column is absent from this table.
    needed = list(dict.fromkeys(shown + [objective, secondary, "root_fraction"]))
    needed = [m for m in needed if m in fields]
    cells = _aggregate(rows, needed)
    fp, fn = _global_rates(cells)
    for cell in cells:
        cell["fp"] = fp
        cell["fn"] = fn

    best_primary = _extreme(cells, objective, want_best=True)
    worst_primary = _extreme(cells, objective, want_best=False)
    flat = _within(best_primary, worst_primary, objective)
    tied_primary = _tied_with(cells, best_primary, objective)
    global_cell = min(cells, key=lambda c: _closeness_key(c, fp, fn))

    chosen = None
    tied_secondary: list[dict] = []
    if not flat:
        sec_best = _extreme(tied_primary, secondary, want_best=True)
        tied_secondary = _tied_with(tied_primary, sec_best, secondary)
        chosen = min(tied_secondary, key=lambda c: _closeness_key(c, fp, fn))

    _write_summary(summary, cells, shown)

    if flat:
        result = {
            "best": None,
            "hp_fp_rate": None,
            "hp_fn_rate": None,
            "objective": objective,
            "value": None,
            "se": None,
            "n_seeds": None,
            "flat": True,
            "root_fraction_best": None,
            "root_fraction_global": _stat_mean(global_cell, "root_fraction"),
            "reason": (
                "The grid is flat: the primary objective's best and worst cells "
                "are within 2·SE, so no parameter should change. "
                "Leave the default alone."
            ),
        }
    else:
        record = _cell_record(chosen, objective)
        result = {
            "best": record,
            "hp_fp_rate": record["hp_fp_rate"],
            "hp_fn_rate": record["hp_fn_rate"],
            "objective": record["objective"],
            "value": record["value"],
            "se": record["se"],
            "n_seeds": record["n_seeds"],
            "flat": False,
            "root_fraction_best": _stat_mean(chosen, "root_fraction"),
            "root_fraction_global": _stat_mean(global_cell, "root_fraction"),
        }

    text = _report(
        cells, shown, objective, secondary, fp, fn,
        best_primary, worst_primary, tied_primary, tied_secondary,
        chosen, global_cell, flat,
    )
    with open(json_out, "w") as fh:
        json.dump(result, fh, indent=2)
        fh.write("\n")
    with open(report, "w") as fh:
        fh.write(text)
    return result


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--metrics", required=True, help="combined metrics.tsv")
    ap.add_argument("--summary", required=True, help="hp_sweep_summary.tsv")
    ap.add_argument("--json", dest="json_out", required=True, help="best_hp_rates.json")
    ap.add_argument("--report", required=True, help="hp_sweep_report.md")
    ap.add_argument("--objective", default="hp_intrusion")
    ap.add_argument("--secondary", default="ad_f1")
    args = ap.parse_args(argv)
    try:
        pick_rates(
            args.metrics,
            summary=args.summary,
            json_out=args.json_out,
            report=args.report,
            objective=args.objective,
            secondary=args.secondary,
        )
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
