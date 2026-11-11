#!/usr/bin/env python3
"""Draw grouped box plots from lineage simulation metric TSVs."""

import argparse
import csv
import math
from collections.abc import Mapping, Sequence
from pathlib import Path

from matplotlib import pyplot as plt
from matplotlib.figure import Figure
from matplotlib.patches import Patch


SCORE_NAMES = ("Precision", "Recall", "F1")
SCORE_GROUPS = (
    ("Variant detection", "var"),
    ("Ancestor–descendant pairs", "ad"),
    ("Haplotype recovery", "hap"),
)


def read_metrics(path: str | Path, columns: Sequence[str]) -> dict[str, list[float]]:
    return _read_metric_groups(path, columns, by_coverage=False)[None]


def _read_metric_groups(
    path: str | Path, columns: Sequence[str], *, by_coverage: bool,
) -> dict[int | None, dict[str, list[float]]]:
    groups: dict[int | None, dict[str, list[float]]] = {}
    with Path(path).open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = [*columns, "total_depth"] if by_coverage else columns
        missing = [column for column in required if column not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(f"{path}: missing metric columns: {', '.join(missing)}")
        for row in reader:
            depth = None
            if by_coverage:
                try:
                    depth = int(row["total_depth"] or "")
                except ValueError:
                    raise ValueError(f"{path}:{reader.line_num}: invalid total_depth") from None
            if depth not in groups:
                groups[depth] = {column: [] for column in columns}
            for column in columns:
                cell = row[column]
                if cell is None:
                    raise ValueError(f"{path}:{reader.line_num}: truncated row")
                cell = cell.strip()
                if cell.lower() in {"", "na", "nan"}:
                    continue
                try:
                    value = float(cell)
                    if not math.isfinite(value):
                        raise ValueError
                except ValueError:
                    raise ValueError(
                        f"{path}:{reader.line_num}: invalid value {cell!r} for {column}"
                    ) from None
                groups[depth][column].append(value)
    if not any(values for group in groups.values() for values in group.values()):
        raise ValueError(f"{path}: no numeric values for the selected metrics")
    return groups


def plot_simulation_metrics(
    inputs: Mapping[str, str | Path],
    *,
    metrics: Sequence[str] | None = None,
    by_coverage: bool = False,
    title: str = "Lineage simulation metrics",
) -> Figure:
    """Return a figure with one box per input and metric, pooling rows within each TSV.

    Input labels define the legend and box order. The default panels use the nine
    variant, ancestry, and haplotype scores. Pass metrics with compatible units
    to plot other numeric columns on a shared, automatically scaled axis.
    Set by_coverage to draw one subplot per total_depth with all selected metrics.
    """
    if not inputs:
        raise ValueError("at least one labeled input TSV is required")
    if metrics is None:
        panels = [
            (heading, [f"{prefix}_{score.lower()}" for score in SCORE_NAMES], SCORE_NAMES)
            for heading, prefix in SCORE_GROUPS
        ]
    else:
        if not metrics or len(set(metrics)) != len(metrics):
            raise ValueError("metrics must contain distinct column names")
        panels = [("", list(metrics), metrics)]

    columns = [column for _, panel_columns, _ in panels for column in panel_columns]
    data = [_read_metric_groups(path, columns, by_coverage=by_coverage) for path in inputs.values()]
    if metrics is None and any(
        not 0.0 <= value <= 1.0
        for groups in data
        for series in groups.values()
        for values in series.values()
        for value in values
    ):
        raise ValueError("default score columns must be between 0 and 1")

    if by_coverage:
        depths = sorted({depth for groups in data for depth in groups if depth is not None})
        coverage_labels = SCORE_NAMES * len(SCORE_GROUPS) if metrics is None else metrics
        panels = [(f"{depth:,}× coverage", columns, coverage_labels) for depth in depths]
        ncols = min(2, len(panels))
        nrows = math.ceil(len(panels) / ncols)
        figsize = (max(6.5, 0.78 * len(columns)) * ncols, 3.1 * nrows + 0.5)
    else:
        depths = [None] * len(panels)
        nrows, ncols = 1, len(panels)
        figsize = (max(6.5, 1.32 * len(columns)), 4.2)

    figure, axes = plt.subplots(nrows, ncols, squeeze=False, sharey=True, figsize=figsize)
    palette = plt.get_cmap("tab10" if len(inputs) <= 10 else "tab20")
    colors = [palette(i % palette.N) for i in range(len(inputs))]
    spacing = 0.76 / len(inputs)
    for axis, depth, (heading, panel_columns, tick_labels) in zip(axes.flat, depths, panels):
        for i, (groups, color) in enumerate(zip(data, colors)):
            series = groups.get(depth)
            if series is None:
                continue
            offset = (i - (len(inputs) - 1) / 2) * spacing
            present = [(j + offset, series[column]) for j, column in enumerate(panel_columns)
                       if series[column]]
            if not present:
                continue
            positions, samples = zip(*present)
            axis.boxplot(
                samples, positions=positions, widths=spacing * 0.84,
                patch_artist=True, manage_ticks=False,
                boxprops={"facecolor": color, "edgecolor": "#26313b", "alpha": 0.72},
                medianprops={"color": "#18212a", "linewidth": 1.5},
                flierprops={"marker": "o", "markersize": 2.8, "alpha": 0.55,
                            "markerfacecolor": color, "markeredgecolor": color},
            )
        axis.set_title(heading, fontweight="semibold", pad=10)
        axis.set_xticks(range(len(panel_columns)), tick_labels)
        if by_coverage and metrics is None:
            for x, name in [(1, "Variants"), (4, "Ancestor–descendant"), (7, "Haplotypes")]:
                axis.text(x, -0.13, name, transform=axis.get_xaxis_transform(),
                          ha="center", va="top", fontsize=9)
            for x in (2.5, 5.5):
                axis.axvline(x, color="#dbe1e6", linewidth=0.8, zorder=0)
        if metrics is not None:
            axis.tick_params(axis="x", labelrotation=25)
            plt.setp(axis.get_xticklabels(), ha="right")
        axis.set_xlim(-0.55, len(panel_columns) - 0.45)
        axis.set_axisbelow(True)
        axis.grid(axis="y", color="#dbe1e6", linewidth=0.8)
        for spine in ("top", "right"):
            axis.spines[spine].set_visible(False)

    for axis in axes.flat[len(panels):]:
        figure.delaxes(axis)
    figure.supylabel("Score" if metrics is None else "Value")
    if metrics is None:
        axes[0, 0].set_ylim(-0.025, 1.045)
        axes[0, 0].set_yticks([0, 0.25, 0.5, 0.75, 1])
    figure.suptitle(title, fontsize=16, fontweight="semibold")
    handles = [Patch(facecolor=color, edgecolor="#26313b", alpha=0.72, label=label)
               for label, color in zip(inputs, colors)]
    figure.legend(handles=handles, loc="lower center", ncol=min(len(inputs), 4), frameon=False)
    figure.canvas.draw()
    legend = figure.legends[0]
    bottom = legend.get_window_extent().transformed(figure.transFigure.inverted()).y1
    figure.tight_layout(rect=(0, bottom + 0.01, 1, 0.98))
    return figure


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", action="append", nargs=2, metavar=("LABEL", "TSV"), required=True,
                        help="labeled metric TSV, repeat for each group")
    parser.add_argument("--metrics", nargs="+", metavar="COLUMN", help="plot these numeric columns")
    parser.add_argument("--by-coverage", action="store_true", help="draw one subplot per total_depth")
    parser.add_argument("--title", default="Lineage simulation metrics")
    parser.add_argument("--output", "-o", type=Path, required=True, help="figure path, e.g. metrics.png or metrics.pdf")
    args = parser.parse_args()
    inputs = dict(args.input)
    if len(inputs) != len(args.input):
        parser.error("input labels must be unique")
    try:
        figure = plot_simulation_metrics(
            inputs, metrics=args.metrics, by_coverage=args.by_coverage, title=args.title,
        )
        figure.savefig(args.output, dpi=240, bbox_inches="tight")
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    finally:
        if "figure" in locals():
            plt.close(figure)


if __name__ == "__main__":
    main()
