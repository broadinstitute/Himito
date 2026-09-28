# Plan 2026-09-27 T1
"""rCRS homopolymer tract catalogue: runs, background anchors, Rust parity."""
import inspect
import os
import pathlib
import re
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import hp_tracts as hp
import simulate_tree as st

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
REF = os.path.join(_REPO, "rCRS.fasta")
TRACTS_RS = os.path.join(_REPO, "src", "tracts.rs")

# Plan table, 1-based inclusive. anchor is the base before the run.
_EXPECTED = [
    ("D310", 303, 309, "C", 7, 302, True, True),
    ("HVS1", 16184, 16188, "C", 5, 16183, True, True),
    ("C573", 568, 573, "C", 6, 567, False, True),
    ("T455", 452, 455, "T", 4, 451, False, True),
    ("C960", 956, 960, "C", 5, 955, False, True),
    ("C5899", 5895, 5899, "C", 5, 5894, False, True),
    ("A12418", 12418, 12425, "A", 8, 12417, False, True),
    ("CA514", 514, 523, "CA", 5, 513, True, True),
    ("R8281", 8272, 8289, "CCCCCTCTA", 2, 8271, False, False),
]

_RUNS_AT_LEAST_6 = (
    (303, 309, "C", 7),
    (568, 573, "C", 6),
    (12418, 12425, "A", 8),
    (9478, 9484, "T", 7),
)

# Rust KNOWN_TRACTS name -> catalogue id. Windows are wider than the run.
_RUST_TO_CATALOGUE = {"D310": "D310", "HVS1-C": "HVS1", "CA514": "CA514"}
_TRACT_RE = re.compile(
    r'Tract \{\s*name:\s*"([^"]+)",\s*start:\s*(\d+),\s*end:\s*(\d+),\s*ref_seq:\s*"([^"]*)"\s*\}'
)


def _tract(tract_id: str):
    return next(t for t in hp.CATALOGUE if t.id == tract_id)


def _row(tract):
    return (
        tract.id,
        tract.start,
        tract.end,
        tract.unit,
        tract.copies,
        tract.anchor,
        tract.in_himito,
        tract.simulate,
    )


def test_find_homopolymer_runs_length_at_least_6_on_rcrs():
    _, seq = hp.load_reference(REF)
    runs = hp.find_homopolymer_runs(seq, min_len=6)
    assert len(runs) == 35
    for expected in _RUNS_AT_LEAST_6:
        assert expected in runs
    for start, end, base, length in runs:
        assert length >= 6
        assert end - start + 1 == length
        assert seq[start - 1 : end] == base * length


def test_find_homopolymer_runs_default_min_len_is_4():
    assert inspect.signature(hp.find_homopolymer_runs).parameters["min_len"].default == 4
    _, seq = hp.load_reference(REF)
    runs = hp.find_homopolymer_runs(seq)
    assert len(runs) > 35
    assert all(length >= 4 for _, _, _, length in runs)
    for expected in _RUNS_AT_LEAST_6:
        assert expected in runs


def test_catalogue_matches_rcrs_and_plan_table():
    _, seq = hp.load_reference(REF)
    assert [_row(t) for t in hp.CATALOGUE] == _EXPECTED
    for tract in hp.CATALOGUE:
        assert seq[tract.start - 1 : tract.end] == tract.unit * tract.copies
        assert tract.anchor == tract.start - 1
        assert tract.end - tract.start + 1 == len(tract.unit) * tract.copies


def test_apply_background_hvs1_anchors():
    _, seq = hp.load_reference(REF)
    hvs1 = _tract("HVS1")

    both = hp.apply_background(seq, ["T16189C", "A16183C"])
    assert len(both) == len(seq)
    assert seq[16182] == "A" and both[16182] == "C"  # 16183 A>C
    assert seq[16188] == "T" and both[16188] == "C"  # 16189 T>C
    runs_both = hp.find_homopolymer_runs(both, min_len=6)
    assert (16183, 16193, "C", 11) in runs_both
    eff_both = hp.effective_tract(hvs1, both)
    assert (eff_both.start, eff_both.end, eff_both.anchor) == (16183, 16193, 16182)
    assert (hvs1.start, hvs1.end, hvs1.anchor) == (16184, 16188, 16183)

    only = hp.apply_background(seq, ["T16189C"])
    assert seq[16182] == "A" and only[16182] == "A"
    assert only[16188] == "C"
    runs_only = hp.find_homopolymer_runs(only, min_len=6)
    assert (16184, 16193, "C", 10) in runs_only
    eff_only = hp.effective_tract(hvs1, only)
    assert (eff_only.start, eff_only.end, eff_only.anchor) == (16184, 16193, 16183)


def test_effective_tract_extends_only_over_run_base():
    _, seq = hp.load_reference(REF)
    hvs1 = hp.effective_tract(_tract("HVS1"), seq)
    assert (hvs1.start, hvs1.end, hvs1.anchor, hvs1.copies) == (16184, 16188, 16183, 5)
    d310 = hp.effective_tract(_tract("D310"), seq)
    assert (d310.start, d310.end, d310.anchor, d310.copies) == (303, 309, 302, 7)
    # Period-2 unit: do not grow by the flanking single C's at 524-525.
    ca = hp.effective_tract(_tract("CA514"), seq)
    assert (ca.start, ca.end, ca.anchor, ca.copies) == (514, 523, 513, 5)


def test_rust_known_tracts_parity():
    _, seq = hp.load_reference(REF)
    text = open(TRACTS_RS).read()
    rust = {
        name: (int(start), int(end), ref_seq)
        for name, start, end, ref_seq in _TRACT_RE.findall(text)
    }
    by_id = {t.id: t for t in hp.CATALOGUE}
    for rust_name, cat_id in _RUST_TO_CATALOGUE.items():
        start, end, ref_seq = rust[rust_name]
        assert seq[start - 1 : end] == ref_seq
        tract = by_id[cat_id]
        assert tract.in_himito is True
    rust_names = set(rust)
    for tract in hp.CATALOGUE:
        if not tract.in_himito:
            assert tract.id not in rust_names


def test_load_reference_reexported_from_simulate_tree():
    assert st.load_reference is hp.load_reference
