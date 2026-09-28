# Plan 2026-09-27 T1
"""rCRS homopolymer tract catalogue and background-aware run anchors.

Coordinates are 1-based inclusive. ``anchor`` is the base immediately before
the run (the position Himito uses for a left-aligned insertion).
"""
import re
from dataclasses import dataclass

_BACKGROUND_TOKEN = re.compile(r"^([ACGT])(\d+)([ACGT])$")


@dataclass(frozen=True)
class Tract:
    id: str
    start: int
    end: int
    unit: str
    copies: int
    anchor: int
    in_himito: bool
    simulate: bool


# Run intervals, not the wider windows in src/tracts.rs. in_himito is True
# only for the three tracts already in KNOWN_TRACTS (HVS1-C maps to HVS1).
CATALOGUE: list[Tract] = [
    Tract("D310", 303, 309, "C", 7, 302, True, True),
    Tract("HVS1", 16184, 16188, "C", 5, 16183, True, True),
    Tract("C573", 568, 573, "C", 6, 567, False, True),
    Tract("T455", 452, 455, "T", 4, 451, False, True),
    Tract("C960", 956, 960, "C", 5, 955, False, True),
    Tract("C5899", 5895, 5899, "C", 5, 5894, False, True),
    Tract("A12418", 12418, 12425, "A", 8, 12417, False, True),
    Tract("CA514", 514, 523, "CA", 5, 513, True, True),
    Tract("R8281", 8272, 8289, "CCCCCTCTA", 2, 8271, False, False),
]


def load_reference(path: str) -> tuple[str, str]:
    """Return ``(name, sequence)`` from a FASTA file. Sequence is uppercased."""
    name, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                name = line[1:].split()[0]
            else:
                chunks.append(line.strip())
    return name, "".join(chunks).upper()


def find_homopolymer_runs(
    seq: str, min_len: int = 4
) -> list[tuple[int, int, str, int]]:
    """Homopolymer runs as ``(start, end, base, length)``, 1-based inclusive."""
    runs: list[tuple[int, int, str, int]] = []
    i = 0
    n = len(seq)
    while i < n:
        j = i + 1
        while j < n and seq[j] == seq[i]:
            j += 1
        length = j - i
        if length >= min_len:
            runs.append((i + 1, j, seq[i], length))
        i = j
    return runs


def apply_background(seq: str, variants: list[str]) -> str:
    """Apply ``<ref><pos><alt>`` SNVs (e.g. ``T16189C``) and return the sequence.

    Tokens are substitutions at 1-based positions. The reference base must
    already be ``ref``.
    """
    bases = list(seq)
    for token in variants:
        match = _BACKGROUND_TOKEN.fullmatch(token)
        if match is None:
            raise ValueError(f"background token {token!r} is not <ref><pos><alt>")
        ref, pos_s, alt = match.group(1), match.group(2), match.group(3)
        pos = int(pos_s)
        if pos < 1 or pos > len(bases):
            raise ValueError(f"background position {pos} is outside the sequence")
        found = bases[pos - 1]
        if found != ref:
            raise ValueError(f"{token}: expected {ref} at {pos}, found {found}")
        bases[pos - 1] = alt
    return "".join(bases)


def effective_tract(tract: Tract, background_seq: str) -> Tract:
    """Run boundaries and anchor of ``tract`` on ``background_seq``.

    Homopolymer units (length 1) extend left and right while the neighbouring
    base equals the run base. The anchor is the base before that extended run.
    Tandem units are left at their catalogue coordinates; step 1(d) specifies
    the homopolymer case.
    """
    if len(tract.unit) != 1:
        return tract
    base = tract.unit
    start, end = tract.start, tract.end
    while start > 1 and background_seq[start - 2] == base:
        start -= 1
    while end < len(background_seq) and background_seq[end] == base:
        end += 1
    return Tract(
        id=tract.id,
        start=start,
        end=end,
        unit=tract.unit,
        copies=end - start + 1,
        anchor=start - 1,
        in_himito=tract.in_himito,
        simulate=tract.simulate,
    )
