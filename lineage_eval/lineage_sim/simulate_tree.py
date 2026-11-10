#!/usr/bin/env python3
"""Simulate a clonal mitochondrial mutation tree (SCITE-style truth).

Root = rCRS, plus any fixed ``--background`` homoplasmies. Each edge
introduces one heteroplasmic variant: an SNV, or, when
``--hp-lineage-indels`` is set, a homopolymer length change on a catalogue
tract. ``--hp-unstable`` additionally splits every clone, including ref, into
sub-clones by tract length; those alleles are not tree edges. Each node
carries a cumulative heteroplasmy frequency (fraction of molecules whose
lineage passes through it), and that determines how many reads the clone
contributes.
"""
import argparse
import itertools
import os
import random
import re
import sys
from dataclasses import dataclass, field

# Plan 2026-09-27 T1. Shared with the tract catalogue; sibling import matches
# how this file is run from lineage_eval/lineage_sim/.
from hp_tracts import CATALOGUE, apply_background, effective_tract, load_reference

# Regions to avoid when placing SNVs: rCRS control-region homopolymers / the
# 3107 'N' placeholder, where Himito's indel-prone calls would confound a
# clean SNV benchmark.
AVOID_RANGES = [(1, 40), (295, 320), (3105, 3110), (16180, 16200), (16560, 16569)]

# Fraction of an internal node's mass held back for its own terminal clone,
# rather than passed down to its children.
#
# This is the identifiability knob for *mutation order*. On a unary path the
# parent's own clone is the only thing that produces "parent without child"
# reads, at mass `keep`, while the child sits at `1 - keep`. Sequencing dropout
# at rate `fn` fakes "child without parent" reads at mass `(1 - keep) * fn`, so
# the true orientation only outvotes the artefact when
#
#     keep / (1 - keep)  >>  fn
#
# At keep=0.10 that ratio is 0.111, which does not clear a realistic ONT dropout
# of ~0.10 — the two signals cancel and the order of a 3-mutation chain becomes a
# coin flip (observed on the n=10 ont-r10 eval: 8 vs 9 reads, sign-test p ~ 0.5,
# with intermediate clones near 1% of molecules).
#
# Raising `keep` costs depth, though: children sit at `1 - keep` of the parent, so
# a larger keep drives the deepest mutation of a big tree toward the caller's
# detection floor. Measured worst-case truth HF at n=10 over seeds 1-5:
#
#     keep   order margin   min truth HF
#     0.10       1.11x          0.037     <- order is a coin flip
#     0.20       2.50x          0.020     <- resolvable, still 2x the 0.01 floor
#     0.30       4.29x          0.010     <- at the floor; costs variant recall
#
# 0.20 is the balance: it takes the order from unresolvable to resolvable while
# keeping every truth variant at twice the detection limit.
DEFAULT_INTERNAL_KEEP = 0.20


# Plan 2026-09-27 T2. 1-based VCF POS: the changed base for a SNV, the anchor
# base for an indel. ref/alt include that anchor, matching Himito's VCF alleles.
@dataclass(frozen=True)
class Edit:
    pos: int
    ref: str
    alt: str

    def is_indel(self) -> bool:
        return len(self.ref) != len(self.alt)


@dataclass(init=False)
class Node:
    """Mutation-tree node. ``pos``/``ref``/``alt`` mirror ``edit`` for callers
    that still pass the SNV kwargs."""

    id: int
    variant: str | None
    parent: int
    children: list[int] = field(default_factory=list)
    cum_freq: float = 0.0
    edit: Edit | None = None

    def __init__(
        self,
        id: int,
        variant: str | None,
        parent: int,
        pos: int | None = None,
        ref: str | None = None,
        alt: str | None = None,
        children: list[int] | None = None,
        cum_freq: float = 0.0,
        edit: Edit | None = None,
        tract_id: str | None = None,
        hp_delta: int | None = None,
    ) -> None:
        self.id = id
        self.variant = variant
        self.parent = parent
        self.children = [] if children is None else list(children)
        self.cum_freq = cum_freq
        # Plan 2026-09-27 T3. None on SNVs and on ROOT.
        self.tract_id = tract_id
        self.hp_delta = hp_delta
        if edit is not None:
            self.edit = edit
        elif pos is None:
            self.edit = None
        else:
            self.edit = Edit(pos, ref, alt)

    @property
    def pos(self) -> int | None:
        return None if self.edit is None else self.edit.pos

    @pos.setter
    def pos(self, value: int | None) -> None:
        if value is None:
            self.edit = None
            return
        ref = None if self.edit is None else self.edit.ref
        alt = None if self.edit is None else self.edit.alt
        self.edit = Edit(value, ref, alt)

    @property
    def ref(self) -> str | None:
        return None if self.edit is None else self.edit.ref

    @ref.setter
    def ref(self, value: str | None) -> None:
        if self.edit is None:
            return
        self.edit = Edit(self.edit.pos, value, self.edit.alt)

    @property
    def alt(self) -> str | None:
        return None if self.edit is None else self.edit.alt

    @alt.setter
    def alt(self, value: str | None) -> None:
        if self.edit is None:
            return
        self.edit = Edit(self.edit.pos, self.edit.ref, value)


@dataclass
class Tree:
    nodes: list[Node]


def himito_name(edit: Edit) -> str:
    """``m.{pos}{ref}>{alt}`` in 1-based VCF coordinates.

    Matches ``call::generate_variant_name`` once a SNP's 0-based ``pos+1`` and
    an indel's already-1-based anchor are both expressed as VCF POS.
    """
    return f"m.{edit.pos}{edit.ref}>{edit.alt}"


def variant_id(pos: int, ref: str, alt: str) -> str:
    return himito_name(Edit(pos, ref, alt))


def _coerce_edit(item: Edit | tuple[int, str, str]) -> Edit:
    if isinstance(item, Edit):
        return item
    pos, ref, alt = item
    return Edit(pos, ref, alt)


def apply_edits(seq: str, edits: list[Edit | tuple[int, str, str]]) -> str:
    """Apply SNVs and indels, highest anchor first.

    An upstream insertion then cannot shift a later coordinate. The result
    matches applying the edits in ascending position order while shifting each
    later edit by the net length change of the edits already applied.
    """
    ordered = sorted((_coerce_edit(e) for e in edits), key=lambda e: e.pos, reverse=True)
    g = list(seq)
    for edit in ordered:
        i = edit.pos - 1
        got = "".join(g[i:i + len(edit.ref)])
        assert got == edit.ref, f"ref mismatch at {edit.pos}: {got} != {edit.ref}"
        g[i:i + len(edit.ref)] = list(edit.alt)
    return "".join(g)


def left_align(edit: Edit, seq: str) -> Edit:
    """Python twin of ``call::left_align_variant`` + ``normalize_left``.

    Walk the anchor left while the inserted or deleted base equals the
    reference base at the anchor, then rebuild both alleles from the reference
    at the new anchor. SNVs are unchanged. ``edit.pos`` is already 1-based.
    """
    if (
        edit.pos < 1
        or not edit.ref
        or not edit.alt
        or len(edit.ref) == len(edit.alt)
    ):
        return edit

    if len(edit.ref) == 1 and len(edit.alt) > 1:
        a = edit.pos - 1
        inserted = list(edit.alt[1:])
        while a > 0 and inserted and seq[a] == inserted[-1]:
            inserted = [inserted[-1]] + inserted[:-1]
            a -= 1
        if a >= len(seq):
            return edit
        anchor = seq[a]
        return Edit(a + 1, anchor, anchor + "".join(inserted))

    if len(edit.alt) == 1 and len(edit.ref) > 1:
        n = len(edit.ref) - 1
        a = edit.pos - 1
        if a + n >= len(seq):
            return edit
        while a > 0 and a + n < len(seq) and seq[a] == seq[a + n]:
            a -= 1
        if a + n >= len(seq):
            return edit
        anchor = seq[a]
        return Edit(a + 1, seq[a:a + 1 + n], anchor)

    return edit


def _is_avoided(pos: int) -> bool:
    return any(lo <= pos <= hi for lo, hi in AVOID_RANGES)


def pick_snv(seq: str, pos: int, rng: random.Random) -> tuple[str, str]:
    ref = seq[pos - 1]
    alt = rng.choice([b for b in "ACGT" if b != ref])
    return ref, alt


def _sample_position(seq: str, used: set[int], rng: random.Random) -> int:
    while True:
        pos = rng.randint(1, len(seq))
        if pos in used or _is_avoided(pos) or seq[pos - 1] == "N":
            continue
        return pos


def _drawable_tracts():
    # Plan 2026-09-27 T3. simulate=True is the drawable set (everything but R8281).
    return [t for t in CATALOGUE if t.simulate]


def _draw_hp_delta(tract, hp_rng: random.Random) -> int:
    """Signed repeat-unit change. Homopolymers +1/+2/-1; CA514 ±1 unit."""
    # Plan 2026-09-27 T3.
    if len(tract.unit) == 1:
        roll = hp_rng.random()
        if roll < 0.6:
            return 1
        if roll < 0.8:
            return 2
        return -1
    if tract.id == "CA514":
        return 1 if hp_rng.random() < 0.5 else -1
    raise ValueError(f"no length distribution for tract {tract.id}")


def _raw_length_edit(seq: str, tract, delta: int) -> Edit:
    """Anchor-based length change of ``abs(delta)`` repeat units, before left-align."""
    anchor = tract.anchor
    base = seq[anchor - 1]
    n_units = abs(delta)
    if delta > 0:
        return Edit(anchor, base, base + tract.unit * n_units)
    n_bp = n_units * len(tract.unit)
    deleted = seq[anchor:anchor + n_bp]
    if deleted != tract.unit * n_units:
        raise ValueError(
            f"{tract.id}: expected {tract.unit * n_units!r} after {anchor}, got {deleted!r}"
        )
    return Edit(anchor, base + deleted, base)


def _tract_blocks_snv(pos: int, tract) -> bool:
    """True when ``pos`` is inside the run or is the anchor the indel consumes."""
    return tract.anchor <= pos <= tract.end


def _apply_hp_lineage_indels(seq: str, nodes: list[Node], k: int,
                             hp_rng: random.Random | None,
                             background_positions: set[int] | None = None) -> None:
    """Replace K SNV edges with distinct catalogue length changes.

    Draws come only from ``hp_rng``. Called after SNV placement so ``rng`` is
    untouched. A surviving SNV that landed in a chosen tract (or on its
    anchor) is redrawn from ``hp_rng``.
    """
    # Plan 2026-09-27 T3.
    if hp_rng is None:
        raise ValueError("hp_rng is required when hp_lineage_indels > 0")
    n_mutations = len(nodes) - 1
    # Plan 2026-09-27 T5. Anchors are the run on this root, which may already
    # include fixed background homoplasmies.
    drawable = [effective_tract(t, seq) for t in _drawable_tracts()]
    if k > n_mutations or k > len(drawable):
        raise ValueError(
            f"hp_lineage_indels={k} exceeds n_mutations={n_mutations} "
            f"or the {len(drawable)} simulate=True tracts"
        )
    node_ids = hp_rng.sample(range(1, n_mutations + 1), k)
    tracts = hp_rng.sample(drawable, k)
    chosen = []
    hp_ids = set(node_ids)
    for nid, tract in zip(node_ids, tracts):
        delta = _draw_hp_delta(tract, hp_rng)
        edit = left_align(_raw_length_edit(seq, tract, delta), seq)
        node = nodes[nid]
        node.edit = edit
        node.variant = himito_name(edit)
        node.tract_id = tract.id
        node.hp_delta = delta
        chosen.append(tract)

    def blocked(pos: int) -> bool:
        return any(_tract_blocks_snv(pos, tract) for tract in chosen)

    used: set[int] = set()
    # Plan 2026-09-27 T5. Same rule as the SNV loop: only a non-empty set is
    # occupied, and the redraw stays on hp_rng.
    if background_positions:
        used.update(background_positions)
    moving: list[int] = []
    for i in range(1, n_mutations + 1):
        if i in hp_ids:
            continue
        if blocked(nodes[i].pos):
            moving.append(i)
        else:
            used.add(nodes[i].pos)
    for i in moving:
        while True:
            pos = hp_rng.randint(1, len(seq))
            if pos in used or _is_avoided(pos) or seq[pos - 1] == "N" or blocked(pos):
                continue
            break
        ref, alt = pick_snv(seq, pos, hp_rng)
        nodes[i].edit = Edit(pos, ref, alt)
        nodes[i].variant = himito_name(nodes[i].edit)
        nodes[i].tract_id = None
        nodes[i].hp_delta = None
        used.add(pos)


TOPOLOGIES = ("random", "chain", "star")


def build_tree(seq: str, n_mutations: int, rng: random.Random,
               topology: str = "random", hp_lineage_indels: int = 0,
               hp_rng: random.Random | None = None,
               background_positions: set[int] | None = None) -> Tree:
    """Rooted mutation tree in one of three shapes.

    ``random`` (default) attaches each new mutation to a uniformly chosen existing
    node, giving a mix of chains and branch points. ``chain`` makes one unary path
    (every mutation is the child of the previous one); ``star`` hangs every
    mutation directly off ROOT.

    Shape is the knob that varies *ordering* difficulty independently of
    *detection* difficulty, which `--n-mutations` cannot do -- raising the mutation
    count shrinks every clone, drops the minimum HF and just costs variant recall.
    Ordering is only ever hard along a unary path, because that is where the
    evidence is a handful of "parent without child" reads (see
    `DEFAULT_INTERNAL_KEEP`); a star has no ancestral pairs to get wrong at all.

    A chain is also the *gentlest* shape on frequency: with one child per node no
    mass is split between siblings, so cum_freq decays only by `internal_keep` per
    level. At n=10 the deepest chain clone sits at 0.85 * 0.8^10 = 0.091, above the
    0.077 a random tree produced at the same settings. So `--topology chain` raises
    ordering difficulty while *improving* detectability -- the decoupling that
    raising --n-mutations failed to deliver.

    Only the ``random`` branch draws from `rng` for parent choice, so existing
    seeds reproduce byte-identically.

    ``hp_lineage_indels`` (plan 2026-09-27 T3) then replaces that many edges
    with catalogue length changes drawn from ``hp_rng``. Zero leaves the SNV
    tree alone and does not draw from ``hp_rng``.
    """
    if topology not in TOPOLOGIES:
        raise ValueError(f"unknown topology {topology!r}; expected one of {TOPOLOGIES}")
    root = Node(id=0, variant=None, pos=None, ref=None, alt=None, parent=0)
    nodes = [root]
    used: set[int] = set()
    # Plan 2026-09-27 T5. Occupied only when the caller passed a non-empty
    # set, so an omitted background does not touch this RNG.
    if background_positions:
        used.update(background_positions)
    for i in range(1, n_mutations + 1):
        if topology == "chain":
            parent = i - 1          # extend the single path
        elif topology == "star":
            parent = 0              # every mutation hangs off ROOT
        else:
            parent = rng.randrange(len(nodes))
        pos = _sample_position(seq, used, rng)
        used.add(pos)
        ref, alt = pick_snv(seq, pos, rng)
        node = Node(id=i, variant=variant_id(pos, ref, alt), pos=pos,
                    ref=ref, alt=alt, parent=parent)
        nodes.append(node)
        nodes[parent].children.append(i)
    # Plan 2026-09-27 T3. After the SNV loop so `rng` draws stay put.
    if hp_lineage_indels:
        _apply_hp_lineage_indels(
            seq, nodes, hp_lineage_indels, hp_rng,
            background_positions=background_positions,
        )
    return Tree(nodes=nodes)


# Plan 2026-09-27 T5. Preset order is the order written to truth/background.txt.
_BACKGROUND_TOKEN = re.compile(r"^([ACGT])(\d+)([ACGT])$")
_BACKGROUND_PRESETS = {
    "hvs1-polyC": ["T16189C", "A16183C"],
}


def parse_background(spec: str) -> list[str]:
    """``none`` applies nothing. ``hvs1-polyC`` is ``T16189C,A16183C``.

    Any other spec is a comma-separated list of ``<REF><pos><ALT>`` tokens.
    ``apply_background`` checks each REF base against rCRS.
    """
    # Plan 2026-09-27 T5.
    if spec == "none":
        return []
    if spec in _BACKGROUND_PRESETS:
        return list(_BACKGROUND_PRESETS[spec])
    tokens = spec.split(",")
    if not tokens or any(_BACKGROUND_TOKEN.fullmatch(token) is None for token in tokens):
        raise ValueError(
            f"background spec {spec!r} is not none, hvs1-polyC, "
            "or comma-separated <REF><pos><ALT> tokens"
        )
    return tokens


def _background_positions(tokens: list[str]) -> set[int]:
    return {int(_BACKGROUND_TOKEN.fullmatch(token).group(2)) for token in tokens}


def _background_himito_name(token: str) -> str:
    """Himito SNV name of a background token, on rCRS (the token's own REF)."""
    match = _BACKGROUND_TOKEN.fullmatch(token)
    ref, pos, alt = match.group(1), match.group(2), match.group(3)
    return f"m.{pos}{ref}>{alt}"


def deepest_chain_freq(n_mutations: int, ref_fraction: float,
                       internal_keep: float) -> float:
    """cum_freq of the last node of a `--topology chain` tree.

    Closed form: a unary path never splits mass between siblings, so each level
    just multiplies by (1 - internal_keep). Used to fail fast, because chain
    frequencies are fully deterministic -- retrying with another seed cannot help.
    """
    return (1.0 - ref_fraction) * (1.0 - internal_keep) ** n_mutations


def assign_frequencies(
    tree: Tree,
    ref_fraction: float,
    min_hf: float,
    max_hf: float,
    rng: random.Random,
    internal_keep: float = DEFAULT_INTERNAL_KEEP,
) -> None:
    """Assign each node a *stay* fraction, then propagate cumulative frequency
    top-down. The root reserves `ref_fraction` of molecules as mutation-free
    reference reads (so no mutation reaches 100% -> everything stays < max_hf).
    Every node keeps `internal_keep` of its incoming frequency for its own clone,
    so every mutation ends strictly inside (0.01, 0.95) and — more importantly —
    every internal node keeps enough terminal reads for its position in the chain
    to be recoverable. See `DEFAULT_INTERNAL_KEEP` for why that floor matters."""
    # Root gets total mass 1.0; it keeps ref_fraction for the reference clone.
    tree.nodes[0].cum_freq = 1.0

    def distribute(node_id: int, incoming: float) -> None:
        node = tree.nodes[node_id]
        node.cum_freq = incoming
        kids = node.children
        # Mass this node passes down vs. keeps for its own terminal clone.
        reserve = ref_fraction if node_id == 0 else 0.0
        available = incoming * (1.0 - reserve)
        if not kids:
            return
        # Give each child a share; hold back `internal_keep` of `available` for
        # this node's own clone so internal nodes still have terminal reads.
        share_pool = available * (1.0 - internal_keep)
        weights = [rng.uniform(0.5, 1.0) for _ in kids]
        wsum = sum(weights)
        for k, w in zip(kids, weights):
            distribute(k, share_pool * w / wsum)

    distribute(0, 1.0)
    # Root's cum_freq is definitional (1.0); mutations must be in-band.
    for node in tree.nodes[1:]:
        if not (min_hf < node.cum_freq < max_hf):
            raise ValueError(
                f"node {node.id} has cum_freq={node.cum_freq:.6f} outside "
                f"({min_hf}, {max_hf}); try a smaller --n-mutations, a lower "
                "--internal-keep, or a different --seed"
            )


def path_to_root(tree: Tree, node_id: int) -> list[Edit]:
    """Edits from root down to ``node_id`` (root excluded)."""
    chain = []
    cur = node_id
    while cur != 0:
        n = tree.nodes[cur]
        chain.append(n.edit)
        cur = n.parent
    chain.reverse()
    return chain


def clone_genome(seq: str, path: list[Edit | tuple[int, str, str]]) -> str:
    """Apply a root-to-node path. Tuples stay accepted for SNV callers."""
    return apply_edits(seq, path)


# Plan 2026-09-27 T4. Deltas are relative to the catalogue run. "common" is the
# only preset; an empty string is not a spec (the CLI omits the flag).
_COMMON_HP_UNSTABLE: dict[str, dict[int, float]] = {
    "D310": {0: 0.55, 1: 0.35, 2: 0.10},
    "HVS1": {0: 0.6, 1: 0.3, 2: 0.1},
    "C573": {0: 0.7, 1: 0.2, 2: 0.1},
}


def parse_hp_unstable(spec: str) -> dict[str, dict[int, float]]:
    """Parse ``common`` or ``TRACT:delta=prob,delta=prob;TRACT:...``.

    Every tract id must be in ``CATALOGUE``. Each distribution must sum to
    1 within 1e-9. Deltas are integers.
    """
    # Plan 2026-09-27 T4.
    if spec == "common":
        return {tid: dict(dist) for tid, dist in _COMMON_HP_UNSTABLE.items()}
    if not spec:
        raise ValueError("empty --hp-unstable spec; omit the flag instead")
    known = {t.id for t in CATALOGUE}
    out: dict[str, dict[int, float]] = {}
    for chunk in spec.split(";"):
        if ":" not in chunk:
            raise ValueError(f"hp-unstable term {chunk!r} is not TRACT:delta=prob")
        tid, body = chunk.split(":", 1)
        if tid not in known:
            raise ValueError(f"hp-unstable tract {tid!r} is not in CATALOGUE")
        if tid in out:
            raise ValueError(f"hp-unstable tract {tid!r} is repeated")
        if not body:
            raise ValueError(f"hp-unstable tract {tid} has no distribution")
        dist: dict[int, float] = {}
        for pair in body.split(","):
            if "=" not in pair:
                raise ValueError(f"hp-unstable pair {pair!r} is not delta=prob")
            delta_s, prob_s = pair.split("=", 1)
            try:
                delta = int(delta_s)
                prob = float(prob_s)
            except ValueError as exc:
                raise ValueError(f"hp-unstable pair {pair!r} is not delta=prob") from exc
            if prob < 0 or prob != prob:
                raise ValueError(f"hp-unstable probability {prob} for {tid} is negative")
            if delta in dist:
                raise ValueError(f"hp-unstable delta {delta} repeated for {tid}")
            dist[delta] = prob
        total = sum(dist.values())
        if abs(total - 1.0) > 1e-9:
            raise ValueError(
                f"hp-unstable distribution for {tid} sums to {total}, not 1"
            )
        out[tid] = dist
    if not out:
        raise ValueError("empty --hp-unstable spec; omit the flag instead")
    return out


@dataclass
class SubClone:
    """One terminal clone, or one length-state of a terminal clone."""

    id: str
    frequency: float
    hp_state: dict = field(default_factory=dict)
    node_variant: str = ""
    variant_path: str = ""
    edits: list = field(default_factory=list)
    lineage_tracts: frozenset = field(default_factory=frozenset)


def _clone_field(clone, name: str, default=None):
    if isinstance(clone, dict):
        return clone.get(name, default)
    return getattr(clone, name, default)


def _hp_product(spec: dict[str, dict[int, float]]):
    """Cartesian product. Tracts stay in spec order; deltas ascend."""
    tracts = list(spec)
    delta_lists = [sorted(spec[tid]) for tid in tracts]
    for index, deltas in enumerate(itertools.product(*delta_lists)):
        state = {tid: delta for tid, delta in zip(tracts, deltas)}
        prob = 1.0
        for tid, delta in zip(tracts, deltas):
            prob *= spec[tid][delta]
        yield index, state, prob


def _expand_one(clone, spec, min_freq: float, depth_drop=None):
    """Split one clone. ``depth_drop(raw_freq)`` removes combos before renormalising.

    Returns ``(subclones, depth_dropped_mass)``. Mass removed only by
    ``min_freq`` is not included in the dropped total. If every combination
    is removed, the modal state (lowest product index on a tie) is kept at
    the parent's frequency so that mass is not lost.
    """
    # Plan 2026-09-27 T4.
    parent_id = str(_clone_field(clone, "id"))
    parent_freq = float(_clone_field(clone, "frequency"))
    node_variant = _clone_field(clone, "node_variant", "") or ""
    variant_path = _clone_field(clone, "variant_path", "") or ""
    edits = list(_clone_field(clone, "edits", []) or [])
    lineage = _clone_field(clone, "lineage_tracts", ()) or ()
    raw = []
    for index, state, prob in _hp_product(spec):
        raw.append((index, dict(state), parent_freq * prob))
    depth_dropped = 0.0
    kept = []
    for index, state, freq in raw:
        if depth_drop is not None and depth_drop(freq):
            depth_dropped += freq
            continue
        if freq < min_freq:
            continue
        kept.append((index, state, freq))
    if not kept and raw and parent_freq > 0:
        best = max(raw, key=lambda item: (item[2], -item[0]))
        kept = [best]
        if depth_drop is not None and depth_drop(best[2]):
            depth_dropped -= best[2]
    kept_sum = sum(freq for _, _, freq in kept)
    if kept_sum <= 0:
        return [], depth_dropped
    scale = parent_freq / kept_sum
    out = []
    for index, state, freq in kept:
        out.append(SubClone(
            id=f"{parent_id}.h{index}",
            frequency=freq * scale,
            hp_state=state,
            node_variant=node_variant,
            variant_path=variant_path,
            edits=list(edits),
            lineage_tracts=frozenset(lineage),
        ))
    return out, depth_dropped


def expand_unstable(clones, spec, min_freq: float) -> list:
    """Split each clone across tract-length states. Depth is not consulted.

    Combinations with frequency < ``min_freq`` are dropped, then the survivors
    are rescaled so the parent clone's mass is unchanged. Output ids are
    ``{clone_id}.h{i}`` with ``i`` the cartesian-product index.
    """
    # Plan 2026-09-27 T4.
    out = []
    for clone in clones:
        subs, _dropped = _expand_one(clone, spec, min_freq, depth_drop=None)
        out.extend(subs)
    return out


def scheduled_read_count(frequencies, total_depth: int) -> int:
    """Reads ``simulate_reads.sh`` would schedule, including the per-rotation floor.

    Frequencies are formatted to 6 decimals first, matching ``clones.tsv``.
    """
    # Plan 2026-09-27 T4.
    total = 0
    for freq in frequencies:
        shown = float(f"{freq:.6f}")
        per_rotation = max(1, int(round(total_depth * shown / 2)))
        total += per_rotation * 2
    return total


def _per_rotation_rounds_to_zero(freq: float, total_depth: int) -> bool:
    return round(total_depth * freq / 2) == 0


def resolve_hp_subclones(parents, spec, min_freq: float, total_depth: int | None):
    """Expand, and apply the depth floor only when it would inflate reads by >5%.

    Returns ``(subclones, dropped_mass_or_None)``. ``dropped_mass`` is the raw
    frequency removed by the round-to-zero rule before renormalising. ``None``
    means that rule was not applied. ``expand_unstable`` itself stays depth-free.
    """
    # Plan 2026-09-27 T4.
    subs = expand_unstable(parents, spec, min_freq)
    if total_depth is None:
        return subs, None
    no_hp = scheduled_read_count((p.frequency for p in parents), total_depth)
    with_hp = scheduled_read_count((s.frequency for s in subs), total_depth)
    if not with_hp > no_hp * 1.05:
        return subs, None
    kept = []
    dropped = 0.0

    def drop(freq, depth=total_depth):
        return _per_rotation_rounds_to_zero(freq, depth)

    for parent in parents:
        part, mass = _expand_one(parent, spec, min_freq, depth_drop=drop)
        kept.extend(part)
        dropped += mass
    return kept, dropped


def format_hp_state(hp_state: dict, tract_order) -> str:
    """``D310:+1;C573:0``. Positive deltas carry a sign; zero does not."""
    parts = []
    for tid in tract_order:
        delta = hp_state[tid]
        text = f"+{delta}" if delta > 0 else str(int(delta))
        parts.append(f"{tid}:{text}")
    return ";".join(parts)


def apply_unstable_deltas(seq: str, path: list, hp_state: dict) -> str:
    """Path edits first, then each non-zero δ at that tract's anchor on the result.

    Higher anchors are applied first so an insertion does not move a later
    (leftward) coordinate. Anchors come from ``effective_tract`` on ``seq``.
    A zero delta adds no edit.
    """
    # Plan 2026-09-27 T4. T5: the anchor is the run on this root sequence.
    genome = apply_edits(seq, path) if path else seq
    by_id = {t.id: effective_tract(t, seq) for t in CATALOGUE}
    pending = [(by_id[tid], delta) for tid, delta in hp_state.items() if delta]
    pending.sort(key=lambda item: item[0].anchor, reverse=True)
    for tract, delta in pending:
        shift = 0
        for edit in path:
            if edit.pos < tract.anchor:
                shift += len(edit.alt) - len(edit.ref)
        anchor = tract.anchor + shift
        base = genome[anchor - 1]
        n_units = abs(delta)
        if delta > 0:
            raw = Edit(anchor, base, base + tract.unit * n_units)
        else:
            n_bp = n_units * len(tract.unit)
            deleted = genome[anchor:anchor + n_bp]
            if deleted != tract.unit * n_units:
                raise ValueError(
                    f"{tract.id}: expected {tract.unit * n_units!r} after {anchor}, "
                    f"got {deleted!r}"
                )
            raw = Edit(anchor, base + deleted, base)
        genome = apply_edits(genome, [left_align(raw, genome)])
    return genome


def _parents_from_tree(tree: Tree, seq: str) -> list:
    parents = [SubClone(
        id="ref",
        frequency=_terminal_freq(tree, 0),
        node_variant="ROOT",
        variant_path="",
        edits=[],
        lineage_tracts=frozenset(),
    )]
    for n in tree.nodes[1:]:
        path = path_to_root(tree, n.id)
        vpath = ",".join(himito_name(left_align(e, seq)) for e in path)
        tracts = []
        cur = n.id
        while cur != 0:
            node = tree.nodes[cur]
            if node.tract_id:
                tracts.append(node.tract_id)
            cur = node.parent
        parents.append(SubClone(
            id=str(n.id),
            frequency=_terminal_freq(tree, n.id),
            node_variant=n.variant or "",
            variant_path=vpath,
            edits=list(path),
            lineage_tracts=frozenset(tracts),
        ))
    return parents


def _write_hp_tracts(path: str, subclones, spec, seq: str) -> None:
    """One row per (tract, δ ≠ 0, origin). ``origin`` is ``mixed`` when the
    lineage path already carries an indel in that tract."""
    # Plan 2026-09-27 T4.
    by_id = {t.id: t for t in CATALOGUE}
    mass: dict[tuple, float] = {}
    for sub in subclones:
        for tid, delta in sub.hp_state.items():
            if delta == 0:
                continue
            origin = "mixed" if tid in sub.lineage_tracts else "unstable"
            key = (tid, delta, origin)
            mass[key] = mass.get(key, 0.0) + sub.frequency
    with open(path, "w") as fh:
        fh.write(
            "tract\tanchor\tref_len\tbackground_len\tdelta\ttruth_freq\t"
            "canonical_name\torigin\n"
        )
        for tid in spec:
            tract = by_id[tid]
            # Plan 2026-09-27 T5. ref_len stays the catalogue run. background_len
            # and anchor follow the run on this root sequence.
            eff = effective_tract(tract, seq)
            ref_len = tract.copies * len(tract.unit)
            background_len = eff.copies * len(eff.unit)
            for delta in sorted(d for d in spec[tid] if d != 0):
                name = himito_name(left_align(_raw_length_edit(seq, eff, delta), seq))
                for origin in ("unstable", "mixed"):
                    freq = mass.get((tid, delta, origin), 0.0)
                    if freq <= 0:
                        continue
                    fh.write(
                        f"{tid}\t{eff.anchor}\t{ref_len}\t{background_len}\t{delta:+d}\t"
                        f"{freq:.6f}\t{name}\t{origin}\n"
                    )


def write_truth(tree: Tree, seq: str, outdir: str, *, hp_spec=None,
                min_subclone_freq: float = 0.002,
                total_depth: int | None = None,
                background: list[str] | None = None) -> None:
    tdir = os.path.join(outdir, "truth")
    os.makedirs(tdir, exist_ok=True)

    # Plan 2026-09-27 T5. Omitted when the background list is empty, so a
    # no-background run does not gain a file.
    if background:
        with open(os.path.join(tdir, "background.txt"), "w") as fh:
            for token in background:
                fh.write(_background_himito_name(token) + "\n")

    with open(os.path.join(tdir, "truth_mutation_tree.tsv"), "w") as fh:
        fh.write("node_id\tvariant\tparent_id\tparent_variant\n")
        for n in tree.nodes:
            var = "ROOT" if n.variant is None else n.variant
            pvar = "ROOT" if tree.nodes[n.parent].variant is None else tree.nodes[n.parent].variant
            fh.write(f"{n.id}\t{var}\t{n.parent}\t{pvar}\n")

    with open(os.path.join(tdir, "truth_variants.txt"), "w") as fh:
        for n in tree.nodes[1:]:
            fh.write(n.variant + "\n")

    # Clones: the reference clone + one terminal clone per node that keeps mass.
    # Frequency must be *terminal* mass (molecules that stop at this node), not
    # cum_freq — otherwise ancestral mass is double-counted into parent clones
    # and read simulation depth blows up (especially costly for HiFi/ccs).
    # Plan 2026-09-27 T4. The hp_state column and sub-clone ids exist only when
    # --hp-unstable is set, so a run without the flag stays byte-identical.
    with open(os.path.join(tdir, "clones.tsv"), "w") as fh, \
         open(os.path.join(tdir, "clone_genomes.fa"), "w") as fa:
        if hp_spec is None:
            fh.write("clone_id\tnode_variant\tvariant_path\tfrequency\n")
            # reference clone (mutation-free)
            fh.write(f"ref\tROOT\t\t{_terminal_freq(tree, 0):.6f}\n")
            fa.write(">clone_ref\n")
            _write_fasta_seq(fa, seq)
            for n in tree.nodes[1:]:
                path = path_to_root(tree, n.id)
                genome = clone_genome(seq, path)
                # Left-aligned Himito names. Canonicalisation of HP alleles is T6;
                # until then this is the name `Himito call` writes for the same edit.
                vpath = ",".join(himito_name(left_align(e, seq)) for e in path)
                fh.write(f"{n.id}\t{n.variant}\t{vpath}\t{_terminal_freq(tree, n.id):.6f}\n")
                fa.write(f">clone_{n.id}\n")
                _write_fasta_seq(fa, genome)
        else:
            parents = _parents_from_tree(tree, seq)
            subclones, dropped = resolve_hp_subclones(
                parents, hp_spec, min_subclone_freq, total_depth,
            )
            if dropped is not None:
                print(
                    f"hp-unstable: dropped {dropped:.6f} mass with per-rotation "
                    "depth rounding to 0 before renormalising",
                    file=sys.stderr,
                )
            fh.write("clone_id\tnode_variant\tvariant_path\tfrequency\thp_state\n")
            order = list(hp_spec)
            for sub in subclones:
                state = format_hp_state(sub.hp_state, order)
                fh.write(
                    f"{sub.id}\t{sub.node_variant}\t{sub.variant_path}\t"
                    f"{sub.frequency:.6f}\t{state}\n"
                )
                fa.write(f">clone_{sub.id}\n")
                _write_fasta_seq(fa, apply_unstable_deltas(seq, sub.edits, sub.hp_state))
            _write_hp_tracts(
                os.path.join(tdir, "hp_tracts.tsv"), subclones, hp_spec, seq,
            )

    # Plan 2026-09-27 T3. Omitted when no lineage HP indel was placed (K=0).
    hp_nodes = [n for n in tree.nodes if n.tract_id]
    if hp_nodes:
        with open(os.path.join(tdir, "hp_lineage.tsv"), "w") as fh:
            fh.write("node_id\tvariant\ttract\tdelta\n")
            for n in hp_nodes:
                fh.write(
                    f"{n.id}\t{n.variant}\t{n.tract_id}\t{n.hp_delta:+d}\n"
                )


def _terminal_freq(tree: Tree, node_id: int) -> float:
    """Mass kept at `node_id` (not passed to children). Sums to 1 over all nodes."""
    node = tree.nodes[node_id]
    child_sum = sum(tree.nodes[c].cum_freq for c in node.children)
    return max(0.0, node.cum_freq - child_sum)


def _write_fasta_seq(fh, seq: str, width: int = 60) -> None:
    for i in range(0, len(seq), width):
        fh.write(seq[i:i + width] + "\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    _repo = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    ap.add_argument("--reference", default=os.path.join(_repo, "rCRS.fasta"))
    ap.add_argument("--n-mutations", type=int, default=10)
    ap.add_argument("--ref-fraction", type=float, default=0.15)
    ap.add_argument(
        "--internal-keep",
        type=float,
        default=DEFAULT_INTERNAL_KEEP,
        help="fraction of an internal node's mass kept for its own terminal "
        "clone; controls whether mutation order along a chain is recoverable "
        f"under sequencing dropout (default {DEFAULT_INTERNAL_KEEP})",
    )
    ap.add_argument(
        "--topology",
        choices=TOPOLOGIES,
        default="random",
        help="tree shape. 'random' (default) mixes chains and branch points; "
        "'chain' makes one unary path, the hardest shape for mutation *ordering* "
        "and the gentlest on clone frequency; 'star' hangs every mutation off ROOT, "
        "which has no ancestral pairs at all (a control, not a difficulty setting)",
    )
    ap.add_argument("--min-hf", type=float, default=0.01)
    ap.add_argument("--max-hf", type=float, default=0.99)
    ap.add_argument(
        "--hp-lineage-indels",
        type=int,
        default=0,
        help="replace K mutation edges with homopolymer length changes "
        "drawn from a separate RNG (plan 2026-09-27 T3)",
    )
    ap.add_argument(
        "--hp-unstable",
        default=None,
        help="lineage-independent tract-length distribution: 'common' or "
        "TRACT:delta=prob,delta=prob;TRACT:... (plan 2026-09-27 T4)",
    )
    ap.add_argument(
        "--hp-min-subclone-freq",
        type=float,
        default=0.002,
        help="drop unstable sub-clones below this frequency and renormalise "
        "within the parent clone (default 0.002)",
    )
    ap.add_argument(
        "--total-depth",
        type=int,
        default=None,
        help="with --hp-unstable, drop sub-clones whose per-rotation depth "
        "rounds to 0 when that floor would inflate scheduled reads by more "
        "than 5 percent (plan 2026-09-27 T4)",
    )
    ap.add_argument(
        "--background",
        default="none",
        help="fixed homoplasmies applied to the root before the tree: "
        "none (default), hvs1-polyC, or comma-separated <REF><pos><ALT> "
        "tokens such as T16189C,A16183C (plan 2026-09-27 T5)",
    )
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    _, seq = load_reference(args.reference)
    # Plan 2026-09-27 T5. Apply before the tree and before unstable deltas.
    # An empty list leaves the sequence and the SNV RNG alone.
    background = parse_background(args.background)
    if background:
        root_seq = apply_background(seq, background)
        occupied = _background_positions(background)
    else:
        root_seq = seq
        occupied = None
    # A chain's frequencies are deterministic (no sibling split, so no rng in the
    # mass assignment), which makes the retry loop below pointless for it: 100
    # reseeds all produce the same cum_freq. Check the closed form up front and
    # say what would actually help.
    if args.topology == "chain":
        deepest = deepest_chain_freq(args.n_mutations, args.ref_fraction,
                                     args.internal_keep)
        if not deepest > args.min_hf:
            import math
            feasible = math.floor(
                math.log(args.min_hf / (1.0 - args.ref_fraction))
                / math.log(1.0 - args.internal_keep)
            )
            raise SystemExit(
                f"--topology chain with --n-mutations {args.n_mutations} puts the "
                f"deepest clone at cum_freq={deepest:.6f}, at or below --min-hf "
                f"{args.min_hf}. Chain frequencies are deterministic, so a different "
                f"--seed cannot help. At these --ref-fraction/--internal-keep the "
                f"limit is --n-mutations {feasible}; otherwise lower --min-hf or "
                f"--internal-keep."
            )

    max_attempts = 100
    hp_k = args.hp_lineage_indels
    for attempt in range(max_attempts):
        seed = args.seed + attempt
        rng = random.Random(seed)
        # Plan 2026-09-27 T3. Not constructed at all when K is 0.
        hp_rng = random.Random(seed + 1_000_003) if hp_k else None
        tree = build_tree(
            root_seq, args.n_mutations, rng, topology=args.topology,
            hp_lineage_indels=hp_k, hp_rng=hp_rng,
            background_positions=occupied,
        )
        try:
            assign_frequencies(
                tree,
                args.ref_fraction,
                args.min_hf,
                args.max_hf,
                rng,
                internal_keep=args.internal_keep,
            )
            break
        except ValueError as exc:
            if attempt == max_attempts - 1:
                raise SystemExit(
                    f"Could not find valid frequency assignment after {max_attempts} attempts. "
                    "Try a smaller --n-mutations or a different --seed."
                ) from exc
    # Plan 2026-09-27 T4. --total-depth is consulted only with --hp-unstable.
    hp_spec = parse_hp_unstable(args.hp_unstable) if args.hp_unstable is not None else None
    write_truth(
        tree, root_seq, args.outdir,
        hp_spec=hp_spec,
        min_subclone_freq=args.hp_min_subclone_freq,
        total_depth=args.total_depth if hp_spec is not None else None,
        background=background or None,
    )
    print(f"wrote truth for {args.n_mutations} mutations "
          f"(topology={args.topology}) to {args.outdir}/truth")


if __name__ == "__main__":
    main()
