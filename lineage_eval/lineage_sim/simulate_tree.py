#!/usr/bin/env python3
"""Simulate a clonal mitochondrial mutation tree (SCITE-style truth).

Root = rCRS. Each edge introduces one heteroplasmic variant: an SNV, or, when
``--hp-lineage-indels`` is set, a homopolymer length change on a catalogue
tract. Each node carries a cumulative heteroplasmy frequency (fraction of
molecules whose lineage passes through it), and that determines how many
reads the clone contributes.
"""
import argparse
import os
import random
from dataclasses import dataclass, field

# Plan 2026-09-27 T1. Shared with the tract catalogue; sibling import matches
# how this file is run from lineage_eval/lineage_sim/.
from hp_tracts import CATALOGUE, load_reference

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
                             hp_rng: random.Random | None) -> None:
    """Replace K SNV edges with distinct catalogue length changes.

    Draws come only from ``hp_rng``. Called after SNV placement so ``rng`` is
    untouched. A surviving SNV that landed in a chosen tract (or on its
    anchor) is redrawn from ``hp_rng``.
    """
    # Plan 2026-09-27 T3.
    if hp_rng is None:
        raise ValueError("hp_rng is required when hp_lineage_indels > 0")
    n_mutations = len(nodes) - 1
    drawable = _drawable_tracts()
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
               hp_rng: random.Random | None = None) -> Tree:
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
        _apply_hp_lineage_indels(seq, nodes, hp_lineage_indels, hp_rng)
    return Tree(nodes=nodes)


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


def write_truth(tree: Tree, seq: str, outdir: str) -> None:
    tdir = os.path.join(outdir, "truth")
    os.makedirs(tdir, exist_ok=True)

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
    with open(os.path.join(tdir, "clones.tsv"), "w") as fh, \
         open(os.path.join(tdir, "clone_genomes.fa"), "w") as fa:
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
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    _, seq = load_reference(args.reference)
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
            seq, args.n_mutations, rng, topology=args.topology,
            hp_lineage_indels=hp_k, hp_rng=hp_rng,
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
    write_truth(tree, seq, args.outdir)
    print(f"wrote truth for {args.n_mutations} mutations "
          f"(topology={args.topology}) to {args.outdir}/truth")


if __name__ == "__main__":
    main()
