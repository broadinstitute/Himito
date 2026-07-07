#!/usr/bin/env python3
"""Simulate a clonal mitochondrial mutation tree (SCITE-style truth).

Root = rCRS. Each edge introduces one heteroplasmic SNV. Each node carries a
cumulative heteroplasmy frequency (fraction of molecules whose lineage passes
through it), and that determines how many reads the clone contributes.
"""
import argparse
import os
import random
from dataclasses import dataclass, field

# Regions to avoid when placing SNVs: rCRS control-region homopolymers / the
# 3107 'N' placeholder, where Himito's indel-prone calls would confound a
# clean SNV benchmark.
AVOID_RANGES = [(1, 40), (295, 320), (3105, 3110), (16180, 16200), (16560, 16569)]


@dataclass
class Node:
    id: int
    variant: str | None
    pos: int | None
    ref: str | None
    alt: str | None
    parent: int
    children: list[int] = field(default_factory=list)
    cum_freq: float = 0.0


@dataclass
class Tree:
    nodes: list[Node]


def variant_id(pos: int, ref: str, alt: str) -> str:
    return f"m.{pos}{ref}>{alt}"


def load_reference(path: str) -> tuple[str, str]:
    name, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                name = line[1:].split()[0]
            else:
                chunks.append(line.strip())
    return name, "".join(chunks).upper()


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


def build_tree(seq: str, n_mutations: int, rng: random.Random) -> Tree:
    """Random rooted tree: each new mutation attaches to a uniformly chosen
    existing node (root or a prior mutation), giving a mix of chains and
    branch points."""
    root = Node(id=0, variant=None, pos=None, ref=None, alt=None, parent=0)
    nodes = [root]
    used: set[int] = set()
    for i in range(1, n_mutations + 1):
        parent = rng.randrange(len(nodes))
        pos = _sample_position(seq, used, rng)
        used.add(pos)
        ref, alt = pick_snv(seq, pos, rng)
        node = Node(id=i, variant=variant_id(pos, ref, alt), pos=pos,
                    ref=ref, alt=alt, parent=parent)
        nodes.append(node)
        nodes[parent].children.append(i)
    return Tree(nodes=nodes)


def assign_frequencies(tree: Tree, ref_fraction: float, rng: random.Random) -> None:
    """Assign each node a *stay* fraction, then propagate cumulative frequency
    top-down. The root reserves `ref_fraction` of molecules as mutation-free
    reference reads (so no mutation reaches 100% -> everything stays < max_hf).
    Every node keeps a floor of its incoming frequency for its own clone, so
    every mutation ends strictly inside (0.01, 0.95)."""
    n = len(tree.nodes)
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
        # Give each child a share; keep >= 10% of `available` for this node's
        # own clone so internal nodes still have terminal reads.
        keep = 0.10
        share_pool = available * (1.0 - keep)
        weights = [rng.uniform(0.5, 1.0) for _ in kids]
        wsum = sum(weights)
        for k, w in zip(kids, weights):
            distribute(k, share_pool * w / wsum)

    distribute(0, 1.0)
    # Root's cum_freq is definitional (1.0); mutations must be in-band.
    for node in tree.nodes[1:]:
        assert 0.01 < node.cum_freq < 0.95, (node.id, node.cum_freq)


def path_to_root(tree: Tree, node_id: int) -> list[tuple[int, str, str]]:
    """Ordered (pos, ref, alt) SNVs from root down to `node_id`."""
    chain = []
    cur = node_id
    while cur != 0:
        n = tree.nodes[cur]
        chain.append((n.pos, n.ref, n.alt))
        cur = n.parent
    chain.reverse()
    return chain


def clone_genome(seq: str, path: list[tuple[int, str, str]]) -> str:
    g = list(seq)
    for pos, ref, alt in path:
        assert g[pos - 1] == ref, f"ref mismatch at {pos}: {g[pos-1]} != {ref}"
        g[pos - 1] = alt
    return "".join(g)


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
    with open(os.path.join(tdir, "clones.tsv"), "w") as fh, \
         open(os.path.join(tdir, "clone_genomes.fa"), "w") as fa:
        fh.write("clone_id\tnode_variant\tvariant_path\tfrequency\tn_reads\n")
        # reference clone (mutation-free)
        fh.write(f"ref\tROOT\t\t{_ref_freq(tree):.6f}\t0\n")
        fa.write(">clone_ref\n")
        _write_fasta_seq(fa, seq)
        for n in tree.nodes[1:]:
            path = path_to_root(tree, n.id)
            genome = clone_genome(seq, path)
            vpath = ",".join(variant_id(p, r, a) for p, r, a in path)
            fh.write(f"{n.id}\t{n.variant}\t{vpath}\t{n.cum_freq:.6f}\t0\n")
            fa.write(f">clone_{n.id}\n")
            _write_fasta_seq(fa, genome)


def _ref_freq(tree: Tree) -> float:
    # The reference clone's share is whatever mass the root reserved.
    child_sum = sum(tree.nodes[c].cum_freq for c in tree.nodes[0].children)
    return max(0.0, 1.0 - child_sum)


def _write_fasta_seq(fh, seq: str, width: int = 60) -> None:
    for i in range(0, len(seq), width):
        fh.write(seq[i:i + width] + "\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference", default="/Users/suhang/Analysis/Himito/rCRS.fasta")
    ap.add_argument("--n-mutations", type=int, default=10)
    ap.add_argument("--ref-fraction", type=float, default=0.15)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    _, seq = load_reference(args.reference)
    tree = build_tree(seq, args.n_mutations, rng)
    assign_frequencies(tree, args.ref_fraction, rng)
    write_truth(tree, seq, args.outdir)
    print(f"wrote truth for {args.n_mutations} mutations to {args.outdir}/truth")


if __name__ == "__main__":
    main()
