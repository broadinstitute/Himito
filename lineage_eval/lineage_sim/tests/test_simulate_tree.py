# eval/lineage_sim/tests/test_simulate_tree.py
import os
import random
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import simulate_tree as st

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
REF = os.path.join(_REPO, "rCRS.fasta")


def test_variant_id_format():
    assert st.variant_id(3144, "A", "G") == "m.3144A>G"


def test_load_reference_reads_chrm_16569():
    name, seq = st.load_reference(REF)
    assert name == "chrM"
    assert len(seq) == 16569
    assert set(seq) <= set("ACGTN")


def test_pick_snv_returns_alt_differing_from_ref():
    _, seq = st.load_reference(REF)
    rng = random.Random(0)
    ref, alt = st.pick_snv(seq, 3144, rng)
    assert ref == seq[3143]  # 1-based pos -> 0-based index
    assert alt in "ACGT" and alt != ref


def test_build_tree_is_connected_and_sized():
    _, seq = st.load_reference(REF)
    rng = random.Random(42)
    tree = st.build_tree(seq, n_mutations=8, rng=rng)
    # 8 mutation nodes + 1 root
    assert len(tree.nodes) == 9
    assert tree.nodes[0].variant is None  # root
    # every non-root node reaches root by following parents
    for n in tree.nodes[1:]:
        cur, steps = n.id, 0
        while tree.nodes[cur].parent != cur:
            cur = tree.nodes[cur].parent
            steps += 1
            assert steps <= len(tree.nodes)
        assert cur == 0
    # all mutation positions distinct
    positions = [n.pos for n in tree.nodes[1:]]
    assert len(positions) == len(set(positions))


def test_frequencies_stay_inside_heteroplasmy_band():
    _, seq = st.load_reference(REF)
    rng = random.Random(7)
    tree = st.build_tree(seq, n_mutations=10, rng=rng)
    st.assign_frequencies(tree, ref_fraction=0.15, min_hf=0.01, max_hf=0.95, rng=rng)
    for n in tree.nodes[1:]:  # every mutation
        assert 0.01 < n.cum_freq < 0.95, (n.id, n.cum_freq)


def test_assign_frequencies_meets_tight_min_hf_when_topology_allows():
    """n=10 / min_hf=0.1 is the WDL high-HF setting. Seed 5's first tree is
    topologically feasible (subtree minima fit in the root share pool), but
    drawing child weights uniformly often puts a leaf below 0.1 and the 100
    retry loop never recovers. Assignment must give each subtree at least its
    minimum mass rather than hoping random splits land in-band."""
    _, seq = st.load_reference(REF)
    rng = random.Random(5)
    tree = st.build_tree(seq, n_mutations=10, rng=rng)
    st.assign_frequencies(tree, ref_fraction=0.15, min_hf=0.1, max_hf=0.95, rng=rng)
    for n in tree.nodes[1:]:
        assert 0.1 < n.cum_freq < 0.95, (n.id, n.cum_freq)


def test_unary_chain_of_10_is_infeasible_at_min_hf_0_1():
    """A 10-step chain with keep=0.20 bottoms out at 0.85*0.8^10 ≈ 0.091, so
    no frequency split can put every mutation above 0.1. The allocator must
    reject the topology rather than emit an out-of-band tree."""
    _, seq = st.load_reference(REF)
    rng = random.Random(0)
    root = st.Node(id=0, variant=None, pos=None, ref=None, alt=None, parent=0)
    nodes = [root]
    used: set[int] = set()
    for i in range(1, 11):
        pos = st._sample_position(seq, used, rng)
        used.add(pos)
        ref, alt = st.pick_snv(seq, pos, rng)
        nodes.append(st.Node(
            id=i, variant=st.variant_id(pos, ref, alt), pos=pos,
            ref=ref, alt=alt, parent=i - 1,
        ))
        nodes[i - 1].children.append(i)
    try:
        st.assign_frequencies(st.Tree(nodes=nodes), 0.15, 0.1, 0.95, rng)
    except ValueError as exc:
        assert "topology cannot keep" in str(exc)
    else:
        raise AssertionError("unary n=10 at min_hf=0.1 should be infeasible")


def test_apply_edits_snv_plus_upstream_insertion_shifts_c_run():
    """Plan 2026-09-27 T2. T16189C and an A>AC insertion at 16183.

    Right-to-left application matches ascending order with later positions
    shifted by the net length change: C11 over 1-based 16184-16194, length 16570.
    """
    _, seq = st.load_reference(REF)
    edits = [st.Edit(16189, "T", "C"), st.Edit(16183, "A", "AC")]
    out = st.apply_edits(seq, edits)
    assert len(out) == 16570
    assert out[16183:16194] == "C" * 11

    bases = list(seq)
    delta = 0
    for edit in sorted(edits, key=lambda e: e.pos):
        i = edit.pos - 1 + delta
        assert "".join(bases[i:i + len(edit.ref)]) == edit.ref
        bases[i:i + len(edit.ref)] = list(edit.alt)
        delta += len(edit.alt) - len(edit.ref)
    assert out == "".join(bases)


def test_himito_name_snv_and_anchor_indels():
    """Plan 2026-09-27 T2. Names use 1-based VCF POS for SNVs and indel anchors."""
    assert st.himito_name(st.Edit(3144, "A", "G")) == "m.3144A>G"
    assert st.himito_name(st.Edit(302, "A", "AC")) == "m.302A>AC"
    assert st.himito_name(st.Edit(302, "ACC", "A")) == "m.302ACC>A"
    assert st.variant_id(3144, "A", "G") == st.himito_name(st.Edit(3144, "A", "G"))


def test_himito_name_matches_five_test_vcf_records():
    """Plan 2026-09-27 T2. test.vcf has no ID; names come from POS/REF/ALT."""
    vcf = pathlib.Path(_REPO) / "test.vcf"
    wanted = {
        (263, "A", "G"): "m.263A>G",
        (302, "A", "AC"): "m.302A>AC",
        (302, "A", "ACC"): "m.302A>ACC",
        (310, "T", "TC"): "m.310T>TC",
        (456, "C", "T"): "m.456C>T",
    }
    found = set()
    for line in vcf.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        _chrom, pos, _id, ref, alt, *_rest = line.split("\t")
        key = (int(pos), ref, alt)
        if key in wanted:
            found.add(key)
            assert st.himito_name(st.Edit(int(pos), ref, alt)) == wanted[key]
    assert found == set(wanted)


def test_left_align_poly_c_insertion_and_d310_deletion():
    """Plan 2026-09-27 T2. Python twin of call::left_align_variant."""
    from hp_tracts import apply_background

    _, seq = st.load_reference(REF)
    ins = st.Edit(16188, "C", "CC")
    aligned = st.left_align(ins, seq)
    assert aligned == st.Edit(16183, "A", "AC")
    assert st.himito_name(aligned) == "m.16183A>AC"

    bg = apply_background(seq, ["T16189C"])
    assert st.himito_name(st.left_align(ins, bg)) == "m.16183A>AC"

    bg2 = apply_background(seq, ["T16189C", "A16183C"])
    aligned2 = st.left_align(ins, bg2)
    assert aligned2 == st.Edit(16182, "A", "AC")
    assert st.himito_name(aligned2) == "m.16182A>AC"

    deletion = st.Edit(308, "CC", "C")
    assert st.himito_name(st.left_align(deletion, seq)) == "m.302AC>A"

    snv = st.Edit(3144, "A", "G")
    assert st.left_align(snv, seq) == snv


def test_node_edit_round_trip_and_path_to_root():
    """Plan 2026-09-27 T2. pos/ref/alt stay usable; the stored allele is an Edit."""
    root = st.Node(id=0, variant=None, pos=None, ref=None, alt=None, parent=0)
    assert root.edit is None
    node = st.Node(id=1, variant="m.3144A>G", pos=3144, ref="A", alt="G", parent=0)
    assert node.edit == st.Edit(3144, "A", "G")
    assert (node.pos, node.ref, node.alt) == (3144, "A", "G")
    node.alt = "T"
    assert node.edit == st.Edit(3144, "A", "T")
    root.children.append(1)
    tree = st.Tree(nodes=[root, node])
    assert st.path_to_root(tree, 1) == [st.Edit(3144, "A", "T")]


def test_clone_genome_accepts_edits_and_tuples():
    """Plan 2026-09-27 T2. Indels apply right-to-left; SNV tuples still work."""
    _, seq = st.load_reference(REF)
    edits = [st.Edit(16189, "T", "C"), st.Edit(16183, "A", "AC")]
    from_edits = st.clone_genome(seq, edits)
    from_tuples = st.clone_genome(seq, [(16189, "T", "C"), (16183, "A", "AC")])
    assert from_edits == from_tuples
    assert len(from_edits) == 16570
    assert from_edits[16183:16194] == "C" * 11


def test_write_truth_left_aligns_variant_path(tmp_path):
    """Plan 2026-09-27 T2. variant_path is the left-aligned Himito name."""
    _, seq = st.load_reference(REF)
    root = st.Node(id=0, variant=None, pos=None, ref=None, alt=None, parent=0)
    node = st.Node(
        id=1, variant="m.16188C>CC", pos=16188, ref="C", alt="CC", parent=0,
    )
    root.children.append(1)
    root.cum_freq = 1.0
    node.cum_freq = 0.5
    st.write_truth(st.Tree(nodes=[root, node]), seq, str(tmp_path))
    rows = (tmp_path / "truth" / "clones.tsv").read_text().splitlines()
    assert rows[0] == "clone_id\tnode_variant\tvariant_path\tfrequency"
    assert rows[2].split("\t")[2] == "m.16183A>AC"
    assert len(rows[2].split("\t")) == 4


def test_clone_genome_applies_only_path_snvs():
    _, seq = st.load_reference(REF)
    # a two-SNV path
    v1 = st.variant_id(100, seq[99], "A" if seq[99] != "A" else "C")
    alt1 = v1.split(">")[1]
    genome = st.clone_genome(seq, [(100, seq[99], alt1)])
    assert genome[99] == alt1
    assert genome[:99] == seq[:99]
    assert genome[100:] == seq[100:]
    assert len(genome) == len(seq)


def test_frequencies_decrease_from_parent_to_child():
    """cum_freq must strictly decrease along every edge (root excluded)."""
    _, seq = st.load_reference(REF)
    rng = random.Random(13)
    tree = st.build_tree(seq, n_mutations=10, rng=rng)
    st.assign_frequencies(tree, ref_fraction=0.15, min_hf=0.01, max_hf=0.95, rng=rng)
    for node in tree.nodes[1:]:  # skip root
        parent = tree.nodes[node.parent]
        assert node.cum_freq < parent.cum_freq, (
            f"node {node.id} cum_freq={node.cum_freq:.6f} >= "
            f"parent {parent.id} cum_freq={parent.cum_freq:.6f}"
        )


def test_internal_clones_stay_big_enough_for_order_to_be_recoverable():
    """Mutation *order* along a chain is only recoverable if each internal node
    keeps enough of its own terminal clone to outvote dropout artefacts.

    On a parent->child edge the parent's own clone is the sole source of genuine
    "parent without child" reads, while dropout at rate `fn` fabricates "child
    without parent" reads from the child's clone. At the old keep=0.10 the two
    were the same size for ONT (~0.10), which is what made n=10 chain order a
    coin flip. The binding case is a unary edge, where the ratio is exactly
    keep/(1-keep) = 0.25 against a 0.10 artefact, i.e. 2.5x.
    """
    ont_dropout = 0.10
    _, seq = st.load_reference(REF)
    for seed in (7, 13, 42, 99):
        rng = random.Random(seed)
        tree = st.build_tree(seq, n_mutations=10, rng=rng)
        st.assign_frequencies(tree, ref_fraction=0.15, min_hf=0.01, max_hf=0.95, rng=rng)
        for node in tree.nodes[1:]:
            parent = tree.nodes[node.parent]
            if parent.id == 0:
                continue  # root is the reference, not a mutation to order against
            own = parent.cum_freq - sum(tree.nodes[k].cum_freq for k in parent.children)
            assert own > 2 * ont_dropout * node.cum_freq, (
                f"seed {seed}: node {parent.id}'s own clone ({own:.4f}) is not "
                f"clearly above the dropout artefact from its child "
                f"({ont_dropout * node.cum_freq:.4f}) — order is not identifiable"
            )


def test_write_truth_emits_expected_files_and_headers(tmp_path):
    """write_truth produces all four output files with correct headers/content."""
    _, seq = st.load_reference(REF)
    n_mutations = 5
    rng = random.Random(99)
    tree = st.build_tree(seq, n_mutations=n_mutations, rng=rng)
    st.assign_frequencies(tree, ref_fraction=0.15, min_hf=0.01, max_hf=0.95, rng=rng)
    st.write_truth(tree, seq, str(tmp_path))

    truth_dir = tmp_path / "truth"

    # --- truth_mutation_tree.tsv ---
    tree_tsv = truth_dir / "truth_mutation_tree.tsv"
    assert tree_tsv.exists()
    lines = tree_tsv.read_text().splitlines()
    assert lines[0] == "node_id\tvariant\tparent_id\tparent_variant"
    variants_col = [line.split("\t")[1] for line in lines[1:]]
    assert "ROOT" in variants_col, "no ROOT row found in truth_mutation_tree.tsv"

    # --- truth_variants.txt ---
    variants_txt = truth_dir / "truth_variants.txt"
    assert variants_txt.exists()
    variant_lines = [l for l in variants_txt.read_text().splitlines() if l]
    assert len(variant_lines) == n_mutations
    import re
    pattern = re.compile(r"^m\.\d+[ACGT]>[ACGT]$")
    for v in variant_lines:
        assert pattern.match(v), f"variant line {v!r} does not match m.<pos><ref>><alt>"

    # --- clones.tsv ---
    clones_tsv = truth_dir / "clones.tsv"
    assert clones_tsv.exists()
    clone_lines = clones_tsv.read_text().splitlines()
    assert clone_lines[0] == "clone_id\tnode_variant\tvariant_path\tfrequency"
    clone_ids = [line.split("\t")[0] for line in clone_lines[1:]]
    assert "ref" in clone_ids, "no 'ref' clone row in clones.tsv"
    # Terminal clone frequencies must partition the molecule pool.
    freqs = [float(line.split("\t")[3]) for line in clone_lines[1:]]
    assert abs(sum(freqs) - 1.0) < 1e-5, f"clone frequencies sum to {sum(freqs)}, not 1"

    # --- clone_genomes.fa ---
    fa_file = truth_dir / "clone_genomes.fa"
    assert fa_file.exists()
    headers = [l for l in fa_file.read_text().splitlines() if l.startswith(">")]
    assert ">clone_ref" in headers
    assert len(headers) == n_mutations + 1, (
        f"expected {n_mutations + 1} FASTA records, got {len(headers)}"
    )
