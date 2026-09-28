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


# Plan 2026-09-27 T3. Frozen truth_mutation_tree.tsv for build_tree(n=10)
# before --hp-lineage-indels existed. Seeds 1..5.
_SNV_ONLY_TRUTH = {
    1: (
        "node_id\tvariant\tparent_id\tparent_variant\n"
        "0\tROOT\t0\tROOT\n"
        "1\tm.2068C>G\t0\tROOT\n"
        "2\tm.16235A>G\t0\tROOT\n"
        "3\tm.12440A>C\t1\tm.2068C>G\n"
        "4\tm.15987C>A\t0\tROOT\n"
        "5\tm.14181A>T\t3\tm.12440A>C\n"
        "6\tm.14595A>G\t0\tROOT\n"
        "7\tm.7497G>T\t5\tm.14181A>T\n"
        "8\tm.10402A>C\t1\tm.2068C>G\n"
        "9\tm.834G>T\t0\tROOT\n"
        "10\tm.12492A>T\t8\tm.10402A>C\n"
    ),
    2: (
        "node_id\tvariant\tparent_id\tparent_variant\n"
        "0\tROOT\t0\tROOT\n"
        "1\tm.3002G>A\t0\tROOT\n"
        "2\tm.5541C>T\t1\tm.3002G>A\n"
        "3\tm.10098G>C\t2\tm.5541C>T\n"
        "4\tm.1171A>T\t1\tm.3002G>A\n"
        "5\tm.14113T>G\t1\tm.3002G>A\n"
        "6\tm.12192G>T\t3\tm.10098G>C\n"
        "7\tm.16452T>C\t3\tm.10098G>C\n"
        "8\tm.900C>G\t0\tROOT\n"
        "9\tm.10436C>G\t7\tm.16452T>C\n"
        "10\tm.5390A>T\t6\tm.12192G>T\n"
    ),
    3: (
        "node_id\tvariant\tparent_id\tparent_variant\n"
        "0\tROOT\t0\tROOT\n"
        "1\tm.4274T>C\t0\tROOT\n"
        "2\tm.2148A>T\t1\tm.4274T>C\n"
        "3\tm.15376A>G\t0\tROOT\n"
        "4\tm.6284A>T\t1\tm.4274T>C\n"
        "5\tm.15610A>G\t3\tm.15376A>G\n"
        "6\tm.4936C>A\t5\tm.15610A>G\n"
        "7\tm.4969G>T\t5\tm.15610A>G\n"
        "8\tm.497C>T\t6\tm.4936C>A\n"
        "9\tm.5224G>T\t1\tm.4274T>C\n"
        "10\tm.9872A>C\t0\tROOT\n"
    ),
    4: (
        "node_id\tvariant\tparent_id\tparent_variant\n"
        "0\tROOT\t0\tROOT\n"
        "1\tm.9939G>A\t0\tROOT\n"
        "2\tm.15692A>C\t1\tm.9939G>A\n"
        "3\tm.2180A>C\t0\tROOT\n"
        "4\tm.9483T>A\t3\tm.2180A>C\n"
        "5\tm.11805A>G\t1\tm.9939G>A\n"
        "6\tm.3480A>G\t1\tm.9939G>A\n"
        "7\tm.841A>T\t1\tm.9939G>A\n"
        "8\tm.8905C>A\t4\tm.9483T>A\n"
        "9\tm.10154A>G\t2\tm.15692A>C\n"
        "10\tm.2842C>T\t5\tm.11805A>G\n"
    ),
    5: (
        "node_id\tvariant\tparent_id\tparent_variant\n"
        "0\tROOT\t0\tROOT\n"
        "1\tm.11749A>T\t0\tROOT\n"
        "2\tm.15258A>C\t0\tROOT\n"
        "3\tm.1700T>A\t2\tm.15258A>C\n"
        "4\tm.12183G>C\t0\tROOT\n"
        "5\tm.12477T>G\t1\tm.11749A>T\n"
        "6\tm.8171T>A\t0\tROOT\n"
        "7\tm.7102T>C\t5\tm.12477T>G\n"
        "8\tm.5967T>C\t4\tm.12183G>C\n"
        "9\tm.2358A>C\t2\tm.15258A>C\n"
        "10\tm.14577T>A\t9\tm.2358A>C\n"
    ),
}

_SIM_DIR = pathlib.Path(__file__).resolve().parents[1]


def _catalogue():
    from hp_tracts import CATALOGUE
    return CATALOGUE


def _drawable():
    return [t for t in _catalogue() if t.simulate]


def _tract_blocks_snv(pos, tract) -> bool:
    """Run interval plus the anchor base the indel edit consumes."""
    return tract.anchor <= pos <= tract.end


def _unit_delta(edit, tract) -> int:
    bp = len(edit.alt) - len(edit.ref)
    assert bp % len(tract.unit) == 0
    return bp // len(tract.unit)


def _expected_length_edit(seq, tract, delta):
    anchor = tract.anchor
    base = seq[anchor - 1]
    if delta > 0:
        raw = st.Edit(anchor, base, base + tract.unit * delta)
    else:
        n_bp = abs(delta) * len(tract.unit)
        deleted = seq[anchor:anchor + n_bp]
        raw = st.Edit(anchor, base + deleted, base)
    return st.left_align(raw, seq)


def test_build_tree_hp_lineage_indels_are_distinct_tract_edits():
    """Plan 2026-09-27 T3. K edges become length changes on distinct tracts."""
    _, seq = st.load_reference(REF)
    by_id = {t.id: t for t in _drawable()}
    assert "R8281" not in by_id
    tree = st.build_tree(
        seq, n_mutations=10, rng=random.Random(11),
        hp_lineage_indels=3, hp_rng=random.Random(99),
    )
    assert len(tree.nodes) == 11
    indels = [n for n in tree.nodes[1:] if n.edit.is_indel()]
    snvs = [n for n in tree.nodes[1:] if not n.edit.is_indel()]
    assert len(indels) == 3
    assert len(snvs) == 7
    assert len({n.tract_id for n in indels}) == 3
    for n in indels:
        tract = by_id[n.tract_id]
        assert tract.simulate
        delta = _unit_delta(n.edit, tract)
        if len(tract.unit) == 1:
            assert delta in (1, 2, -1)
        else:
            assert tract.id == "CA514"
            assert delta in (1, -1)
        assert n.edit == _expected_length_edit(seq, tract, delta)
        assert n.variant == st.himito_name(n.edit)
        assert seq[n.pos - 1:n.pos - 1 + len(n.ref)] == n.ref
        genome = st.clone_genome(seq, st.path_to_root(tree, n.id))
        assert isinstance(genome, str) and genome
    for n in snvs:
        assert n.tract_id is None
        assert not st._is_avoided(n.pos)
        assert seq[n.pos - 1] != "N"
        for indel in indels:
            assert not _tract_blocks_snv(n.pos, by_id[indel.tract_id])

    # Without replacement: K equal to the drawable set uses each tract once.
    k_all = len(by_id)
    for seed in range(5):
        full = st.build_tree(
            seq, n_mutations=k_all, rng=random.Random(seed),
            hp_lineage_indels=k_all, hp_rng=random.Random(1_000 + seed),
        )
        got = [n.tract_id for n in full.nodes[1:] if n.edit.is_indel()]
        assert sorted(got) == sorted(by_id)


def test_d310_length_spellings_are_left_aligned_on_rcrs():
    """Plan 2026-09-27 T3. +1 m.302A>AC, -1 m.302AC>A, +2 m.302A>ACC."""
    _, seq = st.load_reference(REF)
    d310 = next(t for t in _catalogue() if t.id == "D310")
    assert st.himito_name(_expected_length_edit(seq, d310, 1)) == "m.302A>AC"
    assert st.himito_name(_expected_length_edit(seq, d310, -1)) == "m.302AC>A"
    assert st.himito_name(_expected_length_edit(seq, d310, 2)) == "m.302A>ACC"
    # The tree must actually be able to emit each of those three.
    seen = set()
    for seed in range(400):
        tree = st.build_tree(
            seq, n_mutations=8, rng=random.Random(seed),
            hp_lineage_indels=1, hp_rng=random.Random(50_000 + seed),
        )
        node = next(n for n in tree.nodes[1:] if n.edit.is_indel())
        if node.tract_id == "D310":
            seen.add(node.variant)
    assert {"m.302A>AC", "m.302AC>A", "m.302A>ACC"} <= seen


def test_hp_lineage_snvs_avoid_ranges_and_chosen_tracts():
    """Plan 2026-09-27 T3. Surviving SNVs stay out of AVOID_RANGES and HP tracts.

    A SNV that the SNV-only placement put inside a chosen tract (or on its
    anchor) is redrawn from hp_rng. Every other SNV keeps the SNV-only allele,
    which is how K=0 stays byte-identical: hp_rng is never mixed into rng.
    """
    _, seq = st.load_reference(REF)
    by_id = {t.id: t for t in _drawable()}
    moved = 0
    for seed in range(400):
        rng_hp = random.Random(seed)
        tree = st.build_tree(
            seq, n_mutations=40, rng=rng_hp,
            hp_lineage_indels=6, hp_rng=random.Random(seed + 1_000_003),
        )
        rng_plain = random.Random(seed)
        plain = st.build_tree(seq, n_mutations=40, rng=rng_plain)
        assert [n.parent for n in tree.nodes] == [n.parent for n in plain.nodes]
        chosen = {
            n.tract_id: by_id[n.tract_id]
            for n in tree.nodes[1:] if n.edit.is_indel()
        }
        assert len(chosen) == 6
        snv_pos = []
        for node, orig in zip(tree.nodes[1:], plain.nodes[1:]):
            if node.edit.is_indel():
                continue
            assert not node.edit.is_indel()
            assert not st._is_avoided(node.pos)
            assert seq[node.pos - 1] != "N"
            assert node.tract_id is None
            for tract in chosen.values():
                assert not _tract_blocks_snv(node.pos, tract)
            hit = any(_tract_blocks_snv(orig.pos, t) for t in chosen.values())
            if hit:
                assert (node.pos, node.ref, node.alt) != (orig.pos, orig.ref, orig.alt)
                moved += 1
            else:
                assert (node.pos, node.ref, node.alt, node.variant) == (
                    orig.pos, orig.ref, orig.alt, orig.variant,
                )
            snv_pos.append(node.pos)
        assert len(snv_pos) == len(set(snv_pos))
    assert moved > 0


def test_hp_delta_distribution_matches_catalogue():
    """Plan 2026-09-27 T3. HP tracts {+1: 0.6, +2: 0.2, -1: 0.2}; CA514 ±1 unit."""
    _, seq = st.load_reference(REF)
    by_id = {t.id: t for t in _drawable()}
    hp_counts = {1: 0, 2: 0, -1: 0}
    ca_counts = {1: 0, -1: 0}
    for seed in range(600):
        tree = st.build_tree(
            seq, n_mutations=4, rng=random.Random(seed),
            hp_lineage_indels=1, hp_rng=random.Random(70_000 + seed),
        )
        node = next(n for n in tree.nodes[1:] if n.edit.is_indel())
        tract = by_id[node.tract_id]
        delta = _unit_delta(node.edit, tract)
        if tract.id == "CA514":
            assert delta in ca_counts
            assert abs(len(node.edit.alt) - len(node.edit.ref)) == 2
            ca_counts[delta] += 1
        else:
            assert len(tract.unit) == 1
            assert delta in hp_counts
            hp_counts[delta] += 1
    hp_n = sum(hp_counts.values())
    assert hp_n > 400
    assert 0.50 <= hp_counts[1] / hp_n <= 0.70
    assert 0.10 <= hp_counts[2] / hp_n <= 0.32
    assert 0.10 <= hp_counts[-1] / hp_n <= 0.32
    assert ca_counts[1] > 0 and ca_counts[-1] > 0
    ca_n = sum(ca_counts.values())
    assert abs(ca_counts[1] / ca_n - 0.5) < 0.25


def test_hp_lineage_indels_zero_matches_snv_only_truth(tmp_path):
    """Plan 2026-09-27 T3. K=0 does not draw from hp_rng and matches today's trees."""
    _, seq = st.load_reference(REF)
    for seed, expected in _SNV_ONLY_TRUTH.items():
        rng = random.Random(seed)
        hp_rng = random.Random(seed + 1_000_003)
        hp_state = hp_rng.getstate()
        tree = st.build_tree(
            seq, n_mutations=10, rng=rng,
            hp_lineage_indels=0, hp_rng=hp_rng,
        )
        assert hp_rng.getstate() == hp_state
        st.write_truth(tree, seq, str(tmp_path / f"k0_{seed}"))
        got = (tmp_path / f"k0_{seed}" / "truth" / "truth_mutation_tree.tsv").read_text()
        assert got == expected
        assert not (tmp_path / f"k0_{seed}" / "truth" / "hp_lineage.tsv").exists()

        rng_default = random.Random(seed)
        default = st.build_tree(seq, n_mutations=10, rng=rng_default)
        st.write_truth(default, seq, str(tmp_path / f"default_{seed}"))
        default_txt = (
            tmp_path / f"default_{seed}" / "truth" / "truth_mutation_tree.tsv"
        ).read_text()
        assert default_txt == expected
        # The SNV generator is consumed exactly as far as the default call.
        rng_again = random.Random(seed)
        st.build_tree(seq, n_mutations=10, rng=rng_again, hp_lineage_indels=0)
        assert rng_again.getstate() == rng_default.getstate()


def test_hvs1_plus1_with_background_snvs_named_on_background():
    """Plan 2026-09-27 T3. Direct call, not the CLI.

    An HVS1 +1 C together with T16189C and A16183C is named by left_align on
    apply_background output. Default write_truth still left-aligns on plain
    rCRS (no --background until T5), so the same edit stays m.16183A>AC there.
    """
    from hp_tracts import apply_background

    _, seq = st.load_reference(REF)
    ins = st.Edit(16183, "A", "AC")
    assert st.himito_name(st.left_align(ins, seq)) == "m.16183A>AC"
    bg = apply_background(seq, ["T16189C", "A16183C"])
    named = st.himito_name(st.left_align(ins, bg))
    assert named == "m.16182A>AC"


def test_write_truth_hp_lineage_tsv_one_row_per_indel(tmp_path):
    """Plan 2026-09-27 T3. truth/hp_lineage.tsv is node_id, variant, tract, delta."""
    _, seq = st.load_reference(REF)
    by_id = {t.id: t for t in _drawable()}
    tree = st.build_tree(
        seq, n_mutations=10, rng=random.Random(3),
        hp_lineage_indels=3, hp_rng=random.Random(4),
    )
    st.write_truth(tree, seq, str(tmp_path))
    text = (tmp_path / "truth" / "hp_lineage.tsv").read_text()
    rows = text.splitlines()
    assert rows[0] == "node_id\tvariant\ttract\tdelta"
    assert len(rows) == 4
    indels = [n for n in tree.nodes[1:] if n.edit.is_indel()]
    assert [int(r.split("\t")[0]) for r in rows[1:]] == [n.id for n in indels]
    for line, node in zip(rows[1:], indels):
        node_id, variant, tract_id, delta_s = line.split("\t")
        tract = by_id[tract_id]
        delta = _unit_delta(node.edit, tract)
        assert node_id == str(node.id)
        assert variant == node.variant == st.himito_name(node.edit)
        assert tract_id == node.tract_id
        assert delta_s == f"{delta:+d}"
        assert delta_s in {"+1", "+2", "-1"}
    # clones.tsv gains no column.
    header = (tmp_path / "truth" / "clones.tsv").read_text().splitlines()[0]
    assert header == "clone_id\tnode_variant\tvariant_path\tfrequency"


def test_cli_hp_lineage_indels_writes_hp_lineage_tsv(tmp_path):
    """Plan 2026-09-27 T3. --hp-lineage-indels K; K=0 omits hp_lineage.tsv."""
    import subprocess
    script = _SIM_DIR / "simulate_tree.py"
    py = sys.executable
    zero = tmp_path / "zero"
    r0 = subprocess.run(
        [py, str(script), "--reference", REF, "--n-mutations", "8",
         "--seed", "1", "--hp-lineage-indels", "0", "--outdir", str(zero)],
        capture_output=True, text=True,
    )
    assert r0.returncode == 0, r0.stderr
    assert not (zero / "truth" / "hp_lineage.tsv").exists()
    plain = tmp_path / "plain"
    rp = subprocess.run(
        [py, str(script), "--reference", REF, "--n-mutations", "8",
         "--seed", "1", "--outdir", str(plain)],
        capture_output=True, text=True,
    )
    assert rp.returncode == 0, rp.stderr
    assert (zero / "truth" / "truth_mutation_tree.tsv").read_text() == (
        plain / "truth" / "truth_mutation_tree.tsv"
    ).read_text()

    out = tmp_path / "k2"
    r = subprocess.run(
        [py, str(script), "--reference", REF, "--n-mutations", "8",
         "--seed", "1", "--hp-lineage-indels", "2", "--outdir", str(out)],
        capture_output=True, text=True,
    )
    assert r.returncode == 0, r.stderr
    rows = (out / "truth" / "hp_lineage.tsv").read_text().splitlines()
    assert rows[0] == "node_id\tvariant\ttract\tdelta"
    assert len(rows) == 3


def test_run_eval_forwards_hp_lineage_indels_not_rates(tmp_path):
    """Plan 2026-09-27 T3. Empty omits the flag; rates are not for simulate_tree."""
    import subprocess
    bindir = tmp_path / "bin"
    bindir.mkdir()
    log = tmp_path / "argv.txt"
    py = bindir / "python"
    py.write_text(
        "#!/bin/bash\n"
        f'printf "%s\\n" "$*" >> "{log}"\n'
        "exit 42\n"
    )
    py.chmod(0o755)
    env = os.environ.copy()
    env["PATH"] = str(bindir) + os.pathsep + env.get("PATH", "")
    r = subprocess.run(
        ["bash", str(_SIM_DIR / "run_eval.sh"),
         "--outdir", str(tmp_path / "out"),
         "--hp-lineage-indels", "3",
         "--hp-fp-rate", "0.2",
         "--hp-fn-rate", "0.4",
         "--n-mutations", "6", "--seed", "1",
         "--ref", REF],
        capture_output=True, text=True, env=env,
    )
    assert r.returncode != 0
    argv = log.read_text()
    assert "simulate_tree.py" in argv
    assert "--hp-lineage-indels 3" in argv
    assert "--hp-fp-rate" not in argv
    assert "--hp-fn-rate" not in argv


def test_run_himito_forwards_hp_rates_only_to_lineage(tmp_path):
    """Plan 2026-09-27 T3. --hp-fp-rate / --hp-fn-rate reach Himito lineage only."""
    import subprocess
    bindir = tmp_path / "bin"
    bindir.mkdir()
    log = tmp_path / "calls.txt"
    minimap = bindir / "minimap2"
    minimap.write_text("#!/bin/bash\nexit 0\n")
    samtools = bindir / "samtools"
    samtools.write_text(
        "#!/bin/bash\n"
        f'echo "samtools $*" >> "{log}"\n'
        "if [[ \"$1\" == \"sort\" ]]; then\n"
        "  prev=\n"
        "  for a in \"$@\"; do\n"
        "    if [[ \"$prev\" == \"-o\" ]]; then : > \"$a\"; fi\n"
        "    prev=\"$a\"\n"
        "  done\n"
        "fi\n"
        "exit 0\n"
    )
    himito = tmp_path / "Himito"
    himito.write_text(
        "#!/bin/bash\n"
        f'echo "himito $*" >> "{log}"\n'
        "exit 0\n"
    )
    for path in (minimap, samtools, himito):
        path.chmod(0o755)
    out = tmp_path / "cell"
    (out / "reads").mkdir(parents=True)
    (out / "reads" / "reads.fastq.gz").write_bytes(b"@r\nACGT\n+\n!!!!\n")
    env = os.environ.copy()
    env["PATH"] = str(bindir) + os.pathsep + env.get("PATH", "")
    r = subprocess.run(
        ["bash", str(_SIM_DIR / "run_himito.sh"),
         "--outdir", str(out), "--profile", "hifi",
         "--hp-fp-rate", "0.25", "--hp-fn-rate", "0.5",
         "--himito", str(himito), "--ref", REF],
        capture_output=True, text=True, env=env,
    )
    assert r.returncode == 0, r.stderr + r.stdout
    himito_lines = [
        ln for ln in log.read_text().splitlines() if ln.startswith("himito ")
    ]
    call = [ln for ln in himito_lines if ln.split()[1] == "call"]
    lineage = [ln for ln in himito_lines if ln.split()[1] == "lineage"]
    build = [ln for ln in himito_lines if ln.split()[1] == "build"]
    assert len(call) == 1 and len(lineage) == 1 and len(build) == 1
    assert "--hp-fp-rate" not in call[0] and "--hp-fn-rate" not in call[0]
    assert "--hp-fp-rate" not in build[0] and "--hp-fn-rate" not in build[0]
    assert "--hp-fp-rate 0.25" in lineage[0]
    assert "--hp-fn-rate 0.5" in lineage[0]


def test_sweep_seeds_forwards_hp_lineage_indels():
    """Plan 2026-09-27 T3. sweep_seeds.sh accepts the flag and forwards it."""
    import re
    import subprocess
    script = (_SIM_DIR / "sweep_seeds.sh").read_text()
    assert re.search(r"--hp-lineage-indels\)", script)
    assert re.search(
        r'EVAL_ARGS\+=\(--hp-lineage-indels "\$HP_LINEAGE_INDELS"\)', script,
    )
    r = subprocess.run(
        ["bash", str(_SIM_DIR / "sweep_seeds.sh"), "--hp-lineage-indels", "2"],
        capture_output=True, text=True,
    )
    assert "unknown arg" not in r.stderr
    assert "usage:" in r.stderr
    assert r.returncode != 0


# Plan 2026-09-27 T4. sha256 of clones.tsv and clone_genomes.fa from
# simulate_tree.py --n-mutations 10 --seed 1 at HEAD 8f649cd (no --hp-unstable).
_NO_HP_CLONES_SHA256 = "f1032170b13619779c2d39aaaf389c6aa0fee3ab11aa6d16c16be019a7ca900a"
_NO_HP_GENOMES_SHA256 = "7d066368b703d783dd3fb8e492fd8c21a54b2bdc45b2c8981a89bc4d3a4b1c39"


class _FreqClone:
    def __init__(self, id, frequency):
        self.id = id
        self.frequency = frequency


def _fasta_records(path):
    records = {}
    header, chunks = None, []
    for line in pathlib.Path(path).read_text().splitlines():
        if line.startswith(">"):
            if header is not None:
                records[header] = "".join(chunks)
            header, chunks = line[1:], []
        else:
            chunks.append(line.strip())
    if header is not None:
        records[header] = "".join(chunks)
    return records


def _scheduled_reads(freq_fields, total_depth):
    """Per clone, per rotation max(1, round(total_depth * freq / 2)), times 2.

    ``freq_fields`` are the clones.tsv strings simulate_reads.sh parses.
    """
    total = 0
    for field in freq_fields:
        per_rotation = max(1, int(round(total_depth * float(field) / 2)))
        total += per_rotation * 2
    return total


def _clone_freq_fields(clones_tsv):
    lines = pathlib.Path(clones_tsv).read_text().splitlines()[1:]
    return [line.split("\t")[3] for line in lines if line]


def test_parse_hp_unstable_common_and_explicit():
    """Plan 2026-09-27 T4. Preset and TRACT:delta=prob lists; sums are 1."""
    common = st.parse_hp_unstable("common")
    assert common == {
        "D310": {0: 0.55, 1: 0.35, 2: 0.10},
        "HVS1": {0: 0.6, 1: 0.3, 2: 0.1},
        "C573": {0: 0.7, 1: 0.2, 2: 0.1},
    }
    explicit = st.parse_hp_unstable("D310:0=0.5,1=0.5;C573:0=0.8,-1=0.2")
    assert explicit == {
        "D310": {0: 0.5, 1: 0.5},
        "C573": {0: 0.8, -1: 0.2},
    }
    catalogue_ids = {t.id for t in _catalogue()}
    for spec in (common, explicit):
        assert set(spec) <= catalogue_ids
        for dist in spec.values():
            assert abs(sum(dist.values()) - 1.0) <= 1e-9
            assert all(isinstance(delta, int) for delta in dist)


def test_parse_hp_unstable_rejects_unknown_sum_and_empty():
    """Plan 2026-09-27 T4. Unknown ids, sums off by >1e-9, and '' are errors."""
    for spec in ("", "NOT_A_TRACT:0=1", "D310:0=0.5,1=0.4", "D310:0=1;D310:0=1"):
        try:
            st.parse_hp_unstable(spec)
        except ValueError:
            continue
        raise AssertionError(f"expected ValueError for {spec!r}")


def test_expand_unstable_splits_clones_by_tract_state():
    """Plan 2026-09-27 T4. One tract, two clones, no drops. Product index in the id."""
    spec = {"D310": {0: 0.5, 1: 0.5}}
    out = st.expand_unstable(
        [_FreqClone("1", 0.6), _FreqClone("2", 0.4)], spec, 0.002,
    )
    assert [s.id for s in out] == ["1.h0", "1.h1", "2.h0", "2.h1"]
    assert [s.hp_state for s in out] == [
        {"D310": 0}, {"D310": 1}, {"D310": 0}, {"D310": 1},
    ]
    got = [s.frequency for s in out]
    for freq, want in zip(got, (0.3, 0.3, 0.2, 0.2)):
        assert abs(freq - want) <= 1e-9
    assert abs(sum(got) - 1.0) <= 1e-9


def test_expand_unstable_product_order_and_min_freq_renormalise():
    """Plan 2026-09-27 T4. Spec order, deltas ascending; drop then renormalise."""
    spec = {"C573": {-1: 0.2, 0: 0.8}, "D310": {1: 0.5, 0: 0.5}}
    out = st.expand_unstable([{"id": "1", "frequency": 1.0}], spec, 0)
    assert [s.id for s in out] == ["1.h0", "1.h1", "1.h2", "1.h3"]
    assert [s.hp_state for s in out] == [
        {"C573": -1, "D310": 0},
        {"C573": -1, "D310": 1},
        {"C573": 0, "D310": 0},
        {"C573": 0, "D310": 1},
    ]
    # 0.6 * 0.01 = 0.006 is below 0.01; index 0 is dropped, index 1 kept.
    dropped = st.expand_unstable(
        [{"id": "1", "frequency": 0.6}], {"D310": {0: 0.01, 1: 0.99}}, 0.01,
    )
    assert len(dropped) == 1
    assert dropped[0].id == "1.h1"
    assert dropped[0].hp_state == {"D310": 1}
    assert abs(dropped[0].frequency - 0.6) <= 1e-9


def test_expand_unstable_keeps_modal_state_when_filter_empties_clone():
    """Plan 2026-09-27 T4. Mass stays with the clone when every combo is below the floor."""
    out = st.expand_unstable(
        [{"id": "1", "frequency": 0.005}], {"D310": {0: 0.2, 1: 0.8}}, 0.01,
    )
    assert len(out) == 1
    assert out[0].id == "1.h1"
    assert out[0].hp_state == {"D310": 1}
    assert abs(out[0].frequency - 0.005) <= 1e-12
    tied = st.expand_unstable(
        [{"id": "1", "frequency": 0.005}], {"D310": {0: 0.5, 1: 0.5}}, 0.003,
    )
    assert len(tied) == 1
    assert tied[0].id == "1.h0"
    assert tied[0].hp_state == {"D310": 0}
    assert abs(tied[0].frequency - 0.005) <= 1e-12


def test_expand_unstable_is_depth_free():
    """Plan 2026-09-27 T4. round-to-zero combinations stay; depth filtering is not here."""
    out = st.expand_unstable(
        [{"id": "1", "frequency": 1.0}], {"D310": {0: 0.999, 1: 0.001}}, 0,
    )
    assert len(out) == 2
    assert abs(out[1].frequency - 0.001) <= 1e-12
    assert round(1000 * out[1].frequency / 2) == 0


def test_subclone_h1_genome_is_parent_plus_one_c(tmp_path):
    """Plan 2026-09-27 T4. Unstable +1 is one extra C at the D310 anchor, after path edits."""
    _, seq = st.load_reference(REF)
    alt = "A" if seq[999] != "A" else "C"
    root = st.Node(id=0, variant=None, pos=None, ref=None, alt=None, parent=0)
    node = st.Node(
        id=1, variant=st.variant_id(1000, seq[999], alt),
        pos=1000, ref=seq[999], alt=alt, parent=0,
    )
    root.children.append(1)
    root.cum_freq = 1.0
    node.cum_freq = 0.5
    st.write_truth(
        st.Tree(nodes=[root, node]), seq, str(tmp_path),
        hp_spec={"D310": {0: 0.5, 1: 0.5}},
    )
    clones = (tmp_path / "truth" / "clones.tsv").read_text().splitlines()
    assert clones[0] == "clone_id\tnode_variant\tvariant_path\tfrequency\thp_state"
    by_id = {row.split("\t")[0]: row.split("\t") for row in clones[1:]}
    assert by_id["1.h0"][4] == "D310:0"
    assert by_id["1.h1"][4] == "D310:+1"
    assert by_id["1.h1"][2] == st.variant_id(1000, seq[999], alt)
    records = _fasta_records(tmp_path / "truth" / "clone_genomes.fa")
    parent = st.clone_genome(seq, st.path_to_root(st.Tree(nodes=[root, node]), 1))
    assert records["clone_1.h0"] == parent
    extra = st.left_align(st.Edit(302, "A", "AC"), parent)
    assert records["clone_1.h1"] == st.apply_edits(parent, [extra])
    assert len(records["clone_1.h1"]) == len(parent) + 1


def test_mixed_origin_when_lineage_indel_in_same_tract(tmp_path):
    """Plan 2026-09-27 T4. Lineage indel stays; unstable δ stacks; origin is mixed."""
    _, seq = st.load_reference(REF)
    d310 = next(t for t in _catalogue() if t.id == "D310")
    edit = _expected_length_edit(seq, d310, 1)
    root = st.Node(id=0, variant=None, pos=None, ref=None, alt=None, parent=0)
    node = st.Node(
        id=1, variant=st.himito_name(edit), parent=0, edit=edit,
        tract_id="D310", hp_delta=1,
    )
    root.children.append(1)
    root.cum_freq = 1.0
    node.cum_freq = 0.5
    tree = st.Tree(nodes=[root, node])
    st.write_truth(tree, seq, str(tmp_path), hp_spec={"D310": {0: 0.5, 1: 0.5}})
    lineage = st.clone_genome(seq, st.path_to_root(tree, 1))
    records = _fasta_records(tmp_path / "truth" / "clone_genomes.fa")
    assert records["clone_1.h0"] == lineage
    stacked = st.apply_edits(lineage, [st.left_align(st.Edit(302, "A", "AC"), lineage)])
    assert records["clone_1.h1"] == stacked
    assert len(stacked) == len(seq) + 2

    hp_rows = (tmp_path / "truth" / "hp_tracts.tsv").read_text().splitlines()
    assert hp_rows[0] == (
        "tract\tanchor\tref_len\tbackground_len\tdelta\ttruth_freq\tcanonical_name\torigin"
    )
    parsed = [line.split("\t") for line in hp_rows[1:]]
    assert parsed
    assert all(row[4] != "0" and row[4] != "+0" for row in parsed)
    by_origin = {row[7]: row for row in parsed}
    assert set(by_origin) == {"unstable", "mixed"}
    for origin, freq in (("unstable", 0.25), ("mixed", 0.25)):
        row = by_origin[origin]
        assert row[0] == "D310"
        assert row[1] == "302"
        assert row[2] == str(d310.copies * len(d310.unit))
        assert row[3] == row[2]
        assert row[4] == "+1"
        assert abs(float(row[5]) - freq) <= 1e-9
        assert row[6] == "m.302A>AC"
    variants = (tmp_path / "truth" / "truth_variants.txt").read_text().splitlines()
    assert variants == [st.himito_name(edit)]
    tree_rows = (tmp_path / "truth" / "truth_mutation_tree.tsv").read_text().splitlines()
    assert len(tree_rows) == 3
    assert "h0" not in (tmp_path / "truth" / "truth_mutation_tree.tsv").read_text()


def test_no_hp_unstable_clones_match_head_bytes(tmp_path):
    """Plan 2026-09-27 T4. No flag: clones.tsv and clone_genomes.fa match HEAD 8f649cd."""
    import hashlib
    import subprocess
    out = tmp_path / "plain"
    r = subprocess.run(
        [sys.executable, str(_SIM_DIR / "simulate_tree.py"),
         "--reference", REF, "--n-mutations", "10", "--seed", "1",
         "--outdir", str(out)],
        capture_output=True, text=True,
    )
    assert r.returncode == 0, r.stderr
    clones = (out / "truth" / "clones.tsv").read_bytes()
    genomes = (out / "truth" / "clone_genomes.fa").read_bytes()
    assert hashlib.sha256(clones).hexdigest() == _NO_HP_CLONES_SHA256
    assert hashlib.sha256(genomes).hexdigest() == _NO_HP_GENOMES_SHA256
    assert clones.splitlines()[0] == b"clone_id\tnode_variant\tvariant_path\tfrequency"
    assert not (out / "truth" / "hp_tracts.tsv").exists()


def test_awk_split_keeps_dotted_subclone_id(tmp_path):
    """Plan 2026-09-27 T4. simulate_reads.sh splits >clone_3.h2 into clone_3.h2.fa."""
    import subprocess
    script = (_SIM_DIR / "simulate_reads.sh").read_text()
    awk_prog = '/^>/{id=substr($0,8); f=dir"/clone_"id".fa"; print > f; next}{print >> f}'
    assert awk_prog in script
    genomes = tmp_path / "clone_genomes.fa"
    genomes.write_text(">clone_3.h2\nACGTACGT\nTTTTAAAA\n")
    work = tmp_path / "work"
    work.mkdir()
    r = subprocess.run(
        ["awk", "-v", f"dir={work}", awk_prog, str(genomes)],
        capture_output=True, text=True,
    )
    assert r.returncode == 0, r.stderr
    out = work / "clone_3.h2.fa"
    assert out.is_file()
    assert out.read_text().splitlines()[0] == ">clone_3.h2"


def test_simulate_reads_frequency_field_stays_column_four():
    """Plan 2026-09-27 T4. $1 is the clone id and $4 is the frequency."""
    script = (_SIM_DIR / "simulate_reads.sh").read_text()
    assert """awk -F'\\t' '{print $1 \"\\t\" $4}'""" in script


def test_hp_unstable_common_scheduled_depth_within_5_percent(tmp_path):
    """Plan 2026-09-27 T4. common at depth 1000 stays within 5% of the no-HP schedule."""
    import subprocess
    script = str(_SIM_DIR / "simulate_tree.py")
    py = sys.executable
    plain = tmp_path / "plain"
    hp = tmp_path / "hp"
    base = [py, script, "--reference", REF, "--n-mutations", "10", "--seed", "1"]
    r0 = subprocess.run(base + ["--outdir", str(plain)], capture_output=True, text=True)
    assert r0.returncode == 0, r0.stderr
    r1 = subprocess.run(
        base + ["--total-depth", "1000", "--hp-unstable", "common", "--outdir", str(hp)],
        capture_output=True, text=True,
    )
    assert r1.returncode == 0, r1.stderr
    no_hp = _scheduled_reads(_clone_freq_fields(plain / "truth" / "clones.tsv"), 1000)
    with_hp = _scheduled_reads(_clone_freq_fields(hp / "truth" / "clones.tsv"), 1000)
    assert with_hp <= no_hp * 1.05, f"scheduled {with_hp} vs no-HP {no_hp}"
    assert "hp-unstable: dropped" not in r1.stderr
    header = (hp / "truth" / "clones.tsv").read_text().splitlines()[0]
    assert header == "clone_id\tnode_variant\tvariant_path\tfrequency\thp_state"
    assert (hp / "truth" / "hp_tracts.tsv").is_file()


def test_hp_unstable_depth_floor_when_min_freq_would_inflate(tmp_path):
    """Plan 2026-09-27 T4. Drop round-to-zero sub-clones before renormalising."""
    import subprocess
    hp = tmp_path / "hp"
    r = subprocess.run(
        [sys.executable, str(_SIM_DIR / "simulate_tree.py"),
         "--reference", REF, "--n-mutations", "10", "--seed", "1",
         "--total-depth", "1000", "--hp-unstable", "common",
         "--hp-min-subclone-freq", "0", "--outdir", str(hp)],
        capture_output=True, text=True,
    )
    assert r.returncode == 0, r.stderr
    assert "hp-unstable: dropped" in r.stderr
    plain = tmp_path / "plain"
    r0 = subprocess.run(
        [sys.executable, str(_SIM_DIR / "simulate_tree.py"),
         "--reference", REF, "--n-mutations", "10", "--seed", "1",
         "--outdir", str(plain)],
        capture_output=True, text=True,
    )
    assert r0.returncode == 0, r0.stderr
    no_hp = _scheduled_reads(_clone_freq_fields(plain / "truth" / "clones.tsv"), 1000)
    with_hp = _scheduled_reads(_clone_freq_fields(hp / "truth" / "clones.tsv"), 1000)
    assert with_hp <= no_hp * 1.05, f"scheduled {with_hp} vs no-HP {no_hp}"
    # Parent mass is unchanged: written frequencies still partition 1.
    freqs = [float(f) for f in _clone_freq_fields(hp / "truth" / "clones.tsv")]
    assert abs(sum(freqs) - 1.0) < 1e-4


def test_unstable_alleles_omitted_from_mutation_truth(tmp_path):
    """Plan 2026-09-27 T4. hp_tracts names are not tree variants."""
    import subprocess
    out = tmp_path / "hp"
    r = subprocess.run(
        [sys.executable, str(_SIM_DIR / "simulate_tree.py"),
         "--reference", REF, "--n-mutations", "10", "--seed", "1",
         "--hp-unstable", "common", "--outdir", str(out)],
        capture_output=True, text=True,
    )
    assert r.returncode == 0, r.stderr
    variants = (out / "truth" / "truth_variants.txt").read_text()
    tree = (out / "truth" / "truth_mutation_tree.tsv").read_text()
    names = []
    for line in (out / "truth" / "hp_tracts.tsv").read_text().splitlines()[1:]:
        names.append(line.split("\t")[6])
    assert names
    for name in names:
        assert name not in variants
        assert name not in tree


def test_run_eval_forwards_hp_unstable(tmp_path):
    """Plan 2026-09-27 T4. SPEC and min freq reach simulate_tree; empty SPEC does not."""
    import subprocess
    bindir = tmp_path / "bin"
    bindir.mkdir()
    log = tmp_path / "argv.txt"
    py = bindir / "python"
    py.write_text(
        "#!/bin/bash\n"
        f'printf "%s\\n" "$*" >> "{log}"\n'
        "exit 42\n"
    )
    py.chmod(0o755)
    env = os.environ.copy()
    env["PATH"] = str(bindir) + os.pathsep + env.get("PATH", "")
    r = subprocess.run(
        ["bash", str(_SIM_DIR / "run_eval.sh"),
         "--outdir", str(tmp_path / "out"),
         "--hp-unstable", "common",
         "--hp-min-subclone-freq", "0.01",
         "--total-depth", "1000",
         "--n-mutations", "6", "--seed", "1",
         "--ref", REF],
        capture_output=True, text=True, env=env,
    )
    assert r.returncode != 0
    argv = log.read_text()
    assert "simulate_tree.py" in argv
    assert "--hp-unstable common" in argv
    assert "--hp-min-subclone-freq 0.01" in argv
    assert "--total-depth 1000" in argv

    log.unlink()
    r2 = subprocess.run(
        ["bash", str(_SIM_DIR / "run_eval.sh"),
         "--outdir", str(tmp_path / "out2"),
         "--n-mutations", "6", "--seed", "1",
         "--ref", REF],
        capture_output=True, text=True, env=env,
    )
    assert r2.returncode != 0
    argv2 = log.read_text()
    assert "--hp-unstable" not in argv2
    assert "--hp-min-subclone-freq" not in argv2
    assert "--total-depth" not in argv2


def test_sweep_seeds_forwards_hp_unstable():
    """Plan 2026-09-27 T4. sweep_seeds.sh accepts both flags and forwards them."""
    import re
    import subprocess
    script = (_SIM_DIR / "sweep_seeds.sh").read_text()
    assert re.search(r"--hp-unstable\)", script)
    assert re.search(r"--hp-min-subclone-freq\)", script)
    assert re.search(
        r'EVAL_ARGS\+=\(--hp-unstable "\$HP_UNSTABLE"\)', script,
    )
    assert re.search(
        r'EVAL_ARGS\+=\(--hp-min-subclone-freq "\$HP_MIN_SUBCLONE_FREQ"\)', script,
    )
    r = subprocess.run(
        ["bash", str(_SIM_DIR / "sweep_seeds.sh"),
         "--hp-unstable", "D310:0=0.5,1=0.5",
         "--hp-min-subclone-freq", "0.01"],
        capture_output=True, text=True,
    )
    assert "unknown arg" not in r.stderr
    assert "usage:" in r.stderr
    assert r.returncode != 0


def _run_simulate(out, *extra):
    import subprocess
    r = subprocess.run(
        [sys.executable, str(_SIM_DIR / "simulate_tree.py"),
         "--reference", REF, "--outdir", str(out), *extra],
        capture_output=True, text=True,
    )
    return r


def _variant_pos(name: str) -> int:
    import re
    match = re.match(r"m\.(\d+)", name)
    assert match, name
    return int(match.group(1))


def test_parse_background_none_and_hvs1_preset():
    """Plan 2026-09-27 T5. none is empty; hvs1-polyC expands in that order."""
    assert st.parse_background("none") == []
    assert st.parse_background("hvs1-polyC") == ["T16189C", "A16183C"]
    assert st.parse_background("T16189C,A16183C") == ["T16189C", "A16183C"]


def test_background_none_truth_matches_omitted_flag(tmp_path):
    """Plan 2026-09-27 T5. none and an omitted flag leave seed-1 truth untouched."""
    import hashlib
    plain = tmp_path / "plain"
    none = tmp_path / "none"
    r0 = _run_simulate(plain, "--n-mutations", "10", "--seed", "1")
    assert r0.returncode == 0, r0.stderr
    r1 = _run_simulate(none, "--n-mutations", "10", "--seed", "1", "--background", "none")
    assert r1.returncode == 0, r1.stderr
    plain_truth = plain / "truth"
    none_truth = none / "truth"
    assert not (plain_truth / "background.txt").exists()
    assert not (none_truth / "background.txt").exists()
    assert sorted(p.name for p in plain_truth.iterdir()) == sorted(
        p.name for p in none_truth.iterdir()
    )
    for name in sorted(p.name for p in plain_truth.iterdir()):
        assert (plain_truth / name).read_bytes() == (none_truth / name).read_bytes()
    clones = (plain_truth / "clones.tsv").read_bytes()
    genomes = (plain_truth / "clone_genomes.fa").read_bytes()
    assert hashlib.sha256(clones).hexdigest() == _NO_HP_CLONES_SHA256
    assert hashlib.sha256(genomes).hexdigest() == _NO_HP_GENOMES_SHA256


def test_background_snvs_edit_root_and_write_himito_names(tmp_path):
    """Plan 2026-09-27 T5. Root differs from rCRS only at the background sites.

    Names in truth/background.txt are Himito SNVs on rCRS, in CLI order.
    Those sites are not tree variants. AVOID_RANGES still excludes SNVs.
    """
    from hp_tracts import apply_background

    out = tmp_path / "bg"
    r = _run_simulate(
        out, "--n-mutations", "10", "--seed", "1",
        "--background", "T16189C,A16183C",
    )
    assert r.returncode == 0, r.stderr
    _, rcrs = st.load_reference(REF)
    records = _fasta_records(out / "truth" / "clone_genomes.fa")
    root = records["clone_ref"]
    assert root == apply_background(rcrs, ["T16189C", "A16183C"])
    diffs = [i + 1 for i, (a, b) in enumerate(zip(rcrs, root)) if a != b]
    assert diffs == [16183, 16189]
    others = [genome for name, genome in records.items() if name != "clone_ref"]
    assert others
    assert any(genome != root for genome in others)
    for genome in records.values():
        assert len(genome) == len(rcrs)
        assert genome[16182] == "C" and genome[16188] == "C"
    assert (out / "truth" / "background.txt").read_text() == (
        "m.16189T>C\nm.16183A>C\n"
    )
    variants = (out / "truth" / "truth_variants.txt").read_text().splitlines()
    assert "m.16189T>C" not in variants
    assert "m.16183A>C" not in variants
    positions = [_variant_pos(name) for name in variants]
    assert 16183 not in positions and 16189 not in positions
    assert all(not st._is_avoided(pos) for pos in positions)


def test_hvs1_polyc_preset_matches_explicit_tokens(tmp_path):
    """Plan 2026-09-27 T5. hvs1-polyC is T16189C,A16183C in that order."""
    preset = tmp_path / "preset"
    explicit = tmp_path / "explicit"
    r0 = _run_simulate(
        preset, "--n-mutations", "10", "--seed", "1",
        "--background", "hvs1-polyC",
    )
    assert r0.returncode == 0, r0.stderr
    r1 = _run_simulate(
        explicit, "--n-mutations", "10", "--seed", "1",
        "--background", "T16189C,A16183C",
    )
    assert r1.returncode == 0, r1.stderr
    assert (preset / "truth" / "background.txt").read_text() == (
        explicit / "truth" / "background.txt"
    ).read_text() == "m.16189T>C\nm.16183A>C\n"
    assert (preset / "truth" / "clone_genomes.fa").read_bytes() == (
        explicit / "truth" / "clone_genomes.fa"
    ).read_bytes()


def test_background_position_outside_avoid_ranges_is_not_placed():
    """Plan 2026-09-27 T5. Occupied sites are resampled on the same RNG.

    An empty or omitted set does not change the no-background draw.
    AVOID_RANGES still apply to the resampled positions.
    """
    _, seq = st.load_reference(REF)
    pos = 1000
    assert not st._is_avoided(pos) and seq[pos - 1] != "N"
    hit = None
    for seed in range(1, 20001):
        tree = st.build_tree(seq, n_mutations=40, rng=random.Random(seed))
        if any(n.pos == pos for n in tree.nodes[1:]):
            hit = seed
            break
    assert hit is not None
    plain = st.build_tree(seq, n_mutations=40, rng=random.Random(hit))
    assert any(n.pos == pos for n in plain.nodes[1:])
    blocked = st.build_tree(
        seq, n_mutations=40, rng=random.Random(hit),
        background_positions={pos},
    )
    blocked_pos = [n.pos for n in blocked.nodes[1:]]
    assert pos not in blocked_pos
    assert len(blocked_pos) == len(set(blocked_pos)) == 40
    assert all(not st._is_avoided(p) for p in blocked_pos)
    assert all(seq[p - 1] != "N" for p in blocked_pos)
    again = st.build_tree(
        seq, n_mutations=40, rng=random.Random(hit),
        background_positions=None,
    )
    empty = st.build_tree(
        seq, n_mutations=40, rng=random.Random(hit),
        background_positions=set(),
    )
    plain_pos = [n.pos for n in plain.nodes[1:]]
    assert [n.pos for n in again.nodes[1:]] == plain_pos
    assert [n.pos for n in empty.nodes[1:]] == plain_pos


def _root_only_tree():
    root = st.Node(id=0, variant=None, pos=None, ref=None, alt=None, parent=0)
    root.cum_freq = 1.0
    return st.Tree(nodes=[root])


def test_hp_tracts_use_effective_anchor_on_background(tmp_path):
    """Plan 2026-09-27 T5. HVS1 +1 on the polyC background is m.16182A>AC.

    background_len is the extended run (11 with both SNVs, 10 with T16189C).
    ref_len stays the catalogue length. The anchor is the effective anchor.
    """
    from hp_tracts import apply_background

    _, seq = st.load_reference(REF)
    both = tmp_path / "both"
    st.write_truth(
        _root_only_tree(),
        apply_background(seq, ["T16189C", "A16183C"]),
        str(both),
        hp_spec={"HVS1": {0: 0.5, 1: 0.5}},
    )
    rows = (both / "truth" / "hp_tracts.tsv").read_text().splitlines()
    assert rows[1].split("\t")[:7] == [
        "HVS1", "16182", "5", "11", "+1", "0.500000", "m.16182A>AC",
    ]

    only = tmp_path / "only"
    st.write_truth(
        _root_only_tree(),
        apply_background(seq, ["T16189C"]),
        str(only),
        hp_spec={"HVS1": {0: 0.5, 1: 0.5}},
    )
    only_row = (only / "truth" / "hp_tracts.tsv").read_text().splitlines()[1].split("\t")
    assert only_row[:7] == [
        "HVS1", "16183", "5", "10", "+1", "0.500000", "m.16183A>AC",
    ]


def test_cli_background_hp_tracts_follow_effective_anchor(tmp_path):
    """Plan 2026-09-27 T5. The CLI applies the preset before unstable deltas."""
    out = tmp_path / "hp"
    r = _run_simulate(
        out, "--n-mutations", "4", "--seed", "1",
        "--background", "hvs1-polyC",
        "--hp-unstable", "HVS1:0=0.5,1=0.5",
    )
    assert r.returncode == 0, r.stderr
    hvs = [
        line.split("\t")
        for line in (out / "truth" / "hp_tracts.tsv").read_text().splitlines()[1:]
        if line.startswith("HVS1\t")
    ]
    assert hvs
    assert all(row[1] == "16182" and row[2] == "5" and row[3] == "11" for row in hvs)
    plus = [row for row in hvs if row[4] == "+1"]
    assert plus and all(row[6] == "m.16182A>AC" for row in plus)
    assert (out / "truth" / "background.txt").read_text() == "m.16189T>C\nm.16183A>C\n"


def test_run_eval_forwards_background_not_truth_background(tmp_path):
    """Plan 2026-09-27 T6. simulate_tree gets --background; the scorer gets truth files."""
    import subprocess
    script = (_SIM_DIR / "run_eval.sh").read_text()
    assert '--truth-background "$OUTDIR/truth/background.txt"' in script
    assert '--truth-hp-tracts "$OUTDIR/truth/hp_tracts.tsv"' in script
    assert '--truth-hp-lineage "$OUTDIR/truth/hp_lineage.tsv"' in script
    assert 'SCORE_ARGS+=(--background "$BACKGROUND")' in script
    bindir = tmp_path / "bin"
    bindir.mkdir()
    log = tmp_path / "argv.txt"
    py = bindir / "python"
    py.write_text(
        "#!/bin/bash\n"
        f'printf "%s\\n" "$*" >> "{log}"\n'
        "exit 42\n"
    )
    py.chmod(0o755)
    env = os.environ.copy()
    env["PATH"] = str(bindir) + os.pathsep + env.get("PATH", "")
    r = subprocess.run(
        ["bash", str(_SIM_DIR / "run_eval.sh"),
         "--outdir", str(tmp_path / "out"),
         "--background", "hvs1-polyC",
         "--n-mutations", "6", "--seed", "1",
         "--ref", REF],
        capture_output=True, text=True, env=env,
    )
    assert r.returncode != 0
    argv = log.read_text()
    assert "simulate_tree.py" in argv
    assert "--background hvs1-polyC" in argv
    assert "--truth-background" not in argv

    log.unlink()
    r2 = subprocess.run(
        ["bash", str(_SIM_DIR / "run_eval.sh"),
         "--outdir", str(tmp_path / "out2"),
         "--n-mutations", "6", "--seed", "1",
         "--ref", REF],
        capture_output=True, text=True, env=env,
    )
    assert r2.returncode != 0
    assert "--background" not in log.read_text()


def test_sweep_seeds_forwards_background():
    """Plan 2026-09-27 T5. sweep_seeds.sh accepts --background and forwards it."""
    import re
    import subprocess
    script = (_SIM_DIR / "sweep_seeds.sh").read_text()
    assert "--truth-background" not in script
    assert re.search(r"--background\)", script)
    assert re.search(r'EVAL_ARGS\+=\(--background "\$BACKGROUND"\)', script)
    r = subprocess.run(
        ["bash", str(_SIM_DIR / "sweep_seeds.sh"), "--background", "hvs1-polyC"],
        capture_output=True, text=True,
    )
    assert "unknown arg" not in r.stderr
    assert "usage:" in r.stderr
    assert r.returncode != 0
