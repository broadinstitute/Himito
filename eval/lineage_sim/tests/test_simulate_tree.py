# eval/lineage_sim/tests/test_simulate_tree.py
import random
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import simulate_tree as st

REF = "/Users/suhang/Analysis/Himito/rCRS.fasta"


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
    st.assign_frequencies(tree, ref_fraction=0.15, rng=rng)
    for n in tree.nodes[1:]:  # every mutation
        assert 0.01 < n.cum_freq < 0.95, (n.id, n.cum_freq)


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
