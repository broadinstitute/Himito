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


def test_frequencies_decrease_from_parent_to_child():
    """cum_freq must strictly decrease along every edge (root excluded)."""
    _, seq = st.load_reference(REF)
    rng = random.Random(13)
    tree = st.build_tree(seq, n_mutations=10, rng=rng)
    st.assign_frequencies(tree, ref_fraction=0.15, rng=rng)
    for node in tree.nodes[1:]:  # skip root
        parent = tree.nodes[node.parent]
        assert node.cum_freq < parent.cum_freq, (
            f"node {node.id} cum_freq={node.cum_freq:.6f} >= "
            f"parent {parent.id} cum_freq={parent.cum_freq:.6f}"
        )


def test_write_truth_emits_expected_files_and_headers(tmp_path):
    """write_truth produces all four output files with correct headers/content."""
    import csv

    _, seq = st.load_reference(REF)
    n_mutations = 5
    rng = random.Random(99)
    tree = st.build_tree(seq, n_mutations=n_mutations, rng=rng)
    st.assign_frequencies(tree, ref_fraction=0.15, rng=rng)
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
    assert clone_lines[0] == "clone_id\tnode_variant\tvariant_path\tfrequency\tn_reads"
    clone_ids = [line.split("\t")[0] for line in clone_lines[1:]]
    assert "ref" in clone_ids, "no 'ref' clone row in clones.tsv"

    # --- clone_genomes.fa ---
    fa_file = truth_dir / "clone_genomes.fa"
    assert fa_file.exists()
    headers = [l for l in fa_file.read_text().splitlines() if l.startswith(">")]
    assert ">clone_ref" in headers
    assert len(headers) == n_mutations + 1, (
        f"expected {n_mutations + 1} FASTA records, got {len(headers)}"
    )
