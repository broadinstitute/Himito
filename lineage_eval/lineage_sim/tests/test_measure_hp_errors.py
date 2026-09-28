# Plan 2026-09-27 T7
"""Homopolymer length confusion measured from a BAM, mirroring hp_allele."""
import pathlib
import subprocess
import sys

import pysam

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import hp_tracts as hp
import measure_hp_errors as mhe

REPO = HERE.parents[2]
RCRS = REPO / "rCRS.fasta"
SCRIPT = HERE.parent / "measure_hp_errors.py"

# Tiny reference: C6 at 1-based [9, 14]. 0-based anchor 7, exclusive end 14.
# hp_allele requires the alignment to reach end + INDEL_FLANK (5).
REF_SEQ = "ACGTACGTCCCCCCGTACGTACGTAC"

# Rust indel_tests::REF29. A-run is 0-based 8..=15; anchor 7; exclusive end 16.
REF29 = "ACGTACGTAAAAAAAACGTACGTACGTAC"
HP_A, HP_E = 7, 16


def _write_fasta(path, seq, name="chrM"):
    path.write_text(f">{name}\n{seq}\n")


def _write_bam(path, ref_seq, reads, ref_name="chrM"):
    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "unsorted"},
        "SQ": [{"SN": ref_name, "LN": len(ref_seq)}],
    })
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for i, read in enumerate(reads):
            aln = pysam.AlignedSegment(header)
            aln.query_name = read.get("name", f"r{i}")
            aln.query_sequence = read["seq"]
            aln.query_qualities = pysam.qualitystring_to_array("I" * len(read["seq"]))
            aln.reference_id = 0
            aln.reference_start = read.get("start", 0)
            aln.cigartuples = read["cigar"]
            aln.flag = read.get("flag", 0)
            aln.mapping_quality = 60
            bam.write(aln)
    pysam.index(str(path))


def _rows(text):
    lines = [ln.split("\t") for ln in text.splitlines() if ln]
    assert lines[0] == ["base", "ref_len", "obs_len", "strand", "n"]
    return [
        (base, int(ref_len), int(obs_len), strand, int(n))
        for base, ref_len, obs_len, strand, n in lines[1:]
    ]


def test_three_read_synthetic_bam_confusion(tmp_path):
    """Plan 2026-09-27 T7. Match, reverse match, and a 1-bp deletion inside the run."""
    ref = tmp_path / "ref.fa"
    bam = tmp_path / "reads.bam"
    out = tmp_path / "hp_confusion.tsv"
    _write_fasta(ref, REF_SEQ)
    deleted = REF_SEQ[:10] + REF_SEQ[11:]
    _write_bam(bam, REF_SEQ, [
        {"name": "eq", "seq": REF_SEQ, "cigar": [(7, len(REF_SEQ))], "flag": 0},
        {"name": "rev", "seq": REF_SEQ, "cigar": [(0, len(REF_SEQ))], "flag": 16},
        {
            "name": "del",
            "seq": deleted,
            "cigar": [(0, 10), (2, 1), (0, 15)],
            "flag": 0,
        },
    ])
    subprocess.run(
        [sys.executable, str(SCRIPT), "--bam", str(bam), "--ref", str(ref),
         "--output", str(out)],
        check=True,
    )
    assert _rows(out.read_text()) == [
        ("C", 6, 5, "+", 1),
        ("C", 6, 6, "+", 1),
        ("C", 6, 6, "-", 1),
    ]


def test_hp_allele_of_a_pure_match_read_is_zero():
    assert mhe.hp_allele(0, [(0, 29)], REF29, REF29, HP_A, HP_E) == 0


def test_equal_and_diff_ops_do_not_change_length():
    assert mhe.hp_allele(0, [(7, 29)], REF29, REF29, HP_A, HP_E) == 0
    assert mhe.hp_allele(0, [(8, 29)], REF29, REF29, HP_A, HP_E) == 0


def test_soft_clip_does_not_shift_the_run():
    seq = "NNNNN" + REF29
    assert mhe.hp_allele(0, [(4, 5), (0, 29)], seq, REF29, HP_A, HP_E) == 0


def test_hp_allele_counts_an_insertion_of_the_run_base():
    seq = REF29[:16] + "A" + REF29[16:]
    cig = [(0, 16), (1, 1), (0, 13)]
    assert mhe.hp_allele(0, cig, seq, REF29, HP_A, HP_E) == 1


def test_hp_allele_counts_a_deletion_inside_the_run():
    seq = REF29[:10] + REF29[12:]
    cig = [(0, 10), (2, 2), (0, 17)]
    assert mhe.hp_allele(0, cig, seq, REF29, HP_A, HP_E) == -2


def test_hp_allele_of_an_offsetting_insertion_and_deletion_is_zero():
    seq = REF29[:10] + "A" + REF29[10:12] + REF29[13:]
    cig = [(0, 10), (1, 1), (0, 2), (2, 1), (0, 16)]
    assert mhe.hp_allele(0, cig, seq, REF29, HP_A, HP_E) == 0


def test_hp_allele_is_undefined_when_the_inserted_base_is_not_the_run_base():
    seq = REF29[:16] + "G" + REF29[16:]
    cig = [(0, 16), (1, 1), (0, 13)]
    assert mhe.hp_allele(0, cig, seq, REF29, HP_A, HP_E) is None


def test_hp_allele_is_undefined_for_a_read_ending_inside_the_run():
    assert mhe.hp_allele(0, [(0, 13)], REF29[:13], REF29, HP_A, HP_E) is None


def test_hp_allele_is_undefined_without_the_right_hand_flank():
    need = HP_E + mhe.INDEL_FLANK
    assert mhe.hp_allele(0, [(0, need - 1)], REF29[: need - 1], REF29, HP_A, HP_E) is None
    assert mhe.hp_allele(0, [(0, need)], REF29[:need], REF29, HP_A, HP_E) == 0


def test_hp_allele_is_undefined_when_the_anchor_base_is_not_aligned():
    seq = REF29[:5] + REF29[10:]
    cig = [(0, 5), (2, 5), (0, 19)]
    assert mhe.hp_allele(0, cig, seq, REF29, HP_A, HP_E) is None


def test_run_index_matches_hp_allele():
    """The one-pass walk used on real BAMs agrees with hp_allele."""
    runs = [(HP_A, HP_E, "A", HP_E - HP_A - 1, HP_A + 2)]
    index = mhe.RunIndex(runs)
    cases = [
        (0, [(0, 29)], REF29),
        (0, [(7, 29)], REF29),
        (0, [(8, 29)], REF29),
        (0, [(4, 5), (0, 29)], "NNNNN" + REF29),
        (0, [(0, 16), (1, 1), (0, 13)], REF29[:16] + "A" + REF29[16:]),
        (0, [(0, 10), (2, 2), (0, 17)], REF29[:10] + REF29[12:]),
        (0, [(0, 10), (1, 1), (0, 2), (2, 1), (0, 16)],
         REF29[:10] + "A" + REF29[10:12] + REF29[13:]),
        (0, [(0, 16), (1, 1), (0, 13)], REF29[:16] + "G" + REF29[16:]),
        (0, [(0, 13)], REF29[:13]),
        (0, [(0, HP_E + mhe.INDEL_FLANK - 1)], REF29[: HP_E + mhe.INDEL_FLANK - 1]),
        (0, [(0, HP_E + mhe.INDEL_FLANK)], REF29[: HP_E + mhe.INDEL_FLANK]),
        (0, [(0, 5), (2, 5), (0, 19)], REF29[:5] + REF29[10:]),
        (0, [(0, 10), (3, 2), (0, 17)], REF29[:10] + REF29[12:]),
    ]
    for pos, cigar, seq in cases:
        allele = mhe.hp_allele(pos, cigar, seq, REF29, HP_A, HP_E)
        got = index.alleles(pos, cigar, seq)
        if allele is None:
            assert got == {}, (cigar, got)
        else:
            assert got == {0: allele}, (cigar, got, allele)


def test_hp_allele_is_undefined_when_an_n_op_crosses_the_run():
    seq = REF29[:10] + REF29[12:]
    cig = [(0, 10), (3, 2), (0, 17)]
    assert mhe.hp_allele(0, cig, seq, REF29, HP_A, HP_E) is None


def test_foreign_insertion_is_uninformative_and_run_base_insertion_lengthens(tmp_path):
    ref = tmp_path / "ref.fa"
    bam = tmp_path / "reads.bam"
    _write_fasta(ref, REF_SEQ)
    ins_c = REF_SEQ[:14] + "C" + REF_SEQ[14:]
    ins_g = REF_SEQ[:14] + "G" + REF_SEQ[14:]
    _write_bam(bam, REF_SEQ, [
        {"name": "m1", "seq": REF_SEQ, "cigar": [(0, len(REF_SEQ))]},
        {"name": "m2", "seq": REF_SEQ, "cigar": [(0, len(REF_SEQ))]},
        {"name": "insC", "seq": ins_c, "cigar": [(0, 14), (1, 1), (0, 12)]},
        {"name": "insG", "seq": ins_g, "cigar": [(0, 14), (1, 1), (0, 12)]},
    ])
    assert mhe.measure(str(bam), str(ref)) == [
        ("C", 6, 6, "+", 2),
        ("C", 6, 7, "+", 1),
    ]


def test_run_dropped_when_modal_length_is_not_reference(tmp_path):
    ref = tmp_path / "ref.fa"
    bam = tmp_path / "reads.bam"
    _write_fasta(ref, REF_SEQ)
    deleted = REF_SEQ[:10] + REF_SEQ[11:]
    cig = [(0, 10), (2, 1), (0, 15)]
    _write_bam(bam, REF_SEQ, [
        {"name": "d1", "seq": deleted, "cigar": cig},
        {"name": "d2", "seq": deleted, "cigar": cig},
        {"name": "d3", "seq": deleted, "cigar": cig},
    ])
    assert mhe.measure(str(bam), str(ref)) == []


def test_tied_modal_length_drops_the_run(tmp_path):
    ref = tmp_path / "ref.fa"
    bam = tmp_path / "reads.bam"
    _write_fasta(ref, REF_SEQ)
    deleted = REF_SEQ[:10] + REF_SEQ[11:]
    _write_bam(bam, REF_SEQ, [
        {"name": "m", "seq": REF_SEQ, "cigar": [(0, len(REF_SEQ))]},
        {"name": "d", "seq": deleted, "cigar": [(0, 10), (2, 1), (0, 15)]},
    ])
    assert mhe.measure(str(bam), str(ref)) == []


def test_secondary_supplementary_and_unmapped_are_skipped(tmp_path):
    ref = tmp_path / "ref.fa"
    bam = tmp_path / "reads.bam"
    _write_fasta(ref, REF_SEQ)
    deleted = REF_SEQ[:10] + REF_SEQ[11:]
    _write_bam(bam, REF_SEQ, [
        {"name": "m", "seq": REF_SEQ, "cigar": [(0, len(REF_SEQ))], "flag": 0},
        {"name": "sec", "seq": deleted, "cigar": [(0, 10), (2, 1), (0, 15)], "flag": 256},
        {"name": "sup", "seq": deleted, "cigar": [(0, 10), (2, 1), (0, 15)], "flag": 2048},
        {"name": "un", "seq": deleted, "cigar": [(0, len(deleted))], "flag": 4, "start": 0},
    ])
    assert mhe.measure(str(bam), str(ref)) == [("C", 6, 6, "+", 1)]


def test_min_len_drops_shorter_runs(tmp_path):
    ref = tmp_path / "ref.fa"
    bam = tmp_path / "reads.bam"
    out = tmp_path / "hp_confusion.tsv"
    _write_fasta(ref, REF_SEQ)
    _write_bam(bam, REF_SEQ, [
        {"name": "m", "seq": REF_SEQ, "cigar": [(0, len(REF_SEQ))]},
    ])
    subprocess.run(
        [sys.executable, str(SCRIPT), "--bam", str(bam), "--ref", str(ref),
         "--min-len", "7", "--output", str(out)],
        check=True,
    )
    assert _rows(out.read_text()) == []


def test_tract_flags_skip_catalogue_runs(tmp_path):
    """`--tracts hp_tracts` skips catalogue ids; a comma list and a file do too."""
    _, seq = hp.load_reference(str(RCRS))
    bam = tmp_path / "full.bam"
    _write_bam(bam, seq, [
        {"name": "full", "seq": seq, "cigar": [(0, len(seq))]},
    ])
    kept = mhe.measure(str(bam), str(RCRS))
    skipped = mhe.measure(str(bam), str(RCRS), tracts="hp_tracts")
    homo = [t for t in hp.CATALOGUE if len(t.unit) == 1]

    def counts(rows):
        return {(b, rl, ol, st): n for b, rl, ol, st, n in rows}

    delta_n = 0
    kept_c, skip_c = counts(kept), counts(skipped)
    for key, n in kept_c.items():
        delta_n += n - skip_c.get(key, 0)
    assert delta_n == len(homo)
    for tract in homo:
        key = (tract.unit, tract.copies, tract.copies, "+")
        assert kept_c[key] == skip_c.get(key, 0) + sum(
            1 for t in homo if (t.unit, t.copies) == (tract.unit, tract.copies)
        )

    only = mhe.measure(str(bam), str(RCRS), tracts="D310,HVS1")
    only_c = counts(only)
    assert kept_c[("C", 7, 7, "+")] == only_c.get(("C", 7, 7, "+"), 0) + 1
    assert kept_c[("C", 5, 5, "+")] == only_c[("C", 5, 5, "+")] + 1

    tract_file = tmp_path / "tracts.txt"
    tract_file.write_text("tract\nD310\n")
    from_file = mhe.measure(str(bam), str(RCRS), tracts=str(tract_file))
    file_c = counts(from_file)
    assert kept_c[("C", 7, 7, "+")] == file_c.get(("C", 7, 7, "+"), 0) + 1
    assert kept_c[("C", 5, 5, "+")] == file_c[("C", 5, 5, "+")]
