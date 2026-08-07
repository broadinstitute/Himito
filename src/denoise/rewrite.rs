//! Length-changing read rewrite for the ONT denoiser.
//!
//! Every transition here preserves the read's total reference consumption:
//! `Ref->Ins` adds query only; `Ref->Del` converts M(n) to D(n) (both consume n
//! reference bases); reverting `Ins` removes query only; reverting `Del` converts
//! D(n) back to M(n). POS, reference end, and coordinate-sort order are therefore
//! invariant, so the output never needs re-sorting or re-indexing.

use std::collections::HashMap;

use rust_htslib::bam::record::{Cigar, CigarString};

use super::indel::Allele;

/// One read's reassignment at one normalized indel site.
#[derive(Debug, Clone)]
pub struct IndelEdit {
    /// Reference position of the last base BEFORE the event (the pileup anchor).
    pub ref_pos: u32,
    pub from: Allele,
    pub to: Allele,
}

/// All edits for one read, keyed by REFERENCE coordinate so substitutions and
/// length-changing indel edits can be applied in a single ordered walk.
#[derive(Debug, Default, Clone)]
pub struct ReadEdits {
    pub subs: Vec<(u32, u8)>,
    pub indels: Vec<IndelEdit>,
}

impl ReadEdits {
    pub fn is_empty(&self) -> bool {
        self.subs.is_empty() && self.indels.is_empty()
    }
}

pub struct RewriteResult {
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
    pub cigar: CigarString,
    pub len_changed: bool,
}

/// Total reference bases consumed by a CIGAR (M/=/X/D/N). The rewrite invariant.
pub fn ref_consumed(c: &CigarString) -> u32 {
    c.iter()
        .map(|op| match op {
            Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) | Cigar::Del(n) | Cigar::RefSkip(n) => *n,
            _ => 0,
        })
        .sum()
}

fn push_op(ops: &mut Vec<Cigar>, op: Cigar) {
    let len = match op {
        Cigar::Match(n) | Cigar::Ins(n) | Cigar::Del(n) | Cigar::RefSkip(n)
        | Cigar::SoftClip(n) | Cigar::HardClip(n) | Cigar::Pad(n) | Cigar::Equal(n)
        | Cigar::Diff(n) => n,
    };
    if len == 0 {
        return;
    }
    if let Some(last) = ops.last_mut() {
        match (&last, &op) {
            (Cigar::Match(a), Cigar::Match(b)) => { *last = Cigar::Match(a + b); return; }
            (Cigar::Ins(a), Cigar::Ins(b)) => { *last = Cigar::Ins(a + b); return; }
            (Cigar::Del(a), Cigar::Del(b)) => { *last = Cigar::Del(a + b); return; }
            (Cigar::SoftClip(a), Cigar::SoftClip(b)) => { *last = Cigar::SoftClip(a + b); return; }
            _ => {}
        }
    }
    ops.push(op);
}

pub fn rewrite_read(
    pos: i64,
    cigar: &CigarString,
    seq: &[u8],
    qual: &[u8],
    edits: &ReadEdits,
    refseq: &[u8],
) -> RewriteResult {
    let subs: HashMap<u32, u8> = edits.subs.iter().copied().collect();
    let indels: HashMap<u32, &IndelEdit> =
        edits.indels.iter().map(|e| (e.ref_pos, e)).collect();

    let mut out_seq: Vec<u8> = Vec::with_capacity(seq.len() + 8);
    let mut out_qual: Vec<u8> = Vec::with_capacity(qual.len() + 8);
    let mut ops: Vec<Cigar> = Vec::new();

    let mut rp: u32 = pos.max(0) as u32; // reference position of the next ref-consuming base
    let mut qp: usize = 0; // query position of the next original base
    let mut last_ref: u32 = rp; // ref position of the most recent aligned base

    // Quality assigned to synthesized bases: never more confident than the flanks.
    let fill_qual = |out_qual: &Vec<u8>, next: Option<u8>| -> u8 {
        match (out_qual.last().copied(), next) {
            (Some(l), Some(r)) => l.min(r),
            (Some(l), None) => l,
            (None, Some(r)) => r,
            (None, None) => 0,
        }
    };

    for op in cigar.iter() {
        match *op {
            Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                let mut k: u32 = 0;
                while k < n {
                    let cur = rp + k;
                    let qi = qp + k as usize;
                    out_seq.push(*subs.get(&cur).unwrap_or(&seq[qi]));
                    out_qual.push(qual[qi]);
                    push_op(&mut ops, Cigar::Match(1));
                    last_ref = cur;
                    k += 1;

                    // An indel anchored on this aligned base?
                    if let Some(e) = indels.get(&cur) {
                        match (&e.from, &e.to) {
                            (Allele::Ref, Allele::Ins(s)) => {
                                let next = qual.get(qp + k as usize).copied();
                                let f = fill_qual(&out_qual, next);
                                out_seq.extend_from_slice(s);
                                out_qual.extend(std::iter::repeat(f).take(s.len()));
                                push_op(&mut ops, Cigar::Ins(s.len() as u32));
                            }
                            (Allele::Ref, Allele::Del(m)) => {
                                // The following m reference bases become a deletion.
                                // The flank guard in denoise.rs guarantees they lie
                                // inside this same match run.
                                debug_assert!(k + m <= n, "Ref->Del ran past its match op");
                                let m = (*m).min(n - k);
                                push_op(&mut ops, Cigar::Del(m));
                                k += m; // skip both the ref bases and their query partners
                            }
                            _ => {}
                        }
                    }
                }
                qp += n as usize;
                rp += n;
            }

            Cigar::Ins(n) => {
                let e = indels
                    .get(&last_ref)
                    .filter(|e| matches!(e.from, Allele::Ins(_)));
                match e.map(|e| &e.to) {
                    // Reverted: drop the inserted bases entirely.
                    Some(Allele::Ref) => {}
                    // Replaced: emit the site's consensus inserted bases instead.
                    Some(Allele::Ins(s)) => {
                        let next = qual.get(qp + n as usize).copied();
                        let f = fill_qual(&out_qual, next);
                        out_seq.extend_from_slice(s);
                        out_qual.extend(std::iter::repeat(f).take(s.len()));
                        push_op(&mut ops, Cigar::Ins(s.len() as u32));
                    }
                    // No edit, or an unsupported cross-kind target: copy verbatim.
                    _ => {
                        out_seq.extend_from_slice(&seq[qp..qp + n as usize]);
                        out_qual.extend_from_slice(&qual[qp..qp + n as usize]);
                        push_op(&mut ops, Cigar::Ins(n));
                    }
                }
                qp += n as usize;
            }

            Cigar::Del(n) => {
                let e = indels
                    .get(&last_ref)
                    .filter(|e| matches!(e.from, Allele::Del(_)));
                let splice = |out_seq: &mut Vec<u8>, out_qual: &mut Vec<u8>, from: u32, count: u32| {
                    let f = match (out_qual.last().copied(), qual.get(qp).copied()) {
                        (Some(l), Some(r)) => l.min(r),
                        (Some(l), None) => l,
                        (None, Some(r)) => r,
                        (None, None) => 0,
                    };
                    for j in 0..count {
                        let idx = (from + j) as usize;
                        out_seq.push(refseq.get(idx).copied().unwrap_or(b'N'));
                        out_qual.push(f);
                    }
                };
                match e.map(|e| &e.to) {
                    // Reverted: the deleted reference bases come back as matches.
                    Some(Allele::Ref) => {
                        splice(&mut out_seq, &mut out_qual, rp, n);
                        push_op(&mut ops, Cigar::Match(n));
                    }
                    // Shrunk: left-aligned, so the shorter deletion comes first and
                    // the freed reference bases follow as matches.
                    Some(Allele::Del(m)) if *m <= n => {
                        push_op(&mut ops, Cigar::Del(*m));
                        let back = n - *m;
                        if back > 0 {
                            splice(&mut out_seq, &mut out_qual, rp + *m, back);
                            push_op(&mut ops, Cigar::Match(back));
                        }
                    }
                    // No edit, or an unsupported target: copy verbatim.
                    _ => push_op(&mut ops, Cigar::Del(n)),
                }
                rp += n;
            }

            Cigar::SoftClip(n) => {
                out_seq.extend_from_slice(&seq[qp..qp + n as usize]);
                out_qual.extend_from_slice(&qual[qp..qp + n as usize]);
                push_op(&mut ops, Cigar::SoftClip(n));
                qp += n as usize;
            }
            Cigar::RefSkip(n) => {
                push_op(&mut ops, Cigar::RefSkip(n));
                rp += n;
            }
            Cigar::HardClip(n) => push_op(&mut ops, Cigar::HardClip(n)),
            Cigar::Pad(n) => push_op(&mut ops, Cigar::Pad(n)),
        }
    }

    let len_changed = out_seq.len() != seq.len();
    RewriteResult {
        seq: out_seq,
        qual: out_qual,
        cigar: CigarString(ops),
        len_changed,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Reference used by every test below (index 0-based):
    //   0123456789...
    //   A C G T A A A A A C G T A C G T
    const REF: &[u8] = b"ACGTAAAAACGTACGT";

    fn cig(ops: Vec<Cigar>) -> CigarString {
        CigarString(ops)
    }

    /// Assert the design's core invariant plus SEQ/QUAL length agreement.
    fn assert_invariants(before: &CigarString, r: &RewriteResult) {
        assert_eq!(
            ref_consumed(before),
            ref_consumed(&r.cigar),
            "reference consumption changed: {before:?} -> {:?}",
            r.cigar
        );
        assert_eq!(r.seq.len(), r.qual.len(), "SEQ and QUAL lengths disagree");
        let query: u32 = r
            .cigar
            .iter()
            .map(|op| match op {
                Cigar::Match(n) | Cigar::Ins(n) | Cigar::SoftClip(n) | Cigar::Equal(n)
                | Cigar::Diff(n) => *n,
                _ => 0,
            })
            .sum();
        assert_eq!(query as usize, r.seq.len(), "CIGAR query length != SEQ length");
    }

    #[test]
    fn no_edits_is_a_faithful_copy() {
        let c = cig(vec![Cigar::Match(8)]);
        let seq = b"ACGTAAAA".to_vec();
        let qual = vec![30u8; 8];
        let r = rewrite_read(0, &c, &seq, &qual, &ReadEdits::default(), REF);
        assert_eq!(r.seq, seq);
        assert_eq!(r.qual, qual);
        assert_eq!(r.cigar, c);
        assert!(!r.len_changed);
        assert_invariants(&c, &r);
    }

    #[test]
    fn substitution_only_preserves_length_and_cigar() {
        let c = cig(vec![Cigar::Match(8)]);
        let seq = b"ACGTAAAA".to_vec();
        let qual = vec![30u8; 8];
        let edits = ReadEdits { subs: vec![(2, b'T')], indels: vec![] };
        let r = rewrite_read(0, &c, &seq, &qual, &edits, REF);
        assert_eq!(r.seq, b"ACTTAAAA".to_vec());
        assert_eq!(r.qual, qual);
        assert!(!r.len_changed);
        assert_invariants(&c, &r);
    }

    #[test]
    fn reverting_an_insertion_removes_query_bases() {
        // Read = ACGT + "GG" inserted after ref index 3 + AAAA.
        let c = cig(vec![Cigar::Match(4), Cigar::Ins(2), Cigar::Match(4)]);
        let seq = b"ACGTGGAAAA".to_vec();
        let qual = vec![30u8; 10];
        let edits = ReadEdits {
            subs: vec![],
            indels: vec![IndelEdit { ref_pos: 3, from: Allele::Ins(b"GG".to_vec()), to: Allele::Ref }],
        };
        let r = rewrite_read(0, &c, &seq, &qual, &edits, REF);
        assert_eq!(r.seq, b"ACGTAAAA".to_vec());
        assert_eq!(r.cigar, cig(vec![Cigar::Match(8)]));
        assert!(r.len_changed);
        assert_invariants(&c, &r);
    }

    #[test]
    fn reverting_a_deletion_splices_reference_bases_back() {
        // Read skips ref indices 4..6 (AA), i.e. 4M 2D 3M over ACGT|AA|AAA...
        let c = cig(vec![Cigar::Match(4), Cigar::Del(2), Cigar::Match(3)]);
        let seq = b"ACGTAAA".to_vec();
        let qual = vec![20u8, 21, 22, 23, 24, 25, 26];
        let edits = ReadEdits {
            subs: vec![],
            indels: vec![IndelEdit { ref_pos: 3, from: Allele::Del(2), to: Allele::Ref }],
        };
        let r = rewrite_read(0, &c, &seq, &qual, &edits, REF);
        assert_eq!(r.seq, b"ACGTAAAAA".to_vec());
        assert_eq!(r.cigar, cig(vec![Cigar::Match(9)]));
        // Spliced bases take min(left flank, right flank) = min(23, 24) = 23.
        assert_eq!(r.qual, vec![20, 21, 22, 23, 23, 23, 24, 25, 26]);
        assert!(r.len_changed);
        assert_invariants(&c, &r);
    }

    #[test]
    fn gaining_an_insertion_splits_a_match_run() {
        // Plain 8M read gains a consensus "A" insertion after ref index 3.
        let c = cig(vec![Cigar::Match(8)]);
        let seq = b"ACGTAAAA".to_vec();
        let qual = vec![20u8, 21, 22, 23, 24, 25, 26, 27];
        let edits = ReadEdits {
            subs: vec![],
            indels: vec![IndelEdit { ref_pos: 3, from: Allele::Ref, to: Allele::Ins(b"A".to_vec()) }],
        };
        let r = rewrite_read(0, &c, &seq, &qual, &edits, REF);
        assert_eq!(r.seq, b"ACGTAAAAA".to_vec());
        assert_eq!(r.cigar, cig(vec![Cigar::Match(4), Cigar::Ins(1), Cigar::Match(4)]));
        // Inserted base takes min(23, 24) = 23.
        assert_eq!(r.qual, vec![20, 21, 22, 23, 23, 24, 25, 26, 27]);
        assert!(r.len_changed);
        assert_invariants(&c, &r);
    }

    #[test]
    fn gaining_a_deletion_converts_match_to_del() {
        // Plain 8M read gains a 2bp deletion after ref index 3.
        let c = cig(vec![Cigar::Match(8)]);
        let seq = b"ACGTAAAA".to_vec();
        let qual = vec![30u8; 8];
        let edits = ReadEdits {
            subs: vec![],
            indels: vec![IndelEdit { ref_pos: 3, from: Allele::Ref, to: Allele::Del(2) }],
        };
        let r = rewrite_read(0, &c, &seq, &qual, &edits, REF);
        assert_eq!(r.seq, b"ACGTAA".to_vec());
        assert_eq!(r.cigar, cig(vec![Cigar::Match(4), Cigar::Del(2), Cigar::Match(2)]));
        assert!(r.len_changed);
        assert_invariants(&c, &r);
    }

    #[test]
    fn replacing_an_insertion_swaps_the_inserted_bases() {
        let c = cig(vec![Cigar::Match(4), Cigar::Ins(2), Cigar::Match(4)]);
        let seq = b"ACGTGGAAAA".to_vec();
        let qual = vec![30u8; 10];
        let edits = ReadEdits {
            subs: vec![],
            indels: vec![IndelEdit {
                ref_pos: 3,
                from: Allele::Ins(b"GG".to_vec()),
                to: Allele::Ins(b"A".to_vec()),
            }],
        };
        let r = rewrite_read(0, &c, &seq, &qual, &edits, REF);
        assert_eq!(r.seq, b"ACGTAAAAA".to_vec());
        assert_eq!(r.cigar, cig(vec![Cigar::Match(4), Cigar::Ins(1), Cigar::Match(4)]));
        assert_invariants(&c, &r);
    }

    #[test]
    fn shrinking_a_deletion_emits_del_then_match() {
        // 3bp deletion shrinks to 1bp: Del(1) first (left-aligned), then the two
        // reference bases that stop being deleted come back as matches.
        let c = cig(vec![Cigar::Match(4), Cigar::Del(3), Cigar::Match(2)]);
        let seq = b"ACGTAC".to_vec();
        let qual = vec![30u8; 6];
        let edits = ReadEdits {
            subs: vec![],
            indels: vec![IndelEdit { ref_pos: 3, from: Allele::Del(3), to: Allele::Del(1) }],
        };
        let r = rewrite_read(0, &c, &seq, &qual, &edits, REF);
        // ref[5..7] = "AA" return to the read.
        assert_eq!(r.seq, b"ACGTAAAC".to_vec());
        assert_eq!(
            r.cigar,
            cig(vec![Cigar::Match(4), Cigar::Del(1), Cigar::Match(4)])
        );
        assert_invariants(&c, &r);
    }

    #[test]
    fn one_read_carrying_sub_insertion_and_deletion_edits() {
        // 12M read: substitution at ref 1, gained insertion after ref 3, gained
        // deletion after ref 8. All applied in a single walk, left to right.
        let c = cig(vec![Cigar::Match(12)]);
        let seq = b"ACGTAAAAACGT".to_vec();
        let qual = vec![30u8; 12];
        let edits = ReadEdits {
            subs: vec![(1, b'T')],
            indels: vec![
                IndelEdit { ref_pos: 3, from: Allele::Ref, to: Allele::Ins(b"A".to_vec()) },
                IndelEdit { ref_pos: 8, from: Allele::Ref, to: Allele::Del(2) },
            ],
        };
        let r = rewrite_read(0, &c, &seq, &qual, &edits, REF);
        assert_eq!(r.seq, b"ATGTAAAAAAT".to_vec());
        assert_eq!(
            r.cigar,
            cig(vec![Cigar::Match(4), Cigar::Ins(1), Cigar::Match(5), Cigar::Del(2), Cigar::Match(1)])
        );
        assert_invariants(&c, &r);
    }

    #[test]
    fn soft_clips_are_copied_verbatim() {
        let c = cig(vec![Cigar::SoftClip(3), Cigar::Match(8), Cigar::SoftClip(2)]);
        let seq = b"TTTACGTAAAAGG".to_vec();
        let qual = vec![30u8; 13];
        let edits = ReadEdits {
            subs: vec![],
            indels: vec![IndelEdit { ref_pos: 3, from: Allele::Ref, to: Allele::Ins(b"A".to_vec()) }],
        };
        let r = rewrite_read(0, &c, &seq, &qual, &edits, REF);
        assert_eq!(&r.seq[..3], b"TTT");
        assert_eq!(&r.seq[r.seq.len() - 2..], b"GG");
        assert_eq!(
            r.cigar,
            cig(vec![Cigar::SoftClip(3), Cigar::Match(4), Cigar::Ins(1), Cigar::Match(4), Cigar::SoftClip(2)])
        );
        assert_invariants(&c, &r);
    }

    #[test]
    fn read_starting_at_nonzero_pos_uses_absolute_reference_coordinates() {
        // Read aligned at POS 4 (the AAAAA run). Edit anchored at absolute ref 5.
        let c = cig(vec![Cigar::Match(6)]);
        let seq = b"AAAAAC".to_vec();
        let qual = vec![30u8; 6];
        let edits = ReadEdits {
            subs: vec![],
            indels: vec![IndelEdit { ref_pos: 5, from: Allele::Ref, to: Allele::Del(1) }],
        };
        let r = rewrite_read(4, &c, &seq, &qual, &edits, REF);
        assert_eq!(r.seq, b"AAAAC".to_vec());
        assert_eq!(r.cigar, cig(vec![Cigar::Match(2), Cigar::Del(1), Cigar::Match(3)]));
        assert_invariants(&c, &r);
    }
}
