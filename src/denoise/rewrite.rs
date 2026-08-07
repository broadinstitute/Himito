//! Length-changing read rewrite for the ONT denoiser.
//!
//! Every transition here preserves the read's total reference consumption:
//! `Ref->Ins` adds query only; `Ref->Del` converts M(n) to D(n) (both consume n
//! reference bases); reverting `Ins` removes query only; reverting `Del` converts
//! D(n) back to M(n). POS, reference end, and coordinate-sort order are therefore
//! invariant, so the output never needs re-sorting or re-indexing.

use std::collections::{HashMap, HashSet};

use rust_htslib::bam::record::{Cigar, CigarString};

use super::indel::Allele;

/// One read's reassignment at one normalized indel site.
#[derive(Debug, Clone)]
pub struct IndelEdit {
    /// The coordinate `rewrite_read` looks this edit up by. Direction-dependent,
    /// NOT a single anchor convention: the left-normalized SITE position when
    /// `from` is `Ref` (a gained indel has no CIGAR anchor of its own to key on
    /// -- it is keyed to wherever the walk should splice it into plain match
    /// context); the read's OWN CIGAR anchor (the last aligned base before its
    /// own indel, exactly as htslib's pileup reports it) when `from` is `Ins` or
    /// `Del` (reverting or replacing an indel the read already carries).
    ///
    /// Two edits for the SAME read can legitimately share a `ref_pos` while
    /// differing in `from`'s kind -- e.g. a read's own insertion anchored at
    /// position P that gets reverted, plus an unrelated gained deletion whose
    /// normalized site also happens to be P. `rewrite_read` disambiguates by
    /// keying its internal lookup on `(ref_pos, kind_of(from))`, not `ref_pos`
    /// alone (see FIX 1).
    pub ref_pos: u32,
    pub from: Allele,
    pub to: Allele,
}

/// Which of the three edit varieties an `IndelEdit` is, discriminated by `from`.
/// Used only to disambiguate `IndelEdit`s that share a `ref_pos` (see FIX 1):
/// `rewrite_read`'s three lookup sites each want their OWN kind independently,
/// so keying by `(ref_pos, EditKind)` instead of `ref_pos` alone makes those
/// lookups mutually exclusive rather than ambiguous.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub(crate) enum EditKind {
    Ref,
    Ins,
    Del,
}

pub(crate) fn edit_kind(a: &Allele) -> EditKind {
    match a {
        Allele::Ref => EditKind::Ref,
        Allele::Ins(_) => EditKind::Ins,
        Allele::Del(_) => EditKind::Del,
    }
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
    /// `true` iff at least one indel edit (gain/revert/replace `Ins`, or
    /// gain/revert/shrink `Del`) was actually applied during the walk.
    /// A substitution alone never sets this.
    ///
    /// Callers use this to decide whether to strip base-modification tags
    /// (`MM`/`ML`/`MN`) and stale alignment tags (`NM`/`MD`/`cs`): any indel
    /// edit shifts every query offset downstream of it and rewrites the
    /// CIGAR, invalidating those tags even when they key off query position
    /// rather than SEQ length.
    ///
    /// Comparing output SEQ length to input SEQ length is NOT a sound proxy:
    /// a read with an offsetting 1bp insertion at one site and a 1bp deletion
    /// at another has unchanged SEQ length but a changed CIGAR and shifted
    /// offsets, so a length-based check would be a false negative.
    pub structure_changed: bool,
    /// Keys of the edits this walk ACTUALLY applied -- gain/revert/replace
    /// `Ins`, or gain/revert/shrink `Del` -- keyed exactly like the internal
    /// lookup map (`ref_pos`, kind of `from`). A verbatim copy or an
    /// unsupported-target fallthrough never adds its key here.
    ///
    /// Callers must drive reassignment statistics from THIS set, not from
    /// the edit list they handed in: pass 2 predicts what this walk will do
    /// with each edit, and that prediction can be wrong (a normalized site
    /// vs. CIGAR-anchor mismatch, an insertion anchored inside the read's own
    /// deletion, a site swallowed by an earlier gained deletion, ...). When
    /// it is wrong this walk silently no-ops -- the reference-consumption
    /// invariant still holds -- but counting the edit as applied anyway
    /// would report a correction that never happened.
    pub applied: HashSet<(u32, EditKind)>,
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

/// Rewrite one read's SEQ/QUAL/CIGAR against `edits`, which are keyed by
/// REFERENCE coordinate.
///
/// Contract: `refseq` MUST be the full contig the read is aligned to,
/// indexed by absolute reference coordinate (not a windowed slice around the
/// read). `pos` MUST be the read's absolute POS (0-based) into that same
/// contig. Reverted deletions splice bases out of `refseq` at absolute
/// coordinates; a windowed slice would silently splice in the wrong bases
/// with no error.
pub fn rewrite_read(
    pos: i64,
    cigar: &CigarString,
    seq: &[u8],
    qual: &[u8],
    edits: &ReadEdits,
    refseq: &[u8],
) -> RewriteResult {
    if edits.is_empty() {
        // Common fast path: no edits, so nothing about the read changes.
        // Preserve the original CIGAR verbatim (including any `=`/`X` ops)
        // rather than running it through the walk, which folds `=`/`X` into
        // `M`.
        return RewriteResult {
            seq: seq.to_vec(),
            qual: qual.to_vec(),
            cigar: cigar.clone(),
            structure_changed: false,
            applied: HashSet::new(),
        };
    }

    let subs: HashMap<u32, u8> = edits.subs.iter().copied().collect();
    // Keyed by (ref_pos, kind_of(from)), not ref_pos alone (see FIX 1): two edits
    // for the same read can share a ref_pos while differing in `from`'s kind (a
    // gained indel keyed at a normalized site can coincide with this read's own
    // indel anchored at the same coordinate), and each of the three lookup sites
    // below wants its OWN kind specifically.
    let mut indels: HashMap<(u32, EditKind), &IndelEdit> = HashMap::new();
    for e in &edits.indels {
        let key = (e.ref_pos, edit_kind(&e.from));
        debug_assert!(
            !indels.contains_key(&key),
            "duplicate IndelEdit at ref_pos {}: two edits with the same `from` kind \
             ({:?}) collide -- each (ref_pos, kind) pair must be unique per read",
            e.ref_pos,
            key.1
        );
        indels.insert(key, e);
    }

    let mut out_seq: Vec<u8> = Vec::with_capacity(seq.len() + 8);
    let mut out_qual: Vec<u8> = Vec::with_capacity(qual.len() + 8);
    let mut ops: Vec<Cigar> = Vec::new();
    let mut structure_changed = false;
    // Edits genuinely applied, keyed exactly like `indels` above (see FIX 1
    // in the task-9 review round). Populated at each of the six points below
    // where an edit is actually acted on; a verbatim copy never touches this.
    let mut applied: HashSet<(u32, EditKind)> = HashSet::new();

    let mut rp: u32 = pos.max(0) as u32; // reference position of the next ref-consuming base
    let mut qp: usize = 0; // query position of the next original base
    // Ref position of the most recent *aligned* base emitted so far — the
    // pileup anchor convention. `None` until the first aligned base is
    // emitted, so a read whose first ref-relevant op is an `I` or `D` (no
    // aligned base precedes it) cannot alias an edit anchored at some other,
    // later site.
    let mut last_ref: Option<u32> = None;

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
                    last_ref = Some(cur);
                    k += 1;

                    // An indel anchored on this aligned base? Look up the `Ref`
                    // kind specifically: a gained indel is always keyed at the
                    // normalized site with `from: Ref` (see FIX 1), independent
                    // of any `Ins`/`Del`-kind edit that might share this same
                    // `ref_pos` via the read's own CIGAR anchor.
                    if let Some(e) = indels.get(&(cur, EditKind::Ref)) {
                        match (&e.from, &e.to) {
                            (Allele::Ref, Allele::Ins(s)) => {
                                let next = qual.get(qp + k as usize).copied();
                                let f = fill_qual(&out_qual, next);
                                out_seq.extend_from_slice(s);
                                out_qual.extend(std::iter::repeat(f).take(s.len()));
                                push_op(&mut ops, Cigar::Ins(s.len() as u32));
                                applied.insert((cur, EditKind::Ref));
                                structure_changed = true;
                            }
                            (Allele::Ref, Allele::Del(m)) => {
                                // The following m reference bases become a deletion.
                                // The flank guard in denoise.rs guarantees they lie
                                // inside this same match run.
                                debug_assert!(k + m <= n, "Ref->Del ran past its match op");
                                let m = (*m).min(n - k);
                                push_op(&mut ops, Cigar::Del(m));
                                k += m; // skip both the ref bases and their query partners
                                applied.insert((cur, EditKind::Ref));
                                structure_changed = true;
                            }
                            _ => {}
                        }
                    }
                }
                qp += n as usize;
                rp += n;
            }

            Cigar::Ins(n) => {
                // While `last_ref` is `None` there is no aligned base before
                // this insertion to anchor an edit to, so skip the lookup
                // entirely and copy the op verbatim (see FIX 3 above). The
                // `Ins` kind is baked into the key (see FIX 1), so this can
                // only ever find an edit whose `from` is this read's own
                // carried insertion -- no separate `filter` needed.
                let key = last_ref.map(|lr| (lr, EditKind::Ins));
                let e = key.and_then(|k| indels.get(&k));
                match e.map(|e| &e.to) {
                    // Reverted: drop the inserted bases entirely.
                    Some(Allele::Ref) => {
                        applied.insert(key.expect("key set whenever `e` matched"));
                        structure_changed = true;
                    }
                    // Replaced: emit the site's consensus inserted bases instead.
                    Some(Allele::Ins(s)) => {
                        let next = qual.get(qp + n as usize).copied();
                        let f = fill_qual(&out_qual, next);
                        out_seq.extend_from_slice(s);
                        out_qual.extend(std::iter::repeat(f).take(s.len()));
                        push_op(&mut ops, Cigar::Ins(s.len() as u32));
                        applied.insert(key.expect("key set whenever `e` matched"));
                        structure_changed = true;
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
                // See the `Ins` arm above: skip the lookup entirely while
                // there is no preceding aligned base to anchor to. The `Del`
                // kind is baked into the key (see FIX 1), so no separate
                // `filter` is needed here either.
                let key = last_ref.map(|lr| (lr, EditKind::Del));
                let e = key.and_then(|k| indels.get(&k));
                let splice = |out_seq: &mut Vec<u8>, out_qual: &mut Vec<u8>, from: u32, count: u32| {
                    let f = match (out_qual.last().copied(), qual.get(qp).copied()) {
                        (Some(l), Some(r)) => l.min(r),
                        (Some(l), None) => l,
                        (None, Some(r)) => r,
                        (None, None) => 0,
                    };
                    for j in 0..count {
                        let idx = (from + j) as usize;
                        debug_assert!(
                            idx < refseq.len(),
                            "refseq index {idx} out of range (refseq.len() = {}); \
                             refseq must be the full contig indexed by absolute \
                             reference coordinate",
                            refseq.len()
                        );
                        out_seq.push(refseq.get(idx).copied().unwrap_or(b'N'));
                        out_qual.push(f);
                    }
                };
                match e.map(|e| &e.to) {
                    // Reverted: the deleted reference bases come back as matches.
                    Some(Allele::Ref) => {
                        splice(&mut out_seq, &mut out_qual, rp, n);
                        push_op(&mut ops, Cigar::Match(n));
                        applied.insert(key.expect("key set whenever `e` matched"));
                        structure_changed = true;
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
                        applied.insert(key.expect("key set whenever `e` matched"));
                        structure_changed = true;
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

    RewriteResult {
        seq: out_seq,
        qual: out_qual,
        cigar: CigarString(ops),
        structure_changed,
        applied,
    }
}

/// Remove tags invalidated by a structural change.
///
/// `MM`/`ML`/`MN` are base-modification tags indexed by the read's own base
/// positions; any applied indel edit (see `RewriteResult::structure_changed`)
/// silently desynchronizes them. Output length is not a sound proxy for this:
/// an insertion at one site and a deletion at another can offset each other,
/// leaving SEQ length unchanged while the CIGAR and every downstream query
/// offset still shift. Stripping is safe here because methylation aggregation
/// runs on the ORIGINAL mt BAM (see main.rs QuickStart), not the denoised
/// one. `NM`/`MD`/`cs` describe the pre-edit alignment and are simply stale.
pub fn strip_stale_tags(rec: &mut rust_htslib::bam::Record) {
    for tag in [
        b"MM".as_slice(),
        b"ML".as_slice(),
        b"MN".as_slice(),
        b"NM".as_slice(),
        b"MD".as_slice(),
        b"cs".as_slice(),
    ] {
        // Absent tags are expected; removal failure just means "wasn't there".
        let _ = rec.remove_aux(tag);
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
        assert!(!r.structure_changed);
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
        assert!(!r.structure_changed);
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
        assert!(r.structure_changed);
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
        assert!(r.structure_changed);
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
        assert!(r.structure_changed);
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
        assert!(r.structure_changed);
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

    #[test]
    fn two_edits_sharing_a_ref_pos_but_differing_in_kind_are_both_applied() {
        // A single read can legitimately carry two IndelEdits with the SAME
        // ref_pos: one keyed at the read's own CIGAR anchor (from: Ins,
        // reverting an insertion the read itself carries), and one keyed at
        // the normalized site (from: Ref, a gained insertion into plain match
        // context) -- which happen to coincide at position 15 here. Before
        // keying by (ref_pos, kind), a plain `HashMap<u32, &IndelEdit>` could
        // only hold ONE of these two per read, so whichever was inserted
        // second would silently overwrite the first.
        //
        // 16M read that also carries its own 1bp insertion right at the end
        // (anchor 15, the read's last aligned base): REF[0..16] + inserted "A".
        let c = cig(vec![Cigar::Match(16), Cigar::Ins(1)]);
        let mut seq = REF.to_vec();
        seq.extend_from_slice(b"A");
        let qual = vec![30u8; 17];
        let edits = ReadEdits {
            subs: vec![],
            indels: vec![
                // Gained: at ref_pos 15 (from: Ref), splice in a 2bp consensus
                // insertion "AT" right after the last matched base.
                IndelEdit {
                    ref_pos: 15,
                    from: Allele::Ref,
                    to: Allele::Ins(b"AT".to_vec()),
                },
                // Reverted: this read's OWN carried insertion (from: Ins),
                // also anchored at ref_pos 15, is dropped entirely.
                IndelEdit {
                    ref_pos: 15,
                    from: Allele::Ins(b"A".to_vec()),
                    to: Allele::Ref,
                },
            ],
        };
        let r = rewrite_read(0, &c, &seq, &qual, &edits, REF);
        // The gained "AT" is spliced in (proving the Ref-kind lookup fired),
        // AND the read's own inserted "A" is dropped (proving the Ins-kind
        // lookup ALSO fired independently) -- not one overwriting the other.
        let mut expected_seq = REF.to_vec();
        expected_seq.extend_from_slice(b"AT");
        assert_eq!(r.seq, expected_seq);
        assert_eq!(r.cigar, cig(vec![Cigar::Match(16), Cigar::Ins(2)]));
        assert!(r.structure_changed);
        assert_invariants(&c, &r);
    }

    #[test]
    #[cfg(debug_assertions)]
    fn duplicate_ref_pos_and_kind_trips_the_debug_assert() {
        // Two edits with the IDENTICAL (ref_pos, kind) pair -- both `from: Ref`
        // at position 3 -- are a genuine construction bug (pass 2 should never
        // emit two gained-indel edits at the same site for one read); the
        // dedup guard must catch it rather than silently keeping whichever one
        // happened to be inserted last. `debug_assert!` is compiled out in
        // release, so this test is debug-only (see `#[cfg(debug_assertions)]`
        // above) -- there is nothing for it to catch in a release build.
        let c = cig(vec![Cigar::Match(8)]);
        let seq = b"ACGTAAAA".to_vec();
        let qual = vec![30u8; 8];
        let edits = ReadEdits {
            subs: vec![],
            indels: vec![
                IndelEdit { ref_pos: 3, from: Allele::Ref, to: Allele::Ins(b"A".to_vec()) },
                IndelEdit { ref_pos: 3, from: Allele::Ref, to: Allele::Del(1) },
            ],
        };
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            rewrite_read(0, &c, &seq, &qual, &edits, REF)
        }));
        assert!(result.is_err(), "a duplicate (ref_pos, kind) pair must trip the debug_assert");
    }

    #[test]
    fn offsetting_insertion_and_deletion_leave_length_unchanged_but_set_structure_changed() {
        // 12M read gains a 1bp insertion after ref index 3 and a 1bp
        // deletion after ref index 8. The two length changes cancel out, so
        // output SEQ length equals input SEQ length -- but the CIGAR changed
        // and every offset after ref index 3 shifted. `structure_changed`
        // must still be true: this is exactly the false-negative FIX 1
        // closes (output-length comparison would wrongly say "unchanged").
        let c = cig(vec![Cigar::Match(12)]);
        let seq = b"ACGTAAAAACGT".to_vec();
        let qual = vec![30u8; 12];
        let edits = ReadEdits {
            subs: vec![],
            indels: vec![
                IndelEdit { ref_pos: 3, from: Allele::Ref, to: Allele::Ins(b"A".to_vec()) },
                IndelEdit { ref_pos: 8, from: Allele::Ref, to: Allele::Del(1) },
            ],
        };
        let r = rewrite_read(0, &c, &seq, &qual, &edits, REF);
        assert_eq!(
            r.seq.len(),
            seq.len(),
            "offsetting insertion+deletion should leave SEQ length unchanged"
        );
        assert_eq!(r.seq, b"ACGTAAAAAAGT".to_vec());
        assert!(
            r.structure_changed,
            "an indel edit was applied even though output length is unchanged"
        );
        assert_invariants(&c, &r);
    }

    #[test]
    fn leading_insertion_is_unaffected_by_an_edit_anchored_before_any_aligned_base() {
        // `last_ref` starts as `None`: a read whose first ref-relevant op is
        // an insertion has no aligned base preceding it, so an edit
        // anchored at the read's POS must not be mistaken for one anchored
        // on this leading insertion (FIX 3).
        let c = cig(vec![Cigar::Ins(2), Cigar::Match(10)]);
        let seq = b"GGACGTAAAAAC".to_vec();
        let qual = vec![30u8; 12];
        let edits = ReadEdits {
            subs: vec![],
            indels: vec![IndelEdit {
                ref_pos: 0,
                from: Allele::Ins(b"GG".to_vec()),
                to: Allele::Ref,
            }],
        };
        let r = rewrite_read(0, &c, &seq, &qual, &edits, REF);
        assert_eq!(r.seq, seq, "leading insertion must be copied verbatim");
        assert_eq!(r.cigar, c);
        assert!(!r.structure_changed);
        assert_invariants(&c, &r);
    }

    #[test]
    fn no_edits_preserves_equal_and_diff_ops_verbatim() {
        // FIX 4's fast path: with zero edits, the original CIGAR (including
        // `=`/`X` ops) must come back byte-identical, not folded into `M`.
        let c = cig(vec![Cigar::Equal(4), Cigar::Diff(1), Cigar::Equal(3)]);
        let seq = b"ACGTAAAA".to_vec();
        let qual = vec![30u8; 8];
        let r = rewrite_read(0, &c, &seq, &qual, &ReadEdits::default(), REF);
        assert_eq!(r.cigar, c, "CIGAR must be byte-identical when there are no edits");
        assert_eq!(r.seq, seq);
        assert_eq!(r.qual, qual);
        assert!(!r.structure_changed);
        assert_invariants(&c, &r);
    }

    #[test]
    fn stale_tags_are_stripped_and_others_survive() {
        use rust_htslib::bam::record::{Aux, Record};

        let mut rec = Record::new();
        rec.set(
            b"r1",
            Some(&cig(vec![Cigar::Match(4)])),
            b"ACGT",
            &[30u8; 4],
        );
        rec.push_aux(b"MM", Aux::String("C+m,0;")).unwrap();
        rec.push_aux(b"ML", Aux::U8(200)).unwrap();
        rec.push_aux(b"MN", Aux::I32(4)).unwrap();
        rec.push_aux(b"NM", Aux::I32(2)).unwrap();
        rec.push_aux(b"MD", Aux::String("4")).unwrap();
        rec.push_aux(b"cs", Aux::String(":4")).unwrap();
        rec.push_aux(b"XY", Aux::I32(7)).unwrap();

        strip_stale_tags(&mut rec);

        assert!(rec.aux(b"MM").is_err(), "MM must be stripped");
        assert!(rec.aux(b"ML").is_err(), "ML must be stripped");
        assert!(rec.aux(b"MN").is_err(), "MN must be stripped");
        assert!(rec.aux(b"NM").is_err(), "NM must be stripped");
        assert!(rec.aux(b"MD").is_err(), "MD must be stripped");
        assert!(rec.aux(b"cs").is_err(), "cs must be stripped");
        assert!(matches!(rec.aux(b"XY").unwrap(), Aux::I32(7)), "unrelated tags must survive");
    }

    #[test]
    fn stripping_absent_tags_is_a_no_op() {
        use rust_htslib::bam::record::{Aux, Record};

        let mut rec = Record::new();
        rec.set(b"r1", Some(&cig(vec![Cigar::Match(4)])), b"ACGT", &[30u8; 4]);
        rec.push_aux(b"XY", Aux::I32(7)).unwrap();
        strip_stale_tags(&mut rec); // must not panic on missing tags
        assert!(matches!(rec.aux(b"XY").unwrap(), Aux::I32(7)));
    }
}
