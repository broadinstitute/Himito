//! Read-linkage rescue of low-frequency SNVs: a sub-2% lineage marker rides on the same molecules as the
//! other markers; a sequencing error does not.

use std::path::Path;
use std::time::Instant;

use anyhow::{Context, Result};
use rayon::prelude::*;
use rust_htslib::bam::{self, ext::BamRecordExtensions, record::{Aux, Cigar}, Read as HtsRead};
use rustc_hash::FxHashMap;
use crate::denoise::base_to_idx;
use crate::scite::lchoose;

// Candidate alleles: HF in [MIN_HF, MAX_HF) and at least MIN_AC carriers among the good reads.
const MIN_HF: f64 = 0.004;
const MAX_HF: f64 = 0.20;
const MIN_AC: u32 = 10;
/// Pairs at most this many bp apart are not tested.
const MIN_DIST: u32 = 20;
/// Reads whose error rate is above this quantile are dropped.
const ERR_QUANTILE: f64 = 0.95;
// Pair rule thresholds, on both strands together and on each strand.
const MIN_SHARED: u32 = 5;
const MIN_SHARED_STRAND: u32 = 3;
const MIN_COCOVER: u32 = 50;
const MIN_COCOVER_STRAND: u32 = 30;
const MIN_ENRICH: f64 = 10.0;
const MAX_LOG10P: f64 = -10.0;
const MAX_STRAND_LOG10P: f64 = -5.0;

/// An accepted allele with its most significant partner; `pos` is 0-based like `call::Variant::pos` for SNPs.
#[derive(Debug, Clone, PartialEq)]
pub struct LinkCall {
    pub pos: usize,
    pub ref_allele: String,
    pub alt_allele: String,
    pub ac: u32,
    pub dp: u32,
    pub fwd: u32,
    pub rev: u32,
    pub partner_name: String,
    pub shared: u32,
    pub enrich: f64,
    pub log10p: f64,
}

#[derive(Debug, Default)]
pub struct LinkSummary {
    pub reads: usize,
    pub bad_reads: usize,
    pub candidates: usize,
    pub pairs_tested: usize,
    pub pairs_passing: usize,
    pub accepted: usize,
}

/// log10 P(X >= k) for X ~ Hypergeometric(population `n_total`, `k_succ` successes, `n_draw` draws).
pub fn log10_hyper_sf(k: u32, n_total: u32, k_succ: u32, n_draw: u32) -> f64 {
    if k == 0 {
        return 0.0;
    }
    let k_succ = k_succ.min(n_total);
    let n_draw = n_draw.min(n_total);
    let lo = k.max(n_draw.saturating_sub(n_total - k_succ));
    let hi = k_succ.min(n_draw);
    if lo > hi {
        return f64::NEG_INFINITY;
    }
    let ln_den = lchoose(n_total as f64, n_draw as f64);
    let terms: Vec<f64> = (lo..=hi)
        .map(|x| {
            lchoose(k_succ as f64, x as f64) + lchoose((n_total - k_succ) as f64, (n_draw - x) as f64) - ln_den
        })
        .collect();
    let m = terms.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let ln_sf = m + terms.iter().map(|t| (t - m).exp()).sum::<f64>().ln();
    (ln_sf / std::f64::consts::LN_10).min(0.0)
}

struct Fenwick {
    t: Vec<u32>,
}

impl Fenwick {
    fn new(n: usize) -> Self {
        Fenwick { t: vec![0; n + 1] }
    }

    fn add(&mut self, i: usize) {
        let mut i = i + 1;
        while i < self.t.len() {
            self.t[i] += 1;
            i += i & i.wrapping_neg();
        }
    }

    /// Sum over indices `0..=i`.
    fn prefix(&self, i: usize) -> u32 {
        let mut i = (i + 1).min(self.t.len() - 1);
        let mut s = 0;
        while i > 0 {
            s += self.t[i];
            i &= i - 1;
        }
        s
    }
}

/// For each query `(lo, hi)` with `lo <= hi`: the (forward, reverse) number of half-open spans with
/// `start <= lo` and `end > hi`.
fn cocover_counts(spans: &[(u32, u32, bool)], queries: &[(u32, u32)], len: u32) -> Vec<[u32; 2]> {
    let mut by_start: Vec<usize> = (0..spans.len()).collect();
    by_start.sort_unstable_by_key(|&i| spans[i].0);
    let mut by_lo: Vec<usize> = (0..queries.len()).collect();
    by_lo.sort_unstable_by_key(|&i| queries[i].0);

    let mut trees = [Fenwick::new(len as usize + 1), Fenwick::new(len as usize + 1)];
    let mut inserted = [0u32; 2];
    let mut out = vec![[0, 0]; queries.len()];
    let mut si = 0;
    for qi in by_lo {
        let (lo, hi) = queries[qi];
        while si < by_start.len() && spans[by_start[si]].0 <= lo {
            let (_, end, rev) = spans[by_start[si]];
            trees[rev as usize].add(end.min(len) as usize);
            inserted[rev as usize] += 1;
            si += 1;
        }
        out[qi] = [0, 1].map(|s| inserted[s] - trees[s].prefix(hi as usize));
    }
    out
}

struct ReadInfo {
    reverse: bool,
    has_primary: bool,
    err: f64,
    spans: Vec<(u32, u32)>,
}

struct Candidate {
    pos: u32,
    ref_base: u8,
    alt_base: u8,
    ac: u32,
    dp: u32,
    fwd: u32,
    rev: u32,
    reads: Vec<u32>,
}

impl Candidate {
    fn name(&self) -> String {
        format!("m.{}{}>{}", self.pos + 1, self.ref_base as char, self.alt_base as char)
    }
}

fn nm_tag(rec: &bam::Record) -> Option<f64> {
    match rec.aux(b"NM") {
        Ok(Aux::I8(v)) => Some(v as f64),
        Ok(Aux::U8(v)) => Some(v as f64),
        Ok(Aux::I16(v)) => Some(v as f64),
        Ok(Aux::U16(v)) => Some(v as f64),
        Ok(Aux::I32(v)) => Some(v as f64),
        Ok(Aux::U32(v)) => Some(v as f64),
        _ => None,
    }
}

/// Query bases in the alignment, soft clips excluded (pysam's `query_alignment_length`).
fn query_alignment_length(rec: &bam::Record) -> u32 {
    rec.cigar()
        .iter()
        .map(|c| match c {
            Cigar::Match(l) | Cigar::Ins(l) | Cigar::Equal(l) | Cigar::Diff(l) => *l,
            _ => 0,
        })
        .sum()
}

fn span_index(spans: &[(u32, u32)], p: u32) -> Option<usize> {
    spans.iter().position(|&(s, e)| s <= p && p < e)
}

/// Pair rule on both strands together or on one strand (`per_strand`): of `n_total` reads co-covering
/// alleles A and B, `k_succ` carry A, `n_draw` carry B and `k` carry both. `Some((enrichment, log10 p))` if it passes.
fn pair_rule(k: u32, n_total: u32, k_succ: u32, n_draw: u32, per_strand: bool) -> Option<(f64, f64)> {
    let (min_cover, min_shared, max_log10p) = if per_strand {
        (MIN_COCOVER_STRAND, MIN_SHARED_STRAND, MAX_STRAND_LOG10P)
    } else {
        (MIN_COCOVER, MIN_SHARED, MAX_LOG10P)
    };
    if n_total < min_cover || k < min_shared {
        return None;
    }
    let enrich = k as f64 * n_total as f64 / (k_succ as f64 * n_draw as f64);
    if enrich < MIN_ENRICH {
        return None;
    }
    let lp = log10_hyper_sf(k, n_total, k_succ, n_draw);
    (lp <= max_log10p).then_some((enrich, lp))
}

/// Read-linkage rescue on `ref_name` of a coordinate-sorted BAM (the raw, pre-denoise reads).
/// Returns every accepted allele, called already or not, sorted by position.
pub fn rescue(bam_path: &Path, refseq: &[u8], ref_name: &str) -> Result<(Vec<LinkCall>, LinkSummary)> {
    let t0 = Instant::now();
    let mut summary = LinkSummary::default();
    let open = || bam::Reader::from_path(bam_path).with_context(|| format!("cannot open linkage BAM {bam_path:?}"));

    // Pass 1, per read: strand of the primary, error rate (max NM / aligned query length over its
    // records) and merged reference spans.
    let mut reader = open()?;
    let tid = match reader.header().tid(ref_name.as_bytes()) {
        Some(t) => t as i32,
        None => {
            eprintln!("Warning: linkage rescue: contig {ref_name} not in {}; skipped", bam_path.display());
            return Ok((vec![], summary));
        }
    };
    let mut ids: FxHashMap<Vec<u8>, u32> = FxHashMap::default();
    let mut reads: Vec<ReadInfo> = Vec::new();
    let mut warned_nm = false;
    let mut rec = bam::Record::new();
    while let Some(r) = reader.read(&mut rec) {
        r.context("linkage BAM read error")?;
        if rec.is_unmapped() || rec.is_secondary() || rec.tid() != tid {
            continue;
        }
        let id = match ids.get(rec.qname()) {
            Some(&id) => id,
            None => {
                let id = reads.len() as u32;
                ids.insert(rec.qname().to_vec(), id);
                reads.push(ReadInfo { reverse: rec.is_reverse(), has_primary: false, err: 0.0, spans: vec![] });
                id
            }
        };
        let info = &mut reads[id as usize];
        if !rec.is_supplementary() && !info.has_primary {
            info.reverse = rec.is_reverse();
            info.has_primary = true;
        }
        let nm = nm_tag(&rec);
        if nm.is_none() && !warned_nm {
            eprintln!("Warning: linkage rescue: records without an NM tag count as error rate 0");
            warned_nm = true;
        }
        let qlen = query_alignment_length(&rec);
        if qlen > 0 {
            info.err = info.err.max(nm.unwrap_or(0.0) / qlen as f64);
        }
        let mut blocks = rec.aligned_blocks();
        if let Some(first) = blocks.next() {
            let last = blocks.last().unwrap_or(first);
            info.spans.push((first[0].max(0) as u32, last[1].max(0) as u32));
        }
    }
    for info in reads.iter_mut() {
        info.spans.sort_unstable();
        let mut merged: Vec<(u32, u32)> = Vec::with_capacity(info.spans.len());
        for &(s, e) in &info.spans {
            match merged.last_mut() {
                Some(m) if s <= m.1 => m.1 = m.1.max(e),
                _ => merged.push((s, e)),
            }
        }
        info.spans = merged;
    }
    summary.reads = reads.len();
    if reads.is_empty() {
        return Ok((vec![], summary));
    }
    let mut errs: Vec<f64> = reads.iter().map(|r| r.err).collect();
    errs.sort_unstable_by(|a, b| a.total_cmp(b));
    let err_cut = errs[(ERR_QUANTILE * (errs.len() - 1) as f64) as usize];
    let bad: Vec<bool> = reads.iter().map(|r| r.err > err_cut).collect();
    summary.bad_reads = bad.iter().filter(|&&b| b).count();

    // Pass 2: candidate alleles from a pileup of the good reads.
    let mut reader = open()?;
    let mut plp = reader.pileup();
    plp.set_max_depth(i32::MAX as u32);
    let mut last_col = vec![u32::MAX; reads.len()];
    let mut cands: Vec<Candidate> = Vec::new();
    for p in plp {
        let col = p.context("linkage pileup error")?;
        if col.tid() as i32 != tid {
            continue;
        }
        let pos = col.pos();
        let Some(ref_base) = refseq.get(pos as usize).map(|b| b.to_ascii_uppercase()) else { continue };
        let Some(ref_idx) = base_to_idx(ref_base) else { continue };
        let mut counts = [0u32; 4];
        let mut depth = 0u32;
        let mut alt_obs: Vec<(u32, usize)> = Vec::new();
        for a in col.alignments() {
            if a.is_del() || a.is_refskip() {
                continue;
            }
            let Some(qpos) = a.qpos() else { continue };
            let rec = a.record();
            if rec.is_secondary() {
                continue;
            }
            let Some(&id) = ids.get(rec.qname()) else { continue };
            // One observation per read per column, even with overlapping records.
            if bad[id as usize] || last_col[id as usize] == pos {
                continue;
            }
            last_col[id as usize] = pos;
            depth += 1;
            if let Some(bi) = base_to_idx(rec.seq()[qpos]) {
                counts[bi] += 1;
                if bi != ref_idx {
                    alt_obs.push((id, bi));
                }
            }
        }
        if depth == 0 {
            continue;
        }
        for bi in 0..4 {
            let ac = counts[bi];
            let hf = ac as f64 / depth as f64;
            if bi == ref_idx || ac < MIN_AC || hf < MIN_HF || hf >= MAX_HF {
                continue;
            }
            let mut ids_here: Vec<u32> = alt_obs.iter().filter(|o| o.1 == bi).map(|o| o.0).collect();
            ids_here.sort_unstable();
            let rev = ids_here.iter().filter(|&&r| reads[r as usize].reverse).count() as u32;
            cands.push(Candidate {
                pos,
                ref_base,
                alt_base: b"ACGT"[bi],
                ac,
                dp: depth,
                fwd: ac - rev,
                rev,
                reads: ids_here,
            });
        }
    }
    summary.candidates = cands.len();

    // Candidate pairs more than MIN_DIST apart sharing at least MIN_SHARED carriers.
    let n_c = cands.len();
    let pos_of: Vec<u32> = cands.iter().map(|c| c.pos).collect();
    let mut index: Vec<Vec<u32>> = vec![Vec::new(); reads.len()];
    for (ci, c) in cands.iter().enumerate() {
        for &r in &c.reads {
            index[r as usize].push(ci as u32);
        }
    }
    let mut pairs: Vec<(u32, u32)> = (0..n_c)
        .into_par_iter()
        .map_init(
            || (vec![0u32; n_c], Vec::<u32>::new()),
            |(cnt, touched), a| {
                let first_b = pos_of.partition_point(|&p| p <= pos_of[a] + MIN_DIST) as u32;
                for &r in &cands[a].reads {
                    let lst = &index[r as usize];
                    let from = lst.partition_point(|&b| b < first_b);
                    for &b in &lst[from..] {
                        if cnt[b as usize] == 0 {
                            touched.push(b);
                        }
                        cnt[b as usize] += 1;
                    }
                }
                let mut out = Vec::new();
                for &b in touched.iter() {
                    if cnt[b as usize] >= MIN_SHARED {
                        out.push((a as u32, b));
                    }
                    cnt[b as usize] = 0;
                }
                touched.clear();
                out
            },
        )
        .flatten()
        .collect();
    pairs.sort_unstable();
    summary.pairs_tested = pairs.len();

    // Good reads co-covering each pair, by strand: one span covering both positions (Fenwick sweep),
    // or two different spans of a split read.
    let good: Vec<&ReadInfo> = reads.iter().zip(&bad).filter(|(_, &b)| !b).map(|(r, _)| r).collect();
    let spans: Vec<(u32, u32, bool)> = good
        .iter()
        .flat_map(|r| r.spans.iter().map(move |&(s, e)| (s, e, r.reverse)))
        .collect();
    let multi: Vec<&ReadInfo> = good.iter().copied().filter(|r| r.spans.len() > 1).collect();
    let queries: Vec<(u32, u32)> = pairs.iter().map(|&(a, b)| (pos_of[a as usize], pos_of[b as usize])).collect();
    let mut cocover = cocover_counts(&spans, &queries, refseq.len() as u32);
    cocover.par_iter_mut().zip(queries.par_iter()).for_each(|(c, &(lo, hi))| {
        for r in &multi {
            if let (Some(i), Some(j)) = (span_index(&r.spans, lo), span_index(&r.spans, hi)) {
                if i != j {
                    c[r.reverse as usize] += 1;
                }
            }
        }
    });

    // Pair rule overall and on each strand, counting only carriers that co-cover both positions.
    let passing: Vec<(u32, u32, u32, f64, f64)> = pairs
        .par_iter()
        .zip(cocover.par_iter())
        .filter_map(|(&(ia, ib), &n)| {
            let (a, b) = (&cands[ia as usize], &cands[ib as usize]);
            let (lo, hi) = (a.pos, b.pos);
            let strand_counts = |ids: &[u32]| {
                let mut c = [0u32; 2];
                for &r in ids {
                    let info = &reads[r as usize];
                    if span_index(&info.spans, lo).is_some() && span_index(&info.spans, hi).is_some() {
                        c[info.reverse as usize] += 1;
                    }
                }
                c
            };
            let big_k = strand_counts(&a.reads);
            let small_n = strand_counts(&b.reads);
            let shared: Vec<u32> = a.reads.iter().copied().filter(|r| b.reads.binary_search(r).is_ok()).collect();
            let k = strand_counts(&shared);
            let both = |v: [u32; 2]| v[0] + v[1];
            let (enrich, lp) = pair_rule(both(k), both(n), both(big_k), both(small_n), false)?;
            for s in 0..2 {
                pair_rule(k[s], n[s], big_k[s], small_n[s], true)?;
            }
            Some((ia, ib, both(k), enrich, lp))
        })
        .collect();
    summary.pairs_passing = passing.len();

    // Each allele keeps its most significant partner.
    let mut best: Vec<Option<(u32, u32, f64, f64)>> = vec![None; n_c];
    for &(ia, ib, k, enrich, lp) in &passing {
        for (me, other) in [(ia, ib), (ib, ia)] {
            let slot = &mut best[me as usize];
            if slot.map_or(true, |b| lp < b.3) {
                *slot = Some((other, k, enrich, lp));
            }
        }
    }
    let calls: Vec<LinkCall> = cands
        .iter()
        .zip(&best)
        .filter_map(|(c, b)| {
            let (other, shared, enrich, log10p) = (*b)?;
            Some(LinkCall {
                pos: c.pos as usize,
                ref_allele: (c.ref_base as char).to_string(),
                alt_allele: (c.alt_base as char).to_string(),
                ac: c.ac,
                dp: c.dp,
                fwd: c.fwd,
                rev: c.rev,
                partner_name: cands[other as usize].name(),
                shared,
                enrich,
                log10p,
            })
        })
        .collect();
    summary.accepted = calls.len();
    println!(
        "Linkage rescue: {} reads ({} above error {:.3} dropped), {} candidates, {} pairs tested, {} pairs passing, {} accepted, {:.1} s",
        summary.reads,
        summary.bad_reads,
        err_cut,
        summary.candidates,
        summary.pairs_tested,
        summary.pairs_passing,
        summary.accepted,
        t0.elapsed().as_secs_f64()
    );
    Ok((calls, summary))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::call::{merge_link_calls, write_vcf, Variant};
    use rand::rngs::StdRng;
    use rand::{Rng, SeedableRng};
    use rust_htslib::bam::header::{Header, HeaderRecord};
    use rust_htslib::bam::record::{CigarString, Record};
    use rust_htslib::bam::{Format, Writer};
    use std::collections::HashMap;
    use std::path::PathBuf;

    const REF_LEN: usize = 1000;
    const N_READS: usize = 2000;

    fn fixture_ref() -> Vec<u8> {
        let mut rng = StdRng::seed_from_u64(11);
        (0..REF_LEN).map(|_| b"ACGT"[rng.random_range(0..4)]).collect()
    }

    fn alt_of(b: u8) -> u8 {
        match b {
            b'A' => b'G',
            b'G' => b'A',
            b'C' => b'T',
            _ => b'C',
        }
    }

    fn tmp_dir(tag: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!("himito_linkage_{}_{}", tag, std::process::id()));
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    fn rec(name: &str, tid: i32, pos: usize, seq: &[u8], flags: u16, nm: i32) -> Record {
        let mut rec = Record::new();
        rec.set(name.as_bytes(), Some(&CigarString(vec![Cigar::Match(seq.len() as u32)])), seq, &vec![30u8; seq.len()]);
        rec.set_tid(tid);
        rec.set_pos(pos as i64);
        rec.set_mapq(60);
        rec.set_flags(flags);
        rec.push_aux(b"NM", Aux::I32(nm)).unwrap();
        rec
    }

    /// Writes `chrm` after a few chr1 reads under header `chr1`, `chrM`: chrM is tid 1, so a
    /// `tid() == 0` assumption sees nothing.
    fn write_bam(tag: &str, chrm: &[Record]) -> PathBuf {
        let mut header = Header::new();
        for (name, len) in [("chr1", 100usize), ("chrM", REF_LEN)] {
            let mut sq = HeaderRecord::new(b"SQ");
            sq.push_tag(b"SN", name);
            sq.push_tag(b"LN", &len.to_string());
            header.push_record(&sq);
        }
        let path = tmp_dir(tag).join("in.bam");
        let mut w = Writer::from_path(&path, &header, Format::Bam).unwrap();
        for i in 0..5 {
            w.write(&rec(&format!("c1_{i}"), 0, 0, &[b'T'; 50], 0, 0)).unwrap();
        }
        for r in chrm {
            w.write(r).unwrap();
        }
        path
    }

    /// `N_READS` full-length reads alternating forward/reverse; read `i` carries allele A (at `pa`) when `a(i)`.
    fn run_scenario(
        tag: &str,
        pa: usize,
        pb: usize,
        a: impl Fn(usize) -> bool,
        b: impl Fn(usize) -> bool,
        nm: impl Fn(usize) -> i32,
    ) -> (Vec<LinkCall>, LinkSummary) {
        let refseq = fixture_ref();
        let reads: Vec<Record> = (0..N_READS)
            .map(|i| {
                let mut seq = refseq.clone();
                for (p, on) in [(pa, a(i)), (pb, b(i))] {
                    if on {
                        seq[p] = alt_of(refseq[p]);
                    }
                }
                rec(&format!("r{i}"), 1, 0, &seq, if i % 2 == 1 { 16 } else { 0 }, nm(i))
            })
            .collect();
        rescue(&write_bam(tag, &reads), &refseq, "chrM").unwrap()
    }

    #[test]
    fn log10_hyper_sf_matches_reference_vectors() {
        for (k, n_total, k_succ, n_draw, want) in [
            (53, 2635, 68, 65, -84.22),
            (17, 4040, 30, 25, -32.64),
            (11, 4000, 129, 117, -2.92),
            (2, 1000, 20, 20, -1.24),
            (1, 1000, 20, 20, -0.48),
        ] {
            let got = log10_hyper_sf(k, n_total, k_succ, n_draw);
            assert!((got - want).abs() <= 0.01, "({k},{n_total},{k_succ},{n_draw}) -> {got}, want {want}");
        }
        assert_eq!(log10_hyper_sf(0, 1000, 20, 20), 0.0);
    }

    #[test]
    fn fenwick_cocover_matches_brute_force() {
        let mut rng = StdRng::seed_from_u64(5);
        let len = 1000u32;
        let spans: Vec<(u32, u32, bool)> = (0..2000)
            .map(|_| {
                let s = rng.random_range(0..len);
                let e = rng.random_range(s + 1..=len);
                (s, e, rng.random_bool(0.5))
            })
            .collect();
        let queries: Vec<(u32, u32)> = (0..500)
            .map(|_| {
                let a = rng.random_range(0..len);
                let b = rng.random_range(0..len);
                (a.min(b), a.max(b))
            })
            .collect();
        let got = cocover_counts(&spans, &queries, len);
        for (q, &(lo, hi)) in queries.iter().enumerate() {
            let (mut f, mut r) = (0, 0);
            for &(s, e, rev) in &spans {
                if s <= lo && e > hi {
                    if rev { r += 1 } else { f += 1 }
                }
            }
            assert_eq!(got[q], [f, r], "query {lo}-{hi}");
        }
    }

    #[test]
    fn linked_on_both_strands_is_accepted() {
        // Reads 0..40: 20 forward + 20 reverse carry both alleles.
        let (calls, summary) = run_scenario("both", 200, 700, |i| i < 40, |i| i < 40, |_| 10);
        assert_eq!(summary.candidates, 2);
        assert_eq!(summary.pairs_tested, 1);
        assert_eq!(summary.pairs_passing, 1);
        assert_eq!(calls.len(), 2, "{calls:?}");
        let refseq = fixture_ref();
        let a = &calls[0];
        assert_eq!(a.pos, 200);
        assert_eq!(a.ref_allele, (refseq[200] as char).to_string());
        assert_eq!(a.alt_allele, (alt_of(refseq[200]) as char).to_string());
        assert_eq!((a.ac, a.dp, a.fwd, a.rev, a.shared), (40, 2000, 20, 20, 40));
        assert_eq!(a.partner_name, format!("m.701{}>{}", refseq[700] as char, alt_of(refseq[700]) as char));
        assert!(a.log10p <= -10.0 && a.enrich >= 10.0);
        assert_eq!(calls[1].pos, 700);
    }

    #[test]
    fn linked_on_one_strand_only_is_rejected() {
        // The 2542A/2567A case: all 40 shared carriers are forward (even indices).
        let fwd40 = |i: usize| i % 2 == 0 && i < 80;
        let (calls, summary) = run_scenario("one", 200, 700, fwd40, fwd40, |_| 10);
        assert_eq!(summary.candidates, 2);
        assert_eq!(summary.pairs_tested, 1);
        assert_eq!(summary.pairs_passing, 0);
        assert!(calls.is_empty(), "{calls:?}");
    }

    #[test]
    fn independent_alleles_are_rejected() {
        let mut rng = StdRng::seed_from_u64(3);
        let a: Vec<usize> = rand::seq::index::sample(&mut rng, N_READS, 40).into_vec();
        let b: Vec<usize> = rand::seq::index::sample(&mut rng, N_READS, 40).into_vec();
        let (calls, summary) =
            run_scenario("indep", 200, 700, |i| a.contains(&i), |i| b.contains(&i), |_| 10);
        assert_eq!(summary.candidates, 2);
        assert_eq!(summary.pairs_passing, 0);
        assert!(calls.is_empty(), "{calls:?}");
    }

    #[test]
    fn linkage_only_on_bad_reads_is_rejected() {
        // The 40 carriers are the highest-error reads, above the 95th percentile.
        let (calls, summary) =
            run_scenario("bad", 200, 700, |i| i < 40, |i| i < 40, |i| if i < 40 { 200 } else { 10 });
        assert_eq!(summary.bad_reads, 40);
        assert_eq!(summary.candidates, 0);
        assert!(calls.is_empty(), "{calls:?}");
    }

    #[test]
    fn reads_split_into_two_records_co_cover_both_alleles() {
        // Like a read split at the circular origin.
        let refseq = fixture_ref();
        let mut reads = Vec::new();
        for (start, end, alt_at, flag) in [(0usize, 300usize, 100usize, 0u16), (600, 1000, 900, 2048)] {
            for i in 0..N_READS {
                let mut seq = refseq[start..end].to_vec();
                if i < 40 {
                    seq[alt_at - start] = alt_of(refseq[alt_at]);
                }
                reads.push(rec(&format!("r{i}"), 1, start, &seq, flag | if i % 2 == 1 { 16 } else { 0 }, 3));
            }
        }
        let (calls, summary) = rescue(&write_bam("split", &reads), &refseq, "chrM").unwrap();
        assert_eq!(summary.reads, N_READS);
        assert_eq!(summary.pairs_passing, 1, "{summary:?}");
        assert_eq!(calls.len(), 2, "{calls:?}");
        assert_eq!(calls[0].dp, N_READS as u32);
    }

    #[test]
    fn alleles_within_min_dist_are_never_tested() {
        // Exactly MIN_DIST apart and perfectly linked.
        let (calls, summary) = run_scenario("dist", 200, 200 + MIN_DIST as usize, |i| i < 40, |i| i < 40, |_| 10);
        assert_eq!(summary.candidates, 2);
        assert_eq!(summary.pairs_tested, 0);
        assert!(calls.is_empty(), "{calls:?}");
    }

    fn snp(pos: usize, r: &str, a: &str, filter: Option<&str>) -> Variant {
        Variant {
            pos,
            ref_allele: r.into(),
            alt_allele: a.into(),
            variant_type: "SNP".into(),
            allele_count: 30,
            filter: filter.map(String::from),
        }
    }

    fn link(pos: usize, r: &str, a: &str) -> LinkCall {
        LinkCall {
            pos,
            ref_allele: r.into(),
            alt_allele: a.into(),
            ac: 12,
            dp: 1000,
            fwd: 6,
            rev: 6,
            partner_name: "m.901C>T".into(),
            shared: 9,
            enrich: 40.0,
            log10p: -14.5,
        }
    }

    #[test]
    fn already_pass_allele_is_not_reemitted_and_filtered_graph_snv_is_cleared() {
        let graph = vec![
            snp(199, "A", "G", None),
            snp(699, "C", "T", Some("Strand_bias")),
            snp(799, "G", "A", Some("Strand_bias")),
        ];
        let links = vec![link(199, "A", "G"), link(400, "T", "C"), link(699, "C", "T")];
        let (merged, new_calls) = merge_link_calls(&graph, &links);
        assert_eq!(new_calls, vec![link(400, "T", "C")]);
        assert_eq!(merged.len(), 3);
        assert_eq!(merged[0].filter, None);
        assert_eq!(merged[1].filter, None, "linked filtered SNV gets FILTER cleared");
        assert_eq!(merged[2].filter.as_deref(), Some("Strand_bias"), "unlinked stays flagged");

        let path = tmp_dir("vcf").join("out.vcf");
        let coverage: HashMap<usize, usize> = [(199, 100), (699, 100), (799, 100)].into_iter().collect();
        write_vcf(&merged, &coverage, &path, "s1", "chrM", REF_LEN, 0.015, true, Some(&new_calls), &std::collections::HashSet::new()).unwrap();
        let vcf = std::fs::read_to_string(&path).unwrap();
        assert!(vcf.contains("##INFO=<ID=LINK,"), "{vcf}");
        let body: Vec<&str> = vcf.lines().filter(|l| !l.starts_with('#')).collect();
        let pos: Vec<&str> = body.iter().map(|l| l.split('\t').nth(1).unwrap()).collect();
        assert_eq!(pos, vec!["200", "401", "700", "800"], "{vcf}");
        assert_eq!(body[1], "chrM\t401\t.\tT\tC\t.\tPASS\tDP=1000;LINK=m.901C>T,-14.50\tGT:AD:HF\t1:12:0.012");
        assert!(body[2].contains("\tPASS\t"), "{}", body[2]);

        // Without a linkage run the header has no LINK line.
        write_vcf(&merged, &coverage, &path, "s1", "chrM", REF_LEN, 0.015, true, None, &std::collections::HashSet::new()).unwrap();
        assert!(!std::fs::read_to_string(&path).unwrap().contains("ID=LINK"));
    }
}
