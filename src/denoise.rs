use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};
use anyhow::{Context, Result};
use bio::io::fasta;
use log::{info, warn};
use rust_htslib::bam::{self, Read, Reader};
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::record::CigarString;

use crate::denoise_indel::{rewrite_read, strip_stale_tags};

use crate::denoise_indel as indel;

pub(crate) const NQ: usize = 94; // Phred qualities 0..=93

/// Ceiling on the per-base error probability.
pub(crate) const E_MAX: f64 = 0.75;

#[inline]
pub(crate) fn phred_to_e(q: u8) -> f64 {
    let q = q.min((NQ - 1) as u8);
    10f64.powf(-(q as f64) / 10.0).min(E_MAX)
}

#[inline]
fn base_to_idx(b: u8) -> Option<usize> {
    match b.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

#[inline]
fn idx_to_base(i: usize) -> u8 {
    [b'A', b'C', b'G', b'T'][i]
}

pub(crate) struct ColumnCounts {
    n: [[u32; NQ]; 4],     // n[allele][quality]
    strand: [[u32; 2]; 4], // strand[allele][0=fwd, 1=rev]
    ref_idx: Option<usize>,
}

impl ColumnCounts {
    fn new(ref_idx: Option<usize>) -> Self {
        ColumnCounts { n: [[0; NQ]; 4], strand: [[0; 2]; 4], ref_idx }
    }

    fn add(&mut self, allele: usize, q: u8, reverse: bool) {
        let q = q.min((NQ - 1) as u8) as usize;
        self.n[allele][q] += 1;
        self.strand[allele][if reverse { 1 } else { 0 }] += 1;
    }

    fn non_ref(&self) -> u32 {
        (0..4)
            .filter(|&a| Some(a) != self.ref_idx)
            .map(|a| self.strand[a][0] + self.strand[a][1])
            .sum()
    }

    /// (fwd, rev) observation totals across all four alleles in this column.
    fn strand_totals(&self) -> (u32, u32) {
        (0..4).fold((0, 0), |(f, r), a| (f + self.strand[a][0], r + self.strand[a][1]))
    }
}

pub(crate) struct SiteModel {
    freq: [f64; 4],
    kept: [bool; 4],
    /// Not kept; `compute_corrections` defers their correction until linkage to a kept minor allele is tested.
    pending: [bool; 4],
    /// Non-ref alleles that cleared the per-strand count gate but were rejected as strand-biased.
    strand_biased_alts: usize,
}

impl SiteModel {
    fn kept_alt_count(&self, ref_idx: Option<usize>) -> usize {
        (0..4)
            .filter(|&a| self.kept[a] && Some(a) != ref_idx)
            .count()
    }
}

fn fit_site(
    c: &ColumnCounts,
    min_strand: u32,
    vaf: f64,
    strand_bias_p: f64,
    homoplasmic_vaf: f64,
) -> SiteModel {
    // The null for the strand-bias test is this column's OWN forward fraction, not a library-wide constant.
    let (col_fwd, col_rev) = c.strand_totals();
    let col_total = col_fwd + col_rev;
    let expected_fwd_frac = if col_total == 0 {
        0.5
    } else {
        col_fwd as f64 / col_total as f64
    };

    // Eligibility: reference always eligible.
    let mut eligible = [false; 4];
    let mut strand_biased_alts = 0usize;
    for a in 0..4 {
        if Some(a) == c.ref_idx {
            eligible[a] = true;
            continue;
        }
        let (fwd, rev) = (c.strand[a][0], c.strand[a][1]);
        if !(fwd >= min_strand && rev >= min_strand) {
            continue;
        }
        // Mirrors call::filter_strand_bias's near-homoplasmic exemption.
        let prop = if col_total == 0 {
            0.0
        } else {
            (fwd + rev) as f64 / col_total as f64
        };
        if strand_bias_p > 0.0
            && prop < homoplasmic_vaf
            && crate::call::strand_bias_pvalue(fwd as usize, rev as usize, expected_fwd_frac)
                < strand_bias_p
            && (fwd as f64 / (fwd + rev) as f64 - expected_fwd_frac).abs()
                >= crate::STRAND_BIAS_MIN_EFFECT
        {
            strand_biased_alts += 1;
            continue;
        }
        eligible[a] = true;
    }

    // Precompute per-quality error and total observations.
    let e: [f64; NQ] = std::array::from_fn(|q| phred_to_e(q as u8));
    let total_obs: f64 = c.n.iter().flatten().map(|&x| x as f64).sum();

    // Initialize f from raw proportions restricted to eligible alleles.
    let mut f = [0.0f64; 4];
    let mut init_sum = 0.0;
    for a in 0..4 {
        if eligible[a] {
            let t: u32 = c.n[a].iter().sum();
            f[a] = t as f64;
            init_sum += t as f64;
        }
    }
    if init_sum <= 0.0 {
        // No eligible observations: fall back to reference-only.
        if let Some(r) = c.ref_idx {
            f[r] = 1.0;
        }
    } else {
        for a in 0..4 {
            f[a] /= init_sum;
        }
    }

    // EM.
    if total_obs > 0.0 {
        for _iter in 0..100 {
            let mut newf = [0.0f64; 4];
            for a in 0..4 {
                for qi in 0..NQ {
                    let cnt = c.n[a][qi];
                    if cnt == 0 {
                        continue;
                    }
                    let eq = e[qi];
                    // p_k ∝ f_k * P(obs=a | true=k)
                    let mut p = [0.0f64; 4];
                    let mut psum = 0.0;
                    for k in 0..4 {
                        if !eligible[k] || f[k] <= 0.0 {
                            continue;
                        }
                        let like = if k == a { 1.0 - eq } else { eq / 3.0 };
                        p[k] = f[k] * like;
                        psum += p[k];
                    }
                    if psum <= 0.0 {
                        continue;
                    }
                    let cntf = cnt as f64;
                    for k in 0..4 {
                        if p[k] > 0.0 {
                            newf[k] += cntf * p[k] / psum;
                        }
                    }
                }
            }
            for k in 0..4 {
                newf[k] /= total_obs;
            }
            let delta: f64 = (0..4).map(|k| (newf[k] - f[k]).abs()).fold(0.0, f64::max);
            f = newf;
            if delta < 1e-6 {
                break;
            }
        }
    }

    // Keep rule: reference always kept; a non-ref eligible allele kept iff f >= vaf.
    let mut kept = [false; 4];
    let mut pending = [false; 4];
    for a in 0..4 {
        if Some(a) == c.ref_idx {
            kept[a] = true;
        } else if eligible[a] && f[a] >= vaf {
            kept[a] = true;
        } else if eligible[a] && f[a] >= crate::PENDING_FLOOR {
            pending[a] = true;
        }
    }
    // Renormalize freq over kept alleles (guards map_allele against dropped mass).
    let ksum: f64 = (0..4).filter(|&a| kept[a]).map(|a| f[a]).sum();
    let mut freq = [0.0f64; 4];
    if ksum > 0.0 {
        for a in 0..4 {
            if kept[a] {
                freq[a] = f[a] / ksum;
            }
        }
    } else if let Some(r) = c.ref_idx {
        freq[r] = 1.0;
    }

    SiteModel { freq, kept, pending, strand_biased_alts }
}

fn map_allele(observed: usize, q: u8, m: &SiteModel) -> usize {
    // An observation of an allele the site model already kept is evidence FOR that allele, so it stands as observed.
    if m.kept.get(observed).copied().unwrap_or(false) {
        return observed;
    }

    let eq = phred_to_e(q);
    let mut best = observed;
    let mut best_val = -1.0f64;
    for k in 0..4 {
        if !m.kept[k] {
            continue;
        }
        let like = if k == observed { 1.0 - eq } else { eq / 3.0 };
        let val = m.freq[k] * like;
        if val > best_val {
            best_val = val;
            best = k;
        }
    }
    best
}

pub type Corrections = HashMap<Vec<u8>, indel::ReadEdits>;

/// Decided indel sites, per contig.
pub type IndelModels = HashMap<Vec<u8>, Vec<(u32, indel::IndelSiteModel)>>;

#[derive(Default, Debug, serde::Serialize)]
pub struct DenoiseStats {
    pub reads_processed: u64,
    pub reads_modified: u64,
    pub bases_examined: u64,
    pub bases_modified: u64,
    pub columns_processed: u64,
    pub columns_corrected: u64,
    pub columns_skipped: u64,
    pub alt_sites_preserved: u64,
    pub alt_alleles_strand_biased: u64,
    pub alt_alleles_rescued_by_linkage: u64,
    pub substitution_matrix: [[u64; 4]; 4], // [from][to]
    // --- indel denoising (all zero unless --indels is enabled) ---
    pub indel_sites_examined: u64,
    /// Non-reference ALLELES (not sites) that survived candidacy across all sites.
    pub indel_alt_alleles_kept: u64,
    /// Non-reference ALLELES (not sites) rejected specifically as strand-biased.
    pub indel_alt_alleles_strand_biased: u64,
    pub indel_events_examined: u64,
    /// Votes dropped in pass 1: the read has no aligned query base at the normalized column.
    pub indel_events_out_of_span: u64,
    /// The gate is applied at the LEFT-NORMALIZED position.
    /// Tandem repeats whose unit is longer than one base are OUT OF SCOPE here and get no span
    /// requirement beyond the event position.
    pub indel_events_not_spanning_repeat: u64,
    /// Votes dropped in pass 1: the same read already cast an ALT vote at this normalized site.
    pub indel_events_duplicate_site: u64,
    /// A column no such event reached is ordinary missing coverage and is not counted here.
    /// Booked after the pileup from keys remembered at the gate, not from `sites`: a site is inserted
    /// there only after a vote clears the same gate depth uses, so this state never appears in `sites`.
    pub indel_sites_dropped_no_spanning_depth: u64,
    pub reads_reassigned_ref_to_ins: u64,
    pub reads_reassigned_ref_to_del: u64,
    pub reads_reassigned_indel_to_ref: u64,
    pub reads_reassigned_indel_to_indel: u64,
    pub indel_bases_inserted: u64,
    pub indel_bases_removed: u64,
    /// Index = repeat-context length L; audits how much correction fired in long runs.
    pub reassignments_by_hp_length: Vec<u64>,
    /// Zero unless `--indel-hp-length-model` is on.
    pub subs_dropped_in_hp_rewrite: u64,
    /// (read, site) pairs skipped: the read lacks the match-context flank the rewrite walk needs.
    pub assignments_skipped_span_guard: u64,
    /// (read, site) pairs skipped: the MAP target needs an unsupported transition.
    pub assignments_skipped_unsupported_transition: u64,
    /// (read, site) pairs skipped: an untouchable oversized indel of the read lies too close to the site.
    pub assignments_skipped_oversized_neighbor: u64,
    /// (read, site) pairs skipped: two CIGAR events of the read normalize to this site (ambiguous).
    pub assignments_skipped_ambiguous_carry: u64,
    /// (read, site) PAIRS skipped because an earlier.
    pub assignments_skipped_swallowed_by_gained_deletion: u64,
    /// DIAGNOSTIC, should normally be zero.
    pub indel_edits_emitted_but_not_applied: u64,
}

fn compute_corrections(
    bam_path: &Path,
    refs: &HashMap<Vec<u8>, Vec<u8>>,
    min_strand: u32,
    vaf: f64,
    strand_bias_p: f64,
    homoplasmic_vaf: f64,
    iopts: &indel::IndelOpts,
) -> Result<(Corrections, IndelModels, DenoiseStats)> {
    let mut reader = Reader::from_path(bam_path)
        .with_context(|| format!("cannot open BAM {bam_path:?}"))?;
    // Own the tid->contig-name map so we can borrow the reader mutably for pileup.
    let header = reader.header().to_owned();

    let mut corrections: Corrections = HashMap::new();
    let mut stats = DenoiseStats::default();

    // ALT event counts per left-normalized site.
    let mut sites: HashMap<(Vec<u8>, u32), indel::IndelSite> = HashMap::new();
    // Per-column (depth, fwd, rev).
    let mut depths: HashMap<(Vec<u8>, u32), (u32, u32, u32)> = HashMap::new();
    // Remembered here because `sites` only receives votes that cleared that gate.
    let mut span_gated_sites: HashSet<(Vec<u8>, u32)> = HashSet::new();

    // `deferred` holds each pending allele's would-be corrections.
    type AlleleKey = (u32, u32, usize);
    let mut deferred: HashMap<AlleleKey, Vec<(Vec<u8>, u32, u8)>> = HashMap::new();
    let mut carriers: HashMap<AlleleKey, HashSet<usize>> = HashMap::new();
    let mut kept_minors: HashSet<AlleleKey> = HashSet::new();
    let mut read_idx: HashMap<Vec<u8>, usize> = HashMap::new();
    let mut spans: Vec<(u32, u32, u32)> = Vec::new();
    let mut corrected_cols: HashSet<(u32, u32)> = HashSet::new();

    let mut plp = reader.pileup();
    plp.set_max_depth(1_000_000);
    for p in plp {
        let col = p.context("pileup error")?;
        let name = header.tid2name(col.tid());
        let refseq = match refs.get(name) {
            Some(s) => s,
            None => continue,
        };
        let refpos = col.pos() as usize;
        let ref_idx = refseq.get(refpos).and_then(|&b| base_to_idx(b));

        let mut counts = ColumnCounts::new(ref_idx);
        // Remember each usable observation to correct after the fit.
        let mut obs: Vec<(Vec<u8>, usize, u8, usize)> = Vec::new(); // (qname, allele, qual, read idx)
        let tid = col.tid();

        // Indel collection MUST happen here, before the `non_ref() == 0` skip below.
        let mut col_depth = 0u32;
        let mut col_fwd = 0u32;
        let mut col_rev = 0u32;
        // (qname, allele, reverse, read's own POS, read's own CIGAR, read's exclusive reference end).
        let mut col_events: Vec<(Vec<u8>, indel::Allele, bool, i64, CigarString, i64)> =
            Vec::new();
        // A read that ends inside the run cannot show the run's length.
        let tract_need_end = if iopts.enabled {
            indel::hp_tract_end(refseq, refpos as u32) as i64 + 1
        } else {
            0
        };
        // Then the column casts one LENGTH vote per spanning read.
        let hp_votes_here = if iopts.enabled && iopts.hp_length_model {
            indel::hp_site(refseq, refpos as u32)
        } else {
            None
        };
        // (qname, length allele, reverse) for this column's HP-length ALT votes.
        let mut hp_col_votes: Vec<(Vec<u8>, indel::Allele, bool)> = Vec::new();
        // Kept apart from the ALT votes because reference support is not an entry in `IndelSite::obs`.
        let mut hp_col_ref: Vec<Vec<u8>> = Vec::new();

        for a in col.alignments() {
            let rec = a.record();
            if rec.is_secondary() || rec.is_supplementary() || rec.is_unmapped() {
                continue;
            }

            if iopts.enabled && !a.is_refskip() {
                if let Some(q) = a.qpos() {
                    // This decides only whether the read joins THIS column's depth.
                    let reference_end = rec.reference_end();
                    if reference_end >= tract_need_end {
                        col_depth += 1;
                        if rec.is_reverse() { col_rev += 1 } else { col_fwd += 1 }

                        // `hp_allele`'s span gate is one flank STRICTER than `tract_need_end`.
                        if let Some((hp_e, b)) = hp_votes_here {
                            let cigar = rec.cigar().take();
                            let seq = rec.seq().as_bytes();
                            match indel::hp_allele(
                                rec.pos(), &cigar, &seq, refseq, refpos as u32, hp_e,
                            ) {
                                Some(0) => hp_col_ref.push(rec.qname().to_vec()),
                                Some(k) => hp_col_votes.push((
                                    rec.qname().to_vec(),
                                    indel::hp_length_allele(k, b),
                                    rec.is_reverse(),
                                )),
                                None => {}
                            }
                        }
                    }
                    // Decide the event's allele first, WITHOUT touching the CIGAR.
                    let ev_allele = match a.indel() {
                        bam::pileup::Indel::Ins(n) if n < iopts.max_len => {
                            let s: Vec<u8> =
                                (0..n as usize).map(|i| rec.seq()[q + 1 + i]).collect();
                            Some(indel::Allele::Ins(s))
                        }
                        bam::pileup::Indel::Del(n) if n < iopts.max_len => {
                            Some(indel::Allele::Del(n))
                        }
                        // Indel::None contributes depth only; an indel at or above max_len contributes depth and is never corrected.
                        _ => None,
                    };
                    if let Some(allele) = ev_allele {
                        stats.indel_events_examined += 1;
                        col_events.push((
                            rec.qname().to_vec(),
                            allele,
                            rec.is_reverse(),
                            rec.pos(),
                            rec.cigar().take(),
                            reference_end,
                        ));
                    }
                }
            }

            if a.is_del() || a.is_refskip() {
                continue;
            }
            let qpos = match a.qpos() {
                Some(x) => x,
                None => continue,
            };
            let base = rec.seq()[qpos];
            let allele = match base_to_idx(base) {
                Some(i) => i,
                None => continue,
            };
            let qual = rec.qual()[qpos];
            counts.add(allele, qual, rec.is_reverse());
            let rp = refpos as u32;
            let ri = match read_idx.get(rec.qname()) {
                Some(&i) => {
                    let s = &mut spans[i];
                    if s.0 == tid {
                        s.1 = s.1.min(rp);
                        s.2 = s.2.max(rp);
                    }
                    i
                }
                None => {
                    read_idx.insert(rec.qname().to_vec(), spans.len());
                    spans.push((tid, rp, rp));
                    spans.len() - 1
                }
            };
            obs.push((rec.qname().to_vec(), allele, qual, ri));
        }

        if iopts.enabled {
            depths.insert((name.to_vec(), refpos as u32), (col_depth, col_fwd, col_rev));
            // `add_from_read` still enforces one ALT vote per read per site.
            if !hp_col_votes.is_empty() {
                let site = sites.entry((name.to_vec(), refpos as u32)).or_default();
                for (qname, allele, reverse) in hp_col_votes {
                    if !site.add_from_read(&qname, allele, reverse) {
                        stats.indel_events_duplicate_site += 1;
                    }
                }
                for qname in hp_col_ref {
                    site.add_hp_ref(&qname);
                }
            }
            for (qname, allele, reverse, rpos, cigar, reference_end) in col_events {
                let (norm_pos, norm_allele) =
                    indel::normalize_left(refseq, refpos as u32, &allele);
                // Dropped silently: the read has already spoken at this site.
                if iopts.hp_length_model && indel::hp_site(refseq, norm_pos).is_some() {
                    continue;
                }
                // `normalize_left` walks the reference with no knowledge of read boundaries.
                if indel::match_run_len(rpos, &cigar, norm_pos) == 0 {
                    stats.indel_events_out_of_span += 1;
                    continue;
                }
                // The two anchors really do disagree: `hp_tract_end` reports `anchor + 1` when only a single base follows.
                if reference_end < indel::hp_tract_end(refseq, norm_pos) as i64 + 1 {
                    stats.indel_events_not_spanning_repeat += 1;
                    // Counted after the pileup, once depth is final.
                    span_gated_sites.insert((name.to_vec(), norm_pos));
                    continue;
                }
                // A read with two events normalizing to this site contributes at most one ALT vote.
                let voted = sites
                    .entry((name.to_vec(), norm_pos))
                    .or_default()
                    .add_from_read(&qname, norm_allele, reverse);
                if !voted {
                    stats.indel_events_duplicate_site += 1;
                }
            }
        }

        stats.columns_processed += 1;
        stats.bases_examined += obs.len() as u64;

        if counts.non_ref() == 0 {
            stats.columns_skipped += 1;
            continue;
        }

        let model = fit_site(&counts, min_strand, vaf, strand_bias_p, homoplasmic_vaf);
        stats.alt_sites_preserved += model.kept_alt_count(ref_idx) as u64;
        stats.alt_alleles_strand_biased += model.strand_biased_alts as u64;

        for a in 0..4 {
            if model.kept[a] && Some(a) != ref_idx && model.freq[a] < homoplasmic_vaf {
                kept_minors.insert((tid, refpos as u32, a));
            }
        }

        let mut changed = false;
        for (qname, allele, qual, ri) in obs {
            let corr = map_allele(allele, qual, &model);
            let key = (tid, refpos as u32, allele);
            if model.pending[allele] || kept_minors.contains(&key) {
                carriers.entry(key).or_default().insert(ri);
            }
            if model.pending[allele] {
                deferred.entry(key).or_default().push((qname, refpos as u32, idx_to_base(corr)));
                continue;
            }
            if corr != allele {
                corrections
                    .entry(qname)
                    .or_default()
                    .subs
                    .push((refpos as u32, idx_to_base(corr)));
                // KNOWN LIMITATION: these two counters are booked eagerly here, at decision time.
                stats.bases_modified += 1;
                stats.substitution_matrix[allele][corr] += 1;
                changed = true;
            }
        }
        if changed {
            corrected_cols.insert((tid, refpos as u32));
        }
    }

    // ---- ONT-recall Phase2 T6b: linkage rescue of pending alleles ----
    let rescued = rescue_linked_pending(&deferred, &carriers, &kept_minors, &spans);
    stats.alt_alleles_rescued_by_linkage = rescued.len() as u64;
    let mut pending_keys: Vec<&AlleleKey> = deferred.keys().collect();
    pending_keys.sort();
    for key in pending_keys {
        if rescued.contains(key) {
            continue;
        }
        let (tid, pos, allele) = *key;
        for (qname, refpos, base) in &deferred[key] {
            corrections.entry(qname.clone()).or_default().subs.push((*refpos, *base));
            stats.bases_modified += 1;
            stats.substitution_matrix[allele][base_to_idx(*base).unwrap()] += 1;
        }
        if !deferred[key].is_empty() {
            corrected_cols.insert((tid, pos));
        }
    }
    stats.columns_corrected = corrected_cols.len() as u64;
    // Deferred substitutions land after the pileup, so subs must be re-sorted whether or not indels are enabled.
    for edits in corrections.values_mut() {
        edits.subs.sort_by_key(|&(p, _)| p);
    }

    // Missing `depths` and an explicit zero are the same outcome: no read spanned the run.
    for key in &span_gated_sites {
        let depth = depths.get(key).map(|(d, _, _)| *d).unwrap_or(0);
        if depth == 0 {
            stats.indel_sites_dropped_no_spanning_depth += 1;
        }
    }

    // ---- deferred indel site decisions ----
    let mut models: IndelModels = HashMap::new();
    if iopts.enabled {
        let n_sites = sites.len();
        for ((contig, norm_pos), mut site) in sites.into_iter() {
            let refseq = match refs.get(&contig) {
                Some(s) => s,
                None => continue,
            };
            // Fill the denominator from the column totals: each spanning read once.
            let (depth, fwd, rev) = depths
                .get(&(contig.clone(), norm_pos))
                .copied()
                .unwrap_or((0, 0, 0));
            if depth == 0 {
                // Nothing spans this site, so there is nothing to decide.
                continue;
            }
            site.depth = depth;
            site.col_fwd = fwd;
            site.col_rev = rev;

            stats.indel_sites_examined += 1;
            let model = indel::fit_indel_site(
                &site, refseq, norm_pos, min_strand, strand_bias_p, homoplasmic_vaf, n_sites,
                iopts,
            );
            stats.indel_alt_alleles_strand_biased += model.strand_biased as u64;
            // `kept` includes REF; count only surviving non-reference alleles.
            stats.indel_alt_alleles_kept +=
                model.kept.iter().filter(|(a, _)| *a != indel::Allele::Ref).count() as u64;
            models.entry(contig).or_default().push((norm_pos, model));
        }
        // The write-back pass binary-searches each contig's list.
        for list in models.values_mut() {
            list.sort_by_key(|(p, _)| *p);
        }
    }

    Ok((corrections, models, stats))
}

const RESCUE_MIN_DIST: u32 = 20;
const RESCUE_MIN_SHARED: u32 = 3;
const RESCUE_ALPHA: f64 = 1e-4;
const RESCUE_ROUNDS: usize = 3;

/// Partners start as the kept minor alleles; each round adds what the previous one rescued.
fn rescue_linked_pending<V>(
    pending: &HashMap<(u32, u32, usize), V>,
    carriers: &HashMap<(u32, u32, usize), HashSet<usize>>,
    kept_minors: &HashSet<(u32, u32, usize)>,
    spans: &[(u32, u32, u32)],
) -> HashSet<(u32, u32, usize)> {
    let empty = HashSet::new();
    let mut partners: Vec<(u32, u32, usize)> = kept_minors.iter().copied().collect();
    partners.sort();
    let mut open: Vec<(u32, u32, usize)> = pending.keys().copied().collect();
    open.sort();
    let mut rescued = HashSet::new();

    for _ in 0..RESCUE_ROUNDS {
        let mut newly = Vec::new();
        for &p in &open {
            let (ptid, ppos, _) = p;
            let pc = carriers.get(&p).unwrap_or(&empty);
            let over_p: Vec<usize> = (0..spans.len())
                .filter(|&r| spans[r].0 == ptid && spans[r].1 <= ppos && spans[r].2 >= ppos)
                .collect();
            let mut n_partners = 0u32;
            let mut min_log10_p = f64::INFINITY;
            for &m in &partners {
                let (mtid, mpos, _) = m;
                if mtid != ptid || mpos.abs_diff(ppos) <= RESCUE_MIN_DIST {
                    continue;
                }
                let mc = carriers.get(&m).unwrap_or(&empty);
                let (lo, hi) = (ppos.min(mpos), ppos.max(mpos));
                let (mut n, mut k_p, mut k_m, mut a) = (0u32, 0u32, 0u32, 0u32);
                for &r in &over_p {
                    if spans[r].1 > lo || spans[r].2 < hi {
                        continue;
                    }
                    n += 1;
                    let (in_p, in_m) = (pc.contains(&r), mc.contains(&r));
                    k_p += in_p as u32;
                    k_m += in_m as u32;
                    a += (in_p && in_m) as u32;
                }
                if n == 0 {
                    continue;
                }
                n_partners += 1;
                if a >= RESCUE_MIN_SHARED {
                    let lp = crate::linkage::log10_hyper_sf(a, n, k_p, k_m);
                    min_log10_p = min_log10_p.min(lp);
                }
            }
            if n_partners > 0
                && min_log10_p + (n_partners as f64).log10() <= RESCUE_ALPHA.log10()
            {
                newly.push(p);
            }
        }
        if newly.is_empty() {
            break;
        }
        open.retain(|p| !newly.contains(p));
        for p in newly {
            rescued.insert(p);
            partners.push(p);
        }
        partners.sort();
    }
    rescued
}

/// Book the read and base counters for ONE reassignment that was actually applied.
/// Callers must only reach this for an edit the rewrite reports it really made; see
/// `RewriteResult::applied`.
fn book_reassignment(stats: &mut DenoiseStats, from: &indel::Allele, to: &indel::Allele) {
    match (from, to) {
        (indel::Allele::Ref, indel::Allele::Ins(s)) => {
            stats.reads_reassigned_ref_to_ins += 1;
            stats.indel_bases_inserted += s.len() as u64;
        }
        (indel::Allele::Ref, indel::Allele::Del(n)) => {
            stats.reads_reassigned_ref_to_del += 1;
            stats.indel_bases_removed += *n as u64;
        }
        (indel::Allele::Ins(from_s), indel::Allele::Ref) => {
            // Reverting a spurious insertion removes the bases it added -- the single most common correction this feature makes.
            stats.reads_reassigned_indel_to_ref += 1;
            stats.indel_bases_removed += from_s.len() as u64;
        }
        (indel::Allele::Del(from_n), indel::Allele::Ref) => {
            // Reverting a spurious deletion restores the reference bases it dropped.
            stats.reads_reassigned_indel_to_ref += 1;
            stats.indel_bases_inserted += *from_n as u64;
        }
        (indel::Allele::Ins(from_s), indel::Allele::Ins(to_s)) => {
            stats.reads_reassigned_indel_to_indel += 1;
            if to_s.len() > from_s.len() {
                stats.indel_bases_inserted += (to_s.len() - from_s.len()) as u64;
            } else if from_s.len() > to_s.len() {
                stats.indel_bases_removed += (from_s.len() - to_s.len()) as u64;
            }
        }
        (indel::Allele::Del(from_n), indel::Allele::Del(to_n)) => {
            stats.reads_reassigned_indel_to_indel += 1;
            // shrinking a deletion splices reference bases BACK INTO the read, e.g. Del(3) -> Del(1) restores 2 bases.
            if to_n < from_n {
                stats.indel_bases_inserted += (from_n - to_n) as u64;
            } else if to_n > from_n {
                // `rewrite_hp_run` replaces the whole run at once.
                stats.indel_bases_removed += (to_n - from_n) as u64;
            }
        }
        // Both the read's own indel and the target's are charged.
        (indel::Allele::Ins(from_s), indel::Allele::Del(to_n)) => {
            stats.reads_reassigned_indel_to_indel += 1;
            stats.indel_bases_removed += from_s.len() as u64 + *to_n as u64;
        }
        (indel::Allele::Del(from_n), indel::Allele::Ins(to_s)) => {
            stats.reads_reassigned_indel_to_indel += 1;
            stats.indel_bases_inserted += *from_n as u64 + to_s.len() as u64;
        }
        // Ref -> Ref: not a reassignment at all.
        (indel::Allele::Ref, indel::Allele::Ref) => {}
    }
}

fn apply_corrections(
    in_bam: &Path,
    out_bam: &Path,
    corrections: &Corrections,
    refs: &HashMap<Vec<u8>, Vec<u8>>,
    models: &IndelModels,
    iopts: &indel::IndelOpts,
    stats: &mut DenoiseStats,
) -> Result<()> {
    let mut reader = Reader::from_path(in_bam)
        .with_context(|| format!("cannot open BAM {in_bam:?}"))?;
    let header_view = reader.header().to_owned();
    let header = bam::Header::from_template(reader.header());
    let mut writer = bam::Writer::from_path(out_bam, &header, bam::Format::Bam)
        .with_context(|| format!("cannot create BAM {out_bam:?}"))?;

    let mut reads_processed = 0u64;
    let mut reads_modified = 0u64;
    // Contigs referenced by a correctable, edit-bearing record but absent from `refs`.
    let mut missing_ref_contigs: std::collections::BTreeSet<Vec<u8>> = std::collections::BTreeSet::new();
    let mut missing_ref_reads = 0u64;

    for result in reader.records() {
        let mut rec = result.context("read BAM record")?;
        // Only primary mapped reads are candidates for correction.
        let correctable =
            !(rec.is_secondary() || rec.is_supplementary() || rec.is_unmapped());
        if correctable {
            reads_processed += 1;
            // A mapped-but-contig-less record (tid == -1) is a degenerate BAM state.
            if rec.tid() >= 0 {
                let contig = header_view.tid2name(rec.tid() as u32).to_vec();
                match refs.get(&contig) {
                    Some(refseq) => {
                        let mut edits = corrections.get(rec.qname()).cloned().unwrap_or_default();
                        // (ref_pos, kind_of(from)) -> the repeat-context length that decided this edit.
                        let mut edit_context_len: HashMap<(u32, indel::EditKind), usize> = HashMap::new();
                        // They are deliberately NOT `IndelEdit`s: `rewrite_read` edits one.
                        let mut hp_moves: Vec<(u32, u32, u8, indel::Allele, indel::Allele)> =
                            Vec::new();

                        if iopts.enabled {
                            // `end_pos()` is only defined on `CigarStringView`, which carries the record's POS.
                            let view = rec.cigar();
                            let rstart = rec.pos();
                            let rend = view.end_pos();
                            let cigar = view.take();
                            let seq = rec.seq().as_bytes();
                            let evs = indel::read_events(rstart, &cigar, &seq, refseq, iopts.max_len);

                            // What this read carries at each normalized site -- the allele.
                            let mut carried: HashMap<u32, (indel::Allele, u32)> = HashMap::new();
                            let mut blocked: Vec<i64> = Vec::new();
                            // Sites where this read's OWN CIGAR carries two DIFFERENT events that both normalize here.
                            let mut ambiguous: std::collections::HashSet<u32> = std::collections::HashSet::new();
                            for e in &evs {
                                match &e.allele {
                                    Some(a) => {
                                        if carried.insert(e.norm_pos, (a.clone(), e.anchor)).is_some() {
                                            ambiguous.insert(e.norm_pos);
                                        }
                                    }
                                    None => blocked.push(e.norm_pos as i64),
                                }
                            }

                            // Fixed flank (>= 1 by construction; see `INDEL_FLANK`).
                            let flank = indel::INDEL_FLANK as i64;
                            // A gained deletion consumes up to max_len reference bases to the right.
                            let reach = flank + iopts.max_len as i64;
                            // Furthest reference position already consumed by a GAINED deletion emitted earlier in this loop for this same read.
                            let mut consumed_until: i64 = -1;

                            if let Some(list) = models.get(&contig) {
                                let lo = list.partition_point(|(p, _)| (*p as i64) < rstart);
                                for (pos, model) in &list[lo..] {
                                    let s = *pos as i64;
                                    // `rend` is EXCLUSIVE, so `s == rend` is a site the read does not span at all.
                                    if s >= rend {
                                        break;
                                    }
                                    // The tract is keyed by LENGTH, so the read's allele there comes from `hp_allele`, as in pass 1.
                                    if let Some((hp_e, hp_b)) = model.hp_run {
                                        let observed_k = match indel::hp_allele(
                                            rstart, &cigar, &seq, refseq, *pos, hp_e,
                                        ) {
                                            Some(k) => k,
                                            None => {
                                                stats.assignments_skipped_span_guard += 1;
                                                continue;
                                            }
                                        };
                                        if s - flank < rstart {
                                            stats.assignments_skipped_span_guard += 1;
                                            continue;
                                        }
                                        // An earlier gained deletion for this read already consumed the anchor or part of the run.
                                        if s <= consumed_until {
                                            stats
                                                .assignments_skipped_swallowed_by_gained_deletion += 1;
                                            continue;
                                        }
                                        if blocked.iter().any(|b| (b - s).abs() <= reach) {
                                            stats.assignments_skipped_oversized_neighbor += 1;
                                            continue;
                                        }
                                        let observed = indel::hp_length_allele(observed_k, hp_b);
                                        let target = match indel::assign_allele(
                                            &observed, model, iopts,
                                        ) {
                                            indel::Assignment::Keep => continue,
                                            // Unreachable at an HP site (T8c), but counted rather than silently dropped if that ever changes.
                                            indel::Assignment::Refused => {
                                                stats
                                                    .assignments_skipped_unsupported_transition += 1;
                                                continue;
                                            }
                                            indel::Assignment::Move(t) => t,
                                        };
                                        hp_moves.push((*pos, hp_e, hp_b, observed, target));
                                        continue;
                                    }
                                    // this read's own CIGAR carries two DIFFERENT events that both normalize to this site.
                                    if ambiguous.contains(pos) {
                                        stats.assignments_skipped_ambiguous_carry += 1;
                                        continue;
                                    }
                                    let (observed, anchor) = carried
                                        .get(pos)
                                        .cloned()
                                        .unwrap_or((indel::Allele::Ref, *pos));
                                    // an earlier gained deletion already emitted for this read consumed this site's reference position.
                                    if observed == indel::Allele::Ref && s <= consumed_until {
                                        stats.assignments_skipped_swallowed_by_gained_deletion += 1;
                                        continue;
                                    }
                                    let target = match indel::assign_allele(&observed, model, iopts) {
                                        indel::Assignment::Keep => continue,
                                        indel::Assignment::Refused => {
                                            stats.assignments_skipped_unsupported_transition += 1;
                                            continue;
                                        }
                                        indel::Assignment::Move(t) => t,
                                    };
                                    // Left-normalization only ever moves an event LEFT.
                                    let far = (anchor as i64).max(s);
                                    // Require the whole tract plus the usual flank of match context beyond it.
                                    let tract_need =
                                        indel::hp_tract_end(refseq, s as u32) as i64 + flank;
                                    if s - flank < rstart || (far + reach).max(tract_need) > rend {
                                        stats.assignments_skipped_span_guard += 1;
                                        continue;
                                    }
                                    // the check above bounds only the read's.
                                    if observed == indel::Allele::Ref {
                                        let run = indel::match_run_len(rstart, &cigar, s as u32) as i64;
                                        let need = match &target {
                                            indel::Allele::Del(m) => *m as i64 + flank + 1,
                                            _ => 1,
                                        };
                                        if run < need {
                                            stats.assignments_skipped_span_guard += 1;
                                            continue;
                                        }
                                    }
                                    // An oversized indel nearby means this read's alignment in this window cannot be represented by the site's alleles.
                                    if blocked.iter().any(|b| (b - s).abs() <= reach) {
                                        stats.assignments_skipped_oversized_neighbor += 1;
                                        continue;
                                    }
                                    // A gained deletion consumes reference positions (s, s+m]; later sites inside that span are skipped.
                                    if let (indel::Allele::Ref, indel::Allele::Del(m)) = (&observed, &target) {
                                        consumed_until = consumed_until.max(s + *m as i64);
                                    }

                                    // Direction determines which coordinate `rewrite_read` will actually look this edit up by.
                                    let ref_pos = match &observed {
                                        indel::Allele::Ref => *pos,
                                        _ => anchor,
                                    };
                                    // Reassignment/base-count/hp-length statistics are NOT booked here.
                                    edit_context_len.insert(
                                        (ref_pos, indel::edit_kind(&observed)),
                                        model.context_len as usize,
                                    );
                                    edits.indels.push(indel::IndelEdit {
                                        ref_pos,
                                        from: observed,
                                        to: target,
                                    });
                                }
                            }
                            edits.indels.sort_by_key(|e| e.ref_pos);
                        }

                        if !edits.is_empty() || !hp_moves.is_empty() {
                            // Extract owned parts, rewrite, rebuild (aux preserved by set()).
                            let qname = rec.qname().to_vec();
                            let cigar = rec.cigar().take();
                            let seq = rec.seq().as_bytes();
                            let qual = rec.qual().to_vec();

                            // Each `rewrite_hp_run` preserves reference consumption.
                            let mut hp_cigar: Option<CigarString> = None;
                            let mut hp_seq: Option<Vec<u8>> = None;
                            let mut hp_qual: Option<Vec<u8>> = None;
                            // Half-open reference spans of the runs actually rewritten.
                            let mut hp_runs: Vec<(u32, u32)> = Vec::new();
                            for (a, e, b, from, to) in &hp_moves {
                                let rewritten = {
                                    let c = hp_cigar.as_ref().unwrap_or(&cigar);
                                    let s = hp_seq.as_ref().unwrap_or(&seq);
                                    let q = hp_qual.as_ref().unwrap_or(&qual);
                                    indel::rewrite_hp_run(
                                        rec.pos(), c, s, q, *a, *e, *b, to.net_len(),
                                    )
                                };
                                match rewritten {
                                    Some((c, s, q)) => {
                                        hp_cigar = Some(c);
                                        hp_seq = Some(s);
                                        hp_qual = Some(q);
                                        hp_runs.push((*a + 1, *e));
                                        book_reassignment(stats, from, to);
                                        let l = (*e - *a - 1) as usize;
                                        if stats.reassignments_by_hp_length.len() <= l {
                                            stats.reassignments_by_hp_length.resize(l + 1, 0);
                                        }
                                        stats.reassignments_by_hp_length[l] += 1;
                                    }
                                    // Same accounting as an event edit the walk never reached: a harmless missed correction.
                                    None => stats.indel_edits_emitted_but_not_applied += 1,
                                }
                            }
                            if !hp_runs.is_empty() {
                                let before = edits.subs.len();
                                edits.subs.retain(|&(p, _)| {
                                    !hp_runs.iter().any(|&(lo, hi)| p >= lo && p < hi)
                                });
                                stats.subs_dropped_in_hp_rewrite +=
                                    (before - edits.subs.len()) as u64;
                            }

                            let in_cigar = hp_cigar.as_ref().unwrap_or(&cigar);
                            let in_seq = hp_seq.as_ref().unwrap_or(&seq);
                            let in_qual = hp_qual.as_ref().unwrap_or(&qual);
                            let out =
                                rewrite_read(rec.pos(), in_cigar, in_seq, in_qual, &edits, refseq);
                            // A non-empty edit list does not guarantee `rewrite_read` changed anything.
                            let changed = out.seq != seq || out.cigar != cigar;

                            // the reassignment, base-count, and hp-length counters are booked HERE, and only here.
                            let mut applied = out.applied;
                            for edit in &edits.indels {
                                let key = (edit.ref_pos, indel::edit_kind(&edit.from));
                                if !applied.remove(&key) {
                                    stats.indel_edits_emitted_but_not_applied += 1;
                                    continue;
                                }
                                book_reassignment(stats, &edit.from, &edit.to);
                                let l = *edit_context_len.get(&key).unwrap_or(&0);
                                if stats.reassignments_by_hp_length.len() <= l {
                                    stats.reassignments_by_hp_length.resize(l + 1, 0);
                                }
                                stats.reassignments_by_hp_length[l] += 1;
                            }

                            rec.set(&qname, Some(&out.cigar), &out.seq, &out.qual);
                            // An HP-length move shifts every query offset downstream of the run just as an event edit does.
                            if out.structure_changed || !hp_runs.is_empty() {
                                strip_stale_tags(&mut rec);
                            }
                            if changed {
                                reads_modified += 1;
                            }
                        }
                    }
                    None => {
                        // Contig not in the supplied reference: do not rewrite this record.
                        if corrections.get(rec.qname()).is_some_and(|e| !e.is_empty()) {
                            missing_ref_contigs.insert(contig.clone());
                            missing_ref_reads += 1;
                        }
                    }
                }
            }
        }
        writer.write(&rec).context("write BAM record")?;
    }
    if missing_ref_reads > 0 {
        let contigs = missing_ref_contigs
            .iter()
            .map(|c| String::from_utf8_lossy(c).into_owned())
            .collect::<Vec<_>>()
            .join(", ");
        warn!(
            "[denoise] {missing_ref_reads} correctable read(s) on contig(s) [{contigs}] \
             not found in the supplied reference; written through unchanged."
        );
    }
    stats.reads_processed = reads_processed;
    stats.reads_modified = reads_modified;
    Ok(())
}

pub fn start(
    input: &PathBuf,
    output: &PathBuf,
    reference: &PathBuf,
    data_type: &str,
    vaf: f64,
    homoplasmic_vaf: f64,
    indels: &indel::IndelOpts,
    stats_out: Option<&PathBuf>,
) -> Result<()> {
    // The per-strand count gate and the strand-bias p are fixed crate constants.
    let (min_strand, strand_bias_p) = (crate::DENOISE_MIN_STRAND, crate::DENOISE_STRAND_BIAS_P);
    // PacBio (and anything non-ONT) is a byte-identical passthrough.
    if !data_type.starts_with("ont") {
        std::fs::copy(input, output)
            .with_context(|| format!("passthrough copy {input:?} -> {output:?}"))?;
        info!("[denoise] data-type '{data_type}' is not ONT; passthrough copy.");
        return Ok(());
    }

    // Load reference contigs into memory (uppercased).
    let mut refs: HashMap<Vec<u8>, Vec<u8>> = HashMap::new();
    let rdr = fasta::Reader::from_file(reference)
        .with_context(|| format!("cannot open reference {reference:?}"))?;
    for rec in rdr.records() {
        let rec = rec.context("read reference record")?;
        refs.insert(rec.id().as_bytes().to_vec(), rec.seq().to_ascii_uppercase());
    }

    let (corrections, models, mut stats) = compute_corrections(
        input, &refs, min_strand, vaf, strand_bias_p, homoplasmic_vaf, indels,
    )?;
    apply_corrections(input, output, &corrections, &refs, &models, indels, &mut stats)?;

    let rate = if stats.bases_examined > 0 {
        stats.bases_modified as f64 / stats.bases_examined as f64
    } else {
        0.0
    };
    info!(
        "[denoise] reads {}/{} modified; bases {}/{} modified ({:.2}%); \
         columns {} processed / {} corrected / {} skipped; alt sites preserved {}; \
         alt alleles rejected as strand-biased {}",
        stats.reads_modified, stats.reads_processed,
        stats.bases_modified, stats.bases_examined, rate * 100.0,
        stats.columns_processed, stats.columns_corrected, stats.columns_skipped,
        stats.alt_sites_preserved, stats.alt_alleles_strand_biased,
    );

    if indels.enabled {
        info!(
            "[denoise] indel sites {} examined / {} alt alleles kept / {} strand-biased; \
             reads reassigned ref->ins {} ref->del {} indel->ref {} indel->indel {}; \
             skipped span-guard {} unsupported-transition {}",
            stats.indel_sites_examined, stats.indel_alt_alleles_kept, stats.indel_alt_alleles_strand_biased,
            stats.reads_reassigned_ref_to_ins, stats.reads_reassigned_ref_to_del,
            stats.reads_reassigned_indel_to_ref, stats.reads_reassigned_indel_to_indel,
            stats.assignments_skipped_span_guard, stats.assignments_skipped_unsupported_transition,
        );
    }

    if let Some(path) = stats_out {
        let json = serde_json::to_string_pretty(&stats).context("serialize stats")?;
        std::fs::write(path, json).with_context(|| format!("write stats {path:?}"))?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::{self, Format, Header, Read, Reader, Record, Writer};
    use rust_htslib::bam::header::HeaderRecord;
    use rust_htslib::bam::record::{Cigar, CigarString};

    fn q(v: u8) -> u8 { v } // readability helper for quality values

    // Defaults mirrored from the CLI so tests exercise the shipped configuration.
    const SB_P: f64 = 0.01;
    const HOM_VAF: f64 = 0.7;

    #[test]
    fn phred_to_e_matches_definition() {
        assert!((phred_to_e(10) - 0.1).abs() < 1e-12);
        assert!((phred_to_e(20) - 0.01).abs() < 1e-12);
        assert!((phred_to_e(2) - 0.6309573444801932).abs() < 1e-12); // below the clamp
    }

    #[test]
    fn phred_to_e_clamps_degenerate_qualities() {
        // Q0 nominally means P(error) = 1, which would make P(obs | truth = obs) zero and invert every base.
        assert!((phred_to_e(0) - E_MAX).abs() < 1e-12);
        assert!((phred_to_e(1) - E_MAX).abs() < 1e-12);
        for q in 0..=93u8 {
            assert!(phred_to_e(q) <= E_MAX, "q={q} exceeds the clamp");
        }
    }

    #[test]
    fn zero_quality_observations_are_never_inverted() {
        // ref A(0) at 60%, alt G(2) at 40%, both strand-balanced, ALL bases at Q0 -- the shape pbsim's errhmm models emit.
        let mut obs = vec![];
        for i in 0..60 { obs.push((0usize, q(0), i % 2 == 0)); }
        for i in 0..40 { obs.push((2usize, q(0), i % 2 == 0)); }
        let c = column(0, &obs);
        let m = fit_site(&c, 2, 0.01, SB_P, HOM_VAF);

        // Before the clamp, like(k == observed) was 0 while like(k != observed) was e/3 > 0.
        assert_eq!(map_allele(0, q(0), &m), 0);
        // At the clamp the likelihood is flat (1 - e == e/3 == 0.25), so a Q0 base carries no evidence of its own.
        assert_eq!(map_allele(2, q(0), &m), 2);
        // The alleles the model rejected have no evidence and no prior, so they fall back to the most frequent kept allele.
        for observed in 0..4 {
            let got = map_allele(observed, q(0), &m);
            let want = if m.kept[observed] { observed } else { 0 };
            assert_eq!(got, want, "observed={observed}");
        }
    }

    #[test]
    fn low_quality_observations_of_a_kept_alt_are_not_reverted_to_ref() {
        // A ~8% heteroplasmy in ONT-quality reads: the alt clears the keep rule, so the site model has already ruled it real.
        let quals = [3u8, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 20];
        let mut obs = vec![];
        for i in 0..736 { obs.push((0usize, q(quals[i % quals.len()]), i % 2 == 0)); }
        for i in 0..64 { obs.push((2usize, q(quals[i % quals.len()]), i % 2 == 0)); }
        let c = column(0, &obs);
        let m = fit_site(&c, 2, 0.01, SB_P, HOM_VAF);
        assert!(m.kept[0] && m.kept[2], "both ref and the 8% alt must be kept");

        // The alt is kept, so an observation of it stands regardless of quality.
        for ql in [q(0), q(4), q(6), q(7), q(30)] {
            assert_eq!(map_allele(2, ql, &m), 2, "kept alt reverted at Q{ql}");
            assert_eq!(map_allele(0, ql, &m), 0, "ref reverted at Q{ql}");
        }
        // Alleles the model rejected are still errors and still correct to ref.
        assert_eq!(map_allele(1, q(6), &m), 0);
        assert_eq!(map_allele(3, q(6), &m), 0);
    }

    #[test]
    fn base_index_round_trips_and_rejects_n() {
        for (b, i) in [(b'A', 0usize), (b'C', 1), (b'G', 2), (b'T', 3)] {
            assert_eq!(base_to_idx(b), Some(i));
            assert_eq!(idx_to_base(i), b);
        }
        assert_eq!(base_to_idx(b'N'), None);
        assert_eq!(base_to_idx(b'a'), Some(0)); // lowercase accepted
    }

    // Build a column: `obs` is (allele_idx, qual, reverse) repeated.
    fn column(ref_idx: usize, obs: &[(usize, u8, bool)]) -> ColumnCounts {
        let mut c = ColumnCounts::new(Some(ref_idx));
        for &(a, ql, rev) in obs { c.add(a, ql, rev); }
        c
    }

    #[test]
    fn pure_ref_column_with_sparse_errors_keeps_only_ref() {
        // ref = A(0); 30 clean A plus 2 stray G/C at low quality, not strand-balanced.
        let mut obs = vec![];
        for i in 0..30 { obs.push((0usize, q(30), i % 2 == 0)); }
        obs.push((2, q(8), false)); // a single G error (fwd only)
        obs.push((1, q(8), true));  // a single C error (rev only)
        let c = column(0, &obs);
        let m = fit_site(&c, 2, 0.01, SB_P, HOM_VAF);
        assert!(m.kept[0]); // ref kept
        assert_eq!(m.kept_alt_count(Some(0)), 0); // no alt survives
        // the two errors correct to ref A
        assert_eq!(map_allele(2, q(8), &m), 0);
        assert_eq!(map_allele(1, q(8), &m), 0);
    }

    #[test]
    fn strand_balanced_heteroplasmy_is_preserved() {
        // ref A(0) ~60%, alt G(2) ~40%, both strand-balanced, decent quality.
        let mut obs = vec![];
        for i in 0..60 { obs.push((0usize, q(30), i % 2 == 0)); }
        for i in 0..40 { obs.push((2usize, q(30), i % 2 == 0)); }
        let c = column(0, &obs);
        let m = fit_site(&c, 2, 0.01, SB_P, HOM_VAF);
        assert!(m.kept[0] && m.kept[2]);
        assert_eq!(m.kept_alt_count(Some(0)), 1);
        // confident observations of each allele are preserved
        assert_eq!(map_allele(0, q(30), &m), 0);
        assert_eq!(map_allele(2, q(30), &m), 2);
        // a third-base error (T) collapses to the higher-frequency allele (ref A)
        assert_eq!(map_allele(3, q(8), &m), 0);
    }

    #[test]
    fn strand_skewed_candidate_is_corrected_away() {
        // ref A(0); an alt G(2) at 10% but ALL on the forward strand -> fails candidacy.
        let mut obs = vec![];
        for i in 0..90 { obs.push((0usize, q(30), i % 2 == 0)); }
        for _ in 0..10 { obs.push((2usize, q(30), false)); } // fwd only
        let c = column(0, &obs);
        let m = fit_site(&c, 2, 0.01, SB_P, HOM_VAF);
        assert!(!m.kept[2]);
        assert_eq!(map_allele(2, q(30), &m), 0); // corrected to ref
    }

    #[test]
    fn low_vaf_but_balanced_allele_below_threshold_is_dropped() {
        // alt G(2) present twice per strand (passes candidacy) but VAF ~0.02 < 0.05 threshold.
        let mut obs = vec![];
        for i in 0..196 { obs.push((0usize, q(30), i % 2 == 0)); }
        for _ in 0..2 { obs.push((2usize, q(30), false)); }
        for _ in 0..2 { obs.push((2usize, q(30), true)); }
        let c = column(0, &obs);
        let m = fit_site(&c, 2, 0.05, SB_P, HOM_VAF); // higher vaf threshold
        assert!(!m.kept[2]);
    }

    #[test]
    fn strand_biased_candidate_passing_min_strand_is_rejected() {
        // ref A(0) 40fwd/40rev; alt G(2) 18fwd/2rev. min_strand=2 is satisfied on BOTH strands.
        let mut obs = vec![];
        for i in 0..80 { obs.push((0usize, q(30), i % 2 == 0)); }
        for _ in 0..18 { obs.push((2usize, q(30), false)); }
        for _ in 0..2 { obs.push((2usize, q(30), true)); }
        let c = column(0, &obs);

        // With the test disabled (p = 0), the count gate keeps it: VAF 0.20 >> 0.01.
        let m_off = fit_site(&c, 2, 0.01, 0.0, HOM_VAF);
        assert!(m_off.kept[2], "min_strand alone should not reject this allele");
        assert_eq!(m_off.strand_biased_alts, 0);

        // With the test on, the allele loses candidacy and its reads correct to ref.
        let m_on = fit_site(&c, 2, 0.01, SB_P, HOM_VAF);
        assert!(!m_on.kept[2]);
        assert_eq!(m_on.kept_alt_count(Some(0)), 0);
        assert_eq!(m_on.strand_biased_alts, 1);
        assert_eq!(map_allele(2, q(30), &m_on), 0);
    }

    #[test]
    fn near_homoplasmic_allele_is_exempt_from_strand_bias_test() {
        // ref A(0) 0fwd/25rev; alt G(2) 73fwd/2rev -> alt is 75% of the column.
        let mut obs = vec![];
        for _ in 0..25 { obs.push((0usize, q(30), true)); }
        for _ in 0..73 { obs.push((2usize, q(30), false)); }
        for _ in 0..2 { obs.push((2usize, q(30), true)); }
        let c = column(0, &obs);

        let m = fit_site(&c, 2, 0.01, SB_P, 0.7);
        assert!(m.kept[2], "near-homoplasmic alt must survive the strand-bias test");
        assert_eq!(m.strand_biased_alts, 0);
        assert_eq!(map_allele(2, q(30), &m), 2);

        // Raising the exemption bar above the allele's frequency re-enables the test.
        let m_strict = fit_site(&c, 2, 0.01, SB_P, 0.95);
        assert!(!m_strict.kept[2]);
        assert_eq!(m_strict.strand_biased_alts, 1);
    }

    fn strand_column(ref_f: u32, ref_r: u32, alt_f: u32, alt_r: u32) -> ColumnCounts {
        let mut obs = vec![];
        for _ in 0..ref_f { obs.push((0usize, q(30), false)); }
        for _ in 0..ref_r { obs.push((0usize, q(30), true)); }
        for _ in 0..alt_f { obs.push((2usize, q(30), false)); }
        for _ in 0..alt_r { obs.push((2usize, q(30), true)); }
        column(0, &obs)
    }

    #[test]
    fn deep_allele_with_small_strand_skew_stays_eligible() {
        // Alt 328/196 at 10% in a column whose forward fraction is 2693/5240 = 0.514.
        let c = strand_column(2365, 2351, 328, 196);
        assert!(
            crate::call::strand_bias_pvalue(328, 196, 2693.0 / 5240.0) < SB_P,
            "test premise: the p-value alone would reject this allele"
        );
        let m = fit_site(&c, 2, 0.01, SB_P, HOM_VAF);
        assert!(m.kept[2], "a 0.11 strand effect must not erase a real allele");
        assert_eq!(m.strand_biased_alts, 0);
        assert_eq!(map_allele(2, q(30), &m), 2);
    }

    #[test]
    fn deep_single_strand_artifact_is_still_rejected() {
        // Alt 449/15 in a column whose forward fraction is 2171/4640 = 0.468, so the effect is 0.968 - 0.468 = 0.50.
        let c = strand_column(1722, 2454, 449, 15);
        let m = fit_site(&c, 2, 0.01, SB_P, HOM_VAF);
        assert!(!m.kept[2], "a 0.50 strand effect is still an artifact");
        assert_eq!(m.strand_biased_alts, 1);
        assert_eq!(map_allele(2, q(30), &m), 0);
    }

    #[test]
    fn empty_column_falls_back_to_reference_only() {
        // No observations at all: expected-fwd-fraction fallback must not panic and the model must degenerate to reference-only.
        let c = column(0, &[]);
        let m = fit_site(&c, 2, 0.01, SB_P, HOM_VAF);
        assert!(m.kept[0]);
        assert_eq!(m.kept_alt_count(Some(0)), 0);
        assert!((m.freq[0] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn eligible_allele_between_floor_and_keep_vaf_is_pending() {
        let c = strand_column(486, 486, 14, 14); // 2.8%
        let m = fit_site(&c, 2, 0.03, SB_P, HOM_VAF);
        assert!(!m.kept[2]);
        assert!(m.pending[2]);
        assert!(!m.pending[0], "reference is never pending");

        let kept = fit_site(&strand_column(482, 482, 18, 18), 2, 0.03, SB_P, HOM_VAF);
        assert!(kept.kept[2] && !kept.pending[2]);

        let below = fit_site(&strand_column(998, 998, 2, 2), 2, 0.03, SB_P, HOM_VAF); // 0.2%
        assert!(!below.kept[2] && !below.pending[2]);

        let biased = fit_site(&strand_column(486, 486, 26, 2), 2, 0.03, SB_P, HOM_VAF);
        assert_eq!(biased.strand_biased_alts, 1);
        assert!(!biased.pending[2]);
    }

    // Reads in `m_reads` carry C at `m_pos`, reads in `p_reads` carry G at `p_pos`; `extra` adds substitutions.
    fn linkage_bam(
        path: &Path,
        m_pos: usize,
        m_reads: std::ops::Range<usize>,
        p_pos: usize,
        p_reads: std::ops::Range<usize>,
        extra: &[(usize, usize, u8)],
    ) {
        let names: Vec<String> = (0..400).map(|i| format!("r{i:03}")).collect();
        let seqs: Vec<Vec<u8>> = (0..400)
            .map(|i| {
                let mut s = vec![b'A'; 100];
                if m_reads.contains(&i) { s[m_pos] = b'C'; }
                if p_reads.contains(&i) { s[p_pos] = b'G'; }
                for &(r, pos, b) in extra {
                    if r == i { s[pos] = b; }
                }
                s
            })
            .collect();
        let reads: Vec<(&str, i64, &[u8], bool)> = (0..400)
            .map(|i| (names[i].as_str(), 0i64, seqs[i].as_slice(), i % 2 == 1))
            .collect();
        write_test_bam(path, "chrM", 100, &reads);
    }

    fn linkage_refs() -> HashMap<Vec<u8>, Vec<u8>> {
        HashMap::from([(b"chrM".to_vec(), vec![b'A'; 100])])
    }

    fn subs_of(corr: &Corrections, name: &str) -> Vec<(u32, u8)> {
        corr.get(name.as_bytes()).map(|e| e.subs.clone()).unwrap_or_default()
    }

    #[test]
    fn pending_allele_riding_on_a_kept_minor_is_preserved() {
        let dir = std::env::temp_dir().join("himito_denoise_t6b_linked");
        std::fs::create_dir_all(&dir).unwrap();
        let bam = dir.join("in.bam");
        linkage_bam(&bam, 10, 0..14, 60, 0..11, &[]);
        let (corr, _, stats) = compute_corrections(
            &bam, &linkage_refs(), 2, 0.03, SB_P, HOM_VAF, &indel::IndelOpts::default(),
        )
        .unwrap();
        for i in 0..14 {
            assert!(subs_of(&corr, &format!("r{i:03}")).is_empty(), "r{i:03} was corrected");
        }
        assert_eq!(stats.alt_alleles_rescued_by_linkage, 1);
        assert_eq!(stats.bases_modified, 0);
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn unlinked_pending_allele_is_still_corrected_and_subs_stay_sorted() {
        let dir = std::env::temp_dir().join("himito_denoise_t6b_unlinked");
        std::fs::create_dir_all(&dir).unwrap();
        let bam = dir.join("in.bam");
        linkage_bam(&bam, 10, 0..14, 60, 100..111, &[(100, 80, b'T')]);
        let (corr, _, stats) = compute_corrections(
            &bam, &linkage_refs(), 2, 0.03, SB_P, HOM_VAF, &indel::IndelOpts::default(),
        )
        .unwrap();
        assert_eq!(subs_of(&corr, "r100"), vec![(60, b'A'), (80, b'A')]);
        for i in 101..111 {
            assert_eq!(subs_of(&corr, &format!("r{i:03}")), vec![(60, b'A')]);
        }
        for i in 0..14 {
            assert!(subs_of(&corr, &format!("r{i:03}")).is_empty());
        }
        assert_eq!(stats.alt_alleles_rescued_by_linkage, 0);
        assert_eq!(stats.bases_modified, 12);
        assert_eq!(stats.substitution_matrix[2][0], 11);
        assert_eq!(stats.columns_corrected, 2);
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn linked_partner_within_20bp_does_not_rescue() {
        let dir = std::env::temp_dir().join("himito_denoise_t6b_near");
        std::fs::create_dir_all(&dir).unwrap();
        let bam = dir.join("in.bam");
        linkage_bam(&bam, 10, 0..14, 25, 0..11, &[]);
        let (corr, _, stats) = compute_corrections(
            &bam, &linkage_refs(), 2, 0.03, SB_P, HOM_VAF, &indel::IndelOpts::default(),
        )
        .unwrap();
        for i in 0..11 {
            assert_eq!(subs_of(&corr, &format!("r{i:03}")), vec![(25, b'A')]);
        }
        assert_eq!(stats.alt_alleles_rescued_by_linkage, 0);
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn linkage_rescue_is_deterministic() {
        let dir = std::env::temp_dir().join("himito_denoise_t6b_det");
        std::fs::create_dir_all(&dir).unwrap();
        let bam = dir.join("in.bam");
        let extra: Vec<(usize, usize, u8)> = (200..211).map(|r| (r, 85, b'T')).collect();
        linkage_bam(&bam, 10, 0..14, 60, 0..11, &extra);
        let run = || {
            let (corr, _, stats) = compute_corrections(
                &bam, &linkage_refs(), 2, 0.03, SB_P, HOM_VAF, &indel::IndelOpts::default(),
            )
            .unwrap();
            let mut v: Vec<(Vec<u8>, Vec<(u32, u8)>)> =
                corr.into_iter().map(|(k, e)| (k, e.subs)).collect();
            v.sort();
            (v, serde_json::to_string(&stats).unwrap())
        };
        let (a, sa) = run();
        let (b, sb) = run();
        assert_eq!(a, b);
        assert_eq!(sa, sb);
        assert_eq!(a.len(), 11, "only the unlinked T@85 carriers are corrected");
        assert!(a.iter().all(|(_, s)| s == &vec![(85u32, b'A')]));
        std::fs::remove_dir_all(&dir).ok();
    }

    // Build a tiny coordinate-sorted BAM.
    fn write_test_bam(path: &Path, contig: &str, contig_len: usize, reads: &[(&str, i64, &[u8], bool)]) {
        let mut header = Header::new();
        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", contig);
        sq.push_tag(b"LN", &contig_len.to_string());
        header.push_record(&sq);
        let mut w = Writer::from_path(path, &header, Format::Bam).unwrap();
        for (name, pos, seq, reverse) in reads {
            let mut rec = Record::new();
            let quals = vec![30u8; seq.len()];
            let cigar = CigarString(vec![Cigar::Match(seq.len() as u32)]);
            rec.set(name.as_bytes(), Some(&cigar), seq, &quals);
            rec.set_tid(0);
            rec.set_pos(*pos);
            rec.set_mapq(60);
            rec.set_flags(if *reverse { 16 } else { 0 });
            w.write(&rec).unwrap();
        }
    }

    #[test]
    fn compute_corrections_flags_a_single_base_error() {
        let dir = std::env::temp_dir().join("himito_denoise_t2");
        std::fs::create_dir_all(&dir).unwrap();
        let bam = dir.join("in.bam");
        // Reference AAAAAA at chrM.
        let refseq = b"AAAAAA".to_vec();
        let mut reads: Vec<(&str, i64, &[u8], bool)> = vec![];
        let clean = b"AAAAAA";
        let errd = b"AAGAAA"; // error at index 2
        for i in 0..9 { reads.push((["r0","r1","r2","r3","r4","r5","r6","r7","r8"][i], 0, clean, i % 2 == 0)); }
        reads.push(("rerr", 0, errd, false));
        write_test_bam(&bam, "chrM", 6, &reads);

        let mut refs = HashMap::new();
        refs.insert(b"chrM".to_vec(), refseq);
        let (corr, _models, stats) =
            compute_corrections(&bam, &refs, 2, 0.01, SB_P, HOM_VAF, &indel::IndelOpts::default())
                .unwrap();

        // The single G (not strand-balanced, VAF 10% but fwd-only) is corrected to A.
        let edits = corr.get(&b"rerr".to_vec()).expect("rerr must be corrected");
        assert_eq!(edits.subs, vec![(2u32, b'A')]);
        assert_eq!(stats.bases_modified, 1);
        assert_eq!(stats.substitution_matrix[2][0], 1); // G(2) -> A(0)
        assert!(stats.columns_processed >= 1);
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn apply_corrections_rewrites_seq_and_preserves_tags() {
        use rust_htslib::bam::record::Aux;
        let dir = std::env::temp_dir().join("himito_denoise_t3");
        std::fs::create_dir_all(&dir).unwrap();
        let inb = dir.join("in.bam");
        let outb = dir.join("out.bam");

        // One read "AAGAAA" with an aux tag we must preserve.
        {
            let mut header = Header::new();
            let mut sq = HeaderRecord::new(b"SQ");
            sq.push_tag(b"SN", "chrM");
            sq.push_tag(b"LN", &6usize.to_string());
            header.push_record(&sq);
            let mut w = Writer::from_path(&inb, &header, Format::Bam).unwrap();
            let mut rec = Record::new();
            let seq = b"AAGAAA";
            let quals = vec![30u8; 6];
            rec.set(b"rerr", Some(&CigarString(vec![Cigar::Match(6)])), seq, &quals);
            rec.set_tid(0);
            rec.set_pos(0);
            // Record::new() initializes as unmapped (see rust-htslib set_unmapped() in Record::new()).
            rec.unset_unmapped();
            rec.push_aux(b"XY", Aux::I32(7)).unwrap();
            w.write(&rec).unwrap();
        }

        let mut corr: Corrections = HashMap::new();
        corr.insert(
            b"rerr".to_vec(),
            indel::ReadEdits { subs: vec![(2u32, b'A')], indels: vec![] },
        );
        let mut refs = HashMap::new();
        refs.insert(b"chrM".to_vec(), b"AAAAAA".to_vec());
        let models = IndelModels::new();
        let mut stats = DenoiseStats::default();
        apply_corrections(
            &inb, &outb, &corr, &refs, &models, &indel::IndelOpts::default(), &mut stats,
        )
        .unwrap();
        let (proc, modified) = (stats.reads_processed, stats.reads_modified);
        assert_eq!(proc, 1);
        assert_eq!(modified, 1);

        // Verify SEQ was corrected and the aux tag survived.
        let mut r = Reader::from_path(&outb).unwrap();
        let rec = r.records().next().unwrap().unwrap();
        assert_eq!(rec.seq().as_bytes(), b"AAAAAA");
        assert!(matches!(rec.aux(b"XY").unwrap(), Aux::I32(7)));
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn start_pacbio_is_byte_identical_passthrough() {
        let dir = std::env::temp_dir().join("himito_denoise_t4pb");
        std::fs::create_dir_all(&dir).unwrap();
        let inb = dir.join("in.bam");
        let outb = dir.join("out.bam");
        let reff = dir.join("ref.fa");
        std::fs::write(&reff, ">chrM\nAAAAAA\n").unwrap();
        write_test_bam(&inb, "chrM", 6, &[("r0", 0, b"AAGAAA", false)]);

        start(&inb, &outb, &reff, "pacbio", 0.01, HOM_VAF, &indel::IndelOpts::default(), None).unwrap();
        assert_eq!(std::fs::read(&inb).unwrap(), std::fs::read(&outb).unwrap());
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn start_ont_denoises_and_writes_stats_json() {
        let dir = std::env::temp_dir().join("himito_denoise_t4ont");
        std::fs::create_dir_all(&dir).unwrap();
        let inb = dir.join("in.bam");
        let outb = dir.join("out.bam");
        let reff = dir.join("ref.fa");
        let statsf = dir.join("stats.json");
        std::fs::write(&reff, ">chrM\nAAAAAA\n").unwrap();
        let clean = b"AAAAAA";
        let mut reads: Vec<(&str, i64, &[u8], bool)> = vec![];
        for i in 0..9 { reads.push((["r0","r1","r2","r3","r4","r5","r6","r7","r8"][i], 0, clean, i % 2 == 0)); }
        reads.push(("rerr", 0, b"AAGAAA", false));
        write_test_bam(&inb, "chrM", 6, &reads);

        start(&inb, &outb, &reff, "ont-r10", 0.01, HOM_VAF, &indel::IndelOpts::default(), Some(&statsf)).unwrap();

        let mut r = Reader::from_path(&outb).unwrap();
        let mut modified_ok = false;
        for rec in r.records() {
            let rec = rec.unwrap();
            if rec.qname() == b"rerr" {
                assert_eq!(rec.seq().as_bytes(), b"AAAAAA");
                modified_ok = true;
            }
        }
        assert!(modified_ok);
        let json = std::fs::read_to_string(&statsf).unwrap();
        // serde_json::to_string_pretty inserts a space after the colon, unlike compact serialization.
        assert!(json.contains("\"bases_modified\": 1"));
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn corrections_are_keyed_by_reference_position() {
        // The single G error sits at reference position 2 (reads start at POS 0).
        let dir = std::env::temp_dir().join("himito_denoise_refcoord");
        std::fs::create_dir_all(&dir).unwrap();
        let bam = dir.join("in.bam");
        let mut reads: Vec<(&str, i64, &[u8], bool)> = vec![];
        let clean = b"AAAAAA";
        for i in 0..9 {
            reads.push((["r0","r1","r2","r3","r4","r5","r6","r7","r8"][i], 0, clean, i % 2 == 0));
        }
        reads.push(("rerr", 0, b"AAGAAA", false));
        write_test_bam(&bam, "chrM", 6, &reads);

        let mut refs = HashMap::new();
        refs.insert(b"chrM".to_vec(), b"AAAAAA".to_vec());
        let (corr, _models, _) =
            compute_corrections(&bam, &refs, 2, 0.01, SB_P, HOM_VAF, &indel::IndelOpts::default())
                .unwrap();

        let edits = corr.get(&b"rerr".to_vec()).expect("rerr must be corrected");
        assert_eq!(edits.subs, vec![(2u32, b'A')]);
        assert!(edits.indels.is_empty());
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn corrections_use_reference_position_not_query_offset() {
        // Reads aligned at POS 10 (not 0): the erroneous base sits at read/query offset 2, which maps to REFERENCE position 12.
        let dir = std::env::temp_dir().join("himito_denoise_refcoord_nonzero_pos");
        std::fs::create_dir_all(&dir).unwrap();
        let bam = dir.join("in.bam");
        // Reference long enough to cover POS 10 + a 6bp read.
        let refseq = vec![b'A'; 20];
        let clean = b"AAAAAA";
        let errd = b"AAGAAA"; // error at read offset 2 -> reference position 10+2=12
        let mut reads: Vec<(&str, i64, &[u8], bool)> = vec![];
        for i in 0..9 {
            reads.push((["r0","r1","r2","r3","r4","r5","r6","r7","r8"][i], 10, clean, i % 2 == 0));
        }
        reads.push(("rerr", 10, errd, false));
        write_test_bam(&bam, "chrM", 20, &reads);

        let mut refs = HashMap::new();
        refs.insert(b"chrM".to_vec(), refseq);
        let (corr, _models, _) =
            compute_corrections(&bam, &refs, 2, 0.01, SB_P, HOM_VAF, &indel::IndelOpts::default())
                .unwrap();

        let edits = corr.get(&b"rerr".to_vec()).expect("rerr must be corrected");
        // The recorded key must be the REFERENCE position.
        assert_eq!(edits.subs, vec![(12u32, b'A')]);
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn apply_corrections_missing_contig_writes_through_unchanged() {
        // A record on a contig absent from `refs`, with a correction registered for its qname, must NOT be rewritten.
        let dir = std::env::temp_dir().join("himito_denoise_missing_contig");
        std::fs::create_dir_all(&dir).unwrap();
        let inb = dir.join("in.bam");
        let outb = dir.join("out.bam");

        write_test_bam(&inb, "chrM", 6, &[("rerr", 0, b"AAGAAA", false)]);

        let mut corr: Corrections = HashMap::new();
        corr.insert(
            b"rerr".to_vec(),
            indel::ReadEdits { subs: vec![(2u32, b'A')], indels: vec![] },
        );
        // Deliberately do NOT insert "chrM" -- simulates a BAM whose contig is missing from the supplied reference FASTA.
        let refs: HashMap<Vec<u8>, Vec<u8>> = HashMap::new();

        let models = IndelModels::new();
        let mut stats = DenoiseStats::default();
        apply_corrections(
            &inb, &outb, &corr, &refs, &models, &indel::IndelOpts::default(), &mut stats,
        )
        .unwrap();
        let (proc, modified) = (stats.reads_processed, stats.reads_modified);
        assert_eq!(proc, 1);
        assert_eq!(modified, 0, "a record on a contig missing from refs must not be rewritten");

        let mut orig_reader = Reader::from_path(&inb).unwrap();
        let orig_rec = orig_reader.records().next().unwrap().unwrap();
        let mut out_reader = Reader::from_path(&outb).unwrap();
        let out_rec = out_reader.records().next().unwrap().unwrap();
        assert_eq!(out_rec.seq().as_bytes(), orig_rec.seq().as_bytes(), "SEQ must be byte-for-byte unchanged");
        assert_eq!(out_rec.cigar().take(), orig_rec.cigar().take(), "CIGAR must be byte-for-byte unchanged");
        std::fs::remove_dir_all(&dir).ok();
    }

    // Reference shared by the indel wiring tests.
    const REF29: &[u8] = b"ACGTACGTAAAAAAAACGTACGTACGTAC";

    // Build a BAM whose reads carry explicit CIGARs, so indel-bearing reads can be constructed directly.
    fn write_cigar_bam(
        path: &Path,
        contig: &str,
        contig_len: usize,
        reads: &[(&str, i64, &[u8], Vec<Cigar>, bool)],
    ) {
        let mut header = Header::new();
        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", contig);
        sq.push_tag(b"LN", &contig_len.to_string());
        header.push_record(&sq);
        let mut w = Writer::from_path(path, &header, Format::Bam).unwrap();
        for (name, pos, seq, ops, reverse) in reads {
            let mut rec = Record::new();
            let quals = vec![30u8; seq.len()];
            rec.set(name.as_bytes(), Some(&CigarString(ops.clone())), seq, &quals);
            rec.set_tid(0);
            rec.set_pos(*pos);
            rec.set_mapq(60);
            rec.set_flags(if *reverse { 16 } else { 0 });
            w.write(&rec).unwrap();
        }
    }

    fn indel_refs() -> HashMap<Vec<u8>, Vec<u8>> {
        let mut refs = HashMap::new();
        refs.insert(b"chrM".to_vec(), REF29.to_vec());
        refs
    }

    fn indels_on() -> indel::IndelOpts {
        indel::IndelOpts { enabled: true, ..Default::default() }
    }

    fn indels_hp_on() -> indel::IndelOpts {
        indel::IndelOpts { enabled: true, hp_length_model: true, ..Default::default() }
    }

    const CLEAN_SEQ: &[u8] = b"ACGTACGTAAAAAAAACGTACGTACGTAC"; // == REF29, 29bp
    // An extra A inside the run: 16M 1I 13M.
    const HP_INS_SEQ: &[u8] = b"ACGTACGTAAAAAAAAACGTACGTACGTAC"; // 30bp
    // A spurious G after index 18, in unique sequence: 19M 1I 10M.
    const UNIQ_INS_SEQ: &[u8] = b"ACGTACGTAAAAAAAACGTGACGTACGTAC"; // 30bp

    // One A missing from the run: 8M 1D 20M.
    const HP_DEL_SEQ: &[u8] = b"ACGTACGTAAAAAAACGTACGTACGTAC"; // 28bp

    fn clean_cigar() -> Vec<Cigar> { vec![Cigar::Match(29)] }
    fn hp_del_cigar() -> Vec<Cigar> {
        vec![Cigar::Match(8), Cigar::Del(1), Cigar::Match(20)]
    }
    fn hp_ins_cigar() -> Vec<Cigar> { vec![Cigar::Match(16), Cigar::Ins(1), Cigar::Match(13)] }
    fn uniq_ins_cigar() -> Vec<Cigar> { vec![Cigar::Match(19), Cigar::Ins(1), Cigar::Match(10)] }

    const NAMES: [&str; 9] = ["r0","r1","r2","r3","r4","r5","r6","r7","r8"];

    #[test]
    fn indels_disabled_by_default_leaves_indel_reads_untouched() {
        let dir = std::env::temp_dir().join("himito_denoise_indel_off");
        std::fs::create_dir_all(&dir).unwrap();
        let bam = dir.join("in.bam");
        let mut reads: Vec<(&str, i64, &[u8], Vec<Cigar>, bool)> = vec![];
        for (i, n) in NAMES.iter().enumerate() {
            reads.push((n, 0, CLEAN_SEQ, clean_cigar(), i % 2 == 0));
        }
        reads.push(("rins", 0, UNIQ_INS_SEQ, uniq_ins_cigar(), false));
        write_cigar_bam(&bam, "chrM", 29, &reads);

        let off = indel::IndelOpts::default();
        assert!(!off.enabled);
        let (_corr, models, stats) =
            compute_corrections(&bam, &indel_refs(), 2, 0.01, SB_P, HOM_VAF, &off).unwrap();
        assert!(models.is_empty(), "no sites may be built while --indels is off");
        assert_eq!(stats.indel_sites_examined, 0);
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn site_depth_counts_each_spanning_read_exactly_once() {
        // Double-count trap: the read reports Indel::None at the normalized column and its insertion at the anchor.
        let dir = std::env::temp_dir().join("himito_denoise_indel_depth");
        std::fs::create_dir_all(&dir).unwrap();
        let bam = dir.join("in.bam");
        let mut reads: Vec<(&str, i64, &[u8], Vec<Cigar>, bool)> = vec![];
        for (i, n) in NAMES.iter().enumerate() {
            reads.push((n, 0, HP_INS_SEQ, hp_ins_cigar(), i % 2 == 0));
        }
        reads.push(("rdrop", 0, CLEAN_SEQ, clean_cigar(), false));
        write_cigar_bam(&bam, "chrM", 29, &reads);

        let (_corr, models, _stats) =
            compute_corrections(&bam, &indel_refs(), 2, 0.01, SB_P, HOM_VAF, &indels_on()).unwrap();

        let list = models.get(&b"chrM".to_vec()).expect("chrM must have sites");
        let (_, model) = list.iter().find(|(p, _)| *p == 7).expect("site at 7");
        let ins = indel::Allele::Ins(b"A".to_vec());
        let f_ins = model.kept.iter().find(|(a, _)| *a == ins).expect("insertion kept").1;
        // 9 of 10 reads carry it.
        assert!(
            (f_ins - 0.9).abs() < 1e-9,
            "insertion frequency should be 0.9, got {f_ins} (depth double-counted?)"
        );
        assert_eq!(model.context_len, 8, "site 7 sits in an 8bp homopolymer run");
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn hp_length_votes_give_the_same_site_frequencies_as_the_event_votes() {
        // Depth is untouched: the model stays Ins("A") at 0.9 in an 8bp context.
        let dir = std::env::temp_dir().join("himito_denoise_indel_hp_votes");
        std::fs::create_dir_all(&dir).unwrap();
        let bam = dir.join("in.bam");
        let mut reads: Vec<(&str, i64, &[u8], Vec<Cigar>, bool)> = vec![];
        for (i, n) in NAMES.iter().enumerate() {
            reads.push((n, 0, HP_INS_SEQ, hp_ins_cigar(), i % 2 == 0));
        }
        reads.push(("rdrop", 0, CLEAN_SEQ, clean_cigar(), false));
        write_cigar_bam(&bam, "chrM", 29, &reads);

        let (_corr, models, _stats) =
            compute_corrections(&bam, &indel_refs(), 2, 0.01, SB_P, HOM_VAF, &indels_hp_on())
                .unwrap();

        let list = models.get(&b"chrM".to_vec()).expect("chrM must have sites");
        let (_, model) = list.iter().find(|(p, _)| *p == 7).expect("site at 7");
        let ins = indel::Allele::Ins(b"A".to_vec());
        let f_ins = model.kept.iter().find(|(a, _)| *a == ins).expect("insertion kept").1;
        assert!(
            (f_ins - 0.9).abs() < 1e-9,
            "insertion frequency must still be 0.9 with the HP length model on, got {f_ins}"
        );
        assert_eq!(model.context_len, 8, "site 7 sits in an 8bp homopolymer run");
        assert_eq!(
            model.hp_run,
            Some((16u32, b'A')),
            "site 7 must be recognized as the HP site of REF29's A-run"
        );
        // Site 7 is the run's only normalized anchor: positions 8..=13 also have a >= 2bp A-run after them.
        let hp_sites: Vec<u32> =
            list.iter().filter(|(_, m)| m.hp_run.is_some()).map(|(p, _)| *p).collect();
        assert_eq!(hp_sites, vec![7], "the A-run must yield exactly one HP site");
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn an_undefined_hp_allele_does_not_count_as_reference_support() {
        // `hp_allele` returns `None` for a read that is UNINFORMATIVE about the run's length.
        let dir = std::env::temp_dir().join("himito_denoise_indel_hp_ref_support");
        std::fs::create_dir_all(&dir).unwrap();
        let bam = dir.join("in.bam");
        let mut reads: Vec<(&str, i64, &[u8], Vec<Cigar>, bool)> = vec![];
        for (i, n) in NAMES.iter().enumerate() {
            reads.push((n, 0, HP_INS_SEQ, hp_ins_cigar(), i % 2 == 0));
        }
        // The ONLY read that observes the reference length.
        reads.push(("rdrop", 0, CLEAN_SEQ, clean_cigar(), false));
        // Reference end 18: past the depth gate (17), short of the vote gate (21).
        reads.push(("rflank0", 0, &CLEAN_SEQ[..18], vec![Cigar::Match(18)], false));
        reads.push(("rflank1", 0, &CLEAN_SEQ[..18], vec![Cigar::Match(18)], true));
        write_cigar_bam(&bam, "chrM", 29, &reads);

        let (_corr, models, _stats) =
            compute_corrections(&bam, &indel_refs(), 2, 0.01, SB_P, HOM_VAF, &indels_hp_on())
                .unwrap();

        let list = models.get(&b"chrM".to_vec()).expect("chrM must have sites");
        let (_, model) = list.iter().find(|(p, _)| *p == 7).expect("site at 7");
        let ins = indel::Allele::Ins(b"A".to_vec());
        let f_ins = model.kept.iter().find(|(a, _)| *a == ins).expect("insertion kept").1;
        assert!(
            (f_ins - 0.9).abs() < 1e-9,
            "9 carriers against 1 reference observer is 0.9; the two uninformative \
             reads must not be counted as reference support, got {f_ins}"
        );
        let f_ref = model
            .kept
            .iter()
            .find(|(a, _)| *a == indel::Allele::Ref)
            .expect("REF is always kept")
            .1;
        assert!(
            (f_ref - 0.1).abs() < 1e-9,
            "exactly one read observes the reference length, got {f_ref}"
        );
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn site_depth_excludes_reads_that_end_inside_the_homopolymer() {
        // Adding 11 such reads must keep the insertion frequency at site 7 at 0.9.
        let dir = std::env::temp_dir().join("himito_denoise_indel_depth_short");
        std::fs::create_dir_all(&dir).unwrap();
        let bam = dir.join("in.bam");
        let mut reads: Vec<(&str, i64, &[u8], Vec<Cigar>, bool)> = vec![];
        for (i, n) in NAMES.iter().enumerate() {
            reads.push((n, 0, HP_INS_SEQ, hp_ins_cigar(), i % 2 == 0));
        }
        reads.push(("rdrop", 0, CLEAN_SEQ, clean_cigar(), false));
        // Reference end 13: five bases into the run, so the run's true length is invisible to these reads.
        const SHORT: [&str; 10] = ["s0","s1","s2","s3","s4","s5","s6","s7","s8","s9"];
        for (i, n) in SHORT.iter().enumerate() {
            reads.push((n, 0, &CLEAN_SEQ[..13], vec![Cigar::Match(13)], i % 2 == 0));
        }
        // One non-spanning read that also CARRIES an extra A (13M 1I 1M, reference end 14).
        let sins_seq: Vec<u8> = {
            let mut s = CLEAN_SEQ[..13].to_vec();
            s.push(b'A');
            s.push(CLEAN_SEQ[13]);
            s
        };
        reads.push((
            "sins",
            0,
            sins_seq.as_slice(),
            vec![Cigar::Match(13), Cigar::Ins(1), Cigar::Match(1)],
            true,
        ));
        write_cigar_bam(&bam, "chrM", 29, &reads);

        let (_corr, models, stats) =
            compute_corrections(&bam, &indel_refs(), 2, 0.01, SB_P, HOM_VAF, &indels_on()).unwrap();

        let list = models.get(&b"chrM".to_vec()).expect("chrM must have sites");
        let (_, model) = list.iter().find(|(p, _)| *p == 7).expect("site at 7");
        let ins = indel::Allele::Ins(b"A".to_vec());
        let f_ins = model
            .kept
            .iter()
            .find(|(a, _)| *a == ins)
            .unwrap_or_else(|| {
                panic!("insertion must survive candidacy; diluted model = {model:?}")
            })
            .1;
        assert!(
            (f_ins - 0.9).abs() < 1e-9,
            "insertion frequency should be 0.9 (only the 10 run-spanning reads count), \
             got {f_ins} -- reads ending inside the run were counted as REF"
        );
        assert_eq!(model.context_len, 8, "site 7 sits in an 8bp homopolymer run");
        assert_eq!(
            stats.indel_events_not_spanning_repeat, 1,
            "`sins`'s vote is the one event dropped for not spanning the run"
        );
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn alt_vote_is_gated_at_the_normalized_site_not_the_observed_column() {
        // `15M 1I 1M` from POS 0 anchors its insertion on REF29[14] -- the PENULTIMATE base of the A-run.
        let dir = std::env::temp_dir().join("himito_denoise_indel_vote_gate_norm_pos");
        std::fs::create_dir_all(&dir).unwrap();
        let bam = dir.join("in.bam");

        // ref[0..15] + an extra 'A' + ref[15]: 17 query bases, 15M 1I 1M, end 16.
        let pen_seq: Vec<u8> = {
            let mut s = REF29[0..15].to_vec();
            s.push(b'A');
            s.push(REF29[15]);
            s
        };
        let pen_cigar = vec![Cigar::Match(15), Cigar::Ins(1), Cigar::Match(1)];

        let mut reads: Vec<(&str, i64, &[u8], Vec<Cigar>, bool)> = vec![];
        // The only read reaching past the whole run plus its anchor: depth 1 at 7.
        reads.push(("rspan", 0, CLEAN_SEQ, clean_cigar(), false));
        for (i, n) in NAMES.iter().enumerate() {
            reads.push((n, 0, pen_seq.as_slice(), pen_cigar.clone(), i % 2 == 0));
        }
        write_cigar_bam(&bam, "chrM", 29, &reads);

        // Must not panic: nine votes against a depth of one is the invariant break.
        let (_corr, models, stats) =
            compute_corrections(&bam, &indel_refs(), 2, 0.01, SB_P, HOM_VAF, &indels_on()).unwrap();

        assert_eq!(
            stats.indel_events_not_spanning_repeat,
            NAMES.len() as u64,
            "every 15M1I1M read's vote must be dropped at the NORMALIZED site 7, \
             whose run it does not span -- not waved through by column 14's laxer bound"
        );
        let kept_ins_at_7 = models.get(&b"chrM".to_vec()).is_some_and(|list| {
            list.iter().any(|(p, m)| {
                *p == 7
                    && m.kept.iter().any(|(a, _)| *a == indel::Allele::Ins(b"A".to_vec()))
            })
        });
        assert!(
            !kept_ins_at_7,
            "site 7 must not keep an insertion cast only by reads that never spanned its run"
        );
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn indel_site_with_only_nonspanning_alt_votes_counts_dropped_depth() {
        // The inserted A left-normalizes onto site 7, and that site's A-run requires end >= 17, so every vote is dropped.
        let dir = std::env::temp_dir().join("himito_denoise_indel_site_no_span");
        std::fs::create_dir_all(&dir).unwrap();
        let bam = dir.join("in.bam");

        let pen_seq: Vec<u8> = {
            let mut s = REF29[0..15].to_vec();
            s.push(b'A');
            s.push(REF29[15]);
            s
        };
        let pen_cigar = vec![Cigar::Match(15), Cigar::Ins(1), Cigar::Match(1)];

        const SHORT: [&str; 4] = ["p0", "p1", "p2", "p3"];
        let mut reads: Vec<(&str, i64, &[u8], Vec<Cigar>, bool)> = vec![];
        for (i, n) in SHORT.iter().enumerate() {
            reads.push((n, 0, pen_seq.as_slice(), pen_cigar.clone(), i % 2 == 0));
        }
        write_cigar_bam(&bam, "chrM", 29, &reads);

        let (_corr, _models, stats) =
            compute_corrections(&bam, &indel_refs(), 2, 0.01, SB_P, HOM_VAF, &indels_on()).unwrap();

        assert!(
            stats.indel_events_not_spanning_repeat >= 1,
            "each short read's insertion must be dropped by the span gate at site 7, \
             got {}",
            stats.indel_events_not_spanning_repeat
        );
        assert!(
            stats.indel_sites_dropped_no_spanning_depth >= 1,
            "site 7 had ALT votes and no spanning depth, so it must be counted, got {}",
            stats.indel_sites_dropped_no_spanning_depth
        );
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn models_are_sorted_by_position_within_each_contig() {
        // The record walk binary-searches this list; unsorted input silently drops sites from every read's span.
        let dir = std::env::temp_dir().join("himito_denoise_indel_sorted");
        std::fs::create_dir_all(&dir).unwrap();
        let bam = dir.join("in.bam");
        let both: Vec<Cigar> = vec![
            Cigar::Match(16), Cigar::Ins(1), Cigar::Match(4), Cigar::Ins(1), Cigar::Match(9),
        ];
        // Extra A in the run and a spurious G after index 19.
        let both_seq: &[u8] = b"ACGTACGTAAAAAAAAACGTAGCGTACGTAC"; // 31bp
        let mut reads: Vec<(&str, i64, &[u8], Vec<Cigar>, bool)> = vec![];
        for (i, n) in NAMES.iter().enumerate() {
            reads.push((n, 0, both_seq, both.clone(), i % 2 == 0));
        }
        reads.push(("rclean", 0, CLEAN_SEQ, clean_cigar(), false));
        write_cigar_bam(&bam, "chrM", 29, &reads);

        let (_corr, models, _stats) =
            compute_corrections(&bam, &indel_refs(), 2, 0.01, SB_P, HOM_VAF, &indels_on()).unwrap();
        let list = models.get(&b"chrM".to_vec()).expect("chrM must have sites");
        // Assert the exact expected positions rather than checking ascending order pairwise: with only two elements.
        let positions: Vec<u32> = list.iter().map(|(p, _)| *p).collect();
        assert_eq!(positions, vec![7, 19], "sites must be exactly [7, 19], in that order");
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn lone_spurious_insertion_is_removed_from_a_substitution_clean_column() {
        // Every aligned base matches the reference, so `counts.non_ref() == 0` at every column and the SNV path skips them all.
        let dir = std::env::temp_dir().join("himito_denoise_indel_cleancol");
        std::fs::create_dir_all(&dir).unwrap();
        let bam = dir.join("in.bam");
        let outb = dir.join("out.bam");
        let mut reads: Vec<(&str, i64, &[u8], Vec<Cigar>, bool)> = vec![];
        for (i, n) in NAMES.iter().enumerate() {
            reads.push((n, 0, CLEAN_SEQ, clean_cigar(), i % 2 == 0));
        }
        reads.push(("rins", 0, UNIQ_INS_SEQ, uniq_ins_cigar(), false));
        write_cigar_bam(&bam, "chrM", 29, &reads);

        let refs = indel_refs();
        let iopts = indels_on();
        let (corr, models, mut stats) =
            compute_corrections(&bam, &refs, 2, 0.01, SB_P, HOM_VAF, &iopts).unwrap();
        assert_eq!(stats.columns_skipped, stats.columns_processed,
            "test premise: every column is substitution-clean");
        apply_corrections(&bam, &outb, &corr, &refs, &models, &iopts, &mut stats).unwrap();

        let mut r = Reader::from_path(&outb).unwrap();
        let mut seen = false;
        for rec in r.records() {
            let rec = rec.unwrap();
            if rec.qname() == b"rins" {
                assert_eq!(rec.seq().as_bytes(), CLEAN_SEQ, "spurious insertion must be removed");
                assert_eq!(rec.cigar().to_string(), "29M");
                seen = true;
            } else {
                assert_eq!(rec.seq().as_bytes(), CLEAN_SEQ, "clean reads must be untouched");
            }
        }
        assert!(seen);
        assert_eq!(stats.reads_reassigned_indel_to_ref, 1);
        // reverting a 1bp insertion removes 1 base from the read.
        assert_eq!(stats.indel_bases_removed, 1);
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn homoplasmic_indel_survives_and_repairs_the_read_that_dropped_it() {
        // 9 of 10 reads carry a 1bp insertion inside the 8bp A-run, so eps is large enough for the repair direction to fire.
        let dir = std::env::temp_dir().join("himito_denoise_indel_repair");
        std::fs::create_dir_all(&dir).unwrap();
        let bam = dir.join("in.bam");
        let outb = dir.join("out.bam");
        let mut reads: Vec<(&str, i64, &[u8], Vec<Cigar>, bool)> = vec![];
        for (i, n) in NAMES.iter().enumerate() {
            reads.push((n, 0, HP_INS_SEQ, hp_ins_cigar(), i % 2 == 0));
        }
        reads.push(("rdrop", 0, CLEAN_SEQ, clean_cigar(), false));
        write_cigar_bam(&bam, "chrM", 29, &reads);

        let refs = indel_refs();
        let iopts = indels_on();
        let (corr, models, mut stats) =
            compute_corrections(&bam, &refs, 2, 0.01, SB_P, HOM_VAF, &iopts).unwrap();
        assert_eq!(stats.indel_alt_alleles_kept, 1, "the 90% insertion must survive candidacy");
        apply_corrections(&bam, &outb, &corr, &refs, &models, &iopts, &mut stats).unwrap();

        let mut r = Reader::from_path(&outb).unwrap();
        let mut repaired = false;
        for rec in r.records() {
            let rec = rec.unwrap();
            assert_eq!(rec.seq().as_bytes(), HP_INS_SEQ, "every read should carry the insertion");
            if rec.qname() == b"rdrop" {
                repaired = true;
            }
        }
        assert!(repaired);
        assert_eq!(stats.reads_reassigned_ref_to_ins, 1);
        assert_eq!(stats.indel_bases_inserted, 1);
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn read_ending_too_close_to_a_site_is_skipped_by_the_flank_guard() {
        // Same site, but the read that dropped the insertion ends before the needed match context.
        let dir = std::env::temp_dir().join("himito_denoise_indel_guard");
        std::fs::create_dir_all(&dir).unwrap();
        let bam = dir.join("in.bam");
        let outb = dir.join("out.bam");
        let mut reads: Vec<(&str, i64, &[u8], Vec<Cigar>, bool)> = vec![];
        for (i, n) in NAMES.iter().enumerate() {
            reads.push((n, 0, HP_INS_SEQ, hp_ins_cigar(), i % 2 == 0));
        }
        let short_seq: &[u8] = b"ACGTACGTAAAAAA"; // ref[0..14]
        reads.push(("rshort", 0, short_seq, vec![Cigar::Match(14)], false));
        write_cigar_bam(&bam, "chrM", 29, &reads);

        let refs = indel_refs();
        let iopts = indels_on();
        let (corr, models, mut stats) =
            compute_corrections(&bam, &refs, 2, 0.01, SB_P, HOM_VAF, &iopts).unwrap();
        apply_corrections(&bam, &outb, &corr, &refs, &models, &iopts, &mut stats).unwrap();

        let mut r = Reader::from_path(&outb).unwrap();
        for rec in r.records() {
            let rec = rec.unwrap();
            if rec.qname() == b"rshort" {
                assert_eq!(rec.seq().as_bytes(), short_seq, "guarded read must be untouched");
            }
        }
        assert!(stats.assignments_skipped_span_guard >= 1, "the skip must be counted");
        assert_eq!(stats.reads_reassigned_ref_to_ins, 0);
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn read_carrying_an_oversized_indel_is_not_assigned_nearby() {
        // IMPORTANT 3: the original version of this fixture placed the site in UNIQUE sequence (context_len = 1, eps = 0.01).
        let dir = std::env::temp_dir().join("himito_denoise_indel_large");
        std::fs::create_dir_all(&dir).unwrap();
        let bam = dir.join("in.bam");
        let outb = dir.join("out.bam");
        // ref[0..9] + 6 spurious G's (oversized, anchor = 8, inside the run) + ref[9..29].
        let mut big_seq: Vec<u8> = Vec::new();
        big_seq.extend_from_slice(&REF29[0..9]);
        big_seq.extend_from_slice(b"GGGGGG");
        big_seq.extend_from_slice(&REF29[9..29]);
        let big_cigar = vec![Cigar::Match(9), Cigar::Ins(6), Cigar::Match(20)];
        let mut reads: Vec<(&str, i64, &[u8], Vec<Cigar>, bool)> = vec![];
        for (i, n) in NAMES.iter().enumerate() {
            reads.push((n, 0, HP_INS_SEQ, hp_ins_cigar(), i % 2 == 0));
        }
        reads.push(("rbig", 0, big_seq.as_slice(), big_cigar.clone(), false));
        write_cigar_bam(&bam, "chrM", 29, &reads);

        let refs = indel_refs();
        let iopts = indels_on();
        let (corr, models, mut stats) =
            compute_corrections(&bam, &refs, 2, 0.01, SB_P, HOM_VAF, &iopts).unwrap();
        assert_eq!(
            stats.indel_alt_alleles_kept, 1,
            "test premise: the 90% insertion at site 7 must clear candidacy"
        );
        apply_corrections(&bam, &outb, &corr, &refs, &models, &iopts, &mut stats).unwrap();

        let mut r = Reader::from_path(&outb).unwrap();
        let mut seen = false;
        for rec in r.records() {
            let rec = rec.unwrap();
            if rec.qname() == b"rbig" {
                assert_eq!(rec.seq().as_bytes(), big_seq.as_slice(), "oversized indel must be untouched");
                assert_eq!(
                    rec.cigar().take(),
                    CigarString(big_cigar.clone()),
                    "CIGAR must be untouched"
                );
                seen = true;
            }
        }
        assert!(seen);
        assert!(
            stats.assignments_skipped_oversized_neighbor >= 1,
            "the mask's skip must be counted, proving it is what fired"
        );
        std::fs::remove_dir_all(&dir).ok();
    }

    // ---- normalized position vs CIGAR anchor mismatch ----

    #[test]
    fn spurious_homopolymer_insertion_is_reverted_using_its_own_cigar_anchor() {
        // 9 clean reads plus one carrying a lone extra A inside the 8bp A-run, anchored by the aligner at reference index 15.
        let dir = std::env::temp_dir().join("himito_denoise_indel_anchor_mismatch");
        std::fs::create_dir_all(&dir).unwrap();
        let bam = dir.join("in.bam");
        let outb = dir.join("out.bam");
        let mut reads: Vec<(&str, i64, &[u8], Vec<Cigar>, bool)> = vec![];
        for (i, n) in NAMES.iter().enumerate() {
            reads.push((n, 0, CLEAN_SEQ, clean_cigar(), i % 2 == 0));
        }
        reads.push(("rspur", 0, HP_INS_SEQ, hp_ins_cigar(), false));
        write_cigar_bam(&bam, "chrM", 29, &reads);

        let refs = indel_refs();
        let iopts = indels_on();
        let (corr, models, mut stats) =
            compute_corrections(&bam, &refs, 2, 0.01, SB_P, HOM_VAF, &iopts).unwrap();
        apply_corrections(&bam, &outb, &corr, &refs, &models, &iopts, &mut stats).unwrap();

        let mut r = Reader::from_path(&outb).unwrap();
        let mut seen = false;
        for rec in r.records() {
            let rec = rec.unwrap();
            if rec.qname() == b"rspur" {
                assert_eq!(
                    rec.seq().as_bytes(), CLEAN_SEQ,
                    "the spurious insertion must actually be removed from SEQ, not left in place"
                );
                assert_eq!(
                    rec.cigar().to_string(), "29M",
                    "CIGAR must actually collapse to a plain full-length match"
                );
                seen = true;
            } else {
                assert_eq!(rec.seq().as_bytes(), CLEAN_SEQ, "clean reads must be untouched");
            }
        }
        assert!(seen);
        assert_eq!(stats.reads_reassigned_indel_to_ref, 1);
        // this revert removes the 1 spurious base from the read.
        assert_eq!(stats.indel_bases_removed, 1);
        std::fs::remove_dir_all(&dir).ok();
    }

    // ---- depth/ALT sourced inconsistently across span ----

    #[test]
    fn indel_event_normalizing_outside_its_own_span_is_not_fabricated_into_a_site() {
        // Four reads start at POS 10 -- inside the 8bp A-run -- and each carries an extra A anchored at ref idx 10.
        let dir = std::env::temp_dir().join("himito_denoise_indel_out_of_span");
        std::fs::create_dir_all(&dir).unwrap();
        let bam = dir.join("in.bam");
        let outb = dir.join("out.bam");

        // seq = ref[10] + inserted 'A' + ref[11..29]; anchor = 10.
        let mut seq: Vec<u8> = Vec::new();
        seq.push(REF29[10]);
        seq.push(b'A');
        seq.extend_from_slice(&REF29[11..29]);
        let cigar = vec![Cigar::Match(1), Cigar::Ins(1), Cigar::Match(18)];

        let mut reads: Vec<(&str, i64, &[u8], Vec<Cigar>, bool)> = vec![];
        // One read spanning the WHOLE contig, clean, reference-matching -- this is the read that must come out unchanged.
        reads.push(("rfull", 0, CLEAN_SEQ, clean_cigar(), false));
        for (i, n) in ["r0", "r1", "r2", "r3"].iter().enumerate() {
            reads.push((n, 10, seq.as_slice(), cigar.clone(), i % 2 == 0));
        }
        write_cigar_bam(&bam, "chrM", 29, &reads);

        let refs = indel_refs();
        let iopts = indels_on();
        let (corr, models, mut stats) =
            compute_corrections(&bam, &refs, 2, 0.01, SB_P, HOM_VAF, &iopts).unwrap();

        assert_eq!(
            stats.indel_events_out_of_span, 4,
            "all four out-of-span votes must be dropped and counted"
        );
        // This fixture's only would-be evidence for site 7 is the four dropped votes above.
        let no_site_7 = models
            .get(&b"chrM".to_vec())
            .map_or(true, |list| list.iter().all(|(p, _)| *p != 7));
        assert!(no_site_7, "site 7 must not be fabricated from reads that never span it");

        apply_corrections(&bam, &outb, &corr, &refs, &models, &iopts, &mut stats).unwrap();
        let mut r = Reader::from_path(&outb).unwrap();
        for rec in r.records() {
            let rec = rec.unwrap();
            if rec.qname() == b"rfull" {
                assert_eq!(
                    rec.seq().as_bytes(), CLEAN_SEQ,
                    "the fully-spanning reference-matching read must be unchanged"
                );
                assert_eq!(rec.cigar().to_string(), "29M");
            }
        }
        std::fs::remove_dir_all(&dir).ok();
    }

    // ---- pass-2 per-read vote cap and edit guards ----

    // A read whose own insertion and deletion BOTH left-normalize to site 7 -- short tandem repeats are exactly.
    fn dup_cigar() -> Vec<Cigar> {
        vec![Cigar::Match(9), Cigar::Ins(1), Cigar::Match(4), Cigar::Del(1), Cigar::Match(15)]
    }

    #[test]
    fn read_casting_two_events_that_normalize_to_the_same_site_votes_once() {
        // without a per-read cap, `rdup` would cast TWO ALT votes at site 7 (one for its Ins("A"), one for its Del(1)).
        let dir = std::env::temp_dir().join("himito_denoise_indel_dup_site_vote");
        std::fs::create_dir_all(&dir).unwrap();
        let bam = dir.join("in.bam");
        let mut reads: Vec<(&str, i64, &[u8], Vec<Cigar>, bool)> = vec![];
        for (i, n) in NAMES.iter().enumerate() {
            reads.push((n, 0, CLEAN_SEQ, clean_cigar(), i % 2 == 0));
        }
        reads.push(("rdup", 0, CLEAN_SEQ, dup_cigar(), false));
        write_cigar_bam(&bam, "chrM", 29, &reads);

        let (_corr, models, stats) =
            compute_corrections(&bam, &indel_refs(), 2, 0.01, SB_P, HOM_VAF, &indels_on()).unwrap();

        assert_eq!(
            stats.indel_events_duplicate_site, 1,
            "the second vote from the same read at the same site must be dropped and counted"
        );
        let list = models.get(&b"chrM".to_vec()).expect("chrM must have a site at 7");
        let (_, model) = list.iter().find(|(p, _)| *p == 7).expect("site at 7 must exist");
        let ins = indel::Allele::Ins(b"A".to_vec());
        let del = indel::Allele::Del(1);
        let ins_vaf = model.obs_vaf.iter().find(|(a, _)| *a == ins).map(|(_, v)| *v);
        assert!(
            (ins_vaf.expect("the FIRST vote (Ins) must be recorded") - 0.1).abs() < 1e-9,
            "the recorded vote must reflect exactly one voting read (1/10 depth), got {ins_vaf:?}"
        );
        assert!(
            model.obs_vaf.iter().all(|(a, _)| *a != del),
            "the SECOND vote (Del) from the same read must never appear at this site at all"
        );
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn read_with_ambiguous_carried_events_at_one_site_is_left_untouched() {
        // In pass 2, `rdup`'s CIGAR yields two different events that both normalize to site 7.
        let dir = std::env::temp_dir().join("himito_denoise_indel_ambiguous_carry");
        std::fs::create_dir_all(&dir).unwrap();
        let bam = dir.join("in.bam");
        let outb = dir.join("out.bam");
        let mut reads: Vec<(&str, i64, &[u8], Vec<Cigar>, bool)> = vec![];
        for (i, n) in NAMES.iter().enumerate() {
            reads.push((n, 0, CLEAN_SEQ, clean_cigar(), i % 2 == 0));
        }
        reads.push(("rdup", 0, CLEAN_SEQ, dup_cigar(), false));
        write_cigar_bam(&bam, "chrM", 29, &reads);

        let refs = indel_refs();
        let iopts = indels_on();
        let (corr, models, mut stats) =
            compute_corrections(&bam, &refs, 2, 0.01, SB_P, HOM_VAF, &iopts).unwrap();
        apply_corrections(&bam, &outb, &corr, &refs, &models, &iopts, &mut stats).unwrap();

        assert!(
            stats.assignments_skipped_ambiguous_carry >= 1,
            "the ambiguous-carry skip must be counted"
        );
        let mut r = Reader::from_path(&outb).unwrap();
        let mut seen = false;
        for rec in r.records() {
            let rec = rec.unwrap();
            if rec.qname() == b"rdup" {
                assert_eq!(rec.seq().as_bytes(), CLEAN_SEQ, "SEQ must be untouched");
                assert_eq!(
                    rec.cigar().take(),
                    CigarString(dup_cigar()),
                    "CIGAR must be untouched -- no arbitrary pick between the two events"
                );
                seen = true;
            }
        }
        assert!(seen);
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn gained_insertion_at_a_column_the_read_deletes_is_rejected() {
        // `rdel3`'s own D(3) is anchored at 9; this fixture's reference is a 30bp homopolymer.
        let dir = std::env::temp_dir().join("himito_denoise_indel_gained_ins_in_own_del");
        std::fs::create_dir_all(&dir).unwrap();
        let inb = dir.join("in.bam");
        let outb = dir.join("out.bam");

        // 10M 3D 10M: ref[0..10) matched, ref[10..13) deleted, ref[13..23) matched.
        let seq = vec![b'A'; 20]; // 10 (M) + 10 (M) query bases; D consumes none
        let cigar = vec![Cigar::Match(10), Cigar::Del(3), Cigar::Match(10)];
        write_cigar_bam(&inb, "chrM", 30, &[("rdel3", 0, seq.as_slice(), cigar.clone(), false)]);

        let mut refs = HashMap::new();
        refs.insert(b"chrM".to_vec(), vec![b'A'; 30]);

        let model = indel::IndelSiteModel {
            kept: vec![(indel::Allele::Ref, 0.1), (indel::Allele::Ins(b"A".to_vec()), 0.9)],
            eps: 0.15,
            context_len: 1,
            strand_biased: 0,
            obs_vaf: vec![(indel::Allele::Ref, 0.1), (indel::Allele::Ins(b"A".to_vec()), 0.9)],
            strand_rejected: vec![],
            hp_run: None,
        };
        let mut models = IndelModels::new();
        models.entry(b"chrM".to_vec()).or_default().push((12u32, model));

        let corr: Corrections = HashMap::new();
        let mut stats = DenoiseStats::default();
        apply_corrections(
            &inb, &outb, &corr, &refs, &models, &indels_on(), &mut stats,
        )
        .unwrap();
        let (proc, modified) = (stats.reads_processed, stats.reads_modified);
        assert_eq!(proc, 1);
        assert_eq!(modified, 0, "the read must not be rewritten -- it has no query base at 12");

        let mut r = Reader::from_path(&outb).unwrap();
        let rec = r.records().next().unwrap().unwrap();
        assert_eq!(rec.seq().as_bytes(), seq.as_slice(), "SEQ must be untouched");
        assert_eq!(rec.cigar().take(), CigarString(cigar), "CIGAR must be untouched");
        assert_eq!(stats.reads_reassigned_ref_to_ins, 0, "no reassignment must be counted for a rejected edit");
        assert_eq!(stats.indel_bases_inserted, 0, "no bases-inserted count for a rejected edit");
        assert!(
            stats.assignments_skipped_span_guard >= 1,
            "the missing-aligned-base rejection must be counted"
        );
        std::fs::remove_dir_all(&dir).ok();
    }

    // Site 7 sits before REF29's 8bp A-run; Ins("A") dominates 90/10 over Ref.
    fn hp_ins_site_model() -> indel::IndelSiteModel {
        indel::IndelSiteModel {
            kept: vec![(indel::Allele::Ref, 0.1), (indel::Allele::Ins(b"A".to_vec()), 0.9)],
            eps: 0.15,
            context_len: 8,
            strand_biased: 0,
            obs_vaf: vec![(indel::Allele::Ref, 0.1), (indel::Allele::Ins(b"A".to_vec()), 0.9)],
            strand_rejected: vec![],
            hp_run: None,
        }
    }

    /// Del(1) dominates 0.89/0.11 over Ref.
    /// The event-based model refuses `Ins -> Del`; the length model must move it.
    fn hp_del_site_model() -> indel::IndelSiteModel {
        let ins = indel::Allele::Ins(b"A".to_vec());
        indel::IndelSiteModel {
            kept: vec![(indel::Allele::Ref, 0.111), (indel::Allele::Del(1), 0.889)],
            eps: indel::error_rate(8, &indels_hp_on()),
            context_len: 8,
            strand_biased: 0,
            obs_vaf: vec![
                (indel::Allele::Ref, 0.1),
                (ins.clone(), 0.1),
                (indel::Allele::Del(1), 0.8),
            ],
            strand_rejected: vec![ins],
            hp_run: Some((16, b'A')),
        }
    }

    #[test]
    fn hp_site_move_rewrites_the_run_and_drops_substitutions_inside_it() {
        // The rewrite is done by `rewrite_hp_run` BEFORE `rewrite_read` and never enters `edits.indels`.
        let dir = std::env::temp_dir().join("himito_denoise_indel_hp_move");
        std::fs::create_dir_all(&dir).unwrap();
        let inb = dir.join("in.bam");
        let outb = dir.join("out.bam");
        write_cigar_bam(
            &inb,
            "chrM",
            29,
            &[("rins", 0, HP_INS_SEQ, hp_ins_cigar(), false)],
        );

        let mut models = IndelModels::new();
        models.entry(b"chrM".to_vec()).or_default().push((7u32, hp_del_site_model()));

        // One substitution inside the run (reference 10) and one outside it (20).
        let mut corr: Corrections = HashMap::new();
        corr.insert(
            b"rins".to_vec(),
            indel::ReadEdits { subs: vec![(10, b'G'), (20, b'T')], indels: vec![] },
        );
        let mut stats = DenoiseStats::default();
        apply_corrections(
            &inb, &outb, &corr, &indel_refs(), &models, &indels_hp_on(), &mut stats,
        )
        .unwrap();

        assert_eq!(stats.reads_modified, 1);
        assert_eq!(
            stats.reads_reassigned_indel_to_indel, 1,
            "Ins(1) -> Del(1) is one indel->indel reassignment"
        );
        assert_eq!(
            stats.indel_bases_removed, 2,
            "the inserted base leaves the read and one reference base becomes a deletion"
        );
        assert_eq!(
            stats.reassignments_by_hp_length.get(8).copied().unwrap_or(0),
            1,
            "booked against the run length R = 8"
        );
        assert_eq!(
            stats.subs_dropped_in_hp_rewrite, 1,
            "only the substitution INSIDE the run is dropped"
        );
        assert_eq!(stats.indel_edits_emitted_but_not_applied, 0);

        let mut r = Reader::from_path(&outb).unwrap();
        let rec = r.records().next().unwrap().unwrap();
        let cigar = rec.cigar().take();
        assert_eq!(
            cigar,
            CigarString(vec![Cigar::Match(8), Cigar::Del(1), Cigar::Match(20)]),
            "the run must come out in canonical left-aligned form"
        );
        // samtools-style validity: SEQ length == the CIGAR's query consumption.
        let query: u32 = cigar
            .iter()
            .map(|op| match op {
                Cigar::Match(n) | Cigar::Ins(n) | Cigar::SoftClip(n) | Cigar::Equal(n)
                | Cigar::Diff(n) => *n,
                _ => 0,
            })
            .sum();
        let seq = rec.seq().as_bytes();
        assert_eq!(seq.len() as u32, query, "SEQ length must match the CIGAR");
        assert_eq!(seq.len(), rec.qual().len(), "SEQ and QUAL must agree");
        assert_eq!(
            seq[19], b'T',
            "the substitution OUTSIDE the run must still be applied, at reference 20"
        );
        std::fs::remove_dir_all(&dir).ok();
    }

    /// The site keeps `Del(1)`, so the lone `Ins(1)` read carries a length nobody else has.
    fn hp_wrong_direction_reads() -> Vec<(&'static str, i64, &'static [u8], Vec<Cigar>, bool)> {
        let mut reads: Vec<(&str, i64, &[u8], Vec<Cigar>, bool)> = vec![];
        for (i, n) in NAMES[..8].iter().enumerate() {
            reads.push((n, 0, HP_DEL_SEQ, hp_del_cigar(), i % 2 == 0));
        }
        reads.push(("rins", 0, HP_INS_SEQ, hp_ins_cigar(), false));
        reads.push(("rref", 0, CLEAN_SEQ, clean_cigar(), true));
        reads
    }

    /// SEQ length must equal the CIGAR's query consumption, as `samtools quickcheck` requires.
    fn assert_samtools_valid(rec: &Record) {
        let query: u32 = rec
            .cigar()
            .iter()
            .map(|op| match op {
                Cigar::Match(n) | Cigar::Ins(n) | Cigar::SoftClip(n) | Cigar::Equal(n)
                | Cigar::Diff(n) => *n,
                _ => 0,
            })
            .sum();
        let seq = rec.seq().as_bytes();
        assert_eq!(
            seq.len() as u32,
            query,
            "SEQ length must equal the CIGAR's query consumption for {}",
            String::from_utf8_lossy(rec.qname())
        );
        assert_eq!(seq.len(), rec.qual().len(), "SEQ and QUAL lengths must agree");
    }

    fn run_denoise_pair(
        dir: &Path,
        iopts: &indel::IndelOpts,
    ) -> (DenoiseStats, HashMap<Vec<u8>, (Vec<u8>, CigarString)>) {
        std::fs::create_dir_all(dir).unwrap();
        let inb = dir.join("in.bam");
        let outb = dir.join("out.bam");
        write_cigar_bam(&inb, "chrM", 29, &hp_wrong_direction_reads());
        let refs = indel_refs();
        let (corr, models, mut stats) =
            compute_corrections(&inb, &refs, 2, 0.01, SB_P, HOM_VAF, iopts).unwrap();
        apply_corrections(&inb, &outb, &corr, &refs, &models, iopts, &mut stats).unwrap();

        let mut out = HashMap::new();
        let mut r = Reader::from_path(&outb).unwrap();
        for rec in r.records() {
            let rec = rec.unwrap();
            assert_samtools_valid(&rec);
            out.insert(rec.qname().to_vec(), (rec.seq().as_bytes(), rec.cigar().take()));
        }
        (stats, out)
    }

    #[test]
    fn hp_wrong_direction_read_is_refused_with_the_model_off() {
        // Today's behaviour, pinned: `Ins(1)` -> `Del(1)` is cross-kind, which `rewrite_read` cannot express.
        let dir = std::env::temp_dir().join("himito_denoise_indel_hp_e2e_off");
        let (stats, out) = run_denoise_pair(&dir, &indels_on());

        let (seq, cigar) = out.get(&b"rins".to_vec()).expect("rins must be written out");
        assert_eq!(seq.as_slice(), HP_INS_SEQ, "SEQ must be untouched");
        assert_eq!(*cigar, CigarString(hp_ins_cigar()), "CIGAR must be untouched");
        assert!(
            stats.assignments_skipped_unsupported_transition >= 1,
            "the refusal must be counted"
        );
        assert_eq!(stats.subs_dropped_in_hp_rewrite, 0, "no HP rewrite may happen");
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn hp_wrong_direction_read_is_corrected_with_the_model_on() {
        // The `Ins(1)` read is rewritten to the site's kept length, in canonical left-aligned form.
        let dir = std::env::temp_dir().join("himito_denoise_indel_hp_e2e_on");
        let (stats, out) = run_denoise_pair(&dir, &indels_hp_on());

        let (seq, cigar) = out.get(&b"rins".to_vec()).expect("rins must be written out");
        assert_eq!(*cigar, CigarString(hp_del_cigar()), "the run must become 8M 1D 20M");
        assert_eq!(seq.as_slice(), HP_DEL_SEQ, "the read must carry a 7bp A-run");
        assert_eq!(
            stats.assignments_skipped_unsupported_transition, 0,
            "nothing may be refused at an HP site under the length model"
        );
        assert_eq!(
            stats.reads_reassigned_indel_to_indel, 1,
            "exactly one Ins(1) -> Del(1) move"
        );
        assert_eq!(
            stats.reassignments_by_hp_length.get(8).copied().unwrap_or(0),
            2,
            "both moves are booked against the 8bp run: rins and rref"
        );
        assert_eq!(stats.indel_edits_emitted_but_not_applied, 0);

        // The eight reads that already carry the kept length must be untouched.
        for n in NAMES[..8].iter() {
            let (seq, cigar) = out.get(&n.as_bytes().to_vec()).expect("read must survive");
            assert_eq!(seq.as_slice(), HP_DEL_SEQ, "{n} must be untouched");
            assert_eq!(*cigar, CigarString(hp_del_cigar()), "{n} must be untouched");
        }
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn gained_insertion_is_rejected_when_the_read_stops_short_of_the_run_plus_flank() {
        // The OLD guard asked only for `flank + max_len` past the site.
        let dir = std::env::temp_dir().join("himito_denoise_indel_short_of_tract");
        std::fs::create_dir_all(&dir).unwrap();
        let inb = dir.join("in.bam");
        let outb = dir.join("out.bam");

        let seq = &CLEAN_SEQ[..18];
        let cigar = vec![Cigar::Match(18)];
        write_cigar_bam(&inb, "chrM", 29, &[("rshort", 0, seq, cigar.clone(), false)]);

        let mut models = IndelModels::new();
        models.entry(b"chrM".to_vec()).or_default().push((7u32, hp_ins_site_model()));

        let corr: Corrections = HashMap::new();
        let mut stats = DenoiseStats::default();
        apply_corrections(
            &inb, &outb, &corr, &indel_refs(), &models, &indels_on(), &mut stats,
        )
        .unwrap();
        let (proc, modified) = (stats.reads_processed, stats.reads_modified);
        assert_eq!(proc, 1);
        assert_eq!(
            modified, 0,
            "the read must not be rewritten -- it does not span the run plus flank"
        );

        let mut r = Reader::from_path(&outb).unwrap();
        let rec = r.records().next().unwrap().unwrap();
        assert_eq!(rec.seq().as_bytes(), seq, "SEQ must be untouched");
        assert_eq!(rec.cigar().take(), CigarString(cigar), "CIGAR must be untouched");
        assert_eq!(
            stats.reads_reassigned_ref_to_ins, 0,
            "no reassignment must be counted for a rejected edit"
        );
        assert!(
            stats.assignments_skipped_span_guard >= 1,
            "the short-of-the-tract rejection must be counted"
        );
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn gained_insertion_is_applied_when_the_read_spans_the_run_plus_flank() {
        // `rspan` ends at reference 22 >= 16 + 5, so the SAME decided site still rewrites it.
        let dir = std::env::temp_dir().join("himito_denoise_indel_spans_tract");
        std::fs::create_dir_all(&dir).unwrap();
        let inb = dir.join("in.bam");
        let outb = dir.join("out.bam");

        let seq = &CLEAN_SEQ[..22];
        let cigar = vec![Cigar::Match(22)];
        write_cigar_bam(&inb, "chrM", 29, &[("rspan", 0, seq, cigar.clone(), false)]);

        let mut models = IndelModels::new();
        models.entry(b"chrM".to_vec()).or_default().push((7u32, hp_ins_site_model()));

        let corr: Corrections = HashMap::new();
        let mut stats = DenoiseStats::default();
        apply_corrections(
            &inb, &outb, &corr, &indel_refs(), &models, &indels_on(), &mut stats,
        )
        .unwrap();
        let (proc, modified) = (stats.reads_processed, stats.reads_modified);
        assert_eq!(proc, 1);
        assert_eq!(modified, 1, "a read spanning the run plus flank must still be corrected");
        assert_eq!(
            stats.reads_reassigned_ref_to_ins, 1,
            "the gained insertion must be booked once"
        );

        let mut r = Reader::from_path(&outb).unwrap();
        let rec = r.records().next().unwrap().unwrap();
        assert_eq!(
            rec.seq().as_bytes().len(),
            seq.len() + 1,
            "one base must be inserted"
        );
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn gained_deletion_needing_exactly_the_old_bound_is_now_rejected() {
        // the match-run bound for a gained Del(m) must be `m + flank + 1` -- the OLD bound `m + flank` was one short.
        let dir = std::env::temp_dir().join("himito_denoise_indel_offbyone_bound");
        std::fs::create_dir_all(&dir).unwrap();
        let inb = dir.join("in.bam");
        let outb = dir.join("out.bam");

        let seq = vec![b'A'; 35]; // 25 (M) + 10 (M) query bases; RefSkip consumes none
        let cigar = vec![Cigar::Match(25), Cigar::RefSkip(3), Cigar::Match(10)];
        write_cigar_bam(&inb, "chrM", 40, &[("rguard", 0, seq.as_slice(), cigar.clone(), false)]);

        let mut refs = HashMap::new();
        refs.insert(b"chrM".to_vec(), vec![b'A'; 40]);

        let model = indel::IndelSiteModel {
            kept: vec![(indel::Allele::Ref, 0.1), (indel::Allele::Del(2), 0.9)],
            eps: 0.15,
            context_len: 1,
            strand_biased: 0,
            obs_vaf: vec![(indel::Allele::Ref, 0.1), (indel::Allele::Del(2), 0.9)],
            strand_rejected: vec![],
            hp_run: None,
        };
        let mut models = IndelModels::new();
        models.entry(b"chrM".to_vec()).or_default().push((18u32, model));

        let corr: Corrections = HashMap::new();
        let mut stats = DenoiseStats::default();
        apply_corrections(
            &inb, &outb, &corr, &refs, &models, &indels_on(), &mut stats,
        )
        .unwrap();
        let (proc, modified) = (stats.reads_processed, stats.reads_modified);
        assert_eq!(proc, 1);
        assert_eq!(modified, 0, "the read must not be rewritten -- its real flank is one base short");

        let mut r = Reader::from_path(&outb).unwrap();
        let rec = r.records().next().unwrap().unwrap();
        assert_eq!(rec.seq().as_bytes(), seq.as_slice(), "SEQ must be untouched");
        assert_eq!(rec.cigar().take(), CigarString(cigar), "CIGAR must be untouched");
        assert!(
            stats.assignments_skipped_span_guard >= 1,
            "the one-short flank must be counted as a span-guard skip"
        );
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn shrinking_a_deletion_is_counted_as_bases_inserted_not_removed() {
        // `Del(3) -> Del(1)` splices 2 reference bases BACK INTO the read.
        let dir = std::env::temp_dir().join("himito_denoise_indel_shrink_counter");
        std::fs::create_dir_all(&dir).unwrap();
        let inb = dir.join("in.bam");
        let outb = dir.join("out.bam");

        // 10M 3D 10M: this read's OWN deletion is anchored at 9.
        let seq = vec![b'A'; 20];
        let cigar = vec![Cigar::Match(10), Cigar::Del(3), Cigar::Match(10)];
        write_cigar_bam(&inb, "chrM", 40, &[("rshrink", 0, seq.as_slice(), cigar.clone(), false)]);

        // NOT a uniform-base reference: a homopolymer run would left-normalize this read's own Del(3) all the way to position 0.
        let refseq: Vec<u8> = (0..40).map(|i| b"ACGT"[i % 4]).collect();
        let mut refs = HashMap::new();
        refs.insert(b"chrM".to_vec(), refseq);

        // Decided site at 9: Del(1) dominates 90/10 over Ref.
        let model = indel::IndelSiteModel {
            kept: vec![(indel::Allele::Ref, 0.1), (indel::Allele::Del(1), 0.9)],
            eps: 0.15,
            context_len: 1,
            strand_biased: 0,
            obs_vaf: vec![(indel::Allele::Ref, 0.1), (indel::Allele::Del(1), 0.9)],
            strand_rejected: vec![],
            hp_run: None,
        };
        let mut models = IndelModels::new();
        models.entry(b"chrM".to_vec()).or_default().push((9u32, model));

        let corr: Corrections = HashMap::new();
        let mut stats = DenoiseStats::default();
        apply_corrections(
            &inb, &outb, &corr, &refs, &models, &indels_on(), &mut stats,
        )
        .unwrap();
        let (proc, modified) = (stats.reads_processed, stats.reads_modified);
        assert_eq!(proc, 1);
        assert_eq!(modified, 1);

        let mut r = Reader::from_path(&outb).unwrap();
        let rec = r.records().next().unwrap().unwrap();
        assert_eq!(
            rec.cigar().take(),
            CigarString(vec![Cigar::Match(10), Cigar::Del(1), Cigar::Match(12)]),
            "the deletion must shrink to 1bp, splicing the other 2 ref bases back as matches"
        );
        assert_eq!(rec.seq().len(), 22, "the 2 restored reference bases must lengthen SEQ");
        assert_eq!(stats.reads_reassigned_indel_to_indel, 1);
        assert_eq!(
            stats.indel_bases_inserted, 2,
            "shrinking Del(3) -> Del(1) restores 2 bases to the read: bases INSERTED, not removed"
        );
        assert_eq!(
            stats.indel_bases_removed, 0,
            "the shrink must not be double-booked as bases removed"
        );
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn gained_deletion_swallows_a_decided_site_inside_its_own_span() {
        // A gained Del(3) at site 5 consumes reference positions 6..=8 once applied.
        let dir = std::env::temp_dir().join("himito_denoise_indel_swallowed_span");
        std::fs::create_dir_all(&dir).unwrap();
        let inb = dir.join("in.bam");
        let outb = dir.join("out.bam");

        let seq = vec![b'A'; 30];
        let cigar = vec![Cigar::Match(30)];
        write_cigar_bam(&inb, "chrM", 40, &[("rboth", 0, seq.as_slice(), cigar.clone(), false)]);

        // The tract span guard refuses to rewrite a read that cannot see the end of the run.
        let refseq: Vec<u8> = (0..40).map(|i| b"ACGT"[i % 4]).collect();
        let mut refs = HashMap::new();
        refs.insert(b"chrM".to_vec(), refseq);

        // Site 5: Ref -> Del(3) dominates 90/10.
        let del_model = indel::IndelSiteModel {
            kept: vec![(indel::Allele::Ref, 0.1), (indel::Allele::Del(3), 0.9)],
            eps: 0.15,
            context_len: 1,
            strand_biased: 0,
            obs_vaf: vec![(indel::Allele::Ref, 0.1), (indel::Allele::Del(3), 0.9)],
            strand_rejected: vec![],
            hp_run: None,
        };
        // Site 7 (inside site 5's Del(3) span, 6..=8): Ref -> Ins("A") dominates 90/10.
        let ins_model = indel::IndelSiteModel {
            kept: vec![(indel::Allele::Ref, 0.1), (indel::Allele::Ins(b"A".to_vec()), 0.9)],
            eps: 0.15,
            context_len: 1,
            strand_biased: 0,
            obs_vaf: vec![(indel::Allele::Ref, 0.1), (indel::Allele::Ins(b"A".to_vec()), 0.9)],
            strand_rejected: vec![],
            hp_run: None,
        };
        let mut models = IndelModels::new();
        models
            .entry(b"chrM".to_vec())
            .or_default()
            .extend([(5u32, del_model), (7u32, ins_model)]);

        let corr: Corrections = HashMap::new();
        let mut stats = DenoiseStats::default();
        apply_corrections(
            &inb, &outb, &corr, &refs, &models, &indels_on(), &mut stats,
        )
        .unwrap();
        let (proc, modified) = (stats.reads_processed, stats.reads_modified);
        assert_eq!(proc, 1);
        assert_eq!(modified, 1, "the site-5 deletion must still be applied");

        let mut r = Reader::from_path(&outb).unwrap();
        let rec = r.records().next().unwrap().unwrap();
        assert_eq!(
            rec.cigar().take(),
            // The Match arm emits the anchor base in the preceding Match run before splicing in the deletion.
            CigarString(vec![Cigar::Match(6), Cigar::Del(3), Cigar::Match(21)]),
            "only the site-5 deletion may appear; no insertion from site 7"
        );
        assert_eq!(rec.seq().len(), 27, "only 3 ref bases lost, matching a single Del(3)");
        assert_eq!(stats.reads_reassigned_ref_to_del, 1);
        assert_eq!(stats.indel_bases_removed, 3);
        assert_eq!(
            stats.reads_reassigned_ref_to_ins, 0,
            "site 7 must never be evaluated for this read, let alone counted as a reassignment"
        );
        assert_eq!(stats.indel_bases_inserted, 0, "no phantom insertion bases may be counted");
        assert!(
            stats.assignments_skipped_swallowed_by_gained_deletion >= 1,
            "the swallowed site must be counted"
        );
        std::fs::remove_dir_all(&dir).ok();
    }

    // ---- a read's own D/N covering the normalized column ----

    #[test]
    fn read_deleting_the_normalized_column_does_not_cast_a_phantom_alt_vote() {
        // Depth is only incremented where a read has qpos().is_some() -- never inside that SAME read's own D/N.
        let dir = std::env::temp_dir().join("himito_denoise_indel_deletes_own_norm_col");
        std::fs::create_dir_all(&dir).unwrap();
        let bam = dir.join("in.bam");
        let mut reads: Vec<(&str, i64, &[u8], Vec<Cigar>, bool)> = vec![];
        for (i, n) in NAMES.iter().enumerate() {
            reads.push((n, 0, CLEAN_SEQ, clean_cigar(), i % 2 == 0));
        }
        // ref[0..7) + (skip ref[7]) + ref[8..16) + inserted "A" + ref[16..29).
        let mut del_seq: Vec<u8> = Vec::new();
        del_seq.extend_from_slice(&REF29[0..7]);
        del_seq.extend_from_slice(&REF29[8..16]);
        del_seq.push(b'A');
        del_seq.extend_from_slice(&REF29[16..29]);
        let del_cigar = vec![
            Cigar::Match(7), Cigar::Del(1), Cigar::Match(8), Cigar::Ins(1), Cigar::Match(13),
        ];
        reads.push(("rdel", 0, del_seq.as_slice(), del_cigar, false));
        write_cigar_bam(&bam, "chrM", 29, &reads);

        let (_corr, models, stats) =
            compute_corrections(&bam, &indel_refs(), 2, 0.01, SB_P, HOM_VAF, &indels_on()).unwrap();

        assert!(
            stats.indel_events_out_of_span >= 1,
            "the phantom vote at column 7 must be rejected and counted"
        );
        let no_site_7 = models
            .get(&b"chrM".to_vec())
            .map_or(true, |list| list.iter().all(|(p, _)| *p != 7));
        assert!(no_site_7, "site 7 must not be fabricated from a read with no depth there");
        std::fs::remove_dir_all(&dir).ok();
    }

    // ---- a gained deletion must fit inside ONE match run ----

    #[test]
    fn gained_deletion_is_skipped_when_its_own_match_run_is_too_short() {
        // The span/flank guard bounds only the read's OVERALL alignment span.
        let dir = std::env::temp_dir().join("himito_denoise_indel_shortrun");
        std::fs::create_dir_all(&dir).unwrap();
        let inb = dir.join("in.bam");
        let outb = dir.join("out.bam");

        let seq = vec![b'A'; 30];
        let cigar = vec![Cigar::Match(20), Cigar::RefSkip(3), Cigar::Match(10)];
        write_cigar_bam(&inb, "chrM", 40, &[("rguard", 0, seq.as_slice(), cigar, false)]);

        let mut refs = HashMap::new();
        refs.insert(b"chrM".to_vec(), vec![b'A'; 40]);

        // Hand-built decided site at position 18: Del(2) dominates 90/10 over Ref.
        let model = indel::IndelSiteModel {
            kept: vec![(indel::Allele::Ref, 0.1), (indel::Allele::Del(2), 0.9)],
            eps: 0.15,
            context_len: 1,
            strand_biased: 0,
            obs_vaf: vec![(indel::Allele::Ref, 0.1), (indel::Allele::Del(2), 0.9)],
            strand_rejected: vec![],
            hp_run: None,
        };
        let mut models = IndelModels::new();
        models.entry(b"chrM".to_vec()).or_default().push((18u32, model));

        let corr: Corrections = HashMap::new();
        let mut stats = DenoiseStats::default();
        apply_corrections(
            &inb, &outb, &corr, &refs, &models, &indels_on(), &mut stats,
        )
        .unwrap();
        let (proc, modified) = (stats.reads_processed, stats.reads_modified);
        assert_eq!(proc, 1);
        assert_eq!(modified, 0, "the guarded read must not be rewritten");

        let mut r = Reader::from_path(&outb).unwrap();
        let rec = r.records().next().unwrap().unwrap();
        assert_eq!(rec.seq().as_bytes(), seq.as_slice(), "SEQ must be untouched");
        assert_eq!(
            rec.cigar().take(),
            CigarString(vec![Cigar::Match(20), Cigar::RefSkip(3), Cigar::Match(10)]),
            "CIGAR must be untouched"
        );
        assert!(
            stats.assignments_skipped_span_guard >= 1,
            "the too-short match run must be counted as a span-guard skip"
        );
        std::fs::remove_dir_all(&dir).ok();
    }

    // ---- revert-direction base counts ----

    #[test]
    fn lone_spurious_deletion_is_reverted_and_counted_as_bases_inserted() {
        // `Del(n)->Ref` restores n bases to the read but was not being counted in `indel_bases_inserted`.
        let dir = std::env::temp_dir().join("himito_denoise_indel_lone_del");
        std::fs::create_dir_all(&dir).unwrap();
        let bam = dir.join("in.bam");
        let outb = dir.join("out.bam");
        let mut reads: Vec<(&str, i64, &[u8], Vec<Cigar>, bool)> = vec![];
        for (i, n) in NAMES.iter().enumerate() {
            reads.push((n, 0, CLEAN_SEQ, clean_cigar(), i % 2 == 0));
        }
        // ref[0..7) + ref[8..29): a 1bp deletion anchored at 6.
        let mut del_seq: Vec<u8> = Vec::new();
        del_seq.extend_from_slice(&REF29[0..7]);
        del_seq.extend_from_slice(&REF29[8..29]);
        let del_cigar = vec![Cigar::Match(7), Cigar::Del(1), Cigar::Match(21)];
        reads.push(("rdel", 0, del_seq.as_slice(), del_cigar, false));
        write_cigar_bam(&bam, "chrM", 29, &reads);

        let refs = indel_refs();
        let iopts = indels_on();
        let (corr, models, mut stats) =
            compute_corrections(&bam, &refs, 2, 0.01, SB_P, HOM_VAF, &iopts).unwrap();
        apply_corrections(&bam, &outb, &corr, &refs, &models, &iopts, &mut stats).unwrap();

        let mut r = Reader::from_path(&outb).unwrap();
        let mut seen = false;
        for rec in r.records() {
            let rec = rec.unwrap();
            if rec.qname() == b"rdel" {
                assert_eq!(rec.seq().as_bytes(), CLEAN_SEQ, "the spurious deletion must be reverted");
                assert_eq!(rec.cigar().to_string(), "29M");
                seen = true;
            } else {
                assert_eq!(rec.seq().as_bytes(), CLEAN_SEQ, "clean reads must be untouched");
            }
        }
        assert!(seen);
        assert_eq!(stats.reads_reassigned_indel_to_ref, 1);
        assert_eq!(stats.indel_bases_inserted, 1, "restoring the deleted base must count as inserted");
        std::fs::remove_dir_all(&dir).ok();
    }

    // ---- IMPORTANT 4 regression: tid2name must not be called for tid == -1 ----

    #[test]
    fn apply_corrections_treats_mapped_record_with_no_contig_like_missing_ref() {
        // A mapped record with tid == -1 is a degenerate BAM state.
        let dir = std::env::temp_dir().join("himito_denoise_tidless");
        std::fs::create_dir_all(&dir).unwrap();
        let inb = dir.join("in.bam");
        let outb = dir.join("out.bam");
        {
            let mut header = Header::new();
            let mut sq = HeaderRecord::new(b"SQ");
            sq.push_tag(b"SN", "chrM");
            sq.push_tag(b"LN", &6usize.to_string());
            header.push_record(&sq);
            let mut w = Writer::from_path(&inb, &header, Format::Bam).unwrap();
            let mut rec = Record::new();
            let seq = b"AAAAAA";
            let quals = vec![30u8; 6];
            rec.set(b"rtidless", Some(&CigarString(vec![Cigar::Match(6)])), seq, &quals);
            // Deliberately do NOT set_tid: tid stays -1 (Record::new()'s default).
            rec.unset_unmapped();
            w.write(&rec).unwrap();
        }

        let refs = indel_refs();
        let models = IndelModels::new();
        let mut stats = DenoiseStats::default();
        let corr: Corrections = HashMap::new();
        apply_corrections(
            &inb, &outb, &corr, &refs, &models, &indel::IndelOpts::default(), &mut stats,
        )
        .unwrap();
        let (proc, modified) = (stats.reads_processed, stats.reads_modified); // must not panic on an unchecked tid2name offset
        assert_eq!(proc, 1);
        assert_eq!(modified, 0);
        std::fs::remove_dir_all(&dir).ok();
    }

    // ---- drive statistics from what was

    #[test]
    fn edit_emitted_for_an_unreachable_ref_pos_is_a_harmless_missed_correction() {
        // Rather than trust that every future instance will be caught by name.
        let dir = std::env::temp_dir().join("himito_denoise_indel_unreachable_edit");
        std::fs::create_dir_all(&dir).unwrap();
        let inb = dir.join("in.bam");
        let outb = dir.join("out.bam");

        let seq = vec![b'A'; 20];
        let cigar = vec![Cigar::Match(10), Cigar::Del(3), Cigar::Match(10)];
        write_cigar_bam(&inb, "chrM", 30, &[("rvictim", 0, seq.as_slice(), cigar.clone(), false)]);

        let mut refs = HashMap::new();
        refs.insert(b"chrM".to_vec(), vec![b'A'; 30]);

        let mut corr: Corrections = HashMap::new();
        corr.insert(
            b"rvictim".to_vec(),
            indel::ReadEdits {
                subs: vec![],
                indels: vec![indel::IndelEdit {
                    ref_pos: 11,
                    from: indel::Allele::Ref,
                    to: indel::Allele::Ins(b"A".to_vec()),
                }],
            },
        );

        let models = IndelModels::new();
        let mut stats = DenoiseStats::default();
        apply_corrections(
            &inb, &outb, &corr, &refs, &models, &indels_on(), &mut stats,
        )
        .unwrap();
        let (proc, modified) = (stats.reads_processed, stats.reads_modified);
        assert_eq!(proc, 1);
        assert_eq!(modified, 0, "the walk never reaches ref_pos 11, so nothing may change");

        let mut r = Reader::from_path(&outb).unwrap();
        let rec = r.records().next().unwrap().unwrap();
        assert_eq!(rec.seq().as_bytes(), seq.as_slice(), "SEQ must be untouched");
        assert_eq!(rec.cigar().take(), CigarString(cigar), "CIGAR must be untouched");

        assert_eq!(stats.reads_reassigned_ref_to_ins, 0);
        assert_eq!(stats.reads_reassigned_ref_to_del, 0);
        assert_eq!(stats.reads_reassigned_indel_to_ref, 0);
        assert_eq!(stats.reads_reassigned_indel_to_indel, 0);
        assert_eq!(stats.indel_bases_inserted, 0);
        assert_eq!(stats.indel_bases_removed, 0);
        // The edit is hand-fed to `apply_corrections`, so the by-length vector stays empty.
        assert_eq!(
            stats.reassignments_by_hp_length.len(),
            0,
            "no hp-length bucket may exist -- none may have been touched"
        );
        // also confirm this edit bypassed EVERY pass-2 guard.
        assert_eq!(stats.assignments_skipped_span_guard, 0);
        assert_eq!(stats.assignments_skipped_unsupported_transition, 0);
        assert_eq!(stats.assignments_skipped_oversized_neighbor, 0);
        assert_eq!(stats.assignments_skipped_ambiguous_carry, 0);
        assert_eq!(stats.assignments_skipped_swallowed_by_gained_deletion, 0);
        assert_eq!(
            stats.indel_edits_emitted_but_not_applied, 1,
            "the emitted-but-unreachable edit must be loudly visible in the diagnostic counter"
        );
        std::fs::remove_dir_all(&dir).ok();
    }

    // ---- duplicate emitted edits must not be

    #[test]
    fn two_emitted_edits_sharing_ref_pos_and_kind_are_not_double_counted() {
        let dir = std::env::temp_dir().join("himito_denoise_indel_duplicate_emitted_edit");
        std::fs::create_dir_all(&dir).unwrap();
        let inb = dir.join("in.bam");
        let outb = dir.join("out.bam");

        let mut seq = vec![b'A'; 10];
        seq.push(b'T'); // first inserted base
        seq.push(b'C'); // second inserted base
        seq.extend(vec![b'A'; 10]);
        let cigar = vec![
            Cigar::Match(10),
            Cigar::Ins(1),
            Cigar::Del(2),
            Cigar::Ins(1),
            Cigar::Match(10),
        ];
        write_cigar_bam(&inb, "chrM", 30, &[("rdup", 0, seq.as_slice(), cigar.clone(), false)]);

        let mut refs = HashMap::new();
        refs.insert(b"chrM".to_vec(), vec![b'A'; 30]);

        let mut corr: Corrections = HashMap::new();
        corr.insert(
            b"rdup".to_vec(),
            indel::ReadEdits {
                subs: vec![],
                indels: vec![
                    indel::IndelEdit {
                        ref_pos: 9,
                        from: indel::Allele::Ins(b"T".to_vec()),
                        to: indel::Allele::Ref,
                    },
                    indel::IndelEdit {
                        ref_pos: 9,
                        from: indel::Allele::Ins(b"C".to_vec()),
                        to: indel::Allele::Ref,
                    },
                ],
            },
        );

        let models = IndelModels::new();
        let mut stats = DenoiseStats::default();
        apply_corrections(
            &inb, &outb, &corr, &refs, &models, &indels_on(), &mut stats,
        )
        .unwrap();
        let (proc, modified) = (stats.reads_processed, stats.reads_modified);
        assert_eq!(proc, 1);
        assert_eq!(modified, 1, "one of the two Ins ops must actually be reverted");

        // Exactly one of the two duplicate-keyed edits is booked as applied.
        assert_eq!(
            stats.reads_reassigned_indel_to_ref, 1,
            "only ONE of the two duplicate-keyed edits may be booked as applied"
        );
        assert_eq!(
            stats.indel_edits_emitted_but_not_applied, 1,
            "the other duplicate-keyed edit must be reported as not applied"
        );
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn start_with_indels_disabled_is_deterministic_across_runs() {
        // NOTE: this only proves DETERMINISM -- two independent runs of the SAME code produce identical bytes.
        let dir = std::env::temp_dir().join("himito_denoise_optin");
        std::fs::create_dir_all(&dir).unwrap();
        let inb = dir.join("in.bam");
        let out_a = dir.join("a.bam");
        let out_b = dir.join("b.bam");
        let reff = dir.join("ref.fa");
        std::fs::write(&reff, ">chrM\nACGTACGTAAAAAAAACGTACGTACGTAC\n").unwrap();

        let mut reads: Vec<(&str, i64, &[u8], Vec<Cigar>, bool)> = vec![];
        for (i, n) in NAMES.iter().enumerate() {
            reads.push((n, 0, CLEAN_SEQ, clean_cigar(), i % 2 == 0));
        }
        reads.push(("rins", 0, UNIQ_INS_SEQ, uniq_ins_cigar(), false));
        write_cigar_bam(&inb, "chrM", 29, &reads);

        let off = indel::IndelOpts::default();
        start(&inb, &out_a, &reff, "ont-r10", 0.01, HOM_VAF, &off, None).unwrap();
        start(&inb, &out_b, &reff, "ont-r10", 0.01, HOM_VAF, &off, None).unwrap();
        assert_eq!(std::fs::read(&out_a).unwrap(), std::fs::read(&out_b).unwrap());

        // And the indel-bearing read keeps its insertion.
        let mut r = Reader::from_path(&out_a).unwrap();
        let mut seen = false;
        for rec in r.records() {
            let rec = rec.unwrap();
            if rec.qname() == b"rins" {
                assert_eq!(rec.seq().as_bytes(), UNIQ_INS_SEQ);
                seen = true;
            }
        }
        assert!(seen);
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn start_with_indels_disabled_preserves_v1_cigar_on_substitution_reads() {
        // a/2b: the opt-in guarantee that actually matters is not "two runs give the same bytes" but the real v1 property.
        let dir = std::env::temp_dir().join("himito_denoise_v1_cigar_identity");
        std::fs::create_dir_all(&dir).unwrap();
        let inb = dir.join("in.bam");
        let outb = dir.join("out.bam");
        let reff = dir.join("ref.fa");
        std::fs::write(&reff, ">chrM\nAAAAAA\n").unwrap();

        // "AAGAAA" against ref "AAAAAA": one substitution at ref pos 2.
        let err_cigar = vec![Cigar::Equal(2), Cigar::Diff(1), Cigar::Equal(3)];
        let clean: &[u8] = b"AAAAAA";
        let clean_names: [&str; 9] = ["r0","r1","r2","r3","r4","r5","r6","r7","r8"];
        let mut reads: Vec<(&str, i64, &[u8], Vec<Cigar>, bool)> = vec![];
        for (i, n) in clean_names.iter().enumerate() {
            reads.push((n, 0, clean, vec![Cigar::Match(6)], i % 2 == 0));
        }
        reads.push(("rerr", 0, b"AAGAAA", err_cigar.clone(), false));
        write_cigar_bam(&inb, "chrM", 6, &reads);

        start(&inb, &outb, &reff, "ont-r10", 0.01, HOM_VAF, &indel::IndelOpts::default(), None).unwrap();

        let mut r = Reader::from_path(&outb).unwrap();
        let mut seen = false;
        for rec in r.records() {
            let rec = rec.unwrap();
            if rec.qname() == b"rerr" {
                assert_eq!(rec.seq().as_bytes(), b"AAAAAA", "the substitution must still be corrected");
                assert_eq!(
                    rec.cigar().take(),
                    CigarString(err_cigar.clone()),
                    "CIGAR must be byte-identical to the input, `=`/`X` ops included"
                );
                seen = true;
            }
        }
        assert!(seen);
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn start_with_indels_enabled_removes_the_spurious_insertion() {
        let dir = std::env::temp_dir().join("himito_denoise_optin_on");
        std::fs::create_dir_all(&dir).unwrap();
        let inb = dir.join("in.bam");
        let outb = dir.join("out.bam");
        let reff = dir.join("ref.fa");
        let statsf = dir.join("stats.json");
        std::fs::write(&reff, ">chrM\nACGTACGTAAAAAAAACGTACGTACGTAC\n").unwrap();

        let mut reads: Vec<(&str, i64, &[u8], Vec<Cigar>, bool)> = vec![];
        for (i, n) in NAMES.iter().enumerate() {
            reads.push((n, 0, CLEAN_SEQ, clean_cigar(), i % 2 == 0));
        }
        reads.push(("rins", 0, UNIQ_INS_SEQ, uniq_ins_cigar(), false));
        write_cigar_bam(&inb, "chrM", 29, &reads);

        let mut on = indel::IndelOpts::default();
        on.enabled = true;
        start(&inb, &outb, &reff, "ont-r10", 0.01, HOM_VAF, &on, Some(&statsf)).unwrap();

        let mut r = Reader::from_path(&outb).unwrap();
        let mut seen = false;
        for rec in r.records() {
            let rec = rec.unwrap();
            if rec.qname() == b"rins" {
                assert_eq!(rec.seq().as_bytes(), CLEAN_SEQ, "insertion must be removed");
                assert_eq!(rec.cigar().to_string(), "29M", "CIGAR must be rebuilt");
                seen = true;
            }
        }
        assert!(seen);
        let json = std::fs::read_to_string(&statsf).unwrap();
        assert!(json.contains("\"reads_reassigned_indel_to_ref\": 1"));
        std::fs::remove_dir_all(&dir).ok();
    }

    // ---- a genuine heteroplasmic indel must

    #[test]
    fn genuine_heteroplasmic_indel_survives_end_to_end_via_candidacy_not_the_protect_vaf_backstop() {
        // End-to-end check that heteroplasmy is preserved; other tests use hand-built `IndelSiteModel`s.
        let dir = std::env::temp_dir().join("himito_denoise_e2e_genuine_heteroplasmy");
        std::fs::create_dir_all(&dir).unwrap();
        let bam = dir.join("in.bam");
        let outb = dir.join("out.bam");

        let names20: Vec<String> = (0..20).map(|i| format!("r{i:02}")).collect();
        let mut reads: Vec<(&str, i64, &[u8], Vec<Cigar>, bool)> = vec![];
        // 12 reads carrying the insertion plus 8 clean reads -- both groups individually strand-balanced.
        for i in 0..12 {
            reads.push((&names20[i], 0, HP_INS_SEQ, hp_ins_cigar(), i % 2 == 0));
        }
        for i in 12..20 {
            reads.push((&names20[i], 0, CLEAN_SEQ, clean_cigar(), i % 2 == 0));
        }
        write_cigar_bam(&bam, "chrM", 29, &reads);

        let refs = indel_refs();
        let iopts = indels_on();
        let (corr, models, mut stats) =
            compute_corrections(&bam, &refs, 2, 0.01, SB_P, HOM_VAF, &iopts).unwrap();

        let list = models.get(&b"chrM".to_vec()).expect("chrM must have a site at 7");
        let (_, model) = list.iter().find(|(p, _)| *p == 7).expect("site at 7");
        assert_eq!(model.context_len, 8, "test premise: site 7 is an 8bp homopolymer");
        let ins = indel::Allele::Ins(b"A".to_vec());
        let ins_vaf = model.obs_vaf.iter().find(|(a, _)| *a == ins).map(|(_, v)| *v);
        assert!(
            (ins_vaf.expect("insertion must be observed") - 0.6).abs() < 1e-9,
            "test premise: the insertion's raw VAF must be exactly 60%, got {ins_vaf:?}"
        );
        assert!(
            model.kept.iter().any(|(a, _)| *a == ins),
            "60% must genuinely clear candidacy at context length 8 -- this must not be \
             a protect_vaf rescue"
        );
        assert!(
            model.kept.iter().any(|(a, _)| *a == indel::Allele::Ref),
            "REF is always kept"
        );

        apply_corrections(&bam, &outb, &corr, &refs, &models, &iopts, &mut stats).unwrap();

        let mut r = Reader::from_path(&outb).unwrap();
        let mut seen = 0;
        for (i, rec) in r.records().enumerate() {
            let rec = rec.unwrap();
            let (expected_seq, expected_cigar): (&[u8], Vec<Cigar>) = if i < 12 {
                (HP_INS_SEQ, hp_ins_cigar())
            } else {
                (CLEAN_SEQ, clean_cigar())
            };
            assert_eq!(rec.seq().as_bytes(), expected_seq, "read {i} SEQ must be unchanged");
            assert_eq!(
                rec.cigar().take(),
                CigarString(expected_cigar),
                "read {i} CIGAR must be unchanged"
            );
            seen += 1;
        }
        assert_eq!(seen, 20);

        assert_eq!(stats.reads_reassigned_ref_to_ins, 0);
        assert_eq!(stats.reads_reassigned_ref_to_del, 0);
        assert_eq!(stats.reads_reassigned_indel_to_ref, 0);
        assert_eq!(stats.reads_reassigned_indel_to_indel, 0);
        std::fs::remove_dir_all(&dir).ok();
    }
}
