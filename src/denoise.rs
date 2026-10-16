use std::collections::{BTreeSet, HashMap, HashSet};
use std::path::{Path, PathBuf};
use anyhow::{Context, Result};
use bio::io::fasta;
use log::{info, warn};
use rust_htslib::bam::{self, Read, Reader};
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::record::CigarString;

use crate::denoise_indel as indel;

const NQ: usize = 94; // Phred qualities 0..=93

/// Ceiling on the per-base error probability.
const E_MAX: f64 = 0.75;

/// Minimum reads per strand for a non-reference allele to be eligible.
pub(crate) const MIN_STRAND: u32 = 2;
/// Strand-bias test p-value threshold.
pub(crate) const STRAND_BIAS_P: f64 = 0.01;
/// Alleles at or above this within-column frequency skip the strand-bias test (a single-strand artifact cannot
/// get there) and are not minor alleles for linkage rescue. Matches `call`'s near-homoplasmic exemption (0.7).
pub(crate) const HOMOPLASMIC_VAF: f64 = 0.7;
/// Non-kept alleles at or above this frequency are pending: corrected unless rescued by linkage.
const PENDING_FLOOR: f64 = 0.005;

#[inline]
fn phred_to_e(q: u8) -> f64 {
    let q = q.min((NQ - 1) as u8);
    10f64.powf(-(q as f64) / 10.0).min(E_MAX)
}

#[inline]
pub(crate) fn base_to_idx(b: u8) -> Option<usize> {
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

struct ColumnCounts {
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

struct SiteModel {
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

/// EM fit of one column's allele frequencies under a per-quality error model, then the keep/pending rule.
fn fit_site(c: &ColumnCounts, vaf: f64) -> SiteModel {
    // The null for the strand-bias test is this column's OWN forward fraction, not a library-wide constant.
    let (col_fwd, col_rev) = c.strand_totals();
    let col_total = col_fwd + col_rev;
    let expected_fwd_frac = if col_total == 0 {
        0.5
    } else {
        col_fwd as f64 / col_total as f64
    };

    // The reference is always eligible.
    let mut eligible = [false; 4];
    let mut strand_biased_alts = 0usize;
    for a in 0..4 {
        if Some(a) == c.ref_idx {
            eligible[a] = true;
            continue;
        }
        let (fwd, rev) = (c.strand[a][0], c.strand[a][1]);
        if fwd < MIN_STRAND || rev < MIN_STRAND {
            continue;
        }
        // Mirrors call::filter_strand_bias's near-homoplasmic exemption.
        let prop = (fwd + rev) as f64 / col_total as f64;
        if prop < HOMOPLASMIC_VAF
            && crate::call::is_strand_skewed(fwd as usize, rev as usize, expected_fwd_frac, STRAND_BIAS_P)
        {
            strand_biased_alts += 1;
            continue;
        }
        eligible[a] = true;
    }

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
        } else if eligible[a] && f[a] >= PENDING_FLOOR {
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

/// MAP allele for one observation of `observed` at quality `q`.
fn map_allele(observed: usize, q: u8, m: &SiteModel) -> usize {
    // An observation of an allele the site model already kept is evidence FOR that allele, so it stands as observed.
    if m.kept[observed] {
        return observed;
    }

    // `observed` is not kept, so every candidate's likelihood is e/3: the most frequent kept allele wins.
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

/// Per-read edits from pass 1, by qname: substitutions only (indel edits are planned in pass 2).
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
    /// Votes dropped in pass 1: the read ends inside the homopolymer run at the left-normalized site
    /// (longer-unit tandem repeats are not gated).
    pub indel_events_not_spanning_repeat: u64,
    /// Votes dropped in pass 1: the same read already cast an ALT vote at this normalized site.
    pub indel_events_duplicate_site: u64,
    /// Sites whose only ALT votes were dropped as not spanning the run, and with no spanning depth.
    pub indel_sites_dropped_no_spanning_depth: u64,
    pub reads_reassigned_ref_to_ins: u64,
    pub reads_reassigned_ref_to_del: u64,
    pub reads_reassigned_indel_to_ref: u64,
    pub reads_reassigned_indel_to_indel: u64,
    pub indel_bases_inserted: u64,
    pub indel_bases_removed: u64,
    /// Index = repeat-context length L; audits how much correction fired in long runs.
    pub reassignments_by_hp_length: Vec<u64>,
    /// Substitutions dropped because they fall inside a rewritten homopolymer run.
    pub subs_dropped_in_hp_rewrite: u64,
    /// (read, site) pairs skipped: the read lacks the match-context flank the rewrite walk needs.
    pub assignments_skipped_span_guard: u64,
    /// (read, site) pairs skipped: the MAP target needs an unsupported transition.
    pub assignments_skipped_unsupported_transition: u64,
    /// (read, site) pairs skipped: an untouchable oversized indel of the read lies too close to the site.
    pub assignments_skipped_oversized_neighbor: u64,
    /// (read, site) pairs skipped: two CIGAR events of the read normalize to this site (ambiguous).
    pub assignments_skipped_ambiguous_carry: u64,
    /// (read, site) pairs skipped: an earlier gained deletion of the same read consumed the site.
    pub assignments_skipped_swallowed_by_gained_deletion: u64,
    /// Planned edits the rewrite did not apply (diagnostic; normally zero).
    pub indel_edits_emitted_but_not_applied: u64,
}

/// (tid, 0-based reference position, allele index) of a substitution allele.
type AlleleKey = (u32, u32, usize);

/// Pass 1: pileup the BAM, fit each column's substitution model and collect indel votes.
fn compute_corrections(
    bam_path: &Path,
    refs: &HashMap<Vec<u8>, Vec<u8>>,
    vaf: f64,
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
    // Remembered here because `sites` only receives votes that cleared the span gate.
    let mut span_gated_sites: HashSet<(Vec<u8>, u32)> = HashSet::new();

    // Each pending allele's would-be corrections (qname, MAP allele), held until the linkage rescue decides.
    let mut deferred: HashMap<AlleleKey, Vec<(Vec<u8>, usize)>> = HashMap::new();
    // Read indices carrying each pending or kept minor allele.
    let mut carriers: HashMap<AlleleKey, HashSet<usize>> = HashMap::new();
    let mut kept_minors: HashSet<AlleleKey> = HashSet::new();
    let mut read_idx: HashMap<Vec<u8>, usize> = HashMap::new();
    // Per read index: (tid, first, last) reference position with an observed base.
    let mut spans: Vec<(u32, u32, u32)> = Vec::new();
    let mut corrected_cols: HashSet<(u32, u32)> = HashSet::new();

    let mut plp = reader.pileup();
    plp.set_max_depth(1_000_000);
    for p in plp {
        let col = p.context("pileup error")?;
        let tid = col.tid();
        let name = header.tid2name(tid);
        let Some(refseq) = refs.get(name) else { continue };
        let pos = col.pos();
        let ref_idx = refseq.get(pos as usize).and_then(|&b| base_to_idx(b));

        let mut counts = ColumnCounts::new(ref_idx);
        // (qname, allele, qual, read idx) of each usable observation, corrected after the fit.
        let mut obs: Vec<(Vec<u8>, usize, u8, usize)> = Vec::new();

        // Indel collection MUST happen here, before the `non_ref() == 0` skip below.
        let mut col_depth = 0u32;
        let mut col_fwd = 0u32;
        let mut col_rev = 0u32;
        // (qname, allele, reverse, read's own POS, read's own CIGAR, read's exclusive reference end).
        let mut col_events: Vec<(Vec<u8>, indel::Allele, bool, i64, CigarString, i64)> =
            Vec::new();
        // A read that ends inside the run cannot show the run's length.
        let tract_need_end = if iopts.enabled {
            indel::hp_tract_end(refseq, pos) as i64 + 1
        } else {
            0
        };
        // At an HP site the column casts one LENGTH vote per spanning read.
        let hp_votes_here = if iopts.enabled && iopts.hp_length_model {
            indel::hp_site(refseq, pos)
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
            // `None` for a deletion or ref-skip: no query base here.
            let Some(qpos) = a.qpos() else { continue };
            let reverse = rec.is_reverse();

            if iopts.enabled {
                // This decides only whether the read joins THIS column's depth.
                let reference_end = rec.reference_end();
                if reference_end >= tract_need_end {
                    col_depth += 1;
                    if reverse { col_rev += 1 } else { col_fwd += 1 }

                    // `hp_allele`'s span gate is one flank STRICTER than `tract_need_end`.
                    if let Some((hp_e, b)) = hp_votes_here {
                        let cigar = rec.cigar().take();
                        let seq = rec.seq().as_bytes();
                        match indel::hp_allele(rec.pos(), &cigar, &seq, refseq, pos, hp_e) {
                            Some(0) => hp_col_ref.push(rec.qname().to_vec()),
                            Some(k) => hp_col_votes.push((
                                rec.qname().to_vec(),
                                indel::hp_length_allele(k, b),
                                reverse,
                            )),
                            None => {}
                        }
                    }
                }
                // An indel of MAX_INDEL_LEN or longer contributes depth only and is never corrected.
                let ev_allele = match a.indel() {
                    bam::pileup::Indel::Ins(n) if n < indel::MAX_INDEL_LEN => {
                        let s: Vec<u8> =
                            (0..n as usize).map(|i| rec.seq()[qpos + 1 + i]).collect();
                        Some(indel::Allele::Ins(s))
                    }
                    bam::pileup::Indel::Del(n) if n < indel::MAX_INDEL_LEN => Some(indel::Allele::Del(n)),
                    _ => None,
                };
                if let Some(allele) = ev_allele {
                    stats.indel_events_examined += 1;
                    col_events.push((
                        rec.qname().to_vec(),
                        allele,
                        reverse,
                        rec.pos(),
                        rec.cigar().take(),
                        reference_end,
                    ));
                }
            }

            let Some(allele) = base_to_idx(rec.seq()[qpos]) else { continue };
            let qual = rec.qual()[qpos];
            counts.add(allele, qual, reverse);
            let ri = match read_idx.get(rec.qname()) {
                Some(&i) => {
                    let s = &mut spans[i];
                    if s.0 == tid {
                        s.1 = s.1.min(pos);
                        s.2 = s.2.max(pos);
                    }
                    i
                }
                None => {
                    read_idx.insert(rec.qname().to_vec(), spans.len());
                    spans.push((tid, pos, pos));
                    spans.len() - 1
                }
            };
            obs.push((rec.qname().to_vec(), allele, qual, ri));
        }

        if iopts.enabled {
            depths.insert((name.to_vec(), pos), (col_depth, col_fwd, col_rev));
            // `add_from_read` still enforces one ALT vote per read per site.
            if !hp_col_votes.is_empty() {
                let site = sites.entry((name.to_vec(), pos)).or_default();
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
                let (norm_pos, norm_allele) = indel::normalize_left(refseq, pos, &allele);
                // At an HP site the column's length votes replace event votes.
                if iopts.hp_length_model && indel::hp_site(refseq, norm_pos).is_some() {
                    continue;
                }
                // `normalize_left` walks the reference with no knowledge of read boundaries.
                if indel::match_run_len(rpos, &cigar, norm_pos) == 0 {
                    stats.indel_events_out_of_span += 1;
                    continue;
                }
                // Gate on the run at the NORMALIZED site, not at the column the event was seen in.
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

        let model = fit_site(&counts, vaf);
        stats.alt_sites_preserved += model.kept_alt_count(ref_idx) as u64;
        stats.alt_alleles_strand_biased += model.strand_biased_alts as u64;

        for a in 0..4 {
            if model.kept[a] && Some(a) != ref_idx && model.freq[a] < HOMOPLASMIC_VAF {
                kept_minors.insert((tid, pos, a));
            }
        }

        let mut changed = false;
        for (qname, allele, qual, ri) in obs {
            let corr = map_allele(allele, qual, &model);
            let key = (tid, pos, allele);
            if model.pending[allele] || kept_minors.contains(&key) {
                carriers.entry(key).or_default().insert(ri);
            }
            if model.pending[allele] {
                deferred.entry(key).or_default().push((qname, corr));
            } else if corr != allele {
                record_sub(&mut corrections, &mut stats, qname, pos, allele, corr);
                changed = true;
            }
        }
        if changed {
            corrected_cols.insert((tid, pos));
        }
    }

    // Linkage rescue of pending alleles; the rest are corrected like any other error.
    let mut pending_keys: Vec<AlleleKey> = deferred.keys().copied().collect();
    pending_keys.sort();
    let rescued = rescue_linked_pending(&pending_keys, &carriers, &kept_minors, &spans);
    stats.alt_alleles_rescued_by_linkage = rescued.len() as u64;
    for key in pending_keys {
        if rescued.contains(&key) {
            continue;
        }
        let (tid, pos, allele) = key;
        for (qname, corr) in deferred.remove(&key).unwrap() {
            record_sub(&mut corrections, &mut stats, qname, pos, allele, corr);
        }
        corrected_cols.insert((tid, pos));
    }
    stats.columns_corrected = corrected_cols.len() as u64;
    // Deferred substitutions land after the pileup, so subs must be re-sorted.
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

    // Indel sites are decided once the pileup has filled in their depth.
    let mut models: IndelModels = HashMap::new();
    if iopts.enabled {
        let n_sites = sites.len();
        for ((contig, norm_pos), mut site) in sites.into_iter() {
            let Some(refseq) = refs.get(&contig) else { continue };
            // Fill the denominator from the column totals: each spanning read once.
            let (depth, fwd, rev) = depths
                .get(&(contig.clone(), norm_pos))
                .copied()
                .unwrap_or((0, 0, 0));
            // Nothing spans this site, so there is nothing to decide.
            if depth == 0 {
                continue;
            }
            site.depth = depth;
            site.col_fwd = fwd;
            site.col_rev = rev;

            stats.indel_sites_examined += 1;
            let model = indel::fit_indel_site(&site, refseq, norm_pos, n_sites, iopts);
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

/// Queue the substitution `from -> to` at `pos` for `qname` and book it.
/// Booked at decision time, so a sub later dropped by an HP-run rewrite is still counted.
fn record_sub(
    corrections: &mut Corrections,
    stats: &mut DenoiseStats,
    qname: Vec<u8>,
    pos: u32,
    from: usize,
    to: usize,
) {
    corrections.entry(qname).or_default().subs.push((pos, idx_to_base(to)));
    stats.bases_modified += 1;
    stats.substitution_matrix[from][to] += 1;
}

const RESCUE_MIN_DIST: u32 = 20;
const RESCUE_MIN_SHARED: u32 = 3;
const RESCUE_ALPHA: f64 = 1e-4;
const RESCUE_ROUNDS: usize = 3;

/// Pending alleles (sorted) co-carried with a partner allele more than chance allows.
/// Partners start as the kept minor alleles; each round adds what the previous one rescued.
fn rescue_linked_pending(
    pending: &[AlleleKey],
    carriers: &HashMap<AlleleKey, HashSet<usize>>,
    kept_minors: &HashSet<AlleleKey>,
    spans: &[(u32, u32, u32)],
) -> HashSet<AlleleKey> {
    let empty = HashSet::new();
    let mut partners: Vec<AlleleKey> = kept_minors.iter().copied().collect();
    partners.sort();
    let mut open = pending.to_vec();
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

/// Book one reassignment the rewrite actually applied (see `RewriteResult::applied`),
/// bucketed by the repeat-context length that decided it.
fn book_reassignment(
    stats: &mut DenoiseStats,
    from: &indel::Allele,
    to: &indel::Allele,
    context_len: usize,
) {
    use indel::Allele::{Del, Ins, Ref};
    match (from, to) {
        (Ref, Ref) => {}
        (Ref, Ins(_)) => stats.reads_reassigned_ref_to_ins += 1,
        (Ref, Del(_)) => stats.reads_reassigned_ref_to_del += 1,
        (_, Ref) => stats.reads_reassigned_indel_to_ref += 1,
        _ => stats.reads_reassigned_indel_to_indel += 1,
    }
    // The read changes length by the net-length difference, e.g. Del(3) -> Del(1) restores 2 bases.
    let d = to.net_len() - from.net_len();
    if d > 0 {
        stats.indel_bases_inserted += d as u64;
    } else {
        stats.indel_bases_removed += d.unsigned_abs() as u64;
    }
    if stats.reassignments_by_hp_length.len() <= context_len {
        stats.reassignments_by_hp_length.resize(context_len + 1, 0);
    }
    stats.reassignments_by_hp_length[context_len] += 1;
}

/// The site's MAP allele for a read observing `observed`, or `None` to leave the read alone.
fn move_target(
    observed: &indel::Allele,
    model: &indel::IndelSiteModel,
    iopts: &indel::IndelOpts,
    stats: &mut DenoiseStats,
) -> Option<indel::Allele> {
    match indel::assign_allele(observed, model, iopts) {
        indel::Assignment::Keep => None,
        indel::Assignment::Refused => {
            stats.assignments_skipped_unsupported_transition += 1;
            None
        }
        indel::Assignment::Move(t) => Some(t),
    }
}

/// Pass 2: write every record, rewriting correctable reads with their substitutions and indel edits.
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

    // Contigs referenced by a correctable, edit-bearing record but absent from `refs`.
    let mut missing_ref_contigs: BTreeSet<Vec<u8>> = BTreeSet::new();
    let mut missing_ref_reads = 0u64;

    for result in reader.records() {
        let mut rec = result.context("read BAM record")?;
        // Only primary mapped reads are candidates; a mapped record with tid == -1 is written through.
        if !(rec.is_secondary() || rec.is_supplementary() || rec.is_unmapped()) {
            stats.reads_processed += 1;
            if rec.tid() >= 0 {
                let contig = header_view.tid2name(rec.tid() as u32);
                match refs.get(contig) {
                    Some(refseq) => {
                        let sites = models.get(contig).map_or(&[][..], |l| l.as_slice());
                        if rewrite_record(&mut rec, refseq, corrections, sites, iopts, stats) {
                            stats.reads_modified += 1;
                        }
                    }
                    None => {
                        if corrections.get(rec.qname()).is_some_and(|e| !e.is_empty()) {
                            missing_ref_contigs.insert(contig.to_vec());
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
    Ok(())
}

/// Plan this read's indel edits against the decided `sites` of its contig, then apply them and its
/// substitutions in place. Returns whether SEQ or CIGAR changed.
fn rewrite_record(
    rec: &mut bam::Record,
    refseq: &[u8],
    corrections: &Corrections,
    sites: &[(u32, indel::IndelSiteModel)],
    iopts: &indel::IndelOpts,
    stats: &mut DenoiseStats,
) -> bool {
    let mut edits = corrections.get(rec.qname()).cloned().unwrap_or_default();
    // (ref_pos, kind_of(from)) -> the repeat-context length that decided this edit.
    let mut edit_context_len: HashMap<(u32, indel::EditKind), usize> = HashMap::new();
    // (anchor, run end, run base, from, to): HP-length moves, applied by `rewrite_hp_run` rather than
    // `rewrite_read`, which edits one CIGAR event at a time.
    let mut hp_moves: Vec<(u32, u32, u8, indel::Allele, indel::Allele)> = Vec::new();

    if iopts.enabled {
        // `end_pos()` is only defined on `CigarStringView`, which carries the record's POS.
        let view = rec.cigar();
        let rstart = rec.pos();
        let rend = view.end_pos();
        let cigar = view.take();
        let seq = rec.seq().as_bytes();
        let evs = indel::read_events(rstart, &cigar, &seq, refseq);

        // What this read carries at each normalized site: (allele, CIGAR anchor).
        let mut carried: HashMap<u32, (indel::Allele, u32)> = HashMap::new();
        // Normalized positions of this read's oversized (untouchable) indels.
        let mut blocked: Vec<i64> = Vec::new();
        // Sites where this read's own CIGAR carries two events that both normalize here.
        let mut ambiguous: HashSet<u32> = HashSet::new();
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

        let flank = indel::INDEL_FLANK as i64;
        // A gained deletion consumes up to MAX_INDEL_LEN reference bases to the right.
        let reach = flank + indel::MAX_INDEL_LEN as i64;
        // Near an oversized indel the read's alignment cannot be represented by the site's alleles.
        let near_oversized = |s: i64| blocked.iter().any(|b| (b - s).abs() <= reach);
        // Furthest reference position consumed by a gained deletion already planned for this read.
        let mut consumed_until: i64 = -1;

        let lo = sites.partition_point(|(p, _)| (*p as i64) < rstart);
        for (pos, model) in &sites[lo..] {
            let s = *pos as i64;
            // `rend` is EXCLUSIVE.
            if s >= rend {
                break;
            }
            // An HP site is keyed by run LENGTH, so the read's allele comes from `hp_allele`, as in pass 1.
            if let Some((hp_e, hp_b)) = model.hp_run {
                let observed_k = match indel::hp_allele(rstart, &cigar, &seq, refseq, *pos, hp_e) {
                    Some(k) if s - flank >= rstart => k,
                    _ => {
                        stats.assignments_skipped_span_guard += 1;
                        continue;
                    }
                };
                if s <= consumed_until {
                    stats.assignments_skipped_swallowed_by_gained_deletion += 1;
                    continue;
                }
                if near_oversized(s) {
                    stats.assignments_skipped_oversized_neighbor += 1;
                    continue;
                }
                let observed = indel::hp_length_allele(observed_k, hp_b);
                let Some(target) = move_target(&observed, model, iopts, stats) else { continue };
                hp_moves.push((*pos, hp_e, hp_b, observed, target));
                continue;
            }
            if ambiguous.contains(pos) {
                stats.assignments_skipped_ambiguous_carry += 1;
                continue;
            }
            let (observed, anchor) =
                carried.get(pos).cloned().unwrap_or((indel::Allele::Ref, *pos));
            if observed == indel::Allele::Ref && s <= consumed_until {
                stats.assignments_skipped_swallowed_by_gained_deletion += 1;
                continue;
            }
            let Some(target) = move_target(&observed, model, iopts, stats) else { continue };
            // Left-normalization only moves an event LEFT, so the CIGAR anchor may lie right of the site.
            let far = (anchor as i64).max(s);
            // Require the whole tract plus the usual flank of match context beyond it.
            let tract_need = indel::hp_tract_end(refseq, *pos) as i64 + flank;
            if s - flank < rstart || (far + reach).max(tract_need) > rend {
                stats.assignments_skipped_span_guard += 1;
                continue;
            }
            // The span check bounds only the whole alignment; a gained edit needs a match run AT the site.
            if observed == indel::Allele::Ref {
                let run = indel::match_run_len(rstart, &cigar, *pos) as i64;
                let need = match &target {
                    indel::Allele::Del(m) => *m as i64 + flank + 1,
                    _ => 1,
                };
                if run < need {
                    stats.assignments_skipped_span_guard += 1;
                    continue;
                }
            }
            if near_oversized(s) {
                stats.assignments_skipped_oversized_neighbor += 1;
                continue;
            }
            // A gained deletion consumes reference positions (s, s+m].
            if let (indel::Allele::Ref, indel::Allele::Del(m)) = (&observed, &target) {
                consumed_until = consumed_until.max(s + *m as i64);
            }
            // `rewrite_read` looks a gained edit up by its site, any other by the read's CIGAR anchor.
            let ref_pos = if observed == indel::Allele::Ref { *pos } else { anchor };
            edit_context_len
                .insert((ref_pos, indel::edit_kind(&observed)), model.context_len as usize);
            edits.indels.push(indel::IndelEdit { ref_pos, from: observed, to: target });
        }
        edits.indels.sort_by_key(|e| e.ref_pos);
    }

    if edits.is_empty() && hp_moves.is_empty() {
        return false;
    }
    let qname = rec.qname().to_vec();
    let mut cigar = rec.cigar().take();
    let mut seq = rec.seq().as_bytes();
    let mut qual = rec.qual().to_vec();

    // Half-open reference spans of the runs actually rewritten; each rewrite preserves reference consumption.
    let mut hp_runs: Vec<(u32, u32)> = Vec::new();
    for (a, e, b, from, to) in &hp_moves {
        match indel::rewrite_hp_run(rec.pos(), &cigar, &seq, &qual, *a, *e, *b, to.net_len()) {
            Some((c, s, q)) => {
                cigar = c;
                seq = s;
                qual = q;
                hp_runs.push((*a + 1, *e));
                book_reassignment(stats, from, to, (*e - *a - 1) as usize);
            }
            // Same accounting as an event edit the walk never reached: a harmless missed correction.
            None => stats.indel_edits_emitted_but_not_applied += 1,
        }
    }
    if !hp_runs.is_empty() {
        let before = edits.subs.len();
        edits.subs.retain(|&(p, _)| !hp_runs.iter().any(|&(lo, hi)| p >= lo && p < hi));
        stats.subs_dropped_in_hp_rewrite += (before - edits.subs.len()) as u64;
    }

    let out = indel::rewrite_read(rec.pos(), &cigar, &seq, &qual, &edits, refseq);
    // A non-empty edit list does not guarantee `rewrite_read` changed anything.
    let changed = out.seq != rec.seq().as_bytes() || out.cigar != rec.cigar().take();

    let mut applied = out.applied;
    for edit in &edits.indels {
        let key = (edit.ref_pos, indel::edit_kind(&edit.from));
        if !applied.remove(&key) {
            stats.indel_edits_emitted_but_not_applied += 1;
            continue;
        }
        let context_len = edit_context_len.get(&key).copied().unwrap_or(0);
        book_reassignment(stats, &edit.from, &edit.to, context_len);
    }

    rec.set(&qname, Some(&out.cigar), &out.seq, &out.qual);
    // An HP-length move shifts query offsets downstream of the run just as an event edit does.
    if out.structure_changed || !hp_runs.is_empty() {
        indel::strip_stale_tags(rec);
    }
    changed
}

/// Denoise an ONT BAM (substitutions, plus indels when enabled); other data types are copied through.
pub fn start(
    input: &PathBuf,
    output: &PathBuf,
    reference: &PathBuf,
    data_type: &str,
    vaf: f64,
    indels: &indel::IndelOpts,
    stats_out: Option<&PathBuf>,
) -> Result<()> {
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

    let (corrections, models, mut stats) =
        compute_corrections(input, &refs, vaf, indels)?;
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
    use crate::denoise_indel::Allele;
    use rust_htslib::bam::header::HeaderRecord;
    use rust_htslib::bam::record::{Aux, Cigar};
    use rust_htslib::bam::{Format, Header, Record, Writer};

    // QuickStart's homoplasmic exemption.

    /// A column with reference A (0); each `(allele, qual, fwd, rev)` adds `fwd` + `rev` observations.
    fn column(groups: &[(usize, u8, u32, u32)]) -> ColumnCounts {
        let mut c = ColumnCounts::new(Some(0));
        for &(a, q, fwd, rev) in groups {
            (0..fwd).for_each(|_| c.add(a, q, false));
            (0..rev).for_each(|_| c.add(a, q, true));
        }
        c
    }

    #[test]
    fn phred_to_e_matches_definition_and_clamps() {
        // Q0 nominally means P(error) = 1, which would make P(obs | truth = obs) zero and invert every base.
        for (q, e) in [(20u8, 0.01), (10, 0.1), (2, 0.6309573444801932), (1, E_MAX), (0, E_MAX)] {
            assert!((phred_to_e(q) - e).abs() < 1e-12, "Q{q}");
        }
    }

    #[test]
    fn zero_quality_observations_are_never_inverted() {
        // ref A 60%, alt G 40%, all at Q0 -- the shape pbsim's errhmm models emit.
        let m = fit_site(&column(&[(0, 0, 30, 30), (2, 0, 20, 20)]), 0.01);
        assert!(m.kept[0] && m.kept[2]);
        // At the clamp the likelihood is flat (1 - e == e/3), so kept alleles stand and rejected ones,
        // with no evidence and no prior, fall back to the most frequent kept allele.
        for observed in 0..4 {
            let want = if m.kept[observed] { observed } else { 0 };
            assert_eq!(map_allele(observed, 0, &m), want, "observed={observed}");
        }
    }

    #[test]
    fn low_quality_observations_of_a_kept_alt_are_not_reverted_to_ref() {
        // An 8% heteroplasmy in ONT-quality reads: the site model already ruled the alt real.
        let groups: Vec<_> = [3u8, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 20]
            .iter()
            .flat_map(|&q| [(0, q, 23, 23), (2, q, 2, 2)])
            .collect();
        let m = fit_site(&column(&groups), 0.01);
        assert!(m.kept[0] && m.kept[2]);
        for q in [0, 4, 6, 7, 30] {
            assert_eq!(map_allele(2, q, &m), 2, "kept alt reverted at Q{q}");
            assert_eq!(map_allele(0, q, &m), 0, "ref reverted at Q{q}");
        }
        // Rejected alleles are still errors.
        assert_eq!(map_allele(1, 6, &m), 0);
        assert_eq!(map_allele(3, 6, &m), 0);
    }

    #[test]
    fn pure_ref_column_with_sparse_errors_keeps_only_ref() {
        // 30 clean A plus a stray G (fwd only) and C (rev only) at low quality.
        let m = fit_site(&column(&[(0, 30, 15, 15), (2, 8, 1, 0), (1, 8, 0, 1)]), 0.01);
        assert_eq!(m.kept, [true, false, false, false]);
        assert_eq!(map_allele(2, 8, &m), 0);
        assert_eq!(map_allele(1, 8, &m), 0);
    }

    #[test]
    fn strand_balanced_heteroplasmy_is_preserved() {
        let m = fit_site(&column(&[(0, 30, 30, 30), (2, 30, 20, 20)]), 0.01);
        assert_eq!(m.kept, [true, false, true, false]);
        assert_eq!(map_allele(0, 30, &m), 0);
        assert_eq!(map_allele(2, 30, &m), 2);
        // A third-base error collapses to the more frequent allele.
        assert_eq!(map_allele(3, 8, &m), 0);
    }

    #[test]
    fn single_strand_candidate_fails_the_per_strand_gate() {
        // A 10% alt seen only on the forward strand.
        let m = fit_site(&column(&[(0, 30, 45, 45), (2, 30, 10, 0)]), 0.01);
        assert!(!m.kept[2] && !m.pending[2]);
        assert_eq!(m.strand_biased_alts, 0, "rejected by the count gate, not the strand-bias test");
        assert_eq!(map_allele(2, 30, &m), 0);
    }

    #[test]
    fn strand_bias_test_rejects_only_large_strand_effects() {
        // (ref fwd, ref rev, alt fwd, alt rev, alt kept)
        for (rf, rr, af, ar, kept) in [
            // Passes the per-strand gate, fails the bias test.
            (40, 40, 18, 2, false),
            // 75% of the column: near-homoplasmic, so exempt.
            (0, 25, 73, 2, true),
            // Deep and significant, but a strand effect of only 0.11 vs the column's 0.514.
            (2365, 2351, 328, 196, true),
            // Deep single-strand artifact: effect 0.968 - 0.468 = 0.50.
            (1722, 2454, 449, 15, false),
        ] {
            let m = fit_site(&column(&[(0, 30, rf, rr), (2, 30, af, ar)]), 0.01);
            assert_eq!(m.kept[2], kept, "alt {af}/{ar} vs ref {rf}/{rr}");
            assert_eq!(m.strand_biased_alts, !kept as usize);
            assert_eq!(map_allele(2, 30, &m), if kept { 2 } else { 0 });
        }
        // Premise of the deep kept row: only the effect-size floor saves it.
        assert!(crate::call::strand_bias_pvalue(328, 196, 2693.0 / 5240.0) < STRAND_BIAS_P);
    }

    #[test]
    fn empty_column_falls_back_to_reference_only() {
        let m = fit_site(&column(&[]), 0.01);
        assert_eq!(m.kept, [true, false, false, false]);
        assert_eq!(m.freq, [1.0, 0.0, 0.0, 0.0]);
    }

    #[test]
    fn eligible_allele_between_floor_and_keep_vaf_is_pending() {
        // (ref per strand, alt fwd, alt rev) at keep vaf 0.03 -> (kept, pending, strand-biased)
        for (r, af, ar, want) in [
            (486, 14, 14, (false, true, 0)),  // 2.8%
            (482, 18, 18, (true, false, 0)),  // 3.6%
            (998, 2, 2, (false, false, 0)),   // 0.2%: below the pending floor
            (486, 26, 2, (false, false, 1)),  // strand-biased alleles are never pending
        ] {
            let m = fit_site(&column(&[(0, 30, r, r), (2, 30, af, ar)]), 0.03);
            assert_eq!((m.kept[2], m.pending[2], m.strand_biased_alts), want, "alt {af}/{ar}");
            assert!(!m.pending[0], "reference is never pending");
        }
    }

    // ---- BAM fixtures ----

    /// Scratch directory removed on drop.
    struct Scratch(PathBuf);

    impl Scratch {
        fn new(tag: &str) -> Self {
            let dir = std::env::temp_dir()
                .join(format!("himito_denoise_{tag}_{}", std::process::id()));
            std::fs::create_dir_all(&dir).unwrap();
            Scratch(dir)
        }

        fn path(&self, file: &str) -> PathBuf {
            self.0.join(file)
        }
    }

    impl Drop for Scratch {
        fn drop(&mut self) {
            std::fs::remove_dir_all(&self.0).ok();
        }
    }

    fn read(name: &str, pos: i64, seq: &[u8], cigar: &str, reverse: bool) -> Record {
        let mut rec = Record::new();
        let cigar = CigarString::try_from(cigar).unwrap();
        rec.set(name.as_bytes(), Some(&cigar), seq, &vec![30; seq.len()]);
        rec.set_tid(0);
        rec.set_pos(pos);
        rec.set_mapq(60);
        rec.set_flags(if reverse { 16 } else { 0 });
        rec
    }

    /// Nine strand-alternating reads `r0`..`r8` with the same alignment.
    fn nine(pos: i64, seq: &[u8], cigar: &str) -> Vec<Record> {
        (0..9).map(|i| read(&format!("r{i}"), pos, seq, cigar, i % 2 == 0)).collect()
    }

    /// Write `recs` (coordinate-sorted) on the single contig `chrM`.
    fn write_bam(path: &Path, contig_len: usize, recs: &[Record]) {
        let mut header = Header::new();
        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chrM");
        sq.push_tag(b"LN", &contig_len.to_string());
        header.push_record(&sq);
        let mut w = Writer::from_path(path, &header, Format::Bam).unwrap();
        for rec in recs {
            w.write(rec).unwrap();
        }
    }

    fn refs_of(refseq: &[u8]) -> HashMap<Vec<u8>, Vec<u8>> {
        HashMap::from([(b"chrM".to_vec(), refseq.to_vec())])
    }

    /// (SEQ, CIGAR) of a record.
    fn aln(rec: &Record) -> (Vec<u8>, String) {
        (rec.seq().as_bytes(), rec.cigar().to_string())
    }

    fn al(seq: &[u8], cigar: &str) -> (Vec<u8>, String) {
        (seq.to_vec(), cigar.to_string())
    }

    /// Output records by name, each checked for SEQ/QUAL/CIGAR agreement as `samtools quickcheck` requires.
    fn read_bam(path: &Path) -> HashMap<String, Record> {
        let mut out = HashMap::new();
        for rec in Reader::from_path(path).unwrap().records() {
            let rec = rec.unwrap();
            let query: u32 = rec
                .cigar()
                .iter()
                .map(|op| match op {
                    Cigar::Match(n) | Cigar::Ins(n) | Cigar::SoftClip(n) | Cigar::Equal(n)
                    | Cigar::Diff(n) => *n,
                    _ => 0,
                })
                .sum();
            assert_eq!(rec.seq_len() as u32, query, "SEQ length must match the CIGAR");
            assert_eq!(rec.seq_len(), rec.qual().len(), "SEQ and QUAL must agree");
            out.insert(String::from_utf8(rec.qname().to_vec()).unwrap(), rec);
        }
        out
    }

    struct Run {
        corr: Corrections,
        models: IndelModels,
        stats: DenoiseStats,
        out: HashMap<String, Record>,
    }

    impl Run {
        fn site(&self, pos: u32) -> Option<&indel::IndelSiteModel> {
            let list = self.models.get(b"chrM".as_slice())?;
            list.iter().find(|(p, _)| *p == pos).map(|(_, m)| m)
        }

        fn subs(&self, name: &str) -> Vec<(u32, u8)> {
            self.corr.get(name.as_bytes()).map(|e| e.subs.clone()).unwrap_or_default()
        }

        fn aln(&self, name: &str) -> (Vec<u8>, String) {
            aln(&self.out[name])
        }

        fn assert_untouched(&self, recs: &[Record]) {
            for rec in recs {
                let name = std::str::from_utf8(rec.qname()).unwrap();
                assert_eq!(self.aln(name), aln(rec), "{name} must be untouched");
            }
        }
    }

    /// Both passes over `recs` against `refseq`.
    fn run(tag: &str, recs: &[Record], refseq: &[u8], vaf: f64, iopts: indel::IndelOpts) -> Run {
        let s = Scratch::new(tag);
        let (inb, outb) = (s.path("in.bam"), s.path("out.bam"));
        write_bam(&inb, refseq.len(), recs);
        let refs = refs_of(refseq);
        let (corr, models, mut stats) =
            compute_corrections(&inb, &refs, vaf, &iopts).unwrap();
        apply_corrections(&inb, &outb, &corr, &refs, &models, &iopts, &mut stats).unwrap();
        Run { corr, models, stats, out: read_bam(&outb) }
    }

    /// Pass 2 alone, with hand-built corrections and decided indel sites.
    fn run_apply(
        tag: &str,
        recs: &[Record],
        refseq: &[u8],
        corr: Corrections,
        sites: Vec<(u32, indel::IndelSiteModel)>,
        iopts: indel::IndelOpts,
    ) -> Run {
        let s = Scratch::new(tag);
        let (inb, outb) = (s.path("in.bam"), s.path("out.bam"));
        write_bam(&inb, refseq.len(), recs);
        let models = IndelModels::from([(b"chrM".to_vec(), sites)]);
        let mut stats = DenoiseStats::default();
        apply_corrections(&inb, &outb, &corr, &refs_of(refseq), &models, &iopts, &mut stats)
            .unwrap();
        Run { corr, models, stats, out: read_bam(&outb) }
    }

    /// `start` on ONT data; returns the stats JSON and the output records.
    fn run_start(
        tag: &str,
        recs: &[Record],
        refseq: &[u8],
        iopts: indel::IndelOpts,
    ) -> (serde_json::Value, HashMap<String, Record>) {
        let s = Scratch::new(tag);
        let (inb, outb, reff, statsf) =
            (s.path("in.bam"), s.path("out.bam"), s.path("ref.fa"), s.path("stats.json"));
        write_bam(&inb, refseq.len(), recs);
        std::fs::write(&reff, [&b">chrM\n"[..], refseq, b"\n"].concat()).unwrap();
        start(&inb, &outb, &reff, "ont-r10", 0.01, &iopts, Some(&statsf)).unwrap();
        let json = serde_json::from_str(&std::fs::read_to_string(&statsf).unwrap()).unwrap();
        (json, read_bam(&outb))
    }

    // ---- substitutions ----

    /// 400 reads of 100 A's on alternating strands; reads in `m.1` carry C at `m.0`, reads in `p.1` G at `p.0`.
    fn linkage_run(
        tag: &str,
        m: (usize, std::ops::Range<usize>),
        p: (usize, std::ops::Range<usize>),
        extra: &[(usize, usize, u8)],
    ) -> Run {
        let recs: Vec<Record> = (0..400)
            .map(|i| {
                let mut s = vec![b'A'; 100];
                if m.1.contains(&i) { s[m.0] = b'C'; }
                if p.1.contains(&i) { s[p.0] = b'G'; }
                for &(r, pos, b) in extra {
                    if r == i { s[pos] = b; }
                }
                read(&format!("r{i:03}"), 0, &s, "100M", i % 2 == 1)
            })
            .collect();
        run(tag, &recs, &[b'A'; 100], 0.03, indel::IndelOpts::default())
    }

    #[test]
    fn pending_allele_riding_on_a_kept_minor_is_preserved() {
        // C@10 is kept (3.5%); G@60 is pending (2.75%) and carried only by C@10 reads.
        let r = linkage_run("link_rescued", (10, 0..14), (60, 0..11), &[]);
        assert_eq!(r.stats.alt_alleles_rescued_by_linkage, 1);
        assert_eq!(r.stats.bases_modified, 0);
    }

    #[test]
    fn unlinked_pending_allele_is_still_corrected_and_subs_stay_sorted() {
        // r100 also carries a lone T@80, corrected at once; its deferred G@60 lands after it and must be re-sorted.
        let r = linkage_run("link_unlinked", (10, 0..14), (60, 100..111), &[(100, 80, b'T')]);
        assert_eq!(r.subs("r100"), [(60, b'A'), (80, b'A')]);
        for i in 101..111 {
            assert_eq!(r.subs(&format!("r{i}")), [(60, b'A')]);
        }
        assert_eq!(r.stats.alt_alleles_rescued_by_linkage, 0);
        assert_eq!(r.stats.bases_modified, 12);
        assert_eq!(r.stats.substitution_matrix[2][0], 11);
        assert_eq!(r.stats.columns_corrected, 2);
    }

    #[test]
    fn linked_partner_within_20bp_does_not_rescue() {
        let r = linkage_run("link_near", (10, 0..14), (25, 0..11), &[]);
        assert_eq!(r.stats.alt_alleles_rescued_by_linkage, 0);
        for i in 0..11 {
            assert_eq!(r.subs(&format!("r{i:03}")), [(25, b'A')]);
        }
    }

    #[test]
    fn linkage_rescue_is_deterministic() {
        let extra: Vec<_> = (200..211).map(|r| (r, 85, b'T')).collect();
        let summary = || {
            let r = linkage_run("link_det", (10, 0..14), (60, 0..11), &extra);
            let mut subs: Vec<_> = r.corr.into_iter().map(|(k, e)| (k, e.subs)).collect();
            subs.sort();
            (subs, serde_json::to_string(&r.stats).unwrap())
        };
        let first = summary();
        assert_eq!(first, summary());
        // Only the unlinked T@85 carriers are corrected.
        assert_eq!(first.0.len(), 11);
        assert!(first.0.iter().all(|(_, s)| s == &[(85u32, b'A')]));
    }

    #[test]
    fn single_base_error_is_corrected_at_its_reference_position() {
        // Reads start at POS 10, so the G at query offset 2 is reference position 12. It is fwd-only,
        // so it fails the per-strand gate.
        let mut recs = nine(10, b"AAAAAA", "6M");
        let mut err = read("rerr", 10, b"AAGAAA", "6M", false);
        err.push_aux(b"XY", Aux::I32(7)).unwrap();
        recs.push(err);
        let r = run("sub_refpos", &recs, &[b'A'; 20], 0.01, indel::IndelOpts::default());
        assert_eq!(r.subs("rerr"), [(12, b'A')]);
        assert_eq!(r.stats.bases_modified, 1);
        assert_eq!(r.stats.substitution_matrix[2][0], 1);
        assert_eq!((r.stats.reads_processed, r.stats.reads_modified), (10, 1));
        assert_eq!(r.aln("rerr"), al(b"AAAAAA", "6M"));
        assert!(matches!(r.out["rerr"].aux(b"XY").unwrap(), Aux::I32(7)), "aux tags survive");
    }

    #[test]
    fn start_pacbio_is_byte_identical_passthrough() {
        let s = Scratch::new("pacbio");
        let (inb, outb, reff) = (s.path("in.bam"), s.path("out.bam"), s.path("ref.fa"));
        std::fs::write(&reff, ">chrM\nAAAAAA\n").unwrap();
        write_bam(&inb, 6, &[read("r0", 0, b"AAGAAA", "6M", false)]);
        start(&inb, &outb, &reff, "pacbio", 0.01, &indel::IndelOpts::default(), None)
            .unwrap();
        assert_eq!(std::fs::read(&inb).unwrap(), std::fs::read(&outb).unwrap());
    }

    #[test]
    fn start_ont_corrects_substitutions_without_touching_the_cigar() {
        // With indel correction off, `=`/`X` ops must come out byte-identical.
        let mut recs = nine(0, b"AAAAAA", "6M");
        recs.push(read("rerr", 0, b"AAGAAA", "2=1X3=", false));
        let (stats, out) = run_start("start_ont", &recs, b"AAAAAA", indel::IndelOpts::default());
        assert_eq!(aln(&out["rerr"]), al(b"AAAAAA", "2=1X3="));
        assert_eq!(stats["bases_modified"], 1);
    }

    #[test]
    fn records_without_a_reference_contig_are_written_through() {
        // `rerr`'s contig is absent from the reference; `rtidless` is mapped with tid == -1, which must
        // never reach `tid2name`.
        let s = Scratch::new("no_contig");
        let (inb, outb) = (s.path("in.bam"), s.path("out.bam"));
        let mut tidless = read("rtidless", 0, b"AAAAAA", "6M", false);
        tidless.set_tid(-1);
        write_bam(&inb, 6, &[read("rerr", 0, b"AAGAAA", "6M", false), tidless]);
        let edits = indel::ReadEdits { subs: vec![(2, b'A')], indels: vec![] };
        let corr = Corrections::from([(b"rerr".to_vec(), edits)]);
        let mut stats = DenoiseStats::default();
        let iopts = indel::IndelOpts::default();
        apply_corrections(&inb, &outb, &corr, &HashMap::new(), &IndelModels::new(), &iopts, &mut stats)
            .unwrap();
        assert_eq!((stats.reads_processed, stats.reads_modified), (2, 0));
        let out = read_bam(&outb);
        assert_eq!(aln(&out["rerr"]), al(b"AAGAAA", "6M"));
        assert_eq!(aln(&out["rtidless"]), al(b"AAAAAA", "6M"));
    }

    // ---- indels ----

    // An 8bp A-run at [8, 16) after the anchor T at 7: its single HP site is 7.
    const REF29: &[u8] = b"ACGTACGTAAAAAAAACGTACGTACGTAC";
    const HP_INS: &[u8] = b"ACGTACGTAAAAAAAAACGTACGTACGTAC"; // 16M1I13M: one extra A in the run
    const HP_DEL: &[u8] = b"ACGTACGTAAAAAAACGTACGTACGTAC"; // 8M1D20M: one A missing from the run
    const UNIQ_INS: &[u8] = b"ACGTACGTAAAAAAAACGTGACGTACGTAC"; // 19M1I10M: a G in unique sequence

    fn indels(hp_length_model: bool) -> indel::IndelOpts {
        indel::IndelOpts { enabled: true, hp_length_model, ..Default::default() }
    }

    fn ins_a() -> Allele {
        Allele::Ins(b"A".to_vec())
    }

    fn freq(list: &[(Allele, f64)], a: &Allele) -> Option<f64> {
        list.iter().find(|(k, _)| k == a).map(|(_, f)| *f)
    }

    /// Non-repetitive reference, so no hand-placed site sits in a run.
    fn acgt(len: usize) -> Vec<u8> {
        (0..len).map(|i| b"ACGT"[i % 4]).collect()
    }

    /// A hand-decided site keeping `alt` at 0.9 over Ref.
    fn site_model(alt: Allele, context_len: u32) -> indel::IndelSiteModel {
        let kept = vec![(Allele::Ref, 0.1), (alt, 0.9)];
        indel::IndelSiteModel {
            obs_vaf: kept.clone(),
            kept,
            eps: 0.15,
            context_len,
            strand_biased: 0,
            strand_rejected: vec![],
            hp_run: None,
        }
    }

    #[test]
    fn homoplasmic_insertion_is_kept_and_repairs_the_read_that_dropped_it() {
        // 9 of 10 reads carry an extra A in the run. Depth counts each spanning read once (the pileup also
        // reports these reads at the anchor column), and the HP length model reaches the same frequencies.
        let mut recs = nine(0, HP_INS, "16M1I13M");
        recs.push(read("rdrop", 0, REF29, "29M", false));
        for hp in [false, true] {
            let r = run(&format!("hp_ins_kept_{hp}"), &recs, REF29, 0.01, indels(hp));
            let m = r.site(7).expect("site 7");
            assert_eq!(freq(&m.kept, &ins_a()), Some(0.9), "hp={hp}");
            assert_eq!(m.context_len, 8);
            assert_eq!(m.hp_run, hp.then_some((16, b'A')));
            let hp_sites: Vec<u32> = r.models[b"chrM".as_slice()]
                .iter()
                .filter(|(_, m)| m.hp_run.is_some())
                .map(|(p, _)| *p)
                .collect();
            assert_eq!(hp_sites, if hp { vec![7] } else { vec![] }, "one HP site per run");
            assert_eq!(r.stats.indel_alt_alleles_kept, 1);
            assert!(r.out.values().all(|rec| rec.seq().as_bytes() == HP_INS), "hp={hp}");
            assert_eq!((r.stats.reads_reassigned_ref_to_ins, r.stats.indel_bases_inserted), (1, 1));
        }
    }

    #[test]
    fn an_undefined_hp_allele_does_not_count_as_reference_support() {
        // Two reads ending at 18 pass the depth gate (17) but not the HP vote gate (21): they cannot show
        // the run's length, so `rdrop` is the only reference observer.
        let mut recs = nine(0, HP_INS, "16M1I13M");
        recs.push(read("rdrop", 0, REF29, "29M", false));
        recs.push(read("rflank0", 0, &REF29[..18], "18M", false));
        recs.push(read("rflank1", 0, &REF29[..18], "18M", true));
        let r = run("hp_undefined", &recs, REF29, 0.01, indels(true));
        let m = r.site(7).expect("site 7");
        assert_eq!(freq(&m.kept, &ins_a()), Some(0.9));
        assert_eq!(freq(&m.kept, &Allele::Ref), Some(0.1));
    }

    #[test]
    fn site_depth_excludes_reads_that_end_inside_the_homopolymer() {
        // Ten reads end five bases into the run, and `sins` also carries an extra A there (13M1I1M):
        // none can show the run's length, so the insertion stays at 0.9.
        let mut recs = nine(0, HP_INS, "16M1I13M");
        recs.push(read("rdrop", 0, REF29, "29M", false));
        recs.extend((0..10).map(|i| read(&format!("s{i}"), 0, &REF29[..13], "13M", i % 2 == 0)));
        recs.push(read("sins", 0, b"ACGTACGTAAAAAAA", "13M1I1M", true));
        let r = run("depth_short", &recs, REF29, 0.01, indels(false));
        let m = r.site(7).expect("site 7");
        assert_eq!(freq(&m.kept, &ins_a()), Some(0.9), "reads ending in the run counted as REF: {m:?}");
        assert_eq!(m.context_len, 8);
        assert_eq!(r.stats.indel_events_not_spanning_repeat, 1, "only `sins`'s vote is dropped");
    }

    #[test]
    fn alt_votes_are_span_gated_at_the_normalized_site() {
        // 15M1I1M anchors its insertion on the run's penultimate base (14), where the reads (end 16) look
        // spanning, but it normalizes to site 7, whose run they do not span. Without a spanning read the
        // site is booked as dropped for lack of depth.
        for spanning in [true, false] {
            let mut recs = nine(0, b"ACGTACGTAAAAAAAAA", "15M1I1M");
            if spanning {
                recs.push(read("rspan", 0, REF29, "29M", false));
            }
            let r = run(&format!("norm_gate_{spanning}"), &recs, REF29, 0.01, indels(false));
            assert_eq!(r.stats.indel_events_not_spanning_repeat, 9);
            assert!(r.site(7).is_none(), "no site from reads that never spanned its run");
            assert_eq!(r.stats.indel_sites_dropped_no_spanning_depth, !spanning as u64);
        }
    }

    #[test]
    fn models_are_sorted_by_position_within_each_contig() {
        // The record walk binary-searches each contig's list. Reads carry the HP insertion and a G after 19.
        let mut recs = nine(0, b"ACGTACGTAAAAAAAAACGTAGCGTACGTAC", "16M1I4M1I9M");
        recs.push(read("rclean", 0, REF29, "29M", false));
        let r = run("sorted", &recs, REF29, 0.01, indels(false));
        let positions: Vec<u32> = r.models[b"chrM".as_slice()].iter().map(|(p, _)| *p).collect();
        assert_eq!(positions, [7, 19]);
    }

    #[test]
    fn indel_event_normalizing_outside_its_own_span_is_not_fabricated_into_a_site() {
        // Four reads start at 10, inside the run, with an extra A anchored there: it normalizes to 7, before
        // the reads begin.
        let mut recs = vec![read("rfull", 0, REF29, "29M", false)];
        let seq = [&REF29[10..11], b"A", &REF29[11..]].concat();
        recs.extend((0..4).map(|i| read(&format!("r{i}"), 10, &seq, "1M1I18M", i % 2 == 0)));
        let r = run("out_of_span", &recs, REF29, 0.01, indels(false));
        assert_eq!(r.stats.indel_events_out_of_span, 4);
        assert!(r.site(7).is_none());
        r.assert_untouched(&recs[..1]);
    }

    #[test]
    fn read_deleting_the_normalized_column_does_not_cast_a_phantom_alt_vote() {
        // `rdel` deletes REF29[7] and inserts an A after 15; the insertion normalizes to 7, which the read
        // itself deletes.
        let mut recs = nine(0, REF29, "29M");
        let seq = [&REF29[..7], &REF29[8..16], b"A", &REF29[16..]].concat();
        recs.push(read("rdel", 0, &seq, "7M1D8M1I13M", false));
        let r = run("phantom", &recs, REF29, 0.01, indels(false));
        assert_eq!(r.stats.indel_events_out_of_span, 1);
        assert!(r.site(7).is_none());
    }

    #[test]
    fn read_with_two_events_at_one_site_votes_once_and_is_left_untouched() {
        // `rdup`'s insertion and deletion both left-normalize to site 7: pass 1 records only its first vote,
        // and pass 2 refuses to pick between them.
        let mut recs = nine(0, REF29, "29M");
        recs.push(read("rdup", 0, REF29, "9M1I4M1D15M", false));
        let r = run("dup_site", &recs, REF29, 0.01, indels(false));
        assert_eq!(r.stats.indel_events_duplicate_site, 1);
        let m = r.site(7).expect("site 7");
        assert_eq!(freq(&m.obs_vaf, &ins_a()), Some(0.1));
        assert_eq!(freq(&m.obs_vaf, &Allele::Del(1)), None);
        assert!(r.stats.assignments_skipped_ambiguous_carry >= 1);
        r.assert_untouched(&recs);
    }

    #[test]
    fn lone_spurious_indel_is_reverted_to_reference() {
        // The HP insertion's CIGAR anchor (15) differs from its normalized site (7); reverting the deletion
        // of REF29[7] restores a base.
        for (seq, cigar, removed, inserted) in [
            (HP_INS.to_vec(), "16M1I13M", 1, 0),
            ([&REF29[..7], &REF29[8..]].concat(), "7M1D21M", 0, 1),
        ] {
            let mut recs = nine(0, REF29, "29M");
            recs.push(read("rodd", 0, &seq, cigar, false));
            let r = run(&format!("lone_{cigar}"), &recs, REF29, 0.01, indels(false));
            assert!(r.out.values().all(|rec| aln(rec) == al(REF29, "29M")), "{cigar}");
            assert_eq!(r.stats.reads_reassigned_indel_to_ref, 1);
            assert_eq!((r.stats.indel_bases_removed, r.stats.indel_bases_inserted), (removed, inserted));
        }
    }

    #[test]
    fn start_removes_a_spurious_insertion_only_with_indels_on() {
        // Every column is substitution-clean, so only the indel path can touch `rins`.
        let mut recs = nine(0, REF29, "29M");
        recs.push(read("rins", 0, UNIQ_INS, "19M1I10M", false));
        for (iopts, want) in [
            (indel::IndelOpts::default(), al(UNIQ_INS, "19M1I10M")),
            (indels(false), al(REF29, "29M")),
        ] {
            let on = iopts.enabled;
            let (stats, out) = run_start(&format!("start_indels_{on}"), &recs, REF29, iopts);
            assert_eq!(aln(&out["rins"]), want);
            assert_eq!(stats["columns_skipped"], stats["columns_processed"]);
            assert_eq!(stats["indel_sites_examined"] != 0, on);
            assert_eq!(stats["reads_reassigned_indel_to_ref"], on as u64);
            assert_eq!(stats["indel_bases_removed"], on as u64);
        }
    }

    #[test]
    fn read_carrying_an_oversized_indel_is_not_assigned_nearby() {
        // `rbig`'s 6bp insertion (>= MAX_INDEL_LEN) is anchored at 8, inside the run whose site 7 keeps a 90% insertion.
        let mut recs = nine(0, HP_INS, "16M1I13M");
        recs.push(read("rbig", 0, &[&REF29[..9], b"GGGGGG", &REF29[9..]].concat(), "9M6I20M", false));
        let r = run("oversized", &recs, REF29, 0.01, indels(false));
        assert_eq!(r.stats.indel_alt_alleles_kept, 1, "premise: the insertion at 7 is kept");
        assert!(r.stats.assignments_skipped_oversized_neighbor >= 1);
        r.assert_untouched(&recs[9..]);
    }

    #[test]
    fn hp_wrong_direction_read_is_corrected_only_by_the_length_model() {
        // The site keeps Del(1); moving the lone Ins(1) read there is a cross-kind edit `rewrite_read` refuses.
        let mut recs: Vec<Record> =
            (0..8).map(|i| read(&format!("r{i}"), 0, HP_DEL, "8M1D20M", i % 2 == 0)).collect();
        recs.push(read("rins", 0, HP_INS, "16M1I13M", false));
        recs.push(read("rref", 0, REF29, "29M", true));

        let off = run("hp_dir_off", &recs, REF29, 0.01, indels(false));
        assert_eq!(off.aln("rins"), al(HP_INS, "16M1I13M"));
        assert!(off.stats.assignments_skipped_unsupported_transition >= 1);

        let on = run("hp_dir_on", &recs, REF29, 0.01, indels(true));
        // Rewritten to the kept length, in canonical left-aligned form.
        assert_eq!(on.aln("rins"), al(HP_DEL, "8M1D20M"));
        assert_eq!(on.stats.assignments_skipped_unsupported_transition, 0);
        assert_eq!(on.stats.reads_reassigned_indel_to_indel, 1);
        // Both moves, rins and rref, are booked against the 8bp run.
        assert_eq!(on.stats.reassignments_by_hp_length.get(8), Some(&2));
        assert_eq!(on.stats.indel_edits_emitted_but_not_applied, 0);
        on.assert_untouched(&recs[..8]);
    }

    #[test]
    fn genuine_heteroplasmic_indel_survives_candidacy_end_to_end() {
        // 12 of 20 reads carry the HP insertion: 60% must clear candidacy at context 8 on its own, not just
        // be shielded by the `protect_vaf` backstop.
        let recs: Vec<Record> = (0..20)
            .map(|i| {
                let (seq, cigar) = if i < 12 { (HP_INS, "16M1I13M") } else { (REF29, "29M") };
                read(&format!("r{i:02}"), 0, seq, cigar, i % 2 == 0)
            })
            .collect();
        let r = run("genuine_het", &recs, REF29, 0.01, indels(false));
        let m = r.site(7).expect("site 7");
        assert_eq!(m.context_len, 8);
        assert_eq!(freq(&m.obs_vaf, &ins_a()), Some(0.6));
        assert!(freq(&m.kept, &ins_a()).is_some());
        assert_eq!(r.stats.reads_modified, 0);
        r.assert_untouched(&recs);
    }

    /// Site 7 of REF29 under the HP length model: Del(1) dominates; Ins(1) was strand-rejected.
    fn hp_del_site_model() -> indel::IndelSiteModel {
        indel::IndelSiteModel {
            kept: vec![(Allele::Ref, 0.111), (Allele::Del(1), 0.889)],
            eps: indel::error_rate(8),
            context_len: 8,
            strand_biased: 0,
            obs_vaf: vec![(Allele::Ref, 0.1), (ins_a(), 0.1), (Allele::Del(1), 0.8)],
            strand_rejected: vec![ins_a()],
            hp_run: Some((16, b'A')),
        }
    }

    #[test]
    fn hp_site_move_rewrites_the_run_and_drops_substitutions_inside_it() {
        // `rewrite_hp_run` moves Ins(1) -> Del(1) before `rewrite_read`; the sub inside the run (10) is
        // dropped, the one outside it (20) still lands.
        let edits = indel::ReadEdits { subs: vec![(10, b'G'), (20, b'T')], indels: vec![] };
        let corr = Corrections::from([(b"rins".to_vec(), edits)]);
        let recs = [read("rins", 0, HP_INS, "16M1I13M", false)];
        let r = run_apply("hp_move", &recs, REF29, corr, vec![(7, hp_del_site_model())], indels(true));
        assert_eq!(r.stats.reads_modified, 1);
        assert_eq!(r.stats.reads_reassigned_indel_to_indel, 1);
        // The inserted base leaves the read and one reference base becomes a deletion.
        assert_eq!(r.stats.indel_bases_removed, 2);
        assert_eq!(r.stats.reassignments_by_hp_length.get(8), Some(&1));
        assert_eq!(r.stats.subs_dropped_in_hp_rewrite, 1);
        assert_eq!(r.stats.indel_edits_emitted_but_not_applied, 0);
        let (seq, cigar) = r.aln("rins");
        assert_eq!(cigar, "8M1D20M", "canonical left-aligned form");
        assert_eq!(seq[19], b'T');
    }

    #[test]
    fn gained_insertion_needs_the_whole_run_plus_flank() {
        // Site 7 keeps an insertion in the run [8, 16): the read must reach 16 + flank = 21.
        for (len, applied) in [(18usize, false), (22, true)] {
            let recs = [read("r", 0, &REF29[..len], &format!("{len}M"), false)];
            let sites = vec![(7, site_model(ins_a(), 8))];
            let r = run_apply(&format!("tract_{len}"), &recs, REF29, Corrections::new(), sites, indels(false));
            assert_eq!(r.stats.reads_modified, applied as u64, "len {len}");
            assert_eq!(r.stats.reads_reassigned_ref_to_ins, applied as u64);
            assert_eq!(r.stats.assignments_skipped_span_guard, !applied as u64);
            assert_eq!(r.out["r"].seq_len(), len + applied as usize);
        }
    }

    #[test]
    fn gained_edit_needs_a_match_run_at_the_site() {
        // An insertion needs an aligned base at the site; a gained Del(m) needs m + flank + 1 matched bases.
        for (cigar, qlen, site, alt, applied) in [
            ("10M3D10M", 20, 12, ins_a(), false),      // 12 lies in the read's own deletion
            ("25M3N10M", 35, 18, Allele::Del(2), false), // 7 matched bases from 18
            ("26M3N10M", 36, 18, Allele::Del(2), true),  // 8
        ] {
            let recs = [read("r", 0, &vec![b'A'; qlen], cigar, false)];
            let sites = vec![(site, site_model(alt, 1))];
            let r = run_apply(&format!("match_run_{cigar}"), &recs, &acgt(40), Corrections::new(), sites, indels(false));
            assert_eq!(r.stats.reads_modified, applied as u64, "{cigar}");
            assert_eq!(r.stats.assignments_skipped_span_guard, !applied as u64, "{cigar}");
        }
    }

    #[test]
    fn shrinking_a_deletion_is_counted_as_bases_inserted_not_removed() {
        // Del(3) -> Del(1) splices 2 reference bases back into the read.
        let recs = [read("r", 0, &[b'A'; 20], "10M3D10M", false)];
        let sites = vec![(9, site_model(Allele::Del(1), 1))];
        let r = run_apply("shrink", &recs, &acgt(40), Corrections::new(), sites, indels(false));
        assert_eq!(r.stats.reads_modified, 1);
        assert_eq!(r.aln("r").1, "10M1D12M");
        assert_eq!(r.stats.reads_reassigned_indel_to_indel, 1);
        assert_eq!((r.stats.indel_bases_inserted, r.stats.indel_bases_removed), (2, 0));
    }

    #[test]
    fn gained_deletion_swallows_a_decided_site_inside_its_own_span() {
        // A gained Del(3) at site 5 consumes reference 6..=8, so site 7's insertion is never evaluated.
        let recs = [read("r", 0, &[b'A'; 30], "30M", false)];
        let sites = vec![(5, site_model(Allele::Del(3), 1)), (7, site_model(ins_a(), 1))];
        let r = run_apply("swallowed", &recs, &acgt(40), Corrections::new(), sites, indels(false));
        assert_eq!(r.stats.reads_modified, 1);
        assert_eq!(r.aln("r").1, "6M3D21M");
        assert_eq!((r.stats.reads_reassigned_ref_to_del, r.stats.indel_bases_removed), (1, 3));
        assert_eq!((r.stats.reads_reassigned_ref_to_ins, r.stats.indel_bases_inserted), (0, 0));
        assert!(r.stats.assignments_skipped_swallowed_by_gained_deletion >= 1);
    }

    #[test]
    fn edit_emitted_for_an_unreachable_ref_pos_is_a_harmless_missed_correction() {
        // ref_pos 11 lies in the read's own deletion, so the rewrite walk never reaches it.
        let edit = indel::IndelEdit { ref_pos: 11, from: Allele::Ref, to: ins_a() };
        let corr = Corrections::from([(b"r".to_vec(), indel::ReadEdits { subs: vec![], indels: vec![edit] })]);
        let recs = [read("r", 0, &[b'A'; 20], "10M3D10M", false)];
        let r = run_apply("unreachable", &recs, &[b'A'; 30], corr, vec![], indels(false));
        assert_eq!(r.stats.reads_modified, 0);
        r.assert_untouched(&recs);
        assert_eq!(r.stats.indel_edits_emitted_but_not_applied, 1);
        assert_eq!((r.stats.reads_reassigned_ref_to_ins, r.stats.indel_bases_inserted), (0, 0));
        assert!(r.stats.reassignments_by_hp_length.is_empty());
    }

    #[test]
    fn two_emitted_edits_sharing_ref_pos_and_kind_are_not_double_counted() {
        // Both reverts key on (9, Ins): the rewrite applies one, and only one may be booked.
        let edits = [b"T", b"C"]
            .map(|b| indel::IndelEdit { ref_pos: 9, from: Allele::Ins(b.to_vec()), to: Allele::Ref })
            .to_vec();
        let corr = Corrections::from([(b"r".to_vec(), indel::ReadEdits { subs: vec![], indels: edits })]);
        let seq = [&[b'A'; 10][..], b"TC", &[b'A'; 10]].concat();
        let recs = [read("r", 0, &seq, "10M1I2D1I10M", false)];
        let r = run_apply("dup_edit", &recs, &[b'A'; 30], corr, vec![], indels(false));
        assert_eq!(r.stats.reads_modified, 1);
        assert_eq!(r.stats.reads_reassigned_indel_to_ref, 1);
        assert_eq!(r.stats.indel_edits_emitted_but_not_applied, 1);
    }
}
