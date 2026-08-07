use std::collections::HashMap;
use std::path::{Path, PathBuf};
use anyhow::{Context, Result};
use bio::io::fasta;
use log::{info, warn};
use rust_htslib::bam::{self, Read, Reader};
use rust_htslib::bam::record::CigarString;

use crate::denoise::rewrite::{rewrite_read, strip_stale_tags};

pub(crate) mod indel;
pub(crate) mod rewrite;

pub(crate) const NQ: usize = 94; // Phred qualities 0..=93

/// Ceiling on the per-base error probability. Phred 0 literally means P(error) = 1,
/// which would make P(observed | truth = observed) = 1 - e = 0: the model would assert
/// the observed base is impossible and `map_allele` would rewrite EVERY base to some
/// other allele, corrupting reads instead of correcting them. With four alleles the
/// most error a base can carry is uniform-over-4, so e is capped at 3/4. At the cap the
/// likelihood is flat (1 - e == e/3 == 0.25), i.e. the base carries no evidence and
/// correction falls back to the site's frequency prior.
///
/// This is reachable in practice: pbsim3's errhmm models (and any BAM with QUAL '*')
/// emit an all-Q0 quality string.
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

    // Used by the pileup pass (next task); harmless-until-then dead-code allow.
    #[allow(dead_code)]
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
    /// Non-ref alleles that cleared the per-strand count gate but were rejected as
    /// strand-biased; surfaced in DenoiseStats so the gate's impact is auditable.
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
    // The null for the strand-bias test is this column's OWN forward fraction, not a
    // library-wide constant: mito coverage is frequently asymmetric, and a locally
    // skewed column would otherwise make every allele in it look biased. Empty columns
    // fall back to a balanced null (the test is vacuous there anyway).
    let (col_fwd, col_rev) = c.strand_totals();
    let col_total = col_fwd + col_rev;
    let expected_fwd_frac = if col_total == 0 {
        0.5
    } else {
        col_fwd as f64 / col_total as f64
    };

    // Eligibility: reference always eligible; a non-ref allele needs >= min_strand
    // observations on EACH strand (baldur's strand-balance candidacy) AND must not be
    // significantly strand-skewed by the same binomial test call.rs applies to variants
    // (call::strand_bias_pvalue), which catches skews the raw count gate lets through.
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
        // Mirrors call::filter_strand_bias's near-homoplasmic exemption: a single-strand
        // artifact is absent from the other strand's reads, so it can never reach a high
        // apparent frequency. Skew in an allele carried by most reads at the locus is not
        // evidence of artifact, and testing it would correct away a real variant.
        let prop = if col_total == 0 {
            0.0
        } else {
            (fwd + rev) as f64 / col_total as f64
        };
        if strand_bias_p > 0.0
            && prop < homoplasmic_vaf
            && crate::call::strand_bias_pvalue(fwd as usize, rev as usize, expected_fwd_frac)
                < strand_bias_p
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

    // EM. Observations of non-eligible alleles (errors) contribute their
    // responsibility to the eligible alleles; frequencies normalize over all obs.
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
    for a in 0..4 {
        if Some(a) == c.ref_idx {
            kept[a] = true;
        } else if eligible[a] && f[a] >= vaf {
            kept[a] = true;
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

    SiteModel { freq, kept, strand_biased_alts }
}

fn map_allele(observed: usize, q: u8, m: &SiteModel) -> usize {
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

pub type Corrections = HashMap<Vec<u8>, rewrite::ReadEdits>;

/// Decided indel sites, per contig, SORTED ASCENDING BY POSITION so the write-back
/// pass can binary-search the sites a read spans.
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
    pub substitution_matrix: [[u64; 4]; 4], // [from][to]
    // --- indel denoising (all zero unless --indels is enabled) ---
    pub indel_sites_examined: u64,
    /// Non-reference ALLELES (not sites) that survived candidacy across all sites.
    pub indel_alt_alleles_kept: u64,
    /// Non-reference ALLELES (not sites) rejected specifically as strand-biased.
    pub indel_alt_alleles_strand_biased: u64,
    pub indel_events_examined: u64,
    /// Votes dropped in pass 1 because the read that cast them has no ALIGNED QUERY
    /// BASE at the column its event left-normalized into -- not merely because that
    /// column falls outside the read's overall `[pos, end_pos())` span (see
    /// CRITICAL 2), but also when the read's span covers the column only through
    /// its OWN `D`/`N` there, where depth was never incremented for this read
    /// either (see FIX 2). Depth already excludes both cases, so this only tracks
    /// the ALT side.
    pub indel_events_out_of_span: u64,
    pub reads_reassigned_ref_to_ins: u64,
    pub reads_reassigned_ref_to_del: u64,
    pub reads_reassigned_indel_to_ref: u64,
    pub reads_reassigned_indel_to_indel: u64,
    pub indel_bases_inserted: u64,
    pub indel_bases_removed: u64,
    /// Index = repeat-context length L; audits how much correction fired in long runs.
    pub reassignments_by_hp_length: Vec<u64>,
    /// (read, site) PAIRS skipped because the read lacks the match-context flank
    /// the rewrite walk needs on one side or the other of the edit.
    pub assignments_skipped_span_guard: u64,
    /// (read, site) PAIRS skipped because the MAP target is only reachable through
    /// a transition the rewrite walk does not support (cross-kind, or growing a
    /// deletion).
    pub assignments_skipped_unsupported_transition: u64,
    /// (read, site) PAIRS skipped because the read carries its own untouchable
    /// oversized indel close enough to this site that its alignment in this
    /// window cannot be represented by the site's alleles. Filed separately from
    /// the span guard: a distinct cause, not a coordinate-flank failure.
    pub assignments_skipped_oversized_neighbor: u64,
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

    // ALT event counts per left-normalized site. Normalization can move an event into
    // a column the loop has already passed, so sites are decided after it closes.
    let mut sites: HashMap<(Vec<u8>, u32), indel::IndelSite> = HashMap::new();
    // Per-column (depth, fwd, rev). This -- NOT a tally of observations -- is the VAF
    // denominator: a read carrying an event reports Indel::None at the normalized
    // column as well as its event at the anchor column, so counting observations
    // would count that read twice and halve every real indel's apparent frequency.
    let mut depths: HashMap<(Vec<u8>, u32), (u32, u32, u32)> = HashMap::new();

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
        let mut obs: Vec<(Vec<u8>, usize, usize, u8)> = Vec::new(); // (qname, qpos, allele, qual)

        // Indel collection MUST happen here, before the `non_ref() == 0` skip below:
        // that gate is about non-reference SUBSTITUTIONS, and a column with a
        // perfectly clean base composition can still carry indel events. Collecting
        // after it would silently lose every indel at a substitution-clean column.
        let mut col_depth = 0u32;
        let mut col_fwd = 0u32;
        let mut col_rev = 0u32;
        // (allele, reverse, read's own POS, read's own CIGAR). The read's full
        // CIGAR travels with the event so the ALT-vote gate can be checked AFTER
        // normalization, once we know which column the event actually folds into
        // (see CRITICAL 2 / FIX 2), via `indel::match_run_len`.
        let mut col_events: Vec<(indel::Allele, bool, i64, CigarString)> = Vec::new();

        for a in col.alignments() {
            let rec = a.record();
            if rec.is_secondary() || rec.is_supplementary() || rec.is_unmapped() {
                continue;
            }

            if iopts.enabled && !a.is_refskip() {
                if let Some(q) = a.qpos() {
                    col_depth += 1;
                    if rec.is_reverse() { col_rev += 1 } else { col_fwd += 1 }
                    // Decide the event's allele first, WITHOUT touching the CIGAR;
                    // `rec.cigar()` re-parses and allocates the whole CIGAR from raw
                    // bytes on every call, so it must be hoisted to a single call
                    // site reached at most once per alignment (only when this
                    // alignment actually carries a correctable event), not called
                    // separately from each of the two (mutually exclusive) arms.
                    let ev_allele = match a.indel() {
                        bam::pileup::Indel::Ins(n) if n < iopts.max_len => {
                            let s: Vec<u8> =
                                (0..n as usize).map(|i| rec.seq()[q + 1 + i]).collect();
                            Some(indel::Allele::Ins(s))
                        }
                        bam::pileup::Indel::Del(n) if n < iopts.max_len => {
                            Some(indel::Allele::Del(n))
                        }
                        // Indel::None contributes depth only; an indel at or above
                        // max_len contributes depth and is never corrected.
                        _ => None,
                    };
                    if let Some(allele) = ev_allele {
                        col_events.push((allele, rec.is_reverse(), rec.pos(), rec.cigar().take()));
                        stats.indel_events_examined += 1;
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
            obs.push((rec.qname().to_vec(), qpos, allele, qual));
        }

        if iopts.enabled {
            depths.insert((name.to_vec(), refpos as u32), (col_depth, col_fwd, col_rev));
            for (allele, reverse, rpos, cigar) in col_events {
                let (norm_pos, norm_allele) =
                    indel::normalize_left(refseq, refpos as u32, &allele);
                // `normalize_left` walks the reference with no knowledge of read
                // boundaries: a read whose alignment starts (or ends) inside a
                // repeat can normalize an event to a column outside its own
                // aligned span -- or even to a column its span covers only through
                // this SAME read's own `D`/`N`, where depth was never incremented
                // for this read either (depth only counts
                // `!is_refskip() && qpos().is_some()`; see FIX 2). Either way,
                // counting the ALT vote would let alt_total exceed depth at that
                // column -- see CRITICAL 2. The correct test is whether the read
                // has an actual ALIGNED QUERY BASE at norm_pos, i.e. a nonzero
                // contiguous M/=/X run starting there.
                if indel::match_run_len(rpos, &cigar, norm_pos) == 0 {
                    stats.indel_events_out_of_span += 1;
                    continue;
                }
                sites
                    .entry((name.to_vec(), norm_pos))
                    .or_default()
                    .add(norm_allele, reverse);
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

        let mut changed = false;
        for (qname, _qpos, allele, qual) in obs {
            let corr = map_allele(allele, qual, &model);
            if corr != allele {
                corrections
                    .entry(qname)
                    .or_default()
                    .subs
                    .push((refpos as u32, idx_to_base(corr)));
                stats.bases_modified += 1;
                stats.substitution_matrix[allele][corr] += 1;
                changed = true;
            }
        }
        if changed {
            stats.columns_corrected += 1;
        }
    }

    // ---- deferred indel site decisions ----
    let mut models: IndelModels = HashMap::new();
    if iopts.enabled {
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
                continue; // nothing spans this site; nothing to decide
            }
            site.depth = depth;
            site.col_fwd = fwd;
            site.col_rev = rev;

            stats.indel_sites_examined += 1;
            let model = indel::fit_indel_site(
                &site, refseq, norm_pos, min_strand, strand_bias_p, homoplasmic_vaf, iopts,
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
        for edits in corrections.values_mut() {
            edits.subs.sort_by_key(|&(p, _)| p);
        }
    }

    Ok((corrections, models, stats))
}

fn apply_corrections(
    in_bam: &Path,
    out_bam: &Path,
    corrections: &Corrections,
    refs: &HashMap<Vec<u8>, Vec<u8>>,
    models: &IndelModels,
    iopts: &indel::IndelOpts,
    stats: &mut DenoiseStats,
) -> Result<(u64, u64)> {
    let mut reader = Reader::from_path(in_bam)
        .with_context(|| format!("cannot open BAM {in_bam:?}"))?;
    let header_view = reader.header().to_owned();
    let header = bam::Header::from_template(reader.header());
    let mut writer = bam::Writer::from_path(out_bam, &header, bam::Format::Bam)
        .with_context(|| format!("cannot create BAM {out_bam:?}"))?;

    let mut reads_processed = 0u64;
    let mut reads_modified = 0u64;
    // Contigs referenced by a correctable, edit-bearing record but absent from
    // `refs`, and how many such records were affected. Tracked so a missing
    // reference contig is loud (one warning) rather than silently splicing
    // filler bases once indel edits exist (see FIX 1).
    let mut missing_ref_contigs: std::collections::BTreeSet<Vec<u8>> = std::collections::BTreeSet::new();
    let mut missing_ref_reads = 0u64;

    for result in reader.records() {
        let mut rec = result.context("read BAM record")?;
        // Only primary mapped reads are candidates for correction.
        let correctable =
            !(rec.is_secondary() || rec.is_supplementary() || rec.is_unmapped());
        if correctable {
            reads_processed += 1;
            // A mapped-but-contig-less record (tid == -1) is a degenerate BAM
            // state. Casting -1 to u32 before calling `tid2name` is an unchecked
            // pointer offset in rust-htslib (see IMPORTANT 4); guard it exactly
            // like the missing-reference-contig case below and skip straight to
            // the unconditional write at the bottom of the loop.
            if rec.tid() >= 0 {
                let contig = header_view.tid2name(rec.tid() as u32).to_vec();
                match refs.get(&contig) {
                    Some(refseq) => {
                        let mut edits = corrections.get(rec.qname()).cloned().unwrap_or_default();

                        if iopts.enabled {
                            // `end_pos()` (exclusive alignment end on the reference) is only
                            // defined on `CigarStringView`, which carries the record's POS;
                            // it must be read before `.take()` discards that view down to a
                            // plain `CigarString`.
                            let view = rec.cigar();
                            let rstart = rec.pos();
                            let rend = view.end_pos();
                            let cigar = view.take();
                            let seq = rec.seq().as_bytes();
                            let evs = indel::read_events(rstart, &cigar, &seq, refseq, iopts.max_len);

                            // What this read carries at each normalized site -- the allele
                            // AND the CIGAR anchor where that allele actually lives, since
                            // `rewrite_read` keys its Ins/Del lookups on the anchor, not the
                            // normalized site (see CRITICAL 1) -- and where it has an indel
                            // too large to represent.
                            let mut carried: HashMap<u32, (indel::Allele, u32)> = HashMap::new();
                            let mut blocked: Vec<i64> = Vec::new();
                            for e in &evs {
                                match &e.allele {
                                    Some(a) => { carried.insert(e.norm_pos, (a.clone(), e.anchor)); }
                                    None => blocked.push(e.norm_pos as i64),
                                }
                            }

                            let flank = iopts.flank as i64;
                            // A gained deletion consumes up to max_len reference bases to the
                            // right, so the right-hand requirement is larger than the left.
                            let reach = flank + iopts.max_len as i64;

                            if let Some(list) = models.get(&contig) {
                                let lo = list.partition_point(|(p, _)| (*p as i64) < rstart);
                                for (pos, model) in &list[lo..] {
                                    let s = *pos as i64;
                                    // `rend` is EXCLUSIVE, so `s == rend` is a site the
                                    // read does not span at all; admitting it here only
                                    // inflates the span-guard counter below with sites
                                    // that were never real candidates for this read.
                                    if s >= rend {
                                        break;
                                    }
                                    let (observed, anchor) = carried
                                        .get(pos)
                                        .cloned()
                                        .unwrap_or((indel::Allele::Ref, *pos));
                                    let target = match indel::assign_allele(&observed, model, iopts) {
                                        indel::Assignment::Keep => continue,
                                        indel::Assignment::Refused => {
                                            stats.assignments_skipped_unsupported_transition += 1;
                                            continue;
                                        }
                                        indel::Assignment::Move(t) => t,
                                    };
                                    // Left-normalization only ever moves an event LEFT, so the
                                    // read's own CIGAR anchor (where an observed indel actually
                                    // lives) is always >= the normalized site position. The
                                    // rewrite walk needs plain match context flanking whichever
                                    // coordinate it will actually touch.
                                    let far = (anchor as i64).max(s);
                                    if s - flank < rstart || far + reach > rend {
                                        stats.assignments_skipped_span_guard += 1;
                                        continue;
                                    }
                                    // FIX 3: the check above bounds only the read's
                                    // OVERALL span, which can pass even when the
                                    // read's own alignment breaks up into several
                                    // Match/=/X runs -- e.g. an unrelated `D`/`N`
                                    // shortly after this site. A GAINED deletion
                                    // (Ref -> Del(m)) is spliced into `rewrite_read`'s
                                    // Match arm at the normalized site, so it and its
                                    // right-hand flank must both fit inside that ONE
                                    // contiguous run, or the walk truncates it to the
                                    // wrong length (debug_assert in debug, silent
                                    // truncation in release).
                                    if observed == indel::Allele::Ref {
                                        if let indel::Allele::Del(m) = &target {
                                            let need = *m as i64 + flank;
                                            if (indel::match_run_len(rstart, &cigar, s as u32) as i64)
                                                < need
                                            {
                                                stats.assignments_skipped_span_guard += 1;
                                                continue;
                                            }
                                        }
                                    }
                                    // An oversized indel nearby means this read's alignment in
                                    // this window cannot be represented by the site's alleles.
                                    if blocked.iter().any(|b| (b - s).abs() <= reach) {
                                        stats.assignments_skipped_oversized_neighbor += 1;
                                        continue;
                                    }
                                    match (&observed, &target) {
                                        (indel::Allele::Ref, indel::Allele::Ins(s2)) => {
                                            stats.reads_reassigned_ref_to_ins += 1;
                                            stats.indel_bases_inserted += s2.len() as u64;
                                        }
                                        (indel::Allele::Ref, indel::Allele::Del(n)) => {
                                            stats.reads_reassigned_ref_to_del += 1;
                                            stats.indel_bases_removed += *n as u64;
                                        }
                                        (indel::Allele::Ins(from_s), indel::Allele::Ref) => {
                                            // Reverting a spurious insertion removes
                                            // the bases it added -- the single most
                                            // common correction this feature makes.
                                            stats.reads_reassigned_indel_to_ref += 1;
                                            stats.indel_bases_removed += from_s.len() as u64;
                                        }
                                        (indel::Allele::Del(from_n), indel::Allele::Ref) => {
                                            // Reverting a spurious deletion restores
                                            // the reference bases it dropped.
                                            stats.reads_reassigned_indel_to_ref += 1;
                                            stats.indel_bases_inserted += *from_n as u64;
                                        }
                                        (indel::Allele::Ins(from_s), indel::Allele::Ins(to_s)) => {
                                            stats.reads_reassigned_indel_to_indel += 1;
                                            if to_s.len() > from_s.len() {
                                                stats.indel_bases_inserted +=
                                                    (to_s.len() - from_s.len()) as u64;
                                            } else if from_s.len() > to_s.len() {
                                                stats.indel_bases_removed +=
                                                    (from_s.len() - to_s.len()) as u64;
                                            }
                                        }
                                        (indel::Allele::Del(from_n), indel::Allele::Del(to_n)) => {
                                            stats.reads_reassigned_indel_to_indel += 1;
                                            if to_n < from_n {
                                                stats.indel_bases_removed += (from_n - to_n) as u64;
                                            } else if to_n > from_n {
                                                stats.indel_bases_inserted += (to_n - from_n) as u64;
                                            }
                                        }
                                        _ => stats.reads_reassigned_indel_to_indel += 1,
                                    }
                                    let l = model.context_len as usize;
                                    if stats.reassignments_by_hp_length.len() <= l {
                                        stats.reassignments_by_hp_length.resize(l + 1, 0);
                                    }
                                    stats.reassignments_by_hp_length[l] += 1;

                                    // Direction determines which coordinate `rewrite_read`
                                    // will actually look this edit up by (see CRITICAL 1):
                                    // reverting/replacing an indel the read itself carries is
                                    // found via `last_ref`, which is the CIGAR anchor; gaining
                                    // an indel where the read has plain match context is found
                                    // by walking that Match run at the normalized site.
                                    let ref_pos = match &observed {
                                        indel::Allele::Ref => *pos,
                                        _ => anchor,
                                    };
                                    edits.indels.push(rewrite::IndelEdit {
                                        ref_pos,
                                        from: observed,
                                        to: target,
                                    });
                                }
                            }
                            edits.indels.sort_by_key(|e| e.ref_pos);
                        }

                        if !edits.is_empty() {
                            // Extract owned parts, rewrite, rebuild (aux preserved by set()).
                            let qname = rec.qname().to_vec();
                            let cigar = rec.cigar().take();
                            let seq = rec.seq().as_bytes();
                            let qual = rec.qual().to_vec();
                            let out =
                                rewrite_read(rec.pos(), &cigar, &seq, &qual, &edits, refseq);
                            // A non-empty edit list does not guarantee `rewrite_read`
                            // actually changed anything -- CRITICAL 1's bug is exactly a
                            // case where it silently didn't. Compare the real output.
                            let changed = out.seq != seq || out.cigar != cigar;
                            rec.set(&qname, Some(&out.cigar), &out.seq, &out.qual);
                            if out.structure_changed {
                                strip_stale_tags(&mut rec);
                            }
                            if changed {
                                reads_modified += 1;
                            }
                        }
                    }
                    None => {
                        // Contig not in the supplied reference: do not rewrite
                        // this record. Writing it through unchanged is safe;
                        // fabricating a reference slice is not (it would let
                        // the deletion-reversion path splice `N` filler into
                        // real reads).
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
    Ok((reads_processed, reads_modified))
}

pub fn start(
    input: &PathBuf,
    output: &PathBuf,
    reference: &PathBuf,
    data_type: &str,
    vaf: f64,
    min_strand: u32,
    strand_bias_p: f64,
    homoplasmic_vaf: f64,
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

    // The `indels` parameter arrives in Task 9; until then indel correction stays
    // disabled here, matching `IndelOpts::default()`.
    let default_iopts = indel::IndelOpts::default();
    let (corrections, models, mut stats) = compute_corrections(
        input, &refs, min_strand, vaf, strand_bias_p, homoplasmic_vaf, &default_iopts,
    )?;
    let (reads_processed, reads_modified) = apply_corrections(
        input, output, &corrections, &refs, &models, &default_iopts, &mut stats,
    )?;
    stats.reads_processed = reads_processed;
    stats.reads_modified = reads_modified;

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
        // Q0 nominally means P(error) = 1, which would make P(obs | truth = obs) zero
        // and invert every base. With 4 alleles the most error a base can carry is
        // uniform-over-4, i.e. 0.75, so everything at or below Q1 clamps there.
        assert!((phred_to_e(0) - E_MAX).abs() < 1e-12);
        assert!((phred_to_e(1) - E_MAX).abs() < 1e-12);
        for q in 0..=93u8 {
            assert!(phred_to_e(q) <= E_MAX, "q={q} exceeds the clamp");
        }
    }

    #[test]
    fn zero_quality_observations_are_never_inverted() {
        // ref A(0) at 60%, alt G(2) at 40%, both strand-balanced, ALL bases at Q0 --
        // the shape pbsim's errhmm models emit (quality string of all '!').
        let mut obs = vec![];
        for i in 0..60 { obs.push((0usize, q(0), i % 2 == 0)); }
        for i in 0..40 { obs.push((2usize, q(0), i % 2 == 0)); }
        let c = column(0, &obs);
        let m = fit_site(&c, 2, 0.01, SB_P, HOM_VAF);

        // Before the clamp, like(k == observed) was 0 while like(k != observed) was
        // e/3 > 0, so an A read was rewritten to some OTHER base -- corrupting the
        // read instead of correcting it. A ref observation must stay ref.
        assert_eq!(map_allele(0, q(0), &m), 0);
        // At the clamp the likelihood is flat (1 - e == e/3 == 0.25), so a Q0 base
        // carries no evidence and correction falls back to the frequency prior:
        // every observation collapses to the most frequent kept allele (here ref).
        // That is degenerate-input behavior, but it never invents a third base.
        for observed in 0..4 {
            assert_eq!(map_allele(observed, q(0), &m), 0, "observed={observed}");
        }
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
        // ref A(0) 40fwd/40rev; alt G(2) 18fwd/2rev. min_strand=2 is satisfied on BOTH
        // strands, so the count gate alone lets this through -- but against the column's
        // own composition (58/100 fwd) an 18/2 split has p = 4.3e-3 < 0.01.
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
        // ref A(0) 2fwd/18rev; alt G(2) 75fwd/5rev -> alt is 80% of the column. Against
        // the column composition (77/100 fwd) that split gives p = 1.2e-4, but an allele
        // in 80% of reads cannot be a single-strand artifact, so it is kept untested.
        let mut obs = vec![];
        for _ in 0..2 { obs.push((0usize, q(30), false)); }
        for _ in 0..18 { obs.push((0usize, q(30), true)); }
        for _ in 0..75 { obs.push((2usize, q(30), false)); }
        for _ in 0..5 { obs.push((2usize, q(30), true)); }
        let c = column(0, &obs);

        let m = fit_site(&c, 2, 0.01, SB_P, 0.7);
        assert!(m.kept[2], "near-homoplasmic alt must survive the strand-bias test");
        assert_eq!(m.strand_biased_alts, 0);
        assert_eq!(map_allele(2, q(30), &m), 2);

        // Raising the exemption bar above the allele's frequency re-enables the test,
        // which then rejects it -- confirming the exemption is what saved it above.
        let m_strict = fit_site(&c, 2, 0.01, SB_P, 0.95);
        assert!(!m_strict.kept[2]);
        assert_eq!(m_strict.strand_biased_alts, 1);
    }

    #[test]
    fn empty_column_falls_back_to_reference_only() {
        // No observations at all: expected-fwd-fraction fallback must not panic and the
        // model must degenerate to reference-only.
        let c = column(0, &[]);
        let m = fit_site(&c, 2, 0.01, SB_P, HOM_VAF);
        assert!(m.kept[0]);
        assert_eq!(m.kept_alt_count(Some(0)), 0);
        assert!((m.freq[0] - 1.0).abs() < 1e-12);
    }

    // Build a tiny coordinate-sorted BAM. Reads are full-length matches at `pos`.
    // Each read: (qname, pos, seq bytes, reverse). Quality fixed at Q30.
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
        // Reference AAAAAA at chrM. Ten reads all "AAAAAA" except one with a G at pos 2.
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
            // Record::new() initializes as unmapped (see rust-htslib set_unmapped() in
            // Record::new()); set_tid/set_pos alone do not clear that flag, so we must
            // explicitly mark this record mapped to make it a candidate primary read.
            rec.unset_unmapped();
            rec.push_aux(b"XY", Aux::I32(7)).unwrap();
            w.write(&rec).unwrap();
        }

        let mut corr: Corrections = HashMap::new();
        corr.insert(
            b"rerr".to_vec(),
            rewrite::ReadEdits { subs: vec![(2u32, b'A')], indels: vec![] },
        );
        let mut refs = HashMap::new();
        refs.insert(b"chrM".to_vec(), b"AAAAAA".to_vec());
        let models = IndelModels::new();
        let mut stats = DenoiseStats::default();
        let (proc, modified) = apply_corrections(
            &inb, &outb, &corr, &refs, &models, &indel::IndelOpts::default(), &mut stats,
        )
        .unwrap();
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

        start(&inb, &outb, &reff, "pacbio", 0.01, 2, SB_P, HOM_VAF, None).unwrap();
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

        start(&inb, &outb, &reff, "ont-r10", 0.01, 2, SB_P, HOM_VAF, Some(&statsf)).unwrap();

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
        // serde_json::to_string_pretty (used by `start`, per spec) inserts a space
        // after the colon, unlike compact serialization; the brief's literal
        // "\"bases_modified\":1" assumes compact output and never matches pretty
        // JSON, so this assertion matches the actual (and required) pretty format.
        assert!(json.contains("\"bases_modified\": 1"));
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn corrections_are_keyed_by_reference_position() {
        // The single G error sits at reference position 2 (reads start at POS 0).
        // Recording it in reference coordinates is what lets Task 7 mix
        // length-changing indel edits into the same per-read payload (Task 8).
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
        // Reads aligned at POS 10 (not 0): the erroneous base sits at read/query
        // offset 2, which maps to REFERENCE position 12 (10 + 2). At POS 0 those
        // two numbers coincide (2 == 2), so a fixture like
        // `corrections_are_keyed_by_reference_position` above cannot distinguish
        // a reference-keyed implementation from a query-offset-keyed one -- both
        // would record the key as 2. With a nonzero POS the two coordinate
        // systems disagree (12 != 2), so asserting the recorded key is 12 (not
        // 2) actually exercises the reference-coordinate semantics.
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
        // The recorded key must be the REFERENCE position (12), not the query
        // offset (2) at which the erroneous base actually occurs in the read.
        assert_eq!(edits.subs, vec![(12u32, b'A')]);
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn apply_corrections_missing_contig_writes_through_unchanged() {
        // A record on a contig absent from `refs`, with a correction registered
        // for its qname, must NOT be rewritten: fabricating a reference slice
        // for a missing contig is unsafe once indel edits exist (the
        // deletion-reversion path would splice `N` filler into real reads).
        let dir = std::env::temp_dir().join("himito_denoise_missing_contig");
        std::fs::create_dir_all(&dir).unwrap();
        let inb = dir.join("in.bam");
        let outb = dir.join("out.bam");

        write_test_bam(&inb, "chrM", 6, &[("rerr", 0, b"AAGAAA", false)]);

        let mut corr: Corrections = HashMap::new();
        corr.insert(
            b"rerr".to_vec(),
            rewrite::ReadEdits { subs: vec![(2u32, b'A')], indels: vec![] },
        );
        // Deliberately do NOT insert "chrM" -- simulates a BAM whose contig is
        // missing from the supplied reference FASTA.
        let refs: HashMap<Vec<u8>, Vec<u8>> = HashMap::new();

        let models = IndelModels::new();
        let mut stats = DenoiseStats::default();
        let (proc, modified) = apply_corrections(
            &inb, &outb, &corr, &refs, &models, &indel::IndelOpts::default(), &mut stats,
        )
        .unwrap();
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

    // Reference shared by the indel wiring tests. Index map:
    //   0-7   ACGTACGT
    //   8-15  AAAAAAAA   (an 8bp homopolymer run)
    //   16-28 CGTACGTACGTAC
    // Long enough that the shipped default --indel-flank 5 is satisfiable at the
    // sites used below (normalized positions 7 and 19).
    const REF29: &[u8] = b"ACGTACGTAAAAAAAACGTACGTACGTAC";

    // Build a BAM whose reads carry explicit CIGARs, so indel-bearing reads can be
    // constructed directly. Each read: (qname, pos, seq, cigar ops, reverse).
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

    const CLEAN_SEQ: &[u8] = b"ACGTACGTAAAAAAAACGTACGTACGTAC"; // == REF29, 29bp
    // An extra A inside the run: 16M 1I 13M. Aligner anchors at index 15;
    // normalization pulls the site to index 7.
    const HP_INS_SEQ: &[u8] = b"ACGTACGTAAAAAAAAACGTACGTACGTAC"; // 30bp
    // A spurious G after index 18, in unique sequence: 19M 1I 10M.
    const UNIQ_INS_SEQ: &[u8] = b"ACGTACGTAAAAAAAACGTGACGTACGTAC"; // 30bp

    fn clean_cigar() -> Vec<Cigar> { vec![Cigar::Match(29)] }
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
        // The double-count trap: a read carrying the insertion reports Indel::None at
        // the normalized column (7) AND its insertion at the anchor column (15).
        // Depth at site 7 must be 10 -- the number of spanning reads -- not 19.
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
        // 9 of 10 reads carry it. A double-counted denominator would give ~0.47.
        assert!(
            (f_ins - 0.9).abs() < 1e-9,
            "insertion frequency should be 0.9, got {f_ins} (depth double-counted?)"
        );
        assert_eq!(model.context_len, 8, "site 7 sits in an 8bp homopolymer run");
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn models_are_sorted_by_position_within_each_contig() {
        // The record walk binary-searches this list; unsorted input silently drops
        // sites from every read's span.
        let dir = std::env::temp_dir().join("himito_denoise_indel_sorted");
        std::fs::create_dir_all(&dir).unwrap();
        let bam = dir.join("in.bam");
        let both: Vec<Cigar> = vec![
            Cigar::Match(16), Cigar::Ins(1), Cigar::Match(4), Cigar::Ins(1), Cigar::Match(9),
        ];
        // Extra A in the run (site 7) AND a spurious G after index 19 (site 19).
        // 16M consumes ref[0..16] (query "ACGTACGT"+"AAAAAAAA", 16bp) then Ins(1)
        // inserts the extra A (anchor 15 -> normalizes to 7); 4M consumes ref[16..20]
        // = "CGTA" (query); Ins(1) inserts the spurious G right after ref index 19
        // (anchor 19, no normalization: refseq[19] == 'A' != 'G'); 9M consumes
        // ref[20..29] = "CGTACGTAC". Query length = 16+1+4+1+9 = 31.
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
        // Assert the exact expected positions rather than checking ascending order
        // pairwise (the old form: `list.windows(2).all(|w| w[0].0 <= w[1].0)`):
        // with only two elements, an unsorted 2-element list still satisfies
        // pairwise ascending order about half the time by pure chance (whichever
        // of the two orderings HashMap iteration happens to produce), so that
        // check would pass about half the time even with the `sort_by_key`
        // deleted.
        let positions: Vec<u32> = list.iter().map(|(p, _)| *p).collect();
        assert_eq!(positions, vec![7, 19], "sites must be exactly [7, 19], in that order");
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn lone_spurious_insertion_is_removed_from_a_substitution_clean_column() {
        // Every aligned base matches the reference, so `counts.non_ref() == 0` at
        // every column and the SNV path skips them all. The lone spurious insertion
        // must STILL be collected and corrected -- this is the early-continue trap.
        // It is also the REF-only-site case: no alt allele clears candidacy, so the
        // site keeps REF alone and the stray read must still be reverted.
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
        // FIX 4: reverting a 1bp insertion removes 1 base from the read.
        assert_eq!(stats.indel_bases_removed, 1);
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn homoplasmic_indel_survives_and_repairs_the_read_that_dropped_it() {
        // 9 of 10 reads carry a 1bp insertion inside the 8bp A-run, so eps is large
        // enough for the repair direction to fire (in unique sequence eps is 0.01 and
        // a 10% REF allele correctly would NOT flip).
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
        // Same 90%-insertion site (normalized position 7), but the read that dropped
        // the insertion is aligned 0..14 -- it ends before 7 + flank(5) + max_len(5)
        // = 17, so the rewrite walk would lack the match context it needs.
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
        // IMPORTANT 3: the original version of this fixture placed the site in
        // UNIQUE sequence (context_len = 1, eps = 0.01). There, for a Ref-observing
        // read MAP already returns Keep on the arithmetic alone (0.1 * 0.99 = 0.099
        // beats 0.9 * 0.01 = 0.009), so the `blocked` mask is never even reached --
        // the assertion held identically with the mask deleted.
        //
        // This version places the site inside the 8bp A-run (site 7), where eps is
        // large enough that a Ref-observing read WOULD genuinely be repaired into
        // the kept 90% insertion without the mask: eps(8) ~= 0.171, and
        // 0.9 * 0.171 = 0.154 beats 0.1 * 0.829 = 0.083 (identical arithmetic to
        // `near_homoplasmic_site_repairs_reads_that_dropped_the_indel` in
        // indel.rs). `rbig` observes plain Ref at site 7 (its own event is
        // elsewhere), so only the oversized-neighbor mask -- not the assignment
        // math -- can be what leaves it untouched.
        let dir = std::env::temp_dir().join("himito_denoise_indel_large");
        std::fs::create_dir_all(&dir).unwrap();
        let bam = dir.join("in.bam");
        let outb = dir.join("out.bam");
        // ref[0..9] + 6 spurious G's (oversized, anchor = 8, inside the run) +
        // ref[9..29]. `rbig` carries no event of its own at site 7.
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

    // ---- CRITICAL 1 regression: normalized position vs CIGAR anchor mismatch ----

    #[test]
    fn spurious_homopolymer_insertion_is_reverted_using_its_own_cigar_anchor() {
        // 9 clean reads plus one (`rspur`) carrying a lone extra A inside the 8bp
        // A-run, anchored by the aligner at reference index 15 (`HP_INS_SEQ`),
        // which left-normalizes to site 7. The aligner's own CIGAR anchor (15)
        // and the normalized site (7) genuinely differ here -- unlike every other
        // pre-existing rewrite fixture, which sits in unique sequence where the
        // two coincide.
        //
        // Before CRITICAL 1's fix, `apply_corrections` emitted the `IndelEdit`
        // keyed at the NORMALIZED position (7). `rewrite_read` looks an `Ins`
        // reversion up via `last_ref` -- the read's own CIGAR anchor (15) -- so
        // the lookup at 7 found nothing and the `Ins` op was copied verbatim: a
        // silent no-op that still incremented `reads_reassigned_indel_to_ref`.
        // This test asserts the actual SEQ/CIGAR, not just the counter.
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
        // FIX 4: this revert removes the 1 spurious base from the read; before
        // FIX 4 `indel_bases_removed` stayed 0 for exactly this transition.
        assert_eq!(stats.indel_bases_removed, 1);
        std::fs::remove_dir_all(&dir).ok();
    }

    // ---- CRITICAL 2 regression: depth/ALT sourced inconsistently across span ----

    #[test]
    fn indel_event_normalizing_outside_its_own_span_is_not_fabricated_into_a_site() {
        // Four reads start at POS 10 -- inside the 8bp A-run (ref idx 8..=15) --
        // and each carries an extra A anchored at ref idx 10 (inside their own
        // aligned span). `normalize_left` walks the reference with no knowledge
        // of read boundaries and pulls the event back to site 7, which lies
        // BEFORE these reads' own start (10). Without the span check, this
        // fabricates ALT support at site 7 for reads that never covered it: the
        // per-column depth at site 7 (contributed only by the one fully-spanning
        // read below) stays small while alt_total grows unbounded, letting
        // alt_total exceed depth (the `debug_assert!` CRITICAL 2's probe tripped)
        // and, in release, fabricating a homoplasmic insertion that would splice
        // into the one genuinely reference-matching read.
        let dir = std::env::temp_dir().join("himito_denoise_indel_out_of_span");
        std::fs::create_dir_all(&dir).unwrap();
        let bam = dir.join("in.bam");
        let outb = dir.join("out.bam");

        // seq = ref[10] + inserted 'A' + ref[11..29]; anchor = 10, inside this
        // read's own span [10, 29) -- but normalizes to site 7, outside it.
        let mut seq: Vec<u8> = Vec::new();
        seq.push(REF29[10]);
        seq.push(b'A');
        seq.extend_from_slice(&REF29[11..29]);
        let cigar = vec![Cigar::Match(1), Cigar::Ins(1), Cigar::Match(18)];

        let mut reads: Vec<(&str, i64, &[u8], Vec<Cigar>, bool)> = vec![];
        // One read spanning the WHOLE contig, clean, reference-matching -- this
        // is the read that must come out unchanged. Written first: the pileup
        // requires coordinate-sorted input, and this read's POS (0) precedes
        // the others' (10).
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
        // This fixture's only would-be evidence for site 7 is the four dropped
        // votes above, so with them correctly rejected "chrM" has no entry in
        // `models` at all -- assert that absence explicitly instead of nesting
        // the check inside an `if let Some(list)` that never runs (and so never
        // exercises this assertion) on this fixture.
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

    // ---- FIX 2 regression: a read's own D/N covering the normalized column ----

    #[test]
    fn read_deleting_the_normalized_column_does_not_cast_a_phantom_alt_vote() {
        // Depth is only incremented where a read has qpos().is_some() -- never
        // inside that SAME read's own D/N. `rdel` deletes reference column 7
        // (Match(7) + Del(1)) and, separately, carries a spurious insertion
        // inside the downstream A-run (anchor 15) that left-normalizes to
        // that very column 7. `rdel`'s own [pos, end_pos()) span covers
        // column 7 (via the deletion), so the OLD span check would have
        // accepted this as an ALT vote there with no matching depth --
        // alt_total could exceed depth in miniature, CRITICAL 2's failure
        // mode reappearing. The new contiguous-match-run check must reject
        // it, because `rdel` has no aligned M/=/X base at column 7 at all.
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

    // ---- FIX 3 regression: a gained deletion must fit inside ONE match run ----

    #[test]
    fn gained_deletion_is_skipped_when_its_own_match_run_is_too_short() {
        // The span/flank guard bounds only the read's OVERALL alignment span.
        // Read `rguard`'s alignment is 0..33 overall -- comfortably inside
        // flank(5)+reach(10) of the decided site at 18 -- but its FIRST
        // Match run ends at ref pos 20, only 2 bases past the site, because
        // of an unrelated RefSkip further along. A gained Del(2) at 18 needs
        // a contiguous match run of at least m(2)+flank(5)=7 bases from the
        // site; this read only has 2, so it must be skipped, not truncated
        // (the old guard alone would have let `rewrite_read`'s Match arm
        // walk `k + m` past its op length).
        let dir = std::env::temp_dir().join("himito_denoise_indel_shortrun");
        std::fs::create_dir_all(&dir).unwrap();
        let inb = dir.join("in.bam");
        let outb = dir.join("out.bam");

        let seq = vec![b'A'; 30];
        let cigar = vec![Cigar::Match(20), Cigar::RefSkip(3), Cigar::Match(10)];
        write_cigar_bam(&inb, "chrM", 40, &[("rguard", 0, seq.as_slice(), cigar, false)]);

        let mut refs = HashMap::new();
        refs.insert(b"chrM".to_vec(), vec![b'A'; 40]);

        // Hand-built decided site at position 18: Del(2) dominates 90/10 over
        // Ref, with eps large enough that MAP assignment genuinely wants to
        // move a Ref-observing read to Del(2) (same arithmetic shape as
        // `near_homoplasmic_site_repairs_reads_that_dropped_the_indel` in
        // indel.rs). Hand-building this model bypasses the pileup/candidacy
        // machinery entirely, which is irrelevant to what this test checks.
        let model = indel::IndelSiteModel {
            kept: vec![(indel::Allele::Ref, 0.1), (indel::Allele::Del(2), 0.9)],
            eps: 0.15,
            context_len: 1,
            strand_biased: 0,
            obs_vaf: vec![(indel::Allele::Ref, 0.1), (indel::Allele::Del(2), 0.9)],
            strand_rejected: vec![],
        };
        let mut models = IndelModels::new();
        models.entry(b"chrM".to_vec()).or_default().push((18u32, model));

        let corr: Corrections = HashMap::new();
        let mut stats = DenoiseStats::default();
        let (proc, modified) = apply_corrections(
            &inb, &outb, &corr, &refs, &models, &indels_on(), &mut stats,
        )
        .unwrap();
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

    // ---- FIX 4 regression: revert-direction base counts ----

    #[test]
    fn lone_spurious_deletion_is_reverted_and_counted_as_bases_inserted() {
        // `Del(n)->Ref` restores n bases to the read but was not being counted
        // in `indel_bases_inserted`. Mirrors
        // `lone_spurious_insertion_is_removed_from_a_substitution_clean_column`
        // but with a deletion: a single read carries a spurious 1bp deletion in
        // unique sequence; candidacy rejects it (a lone single-strand
        // observation fails min_strand), but the stray read must still be
        // repaired back to REF, and that repair must show up as 1 restored
        // (inserted) base.
        let dir = std::env::temp_dir().join("himito_denoise_indel_lone_del");
        std::fs::create_dir_all(&dir).unwrap();
        let bam = dir.join("in.bam");
        let outb = dir.join("out.bam");
        let mut reads: Vec<(&str, i64, &[u8], Vec<Cigar>, bool)> = vec![];
        for (i, n) in NAMES.iter().enumerate() {
            reads.push((n, 0, CLEAN_SEQ, clean_cigar(), i % 2 == 0));
        }
        // ref[0..7) + (skip ref[7]) + ref[8..29): a 1bp deletion anchored at 6,
        // in the unique prefix, far enough from both read ends to clear the
        // span/flank guard.
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
        // A mapped record with tid == -1 is a degenerate BAM state. Casting -1 to
        // u32 and calling `tid2name` on it is an unchecked pointer offset in
        // rust-htslib; it must be guarded and the record written through
        // unchanged, exactly like the missing-reference-contig path.
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
            // Deliberately do NOT set_tid: tid stays -1 (Record::new()'s default),
            // but clearing the unmapped flag makes this record "correctable".
            rec.unset_unmapped();
            w.write(&rec).unwrap();
        }

        let refs = indel_refs();
        let models = IndelModels::new();
        let mut stats = DenoiseStats::default();
        let corr: Corrections = HashMap::new();
        let (proc, modified) = apply_corrections(
            &inb, &outb, &corr, &refs, &models, &indel::IndelOpts::default(), &mut stats,
        )
        .unwrap(); // must not panic on an unchecked tid2name offset
        assert_eq!(proc, 1);
        assert_eq!(modified, 0);
        std::fs::remove_dir_all(&dir).ok();
    }
}
