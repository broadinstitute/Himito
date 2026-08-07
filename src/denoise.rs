use std::collections::HashMap;
use std::path::{Path, PathBuf};
use anyhow::{Context, Result};
use bio::io::fasta;
use log::{info, warn};
use rust_htslib::bam::{self, Read, Reader};

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
}

fn compute_corrections(
    bam_path: &Path,
    refs: &HashMap<Vec<u8>, Vec<u8>>,
    min_strand: u32,
    vaf: f64,
    strand_bias_p: f64,
    homoplasmic_vaf: f64,
) -> Result<(Corrections, DenoiseStats)> {
    let mut reader = Reader::from_path(bam_path)
        .with_context(|| format!("cannot open BAM {bam_path:?}"))?;
    // Own the tid->contig-name map so we can borrow the reader mutably for pileup.
    let header = reader.header().to_owned();

    let mut corrections: Corrections = HashMap::new();
    let mut stats = DenoiseStats::default();

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

        for a in col.alignments() {
            let rec = a.record();
            if rec.is_secondary() || rec.is_supplementary() || rec.is_unmapped() {
                continue;
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

    Ok((corrections, stats))
}

fn apply_corrections(
    in_bam: &Path,
    out_bam: &Path,
    corrections: &Corrections,
    refs: &HashMap<Vec<u8>, Vec<u8>>,
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
            if let Some(edits) = corrections.get(rec.qname()) {
                if !edits.is_empty() {
                    let contig = header_view.tid2name(rec.tid() as u32);
                    match refs.get(contig) {
                        Some(refseq) => {
                            // Extract owned parts, rewrite, rebuild (aux preserved by set()).
                            let qname = rec.qname().to_vec();
                            let cigar = rec.cigar().take();
                            let seq = rec.seq().as_bytes();
                            let qual = rec.qual().to_vec();
                            let out =
                                rewrite_read(rec.pos(), &cigar, &seq, &qual, edits, refseq);
                            rec.set(&qname, Some(&out.cigar), &out.seq, &out.qual);
                            if out.structure_changed {
                                strip_stale_tags(&mut rec);
                            }
                            reads_modified += 1;
                        }
                        None => {
                            // Contig not in the supplied reference: do not rewrite
                            // this record. Writing it through unchanged is safe;
                            // fabricating a reference slice is not (it would let
                            // the deletion-reversion path splice `N` filler into
                            // real reads).
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

    let (corrections, mut stats) =
        compute_corrections(input, &refs, min_strand, vaf, strand_bias_p, homoplasmic_vaf)?;
    let (reads_processed, reads_modified) = apply_corrections(input, output, &corrections, &refs)?;
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
        let (corr, stats) = compute_corrections(&bam, &refs, 2, 0.01, SB_P, HOM_VAF).unwrap();

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
        let (proc, modified) = apply_corrections(&inb, &outb, &corr, &refs).unwrap();
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
        let (corr, _) = compute_corrections(&bam, &refs, 2, 0.01, SB_P, HOM_VAF).unwrap();

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
        let (corr, _) = compute_corrections(&bam, &refs, 2, 0.01, SB_P, HOM_VAF).unwrap();

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

        let (proc, modified) = apply_corrections(&inb, &outb, &corr, &refs).unwrap();
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
}
