use std::collections::HashMap;
use std::path::{Path, PathBuf};
use anyhow::{Context, Result};
use rust_htslib::bam::{self, Read, Reader};

pub(crate) const NQ: usize = 94; // Phred qualities 0..=93

#[inline]
pub(crate) fn phred_to_e(q: u8) -> f64 {
    let q = q.min((NQ - 1) as u8);
    10f64.powf(-(q as f64) / 10.0)
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
}

pub(crate) struct SiteModel {
    freq: [f64; 4],
    kept: [bool; 4],
}

impl SiteModel {
    fn kept_alt_count(&self, ref_idx: Option<usize>) -> usize {
        (0..4)
            .filter(|&a| self.kept[a] && Some(a) != ref_idx)
            .count()
    }
}

fn fit_site(c: &ColumnCounts, min_strand: u32, vaf: f64) -> SiteModel {
    // Eligibility: reference always eligible; a non-ref allele needs >= min_strand
    // observations on EACH strand (baldur's strand-balance candidacy).
    let mut eligible = [false; 4];
    for a in 0..4 {
        eligible[a] = Some(a) == c.ref_idx
            || (c.strand[a][0] >= min_strand && c.strand[a][1] >= min_strand);
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

    SiteModel { freq, kept }
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

pub type Corrections = HashMap<Vec<u8>, Vec<(usize, u8)>>;

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
    pub substitution_matrix: [[u64; 4]; 4], // [from][to]
}

fn compute_corrections(
    bam_path: &Path,
    refs: &HashMap<Vec<u8>, Vec<u8>>,
    min_strand: u32,
    vaf: f64,
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

        let model = fit_site(&counts, min_strand, vaf);
        stats.alt_sites_preserved += model.kept_alt_count(ref_idx) as u64;

        let mut changed = false;
        for (qname, qpos, allele, qual) in obs {
            let corr = map_allele(allele, qual, &model);
            if corr != allele {
                corrections
                    .entry(qname)
                    .or_default()
                    .push((qpos, idx_to_base(corr)));
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

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::{self, Format, Header, Read, Reader, Record, Writer};
    use rust_htslib::bam::header::HeaderRecord;
    use rust_htslib::bam::record::{Cigar, CigarString};

    fn q(v: u8) -> u8 { v } // readability helper for quality values

    #[test]
    fn phred_to_e_matches_definition() {
        assert!((phred_to_e(10) - 0.1).abs() < 1e-12);
        assert!((phred_to_e(20) - 0.01).abs() < 1e-12);
        assert!((phred_to_e(0) - 1.0).abs() < 1e-12);
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
        let m = fit_site(&c, 2, 0.01);
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
        let m = fit_site(&c, 2, 0.01);
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
        let m = fit_site(&c, 2, 0.01);
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
        let m = fit_site(&c, 2, 0.05); // higher vaf threshold
        assert!(!m.kept[2]);
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
        let (corr, stats) = compute_corrections(&bam, &refs, 2, 0.01).unwrap();

        // The single G (not strand-balanced, VAF 10% but fwd-only) is corrected to A.
        assert_eq!(corr.get(&b"rerr".to_vec()), Some(&vec![(2usize, b'A')]));
        assert_eq!(stats.bases_modified, 1);
        assert_eq!(stats.substitution_matrix[2][0], 1); // G(2) -> A(0)
        assert!(stats.columns_processed >= 1);
        std::fs::remove_dir_all(&dir).ok();
    }
}
