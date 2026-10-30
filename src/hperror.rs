//! Homopolymer length confusion of a sample's own reads, for `lineage`'s tract-allele error rates.
//!
//! A Rust port of `lineage_eval/lineage_sim/measure_hp_errors.py` (plan 2026-09-27 T7): every reference
//! homopolymer run of length >= `min_len` is scored on every primary alignment with
//! [`crate::denoise_indel::hp_allele`], runs whose modal observed length is not the reference length are dropped
//! (the sample is not reference there), and the counts are written as `base ref_len obs_len strand n`.
//! The output must stay byte-identical to the Python script's.

use anyhow::{bail, Context, Result};
use rust_htslib::bam::{self, ext::BamRecordExtensions, Read as BamRead};
use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

/// `hp_tracts.CATALOGUE` run intervals, 1-based inclusive, by id.
pub const CATALOGUE: &[(&str, u32, u32)] = &[
    ("D310", 303, 309),
    ("HVS1", 16184, 16188),
    ("C573", 568, 573),
    ("T455", 452, 455),
    ("C960", 956, 960),
    ("C5899", 5895, 5899),
    ("A12418", 12418, 12425),
    ("CA514", 514, 523),
    ("R8281", 8272, 8289),
];

/// The tracts `lineage` already treats as variable; skipped unless `--tracts` says otherwise.
pub const DEFAULT_SKIP: &[&str] = &["D310", "HVS1", "CA514"];

/// Shortest run measured, as in `measure_hp_errors.py --min-len`.
pub const DEFAULT_MIN_LEN: u32 = 4;

/// One reference homopolymer run.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HpRun {
    /// 0-based base immediately before the run (the `a` of `hp_allele`).
    pub anchor: u32,
    /// 0-based exclusive end of the run (the `e` of `hp_allele`).
    pub end: u32,
    pub base: u8,
    pub len: u32,
    /// 1-based first base of the run.
    pub start1: u32,
}

/// `(base, ref_len, obs_len, strand, n)`, one row of `hp_confusion.tsv`.
pub type ConfusionRow = (char, u32, u32, char, u64);

/// `find_homopolymer_runs(seq, min_len)` of `hp_tracts.py`, minus runs whose 1-based `(start, end)` is in `skip`
/// and runs with no base before them.
pub fn homopolymer_runs(seq: &[u8], min_len: u32, skip: &[(u32, u32)]) -> Vec<HpRun> {
    let mut runs = Vec::new();
    let n = seq.len();
    let mut i = 0;
    while i < n {
        let mut j = i + 1;
        while j < n && seq[j] == seq[i] {
            j += 1;
        }
        let len = (j - i) as u32;
        let (start1, end1) = (i as u32 + 1, j as u32);
        if len >= min_len && !skip.contains(&(start1, end1)) && start1 >= 2 {
            runs.push(HpRun { anchor: start1 - 2, end: end1, base: seq[i], len, start1 });
        }
        i = j;
    }
    runs
}

/// Catalogue ids from a comma list, `hp_tracts` (every id), or `none`/empty.
fn parse_skip(spec: &str) -> Result<Vec<&'static str>> {
    let spec = spec.trim();
    if spec.is_empty() || spec.eq_ignore_ascii_case("none") {
        return Ok(Vec::new());
    }
    if spec == "hp_tracts" {
        return Ok(CATALOGUE.iter().map(|&(id, _, _)| id).collect());
    }
    spec.split(',')
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(|id| {
            CATALOGUE
                .iter()
                .find(|&&(c, _, _)| c == id)
                .map(|&(c, _, _)| c)
                .with_context(|| format!("unknown tract id {id:?} (known: {})", known_ids()))
        })
        .collect()
}

fn known_ids() -> String {
    CATALOGUE.iter().map(|&(id, _, _)| id).collect::<Vec<_>>().join(", ")
}

/// `{name: uppercased sequence}` in file order; the name is the first word of the header.
fn load_fasta(path: &Path) -> Result<Vec<(String, Vec<u8>)>> {
    let text = std::fs::read_to_string(path).with_context(|| format!("Cannot read FASTA {}", path.display()))?;
    let mut seqs: Vec<(String, Vec<u8>)> = Vec::new();
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        if let Some(header) = line.strip_prefix('>') {
            seqs.push((header.split_whitespace().next().unwrap_or("").to_string(), Vec::new()));
        } else if let Some((_, seq)) = seqs.last_mut() {
            seq.extend(line.bytes().map(|b| b.to_ascii_uppercase()));
        }
    }
    if seqs.is_empty() {
        bail!("no sequences in {}", path.display());
    }
    Ok(seqs)
}

/// Confusion rows for the reference homopolymer runs of `fasta`, from the primary alignments in `bam_path`.
pub fn measure(bam_path: &Path, fasta: &Path, min_len: u32, skip_ids: &[&str]) -> Result<Vec<ConfusionRow>> {
    let seqs = load_fasta(fasta)?;
    let skip: Vec<(u32, u32)> =
        CATALOGUE.iter().filter(|(id, _, _)| skip_ids.contains(id)).map(|&(_, s, e)| (s, e)).collect();
    let runs: Vec<Vec<HpRun>> = seqs.iter().map(|(_, s)| homopolymer_runs(s, min_len, &skip)).collect();
    let anchors: Vec<Vec<u32>> = runs.iter().map(|r| r.iter().map(|x| x.anchor).collect()).collect();

    let mut reader =
        bam::Reader::from_path(bam_path).with_context(|| format!("Cannot open BAM {}", bam_path.display()))?;
    let header = reader.header().clone();
    let by_name: HashMap<&str, usize> = seqs.iter().enumerate().map(|(i, (n, _))| (n.as_str(), i)).collect();
    let mut tid_to_seq: HashMap<i32, usize> = HashMap::new();
    for tid in 0..header.target_count() {
        let name = String::from_utf8_lossy(header.tid2name(tid)).to_string();
        if let Some(&i) = by_name.get(name.as_str()) {
            tid_to_seq.insert(tid as i32, i);
        }
    }
    // One reference on each side: map them whatever their names, as the Python script does.
    if tid_to_seq.is_empty() && seqs.len() == 1 && header.target_count() == 1 {
        tid_to_seq.insert(0, 0);
    }
    // The Python script silently writes an empty table here; an empty table only hides the mismatch.
    if tid_to_seq.is_empty() {
        bail!(
            "no BAM reference of {} matches a FASTA record of {} ({})",
            bam_path.display(),
            fasta.display(),
            seqs.iter().map(|(n, _)| n.as_str()).collect::<Vec<_>>().join(", ")
        );
    }

    // (sequence, run) -> (obs_len, strand) -> n
    let mut observed: HashMap<(usize, usize), HashMap<(u32, char), u64>> = HashMap::new();
    for rec in reader.records() {
        let rec = rec.context("failed to read BAM record")?;
        if rec.is_unmapped() || rec.is_secondary() || rec.is_supplementary() || rec.pos() < 0 {
            continue;
        }
        let Some(&si) = tid_to_seq.get(&rec.tid()) else { continue };
        let cigar = rec.cigar().take();
        if rec.seq_len() == 0 || cigar.is_empty() {
            continue;
        }
        let seq = rec.seq().as_bytes().to_ascii_uppercase();
        let (start, end) = (rec.pos() as u32, rec.reference_end().max(0) as u32);
        let strand = if rec.is_reverse() { '-' } else { '+' };
        let refseq = &seqs[si].1;
        let first = anchors[si].partition_point(|&a| a < start);
        for (ri, run) in runs[si].iter().enumerate().skip(first) {
            if run.anchor >= end {
                break;
            }
            if let Some(k) = crate::denoise_indel::hp_allele(rec.pos(), &cigar, &seq, refseq, run.anchor, run.end) {
                let obs = (run.len as i64 + k as i64).max(0) as u32;
                *observed.entry((si, ri)).or_default().entry((obs, strand)).or_insert(0) += 1;
            }
        }
    }

    let mut pooled: BTreeMap<(char, u32, u32, char), u64> = BTreeMap::new();
    for ((si, ri), items) in &observed {
        let run = &runs[*si][*ri];
        let mut by_len: HashMap<u32, u64> = HashMap::new();
        for (&(obs, _), &n) in items {
            *by_len.entry(obs).or_insert(0) += n;
        }
        let best = by_len.values().copied().max().unwrap_or(0);
        let modes: Vec<u32> = by_len.iter().filter(|&(_, &n)| n == best).map(|(&l, _)| l).collect();
        if modes != [run.len] {
            continue;
        }
        for (&(obs, strand), &n) in items {
            *pooled.entry((run.base as char, run.len, obs, strand)).or_insert(0) += n;
        }
    }
    Ok(pooled.into_iter().map(|((b, rl, ol, st), n)| (b, rl, ol, st, n)).collect())
}

const HEADER: &str = "base\tref_len\tobs_len\tstrand\tn";

/// Measure this sample and write `hp_confusion.tsv` at `path`. Runs shorter than [`DEFAULT_MIN_LEN`] and the
/// tracts in [`DEFAULT_SKIP`] are left out.
pub fn write_table(bam: &Path, fasta: &Path, path: &Path) -> Result<()> {
    let rows = measure(bam, fasta, DEFAULT_MIN_LEN, DEFAULT_SKIP)?;
    write_confusion(&rows, path)
}

pub fn write_confusion(rows: &[ConfusionRow], path: &Path) -> Result<()> {
    let mut w = BufWriter::new(File::create(path).with_context(|| format!("Cannot create {}", path.display()))?);
    writeln!(w, "{HEADER}")?;
    for &(b, rl, ol, st, n) in rows {
        writeln!(w, "{b}\t{rl}\t{ol}\t{st}\t{n}")?;
    }
    Ok(())
}

/// Parse `hp_confusion.tsv` text (taking text lets embedded defaults use `include_str!`).
pub fn read_confusion(text: &str) -> Result<Vec<ConfusionRow>> {
    let mut lines = text.lines();
    match lines.next() {
        Some(h) if h.trim_end() == HEADER => {}
        other => bail!("not an hp_confusion.tsv: header {other:?}, expected {HEADER:?}"),
    }
    let mut rows = Vec::new();
    for (i, line) in lines.enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        let bad = || format!("hp_confusion.tsv line {}: {line:?}", i + 2);
        if f.len() != 5 {
            bail!(bad());
        }
        let base = f[0].chars().next().with_context(bad)?;
        let strand = f[3].chars().next().with_context(bad)?;
        rows.push((
            base,
            f[1].parse().with_context(bad)?,
            f[2].parse().with_context(bad)?,
            strand,
            f[4].parse().with_context(bad)?,
        ));
    }
    Ok(rows)
}

#[cfg(test)]
mod tests {
    use super::*;

    //                          0         1         2         3
    //                          0123456789012345678901234567890 1
    const REF: &str = "GATCGATCGTAAAAAACTGACTGACTGACTGA";
    // AAAAAA is 0-based [10, 16) (1-based 11..16); anchor T at 9; hp_allele needs the read to reach 16 + 5 = 21.

    #[test]
    fn homopolymer_runs_matches_find_homopolymer_runs() {
        let seq = b"ACGTTTTACCCCCGAAAAT";
        let runs = homopolymer_runs(seq, 4, &[]);
        let got: Vec<(u32, u32, u8, u32)> = runs.iter().map(|r| (r.start1, r.end, r.base, r.len)).collect();
        assert_eq!(got, vec![(4, 7, b'T', 4), (9, 13, b'C', 5), (15, 18, b'A', 4)]);
        assert_eq!(runs[1].anchor, 7, "the base before the run, 0-based");
        assert_eq!(homopolymer_runs(seq, 5, &[]).len(), 1);
        assert_eq!(homopolymer_runs(seq, 4, &[(9, 13)]).len(), 2);
        // A run at the very start has no anchor base.
        assert!(homopolymer_runs(b"AAAACG", 4, &[]).is_empty());
    }

    #[test]
    fn parse_skip_accepts_ids_hp_tracts_and_none() {
        assert_eq!(parse_skip("D310,HVS1,CA514").unwrap(), DEFAULT_SKIP);
        assert_eq!(parse_skip("hp_tracts").unwrap().len(), CATALOGUE.len());
        assert!(parse_skip("none").unwrap().is_empty());
        assert!(parse_skip("D311").is_err());
    }

    /// One SAM record on chrM at 1-based POS 1.
    fn sam_line(name: &str, flag: u16, cigar: &str, seq: &str) -> String {
        format!("{name}\t{flag}\tchrM\t1\t60\t{cigar}\t*\t0\t0\t{seq}\t*\n")
    }

    fn write_inputs(tag: &str, reads: &[String]) -> (std::path::PathBuf, std::path::PathBuf) {
        let dir = std::env::temp_dir().join(format!("himito_hperror_{tag}"));
        std::fs::create_dir_all(&dir).unwrap();
        let fasta = dir.join("ref.fa");
        std::fs::write(&fasta, format!(">chrM test\n{REF}\n")).unwrap();
        let sam = dir.join("reads.sam");
        let mut text = format!("@HD\tVN:1.6\n@SQ\tSN:chrM\tLN:{}\n", REF.len());
        for r in reads {
            text.push_str(r);
        }
        std::fs::write(&sam, text).unwrap();
        (sam, fasta)
    }

    fn with_base(pos: usize, ins: &str) -> String {
        format!("{}{ins}{}", &REF[..pos], &REF[pos..])
    }

    #[test]
    fn measure_counts_hp_allele_per_primary_alignment() {
        let del = format!("{}{}", &REF[..12], &REF[13..]);
        let reads = vec![
            sam_line("exact", 0, "32M", REF),
            sam_line("exact_rev", 16, "32M", REF),
            sam_line("del", 16, "12M1D19M", &del),
            sam_line("ins", 0, "16M1I16M", &with_base(16, "A")),
            sam_line("bad_ins", 0, "16M1I16M", &with_base(16, "G")), // not the run base: uninformative
            sam_line("short", 0, "19M", &REF[..19]),                  // stops inside the flank
            sam_line("supp", 2048, "32M", &with_base(16, "A")[..32]), // supplementary: skipped
            sam_line("sec", 256, "32M", REF),                          // secondary: skipped
        ];
        let (sam, fasta) = write_inputs("counts", &reads);
        let rows = measure(&sam, &fasta, 4, &[]).unwrap();
        assert_eq!(
            rows,
            vec![('A', 6, 5, '-', 1), ('A', 6, 6, '+', 1), ('A', 6, 6, '-', 1), ('A', 6, 7, '+', 1)]
        );
    }

    #[test]
    fn modal_rule_drops_non_reference_runs() {
        // Every informative read is one A short: the sample is not reference here.
        let del = format!("{}{}", &REF[..12], &REF[13..]);
        let reads = vec![sam_line("d1", 0, "12M1D19M", &del), sam_line("d2", 16, "12M1D19M", &del)];
        let (sam, fasta) = write_inputs("modal", &reads);
        assert!(measure(&sam, &fasta, 4, &[]).unwrap().is_empty());
        // A tie between two lengths is dropped too.
        let reads = vec![sam_line("d1", 0, "12M1D19M", &del), sam_line("e1", 0, "32M", REF)];
        let (sam, fasta) = write_inputs("tie", &reads);
        assert!(measure(&sam, &fasta, 4, &[]).unwrap().is_empty());
    }

    #[test]
    fn write_then_read_confusion_round_trips() {
        let rows = vec![('A', 6, 5, '-', 1), ('C', 4, 4, '+', 12)];
        let path = std::env::temp_dir().join("himito_hperror_roundtrip.tsv");
        write_confusion(&rows, &path).unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        assert!(text.starts_with("base\tref_len\tobs_len\tstrand\tn\n"));
        assert_eq!(read_confusion(&text).unwrap(), rows);
        assert!(read_confusion("x\ty\n").is_err());
    }
}
