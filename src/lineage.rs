use std::collections::HashMap;
use anyhow::{Context, Result};
use rust_htslib::bcf::{Read, Reader}; // Read trait must be in scope for .records()
use std::path::PathBuf;

/// One unique binary profile and the set of reads that share it.
pub struct Haplotype {
    pub id: String,
    /// Binary state over filtered variants (`profile[v]` = 0 or 1)
    pub profile: Vec<u8>,
    pub reads: Vec<String>,
    pub count: usize,
}

/// Haplotypes together with the (shared) variant list.
pub struct HaplotypeMatrix {
    pub variants: Vec<String>,
    pub haplotypes: Vec<Haplotype>,
}

/// HF (heteroplasmic fraction) keyed by Himito variant ID (e.g. "m.13376T>C").
pub type HfMap = HashMap<String, f64>;

/// Binary matrix after filtering.
/// `data[variant_idx][read_idx]` is 0 or 1.
pub struct BinaryMatrix {
    pub variants: Vec<String>,
    pub reads: Vec<String>,
    pub data: Vec<Vec<u8>>,
}

/// Parse a Himito VCF file and return a map from variant ID to HF value.
///
/// Variant IDs are reconstructed as `m.<POS><REF>><ALT>` (the format Himito
/// uses in the matrix CSV), e.g. `m.13376T>C`.
pub fn parse_vcf(vcf_path: &str, min_hf: f64, max_hf: f64) -> Result<HfMap> {
    let mut reader = Reader::from_path(vcf_path)
        .with_context(|| format!("Cannot read VCF: {vcf_path}"))?;
    let mut map = HfMap::new();

    for result in reader.records() {
        #[allow(unused_mut)] // rec.format(b"HF") requires &mut self in rust_htslib
        let mut rec = result?;

        // Extract position (0-based in htslib → 1-based in VCF/Himito IDs)
        let pos = rec.pos() + 1;

        // rec.alleles() borrows rec immutably; collect into owned Strings so
        // we can release the borrow before calling rec.format() (mut borrow).
        let (ref_allele, alt_allele) = {
            let alleles = rec.alleles();
            if alleles.len() < 2 {
                continue;
            }
            (
                std::str::from_utf8(alleles[0]).unwrap_or("").to_owned(),
                std::str::from_utf8(alleles[1]).unwrap_or("").to_owned(),
            )
        };

        let vid = format!("m.{pos}{ref_allele}>{alt_allele}");

        // Read the HF FORMAT field (float, first sample, first value)
        if let Ok(hf_data) = rec.format(b"HF").float() {
            if let Some(&hf_val) = hf_data.get(0).and_then(|s| s.first()) {
                let hf = hf_val as f64;
                if hf >= min_hf && hf < max_hf {
                    map.insert(vid, hf);
                }
            }
        }
    }
    Ok(map)
}

/// Read `matrix_path` (Himito `.matrix.csv`) and apply filters:
///
/// * **HF bounds** — keep rows where `min_hf ≤ HF < max_hf`.
/// * **Prevalence** — keep rows where the variant is present in ≥ `min_presence`
///   reads AND absent from ≥ `min_absence` reads (guarantees a bifurcation).
///
/// Counts are binarised: any count ≥ 1 becomes 1.
pub fn load_and_filter_matrix(
    matrix_path: &str,
    hf_map: &HfMap,   // contains only variants that already passed HF filtering
    min_presence: usize,
    min_absence: usize,
) -> Result<BinaryMatrix> {
    let mut rdr = csv::Reader::from_path(matrix_path)
        .with_context(|| format!("Cannot read matrix CSV: {matrix_path}"))?;

    // Column 0 is "variant", then one column per read
    let headers = rdr.headers()?.clone();
    let reads: Vec<String> = headers.iter().skip(1).map(String::from).collect();

    let mut variants: Vec<String> = Vec::new();
    let mut data: Vec<Vec<u8>>   = Vec::new();

    for result in rdr.records() {
        let rec = result?;
        let vid = rec[0].to_string();

        // Keep only variants that passed HF filtering in parse_vcf.
        // If hf_map is empty (no VCF provided), terminate the program.
        if !hf_map.is_empty() && !hf_map.contains_key(&vid) {
            continue;
        }

        // Binarise counts
        let row: Vec<u8> = rec
            .iter()
            .skip(1)
            .map(|v| u8::from(v.parse::<u32>().unwrap_or(0) >= 1))
            .collect();

        // Prevalence filter
        let present = row.iter().filter(|&&b| b == 1).count();
        let absent  = row.iter().filter(|&&b| b == 0).count();
        if present < min_presence || absent < min_absence {
            continue;
        }

        variants.push(vid);
        data.push(row);
    }

    Ok(BinaryMatrix { variants, reads, data })
}

/// Return a matrix containing only the variant rows at `indices` (in order).
pub fn subset_variants(matrix: &BinaryMatrix, indices: &[usize]) -> BinaryMatrix {
    BinaryMatrix {
        variants: indices
            .iter()
            .map(|&i| matrix.variants[i].clone())
            .collect(),
        reads: matrix.reads.clone(),
        data: indices.iter().map(|&i| matrix.data[i].clone()).collect(),
    }
}

/// Collapse reads with identical binary profiles into unique haplotypes.
///
/// Haplotypes are sorted by ascending mutation count (number of 1-bits),
/// then lexicographically by profile, so that the all-zero haplotype (if
/// present) is always `H0000`.
pub fn deduplicate(matrix: &BinaryMatrix, min_reads: usize) -> HaplotypeMatrix {
    let n_reads    = matrix.reads.len();
    let n_variants = matrix.variants.len();

    // profile (per-read, across variants) → list of read indices
    let mut profile_map: HashMap<Vec<u8>, Vec<usize>> = HashMap::new();

    for read_idx in 0..n_reads {
        let profile: Vec<u8> = (0..n_variants)
            .map(|v| matrix.data[v][read_idx])
            .collect();
        profile_map.entry(profile).or_default().push(read_idx);
    }

    let mut haplotypes: Vec<Haplotype> = profile_map
        .into_iter()
        .map(|(profile, indices)| {
            let reads = indices.iter().map(|&r| matrix.reads[r].clone()).collect();
            let count = indices.len();
            Haplotype {
                id: String::new(), // assigned after sort
                profile,
                reads,
                count,
            }
        })
        .collect();
    // filter out haplotypes with less than min_reads
    haplotypes = haplotypes.into_iter().filter(|h| h.count >= min_reads).collect();
    // Deterministic ordering: fewest mutations first, then lexicographic profile
    haplotypes.sort_by(|a, b| {
        let ma: usize = a.profile.iter().map(|&x| x as usize).sum();
        let mb: usize = b.profile.iter().map(|&x| x as usize).sum();
        ma.cmp(&mb).then_with(|| a.profile.cmp(&b.profile))
    });

    for (i, h) in haplotypes.iter_mut().enumerate() {
        h.id = format!("H{i:04}");
    }

    HaplotypeMatrix {
        variants: matrix.variants.clone(),
        haplotypes,
    }
}


/// Run the four-gamete test on every pair of variants.
///
/// Returns a list of `(variant_a, variant_b)` pairs where all four
/// combinations {00, 01, 10, 11} are observed across haplotypes — a
/// violation of the Infinite Sites Assumption (ISA).
pub fn four_gamete_test(matrix: &HaplotypeMatrix) -> Vec<(String, String)> {
    let n_var = matrix.variants.len();
    let mut violations = Vec::new();

    for i in 0..n_var {
        for j in (i + 1)..n_var {
            let mut seen = [false; 4]; // index = a*2 + b for states a,b ∈ {0,1}
            for h in &matrix.haplotypes {
                let idx = (h.profile[i] * 2 + h.profile[j]) as usize;
                seen[idx] = true;
                if seen.iter().all(|&s| s) {
                    break;
                }
            }
            if seen.iter().all(|&s| s) {
                violations.push((
                    matrix.variants[i].clone(),
                    matrix.variants[j].clone(),
                ));
            }
        }
    }

    violations
}


pub fn start (matrix_file: &str, vcf_file: Option<&str>, min_hf: f64, max_hf: f64, min_presence: usize, min_absence: usize, min_reads: usize) -> Result<()> {
    println!("Starting lineage analysis...");

    // ── Step 1: load VCF HF values, then filter the binary matrix ─────────────
    println!("[1/6] Parsing VCF: {}", vcf_file.as_ref().unwrap_or(&""));
    let hf_map = parse_vcf(vcf_file.as_ref().unwrap_or(&""), min_hf, max_hf)?;

    eprintln!("[1/6] Loading and filtering matrix: {}", matrix_file);
    let binary = load_and_filter_matrix(
        &matrix_file,
        &hf_map,
        min_presence,
        min_absence,
    )?;

    // ── Step 2: deduplicate reads into haplotypes ───────────────────────────────
    eprintln!("[2/6] Deduplicating reads into haplotypes...");
    let hap_matrix = deduplicate(&binary, min_reads);
    // write how many haplotypes and how many heteroplasmic variants on each haplotype
    eprintln!("[2/6] Found {} haplotypes across {} variants.", hap_matrix.haplotypes.len(), hap_matrix.variants.len());


    Ok(())

}