//! Describe matrix changes without treating inferred genotypes as ground truth.

use std::collections::{BTreeSet, HashMap};
use std::fs::File;
use std::io::{BufWriter, Write};

use anyhow::{Context, Result};
use serde_json::{json, Value};

use crate::lineage::{BinaryMatrix, RawMatrix};
use crate::scite::{CleanedMatrix, ErrorRates};

#[derive(Default)]
struct Changes {
    rows: usize,
    // 0→0, 0→1, 1→0, 1→1, missing→0, missing→1.
    transitions: [usize; 6],
    added_cells_without_input_row: usize,
}

impl Changes {
    fn compare(input: Option<&[Option<u8>]>, output: &[u8]) -> Self {
        let mut changes = Self { rows: 1, ..Self::default() };
        if let Some(input) = input {
            for (&before, &after) in input.iter().zip(output) {
                let index = match before {
                    Some(v) => usize::from(v) * 2 + usize::from(after),
                    None => 4 + usize::from(after),
                };
                changes.transitions[index] += 1;
            }
        } else {
            changes.added_cells_without_input_row = output.len();
        }
        changes
    }

    fn add(&mut self, other: &Self) {
        self.rows += other.rows;
        for (a, b) in self.transitions.iter_mut().zip(other.transitions) {
            *a += b;
        }
        self.added_cells_without_input_row += other.added_cells_without_input_row;
    }

    fn observed(&self) -> usize {
        self.transitions[..4].iter().sum()
    }

    fn json(&self) -> Value {
        let [ref_ref, ref_alt, alt_ref, alt_alt, missing_ref, missing_alt] = self.transitions;
        let observed = self.observed();
        let changed = ref_alt + alt_ref;
        json!({
            "variant_rows": self.rows,
            "observed_calls": observed,
            "observed_ref_to_ref": ref_ref,
            "observed_ref_to_alt": ref_alt,
            "observed_alt_to_ref": alt_ref,
            "observed_alt_to_alt": alt_alt,
            "observed_calls_changed": changed,
            "observed_call_change_fraction": (observed > 0).then(|| changed as f64 / observed as f64),
            "missing_to_ref": missing_ref,
            "missing_to_alt": missing_alt,
            "missing_calls_filled": missing_ref + missing_alt,
            "added_cells_without_input_row": self.added_cells_without_input_row,
        })
    }
}

/// Compare matched read/variant cells, keeping inferred calls and annotations separate.
pub fn write_imputation_summary(
    raw: &RawMatrix,
    search: &BinaryMatrix,
    cleaned: &CleanedMatrix,
    homoplasmic: &BTreeSet<String>,
    rates: &ErrorRates,
    log_likelihood: f64,
    path: &str,
) -> Result<()> {
    let raw_rows: HashMap<_, _> = raw.variants.iter().zip(&raw.data)
        .map(|(name, row)| (name.as_str(), row.as_slice())).collect();
    let search_rows: HashMap<_, _> = search.variants.iter().enumerate()
        .map(|(i, name)| (name.as_str(), i)).collect();
    let exported: BTreeSet<_> = cleaned.variants.iter().map(String::as_str).collect();
    let mut total = Changes::default();
    let mut modeled = Changes::default();
    let mut fixed = Changes::default();
    let mut fallback = Changes::default();
    let mut variants = Vec::with_capacity(cleaned.variants.len());

    for (name, row) in cleaned.variants.iter().zip(&cleaned.data) {
        let search_idx = search_rows.get(name.as_str()).copied();
        let before = search_idx.map(|i| search.data[i].as_slice())
            .or_else(|| raw_rows.get(name.as_str()).copied());
        let changes = Changes::compare(before, row);
        let (source, category) = if search_idx.is_some() {
            ("scite", &mut modeled)
        } else if homoplasmic.contains(name) {
            ("homoplasmic", &mut fixed)
        } else {
            ("zero_fallback", &mut fallback)
        };
        category.add(&changes);
        total.add(&changes);
        let model_rates = search_idx.map(|i| {
            let (fp, fnr) = rates.clamped_for(i);
            json!({ "fp": fp, "fn": fnr })
        });
        variants.push(json!({
            "variant": name,
            "source": source,
            "model_error_rates": model_rates,
            "changes": changes.json(),
        }));
    }

    let report = json!({
        "n_reads": cleaned.reads.len(),
        "n_input_variant_rows": raw.variants.len(),
        "n_output_variant_rows": cleaned.variants.len(),
        "input_variant_rows_absent_from_export": raw.variants.iter()
            .filter(|name| !exported.contains(name.as_str())).count(),
        "model": {
            "status": if search.variants.is_empty() { "no_informative_variants" } else { "fitted" },
            "global_fp": rates.fp_rate,
            "global_fn": rates.fn_rate,
            "final_log_likelihood": log_likelihood,
            "log_likelihood_per_observed_call": (modeled.observed() > 0)
                .then(|| log_likelihood / modeled.observed() as f64),
        },
        "total_changes": total.json(),
        "scite_changes": modeled.json(),
        "homoplasmic_changes": fixed.json(),
        "zero_fallback_changes": fallback.json(),
        "variants": variants,
    });
    let mut writer = BufWriter::new(File::create(path).with_context(|| format!("Cannot create {path}"))?);
    serde_json::to_writer_pretty(&mut writer, &report)?;
    writeln!(writer)?;
    writer.flush()?;
    Ok(())
}
