// Sequence-resolved NUMT caller.
//
// Parses read CIGARs (primary + SA) to locate nuclear<->mitochondrial junctions,
// extracts the mt-derived sequence carried by each read, clusters junctions by
// signature, applies false-positive filters, and writes a sequence-resolved VCF
// (explicit-sequence INS when a single read spans both boundaries, BND fallback
// otherwise). Everything lives in this one file, organized into sections:
//   - align:    AlignSeg model + CIGAR arithmetic (pure)
//   - junction: per-read junction extraction + sequence resolution (pure)
//   - cluster:  junction-signature single-linkage clustering + AC (pure)
//   - filters:  FP metrics -> FILTER tags (pure)
//   - vcf:      VCF record construction + writer
//   - config + start(): orchestrator + BAM reading
//
// Grounding (2026-07-16 field survey): thresholds follow the split-read SV / NUMT
// caller consensus (Sniffles2/cuteSV MAPQ 20, ~100bp clustering; dinumt ref-NUMT
// masking; Sniffles2 3+0.1/kb max-splits chimera rule). PALMER's TSD+poly(A)
// hallmarks are TPRT-specific and deliberately NOT used (NUMTs do not integrate
// by target-primed reverse transcription).

use std::{collections::{HashMap, HashSet}, path::PathBuf, io::Write, fs::File};
use rust_htslib::bam::{self, Read, record::{Aux, Cigar}};
use rust_htslib::faidx;

// ---------------------------------------------------------------------------
// config
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
pub struct CallNumtsConfig {
    pub tolerance: i64,        // clustering epsilon (bp)
    pub min_mapq: u8,          // segment MAPQ floor
    pub min_seg_len: i64,      // min aligned segment length (bp)
    pub ac_threshold: usize,   // min distinct supporting reads
    pub max_splits: i64,       // base allowance for supplementary alignments per read
    pub max_splits_per_kb: f64,// extra split allowance per read-kb
    pub emit_bnd: bool,        // force breakend output for all calls
    pub ref_numt_bed: Option<PathBuf>,
}

impl Default for CallNumtsConfig {
    fn default() -> Self {
        CallNumtsConfig {
            tolerance: 100,
            min_mapq: 20,
            min_seg_len: 100,
            ac_threshold: 2,
            max_splits: 3,
            max_splits_per_kb: 0.1,
            emit_bnd: false,
            ref_numt_bed: None,
        }
    }
}

// ---------------------------------------------------------------------------
// align section
// ---------------------------------------------------------------------------

/// One alignment of a read (primary record OR one SA entry), with query spans
/// normalized to the PRIMARY record's stored-SEQ frame ("record frame"), so that
/// spans from the primary and from SA entries are directly comparable and can
/// index the primary's `seq()`.
#[derive(Debug, Clone)]
pub struct AlignSeg {
    pub contig: String,
    pub is_mt: bool,
    pub ref_start: i64,       // 0-based, ref-forward
    pub ref_end: i64,         // exclusive
    pub q_start: i64,         // 0-based in record frame
    pub q_end: i64,           // exclusive
    pub strand: char,         // '+' | '-' (absolute, vs reference forward)
    pub mapq: u8,
    pub seg_len: i64,         // aligned reference length
    pub same_as_primary: bool,// segment orientation matches the primary's
}

/// Parse a CIGAR string ("10S90M") into (op_byte, length) pairs.
fn parse_cigar_str(s: &str) -> Vec<(u8, i64)> {
    let mut out = Vec::new();
    let mut num = String::new();
    for ch in s.chars() {
        if ch.is_ascii_digit() {
            num.push(ch);
        } else {
            if let Ok(n) = num.parse::<i64>() {
                out.push((ch as u8, n));
            }
            num.clear();
        }
    }
    out
}

/// Convert an htslib Cigar op into a (op_byte, length) pair.
fn cigar_to_pair(c: &Cigar) -> (u8, i64) {
    match c {
        Cigar::Match(l) => (b'M', *l as i64),
        Cigar::Ins(l) => (b'I', *l as i64),
        Cigar::Del(l) => (b'D', *l as i64),
        Cigar::RefSkip(l) => (b'N', *l as i64),
        Cigar::SoftClip(l) => (b'S', *l as i64),
        Cigar::HardClip(l) => (b'H', *l as i64),
        Cigar::Pad(l) => (b'P', *l as i64),
        Cigar::Equal(l) => (b'=', *l as i64),
        Cigar::Diff(l) => (b'X', *l as i64),
    }
}

/// From CIGAR pairs compute (ref_len, leading clip, trailing clip, aligned query len).
/// Leading/trailing clip counts both soft (S) and hard (H) clips so that spans are
/// on a consistent full-read coordinate axis across primary (soft-clipped) and
/// supplementary (often hard-clipped) alignments.
fn cigar_metrics(ops: &[(u8, i64)]) -> (i64, i64, i64, i64) {
    let mut ref_len = 0i64;
    let mut qlen = 0i64;
    for (op, len) in ops {
        match op {
            b'M' | b'D' | b'N' | b'=' | b'X' => ref_len += len,
            _ => {}
        }
        match op {
            b'M' | b'I' | b'=' | b'X' => qlen += len,
            _ => {}
        }
    }
    // Sum consecutive clip ops at each end (a valid SAM clip run may be H then S,
    // e.g. "10H20S..."), so spans land on a consistent full-read coordinate axis
    // across soft-clipped primaries and hard-clipped supplementaries. A valid
    // alignment always has >=1 consuming op between the runs, so they never overlap.
    let is_clip = |op: u8| op == b'S' || op == b'H';
    let mut lead = 0i64;
    for (op, len) in ops {
        if is_clip(*op) {
            lead += len;
        } else {
            break;
        }
    }
    let mut trail = 0i64;
    for (op, len) in ops.iter().rev() {
        if is_clip(*op) {
            trail += len;
        } else {
            break;
        }
    }
    (ref_len, lead, trail, qlen)
}

/// Build an AlignSeg from CIGAR pairs and orientation info.
/// `ref_start` is 0-based. `seg_is_reverse` / `primary_is_reverse` are absolute
/// (vs reference forward). Query spans are placed in the primary's record frame.
fn build_seg(
    contig: String,
    is_mt: bool,
    ref_start: i64,
    ops: &[(u8, i64)],
    seg_is_reverse: bool,
    primary_is_reverse: bool,
    mapq: u8,
) -> AlignSeg {
    let (ref_len, lead, trail, qlen) = cigar_metrics(ops);
    let same = seg_is_reverse == primary_is_reverse;
    let (q_start, q_end) = if same {
        (lead, lead + qlen)
    } else {
        (trail, trail + qlen)
    };
    AlignSeg {
        contig,
        is_mt,
        ref_start,
        ref_end: ref_start + ref_len,
        q_start,
        q_end,
        strand: if seg_is_reverse { '-' } else { '+' },
        mapq,
        seg_len: ref_len,
        same_as_primary: same,
    }
}

/// Reference base index (0-based) of the junction-facing anchor of a segment.
/// `is_left` = this segment is the upstream (smaller record-frame q) side of the
/// junction, so its junction-facing end is the q_end side.
fn anchor_ref(seg: &AlignSeg, is_left: bool) -> i64 {
    match (is_left, seg.same_as_primary) {
        (true, true) => seg.ref_end - 1,   // last aligned base, upstream in ref
        (true, false) => seg.ref_start,    // reversed: read-3' maps to ref_start
        (false, true) => seg.ref_start,    // downstream segment's first aligned base
        (false, false) => seg.ref_end - 1,
    }
}

fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            _ => b'N',
        })
        .collect()
}

// ---------------------------------------------------------------------------
// junction section
// ---------------------------------------------------------------------------

/// A single nuclear<->mt junction extracted from one read.
#[derive(Debug, Clone)]
pub struct Junction {
    pub qname: String,
    pub nuclear_chrom: String,
    pub nuclear_bp: i64,       // 0-based anchor reference base
    pub mt_chrom: String,
    pub mt_bp: i64,            // 0-based mt breakpoint base
    pub mt_ref_start: i64,     // mt segment span (for MTREGION)
    pub mt_ref_end: i64,
    pub orientation: char,     // '+' same strand, '-' opposite
    pub microhom: i64,         // >=0 read-coord overlap (NHEJ junction signal, not a TSD)
    pub gap: i64,              // >=0 unaligned read bases between segments
    pub ins_seq: Vec<u8>,      // mt-segment read bases, nuclear-forward oriented
    pub is_full_ins: bool,     // nuclear-mt-nuclear spanned by this read
    pub n_splits: usize,       // total alignment segments for this read
    pub read_len: i64,         // full read length (primary stored SEQ), for chimera rule
    pub min_seg_mapq: u8,
    pub min_seg_len: i64,
}

/// Extract junctions from one read's segment set (already in record frame) and the
/// primary's stored SEQ (record frame). Emits at most one junction per mt segment.
pub fn extract_junctions(qname: &str, segs: &[AlignSeg], seq: &[u8]) -> Vec<Junction> {
    let mut sorted: Vec<AlignSeg> = segs.to_vec();
    sorted.sort_by_key(|s| s.q_start);
    let n = sorted.len();
    let mut out = Vec::new();

    for i in 0..n {
        if !sorted[i].is_mt {
            continue;
        }
        let mt = &sorted[i];
        let left = if i > 0 { Some(&sorted[i - 1]) } else { None };
        let right = if i + 1 < n { Some(&sorted[i + 1]) } else { None };
        let left_nuc = left.map_or(false, |s| !s.is_mt);
        let right_nuc = right.map_or(false, |s| !s.is_mt);

        // Choose the nuclear anchor and whether this is a full insertion.
        let (nuclear, nuclear_is_left, is_full) = if left_nuc
            && right_nuc
            && left.unwrap().contig == right.unwrap().contig
        {
            (left.unwrap(), true, true) // nuclear-mt-nuclear: anchor on upstream nuclear
        } else if left_nuc {
            (left.unwrap(), true, false)
        } else if right_nuc {
            (right.unwrap(), false, false)
        } else {
            continue; // mt segment with no nuclear neighbor
        };

        let nuclear_bp = anchor_ref(nuclear, nuclear_is_left);
        let mt_bp = anchor_ref(mt, !nuclear_is_left);

        // microhomology / gap in record-frame read coordinates at the junction.
        let (a, b) = if nuclear_is_left { (nuclear, mt) } else { (mt, nuclear) };
        let microhom = (a.q_end - b.q_start).max(0);
        let gap = (b.q_start - a.q_end).max(0);

        // mt-derived sequence, sliced from the record frame, oriented nuclear-forward.
        let lo = mt.q_start.max(0) as usize;
        let hi = (mt.q_end.max(0) as usize).min(seq.len());
        let mut ins_seq = if lo < hi { seq[lo..hi].to_vec() } else { Vec::new() };
        // record frame == nuclear-forward iff the nuclear segment matches the
        // primary orientation; otherwise reverse-complement to nuclear-forward.
        if !nuclear.same_as_primary {
            ins_seq = revcomp(&ins_seq);
        }

        let orientation = if nuclear.strand == mt.strand { '+' } else { '-' };
        let min_seg_mapq = nuclear.mapq.min(mt.mapq);
        let min_seg_len = nuclear.seg_len.min(mt.seg_len);

        out.push(Junction {
            qname: qname.to_string(),
            nuclear_chrom: nuclear.contig.clone(),
            nuclear_bp,
            mt_chrom: mt.contig.clone(),
            mt_bp,
            mt_ref_start: mt.ref_start,
            mt_ref_end: mt.ref_end,
            orientation,
            microhom,
            gap,
            ins_seq,
            is_full_ins: is_full,
            n_splits: n,
            read_len: seq.len() as i64,
            min_seg_mapq,
            min_seg_len,
        });
    }
    out
}

// ---------------------------------------------------------------------------
// cluster section
// ---------------------------------------------------------------------------

/// A clustered NUMT call.
#[derive(Debug, Clone)]
pub struct NumtCall {
    pub nuclear_chrom: String,
    pub nuclear_bp: i64,        // representative (median) 0-based anchor
    pub mt_chrom: String,
    pub mt_bp: i64,
    pub mt_region_start: i64,
    pub mt_region_end: i64,
    pub orientation: char,
    pub ac: usize,             // distinct qnames
    pub ins_seq: Vec<u8>,      // longest/consensus resolved sequence
    pub is_full_ins: bool,
    pub microhom: i64,
    pub min_seg_mapq: u8,
    pub min_seg_len: i64,
    pub mt_bp_spread: i64,
    pub max_n_splits: usize,
    pub max_split_read_len: i64, // read length of the most-fragmented supporting read
    pub in_ref_numt: bool,
    pub filters: Vec<String>,  // FILTER tags, empty => PASS
}

/// Union-find helper.
fn uf_find(parent: &mut [usize], x: usize) -> usize {
    let mut r = x;
    while parent[r] != r {
        r = parent[r];
    }
    let mut cur = x;
    while parent[cur] != r {
        let next = parent[cur];
        parent[cur] = r;
        cur = next;
    }
    r
}

/// Single-linkage cluster junctions sharing (nuclear_chrom, nuclear_bp +/- eps,
/// mt_chrom, mt_bp +/- eps, orientation). Returns one NumtCall per cluster.
pub fn cluster_junctions(junctions: &[Junction], eps: i64) -> Vec<NumtCall> {
    let n = junctions.len();
    let mut parent: Vec<usize> = (0..n).collect();

    for i in 0..n {
        for j in (i + 1)..n {
            let a = &junctions[i];
            let b = &junctions[j];
            if a.nuclear_chrom == b.nuclear_chrom
                && a.mt_chrom == b.mt_chrom
                && a.orientation == b.orientation
                && (a.nuclear_bp - b.nuclear_bp).abs() <= eps
                && (a.mt_bp - b.mt_bp).abs() <= eps
            {
                let ra = uf_find(&mut parent, i);
                let rb = uf_find(&mut parent, j);
                if ra != rb {
                    parent[ra] = rb;
                }
            }
        }
    }

    let mut groups: HashMap<usize, Vec<usize>> = HashMap::new();
    for i in 0..n {
        let r = uf_find(&mut parent, i);
        groups.entry(r).or_default().push(i);
    }

    let mut calls = Vec::new();
    for (_root, idxs) in groups {
        let members: Vec<&Junction> = idxs.iter().map(|&i| &junctions[i]).collect();

        let mut nuc_bps: Vec<i64> = members.iter().map(|m| m.nuclear_bp).collect();
        let mut mt_bps: Vec<i64> = members.iter().map(|m| m.mt_bp).collect();
        nuc_bps.sort();
        mt_bps.sort();
        let median_nuclear_bp = nuc_bps[nuc_bps.len() / 2];
        let median_mt_bp = mt_bps[mt_bps.len() / 2];
        let mt_bp_spread = mt_bps[mt_bps.len() - 1] - mt_bps[0];

        // distinct supporting reads
        let mut qnames: Vec<&str> = members.iter().map(|m| m.qname.as_str()).collect();
        qnames.sort();
        qnames.dedup();
        let ac = qnames.len();

        // representative sequence: prefer full-ins, then longest.
        let rep = members
            .iter()
            .max_by_key(|m| (m.is_full_ins as usize, m.ins_seq.len()))
            .unwrap();

        let is_full_ins = members.iter().any(|m| m.is_full_ins);
        // For an explicit-sequence INS the emitted ALT comes from `rep`; anchor POS
        // to the SAME read so POS/REF/ins_seq/MTREGION are internally consistent.
        // Otherwise (BND fallback) use the robust cluster median.
        let (nuclear_bp, mt_bp) = if is_full_ins {
            (rep.nuclear_bp, rep.mt_bp)
        } else {
            (median_nuclear_bp, median_mt_bp)
        };
        let min_seg_mapq = members.iter().map(|m| m.min_seg_mapq).min().unwrap();
        let min_seg_len = members.iter().map(|m| m.min_seg_len).min().unwrap();
        // Chimera rule is per-read: evaluate the most-fragmented supporting read
        // against the allowance for ITS OWN length.
        let most_split = members.iter().max_by_key(|m| m.n_splits).unwrap();
        let max_n_splits = most_split.n_splits;
        let max_split_read_len = most_split.read_len;

        calls.push(NumtCall {
            nuclear_chrom: members[0].nuclear_chrom.clone(),
            nuclear_bp,
            mt_chrom: members[0].mt_chrom.clone(),
            mt_bp,
            mt_region_start: rep.mt_ref_start,
            mt_region_end: rep.mt_ref_end,
            orientation: members[0].orientation,
            ac,
            ins_seq: rep.ins_seq.clone(),
            is_full_ins,
            microhom: rep.microhom,
            min_seg_mapq,
            min_seg_len,
            mt_bp_spread,
            max_n_splits,
            max_split_read_len,
            in_ref_numt: false,
            filters: Vec::new(),
        });
    }

    calls.sort_by(|a, b| {
        a.nuclear_chrom
            .cmp(&b.nuclear_chrom)
            .then(a.nuclear_bp.cmp(&b.nuclear_bp))
    });
    calls
}

// ---------------------------------------------------------------------------
// filters section
// ---------------------------------------------------------------------------

/// Maximum allowed alignment segments for a read of `read_len_bp` bases
/// (Sniffles2 "3 + 0.1/kb" rule). Split count above this flags a chimera.
fn max_allowed_splits(cfg: &CallNumtsConfig, read_len_bp: i64) -> i64 {
    cfg.max_splits + (cfg.max_splits_per_kb * (read_len_bp as f64 / 1000.0)).round() as i64
}

/// Populate `filters` for a call. Empty => PASS.
pub fn apply_filters(call: &mut NumtCall, cfg: &CallNumtsConfig) {
    call.filters.clear();
    if call.min_seg_mapq < cfg.min_mapq {
        call.filters.push("LOWMAPQ".to_string());
    }
    if call.min_seg_len < cfg.min_seg_len {
        call.filters.push("SHORTSEG".to_string());
    }
    if call.ac < cfg.ac_threshold {
        call.filters.push("LOW_SUPPORT".to_string());
    }
    if call.mt_bp_spread > cfg.tolerance {
        call.filters.push("DISCORDANT".to_string());
    }
    // Chimera: compare the most-fragmented supporting read's segment count against
    // the "3 + 0.1/kb" allowance for THAT read's actual length.
    if (call.max_n_splits as i64) > max_allowed_splits(cfg, call.max_split_read_len) {
        call.filters.push("MULTISPLIT".to_string());
    }
    if call.in_ref_numt {
        call.filters.push("REFNUMT".to_string());
    }
}

/// Parse a BED file into per-chrom sorted (start,end) intervals (0-based, half-open).
fn load_bed(path: &PathBuf) -> Result<HashMap<String, Vec<(i64, i64)>>, Box<dyn std::error::Error>> {
    let content = std::fs::read_to_string(path)?;
    let mut map: HashMap<String, Vec<(i64, i64)>> = HashMap::new();
    for line in content.lines() {
        if line.is_empty() || line.starts_with('#') || line.starts_with("track") {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 3 {
            continue;
        }
        if let (Ok(s), Ok(e)) = (f[1].parse::<i64>(), f[2].parse::<i64>()) {
            map.entry(f[0].to_string()).or_default().push((s, e));
        }
    }
    for v in map.values_mut() {
        v.sort();
    }
    Ok(map)
}

fn bed_overlaps(bed: &HashMap<String, Vec<(i64, i64)>>, chrom: &str, pos: i64) -> bool {
    if let Some(intervals) = bed.get(chrom) {
        intervals.iter().any(|(s, e)| pos >= *s && pos < *e)
    } else {
        false
    }
}

// ---------------------------------------------------------------------------
// vcf section
// ---------------------------------------------------------------------------

/// BND ALT string for a nuclear-anchored breakend pointing at the mt mate.
fn bnd_alt(ref_base: &str, mt_chrom: &str, mt_pos_1based: i64, orientation: char) -> String {
    if orientation == '+' {
        format!("{}[{}:{}[", ref_base, mt_chrom, mt_pos_1based)
    } else {
        format!("{}]{}:{}]", ref_base, mt_chrom, mt_pos_1based)
    }
}

fn write_vcf(
    output_vcf: &PathBuf,
    calls: &[NumtCall],
    fasta: &faidx::Reader,
    sample_name: &str,
    emit_bnd: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut vcf = File::create(output_vcf)?;
    writeln!(vcf, "##fileformat=VCFv4.3")?;
    for i in 0..fasta.n_seqs() as i32 {
        let name = fasta.seq_name(i)?;
        let len = fasta.fetch_seq_len(&name);
        writeln!(vcf, "##contig=<ID={},length={}>", name, len)?;
    }
    writeln!(vcf, "##ALT=<ID=INS,Description=\"NUMT insertion (mt-derived sequence)\">")?;
    writeln!(vcf, "##ALT=<ID=BND,Description=\"NUMT junction breakend\">")?;
    writeln!(vcf, "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">")?;
    writeln!(vcf, "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of inserted mt sequence\">")?;
    writeln!(vcf, "##INFO=<ID=MTREGION,Number=1,Type=String,Description=\"mt origin chrom:start-end (1-based)\">")?;
    writeln!(vcf, "##INFO=<ID=MTSTRAND,Number=1,Type=String,Description=\"Relative orientation nuclear->mt\">")?;
    writeln!(vcf, "##INFO=<ID=MICROHOM,Number=1,Type=Integer,Description=\"Junction microhomology length (NHEJ signal, not a TSD)\">")?;
    writeln!(vcf, "##INFO=<ID=SVINSSEQ,Number=1,Type=String,Description=\"Inserted mt sequence (partial for BND)\">")?;
    writeln!(vcf, "##INFO=<ID=SEGMAPQ,Number=1,Type=Integer,Description=\"Min segment MAPQ\">")?;
    writeln!(vcf, "##INFO=<ID=SEGLEN,Number=1,Type=Integer,Description=\"Min aligned segment length\">")?;
    writeln!(vcf, "##INFO=<ID=MTBPSPREAD,Number=1,Type=Integer,Description=\"mt breakpoint spread across supporting reads\">")?;
    writeln!(vcf, "##INFO=<ID=NSPLITS,Number=1,Type=Integer,Description=\"Max alignment-segment count among supporting reads\">")?;
    writeln!(vcf, "##FILTER=<ID=LOWMAPQ,Description=\"segment MAPQ below threshold\">")?;
    writeln!(vcf, "##FILTER=<ID=SHORTSEG,Description=\"aligned segment shorter than threshold\">")?;
    writeln!(vcf, "##FILTER=<ID=LOW_SUPPORT,Description=\"AC below threshold\">")?;
    writeln!(vcf, "##FILTER=<ID=DISCORDANT,Description=\"inconsistent mt breakpoint across supporting reads\">")?;
    writeln!(vcf, "##FILTER=<ID=MULTISPLIT,Description=\"supporting read shredded into too many supplementary alignments (chimera)\">")?;
    writeln!(vcf, "##FILTER=<ID=REFNUMT,Description=\"overlaps a known reference NUMT (from --ref-numt-bed)\">")?;
    writeln!(vcf, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")?;
    writeln!(vcf, "##FORMAT=<ID=SR,Number=1,Type=Integer,Description=\"Distinct supporting read count (AC)\">")?;
    writeln!(vcf, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}", sample_name)?;

    for call in calls {
        let contig_len = fasta.fetch_seq_len(&call.nuclear_chrom);
        if contig_len == 0 {
            return Err(format!(
                "contig {} not found in reference FASTA",
                call.nuclear_chrom
            )
            .into());
        }
        if call.nuclear_bp < 0 || call.nuclear_bp as u64 >= contig_len {
            return Err(format!(
                "breakpoint {}:{} out of range for reference",
                call.nuclear_chrom, call.nuclear_bp
            )
            .into());
        }
        let pos = call.nuclear_bp + 1; // 1-based
        // faidx fetch_seq uses 0-based, INCLUSIVE end: [p, p] returns a single base.
        let bp = call.nuclear_bp as usize;
        let base_bytes = fasta.fetch_seq(&call.nuclear_chrom, bp, bp)?;
        let ref_base = base_bytes
            .first()
            .map(|b| (b.to_ascii_uppercase() as char).to_string())
            .ok_or_else(|| {
                format!(
                    "empty reference base at {}:{}",
                    call.nuclear_chrom, call.nuclear_bp
                )
            })?;
        let ins_seq_str = String::from_utf8_lossy(&call.ins_seq).to_string();
        let variant_id = format!(
            "NUMT_{}_{}_{}_{}_{}",
            call.nuclear_chrom, call.nuclear_bp, call.orientation, call.mt_chrom, call.mt_bp
        );
        let filter = if call.filters.is_empty() {
            "PASS".to_string()
        } else {
            call.filters.join(";")
        };

        let common_info = format!(
            "MTREGION={}:{}-{};MTSTRAND={};MICROHOM={};SEGMAPQ={};SEGLEN={};MTBPSPREAD={};NSPLITS={}",
            call.mt_chrom,
            call.mt_region_start + 1,
            call.mt_region_end,
            call.orientation,
            call.microhom,
            call.min_seg_mapq,
            call.min_seg_len,
            call.mt_bp_spread,
            call.max_n_splits
        );

        let (alt, info) = if call.is_full_ins && !emit_bnd {
            let alt = format!("{}{}", ref_base, ins_seq_str);
            let info = format!("SVTYPE=INS;SVLEN={};{}", call.ins_seq.len(), common_info);
            (alt, info)
        } else {
            let alt = bnd_alt(&ref_base, &call.mt_chrom, call.mt_bp + 1, call.orientation);
            let info = format!("SVTYPE=BND;SVINSSEQ={};{}", ins_seq_str, common_info);
            (alt, info)
        };

        writeln!(
            vcf,
            "{}\t{}\t{}\t{}\t{}\t.\t{}\t{}\tGT:SR\t0/1:{}",
            call.nuclear_chrom, pos, variant_id, ref_base, alt, filter, info, call.ac
        )?;
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// orchestrator
// ---------------------------------------------------------------------------

/// True if the CIGAR contains any hard clip (`H`). A hard-clipped record's stored
/// `seq()` is missing the clipped bases, so it cannot anchor `ins_seq` extraction.
fn has_hard_clip(ops: &[(u8, i64)]) -> bool {
    ops.iter().any(|(op, _)| *op == b'H')
}

/// Reconstruct a read's full alignment-segment set from a SINGLE record's own
/// alignment plus its `SA` tag. Coordinates land in that record's stored-SEQ frame
/// (`self_is_reverse` is the reference orientation), so it works whether the
/// anchoring record is the primary or a supplementary — the latter matters for
/// mt-only BAMs, where a NUMT read's nuclear primary was dropped by the chrM
/// extraction and only the chrM supplementary (with a nuclear `SA` entry) remains.
fn build_read_segments(
    self_contig: &str,
    self_is_reverse: bool,
    self_pos: i64,
    self_ops: &[(u8, i64)],
    self_mapq: u8,
    sa_str: &str,
    mt_contig: &str,
) -> Vec<AlignSeg> {
    let mut segs = vec![build_seg(
        self_contig.to_string(),
        self_contig == mt_contig,
        self_pos,
        self_ops,
        self_is_reverse,
        self_is_reverse,
        self_mapq,
    )];

    for sa in sa_str.split(';').filter(|s| !s.is_empty()) {
        let parts: Vec<&str> = sa.split(',').collect();
        if parts.len() < 6 {
            continue;
        }
        let chrom = parts[0];
        let pos_1based = match parts[1].parse::<i64>() {
            Ok(p) => p,
            Err(_) => continue,
        };
        let seg_is_reverse = parts[2] == "-";
        let ops = parse_cigar_str(parts[3]);
        let sa_mapq = parts[4].parse::<u8>().unwrap_or(0);
        segs.push(build_seg(
            chrom.to_string(),
            chrom == mt_contig,
            pos_1based - 1, // SA tag POS is 1-based
            &ops,
            seg_is_reverse,
            self_is_reverse,
            sa_mapq,
        ));
    }
    segs
}

/// Reconstruct each read's segment set (from that read's own CIGAR + SA tag) and
/// collect nuclear<->mt junctions. Iterates the whole BAM linearly so junctions are
/// found in BOTH directions (mt->nuclear and nuclear->mt) and regardless of which
/// alignment of the read is the primary. Each read is emitted once (deduped by
/// qname); the anchoring record must carry the full read SEQ (no hard clip) so the
/// inserted sequence can be resolved.
fn find_junctions(bam_file: &PathBuf, mt_contig: &str) -> Result<Vec<Junction>, Box<dyn std::error::Error>> {
    let mut bam = bam::Reader::from_path(bam_file)?;
    let header = bam.header().clone();
    let mut junctions = Vec::new();
    let mut seen: HashSet<String> = HashSet::new();
    let mut hardclip_only: HashSet<String> = HashSet::new();

    for r in bam.records() {
        let record = r?;
        // Secondary alignments carry no reliable SEQ/SA; skip. Supplementary records
        // ARE considered here (unlike before) so mt-only BAMs still resolve.
        if record.is_unmapped() || record.is_secondary() {
            continue;
        }
        // SA tag required for a split (junction-bearing) read.
        let sa_str = match record.aux(b"SA") {
            Ok(Aux::String(s)) => s.to_string(),
            _ => continue,
        };

        let qname = String::from_utf8_lossy(record.qname()).to_string();
        if seen.contains(&qname) {
            continue; // already resolved this read from another of its alignments
        }

        let self_is_reverse = record.is_reverse();
        let self_contig = String::from_utf8_lossy(header.tid2name(record.tid() as u32)).to_string();
        let self_ops: Vec<(u8, i64)> = record.cigar().iter().map(cigar_to_pair).collect();

        let segs = build_read_segments(
            &self_contig,
            self_is_reverse,
            record.pos(),
            &self_ops,
            record.mapq(),
            &sa_str,
            mt_contig,
        );

        // Need at least one mt and one non-mt segment.
        let has_mt = segs.iter().any(|s| s.is_mt);
        let has_nuc = segs.iter().any(|s| !s.is_mt);
        if !has_mt || !has_nuc {
            continue;
        }

        // ins_seq is sliced from this record's SEQ, which is complete only when the
        // record is not hard-clipped. If this alignment is hard-clipped, defer: a
        // soft-clipped alignment of the same read (e.g. the primary) may follow.
        if has_hard_clip(&self_ops) {
            hardclip_only.insert(qname);
            continue;
        }

        let seq = record.seq().as_bytes();
        junctions.extend(extract_junctions(&qname, &segs, &seq));
        seen.insert(qname);
    }

    let unresolved = hardclip_only.iter().filter(|q| !seen.contains(*q)).count();
    if unresolved > 0 {
        eprintln!(
            "[warn] {} junction-bearing read(s) had only hard-clipped alignments; \
             skipped (no full SEQ to resolve inserted sequence)",
            unresolved
        );
    }

    Ok(junctions)
}

pub fn start(
    input_bam: &PathBuf,
    chromo: &str,
    output_vcf: &PathBuf,
    reference_file: &PathBuf,
    sample_name: &str,
    config: &CallNumtsConfig,
) -> Result<(), Box<dyn std::error::Error>> {
    let junctions = find_junctions(input_bam, chromo)?;
    println!("found {} nuclear<->mt junctions", junctions.len());

    let mut calls = cluster_junctions(&junctions, config.tolerance);

    // reference-NUMT masking (optional).
    if let Some(bed_path) = &config.ref_numt_bed {
        let bed = load_bed(bed_path)?;
        for call in calls.iter_mut() {
            call.in_ref_numt = bed_overlaps(&bed, &call.nuclear_chrom, call.nuclear_bp);
        }
    }

    for call in calls.iter_mut() {
        apply_filters(call, config);
    }

    // Indexed random access (uses the .fai sidecar): O(1) memory instead of loading
    // the whole combined mt+nuclear reference into RAM, and faidx honors htslib's
    // remote/GCS backends when the crate is built with those features.
    let fasta = faidx::Reader::from_path(reference_file)?;

    write_vcf(output_vcf, &calls, &fasta, sample_name, config.emit_bnd)?;
    println!("wrote {} NUMT calls to {}", calls.len(), output_vcf.display());
    Ok(())
}

// ---------------------------------------------------------------------------
// tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn seg(contig: &str, is_mt: bool, ref_start: i64, cigar: &str, rev: bool, prim_rev: bool, mapq: u8) -> AlignSeg {
        build_seg(contig.to_string(), is_mt, ref_start, &parse_cigar_str(cigar), rev, prim_rev, mapq)
    }

    // --- align: CIGAR metrics & spans ---

    #[test]
    fn cigar_metrics_basic() {
        // 10S 90M -> ref 90, lead 10, trail 0, qlen 90
        let (rl, lead, trail, qlen) = cigar_metrics(&parse_cigar_str("10S90M"));
        assert_eq!((rl, lead, trail, qlen), (90, 10, 0, 90));
        // 10M2D5M3I4S -> ref 10+2+5=17, qlen 10+5+3=18, lead 0, trail 4
        let (rl, lead, trail, qlen) = cigar_metrics(&parse_cigar_str("10M2D5M3I4S"));
        assert_eq!((rl, lead, trail, qlen), (17, 0, 4, 18));
    }

    #[test]
    fn cigar_metrics_sums_combined_hard_and_soft_clips() {
        // Valid SAM clip runs (H outermost, S inner) must be summed, not truncated
        // to the first op: "10H20S30M40S5H" -> lead 30, trail 45.
        let (rl, lead, trail, qlen) = cigar_metrics(&parse_cigar_str("10H20S30M40S5H"));
        assert_eq!((rl, lead, trail, qlen), (30, 30, 45, 30));
    }

    #[test]
    fn forward_primary_query_span_uses_leading_clip() {
        // primary forward, 30S70M: aligned query bases 30..100 in record frame.
        let s = seg("chr1", false, 1000, "30S70M", false, false, 60);
        assert_eq!((s.q_start, s.q_end), (30, 100));
        assert_eq!((s.ref_start, s.ref_end), (1000, 1070));
    }

    #[test]
    fn sa_opposite_strand_uses_trailing_clip() {
        // primary is forward; an SA on the reverse strand (opposite) must flip to
        // the trailing clip when placed in the record frame.
        // SA cigar 70M30S (reverse): lead 0, trail 30, qlen 70 -> q_start = trail = 30? no:
        // opposite uses trail as q_start: trail=30 -> q_start=30, q_end=100.
        let s = seg("chr1", false, 2000, "70M30S", true, false, 60);
        assert_eq!((s.q_start, s.q_end), (30, 100));
    }

    #[test]
    fn hard_and_soft_clip_equivalent_for_spans() {
        let soft = seg("chrM", true, 5, "20S40M", false, false, 60);
        let hard = seg("chrM", true, 5, "20H40M", false, false, 60);
        assert_eq!((soft.q_start, soft.q_end), (hard.q_start, hard.q_end));
    }

    // --- junction: forward nuclear->mt, sequence resolution ---

    #[test]
    fn forward_nuclear_then_mt_partial_junction() {
        // Read (100bp): 0..60 aligns to chr1 (nuclear), 60..100 aligns to chrM (mt).
        // primary = nuclear forward (60M40S); SA = mt forward (60S40M).
        let nuc = seg("chr1", false, 1000, "60M40S", false, false, 60);
        let mt = seg("chrM", true, 200, "60S40M", false, false, 60);
        let seq: Vec<u8> = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\
                             CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG".to_vec();
        assert_eq!(seq.len(), 100);
        let js = extract_junctions("read1", &[nuc, mt], &seq);
        assert_eq!(js.len(), 1);
        let j = &js[0];
        assert_eq!(j.nuclear_chrom, "chr1");
        assert_eq!(j.mt_chrom, "chrM");
        // nuclear anchor = last aligned nuclear base = ref_end-1 = 1059.
        assert_eq!(j.nuclear_bp, 1059);
        // mt anchor: mt is the right/downstream, same_as_primary -> ref_start = 200.
        assert_eq!(j.mt_bp, 200);
        assert_eq!(j.orientation, '+');
        assert_eq!(j.microhom, 0);
        assert_eq!(j.gap, 0);
        assert!(!j.is_full_ins);
        // inserted mt sequence = read bases 60..100, nuclear-forward (no revcomp).
        assert_eq!(j.ins_seq, b"CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG".to_vec());
    }

    #[test]
    fn full_insertion_nuclear_mt_nuclear() {
        // Read: chr1[0..30], chrM[30..70], chr1[70..100] -> full INS.
        let nuc_l = seg("chr1", false, 1000, "30M70S", false, false, 60);
        let mt = seg("chrM", true, 500, "30S40M30S", false, false, 60);
        let nuc_r = seg("chr1", false, 1030, "70S30M", false, false, 60);
        let seq = vec![b'A'; 100];
        let js = extract_junctions("read2", &[nuc_l, mt, nuc_r], &seq);
        assert_eq!(js.len(), 1);
        assert!(js[0].is_full_ins);
        assert_eq!(js[0].ins_seq.len(), 40);
        assert_eq!(js[0].n_splits, 3);
    }

    #[test]
    fn reverse_nuclear_revcomps_inserted_sequence() {
        // primary = mt forward; nuclear SA on opposite strand -> ins_seq revcomp'd.
        // Read frame (record frame = mt primary): mt 0..40 (40M60S).
        let mt = seg("chrM", true, 300, "40M60S", false, false, 60);
        // nuclear SA opposite strand: cigar 40S60M reverse -> opposite uses trail(=0)
        // as q_start -> q_start=0, q_end=60. That collides with mt (0..40); to place
        // nuclear after mt in record frame we instead give it a leading clip of 40
        // on the reverse read: cigar 60M40S reverse -> trail=40 -> q_start=40.
        let nuc = seg("chr2", false, 8000, "60M40S", true, false, 60);
        let seq: Vec<u8> = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".iter().cloned()
            .chain(std::iter::repeat(b'T').take(60)).collect();
        assert_eq!(seq.len(), 100);
        let js = extract_junctions("read3", &[mt, nuc], &seq);
        assert_eq!(js.len(), 1);
        // mt segment record-frame bases 0..40; nuclear is opposite -> revcomp applied.
        let expected = revcomp(b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");
        assert_eq!(js[0].ins_seq, expected);
        assert_eq!(js[0].orientation, '-');
    }

    // --- find_junctions helpers: mt-only-BAM recovery ---

    #[test]
    fn has_hard_clip_detects_h() {
        assert!(has_hard_clip(&parse_cigar_str("11831H925M")));
        assert!(!has_hard_clip(&parse_cigar_str("11848S459=1D466=")));
        assert!(!has_hard_clip(&parse_cigar_str("100M")));
    }

    #[test]
    fn build_read_segments_recovers_junction_from_mt_supplementary() {
        // Real mt-only-BAM shape (HG002 chrM_100x): the read's nuclear PRIMARY on
        // chr3 was dropped by the chrM extraction; only the chrM SUPPLEMENTARY
        // remains, and it is soft-clipped (full SEQ) with the nuclear primary in SA.
        //   chrM self:  11848S459=1D466=   (mt aligned, forward)
        //   SA:         chr3,171111165,+,11831M17I925S,60,1
        let self_ops = parse_cigar_str("11848S459=1D466=");
        let sa = "chr3,171111165,+,11831M17I925S,60,1;";
        let segs = build_read_segments("chrM", false, 5000, &self_ops, 60, sa, "chrM");
        assert_eq!(segs.len(), 2);
        assert!(segs.iter().any(|s| s.is_mt && s.contig == "chrM"));
        assert!(segs.iter().any(|s| !s.is_mt && s.contig == "chr3"));

        // Full read length = 11848 (nuclear query) + 925 (mt query) = 12773.
        let seq = vec![b'A'; 12773];
        let js = extract_junctions("numt_read", &segs, &seq);
        assert_eq!(js.len(), 1);
        assert_eq!(js[0].nuclear_chrom, "chr3");
        assert_eq!(js[0].mt_chrom, "chrM");
        assert!(!js[0].is_full_ins);
        assert_eq!(js[0].ins_seq.len(), 925); // mt-derived bases resolved
    }

    // --- cluster ---

    fn mkj(qname: &str, nuc_bp: i64, mt_bp: i64) -> Junction {
        Junction {
            qname: qname.to_string(),
            nuclear_chrom: "chr1".to_string(),
            nuclear_bp: nuc_bp,
            mt_chrom: "chrM".to_string(),
            mt_bp,
            mt_ref_start: mt_bp,
            mt_ref_end: mt_bp + 40,
            orientation: '+',
            microhom: 0,
            gap: 0,
            ins_seq: vec![b'A'; 40],
            is_full_ins: false,
            n_splits: 2,
            read_len: 100,
            min_seg_mapq: 60,
            min_seg_len: 40,
        }
    }

    #[test]
    fn distinct_nearby_insertions_stay_separate() {
        // Two insertions 5kb apart must NOT merge (regression vs old 10kb gap merge).
        let js = vec![
            mkj("r1", 1000, 200),
            mkj("r2", 1005, 202),
            mkj("r3", 6000, 900),
        ];
        let calls = cluster_junctions(&js, 100);
        assert_eq!(calls.len(), 2);
    }

    #[test]
    fn jittered_breakpoints_merge_and_count_distinct_reads() {
        let js = vec![
            mkj("r1", 1000, 200),
            mkj("r2", 1050, 210), // within eps=100
            mkj("r1", 1010, 205), // same read name r1 -> AC must not double count
        ];
        let calls = cluster_junctions(&js, 100);
        assert_eq!(calls.len(), 1);
        assert_eq!(calls[0].ac, 2); // r1, r2
    }

    // --- filters ---

    fn mkcall(min_mapq: u8, min_seg_len: i64, ac: usize, spread: i64, in_ref: bool) -> NumtCall {
        NumtCall {
            nuclear_chrom: "chr1".to_string(),
            nuclear_bp: 1000,
            mt_chrom: "chrM".to_string(),
            mt_bp: 200,
            mt_region_start: 200,
            mt_region_end: 240,
            orientation: '+',
            ac,
            ins_seq: vec![b'A'; 40],
            is_full_ins: true,
            microhom: 0,
            min_seg_mapq: min_mapq,
            min_seg_len,
            mt_bp_spread: spread,
            max_n_splits: 2,
            max_split_read_len: 10000,
            in_ref_numt: in_ref,
            filters: Vec::new(),
        }
    }

    #[test]
    fn filters_pass_and_tags() {
        let cfg = CallNumtsConfig::default();
        let mut good = mkcall(60, 200, 3, 0, false);
        apply_filters(&mut good, &cfg);
        assert!(good.filters.is_empty());

        let mut lowmapq = mkcall(5, 200, 3, 0, false);
        apply_filters(&mut lowmapq, &cfg);
        assert!(lowmapq.filters.contains(&"LOWMAPQ".to_string()));

        let mut shortseg = mkcall(60, 50, 3, 0, false);
        apply_filters(&mut shortseg, &cfg);
        assert!(shortseg.filters.contains(&"SHORTSEG".to_string()));

        let mut discordant = mkcall(60, 200, 3, 500, false);
        apply_filters(&mut discordant, &cfg);
        assert!(discordant.filters.contains(&"DISCORDANT".to_string()));

        let mut refnumt = mkcall(60, 200, 3, 0, true);
        apply_filters(&mut refnumt, &cfg);
        assert!(refnumt.filters.contains(&"REFNUMT".to_string()));
    }

    #[test]
    fn multisplit_uses_real_read_length_not_insseq_len() {
        let cfg = CallNumtsConfig::default(); // max_splits=3, per_kb=0.1
        // A 10 kb read tiled into 4 segments is within allowance (3 + round(0.1*10) = 4).
        // The old proxy used ins_seq.len() (40 bp) -> allowance 3 -> wrongly flagged.
        let mut c = mkcall(60, 200, 3, 0, false);
        c.max_n_splits = 4;
        c.max_split_read_len = 10_000;
        apply_filters(&mut c, &cfg);
        assert!(!c.filters.contains(&"MULTISPLIT".to_string()));

        // A short 500 bp read shredded into 5 segments exceeds allowance (3 + 0 = 3).
        let mut chimera = mkcall(60, 200, 3, 0, false);
        chimera.max_n_splits = 5;
        chimera.max_split_read_len = 500;
        apply_filters(&mut chimera, &cfg);
        assert!(chimera.filters.contains(&"MULTISPLIT".to_string()));
    }

    #[test]
    fn low_support_triggers_below_threshold() {
        let cfg = CallNumtsConfig::default(); // ac_threshold = 2
        let mut c = mkcall(60, 200, 3, 0, false);
        c.ac = 1;
        apply_filters(&mut c, &cfg);
        assert!(c.filters.contains(&"LOW_SUPPORT".to_string()));
    }
}
