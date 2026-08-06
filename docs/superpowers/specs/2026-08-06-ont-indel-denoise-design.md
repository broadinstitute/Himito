# ONT Small-Indel Denoiser (`Himito denoise --indels`) — Design

**Date:** 2026-08-06
**Status:** Approved (design), pending implementation plan
**Scope:** ONT-only, small-indel (< 5 bp) read error correction, added to the existing `denoise` pass. Opt-in; PacBio and default ONT behavior unchanged.
**Supersedes:** the "Indels / homopolymer correction" non-goal in [`2026-08-03-ont-bam-denoise-design.md`](2026-08-03-ont-bam-denoise-design.md) §2 / §11.

## 1. Motivation

The v1 denoiser corrects substitutions only. ONT's *dominant* residual error mode after that is small indels — homopolymer run-length errors and short tandem-repeat length errors — and those errors are not neutral to Himito's architecture:

- `build.rs` consolidates reads into anchor edges by **exact inter-anchor sequence identity**, and it k-merizes **raw read `SEQ`, ignoring `CIGAR` entirely** (`src/build.rs:920–930`). A 1 bp indel error therefore changes the read's inter-anchor sequence exactly as destructively as a substitution did before v1: the read shatters off into its own edge, and can additionally break the anchor k-mers themselves.
- The resulting spurious edges become spurious columns in the read × variant matrix, degrading haplotype clustering and the `lineage` / SCITE tree.

The existing mitigation is **post-hoc and downstream**: `call.rs:459–469` flags small indels (< 5 bp) whose frequency falls below `indel_false_threshold` as `Potential_Artifact` in the VCF. That labels a bad call after the fact; it does nothing for the graph or the matrix, which are already polluted by the time `call` runs.

**Target failure modes for this work (explicitly, per requirements):** messy graph / spurious bubbles, and a noisy lineage matrix. Reducing VCF false-indel calls is a welcome side effect, not the objective.

## 2. Goal & scope

Extend the existing `denoise` pileup pass to also model **small indel events per reference site** and rewrite reads to their maximum-a-posteriori indel allele, producing a valid denoised BAM that `build` → `call` → `lineage` consume unchanged.

**In scope:**
- ONT data types only; PacBio remains a no-op passthrough.
- Indel events with length < `--indel-max-len` (default 5).
- **Two-directional (symmetric) correction**: a read may lose a spurious indel, gain a consensus indel it dropped, or be snapped from one indel allele to another.
- Same single pileup pass as SNV denoising — no second BAM read.
- Opt-in behind `--indels` (default **off**), so no existing result changes on landing.

**Non-goals:**
- Indels ≥ `--indel-max-len` (treated as untouchable; they may be real SVs or NUMT signal).
- Local realignment / POA consensus (see §11).
- Full EM over the indel allele set (frequencies are raw proportions; see §3.4).
- Any change to methylation, assembly, `call`, `lineage`, or PacBio behavior.
- Re-sorting or re-indexing the output (proven unnecessary — see the invariant in §4).

## 3. Method

### 3.1 Why not the alternatives

| Approach | Verdict |
|---|---|
| **A. Read-level revert/reassign in `denoise.rs`** | **Chosen.** Fixes the graph at its source, since `build.rs` k-merizes raw `SEQ`. Cost: a length-changing BAM writer, the only delicate part. |
| B. Graph-level bubble collapse in `build.rs` | Rejected. Far less code (an extension of `min_edge_reads` + `normalize_homopolymer_indels`), but it fires *after* anchors have been fragmented — it cleans the symptom, and cannot help the matrix, where the read's identity across the bubble is already lost. |
| C. Local realignment / POA | Rejected for now. Most powerful and the standard answer, but a large new subsystem with an aligner in the loop, and its aggressiveness is hard to bound — real heteroplasmic indels are exactly what it tends to smooth away. |

### 3.2 Observation collection

Inside the existing pileup loop, alongside the SNV tally, read `Alignment::indel()` (`rust_htslib` 0.49 `bam::pileup::Indel`) for every primary, mapped alignment at column `p`:

| `indel()` | Recorded as |
|---|---|
| `Ins(n)`, `n < indel_max_len` | event `INS(seq = read[q+1 .. q+1+n])`, anchored at `p` |
| `Del(n)`, `n < indel_max_len` | event `DEL(n, ref[p+1 .. p+1+n])`, anchored at `p` |
| either, `n ≥ indel_max_len` | **untouchable** — the read counts toward depth but is never corrected at this site |
| `None` | depth only (this read votes REF at this junction) |

`indel_max_len` is an **exclusive** bound: with the default of 5, events of length 1–4 are correctable and length ≥ 5 is untouchable.

`indel()` is reported at the anchor column — the last reference base *before* the event — where `qpos()` is always `Some`. Column depth is already computed by the SNV pass (`obs.len()`) and is stored into a `pos → depth` map for use as the VAF denominator.

**Required change to the existing column skip.** `compute_corrections` currently does `if counts.non_ref() == 0 { columns_skipped += 1; continue; }` — that gate is about non-reference *substitutions*, and a column with a perfectly clean base composition can still carry indel events. When `--indels` is on, indel observation collection must happen **before** that early-`continue`, or every indel at a substitution-clean column is silently lost. This is the single easiest way to get a quietly under-corrected result, so it is called out here explicitly.

### 3.3 Left-normalization (and why decisions move out of the loop)

An indel inside a homopolymer or tandem repeat can be legally anchored at several reference positions, so the same biological event arrives from different reads at different columns. Each event is shifted maximally left — *while the base immediately preceding the event equals the event's last base, shift by one* — yielding a canonical key:

```
(tid, norm_pos, kind, len, ins_seq)
```

A normalized event can land in a column the pileup loop has already passed. Sites are therefore accumulated across the whole pass and **decided after the pileup loop closes**, then merged with the in-loop SNV corrections before the write-back pass. This is the one structural difference from the SNV design, which decides per column in-loop.

Under-normalization splits one site's support across positions, deflating its VAF and causing over-correction. The normalization unit tests (§8) are the guard.

### 3.4 Site model

Alleles at a normalized site are indexed by net length change `ℓ` — `REF = 0`, `INS(n) = +n`, `DEL(n) = −n` — plus inserted sequence for insertions.

**Context.** `L` is the reference repeat context at `norm_pos`: the homopolymer run length for 1 bp events, the tandem-copy count of the event's own sequence for longer ones, and `L = 1` in unique sequence.

**One context function drives both the gate and the assignment:**

```
ε(L)     = min(indel_err0 · indel_err_scale^(L−1), indel_err_cap)   # per-junction error prob
floor(L) = max(indel_vaf, indel_floor_mult · ε(L))                  # candidacy VAF floor
```

**Candidacy** (an allele is *kept*, i.e. eligible to receive reads, iff):

1. `fwd ≥ min_strand && rev ≥ min_strand` — unchanged from the SNV design.
2. `vaf = support / depth(norm_pos) ≥ floor(L)`. If `norm_pos` has no recorded column depth (it normalized past the first pileup column of the contig), the depth falls back to the maximum depth among the event's original anchor columns; if that is also zero, the site is left undecided and every read at it is untouched.
3. `call::strand_bias_pvalue(fwd, rev, column_fwd_frac) ≥ strand_bias_p`, with the same near-homoplasmic exemption at `homoplasmic_vaf` — a single-strand artifact cannot reach a near-homoplasmic frequency, so testing one would correct away a real variant.

**REF is always kept**, mirroring "the reference allele is always eligible" in the SNV model.

Allele frequencies `f` are raw proportions renormalized over the kept set — exactly what `fit_site` falls back to. This is deliberately *not* the full EM: upgrading to EM over the dynamic allele set is a drop-in later (§11).

So a 1 bp deletion inside a 7-mer A-run must clear a far higher bar than the same deletion in unique sequence, which is where ONT's error mass actually sits.

### 3.5 Per-read MAP assignment (two-directional)

For a read observing allele `a` at the site:

```
assign(a) = argmax over kept k of  f_k · P(a | k)

P(a | k) = 1 − ε(L)                            if a == k
         = ε(L) · indel_delta^|ℓ_a − ℓ_k| / Z  otherwise

Z = Σ over kept a' ≠ k of indel_delta^|ℓ_a' − ℓ_k|
```

`Z` normalizes the error mass over the kept alleles other than `k`, so `P(· | k)` sums to 1 over the kept set. When `k` is the only kept allele the site is degenerate and no read is reassigned.

**The per-read likelihood term is mandatory, not decorative.** Frequency-only MAP (`argmax_k f_k`) is identical for every read at the site, so *every* read collapses to the single most frequent allele — at a 40/60 heteroplasmic indel site the 40% allele is annihilated. That is structural, not a tuning problem. The SNV `map_allele` escapes it only because base quality supplies a strong per-observation likelihood; indels have no per-base quality, so `ε(L)` supplies it instead.

The crossover is analytic and conservative: a read leaves its own allele only when `f_own / f_other < ε / (1 − ε)`. At a long-homopolymer site (`ε = 0.2`) that means its allele is below ~20% of the competing one; in unique sequence (`ε = 0.01`), below ~1%.

| Site | Outcome |
|---|---|
| 40% INS / 60% REF, `ε = 0.2` | nothing flips — heteroplasmy preserved |
| 90% INS / 10% REF, `ε = 0.2` | REF reads flip to INS — "read dropped a real indel" repaired |
| failed 2 bp DEL where the site keeps a 1 bp DEL | read snaps to the 1 bp allele, not all the way to REF |

The third row is why this design is two-directional rather than revert-only: a one-directional model reverts the 2 bp read all the way to reference, leaving it disagreeing with consensus in the opposite direction.

**Guards.** A read is reassignable at a site only if it spans the site's full footprint with plain `M`/`=`/`X` ops plus `--indel-flank` bases of margin on each side (default 5). Reads soft-clipped at the site are skipped. Skipped reads are counted (`reads_skipped_span_guard`).

## 4. Write-back

### 4.1 Edit representation

Mixing position-preserving edits (substitutions, keyed by `qpos`) with length-changing ones is where this gets bug-prone. Both edit kinds already know their `ref_pos` at record time, so the correction payload is unified on **reference coordinates** and the write-back becomes a single left-to-right walk with no `qpos` arithmetic:

```rust
enum Allele { Ref, Ins(Vec<u8>), Del(u32) }

struct IndelEdit { ref_pos: u32, from: Allele, to: Allele }

struct ReadEdits { subs: Vec<(u32, u8)>, indels: Vec<IndelEdit> }  // both keyed by ref_pos
pub type Corrections = HashMap<Vec<u8>, ReadEdits>;
```

Substitutions are only ever recorded at `M`/`=`/`X` columns, so `ref_pos → qpos` is unique and the walk resolves it. This retires the current `if qpos < seq.len()` guard in `apply_corrections`, which exists only because `qpos` was carried unanchored.

Changing the `Corrections` type touches existing v1 tests; they are updated as part of this work.

### 4.2 The walk

One pass per modified read over the original CIGAR, maintaining `(ref_pos, qpos)` and emitting new `SEQ`, `QUAL`, and `CIGAR` together:

| Op | Action |
|---|---|
| `M` / `=` / `X` | copy; apply substitution at matching `ref_pos`; `REF→INS(s)` appends `s` after the anchor base and emits `I(len(s))`; `REF→DEL(n)` skips `n` query bases and emits `D(n)` |
| `I(n)` | reverted → drop; replaced → emit the new inserted bases; else copy |
| `D(n)` | reverted → splice `n` reference bases into `SEQ`/`QUAL`, emit `M(n)`; replaced → splice/trim to the target length |
| `S` / `H` / `N` / `P` | copied verbatim; no edit is ever anchored inside one |

Adjacent same-kind ops coalesce.

### 4.3 The invariant

**Every transition preserves the read's reference consumption:**

- `REF→INS(n)` adds `I(n)`, which consumes query but not reference.
- `REF→DEL(n)` converts `M(n)→D(n)`; both consume `n` reference bases.
- Reverting `I(n)` removes query only.
- Reverting `D(n)→M(n)`; both consume `n` reference bases.

Therefore `POS`, reference end, and coordinate-sort order are invariant in all four transition types — **no re-sorting, no index invalidation**, and the CIGAR rewrite stays purely local. This becomes a debug assertion inside the walk and an explicit assertion in every rewrite test.

Inserted sequence is always the site's consensus inserted bases, never a per-read guess.

### 4.4 Synthesized base qualities

Bases added by the walk (reverting a deletion, or adding consensus `INS` bases to a REF read) receive `min(qual[left flank], qual[right flank])` — never claiming more confidence than the surrounding read. At a read edge, the single available flank is used.

### 4.5 Auxiliary tags

For any read whose **length changed**:

- Strip `MM`, `ML`, `MN`. Base-modification tags are indexed by the read's own base positions and their semantics are broken by a length change. This is safe because `methyl::start` runs on the **original** mt BAM (`src/main.rs:594`), consistent with the v1 methylation-safety decision.
- Strip `NM`, `MD`, `cs` as stale.

All other tags are preserved, as in v1.

**Pre-existing issue, flagged separately:** substitution-only reads already leave `NM`/`MD`/`cs` stale today, and a substitution that creates or destroys a `C` also breaks `MM`. That is a latent v1 issue, not one this change introduces. Stripping those tags on substitution-only modified reads too is a small hardening; it is **out of scope here unless explicitly requested**.

## 5. CLI & pipeline integration

`start()` already takes 9 arguments; six more would be unreadable. Indel configuration is passed as one struct:

```rust
pub struct IndelOpts {
    pub enabled: bool,      // --indels                 default false
    pub max_len: u32,       // --indel-max-len          default 5
    pub vaf: f64,           // --indel-vaf              default 0.05
    pub err0: f64,          // --indel-err0             default 0.01
    pub err_scale: f64,     // --indel-err-scale        default 1.5
    pub err_cap: f64,       // --indel-err-cap          default 0.4
    pub floor_mult: f64,    // --indel-floor-mult       default 3.0
    pub delta: f64,         // --indel-delta            default 0.3
    pub flank: usize,       // --indel-flank            default 5
}
impl Default for IndelOpts { /* enabled: false */ }
```

`min_strand`, `strand_bias_p`, and `homoplasmic_vaf` are reused from the existing SNV parameters.

`QuickStart` passes `IndelOpts::default()` (disabled), so **no existing QuickStart result shifts when this lands**. The default is flipped only after the validation gate in §9 passes.

**The defaults above are starting points to be tuned on `lineage_sim`, not claims.** `ε` governs how aggressive assignment is and is the one knob that can erase real heteroplasmy.

## 6. Statistics

`DenoiseStats` gains additive fields (JSON stays backward-compatible):

- `indel_sites_examined`, `indel_sites_kept`, `indel_sites_strand_biased`
- `indel_events_examined`
- `reads_reassigned_ref_to_ins`, `reads_reassigned_ref_to_del`, `reads_reassigned_indel_to_ref`, `reads_reassigned_indel_to_indel`
- `indel_bases_inserted`, `indel_bases_removed`
- `reassignments_by_hp_length: Vec<u64>` — audit of how much correction fired inside long homopolymers
- `reads_skipped_span_guard`

## 7. Performance

Unchanged asymptotics: still `O(B)` in total aligned bases, with no pairwise read comparison. The indel site map adds memory proportional to the number of distinct normalized events (bounded in practice by a small multiple of reference length for mito), and the deferred site decision adds one pass over that map after the pileup loop. The write-back walk is `O(read length)` per modified read.

## 8. Testing

In-file unit tests, matching the existing `denoise.rs` style:

- **Normalization** — one homopolymer insertion anchored at three different columns collapses to a single key.
- **Context detection** — `L` for a homopolymer run, a tandem repeat, and unique sequence.
- **`ε(L)` / `floor(L)`** — monotonic in `L` and correctly capped.
- **Gate** — a strand-balanced 20% 1 bp `DEL` is kept in unique sequence and rejected inside a 7-mer run; a strand-skewed candidate passing `min_strand` is rejected by the binomial test; the near-homoplasmic exemption spares a high-frequency allele.
- **Assignment** — the three rows of the §3.5 table: 40/60 preserves both, 90/10 flips REF reads, and a failed 2 bp `DEL` snaps to the kept 1 bp `DEL`.
- **Rewrite, all four transitions** — `SEQ`/`QUAL`/`CIGAR` correct, **with the reference-consumption invariant asserted in each**.
- **Multi-edit read** — one read carrying a substitution, an insertion edit, and a deletion edit applied in a single walk.
- **Span/flank guard** — a read ending 2 bp past a site is skipped and counted.
- **Tags** — `MM`/`ML` stripped only on length-changed reads; other tags survive.
- **Substitution-clean column carrying an indel** — a column whose bases are 100% reference but where reads carry a spurious insertion; the event must still be collected and corrected, proving the `non_ref() == 0` early-`continue` no longer swallows indels (§3.2).
- **End-to-end** — `compute_corrections` + `apply_corrections` on a synthetic BAM, output re-readable by htslib and still coordinate-sorted.
- **Off by default** — with `IndelOpts::default()`, output is byte-identical to current v1 behavior.

## 9. Validation

The stated goals are graph and matrix cleanliness, so those are measured directly on `lineage_eval/lineage_sim` (pbsim3 ONT profiles, known truth tree), `--indels` off vs on:

1. **GFA edge and anchor count** — fragmentation proxy.
2. **Matrix column count vs. true variant count** — false columns are the noise being targeted.
3. **`score_lineage.py` tree score** — the end metric.
4. **Small-indel precision / recall vs. truth** — confirms real indels survived.

**Ship gate:** a measurable drop in (1) and (2) with **no regression** in (3) or (4). Only then is the `--indels` default flipped on for ONT.

## 10. Risks

- **Parameter defaults are unvalidated guesses** until the §9 sweep runs. `ε` is the knob that can erase real heteroplasmy.
- **Under-normalization splits support** across positions, deflating VAF and over-correcting. Guarded by the normalization tests.
- **Single-pass by construction.** Frequencies are computed from the original BAM and decisions are applied once. Re-running the denoiser on its own output would cascade (each pass shifts the frequencies the next pass sees). This must be documented in the user guide.
- **Two-directional correction touches reads that carried no error signal at all** (a REF-observing read gaining an indel). The `ε`-derived crossover bounds this, and `reassignments_by_hp_length` makes it auditable.

## 11. Future work

- Full EM over the dynamic indel allele set, replacing raw renormalized proportions (drop-in at §3.4).
- Empirical per-run estimation of `ε(L)` from the BAM itself (clearly homoplasmic-reference sites), self-calibrating across chemistries and basecallers.
- Local realignment / POA consensus (approach C) if the site-wise model proves insufficient in low-complexity regions.
- Extend `MM`/`ML` remapping across length-changing edits, if a downstream consumer of the denoised BAM ever needs methylation.
- Fold the deferred v1 tag-staleness hardening (§4.5) into the substitution path.
