version 1.0

# WDL port of run_simulation.sh: a Cartesian sweep of the simulate -> reconstruct
# -> score cycle (run_eval.sh) over seeds x n-mutations x depths. Each grid cell
# runs one full RunEvalCell task; CombineMetrics concatenates the per-cell metrics
# (with identifying seed/n_mutations/total_depth columns) into a single metrics.tsv.
#
# ---------------------------------------------------------------------------
# THE THREE FREQUENCY GATES
# ---------------------------------------------------------------------------
# The pipeline applies three independent frequency floors, and a truth variant is
# destroyed by whichever is highest:
#
#   denoise --vaf   (before the graph; follows --vaf when unset)
#   call -v         (= min_hf)
#   call -f         (= frequency_threshold; Himito's own default is 0.2)
#
# EVERY one of them must sit below the LOWEST TRUTH HF, or truth variants are
# deleted before they can be called and var_recall is capped no matter how much
# depth is used -- these are frequency floors, not coverage floors.
#
# Note that is the lowest *truth* HF, not sim_min_hf. sim_min_hf is only the floor
# the simulator may descend to; the achieved minimum is usually well above it (the
# standard chain config floors at 0.296 under sim_min_hf=0.05), so comparing gates
# against sim_min_hf rejects perfectly good configs. ValidateConfig computes the
# exact minimum for topology=chain and falls back to sim_min_hf as a lower bound
# otherwise, then fails the workflow before a single expensive cell runs.
#
# NOTE: topology/denoise_vaf and the HP modes require an image built from the
# current lineage_sim/ scripts, including hp_tracts.py.

workflow HimitoLineageSim {
    input {
        Array[Int] seeds
        Array[Int] n_mutations
        Array[Int] depths
        File reference_fa

        # ont-denoised | ont-r10 | hifi. ont-r10 SKIPS the denoise step, which
        # leaves indel artifacts crowding the graph and yields an empty VCF;
        # hifi needs ccs -> linux/amd64.
        String profile = "ont-denoised"

        # Himito call -f. Himito's own default is 0.2, which is above essentially
        # every simulated truth variant and empties the VCF. Must be < sim_min_hf.
        Float frequency_threshold = 0.05

        # Himito denoise --vaf. Unset follows --vaf (= min_hf). Must be
        # < sim_min_hf when set.
        Float? denoise_vaf

        Float? fp
        Float? fn

        Float min_hf = 0.01
        Float max_hf = 0.95

        # Band for the SIMULATED truth frequencies -- a different thing from
        # min_hf/max_hf, which filter what was called.
        Float sim_min_hf = 0.05
        Float sim_max_hf = 0.99

        # Tree shape. "chain" is the knob that raises mutation-ORDERING difficulty
        # without costing detection (a unary path never splits clone mass between
        # siblings, so its minimum HF is higher than a random tree of the same
        # size). "star" has no ancestral pairs and is a control. Raising
        # n_mutations is NOT a difficulty knob: it only shrinks clones and costs
        # variant recall.
        String topology = "random"

        # Terminal mass each internal node keeps. THIS is what makes an individual
        # ordering call hard; lower is harder, and it costs no detection.
        Float internal_keep = 0.20

        # Homopolymer length heteroplasmy; all off by default so an
        # unchanged inputs JSON reproduces the SNV-only benchmark.
        String background = "none"          # none | hvs1-polyC | comma list e.g. "T16189C,A16183C"
        Int hp_lineage_indels = 0           # K of the n_mutations edges are HP length changes
        String hp_unstable = ""             # "" | common | "D310:0=0.5,1=0.5;C573:..."
        Float hp_min_subclone_freq = 0.002
        Boolean collapse_hp = false         # Himito lineage --collapse-hp
        Float? hp_fp_rate                   # Himito lineage --hp-fp-rate (unset = 2x fp)
        Float? hp_fn_rate                   # Himito lineage --hp-fn-rate (unset = 2x fn)

        String docker = "us.gcr.io/broad-dsp-lrma/hangsuunc/himito-lineage-sim:dev"
        RuntimeAttr? runtime_attr_override
    }

    String fp_arg = if defined(fp) then "--fp " + select_first([fp]) else ""
    String fn_arg = if defined(fn) then "--fn " + select_first([fn]) else ""
    String denoise_vaf_arg = if defined(denoise_vaf) then "--denoise-vaf " + select_first([denoise_vaf]) else ""
    # Empty when the mode is off, same pattern as fp_arg.
    # Quote specs that may contain ';' so the shell does not split them.
    String background_arg = if background == "none" then "" else "--background \"" + background + "\""
    String hp_lineage_arg = if hp_lineage_indels > 0 then "--hp-lineage-indels " + hp_lineage_indels else ""
    String hp_unstable_arg = if hp_unstable != "" then "--hp-unstable \"" + hp_unstable + "\"" else ""
    String hp_min_subclone_arg = if hp_unstable != "" then "--hp-min-subclone-freq " + hp_min_subclone_freq else ""
    String collapse_hp_arg = if collapse_hp then "--collapse-hp" else ""
    String hp_fp_arg = if defined(hp_fp_rate) then "--hp-fp-rate " + select_first([hp_fp_rate]) else ""
    String hp_fn_arg = if defined(hp_fn_rate) then "--hp-fn-rate " + select_first([hp_fn_rate]) else ""
    String hp_unstable_txt = if hp_unstable == "" then "NA" else hp_unstable
    String hp_fp_txt = if defined(hp_fp_rate) then "" + select_first([hp_fp_rate]) else "NA"
    String hp_fn_txt = if defined(hp_fn_rate) then "" + select_first([hp_fn_rate]) else "NA"
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, object {}])

    # Fail the whole workflow on a misconfigured gate before any cell burns a VM.
    call ValidateConfig {
        input:
            profile = profile,
            frequency_threshold = frequency_threshold,
            denoise_vaf = denoise_vaf,
            min_hf = min_hf,
            sim_min_hf = sim_min_hf,
            topology = topology,
            internal_keep = internal_keep,
            n_mutations = n_mutations,
            background = background,
            hp_unstable = hp_unstable,
            hp_lineage_indels = hp_lineage_indels,
            collapse_hp = collapse_hp,
            hp_fp_set = defined(hp_fp_rate),
            hp_fn_set = defined(hp_fn_rate),
            docker = docker
    }

    scatter (seed in seeds) {
        scatter (nmut in n_mutations) {
            scatter (depth in depths) {
                call RunEvalCell {
                    input:
                        seed = seed,
                        n_mutations = nmut,
                        total_depth = depth,
                        profile = profile,
                        fp_arg = fp_arg,
                        fn_arg = fn_arg,
                        denoise_vaf_arg = denoise_vaf_arg,
                        frequency_threshold = frequency_threshold,
                        min_hf = min_hf,
                        max_hf = max_hf,
                        sim_min_hf = sim_min_hf,
                        sim_max_hf = sim_max_hf,
                        topology = topology,
                        internal_keep = internal_keep,
                        background = background,
                        hp_lineage_indels = hp_lineage_indels,
                        hp_unstable = hp_unstable,
                        collapse_hp = collapse_hp,
                        hp_unstable_txt = hp_unstable_txt,
                        hp_fp_txt = hp_fp_txt,
                        hp_fn_txt = hp_fn_txt,
                        background_arg = background_arg,
                        hp_lineage_arg = hp_lineage_arg,
                        hp_unstable_arg = hp_unstable_arg,
                        hp_min_subclone_arg = hp_min_subclone_arg,
                        collapse_hp_arg = collapse_hp_arg,
                        hp_fp_arg = hp_fp_arg,
                        hp_fn_arg = hp_fn_arg,
                        reference_fa = reference_fa,
                        validated = ValidateConfig.ok,
                        docker = docker,
                        runtime_attr_override = runtime_attr
                }
            }
        }
    }

    Array[File] cell_metrics = flatten(flatten(RunEvalCell.cell_metrics))
    Array[File] cell_bundles = flatten(flatten(RunEvalCell.cell_bundle))
    Array[File] cell_truth = flatten(flatten(RunEvalCell.truth_dir_tar))
    Array[File] cell_matrices = flatten(flatten(RunEvalCell.matrix))
    Array[File] cell_vcfs = flatten(flatten(RunEvalCell.vcf))
    Array[File] cell_cleaned = flatten(flatten(RunEvalCell.cleaned_matrix))

    call CombineMetrics {
        input:
            cell_metrics = cell_metrics,
            docker = docker
    }

    output {
        File metrics = CombineMetrics.combined
        Array[File] per_cell_metrics = cell_metrics
        Array[File] per_cell_bundles = cell_bundles
        Array[File] per_cell_truth = cell_truth
        Array[File] per_cell_matrix = cell_matrices
        Array[File] per_cell_vcf = cell_vcfs
        Array[File] per_cell_cleaned_matrix = cell_cleaned
    }
}

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
    String? docker
}

task ValidateConfig {

    meta {
        description: "Check the frequency-gate invariant and topology feasibility before the grid runs. Cheap, and it turns a whole sweep of empty VCFs into one clear error."
    }

    input {
        String profile
        Float frequency_threshold
        Float? denoise_vaf
        Float min_hf
        Float sim_min_hf
        String topology
        Float internal_keep
        Array[Int] n_mutations
        String background
        String hp_unstable
        Int hp_lineage_indels
        Boolean collapse_hp
        Boolean hp_fp_set
        Boolean hp_fn_set
        # The HP-rate sweep sets this. The sim workflow leaves it false, so an
        # SNV-only inputs JSON still validates.
        Boolean require_hp_mode = false
        String docker
    }

    # Unset denoise_vaf follows --vaf, which run_eval.sh sets from min_hf.
    Float effective_denoise_vaf = select_first([denoise_vaf, min_hf])

    command <<<
        set -euo pipefail
        python3 <<'PY'
import sys
sys.path.insert(0, "/opt/lineage_sim")
try:
    import hp_tracts
    import simulate_tree
except ImportError as exc:
    print(
        "ERROR: cannot import hp_tracts/simulate_tree from /opt/lineage_sim "
        f"({exc}). Rebuild the himito-lineage-sim image.",
        file=sys.stderr,
    )
    sys.exit(1)

profile   = "~{profile}"
freq      = float("~{frequency_threshold}")
denoise   = float("~{effective_denoise_vaf}")
min_hf    = float("~{min_hf}")
sim_min   = float("~{sim_min_hf}")
topology  = "~{topology}"
keep      = float("~{internal_keep}")
nmuts     = [int(x) for x in "~{sep=' ' n_mutations}".split()]
# simulate_tree.parse_hp_unstable checks tract ids against hp_tracts.CATALOGUE.
background = "~{background}"
hp_unstable = "~{hp_unstable}"
hp_k = int("~{hp_lineage_indels}")
collapse_hp = "~{collapse_hp}".lower() == "true"
hp_fp_set = "~{hp_fp_set}".lower() == "true"
hp_fn_set = "~{hp_fn_set}".lower() == "true"
# Background alone, or collapse_hp alone, does not give the HP rates a variant
# to act on. The sweep is then guaranteed flat.
require_hp_mode = "~{require_hp_mode}".lower() == "true"
REF_FRACTION = 0.15   # simulate_tree.py default; not exposed by run_eval.sh

errors, warnings = [], []

if topology not in ("random", "chain", "star"):
    errors.append(f"topology must be random|chain|star, got {topology!r}")

# --- the gate invariant -------------------------------------------------------
# Compare gates against the LOWEST TRUTH HF, not against sim_min_hf. sim_min_hf is
# only the floor the simulator is permitted to descend to; the frequency actually
# achieved is usually well above it (the standard chain config floors at 0.296
# under sim_min_hf=0.05), so checking against sim_min_hf rejects working configs.
#
# For a chain the minimum is exact and deterministic. For random/star it is
# seed-dependent and unknowable here, so sim_min_hf is the only guaranteed bound
# and the check is necessarily weaker.
GATES = (("call -f (frequency_threshold)", freq),
         ("denoise --vaf", denoise),
         ("call -v (min_hf)", min_hf))

def check_gates(min_truth_hf, exact, ctx):
    for label, value in GATES:
        if value >= min_truth_hf:
            errors.append(
                f"{label} = {value} is at or above the {'' if exact else 'lowest possible '}"
                f"truth HF {min_truth_hf:.4f}{ctx}. This gate deletes truth variants "
                "before they can be called, and depth cannot compensate."
            )
        elif value > min_truth_hf / 2.0:
            warnings.append(
                f"{label} = {value} has under 2x margin on the {'' if exact else 'lowest possible '}"
                f"truth HF {min_truth_hf:.4f}{ctx}. Measured: denoise --vaf 0.03 against a "
                "truth HF of 0.033 cut var_recall to 0.533, and 0.01 restored it to 0.867 "
                "on the same reads. Check var_recall in the output."
            )

if topology == "chain":
    # A chain never splits mass between siblings, so the deepest clone is exactly
    # (1 - ref_fraction) * (1 - internal_keep)^n -- deterministic, so an
    # infeasible chain cannot be rescued by reseeding.
    for n in sorted(set(nmuts)):
        deepest = (1.0 - REF_FRACTION) * (1.0 - keep) ** n
        if deepest <= sim_min:
            errors.append(
                f"topology=chain with n_mutations={n} puts the deepest clone at "
                f"{deepest:.4f}, at or below sim_min_hf={sim_min}. Chain "
                "frequencies are deterministic, so no seed fixes this: lower "
                "sim_min_hf or internal_keep, or reduce n_mutations."
            )
        else:
            check_gates(deepest, True, f" for n_mutations={n}")
else:
    check_gates(sim_min, False, "")

# --- homopolymer modes ---------------------------------------------------------
try:
    bg_tokens = simulate_tree.parse_background(background)
except ValueError as exc:
    errors.append(str(exc))
    bg_tokens = None
if bg_tokens:
    ref_path = "/opt/lineage_sim/rCRS.fasta"
    try:
        _, rcrs = hp_tracts.load_reference(ref_path)
        hp_tracts.apply_background(rcrs, bg_tokens)
    except Exception as exc:
        errors.append(f"background {background!r} does not match rCRS: {exc}")
if hp_unstable:
    try:
        simulate_tree.parse_hp_unstable(hp_unstable)
    except ValueError as exc:
        errors.append(f"hp_unstable: {exc}")
n_drawable = sum(1 for tract in hp_tracts.CATALOGUE if tract.simulate)
if hp_k < 0:
    errors.append(f"hp_lineage_indels={hp_k} is negative")
if hp_k > n_drawable:
    errors.append(
        f"hp_lineage_indels={hp_k} exceeds the {n_drawable} simulate=True "
        "tracts (draws are without replacement)"
    )
for n in sorted(set(nmuts)):
    if hp_k > n:
        errors.append(
            f"hp_lineage_indels={hp_k} exceeds n_mutations={n}"
        )
if require_hp_mode and not (hp_k > 0 or hp_unstable != ""):
    errors.append(
        "the rates have nothing to act on and the sweep is guaranteed flat"
    )
if hp_unstable and not collapse_hp:
    warnings.append(
        "hp_unstable is set and collapse_hp is false, so unstable alleles "
        "enter the tree under several spellings"
    )
hp_mode = background != "none" or hp_k > 0 or bool(hp_unstable) or collapse_hp
if (hp_fp_set or hp_fn_set) and not hp_mode:
    warnings.append(
        "hp_fp_rate/hp_fn_rate are set and no HP mode is on; they have nothing to act on"
    )

# --- non-fatal ----------------------------------------------------------------
if profile == "ont-r10":
    warnings.append(
        "profile=ont-r10 SKIPS Himito denoise, so indel artifacts crowd the graph "
        "and the VCF usually comes back empty. Use ont-denoised unless you are "
        "deliberately measuring the un-denoised path."
    )
if topology == "star":
    warnings.append(
        "topology=star produces no ancestral pairs, so ad_recall has an empty "
        "denominator and ad_* is uninformative. It is a control, not a config."
    )

for w in warnings:
    print(f"WARNING: {w}", file=sys.stderr)
if errors:
    print()
    for e in errors:
        print(f"ERROR: {e}")
    sys.exit(1)

print(f"config OK: profile={profile} topology={topology} internal_keep={keep}")
if topology == "chain":
    basis = ", ".join(
        f"n={n}:{(1.0 - REF_FRACTION) * (1.0 - keep) ** n:.4f}" for n in sorted(set(nmuts)))
    print(f"  gates checked against the exact chain minimum ({basis})")
else:
    print(f"  gates checked against sim_min_hf={sim_min} (lowest possible truth HF; "
          "the achieved minimum is seed-dependent and usually higher)")
print(f"  call -f={freq}  denoise --vaf={denoise}  call -v={min_hf}")
PY
        echo ok > ok.txt
    >>>

    output {
        String ok = read_string("ok.txt")
    }

    runtime {
        cpu:    1
        memory: "1 GiB"
        disks:  "local-disk 10 HDD"
        docker: docker
    }
}

task RunEvalCell {

    meta {
        description: "One simulate -> reconstruct -> score cycle (run_eval.sh) for a single (seed, n_mutations, total_depth) grid cell."
    }

    parameter_meta {
        seed:          "RNG seed for simulate_tree.py / simulate_reads.sh"
        n_mutations:   "number of heteroplasmic SNVs in the truth tree. NOT a difficulty knob: clone mass is conserved, so raising it only costs variant recall"
        total_depth:   "total simulated read depth across all clones"
        profile:       "read profile: ont-denoised (recommended), ont-r10 (skips denoise; usually empty VCF), or hifi"
        fp_arg:        "optional '--fp <rate>' override; empty uses the profile's Himito default"
        fn_arg:        "optional '--fn <rate>' override; empty uses the profile's Himito default"
        denoise_vaf_arg: "optional '--denoise-vaf <F>' override; empty follows --vaf. Must stay below sim_min_hf"
        frequency_threshold: "Himito call -f. Must stay below sim_min_hf; Himito's own default (0.2) empties the VCF"
        min_hf:        "shared HF floor for Himito call -v and lineage --min-hf"
        sim_min_hf:    "floor for the SIMULATED truth frequencies -- distinct from min_hf, which filters what was called"
        topology:      "truth tree shape: random | chain | star. chain raises ordering difficulty without costing detection"
        internal_keep: "terminal mass kept by each internal node; lower is harder to order and costs no detection"
        background:    "fixed homoplasmies: none | hvs1-polyC | comma-separated <REF><pos><ALT>"
        hp_lineage_indels: "K of the n_mutations edges replaced by HP length changes; 0 leaves the tree SNV-only"
        hp_unstable:   "lineage-independent tract-length spec; empty omits --hp-unstable"
        collapse_hp:   "forward Himito lineage --collapse-hp when true"
        hp_unstable_txt: "hp_unstable written into cell_metrics; NA when unset"
        hp_fp_txt:     "hp_fp_rate written into cell_metrics; NA when unset"
        hp_fn_txt:     "hp_fn_rate written into cell_metrics; NA when unset"
        background_arg: "optional '--background <spec>'; empty when background is none"
        hp_lineage_arg: "optional '--hp-lineage-indels K'; empty when K is 0"
        hp_unstable_arg: "optional '--hp-unstable <spec>'; empty when unset. Quoted because specs contain ';'"
        hp_min_subclone_arg: "optional '--hp-min-subclone-freq F'; empty unless hp_unstable is set"
        collapse_hp_arg: "optional '--collapse-hp'; empty when false"
        hp_fp_arg:     "optional '--hp-fp-rate <rate>'; empty uses Himito's 2x fp default"
        hp_fn_arg:     "optional '--hp-fn-rate <rate>'; empty uses Himito's 2x fn default"
        reference_fa:  "mitochondrial reference FASTA (required)"
        validated:     "sentinel from ValidateConfig; forces gate validation to run before any cell"
    }

    input {
        File reference_fa
        Int seed
        Int n_mutations
        Int total_depth

        String profile
        String fp_arg
        String fn_arg
        String denoise_vaf_arg
        Float frequency_threshold

        Float min_hf
        Float max_hf
        Float sim_min_hf
        Float sim_max_hf
        String topology
        Float internal_keep

        String background
        Int hp_lineage_indels
        String hp_unstable
        Boolean collapse_hp
        String hp_unstable_txt
        String hp_fp_txt
        String hp_fn_txt
        String background_arg
        String hp_lineage_arg
        String hp_unstable_arg
        String hp_min_subclone_arg
        String collapse_hp_arg
        String hp_fp_arg
        String hp_fn_arg

        String validated

        String docker
        RuntimeAttr? runtime_attr_override
    }

    String cell = "seed~{seed}_mut~{n_mutations}_depth~{total_depth}" + (if hp_lineage_indels > 0 || hp_unstable != "" then "_hp" else "")
    command <<<
        set -euxo pipefail

        # Tools + scripts + models + Himito binary are baked into the image;
        # export the paths the shell scripts look for so they run standalone.
        export HIMITO="/Himito/target/release/Himito"
        export PBSIM_MODEL_DIR="/opt/lineage_sim/pbsim3_models"
        export REF="~{reference_fa}"

        /opt/lineage_sim/run_eval.sh \
            --outdir "~{cell}" \
            --profile "~{profile}" \
            --n-mutations ~{n_mutations} \
            --total-depth ~{total_depth} \
            --seed ~{seed} \
            --ref "$REF" \
            --min-hf ~{min_hf} \
            --max-hf ~{max_hf} \
            --sim-min-hf ~{sim_min_hf} \
            --sim-max-hf ~{sim_max_hf} \
            --topology ~{topology} \
            --internal-keep ~{internal_keep} \
            --frequency-threshold ~{frequency_threshold} \
            ~{fp_arg} \
            ~{fn_arg} \
            ~{denoise_vaf_arg} \
            ~{background_arg} \
            ~{hp_lineage_arg} \
            ~{hp_unstable_arg} \
            ~{hp_min_subclone_arg} \
            ~{collapse_hp_arg} \
            ~{hp_fp_arg} \
            ~{hp_fn_arg}

        # Prepend the identifying columns run_simulation.sh adds to the combined
        # table. Read the header from the file rather than assuming
        # score_lineage.py's column order.
        metrics="~{cell}/metrics.tsv"
        header=$(head -n 1 "$metrics")
        row=$(tail -n 1 "$metrics")
        # NA for unset hp_unstable / hp_fp_rate / hp_fn_rate.
        printf 'seed\tn_mutations\ttotal_depth\ttopology\tinternal_keep\tbackground\thp_lineage_indels\thp_unstable\tcollapse_hp\thp_fp_rate\thp_fn_rate\t%s\n' "$header" > cell_metrics.tsv
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "~{seed}" "~{n_mutations}" "~{total_depth}" "~{topology}" "~{internal_keep}" \
            "~{background}" "~{hp_lineage_indels}" "~{hp_unstable_txt}" "~{collapse_hp}" \
            "~{hp_fp_txt}" "~{hp_fn_txt}" "$row" >> cell_metrics.tsv

        # Bundle the full cell (truth/, reads/, himito/) for debugging/inspection.
        tar czf "~{cell}.tar.gz" "~{cell}"
        # truth/ alone, including hp_tracts.tsv, hp_lineage.tsv, background.txt.
        tar czf "~{cell}_truth.tar.gz" -C "~{cell}" truth
    >>>

    output {
        File cell_metrics = "cell_metrics.tsv"
        File cell_bundle = "~{cell}.tar.gz"
        File truth_dir_tar = "~{cell}_truth.tar.gz"
        File matrix = "~{cell}/himito/sim.matrix.csv"
        File vcf = "~{cell}/himito/sim.vcf"
        File cleaned_matrix = "~{cell}/himito/sim_lineage.cleaned_matrix_all_variants.csv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            50,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             docker
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task CombineMetrics {

    meta {
        description: "Concatenate per-cell metrics.tsv files into one table, keeping a single header."
    }

    input {
        Array[File] cell_metrics
        String docker
    }

    command <<<
        set -euxo pipefail
        first=1
        expected=""
        : > metrics.tsv
        for f in ~{sep=' ' cell_metrics}; do
            [[ -s "$f" ]] || continue
            header=$(head -n 1 "$f")
            if [[ $first -eq 1 ]]; then
                expected="$header"
                cat "$f" >> metrics.tsv
                first=0
            else
                # Position-wise concat misaligns mixed schemas.
                if [[ "$header" != "$expected" ]]; then
                    echo "ERROR: header mismatch: $f" >&2
                    exit 1
                fi
                tail -n +2 "$f" >> metrics.tsv
            fi
        done
        column -t metrics.tsv || true
    >>>

    output {
        File combined = "metrics.tsv"
    }

    runtime {
        cpu:    1
        memory: "2 GiB"
        disks:  "local-disk 10 HDD"
        docker: docker
    }
}
