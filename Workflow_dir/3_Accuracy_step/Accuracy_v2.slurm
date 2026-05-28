#!/bin/bash
#SBATCH --job-name=transrate_coverage
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00
#SBATCH -o out.%J
#SBATCH -e err.%J
# =============================================================================
# Accuracy_v2.slurm
#
# Consumes the assemblies produced by Assemblers_v2.sh:
#     <stem>_FASTA_DIR/<stem>_<tool>.fa
# Pairs each with its LOSO manifest:
#     <stem>_samples.txt    (preferred)
#     <stem>.txt            (fallback for legacy non-LOSO manifests)
# and runs Accuracy_v2.sh in read-based mode by default — because
# Assemblers_v2.sh already covers reference-based TransRate against the
# held-out ground truth.
#
# Override mode by exporting before sbatch, e.g.:
#     MODE=both   sbatch Accuracy_v2.slurm
#     MODE=reads  sbatch Accuracy_v2.slurm    # default
# =============================================================================

set -eo pipefail

MODE="${MODE:-reads}"                 # reads | reference | both
SCRIPT="${SCRIPT:-./Accuracy_v2.sh}"  # path to Accuracy_v2.sh

[[ -x "$SCRIPT" ]] || { echo "Error: $SCRIPT not found / not executable." >&2; exit 1; }

mkdir -p chkp_dir

# -----------------------------------------------------------------------------
# Discover assemblies.  Prefer the Assemblers_v2.sh layout (per-fold FASTA_DIR);
# fall back to a flat layout for backward compat with the original pipeline.
# -----------------------------------------------------------------------------
shopt -s nullglob
assemblies=( *_FASTA_DIR/*.fa )
layout="loso"

if (( ${#assemblies[@]} == 0 )); then
    assemblies=( *.fa )
    layout="flat"
fi

if (( ${#assemblies[@]} == 0 )); then
    echo "Error: no assemblies found (looked in *_FASTA_DIR/*.fa, then *.fa)." >&2
    exit 1
fi

echo "Layout=$layout   mode=$MODE   assemblies=${#assemblies[@]}"

# -----------------------------------------------------------------------------
# Resolve assembly -> manifest:
#   loso layout:  parent dir is <stem>_FASTA_DIR  ->  stem is unambiguous
#   flat layout:  strip trailing _<tool> off the file basename (best effort;
#                 only used when the user hasn't migrated to Assemblers_v2.sh)
# -----------------------------------------------------------------------------
locate_manifest() {
    local assembly="$1"
    local d stem c
    d=$(dirname "$assembly")

    if [[ "$(basename "$d")" == *_FASTA_DIR ]]; then
        stem=$(basename "$d"); stem="${stem%_FASTA_DIR}"
    else
        # Legacy: filename is <stem>_<tool>.fa; strip last "_token".
        local fn
        fn=$(basename "${assembly%.*}")
        stem="${fn%_*}"
    fi

    for c in "${stem}_samples.txt" \
             "${stem}.txt" \
             "${d}/../${stem}_samples.txt" \
             "${d}/../${stem}.txt"; do
        [[ -f "$c" ]] && { printf '%s\n' "$c"; return 0; }
    done
    return 1
}

# -----------------------------------------------------------------------------
# Process every assembly with checkpointing.
# -----------------------------------------------------------------------------
n_done=0; n_skip=0; n_fail=0; n_nomf=0

for a in "${assemblies[@]}"; do
    if [[ ! -s "$a" ]]; then
        echo "  [skip-empty] $a"
        continue
    fi

    if ! manifest=$(locate_manifest "$a"); then
        echo "  [WARN] no manifest found for $a — skipping" >&2
        n_nomf=$(( n_nomf + 1 ))
        continue
    fi

    chkp_file="$(basename "${a%.fa}")_${MODE}.chkp"
    if [[ -f "chkp_dir/$chkp_file" ]]; then
        echo "  [skip] $a — checkpoint $chkp_file already exists"
        n_skip=$(( n_skip + 1 ))
        continue
    fi

    echo "------------------------------------------------------------"
    echo "Scoring  assembly = $a"
    echo "         manifest = $manifest"
    echo "         mode     = $MODE"

    if "$SCRIPT" -s "$manifest" -d "$a" -m "$MODE"; then
        touch "chkp_dir/$chkp_file"
        n_done=$(( n_done + 1 ))
    else
        echo "  [FAIL] $a (see Accuracy_v2.sh output above)" >&2
        n_fail=$(( n_fail + 1 ))
    fi
done

echo "============================================================"
echo "Accuracy batch complete:"
echo "    done       = $n_done"
echo "    skipped    = $n_skip"
echo "    failed     = $n_fail"
echo "    no-manifest= $n_nomf"
exit 0
