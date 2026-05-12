#!/bin/bash
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00
# =============================================================================
# Assemblers.sh  (optimized; LOSO-aware)
#
#
# Compatible with:
#   - Original 4-column manifests:  sample_id  fwd_fq  rev_fq  ref_fasta
#   - LOSO manifests from Simulate_rnaseq_loso.sh
#     (e.g. "test_M_superfamily_200x_PE_samples.txt", a single data row whose
#     4th column is the held-out test_<sf>.fasta used as TransRate ground truth).
#
# Key changes vs. the original:
#   * The trailing "_samples" suffix in LOSO manifest names is stripped when
#     deriving paths -> clean fold stems like  test_M_superfamily_200x_PE
#   * Single-row manifests symlink the reads instead of `cat`-ing them
#     (large IO win for LOSO folds, where there's exactly one sample per file).
#   * SLURM_CPUS_ON_NODE / SLURM_MEM_PER_NODE are respected; -t/-m override.
#   * `set -eo pipefail`, robust quoting, manifest comment/blank-line handling,
#     per-fold transrate log instead of >/dev/null, and missing-assembly cases
#     are diverted to issues_dir/ without aborting the whole batch.
#   * Config-template expansion (`eval echo $run_tool ...`) and the variable
#     names it references (forward_fq, reverse_fq, OUTDIR, CPU, MEM,
#     FASTA_DIR, REFERENCE) are preserved verbatim for backward compat.
# =============================================================================

set -eo pipefail

# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------
print_help() {
    cat <<EOF
Usage: $(basename "$0") -s <Manifest> -c <Config> [-t <CPUs>] [-m <MEM_GB>]

  -s   Manifest TSV (sample_id, fwd_fq, rev_fq, ref_fasta).
       Use 'all' to process every *.txt manifest in the current directory.
  -c   Assembler config file (one 'tool=command_template' per line).
  -t   CPU threads.  Default: \$SLURM_CPUS_ON_NODE, else 20.
  -m   Memory in GB. Default: \$SLURM_MEM_PER_NODE/1024, else 100.
  -h   This message.

The config templates are eval-expanded with the in-scope variables:
  \$forward_fq  \$reverse_fq  \$OUTDIR  \$CPU  \$MEM  \$FASTA_DIR  \$REFERENCE
EOF
}

MANIFEST_FILE=""
CONFIG_FILE=""
CPU_OPT=""
MEM_OPT=""

while getopts "s:c:t:m:h" opt; do
    case $opt in
        s) MANIFEST_FILE="$OPTARG";;
        c) CONFIG_FILE="$OPTARG";;
        t) CPU_OPT="$OPTARG";;
        m) MEM_OPT="$OPTARG";;
        h) print_help; exit 0;;
        ?) print_help >&2; exit 1;;
    esac
done
shift $((OPTIND-1))

if [[ -z "$MANIFEST_FILE" || -z "$CONFIG_FILE" ]]; then
    echo "Error: -s and -c are required." >&2
    print_help >&2
    exit 1
fi
if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "Error: config file not found: $CONFIG_FILE" >&2
    exit 1
fi

# -----------------------------------------------------------------------------
# Resource discovery (SLURM-aware, with CLI override)
# -----------------------------------------------------------------------------
CPU="${CPU_OPT:-${SLURM_CPUS_ON_NODE:-20}}"

if [[ -n "$MEM_OPT" ]]; then
    MEM="$MEM_OPT"
elif [[ -n "${SLURM_MEM_PER_NODE:-}" ]]; then
    # SLURM_MEM_PER_NODE is in MB
    MEM=$(( SLURM_MEM_PER_NODE / 1024 ))
    (( MEM < 1 )) && MEM=1
else
    MEM=100
fi
export CPU MEM   # so eval'd templates can read them

# -----------------------------------------------------------------------------
# Make assembler wrappers visible on PATH (if the dir exists)
# -----------------------------------------------------------------------------
if [[ -d "$PWD/run_assemblers_dir" ]]; then
    chmod +x "$PWD"/run_assemblers_dir/*.sh 2>/dev/null || true
    export PATH="$PWD/run_assemblers_dir:$PATH"
else
    echo "Note: $PWD/run_assemblers_dir not found; assuming wrappers are on PATH." >&2
fi

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

# Derive a clean fold stem from a manifest path.
# Strips the trailing "_samples" produced by Simulate_rnaseq_loso.sh, so that
# downstream artifacts use e.g. test_M_superfamily_200x_PE rather than
# test_M_superfamily_200x_PE_samples.
manifest_stem() {
    local bs
    bs=$(basename "${1%.*}")
    printf '%s\n' "${bs%_samples}"
}

# Stream manifest rows with comments + blank lines removed.
read_manifest_rows() {
    awk 'BEGIN{FS=OFS="\t"}
         /^[[:space:]]*$/  {next}
         /^[[:space:]]*#/  {next}
         {print}' "$1"
}

run_reference_transrate() {
    # Args:  $1 assembly FASTA, $2 reference FASTA
    export PATH=/LUSTRE/apps/bioinformatica/.local/bin:$PATH
    export PATH=/LUSTRE/apps/bioinformatica/ruby2/ruby-2.2.0/bin:$PATH
    export LD_LIBRARY_PATH=/LUSTRE/apps/bioinformatica/ruby2/ruby-2.2.0/lib:$LD_LIBRARY_PATH

    local query="$1"
    local target="$2"

    mkdir -p transrate_tmp_dir transrate_contigs_dir

    if [[ ! -s "$query" ]]; then
        echo "  [transrate] skipping: empty/missing assembly $query" >&2
        return 1
    fi
    if [[ ! -s "$target" ]]; then
        echo "  [transrate] skipping: empty/missing reference $target" >&2
        return 1
    fi

    local bs
    bs=$(basename "${query%.*}")
    local TRANSRATE_DIR="${bs}_dir"
    local tr_path="transrate_tmp_dir/$TRANSRATE_DIR"
    local tr_log="transrate_tmp_dir/${bs}.log"

    echo "  [transrate] $bs  vs  $(basename "$target")"
    local call="transrate --assembly $query --reference $target --output $tr_path --threads $CPU"
    echo "    $call"

    if ! eval "$call" > "$tr_log" 2>&1; then
        echo "  [transrate] failed for $bs (log: $tr_log)" >&2
        return 1
    fi

    # Move every contigs.csv (preserving sub-directory structure under
    # transrate_contigs_dir) so multiple folds don't collide.
    local n_moved=0
    while IFS= read -r contig_file; do
        local dst="transrate_contigs_dir/${contig_file#transrate_tmp_dir/}"
        mkdir -p "$(dirname "$dst")"
        mv "$contig_file" "$dst"
        n_moved=$(( n_moved + 1 ))
    done < <(find "$tr_path" -name "contigs.csv")

    if (( n_moved == 0 )); then
        echo "  [transrate] warning: no contigs.csv produced for $bs" >&2
    fi

    rm -rf "$tr_path"
}

run_assembly_batches() {
    local manifest="$1"
    local config="$2"

    if [[ ! -f "$manifest" ]]; then
        echo "Error: manifest not found: $manifest" >&2
        return 1
    fi

    local stem
    stem=$(manifest_stem "$manifest")

    local FASTA_DIR="${stem}_FASTA_DIR"
    mkdir -p "$FASTA_DIR" chkp_dir

    # ---------------------------------------------------------------
    # Slurp manifest rows (sample_id, fwd, rev, ref)
    # ---------------------------------------------------------------
    local -a sids=() fwds=() revs=() refs=()
    while IFS=$'\t' read -r sid fwd rev ref _rest; do
        [[ -z "${sid:-}" ]] && continue
        sids+=("$sid"); fwds+=("$fwd"); revs+=("$rev"); refs+=("$ref")
    done < <(read_manifest_rows "$manifest")

    local n=${#sids[@]}
    if (( n == 0 )); then
        echo "Warning: $manifest has no data rows; skipping." >&2
        return 0
    fi
    echo "============================================================"
    echo "Manifest: $manifest"
    echo "  Stem:        $stem"
    echo "  Samples:     $n"
    echo "  Resources:   CPU=$CPU  MEM=${MEM}GB"

    # All rows within a LOSO fold share the same reference (test_<sf>.fasta).
    # For legacy multi-row manifests we also use row 1's ref, warning if rows
    # disagree.
    local REFERENCE="${refs[0]}"
    for r in "${refs[@]}"; do
        if [[ "$r" != "$REFERENCE" ]]; then
            echo "  Warning: manifest has mixed reference columns; using $REFERENCE" >&2
            break
        fi
    done
    echo "  Reference:   $REFERENCE"

    # ---------------------------------------------------------------
    # Prepare reads: symlink for n==1 (LOSO), concatenate for n>1
    # ---------------------------------------------------------------
    local forward_fq="${stem}_concat_PE1.fq"
    local reverse_fq="${stem}_concat_PE2.fq"
    # Names match the original so config templates that reference $forward_fq
    # / $reverse_fq keep working.

    rm -f "$forward_fq" "$reverse_fq"
    if (( n == 1 )); then
        ln -s "${fwds[0]}" "$forward_fq"
        ln -s "${revs[0]}" "$reverse_fq"
        echo "  Reads:       symlinked (single-sample manifest)"
    else
        : > "$forward_fq"
        : > "$reverse_fq"
        local f
        for f in "${fwds[@]}"; do [[ -s "$f" ]] && cat "$f" >> "$forward_fq"; done
        for f in "${revs[@]}"; do [[ -s "$f" ]] && cat "$f" >> "$reverse_fq"; done
        echo "  Reads:       concatenated ($n samples)"
    fi
    export forward_fq reverse_fq FASTA_DIR REFERENCE

    # ---------------------------------------------------------------
    # Run each tool from the config
    # ---------------------------------------------------------------
    local line tool run_tool OUTDIR chkp_file final_fasta call
    while IFS= read -r line || [[ -n "$line" ]]; do
        [[ "$line" =~ ^[[:space:]]*# ]] && continue
        [[ -z "${line// }" ]] && continue

        tool="${line%%=*}"
        run_tool="${line#*=}"
        # trim leading/trailing whitespace
        tool="$(echo -n "$tool" | awk '{$1=$1;print}')"
        run_tool="$(echo -n "$run_tool" | awk '{$1=$1;print}')"
        OUTDIR="${stem}_${tool}_dir"
        chkp_file="1_${tool}_${stem}.chkp"

        if [[ -f "chkp_dir/$chkp_file" ]]; then
            echo "[$tool @ $stem] checkpoint exists; skipping."
            continue
        fi

        echo "[$tool @ $stem] Step 1: assembling..."
        mkdir -p "$OUTDIR"
        export OUTDIR

        # Preserve the original eval-echo template-expansion behavior.
        call=$(eval echo "$run_tool" \
                    "$forward_fq" "$reverse_fq" "$OUTDIR" \
                    "$CPU" "$MEM" "$FASTA_DIR" "$REFERENCE")
        echo "  $call"

        if ! eval "$call"; then
            echo "  [$tool @ $stem] command exited non-zero; saving to issues_dir/" >&2
            mkdir -p issues_dir
            mv "$OUTDIR" "issues_dir/" 2>/dev/null || true
            continue
        fi

        final_fasta="$FASTA_DIR/${OUTDIR%_dir}.fa"
        if [[ ! -s "$final_fasta" ]]; then
            echo "  [$tool @ $stem] no/empty assembly at $final_fasta; saving to issues_dir/" >&2
            mkdir -p issues_dir
            mv "$OUTDIR" "issues_dir/" 2>/dev/null || true
            continue
        fi

        touch "chkp_dir/$chkp_file"
        rm -rf "$OUTDIR"
        echo "  [$tool @ $stem] assembly -> $final_fasta"

        echo "[$tool @ $stem] Step 2: TransRate (reference = ground truth)"
        run_reference_transrate "$final_fasta" "$REFERENCE" || \
            echo "  [$tool @ $stem] TransRate scoring did not complete (see log)." >&2

    done < "$config"

    # Clean up the staged read pair (symlinks or temp files)
    rm -f "$forward_fq" "$reverse_fq"
}

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
if [[ "$MANIFEST_FILE" == "all" ]]; then
    echo "Batch mode: every *.txt manifest under $PWD"
    shopt -s nullglob
    manifests=( "$PWD"/*.txt )
    if (( ${#manifests[@]} == 0 )); then
        echo "Error: no *.txt manifests in $PWD" >&2
        exit 1
    fi
    for m in "${manifests[@]}"; do
        run_assembly_batches "$m" "$CONFIG_FILE"
    done
else
    run_assembly_batches "$MANIFEST_FILE" "$CONFIG_FILE"
fi

echo "============================================================"
echo "Finished assembly batch."
exit 0