#!/bin/bash
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00
# =============================================================================
# Assemblers_v3.sh  (optimized; LOSO-aware; mode-aware TransRate)
#
# Builds on Assemblers_v2.sh and folds in the read-handling machinery from
# Accuracy_v2.sh so that Step 2 (TransRate) can be run in any of four modes:
#
#     -m none        assemble only; SKIP TransRate entirely
#     -m reference   reference-only metrics  (Assemblers_v2.sh behaviour)
#     -m reads       read-only metrics       (left/right vs assembly)
#     -m both        reference + reads       (legacy "longer" mode)
#
# For -m reads / -m both the staged read pair is sanitized before scoring
# (seqkit sana | seqkit pair) exactly as Accuracy_v2.sh does, preventing
# TransRate's "Unmatched read IDs ... Use the -I option to ignore this"
# failure.  Use -S to skip sanitization.
#
# Compatible with:
#   - Original 4-column manifests:  sample_id  fwd_fq  rev_fq  ref_fasta
#   - LOSO manifests from Simulate_rnaseq_loso.sh
#     (e.g. "test_M_superfamily_200x_PE_samples.txt", a single data row whose
#     4th column is the held-out test_<sf>.fasta used as TransRate ground truth).
#
# Key properties carried over from v2:
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
#
# New in v3:
#   * -m { none | reference | reads | both }   TransRate mode selector.
#   * -S                                       skip seqkit read sanitization.
#   * stage_reads() / sanitize_reads()          ported from Accuracy_v2.sh.
#   * The assembly read pair (forward_fq/reverse_fq) is reused for read-mode
#     TransRate, so reads are never staged twice.
#
# NOTE: the meaning of -m changed between v2 and v3.
#   v2:  -m  ==  memory in GB
#   v3:  -m  ==  TransRate mode  (memory is now -c, mirroring Accuracy_v2.sh)
# =============================================================================

set -eo pipefail

# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------
print_help() {
    cat <<EOF
Usage: $(basename "$0") -s <Manifest> -c <Config>
                       [-m none|reference|reads|both] [-t <CPUs>] [-C <MEM_GB>] [-S]

  -s   Manifest TSV (sample_id, fwd_fq, rev_fq, ref_fasta).
       Use 'all' to process every *.txt manifest in the current directory.
  -c   Assembler config file (one 'tool=command_template' per line).
  -m   TransRate mode for Step 2.  Default: reference
         none       assemble only; do NOT run TransRate
         reference  reference-only metrics (assembly vs held-out ground truth)
         reads      read-only metrics      (left/right vs assembly)
         both       reference + reads
  -t   CPU threads.  Default: \$SLURM_CPUS_ON_NODE, else 20.
  -C   Memory in GB. Default: \$SLURM_MEM_PER_NODE/1024, else 100.
       (Informational; TransRate has no memory cap, but MEM is exported so
        config templates / user hooks can read it.)
  -S   Skip the seqkit sana + seqkit pair sanitization step for reads/both
       modes.  Sanitization (default ON) prevents TransRate's
       "Unmatched read IDs ... Use the -I option to ignore this" failure.
       Requires seqkit on PATH (https://bioinf.shenwei.me/seqkit/).
  -h   This message.

The config templates are eval-expanded with the in-scope variables:
  \$forward_fq  \$reverse_fq  \$OUTDIR  \$CPU  \$MEM  \$FASTA_DIR  \$REFERENCE
EOF
}

MANIFEST_FILE=""
CONFIG_FILE=""
MODE="reference"
CPU_OPT=""
MEM_OPT=""
SANITIZE=1

while getopts "s:c:m:t:C:Sh" opt; do
    case $opt in
        s) MANIFEST_FILE="$OPTARG";;
        c) CONFIG_FILE="$OPTARG";;
        m) MODE="$OPTARG";;
        t) CPU_OPT="$OPTARG";;
        C) MEM_OPT="$OPTARG";;
        S) SANITIZE=0;;
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

case "$MODE" in
    none|reference|reads|both) ;;
    *) echo "Error: invalid -m '$MODE' (use none|reference|reads|both)" >&2
       print_help >&2
       exit 1;;
esac

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

setup_transrate_env() {
    export PATH=/LUSTRE/apps/bioinformatica/.local/bin:$PATH
    export PATH=/LUSTRE/apps/bioinformatica/ruby2/ruby-2.2.0/bin:$PATH
    export LD_LIBRARY_PATH=/LUSTRE/apps/bioinformatica/ruby2/ruby-2.2.0/lib:${LD_LIBRARY_PATH:-}
}

# ---------------------------------------------------------------------------
# stage_reads  (ported from Accuracy_v2.sh)
#
# Creates ${stem}_concat_PE1.fq and ${stem}_concat_PE2.fq.
# Symlinks for n==1 (LOSO folds), concatenates otherwise.  Single-row
# manifests never copy raw reads.
#
# Expects globals: n, fwds[], revs[]  (populated by run_assembly_batches).
# ---------------------------------------------------------------------------
stage_reads() {
    local stem="$1"
    local fwd="${stem}_concat_PE1.fq"
    local rev="${stem}_concat_PE2.fq"
    rm -f "$fwd" "$rev"
    if (( n == 1 )); then
        ln -s "${fwds[0]}" "$fwd"
        ln -s "${revs[0]}" "$rev"
        echo "  Reads:       symlinked (single-sample manifest)"
    else
        : > "$fwd"; : > "$rev"
        local f
        for f in "${fwds[@]}"; do [[ -s "$f" ]] && cat "$f" >> "$fwd"; done
        for f in "${revs[@]}"; do [[ -s "$f" ]] && cat "$f" >> "$rev"; done
        echo "  Reads:       concatenated ($n samples)"
    fi
}

# ---------------------------------------------------------------------------
# sanitize_reads  (ported from Accuracy_v2.sh)
#
# Repair broken FASTQ records and enforce matched R1/R2 read IDs.
# Prevents TransRate's "Unmatched read IDs ... Use the -I option to ignore
# this" failure (snap.cpp).  Strategy:
#   1) seqkit sana   : drop or repair broken single-line FASTQ records
#   2) seqkit pair   : keep only reads whose IDs appear in BOTH R1 and R2
# Operates on the staged ${stem}_concat_PE{1,2}.fq files in-place.  If the
# inputs were symlinks (LOSO single-sample case) the original source files
# under fwds[0]/revs[0] are NEVER mutated — we delete the symlink and replace
# it with a fresh sanitized file.
# ---------------------------------------------------------------------------
sanitize_reads() {

    EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/
    export PATH=$PATH:$EXPORT

    local stem="$1"
    local fwd="${stem}_concat_PE1.fq"
    local rev="${stem}_concat_PE2.fq"

    if ! command -v seqkit >/dev/null 2>&1; then
        echo "  WARNING: seqkit not on PATH; skipping pair sanitization." >&2
        echo "           If TransRate fails with 'Unmatched read IDs', install seqkit" >&2
        echo "           from https://bioinf.shenwei.me/seqkit/ and re-run." >&2
        return 0
    fi

    local sana_fwd="${stem}_sana_PE1.fq"
    local sana_rev="${stem}_sana_PE2.fq"
    local pair_dir="${stem}_paired_tmp"
    local sana_log="transrate_tmp_dir/${stem}_seqkit.log"
    mkdir -p transrate_tmp_dir
    rm -f "$sana_fwd" "$sana_rev"
    rm -rf "$pair_dir"

    echo "  Sanitizing reads (seqkit sana | seqkit pair) ..." >&2

    # ---- Step 1: seqkit sana — repair broken FASTQ records --------------
    if ! seqkit sana --quiet -j "$CPU" "$fwd" -o "$sana_fwd" >"$sana_log"  2>&1 || \
       ! seqkit sana --quiet -j "$CPU" "$rev" -o "$sana_rev" >>"$sana_log" 2>&1
    then
        echo "  WARNING: seqkit sana failed (log: $sana_log); using originals." >&2
        rm -f "$sana_fwd" "$sana_rev"
        return 0
    fi

    # ---- Step 2: seqkit pair — drop orphan reads ------------------------
    if ! seqkit pair --quiet -j "$CPU" \
            -1 "$sana_fwd" -2 "$sana_rev" -O "$pair_dir" >>"$sana_log" 2>&1
    then
        echo "  WARNING: seqkit pair failed (log: $sana_log); using originals." >&2
        rm -f "$sana_fwd" "$sana_rev"
        rm -rf "$pair_dir"
        return 0
    fi

    local paired_fwd="${pair_dir}/$(basename "${sana_fwd%.fq}").paired.fq"
    local paired_rev="${pair_dir}/$(basename "${sana_rev%.fq}").paired.fq"

    if [[ ! -s "$paired_fwd" || ! -s "$paired_rev" ]]; then
        echo "  WARNING: seqkit pair produced empty output; using originals." >&2
        rm -f "$sana_fwd" "$sana_rev"
        rm -rf "$pair_dir"
        return 0
    fi

    # Replace staged files.  `rm -f $fwd` removes the symlink (not its target)
    # when n==1, so original LOSO reads remain untouched on disk.
    rm -f "$fwd" "$rev" "$sana_fwd" "$sana_rev"
    mv "$paired_fwd" "$fwd"
    mv "$paired_rev" "$rev"
    rm -rf "$pair_dir"

    local n_pairs
    n_pairs=$(( $(wc -l < "$fwd") / 4 ))
    echo "  Sanitized pairs retained: $n_pairs   (log: $sana_log)" >&2
}

# ---------------------------------------------------------------------------
# run_transrate  (mode-aware; replaces v2's run_reference_transrate)
#
# Args:
#   $1  query assembly FASTA
#   $2  reference FASTA  (used for reference/both)
#   $3  mode             (reference | reads | both)
#   $4  forward reads    (used for reads/both)
#   $5  reverse reads    (used for reads/both)
#
# The caller is responsible for staging/sanitizing the read pair; this
# function only consumes it.  Mode 'none' is filtered out before we get here.
# ---------------------------------------------------------------------------
run_transrate() {
    local query="$1"
    local target="$2"
    local mode="$3"
    local fwd="$4"
    local rev="$5"

    setup_transrate_env
    mkdir -p transrate_tmp_dir transrate_contigs_dir

    if [[ ! -s "$query" ]]; then
        echo "  [transrate] skipping: empty/missing assembly $query" >&2
        return 1
    fi
    if [[ "$mode" == "reference" || "$mode" == "both" ]] && [[ ! -s "$target" ]]; then
        echo "  [transrate] skipping: empty/missing reference $target" >&2
        return 1
    fi
    if [[ "$mode" == "reads" || "$mode" == "both" ]] && [[ ! -s "$fwd" || ! -s "$rev" ]]; then
        echo "  [transrate] skipping: empty/missing read pair ($fwd / $rev)" >&2
        return 1
    fi

    local bs
    bs=$(basename "${query%.*}")
    # Per-query, per-mode scratch dir keeps assemblies + folds + modes from
    # clobbering each other.
    local TRANSRATE_DIR="${bs}_${mode}_dir"
    local tr_path="transrate_tmp_dir/$TRANSRATE_DIR"
    local tr_log="transrate_tmp_dir/${bs}_${mode}.log"

    local call
    case "$mode" in
        reference) call="transrate --assembly $query --reference $target --output $tr_path --threads $CPU";;
        reads)     call="transrate --left $fwd --right $rev --assembly $query --output $tr_path --threads $CPU";;
        both)      call="transrate --left $fwd --right $rev --assembly $query --reference $target --output $tr_path --threads $CPU";;
    esac

    echo "  [transrate $mode] $bs"
    [[ "$mode" != "reads" ]] && echo "    REFERENCE:  $(basename "$target")"
    [[ "$mode" != "reference" ]] && echo "    LEFT/RIGHT: $fwd / $rev"
    echo "    $call"
    echo "    log: $tr_log"

    local rc=0
    if ! eval "$call" > "$tr_log" 2>&1; then
        echo "  [transrate $mode] failed for $bs (log: $tr_log)" >&2
        rc=1
    else
        # Move every contigs.csv (preserving sub-directory structure under
        # transrate_contigs_dir) so multiple folds/modes don't collide.
        local n_moved=0 contig_file dst
        while IFS= read -r contig_file; do
            dst="transrate_contigs_dir/${contig_file#transrate_tmp_dir/}"
            mkdir -p "$(dirname "$dst")"
            mv "$contig_file" "$dst"
            n_moved=$(( n_moved + 1 ))
        done < <(find "$tr_path" -name "contigs.csv")

        if (( n_moved == 0 )); then
            echo "  [transrate $mode] warning: no contigs.csv produced for $bs" >&2
        fi
    fi

    rm -rf "$tr_path"
    return $rc
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
    echo "  Mode:        $MODE"
    echo "  Resources:   CPU=$CPU  MEM=${MEM}GB  sanitize=$(( SANITIZE ))"

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
    # Names match the original so config templates that reference
    # $forward_fq / $reverse_fq keep working.
    # ---------------------------------------------------------------
    local forward_fq="${stem}_concat_PE1.fq"
    local reverse_fq="${stem}_concat_PE2.fq"

    stage_reads "$stem"

    # For reads/both modes the SAME staged pair is reused by TransRate, so we
    # sanitize it once, up front.  Assemblers are happy with sanitized reads
    # too (matched, repaired FASTQ records).
    if [[ "$MODE" == "reads" || "$MODE" == "both" ]] && (( SANITIZE )); then
        sanitize_reads "$stem"
    elif [[ "$MODE" == "reads" || "$MODE" == "both" ]]; then
        echo "  Skipping read sanitization (-S given)." >&2
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

        # -----------------------------------------------------------
        # Step 2: TransRate — only when MODE != none
        # -----------------------------------------------------------
        if [[ "$MODE" == "none" ]]; then
            echo "[$tool @ $stem] Step 2: skipped (-m none)."
        else
            echo "[$tool @ $stem] Step 2: TransRate (mode=$MODE)"
            run_transrate "$final_fasta" "$REFERENCE" "$MODE" \
                          "$forward_fq" "$reverse_fq" || \
                echo "  [$tool @ $stem] TransRate scoring did not complete (see log)." >&2
        fi

    done < "$config"

    # Clean up the staged read pair (symlinks or temp/sanitized files)
    rm -f "$forward_fq" "$reverse_fq"
}

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
echo "============================================================"
echo "Assemblers_v3.sh   mode=$MODE   CPU=$CPU   MEM=${MEM}GB   sanitize=$(( SANITIZE ))"

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
