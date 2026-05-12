#!/bin/bash
# =============================================================================
# Accuracy_v2.sh  — read-based TransRate scoring for Assemblers_v2.sh output
#
# Pipes directly off Assemblers_v2.sh:
#   * Same manifest format (sample_id, fwd_fq, rev_fq, ref_fasta)
#   * Same stem rule (trailing "_samples" on manifest basenames is stripped, so
#     test_M_superfamily_200x_PE_samples.txt -> stem test_M_superfamily_200x_PE)
#   * Same single-row LOSO optimization: symlink reads rather than cat
#
# Assemblers_v2.sh already produces *reference-based* TransRate scores against
# the held-out LOSO ground truth.  This script therefore defaults to *reads*
# mode (the part Assemblers_v2 does NOT do).  -m reference / -m both are kept
# for back-compat with the legacy Accuracy.sh.
# =============================================================================

set -eo pipefail

print_help() {
    cat <<EOF
Usage: $(basename "$0") -s <Manifest> -d <Assembly FASTA> [-m reference|reads|both]
                       [-t CPU] [-c MEM_GB] [-S]

  -s   Manifest TSV (sample_id, fwd_fq, rev_fq, ref_fasta).  Same format
       Assemblers_v2.sh consumes; trailing "_samples" suffix is auto-stripped.
  -d   Query assembly FASTA, e.g.  <stem>_FASTA_DIR/<stem>_<tool>.fa
  -m   reference | reads | both.  Default: reads
         reference  reference-only metrics (already done by Assemblers_v2.sh;
                    here for one-off re-runs)
         reads      read-only metrics — RECOMMENDED default, complements
                    Assemblers_v2.sh
         both       reference + reads (legacy "longer" mode)
  -t   CPU threads.  Default: \$SLURM_CPUS_ON_NODE, else 20.
  -c   Memory GB.   Default: \$SLURM_MEM_PER_NODE/1024, else 100.  (Informational;
       TransRate doesn't accept a memory cap, but we export MEM for any user hook.)
  -S   Skip the seqkit sana + seqkit pair sanitization step (default: ON for
       reads/both modes).  Sanitization prevents TransRate's
       "Unmatched read IDs ... Use the -I option to ignore this" failure by
       (1) repairing broken single-line FASTQ records via seqkit sana, then
       (2) discarding orphans via seqkit pair so R1 and R2 share identical IDs.
       Requires seqkit on PATH (https://bioinf.shenwei.me/seqkit/).
  -h   This message.
EOF
}

MANIFEST_FILE=""; QUERY=""; MODE="reads"; CPU_OPT=""; MEM_OPT=""; SANITIZE=1

while getopts "s:d:m:t:c:Sh" opt; do
    case $opt in
        s) MANIFEST_FILE="$OPTARG";;
        d) QUERY="$OPTARG";;
        m) MODE="$OPTARG";;
        t) CPU_OPT="$OPTARG";;
        c) MEM_OPT="$OPTARG";;
        S) SANITIZE=0;;
        h) print_help; exit 0;;
        ?) print_help >&2; exit 1;;
    esac
done
shift $((OPTIND-1))

if [[ -z "$MANIFEST_FILE" || -z "$QUERY" ]]; then
    echo "Error: -s and -d are required." >&2
    print_help >&2
    exit 1
fi
[[ -f "$MANIFEST_FILE" ]] || { echo "Error: manifest not found: $MANIFEST_FILE" >&2; exit 1; }
[[ -s "$QUERY"          ]] || { echo "Error: query FASTA empty/missing: $QUERY" >&2; exit 1; }

case "$MODE" in
    reference|reads|both) ;;
    *) echo "Error: invalid -m '$MODE' (use reference|reads|both)" >&2; exit 1;;
esac

# -----------------------------------------------------------------------------
# Resources (SLURM-aware, CLI override)
# -----------------------------------------------------------------------------
CPU="${CPU_OPT:-${SLURM_CPUS_ON_NODE:-20}}"
if [[ -n "$MEM_OPT" ]]; then
    MEM="$MEM_OPT"
elif [[ -n "${SLURM_MEM_PER_NODE:-}" ]]; then
    MEM=$(( SLURM_MEM_PER_NODE / 1024 )); (( MEM < 1 )) && MEM=1
else
    MEM=100
fi
export CPU MEM

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

# Mirrors Assemblers_v2.sh exactly so stems line up between the two scripts.
manifest_stem() {
    local bs
    bs=$(basename "${1%.*}")
    printf '%s\n' "${bs%_samples}"
}

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

# Populates globals: sids fwds revs refs n REFERENCE
load_manifest() {
    local manifest="$1"
    sids=(); fwds=(); revs=(); refs=()
    while IFS=$'\t' read -r sid fwd rev ref _rest; do
        [[ -z "${sid:-}" ]] && continue
        sids+=("$sid"); fwds+=("$fwd"); revs+=("$rev"); refs+=("$ref")
    done < <(read_manifest_rows "$manifest")
    n=${#sids[@]}
    if (( n == 0 )); then
        echo "Error: manifest $manifest has no data rows." >&2
        exit 1
    fi
    REFERENCE="${refs[0]}"
    local r
    for r in "${refs[@]}"; do
        if [[ "$r" != "$REFERENCE" ]]; then
            echo "  Warning: manifest has mixed reference columns; using $REFERENCE" >&2
            break
        fi
    done
}

# Creates ${stem}_concat_PE1.fq and ${stem}_concat_PE2.fq.
# Symlinks for n==1 (LOSO folds), concatenates otherwise. Mirrors the I/O win
# from Assemblers_v2.sh — single-row manifests never copy raw reads.
stage_reads() {
    local stem="$1"
    local fwd="${stem}_concat_PE1.fq"
    local rev="${stem}_concat_PE2.fq"
    rm -f "$fwd" "$rev"
    if (( n == 1 )); then
        ln -s "${fwds[0]}" "$fwd"
        ln -s "${revs[0]}" "$rev"
        echo "    Reads: symlinked (single-sample LOSO manifest)" >&2
    else
        : > "$fwd"; : > "$rev"
        local f
        for f in "${fwds[@]}"; do [[ -s "$f" ]] && cat "$f" >> "$fwd"; done
        for f in "${revs[@]}"; do [[ -s "$f" ]] && cat "$f" >> "$rev"; done
        echo "    Reads: concatenated ($n samples)" >&2
    fi
}

# Repair broken FASTQ records and enforce matched R1/R2 read IDs.
# Prevents TransRate's "Unmatched read IDs ... Use the -I option to ignore
# this" failure (snap.cpp).  Strategy:
#   1) seqkit sana   : drop or repair broken single-line FASTQ records
#   2) seqkit pair   : keep only reads whose IDs appear in BOTH R1 and R2
# Operates on the staged ${stem}_concat_PE{1,2}.fq files in-place.  If the
# inputs were symlinks (LOSO single-sample case) the original source files
# under fwds[0]/revs[0] are NEVER mutated — we delete the symlink and replace
# it with a fresh sanitized file.
sanitize_reads() {

    EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/
    export PATH=$PATH:$EXPORT

    local stem="$1"
    local fwd="${stem}_concat_PE1.fq"
    local rev="${stem}_concat_PE2.fq"

    if ! command -v seqkit >/dev/null 2>&1; then
        echo "    WARNING: seqkit not on PATH; skipping pair sanitization." >&2
        echo "             If TransRate fails with 'Unmatched read IDs', install seqkit" >&2
        echo "             from https://bioinf.shenwei.me/seqkit/ and re-run." >&2
        return 0
    fi

    local sana_fwd="${stem}_sana_PE1.fq"
    local sana_rev="${stem}_sana_PE2.fq"
    local pair_dir="${stem}_paired_tmp"
    local sana_log="transrate_tmp_dir/${stem}_seqkit.log"
    mkdir -p transrate_tmp_dir
    rm -f "$sana_fwd" "$sana_rev"
    rm -rf "$pair_dir"

    echo "    Sanitizing reads (seqkit sana | seqkit pair) ..." >&2

    # ---- Step 1: seqkit sana — repair broken FASTQ records --------------
    # `sana` reads through the (possibly symlinked) input and writes a clean
    # copy.  -j threads it; --quiet keeps the log tidy.
    if ! seqkit sana --quiet -j "$CPU" "$fwd" -o "$sana_fwd" >"$sana_log"  2>&1 || \
       ! seqkit sana --quiet -j "$CPU" "$rev" -o "$sana_rev" >>"$sana_log" 2>&1
    then
        echo "    WARNING: seqkit sana failed (log: $sana_log); using originals." >&2
        rm -f "$sana_fwd" "$sana_rev"
        return 0
    fi

    # ---- Step 2: seqkit pair — drop orphan reads ------------------------
    # Outputs go to $pair_dir/<basename>.paired.fq.
    if ! seqkit pair --quiet -j "$CPU" \
            -1 "$sana_fwd" -2 "$sana_rev" -O "$pair_dir" >>"$sana_log" 2>&1
    then
        echo "    WARNING: seqkit pair failed (log: $sana_log); using originals." >&2
        rm -f "$sana_fwd" "$sana_rev"
        rm -rf "$pair_dir"
        return 0
    fi

    local paired_fwd="${pair_dir}/$(basename "${sana_fwd%.fq}").paired.fq"
    local paired_rev="${pair_dir}/$(basename "${sana_rev%.fq}").paired.fq"

    if [[ ! -s "$paired_fwd" || ! -s "$paired_rev" ]]; then
        echo "    WARNING: seqkit pair produced empty output; using originals." >&2
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
    echo "    Sanitized pairs retained: $n_pairs   (log: $sana_log)" >&2
}

# Moves contigs.csv and *.blast out of the per-run scratch dir into a stable
# transrate_contigs_dir/ tree, preserving sub-directory layout to prevent
# cross-fold collisions.
collect_outputs() {
    local tr_path="$1"
    local n_csv=0 n_blast=0 contig_file blast_file dst

    while IFS= read -r contig_file; do
        dst="transrate_contigs_dir/${contig_file#transrate_tmp_dir/}"
        mkdir -p "$(dirname "$dst")"
        mv "$contig_file" "$dst"
        n_csv=$(( n_csv + 1 ))
    done < <(find "$tr_path" -name "contigs.csv")

    mkdir -p transrate_contigs_dir/blast_outputs
    while IFS= read -r blast_file; do
        dst="transrate_contigs_dir/blast_outputs/${blast_file#transrate_tmp_dir/}"
        mkdir -p "$(dirname "$dst")"
        mv "$blast_file" "$dst"
        n_blast=$(( n_blast + 1 ))
    done < <(find "$tr_path" -name "*[0-9].blast" 2>/dev/null)

    echo "    Saved: $n_csv contigs.csv, $n_blast blast file(s)."
}

run_transrate() {
    local query="$1"
    local manifest="$2"
    local mode="$3"

    setup_transrate_env
    mkdir -p transrate_tmp_dir transrate_contigs_dir

    local mstem qstem
    mstem=$(manifest_stem "$manifest")
    qstem=$(basename "${query%.*}")

    load_manifest "$manifest"

    # Per-query, per-mode scratch dir keeps multiple assemblies + folds + modes
    # from clobbering each other.
    local TRANSRATE_DIR="${qstem}_${mode}_dir"
    local tr_path="transrate_tmp_dir/$TRANSRATE_DIR"
    local tr_log="transrate_tmp_dir/${qstem}_${mode}.log"

    local fwd="${mstem}_concat_PE1.fq"
    local rev="${mstem}_concat_PE2.fq"
    local need_reads=0
    [[ "$mode" == "reads" || "$mode" == "both" ]] && need_reads=1

    if (( need_reads )); then
        stage_reads "$mstem"
        if (( SANITIZE )); then
            sanitize_reads "$mstem"
        else
            echo "    Skipping read sanitization (-S given)." >&2
        fi
    fi

    local call
    case "$mode" in
        reference) call="transrate --assembly $query --reference $REFERENCE --output $tr_path --threads $CPU";;
        reads)     call="transrate --left $fwd --right $rev --assembly $query --output $tr_path --threads $CPU";;
        both)      call="transrate --left $fwd --right $rev --assembly $query --reference $REFERENCE --output $tr_path --threads $CPU";;
    esac

    echo "  [transrate $mode]  query=$qstem  manifest_stem=$mstem"
    (( need_reads ))            && echo "    LEFT/RIGHT: $fwd / $rev"
    [[ "$mode" != "reads" ]]    && echo "    REFERENCE:  $REFERENCE"
    echo "    $call"
    echo "    log: $tr_log"

    local rc=0
    if ! eval "$call" > "$tr_log" 2>&1; then
        echo "  [transrate $mode] FAILED for $qstem (log: $tr_log)" >&2
        rc=1
    else
        collect_outputs "$tr_path"
    fi

    rm -rf "$tr_path"
    (( need_reads )) && rm -f "$fwd" "$rev"
    return $rc
}

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
echo "============================================================"
echo "Accuracy_v2.sh   mode=$MODE   CPU=$CPU   MEM=${MEM}GB   sanitize=$(( SANITIZE ))"
echo "  Manifest: $MANIFEST_FILE"
echo "  Query:    $QUERY"

run_transrate "$QUERY" "$MANIFEST_FILE" "$MODE"
rc=$?

echo "Done (exit=$rc)."
exit $rc
