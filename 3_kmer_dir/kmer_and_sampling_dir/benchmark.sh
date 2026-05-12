#!/bin/bash
# ============================================================================
# benchmark.sh
# Unified De Novo transcriptome assembly benchmarking pipeline
# Merges the logic of kmer.sh + subsampling.sh into a single configurable tool.
#
# Pipeline (per manifest, per parameter combination):
#   1. Concatenate paired-end FASTQ files listed in the manifest
#   2. (optional) Subsample with seqkit at proportion -p
#   3. Run the chosen assembler (-a)  with the chosen k-mer (-k, spades only)
#   4. Run reference-based transrate; keep contigs.csv in transrate_contigs_dir/
#
# Manifest format (TAB or whitespace-separated, 4 columns):
#   factor    forward.fq[.gz]    reverse.fq[.gz]    reference.fa
#
# Modes:
#   kmer        : sweep k-mer sizes,        assembler defaults to spades
#   subsampling : sweep subsampling props,  pick any supported assembler
#   both        : subsample AND vary k-mer (spades), full factorial sweep
#
# Output directory naming (kept human-readable AND machine-parseable):
#   transrate_contigs_dir/<factor>_<assembler>_k<kmer>_p<prop>_dir/contigs.csv
#
# Examples:
#   ./benchmark.sh -s Fold01.txt -m kmer        -k 49
#   ./benchmark.sh -s all         -m kmer        -k 21
#   ./benchmark.sh -s Fold01.txt -m subsampling -p 0.5 -a trinity
#   ./benchmark.sh -s Fold01.txt -m both        -k 49  -p 0.5
# ============================================================================

set -u
set -o pipefail

# ----- defaults -----
MODE=""
MANIFEST_FILE=""
KMERSIZE=""
PROPORTION=""
ASSEMBLER="spades"
OUTPUT_DIR=""
THREADS=24
MEM_GB=100
RAND_SEED=123

usage() {
    cat <<EOF
Usage: $(basename "$0") -s <Manifest|all> -m <mode> [options]

Required:
  -s <Manifest>   Manifest file (4 cols: factor fwd rev reference).
                  Use 'all' to loop over every *.txt manifest in PWD.
  -m <mode>       One of: kmer | subsampling | both

Mode-dependent:
  -k <int|list>   K-mer size(s) for rnaSPAdes (must be odd, <128).
                  Required for modes 'kmer' and 'both'.
  -p <float>      Subsampling proportion (0,1].
                  Required for modes 'subsampling' and 'both'.
  -a <assembler>  spades | trinity | idba | megahit | rnabloom
                  (Defaults to spades. Only spades honors -k.)

Optional:
  -o <dir>        Output root (default: <manifest_basename>_dir or
                  multiple_analysis_dir when -s all).
  -t <int>        Threads (default: $THREADS).
  -M <int>        Max memory GB (default: $MEM_GB).
  -h              Show this help.
EOF
}

while getopts "s:m:k:p:a:o:t:M:h" opt; do
    case $opt in
        s) MANIFEST_FILE="$OPTARG" ;;
        m) MODE="$OPTARG" ;;
        k) KMERSIZE="$OPTARG" ;;
        p) PROPORTION="$OPTARG" ;;
        a) ASSEMBLER="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        M) MEM_GB="$OPTARG" ;;
        h) usage; exit 0 ;;
        ?) usage; exit 1 ;;
    esac
done

# ----- argument validation -----
if [[ -z "$MANIFEST_FILE" || -z "$MODE" ]]; then
    echo "ERROR: -s and -m are required." >&2; usage; exit 1
fi

case "$MODE" in
    kmer)
        [[ -z "$KMERSIZE" ]] && { echo "ERROR: -k required for mode kmer" >&2; exit 1; }
        ASSEMBLER="spades"
        ;;
    subsampling)
        [[ -z "$PROPORTION" ]] && { echo "ERROR: -p required for mode subsampling" >&2; exit 1; }
        ;;
    both)
        [[ -z "$KMERSIZE" || -z "$PROPORTION" ]] && {
            echo "ERROR: both -k and -p required for mode both" >&2; exit 1; }
        ASSEMBLER="spades"
        ;;
    *) echo "ERROR: unknown mode '$MODE'" >&2; usage; exit 1 ;;
esac

case "$ASSEMBLER" in
    spades|trinity|idba|megahit|rnabloom) ;;
    *) echo "ERROR: unknown assembler '$ASSEMBLER'" >&2; exit 1 ;;
esac

# ============================================================================
# Assembler wrappers
# ============================================================================

run_spades() {
    local forward_fq="$1" reverse_fq="$2" output="$3" final_fasta="$4" kmer_sizes="$5"
    export PATH=/LUSTRE/apps/bioinformatica/SPAdes-3.15.5-Linux/bin/:$PATH

    local k_flag=""
    [[ -n "$kmer_sizes" ]] && k_flag="-k $kmer_sizes"

    local call="rnaspades.py $k_flag -1 $forward_fq -2 $reverse_fq -o $output -t $THREADS -m $MEM_GB"
    echo "$call"; eval "$call"

    local f1
    f1=$(find "$output" -maxdepth 1 -type f -name 'transcripts.fasta')
    echo "Assembly: $f1"
    mv "$f1" "$final_fasta"
}

run_trinity() {
    local forward_fq="$1" reverse_fq="$2" output="$3" final_fasta="$4"
    module load trinityrnaseq-v2.15.1
    local output_dir="${output}/Trinity_out_dir"

    local call="Trinity --seqType fq --max_memory ${MEM_GB}G --left $forward_fq --right $reverse_fq \
        --no_normalize_reads --CPU $THREADS --output $output_dir --full_cleanup"
    echo "$call"; eval "$call"

    local f1
    f1=$(find "$output" -type f -name 'Trinity_out_dir.Trinity.fasta')
    mv "$f1" "$final_fasta"
}

run_idba() {
    local forward_fq="$1" reverse_fq="$2" output="$3" final_fasta="$4"
    export PATH=/LUSTRE/apps/bioinformatica/idba-1.1.3/bin:$PATH

    local merged="${output}/idba_merged.fa"
    fq2fa --merge "$forward_fq" "$reverse_fq" "$merged"
    local call="idba_tran -r $merged -o $output --num_threads $THREADS"
    echo "$call"; eval "$call"

    mv "${output}/contig.fa" "$final_fasta"
    rm -f "$merged"
}

run_megahit() {
    local forward_fq="$1" reverse_fq="$2" output="$3" final_fasta="$4"
    module load megahit || true

    rm -rf "$output/megahit_out"
    local call="megahit -1 $forward_fq -2 $reverse_fq -o $output/megahit_out -t $THREADS"
    echo "$call"; eval "$call"

    mv "$output/megahit_out/final.contigs.fa" "$final_fasta"
}

run_rnabloom() {
    local forward_fq="$1" reverse_fq="$2" output="$3" final_fasta="$4"
    module load conda-2025
    source activate base
    conda activate rnabloom

    local call="rnabloom -l $forward_fq -r $reverse_fq -t $THREADS -outdir $output -mem $MEM_GB"
    echo "$call"; eval "$call"

    local f1
    f1=$(find "$output" -maxdepth 1 -type f -name 'rnabloom.transcripts.fa')
    mv "$f1" "$final_fasta"
}

dispatch_assembler() {
    # $1 fwd, $2 rev, $3 outdir, $4 final_fasta, $5 kmer_sizes (spades only)
    case "$ASSEMBLER" in
        spades)   run_spades   "$1" "$2" "$3" "$4" "$5" ;;
        trinity)  run_trinity  "$1" "$2" "$3" "$4" ;;
        idba)     run_idba     "$1" "$2" "$3" "$4" ;;
        megahit)  run_megahit  "$1" "$2" "$3" "$4" ;;
        rnabloom) run_rnabloom "$1" "$2" "$3" "$4" ;;
    esac
}

# ============================================================================
# Shared transrate (kept identical to the original implementation)
# ============================================================================

run_reference_transrate() {
    export PATH=/LUSTRE/apps/bioinformatica/.local/bin:$PATH
    export PATH=/LUSTRE/apps/bioinformatica/ruby2/ruby-2.2.0/bin:$PATH
    export LD_LIBRARY_PATH=/LUSTRE/apps/bioinformatica/ruby2/ruby-2.2.0/lib:${LD_LIBRARY_PATH:-}

    mkdir -p transrate_tmp_dir transrate_contigs_dir

    local QUERY="$1"
    local TARGET; TARGET=$(awk '{print $4}' "$2" | head -n1)
    local BS; BS=$(basename "${QUERY%.*}")
    local TRANSRATE_DIR="${BS}_dir"

    echo "Transrate: $QUERY vs $TARGET -> transrate_tmp_dir/$TRANSRATE_DIR"
    transrate --assembly "$QUERY" --reference "$TARGET" \
              --output "transrate_tmp_dir/$TRANSRATE_DIR" --threads "$THREADS" &> /dev/null

    local contig_file dst
    contig_file=$(find "transrate_tmp_dir/$TRANSRATE_DIR" -name "contigs.csv")
    dst="transrate_contigs_dir/${contig_file#transrate_tmp_dir/}"
    mkdir -p "$(dirname "$dst")"
    mv "$contig_file" "$dst"
    rm -rf "transrate_tmp_dir/$TRANSRATE_DIR"
}

# ============================================================================
# Core per-manifest workflow
# ============================================================================

# Build a stable run-id used in directory names (so EDA can parse all params)
#   <factor>_<assembler>_k<kmer-or-NA>_p<prop-or-NA>
build_runid() {
    local manifest="$1" k="$2" p="$3"
    local factor; factor=$(basename "${manifest%.*}")
    local kpart="kNA" ppart="pNA"
    [[ -n "$k" ]] && kpart="k${k//,/-}"   # keep multi-k lists encoded with '-'
    [[ -n "$p" ]] && ppart="p${p}"
    echo "${factor}_${ASSEMBLER}_${kpart}_${ppart}"
}

# One assembly + transrate cycle. $1 manifest, $2 kmer (or ""), $3 proportion (or "")
run_one() {
    local manifest="$1" kmer="$2" prop="$3" outroot="$4"
    local runid; runid=$(build_runid "$manifest" "$kmer" "$prop")

    local WORK="${outroot}/${runid}_work"
    local FASTA_DIR="${outroot}/${ASSEMBLER}_fasta_dir"
    mkdir -p "$WORK" "$FASTA_DIR"

    local fwd="${WORK}/${runid}_PE1.fq"
    local rev="${WORK}/${runid}_PE2.fq"
    local final_fasta="${FASTA_DIR}/${runid}.fa"

    echo "==> [${runid}] concatenating reads"
    awk '{print $2}' "$manifest" | xargs cat > "$fwd"
    awk '{print $3}' "$manifest" | xargs cat > "$rev"

    # optional subsampling
    if [[ -n "$prop" ]]; then
        echo "==> [${runid}] seqkit sample -p $prop"
        local fwd_s="${WORK}/${runid}_sampled_PE1.fq"
        local rev_s="${WORK}/${runid}_sampled_PE2.fq"
        seqkit sample --proportion "$prop" --rand-seed "$RAND_SEED" -o "$fwd_s" "$fwd"
        seqkit sample --proportion "$prop" --rand-seed "$RAND_SEED" -o "$rev_s" "$rev"
        rm -f "$fwd" "$rev"
        fwd="$fwd_s"; rev="$rev_s"
    fi

    echo "==> [${runid}] running ${ASSEMBLER}"
    dispatch_assembler "$fwd" "$rev" "$WORK" "$final_fasta" "$kmer" &> /dev/null

    echo "==> [${runid}] transrate"
    run_reference_transrate "$final_fasta" "$manifest" &> /dev/null

    rm -rf "$WORK"
    echo "==> [${runid}] DONE"
}

# ============================================================================
# Mode dispatch
# ============================================================================

run_for_manifest() {
    local manifest="$1"
    local outroot="$2"
    mkdir -p "$outroot"

    case "$MODE" in
        kmer)
            for k in ${KMERSIZE//,/ }; do
                run_one "$manifest" "$k" "" "$outroot"
            done
            ;;
        subsampling)
            for p in ${PROPORTION//,/ }; do
                run_one "$manifest" "" "$p" "$outroot"
            done
            ;;
        both)
            for k in ${KMERSIZE//,/ }; do
                for p in ${PROPORTION//,/ }; do
                    run_one "$manifest" "$k" "$p" "$outroot"
                done
            done
            ;;
    esac
}

# ----- main -----

if [[ "$MANIFEST_FILE" == "all" ]]; then
    OUTDIR="${OUTPUT_DIR:-multiple_analysis_dir}"
    echo "Running batch mode over $PWD/*.txt -> $OUTDIR"
    find "$PWD" -maxdepth 1 -name '*.txt' | while read -r m; do
        run_for_manifest "$m" "$OUTDIR"
    done
else
    OUTDIR="${OUTPUT_DIR:-${MANIFEST_FILE%.*}_dir}"
    run_for_manifest "$MANIFEST_FILE" "$OUTDIR"
fi

rm -rf transrate_tmp_dir
echo "All runs finished."
exit 0
