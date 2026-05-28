#!/bin/bash
#SBATCH --job-name=conotoxin_strat3
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00
#SBATCH -o out.%J
#SBATCH -e err.%J
# =============================================================================
# run_clustering_batch.slurm  (Strategy 3 driver)
#
# Drives the Sequence-Identity-Stratified Hold-Out workflow:
#   1. cluster_toxins_mmseqs.sh  -> reference_*/test_* FASTAs + nucleotide NRI
#   2. (optional) protein-level NRI as a robustness axis
#   3. prints the hand-off command for simulation + assembly
#
# Unlike the original batch script (which looped over many FASTA files and
# only counted clusters), Strategy 3 takes ONE curated ConoServer FASTA and
# produces stratified train/test splits. The "batch" dimension here is the
# set of identity thresholds, handled inside cluster_toxins_mmseqs.sh.
#
# Usage:
#   sbatch run_clustering_batch.slurm <curated_conoserver_nuc.fasta>
#   sbatch run_clustering_batch.slurm <nuc.fasta> <protein.fasta>   # adds protein NRI
# =============================================================================

set -euo pipefail

# ---- environment ----
EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/mmseqs/bin
export PATH=$PATH:$EXPORT

WORK_DIR=$PWD
cd "$WORK_DIR"

CLUSTER_SCRIPT="./cluster_toxins_mmseqs.sh"
THREADS="${SLURM_CPUS_ON_NODE:-24}"
SEED=20260520
THRESHOLDS="1.00 0.99 0.98 0.97 0.96 0.95 0.94 0.93 0.92 0.91 0.90 0.80 0.70 0.60 0.50 0.40 0.3"

# ---- arguments ----
if [[ $# -lt 1 ]]; then
    echo "Usage: sbatch $0 <curated_conoserver_nuc.fasta> [protein.fasta]" >&2
    exit 1
fi
NUC_FASTA="$1"
PROT_FASTA="${2:-}"

if [[ ! -f "$CLUSTER_SCRIPT" ]]; then
    echo "ERROR: $CLUSTER_SCRIPT not found in $WORK_DIR" >&2
    exit 1
fi
if [[ ! -f "$NUC_FASTA" ]]; then
    echo "ERROR: nucleotide FASTA not found: $NUC_FASTA" >&2
    exit 1
fi
chmod +x "$CLUSTER_SCRIPT"

BASE_NAME=$(basename "$NUC_FASTA" | sed 's/\.[^.]*$//')
OUTDIR="${BASE_NAME}_strat3_dir"

echo "=========================================="
echo "Strategy 3 driver"
echo "  nucleotide FASTA : $NUC_FASTA"
echo "  protein FASTA    : ${PROT_FASTA:-<none>}"
echo "  thresholds       : $THRESHOLDS"
echo "  output dir       : $OUTDIR"
echo "  started          : $(date)"
echo "=========================================="

# -----------------------------------------------------------------------------
# Step 1: stratified hold-out splitting + nucleotide NRI
# -----------------------------------------------------------------------------
echo "[1] running stratified splitter..."
"$CLUSTER_SCRIPT" \
    -i "$NUC_FASTA" \
    -t "$THRESHOLDS" \
    -o "$OUTDIR" \
    -p "$THREADS" \
    -s "$SEED"

# -----------------------------------------------------------------------------
# Step 2 (optional): protein-level NRI for the held-out test sequences
# -----------------------------------------------------------------------------
# For a hypervariable gene family, nucleotide identity can be depressed by
# synonymous-site saturation. Recomputing NRI at the protein level provides a
# robustness axis. We translate test/reference precursors and re-search.
if [[ -n "$PROT_FASTA" && -f "$PROT_FASTA" ]]; then
    echo "[2] computing protein-level NRI..."
    PWORK="$OUTDIR/protein_nri_work"
    mkdir -p "$PWORK"
    MMLOG="$OUTDIR/protein_nri_mmseqs.log"
    > "$MMLOG"

    PROT_NRI="$OUTDIR/nearest_reference_identity_protein.tsv"
    printf "threshold\ttest_id\tnearest_ref_id\tnri_prot\taln_len\n" > "$PROT_NRI"

    for thr in $THRESHOLDS; do
        thr_str=$(printf "%.2f" "$thr")
        test_fa="$OUTDIR/test_${thr_str}.fasta"
        ref_fa="$OUTDIR/reference_${thr_str}.fasta"

        # report which files are present/absent instead of silently continuing
        if [[ ! -s "$test_fa" || ! -s "$ref_fa" ]]; then
            echo "  [thr=$thr_str] SKIP — missing or empty split files:" \
                 "(test: $([ -s "$test_fa" ] && echo ok || echo MISSING)," \
                 "ref: $([ -s "$ref_fa"  ] && echo ok || echo MISSING))"
            continue
        fi

        # IDs present in this threshold's test / reference splits
        grep '^>' "$test_fa" | sed 's/^>//; s/ .*//' > "$PWORK/test_ids_${thr_str}.txt"
        grep '^>' "$ref_fa"  | sed 's/^>//; s/ .*//' > "$PWORK/ref_ids_${thr_str}.txt"

        # subset the protein FASTA to those IDs
        for role in test ref; do
            awk 'BEGIN{while((getline l < ARGV[1])>0) keep[l]=1; ARGV[1]="";p=0}
                 /^>/{id=substr($1,2); p=(id in keep)}
                 {if(p) print}' \
                 "$PWORK/${role}_ids_${thr_str}.txt" "$PROT_FASTA" \
                 > "$PWORK/${role}_${thr_str}.faa"
        done

        n_test=$(grep -c '^>' "$PWORK/test_${thr_str}.faa" 2>/dev/null || echo 0)
        n_ref=$(grep -c  '^>' "$PWORK/ref_${thr_str}.faa"  2>/dev/null || echo 0)
        echo "  [thr=$thr_str] protein subsets: test=$n_test ref=$n_ref sequences"

        # guard: skip mmseqs if protein IDs didn't match nucleotide IDs
        if [[ "$n_test" -eq 0 || "$n_ref" -eq 0 ]]; then
            echo "  [thr=$thr_str] SKIP — no matching protein IDs (check that nuc/prot FASTAs share IDs)"
            continue
        fi

        # protein dbs (--dbtype 1 = amino acid) and search (--search-type 1 = amino acid)
        mmseqs createdb "$PWORK/ref_${thr_str}.faa"  "$PWORK/refdb_${thr_str}"  --dbtype 1 >>"$MMLOG" 2>&1
        mmseqs createdb "$PWORK/test_${thr_str}.faa" "$PWORK/testdb_${thr_str}" --dbtype 1 >>"$MMLOG" 2>&1
        mmseqs search "$PWORK/testdb_${thr_str}" "$PWORK/refdb_${thr_str}" \
            "$PWORK/res_${thr_str}" "$PWORK/tmp_${thr_str}" \
            --search-type 1 -s 7.5 --threads "$THREADS" >>"$MMLOG" 2>&1
        mmseqs convertalis "$PWORK/testdb_${thr_str}" "$PWORK/refdb_${thr_str}" \
            "$PWORK/res_${thr_str}" "$PWORK/alis_${thr_str}.tsv" \
            --format-output "query,target,fident,alnlen" >>"$MMLOG" 2>&1

        n_hits=$(wc -l < "$PWORK/alis_${thr_str}.tsv" 2>/dev/null || echo 0)
        echo "  [thr=$thr_str] alignment hits: $n_hits"

        # best hit per test sequence
        awk -F'\t' -v thr="$thr_str" '
            { if (!($1 in bestf) || $3>bestf[$1]) { bestf[$1]=$3; bestt[$1]=$2; bestl[$1]=$4 } }
            END { for (q in bestf) printf "%s\t%s\t%s\t%.4f\t%s\n", thr, q, bestt[q], bestf[q], bestl[q] }
        ' "$PWORK/alis_${thr_str}.tsv" >> "$PROT_NRI"
    done
    rm -rf "$PWORK"
    echo "    protein NRI table: $PROT_NRI"
    echo "    mmseqs log:        $MMLOG"
fi

# -----------------------------------------------------------------------------
# Step 3: hand-off
# -----------------------------------------------------------------------------
echo
echo "=========================================="
echo "Strategy 3 splitting complete: $(date)"
echo
echo "Outputs in $OUTDIR/:"
echo "  reference_<thr>.fasta              reference/training sets"
echo "  test_<thr>.fasta                   held-out test sets"
echo "  nearest_reference_identity.tsv     nucleotide NRI  (primary x-axis)"
if [[ -n "$PROT_FASTA" ]]; then
echo "  nearest_reference_identity_protein.tsv  protein NRI (robustness axis)"
fi
echo "  cluster_split_<thr>.tsv            cluster -> rep/test mapping"
echo
echo "Next steps:"
echo "  1. Simulate reads from each test_<thr>.fasta:"
echo "       ./Simulate_rnaseq_loso.sh -d $OUTDIR -p test_"
echo "  2. Assemble each fold:"
echo "       ./Assemblers.sh -s all -c run.config"
echo "  3. Plot recovery vs NRI:"
echo "       Rscript analyze_mmseq_clustering_results.R --strat3_dir $OUTDIR \\"
echo "               --runlog_glob '*_runlog.tsv'"
echo "=========================================="
exit 0
