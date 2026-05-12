#!/bin/bash
# =============================================================================
# Run_LOSO_pipeline.sh
#
# End-to-end LOSO orchestrator for Strategy 1.
#
# Pipeline:
#   1. (Once)  R: LOSO_conoServerDB.R          -> train_*.fasta, test_*.fasta, manifest
#   2.         bash: Simulate_rnaseq_loso.sh   -> simulated reads from test_*.fasta only
#   3.         bash: Assemblers.sh             -> per-fold assemblies + TransRate
#                                                (ground truth = test_*.fasta)
#   4.         bash: leakage check             -> BLAST assembly vs. train_*.fasta to flag
#                                                contigs that map to training (should be rare)
#   5. (Once)  R: LOSO_aggregate.R             -> per-superfamily sensitivity/precision
#
# Usage:
#   ./Run_LOSO_pipeline.sh -d INPUTS/vfolds_loso_resampling_dir -c run.config
#
# Prerequisites:
#   - LOSO_conoServerDB.R has already been run and produced train/test FASTAs
#   - Simulate_rnaseq_loso.sh, Assemblers.sh, and run.config are in PATH or PWD
#   - blastn or makeblastdb available for the leakage check (optional, can be skipped)
# =============================================================================

set -e

CONFIG=""
LOSO_DIR=""
SKIP_LEAKAGE=0

while getopts "d:c:Sh" opt; do
    case $opt in
        d) LOSO_DIR="$OPTARG";;
        c) CONFIG="$OPTARG";;
        S) SKIP_LEAKAGE=1;;
        h)
            echo "Usage: $(basename $0) -d <loso_dir> -c <assembler_config> [-S]"
            echo
            echo "  -d   Path to LOSO output dir from LOSO_conoServerDB.R"
            echo "  -c   Assembler config file (same format as Assemblers.sh expects)"
            echo "  -S   Skip the BLAST-based leakage check (default: run it)"
            exit 0
            ;;
        ?) echo "Invalid option: -$OPTARG" >&2; exit 1;;
    esac
done
shift $((OPTIND-1))

if [[ -z "$LOSO_DIR" || -z "$CONFIG" ]]; then
    echo "Error: missing required arguments. Run with -h for help." >&2
    exit 1
fi

# -----------------------------------------------------------------------------
# Step 2: Simulate reads (only from test_*.fasta)
# -----------------------------------------------------------------------------
echo "==========================================================="
echo "STEP 2/5  Simulating reads from held-out superfamilies..."
echo "==========================================================="

./Simulate_rnaseq_loso.sh -d "$LOSO_DIR"

# Most-recently-created _loso_dir is our simulation output
SIM_DIR=$(ls -td *_loso_dir 2>/dev/null | head -n 1)
if [[ -z "$SIM_DIR" ]]; then
    echo "ERROR: no *_loso_dir found after simulation step" >&2
    exit 1
fi
echo "Simulation outputs in: $SIM_DIR"

# -----------------------------------------------------------------------------
# Step 3: Assemble each LOSO fold and score against held-out ground truth
# -----------------------------------------------------------------------------
echo "==========================================================="
echo "STEP 3/5  Running assemblers on each LOSO fold..."
echo "==========================================================="

# Iterate over each fold's manifest in the simulation output directory
find "$SIM_DIR" -maxdepth 1 -name "*_samples.txt" | while read -r manifest; do
    echo
    echo "--- Assembling fold: $(basename "$manifest") ---"
    ./Assemblers.sh -s "$manifest" -c "$CONFIG"
done

# -----------------------------------------------------------------------------
# Step 4: Leakage check (optional)
# -----------------------------------------------------------------------------
if [[ "$SKIP_LEAKAGE" -eq 0 ]]; then
    echo "==========================================================="
    echo "STEP 4/5  Leakage check: BLAST assemblies vs. training refs"
    echo "==========================================================="

    LEAKAGE_DIR="${SIM_DIR}/leakage_check"
    mkdir -p "$LEAKAGE_DIR"

    META="$SIM_DIR/LOSO_metadata.tsv"
    if [[ ! -f "$META" ]]; then
        echo "WARNING: LOSO_metadata.tsv not found in $SIM_DIR; skipping leakage check" >&2
    else
        # For each (sample, train_ref, sf, fcov) row, BLAST the assembly contigs
        # of that sample against the corresponding training reference. Any high-
        # identity hit is unexpected under LOSO and is recorded for inspection.
        tail -n +2 "$META" | while IFS=$'\t' read -r sample_id train_ref sf_tag fcov; do
            if [[ "$train_ref" == "NA" || ! -f "$train_ref" ]]; then
                continue
            fi

            # Locate this sample's assemblies. Assemblers.sh creates
            # <sample_id>_<tool>_dir/<sample_id>_<tool>.fa under the working dir.
            for asm in ./*_FASTA_DIR/${sample_id}_*.fa; do
                [[ ! -f "$asm" ]] && continue

                asm_bs=$(basename "$asm" .fa)
                db_prefix="$LEAKAGE_DIR/db_${sf_tag}"
                hit_file="$LEAKAGE_DIR/${asm_bs}_vs_train_${sf_tag}.tsv"

                # Build BLAST DB if not yet present
                if [[ ! -f "${db_prefix}.nhr" ]]; then
                    makeblastdb -in "$train_ref" -dbtype nucl -out "$db_prefix" \
                        > /dev/null 2>&1 || { echo "makeblastdb failed for $train_ref"; continue; }
                fi

                # Run BLAST. We flag hits with pident>=90 and qcov>=80 as "leakage candidates":
                # contigs that look like training-set sequences. In a leak-free LOSO run, the
                # number of such hits should be small (only via paralog overlap).
                blastn -query "$asm" -db "$db_prefix" \
                    -outfmt "6 qseqid sseqid pident length qcovs evalue bitscore" \
                    -max_target_seqs 1 -evalue 1e-10 \
                    -out "$hit_file" 2>/dev/null

                # Summary line
                n_total=$(wc -l < "$hit_file" 2>/dev/null || echo 0)
                n_high=$(awk -F'\t' '$3>=90 && $5>=80' "$hit_file" 2>/dev/null | wc -l)
                printf "%s\t%s\t%d\t%d\n" "$asm_bs" "$sf_tag" "$n_total" "$n_high" \
                    >> "$LEAKAGE_DIR/leakage_summary.tsv"
            done
        done

        if [[ -f "$LEAKAGE_DIR/leakage_summary.tsv" ]]; then
            # Prepend header if not already present
            if ! head -n 1 "$LEAKAGE_DIR/leakage_summary.tsv" | grep -q "^assembly"; then
                sed -i '1i assembly\theld_out_sf\ttotal_hits\thits_pident90_qcov80' \
                    "$LEAKAGE_DIR/leakage_summary.tsv"
            fi
            echo "Leakage summary: $LEAKAGE_DIR/leakage_summary.tsv"
        fi
    fi
fi

# -----------------------------------------------------------------------------
# Step 5: Hand off to R for aggregation
# -----------------------------------------------------------------------------
echo "==========================================================="
echo "STEP 5/5  Run LOSO_aggregate.R for per-superfamily summary"
echo "==========================================================="
echo
echo "Next step (manual):"
echo "  Rscript LOSO_aggregate.R \\"
echo "    --loso_dir $LOSO_DIR \\"
echo "    --sim_dir  $SIM_DIR \\"
echo "    --out_dir  $SIM_DIR/LOSO_summary"
echo
echo "Pipeline complete."
exit 0
