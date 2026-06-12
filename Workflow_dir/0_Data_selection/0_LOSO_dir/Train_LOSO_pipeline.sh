#!/bin/bash
# =============================================================================
# Train_LOSO_pipeline.sh
#
# Post-simulation leakage-check step for the LOSO benchmark.
# Consumes the outputs produced by Simulate_rnaseq_loso.sh:
#   - LOSO_metadata.tsv   (sample_id  train_ref  superfamily  coverage)
#   - Assembly FASTAs      produced by Assemblers.sh under <ASM_DIR>/
#
# For each (sample, train_ref) pair in the metadata, BLAST the assembled
# contigs against the corresponding training-set FASTA. High-identity hits
# (pident ≥ 90, qcov ≥ 80) are flagged as leakage candidates and written
# to a leakage_summary.tsv inside <ASM_DIR>/leakage_check/.
#
# Usage:
#   ./Train_LOSO_pipeline.sh \
#       -m <path/to/LOSO_metadata.tsv> \
#       -d <path/to/assemblies_FASTA_dir>
#
#   ./Train_LOSO_pipeline.sh \
#       -m 9f3c_loso_dir/LOSO_metadata.tsv \
#       -d assemblies_dir \
#       --skip-leakage
#
# Dependencies: blastn, makeblastdb (NCBI BLAST+)
# =============================================================================

# Defaults
META=""
ASM_DIR=""
SKIP_LEAKAGE=0

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        -m|--metadata)
            META="$2"; shift 2;;
        -d|--fasta-dir)
            ASM_DIR="$2"; shift 2;;
        --skip-leakage)
            SKIP_LEAKAGE=1; shift;;
        -h|--help)
            echo "Usage: $(basename "$0") -m <LOSO_metadata.tsv> -d <assemblies_dir> [--skip-leakage]"
            echo
            echo "  -m, --metadata       Path to LOSO_metadata.tsv emitted by Simulate_rnaseq_loso.sh"
            echo "  -d, --fasta-dir      Directory that contains assembly FASTAs"
            echo "                       (Assemblers.sh writes <sample_id>_<tool>.fa here)"
            echo "  --skip-leakage       Skip the BLAST-based leakage check (Step 4)"
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2; exit 1;;
    esac
done

# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------
if [[ -z "$META" ]]; then
    echo "Error: -m / --metadata is required. Run with -h for help." >&2
    exit 1
fi
if [[ ! -f "$META" ]]; then
    echo "Error: metadata file not found: $META" >&2
    exit 1
fi

if [[ -z "$ASM_DIR" ]]; then
    echo "Error: -d / --fasta-dir is required. Run with -h for help." >&2
    exit 1
fi
if [[ ! -d "$ASM_DIR" ]]; then
    echo "Error: assembly directory not found: $ASM_DIR" >&2
    exit 1
fi

echo "LOSO metadata:     $META"
echo "Assembly FASTA dir: $ASM_DIR"
echo "Skip leakage check: $SKIP_LEAKAGE"

EXPORT=/LUSTRE/apps/bioinformatica/ncbi-blast-2.14.0+/bin/
export PATH=$PATH:$EXPORT

# ---------------------------------------------------------------------------
# Step 4: Leakage check
# ---------------------------------------------------------------------------
if [[ "$SKIP_LEAKAGE" -eq 0 ]]; then
    echo "==========================================================="
    echo "STEP 4/5  Leakage check: BLAST assemblies vs. training refs"
    echo "==========================================================="

    LEAKAGE_DIR="${ASM_DIR}/leakage_check"
    mkdir -p "$LEAKAGE_DIR"

    # Columns: fold_id  fold_tag  held_out_superfamily  n_test  n_train  test_fasta  train_fasta
    tail -n +2 "$META" | while IFS=$'\t' read -r fold_id fold_tag sf_tag n_test n_train test_fasta train_ref; do
        if [[ "$train_ref" == "NA" || ! -f "$train_ref" ]]; then
            echo "WARNING: skipping fold $fold_id — train_ref not found: $train_ref" >&2
            continue
        fi

        # Derive sample prefix from test FASTA basename (test_B1_superfamily)
        # Assemblers.sh creates <ASM_DIR>/<sample_id>_<tool>.fa
        sample_id=$(basename "$test_fasta" .fasta)
        for asm in "$ASM_DIR"/${sample_id}_*.fa; do
            [[ ! -f "$asm" ]] && continue

            echo "Assembly: $asm \n"
            echo "Train: $train_ref \n"

            asm_bs=$(basename "$asm" .fa)
            db_prefix="$LEAKAGE_DIR/db_${sf_tag//[$'\t\r\n ']}"
            hit_file="$LEAKAGE_DIR/${asm_bs}_vs_train_${sf_tag}.tsv"

            echo "DB prefix: $db_prefix \n"
            echo "Hit file: $hit_file \n"

            # Build BLAST DB once per superfamily
            if [[ ! -f "${db_prefix}.nhr" ]]; then
                echo "CMD: makeblastdb -in $train_ref -dbtype nucl -out $db_prefix"
                makeblastdb -in "$train_ref" -dbtype nucl -out "$db_prefix" \
                    || { echo "makeblastdb failed for $train_ref" >&2; continue; }
            fi

            # Hits with pident>=90 and qcov>=80 are "leakage candidates":
            # contigs that resemble training-set sequences. In a leak-free LOSO
            # run, such hits should be rare (only via genuine paralog overlap).
            echo "CMD: blastn -query $asm -db $db_prefix -outfmt '6 ...' -out $hit_file"
            blastn \
                -query  "$asm" \
                -db     "$db_prefix" \
                -outfmt "6 qseqid sseqid pident length qcovs evalue bitscore" \
                -max_target_seqs 1 \
                -evalue 1e-10 \
                -out    "$hit_file"

            n_total=$(wc -l < "$hit_file" 2>/dev/null || echo 0)
            n_high=$(awk -F'\t' '$3>=90 && $5>=80' "$hit_file" 2>/dev/null | wc -l)
            printf "%s\t%s\t%d\t%d\n" "$asm_bs" "$sf_tag" "$n_total" "$n_high" \
                >> "$LEAKAGE_DIR/leakage_summary.tsv"
        done
    done

    if [[ -f "$LEAKAGE_DIR/leakage_summary.tsv" ]]; then
        if ! head -n 1 "$LEAKAGE_DIR/leakage_summary.tsv" | grep -q "^assembly"; then
            sed -i '1i assembly\theld_out_sf\ttotal_hits\thits_pident90_qcov80' \
                "$LEAKAGE_DIR/leakage_summary.tsv"
        fi
        echo "Leakage summary: $LEAKAGE_DIR/leakage_summary.tsv"
    else
        echo "No leakage results written (no matching assemblies found under $ASM_DIR)." >&2
    fi
fi

echo
echo "Train LOSO pipeline complete."
exit 0
