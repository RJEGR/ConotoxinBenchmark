#!/bin/bash
# =============================================================================
# Simulate_rnaseq_loso.sh
#
# LOSO-aware adaptation of Simulate_rnaseq.sh.
# Simulates reads ONLY from FASTAs whose basename starts with "test_"
# (the held-out superfamilies emitted by LOSO_conoServerDB.R).
# Writes two manifests per simulated dataset:
#   - <output_bs>_samples.txt    : 4-column, backward compatible with Assemblers.sh
#                                  (sample_id  fwd  rev  test_ref_for_transrate)
#   - LOSO_metadata.tsv          : sidecar with (sample_id  train_ref  superfamily  fcov)
#
# Why two manifests?
#   - Assemblers.sh + TransRate score the assembly against the test_*.fasta
#     (because reads were simulated from it = ground truth).
#   - Downstream LOSO scoring (ConoSorter, BLAST DBs, optional StringTie reference)
#     needs the *training* set per fold. That lives in train_<sf>.fasta and is
#     tracked in the sidecar.
#
# Usage:
#   ./Simulate_rnaseq_loso.sh -d INPUTS/vfolds_loso_resampling_dir
#   ./Simulate_rnaseq_loso.sh -d INPUTS/vfolds_loso_resampling_dir -p test_ -c "10 50 200"
#
# =============================================================================

# Note: on omicas/LUSTRE, may require `module load conda-2024` first
EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/art_bin_MountRainier
export PATH=$PATH:$EXPORT

# Defaults
PREFIX_FILTER="test_"
COVERAGES="200"

while getopts "d:p:c:h" opt; do
    case $opt in
        d) FASTA_DIR="$OPTARG";;
        p) PREFIX_FILTER="$OPTARG";;
        c) COVERAGES="$OPTARG";;
        h)
            echo "Usage: $(basename $0) -d <fasta_folder> [-p <prefix>] [-c \"<coverage list>\"]"
            echo
            echo "  -d   Path to LOSO output dir (contains test_*.fasta and train_*.fasta)"
            echo "  -p   Process only FASTAs whose basename starts with this prefix.   Default: test_"
            echo "  -c   Space-separated list of coverages (fold-coverage values).     Default: \"5 10 50 100 200\""
            exit 0
            ;;
        ?) echo "Invalid option: -$OPTARG" >&2; exit 1;;
    esac
done
shift $((OPTIND-1))

WD="$FASTA_DIR"
if [[ -z "$WD" ]]; then
    echo "Error: missing -d argument. Run with -h for help." >&2
    exit 1
fi
if [[ ! -d "$WD" ]]; then
    echo "Error: directory does not exist: $WD" >&2
    exit 1
fi

# -----------------------------------------------------------------------------
# run_depth(): identical to Simulate_rnaseq.sh, kept inline for portability
# -----------------------------------------------------------------------------
run_depth() {
    EXPORT=/LUSTRE/apps/bioinformatica/RSEM/bin/samtools-1.3/
    export PATH=$PATH:$EXPORT

    local reference="$1"
    local sam_file="$2"

    local bs="${sam_file%*.sam}"
    local bam_index_file="${bs}.bam"
    local bam_file="${bs}.sort.bam"
    local touch_file="${bs}.depth.chkpt"

    if [ ! -f "$touch_file" ]; then
        samtools view -b -S "$sam_file" > "$bam_index_file"
        samtools sort -@ 20 "$bam_index_file" -o "$bam_file"
        samtools flagstat "$bam_file" > "${bs}.flagstats.txt"
        samtools depth -a --reference "$reference" "$bam_file" > "${bs}.depth.txt"
        touch "$touch_file"
        rm -f "$bam_index_file" "$bam_file" "$sam_file"
    else
        echo "Depth checkpoint already exists for $bs"
    fi
}

# -----------------------------------------------------------------------------
# Create unique output root
# -----------------------------------------------------------------------------
random_string=$(date +%s%N)
md5sum=$(echo -n "$random_string" | md5sum | awk '{print $1}')
OUTROOT="${md5sum}_loso_dir"
mkdir -p "$OUTROOT"

echo "LOSO simulation output: $OUTROOT"
echo "Filter prefix:          $PREFIX_FILTER"
echo "Coverages:              $COVERAGES"

# ART parameters (held constant for comparability with Simulate_rnaseq.sh)
READ_LEN=100
MEAN_FRAG=350
SD_FRAG=200
SEQ_SYS="HS20"

# Initialize sidecar metadata file with header
META="$OUTROOT/LOSO_metadata.tsv"
printf "sample_id\ttrain_ref\tsuperfamily\tcoverage\n" > "$META"

# -----------------------------------------------------------------------------
# Main loop: only over FASTAs matching the prefix (default: test_*.fasta)
# -----------------------------------------------------------------------------
fasta_count=$(find "$WD" -maxdepth 1 -name "${PREFIX_FILTER}*.fasta" | wc -l)

if [[ "$fasta_count" -eq 0 ]]; then
    echo "WARNING: no FASTAs matched '${PREFIX_FILTER}*.fasta' in $WD" >&2
    exit 1
fi
echo "Found $fasta_count test FASTAs to simulate from."

find "$WD" -maxdepth 1 -name "${PREFIX_FILTER}*.fasta" | while read -r fasta_file; do

    bs=$(basename "$fasta_file" .fasta)
    sf_tag="${bs#${PREFIX_FILTER}}"                 # held-out superfamily label
    train_fasta="$WD/train_${sf_tag}.fasta"

    if [[ ! -f "$train_fasta" ]]; then
        echo "WARNING: missing training FASTA for $fasta_file (expected $train_fasta)" >&2
        train_ref_resolved="NA"
    else
        train_ref_resolved="$(cd "$(dirname "$train_fasta")" && pwd)/$(basename "$train_fasta")"
    fi

    art_dir="$OUTROOT/${bs}_dir"
    mkdir -p "$art_dir"

    echo
    echo "=== Held-out SF: $sf_tag ==="
    echo "Source FASTA: $fasta_file"
    echo "Train FASTA:  $train_ref_resolved"

    for fcov in $COVERAGES; do

        output_bs="${bs}_${fcov}x_PE"
        echo "[SF=$sf_tag, ${fcov}x] simulating..."

        art_illumina \
            -i "$fasta_file" \
            -p \
            -l "$READ_LEN" \
            -f "$fcov" \
            -m "$MEAN_FRAG" \
            -s "$SD_FRAG" \
            -ss "$SEQ_SYS" \
            -na \
            -rs 123 \
            -sam \
            -o "$art_dir/$output_bs"

        # Absolute paths (portable, avoids realpath which is missing on some HPCs)
        forward_file="$(cd "$(dirname "$art_dir/${output_bs}1.fq")" && pwd)/${output_bs}1.fq"
        reverse_file="$(cd "$(dirname "$art_dir/${output_bs}2.fq")" && pwd)/${output_bs}2.fq"
        test_ref="$(cd "$(dirname "$fasta_file")" && pwd)/$(basename "$fasta_file")"

        # 4-column manifest compatible with existing Assemblers.sh
        # Fields: sample_id  fwd  rev  ref_for_transrate(=test_fasta=ground_truth)
        printf "%s\t%s\t%s\t%s\n" \
            "$output_bs" "$forward_file" "$reverse_file" "$test_ref" \
            > "$OUTROOT/${output_bs}_samples.txt"

        # LOSO metadata sidecar (one row per simulated dataset)
        printf "%s\t%s\t%s\t%s\n" \
            "$output_bs" "$train_ref_resolved" "$sf_tag" "$fcov" \
            >> "$META"

        # Depth statistics for this simulated dataset
        run_depth "$fasta_file" "$art_dir/${output_bs}.sam"
    done
done

# Consolidate depth files
mkdir -p "$OUTROOT/depth_stats_dir"
find "$OUTROOT" -name "*.depth.txt" | while read -r depth_file; do
    bs="${depth_file%*.depth.txt}"
    dst="$OUTROOT/depth_stats_dir/${bs##*/}.depth.txt"
    mv "$depth_file" "$dst"
done

echo
echo "LOSO simulation complete."
echo "Output root:       $OUTROOT"
echo "Per-fold manifests:  ${OUTROOT}/*_samples.txt"
echo "Metadata sidecar:    $META"
exit 0
