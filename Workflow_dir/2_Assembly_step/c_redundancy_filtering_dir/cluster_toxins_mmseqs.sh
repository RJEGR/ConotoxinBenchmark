#!/bin/bash

# conda activate mmseq2     
# conda deactivate
# MMseqs2 Clustering Script for DNA Sequences at Multiple Identity Cutoffs
# Usage: ./cluster_toxins_mmseqs.sh input.fasta

set -e  # Exit on any error

# Check if input file is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 <input.fasta>"
    echo "Example: $0 conotoxins.fasta"
    exit 1
fi

INPUT_FASTA="$1"
THREADS=24

# Check if input file exists
if [ ! -f "$INPUT_FASTA" ]; then
    echo "Error: Input file '$INPUT_FASTA' not found!"
    exit 1
fi

# Check if mmseqs is available
if ! command -v mmseqs &> /dev/null; then
    echo "Error: mmseqs command not found. Please install MMseqs2 and ensure it's in your PATH."
    exit 1
fi

# Get base name for output files
BASE_NAME=$(basename "$INPUT_FASTA" | sed 's/\.[^.]*$//')
OUTPUT_CSV="${BASE_NAME}_cluster_counts.csv"
WORK_DIR="${BASE_NAME}_mmseqs_work"

echo "Starting MMseqs2 clustering analysis for: $INPUT_FASTA"
echo "Base name: $BASE_NAME"
echo "Output CSV: $OUTPUT_CSV"
echo "Using $THREADS threads"

# Create working directory
mkdir -p "$WORK_DIR"
cd "$WORK_DIR"

# Create MMseqs2 database
echo "Creating MMseqs2 database..."
DB_NAME="${BASE_NAME}_db"

# First, let's check what type of sequences we have
echo "Analyzing sequence type in input file..."
first_seq=$(grep -A1 "^>" "../$INPUT_FASTA" | tail -1 | head -c 50)
echo "First 50 characters of first sequence: $first_seq"

# Create database with explicit nucleotide type
echo "Creating nucleotide database..."
mmseqs createdb "../$INPUT_FASTA" "$DB_NAME" --dbtype 1

# Verify database type
echo "Database information:"
if [ -f "${DB_NAME}.dbtype" ]; then
    dbtype_value=$(od -t u4 -N 4 "${DB_NAME}.dbtype" | head -1 | awk '{print $2}')
    case $dbtype_value in
        0) echo "  Database type: Amino acid" ;;
        1) echo "  Database type: Nucleotide" ;;
        2) echo "  Database type: Profile" ;;
        *) echo "  Database type: Unknown ($dbtype_value)" ;;
    esac
else
    echo "  Warning: Could not read database type file"
fi

# Initialize CSV file with header
echo "sequence_identity,num_clusters" > "../$OUTPUT_CSV"

# Function to run clustering and count clusters
run_clustering() {
    local cutoff=$1
    local cutoff_str=$(printf "%.2f" $cutoff)
    local tmp_dir="tmp_${cutoff_str}"
    local cluster_db="cluster_${cutoff_str}"
    
    echo "Processing identity cutoff: $cutoff_str"
    
    # Create temporary directory for this cutoff
    mkdir -p "$tmp_dir"
    
    # Run clustering with error output visible
    echo "  Running mmseqs cluster with --min-seq-id $cutoff..."
    if mmseqs cluster "$DB_NAME" "$cluster_db" "$tmp_dir" \
        --min-seq-id "$cutoff" \
        --threads "$THREADS" \
        --cluster-mode 0 \
        --cov-mode 0 \
        -c 0.8 \
        -e 1e-3; then
        
        # Count number of clusters
        if [ -f "${cluster_db}.index" ]; then
            num_clusters=$(wc -l < "${cluster_db}.index")
            echo "  -> Found $num_clusters clusters at $cutoff_str identity"
            
            # Append to CSV
            echo "$cutoff_str,$num_clusters" >> "../$OUTPUT_CSV"
        else
            echo "  -> ERROR: Cluster index file not created for cutoff $cutoff_str"
            echo "$cutoff_str,ERROR" >> "../$OUTPUT_CSV"
        fi
    else
        echo "  -> ERROR: Clustering failed for cutoff $cutoff_str"
        echo "$cutoff_str,ERROR" >> "../$OUTPUT_CSV"
    fi
    
    # Clean up this iteration's files
    rm -rf "$tmp_dir"
    rm -f "${cluster_db}"*
}

echo "Running clustering at different sequence identity cutoffs..."

# Define all cutoffs directly (simpler and more reliable)
declare -a cutoffs=(
    1.00 0.99 0.98 0.97 0.96 0.95 0.94 0.93 0.92 0.91 0.90
    0.80 0.70 0.60 0.50
)

echo "Will test ${#cutoffs[@]} different cutoffs: ${cutoffs[*]}"

# Run clustering for each cutoff
for cutoff in "${cutoffs[@]}"; do
    run_clustering "$cutoff"
done

# Clean up database files and working directory
echo "Cleaning up intermediate files..."
cd ..
rm -rf "$WORK_DIR"

# Display final results
echo ""
echo "=== CLUSTERING ANALYSIS COMPLETE ==="
echo "Results saved to: $OUTPUT_CSV"
echo ""
echo "Summary of results:"
echo "Sequence Identity,Number of Clusters"
cat "$OUTPUT_CSV" | tail -n +2 | head -10  # Show first 10 results
if [ $(wc -l < "$OUTPUT_CSV") -gt 11 ]; then
    echo "... (showing first 10 results)"
fi

echo ""
echo "Total sequence identity cutoffs tested: $(($(wc -l < "$OUTPUT_CSV") - 1))"
echo "Analysis complete!"

# Optional: Display basic statistics
total_seqs=$(grep -c "^>" "../$INPUT_FASTA")
echo "Total sequences in input: $total_seqs"

# Get cluster count at highest identity (most restrictive)
highest_id_clusters=$(head -2 "$OUTPUT_CSV" | tail -1 | cut -d',' -f2)
echo "Clusters at highest identity ($(head -2 "$OUTPUT_CSV" | tail -1 | cut -d',' -f1)): $highest_id_clusters"

# Get cluster count at lowest identity (least restrictive)
lowest_id_clusters=$(tail -1 "$OUTPUT_CSV" | cut -d',' -f2)
lowest_id=$(tail -1 "$OUTPUT_CSV" | cut -d',' -f1)
echo "Clusters at lowest identity ($lowest_id): $lowest_id_clusters"