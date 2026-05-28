#!/bin/bash
#SBATCH --job-name=conotoxin_clustering
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00
#SBATCH -o out.%J
#SBATCH -e err.%J

# Load required modules (adjust based on your cluster setup)

EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/mmseqs/bin
export PATH=$PATH:$EXPORT

# Set working directory (change this to your actual working directory)
WORK_DIR=$PWD
cd "$WORK_DIR"

# Path to the clustering script
CLUSTER_SCRIPT="./cluster_toxins_mmseqs.sh"

# Check if clustering script exists
if [ ! -f "$CLUSTER_SCRIPT" ]; then
    echo "ERROR: Clustering script not found: $CLUSTER_SCRIPT"
    echo "Please ensure cluster_toxins_mmseqs.sh is in the working directory"
    exit 1
fi

# Make sure the script is executable
chmod +x "$CLUSTER_SCRIPT"

# Create output directory for results
RESULTS_DIR="clustering_results_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$RESULTS_DIR"

echo "Starting batch clustering analysis..."
echo "Working directory: $WORK_DIR"
echo "Results directory: $RESULTS_DIR"
echo "Timestamp: $(date)"
echo "=========================================="

# Initialize counters
total_files=0
successful_files=0
failed_files=0

# Create summary log
SUMMARY_LOG="$RESULTS_DIR/batch_summary.log"
echo "Batch Clustering Analysis Summary" > "$SUMMARY_LOG"
echo "Started: $(date)" >> "$SUMMARY_LOG"
echo "Working Directory: $WORK_DIR" >> "$SUMMARY_LOG"
echo "=========================================" >> "$SUMMARY_LOG"

# Find all FASTA files and process them
for fasta_file in *.fa *.fasta; do
    # Check if glob matched any files (avoid processing literal *.fa if no files exist)
    if [ ! -f "$fasta_file" ]; then
        continue
    fi
    
    total_files=$((total_files + 1))
    echo "Processing file $total_files: $fasta_file"
    
    # Get base name for this file
    base_name=$(basename "$fasta_file" | sed 's/\.[^.]*$//')
    
    # Create individual log file for this analysis
    log_file="$RESULTS_DIR/${base_name}_clustering.log"
    
    echo "  -> Starting clustering analysis..."
    echo "  -> Log file: $log_file"
    
    # Run clustering script and capture output
    if timeout 10800 "$CLUSTER_SCRIPT" "$fasta_file" > "$log_file" 2>&1; then
        echo "  -> SUCCESS: Clustering completed"
        successful_files=$((successful_files + 1))
        
        # Move result CSV to results directory
        result_csv="${base_name}_cluster_counts.csv"
        if [ -f "$result_csv" ]; then
            mv "$result_csv" "$RESULTS_DIR/"
            echo "  -> Results saved: $RESULTS_DIR/$result_csv"
        else
            echo "  -> WARNING: Expected CSV file not found: $result_csv"
        fi
        
        # Log success
        echo "SUCCESS: $fasta_file -> $RESULTS_DIR/$result_csv" >> "$SUMMARY_LOG"
        
    else
        exit_code=$?
        echo "  -> FAILED: Clustering failed (exit code: $exit_code)"
        failed_files=$((failed_files + 1))
        
        # Log failure
        echo "FAILED: $fasta_file (exit code: $exit_code)" >> "$SUMMARY_LOG"
        
        # Check if it was a timeout
        if [ $exit_code -eq 124 ]; then
            echo "  -> Reason: Timeout (>3 hours)" | tee -a "$SUMMARY_LOG"
        fi
    fi
    
    echo "  -> Completed processing $fasta_file"
    echo ""
done

# Final summary
echo "=========================================="
echo "Batch clustering analysis completed!"
echo "Total files processed: $total_files"
echo "Successful: $successful_files"
echo "Failed: $failed_files"
echo "Results directory: $RESULTS_DIR"
echo "Completed: $(date)"

# Update summary log
echo "=========================================" >> "$SUMMARY_LOG"
echo "Completed: $(date)" >> "$SUMMARY_LOG"
echo "Total files processed: $total_files" >> "$SUMMARY_LOG"
echo "Successful: $successful_files" >> "$SUMMARY_LOG"
echo "Failed: $failed_files" >> "$SUMMARY_LOG"

# If no files were found, provide helpful message
if [ $total_files -eq 0 ]; then
    echo "WARNING: No FASTA files (*.fa, *.fasta) found in current directory!"
    echo "Please check that:"
    echo "1. You're in the correct working directory"
    echo "2. Your FASTA files have .fa or .fasta extensions"
    echo "3. The files are accessible"
    exit 1
fi

# Create consolidated results file
if [ $successful_files -gt 0 ]; then
    echo "Creating consolidated results..."
    CONSOLIDATED_RESULTS="$RESULTS_DIR/all_clustering_results.csv"
    
    # Header
    echo "dataset,sequence_identity,num_clusters" > "$CONSOLIDATED_RESULTS"
    
    # Combine all individual CSV files
    for csv_file in "$RESULTS_DIR"/*_cluster_counts.csv; do
        if [ -f "$csv_file" ]; then
            dataset=$(basename "$csv_file" | sed 's/_cluster_counts.csv$//')
            # Skip header and add dataset name
            tail -n +2 "$csv_file" | sed "s/^/$dataset,/" >> "$CONSOLIDATED_RESULTS"
        fi
    done
    
    echo "Consolidated results saved: $CONSOLIDATED_RESULTS"
    echo "Consolidated results: $CONSOLIDATED_RESULTS" >> "$SUMMARY_LOG"
fi

echo "Analysis complete. Check $SUMMARY_LOG for detailed summary."
