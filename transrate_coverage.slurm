#!/bin/bash
#SBATCH --job-name=transrate_coverage
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00
#SBATCH -o out.%J
#SBATCH -e err.%J


REFDIR=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/fernando_pub/reads_artificiales/inputs


while getopts "s:h" opt; do
    case $opt in
        s) MANIFEST_FILE="$OPTARG";;
        h)
            echo "Usage: $(basename $0) -s <Manifest>"
            echo
            echo "Arguments:"
            echo "  -s <Manifest>: A text file containing sample information. The first column should contain sample names, and the second column should contain file paths."
            echo
            echo "Description:"
            echo "  This script performs transrate coverage analysis based on the provided Manifest file."
            exit 0
            ;;
        ?)
            echo "Invalid option: -$OPTARG" >&2
            echo "Run '$(basename $0) -h' for usage information."
            exit 1
            ;;
    esac
done
shift $((OPTIND-1))

if [[ -z "$MANIFEST_FILE" ]]; then
    echo "Error: Missing required argument. -s (Manifest) is required."
    echo "Run '$(basename $0) -h' for usage information."
    exit 1
fi

run_transrate() {
    Manifest=$MANIFEST_FILE

    local FASTA_DIR="${Manifest%.*}_FASTA_DIR"
    
    local forward_fq=${Manifest%.*}_concat_PE1.fq
    local reverse_fq=${Manifest%.*}_concat_PE2.fq

    call=$(awk '{print $2}' "$Manifest" | tr "\n" " ")
    cat $call > $forward_fq

    call=$(awk '{print $3}' "$Manifest" | tr "\n" " ")
    cat $call > $reverse_fq

    mkdir -p 2_transrate_dir

    find "$FASTA_DIR" -maxdepth 1 -type f -name '*.fa' | while read -r asm; do
        FAMILYPREFIX=$(basename "${FASTA_DIR%_FASTA_DIR}")
        ref=$REFDIR/${FAMILYPREFIX}.fasta
        BSCONTIG=$(basename "${asm%.fa}")
        TRANSRATE_DIR=${BSCONTIG}_transrate_dir

        call="transrate --left $forward_fq --right $reverse_fq --assembly $asm --reference $ref --output 2_transrate_dir/$TRANSRATE_DIR --threads 20"

        echo "Executing: $call"
        eval $call
    done

    rm "${Manifest%.*}_concat_PE"*.fq
}

export PATH=/LUSTRE/apps/bioinformatica/.local/bin:$PATH
export PATH=/LUSTRE/apps/bioinformatica/ruby2/ruby-2.2.0/bin:$PATH
export LD_LIBRARY_PATH=/LUSTRE/apps/bioinformatica/ruby2/ruby-2.2.0/lib:$LD_LIBRARY_PATH

which transrate

# Loop over every *.txt file and run run_transrate
for Manifest in *.txt; do
    run_transrate "$Manifest"
done

if [[ "$Manifest" == "all" ]]; then
    echo "Running transrate on all batches for all manifests matching the pattern..."
    
    find "$PWD" -maxdepth 1 -type f -name '*.txt' | while read -r MANIFESTBATCH; do 
         call="run_transrate "$MANIFESTBATCH""
         echo "Executing: $call"
         eval $call
    done
else
    call="run_transrate "$Manifest""
    echo "Executing: $call"
    eval $call
   
fi

exit

