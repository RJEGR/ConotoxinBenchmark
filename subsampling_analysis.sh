#!/bin/bash

#SBATCH --job-name=subsampling_analysis
#SBATCH -N 1
#SBATCH --mem=120GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00
#SBATCH -o out.%J
#SBATCH -e err.%J

export PATH=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software:$PATH


EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/fernando_pub/1_assembly_dir/RUN_ARTIFICIAL_PIPELINE_DIR/run_assemblers_dir
export PATH=$PATH:$EXPORT

while getopts "s:h" opt; do
    case $opt in
        s) MANIFEST_FILE="$OPTARG";;
        #p) SAMPLE_FRAC="$OPTARG";;
        h)
            echo "Usage: $(basename $0) -s <Manifest> -p <sample fraction [0.1 - 1]>"
            echo
            echo "Arguments:"
            echo "  -s <Manifest>: A text file containing sample information." Option -s all, allow to run assembly batches for all manifests matching the pattern *.txt.
            #echo "  -p <sample fraction>:  fraction of sample."
            echo
            echo "Description:"
            echo "  This script performs subsampling of reads from a given Manifest file. by concatenating reads from the Manifest file, and executing the specified commands."
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

Manifest=$MANIFEST_FILE 
#Prop=$SAMPLE_FRAC 

#if [[ -z "$Manifest" || -z "$SAMPLE_FRAC" ]]; then
if [[ -z "$Manifest" ]]; then
    echo "Error: Missing required arguments."
    echo "Run '$(basename $0) -h' for usage information."
    exit 1
fi

run_assembly() {

    #local CPU=${SLURM_CPUS_ON_NODE}
    #local MEM=${SLURM_MEM_PER_NODE}

    local manifest="$1"
    local Prop="$2"

    local bs=$(basename "${manifest%.*}")

    local OUTDIR="${bs}_${Prop}_dir"
            
    mkdir -p "$OUTDIR"

    local FASTA_DIR="${manifest%.*}_FASTA_DIR"

    mkdircall="mkdir -p "$FASTA_DIR""
    
    echo $mkdircall
    eval $mkdircall

    local forward_fq="${manifest%.*}_concat_PE1.fq"
    local reverse_fq="${manifest%.*}_concat_PE2.fq"

    # Concatenate forward reads
    local call=$(awk '{print $2}' "$manifest" | tr "\n" " ")
    

    # Filter out empty entries before concatenation
    for fq in $call; do [[ -n "$fq" ]] && echo "$fq"; done | xargs cat > "$forward_fq"

    # Concatenate reverse reads
    call=$(awk '{print $3}' "$manifest" | tr "\n" " ")
    for fq in $call; do [[ -n "$fq" ]] && echo "$fq"; done | xargs cat > "$reverse_fq"

    local forward_sampled_fq="${forward_fq%.fq}_${Prop}_sampled.fq"
    local reverse_sampled_fq="${reverse_fq%.fq}_${Prop}_sampled.fq"



    # Using the same seed ensures SeqKit selects the same corresponding records based on sequence order and IDs. 
    # This keeps paired-end reads in sync. (it working!!)
   
   
    seqkit sample --proportion $Prop --rand-seed 123 -o $forward_sampled_fq $forward_fq
    seqkit sample --proportion $Prop --rand-seed 123 -o $reverse_sampled_fq $reverse_fq

    call="Run_trinity.sh $forward_sampled_fq $reverse_sampled_fq $OUTDIR 20 100"

    echo $call
    
    eval $call

    find "$FASTA_DIR" -maxdepth 1 -type f -name '*.fa' | while read -r ASSEMBLER; do
        REFDIR=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/fernando_pub/reads_artificiales/inputs
        FAMILYPREFIX=$(basename "${FASTA_DIR%_FASTA_DIR}")
        REF="$REFDIR/${FAMILYPREFIX}.fasta"
        BSCONTIG=$(basename "${ASSEMBLER%.fa}")
        TRANSRATE_DIR="${BSCONTIG}_transrate_dir"

      
        transrate_call="Run_transrate.sh $forward_sampled_fq $reverse_sampled_fq $ASSEMBLER $REF $TRANSRATE_DIR"
        echo "Executing: $transrate_call"
        eval $transrate_call

        done

    if [[ -n "${manifest%.*}" ]]; then
        rm -f "${manifest%.*}_concat_PE"*.fq
        echo "files to remove ${manifest%.*}"
        
    else
        echo "No concatenated files to remove."
    fi

    rm -fr "$OUTDIR"

}

make_subsamples() {
    local manifest="$1"
    for P in $(seq 0.1 0.1 1.0); do
        call="run_assembly \"$manifest\" \"$P\""
        #echo $call
        eval $call
    done
}

# run_seqkit $Manifest $Prop

if [[ "$Manifest" == "all" ]]; then
    echo "Running assembly for all manifests matching the pattern..."
    
    find "$PWD" -maxdepth 1 -type f -name '*.txt' | while read -r MANIFESTBATCH; do 
         call="make_subsamples "$MANIFESTBATCH""
         echo "Executing: $call"
         eval $call
    done
else
    echo "Running assembly for the specified manifest..."
    make_subsamples "$Manifest"
   
fi


# run trinity (or other tool)

