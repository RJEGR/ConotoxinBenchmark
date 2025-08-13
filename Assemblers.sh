#!/bin/bash

while getopts "s:c:h" opt; do
    case $opt in
        s) MANIFEST_FILE="$OPTARG";;
        c) CONFIG_FILE="$OPTARG";;
        h)
            echo "Usage: $(basename $0) -s <Manifest> -c <Config>"
            echo
            echo "Arguments:"
            echo "  -s <Manifest>: A text file containing sample information. The first column should contain sample names, and the second column should contain file paths." Option -s all, allow to run assembly batches for all manifests matching the pattern *.txt.
            echo "  -c <Config>:   A configuration file where each line specifies a tool and its command in the format 'tool=command'."
            echo
            echo "Description:"
            echo "  This script performs meta-analysis by concatenating reads from the Manifest file, configuring tools from the Config file, and executing the specified commands."
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
CONFIG=$CONFIG_FILE 

if [[ -z "$Manifest" || -z "$CONFIG" ]]; then
    echo "Error: Missing required arguments."
    echo "Run '$(basename $0) -h' for usage information."
    exit 1
fi

chmod +x run_assemblers_dir/*sh

# set environment variables
EXPORT=$PWD/run_assemblers_dir
export PATH=$PATH:$EXPORT

run_assembly_batches() {
    local manifest="$1"
    local config="$2"

    local forward_fq="${manifest%.*}_concat_PE1.fq"
    local reverse_fq="${manifest%.*}_concat_PE2.fq"

    # Concatenate forward reads
    local call=$(awk '{print $2}' "$manifest" | tr "\n" " ")
    
    #cat $call > "$forward_fq"

    # Filter out empty entries before concatenation
    for fq in $call; do [[ -n "$fq" ]] && echo "$fq"; done | xargs cat > "$forward_fq"

    # Concatenate reverse reads
    local call=$(awk '{print $3}' "$manifest" | tr "\n" " ")
    for fq in $call; do [[ -n "$fq" ]] && echo "$fq"; done | xargs cat > "$reverse_fq"
    

    local CPU=${SLURM_CPUS_ON_NODE}
    local MEM=${SLURM_MEM_PER_NODE}

    local FASTA_DIR="${manifest%.*}_FASTA_DIR"

    mkdircall="mkdir -p "$FASTA_DIR""
    echo $mkdircall
    eval $mkdircall

    mkdir -p chkp_dir

    while read -r line; do 
        # Skip comments and empty lines
        [[ "$line" =~ ^[[:space:]]*# ]] && continue
        [[ -z "${line// }" ]] && continue

        local tool=$(echo "$line" | awk -F "=" '{print $1}')
        local bs=$(basename "${manifest%.*}")
        local call=$(echo "$line" | awk -F "=" '{print $2}') 

        local chkp_file="1_${tool}_${bs}.chkp" 

        if [ ! -f chkp_dir/"$chkp_file" ]; then
        
            local OUTDIR="${bs}_${tool}_dir"
            
            mkdir -p "$OUTDIR"

            call=$(eval echo $call "$forward_fq" "$reverse_fq" "$OUTDIR" "$CPU" "$MEM")         

            echo "Executing: $call"
            
            eval $call

            # write a chunk of code to check if any fasta file was creates in $OUTDIR, if true, touch "$chkp_file", else echo "No fasta file created in $OUTDIR, skipping checkpoint creation."
            if ls "$FASTA_DIR"/${OUTDIR%_dir}.fa 1> /dev/null 2>&1; then
                
                echo "Fasta file created in $OUTDIR, creating checkpoint."

                touch chkp_dir/"$chkp_file"
                
                rm -fr "$OUTDIR"
            else
                echo "No fasta file created in $OUTDIR, skipping checkpoint creation."
                echo "Saving outputs in issues_dir directory for further investigation."
                mkdir -p issues_dir
                mv "$OUTDIR" issues_dir
                continue
            fi
            

        else
            echo "${chkp_file} checkpoint already exists."
        fi
    done < "$config"


    # After running all assemblers, evaluate accuracy of assemblies using transrate or other tool

    # I recomend to write sub-module Accuracy.sh to call this step

     find "$FASTA_DIR" -maxdepth 1 -type f -name '*.fa' | while read -r ASSEMBLER; do
        REFDIR=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/fernando_pub/reads_artificiales/inputs
        FAMILYPREFIX=$(basename "${FASTA_DIR%_FASTA_DIR}")
        REF="$REFDIR/${FAMILYPREFIX}.fasta"
        BSCONTIG=$(basename "${ASSEMBLER%.fa}")
        TRANSRATE_DIR="${BSCONTIG}_transrate_dir"

      
        transrate_call="Run_transrate.sh $forward_fq $reverse_fq $ASSEMBLER $REF $TRANSRATE_DIR"
        echo "Executing: $transrate_call"
        eval $transrate_call

        done

    if [[ -n "${manifest%.*}" ]]; then
        rm -f "${manifest%.*}_concat_PE"*.fq
        echo "files to remove ${manifest%.*}"
        
    else
        echo "No concatenated files to remove."
    fi
}

if [[ "$Manifest" == "all" ]]; then
    echo "Running assembly batches for all manifests matching the pattern..."
    
    find "$PWD" -maxdepth 1 -type f -name '*.txt' | while read -r MANIFESTBATCH; do 
         call="run_assembly_batches "$MANIFESTBATCH" "$CONFIG""
         echo "Executing: $call"
         eval $call
    done
else
    call="run_assembly_batches "$Manifest" "$CONFIG""
    echo "Executing: $call"
    eval $call
   
fi




exit