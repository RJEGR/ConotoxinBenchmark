#!/bin/bash
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00

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


run_reference_transrate() {


    export PATH=/LUSTRE/apps/bioinformatica/.local/bin:$PATH
    export PATH=/LUSTRE/apps/bioinformatica/ruby2/ruby-2.2.0/bin:$PATH
    export LD_LIBRARY_PATH=/LUSTRE/apps/bioinformatica/ruby2/ruby-2.2.0/lib:$LD_LIBRARY_PATH

    mkdir -p transrate_tmp_dir
    mkdir -p transrate_contigs_dir
    
    local QUERY="$1"
    
    #local Manifest="$2"

    local TARGET=$(awk '{print $4}' "$2")
        
    echo "Processing FASTA : $QUERY against $TARGET"

    BS=$(basename "${QUERY%.*}")
        
    TRANSRATE_DIR=${BS}_dir

    echo "Transrate directory is : $TRANSRATE_DIR"

    call="transrate --assembly $QUERY --reference $TARGET --output transrate_tmp_dir/$TRANSRATE_DIR --threads 20"

    echo "Executing: $call"

    eval $call &> /dev/null

    echo "Finished processing $QUERY and $TARGET in transrate_tmp_dir/$TRANSRATE_DIR"

    echo "Moving results to final directory: transrate_contigs_dir"

    # Here is not finding contig_file 

    contig_file=$(find "transrate_tmp_dir/$TRANSRATE_DIR"  -name "contigs.csv")

    dst="transrate_contigs_dir/${contig_file#transrate_tmp_dir}"

    mkdir -p "$(dirname "$dst")"
            
    echo "Moving $contig_file to $dst"

    mv "$contig_file" "$dst"

    rm -fr transrate_tmp_dir/$TRANSRATE_DIR
    

}

run_assembly_batches() {
    local manifest="$1"
    local config="$2"

    local forward_fq="${manifest%.*}_concat_PE1.fq"
    local reverse_fq="${manifest%.*}_concat_PE2.fq"

    # Concatenate forward reads
    local call=$(awk '{print $2}' "$manifest" | tr "\n" " ")
    
    # Filter out empty entries before concatenation
    for fq in $call; do [[ -n "$fq" ]] && echo "$fq"; done | xargs cat > "$forward_fq"

    # Concatenate reverse reads
    local call=$(awk '{print $3}' "$manifest" | tr "\n" " ")
    for fq in $call; do [[ -n "$fq" ]] && echo "$fq"; done | xargs cat > "$reverse_fq"
    

    #local CPU=${SLURM_CPUS_ON_NODE}
    #local MEM=${SLURM_MEM_PER_NODE}

    CPU=20
    MEM=100

    local REFERENCE=$(awk '{print $4}' "$manifest")

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
        local run_tool=$(echo "$line" | awk -F "=" '{print $2}') 

        local chkp_file="1_${tool}_${bs}.chkp" 

        if [ ! -f chkp_dir/"$chkp_file" ]; then

            echo "Step 1: Running assembly with $tool and $bs"
        
            local OUTDIR="${bs}_${tool}_dir"
            
            mkdir -p "$OUTDIR"

            call=$(eval echo $run_tool "$forward_fq" "$reverse_fq" "$OUTDIR" "$CPU" "$MEM" "$FASTA_DIR" "$REFERENCE")         

            echo "Executing: $call"
            
            eval $call #&> /dev/null

            # write a chunk of code to check if any fasta file was creates in $OUTDIR, if true, touch "$chkp_file", else echo "No fasta file created in $OUTDIR, skipping checkpoint creation."
            
            final_fasta="$FASTA_DIR"/${OUTDIR%_dir}.fa
            
            if ls $final_fasta 1> /dev/null 2>&1; then
                
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

        


        echo "Step 2: calculating accuracy of the assemblies..."

        call="run_reference_transrate "$final_fasta" "$manifest"" 
        echo "Executing: $call"
        eval $call &> /dev/null

        else
            echo "${chkp_file} checkpoint already exists."
        fi
    done < "$config"


    if [[ -n "${manifest%.*}" ]]; then
        rm -f "${manifest%.*}_concat_PE"*.fq
        echo "files to remove ${manifest%.*}"
        
    else
        echo "No concatenated files to remove."
    fi
}


if [[ "$Manifest" == "all" ]]; then
    echo "Running assembly batches for all manifests matching the pattern..."
    
    find "$PWD" -maxdepth 1 -name '*.txt' | while read -r MANIFESTBATCH; do 
         call="run_assembly_batches "$MANIFESTBATCH" "$CONFIG""
         echo "Executing: $call"
         eval $call
    done

else

    echo "Step 1: Running assemblies..."

    call="run_assembly_batches "$Manifest" "$CONFIG""
    echo "Executing: $call"
    eval $call

    echo "======================================================================"
    echo "Finished assembly"
   
fi



exit