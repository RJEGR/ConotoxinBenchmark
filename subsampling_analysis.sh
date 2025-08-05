#!/bin/bash

export PATH=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software:$PATH


EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/fernando_pub/1_assembly_dir/RUN_ARTIFICIAL_PIPELINE_DIR/run_assemblers_dir
export PATH=$PATH:$EXPORT

# Add -o flag for output directory, with default as md5sum of manifest+Prop

while getopts "s:o:h" opt; do
    case $opt in
        s) MANIFEST_FILE="$OPTARG";;
        o) OUTPUT_DIR="$OPTARG";;
        h)
            echo "Usage: $(basename $0) -s <Manifest> [-o <OutputDir>]"
            echo
            echo "Arguments:"
            echo "  -s <Manifest>: A text file containing sample information. Option -s all, allow to run assembly batches for all manifests matching the pattern *.txt."
            echo "  -o <OutputDir>: Directory to save all outputs. If not provided, a default folder will be created using md5sum of manifest and proportion."
            echo
            echo "Description:"
            echo "  This script performs subsampling of reads from a given Manifest file by concatenating reads from the Manifest file, and executing the specified commands."
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

if [[ -z "$Manifest" ]]; then
    echo "Error: Missing required arguments."
    echo "Run '$(basename $0) -h' for usage information."
    exit 1
fi

run_trinity() {
    module load trinityrnaseq-v2.15.1
    local forward_fq="$1"
    local reverse_fq="$2"
    local output="$3"

    call="Trinity --seqType fq --max_memory 100G --left $forward_fq --right $reverse_fq --no_normalize_reads --CPU 24 --output ${output}/Trinity_out_dir --full_cleanup"
    
    echo $call
    
    eval $call


  # rm ${output}/Trinity_out_dir/Trinity_out_dir.Trinity.fasta.gene_trans_map*.fasta.gene_trans_map
}

run_assembly() {
    local manifest="$1"
    local Prop="$2"
    local outdir_flag="$3"

    # If output dir not provided, use md5sum of manifest+Prop
    if [[ -z "$outdir_flag" ]]; then
        # local random_string="${manifest}_${Prop}"
        local random_string="${manifest}_$(date +%d%m%y)" 
        md5sum=$(echo -n "$random_string" | md5sum | awk '{print $1}')
        outdir_flag="${md5sum}_dir"
        # outdir="$(echo -n "$random_string" | md5sum | awk '{print $1"_dir"}')"
    fi

    echo "Using output directory: $outdir_flag, and $manifest, and $Prop"
    
    local bs=$(basename "${manifest%.*}")

    local OUTDIR="${outdir_flag}/${bs}_${Prop}_dir"
            
    mkdir -p "$OUTDIR"

    local FASTA_DIR="${outdir_flag}/${manifest%.*}_FASTA_DIR"

    mkdircall="mkdir -p "$FASTA_DIR""
    
    echo $mkdircall
    eval $mkdircall

    local forward_fq="${OUTDIR}/${bs}_concat_PE1.fq"
    local reverse_fq="${OUTDIR}/${bs}_concat_PE2.fq"

    # Concatenate forward reads
    local call=$(awk '{print $2}' "$manifest" | tr "\n" " ")
    for fq in $call; do [[ -n "$fq" ]] && echo "$fq"; done | xargs cat > "$forward_fq"

    # Concatenate reverse reads
    call=$(awk '{print $3}' "$manifest" | tr "\n" " ")
    for fq in $call; do [[ -n "$fq" ]] && echo "$fq"; done | xargs cat > "$reverse_fq"

    local forward_sampled_fq="${OUTDIR}/${bs}_${Prop}_sampled_PE1.fq"
    local reverse_sampled_fq="${OUTDIR}/${bs}_${Prop}_sampled_PE2.fq"

    seqkit sample --proportion $Prop --rand-seed 123 -o $forward_sampled_fq $forward_fq
    seqkit sample --proportion $Prop --rand-seed 123 -o $reverse_sampled_fq $reverse_fq

    #call="Run_trinity.sh $forward_sampled_fq $reverse_sampled_fq $OUTDIR 20 100"

    call="run_trinity $forward_sampled_fq $reverse_sampled_fq $OUTDIR"
    
    eval $call

    f1=$(find ${OUTDIR} -type f -name '*fasta')

    #f2=$(find ${outdir_flag} -type d -name '*FASTA_DIR')

    BS=$(basename "${FASTA_DIR%_FASTA_DIR}")
   
    movecall="mv $f1 $FASTA_DIR/${Prop}.fa"
   
    echo "moving file $f1"

    echo $movecall
   
    eval $movecall


    if [[ -n "${OUTDIR}" ]]; then
        rm -f "${OUTDIR}/${bs}_concat_PE"*.fq
        echo "files to remove ${OUTDIR}/${bs}_concat_PE*.fq"
    else
        echo "No concatenated files to remove."
    fi

}

make_subsamples() {
    local manifest="$1"
    for P in $(seq 0.1 0.1 1); do
        run_assembly "$manifest" "$P" "$OUTPUT_DIR"
    done
}

 
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

# Proccess the FASTA files with transrate


run_transrate() {
    
    REFDIR=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/fernando_pub/reads_artificiales/inputs

    mkdir -p transrate_tmp_dir
    
    find  find */*_FASTA_DIR -maxdepth 1 -type f -name "*.fa" | while read -r asm; do
        
        echo "Processing FASTA : $asm"

        # Substring to extract the sf prefix

        FAMILYPREFIX=$(echo $asm | awk -F'/' '{for(i=1;i<=NF;i++) if($i ~ /_FASTA_DIR/) print $i}')

        FAMILYPREFIX="all_superfamily"

        echo ${FAMILYPREFIX%_FASTA_DIR}

        ref=$REFDIR/${FAMILYPREFIX%_FASTA_DIR}.fasta
        
        #TRANSRATE_DIR=$(dirname "$asm")
        
        TRANSRATE_DIR=${asm%.fa}_dir

        call="transrate --assembly $asm --reference $ref --output transrate_tmp_dir/$TRANSRATE_DIR --threads 20"

        echo "Executing: $call"

        eval $call

        echo "Finished processing FASTA : $asm"
    
    done



}

run_transrate


mkdir -p transrate_contigs_dir

find transrate_tmp_dir -name 'contigs.csv' -exec bash -c '
        for src; do
              dst="transrate_contigs_dir/${src#transrate_tmp_dir}"
              echo $dst
              mkdir -p "$(dirname "$dst")"
              cp "$src" "$dst"

        done
' bash {} +

exit


