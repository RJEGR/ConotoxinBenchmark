#!/bin/bash

while getopts "s:o:k:h" opt; do
    case $opt in
        s) MANIFEST_FILE="$OPTARG";;
        #o) OUTPUT_DIR="$OPTARG";;
        k) KMERSIZE="$OPTARG";;
        h)
            echo "Usage: $(basename $0) -s <Manifest> [-o <OutputDir>]"
            echo
            echo "Arguments:"
            echo "  -s <Manifest>: A text file containing sample information. Option -s all, allow to run assembly batches for all manifests matching the pattern *.txt."
            #echo "  -o <OutputDir>: Directory to save all outputs. If not provided, a default folder will be created using md5sum hash, including {Manifest} prefix and proportion value."
            echo "  -k <Kmer size>:  k-mer sizes (must be odd and less than 128)."
            echo
            echo "Description:"
            echo "  This script performs subsampling with seqkit-tool a given proportion between 0 to 1 of PE-end reads from the input {Manifest} file provided, and executes an assembly of these reads."
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
manifest=$MANIFEST_FILE
kmer_size=$KMERSIZE

if [[ -z "$manifest" ]]; then
    echo "Error: Missing required arguments."
    echo "Run '$(basename $0) -h' for usage information."
    exit 1
fi



run_spades() {
    EXPORT=/LUSTRE/apps/bioinformatica/SPAdes-3.15.5-Linux/bin/
    export PATH=$PATH:$EXPORT

    local forward_fq="$1"
    local reverse_fq="$2"
    local output="$3"
    local final_fasta="$4"
    local kmer_sizes="$5"

   # -k <int> [<int> ...]        list of k-mer sizes (must be odd and less than 128) [default: 'auto']

    call="rnaspades.py -k $kmer_sizes -1 $forward_fq -2 $reverse_fq  -o $output -t 24 -m 100"

    echo $call
    
    eval $call

    f1=$(find "${output}" -maxdepth 1 -type f -name 'transcripts.fasta') 

    echo "Results of the assembly found at: $f1"


    movecall="mv $f1 $final_fasta"

    echo $movecall

    eval $movecall

}

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

    contig_file=$(find "transrate_tmp_dir/$TRANSRATE_DIR"  -name "contigs.csv")

    dst="transrate_contigs_dir/${contig_file#transrate_tmp_dir}"

    mkdir -p "$(dirname "$dst")"
            
    echo "Moving $contig_file to $dst"

    mv "$contig_file" "$dst"

    rm -fr transrate_tmp_dir/$TRANSRATE_DIR
    

}

kmer_assembly() {


    export PATH=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software:$PATH

    
    local manifest="$1"
    local kmer="$2"
    local outdir_flag="$3"

    echo "Using output directory: $outdir_flag, and $manifest, and $kmer"
    
    local bs=$(basename "${manifest%.*}")

    local OUTDIR="${outdir_flag}/${bs}_${kmer}_dir"

    mkdircall="mkdir -p "$OUTDIR""
    
    #echo $mkdircall
    eval $mkdircall

    local FASTA_DIR="${outdir_flag}/${bs}_FASTA_DIR"

    mkdircall="mkdir -p "$FASTA_DIR""
    
    #echo $mkdircall
    
    eval $mkdircall

    final_fasta=$FASTA_DIR/${bs}_${kmer}.fa

    local forward_fq="${OUTDIR}/${bs}_concat_PE1.fq"
    
    local reverse_fq="${OUTDIR}/${bs}_concat_PE2.fq"

    # Concatenate forward reads
    call=$(awk '{print $2}' "$manifest" | tr "\n" " ")
    for fq in $call; do [[ -n "$fq" ]] && echo "$fq"; done | xargs cat > "$forward_fq"

    # Concatenate reverse reads
    call=$(awk '{print $3}' "$manifest" | tr "\n" " ")
    for fq in $call; do [[ -n "$fq" ]] && echo "$fq"; done | xargs cat > "$reverse_fq"

    echo "Running assembly with kmer of $kmer"


    call="run_spades $forward_fq $reverse_fq $OUTDIR $final_fasta $kmer"


    echo "Running $call"
    
    eval $call &> /dev/null

    echo "Finishing $call"


    rm -f $forward_fq
    rm -f $reverse_fq

    rm -f $forward_sampled_fq 
    rm -f $reverse_sampled_fq

    rm -rf $OUTDIR

    echo "..."
    echo "Finished subsampling and assembly for kmer $kmer in directory $OUTDIR"
    echo "Results saved in $FASTA_DIR/${kmer}.fa"


    echo "Step 2: calculating accuracy of the assemblies..."

    run_reference_transrate "$final_fasta" "$manifest" &> /dev/null

}



if [[ "$manifest" == "all" ]]; then
    echo "Running assembly for all manifests matching the pattern. ======================"

    OUTPUT_DIR="multiple_analysis_dir"
    
    find "$PWD" -maxdepth 1 -name '*.txt' | while read -r MANIFESTBATCH; do 
        
         echo "Executing: make_subsamples \"$MANIFESTBATCH\" \"$OUTPUT_DIR\""
         kmer_assembly "$MANIFESTBATCH" "$kmer_size" "$OUTPUT_DIR"
    done
else
    OUTPUT_DIR="${manifest%.*}_dir"

    echo "Running assembly for the specified manifest ======================"
    echo "Executing: make_subsamples \"$manifest\" \"$OUTPUT_DIR\""
    kmer_assembly "$manifest" "$kmer_size" "$OUTPUT_DIR"
   
fi

rm -fr transrate_tmp_dir


exit


