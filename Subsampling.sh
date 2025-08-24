#!/bin/bash

# Read a manifest file (manifest.txt file containing four coulumns of sample information in the order: factor, forward, reverse, and reference)
# For every line in manifest, concatenate the forward and reverse reads, subsample them with seqkit-tool, and perform Trinity (or other) assembly.
# The assembly fasta is used to calculate the accuracy of the assembly using transrate too 

# Import seqkit to the environment

while getopts "s:o:p:h" opt; do
    case $opt in
        s) MANIFEST_FILE="$OPTARG";;
        o) OUTPUT_DIR="$OPTARG";;
        p) PROPORTION="$OPTARG";;
        h)
            echo "Usage: $(basename $0) -s <Manifest> [-o <OutputDir>]"
            echo
            echo "Arguments:"
            echo "  -s <Manifest>: A text file containing sample information. Option -s all, allow to run assembly batches for all manifests matching the pattern *.txt."
            echo "  -o <OutputDir>: Directory to save all outputs. If not provided, a default folder will be created using md5sum hash, including {Manifest} prefix and proportion value."
            echo "  -p <Proportion>: Proportion value."
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
Prop=$PROPORTION

if [[ -z "$manifest" ]]; then
    echo "Error: Missing required arguments."
    echo "Run '$(basename $0) -h' for usage information."
    exit 1
fi

run_trinity() {
    module load trinityrnaseq-v2.15.1
    
    local forward_fq="$1"
    local reverse_fq="$2"
    local output="$3"

    output_dir=${output}/Trinity_out_dir # Trinity_out_dir prefix is mandatory for Trinity

    call="Trinity --seqType fq --max_memory 100G --left $forward_fq --right $reverse_fq --no_normalize_reads --CPU 24 --output $output_dir --full_cleanup"
    
    echo $call
    
    eval $call
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

    # Here is not finding contig_file 

    contig_file=$(find "transrate_tmp_dir/$TRANSRATE_DIR"  -name "contigs.csv")

    dst="transrate_contigs_dir/${contig_file#transrate_tmp_dir}"

    mkdir -p "$(dirname "$dst")"
            
    echo "Moving $contig_file to $dst"

    mv "$contig_file" "$dst"

    rm -fr transrate_tmp_dir/$TRANSRATE_DIR
    

}

subsampling_assembly() {


    export PATH=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software:$PATH

    
    local manifest="$1"
    local Prop="$2"
    local outdir_flag="$3"

    echo "Using output directory: $outdir_flag, and $manifest, and $Prop"
    
    local bs=$(basename "${manifest%.*}")

    local OUTDIR="${outdir_flag}/${bs}_${Prop}_dir"

    mkdircall="mkdir -p "$OUTDIR""
    
    #echo $mkdircall
    eval $mkdircall

    local FASTA_DIR="${outdir_flag}/${bs}_FASTA_DIR"

    mkdircall="mkdir -p "$FASTA_DIR""
    
    #echo $mkdircall
    
    eval $mkdircall

    local forward_fq="${OUTDIR}/${bs}_concat_PE1.fq"
    
    local reverse_fq="${OUTDIR}/${bs}_concat_PE2.fq"

    # Concatenate forward reads
    call=$(awk '{print $2}' "$manifest" | tr "\n" " ")
    for fq in $call; do [[ -n "$fq" ]] && echo "$fq"; done | xargs cat > "$forward_fq"

    # Concatenate reverse reads
    call=$(awk '{print $3}' "$manifest" | tr "\n" " ")
    for fq in $call; do [[ -n "$fq" ]] && echo "$fq"; done | xargs cat > "$reverse_fq"

    local forward_sampled_fq="${OUTDIR}/${bs}_${Prop}_sampled_PE1.fq"
    local reverse_sampled_fq="${OUTDIR}/${bs}_${Prop}_sampled_PE2.fq"

    echo "Running subsampling with seqkit for proportion $Prop"

    seqkit sample --proportion $Prop --rand-seed 123 -o $forward_sampled_fq $forward_fq
    seqkit sample --proportion $Prop --rand-seed 123 -o $reverse_sampled_fq $reverse_fq

    #call="Run_trinity.sh $forward_sampled_fq $reverse_sampled_fq $OUTDIR 20 100"

    call="run_trinity $forward_sampled_fq $reverse_sampled_fq $OUTDIR"

    echo "Running $call"
    
    eval $call &> /dev/null

    echo "Finishing $call"

    f1=$(find "${OUTDIR}" -type f -name 'Trinity_out_dir.Trinity.fasta')
    if [[ ! -f "$f1" ]]; then
        f1=$(find "${OUTDIR}" -type f -name 'Trinity.tmp.fasta')
    fi

    echo "Results of the assembly found at: $f1"
   
    final_fasta=$FASTA_DIR/${bs}_${Prop}.fa

    movecall="mv $f1 $final_fasta"
   
    #echo "moving file $f1"

    echo $movecall
   
    eval $movecall


    rm -f $forward_fq
    rm -f $reverse_fq

    rm -f $forward_sampled_fq 
    rm -f $reverse_sampled_fq

    rm -rf $OUTDIR/Trinity_out_dir

    echo "..."
    echo "Finished subsampling and assembly for proportion $Prop in directory $OUTDIR"
    echo "Results saved in $FASTA_DIR/${Prop}.fa"


    echo "Step 2: calculating accuracy of the assemblies..."

    run_reference_transrate "$final_fasta" "$manifest" &> /dev/null

}


# Using make_subsamples someway $f1 recycle the previous output collapse moving the output to the next runing step
# Therefore, we need to run the subsampling_assembly function for each proportion separately


if [[ "$manifest" == "all" ]]; then
    echo "Running assembly for all manifests matching the pattern. ======================"

    OUTPUT_DIR="multiple_analysis_dir"
    
    find "$PWD" -maxdepth 1 -name '*.txt' | while read -r MANIFESTBATCH; do 
        
         echo "Executing: make_subsamples \"$MANIFESTBATCH\" \"$OUTPUT_DIR\""
         subsampling_assembly "$MANIFESTBATCH" "$Prop" "$OUTPUT_DIR"
    done
else
    OUTPUT_DIR="${manifest%.*}_dir"

    echo "Running assembly for the specified manifest ======================"
    echo "Executing: make_subsamples \"$manifest\" \"$OUTPUT_DIR\""
    subsampling_assembly "$manifest" "$Prop" "$OUTPUT_DIR"
   
fi

rm -fr transrate_tmp_dir


exit

#subsampling_assembly "$manifest" "0.1" "$out_dir"
#subsampling_assembly "$manifest" "0.2" "$out_dir"
#subsampling_assembly "$manifest" "0.3" "$out_dir"
subsampling_assembly "$manifest" "0.3" "$out_dir"
subsampling_assembly "$manifest" "0.4" "$out_dir"
subsampling_assembly "$manifest" "0.5" "$out_dir"
subsampling_assembly "$manifest" "0.6" "$out_dir"
subsampling_assembly "$manifest" "0.7" "$out_dir"
subsampling_assembly "$manifest" "0.8" "$out_dir"
subsampling_assembly "$manifest" "0.9" "$out_dir"
subsampling_assembly "$manifest" "1" "$out_dir"



make_subsamples() {
    
    local manifest="$1"
    local outdir_flag="$2"

    if [[ -z "$outdir_flag" ]]; then
        # local random_string="${manifest}_${Prop}"
        local random_string="${manifest}_$(date +%d%m%y)" 
        md5sum=$(echo -n "$random_string" | md5sum | awk '{print $1}')
        outdir_flag="${md5sum}_dir"
        # outdir="$(echo -n "$random_string" | md5sum | awk '{print $1"_dir"}')"
    fi

    for P in $(seq 0.1 0.1 1); do
        echo "Processing proportion: $P ======================"
        call="subsampling_assembly \"$manifest\" \"$P\" \"$outdir_flag\""
        echo "Executing: $call"
        eval "$call"
    done
}

 
if [[ "$manifest" == "all" ]]; then
    echo "Running assembly for all manifests matching the pattern. ======================"
    
    find "$PWD" -maxdepth 1 -name '*.txt' | while read -r MANIFESTBATCH; do 
         
         echo "Executing: analysis for $MANIFESTBATCH"

         echo "Executing: analysis for $MANIFESTBATCH"
         echo "Executing: make_subsamples \"$MANIFESTBATCH\" \"$OUTPUT_DIR\""
         make_subsamples "$MANIFESTBATCH" "$OUTPUT_DIR"
    done
else
    echo "Running assembly for the specified manifest..."
    echo "Executing: make_subsamples \"$manifest\" \"$OUTPUT_DIR\""
    make_subsamples "$manifest" "$OUTPUT_DIR"
   
fi



# evals example

All commands completed successfully. :-)

** Harvesting all assembled transcripts into a single multi-fasta file...

Saturday, August 23, 2025: 08:29:03	CMD: find /LUSTRE/bioinformatica_data/genomica_funcional/rgomez/fernando_pub/2_subsampling_dir/test_issue_dir/Fold01_200x_PE_samples_dir/Fold01_200x_PE_samples_0.2_dir/Trinity_out_dir/read_partitions/ -name '*inity.fasta'  | /LUSTRE/apps/bioinformatica/trinityrnaseq-v2.15.1/util/support_scripts/partitioned_trinity_aggregator.pl --token_prefix TRINITY_DN --output_prefix /LUSTRE/bioinformatica_data/genomica_funcional/rgomez/fernando_pub/2_subsampling_dir/test_issue_dir/Fold01_200x_PE_samples_dir/Fold01_200x_PE_samples_0.2_dir/Trinity_out_dir/Trinity.tmp
* [Sat Aug 23 08:29:04 2025] Running CMD: /LUSTRE/apps/bioinformatica/trinityrnaseq-v2.15.1/util/support_scripts/salmon_runner.pl Trinity.tmp.fasta /LUSTRE/bioinformatica_data/genomica_funcional/rgomez/fernando_pub/2_subsampling_dir/test_issue_dir/Fold01_200x_PE_samples_dir/Fold01_200x_PE_samples_0.2_dir/Trinity_out_dir/both.fa 24
* [Sat Aug 23 08:29:07 2025] Running CMD: /LUSTRE/apps/bioinformatica/trinityrnaseq-v2.15.1/util/support_scripts/filter_transcripts_require_min_cov.pl Trinity.tmp.fasta /LUSTRE/bioinformatica_data/genomica_funcional/rgomez/fernando_pub/2_subsampling_dir/test_issue_dir/Fold01_200x_PE_samples_dir/Fold01_200x_PE_samples_0.2_dir/Trinity_out_dir/both.fa salmon_outdir/quant.sf 2 > /LUSTRE/bioinformatica_data/genomica_funcional/rgomez/fernando_pub/2_subsampling_dir/test_issue_dir/Fold01_200x_PE_samples_dir/Fold01_200x_PE_samples_0.2_dir/Trinity_out_dir.Trinity.fasta


#############################################################################
Finished.  Final Trinity assemblies are written to /LUSTRE/bioinformatica_data/genomica_funcional/rgomez/fernando_pub/2_subsampling_dir/test_issue_dir/Fold01_200x_PE_samples_dir/Fold01_200x_PE_samples_0.2_dir/Trinity_out_dir.Trinity.fasta
#############################################################################
