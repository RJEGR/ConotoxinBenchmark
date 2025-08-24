#!/bin/sh

# 

while getopts "s:d:m:h" opt; do
    case $opt in
        s) MANIFEST_FILE="$OPTARG";;
        m) MODE="$OPTARG";;
        d) QUERY="$OPTARG";;
        h)
            echo "Usage: $(basename $0) -s <Manifest> [-d <query sequences>]"
            echo
            echo "Arguments:"
            echo "  -s <Manifest>: A text file containing sample information. Option < -s all >, allow to run assembly batches for all manifests matching the pattern *.txt."
            echo "  -d <Query sequences> database and query are each either a fasta format"
            echo "  -m <Mode>: Assess reference-based < -m reference > or read-based < -m reads > metrics. If not provided, defaults to < -m both >."

            echo
            echo "Description:"
            echo "  This script evaluate the full-length transcript reconstruction coverage of a given DNA sequences (example contigs from assembly) based on the reads and reference"
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


if [[ -z "$QUERY" || -z "$MANIFEST_FILE" ]] ; then
    echo "Error: You must provide both QUERY and MANIFEST_FILE."
    echo "Run '$(basename $0) -h' for usage information."
    exit 1
fi


run_longer_transrate() {

    export PATH=/LUSTRE/apps/bioinformatica/.local/bin:$PATH
    export PATH=/LUSTRE/apps/bioinformatica/ruby2/ruby-2.2.0/bin:$PATH
    export LD_LIBRARY_PATH=/LUSTRE/apps/bioinformatica/ruby2/ruby-2.2.0/lib:$LD_LIBRARY_PATH

    #which transrate


    mkdir -p transrate_tmp_dir
    mkdir -p transrate_contigs_dir
    
    #local TARGET="$1"

    local QUERY="$1"

    local Manifest="$2"


    local forward_fq=${Manifest%.*}_concat_PE1.fq
    
    local reverse_fq=${Manifest%.*}_concat_PE2.fq

    #BS=$(awk '{print $1}' "$Manifest")

    #BS=${Manifest%.*}

    BS=${Manifest%.*}_${QUERY%.*}

    local call=$(awk '{print $2}' "$Manifest" | tr "\n" " ")
    for fq in $call; do [[ -n "$fq" ]] && echo "$fq"; done | xargs cat > "$forward_fq"

    local call=$(awk '{print $3}' "$Manifest" | tr "\n" " ")
    for fq in $call; do [[ -n "$fq" ]] && echo "$fq"; done | xargs cat > "$reverse_fq"
    
    local TARGET=$(awk '{print $4}' "$Manifest")

    TRANSRATE_DIR=${BS}_transrate_dir

    call="transrate --left $forward_fq --right $reverse_fq --assembly $QUERY --reference $TARGET --output transrate_tmp_dir/$TRANSRATE_DIR --threads 20"

    echo "Executing: $call"

    eval $call

    rm "$forward_fq"
    rm "$reverse_fq"

    echo "Finished processing FASTA in transrate_tmp_dir/$TRANSRATE_DIR"

    echo "Moving results to final directory: transrate_contigs_dir"

    find "transrate_tmp_dir/$TRANSRATE_DIR" -name "contigs.csv" | while read -r contig_file; do

            dst="transrate_contigs_dir/${contig_file#transrate_tmp_dir}"
            
            echo $dst
            
            mkdir -p "$(dirname "$dst")"
            
            echo "Moving $contig_file to $dst"

            mv "$contig_file" "$dst"
        done

}

run_reference_transrate() {


    export PATH=/LUSTRE/apps/bioinformatica/.local/bin:$PATH
    export PATH=/LUSTRE/apps/bioinformatica/ruby2/ruby-2.2.0/bin:$PATH
    export LD_LIBRARY_PATH=/LUSTRE/apps/bioinformatica/ruby2/ruby-2.2.0/lib:$LD_LIBRARY_PATH

    #which transrate


    mkdir -p transrate_tmp_dir
    mkdir -p transrate_contigs_dir
    
    local QUERY="$1"
    
    #local Manifest="$2"

    local TARGET=$(awk '{print $4}' "$2")


        
    echo "Processing FASTA : $QUERY"

        # Substring to extract the sf prefix

        # FAMILYPREFIX=$(echo $asm | awk -F'/' '{for(i=1;i<=NF;i++) if($i ~ /_QUERY/) print $i}')


        #echo ${FAMILYPREFIX%_QUERY}

        #ref=$TARGET/${FAMILYPREFIX%_QUERY}.fasta
        
        
    TRANSRATE_DIR=${QUERY%.*}_dir

    echo "Transrate directory is : $TRANSRATE_DIR"

    call="transrate --assembly $QUERY --reference $TARGET --output transrate_tmp_dir/$TRANSRATE_DIR --threads 20"

    echo "Executing: $call"

    eval $call

    echo "Finished processing FASTA in transrate_tmp_dir/$TRANSRATE_DIR"

    echo "Moving results to final directory: transrate_contigs_dir"

    find "transrate_tmp_dir/$TRANSRATE_DIR" -name "contigs.csv" | while read -r contig_file; do

            dst="transrate_contigs_dir/${contig_file#transrate_tmp_dir}"
            
            echo $dst
            
            mkdir -p "$(dirname "$dst")"
            
            echo "Moving $contig_file to $dst"

            mv "$contig_file" "$dst"
        done
    

}

run_read_transrate() {

    export PATH=/LUSTRE/apps/bioinformatica/.local/bin:$PATH
    export PATH=/LUSTRE/apps/bioinformatica/ruby2/ruby-2.2.0/bin:$PATH
    export LD_LIBRARY_PATH=/LUSTRE/apps/bioinformatica/ruby2/ruby-2.2.0/lib:$LD_LIBRARY_PATH

    #which transrate


    mkdir -p transrate_tmp_dir
    mkdir -p transrate_contigs_dir
    

    local QUERY="$1"
    local Manifest="$2"


    local forward_fq=${Manifest%.*}_concat_PE1.fq
    
    local reverse_fq=${Manifest%.*}_concat_PE2.fq

    BS=${Manifest%.*}_${QUERY%.*}

    local call=$(awk '{print $2}' "$Manifest" | tr "\n" " ")
    for fq in $call; do [[ -n "$fq" ]] && echo "$fq"; done | xargs cat > "$forward_fq"

    local call=$(awk '{print $3}' "$Manifest" | tr "\n" " ")
    for fq in $call; do [[ -n "$fq" ]] && echo "$fq"; done | xargs cat > "$reverse_fq"

    TRANSRATE_DIR=${BS}_transrate_dir
    
    call="transrate --left $forward_fq --right $reverse_fq --assembly $QUERY  --output transrate_tmp_dir/$TRANSRATE_DIR --threads 20"

    echo "Executing: $call"

    eval $call

    rm "$forward_fq"
    rm "$reverse_fq"

    echo "Finished processing FASTA in transrate_tmp_dir/$TRANSRATE_DIR"

    echo "Moving results to final directory: transrate_contigs_dir"

    find "transrate_tmp_dir/$TRANSRATE_DIR" -name "contigs.csv" | while read -r contig_file; do

            dst="transrate_contigs_dir/${contig_file#transrate_tmp_dir}"
            
            echo $dst
            
            mkdir -p "$(dirname "$dst")"
            
            echo "Moving $contig_file to $dst"

            mv "$contig_file" "$dst"
        done

}


# if $MODE argument is provided, run the corresponding function
if [[ "$MODE" == "reference" ]]; then
    echo "Running reference-based transrate..."
    run_reference_transrate "$QUERY" "$MANIFEST_FILE"
elif [[ "$MODE" == "reads" ]]; then
    echo "Running read-based transrate..."
    run_read_transrate "$QUERY" "$MANIFEST_FILE"
elif [[ "$MODE" == "both" ]]; then
    echo "Running both reference and read-based transrate..."
    run_longer_transrate "$QUERY" "$MANIFEST_FILE"
else
    echo "No mode specified or invalid mode. Running default transrate..."
        run_longer_transrate "$QUERY" "$MANIFEST_FILE"
fi



rm -rf transrate_tmp_dir

exit
