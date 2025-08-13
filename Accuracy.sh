#!/bin/sh

# write help arguments to include REFDIR=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/fernando_pub/reads_artificiales/inputs


while getopts "s:r:d:h" opt; do
    case $opt in
        s) MANIFEST_FILE="$OPTARG";;
        r) REFDIR="$OPTARG";;
        d) FASTA_DIR="$OPTARG";;
        h)
            echo "Usage: $(basename $0) -s <Manifest> [-r <Reference directory>] [-d <sequences Directory>]"
            echo
            echo "Arguments:"
            echo "  -d <Fasta directory> Relative or absolute path to the directory containing fasta files."
            echo "  -s <Manifest>: A text file containing sample information. Option -s all, allow to run assembly batches for all manifests matching the pattern *.txt."
            echo "  -r <Reference>: Relative or absolute path to the directory containing fasta files."

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


if [[ -z "$FASTA_DIR" || -z "$MANIFEST_FILE" || -z "$REFDIR" ]]; then
    echo "Error: Missing required arguments."
    echo "Run '$(basename $0) -h' for usage information."
    exit 1
fi


export PATH=/LUSTRE/apps/bioinformatica/.local/bin:$PATH
export PATH=/LUSTRE/apps/bioinformatica/ruby2/ruby-2.2.0/bin:$PATH
export LD_LIBRARY_PATH=/LUSTRE/apps/bioinformatica/ruby2/ruby-2.2.0/lib:$LD_LIBRARY_PATH

which transrate

run_transrate() {
    
    REFDIR=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/fernando_pub/reads_artificiales/inputs

    mkdir -p transrate_tmp_dir

    local FASTA_DIR="$1"
    
    find $FASTA_DIR -maxdepth 1 -type f -name "*.fa" | while read -r asm; do
        
        echo "Processing FASTA : $asm"

        # Substring to extract the sf prefix

        FAMILYPREFIX=$(echo $asm | awk -F'/' '{for(i=1;i<=NF;i++) if($i ~ /_FASTA_DIR/) print $i}')


        echo ${FAMILYPREFIX%_FASTA_DIR}

        ref=$REFDIR/${FAMILYPREFIX%_FASTA_DIR}.fasta
        
        
        TRANSRATE_DIR=${asm%.fa}_dir

        call="transrate --assembly $asm --reference $ref --output transrate_tmp_dir/$TRANSRATE_DIR --threads 20"

        echo "Executing: $call"

        #eval $call

        echo "Finished processing FASTA : $asm"
    
    done



}

run_transrate "$FASTA_DIR"

exit

mkdir -p transrate_contigs_dir

find transrate_tmp_dir -name 'contigs.csv' -exec bash -c '
        for src; do
              dst="transrate_contigs_dir/${src#transrate_tmp_dir}"
              echo $dst
              mkdir -p "$(dirname "$dst")"
              cp "$src" "$dst"

        done
' bash {} +

# rm -rf transrate_tmp_dir

exit