#!/bin/bash

# Note, if omicas, first load conda
# module load conda-2024

EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/art_bin_MountRainier
export PATH=$PATH:$EXPORT


while getopts "d:h" opt; do
    case $opt in
        d) FASTA_DIR="$OPTARG";;
        h)
            echo "Usage: $(basename $0) -s <Manifest> [-o <OutputDir>]"
            echo
            echo "Arguments:"
            echo "  -d <Fasta directory> Relative or absolute path to the directory containing fasta files."
            echo
            echo "Description:"
            echo "  This script simulate rnaseq data using art_illumina software."
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

WD=$FASTA_DIR

if [[ -z "$WD" ]]; then
    echo "Error: Missing required arguments."
    echo "Run '$(basename $0) -h' for usage information."
    exit 1
fi

run_transrate() {

  mkdir -p transrate_tmp_dir

  export PATH=/LUSTRE/apps/bioinformatica/.local/bin:$PATH
  export PATH=/LUSTRE/apps/bioinformatica/ruby2/ruby-2.2.0/bin:$PATH
  export LD_LIBRARY_PATH=/LUSTRE/apps/bioinformatica/ruby2/ruby-2.2.0/lib:$LD_LIBRARY_PATH

  #which transrate

    local ref=$1
    
    local Manifest=$2


    local forward_fq=${Manifest%.*}_concat_PE1.fq
    
    local reverse_fq=${Manifest%.*}_concat_PE2.fq

    BS=$(awk '{print $1}' "$Manifest")

    local call=$(awk '{print $2}' "$manifest" | tr "\n" " ")
    for fq in $call; do [[ -n "$fq" ]] && echo "$fq"; done | xargs cat > "$forward_fq"

    local call=$(awk '{print $3}' "$manifest" | tr "\n" " ")
    for fq in $call; do [[ -n "$fq" ]] && echo "$fq"; done | xargs cat > "$reverse_fq"

    TRANSRATE_DIR=${BS}_transrate_dir


    call="transrate --left $forward_fq --right $reverse_fq --assembly $ref --output transrate_tmp_dir/$TRANSRATE_DIR --threads 20"

    echo "Executing: $call"
    
    eval $call

    rm "${Manifest%.*}_concat_PE"*.fq
}

run_depth() {

# 1. index the bam file
## 2. sort the bam file
## 3. calculate the flagstat of the bam file
## 4. calculate the stats of the bam fileinput_bam=${bs}.sam

  EXPORT=/LUSTRE/apps/bioinformatica/RSEM/bin/samtools-1.3/
  export PATH=$PATH:$EXPORT

  local reference=$1
  local sam_file=$2
    
  local  bs="${sam_file%*.sam}"
  local  bam_index_file=${bs}.bam
  local  bam_file=${bs}.sort.bam

    touch_file=${bs}.depth.chkpt


    if [ ! -f "$touch_file" ]; then

    
    call="samtools view -b -S $sam_file > $bam_index_file"

    eval $call

    call="samtools sort -@ 20 $bam_index_file -o $bam_file"

    eval $call

    call="samtools flagstat $bam_file > ${bs}.flagstats.txt"

    eval $call

    #call="samtools depth $bam_file > ${bs}.depth.txt"

    call="samtools depth -a --reference $reference $bam_file > ${bs}.depth.txt"

    echo "Executing: $call"

    eval $call

    touch $touch_file

    rm $bam_index_file
    rm $bam_file
    rm $sam_file

    else
    echo "File $touch_file already exists"

    fi

}

# filepath: /Users/cigom/Documents/GitHub/ConotoxinBenchmark/run_artificialset.sh
random_string=$(date +%s%N) # Current timestamp in nanoseconds
md5sum=$(echo -n "$random_string" | md5sum | awk '{print $1}') # Generate MD5 checksum from the random string

# Create a directory using the MD5 checksum as its name
mkdir "${md5sum}_dir"

echo "Directory named '${md5sum}_dir' created successfully."


# Parámetros ART
READ_LEN=100
MEAN_FRAG=350 # the average fragment size that are simulated in base pairs. It may  match with the mean fragment size of the reference
SD_FRAG=200 # : Defines the standard deviation around the mean fragment length. 
SEQ_SYS="HS20"

 find "$WD" -name "*.fasta" | while read -r fasta_file; do

  bs=$(basename "$fasta_file" .fasta)
  
  echo "Running artificial step for $fasta_file."

  art_illumina_dir=${md5sum}_dir/${bs}_dir

  mkdir -p $art_illumina_dir

  for fcov in 10 20 50 100 200 500 700 1000; do

    
    output_bs="${bs}_${fcov}x_PE"
    
    echo "Creating artificial set $fcov"
 
     call="art_illumina \
      -i "$fasta_file" \
      -p \
      -l "$READ_LEN" \
      -f "$fcov" \
      -m "$MEAN_FRAG" \
      -s "$SD_FRAG" \
      -ss "$SEQ_SYS" \
      -na \
      -rs 123 \
      -sam \
      -o $art_illumina_dir/$output_bs"
    
    echo $call
    
    eval $call

    touch ${md5sum}_dir/${output_bs}_samples.txt

    forward_file="$PWD/${art_illumina_dir}/${output_bs}1.fq"
    reverse_file="$PWD/${art_illumina_dir}/${output_bs}2.fq"

    echo "$output_bs" `printf "$forward_file\t$reverse_file"` >> "${md5sum}_dir/${output_bs}_samples.txt"


    echo "Calculating depth for $output_bs"

    call="run_depth $fasta_file ${art_illumina_dir}/${output_bs}.sam"

   # Omit transrate as illumina art may output sam files to calculate coverage stats (sam_stats.sh)
    # call="run_transrate $fasta_file ${md5sum}_dir/${output_bs}_samples.txt"
    
    echo "Executing: $call"
    
    eval $call


  done

done




mkdir -p "${md5sum}_dir"/depth_stats_dir

find "${md5sum}_dir" -name "*.depth.txt" | while read -r depth_file; do

    bs="${depth_file%*.depth.txt}"
    dst="${md5sum}_dir/depth_stats_dir/"${bs##*/}.depth.txt

    echo "Moving $depth_file to $dst"
    #echo ${bs##*/}

    mv "$depth_file" "$dst"
done


exit
