#!/bin/bash

# Note, if omicas, first load conda
# module load conda-2024


# filepath: /Users/cigom/Documents/GitHub/ConotoxinBenchmark/run_artificialset.sh
random_string=$(date +%s%N) # Current timestamp in nanoseconds
md5sum=$(echo -n "$random_string" | md5sum | awk '{print $1}') # Generate MD5 checksum from the random string

# Create a directory using the MD5 checksum as its name
mkdir "${md5sum}_dir"

echo "Directory named '${md5sum}_dir' created successfully."

# WD=/Users/cigom/Documents/GitHub/ConotoxinBenchmark/INPUTS

WD=$1

# Find all FASTA files in the current directory and process them using a while loop
find "$WD" -name "*.fasta" | while read -r fasta_file; do
  # Extract the base name of the file (without path and extension)
  bs=$(basename "$fasta_file" .fasta)
  
  echo "Running Artificialset.py for $fasta_file."

  # Run Artificialset.py 5 times for each FASTA file
  for i in $(seq 1 5); do
    output_fasta="${bs}_output_set${i}.fq"
    
    echo "Creating artificial set $i"
    
    call="python Artificialset.py --insert_size 100 --phred_score 40 --output_format fastq "$fasta_file" "${md5sum}_dir/$output_fasta" 1000000"
    
    eval $call

    touch ${md5sum}_dir/${bs}_samples.txt

    echo "$bs" `printf "$PWD/${md5sum}_dir/$output_fasta"` >> "${md5sum}_dir/${bs}_samples.txt"


  done


done