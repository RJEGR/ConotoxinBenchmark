#!/bin/bash
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00

FASTQC=/LUSTRE/apps/bioinformatica/FastQC_v0.12.1/
export PATH=$PATH:$FASTQC


#Manifest=$1

run_paired_fastqc() {
    
    local manifest="$1"
    local outdir="$2"

    local forward_fq="${manifest%.*}_concat_PE1.fq"
    local reverse_fq="${manifest%.*}_concat_PE2.fq"

    # Concatenate forward reads
    local call=$(awk '{print $2}' "$manifest" | tr "\n" " ")
    
    # Filter out empty entries before concatenation
    for fq in $call; do [[ -n "$fq" ]] && echo "$fq"; done | xargs cat > "$forward_fq"

    # Concatenate reverse reads
    local call=$(awk '{print $3}' "$manifest" | tr "\n" " ")
    for fq in $call; do [[ -n "$fq" ]] && echo "$fq"; done | xargs cat > "$reverse_fq"

    fastqc --nogroup -t 20 -o $outdir $forward_fq 
    
    fastqc --nogroup -t 20 -o $outdir $reverse_fq 

    rm -f "${manifest%.*}_concat_PE"*.fq
    
}

mkdir -p fastqc_dir multiqc_dir

find "$PWD" -maxdepth 1 -name '*.txt' | while read -r Manifest; do 
    call="run_paired_fastqc "$Manifest" "fastqc_dir" "
    echo "Executing: $call"
    eval $call
done

exit

multiqc ./fastqc_dir/*zip -o multiqc_dir
