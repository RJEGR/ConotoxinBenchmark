#!/bin/bash
#SBATCH --job-name=transrate_coverage
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00
#SBATCH -o out.%J
#SBATCH -e err.%J


REFDIR=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/fernando_pub/reads_artificiales/inputs

run_transrate() {
    Manifest=$1

    FASTA_DIR="${Manifest%.*}_FASTA_DIR"
    forward_fq=${Manifest%.*}_concat_PE1.fq
    reverse_fq=${Manifest%.*}_concat_PE2.fq

    call=$(awk '{print $2}' "$Manifest" | tr "\n" " ")
    cat $call > $forward_fq

    call=$(awk '{print $3}' "$Manifest" | tr "\n" " ")
    cat $call > $reverse_fq

    mkdir -p 2_transrate_dir

    find "$FASTA_DIR" -maxdepth 1 -type f -name '*.fa' | while read -r asm; do
        FAMILYPREFIX=$(basename "${FASTA_DIR%_FASTA_DIR}")
        ref=$REFDIR/${FAMILYPREFIX}.fasta
        BSCONTIG=$(basename "${asm%.fa}")
        TRANSRATE_DIR=${BSCONTIG}_transrate_dir

        call="transrate --left $forward_fq --right $reverse_fq --assembly $asm --reference $ref --output 2_transrate_dir/$TRANSRATE_DIR --threads 20"

        echo "Executing: $call"
        eval $call
    done

    rm "${Manifest%.*}_concat_PE"*.fq
}

# Loop over every *.txt file and run run_transrate
for Manifest in *.txt; do
    run_transrate "$Manifest"
done

exit
