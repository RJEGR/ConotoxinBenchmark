#!/bin/bash
#SBATCH --job-name=transrate_coverage
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00
#SBATCH -o out.%J
#SBATCH -e err.%J

REFDIR=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/fernando_pub/reads_artificiales/inputs

FASTA_DIR=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/fernando_pub/1_assembly_dir/RUN_BATCH_DIR/FASTA_DIR_BKP

FASTQ_DIR=

find "$FASTA_DIR" -maxdepth 1 -type f -name '*.fa' | while read -r ASSEMBLER; do

    FAMILYPREFIX=$(basename ${FASTA_DIR%_FASTA_DIR})
    REF=$REFDIR/${FAMILYPREFIX}.fasta
    BSCONTIG=$(basename ${ASSEMBLER%.fa})
    TRANSRATE_DIR=${BSCONTIG}_transrate_dir

    call="transrate --left $forward_fq --right $reverse_fq --assembly $ref --output $outdir --threads $threads"

    echo "Executing: $call"
    eval $call

done



# output=${contigs_a%.*}_to_${true_ref_b%.*}_transrate_dir



