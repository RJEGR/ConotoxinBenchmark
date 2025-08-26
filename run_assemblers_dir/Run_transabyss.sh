#!/bin/bash

forward_fq=$1
reverse_fq=$2
OUTDIR=$3
CPU=$4
MEM=$5
FINAL_DIR=$6

module load conda-2025

source activate base

conda activate transabbys

# transabbys use samtools lower version than trinity, and trinity crash, therefore does not export PATHS in parallel

which transabyss

call="transabyss --pe $forward_fq $reverse_fq --threads $CPU --outdir $OUTDIR"

echo "Executing: $call" 
eval $call


f1=$(find "${OUTDIR}" -maxdepth 1 -type f -name 'transabyss-final.fa')

echo "Results of the assembly found at: $f1"
 
final_fasta=$FINAL_DIR/${OUTDIR%_dir}.fa

movecall="mv $f1 $final_fasta"

echo $movecall

eval $movecall

exit