#!/bin/bash

forward_fq=$1
reverse_fq=$2

OUTDIR=$3

CPU=$4
MEM=$5

module load conda-2025

source activate base

conda activate transabbys

# transabbys use samtools lower version than trinity, and trinity crash, therefore does not export PATHS in parallel

which transabyss

call="transabyss --pe $forward_fq $reverse_fq --threads $CPU --outdir $OUTDIR"

echo $call

eval $call

BS=`echo $OUTDIR | awk -F'_' '{print $1"_"$2}'`

mv $OUTDIR/transabyss-final.fa ${BS}_FASTA_DIR/${OUTDIR%_dir}.fa

exit