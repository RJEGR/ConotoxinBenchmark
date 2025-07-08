#!/bin/bash -ex

single_reads=$1 
CPU=$2
MEM=$3
OUTDIR=$4


module load conda-2025

source activate base

conda activate transabbys

# transabbys use samtools lower version than trinity, and trinity crash, therefore does not export PATHS in parallel

transabyss --se $single_reads --threads $CPU --outdir $OUTDIR

exit 0