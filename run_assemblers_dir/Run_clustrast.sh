#!/bin/bash

forward_fq=$1
reverse_fq=$2
OUTDIR=$3
CPU=$4
MEM=$5
FINAL_DIR=$6

module load conda-2025

source activate base

# conda create --name shannon_cpp_env
# conda activate shannon_cpp_env

conda activate transabbys

which transabyss

EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/minimap2-2.30_x64-linux
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/clustrast-main
export PATH=$PATH:$EXPORT


exit


all="transabyss --pe $forward_fq $reverse_fq --threads $CPU --outdir $OUTDIR"

echo "Executing: $call" 
eval $call


f1=$(find "${OUTDIR}" -maxdepth 1 -type f -name 'transabyss-final.fa')

echo "Results of the assembly found at: $f1"
 
final_fasta=$FINAL_DIR/${OUTDIR%_dir}.fa

movecall="mv $f1 $final_fasta"

echo $movecall

eval $movecall

exit