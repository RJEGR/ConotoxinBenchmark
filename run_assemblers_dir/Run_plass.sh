#!/bin/bash

forward_fq=$1
reverse_fq=$2
OUTDIR=$3
CPU=$4
MEM=$5
FINAL_DIR=$6

module load conda-2025

source activate base

source activate plass


#which plass
which penguin

#call="plass assemble $forward_fq $reverse_fq $OUTDIR/transcripts_plass.fa $OUTDIR/tmp --threads $CPU --min-length 40 --filter-proteins 0"

#echo $call

eval $call

# --rescore-mode 2
# --min-seq-id 0.9

call="penguin nuclassemble $forward_fq $reverse_fq  $OUTDIR/transcripts_penguin.fa $OUTDIR/tmp --threads $CPU --min-contig-len 100"

#call="penguin nuclassemble $forward_fq $reverse_fq  $OUTDIR/transcripts_penguin.fa $OUTDIR/tmp --threads $CPU --min-contig-len 100 --rescore-mode 2 --min-seq-id 0.9"


echo "Executing: $call" 
eval $call


f1=$(find "${OUTDIR}" -maxdepth 1 -type f -name 'transcripts_penguin.fa')

echo "Results of the assembly found at: $f1"

final_fasta=$FINAL_DIR/${OUTDIR%_dir}.fa

movecall="mv $f1 $final_fasta"

echo $movecall
eval $movecall

exit