#!/bin/bash

forward_fq=$1
reverse_fq=$2

OUTDIR=$3

CPU=$4
MEM=$5


module load conda-2025

source activate base

source activate plass


which plass
which penguin

call="plass assemble $forward_fq $reverse_fq $OUTDIR/transcripts_plass.fa $OUTDIR/tmp --threads $CPU --min-length 40 --filter-proteins 0"

echo $call

eval $call

call="penguin nuclassemble $forward_fq $reverse_fq  $OUTDIR/transcripts_penguin.fa $OUTDIR/tmp --threads $CPU"

echo $call

eval $call

BS=`echo $OUTDIR | awk -F'_' '{print $1"_"$2}'`

mv $OUTDIR/transcripts_penguin.fa ${BS}_FASTA_DIR/${OUTDIR%_PLASS_dir}_PENGUIN.fa
mv $OUTDIR/transcripts_plass.fa ${BS}_FASTA_DIR/${OUTDIR%_PLASS_dir}_PLASS.fa

exit