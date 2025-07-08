#!/bin/bash -ex

single_reads=$1 
CPU=$2
MEM=$3
OUTDIR=$4


module load conda-2025

source activate base

source activate plass

 # assemble paired-end reads 
  #plass assemble examples/reads_1.fastq.gz examples/reads_2.fastq.gz assembly.fas tmp

  # assemble single-end reads 
plass assemble $single_reads $OUTDIR/transcripts_plass.fa $OUTDIR/tmp --threads $CPU --min-length 40 --filter-proteins 0

penguin nuclassemble $single_reads $OUTDIR/transcripts_penguin.fa $OUTDIR/tmp


exit 0