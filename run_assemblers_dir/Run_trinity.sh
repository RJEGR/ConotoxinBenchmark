#!/bin/bash -ex

single_reads=$1 
CPU=$2
MEM=$3
OUTDIR=$4

#OUTDIR=${single_reads%.*}_trinity_dir

#mkdir -p $OUTDIR

module load trinityrnaseq-v2.15.1

Trinity --seqType fq --max_memory 100G --single $single_reads --no_normalize_reads --CPU $CPU --output $OUTDIR --full_cleanup

exit 0
