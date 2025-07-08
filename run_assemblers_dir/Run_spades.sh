#!/bin/bash -ex

single_reads=$1 
CPU=$2
MEM=$3
OUTDIR=$4

#OUTDIR=${single_reads%.*}_spades_dir

#mkdir -p $OUTDIR

EXPORT=/LUSTRE/apps/bioinformatica/SPAdes-3.15.5-Linux/bin/
export PATH=$PATH:$EXPORT


rnaspades.py -s $single_reads -t $CPU -m $MEM -o $OUTDIR