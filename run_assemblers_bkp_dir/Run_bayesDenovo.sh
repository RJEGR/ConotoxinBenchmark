#!/bin/bash -ex

single_reads=$1 
CPU=$2
MEM=$3
OUTDIR=$4

export PATH=/LUSTRE/apps/bioinformatica/BayesDenovo/BayesDenovo-v1:$PATH

#OUTDIR=${single_reads%.*}_bayesDenovo_dir

#mkdir -p $OUTDIR

BayesDenovo_v1.pl --seqType fq --single $single_reads --CPU $CPU  --output $OUTDIR --bamfile $OUTDIR/bamfile.bam --fafile $OUTDIR/bayesdenovo.fa