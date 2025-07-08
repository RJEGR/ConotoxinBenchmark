#!/bin/bash -ex

single_reads=$1 
CPU=$2
MEM=$3
OUTDIR=$4

EXPORT=/LUSTRE/apps/bioinformatica/velvet
export PATH=$PATH:$EXPORT


#Short single end reads:

hash_length=31

velveth $OUTDIR $hash_length -short -fastq $single_reads
velvetg $OUTDIR 

#Paired-end short reads (using separate files for the paired reads)
#	velveth Assem 31 -shortPaired -fasta -separate left.fa right.fa

exit