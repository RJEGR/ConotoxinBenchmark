#!/bin/bash

forward_fq=$1
reverse_fq=$2

OUTDIR=$3

CPU=$4
MEM=$5

EXPORT=/LUSTRE/apps/bioinformatica/velvet
export PATH=$PATH:$EXPORT

which velveth
which velvetg

hash_length=31

call="velveth k${hash_length}_${OUTDIR} $hash_length -shortPaired -separate -fastq $forward_fq $reverse_fq"
echo $call

eval $call

call="velvetg k${hash_length}_${OUTDIR} -read_trkg yes"
echo $call

eval $call

BS=`echo $OUTDIR | awk -F'_' '{print $1"_"$2}'`

mv k${hash_length}_${OUTDIR}/contigs.fa ${BS}_FASTA_DIR/${OUTDIR%_dir}.fa 

rm -rf k${hash_length}_${OUTDIR}

exit