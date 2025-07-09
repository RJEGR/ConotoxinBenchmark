#!/bin/bash

forward_fq=$1
reverse_fq=$2

OUTDIR=$3

CPU=$4
MEM=$5


module load trinityrnaseq-v2.15.1

which Trinity

call="Trinity --seqType fq --max_memory 100G --left $forward_fq --right $reverse_fq --no_normalize_reads --CPU $CPU --output trinity_out_dir.${OUTDIR} --full_cleanup"

echo $call

eval $call

BS=`echo $OUTDIR | awk -F'_' '{print $1"_"$2}'`

mv trinity_out_dir.${OUTDIR}.Trinity.fasta ${BS}_FASTA_DIR/${OUTDIR%_dir}.fa

rm *.fasta.gene_trans_map

exit
