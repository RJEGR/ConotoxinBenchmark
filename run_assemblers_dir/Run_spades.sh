#!/bin/bash

forward_fq=$1
reverse_fq=$2

OUTDIR=$3

CPU=$4
MEM=$5


EXPORT=/LUSTRE/apps/bioinformatica/SPAdes-3.15.5-Linux/bin/
export PATH=$PATH:$EXPORT

which rnaspades.py

call="rnaspades.py -1 $forward_fq -2 $reverse_fq  -o $OUTDIR -t $CPU -m $MEM"

echo $call

eval $call

BS=`echo $OUTDIR | awk -F'_' '{print $1"_"$2}'`

mv $OUTDIR/transcripts.fasta ${BS}_FASTA_DIR/${OUTDIR%_dir}.fa

exit