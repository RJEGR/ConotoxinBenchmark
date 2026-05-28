#!/bin/bash

forward_fq=$1
reverse_fq=$2

OUTDIR=$3

CPU=$4
MEM=$5

FINAL_DIR=$6


EXPORT=/LUSTRE/apps/bioinformatica/SPAdes-3.15.5-Linux/bin/
export PATH=$PATH:$EXPORT

which rnaspades.py

call="rnaspades.py -1 $forward_fq -2 $reverse_fq  -o $OUTDIR -t $CPU -m $MEM"

echo "Executing: $call"

eval $call


# Not recursesive search in find. Only the files with the name 'transcripts.fasta' are searched at the $OUTDIR level.

f1=$(find "${OUTDIR}" -maxdepth 1 -type f -name 'transcripts.fasta') 

echo "Results of the assembly found at: $f1"
 
final_fasta=$FINAL_DIR/${OUTDIR%_dir}.fa

movecall="mv $f1 $final_fasta"

echo $movecall
eval $movecall

exit

