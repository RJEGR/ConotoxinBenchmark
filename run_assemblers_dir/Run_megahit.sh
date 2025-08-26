#!/bin/bash

forward_fq=$1
reverse_fq=$2
OUTDIR=$3
CPU=$4
MEM=$5
FINAL_DIR=$6



export PATH=/LUSTRE/apps/bioinformatica/megahit/bin/:$PATH

which megahit

rm -rf $OUTDIR # As megahit does not overwrite existing dirs, we remove it first

call="megahit -1 $forward_fq -2 $reverse_fq -o $OUTDIR -t $CPU --presets meta-sensitive"

echo "Executing: $call" 
eval $call


f1=$(find "${OUTDIR}" -maxdepth 1 -type f -name 'final.contigs.fa')

echo "Results of the assembly found at: $f1"

final_fasta=$FINAL_DIR/${OUTDIR%_dir}.fa

movecall="mv $f1 $final_fasta"

echo $movecall
eval $movecall

exit