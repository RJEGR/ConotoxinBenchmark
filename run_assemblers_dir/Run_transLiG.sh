#!/bin/bash

forward_fq=$1
reverse_fq=$2
OUTDIR=$3
CPU=$4
MEM=$5
FINAL_DIR=$6



export PATH=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/TransLiG_1.3/:$PATH

# -t <int>: minimum length of transcripts, default: 200.

call="TransLiG -s fq -p pair -l $forward_fq -r $reverse_fq -t 150 -o $OUTDIR"


echo "Executing: $call" 

eval $call


f1=$(find "${OUTDIR}" -maxdepth 1 -type f -name 'TransLiG.fa')

echo "Results of the assembly found at: $f1"
 
final_fasta=$FINAL_DIR/${OUTDIR%_dir}.fa

movecall="mv $f1 $final_fasta"

echo $movecall

eval $movecall

exit