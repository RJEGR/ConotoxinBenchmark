#!/bin/bash

forward_fq=$1
reverse_fq=$2
OUTDIR=$3
CPU=$4
MEM=$5
FINAL_DIR=$6



export PATH=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/Bridger_r2014-12-01:$PATH

which Bridger.pl

rm -rf $OUTDIR

call="Bridger.pl --seqType fq --left $forward_fq --right $reverse_fq --output $OUTDIR --CPU $CPU"


echo "Executing: $call" 

eval $call


f1=$(find "${OUTDIR}" -maxdepth 1 -type f -name 'Bridger.fasta')

echo "Results of the assembly found at: $f1"
 
final_fasta=$FINAL_DIR/${OUTDIR%_dir}.fa

movecall="mv $f1 $final_fasta"

echo $movecall

eval $movecall

exit