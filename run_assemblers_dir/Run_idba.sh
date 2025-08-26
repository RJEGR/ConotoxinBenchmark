#!/bin/bash

forward_fq=$1
reverse_fq=$2

OUTDIR=$3

CPU=$4
MEM=$5
FINAL_DIR=$6

export PATH=/LUSTRE/apps/bioinformatica/idba/bin:$PATH

which idba_hybrid

forward_fa=${forward_fq%.fq}.idba.fa
reverse_fa=${reverse_fq%.fq}.idba.fa

awk 'NR%4==1 {print ">" substr($0, 2)} NR%4==2 {print}' $forward_fq > $(basename $forward_fa)
awk 'NR%4==1 {print ">" substr($0, 2)} NR%4==2 {print}' $reverse_fq > $(basename $reverse_fa)


call="idba_hybrid -r $(basename $forward_fa) $(basename $reverse_fa) --num_threads $CPU  -o $OUTDIR"

echo "Executing: $call" 
eval $call


rm $forward_fa $reverse_fa

f1=$(find "${OUTDIR}" -maxdepth 1 -type f -name 'scaffold.fa')

echo "Results of the assembly found at: $f1"
 
final_fasta=$FINAL_DIR/${OUTDIR%_dir}.fa

movecall="mv $f1 $final_fasta"

echo $movecall
eval $movecall

exit

