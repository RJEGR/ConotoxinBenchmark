#!/bin/bash

forward_fq=$1
reverse_fq=$2

OUTDIR=$3

CPU=$4
MEM=$5


export PATH=/LUSTRE/apps/bioinformatica/idba/bin:$PATH

forward_fa=${forward_fq%.fq}.idba.fa
reverse_fa=${basename reverse_fq%.fq}.idba.fa

awk 'NR%4==1 {print ">" substr($0, 2)} NR%4==2 {print}' $forward_fq > $(basename $forward_fa)
awk 'NR%4==1 {print ">" substr($0, 2)} NR%4==2 {print}' $reverse_fq > $(basename $reverse_fa)

idba_hybrid -r $(basename $forward_fa) $(basename $reverse_fa) --num_threads $CPU  -o $OUTDIR

rm $forward_fa $reverse_fa

exit
