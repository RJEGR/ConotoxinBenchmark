#!/bin/bash

forward_fq=$1
reverse_fq=$2

OUTDIR=$3

CPU=$4
MEM=$5


export PATH=/LUSTRE/apps/bioinformatica/idba/bin:$PATH

which idba_hybrid

forward_fa=${forward_fq%.fq}.idba.fa
reverse_fa=${reverse_fq%.fq}.idba.fa

awk 'NR%4==1 {print ">" substr($0, 2)} NR%4==2 {print}' $forward_fq > $(basename $forward_fa)
awk 'NR%4==1 {print ">" substr($0, 2)} NR%4==2 {print}' $reverse_fq > $(basename $reverse_fa)


call="idba_hybrid -r $(basename $forward_fa) $(basename $reverse_fa) --num_threads $CPU  -o $OUTDIR"

echo $call

eval $call

rm $forward_fa $reverse_fa

BS=`echo $OUTDIR | awk -F'_' '{print $1"_"$2}'`

mv $OUTDIR/scaffold.fa ${BS}_FASTA_DIR/${OUTDIR%_dir}.fa

exit
