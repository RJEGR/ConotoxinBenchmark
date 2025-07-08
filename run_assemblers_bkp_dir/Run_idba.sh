#!/bin/bash -ex

single_reads=$1 
CPU=$2
MEM=$3
OUTDIR=$4


#OUTDIR=${single_reads%.*}_idba_dir

#mkdir -p $OUTDIR

export PATH=/LUSTRE/apps/bioinformatica/idba/bin:$PATH

single_reads_fa=${single_reads%.fq}.fa

awk 'NR%4==1 {print ">" substr($0, 2)} NR%4==2 {print}' $single_reads > $single_reads_fa

idba_hybrid -r $single_reads_fa --num_threads $CPU  -o $OUTDIR


exit 0
