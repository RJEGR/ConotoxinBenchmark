#!/bin/bash -ex

single_reads=$1 
CPU=$2
MEM=$3
OUTDIR=$4

module load conda-2025
source activate base
conda activate cap3 

single_reads_fa=${single_reads%.fq}.cap.fa

awk 'NR%4==1 {print ">" substr($0, 2)} NR%4==2 {print}' $single_reads > $single_reads_fa

cap3 $single_reads_fa

mv *.cap.* $OUTDIR

rm $single_reads_fa

exit 0