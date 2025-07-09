#!/bin/bash

forward_fq=$1
reverse_fq=$2

OUTDIR=$3

CPU=$4
MEM=$5

module load conda-2025
source activate base
conda activate cap3 

which cap3

forward_fa=${forward_fq%.fq}_cap.fa
reverse_fa=${reverse_fq%.fq}_cap.fa

awk 'NR%4==1 {print ">" substr($0, 2)} NR%4==2 {print}' $forward_fq > $forward_fa
awk 'NR%4==1 {print ">" substr($0, 2)} NR%4==2 {print}' $reverse_fq > $reverse_fa

call="cap3 $forward_fa $reverse_fa"

echo $call

eval $call

mv *.cap.* $OUTDIR

rm $forward_fa $reverse_fa

exit