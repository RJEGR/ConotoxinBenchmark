#!/bin/bash

forward_fq=$1
reverse_fq=$2

OUTDIR=$3

CPU=$4
MEM=$5


module load conda-2025
source activate base
conda activate rnabloom

which rnabloom

call="rnabloom -l $forward_fq -r $reverse_fq -t $CPU  -outdir $OUTDIR -mem $MEM"

echo $call

eval $call

BS=`echo $OUTDIR | awk -F'_' '{print $1"_"$2}'`

mv $OUTDIR/rnabloom.transcripts.fa ${BS}_FASTA_DIR/${OUTDIR%_dir}.fa

exit