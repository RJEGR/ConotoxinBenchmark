#!/bin/bash

forward_fq=$1
reverse_fq=$2
OUTDIR=$3
CPU=$4




export PATH=/LUSTRE/apps/bioinformatica/megahit/bin/:$PATH

which megahit

# BinPacker -s fq -p pair -l reads.left.fq -r reads.right.fq

rm -rf $OUTDIR

call="megahit -1 $forward_fq -2 $reverse_fq -o $OUTDIR -t $CPU --presets meta-sensitive"

echo $call

eval $call

BS=`echo $(basename $OUTDIR) | awk -F'_' '{print $1"_"$2}'`

# mv $OUTDIR/final.contigs.fa ${BS}_FASTA_DIR/${OUTDIR%_dir}.fa
movecall="mv $OUTDIR/final.contigs.fa ${BS}_FASTA_DIR/${OUTDIR%_dir}.fa"
echo $movecall
eval $movecall

exit