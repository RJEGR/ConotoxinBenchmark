#!/bin/bash

forward_fq=$1
reverse_fq=$2
OUTDIR=$3
CPU=$4




export PATH=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/BinPacker_1.0:$PATH

which BinPacker

rm -rf $OUTDIR

call="BinPacker -s fq -p pair -l $forward_fq -r $reverse_fq -o $OUTDIR"

echo $call

eval $call

BS=`echo $(basename $OUTDIR) | awk -F'_' '{print $1"_"$2}'`

movecall="mv $OUTDIR/BinPacker.fa ${BS}_FASTA_DIR/${OUTDIR%_dir}.fa"
echo $movecall
eval $movecall

exit