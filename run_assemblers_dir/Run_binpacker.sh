#!/bin/bash

forward_fq=$1
reverse_fq=$2
OUTDIR=$3

CPU=$4




export PATH=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/BinPacker_1.0:$PATH

which BinPacker

BinPacker -s fq -p pair -l reads.left.fq -r reads.right.fq

OUTDIR=$3
call="BinPacker -s fq -p pair -l $forward_fq -r $reverse_fq -o $OUTDIR"

echo $call

eval $call

# mv $OUTDIR/bayesdenovo.fa ASSEMBLIES_DIR/BAYES.fa

exit