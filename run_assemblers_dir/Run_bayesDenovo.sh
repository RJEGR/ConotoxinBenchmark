#!/bin/bash

forward_fq=$1
reverse_fq=$2
OUTDIR=$3

CPU=$4
MEM=$5



export PATH=/LUSTRE/apps/bioinformatica/BayesDenovo/BayesDenovo-v1:$PATH

which BayesDenovo_v1.pl

call="BayesDenovo_v1.pl --seqType fq --left $forward_fq --right $reverse_fq --CPU $CPU  --output $OUTDIR --bamfile $OUTDIR/bamfile.bam --fafile $OUTDIR/bayesdenovo.fa"

echo $call

eval $call

# mv $OUTDIR/bayesdenovo.fa ASSEMBLIES_DIR/BAYES.fa

exit