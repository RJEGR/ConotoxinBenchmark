#!/bin/bash

forward_fq=$1
reverse_fq=$2
OUTDIR=$3

CPU=$4




export PATH=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/Bridger_r2014-12-01:$PATH

which Bridger.pl

perl perl Bridger.pl --seqType fq --left reads1.fq --right reads2.fq --CPU 6 

call="Bridger.pl --seqType fq --left $forward_fq --right $reverse_fq --CPU $CPU"

echo $call

eval $call

# mv $OUTDIR/bayesdenovo.fa ASSEMBLIES_DIR/BAYES.fa

exit