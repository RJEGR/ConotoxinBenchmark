#!/bin/bash

forward_fq=$1
reverse_fq=$2
OUTDIR=$3
CPU=$4




export PATH=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/Bridger_r2014-12-01:$PATH

which Bridger.pl
#--output/-o

call="Bridger.pl --seqType fq --left $forward_fq --right $reverse_fq --output $OUTDIR --CPU $CPU"

echo $call

eval $call

BS=`echo $OUTDIR | awk -F'_' '{print $1"_"$2}'`

# mv $OUTDIR/Bridger.fasta ${BS}_FASTA_DIR/${OUTDIR%_dir}.fa
movecall="mv $OUTDIR/Bridger.fasta ${BS}_FASTA_DIR/${OUTDIR%_dir}.fa"
echo $movecall
eval $movecall

exit