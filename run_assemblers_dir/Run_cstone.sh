#!/bin/bash

forward_fq=$1
reverse_fq=$2
OUTDIR=$3

CPU=$4
MEM=$5



#export PATH=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/cstone:$PATH

WD=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/cstone
cstone="java -jar -Xmx64G $WD/cstone.jar"


call="$cstone -r1 $forward_fq -r2 $reverse_fq -p $CPU  -o $OUTDIR"

echo $call

eval $call

BS=`echo $OUTDIR | awk -F'_' '{print $1"_"$2}'

movecall="mv $OUTDIR/contigs.fasta ${BS}_FASTA_DIR/${OUTDIR%_dir}.fa"
echo $movecall
eval $movecall

exit