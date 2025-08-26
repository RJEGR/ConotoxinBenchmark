#!/bin/bash

forward_fq=$1
reverse_fq=$2
OUTDIR=$3

CPU=$4
MEM=$5

FINAL_DIR=$6

#export PATH=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/cstone:$PATH
rm -rf $OUTDIR

WD=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/cstone
cstone="java -jar -Xmx64G $WD/cstone.jar"


call="$cstone -r1 $forward_fq -r2 $reverse_fq -o $OUTDIR -p $CPU"

echo "Executing: $call" 
eval $call


f1=$(find "${OUTDIR}" -maxdepth 1 -type f -name 'contigs.fasta')

echo "Results of the assembly found at: $f1"

final_fasta=$FINAL_DIR/${OUTDIR%_dir}.fa

movecall="mv $f1 $final_fasta"

echo $movecall
eval $movecall

exit