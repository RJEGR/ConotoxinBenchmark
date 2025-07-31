#!/bin/bash

forward_fq=$1
reverse_fq=$2
OUTDIR=$3

CPU=$4
MEM=$5



#export PATH=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/cstone:$PATH
rm -rf $OUTDIR

WD=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/
tool="java -jar -Xmx64G $WD/vt-builder_r8clean.jar"

java -jar -Xmx97152M vt-builder_r8clean.jar -a n -b /LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/./I3_superfamily_PE1.fq -c /LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/./I3_superfamily_PE2.fq -d ./output.txt -e ./log.txt -f 24 -g 150 -h 250 -j 96

-d ./output.txt -e ./log.txt -f 24 -g 150 -h 250 -j 96

call="$tool -a n -b $forward_fq -c $reverse_fq -d $OUTDIR -f $CPU -j $MEM -g 100 -h 200"

echo $call

eval $call

BS=`echo $(basename $OUTDIR) | awk -F'_' '{print $1"_"$2}'`

movecall="mv $OUTDIR/contigs.fasta ${BS}_FASTA_DIR/${OUTDIR%_dir}.fa"
echo $movecall
eval $movecall

exit