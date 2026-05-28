#!/bin/bash

forward_fq=$1
reverse_fq=$2
OUTDIR=$3
CPU=$4
MEM=$5
FINAL_DIR=$6
REFERENCE=$7

BS=${OUTDIR%_dir}

REF_PREFIX=$(basename ${REFERENCE%.f*})


EXPORT=/LUSTRE/apps/bioinformatica/hisat2/
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/apps/bioinformatica/RSEM/bin/samtools-1.3
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/apps/bioinformatica/stringtie/
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/apps/gffread
export PATH=$PATH:$EXPORT


# 1) Create index if does not exist

mkdir -p ${OUTDIR}/INDEX_dir

if [ ! -f "${OUTDIR}/INDEX_dir/${REF_PREFIX}.1.ht2" ]; then
hisat2-build --quiet -p $CPU $REFERENCE ${OUTDIR}/INDEX_dir/$REF_PREFIX

else
    echo "index for '$REFERENCE' already exists."
fi

call="hisat2 --quiet --phred33 --dta -p $CPU -x ${OUTDIR}/INDEX_dir/$REF_PREFIX -1 $forward_fq -2 $reverse_fq -S ${OUTDIR}/${BS}.sam"

echo "Executing: $call"

eval $call

samtools sort -@ $CPU -o ${OUTDIR}/${BS}.sorted.bam ${OUTDIR}/${BS}.sam

rm -rf ${OUTDIR}/${BS}.sam

stringtie --rf -p $CPU -o ${OUTDIR}/${BS}_transcripts.gtf ${OUTDIR}/${BS}.sorted.bam


# 4) Get fasta (This is addit. step to comparative transcriptomic step)

gffread -w ${OUTDIR}/transcripts.fasta -g $REFERENCE ${OUTDIR}/${BS}_transcripts.gtf


f1=$(find "${OUTDIR}" -maxdepth 1 -type f -name 'transcripts.fasta')

echo "Results of the assembly found at: $f1"
 
final_fasta=$FINAL_DIR/${OUTDIR%_dir}.fa

movecall="mv $f1 $final_fasta"

echo $movecall

eval $movecall

exit
