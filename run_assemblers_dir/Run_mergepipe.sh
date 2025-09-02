#!/bin/bash

forward_fq=$1
reverse_fq=$2
OUTDIR=$3
CPU=$4
MEM=$5
FINAL_DIR=$6



module load trinityrnaseq-v2.15.1

which Trinity

#call="Trinity --seqType fq --max_memory 100G --left $forward_fq --right $reverse_fq --no_normalize_reads --CPU $CPU --output trinity_out_dir.${OUTDIR} --full_cleanup"

output_dir=${OUTDIR}/Trinity_out_dir # Trinity_out_dir prefix is mandatory for Trinity

call="Trinity --seqType fq --max_memory 100G --left $forward_fq --right $reverse_fq --no_normalize_reads --CPU 24 --output $output_dir --full_cleanup"
   

echo "Executing: $call"

eval $call

f1=$(find "${OUTDIR}" -maxdepth 1 -type f -name 'Trinity_out_dir.Trinity.fasta')

echo "Results of the assembly found at: $f1"
 
final_fasta=$FINAL_DIR/${OUTDIR%_dir}.fa

movecall="mv $f1 $final_fasta"

echo $movecall
eval $movecall

#rm *.fasta.gene_trans_map

exit
