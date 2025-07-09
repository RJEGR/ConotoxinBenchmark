#!/bin/bash

forward_fq=$1
reverse_fq=$2

OUTDIR=$3

CPU=$4
MEM=$5

configureFile=$(basename ${forward_fq%.*})_SOAPdenovo.config

cat <<EOF > $configureFile
#maximal read length
max_rd_len=155
[LIB]
#maximal read length in this lib
rd_len_cutof=70
#average insert size
avg_ins=300
#if sequence needs to be reversed 
reverse_seq=0
#in which part(s) the reads are used
asm_flags=3
#minimum aligned length to contigs for a reliable read location (at least 32 for short insert size)
map_len=32
#fastq file for paired-end reads
q1=$forward_fq 
q2=$reverse_fq
EOF

module load conda-2024_py3.8
source activate soapdenovo-trans

which SOAPdenovo-Trans-31mer

#export PATH=/LUSTRE/apps/Anaconda/2024/miniconda3-py3.8/envs/soapdenovo-trans/bin/:$PATH

call="SOAPdenovo-Trans-31mer all -s $configureFile -p $CPU  -o $OUTDIR"

echo $call

eval $call

# BS=`echo $OUTDIR | awk -F'_' '{print $1"_"$2}'`

# mv $OUTDIR/rnabloom_assembly.fasta ${BS}_FASTA_DIR/${OUTDIR%_dir}.fa

exit
