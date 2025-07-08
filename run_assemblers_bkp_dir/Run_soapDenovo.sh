#!/bin/bash -ex

single_reads=$1 
CPU=$2
MEM=$3
OUTDIR=$4

configureFile=${single_reads%.*}_SOAPdenovo.config

cat <<EOF > $configureFile
#maximal read length
max_rd_len=50
[LIB]
#maximal read length in this lib
rd_len_cutof=45
#average insert size
avg_ins=100
#if sequence needs to be reversed 
reverse_seq=0
#in which part(s) the reads are used
asm_flags=3
#minimum aligned length to contigs for a reliable read location (at least 32 for short insert size)
map_len=32
#fastq file for single reads
q=$PWD/$single_reads 
EOF

module load conda-2024_py3.8
# source activate soapdenovo-trans

export PATH=/LUSTRE/apps/Anaconda/2024/miniconda3-py3.8/envs/soapdenovo-trans/bin/:$PATH

SOAPdenovo-Trans-31mer all -s $configureFile -p $CPU  -o $OUTDIR"

exit 0
