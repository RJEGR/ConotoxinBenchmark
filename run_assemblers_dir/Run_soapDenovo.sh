#!/bin/bash

forward_fq=$1
reverse_fq=$2
OUTDIR=$3
CPU=$4
MEM=$5
FINAL_DIR=$6

function abspath() {
    if [ -d "$1" ]; then
        (cd "$1"; pwd)
    elif [ -f "$1" ]; then
        if [[ $1 == */* ]]; then
            echo "$(cd "${1%/*}"; pwd)/${1##*/}"
        else
            echo "$(pwd)/$1"
        fi
    fi
}

forward_fq=$(abspath $forward_fq)
reverse_fq=$(abspath $reverse_fq)


echo "Forward reads: $forward_fq"
echo "Reverse reads: $reverse_fq"

forward_fa=${forward_fq%.fq}.fa
reverse_fa=${reverse_fq%.fq}.fa

awk 'NR%4==1 {print ">" substr($0, 2)} NR%4==2 {print}' $forward_fq > $(basename $forward_fa)
awk 'NR%4==1 {print ">" substr($0, 2)} NR%4==2 {print}' $reverse_fq > $(basename $reverse_fa)



configureFile=$(basename ${forward_fa%.*})_SOAPdenovo.config

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
#q1=$forward_fq
#q2=$reverse_fq
#fasta file for paired-end reads
f1=$forward_fa
f2=$reverse_fa
EOF

module load conda-2024_py3.8
source activate soapdenovo-trans

#source  /LUSTRE/bioinformatica_data/abtx/fernando/programas/miniconda3/bin/activate soapdenovo

which SOAPdenovo-Trans-31mer

#export PATH=/LUSTRE/apps/Anaconda/2024/miniconda3-py3.8/envs/soapdenovo-trans/bin/:$PATH

call="SOAPdenovo-Trans-31mer all -s $configureFile -p $CPU  -o $OUTDIR/out"

echo "Executing: $call"

eval $call

rm $forward_fa $reverse_fa

f1=$(find "${OUTDIR}" -maxdepth 1 -type f -name 'out.scafSeq')

echo "Results of the assembly found at: $f1"
 
final_fasta=$FINAL_DIR/${OUTDIR%_dir}.fa

movecall="mv $f1 $final_fasta"

echo $movecall

eval $movecall

exit
