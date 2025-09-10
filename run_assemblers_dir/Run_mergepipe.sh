#!/bin/bash

forward_fq=$1
reverse_fq=$2
OUTDIR=$3
CPU=$4
MEM=$5
FINAL_DIR=$6

EXPORT=$PWD/run_assemblers_dir
export PATH=$PATH:$EXPORT

Run_plass.sh "$forward_fq" "$reverse_fq" "$OUTDIR" "$CPU" "$MEM" "$FINAL_DIR" &> /dev/null

f1=$(find "${FINAL_DIR}" -maxdepth 1 -type f -name '*_MERGEPIPE.fa')

echo "Results of the assembly found at: $f1"

fasta_tmp_1=$FINAL_DIR/assembly_1.fa

movecall="mv $f1 $fasta_tmp_1"

echo $movecall
eval $movecall

rm -rf $OUTDIR/*

Run_spades.sh "$forward_fq" "$reverse_fq" "$OUTDIR" "$CPU" "$MEM" "$FINAL_DIR" &> /dev/null

f1=$(find "${FINAL_DIR}" -maxdepth 1 -type f -name '*_MERGEPIPE.fa')

echo "Results of the assembly found at: $f1"

fasta_tmp_2=$FINAL_DIR/assembly_2.fa

movecall="mv $f1 $fasta_tmp_2"

echo $movecall
eval $movecall

rm -rf $OUTDIR/*

Run_trinity.sh "$forward_fq" "$reverse_fq" "$OUTDIR" "$CPU" "$MEM" "$FINAL_DIR" &> /dev/null

f1=$(find "${FINAL_DIR}" -maxdepth 1 -type f -name '*_MERGEPIPE.fa')

echo "Results of the assembly found at: $f1"

fasta_tmp_3=$FINAL_DIR/assembly_3.fa

movecall="mv $f1 $fasta_tmp_3"

echo $movecall
eval $movecall

#rm -rf $OUTDIR/*

# Concatenate all assemblies into a single file

fasta_file=$FINAL_DIR/${OUTDIR%_dir}.tmp.fa

echo "Concatenating all assemblies into a single file... $fasta_file"

cat $fasta_tmp_1 $fasta_tmp_2 $fasta_tmp_3 > $fasta_file

thread_count=20

INPUT=$fasta_file

ID=0.9

OUTFILE=$FINAL_DIR/${OUTDIR%_dir}.fa

echo "Reduce redundancy filtering sequences with identical length and 100% length overlap in $OUTFILE"


#OUTFILE=$FINAL_DIR/${INPUT%.*}.${ID}.mmseq.fasta


EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/mmseqs/bin
export PATH=$PATH:$EXPORT

mmseqs createdb $INPUT $OUTDIR/sequenceDB 

mmseqs clusthash $OUTDIR/sequenceDB  $OUTDIR/clusthashDB --threads $thread_count --min-seq-id $ID
mmseqs clust $OUTDIR/sequenceDB  $OUTDIR/clusthashDB $OUTDIR/clusterDB --threads $thread_count

mmseqs createsubdb $OUTDIR/clusterDB $OUTDIR/sequenceDB $OUTDIR/createsubDB
mmseqs convert2fasta $OUTDIR/createsubDB $OUTFILE

rm $fasta_tmp_1 $fasta_tmp_2 $fasta_tmp_3 $fasta_file


exit