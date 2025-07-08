#!/bin/bash -ex

single_reads=$1 
CPU=$2
MEM=$3
OUTDIR=$4

EXPORT=/LUSTRE/apps/bioinformatica/oases
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/apps/bioinformatica/velvet
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/apps/bioinformatica/oases/scripts/
export PATH=$PATH:$EXPORT


# oases pipeline

# hash reads for k=31 and k=23

hash_length=31
velveth k${hash_length}_${OUTDIR} $hash_length -short -fastq $single_reads
velvetg k${hash_length}_${OUTDIR} -read_trkg yes
oases k${hash_length}_${OUTDIR}

hash_length=23
velveth k${hash_length}_${OUTDIR} $hash_length -short -fastq $single_reads
velvetg k${hash_length}_${OUTDIR} -read_trkg yes
oases k${hash_length}_${OUTDIR}

# merge the above assemblies with a nominal k=23
velveth mergedAssembly_${OUTDIR} 23 -long k*/*.fa
velvetg mergedAssembly_${OUTDIR} -read_trkg yes -conserveLong yes
oases mergedAssembly_${OUTDIR} -merge yes


exit

#Short single end reads:

hash_length=31

velveth $OUTDIR $hash_length -short -fastq $single_reads
velvetg $OUTDIR -read_trkg yes
oases $OUTDIR


