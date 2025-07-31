#!/bin/bash

forward_fq=$1
reverse_fq=$2

OUTDIR=$3

CPU=$4
MEM=$5


EXPORT=/LUSTRE/apps/bioinformatica/oases
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/apps/bioinformatica/velvet
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/apps/bioinformatica/oases/scripts/
export PATH=$PATH:$EXPORT

which velveth
which velvetg
which oases
# oases pipeline. Using hash reads for k=31 and k=23

hash_length=31
velveth k${hash_length}_${OUTDIR} $hash_length -shortPaired -separate -fastq $forward_fq $reverse_fq
velvetg k${hash_length}_${OUTDIR} -read_trkg yes
oases k${hash_length}_${OUTDIR}

hash_length=23
velveth k${hash_length}_${OUTDIR} $hash_length -shortPaired -separate -fastq $forward_fq $reverse_fq
velvetg k${hash_length}_${OUTDIR} -read_trkg yes
oases k${hash_length}_${OUTDIR}

# merge the above assemblies with a nominal k=23
velveth mergedAssembly_${OUTDIR} 23 -long k*_${OUTDIR}/*.fa
velvetg mergedAssembly_${OUTDIR} -read_trkg yes -conserveLong yes
oases mergedAssembly_${OUTDIR} -merge yes

exit


