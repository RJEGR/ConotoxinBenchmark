#!/bin/bash -ex

single_reads=$1 
CPU=$2
MEM=$3
OUTDIR=$4

#OUTDIR=${single_reads%.*}_rnabloom_dir

#mkdir -p $OUTDIR

module load conda-2025

# conda activate rnabloom

export PATH=/LUSTRE/apps/Anaconda/2025/miniconda3/envs/rnabloom/bin/:$PATH


rnabloom -sef $single_reads -t $CPU  -outdir $OUTDIR