#!/bin/bash

forward_fq=$1
reverse_fq=$2
OUTDIR=$3
CPU=$4
MEM=$5
FINAL_DIR=$6

# 1) Carga Conda y activa el entorno Oases
CONDA_SH="/LUSTRE/bioinformatica_data/abtx/fernando/programas/miniconda3/etc/profile.d/conda.sh"
source "$CONDA_SH"
conda activate oases

# 2) Parámetros de k
KMIN=19; KMAX=31; KSTEP=4; KMERGE=27

# 3) Inserto y cobertura
INS_LEN=300; COV_CUTOFF=1; MIN_TRANS=100

# 4) Directorios
#WORKDIR="/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/fernando_pub/1_assembly_dir/RUN_BATCH_DIR/oases_auto2"
#mkdir -p "$OUTDIR"

echo "=== Ejecutando pipeline Oases para $PREFIX ==="

call="python /LUSTRE/bioinformatica_data/abtx/fernando/programas/miniconda3/envs/oases/bin/oases_pipeline.py \
    --kmin="$KMIN" --kmax="$KMAX" --kstep="$KSTEP" --merge="$KMERGE" \
    --output="$OUTDIR" \
    -d -shortPaired "$forward_fq" "$reverse_fq" \
    -p -ins_length "$INS_LEN" -cov_cutoff "$COV_CUTOFF" -min_trans_lgth "$MIN_TRANS""


echo "Executing: $call"

eval $call

f1=$(find "${OUTDIR}" -maxdepth 1 -type f -name 'transcripts.fasta')

echo "Results of the assembly found at: $f1"
 
final_fasta=$FINAL_DIR/${OUTDIR%_dir}.fa

movecall="mv $f1 $final_fasta"

echo $movecall

eval $movecall
echo "✓ Pipeline Oases completado en $WORKDIR."