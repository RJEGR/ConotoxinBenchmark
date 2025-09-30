#!/bin/bash
#SBATCH --job-name=FASTP
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00

EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software
export PATH=$PATH:$EXPORT

NPROCS=$SLURM_NPROCS


mkdir -p MULTIQC_VIZ_DIR

mkdir -p FASTP_OUT_DIR

mkdir -p CHKPNT_DIR

for f in $(ls *_1.fastq)
do
basename=${f##*/}
bs="${basename%_1.fastq}"
infile="${f%_1.fastq}"

left_file=${infile}_1.fastq
right_file=${infile}_2.fastq

# Consider to
# --disable_length_filtering: length filtering is enabled by default. If this option is specified, length filtering is disabled
# --length_required: reads shorter than length_required will be discarded, default is 15. (int [=15])

if [ ! -f "CHKPNT_DIR/${bs}_fastp.chkpt" ]; then

    call="fastp --thread $NPROCS \
    --detect_adapter_for_pe \
    --disable_length_filtering \
    --trim_poly_x  \
    --json MULTIQC_VIZ_DIR/${bs}_fastp.json \
    --html MULTIQC_VIZ_DIR/${bs}_fastp.html \
    -i $left_file -I $right_file \
    -o FASTP_OUT_DIR/${bs}_R1.fq.gz  -O FASTP_OUT_DIR/${bs}_R2.fq.gz"

    echo "Running fastp in '${bs}' sample."

    echo $call

    eval $call

    # Merge and dedup
    call="fastp --thread $NPROCS \
    --disable_length_filtering \
    --dont_eval_duplication \
    --json MULTIQC_VIZ_DIR/${bs}_fastp.json.tmp.merged \
    --html MULTIQC_VIZ_DIR/${bs}_fastp.html.tmp.merged \
    -i FASTP_OUT_DIR/${bs}_R1.fq.gz -I FASTP_OUT_DIR/${bs}_R2.fq.gz \
    --dedup --dup_calc_accuracy 3 \
    --merge --merged_out FASTP_OUT_DIR/${bs}_dedup_merged.fq.gz"

    echo "Merging and dedup fastp in '${bs}' sample."

    eval $call

    rm -f MULTIQC_VIZ_DIR/${bs}_fastp.*.tmp.merged


    unlink $left_file
    unlink $right_file

    touch  CHKPNT_DIR/${bs}_fastp.chkpt

else
    echo "'CHKPNT_DIR/${bs}_fastp.chkpt' already exists."
    echo "Continue with next sample."
fi


done

WDM=/LUSTRE/apps/Anaconda/2023/miniconda3/bin/

$WDM/multiqc MULTIQC_VIZ_DIR/*_fastp.json  -o MULTIQC_VIZ_DIR --config multiqc_info.conf

exit