#!/bin/bash
#SBATCH --job-name=TrueSet
#SBATCH -N 1
#SBATCH --mem=120GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00

# sbatch run_*.slurm randomize_samples.txt

# Creates the Manifest data

# for file in $(ls *set1.fq); do basename=${file##*/}; base="${basename%*_1.fastq}"; bs="${base}"; fwr=${base}_1.fastq; rev=${base}_2.fastq; echo "$bs" `printf "$PWD/$fwr"`; done > ${sf}_randomize_samples.txt

CPU=$SLURM_NPROCS
MEM=$SLURM_MEM_PER_NODE

Manifest=$1 # samples.txt



single_reads=${Manifest%.txt}.concat.fq


call=`head -n1 $Manifest | awk '{print $2}' |tr "\n" " "`

cat $call > $single_reads


# creates a while to run individual chkpoints per assembly tool

spades="rnaspades.py -s $single_reads -t $CPU -m $MEM -o $OUTDIR"
trinity="Trinity --seqType fq --max_memory 120G --single $single_reads --no_normalize_reads --CPU $CPU --output $OUTDIR"
transabbys="transabbys $single_reads --threads $CPU  --outdir $OUTDIR"
BayesDenovo="BayesDenovo_v1.pl --seqType fq --single $single_reads --CPU $CPU  --output $OUTDIR --bamfile $OUTDIR/bamfile.bam --fafile $OUTDIR/bayesdenovo.fa"
SOAPdenovo="SOAPdenovo-Trans all -s $configureFile -p $CPU  -o $OUTDIR"


tools=("spades" "trinity" "transabbys" "BayesDenovo")
i=0

while [ $i -lt ${#tools[@]} ]; do
    tool=${tools[$i]}
    OUTDIR=${Manifest%.txt}_${tool}_dir
    mkdir -p $OUTDIR

    if [ ! -f "1_${tool}.chkp" ]; then
        call=${!tool}
        eval $call
        touch 1_${tool}.chkp
    else
        echo "${tool} checkpoint already exists."
    fi

    i=$((i + 1))
done











