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

single_reads=`awk '{print "--s" NR, $2}' $Manifest`

OUTDIR=${Manifest%.txt}_spades_dir

mkdir -p $OUTDIR


if [ ! -f "1_spades.chkp" ]; then

    call="rnaspades.py $single_reads -t $CPU -m $MEM -o $OUTDIR"

    eval $call

    touch 1_spades.chkp

else
    echo "..."
fi

# Trinity
module load trinityrnaseq-v2.15.1

OUTDIR=${Manifest%.txt}_trinity_dir

mkdir -p $OUTDIR

single_reads=`awk '{print $2}' $Manifest` # comma-delimited list of single-end reads


if [ ! -f "1_trinity.chkp" ]; then

    #call="Trinity --seqType fq --max_memory 120G --samples_file $Manifest --no_normalize_reads --CPU $CPU --output $OUTDIR"
    call="Trinity --seqType fq --max_memory 120G --single $single_reads --no_normalize_reads --CPU $CPU --output $OUTDIR"
    
    eval $call

    touch 1_trinity.chkp

else
    echo "..."
fi

# SOAPDenovo-trans


module load conda-2024_py3.8
source activate soapdenovo-trans


configureFile=`awk '{print "q="""$2}' $Manifest`


OUTDIR=${Manifest%.txt}_SOAPdenovo_dir

mkdir -p $OUTDIR


if [ ! -f "1_SOAPDenovo.chkp" ]; then

    call="SOAPdenovo-Trans all -s $configureFile -p $CPU  -o $OUTDIR"
    
    eval $call

    touch 1_SOAPDenovo.chkp

else
    echo "..."
fi

# transabbys

module load conda-2025
conda activate transabbys

OUTDIR=${Manifest%.txt}_transabbys_dir

mkdir -p $OUTDIR

single_reads=`awk '{print "--se" NR, $2}' $Manifest`

if [ ! -f "1_transabbys.chkp" ]; then

    call="transabbys $single_reads --threads $CPU  --outdir $OUTDIR"
    
    eval $call

    touch 1_transabbys.chkp

else
    echo "..."
fi

# idba

export PATH=/LUSTRE/apps/bioinformatica/idba/bin:$PATH


if [ ! -f "1_idba.chkp" ]; then

    call="Code $CPU  $OUTDIR"
    
    # eval $call

    #touch 1_idba.chkp

else
    echo "..."
fi

# rnabloom

rnabloom
module load conda-2025
conda activate rnabloom

OUTDIR=${Manifest%.txt}_rnabloom_dir

mkdir -p $OUTDIR


if [ ! -f "1_rnabloom.chkp" ]; then

    # call="java -jar RNA-Bloom.jar -sef $single_end -t $CPU  -outdir $OUTDIR"
    
    eval $call

    touch 1_rnabloom.chkp

else
    echo "..."
fi

# cap3

cap3
module load conda-2025
conda activate cap3 

if [ ! -f "1_cap3.chkp" ]; then

    call="Code $CPU  $OUTDIR"
    
    # eval $call

    # touch 1_cap3.chkp

else
    echo "..."
fi

# Bayesdenovo

export PATH=/LUSTRE/apps/bioinformatica/BayesDenovo/BayesDenovo-v1:$PATH

OUTDIR=${Manifest%.txt}_BayesDenovo_dir

mkdir -p $OUTDIR

if [ ! -f "1_Bayesdenovo.chkp" ]; then

    call="BayesDenovo_v1.pl --seqType fq --single $single_end --CPU $CPU  --output $OUTDIR --bamfile $OUTDIR/bamfile.bam --fafile $OUTDIR/bayesdenovo.fa"
    
    eval $call

    touch 1_Bayesdenovo.chkp

else
    echo "..."
fi











