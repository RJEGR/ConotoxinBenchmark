#!/bin/bash
#SBATCH --job-name=metaanalysis
#SBATCH -N 1
#SBATCH --mem=120GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00

CPU=$SLURM_NPROCS
MEM=$SLURM_MEM_PER_NODE

Manifest=$1 # samples.txt

Config=$2 # Configure file


single_reads=${Manifest%.txt}.concat.fq


call=`head -n1 $Manifest | awk '{print $2}' |tr "\n" " "`

cat $call > $single_reads

#Exporting tools

export PATH=/LUSTRE/apps/bioinformatica/SPAdes-3.15.5-Linux/bin:$PATH
export PATH=/LUSTRE/apps/bioinformatica/idba/bin:$PATH
export PATH=/LUSTRE/apps/bioinformatica/BayesDenovo/BayesDenovo-v1:$PATH


# Tools working under conda-2024

conda deactivate


module load conda-2024_py3.8

export PATH=/LUSTRE/apps/Anaconda/2024/miniconda3-py3.8/envs/soapdenovo-trans/bin/:$PATH

# Tools working under conda-2025

module load conda-2025

export PATH=/LUSTRE/apps/Anaconda/2025/miniconda3/envs/rnabloom/bin/:$PATH
export PATH=/LUSTRE/apps/Anaconda/2025/miniconda3/envs/transabbys/bin/:$PATH

# export PATH=/LUSTRE/apps/Anaconda/2025/miniconda3/envs/cap3/bin/:$PATH

# cap3 needs FASTA FORMAT
# CAP3="cap3 $CPU  $OUTDIR"

configureFile=${Manifest%.txt}_SOAPdenovo.config

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
q=$single_reads 
EOF

while read -r line; do 
    
    tool=$(echo $line | awk -F "=" '{print $1}')
    
    call=$(echo $line | awk -F "=" '{print $2}')     
   
   if [ ! -f "1_${tool}.chkp" ]; then

        OUTDIR=${Manifest%.txt}_${tool}_dir

        echo $OUTDIR

        echo $call

        touch 1_${tool}.chkp

        
    else
        echo "${tool} checkpoint already exists."
    fi
done < $Config

