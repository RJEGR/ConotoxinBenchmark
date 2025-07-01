#!/bin/bash
#SBATCH --job-name=TrueSet
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

module load conda-2024_py3.8
source activate soapdenovo-trans

# Tools working under conda-2025

module load conda-2025
conda activate transabbys
conda activate rnabloom
conda activate cap3 


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

