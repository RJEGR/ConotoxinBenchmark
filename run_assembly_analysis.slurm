#!/bin/bash
#SBATCH --job-name=metaanalysis
#SBATCH -N 1
#SBATCH --mem=120GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00

# Check for --help flag
if [[ "$1" == "--help" ]]; then
    echo "Usage: $(basename $0) <Manifest> <Config>"
    echo
    echo "Arguments:"
    echo "  Manifest: A text file containing sample information. The first column should contain sample names, and the second column should contain file paths."
    echo "  Config:   A configuration file where each line specifies a tool and its command in the format 'tool=command'."
    echo
    echo "Description:"
    echo "  This script performs meta-analysis by concatenating reads from the Manifest file, configuring tools from the Config file, and executing the specified commands."
    exit 0
fi

Manifest=$1 
Config=$2 

if [[ -z "$Manifest" || -z "$Config" ]]; then
    echo "Error: Missing required arguments."
    echo "Run '$(basename $0) --help' for usage information."
    exit 1
fi

single_reads=${Manifest%.txt}.concat.fq

call=`head -n1 $Manifest | awk '{print $2}' |tr "\n" " "`
# call=`awk '{print $2}' $Manifest | tr "\n" " "`

cat $call > $single_reads


CPU=$SLURM_NPROCS
MEM=$SLURM_MEM_PER_NODE


# First run stand-alone tools

export PATH=/LUSTRE/apps/bioinformatica/SPAdes-3.15.5-Linux/bin:$PATH
export PATH=/LUSTRE/apps/bioinformatica/trinityrnaseq-v2.15.1:$PATH
export PATH=/LUSTRE/apps/bioinformatica/idba/bin:$PATH
export PATH=/LUSTRE/apps/bioinformatica/BayesDenovo/BayesDenovo-v1:$PATH

module load conda-2025
export PATH=/LUSTRE/apps/Anaconda/2025/miniconda3/envs/rnabloom/bin/:$PATH


while read -r line; do 
    
    tool=$(echo $line | awk -F "=" '{print $1}')
    
    call=$(echo $line | awk -F "=" '{print $2}') 

    chkp_file="1_${tool}_${Manifest%.txt}.chkp"   
   
   if [ ! -f "$chkp_file" ]; then

        OUTDIR=${Manifest%.txt}_${tool}_dir

        mkdir -p $OUTDIR

        call=$(eval echo $call) # Expands variables in the command to fix issue with quotes         

        echo "Executing: $call"

        eval $call

        touch $chkp_file

        
    else
        echo "${tool} checkpoint already exists."
    fi
done < $Config


exit

# Now run the tools that require a configuration file or python dependencies

# Tools working under conda-2024

module load conda-2024_py3.8

export PATH=/LUSTRE/apps/Anaconda/2024/miniconda3-py3.8/envs/soapdenovo-trans/bin/:$PATH

# Tools working under conda-2025

module load conda-2025
# transabbys use samtools lower version than trinity, and trinity crash, therefore does not export PATHS in parallel

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

# rm *.concat.fq