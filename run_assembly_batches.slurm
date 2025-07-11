#!/bin/bash
#SBATCH --job-name=metaanalysis
#SBATCH -N 1
#SBATCH --mem=120GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00
#SBATCH -o out.%J
#SBATCH -e err.%J


if [[ "$1" == "-h" ]]; then
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
CONFIG=$2 

if [[ -z "$Manifest" || -z "$CONFIG" ]]; then
    echo "Error: Missing required arguments."
    echo "Run '$(basename $0) --help' for usage information."
    exit 1
fi

forward_fq=${Manifest%.*}_concat_PE1.fq
reverse_fq=${Manifest%.*}_concat_PE2.fq

call=`awk '{print $2}' $Manifest | tr "\n" " "`
cat $call > $forward_fq

call=`awk '{print $3}' $Manifest | tr "\n" " "`
cat $call > $reverse_fq

CPU=$SLURM_NPROCS
MEM=$SLURM_MEM_PER_NODE


chmod +x run_assemblers_dir/*sh

FASTA_DIR="${Manifest%.*}_FASTA_DIR"

mkdir -p $FASTA_DIR

while read -r line; do 
    # Skip lines starting with '#'
    [[ "$line" =~ ^[[:space:]]*# ]] && continue
    [[ -z "${line// }" ]] && continue

    tool=$(echo $line | awk -F "=" '{print $1}')
    
    bs=${Manifest%.*}

    call=$(echo $line | awk -F "=" '{print $2}') 

    chkp_file="1_${tool}_${bs}.chkp"   
   
    if [ ! -f "$chkp_file" ]; then

        OUTDIR=${bs}_${tool}_dir
        
        mkdir -p $OUTDIR

        call=$(eval echo $call $forward_fq $reverse_fq  $OUTDIR $CPU $MEM)         

        echo "Executing: $call"

        eval $call

        touch $chkp_file

        rm -fr $OUTDIR

    else
        echo "${chkp_file} checkpoint already exists."
    fi
done < $CONFIG

# Continue with the next steps: transrate

REFDIR=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/fernando_pub/reads_artificiales/inputs

#Trinity/contigs.csv

mkdir -p 2_transrate_dir

chmod +x 2_transrate.sh


find "$FASTA_DIR" -maxdepth 1 -type f -name '*.fa' | while read -r ASSEMBLER; do

    FAMILYPREFIX=$(basename ${FASTA_DIR%_FASTA_DIR})
    REF=$REFDIR/${FAMILYPREFIX}.fasta
    BSCONTIG=$(basename ${ASSEMBLER%.fa})
    TRANSRATE_DIR=${BSCONTIG}_transrate_dir

    call="./2_transrate.sh $forward_fq $reverse_fq $ASSEMBLER $REF $TRANSRATE_DIR"

    #call="./2_transrate.sh $ASSEMBLER $REF $TRANSRATE_DIR"
    
    eval $call
done

rm ${Manifest%.*}_concat_PE*.fq 


exit
