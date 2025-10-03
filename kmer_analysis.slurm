#!/bin/bash
#SBATCH --job-name=kmer_screening
#SBATCH -N 1
#SBATCH --mem=120GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00
#SBATCH -o out.%J
#SBATCH -e err.%J

#Larger k-mers (like 77, 99, and 127) are more unique and are crucial for resolving repetitive regions in the genome, leading to a more contiguous assembly.

# Smaller k-mers (like 21 and 33) are essential for capturing regions of low coverage and ensuring the graph is well-connected, but they perform poorly in repetitive regions.

# Auto mode: SPAdes checks the maximum read length and chooses a pre-defined set of k-mers that are generally well-suited for that data type.
# If max read length is less than 100 bp:

# SPAdes will use: -k 21,33,55

# common for   k: [19, 21, 25, 33, 49, 55, 73]

for k in 19 21 25 33 38 44 49 55 60 73; do
    echo "Processing k: $k ======================"
    ./kmer.sh -s all -k $k
done

# small k-mers
for k in 19,25,33; do
    echo "Processing k: $k ======================"
    #./kmer.sh -s all -k $k
done

exit
