#!/bin/bash
#SBATCH --job-name=sampling_analysis
#SBATCH -N 1
#SBATCH --mem=120GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00
#SBATCH -o out.%J
#SBATCH -e err.%J

#manifest="$1"

for P in $(seq 0.1 0.1 1); do
    echo "Processing proportion: $P ======================"
     ./Subsampling.sh -s all -p $P
done

exit

# Split replicates proportion (not used yet)

for i in $(seq 1 4); do
    echo "running round ${i}"
    ./Subsampling.sh -s $manifest -o boostrap_dir_${i}
done

for i in $(seq 5 8); do
    echo "running round ${i}"
    ./Subsampling.sh -s $manifest -o boostrap_dir_${i}
done

for i in $(seq 9 12); do
    echo "running round ${i}"
    ./Subsampling.sh -s $manifest -o boostrap_dir_${i}
done

exit

find "boostrap_dir_${i}/*_FASTA_DIR" -name "*.fa" | while read -r contig_file; do

            call="Accuracy.sh -d $contig_file -s $manifest"
            
            echo $call
            
            eval $call

        done

exit