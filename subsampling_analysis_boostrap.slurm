#!/bin/bash
#SBATCH --job-name=boostrap_analysis
#SBATCH -N 1
#SBATCH --mem=120GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00
#SBATCH -o out.%J
#SBATCH -e err.%J

manifest="$1"

# Split boostraping proportion

for i in $(seq 1 4); do
    echo "running round ${i}"
    ./subsampling_analysis.sh -s $manifest -o boostrap_dir_${i}
done

for i in $(seq 5 8); do
    echo "running round ${i}"
    ./subsampling_analysis.sh -s $manifest -o boostrap_dir_${i}
done

for i in $(seq 9 12); do
    echo "running round ${i}"
    ./subsampling_analysis.sh -s $manifest -o boostrap_dir_${i}
done

exit