#!/bin/bash
#SBATCH --job-name=boostrap_analysis
#SBATCH -N 1
#SBATCH --mem=120GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00
#SBATCH -o out.%J
#SBATCH -e err.%J

manifest="$1"

for i in {1..10}; do
    echo "running round ${i}"
    ./subsampling_analysis.sh -s $manifest -o boostrap_dir_${i}
done

exit