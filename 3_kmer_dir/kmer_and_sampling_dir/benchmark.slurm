#!/bin/bash
#SBATCH --job-name=benchmark
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00
#SBATCH -o out.%J
#SBATCH -e err.%J

# Import seqkit

EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/
export PATH=$PATH:$EXPORT

for k in 19 21 25 33 38 44 49 55 60 73; do
    for p in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.8 0.9 1.0; do
        call="./benchmark.sh -s all -m both -k $k -p $p -o factorial_dir"
        echo $call
        eval $call
    done
done

exit