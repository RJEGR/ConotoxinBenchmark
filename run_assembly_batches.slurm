#!/bin/bash
#SBATCH --job-name=Assemblers
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00
#SBATCH -o out.%J
#SBATCH -e err.%J

Assemblers.sh -s all -c run.config

exit 0