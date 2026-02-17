#!/bin/bash
#SBATCH --job-name=ConoSorter
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00
#SBATCH -o out.%J
#SBATCH -e err.%J

mkdir -p ConoSorter_outdir

for i in $(ls *.fa); do
  echo "Processing file: $i"
  
  ConoSorter_nuc_pipe.sh $i  

  mv ${i}_Regex.tab ConoSorter_outdir
  mv ${i}_pHMM.tab ConoSorter_outdir
done

exit


# Transdecoder-to-ConoSorter pipeline

mkdir -p ConoSorter_outdir

for i in $(ls *.fa); do
  
  echo "Step 1. Processing file: $i"

  ./TransDecoder_predict.sh $i

  pep_file="${i}.transdecoder.pep"

  echo "Step 2. Processing file: $pep_file"

  ./ConoSorter_pep_pipe.sh $pep_file

  mv ${i%.pep}_Regex.tab ConoSorter_outdir

  mv ${i%.pep}_pHMM.tab ConoSorter_outdir

  rm *.chkp eggnog_mapper.emapper.* *.outfmt6 *.ok *.bed *.gff

done

exit


