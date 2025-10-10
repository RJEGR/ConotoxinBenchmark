#!/bin/bash
## Directivas
#SBATCH --job-name=BUSCOpy
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=20
#SBATCH -t 6-00:00:00

# module load python-2.7-anaconda

module load gcc-12.2.0


IN=$1 # Trinity.longest.fasta
out=${IN%.*}

EXPORT=/home/rgomez/bin/busco-master/scripts/
export PATH=$EXPORT:$PATH

EXPORT=/LUSTRE/apps/bioinformatica/ncbi-blast-2.14.0+/bin/
export PATH=$EXPORT:$PATH

EXPORT=/LUSTRE/apps/bioinformatica/hmmer-3.3.2/bin/
export PATH=$EXPORT:$PATH

BUSCOdb=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/BUSCOdb


eukaryota_odb10=$BUSCOdb/eukaryota_odb10
mollusca_odb10=$BUSCOdb/mollusca_odb10
metazoa_odb10=$BUSCOdb/metazoa_odb10
mammalia_odb10=$BUSCOdb/mammalia_odb10
bacteria_odb10=$BUSCOdb/bacteria_odb10

CPU=$SLURM_NPROCS

run_BUSCO.py -i $IN -l $mollusca_odb10 -m transcriptome -o ${out}_mollusca_odb10 -c $CPU

run_BUSCO.py -i $IN -l $metazoa_odb10 -m transcriptome -o ${out}_metazoa_odb10 -c $CPU

run_BUSCO.py -i $IN -l $eukaryota_odb10 -m transcriptome -o ${out}_eukaryota_odb10 -c $CPU

run_BUSCO.py -i $IN -l $bacteria_odb10 -m transcriptome -o ${out}_bacteria_odb10 -c $CPU

run_BUSCO.py -i $IN -l $mammalia_odb10 -m transcriptome -o ${out}_mammalia_odb10 -c $CPU

mkdir -p full_table_busco_dir/

cp ./full_table*_odb10/*.tsv full_table_busco_dir


exit



mkdir -p summaries

cp ./*_odb10/short_summary* summaries

module load R-4.3.1

generate_plot.py --working_directory ./summaries/

cp summaries/busco_figure.R .
sed -i 's/_odb10//g' summaries/busco_figure.R
Rscript summaries/busco_figure.R
# firefox summaries/busco_figure.png
cp summaries/busco_figure.R .