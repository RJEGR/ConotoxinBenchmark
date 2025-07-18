
EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/hhsuite-3.3.0
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/mmseqs/bin
export PATH=$PATH:$EXPORT

WD=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/DataBases/Conopeptides/Uniclust_dir

INPUT=$WD/Uniclust30_2018_08.fasta
sequenceDB=$WD/uniclust30_2018_08

# uniclust_workflow.sh $INPUT $sequenceDB

# 2. After running uniclust_workflow.sh, creates a database of HHblits compatible MSAs
# For creating an HHblits database from a clustering, the procedure is almost the same, except that
# you have to create symlinks to the ﬀindex _header and _sequence files needed by HHblits:



# Generating a multiple sequence alignment using HHblits

# -i input/query: single sequence or multiple sequence alignment (MSA) in a3m, a2m, or FASTA format, or HMM in hhm format

hhblits -cpu 4 -i data/query.seq -d databases/uniclust30_2018_08/uniclust30_2018_08 -oa3m query.a3m -n 1

hhblits -cpu 4 -i uniprotkb_taxonomy_id_33208_AND_cc_tiss_2025_02_19.fasta -d outDir/2025_07/uniclust30_2025_07 -oa3m query.a3m -n 1

