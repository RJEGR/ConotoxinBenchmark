```bash
git clone https://github.com/soedinglab/hh-suite.git
mkdir -p hh-suite/build && cd hh-suite/build
cmake -DCMAKE_INSTALL_PREFIX=. ..
make -j 4 && make install
export PATH="$(pwd)/bin:$(pwd)/scripts:$PATH"
```

Tools
```bash
EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/hhsuite-3.3.0
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/mmseqs/bin
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/mafft-linux64/
export PATH=$PATH:$EXPORT

```

# Generating a multiple sequence alignment using mmseq or maft

## Local alignment of prefiltered sequence pairs using mmseqs align

https://mmseqs.com/latest/userguide.pdf#page=50.72

```bash
CPU=20 #$SLURM_NPROCS
MEM=100 #$SLURM_MEM_PER_NODE
INPUT=O1_superfamily_californicus.fasta

mmseqs align -h

# mmseqs align <i:queryDB> <i:targetDB> <i:resultDB> <o:alignmentDB> [options]
mmseqs createdb $INPUT sequenceDB 

# mmseqs prefilter <i:queryDB> <i:targetDB> <o:prefilterDB>
#  -s FLOAT           Sensitivity: 1.0 faster; 4.0 fast; 7.5 sensitive [4.000]
# using 7 allow sequences passed prefiltering 
mmseqs prefilter sequenceDB sequenceDB resultDB_pref -s 7 --threads $CPU

# Create a HHblits database:
# mmseqs result2msa <i:queryDB> <i:targetDB> <i:resultDB> <o:msaDB> [options]

mmseqs result2msa sequenceDB sequenceDB resultDB_pref MSA --msa-format-mode 2
mmseqs result2msa sequenceDB sequenceDB resultDB_pref MSA --msa-format-mode 5 --threads $CPU

--msa-format-mode INT        Format MSA as: 0: binary cA3M DB
                              1: binary ca3m w. consensus DB
                              2: aligned FASTA DB
                              3: aligned FASTA w. header summary
                              4: STOCKHOLM flat file
                              5: A3M format
                              6: A3M format w. alignment info [2]

# OR

# mmseqs align <i:queryDB> <i:targetDB> <i:resultDB> <o:alignmentDB>
# mmseqs align sequenceDB sequenceDB resultDB_pref resultDB_aln --thread $CPU

```


```bash
mafft --localpair --maxiterate 1000 $INPUT > ${INPUT%.fasta}.aln

# L-INS-i specializes in generating accurate alignments by:
# Iterative refinement: Repeatedly optimizing the alignment using consistency-based scoring derived from all pairwise comparisons.
# Local alignment focus: Prioritizing conserved local segments rather than forcing full-length matches, making it ideal for sequences with variable domains or non-homologous regions.

```

# Generating a multiple sequence alignment using HHblits
Continue here 
https://github.com/soedinglab/hh-suite/wiki#hhsearchhhblits-model-format-hhm-format

```bash
hhmake -i MSA # issues w/ mmseq format

# Now you can generate a hidden Markov model (HMM) from this MSA:

hhmake -i ${INPUT%.fasta}.aln -M first # this works, see output hhm

hhblits -cpu 4 -i $INPUT -d O1_superfamily_californicus.hhm -oa3m query.a3m -n 1
```


Populate using HHblits
Multiple alignments can be read in A2M, A3M, or aligned FASTA format. 

```bash
hhblits -h all 
# -i input/query: single sequence or multiple sequence alignment (MSA) in a3m, a2m, or FASTA format, or HMM in hhm format

hhblits -cpu 4 -i data/query.seq -d databases/uniclust30_2018_08/uniclust30_2018_08 -oa3m query.a3m -n 1

```
Using HMMs both on the query and the database side greatly enhances the sensitivity/selectivity and alignment quality over sequence-profile based methods such as PSI-BLAST. HHsearch is the first software to employ HMM-HMM comparison and HHblits is the first profile-profile comparison method that is fast enough to do iterative searches to build MSAs.