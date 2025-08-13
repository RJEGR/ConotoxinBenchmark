
```bash
EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/hmmer-3.4/bin
export PATH=$PATH:$EXPORT
EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/mafft-linux64/
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/
export PATH=$PATH:$EXPORT

```

```bash
WD=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/fernando_pub/hhsuite_model_dir/californicus_dir/MSA_dir
cd $WD
```

```bash
# cat all.pep.fa araneomorphae_all.pep.fa mygalomorphae_all.pep.fa | seqkit grep -n -r -p "sodium" > sodium_channel.pep

INPUT=O1_superfamily_californicus.fasta
mafft --localpair --maxiterate 1000 $INPUT > ${INPUT%.fasta}.aln
```

```bash
hmmbuild ${INPUT%.fasta}_hmmfile.out ${INPUT%.fasta}.aln
```

```bash
cat hmmfile1.out hmmfile2.out > combined.hmm

```

```bash
HMM=combined.hmm
INPUT=query
hmmsearch --tblout ${INPUT%.fasta}.tblout --cpu $CPU $HMM $INPUT
```



