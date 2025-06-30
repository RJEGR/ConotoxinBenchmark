
```bash
EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/hmmer-3.4/bin
export PATH=$PATH:$EXPORT
```

```bash
WD=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/fernando_pub/hhsuite_model_dir/californicus_dir/MSA_dir
cd $WD
```

```bash
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
hmmsearch --tblout ${INPUT%.fasta}.tblout --cpu $CPU ${INPUT%.fasta}_hmmfile.out $INPUT
```