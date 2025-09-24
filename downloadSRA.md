Ricardo Gómez-Reyes

### Tools required

- Esearch Installation: [Here](https://www.ncbi.nlm.nih.gov/books/NBK179288/)

- SRA Tool Kit Installation (For prefetch and fastq-dump) : [Here](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit)
- (Optional) SetKit: [Here](https://bioinf.shenwei.me/seqkit/download/)

### 1) Get all SRA runs for a BioProject based on an SRA Run ID

```bash
# PRJNA450372, Single Cell RNA sequencing of Adult Human Breast Epithelial Cells
esearch -db sra -query "PRJNA450372" |  efetch -format docsum | xtract -pattern Runs -ACC @acc  -element "&ACC" > SraAccList.txt
```

### 2) A list of Runs:

Proof-read: [Here](https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump#how-to-use-prefetch-and-fasterq-dump-to-extract-fastq-files-from-sra-accessions)

```bash
# Test
head -n1 SraAccList.txt | xargs -n 1 -P 12 fastq-dump --split-files --gzip --skip-technical --outdir .
# run All downloads (Slurm recommended)
fastq-dump --split-files SraAccList.txt --gzip --skip-technical --outdir .
```

### 3) (Optional)

```bash
seqkit stat *fastq.gz

#file                   format  type   num_seqs      sum_len  min_len  avg_len  max_len
#SRR7008752_1.fastq.gz  FASTQ   DNA   1,859,050  185,905,000      100      100      100
# SRR7008752_2.fastq.gz  FASTQ   DNA   1,859,050  185,905,000      100      100      100
```

### Processing

Samples save at:

/LUSTRE/bioinformatica_data/ecologia_molecular/Francesco/Cancer

Los datos obtenidos mediante fastq-dump, incluía en SRAToolkit; se encuentra en la trayectoria: /LUSTRE/bioinformatica_data/ecologia_molecular/SRAdata

Se descargo (Silvia) unicamente los datos del proyecto: 

- PRJNA450372 (590 librerias con prefijo SRR700*, corresponde a las bibliotecas control)
- pero hace falta descargar las bibliotecas del proyecto de cancer (PRJNA625321), prefilo SRR115.

```bash
# PRJNA450372, Single Cell RNA sequencing of Adult Human Breast Epithelial Cells
esearch -db sra -query "PRJNA625321" |  efetch -format docsum | xtract -pattern Runs -ACC @acc  -element "&ACC" > SraAccList.txt

# 

./fastq-dump -n 1 -P 1 --split-files SraAccList.txt --gzip --skip-technical --outdir .


# Test
head -n1 SraAccList.txt | xargs -n 1 -P 1 ./fastq-dump --split-files --gzip --skip-technical --outdir .

# Then
# for i in $(cat *txt); do echo $i; done

for i in $(cat *txt); do echo $i | xargs -n 1 -P 1 ./fastq-dump --split-files --gzip --skip-technical --outdir .; done

# track the download:
ls -lth

# (optional)
seqkit stat *fastq.gz
```