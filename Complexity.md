# Title

Task: Determine sequence complexity of a given superfamily by computing some sequence complexity measures such as:
- nonsynonymous to synonymous substitution ratio (ka/ks): understand the evolutionary forces shaping gene sequences, revealing whether a gene is under functional constraint or adaptation over time (see discusion at [Zheng et al., 2023](https://doi.org/10.1186/s12864-023-09689-4)). To compute Ka/ks use [orthofinder](https://github.com/OrthoFinder/OrthoFinder) and [ParaAT](https://github.com/wonaya/ParaAT/tree/master) in a pipeline, using FASTA by superfamily as intial input. See introductin by at [Álvarez-Carretero et al., 2023](https://academic.oup.com/mbe/article/40/4/msad041/7140562)

- Phylogeny mining (estimate nucleotide diversity): ...

- information and entropy measures (information content and entropy): quantify the amount of information contained in a sequence, and the degree of uncertainty in the sequence. Here we assume repetitiviness and hypervariability of conotoxins contributes to changes in this measures. Some tools are reviewed by [Orlov and Orlova 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10643780/). DNA sequence complexity as well as protein sequence complexity methods could be broadly classified into large groups—entropy-based and compression-based methods. Entropy based methods of complexity estimates include word frequency (linguistic) approaches, Shannon entropy and its variants. Compression based methods include modifications of Lempel–Ziv compression scheme, could be applied for repeat search in genomes (direct, inverted, generated) and alignment-free genome comparisons (Orlov and Orlova 2023)


- Estimates simpler complexity, using the formula `(base[i] != base[i+1])`, for a sequence of 51-bp, with 3 bases that is different from its next base `AAAATTTTTTTTTTTTTTTTTTTTTGGGGGGGGGGGGGGGGGGGGGGCCCC`, its complexity is `complexity = 4/(51-1) = 6%`. ([Shifu Chen 2023](https://github.com/OpenGene/fastp?tab=readme-ov-file#low-complexity-filter))


## Orthology
After split and output sequences by superfamily, 
```bash
cd /Users/cigom/Documents/GitHub/ConotoxinBenchmark/INPUTS/
orthofinder -S diamond -op -f peptide_dir
orthofinder -d -S diamond -op -f orthofinder_dir
```

Then, test ParaAT and kaKs calculator
```bash
EXPORT=/Users/cigom/Documents/Tools/ParaAT-master
export PATH=$PATH:$EXPORT

EXPORT=/Users/cigom/Documents/Tools/kakscalculator2/bin
export PATH=$PATH:$EXPORT

# generate axt file for KaKs_Calculator
ParaAT.pl -h test.homologs -n test.cds -a test.pep -p proc -o output_axt  -f axt

KaKs_Calculator -i output_axt/NP_000039-NP_598529.cds_aln.axt -o example.axt.kaks -m YN

# ----

ParaAT.pl -h homologs_from_orthofinder -n test.cds -a test.pep -p proc -o output_axt  -f axt
KaKs_Calculator -i output_axt/NP_000039-NP_598529.cds_aln.axt -o example.axt.kaks -m YN

```


Alternatively We can run in r
```r
library(seqinr)
aligned_sequences <- read.alignment("your_alignment_file.fasta", format = "fasta")
result <- kaks(aligned_sequences)
print(result$ka)  # Ka values
print(result$ks)  # Ks values
print(result$ka / result$ks)  # Ka/Ks ratios

# or 

```
[CRBHits](https://github.com/kullrich/CRBHits) R package (https://doi.org/10.21105/joss.02424) is designed to directly handle orthologous coding sequences (from RBH and cluster methods like OrthoFinder), create codon alignments, and calculate Ka/Ks within R using KaKs_Calculator as a backend. This can be a one-stop solution to:

```bash
library(CRBHits)

# Assuming you have CDS files for species 1 and species 2
# and list of orthologous pairs from OrthoFinder

# Convert CDS files to reciprocal best hits with filtering and longest isoform extraction
rbh_pairs <- cdsfile2rbh("species1_cds.fa", "species2_cds.fa", longest.isoform = TRUE, threads = 8)

# Calculate Ka/Ks for these pairs (can use KaKs_Calculator backend)
kaks_results <- rbh2kaks(rbh_pairs)

# View results
head(kaks_results)

```

## Computes Smith-Waterman alignment scores between all sequences
mmseqs clust: Computes a similarity clustering of a sequence database based on Smith
Waterman alignment scores of the sequence pairs computed by mmseqs align.

mmseqs align: Computes Smith-Waterman alignment scores between all sequences in the
query database and the sequences of the target database whose prefiltering scores computed
by mmseqs prefilter pass a minimum threshold.

```bash
cd ~/Documents/GitHub/ConotoxinBenchmark/INPUTS/

RELEASE="${3:-$(date "+%Y_%m")}"
INPUT=conoServerDB.fasta
OUTDIR=`pwd`/$RELEASE
OUTPUT=${INPUT%.*}_SWscores.tsv

TMPDIR=$OUTDIR/tmp
mkdir -p $TMPDIR

#OUTDIR=$(abspath $OUTDIR)
#TMPDIR=$(abspath $TMPDIR)

SEQUENCE_DB="${OUTDIR}/${INPUT%.*}_db"

mmseqs createdb "$INPUT" "${SEQUENCE_DB}"
INPUT="${SEQUENCE_DB}"
mmseqs clusthash $INPUT "$TMPDIR/aln_redundancy" --min-seq-id 1

date ${DATE}
mmseqs clust $INPUT "$TMPDIR/aln_redundancy" "$TMPDIR/clu_redundancy" ${CLUSTER1_PAR}
mmseqs createsubdb "$TMPDIR/clu_redundancy" $INPUT "$TMPDIR/input_step0"

STEP=0
INPUT="$TMPDIR/input_step0" 

mmseqs prefilter "$INPUT" "$INPUT" "$TMPDIR/pref_step$STEP"

ALIGNMENT_COMMON="$COMMON -e 0.001 --max-seq-len 65535 --max-rejected 2147483647"
ALIGNMENT0_PAR="--alignment-mode 0 --min-seq-id 0 --comp-bias-corr 0 ${ALIGNMENT_COMMON}"

# ${ALIGNMENT0_PAR}

mmseqs align "$INPUT" "$INPUT" "$TMPDIR/pref_step$STEP" "$TMPDIR/aln_step$STEP" 
mmseqs convertalis $INPUT $INPUT "$TMPDIR/aln_step$STEP" $OUTPUT


# The file is formatted as a tab-separated list with 12 columns: (1,2) identifiers for query and target
sequences/profiles, (3) sequence identity, (4) alignment length, (5) number of mismatches, (6) number of gap openings, (7-8, 9-10) domain start and end-position in query and in target, (11) E-value,
and (12) bit score.
```