
write a code in python for the follow instruction: Functions

1. read_fasta(file_path)
Reads a (possibly gzipped) FASTA file and returns a list of sequence records.

2. randomize_sequences(sequences, sample_size, kmer_length)
from the every provided sequences randomly (sample_size) selects k-mers (subsequences of length k) of length defined in the input flag kmer_length and sample_size. 

Only includes sequences long enough for the specified by the flag kmer_length.

Each k-mer is randomly chosen from a random position in a random sequence.

For the provided sequences estimates the sampling effort of the subsequences to coverage the 100% of sequence by the sample_size and kmer_length. 

3. write_fasta(sequences, output_path)
By default Writes a list of subsequences to a file in FASTA format.

4. write_fastq(sequences, output_path, phred_score=30)
if option --output_format is fastq Writes a list of subsequences to a file in FASTQ format. 

For this function a phred_score between 10 to 40 may be input. If phred_score is not provided, the function assigns a uniform PHRED quality score of 40 to every base of each subsequences.

5. When run as a script, it uses argparse to provide command-line options for:

Input FASTA file (possibly gzipped), output filename prefix, 
Set the number of random subsequences to sampling (--sample_size)
Set the length of subsequences to (--kmer_length)
Output type (FASTA by default or FASTQ if --output_format is fastq).
PHRED score for FASTQ (optional).

consider use Bio random, and gzip utilities. including SeqIO, Seq and SeqRecord
