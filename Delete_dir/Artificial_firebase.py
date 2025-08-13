import argparse
import gzip
import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def read_fasta(file_path):
 """Reads a (possibly gzipped) FASTA file and returns a list of sequence records."""
    if file_path.endswith('.gz'):
        with gzip.open(file_path, "rt") as handle:
 return list(SeqIO.parse(handle, "fasta"))
 else:
        with open(file_path, "r") as handle:
 return list(SeqIO.parse(handle, "fasta"))

def randomize_sequences(sequences, sample_size, kmer_length):
 """
 From the provided sequences, randomly selects sample_size k-mers (subsequences of length kmer_length).
 Only includes sequences long enough for the specified kmer_length.
 Each k-mer is randomly chosen from a random position in a random sequence.
 Estimates the sampling effort to cover 100% of the sequence with the given sample size and k-mer length.
 """
    valid_sequences = [seq for seq in sequences if len(seq.seq) >= kmer_length]
    randomized_kmers = []
    total_sequence_length = sum(len(seq.seq) for seq in valid_sequences)

    if total_sequence_length == 0:
 print("No sequences long enough for the specified k-mer length.")
 return [], 0

    for _ in range(sample_size):
 # Randomly select a sequence
        random_seq_record = random.choice(valid_sequences)
 random_seq = random_seq_record.seq

 # Randomly select a starting position for the k-mer
        start_pos = random.randint(0, len(random_seq) - kmer_length)

 # Extract the k-mer
        kmer_seq = random_seq[start_pos : start_pos + kmer_length]

 # Create a SeqRecord for the k-mer
 kmer_record = SeqRecord(kmer_seq, id=f"{random_seq_record.id}_kmer_{start_pos}", description="")
 randomized_kmers.append(kmer_record)

 # Estimate sampling effort
 # A rough estimation: the probability of a specific k-mer being sampled is (kmer_length / total_sequence_length)
 # The probability of NOT sampling a specific k-mer in N trials is (1 - (kmer_length / total_sequence_length))^N
 # To cover 100%, we want this to be close to 0.
 # This is a simplified estimation and doesn't account for overlapping k-mers.
 coverage_estimate = (sample_size * kmer_length) / total_sequence_length * 100

 return randomized_kmers, coverage_estimate

def write_fasta(sequences, output_path):
 """Writes a list of subsequences to a file in FASTA format."""
    with open(output_path, "w") as output_handle:
 SeqIO.write(sequences, output_handle, "fasta")

def write_fastq(sequences, output_path, phred_score=30):
 """Writes a list of subsequences to a file in FASTQ format."""
    with open(output_path, "w") as output_handle:
 for record in sequences:
            qualities = [phred_score] * len(record.seq)
            record.letter_annotations["phred_quality"] = qualities
 SeqIO.write(record, output_handle, "fastq")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Randomly sample k-mers from a FASTA file.")
    parser.add_argument("input_fasta", help="Input FASTA file (can be gzipped)")
    parser.add_argument("output_prefix", help="Output filename prefix")
    parser.add_argument("--sample_size", type=int, required=True, help="Number of random subsequences to sample")
    parser.add_argument("--kmer_length", type=int, required=True, help="Length of subsequences to sample")
    parser.add_argument("--output_format", default="fasta", choices=["fasta", "fastq"], help="Output file format (fasta or fastq, default fasta)")
    parser.add_argument("--phred_score", type=int, default=30, choices=range(10, 41), metavar="[10-40]", help="PHRED score for FASTQ (optional, default 30). Must be between 10 and 40.")

 args = parser.parse_args()

    try:
 # Read the input FASTA file
 print(f"Reading FASTA file: {args.input_fasta}")
 sequences = read_fasta(args.input_fasta)
 print(f"Read {len(sequences)} sequences.")

 # Randomize sequences (sample k-mers)
 print(f"Sampling {args.sample_size} k-mers of length {args.kmer_length}...")
        randomized_kmers, coverage_estimate = randomize_sequences(
 sequences, args.sample_size, args.kmer_length
 )
 print(f"Estimated sampling effort for 100% coverage: {coverage_estimate:.2f}%")
 print(f"Generated {len(randomized_kmers)} random k-mers.")

 # Write the output file
 output_filename = f"{args.output_prefix}.{args.output_format}"
        if args.output_format == "fasta":
 print(f"Writing FASTA output to: {output_filename}")
 write_fasta(randomized_kmers, output_filename)
 elif args.output_format == "fastq":
 print(f"Writing FASTQ output to: {output_filename} with PHRED score {args.phred_score}")
 write_fastq(randomized_kmers, output_filename, args.phred_score)

 print("Script finished successfully.")

    except Exception as e:
 print(f"An error occurred: {e}")

Functions
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
