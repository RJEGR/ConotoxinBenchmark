import argparse
import gzip
import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def read_fasta(file_path):
    """
    Reads a (possibly gzipped) FASTA file and returns a list of sequence records.
    """
    if file_path.endswith(".gz"):
        with gzip.open(file_path, "rt") as handle:
            sequences = list(SeqIO.parse(handle, "fasta"))
    else:
        with open(file_path, "r") as handle:
            sequences = list(SeqIO.parse(handle, "fasta"))
    return sequences

def randomize_sequences(sequences, sample_size, kmer_length):
    """
    Randomly selects k-mers from the provided sequences.
    """
    kmer_list = []
    for _ in range(sample_size):
        while True:
            seq_record = random.choice(sequences)
            if len(seq_record.seq) >= kmer_length:
                start_pos = random.randint(0, len(seq_record.seq) - kmer_length)
                kmer = seq_record.seq[start_pos:start_pos + kmer_length]
                kmer_list.append(SeqRecord(kmer, id=f"{seq_record.id}_kmer_{start_pos}", description=""))
                break
    return kmer_list

def write_fasta(sequences, output_path):
    """
    Writes a list of subsequences to a file in FASTA format.
    """
    with open(output_path, "w") as handle:
        SeqIO.write(sequences, handle, "fasta")

def write_fastq(sequences, output_path, phred_score=40):
    """
    Writes a list of subsequences to a file in FASTQ format.
    """
    if not (10 <= phred_score <= 40):
        raise ValueError("PHRED score must be between 10 and 40.")
    for seq_record in sequences:
        seq_record.letter_annotations["phred_quality"] = [phred_score] * len(seq_record.seq)
    with open(output_path, "w") as handle:
        SeqIO.write(sequences, handle, "fastq")

def main():
    parser = argparse.ArgumentParser(description="Process FASTA/FASTQ files and extract random k-mers.")
    parser.add_argument("input_file", help="Input FASTA file (possibly gzipped).")
    parser.add_argument("output_prefix", help="Output filename prefix.")
    parser.add_argument("--sample_size", type=int, required=True, help="Number of random subsequences to sample.")
    parser.add_argument("--kmer_length", type=int, required=True, help="Length of subsequences (k-mers).")
    parser.add_argument("--output_format", choices=["fasta", "fastq"], default="fasta", help="Output format (FASTA or FASTQ).")
    parser.add_argument("--phred_score", type=int, default=40, help="PHRED score for FASTQ output (default: 40).")
    args = parser.parse_args()

    sequences = read_fasta(args.input_file)
    kmer_list = randomize_sequences(sequences, args.sample_size, args.kmer_length)

    output_file = f"{args.output_prefix}.{args.output_format}"
    if args.output_format == "fasta":
        write_fasta(kmer_list, output_file)
    elif args.output_format == "fastq":
        write_fastq(kmer_list, output_file, args.phred_score)

if __name__ == "__main__":
    main()