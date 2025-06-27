import random
import gzip
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def read_fasta(file_path):
    try:
        if file_path.endswith(".gz"):
            with gzip.open(file_path, "rt") as handle:
                return list(SeqIO.parse(handle, "fasta"))
        else:
            return list(SeqIO.parse(file_path, "fasta"))
    except Exception as e:
        print(f"Error reading FASTA file: {e}")
        return []
    
def estimate_full_coverage(sequences, kmer_length):
    """
    Estimates the minimal number of k-mers needed to cover 100% of the sequences.
    """
    
    total_length = sum(len(seq.seq) for seq in sequences)
    
    if total_length == 0:
        print("No sequences found in the input file. Please check the file.")
        return 0
    
    print(f"Total length of sequences: {total_length}")
    
    covered_positions = set()
    iterations = 0


    while len(covered_positions) < total_length:
        iterations += 1
        seq_record = random.choice(sequences)
        if len(seq_record.seq) >= kmer_length:
            start_pos = random.randint(0, len(seq_record.seq) - kmer_length)
            kmer_positions = range(start_pos, start_pos + kmer_length)
            covered_positions.update(kmer_positions)

    
    #print(f"Estimated iterations needed for full coverage: {iterations}")

    return iterations

def main(input_fasta, kmer_length):
    sequences = read_fasta(input_fasta)
    iterations_needed = estimate_full_coverage(sequences, kmer_length)
    print(f"Estimated iterations needed for full coverage: {iterations_needed}")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Randomize sequences in a FASTA file by k-mer length.")
    parser.add_argument("input_fasta", help="Input FASTA (or gzipped FASTA) file")
    parser.add_argument("kmer_length", type=int, help="kmer length to randomize")
    args = parser.parse_args()
    
    main(args.input_fasta, args.kmer_length)