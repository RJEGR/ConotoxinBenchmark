#!/bin/bash

# Quick sequence type checker
# Usage: ./check_sequence_type.sh input.fasta

INPUT_FASTA="$1"

if [ ! -f "$INPUT_FASTA" ]; then
    echo "Usage: $0 <input.fasta>"
    exit 1
fi

echo "Analyzing sequence composition in: $INPUT_FASTA"

# Get first few sequences for analysis
sequences=$(grep -A1 "^>" "$INPUT_FASTA" | grep -v "^>" | grep -v "^--" | head -10)

# Count nucleotides and amino acids
total_chars=0
nucleotides=0
amino_acids=0

while read -r line; do
    if [ ! -z "$line" ]; then
        # Convert to uppercase
        line=$(echo "$line" | tr '[:lower:]' '[:upper:]')
        
        # Count total characters (excluding spaces/newlines)
        line_length=${#line}
        total_chars=$((total_chars + line_length))
        
        # Count nucleotides (A, T, G, C, U, N)
        nt_count=$(echo "$line" | grep -o '[ATGCUN]' | wc -l)
        nucleotides=$((nucleotides + nt_count))
        
        # Count amino acids (standard 20 + some ambiguous)
        aa_count=$(echo "$line" | grep -o '[ARNDCEQGHILKMFPSTWYV]' | wc -l)
        amino_acids=$((amino_acids + aa_count))
    fi
done <<< "$sequences"

echo "Analysis results:"
echo "  Total characters analyzed: $total_chars"
echo "  Nucleotide characters (ATGCUN): $nucleotides"
echo "  Amino acid characters (20 standard): $amino_acids"

if [ $total_chars -gt 0 ]; then
    nt_percent=$((nucleotides * 100 / total_chars))
    aa_percent=$((amino_acids * 100 / total_chars))
    
    echo "  Nucleotide percentage: ${nt_percent}%"
    echo "  Amino acid percentage: ${aa_percent}%"
    
    if [ $nt_percent -gt 80 ]; then
        echo "  CONCLUSION: Likely NUCLEOTIDE sequences"
    elif [ $aa_percent -gt 80 ]; then
        echo "  CONCLUSION: Likely AMINO ACID sequences"
        echo "  WARNING: For amino acid sequences, remove --dbtype 1 from mmseqs createdb"
    else
        echo "  CONCLUSION: Mixed or ambiguous sequence type"
    fi
fi

echo ""
echo "First few sequence lines:"
echo "$sequences" | head -5