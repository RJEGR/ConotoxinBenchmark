# This script requires the Biostrings package for handling DNA sequences.
# You can install it with the following command:
# install.packages("Biostrings")

# Load the required library
library(Biostrings)

#' Calculate the Shannon Entropy for a DNA sequence.
#'
#' This function computes the Shannon Entropy, a measure of the randomness or
#' uncertainty of a DNA sequence based on the frequency of each nucleotide (A, C, G, T).
#'
#' @param sequence A character string representing the DNA sequence.
#' @return The Shannon Entropy value.
#'
#' @examples
#' dna_seq1 <- "ATGCATGC" # High diversity, high entropy
#' calculate_entropy(dna_seq1)
#'
#' dna_seq2 <- "AAAAAAAAAA" # Low diversity, zero entropy
#' calculate_entropy(dna_seq2)
#'
calculate_entropy <- function(sequence) {
  # Convert the sequence to an object that Biostrings can work with.
  dna <- DNAString(sequence)
  
  # Count the frequency of each nucleotide (A, C, G, T).
  counts <- alphabetFrequency(dna)
  
  # Normalize counts to get the probabilities (p_i).
  probabilities <- counts[c("A", "C", "G", "T")] / sum(counts)
  
  # Calculate the entropy using the Shannon formula.
  # We use a small epsilon to avoid log(0) for nucleotides with zero frequency.
  entropy <- -sum(probabilities * log2(probabilities + .Machine$double.eps))
  
  return(entropy)
}

#' Calculate the Information Content for a DNA sequence.
#'
#' This function computes the information content, which is a measure of the
#' non-randomness or conservation of a sequence. It assumes an equal background
#' probability for each nucleotide (0.25).
#'
#' @param sequence A character string representing the DNA sequence.
#' @return The Information Content value.
#'
#' @examples
#' dna_seq1 <- "ATGCATGC" # Lower information content (more random)
#' calculate_information_content(dna_seq1)
#'
#' dna_seq2 <- "AAAAAAAAAA" # High information content (highly conserved)
#' calculate_information_content(dna_seq2)
#'
calculate_information_content <- function(sequence) {
  # Convert the sequence to an object that Biostrings can work with.
  dna <- DNAString(sequence)
  
  # Count the frequency of each nucleotide.
  counts <- alphabetFrequency(dna)
  
  # Normalize counts to get the probabilities (f_b).
  observed_probs <- counts[c("A", "C", "G", "T")] / sum(counts)
  
  # Set the expected background probabilities (p_b).
  # We assume an equal probability for a random sequence (0.25 for each).
  expected_probs <- rep(0.25, 4)
  
  # Calculate the information content using the formula.
  information <- sum(observed_probs * log2(observed_probs / expected_probs + .Machine$double.eps))
  
  return(information)
}

#' Calculate a custom complexity value for a DNA sequence.
#'
#' This function calculates a custom complexity value based on the number of
#' changes between adjacent bases, normalized by the sequence length minus one.
#'
#' @param sequence A character string representing the DNA sequence.
#' @return The custom complexity value as a percentage.
#'
#' @examples
#' dna_seq_custom <- "AAAATTTTTTTTTTTTTTTTTTTTTGGGGGGGGGGGGGGGGGGGGGGCCCC"
#' calculate_custom_complexity(dna_seq_custom)
#'
calculate_custom_complexity <- function(sequence) {
  # Check for an empty or single-base sequence.
  if (nchar(sequence) < 2) {
    return(0)
  }
  
  # Split the sequence into individual bases.
  bases <- strsplit(sequence, "")[[1]]
  
  # Count the number of times a base is different from its next base.
  changes <- 0
  for (i in 1:(length(bases) - 1)) {
    if (bases[i] != bases[i+1]) {
      changes <- changes + 1
    }
  }
  
  # Calculate the complexity value based on the formula.
  complexity <- (changes / (length(bases) - 1)) * 100
  
  return(complexity)
}

# Example usage with your provided sequence:
dna_sequence_example <- "AAAATTTTTTTTTTTTTTTTTTTTTGGGGGGGGGGGGGGGGGGGGGGCCCC"
entropy_value <- calculate_entropy(dna_sequence_example)
info_content_value <- calculate_information_content(dna_sequence_example)
complexity_value <- calculate_custom_complexity(dna_sequence_example)

cat("The Shannon Entropy of the sequence is:", entropy_value, "\n")
cat("The Information Content of the sequence is:", info_content_value, "\n")
cat("The custom complexity of the sequence is:", complexity_value, "%\n")


dna_seq1 <- "ATGCATGC" # Lower information content (more random), but High entropy (diversity)
calculate_entropy(dna_seq1)
calculate_information_content(dna_seq1)
calculate_custom_complexity(dna_seq1)
  
dna_seq2 <- "AAAAAAAAAA" # High information content (highly conserved), but Null entropy
calculate_entropy(dna_seq2)
calculate_information_content(dna_seq2)
calculate_custom_complexity(dna_seq2)
