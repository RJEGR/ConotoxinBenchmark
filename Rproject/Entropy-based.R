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

# Perplexity
# Perplexity is defined as the exponentiation of the cross-entropy loss, which is equivalent to 1 over the probability given to the correct nucleotide

# Our probability in DNA is 1/4 (0.25) per nucleotide

# Example: A vector of 12-mers (strings of length 12 over A,C,G,T)
kmers <- c("ACTGACTGACTA", "ACTAACTTACTA", "ACTGACTAACTA", "ACTGACTGACTG")

# Alphabet order: A, C, G, T
alphabet <- c("A", "C", "G", "T")

k <- 12

# Initialize matrix to hold probabilities: rows=positions, cols=nucleotides
pspm <- matrix(0, nrow = k, ncol = 4, dimnames = list(NULL, alphabet))

# Count nucleotides at each position
for (pos in 1:k) {
  nt_counts <- table(factor(substr(kmers, pos, pos), levels = alphabet))
  pspm[pos, ] <- nt_counts / sum(nt_counts)
}

print(pspm)


# Example usage:
# true_probs and predicted_probs must be numeric vectors of same length and sum to 1
true_probs <- c(0.2, 0.3, 0.5)
predicted_probs <- c(0.1, 0.4, 0.5)

perplexity_value <- calculate_perplexity(true_probs, predicted_probs)
print(perplexity_value)


# s: character vector representing the observed k-mer sequence, e.g. c("A", "C", "T", "G", ...)
# pspm: matrix (k x 4) with predicted probabilities for A,C,G,T at each position
# Columns of pspm in order: A, C, G, T



calculate_perplexity <- function(s) {
  # s: character vector of observed nucleotides in the k-mer (length k)
  # pspm: matrix k x 4, columns: A, C, G, T with predicted probabilities at each position
  
  s <- strsplit(s, "")[[1]]
  
  alphabet <- c("A", "C", "G", "T") # unique(s)
  
  # Initialize matrix to hold probabilities: rows=positions, cols=nucleotides
  pspm <- matrix(0, nrow = k, ncol = 4, dimnames = list(NULL, alphabet))
  
  # Count nucleotides at each position
  for (pos in 1:k) {
    nt_counts <- table(factor(substr(kmers, pos, pos), levels = alphabet))
    pspm[pos, ] <- nt_counts / sum(nt_counts)
  }
  
  print(pspm)
  
  nuc_to_col <- c(A=1, C=2, G=3, T=4)
  
  k <- length(s)
  
  # Construct true_probs vector from observed sequence as one-hot encoding
  true_probs <- numeric(k * 4)   # flattened vector for all positions and nucleotides
  predicted_probs <- numeric(k * 4)
  
  for (i in seq_len(k)) {
    for (nt in names(nuc_to_col)) {
      idx <- (i - 1) * 4 + nuc_to_col[nt]  # position in flattened vector
      
      # True probability: 1 if nt matches s[i], else 0
      true_probs[idx] <- ifelse(nt == s[i], 1, 0)
      
      # Predicted probability from PSPM
      predicted_probs[idx] <- pspm[i, nuc_to_col[nt]]
    }
  }
  
  # Add epsilon to avoid log(0)
  epsilon <- .Machine$double.eps
  cross_entropy <- -sum(true_probs * log(predicted_probs + epsilon))
  perplexity <- exp(cross_entropy / k)  # normalize by sequence length
  
  return(perplexity)
}

# Example:

sequences_list <- lapply(splitted_seqs[10:12], kwindows, k = 20)


# Uniform PSPM example
pspm_uniform <- matrix(0.25, nrow = 12, ncol = 4, dimnames = list(NULL, c("A","C","G","T")))

perplexities <- lapply(sequences_list, calculate_perplexity, pspm = pspm_uniform)

# Optionally convert the list to numeric vector
perplexities_vec <- unlist(perplexities)

# View the perplexities per sequence
print(perplexities_vec)

