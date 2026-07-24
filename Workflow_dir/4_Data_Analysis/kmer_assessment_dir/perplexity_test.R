

calculate_perplexity <- function(sequence) {
  
  calculate_perplexity_ <- function(s, pspm) {
    # s: character vector of observed nucleotides in the k-mer (length k)
    # pspm: matrix k x 4, columns: A, C, G, T with predicted probabilities at each position
    
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
  
  s <- strsplit(sequence, "")[[1]]
  
  k <- length(s)
  
  # Uniform PSPM example
  pspm_uniform <- matrix(0.25, nrow = k, ncol = 4, dimnames = list(NULL, c("A","C","G","T")))
  
  calculate_perplexity_(s, pspm_uniform)
  
}

# Example:

sequence <- "ACTGACTGACTA"
kmers <- c("ACTGACTGACTA", "ACTAACTTACTA", "ACTGACTAACTA", "ACTGACTGACTG")

calculate_perplexity(kmers)

sequences_list <-splitted_seqs

perplexities <- lapply(sequences_list, calculate_perplexity)

# Optionally convert the list to numeric vector
perplexities_vec <- unlist(perplexities)

# View the perplexities per sequence
print(perplexities_vec)
