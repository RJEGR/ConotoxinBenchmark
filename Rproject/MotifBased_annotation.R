
#Smart obj: separate and annotate sequence regions based on query similarity to curated reference or consensus sequence patterns.

# Load necessary libraries
library(DECIPHER)
library(universalmotif)
library(Biostrings)

# Step 1: Read query and reference sequences
query <- readAAStringSet("query_sequences.fasta")
ref_sig <- readAAStringSet("signal_reference.fasta")
ref_pro <- readAAStringSet("propeptide_reference.fasta")
ref_mat <- readAAStringSet("mature_reference.fasta")

# Step 2: Create motifs from reference alignments (using universalmotif)
# https://bioconductor.org/packages/release/bioc/html/universalmotif.html
sig_motif <- create_motif(ref_sig)
pro_motif <- create_motif(ref_pro)
mat_motif <- create_motif(ref_mat)

# Step 3: Scan query for motifs
sig_matches <- scan_sequences(sig_motif, query, threshold = 0.8)
pro_matches <- scan_sequences(pro_motif, query, threshold = 0.8)
mat_matches <- scan_sequences(mat_motif, query, threshold = 0.8)

# Step 4: Compile annotated results
annotations <- data.frame(query = names(query),
  signal = sig_matches,
  propeptide = pro_matches,
  mature = mat_matches)

## Annotation

library(universalmotif)
library(Biostrings)
# If references are PWMs, read as motifs; if FASTA, derive motif models first
ref_signal <- create_motif(readAAStringSet("signal_reference.fasta"))
ref_pro <- create_motif(readAAStringSet("propeptide_reference.fasta"))
ref_mat <- create_motif(readAAStringSet("mature_reference.fasta"))
query <- readAAStringSet("query.fasta")

# Scan each query for each motif and annotate
sig_hits <- scan_sequences(ref_signal, query, threshold=0.8) # threshold: similarity
pro_hits <- scan_sequences(ref_pro, query, threshold=0.8)
mat_hits <- scan_sequences(ref_mat, query, threshold=0.8)

# Combine and format results
df <- data.frame(
  Sequence = names(query),
  SignalStart = sapply(sig_hits, function(x) if(length(x)) x[[1]]$start else NA),
  ProStart = sapply(pro_hits, function(x) if(length(x)) x[[1]]$start else NA),
  MatureStart = sapply(mat_hits, function(x) if(length(x)) x[[1]]$start else NA)
)
