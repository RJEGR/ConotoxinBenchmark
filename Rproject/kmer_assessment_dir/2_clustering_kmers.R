
# Creates a matrix of kmers per superfamily or kmer group
# Clusterng or dim reduction



rm(list = ls())

if(!is.null(dev.list())) dev.off()

library(tidyverse)

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

my_custom_theme <- function(...) {
  base_size = 14
  theme_bw(base_family = "GillSans", base_size = base_size) +
    theme(legend.position = "top",
      strip.placement = "outside",
      strip.background = element_rect(fill = 'white', color = 'white'),
      strip.text = element_text(angle = 0, size = base_size), 
      axis.text = element_text(size = rel(0.7), color = "black"),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.x = element_blank(),
      ...
    )
}

outdir <- "~/Documents/GitHub/ConotoxinBenchmark/INPUTS/"

file_out <- file.path(outdir, "curated_nuc_conoServerDB.rds")

data <- read_rds(file_out)

kmer_merge <- function(sequences, window_vec) {
  
  # Extract k-mers from a sequence for given k
  kwindows <- function(read, k) {
    kmers <- c()
    for (i in 1:(nchar(read) - k + 1)) {
      kmers <- c(kmers, substr(read, i, i + k - 1))
    }
    return(kmers)
  }
  
  # Count k-mers for a list of sequences
  count_kmers <- function(seq_list, k) {
    kmers <- lapply(seq_list, kwindows, k = k)
    counts <- table(unlist(kmers))
    return(counts)
  }
  
  merged_matrices <- lapply(window_vec, function(win) {
    # For each window length, count k-mers per sample group
    count_list <- lapply(sequences, count_kmers, k = win)
    
    # Get all unique k-mers across all groups
    all_kmers <- unique(unlist(lapply(count_list, names)))
    
    # Build count matrix by populating counts per group, missing are zero
    mat <- sapply(count_list, function(cnts) {
      vals <- cnts[all_kmers]
      vals[is.na(vals)] <- 0
      vals
    })
    
    # Name rows and columns for clarity
    rownames(mat) <- all_kmers
    
    colnames(mat) <- names(sequences)
    
    return(mat)
  })
  
  names(merged_matrices) <- paste0("k", window_vec)
  
  return(merged_matrices)
}

kmers_combined <- kmer_merge(splitted_seqs, window_vec = window_vec)

count_kmers <- function(sequences, window_vec) {
  
  require(Biostrings)
  
  sequences <- DNAStringSet(sequences)
  
  kmer_counts_lists <- lapply(sequences, oligonucleotideFrequency, width = window_vec)
  
  return(kmer_counts_lists)
}

# count_list <- lapply(splitted_seqs[10:12], count_kmers, window_vec = 12)

# Estimate k-mer probabilities by normalizing counts to frequencies:

kmer_freq <- sweep(kmer_counts, 1, rowSums(kmer_counts), "/")


# count_mat <- kmers_combined$k73

# count_mat: k-mer count matrix with k-mers as rows and samples as columns
# threshold: minimum total count across all samples to keep a k-mer

filter_low_abundance <- function(count_mat, threshold = 5) {
  # Compute total counts per k-mer across all samples
  total_counts <- rowSums(count_mat)
  
  # Retain only k-mers with counts >= threshold
  filtered_mat <- count_mat[total_counts >= threshold, , drop = FALSE]
  
  return(filtered_mat)
}

# Example usage with one matrix from kmer_combined list:
dim(filtered_mat <- filter_low_abundance(count_mat, threshold = 5))

# Assuming kmers_combined is a list of count matrices, one per window/k-mer size
filtered_mat <- lapply(kmers_combined, function(x) {
  filter_low_abundance(x, threshold = 5)
})

# Then normalize

# Assuming kmers_combined is a list of count matrices, one per window/k-mer size
kmers_relative_freq <- lapply(filtered_mat, function(filtered_mat) {
  # Divide each column by its sum to get relative frequencies
  sweep(filtered_mat, 2, colSums(filtered_mat), FUN = "/")
})


head(kmers_relative_freq$k10)

full_results <- do.call(rbind, lapply(names(kmers_relative_freq), function(nm) {
  df <- kmers_relative_freq[[nm]]
  df$genesuperfamily <- nm
  df
}))


PCA <- function(m) {
  
  
  PCA <- prcomp(m, center = F, scale. = T)
  
  PCAdf <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2])
  
  percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
  
  # percentVar <- round(PCA$sdev/sum(PCA$sdev)*100,1)
  
  PCAvar <- data.frame(
    Eigenvalues = PCA$sdev,
    percentVar = percentVar,
    Varcum = cumsum(percentVar)
    # Method = basename(f)
  )
  
  return(list(PCAdf, PCAvar))
}


dist_mat <- dist(t(kmers_relative_freq$k10))

mds_res <- cmdscale(dist_mat, k = 2)


data.frame(mds_res) %>%
  mutate(genesuperfamily = rownames(.)) %>% 
  # left_join(Manifest, by = "LIBRARY_ID") %>%
  ggplot(., aes(X1, X2)) +
  geom_point() + 
  ggrepel::geom_text_repel(aes(label = genesuperfamily), max.overlaps = 100) +
  my_custom_theme()

PCAdf <- PCA(t(dist_mat))

percentVar1 <- paste0(PCAdf[[2]]$percentVar[1], "% ,",PCAdf[[1]]$percentVar[1])
percentVar1 <- paste0("PC1, VarExp: ", percentVar1, "%")


percentVar2 <- paste0(PCAdf[[2]]$percentVar[2], "% ,",PCAdf[[1]]$percentVar[2])
percentVar2 <- paste0("PC2, VarExp: ", percentVar2, "%")


PCAdf <- PCAdf[[1]]


PCAdf %>%
  mutate(genesuperfamily = rownames(.)) %>% 
  # left_join(Manifest, by = "LIBRARY_ID") %>%
  ggplot(., aes(PC1, PC2)) +
  theme_bw(base_size = 14, base_family = "GillSans") +
  xlab(percentVar1) + ylab(percentVar2) +
  geom_abline(slope = 0, intercept = 0, linetype="dashed", alpha=0.5) +
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.5) +
  # ggrepel::geom_text_repel( family = "GillSans", mapping = aes(label = HPF), size = 3.5, color = "black") +
  # geom_text( family = "GillSans", mapping = aes(label = LIBRARY_ID), size = 3.5, color = "black") +
  theme( legend.position = 'top',
    legend.key.width = unit(0.2, "mm"),
    legend.key.height = unit(0.2, "mm")) -> p

p + geom_point() + 
  ggrepel::geom_text_repel(aes(label = genesuperfamily), max.overlaps = 100) +
  my_custom_theme()


