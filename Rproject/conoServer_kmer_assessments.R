#Smart obj: Per gene_superfamilies, separate and estimates sequence regions complexity based on unsupervised kmer method: k-base wide sliding window. k-mer analysis underpins essential bioinformatics strategies for genome assembly, motif discovery, contamination detection, taxonomic assignment, and comparative genomics, providing a scalable, alignment-free means to dissect and understand sequence datasets in diverse biological contexts.

#' @ Any genomic sequence can be decomposed into a number of consecutive k-mers, and this number will depend on both the length of the sequence (L) and k-mer length (k).
#' @ We estimate the minimum length of k by computing the expected number of occurrences of a given k-mer in a monoploid genome of length G. This is estimated as G/4k, meaning we need to select k to be at least log4(G). For the 3 Gbp human genome (3E9), this would suggest a minimum length of 16; practically, however, the human genome has many 16 bp repeats so it is still advantageous to select even longer k-mers 
#' @ Sequence Composition and Structure Interpretation: Statistics based on k-mers help interpret structural and functional aspects of biological sequences, providing insights into genomic architecture, biomarker identification, and genetic variation detection. They support visualization and analysis of genome composition and facilitate large-scale pattern matching tasks

round(log(3E9, base = 4) )

# for every kmer, evals (Entropy-based.R):
# calculate_entropy
# calculate_information_content
# calculate_custom_complexity
# Conotoxin_protein_fam_diversity.R
# kmer spectra

rm(list = ls())

if(!is.null(dev.list())) dev.off()

library(tidyverse)

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

outdir <- "~/Documents/GitHub/ConotoxinBenchmark/INPUTS/"

file_out <- file.path(outdir, "curated_nuc_conoServerDB.rds")

data <- read_rds(file_out)

data %>%
  mutate(width = nchar(sequence)) %>%
  summarise(width = sum(width))

# 0) choosing kmer according to sequences sizes, and split data in kmers

k <- round(log(591506, base = 4) )

sekmer <- function(read, k = 10) {
  
  kmers <- c()
  
  for (i in 1:(nchar(read) - k + 1)) {
    kmers <- c(kmers, substr(read, i, i + k - 1))
  }
  
  return(kmers)
}

read <- "CAGTCGATT"

sekmer(read, k = 4)

# Here the split by superfamily must be applied

# data %>% count(genesuperfamily)


# for i in 

seqs <- data %>%
  # filter(genesuperfamily == "B1 superfamily") %>% 
  pull(sequence, name = entry_id) # %>% head()

# str(seqs)

kmers <- sapply(seqs, sekmer)

# kmer_counts <- table(unlist(kmers))

# tail(sort(kmer_counts))

# Plot k-mer spectrum
# plot(density(kmer_counts), main = "K-mer spectrum", xlab = "Occurrence of k-mer", ylab = "Frequency of k-mers")

# SPLIT BY SF,  AND COUNT KMER SPECTRUM, richness and global entropy

# head(data.frame(kmer_counts))

# Shannon diversity (entropy)

estimates <- function(kmers) {
  
  kmer_counts <- table(unlist(kmers))
  
  kmer_freqs <- kmer_counts / sum(kmer_counts)
  
  kmer_entropy <- -sum(kmer_freqs * log2(kmer_freqs + .Machine$double.eps))
  
  # kmer_entropy <- -sum(kmer_freqs * log(kmer_freqs))
  
  # Richness (number of unique k-mers)
  kmer_richness <- length(kmer_counts)
  
  # print(paste("Shannon entropy:", entropy))
  # print(paste("Richness:", richness)) 
  
  data.frame(kmer_entropy, kmer_richness)
  
}

# Global (per sf)
estimates(kmers)

# Individual, (per sequences)

OUT <- do.call(rbind, lapply(kmers, estimates)) 

OUT <- OUT %>% mutate(entry_id = rownames(.)) %>% as_tibble()


OUT %>% 
  left_join(distinct(data, sequence, entry_id, genesuperfamily)) %>%
  drop_na(genesuperfamily) %>%
  mutate(
    information_content = sapply(sequence, function(x) sapply(x, calculate_information_content)),
    seq_entropy = sapply(sequence, function(x) sapply(x, calculate_entropy)),
    seq_complexity = sapply(sequence, function(x) sapply(x, calculate_custom_complexity))) %>%
  pivot_longer(cols = c("kmer_entropy", "information_content", "seq_entropy", "seq_complexity")) %>%
  ggplot(aes(y = genesuperfamily, x = value)) +
  facet_grid(~ name, scales = "free_x") +
  geom_boxplot()
  

# 1) estimates number of shared and diverse kmer
#' @ we find that k-mer-based measures of genetic diversity scale consistently with pairwise nucleotide diversity. Some statistics for informing number of shared vs differing k-mers  include cumulative relative entropy (Sims et al. 2009), relative sequence divergence (Sims et al. 2009), average number of common features (Zhang et al. 2017), Shannon’s diversity index (Zhang et al. 2017), and chi-square tests (Bai et al. 2017). https://academic.oup.com/mbe/article/42/3/msaf047/8052716



# Apply the entropy function to every k-mer in each list element

entropy_results <- lapply(kmers, function(kmers) sapply(kmers, calculate_entropy))

info_content_value <- lapply(set_kmer, function(kmers) sapply(kmers, calculate_information_content))

complexity_value <- lapply(set_kmer, function(kmers) sapply(kmers, calculate_custom_complexity))


entropy_results
info_content_value
complexity_value

table(unlist(entropy_results))

LONGERDF <- data.frame(entry_id = rep(names(kmers),
  sapply(kmers, length)),
  kmer = unlist(kmers), row.names = NULL) %>% as_tibble()

LONGERDF <- LONGERDF %>% 
  # count(kmer, sort = T) %>%
  mutate(
    entropy = sapply(kmer, function(x) sapply(x, calculate_entropy)),
    info = sapply(kmer, function(x) sapply(x, calculate_information_content)),
    complexity = sapply(kmer, function(x) sapply(x, calculate_custom_complexity))) %>%
  left_join(data)


# LONGERDF %>%
#   select_if(is.double) %>%
#   mutate_all(~replace(., is.na(.), 0)) %>%
#   cor(method = "spearman") -> M
# # 
# # 
# testRes <- corrplot::cor.mtest(M, conf.level = 0.95)
# # 
# corrplot::corrplot(M, p.mat = testRes$p ,method = "color", type="upper", order = "hclust", insig = "label_sig")
#   

LONGERDF %>%
  ggplot(aes(entropy, group = entry_id)) +
  facet_grid(~ genesuperfamily) +
  geom_density()
  # ggplot2::stat_ecdf()

LONGERDF %>%
  pivot_longer(cols = c("entropy", "info", "complexity")) %>%
  ggplot(aes(value)) +
  facet_wrap(name ~., scales = "free") +
  geom_histogram()


LONGERDF %>% left_join(data) %>%
  group_by(genesuperfamily) %>%
  summarise(entropy = mean(entropy), n = n())

# 2) Estimate genomic complexity using espectra (https://github.com/TGAC/KAT)
#' @ K-mer spectra reveal information not only about the data quality (level of errors, sequencing biases, completeness of sequencing coverage and potential contamination) but also of genomic complexity (size, karyotype, levels of heterozygosity and repeat content; Simpson, 2014). 

# 3) Vizualize as MDS, or t-SNE (https://pmc.ncbi.nlm.nih.gov/articles/PMC12047656/)

