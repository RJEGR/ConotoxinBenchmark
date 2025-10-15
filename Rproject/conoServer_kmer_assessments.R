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

data %>%
  mutate(width = nchar(sequence)) %>%
  summarise(width = sum(width))

# 0) choosing kmer according to sequences sizes, and split data in kmers

k <- round(log(591506, base = 4) )

kmer_assessment <- function(sequences, window_vec = 10) {
  
  
  n_seqs <- length(sequences)
  
  kwindows <- function(read, k = window) {
    
    kmers <- c()
    
    for (i in 1:(nchar(read) - k + 1)) {
      kmers <- c(kmers, substr(read, i, i + k - 1))
    }
    
    return(kmers)
  }
  
  estimates <- function(kmers) {
    
    kmer_counts <- table(unlist(kmers))
    
    kmer_freqs <- kmer_counts / sum(kmer_counts)
    
    # kmer_entropy <- -sum(kmer_freqs * log2(kmer_freqs + .Machine$double.eps))
    
    kmer_entropy <- -sum(kmer_freqs * log(kmer_freqs))
    
    kmer_richness <- length(kmer_counts)
    
    # Calculate evenness: H / log(k)
    # evenness will be between 0 and 1, enabling comparability across samples and window sizes
    evenness <- if (kmer_richness > 1) kmer_entropy / log(kmer_richness) else 0
    
    data.frame(kmer_entropy, kmer_richness, evenness)
    
  }
  
  results <- lapply(window_vec, function(win) {
    kmers <- lapply(sequences, kwindows, k = win)
    data.frame(estimates(kmers), window = win, n_seqs = length(sequences))
  })
  
  # Combine results for all window sizes
  do.call(rbind, results)
  
  # kmers <- kwindows(sequences)
  
  # kmers <- sapply(sequences, kwindows)
  
  # data.frame(estimates(kmers), n_seqs)
  
}


# ==========
# Evals (first) using subsampling, from 0.1 to 1, by = 0.1 ======
# ===========

library(rsample)

# rsample::bootstraps(data, times = 10)

# rsample::initial_split(data, prop = 1)

# Define fractions from 0.1 to 1.0 (step 0.1)
fractions <- seq(0.1, 0.9, by = 0.1)

# Function to create subsample at given fraction
subsample_at_fraction <- function(data, frac) {
  
 
  initial_split(data, prop = frac) %>% analysis() %>% mutate(prop = frac) 
    
}


# Generate list of subsamples
subsamples <- lapply(fractions, function(f) subsample_at_fraction(data, f))

# Optional: name each subsample by fraction
names(subsamples) <- paste0("sample_", fractions)

subsamples <- dplyr::bind_rows(subsamples)

splitted_seqs <- split(strsplit(subsamples$sequence, ";") , subsamples$prop)

splitted_seqs <- lapply(splitted_seqs, unlist)

seq_results <- lapply(splitted_seqs, kmer_assessment, window_vec = k)

# Combine results to a single data.frame if desired:
seq_results <- do.call(rbind, lapply(names(seq_results), function(nm) {
  df <- seq_results[[nm]]
  df$Prop <- nm
  df
}))

head(seq_results)

seq_results <- data %>%
  pull(sequence, name = entry_id) %>%
  kmer_assessment(window_vec = k) %>%
  mutate(Prop = 1.0) %>%
  rbind(seq_results)

seq_results %>% 
  as_tibble() %>%
  ggplot(aes(x = as.factor(Prop), y = kmer_entropy)) +
  # geom_jitter(position = position_jitter(0.1), shape = 1) +
  stat_summary(fun = "mean", geom = "line", aes(group = 1), color="blue") +
  # stat_summary(fun = "mean", geom = "point", aes(group = 1), color="blue") + # , color="blue"
  stat_summary(fun.data=mean_se, geom="pointrange", color="blue") +
  # labs(x = "k-mer size", y = "", caption = "kmer.R (Assembly with spades)", tag = subti) +
  my_custom_theme(plot.tag = element_text(size = 7, hjust = 0))

stratify <- function(data, k = 10, fractions = seq(0.1, 0.9, by = 0.1)) {

  fractions <- seq(0.1, 0.9, by = 0.1)
  
  # Function to create subsample at given fraction
  subsample_at_fraction <- function(data, frac) {
    
    cat("\n Proportion: ", frac, "\n")
    
    initial_split(data, prop = frac) %>% analysis() %>% mutate(prop = frac) 
    
  }
  
  
  # Generate list of subsamples
  subsamples <- lapply(fractions, function(f) subsample_at_fraction(data, f))
  
  # Optional: name each subsample by fraction
  names(subsamples) <- paste0("sample_", fractions)
  
  subsamples <- dplyr::bind_rows(subsamples)
  
  splitted_seqs <- split(strsplit(subsamples$sequence, ";") , subsamples$prop)
  
  splitted_seqs <- lapply(splitted_seqs, unlist)
  
  seq_results <- lapply(splitted_seqs, kmer_assessment, window_vec = k)
  
  # Combine results to a single data.frame if desired:
  seq_results <- do.call(rbind, lapply(names(seq_results), function(nm) {
    df <- seq_results[[nm]]
    df$Prop <- nm
    df
  }))
  
  return(seq_results)
}

stratify(data)

stratify_ <- function(s) {
  
  cat("\n Label ", labels(s)$id, "\n")
  
  analysis(s) %>%
    stratify() %>%
    mutate(strata = labels(s)$id)
}

folds <- rsample::vfold_cv(data, v = 12) # strata = stratified_sampling

# stratify_(folds$splits[[2]])

seq_results <- dplyr::bind_rows(lapply(folds$splits, stratify_)) 

seq_results %>% dplyr::count(Prop)

p <- seq_results %>% 
  as_tibble() %>%
  ggplot(aes(x = as.factor(Prop), y = kmer_entropy)) +
  geom_jitter(position = position_jitter(0.1), shape = 1) +
  stat_summary(fun = "mean", geom = "line", aes(group = 1), color="#56B4EA") +
  # stat_summary(fun = "mean", geom = "point", aes(group = 1), color="blue") + # , color="blue"
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="#56B4EA") +
  # labs(x = "k-mer size", y = "", caption = "kmer.R (Assembly with spades)", tag = subti) +
  my_custom_theme() +
  labs(y = "Shannon Entropy (k-mer diversity)", x = "Sample size")

p

ggsave(p,
    filename = 'conoServer_kmer_assessment_1.png',
    path = outdir, width = 5.5, height = 5, dpi = 1000, device = png)


# ======
# Here the split by superfamily must be applied =====
# ======
# For any (or every) fold-sample, run a screening of kmer assessments?
# Or, for any superfamily, run a screening of kmer assessment?
# data %>% count(genesuperfamily)


# Write a code to split sequences by gene superfamily, then evals kmer assessment

#' @ Note: as some sequence are not well distributed, may be need to subsampling data using rsampling

splitted_seqs <- split(strsplit(data$sequence, ";") , data$genesuperfamily)

splitted_seqs <- lapply(splitted_seqs, unlist)


# kmer_assessment(splitted_seqs$`D superfamily`, window_vec = c(10, 12))

results <- lapply(splitted_seqs, kmer_assessment, window_vec = 10)

# Combine results to a single data.frame if desired:
full_results <- do.call(rbind, lapply(names(results), function(nm) {
  df <- results[[nm]]
  df$genesuperfamily <- nm
  df
}))

head(full_results)


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

kmers_combined <- kmer_merge(splitted_seqs, window_vec = 10)

m <- kmers_combined$k10

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


dist_mat <- dist(t(m))

mds_res <- cmdscale(dist_mat, k = 2)


PCAdf <- PCA(t(m))

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




# -=======
# Evals according to the screening of kmers =====
# =========
window_vec <- c(k, 19, 21, 25, 33, 49, 55, 73)

kmer_assessment(seqs[1:3], window_vec)

# kmer_assessment must accept a sequence set + window_vec
results <- lapply(splitted_seqs, kmer_assessment, window_vec = window_vec)

# Combine results to a single data.frame if desired:
full_results <- do.call(rbind, lapply(names(results), function(nm) {
  df <- results[[nm]]
  df$genesuperfamily <- nm
  df
}))

head(full_results)

# Normalize by the Global 
# kmer_entropy / log(kmer_richness)

seqs <- data %>%
  # filter(genesuperfamily == "B1 superfamily") %>% 
  pull(sequence, name = entry_id) # %>% head()

global_results <- kmer_assessment(seqs, window_vec = window_vec)

ref_vec <- structure(log(global_results$kmer_richness), names = global_results$window)

# Match 'window' values in df to names in ref_vec, getting a vector of divisors
div_vals <- ref_vec[as.character(full_results$window)]

full_results$evenness2 <- full_results$kmer_entropy / div_vals

# 
# Low Entropy: Sequence dataset is dominated by a few k-mers (low diversity, strong bias, or overrepresentation).

# High Entropy: Many different k-mers occur with similar frequency, indicating high genetic/sequence diversity and lower predictability.


full_results %>% 
  as_tibble() %>%
  arrange(desc(n_seqs)) %>% mutate(genesuperfamily = factor(genesuperfamily, levels = unique(genesuperfamily))) %>%
  ggplot(aes(y = genesuperfamily, x = as.factor(window), fill = evenness2)) +
  geom_tile() +
  my_custom_theme(plot.tag = element_text(size = 7, hjust = 0))


full_results %>% 
  as_tibble() %>%
  ggplot(aes(kmer_richness, kmer_entropy, color = as.factor(window))) + geom_point()

full_results %>%
  as_tibble() %>%
  # filter(window == 10) %>%
  arrange(desc(n_seqs)) %>% mutate(genesuperfamily = factor(genesuperfamily, levels = unique(genesuperfamily))) %>%
  ggplot(aes(y = genesuperfamily, x = evenness2, alpha = kmer_entropy)) +
  geom_col() +
  facet_grid(~ window) + labs(x = "Evenness (kmer_entropy / log(kmer_richness)") +
  my_custom_theme(plot.tag = element_text(size = 7, hjust = 0))


full_results %>%
  as_tibble() %>%
  arrange(desc(n_seqs)) %>% mutate(genesuperfamily = factor(genesuperfamily, levels = unique(genesuperfamily))) %>%
  ggplot(aes(y = genesuperfamily, x = kmer_entropy)) +
  # facet_wrap(~ genesuperfamily) +
  geom_jitter(position = position_jitter(0.1), shape = 1) +
  # stat_summary(fun = "mean", geom = "line", aes(group = 1), color="blue") +
  # stat_summary(fun = "mean", geom = "point", aes(group = 1), color="blue") + # , color="blue"
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="blue") +
  # labs(x = "k-mer size", y = "", caption = "kmer.R (Assembly with spades)", tag = subti) +
  my_custom_theme(plot.tag = element_text(size = 7, hjust = 0))


# By sequence, ======

splitted_seqs <- split(strsplit(data$sequence, ";") , data$entry_id)

splitted_seqs <- lapply(splitted_seqs, unlist)

seq_results <- lapply(splitted_seqs, kmer_assessment, window_vec = window_vec)

# Combine results to a single data.frame if desired:
full_seq_results <- do.call(rbind, lapply(names(seq_results), function(nm) {
  df <- seq_results[[nm]]
  df$entry_id <- nm
  df
}))

head(full_seq_results)

gslevels <- full_results %>%
  as_tibble() %>%
  arrange(desc(n_seqs)) %>% mutate(genesuperfamily = factor(genesuperfamily, levels = unique(genesuperfamily))) %>% distinct(genesuperfamily) %>% pull() %>% levels()

full_seq_results %>% 
  as_tibble() %>%
  filter(window == 10) %>%
  left_join(distinct(data, sequence, entry_id, genesuperfamily)) %>%
  as_tibble() %>%
  drop_na(genesuperfamily) %>%
  mutate(genesuperfamily = factor(genesuperfamily, levels = gslevels)) %>%
  ggplot(aes(y = genesuperfamily, x = kmer_entropy)) +
  # facet_wrap(~ genesuperfamily) +
  geom_jitter(position = position_jitter(0.1), shape = 1) +
  # stat_summary(fun = "mean", geom = "line", aes(group = 1), color="blue") +
  # stat_summary(fun = "mean", geom = "point", aes(group = 1), color="blue") + # , color="blue"
  stat_summary(fun.data=mean_se, geom="pointrange", color="blue") +
  # labs(x = "k-mer size", y = "", caption = "kmer.R (Assembly with spades)", tag = subti) +
  my_custom_theme(plot.tag = element_text(size = 7, hjust = 0))

# Assessments by kmer =====

sequences <- splitted_seqs$N00001
  
kwindows <- function(read, k = window) {
  
  kmers <- c()
  
  for (i in 1:(nchar(read) - k + 1)) {
    kmers <- c(kmers, substr(read, i, i + k - 1))
  }
  
  return(kmers)
}

results <- lapply(window_vec, function(win) {
  kmers <- lapply(sequences, kwindows, k = win)
  data.frame(estimates(kmers), window = win, n_seqs = length(sequences))
})

full_seq_results %>% 
  as_tibble() %>%
  filter(window == 10) %>%
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

