#!/usr/bin/env Rscript

################################################################################
# LOSO_SHANNON_ENTROPY_WORKFLOW.R
# 
# Optimized workflow for estimating Shannon entropy-based diversity metrics
# in Leave-One-Superfamily-Out cross-validation, with Grimm Type-2 
# circularity control via Jensen-Shannon divergence.
#
# INTEGRATION POINT: Run this AFTER LOSO_conoServerDB.R
# INPUT: train_*.fasta, test_*.fasta, LOSO_manifest.tsv
# OUTPUT: LOSO_entropy_metrics.csv, LOSO_superfamily_summary.csv, plots
#
# AUTHOR: Based on LLM Council Analysis with 5 independent perspectives
# DATE: May 2026
################################################################################

library(tidyverse)
library(Biostrings)
library(entropy)
library(DECIPHER)
library(patchwork)

# ============================================================================
# CONFIGURATION
# ============================================================================

LOSO_dir <- Sys.getenv('LOSO_DIR', 'INPUTS/vfolds_loso_resampling_dir')
output_dir <- Sys.getenv('OUTPUT_DIR', LOSO_dir)
K_VALUES <- c(3, 5, 7, 9)
MIN_KMER_COUNT <- 10  # Skip k-mers that appear < 10 times (noise filtering)
SEED <- 20260511

set.seed(SEED)

cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("LOSO SHANNON ENTROPY WORKFLOW - Complexity & Diversity Estimation\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("LOSO_dir:", LOSO_dir, "\n")
cat("k-values:", paste(K_VALUES, collapse=", "), "\n")
cat("Seed:", SEED, "\n\n")

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

# Extract k-mers from DNA sequences (vectorized)
extract_kmers <- function(fasta_file, k) {
  tryCatch({
    dna <- readDNAStringSet(fasta_file)
    if (length(dna) == 0) stop("Empty FASTA")
    
    # Use Biostrings oligonucleotideFrequency for speed
    all_kmers <- c()
    for (seq in dna) {
      kmers <- oligonucleotideFrequency(seq, width = k)
      all_kmers <- c(all_kmers, colnames(kmers)[kmers > 0])
    }
    
    # Count frequency
    table(all_kmers)
  }, error = function(e) {
    warning(paste("Error extracting k-mers from", fasta_file, ":", e$message))
    table()
  })
}

# Shannon entropy: H = -Σ(p_i * log₂(p_i))
shannon_entropy <- function(kmer_freq_table) {
  if (length(kmer_freq_table) == 0) return(NA)
  
  p <- as.numeric(kmer_freq_table) / sum(kmer_freq_table)
  p <- p[p > 0]  # Remove zeros to avoid log(0)
  
  -sum(p * log2(p))
}

# Normalize entropy: H_norm = H / H_max where H_max = log2(4^k)
normalize_entropy <- function(H, k) {
  if (is.na(H)) return(NA)
  H / log2(4^k)
}

# Miller-Madow bias correction for finite samples
miller_madow_correction <- function(kmer_freq_table, k) {
  if (length(kmer_freq_table) == 0) return(NA)
  
  n_distinct <- length(kmer_freq_table)
  n_total <- sum(kmer_freq_table)
  
  (n_distinct - 1) / (2 * log(2) * n_total)
}

# Rényi entropy (generalized entropy): H_α
renyi_entropy <- function(kmer_freq_table, alpha) {
  if (length(kmer_freq_table) == 0) return(NA)
  
  p <- as.numeric(kmer_freq_table) / sum(kmer_freq_table)
  p <- p[p > 0]
  
  if (alpha == 0) {
    return(log2(length(p)))  # Support size
  } else if (abs(alpha - 1) < 1e-10) {
    return(shannon_entropy(kmer_freq_table))  # Shannon (α=1)
  } else if (is.infinite(alpha)) {
    return(-log2(max(p)))  # Min-entropy
  } else {
    return((1 / (1 - alpha)) * log2(sum(p^alpha)))
  }
}

# Jensen-Shannon divergence: JS(P||Q) = √(0.5*KL(P||M) + 0.5*KL(Q||M))
jensen_shannon_divergence <- function(freq1, freq2) {
  if (length(freq1) == 0 || length(freq2) == 0) return(NA)
  
  # Get union of all k-mers
  all_kmers <- union(names(freq1), names(freq2))
  p <- numeric(length(all_kmers))
  q <- numeric(length(all_kmers))
  
  names(p) <- all_kmers
  names(q) <- all_kmers
  p[names(freq1)] <- freq1
  q[names(freq2)] <- freq2
  
  p <- p / sum(p)
  q <- q / sum(q)
  
  m <- 0.5 * (p + q)
  
  # KL divergence (with 0*log(0) = 0 convention)
  kl_pm <- sum(p[p > 0] * log2(p[p > 0] / m[p > 0]))
  kl_qm <- sum(q[q > 0] * log2(q[q > 0] / m[q > 0]))
  
  sqrt(0.5 * kl_pm + 0.5 * kl_qm)
}

# K-mer clustering & diversity index
compute_diversity_index <- function(fasta_file, k = 5) {
  tryCatch({
    dna <- readDNAStringSet(fasta_file)
    if (length(dna) == 0) return(NA)
    
    # Use DECIPHER::Clusterize for sequence-level clustering
    clusters <- DECIPHER::Clusterize(dna, cutoff = 0.5, minCoverage = 0.95)
    
    n_clusters <- length(unique(clusters$cluster))
    n_sequences <- length(dna)
    
    # Diversity = 1 - (homogeneity)
    # High diversity (close to 1): sequences are distinct
    # Low diversity (close to 0): sequences are similar
    diversity <- 1 - (n_clusters / n_sequences)
    
    return(diversity)
  }, error = function(e) {
    warning(paste("Error computing diversity for", fasta_file, ":", e$message))
    return(NA)
  })
}

# ============================================================================
# MAIN ANALYSIS LOOP
# ============================================================================

# Load LOSO manifest
manifest_file <- file.path(LOSO_dir, 'LOSO_manifest.tsv')
if (!file.exists(manifest_file)) {
  stop("LOSO_manifest.tsv not found in", LOSO_dir)
}

loso_manifest <- read_tsv(manifest_file, 
                          col_types = cols(.default = col_character(), 
                                           fold_id = col_integer(),
                                           n_test = col_integer(),
                                           n_train = col_integer()))

cat("Loaded", nrow(loso_manifest), "folds from manifest\n\n")

# Initialize results storage
results_list <- list()

# Iterate through folds
for (fold_idx in seq_len(nrow(loso_manifest))) {
  
  fold_row <- loso_manifest[fold_idx, ]
  fold_id <- fold_row$fold_id
  superfamily <- fold_row$held_out_superfamily
  test_fasta <- fold_row$test_fasta
  train_fasta <- fold_row$train_fasta
  n_test <- fold_row$n_test
  n_train <- fold_row$n_train
  
  cat(sprintf("▶ Fold %2d/%d: %s (n_test=%d, n_train=%d)\n", 
              fold_idx, nrow(loso_manifest), 
              substr(superfamily, 1, 40), n_test, n_train))
  
  # Iterate through k-mer lengths
  for (k in K_VALUES) {
    
    cat(sprintf("  ├─ k=%d: ", k))
    
    # Extract k-mers
    test_kmers <- extract_kmers(test_fasta, k)
    train_kmers <- extract_kmers(train_fasta, k)
    
    if (length(test_kmers) == 0 || length(train_kmers) == 0) {
      cat("SKIP (empty k-mer table)\n")
      next
    }
    
    # Shannon entropy
    H_test_raw <- shannon_entropy(test_kmers)
    H_train_raw <- shannon_entropy(train_kmers)
    
    H_test_norm <- normalize_entropy(H_test_raw, k)
    H_train_norm <- normalize_entropy(H_train_raw, k)
    
    # Bias correction (Miller-Madow)
    correction <- miller_madow_correction(test_kmers, k)
    H_test_adj <- H_test_norm + correction
    
    # Rényi spectrum (robustness check)
    H_0 <- renyi_entropy(test_kmers, 0)
    H_1 <- H_test_raw  # Shannon is α=1
    H_2 <- renyi_entropy(test_kmers, 2)
    H_inf <- renyi_entropy(test_kmers, Inf)
    
    # Jensen-Shannon divergence (distributional shift detection)
    JS_div <- jensen_shannon_divergence(test_kmers, train_kmers)
    OOD_flag <- if_else(JS_div > 0.15, 'HIGH_OOD', 
                        if_else(JS_div > 0.10, 'MOD_OOD', 'OK'))
    
    # K-mer diversity index
    diversity <- compute_diversity_index(test_fasta, k = 5)
    
    # Store results
    result_row <- tibble(
      fold_id = fold_id,
      superfamily = superfamily,
      n_test = n_test,
      n_train = n_train,
      k = k,
      n_distinct_kmers = length(test_kmers),
      H_test_raw = round(H_test_raw, 4),
      H_train_raw = round(H_train_raw, 4),
      H_test_norm = round(H_test_norm, 4),
      H_train_norm = round(H_train_norm, 4),
      H_test_adj = round(H_test_adj, 4),
      H_0_renyi = round(H_0, 4),
      H_1_shannon = round(H_1, 4),
      H_2_renyi = round(H_2, 4),
      H_inf_minentropy = round(H_inf, 4),
      JS_divergence = round(JS_div, 4),
      OOD_flag = OOD_flag,
      diversity_index = round(diversity, 4)
    )
    
    results_list[[paste(fold_id, k)]] <- result_row
    
    cat(sprintf("H=%.3f (test) vs H=%.3f (train), JS=%.4f [%s]\n",
                H_test_norm, H_train_norm, JS_div, OOD_flag))
  }
}

# Combine results into single data frame
entropy_results <- bind_rows(results_list)

cat("\n✓ Analysis complete! Processed", nrow(entropy_results), "fold×k combinations\n\n")

# ============================================================================
# AGGREGATE & SUMMARIZE BY SUPERFAMILY
# ============================================================================

cat("Computing per-superfamily summaries...\n\n")

superfamily_summary <- entropy_results %>%
  group_by(superfamily) %>%
  summarize(
    n_folds = n_distinct(fold_id),
    mean_n_test = round(mean(n_test), 0),
    median_H_test_norm = round(median(H_test_norm, na.rm = TRUE), 4),
    sd_H_test_norm = round(sd(H_test_norm, na.rm = TRUE), 4),
    median_H_adj = round(median(H_test_adj, na.rm = TRUE), 4),
    median_diversity = round(median(diversity_index, na.rm = TRUE), 4),
    median_JS_div = round(median(JS_divergence, na.rm = TRUE), 4),
    n_OOD_flags = sum(OOD_flag != 'OK'),
    pct_OOD = round(100 * n_OOD_flags / n(), 1),
    .groups = 'drop'
  ) %>%
  arrange(desc(median_H_test_norm))

# ============================================================================
# SAVE OUTPUT FILES
# ============================================================================

cat("Writing output files to:", output_dir, "\n")

# 1. Full entropy metrics table
output_file_metrics <- file.path(output_dir, 'LOSO_entropy_metrics.tsv')
write_tsv(entropy_results, output_file_metrics)
cat("✓", output_file_metrics, "\n")

# 2. Per-superfamily summary
output_file_summary <- file.path(output_dir, 'LOSO_entropy_superfamily_summary.tsv')
write_tsv(superfamily_summary, output_file_summary)
cat("✓", output_file_summary, "\n")

# ============================================================================
# STATISTICAL TESTS
# ============================================================================

cat("\n--- STATISTICAL TESTS ---\n\n")

# Test 1: Is H_test < H_train? (Wilcoxon signed-rank test)
test1 <- wilcox.test(entropy_results$H_test_norm, 
                     entropy_results$H_train_norm,
                     paired = FALSE, alternative = 'less')

cat("Test 1: Wilcoxon (H_test < H_train)\n")
cat("  W =", test1$statistic, "\n")
cat("  p-value =", format(test1$p.value, scientific = TRUE, digits = 3), "\n")
cat("  Interpretation: Training entropy is", 
    if_else(test1$p.value < 0.05, "SIGNIFICANTLY", "NOT"), 
    "higher than test entropy\n")

# Test 2: Correlation between n_test and entropy
test2 <- cor.test(entropy_results$n_test, entropy_results$H_test_norm,
                  method = 'spearman', use = 'complete.obs')

cat("\nTest 2: Spearman ρ (n_test vs H_test_norm)\n")
cat("  ρ =", round(test2$estimate, 3), "\n")
cat("  p-value =", format(test2$p.value, scientific = TRUE, digits = 3), "\n")
cat("  Interpretation: ", 
    if_else(abs(test2$estimate) > 0.3, "MODERATE", "WEAK"),
    "correlation between sample size and entropy\n")

# Test 3: Friedman test across superfamilies
test3_data <- entropy_results %>%
  filter(k == 5) %>%  # Use k=5 for simplicity
  select(fold_id, superfamily, H_test_norm) %>%
  pivot_wider(names_from = superfamily, values_from = H_test_norm, values_fill = NA)

if (ncol(test3_data) > 2) {
  test3 <- friedmantst(as.matrix(test3_data[, -1]))
  cat("\nTest 3: Friedman (entropy across superfamilies, k=5)\n")
  cat("  χ² =", round(test3$statistic, 3), "\n")
  cat("  p-value =", format(test3$p.value, scientific = TRUE, digits = 3), "\n")
  cat("  Interpretation: Entropy differs", 
      if_else(test3$p.value < 0.05, "SIGNIFICANTLY", "NOT"),
      "across superfamilies\n")
}

# ============================================================================
# VISUALIZATIONS
# ============================================================================

cat("\n--- GENERATING PLOTS ---\n")

# Create plots
theme_loso <- function() {
  theme_bw(base_family = "sans", base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      strip.background = element_rect(fill = 'gray90', color = 'white'),
      strip.text = element_text(size = 9, face = 'bold'),
      legend.position = 'bottom'
    )
}

# Plot 1: Diversity by superfamily size
p1 <- entropy_results %>%
  filter(k == 5) %>%
  ggplot(aes(x = factor(n_test), y = diversity_index, fill = superfamily)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.3, size = 1) +
  scale_fill_viridis_d(guide = 'none') +
  labs(
    title = 'K-mer Diversity Index by Test Set Size (k=5)',
    x = 'n_test (grouped)',
    y = 'Diversity Index (1 - n_clusters/n_seqs)'
  ) +
  theme_loso()

ggsave(p1, filename = file.path(output_dir, 'LOSO_diversity_by_size.png'),
       width = 10, height = 6, dpi = 150)
cat("✓ LOSO_diversity_by_size.png\n")

# Plot 2: Shannon entropy vs test set size
p2 <- entropy_results %>%
  filter(k %in% c(5, 7)) %>%
  ggplot(aes(x = n_test, y = H_test_norm, color = factor(k), shape = factor(k))) +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_smooth(method = 'loess', se = TRUE, alpha = 0.15) +
  scale_color_brewer(palette = 'Set1', name = 'k-mer length') +
  scale_shape_manual(values = c(16, 17), name = 'k-mer length') +
  facet_wrap(~k, ncol = 2, labeller = labeller(k = c('5'='k=5', '7'='k=7'))) +
  labs(
    title = 'Normalized Shannon Entropy vs Test Set Size',
    x = 'n_test (number of sequences)',
    y = 'H_norm (normalized Shannon entropy, 0-1 scale)'
  ) +
  theme_loso() +
  theme(legend.position = 'none')

ggsave(p2, filename = file.path(output_dir, 'LOSO_H_vs_n_test.png'),
       width = 11, height = 6, dpi = 150)
cat("✓ LOSO_H_vs_n_test.png\n")

# Plot 3: Rényi spectrum (robustness check)
p3 <- entropy_results %>%
  filter(k == 7) %>%
  select(fold_id, superfamily, k, H_0_renyi, H_1_shannon, H_2_renyi, H_inf_minentropy) %>%
  pivot_longer(cols = starts_with('H_'), names_to = 'alpha', values_to = 'entropy') %>%
  ggplot(aes(x = alpha, y = entropy, color = superfamily, group = interaction(fold_id, superfamily))) +
  geom_point(alpha = 0.3, size = 1.5) +
  geom_line(alpha = 0.15, linewidth = 0.3) +
  scale_color_viridis_d(name = 'Superfamily', option = 'turbo') +
  labs(
    title = 'Rényi Entropy Spectrum (k=7): Robustness Check',
    x = 'Entropy Order (α)',
    y = 'Rényi Entropy'
  ) +
  theme_loso() +
  theme(legend.position = 'right', legend.text = element_text(size = 7))

ggsave(p3, filename = file.path(output_dir, 'LOSO_renyi_spectrum.png'),
       width = 12, height = 6, dpi = 150)
cat("✓ LOSO_renyi_spectrum.png\n")

# Plot 4: Jensen-Shannon divergence heatmap
p4 <- entropy_results %>%
  filter(k == 7) %>%
  mutate(fold_label = paste0('F', fold_id)) %>%
  ggplot(aes(x = reorder(superfamily, -JS_divergence),
             y = reorder(fold_label, -JS_divergence),
             fill = JS_divergence)) +
  geom_tile(color = 'white', linewidth = 0.2) +
  scale_fill_gradient(low = '#2ecc71', mid = '#f39c12', high = '#e74c3c',
                      midpoint = 0.10,
                      name = 'JS Divergence\n(bits)') +
  labs(
    title = 'Jensen-Shannon Divergence per Fold × Superfamily (k=7)',
    x = 'Test Superfamily',
    y = 'Fold'
  ) +
  theme_loso() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 7))

ggsave(p4, filename = file.path(output_dir, 'LOSO_jensen_shannon_heatmap.png'),
       width = 14, height = 10, dpi = 150)
cat("✓ LOSO_jensen_shannon_heatmap.png\n")

# Plot 5: OOD flag distribution
p5 <- entropy_results %>%
  filter(k == 7) %>%
  group_by(superfamily, OOD_flag) %>%
  summarize(n = n(), .groups = 'drop') %>%
  ggplot(aes(x = reorder(superfamily, -n), y = n, fill = OOD_flag)) +
  geom_col() +
  scale_fill_manual(values = c('HIGH_OOD' = '#e74c3c', 'MOD_OOD' = '#f39c12', 'OK' = '#2ecc71'),
                    name = 'OOD Flag') +
  labs(
    title = 'Out-of-Distribution (OOD) Flags per Superfamily (k=7)',
    x = 'Superfamily',
    y = 'Number of folds'
  ) +
  theme_loso() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))

ggsave(p5, filename = file.path(output_dir, 'LOSO_OOD_flags.png'),
       width = 12, height = 6, dpi = 150)
cat("✓ LOSO_OOD_flags.png\n")

# ============================================================================
# FINAL SUMMARY REPORT
# ============================================================================

cat("\n" %+% rep("━", 70) %+% "\n")
cat("SUMMARY REPORT\n")
cat(rep("━", 70) %+% "\n\n")

cat("Data Summary:\n")
cat("  • Total folds analyzed:", nrow(loso_manifest), "\n")
cat("  • K-mer lengths tested:", paste(K_VALUES, collapse = ", "), "\n")
cat("  • Total fold×k combinations:", nrow(entropy_results), "\n")
cat("  • Superfamilies:", n_distinct(entropy_results$superfamily), "\n\n")

cat("Entropy Results:\n")
cat("  • Mean H_test_norm:", round(mean(entropy_results$H_test_norm, na.rm=T), 3), "\n")
cat("  • Mean H_train_norm:", round(mean(entropy_results$H_train_norm, na.rm=T), 3), "\n")
cat("  • Mean H_delta (train - test):", 
    round(mean(entropy_results$H_train_norm - entropy_results$H_test_norm, na.rm=T), 3), "\n\n")

cat("Circularity Control (Jensen-Shannon Divergence):\n")
cat("  • Mean JS divergence:", round(mean(entropy_results$JS_divergence, na.rm=T), 4), "bits\n")
cat("  • High OOD flags (JS > 0.15):", sum(entropy_results$OOD_flag == 'HIGH_OOD'), "instances\n")
cat("  • Moderate OOD flags (JS > 0.10):", sum(entropy_results$OOD_flag == 'MOD_OOD'), "instances\n")
cat("  • Interpretation:", 
      if_else(mean(entropy_results$JS_divergence, na.rm=T) < 0.10,
              "Good generalization expected",
              "Possible OOD risk; validate with accuracy metrics"),
      "\n\n")

cat("Diversity Index:\n")
cat("  • Mean diversity (1 - n_clusters/n_seqs):", 
      round(mean(entropy_results$diversity_index, na.rm=T), 3), "\n")
cat("  • Range:", round(min(entropy_results$diversity_index, na.rm=T), 3),
      "-", round(max(entropy_results$diversity_index, na.rm=T), 3), "\n\n")

cat("Top 5 Most Diverse Superfamilies (by median H_test_norm):\n")
print(superfamily_summary %>% head(5) %>% select(superfamily, median_H_test_norm, n_OOD_flags))

cat("\n" %+% rep("━", 70) %+% "\n")
cat("✓ WORKFLOW COMPLETE\n")
cat(rep("━", 70) %+% "\n")

cat("\nNext Steps:\n")
cat("1. Review plots for visual patterns\n")
cat("2. Correlate entropy metrics with assembly accuracy (if available)\n")
cat("3. Use H_test_norm and JS_divergence in downstream ML models\n")
cat("4. Report findings in manuscript\n\n")

cat("Output files created:\n")
cat("  • LOSO_entropy_metrics.tsv\n")
cat("  • LOSO_entropy_superfamily_summary.tsv\n")
cat("  • LOSO_diversity_by_size.png\n")
cat("  • LOSO_H_vs_n_test.png\n")
cat("  • LOSO_renyi_spectrum.png\n")
cat("  • LOSO_jensen_shannon_heatmap.png\n")
cat("  • LOSO_OOD_flags.png\n\n")

cat("Session info:\n")
cat("  • R version:", R.version$version.string, "\n")
cat("  • Biostrings version:", packageVersion("Biostrings"), "\n")
cat("  • DECIPHER version:", packageVersion("DECIPHER"), "\n")
cat("  • Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n")
