
# Ricardo Gomez-Reyes

# ==========
# Evals subsampling, from 0.1 to 0.9, by = 0.1 ======
# ===========


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

kmer_assessment <- function(sequences, window_vec = k) {
  
  
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
    
    # Here provide any other sequence based assessment <------
    
    data.frame(kmer_entropy, kmer_richness, evenness)
    
  }
  
  results <- lapply(window_vec, function(win) {
    
    cat("\n Using: ", win, "k\n")
    
    kmers <- lapply(sequences, kwindows, k = win)
    data.frame(estimates(kmers), window = win, n_seqs = length(sequences))
  })
  
  # Combine results for all window sizes
  do.call(rbind, results)

  
}


library(rsample)


# Define fractions from 0.1 to 1.0 (step 0.1)
fractions <- seq(0.1, 0.9, by = 0.1)

# Function to create subsample at given fraction
subsample_at_fraction <- function(data, frac) {
  
  
  initial_split(data, prop = frac) %>% analysis() %>% mutate(prop = frac) 
  
}

# As above test, only represent a single-or-one sampling, lets resampling 12folds and redoing the test 

stratify_data <- function(s, window_vec = c(k, 19, 21, 25, 33, 49, 55, 73)) {
  
  stratify <- function(data, k = window_vec, fractions = seq(0.1, 0.9, by = 0.1)) {
    
    # fractions <- seq(0.1, 0.9, by = 0.1)
    
    # Function to create subsample at given fraction
    subsample_at_fraction <- function(data, frac) {
      
      cat("\n Proportion: ", frac, "...")
      
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
  
  cat("\n Label ", labels(s)$id, "\n")
  
  analysis(s) %>%
    stratify() %>%
    mutate(strata = labels(s)$id)
}


folds <- rsample::vfold_cv(data, v = 12) # strata = stratified_sampling

# take a while ....

seq_results <- dplyr::bind_rows(lapply(folds$splits, stratify_data)) 


seq_results %>% dplyr::count(Prop, window)

window_vec = c(k, 19, 21, 25, 33, 49, 55, 73)

seq_results <- data %>%
  pull(sequence, name = entry_id) %>%
  kmer_assessment(window_vec = window_vec) %>%
  mutate(Prop = 1.0, strata = "data") %>%
  rbind(seq_results)
  


seq_results %>%
  write_tsv(file.path(outdir, "1_entropy_sample_size.tsv"))

discrete_scale <- seq_results %>% distinct(window) %>% pull()

n <- length(discrete_scale)

scale_col <- c(ggsci::pal_startrek()(7), ggsci::pal_cosmic()(n-7))

scale_col <- structure(scale_col, names = sort(discrete_scale))

# scale_fill <- ggsci::pal_uchicago(alpha = 0.5)(n_pallet) 

scale_fill <- c(ggsci::pal_startrek(alpha = 0.5)(7), ggsci::pal_cosmic(alpha = 0.5)(n-7))

scale_fill <- structure(scale_fill, names = sort(discrete_scale))


p <- seq_results %>% 
  as_tibble() %>%
  ggplot(aes(x = as.factor(Prop), y = kmer_entropy, group = as.factor(window), color = as.factor(window))) +
  geom_jitter(position = position_jitter(0.1), shape = 1, size = 0.5, alpha = 0.7) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun.data=mean_sdl, geom="pointrange", size = 0.5, alpha = 0.5) +
  my_custom_theme() +
  scale_color_manual("kmer sizes",values = scale_col) +
  labs(y = "Shannon entropy", x = "Sample size")

# p

ggsave(p,
  filename = 'conoServer_kmer_assessment_1.png',
  path = outdir, width = 4.5, height = 5, dpi = 1000, device = png)


# Provide summary of accuracy per assembly (Subsampling.R), and bind to seq_results, then correlates
# as precision, and entropy are indendent measurment, and unequal sample sizes, compare the relationship must be performed at the aggregate level (mean, etc.)

accuracy_df <- read_tsv(file.path(outdir, "Sumbsampling_accuracy.tsv"))

# accuracy_df <- accuracy_df %>% 
#   mutate(sampling_set = as.character(sampling_set)) %>%
#   select(Precision, sampling_set, Assembly) %>% arrange(sampling_set) %>% 
#   group_by(Assembly, sampling_set) %>% mutate(rowid = row_number()) %>%
#   pivot_wider(names_from = Assembly, values_from = Precision) %>%
#   dplyr::rename("Prop" = sampling_set) %>%
#   ungroup() 

# cordf <- seq_results %>% 
#   as_tibble() %>%
#   select(kmer_entropy, Prop, window) %>% arrange(Prop) %>% 
#   group_by(Prop, window) %>% mutate(rowid = row_number()) %>%
#   pivot_wider(names_from = window, values_from = kmer_entropy) %>%
#   # select(-rowid) %>%
#   ungroup() %>%
#   right_join(accuracy_df) %>%
#   mutate_all(~replace(., is.na(.), 0)) %>%
#   mutate_if(is.character, as.double)


accuracy_df <- accuracy_df %>% 
  mutate(sampling_set = as.character(sampling_set)) %>%
  group_by(Assembly, sampling_set) %>%
  summarise(Precision = mean(Precision)) %>%
  pivot_wider(names_from = Assembly, values_from = Precision) %>%
  dplyr::rename("Prop" = sampling_set) 

cordf <- seq_results %>% 
  as_tibble()  %>%
  group_by(Prop, window) %>%
  summarise(kmer_entropy = mean(kmer_entropy)) %>%
  pivot_wider(names_from = window, values_from = kmer_entropy) %>%
  right_join(accuracy_df) %>% ungroup() 
  # ungroup() %>% mutate_if(is.character, as.double)


cordf %>%
  # select(Prop, `21`, Trinity) %>%
  select_if(is.double) %>%
  cor(method = "pearson") -> M

testRes <- corrplot::cor.mtest(M, conf.level = 0.95)
# 
corrplot::corrplot(M, p.mat = testRes$p ,method = "color", type="upper", order = "hclust", insig = "label_sig")


seq_results %>% 
  as_tibble()  %>%
  group_by(Prop, window) %>%
  summarise(kmer_entropy = mean(kmer_entropy)) %>%
  right_join(accuracy_df)  %>%
  pivot_longer(cols = c("Spades", "Trinity"), names_to = "Assembly", values_to = "Precision") %>%
  ggplot(aes(Precision, kmer_entropy, color = Assembly)) +
  facet_grid(~ Assembly, scales = "free") +
  geom_point(shape = 1) +
  geom_smooth(method = lm, se = FALSE) +
  ggpubr::stat_cor(method = "pearson", 
    cor.coef.name = "R", color = "black", family = "GillSans",
    p.accuracy = 0.01) +
  my_custom_theme() 

