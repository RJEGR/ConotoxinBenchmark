# As we found, as kmer size increase, the shannon (complexity) entropy increaseo, and
# # if sample size increase, also the shannon entropy, in the conoServerDB
# # lets calculate the shannon values for those good contigs assembled in the Subsampling step
# 
# /LUSTRE/bioinformatica_data/genomica_funcional/rgomez/fernando_pub/2_subsampling_dir/multiple_analysis_dir/Trinity_FASTA_DIR
# 
# 
# 
rm(list = ls())

if(!is.null(dev.list())) dev.off()

# pak::pkg_install("tidyverse")

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


calculate_metrics <- function(df, reference_coverage_val = 1) {
  
  is_chimeric_value = 0
  
  calculate_false <- function(df) {
    
    # as many assemblers use width 200 to filter contigs, count number of refseq > 200
 
    count_Nsequences <- function() {
      
      # Temporal directory where vfolds_resampling_dir are found
      
      # outdir <- "~/Documents/GitHub/ConotoxinBenchmark/INPUTS/vfolds_resampling_dir/"
      # 
      outdir <- "C://Users//cinai/OneDrive/Documentos/GitHub/ConotoxinBenchmark/INPUTS/vfolds_resampling_dir/"
      
      
      f <- list.files(path = outdir, pattern = ".fasta", full.names = T)
      
      my_func <- function(x) { 
        
        dna <- Biostrings::readDNAStringSet(x)
        
        structure(
          sum(Biostrings::width(dna) >=200), 
          names = gsub(".fasta", "", basename(x)))
        
        
      }
      
      unlist(lapply(f, my_func))
      
    }
    
    
    InputNsequences <- count_Nsequences()
    
    # False Negatives (FN): Transcripts present in the simulated data but not assembled. 
    ## TP - (N reference sequences in InputNsequences) 
    
    
    data.frame(InputNsequences) %>% 
      as_tibble(rownames = "vfold_set") %>% 
      right_join(df) %>%
      mutate(FN = abs(InputNsequences - TP))
    
  }
  
  
  Totaldf <- df %>% 
    # dplyr::count(vfold_set, sampling_set)  %>% 
    dplyr::count() %>%
    dplyr::rename("rawcontigs" = "n")
  
  # True Positives (TP): Transcripts correctly assembled by the assembler. 
  ## TP = reference_cov >= reference_coverage_val ||reference_cov == 1 (= !is.na(hits))

  
  TP <- df %>%
    filter(reference_coverage >= reference_coverage_val) %>%
    # Important: make distinct hits, to count only unique references (Recall), not over-inflate by N contigs!!!
    distinct(hits) %>% 
    dplyr::count() %>%
    dplyr::rename("TP" = "n")
  
  # False Positives (FP): Transcripts incorrectly assembled by the assembler.
  ## FP = N contig_name where reference_cov < 1 BUT have some identity threshold (reference_coverage > 0)
  
  FP <- df %>%
    mutate(reference_coverage = ifelse(is.na(hits) & is.na(reference_coverage), 0, reference_coverage )) %>%
    filter(reference_coverage > is_chimeric_value &  reference_coverage < reference_coverage_val)  %>%
    # dplyr::count(vfold_set, sampling_set) %>% 
    dplyr::count() %>%
    dplyr::rename("FP" = "n")
  
  
  # Dealing with Overestimate contig number
  # No hay manera de usar estos como TN,
  # contigs where reference_coverage == 0, OR
  # Use Rawcontigs number minus TP + FP to count potential True Negative (TN)
  #     mutate(TN = rawcontigs - (TP+FP)) %>%
  
  TN <- df %>%
    mutate(reference_coverage = ifelse(is.na(hits) & is.na(reference_coverage), 0, reference_coverage )) %>%
    # mutate_all(~replace(., is.na(.), 0)) %>%
    filter(reference_coverage <= is_chimeric_value) %>%
    # filter(is.na(hits)) %>%
    # dplyr::count(vfold_set, sampling_set) %>%
    dplyr::count() %>%
    dplyr::rename("TN" = "n")
  
  
  Totaldf %>% 
    left_join(TP) %>% 
    # left_join(TN) %>% 
    left_join(FP) %>% 
    mutate_all(~replace(., is.na(.), 0)) %>%
    left_join(calculate_false(.)) %>%
    select(-InputNsequences) 
  
  
}


read_files <- function(dir, reference_coverage_val = 0.95, Assembly = ...) {
  
  file_list <- list.files(path = dir, pattern = "contigs.csv", recursive = T, full.names = TRUE)
  
  transratedf <- lapply(file_list, read_transrate_scores)
  
  transratedf <- do.call(rbind,transratedf)
  
  transratedf %>% 
    mutate(rel_path = basename(dirname(rel_path))) %>%
    mutate(vfold_set = sapply(strsplit(rel_path, "_"), `[`, 1)) %>%
    mutate(sampling_set = sapply(strsplit(rel_path, "_"), `[`, 5)) %>%
    mutate(sampling_set = as.double(sampling_set)) %>%
    mutate(Assembly = Assembly) %>%
    filter(length >= 200) %>%
    mutate(reference_coverage = ifelse(is.na(hits) & is.na(reference_coverage), 0, reference_coverage )) %>%
    filter(reference_coverage >=  reference_coverage_val)
}


read_transrate_scores <- function(file_list) {
  
  # dir <- dirname(dir)
  rel_path <- str_remove(file_list, paste0("^", dir, "/"))
  
  cat("\nReading\n")
  cat(rel_path)
  cat("\n")
  
  rel_path <- sapply(strsplit(dirname(rel_path), "/"), `[`, 15)
  vfold_set <- sapply(strsplit(rel_path, "_"), `[`, 1)
  Prop <- sapply(strsplit(rel_path, "_"), `[`, 5)
  
  cat("\n")
  cat(vfold_set, Prop)
  cat("\n")
  
  read_csv(file_list) %>%
    # mutate(file_list = file_list)
    mutate(rel_path, vfold_set, Prop)
}


# A) 
# Starting with TRINITY outputs
# dir <- "C://Users//cinai/OneDrive/Documentos/GitHub/ConotoxinBenchmark/INPUTS/2_subsampling_dir/Trinity_dir/transrate_contigs_dir/"
 

dir <- "C://Users//cinai/OneDrive/Documentos/GitHub/ConotoxinBenchmark/INPUTS/2_subsampling_dir/Trinity_dir/transrate_contigs_dir/"

str(file_list <- list.files(path = dir, pattern = "contigs.csv", recursive = T, full.names = TRUE))

transratedf <- lapply(file_list, read_transrate_scores)

transratedf <- do.call(rbind,transratedf)

transratedf %>% distinct(rel_path, vfold_set, Prop)

Benchmark_results <- transratedf %>%
  filter(length>=200) %>%
  group_by(rel_path, vfold_set, Prop) %>%
  calculate_metrics(reference_coverage_val = 0.95) %>% 
  mutate(
    Ratio = TP/FP,
    # Tell us what percentage of positive classes were correctly identified
    Accuracy = TP / (TP + FN + FP),
    Precision = TP /(TP + FP),
    Sensitivity = TP /(TP + FN),
    Fscore = 2 * (TP) / (2 * (TP) + FP + FN), 
  ) 

Benchmark_results %>% distinct(rel_path)

# Step 2: LOAD DATA from sequences 
# 
# 

outdir <- "C://Users//cinai/OneDrive/Documentos/GitHub/ConotoxinBenchmark/INPUTS/Subsampling_dir/"

f <- list.files(file.path(outdir, "Trinity_FASTA_DIR"), full.names = T)


# k <- round(log(sum(Biostrings::width(stringSet)), base = 4) )

# window_vec <- c(k, 19, 21, 25, 33, 49, 55, 73)

shannon_by_contigs <- function(f, window_vec = 33) {
  
  # this function use by default the transratedf DB to track the good assembled contigs 
  
  stringSet <- Biostrings::readAAStringSet(f)
  
  contig_names_l <- sapply(strsplit(names(stringSet), " "), `[`, 1)
  
  rel_path_f <- gsub(".fa", "",basename(f))
  
  cat("\n Reading: ", rel_path_f, " file\n")
  
  keep_stringSet <- transratedf %>%
    filter(length >= 200) %>%
    mutate(reference_coverage = ifelse(is.na(hits) & is.na(reference_coverage), 0, reference_coverage )) %>%
    filter(reference_coverage >=  0.95) %>%
    filter(rel_path == rel_path_f) %>%
    # The follow step is redundant, but is for sanity check
    filter(contig_name %in% contig_names_l) %>%
    distinct(contig_name) %>% pull()
  
  stringSet <- stringSet[contig_names_l %in% keep_stringSet]  
  
  sequences <- as.character(stringSet)
  
  df <- kmer_assessment(sequences, window_vec = window_vec)
  
  data.frame(df, rel_path = rel_path_f)
  
}

shannon_by_contigs(f[10])

shannon_results <- dplyr::bind_rows(lapply(f, shannon_by_contigs)) 

# Merge shannon with Precision results

DataViz <- Benchmark_results %>%
  ungroup() %>%
  # select(rel_path, sampling_set, Assembly, Accuracy, Precision, Sensitivity) %>%
  left_join(shannon_results, by = "rel_path") 


# DATA VIZ
# 
# 
DataViz %>%
  ungroup() %>% 
  select(-window) %>%
  select_if(is.double) %>%
  cor(method = "pearson") -> M

testRes <- corrplot::cor.mtest(M, conf.level = 0.95)
# 
corrplot::corrplot(M, p.mat = testRes$p ,method = "color", type="upper", order = "hclust", insig = "label_sig")

#

DataViz %>%
  ggplot(aes(x = as.factor(sampling_set ), y = kmer_entropy, group = as.factor(window), color = as.factor(window))) +
  geom_jitter(position = position_jitter(0.1), shape = 1, size = 2, alpha = 0.7) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun.data=mean_sdl, geom="pointrange", size = 0.5, alpha = 0.5) +
  my_custom_theme() +
  scale_color_manual("kmer sizes",values = scale_col) +
  labs(y = "Shannon entropy", x = "Sample size")

discrete_scale <- DataViz %>% distinct(sampling_set) %>% pull()

n <- length(discrete_scale)

scale_col <- c(ggsci::pal_startrek()(7), ggsci::pal_cosmic()(n-7))

scale_col <- structure(scale_col, names = sort(discrete_scale))



DataViz %>%
  ggplot(aes(Ratio, kmer_entropy)) +
  # facet_wrap(~ sampling_set, scales = "free") +
  geom_point(aes(color= as.factor(sampling_set )), shape = 12, size = 2.5) +
  geom_smooth(method = lm, se = FALSE) +
  ggpubr::stat_cor(method = "pearson", 
                   cor.coef.name = "R", color = "black", family = "GillSans",
                   p.accuracy = 0.01) +
  my_custom_theme() +
  scale_color_manual("sample sizes",values = scale_col) 

  


# B)
# DOING for Spades outputs
# 
# 
dir <- "C://Users//cinai/OneDrive/Documentos/GitHub/ConotoxinBenchmark/INPUTS/2_subsampling_dir/Spades_dir/transrate_contigs_dir/"

str(file_list <- list.files(path = dir, pattern = "contigs.csv", recursive = T, full.names = TRUE))

transratedf <- lapply(file_list, read_transrate_scores)

transratedf <- do.call(rbind,transratedf)

transratedf %>% distinct(rel_path, vfold_set, Prop)

Benchmark_results <- transratedf %>%
  filter(length>=200) %>%
  group_by(rel_path, vfold_set, Prop) %>%
  calculate_metrics(reference_coverage_val = 0.95) %>% 
  mutate(
    Ratio = TP/FP,
    # Tell us what percentage of positive classes were correctly identified
    Accuracy = TP / (TP + FN + FP),
    Precision = TP /(TP + FP),
    Sensitivity = TP /(TP + FN),
    Fscore = 2 * (TP) / (2 * (TP) + FP + FN), 
  ) 

Benchmark_results %>% distinct(rel_path)

shannon_by_contigs(f[10])

shannon_results <- dplyr::bind_rows(lapply(f, shannon_by_contigs)) 

outdir <- "C://Users//cinai/OneDrive/Documentos/GitHub/ConotoxinBenchmark/INPUTS/Subsampling_dir/"

f <- list.files(file.path(outdir, "Trinity_FASTA_DIR"), full.names = T)

