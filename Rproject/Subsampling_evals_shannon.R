# As we found, as kmer size increase, the shannon (complexity) entropy increaseo, and
# # if sample size increase, also the shannon entropy, in the conoServerDB
# # lets calculate the shannon values for those good contigs assembled in the Subsampling step
# 
# /LUSTRE/bioinformatica_data/genomica_funcional/rgomez/fernando_pub/2_subsampling_dir/multiple_analysis_dir/Trinity_FASTA_DIR
# 
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
  theme_bw(base_family = "GillSansMT", base_size = base_size) +
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
    # 
    # See entropy-based.R
    # pspm_uniform <- matrix(0.25, nrow = window, ncol = 4, dimnames = list(NULL, c("A","C","G","T")))
    
    # perplexities <- lapply(kmers, calculate_perplexity, pspm = pspm_uniform)
    
    
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


# k <- round(log(sum(Biostrings::width(stringSet)), base = 4) )

kmer_vec <- c(19, 21, 25, 33, 49, 55, 73)

shannon_by_contigs <- function(f, window_vec = kmer_vec) {
  
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

# Step 1
# Starting with TRINITY outputs ======
# dir <- "C://Users//cinai/OneDrive/Documentos/GitHub/ConotoxinBenchmark/INPUTS/2_subsampling_dir/Trinity_dir/transrate_contigs_dir/"
 

dir <- "C://Users//cinai/OneDrive/Documentos/GitHub/ConotoxinBenchmark/INPUTS/2_subsampling_dir/Trinity_dir/transrate_contigs_dir/"

str(file_list <- list.files(path = dir, pattern = "contigs.csv", recursive = T, full.names = TRUE))

transratedf <- lapply(file_list, read_transrate_scores)

transratedf <- do.call(rbind,transratedf)

transratedf %>% distinct(rel_path, vfold_set, Prop)

.transratedf_trinity <- transratedf

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

# Step 2: LOAD DATA from sequences in TRINITY
# 
# 

outdir <- "C://Users//cinai/OneDrive/Documentos/GitHub/ConotoxinBenchmark/INPUTS/Subsampling_dir/"

f <- list.files(file.path(outdir, "Trinity_FASTA_DIR"), full.names = T)

shannon_by_contigs(f[10])

shannon_results <- dplyr::bind_rows(lapply(f, shannon_by_contigs)) 

DataViz <- Benchmark_results %>%
  ungroup() %>%
  # select(rel_path, sampling_set, Assembly, Accuracy, Precision, Sensitivity) %>%
  left_join(shannon_results, by = "rel_path") %>%
  mutate(Assembly = "Trinity") 

# DOING for Spades outputs
# Step 3
# 
dir <- "C://Users//cinai/OneDrive/Documentos/GitHub/ConotoxinBenchmark/INPUTS/2_subsampling_dir/Spades_dir/transrate_contigs_dir/"

str(file_list <- list.files(path = dir, pattern = "contigs.csv", recursive = T, full.names = TRUE))

transratedf <- lapply(file_list, read_transrate_scores)

transratedf <- do.call(rbind,transratedf)

.transratedf <- transratedf %>%
  mutate(Assembler = "Spades") %>%
  rbind(mutate(.transratedf_trinity, Assembler = "Trinity"))


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

# Step 4: LOAD DATA from sequences in SPADES
# 
# 


outdir <- "//wsl.localhost/Debian/home/ricardo/"

f <- list.files(file.path(outdir, "Spades_FASTA_DIR"), full.names = T)

shannon_results_2 <- dplyr::bind_rows(lapply(f, shannon_by_contigs)) 

# Merge shannon with Precision results

DataViz <- Benchmark_results %>%
  ungroup() %>%
  # select(rel_path, sampling_set, Assembly, Accuracy, Precision, Sensitivity) %>%
  left_join(shannon_results_2, by = "rel_path") %>%
  mutate(Assembly = "Spades") %>%
  rbind(DataViz)

outdir <- "C://Users//cinai/OneDrive/Documentos/GitHub/ConotoxinBenchmark/INPUTS/"


# write_tsv(DataViz, file = file.path(outdir, "subsampling_shannon_benchmark.tsv"))

DataViz <- read_tsv(file.path(outdir, "subsampling_shannon_benchmark.tsv"))

# DATA VIZ
# 
# 
# 



DataViz %>%
  ungroup() %>% 
  select(-window) %>%
  select_if(is.double) %>%
  cor(method = "pearson") -> M

testRes <- corrplot::cor.mtest(M, conf.level = 0.95)
 
corrplot::corrplot(M, p.mat = testRes$p ,method = "color", type="upper", order = "hclust", insig = "label_sig")

#



DataViz %>%
  ggplot(aes(Precision , kmer_entropy)) +
  facet_grid(window~ Assembly, scales = "free") +
  geom_point(aes(color= as.factor(Prop   )), shape = 12, size = 2.5) +
  geom_smooth(method = lm, se = FALSE) +
  ggpubr::stat_cor(method = "pearson", 
                   cor.coef.name = "R", color = "black", family = "GillSans",
                   p.accuracy = 0.01) +
  my_custom_theme() 


discrete_scale <- DataViz %>% distinct(window) %>% pull()

n <- length(discrete_scale)

scale_col <- c(ggsci::pal_startrek()(7), ggsci::pal_cosmic()(n-7))

scale_col <- structure(scale_col, names = sort(discrete_scale))

# Accuracy Precision Sensitivity Fscore
# 
DataViz %>%
  ggplot(aes(Sensitivity , kmer_entropy, color= as.factor(window))) +
  facet_grid(~ Assembly, scales = "free") +
  geom_point(shape = 12, size = 2.5) +
  geom_smooth(method = lm, se = FALSE) +
  ggpubr::stat_cor(method = "pearson", 
                   cor.coef.name = "R", color = "black", family = "GillSans",
                   p.accuracy = 0.01) +
  my_custom_theme()  +
  scale_color_manual("kmer sizes",values = scale_col) 

discrete_scale <- DataViz %>% distinct(Prop) %>% pull()

n <- length(discrete_scale)

scale_col <- c(ggsci::pal_startrek()(7), ggsci::pal_cosmic()(n-7))

scale_col <- structure(scale_col, names = sort(discrete_scale))

DataViz %>%
  ggplot(aes(x = as.factor(Prop), y = kmer_entropy, group = as.factor(window), color = as.factor(window))) +
  geom_jitter(position = position_jitter(0.1), shape = 1, size = 2, alpha = 0.7) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun.data=mean_sdl, geom="pointrange", size = 0.5, alpha = 0.5) +
  my_custom_theme() +
  scale_color_manual("kmer sizes",values = scale_col) +
  labs(y = "Shannon entropy", x = "Sample size")

# Now step 5, evaluates entropy of refereces that were not possible to assembly, 
# And figure out how complexity confirm sensitivity of assemblers, 
# 
# 

DataViz %>%
  distinct(rel_path,vfold_set, Prop, kmer_entropy, window) %>%
  filter(window == 33 & Prop == 1) %>% 
  mutate(y = "Assembled") %>%
  ggplot(aes(y = y, x = kmer_entropy)) +
  geom_jitter(position = position_jitter(0.1), shape = 1, size = 2, alpha = 0.7) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun.data=mean_sdl, geom="pointrange", size = 0.5, alpha = 0.5) +
  # scale_color_manual("kmer sizes",values = scale_col) +
  my_custom_theme()


# outdir <- "C://Users//cinai/OneDrive/Documentos/GitHub/ConotoxinBenchmark/INPUTS//"

# f <- list.files(file.path(outdir, "vfolds_resampling_dir"), full.names = T)

outdir <- "//wsl.localhost/Debian/home/ricardo/"

f <- list.files(file.path(outdir, "Spades_FASTA_DIR"), full.names = T)

# shannon_results_2 <- dplyr::bind_rows(lapply(f, shannon_by_contigs)) 


shannon_by_contigs2 <- function(f, window_vec = kmer_vec) {
  
  # Esto funciona perfectamente, pero el n_seqs no esta equilibrado, 
  # y entonces los calculos se sesgan al tamano de la muestra mas
  #  mas que a una razon explicativa de la complejidad

  rel_path_f <- gsub(".fa|.fasta", "",basename(f))
  
  cat("\n Reading: ", rel_path_f, " file\n")
  
  splitted_seqs <- transratedf %>%
    filter(length >= 200) %>%
    filter(rel_path == rel_path_f) %>%
    mutate(summarise = "< 80 % alignment") %>%
    mutate(summarise = ifelse(is.na(hits) & is.na(reference_coverage), "Nohit", summarise )) %>%
    mutate(reference_coverage = ifelse(is.na(hits) & is.na(reference_coverage), 0, reference_coverage )) %>%
    mutate(summarise = ifelse(reference_coverage >= 0.8, ">= 80% alignment", summarise)) %>%
    mutate(summarise = ifelse(reference_coverage >= 0.9, ">= 90% alignment", summarise)) %>%
    mutate(summarise = ifelse(reference_coverage >= 0.95, ">= 95% alignment", summarise)) %>%
    mutate(summarise = ifelse(reference_coverage == 1, "100% alignment", summarise)) %>%
    distinct(contig_name, summarise) 
  
  
  # classify_seqs %>% 
  #   filter(!summarise %in% c(">= 95% alignment", "100% alignment")) %>%
  #   ggplot(aes(linguistic_complexity_6, color = summarise)) + 
  #   geom_density()
  # 
  splitted_seqs <- splitted_seqs %>% 
    filter(!summarise %in% c(">= 95% alignment", "100% alignment")) %>%
    mutate(summarise = "All < 95% alignment") %>%
    rbind(splitted_seqs)
    
    
    
  splitted_seqs <- split(strsplit(splitted_seqs$contig_name, ";") , splitted_seqs$summarise)
  
  splitted_seqs <- lapply(splitted_seqs, unlist)
  
  
  
  kmer_splits <- function(keep_stringSet, window_vec = window_vec) {
    
    # this function use by default the transratedf DB to track the good assembled contigs 
    
    stringSet <- Biostrings::readAAStringSet(f)
    
    contig_names_l <- sapply(strsplit(names(stringSet), " "), `[`, 1)
    
    stringSet <- stringSet[contig_names_l %in% keep_stringSet]  
    
    sequences <- as.character(stringSet)
    
    df <- kmer_assessment(sequences, window_vec = window_vec)
    
    data.frame(df, rel_path = rel_path_f)
    
  }
  
  seq_results <- lapply(splitted_seqs, kmer_splits, window_vec = kmer_vec)
  
  # Combine results to a single data.frame if desired:
  seq_results <- do.call(rbind, lapply(names(seq_results), function(nm) {
    df <- seq_results[[nm]]
    df$splitted_seqs <- nm
    df
  }))
  
  return(seq_results)
}


DataViz <- shannon_by_contigs2(f[1])

# Evenness enabling comparability across samples and window sizes

DataViz %>%
  # distinct(rel_path,vfold_set, Prop, kmer_entropy, window) %>%
  # filter(window == 33 & Prop == 1) %>% 
  # mutate(y = "Assembled") %>%
  ggplot(aes(y = splitted_seqs, x = evenness)) +
  geom_jitter(position = position_jitter(0.1), shape = 1, size = 2, alpha = 0.7) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun.data=mean_sdl, geom="pointrange", size = 0.5, alpha = 0.5) +
  # scale_color_manual("kmer sizes",values = scale_col) +
  my_custom_theme()

seq_results <- dplyr::bind_rows(lapply(f, shannon_by_contigs2))


# Step 6 evals on reference
# 


outdir <- "C://Users//cinai/OneDrive/Documentos/GitHub/ConotoxinBenchmark/INPUTS"

file_out <- file.path(outdir, "curated_nuc_conoServerDB.rds")

data <- read_rds(file_out)


.transratedf %>% 
  ggplot(aes(y = Prop , x = reference_coverage, fill = Assembler, alpha = after_stat(x))) +
  # geom_violin() +
  # facet_grid(~ Assembler) +
  ggridges::geom_density_ridges_gradient(
    jittered_points = F,
    position = ggridges::position_points_jitter(width = 0.05, height = 0),
    point_shape = '|', point_size = 3, point_alpha = 1, alpha = 1)

queries_seqs <- .transratedf %>%
  # drop_na(hits) %>%
  filter(reference_coverage >= 0.95) %>%
  distinct(hits, Assembler) %>%
  dplyr::count(hits, sort = T) %>%
  pull(hits)

length(queries_seqs)

sum(data$entry_id %in% queries_seqs)

splitted_seqs <- data %>% 
  mutate(splitted_seqs = "Not Assembled") %>%
  mutate(splitted_seqs = ifelse(entry_id %in% queries_seqs, ">= 95% alignment", splitted_seqs) ) %>%
  distinct(entry_id, splitted_seqs, sequence)


splitted_seqs <- split(strsplit(splitted_seqs$sequence, ";") , splitted_seqs$splitted_seqs    )

splitted_seqs <- lapply(splitted_seqs, unlist)

str(splitted_seqs)

kmer_vec <- c(19, 21, 25, 33, 49, 55, 73)


split_apply <- function(splitted_seqs) {
  
  seq_results <- lapply(splitted_seqs, kmer_assessment, window_vec = kmer_vec)
  
  seq_results <- do.call(rbind, lapply(names(seq_results), function(nm) {
    df <- seq_results[[nm]]
    df$splitted_seqs <- nm
    df
  }))
  
  return(seq_results)
  
}

seq_results <- split_apply(splitted_seqs)

seq_results %>%
  ggplot(aes(y = splitted_seqs, x = kmer_entropy )) +
  geom_jitter(position = position_jitter(0.1), shape = 1, size = 2, alpha = 0.7) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun.data=mean_sdl, geom="pointrange", size = 0.5, alpha = 0.5) +
  # scale_color_manual("kmer sizes",values = scale_col) +
  my_custom_theme()


####
####
#### 

data %>% 
  mutate(splitted_seqs = "Not Assembled") %>%
  mutate(splitted_seqs = ifelse(entry_id %in% queries_seqs, ">= 95% alignment", splitted_seqs) ) %>%
  distinct(entry_id, splitted_seqs, sequence)

stratify_data <- function(s, window_vec = c(19, 21, 25, 33, 49, 55, 73)) {
  
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

data_notassembled <- data %>%  
  mutate(splitted_seqs = "Not Assembled") %>%
  mutate(splitted_seqs = ifelse(entry_id %in% queries_seqs, ">= 95% alignment", splitted_seqs) ) %>%
  distinct(entry_id, splitted_seqs, sequence) %>% 
  filter(splitted_seqs == "Not Assembled")

folds <- rsample::vfold_cv(data_notassembled, v = 12) # strata = stratified_sampling

require(rsample)

seq_results <- dplyr::bind_rows(lapply(folds$splits, stratify_data)) 

# plot
# 
# 
# 
# 
seq_results %>%
  ggplot(aes(y = kmer_entropy, x = as.factor(Prop) )) +
  geom_jitter(position = position_jitter(0.1), shape = 1, size = 2, alpha = 0.7) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun.data=mean_sdl, geom="pointrange", size = 0.5, alpha = 0.5) +
  # scale_color_manual("kmer sizes",values = scale_col) +
  my_custom_theme()

# Then
# 
# 

data_assembled <- data %>%  
  mutate(splitted_seqs = "Not Assembled") %>%
  mutate(splitted_seqs = ifelse(entry_id %in% queries_seqs, ">= 95% alignment", splitted_seqs) ) %>%
  distinct(entry_id, splitted_seqs, sequence) %>% 
  filter(splitted_seqs != "Not Assembled")

folds <- rsample::vfold_cv(data_assembled, v = 12) # strata = stratified_sampling

require(rsample)

seq_results2 <- dplyr::bind_rows(lapply(folds$splits, stratify_data)) 

 