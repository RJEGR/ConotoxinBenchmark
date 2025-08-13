

# Read transrateDB (to evals accuracy of assembly methods)

# True Positives (TP): Transcripts correctly assembled by the assembler. 
## TP = reference_cov == 1
# False Positives (FP): Transcripts incorrectly assembled by the assembler.
## FP = contig_name where is.na(hits) == TRUE, and reference_cov < 1
# False Negatives (FN): Transcripts present in the simulated data but not assembled. 
## TP - (N reference sequences in InputNsequences) 
# Sensitivity: The proportion of true transcripts that are correctly assembled. 
# Precision: The proportion of assembled transcripts that are actually true. 
# F1-score (Recall): A harmonic mean of precision and sensitivity, providing an overall measure of assembly quality. 


rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)


calculate_metrics <- function(df, reference_coverage_val = 1) {
  
  calculate_false <- function(df) {
    
    # False Negatives (FN): Transcripts present in the simulated data but not assembled. 
    ## TP - (N reference sequences in InputNsequences) 
    
    InputNsequences <- c(all=1835,
      A=251,
      # Conopeptides.fasta:123
      D=18,
      I1=22,
      I2=49,
      I3=10,
      Insulin=23,
      J=28,
      L=11,
      M=445,
      O1=356,
      O2=123,
      O3=36,
      P=15,
      Q=20,
      S=11,
      `T`=215,
      UNDER=78)
    
    data.frame(InputNsequences) %>% 
      as_tibble(rownames = "Superfamily") %>% 
      left_join(df) %>%
      mutate(FN = abs(InputNsequences - TP))
    
  }
  
  Totaldf <- df %>% 
    count(Superfamily, Assembler)  %>% 
    dplyr::rename("rawcontigs" = "n")
  
  # True Positives (TP): Transcripts correctly assembled by the assembler. 
  ## TP = reference_cov >= reference_coverage_val
  
  TP <- df %>%
    filter(reference_coverage >= reference_coverage_val) %>%
    # Use distinct to trim redundancy
    distinct(hits, Superfamily, Assembler) %>%
    count(Superfamily, Assembler) %>% 
    dplyr::rename("TP" = "n")
  
  # False Positives (FP): Transcripts incorrectly assembled by the assembler.
  ## FP = N contig_name where is.na(hits) == TRUE, AND|OR reference_cov < 1
  
  FP <- df %>%
    mutate(reference_coverage = ifelse(is.na(hits) & is.na(reference_coverage), 0,reference_coverage )) %>%
    filter(reference_coverage < reference_coverage_val) %>%
    count(Superfamily, Assembler) %>% 
    dplyr::rename("FP" = "n")
  
  # FPnohit <- df %>%
  #   # filter(reference_coverage < 1) %>%
  #   drop_na(hits) %>%
  #   count(Superfamily, Assembler) %>% 
  #   dplyr::rename("FPnohit" = "n")
  
  Totaldf %>% 
    left_join(TP) %>% 
    left_join(FP) %>% 
    # left_join(FPnohit) %>%
    mutate_all(~replace(., is.na(.), 0)) %>%
    left_join(calculate_false(.))
  
  
  
}

outdir <- "~/Documents/GitHub/ConotoxinBenchmark/INPUTS/"

# read_rds(file.path(outdir, "transrate_assemblers.rds"))

assembly_df <- read_rds(file.path(outdir, "transrate_assemblers_true_contigs.rds"))


metricsdf <- calculate_metrics(assembly_df, reference_coverage_val = 0.7) %>% 
  mutate(
    # RatioCI = TP/FP,
    # Sensitivity: The proportion of true transcripts that are correctly assembled. 
    # Sensitivity = TP / (TP + FP), #<-- validate formula
    Sensitivity = TP / (TP + FN),
    # Precision: The proportion of assembled transcripts that are actually true. 
    # Precision = TP /(TP + FN),  #<-- validate formula
    Precision = TP /(TP + FP),
    # F1-score (aka Recall): A harmonic mean of precision and sensitivity, providing an overall measure of assembly quality. 
    
    Recall = 2 * (Precision * Sensitivity) / (Precision + Sensitivity),  #<-- validate formula
    # F1score = 2 * (TP) / (2 * (TP) + FP + FN) # <-- Recall and F1score formula is the same
  )


# In reference_coverage_val, try a loop to itinerate evaluation of recall


reference_coverage_seq <- seq(0.1, 1, by = 0.1)

metrics_list <- lapply(reference_coverage_seq, function(i) {
  calculate_metrics(assembly_df, reference_coverage_val = i) %>%
    mutate(
      Sensitivity = TP / (TP + FN),
      Precision = TP / (TP + FP),
      Recall = 2 * (Precision * Sensitivity) / (Precision + Sensitivity),
      threshold = i
    )
})

metricsdf <- bind_rows(metrics_list)


assembly_df %>%
  mutate(col = NA) %>%
  filter(Assembler != "VELVET") %>%
  mutate(reference_coverage = ifelse(is.na(hits) & is.na(reference_coverage), 0,reference_coverage )) %>%
  mutate(col = ifelse(reference_coverage >= 0.50, "≥ 50% alignment", col)) %>%
  mutate(col = ifelse(reference_coverage >= 0.70, "≥ 70% alignment", col)) %>%
  mutate(col = ifelse(reference_coverage >= 0.80, "≥ 80% alignment", col)) %>%
  mutate(col = ifelse(reference_coverage >= 0.90, "≥ 90% alignment", col)) %>%
  mutate(col = ifelse(reference_coverage >= 0.95, "≥ 95% alignment", col)) %>%
  mutate(col = ifelse(reference_coverage == 1, "≥ 100% alignment", col)) %>%
  count(Assembler, col) %>%
  drop_na() %>%
  ggplot(aes(y = n, x = Assembler, fill = col)) +
  ggforce::facet_col(col ~., scales = "free_y") +
  geom_col(position = position_dodge2())


metricsdf %>% 
  select(threshold, Recall, Assembler, Superfamily) %>%
  ggplot(aes(x = as.factor(threshold), y = Recall)) +
  facet_grid(~Assembler) +
  geom_boxplot()

