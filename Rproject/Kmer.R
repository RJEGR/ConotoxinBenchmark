
# Read transrate scores (contigs.txt) using read_transrate_scores() (to evals accuracy of assembly methods)
# Output transrateDB containing follow columns
# cd /Users/cigom/Documents/GitHub/ConotoxinBenchmark/3_kmer_dir
# scp -r rgomez@omica.cicese.mx:/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/fernando_pub/3_kmer_dir/transrate_contigs_dir .

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

my_custom_theme <- function(...) {
  base_size = 14
  theme_bw(base_family = "GillSans", base_size = base_size) +
    theme(legend.position = "top",
      strip.placement = "outside",
      strip.background = element_rect(fill = 'gray90', color = 'white'),
      strip.text = element_text(angle = 0, size = base_size, hjust = 0), 
      axis.text = element_text(size = rel(0.7), color = "black"),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.x = element_blank(),
      ...
    )
}

read_transrate_scores <- function(file_list) {
  
  # dir <- dirname(dir)
  rel_path <- str_remove(file_list, paste0("^", dir, "/"))
  
  cat("\nReading\n")
  cat(rel_path)
  cat("\n")
  
  rel_path <- gsub("transrate_contigs_dir/", "", rel_path)
  
  rel_path <- sapply(strsplit(dirname(rel_path), "/"), `[`, 1)
  vfold_set <- sapply(strsplit(rel_path, "_"), `[`, 1)
  kmer <- sapply(strsplit(rel_path, "_"), `[`, 5)
  
  cat("\n")
  cat(vfold_set, kmer)
  cat("\n")
  
  read_csv(file_list) %>%
    # mutate(file_list = file_list)
    mutate(rel_path, vfold_set, kmer)
}

calculate_metrics <- function(df, reference_coverage_val = 1) {
  
  is_chimeric_value = 0
  
  calculate_false <- function(df) {
    
    # as many assemblers use width 200 to filter contigs, count number of refseq > 200
    # InputNsequences <- c(
    #   Fold01=1615,
    #   Fold02=1614,
    #   Fold03=1621,
    #   Fold04=1618,
    #   Fold05=1616,
    #   Fold06=1619,
    #   Fold07=1614,
    #   Fold08=1617,
    #   Fold09=1616,
    #   Fold10=1614,
    #   Fold11=1617,
    #   Fold12=1619)
    
    
    count_Nsequences <- function() {
      
      # Temporal directory where vfolds_resampling_dir are found
      
      outdir <- "~/Documents/GitHub/ConotoxinBenchmark/INPUTS/vfolds_resampling_dir/"
      
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
    # dplyr::count(vfold_set, sampling_set)  %>% 
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


outdir <- "~/Documents/GitHub/ConotoxinBenchmark/INPUTS/"

f <- list.files(path = outdir, pattern = "curated_nuc_conoServerDB.rds", full.names = T)

conoServerDB <- read_rds(f) %>% dplyr::rename("hits" = "entry_id")

dir <- "/Users/cigom/Documents/GitHub/ConotoxinBenchmark/3_kmer_dir/"

str(file_list <- list.files(path = dir, pattern = "contigs.csv", recursive = T, full.names = TRUE))


transratedf <- lapply(file_list, read_transrate_scores)

transratedf <- do.call(rbind,transratedf)

transratedf %>%
  group_by(vfold_set, kmer) %>%
  dplyr::count()  

transratedf %>%
  group_by(kmer) %>%
  summarise(mean = mean(length), min = min(length), max = max(length))


metricsdf <- transratedf %>% 
  filter(length>=200) %>%
  group_by(vfold_set, kmer) %>%
  calculate_metrics(reference_coverage_val = 0.95) %>% 
  mutate(
    Ratio = TP/FP,
    # Tell us what percentage of positive classes were correctly identified
    Accuracy = TP / (TP + FN + FP),
    Precision = TP /(TP + FP),
    Sensitivity = TP /(TP + FN),
    Fscore = 2 * (TP) / (2 * (TP) + FP + FN), 
  ) 


cols_to <- c("Accuracy", "Precision", "Sensitivity")

recode_to <- structure(c("A) Accuracy", "B) Precision", "C) Sensitivity"), names = cols_to)

base_size <- 14

subti <- "Larger k-mers are more unique\nand are crucial for resolving repetitive regions,\nleading to a more contiguous assembly\nSmaller k-mers are essential\nfor capturing regions of low coverage,\nbut they perform poorly in repetitive regions."
  
  # Smaller k-mers are essential for capturing regions of low coverage, but they perform poorly in repetitive regions.
  

p1 <- metricsdf %>%
  select(-TP, -FP, -FN, -rawcontigs) %>%
  pivot_longer(cols = cols_to, values_to = "y", names_to = "facet") %>%
  dplyr::mutate(facet = dplyr::recode_factor(facet, !!!recode_to)) %>%
  drop_na() %>%
  ggplot(aes(y = y, x = as.factor(kmer))) +
  facet_grid(facet ~., scales ="free", switch = "y") +
  geom_jitter(position = position_jitter(0.1), shape = 1) +
  stat_summary(fun = "mean", geom = "line", aes(group = 1), color="blue") +
  stat_summary(fun = "mean", geom = "point", aes(group = 1), color="blue") + # , color="blue"
  # stat_summary(fun.data=mean_sdl, geom="pointrange", color="red") +
  labs(x = "k-mer size", y = "", caption = "kmer.R (Assembly with spades)", tag = subti) +
  my_custom_theme(plot.tag = element_text(size = 7, hjust = 0))

p1

ggsave(p1, filename = 'kmer_analysis.png', 
  path = outdir, width = 4.5, height = 6, dpi = 1000, device = png)


quit()




