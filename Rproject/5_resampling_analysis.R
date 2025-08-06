# Utilizar técnicas de muestreo probabilístico, especialmente el muestreo aleatorio, que da igualdad de oportunidades a todos los individuos de la población.

# Emplear muestreo estratificado, dividiendo la población en subgrupos homogéneos y tomando muestras proporcionales para asegurar la representación adecuada de cada grupo.


# scp -r rgomez@omica.cicese.mx:/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/fernando_pub/2_subsampling_dir/2_transrate_contigs_dir .

# find  boostrap_dir_[0-9]/all_superfamily_FASTA_DIR/ -name "*.fa"
# Rerun transrating, omiting use paired-end evaluation, using only reference metrics

# find  boostrap_dir_[0-9]/all_superfamily_FASTA_DIR/ -name "*.fa"

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

read_transrate_scores <- function(file_list) {
  
  # dir <- dirname(dir)
  
  cat("\nReading\n")
  cat(str_remove(file_list, paste0("^", dir, "/")))
  cat("\n")
  
  read_csv(file_list) %>%
    # mutate(file_list = file_list)
    mutate(rel_path = str_remove(file_list, paste0("^", dirname(dir), "/"))) # %>%
    # separate(rel_path, into = c("subdir1", "subdir2", "filename"), sep = "/", extra = "merge") %>%
    # select(-filename)
}

calculate_metrics <- function(df, reference_coverage_val = 1) {
  
  calculate_false <- function(df) {

    InputNsequences <- c(
      all=1835,
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
    
    
    # False Negatives (FN): Transcripts present in the simulated data but not assembled. 
    ## TP - (N reference sequences in InputNsequences) 
    
    
    data.frame(InputNsequences) %>% 
      as_tibble(rownames = "Superfamily") %>% 
      right_join(df) %>%
      mutate(FN = abs(InputNsequences - TP))
    
  }
  

  Totaldf <- df %>% 
    dplyr::count(subdir1, subdir2, Superfamily)  %>% 
    dplyr::rename("rawcontigs" = "n")
  
  # True Positives (TP): Transcripts correctly assembled by the assembler. 
  ## TP = reference_cov >= reference_coverage_val ||reference_cov == 1 (= !is.na(hits))
  
  TP <- df %>%
    filter(reference_coverage >= reference_coverage_val) %>%
    # Use distinct to trim redundancy between contig hits
    # distinct(hits, Superfamily, subdir1, subdir2) %>%
    dplyr::count(Superfamily, subdir1, subdir2) %>% 
    dplyr::rename("TP" = "n")
  
  # False Positives (FP): Transcripts incorrectly assembled by the assembler.
  ## FP = N contig_name where reference_cov < 1 BUT have some identity threshold (reference_coverage > 0)
  
  FP <- df %>%
    mutate(reference_coverage = ifelse(is.na(hits) & is.na(reference_coverage), 0, reference_coverage )) %>%
    filter(reference_coverage > 0 &  reference_coverage < reference_coverage_val)  %>%
    dplyr::count(Superfamily, subdir1, subdir2) %>% 
    dplyr::rename("FP" = "n")
  

  # Dealing with Overestimate contig number
  # contigs where reference_coverage == 0, OR
  # Use Rawcontigs number minus TP + FP to count potential True Negative (TN)
  #     mutate(TN = rawcontigs - (TP+FP)) %>%
  
  TN <- df %>%
    mutate(reference_coverage = ifelse(is.na(hits) & is.na(reference_coverage), 0, reference_coverage )) %>%
    # mutate_all(~replace(., is.na(.), 0)) %>%
    filter(reference_coverage == 0) %>%
    # filter(is.na(hits)) %>%
    dplyr::count(Superfamily, subdir1, subdir2) %>%
    dplyr::rename("TN" = "n")
  
  
  Totaldf %>% 
    left_join(TP) %>% 
    left_join(TN) %>% 
    left_join(FP) %>% 
    mutate_all(~replace(., is.na(.), 0)) %>%
    left_join(calculate_false(.)) %>%
    select(-InputNsequences) 
  
  
}

# Reading conoServer info

outdir <- "~/Documents/GitHub/ConotoxinBenchmark/INPUTS/"

f <- list.files(path = outdir, pattern = "curated_nuc_conoServerDB.rds", full.names = T)

conoServerDB <- read_rds(f) %>% dplyr::rename("hits" = "entry_id")

# Reading data from 2_subsampling_dir

dir <- "~/Documents/GitHub/ConotoxinBenchmark/2_subsampling_dir/transrate_contigs_dir/"

dir <- dirname(dir)

str(file_list <- list.files(path = dir, pattern = "contigs.csv", recursive = T, full.names = TRUE))


file_list <- file_list[grepl("M_superfamily",file_list)]

# tibble(path = file_list) 

# file_list <- fs::dir_ls(path = dir, recurse = TRUE, regexp = "contigs.csv$")


# command <- paste0("find ", dir, " -name contigs.csv")

# str(file_list <- system(command, intern = TRUE))

transratedf <- lapply(file_list, read_transrate_scores)

transratedf <- do.call(rbind,transratedf)

# Review not duplicates between files <-----

# transratedf %>%  count(rel_path) %>% view()

transratedf <- transratedf %>%
  # select(rel_path) %>%
  mutate(subdir1 = sapply(strsplit(rel_path, "/"), `[`, 3)) %>%
  mutate(Superfamily = sapply(strsplit(rel_path, "/"), `[`, 4)) %>%
  mutate(subdir2 = sapply(strsplit(rel_path, "/"), `[`, 5)) %>%
  # select(-rel_path) %>%
  # mutate(subdir1 = gsub("_superfamily_0.[0-9]_dir","", subdir1)) %>%
  mutate(Superfamily = gsub("_superfamily_FASTA_DIR","", Superfamily)) %>%
  mutate(subdir2 = as.double(gsub("_dir","", subdir2)))


transratedf <- transratedf %>%  
  mutate(hits = sapply(strsplit(hits, "_"), `[`, 1)) 


transratedf %>% dplyr::count(subdir1, subdir2, Superfamily)

# transratedf %>% distinct(hits) 

# Merge with conoServerDB


transratedf %>% distinct(hits) %>% nrow()

# all_contigs_db <- transratedf %>% filter(is.na(hits)) 

sum(sort(conoServerDB$hits) %in% sort(unique(transratedf$hits))) 

# Use right_join to record the reference sequence where assemblers does not support assembly

transratedf <- transratedf %>% 
  # select(contig_name, linguistic_complexity_6, reference_coverage, hits, p_good)
  right_join(conoServerDB,by = "hits")

# Recall test


metricsdf <- calculate_metrics(transratedf, reference_coverage_val = ) %>% 
  mutate(
    Ratio = TP/FP,
    Accuracy = TP / (TP + FN + FP),
    Precision = TP /(TP + FP),
    Recall = TP /(TP + FN),
    Fscore = 2 * (TP) / (2 * (TP) + FP + FN), 
  )


metricsdf %>%
  drop_na() %>%
  ggplot(aes(y = Fscore, x = as.factor(subdir2))) +
  # geom_boxplot() +
  facet_grid(~ Superfamily) +
  geom_jitter(position = position_jitter(0.1), shape = 1) +
  stat_summary(fun = "mean", geom = "line", aes(group = 1), color="red") +
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="red") +
  labs(x = "Subsampling") +
  theme_bw(base_family = "GillSans", base_size = 15) +
  theme(legend.position = "top",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank())


# What is the number of representants per subsampling?
# Split by genesuperfamily

transratedf %>%
  dplyr::count(subdir1, subdir2, genesuperfamily) %>%
  drop_na() %>%
  ggplot(aes(y = n, x = as.factor(subdir2))) +
  facet_grid(genesuperfamily ~.) + geom_boxplot() +
  theme_bw(base_family = "GillSans", base_size = 12) +
  theme(legend.position = "top",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1,size = 7))
  

transratedf %>%
  dplyr::count(subdir1, subdir2, genesuperfamily) %>%
  drop_na() %>%
  ggplot(aes(y = n, x = as.factor(subdir2))) +
  facet_wrap( ~genesuperfamily, scales = "free_y") +
  geom_point()
