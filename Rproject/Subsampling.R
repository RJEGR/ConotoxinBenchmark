#
# After running Subsampling.sh analysis, OPEN contigs.csv files and calculate metrics of accuracy
# EDA of 2_subsampling_dir/2_transrate_contigs_dir


# scp -r rgomez@omica.cicese.mx:/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/fernando_pub/2_subsampling_dir/2_transrate_contigs_dir .

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

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
  
  is_chimeric_value = 0
  
  calculate_false <- function(df) {
    
    InputNsequences <- c(Fold01=1615,
    Fold02=1614,
    Fold03=1621,
    Fold04=1618,
    Fold05=1616,
    Fold06=1619,
    Fold07=1614,
    Fold08=1617,
    Fold09=1616,
    Fold10=1614,
    Fold11=1617,
    Fold12=1619)
    
    # False Negatives (FN): Transcripts present in the simulated data but not assembled. 
    ## TP - (N reference sequences in InputNsequences) 
    
    
    data.frame(InputNsequences) %>% 
      as_tibble(rownames = "vfold_set") %>% 
      right_join(df) %>%
      mutate(FN = abs(InputNsequences - TP))
    
  }
  

  Totaldf <- df %>% 
    dplyr::count(vfold_set, sampling_set)  %>% 
    dplyr::rename("rawcontigs" = "n")
  
  # True Positives (TP): Transcripts correctly assembled by the assembler. 
  ## TP = reference_cov >= reference_coverage_val ||reference_cov == 1 (= !is.na(hits))
  
  TP <- df %>%
    filter(reference_coverage >= reference_coverage_val) %>%
    # Use distinct to trim redundancy between contig hits
    # distinct(hits, Superfamily, subdir1, subdir2) %>%
    dplyr::count(vfold_set, sampling_set)  %>% 
    dplyr::rename("TP" = "n")
  
  # False Positives (FP): Transcripts incorrectly assembled by the assembler.
  ## FP = N contig_name where reference_cov < 1 BUT have some identity threshold (reference_coverage > 0)
  
  FP <- df %>%
    mutate(reference_coverage = ifelse(is.na(hits) & is.na(reference_coverage), 0, reference_coverage )) %>%
    filter(reference_coverage > is_chimeric_value &  reference_coverage < reference_coverage_val)  %>%
    dplyr::count(vfold_set, sampling_set) %>% 
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
    dplyr::count(vfold_set, sampling_set) %>%
    dplyr::rename("TN" = "n")
  
  
  Totaldf %>% 
    left_join(TP) %>% 
    # left_join(TN) %>% 
    left_join(FP) %>% 
    mutate_all(~replace(., is.na(.), 0)) %>%
    left_join(calculate_false(.)) %>%
    select(-InputNsequences) 
  
  
}

# Reading conoServer info

outdir <- "~/Documents/GitHub/ConotoxinBenchmark/INPUTS/"

f <- list.files(path = outdir, pattern = "curated_nuc_conoServerDB.rds", full.names = T)

conoServerDB <- read_rds(f) %>% dplyr::rename("hits" = "entry_id")

dir <- "/Users/cigom/Documents/GitHub/ConotoxinBenchmark/2_subsampling_dir/transrate_contigs_dir/"

str(file_list <- list.files(path = dir, pattern = "contigs.csv", recursive = T, full.names = TRUE))

transratedf <- lapply(file_list, read_transrate_scores)

transratedf <- do.call(rbind,transratedf)

transratedf <- transratedf %>% 
  mutate(rel_path = basename(dirname(rel_path))) %>%
  mutate(vfold_set = sapply(strsplit(rel_path, "_"), `[`, 1)) %>%
  mutate(sampling_set = sapply(strsplit(rel_path, "_"), `[`, 5)) %>%
  mutate(sampling_set = as.double(sampling_set))

transratedf %>%
  dplyr::count(rel_path, vfold_set, sampling_set)  


calculate_metrics(transratedf, reference_coverage_val = 0.95) %>%
  write_tsv(file = file.path(dir, "subsampling_benchmark.tsv"))

metricsdf <- calculate_metrics(transratedf, reference_coverage_val = 0.95) %>% 
  mutate(
    Ratio = TP/FP,
    # Tell us what percentage of positive classes were correctly identified
    Accuracy = TP / (TP + FN + FP),
    Precision = TP /(TP + FP),
    Sensitivity = TP /(TP + FN),
    Fscore = 2 * (TP) / (2 * (TP) + FP + FN), 
  ) 


# (Quantitative): Proxy 1

DataViz <- transratedf %>% 
  drop_na() %>% 
  distinct(sampling_set, vfold_set, hits, reference_coverage) %>%
  # left_join(conoServerDB,by = "hits") %>%
  mutate(summarise = "< 80 % alignment") %>%
  mutate(summarise = ifelse(reference_coverage >= 0.8, ">= 80% alignment", summarise)) %>%
  mutate(summarise = ifelse(reference_coverage >= 0.9, ">= 90% alignment", summarise)) %>%
  mutate(summarise = ifelse(reference_coverage >= 0.95, ">= 95% alignment", summarise)) %>%
  mutate(summarise = ifelse(reference_coverage == 1, "100% alignment", summarise)) %>%
  dplyr::count(sampling_set, vfold_set, summarise)


n_pallet <- length(unique(DataViz$summarise))

scale_col <- ggsci::pal_uchicago()(n_pallet) 

scale_col <- structure(scale_col, names = sort(unique(DataViz$summarise)))

scale_fill <- ggsci::pal_uchicago(alpha = 0.5)(n_pallet) 

scale_fill <- structure(scale_fill, names = sort(unique(DataViz$summarise)))


p2 <- DataViz %>%
  # filter(summarise != "< 80 % alignment") %>%
  ggplot(aes(y = n, x = as.factor(sampling_set), color = summarise, fill = summarise)) +
  # geom_jitter(position = position_jitter(0.1), shape = 1) +
  stat_summary(fun.data=mean_sdl, geom="pointrange", shape = 1) + # position = position_jitter(0.25)
  labs(y = "Number of assembled conotoxins", x = "Sample size (Proportion of the sample)", caption = "3_Subsampling.R") +
  my_custom_theme(legend.text = element_text(size = 5)) +
  scale_color_manual("",values = scale_col ) +
  scale_fill_manual("",values = scale_fill)
  
p2

# ggsave(p2,
#     filename = 'Subsampling_boxplot.png', 
#     path = outdir, width = 5.5, height = 5, dpi = 1000, device = png)

  
# (Quantitative): Proxy 2. What is the distribution in genesuperfamily ? Why some sf are easy to assembly?
# Based on entropy avarege per family, diversity of sequences, and cluster similarity (see )
# 
conoServerDB %>%
  dplyr::count(genesuperfamily)
transratedf %>% 
  # filter(sampling_set == 1) %>%
  distinct(sampling_set, vfold_set, hits, reference_coverage) %>%
  left_join(conoServerDB,by = "hits") %>%   drop_na() %>% 
  mutate(summarise = "< 80 % alignment") %>%
  mutate(summarise = ifelse(reference_coverage >= 0.8, ">= 80% alignment", summarise)) %>%
  mutate(summarise = ifelse(reference_coverage >= 0.9, ">= 90% alignment", summarise)) %>%
  mutate(summarise = ifelse(reference_coverage >= 0.95, ">= 95% alignment", summarise)) %>%
  mutate(summarise = ifelse(reference_coverage == 1, "100% alignment", summarise)) %>%
  dplyr::count(sampling_set, vfold_set, summarise, genesuperfamily) %>% #genesuperfamily, 
  # mutate()
  ggplot(aes(y = n, x = as.factor(genesuperfamily), color = summarise, fill = summarise)) +
  ggforce::facet_col(~ summarise) +
  # geom_jitter(position = position_jitter(0.1), shape = 1) +
  stat_summary(fun.data=mean_sdl, geom="pointrange", shape = 1) + # position = position_jitter(0.25)
  labs(y = "Number of assembled conotoxins", x = "Sample size (Proportion of the sample)", caption = "3_Subsampling.R") +
  my_custom_theme(legend.text = element_text(size = 5), 
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1,size = 12)) +
  ggsci::scale_color_uchicago(name = "") +
  ggsci::scale_fill_uchicago(name = "", alpha = 0.5)



# Facet benchmark

cols_to <- c("Accuracy", "Precision", "Sensitivity")

recode_to <- structure(c("A) Accuracy", "B) Precision", "C) Sensitivity"), names = cols_to)

base_size <- 14

p1 <- metricsdf %>%
  select(-TP, -FP, -FN, -rawcontigs) %>%
  pivot_longer(cols = cols_to, values_to = "y", names_to = "facet") %>%
  dplyr::mutate(facet = dplyr::recode_factor(facet, !!!recode_to)) %>%
  drop_na() %>%
  ggplot(aes(y = y, x = as.factor(sampling_set))) +
  facet_grid(facet ~., scales ="free", switch = "y") +
  geom_jitter(position = position_jitter(0.1), shape = 1) +
  stat_summary(fun = "mean", geom = "line", aes(group = 1), color="blue") +
  stat_summary(fun = "mean", geom = "point", aes(group = 1), color="blue") +
  # stat_summary(fun.data=mean_sdl, geom="pointrange", color="red") +
  labs(x = "Sample size (Proportion of the sample)", y = "", caption = "3_Subsampling.R") +
  # ylim(0,NA) +
  my_custom_theme()

# p1

ggsave(p1, filename = 'Subsampling.png', 
  path = outdir, width = 4, height = 6, dpi = 1000, device = png)


quit()




# (Qualitative) what are the conotoxins not assembled?

transratedf %>% drop_na() %>% distinct(hits) %>% nrow()

sum(sort(conoServerDB$hits) %in% sort(unique(transratedf$hits))) 

# All the sequences !!! but what is the degree?

# (Quantitative) What is the category alignment of dsitribution per superfamily)



# Use right_join to record the reference sequence where assemblers does not support assembly

transratedf %>% right_join(conoServerDB,by = "hits")

# select(contig_name, linguistic_complexity_6, reference_coverage, hits, p_good)

transratedf %>% 
  filter(is.na(hits))

# What is the relation between sequence complexity and precision

transratedf %>%
  ggplot(aes(linguistic_complexity_6, reference_coverage)) + geom_point()

transratedf %>%
  dplyr::count(vfold_set, sampling_set, genesuperfamily) %>%
  drop_na() %>%
  ggplot(aes(y = n, x = as.factor(sampling_set))) +
  # facet_grid(genesuperfamily ~.) + geom_boxplot() +
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


conoServerDB  %>%
  drop_na() %>%
  count(organismlatin, organismdiet, genesuperfamily, sort = T) %>% view()
  mutate(genesuperfamily = factor(genesuperfamily, levels = unique(genesuperfamily))) %>%
  ggplot(aes(organismlatin, genesuperfamily, fill = n)) +
  facet_grid(~ organismdiet, scales = "free", space = "free") +
  geom_tile(color = "white", linewidth = 0.5) +
  theme_bw(base_family = "GillSans", base_size = 14) +
  # scale_x_discrete(position = "top") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

conoServerDB %>% count(organismlatin, organismdiet, sort = T) %>%
  ggplot(aes(y = organismlatin,  x = n)) +
  ggforce::facet_col(organismdiet ~., scales = "free", space = "free") +
  geom_col() +
  theme_bw(base_family = "GillSans", base_size = 14)


# library(taxize)
# tax <- names2wormsdf(query, accepted = TRUE, marine = TRUE)
