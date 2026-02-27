#
# After running Subsampling.sh analysis, OPEN contigs.csv files and calculate metrics of accuracy
# EDA of 2_subsampling_dir/2_transrate_contigs_dir


# scp -r rgomez@omica.cicese.mx:/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/fernando_pub/2_subsampling_dir/2_transrate_contigs_dir .


rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

extrafont::loadfonts(device = "win")


my_custom_theme <- function(base_size = 14, legend_pos = "top", ...) {
  base_size = 14
  theme_bw(base_family = "Gill Sans MT", base_size = base_size) +
    theme(legend.position = legend_pos,
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
# Reading conoServer info


outdir <- "C://Users//cinai/OneDrive/Documentos/GitHub/ConotoxinBenchmark/INPUTS"


# outdir <- "~/Documents/GitHub/ConotoxinBenchmark/INPUTS/"

f <- list.files(path = outdir, pattern = "curated_nuc_conoServerDB.rds", full.names = T)

conoServerDB <- read_rds(f) %>% dplyr::rename("hits" = "entry_id")

# dir <- "/Users/cigom/Documents/GitHub/ConotoxinBenchmark/2_subsampling_dir/Trinity_dir/transrate_contigs_dir/"

dir <- "C://Users//cinai/OneDrive/Documentos/GitHub/ConotoxinBenchmark/INPUTS/2_subsampling_dir/Trinity_dir/transrate_contigs_dir/"


read_files <- function(dir, Assembly = ...) {
  
  file_list <- list.files(path = dir, pattern = "contigs.csv", recursive = T, full.names = TRUE)
  
  transratedf <- lapply(file_list, read_transrate_scores)
  
  transratedf <- do.call(rbind,transratedf)
  
  transratedf %>% 
    mutate(rel_path = basename(dirname(rel_path))) %>%
    mutate(vfold_set = sapply(strsplit(rel_path, "_"), `[`, 1)) %>%
    mutate(sampling_set = sapply(strsplit(rel_path, "_"), `[`, 5)) %>%
    mutate(sampling_set = as.double(sampling_set)) %>%
    mutate(Assembly = Assembly)
}

transratedf <- read_files(dir, Assembly = "Trinity")


# dir <- "/Users/cigom/Documents/GitHub/ConotoxinBenchmark/2_subsampling_dir/Spades_dir/transrate_contigs_dir/"
dir <- "C://Users//cinai/OneDrive/Documentos/GitHub/ConotoxinBenchmark/INPUTS/2_subsampling_dir/Spades_dir/transrate_contigs_dir/"



transratedf2 <- read_files(dir, Assembly = "Spades")

transratedf <- read_files(dir, Assembly = "Spades") %>% rbind(transratedf)

transratedf %>%
  dplyr::count(rel_path, vfold_set, sampling_set, Assembly) 

# calculate_metrics(transratedf, reference_coverage_val = 0.95) %>%
#   write_tsv(file = file.path(dir, "subsampling_benchmark.tsv"))

metricsdf <- transratedf %>% 
  filter(length>=200) %>%
  group_by(vfold_set, sampling_set, Assembly) %>%
  calculate_metrics(reference_coverage_val = 0.95) %>% 
  mutate(
    Ratio = TP/FP,
    # Tell us what percentage of positive classes were correctly identified
    Accuracy = TP / (TP + FN + FP),
    Precision = TP /(TP + FP),
    Sensitivity = TP /(TP + FN),
    Fscore = 2 * (TP) / (2 * (TP) + FP + FN), 
  ) 

metricsdf %>%
  write_tsv(file.path(outdir, "Sumbsampling_accuracy.tsv"))


metricsdf <- read_tsv(file.path(outdir, "Sumbsampling_accuracy.tsv"))

# (Quantitative): Proxy 1

DataViz <- transratedf %>% 
  drop_na() %>% 
  distinct(sampling_set, vfold_set, hits, reference_coverage, Assembly) %>%
  # left_join(conoServerDB,by = "hits") %>%
  mutate(summarise = "< 80 % alignment") %>%
  mutate(summarise = ifelse(reference_coverage >= 0.8, ">= 80% alignment", summarise)) %>%
  mutate(summarise = ifelse(reference_coverage >= 0.9, ">= 90% alignment", summarise)) %>%
  mutate(summarise = ifelse(reference_coverage >= 0.95, ">= 95% alignment", summarise)) %>%
  mutate(summarise = ifelse(reference_coverage == 1, "100% alignment", summarise)) %>%
  dplyr::count(sampling_set, vfold_set, Assembly, summarise)


n_pallet <- length(unique(DataViz$summarise))

scale_col <- ggsci::pal_uchicago()(n_pallet) 

scale_col <- structure(scale_col, names = sort(unique(DataViz$summarise)))

scale_fill <- ggsci::pal_uchicago(alpha = 0.5)(n_pallet) 

scale_fill <- structure(scale_fill, names = sort(unique(DataViz$summarise)))


p2 <- DataViz %>%
  # filter(summarise != "< 80 % alignment") %>%
  mutate(Assembly = ifelse(Assembly %in% "Spades","A) Spades", "B) Trinity")) %>%
  ggplot(aes(y = n, x = as.factor(sampling_set), color = summarise, fill = summarise)) +
  # geom_jitter(position = position_jitter(0.1), shape = 1) +
  stat_summary(fun = "mean", geom = "line", aes(group = summarise)) +
  # stat_summary(fun = "mean", geom = "point", aes(group = summarise)) + 
  stat_summary(fun.data=mean_se, geom="pointrange", shape = 1, size = 0.3) + # position = position_jitter(0.25)
  labs(y = "Number of assembled conotoxins", x = "Sample size (Proportion of the sample)", caption = "3_Subsampling.R") +
  my_custom_theme(legend.text = element_text(size = 5)) +
  scale_color_manual("",values = scale_col ) +
  scale_fill_manual("",values = scale_fill)
  
p2 <- p2 + facet_grid(~ Assembly)

ggsave(p2,
    filename = 'Subsampling_boxplot_spades_trinity.png',
    path = outdir, width = 5.2, height = 4, dpi = 1000, device = png)

  
# (Quantitative): Proxy 2. What is the distribution in genesuperfamily ? Why some sf are easy to assembly?
# Based on entropy avarege per family, diversity of sequences, and cluster similarity (see )
# 
conoServerDB %>%
  dplyr::count(genesuperfamily)

transratedf %>% 
  # filter(sampling_set == 1) %>%
  # distinct(sampling_set, vfold_set, hits, reference_coverage) %>%
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

cols_to <- c("Precision", "Sensitivity", "Accuracy")

recode_to <- structure(c("A) Precision", "B) Sensitivity", "C) Accuracy"), names = cols_to)

base_size <- 14

p1 <- metricsdf %>%
  select(-TP, -FP, -FN, -rawcontigs) %>%
  pivot_longer(cols = cols_to, values_to = "y", names_to = "facet") %>%
  dplyr::mutate(facet = dplyr::recode_factor(facet, !!!recode_to)) %>%
  drop_na() %>%
  mutate(Assembly = ifelse(Assembly %in% "Spades","A) Spades", "B) Trinity")) %>%
  ggplot(aes(y = y, x = as.factor(sampling_set))) +
  facet_grid(facet ~Assembly, scales ="free", switch = "y") +
  # geom_jitter(position = position_jitter(0.1), shape = 1) +
  stat_summary(fun = "mean", geom = "line", aes(group = 1), color="grey20") +
  # stat_summary(fun = "mean", geom = "point", aes(group = 1), color="gray70") + # , color="blue"
  stat_summary(fun.data=mean_se, geom="pointrange", color="grey20", shape = 1) +
  # ylim(0,1) +
  labs(x = "Sample size (Proportion of the sample)", y = "", caption = "Subsampling.R") +
  # ylim(0,NA) +
  my_custom_theme()

p1

# 
# metricsdf %>% 
#   select(-TP, -FP, -FN, -rawcontigs) %>%
#   pivot_longer(cols = cols_to, values_to = "y", names_to = "facet") %>%
#   mutate(Assembler = Assembly) %>%
#   group_by(Assembler, sampling_set, facet) %>% 
#   summarise(sd = sd(y), x = mean(y)) %>%
#   arrange(desc(x)) %>%
#   mutate(label = paste0(Assembler, " (",  round(x, 3), ")")) %>%
#   mutate(Assembler = factor(Assembler, levels = rev(unique(Assembler)))) %>%
#   # mutate(facet = "B) Accuracy") %>%
#   mutate(xmin = x-sd, xmax = x+sd, xlab = xmax+0.1) %>%
#   mutate(xlab = ifelse(Assembler == "STRINGTIE", 0.25, xlab)) %>% 
#   mutate(sampling_set = as.factor(sampling_set)) %>%
#   ggplot(aes(y = sampling_set , x = x, fill= Assembler, color = Assembler)) + 
#   facet_grid(~ facet, scales = "free") +
#   geom_col(position = position_dodge2()) 
#   # geom_text(aes(y = sampling_set      , x = xlab, label = label), hjust = 0.1, size = 3) +
#   geom_errorbar(aes(xmin = xmin, xmax = xmax), width = 0.15, alpha = 0.3, position = position_dodge2()) 
#   # scale_x_continuous("",limits = c(0,1)) +
#   # scale_y_discrete(position = "right")+
#   scale_color_manual("", values = scale_col) +
#   scale_fill_manual("", values = scale_col) +
#   my_custom_theme(legend_pos = "none", axis.ticks.y = element_blank(), 
#                   axis.text.y = element_blank(), axis.title.y = element_blank())  


ggsave(p1, filename = 'Subsampling.png', 
       path = outdir, width = 5.5, height = 5, dpi = 1000, device = png)



# metricsdf %>%
#   drop_na() %>%
#   mutate(Assembly = ifelse(Assembly %in% "Spades","A) Spades", "B) Trinity")) %>%
#   ggplot(aes(y = Precision, x = Sensitivity , color = as.factor(sampling_set))) +
#   facet_grid( ~Assembly, scales ="free", switch = "y") +
#   # geom_jitter(position = position_jitter(0.1), shape = 1) +
#   stat_summary(fun = "mean", geom = "point")  # , color="blue"
#   # stat_summary(fun.data=mean_sdl, geom="pointrange", color="red") +
#   labs(x = "Sample size (Proportion of the sample)", y = "", caption = "Subsampling.R") +
#   # ylim(0,NA) +
#   my_custom_theme()


# Calculate stats
# 
# Priori
# 

data <- metricsdf %>%
  select(-TP, -FP, -FN, -rawcontigs, -Ratio) %>%
  pivot_longer(cols = cols_to, values_to = "y", names_to = "facet") %>%
  # dplyr::mutate(facet = dplyr::recode_factor(facet, !!!recode_to)) %>%
  drop_na()

library(rstatix)

data <- data %>%
  group_by(Assembly , facet, sampling_set) %>%
  rstatix::get_summary_stats() %>% 
  mutate (y = mean )

data %>%
  group_by(Assembly , facet, sampling_set) %>% rstatix::shapiro_test(y) %>%
  mutate(gauss = ifelse(p > 0.05, TRUE, FALSE)) %>%
  count(gauss)

data %>%
  group_by(Assembly , facet) %>% levene_test(y ~ as.factor(sampling_set)) %>%
  mutate(hom_var = ifelse(p > 0.05, TRUE, FALSE))


# ANOVA assumes that the variance of the residuals is equal for all groups.
# By groups Precision and assembly, Evaluate if significatn differences in performance due to samples depth
#  
data %>%
  group_by(facet, Assembly) %>%
  # rstatix::anova_test(y ~ as.factor(sampling_set)) %>%
  kruskal_test(y ~ as.factor(sampling_set)) %>%
  adjust_pvalue(method = "none") %>%
  add_significance("p") -> res_aov

View(res_aov)

# Posteriori
# 
#   # rstatix::pairwise_t_test(Count ~ pH, paired = F) 
data %>%
  group_by(facet, Assembly ) %>%
  # tukey_hsd(y ~ as.factor(sampling_set)) %>% # If ANOVA used as priori
  rstatix::wilcox_test(y ~ sampling_set)
  adjust_pvalue() %>%
  add_significance()-> stat_test


View(stat_test)

stat_test %>%
  count(p.adj.signif)


# 
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
