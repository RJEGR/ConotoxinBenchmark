
# Read transrate scores (to evals accuracy of assembly methods)

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


dir <- "/Users/cigom/Documents/GitHub/ConotoxinBenchmark/2_transrate_dir/2_transrate_dir/"

subdirs <- list.files(dir, pattern = "_transrate_dir")

# Note: PLASS records are aminoacid sequences, therefore transrate results must be empty

subdirs <- subdirs[!grepl("PLASS", subdirs)]

read_transrate_scores <- function(path) {
  
  basedir <- gsub("_transrate_dir","",path)
  
  Superfamily <- sapply(strsplit(basedir, "_"), `[`, 1)
  Assembler <- gsub(paste0(Superfamily, "_superfamily[_|.]|all_superfamilies.fixed[_|.]"),"",basedir)
  
  f <- list.files(file.path(dir,path, basedir), "contigs.csv", full.names = T)
  
  read_csv(f) %>% mutate(Superfamily = Superfamily, Assembler = Assembler)
  
}

transratedf <- lapply(subdirs, read_transrate_scores)

transratedf <- do.call(rbind,transratedf)

transratedf %>%  count(Assembler)

# transratedf %>% 
#   select_if(is.double) %>%
#   mutate_all(~replace(., is.na(.), 0)) %>%
#   cor(method = "spearman") -> M
#   
#   
# testRes <- corrplot::cor.mtest(M, conf.level = 0.95)
# 
# corrplot::corrplot(M, p.mat = testRes$p ,method = "color", type="upper", order = "hclust", insig = "label_sig")
# 
# transratedf %>%
#   group_by(Assembler) %>% sample_frac(size = 0.05) %>% ungroup() %>%
#   select_if(is.double) %>%
#   mutate_all(~replace(., is.na(.), 0)) %>%
#   # pairs(pch = 19, lower.panel = NULL)
#   GGally::ggpairs()
  
transratedf %>% drop_na(hits) %>% count(Assembler)

# Calculate metrics  -----

calculate_metrics <- function(transratedf, reference_coverage_val = 1) {
  
  calculate_false <- function(DF) {
    
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
      left_join(DF) %>%
      mutate(FN = abs(InputNsequences - TP))
    
  }
  
 Totaldf <- transratedf %>% 
    count(Superfamily, Assembler)  %>% 
   dplyr::rename("rawcontigs" = "n")
  
 # True Positives (TP): Transcripts correctly assembled by the assembler. 
 ## TP = reference_cov >= reference_coverage_val
 
 TP <- transratedf %>%
    filter(reference_coverage >= reference_coverage_val) %>%
    drop_na(hits) %>%
   # Use distinct to trim redundancy
   distinct(hits, Superfamily, Assembler) %>%
    count(Superfamily, Assembler) %>% 
   dplyr::rename("TP" = "n")
 
 # False Positives (FP): Transcripts incorrectly assembled by the assembler.
 ## FP = N contig_name where is.na(hits) == TRUE, AND|OR reference_cov < 1
 
 FP <- transratedf %>%
   mutate(reference_coverage = ifelse(is.na(hits) & is.na(reference_coverage),
     0,reference_coverage )) %>%
   filter(reference_coverage < reference_coverage_val) %>%
   # drop_na(hits) %>%
   count(Superfamily, Assembler) %>% 
   dplyr::rename("FP" = "n")
 
 # FPnohit <- transratedf %>%
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

metricsdf <- calculate_metrics(transratedf, reference_coverage_val = 0.9) %>% 
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


# metricsdf %>%
#   select_if(is.double) %>%
#   mutate_all(~replace(., is.na(.), 0)) %>%
#   cor(method = "spearman") -> M
# 
# 
# testRes <- corrplot::cor.mtest(M, conf.level = 0.95)
# 
# corrplot::corrplot(M, p.mat = testRes$p ,method = "color", type="upper", order = "hclust", insig = "label_sig")



plot_val <- "F1score" # Precision, Sensitivity

# subtitle <- "Sensitivity: The proportion of true transcripts that are correctly assembled (TP / (TP + FP))"

# subtitle <- "Precision: The proportion of assembled transcripts that are actually true (TP /(TP + FN))"

subtitle <- "F1-score (aka Recall): A harmonic mean of precision and sensitivity,\nproviding an overall measure of assembly quality."

metricsdf  %>%
  dplyr::rename("fill" = plot_val) %>%
  select(all_of(c("fill", "Superfamily", "Assembler"))) %>%
  mutate(fill = ifelse(fill == 0, NA, fill)) %>%
  ggplot(aes(Superfamily, Assembler, fill = fill)) +
  geom_tile(color = 'white', linewidth = 0.2) +
  scale_fill_gradient2(low = "gray", high = "red", mid = "orange", 
    na.value = "white", midpoint = 0.5, limit = c(0, 1), 
    breaks = c(0, 0.5, 1),
    name = NULL) +
  geom_text(aes(label = scales::percent(fill, accuracy = 1)), color = "white", family = "GillSans") +
  labs(subtitle = subtitle) +
  theme_bw(base_family = "GillSans", base_size = 12) +
  theme(legend.position = "top",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.y = element_text(size = 7),
    # axis.text.x = element_text(angle = 0, hjust = 1, vjust = 1,size = 7),
    strip.text = element_text(
      angle = 0, hjust = 0,
      size = 10)) +
  guides(
    fill = guide_colorbar(
      barwidth = unit(2, "in"),
      barheight = unit(0.1, "in"), 
      label.position = "bottom",
      label.hjust = 0.5,
      title = plot_val,
      title.position  = "top", title.hjust = 0,
      title.theme = element_text(size = 10, family = "GillSans", hjust = 1),
      ticks.colour = "black", ticks.linewidth = 0.35,
      frame.colour = "black", frame.linewidth = 0.35,
      label.theme = element_text(size = 10, family = "GillSans")
    ))


which_cols <- metricsdf %>% select_if(is.double) %>% names()

which_cols <- c("InputNsequences", "rawcontigs", "TP","FP", "FN", 
  "Sensitivity","Precision","Recall")  

metricsdf %>%
  pivot_longer(all_of(which_cols), names_to = "facet", values_to = "x") %>%
  mutate(facet = factor(facet, levels = which_cols)) %>%
  ggplot(aes(y = Superfamily, x = x)) +
  facet_grid(~ facet, scales = "free_x") +
  geom_boxplot() +
  theme_bw(base_family = "GillSans", base_size = 12) +
  theme(legend.position = "top",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.y = element_text(size = 7),
    # axis.text.x = element_text(angle = 0, hjust = 1, vjust = 1,size = 7),
    strip.text = element_text(
      angle = 0, hjust = 0,
      size = 10))

# In reference_coverage_val, try a loop to itinerate evaluation of recall




reference_coverage_seq <- seq(0.1, 1, by = 0.1)

metrics_list <- lapply(reference_coverage_seq, function(i) {
  calculate_metrics(transratedf, reference_coverage_val = i) %>%
    mutate(
      Sensitivity = TP / (TP + FN),
      Precision = TP / (TP + FP),
      Recall = 2 * (Precision * Sensitivity) / (Precision + Sensitivity),
      threshold = i
    )
})

metricsdf <- bind_rows(metrics_list)

metricsdf %>% 
  select(threshold, Recall, Assembler, Superfamily) %>%
  ggplot(aes(x = threshold, y = Recall, color = Assembler)) +
  facet_grid(~Assembler) +
  geom_line(aes(group = Superfamily))



# metricsdf %>%
#   ggplot(aes(color = Assembler, Recall)) +
#   ggplot2::stat_ecdf()

metricsdf %>%
  ggplot(aes(color =Assembler , x = Precision, y = Sensitivity)) +
  labs(x = "False positive rate (Precision)", 
    y = "True positive rate (Sensitivity)") +
  facet_wrap(~ Assembler) +
  geom_line()



transratedf %>%
  distinct(Superfamily, hits) %>%
  count(Superfamily)


library(ggVennDiagram)

transratedf %>%
  drop_na(hits) %>%
  count(Assembler, hits, sort = T)

DF <- transratedf %>% 
  drop_na(hits) %>%
  distinct(Superfamily, Assembler, hits)


DF %>% dplyr::count(Superfamily, Assembler, sort = T) %>%
  mutate(Assembler = factor(Assembler, levels = unique(Assembler))) %>%
  ggplot(aes(y = Assembler, x = n)) +
  # facet_grid(~ Superfamily) +
  geom_col()

transratedf %>%
  drop_na(hits) %>%
  ggplot(aes(y = Assembler, x = p_good, fill = after_stat(x))) +
  # geom_violin() +
  # facet_grid(~ Superfamily) +
  ggridges::geom_density_ridges_gradient(
    jittered_points = T,
    position = ggridges::position_points_jitter(width = 0.05, height = 0),
    point_shape = '|', point_size = 3, point_alpha = 1, alpha = 1) +
  scale_fill_viridis_c(option = "C") +
  labs(y = "", x = "Reference coverage (TransRate)") +
  theme_bw(base_family = "GillSans", base_size = 12) + 
  theme(legend.position = "none", 
    # axis.text.y = element_blank(), 
    # axis.ticks.y = element_blank(), 
    axis.title.x = element_text(size = 7)) +
  scale_x_continuous(position = "top")



barpdf <- transratedf %>% 
  # Label those not found in the reference (Ho: chimera?)
  mutate(col = ifelse(is.na(hits), "No-hit", "Hit")) %>%
  # filter(reference_coverage)
  count(Assembler, col) %>%
  mutate(label = paste0("(", n,")"))

p2 <- barpdf %>%
  ggplot(aes(x = n, y = Assembler, fill = col)) +
  geom_col() + 
  theme_bw(base_family = "GillSans", base_size = 12) +
  labs(y = "Method", x = "Accuracy (N contigs/ N hits to reference)") +
  ggthemes::scale_fill_calc() +
  scale_x_continuous(position = "top") +
  theme(legend.position = "bottom", 
    # axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = 7))

p2 <- p2 + 
  geom_text(data = filter(barpdf, col == "Hit"),
    aes(label= label), hjust= 1.05, vjust = 0.5, size = 2.5, family = "GillSans", color = "white")

p2 <- p2 + guides(fill=guide_legend(title = ""))

p2 

DF <- transratedf %>% 
  drop_na(hits) %>%
  filter(reference_coverage > 0.8) %>%
  distinct(Superfamily, Assembler, hits)


gene2ven <- split(strsplit(DF$hits, ";") , DF$Assembler)

gene2ven <- lapply(gene2ven, unlist)

# ggVennDiagram(gene2ven) + scale_fill_gradient(low="grey90",high = "red")

which_venn <- which(names(gene2ven) %in% c("TRINITY", "SPADES","TRANSABBYS","RNABLOOM"))

ggVennDiagram(gene2ven[which_venn]) + scale_fill_gradient(low="grey90",high = "red")



