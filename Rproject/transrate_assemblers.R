
# Read transrate scores (to evals accuracy of assembly methods)

# True Positives (TP): Transcripts correctly assembled by the assembler. 
# False Positives (FP): Transcripts incorrectly assembled by the assembler. 
# False Negatives (FN): Transcripts present in the simulated data but not assembled. 
# Sensitivity: The proportion of true transcripts that are correctly assembled. 
# Precision: The proportion of assembled transcripts that are actually true. 
# F1-score: A harmonic mean of precision and sensitivity, providing an overall measure of assembly quality. 

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)


dir <- "/Users/cigom/Documents/GitHub/ConotoxinBenchmark/2_transrate_dir/"

subdirs <- list.files(dir, pattern = "_transrate_dir")

read_transrate_scores <- function(path) {
  
  basedir <- gsub("_transrate_dir","",path)
  
  Superfamily <- sapply(strsplit(basedir, "_"), `[`, 1)
  Assembler <- gsub(paste0(Superfamily, "_superfamily_"),"",basedir)
  
  f <- list.files(file.path(dir,path, basedir), "contigs.csv", full.names = T)
  
  read_csv(f) %>% mutate(Superfamily = Superfamily, Assembler = Assembler)
  
}

transratedf <- lapply(subdirs, read_transrate_scores)

transratedf <- do.call(rbind,transratedf)

Totaldf <- transratedf %>% 
  count(Superfamily, Assembler)

Hitsdf <- transratedf %>% 
  drop_na(hits) %>%
  count(Superfamily, Assembler)

transratedf %>% 
  drop_na(hits) %>%
  count(Superfamily, Assembler) %>%
  ggplot(aes(Superfamily, Assembler, fill = n)) +
  geom_tile()

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
  facet_grid(~ Superfamily) +
  geom_col()

transratedf %>%
  drop_na(hits) %>%
  ggplot(aes(y = Assembler, x = reference_coverage, fill = after_stat(x))) +
  # geom_violin() +
  # facet_grid(~ facet) +
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
  ggplot(aes(x = n/123, y = Assembler, fill = col)) +
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

DF <- transratedf %>% 
  drop_na(hits) %>%
  filter(reference_coverage > 0.8) %>%
  distinct(Superfamily, Assembler, hits)


gene2ven <- split(strsplit(DF$hits, ";") , DF$Assembler)

gene2ven <- lapply(gene2ven, unlist)

# ggVennDiagram(gene2ven) + scale_fill_gradient(low="grey90",high = "red")

which_venn <- which(names(gene2ven) %in% c("TRINITY", "SPADES","TRANSABBYS","RNABLOOM"))

ggVennDiagram(gene2ven[which_venn]) + scale_fill_gradient(low="grey90",high = "red")
