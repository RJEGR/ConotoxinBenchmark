
# Computes Smith-Waterman alignment scores between all sequences
# Compute for all conotoxins found in conoServerDB
# (optional) split ER signal from every gene sf
# SAVE conoServerDB to fasta file
# RUN mmseqs (see complexity.md notes)
# LOAD _SWscores.tsv file 
# Estimate sequence identity within (min identity) and between (max. id) gene sf sequences. 
# Diag. heatmap of summary results (see figure 2 from c. marmoreus paper)

library(tidyverse)


rm(list = ls())

if(!is.null(dev.list())) dev.off()

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

f <- list.files(path = outdir, pattern = "curated_nuc_conoServerDB.rds", full.names = T)

conoServerDB <- read_rds(f) 

conoServerDB <- conoServerDB %>% distinct(entry_id, genesuperfamily) %>%  drop_na()


dir <- "/Users/cigom/Documents/GitHub/ConotoxinBenchmark/INPUTS/SmithWaterman_dir/"

f <- list.files(dir, pattern = "_SWscores.tsv", full.names = T)

df <- read_tsv(f, col_names = F)

# The file is formatted as a tab-separated list with 12 columns: (1,2) identifiers for query and target sequences/profiles, (3) sequence identity, (4) alignment length, (5) number of mismatches, (6) number of gap openings, (7-8, 9-10) domain start and end-position in query and in target, (11) E-value, and (12) bit score.


col_names <- c("query", "target", "identity", "alignment_length", "mismatch", "gaps", "qstart", "qend", "tstart", "tend", "evalue", "bit")

names(df) <- col_names

df <- df %>% filter(identity != 1)

# Summarise: The percentage of identity was computed using the length of the smallest sequences
# 1) minimun percentage of sequence identitited computed within sequences belonging to the same gene sf (diagonal)
# 2) the maximum percentage identities measured between sequences of gene sf belonging to different gene sf (rest of diagonal)

# Evals, all sf from conoServerDB found in df ---
# 1,727 sequences fom 43 gene sf

conoServerDB %>% distinct(genesuperfamily)

length(unique(sort(c(df$query, df$target)))) # 1762 (because 35 are na.sf)

df %>% 
  select(c("query", "target", "identity")) %>%
  left_join(conoServerDB, by = c("query" = "entry_id")) %>% dplyr::rename("query_sf" = "genesuperfamily") %>%
  left_join(conoServerDB, by = c("target" = "entry_id")) %>% dplyr::rename("target_sf" = "genesuperfamily") %>%
  filter(query_sf =="A superfamily" | target_sf == "A superfamily") %>% distinct(query_sf, target_sf)

# How to join and split sf to NodeQuery and NodeTarget?

# df %>% ggplot(aes(identity)) + geom_histogram()

long_data <- df %>% 
  select(c("query", "target", "identity")) %>%
  left_join(conoServerDB, by = c("query" = "entry_id")) %>% dplyr::rename("query_sf" = "genesuperfamily") %>%
  left_join(conoServerDB, by = c("target" = "entry_id")) %>% dplyr::rename("target_sf" = "genesuperfamily") %>%
  drop_na() %>%
  group_by(query_sf, target_sf) %>%
  rstatix::get_summary_stats(identity, type = "five_number") %>%
  mutate(group = ifelse(query_sf == target_sf, "Within", "Between"))


lower_triangle_data <- long_data %>%
  filter(group == "Between") %>%
  dplyr::rename("fill" = "max") %>%
  select(query_sf, target_sf, fill, group)

diag_data <- long_data %>%
  filter(group == "Within") %>%
  dplyr::rename("fill" = "min") %>%
  select(query_sf, target_sf, fill, group)


rbind(lower_triangle_data, diag_data) %>%
  pivot_longer(cols = c(query_sf, target_sf), values_to = "gene_sf") %>%
  distinct(fill, group, gene_sf) %>%
  ggplot(aes(x = fill, y = gene_sf)) +
  facet_grid(~ group) +
  geom_point()


matrix <- rbind(lower_triangle_data, diag_data) %>% 
  arrange(query_sf, target_sf) %>%
  mutate_if(is.character, as.factor) %>%
  mutate(query_sf = gsub(" ", ".", query_sf), target_sf = gsub(" ", ".", target_sf)) %>%
  # mutate(target_sf = factor(target_sf, level = unique(query_sf)))
  pivot_wider(names_from = target_sf, values_from = fill, values_fill = 0) %>%
  data.frame(row.names = "query_sf") %>% as("matrix")

matrix[lower.tri(matrix, diag = FALSE)] <- NA # Set diag = TRUE to include diagonal

df_long_upper <- reshape2::melt(matrix, na.rm = TRUE)

ggplot(df_long_upper, aes(Var1, Var2, fill = value)) +
  # geom_tile() +
  geom_text(aes(label = value), family = "GillSans")
  
rbind(lower_triangle_data, diag_data) %>%
  mutate_if(is.character, as.factor) %>%
  # rowwise() %>%
  # mutate(pair = sort(c(query_sf, target_sf)) %>% paste(collapse = ",")) %>%
  # group_by(pair) %>%
  # distinct(pair, .keep_all = T) %>%
  ggplot(aes(x = query_sf, y = target_sf)) +
  geom_tile(color = "white", fill = "gray") + # Add white borders for better visualization
  geom_text(aes(label = fill), family = "GillSans") +
  my_custom_theme(axis.text.x = element_text(angle=90, hjust=TRUE)) +
  coord_fixed() # Ensures square cells

