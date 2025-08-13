

# Read samtools deph data
# Evaluates depth
#


rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)


dir <- "~/Documents/GitHub/ConotoxinBenchmark/0_simulated_data_dir/depth_stats_dir/"

files <- list.files(dir, pattern = "txt",full.names = T)

read_sam_depth_files <- function(f) {

  cat("\nSummarising ", gsub("._PE.depth.txt", "", basename(f)), "\n")  

  read_table(f, col_names = c("chrom", "pos", "depth")) %>%
    group_by(chrom) %>%
    summarise(
      total_depth = sum(depth),
      mean_depth = mean(depth),
      # median_depth = median(depth),
      sd_depth = sd(depth),
      # min_depth = min(depth),
      # max_depth = max(depth),o
      total_bases = n(),
      zero_coverage = sum(depth == 0),
      covered_bases = total_bases-zero_coverage,
      pct_coverage = covered_bases/total_bases
    ) %>%
    mutate(file = gsub("._PE.depth.txt", "", basename(f)))
  
}


df <- lapply(files, read_sam_depth_files)

df <- do.call(rbind, df)

df %>% count(file)

df <- df %>% 
  separate(file, into = c("dataset", "coverage"), sep = "_", remove = F) 

data_text <- df %>% 
  count(dataset, coverage) %>%
  mutate(coverage = factor(coverage, levels = c("10", "20", "50","100","200","500"))) %>%
  filter(coverage == 500)


df %>% 
  # group_by(file) %>% mutate(total_depth = total_depth/sum(total_depth)) %>%
  mutate(coverage = factor(coverage, levels = c("10", "20", "50","100","200","500"))) %>%
  ggplot(aes(y = pct_coverage, x = coverage)) +
  # geom_jitter(position = position_jitter(0.1), shape = 1) +
  stat_summary(fun = "mean", geom = "line", aes(group = 1), color="red") +
  # stat_summary(fun.data=mean_sdl, geom="pointrange", color="red") +
  stat_summary(fun = "mean", geom = "point", aes(group = 1), color="red") +
  labs(y = "Average Coverage (%)", x = "Number of pairs of reads (millions)") +
  theme_bw(base_family = "GillSans", base_size = 15) +
  theme(legend.position = "top",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank())

df %>% 
  # filter(zero_coverage !=0) %>%
  filter(coverage >= 0.95) %>%
  count(dataset, coverage) %>%
  # mutate(coverage = as.integer(coverage)) %>%
  mutate(coverage = factor(coverage, levels = c("10", "20", "50","100","200","500"))) %>%
  # filter(dataset == "Fold05")  %>%
  ggplot(aes(x = coverage, y = n, group = dataset)) + 
  # ggrepel::geom_text_repel(data = data_text, aes(label = dataset), hjust = 2, vjust = -2) +
  geom_point(color = "gray", size = 2, shape = 1) +
  geom_line(color = "gray") +
  stat_summary(fun = "mean", geom = "line", aes(group = 1), color="red") +
  stat_summary(fun = "mean", geom = "point", aes(group = 1), color="red") +
  labs(y = "number of full-length genes", x = "Number of pairs of reads (millions)") +
  theme_bw(base_family = "GillSans", base_size = 15) +
  theme(legend.position = "top",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank())


