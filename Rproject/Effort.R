
# DEREPLICATE_ TO BE REMOVED
# Titration evaluation (sampling effort); 
#  Estimating sampling effort from the produced ground truth paired-end reads. 
#  (A) Phred quality scores of simulated data. 
#  (B) Number of full-length genes (B) and 
#  (C) Average coverage (%) for the obtained number of pairs of reads (millions). 
#  The figure displays the number of full-length conotoxins (y-axis) included at different input read numbers (x-axis) on average. 
#  To assess whether we have simulated sufficient depth to capture all the sequences presented in the groups, we evaluate the sampling effort to 10x, 20x, 50x, 100x, 200x, and 500x coverage.
#  
# Read samtools deph data
# Plot 1: Plot coverage of conotoxins per 
# Evaluates depth
# Dataviz quality read 


rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)
library(patchwork)

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

dir <- "C:/Users/cinai/OneDrive/Escritorio/0_simulated_data_dir/depth_stats_dir/"

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


# df <- lapply(files, read_sam_depth_files)

# df <- do.call(rbind, df)

# df %>% dplyr::count(file)

# df |> write_rds(file = paste0(dir, "depths.rds"))

df <- read_rds(paste0(dir, "depths.rds"))

df <- df %>% 
  separate(file, into = c("dataset", "coverage"), sep = "_", remove = F) |>
  mutate(coverage = factor(coverage, levels = c("10", "20", "50","100","200","500", "700", "1000")))

data_text <- df %>% 
  count(dataset, coverage) %>%
  mutate(coverage = factor(coverage, levels = c("10", "20", "50","100","200","500", "700", "1000"))) %>%
  filter(coverage == 500)

g <- ggplot(df, aes(x = coverage, y = pct_coverage, color = coverage, fill = coverage))

# g + geom_boxplot(alpha = .5, size = 1.5, outlier.size = 5)


df |>
  group_by(dataset, coverage) |>
  sample_n(size = 30) |>
  ggplot(aes(x = coverage, y = pct_coverage)) +
  # ggbeeswarm::geom_quasirandom(
  #   size = 2, width = .33, alpha = .3
  # ) +
  stat_summary(
    fun = mean, geom = "point", 
    shape = 95, size = 12
  ) + 
  ggbeeswarm::geom_quasirandom(
    size = 1, width = .33, shape = 1, color = "black", stroke = .8
  ) +
  my_custom_theme()


library(ggdist)
library(distributional)

g_interval <- df |>
  distinct(dataset, coverage, pct_coverage) |>
  mutate(facet = "B) Total bases covered") |>
  ggplot(aes(x = coverage, y = pct_coverage)) +
  scale_color_viridis_d(
    option = "mako", name = "Level:", direction = -1, 
    begin = .15, end = .9
  ) +
  guides(
    color = guide_legend(reverse = TRUE, title.position = "top")
  ) + my_custom_theme() +
  labs(y = "Covered bases (%)", x = "Number read pairs (millions)") 


p1 <- g_interval +
  # ggdist::stat_interval(
  #   .width = c(.25, .5, .95, 1),
  #   size = 7
  # ) +
  ggdist::stat_gradientinterval(
    width = .25, 
    point_size = 0.5,
    fill = "skyblue"
  ) +
  ggdist::stat_halfeye(
    adjust = .33, ## bandwidth
    width = .5, fill = "grey85",
    interval_colour = NA, point_colour = "black",
    shape = 23, stroke = 1.5, point_size = 0.5, point_fill = "white",
    position = position_nudge(x = .03),
    aes(thickness = stat(f*n))
  ) +
  # scale_color_viridis_d(
  #   option = "mako", name = "Level:", direction = -1, 
  #   begin = .15, end = .9,
  #   labels = function(x) paste0(as.numeric(x)*100, "%")) +
  facet_grid(~ facet)

p1


p3 <- df %>% 
  filter(pct_coverage >= 0.95) %>%
  count(dataset, coverage) %>%
  # mutate(coverage = as.integer(coverage)) %>%
  mutate(coverage = factor(coverage, levels = c("10", "20", "50","100","200","500", "700", "1000"))) %>%
  ggplot(aes(x = coverage, y = n, group = dataset)) + 
  # ggrepel::geom_text_repel(data = data_text, aes(label = dataset), hjust = 2, vjust = -2) +
  # geom_point(color = "gray", size = 2, shape = 1) +
  # geom_line(color = "gray") +
  stat_summary(fun = "mean", geom = "line", aes(group = 1), color="red") +
  stat_summary(fun.data=mean_se, geom="pointrange", aes(group = 1), color="red") +
  # stat_summary(fun = "mean", geom = "point", aes(group = 1), color="red") +
  labs(y = "Number of full-length genes", x = "Number read pairs (millions)") +
  my_custom_theme()

design <- "ABC"


wrap_plots(C = p1, 
                    B = p2, 
                    A = p2, design = design) +
  plot_layout(heights  = c(0.65,0.15,1))

# Step 2: fastqc data
# 
# 

library(fastqcr)

dir <- "//wsl.localhost/Debian/home/ricardo/fastqc_dir"


fileqc <- list.files(dir, pattern = '*zip', full.names = TRUE)

fastq <- data.frame(qc_read_collection(fileqc, modules = "Per base sequence quality", verbose = TRUE)$per_base_sequence_quality)

# qc_plot_collection(fileqc, modules = "Per base sequence quality")

fastq |> ggplot(aes(x = Base, y = Mean)) + geom_boxplot()

fastqStats <- function(fileqc) {
  
  require(fastqcr)
  
  fastq <- data.frame(qc_read_collection(fileqc, modules = "Per base sequence quality", verbose = TRUE)$per_base_sequence_quality)
  
  fastq.tmp <- rbind(data.frame(R=fastq$Base, 
                                Q=fastq$Mean, S=c("Mean"), 
                                E=10^(-fastq$Mean/10), 
                                A=Reduce('+', 10^(-fastq$Mean/10), accumulate = TRUE)),
                     data.frame(R=fastq$Base, 
                                Q=fastq$Median, 
                                S=c("Median"), 
                                E=10^(-fastq$Median/10), 
                                A=Reduce('+', 10^(-fastq$Median/10), accumulate = TRUE)),
                     data.frame(R=fastq$Base, 
                                Q=fastq$Lower.Quartile, 
                                S=c("Lower.Quartile"), 
                                E=10^(-fastq$Lower.Quartile/10), 
                                A=Reduce('+', 10^(-fastq$Lower.Quartile/10), accumulate = TRUE)),
                     data.frame(R=fastq$Base, 
                                Q=fastq$Upper.Quartile, 
                                S=c("Upper.Quartile"), 
                                E=10^(-fastq$Upper.Quartile/10), 
                                A=Reduce('+', 10^(-fastq$Upper.Quartile/10), accumulate = TRUE)),
                     data.frame(R=fastq$Base, 
                                Q=fastq$X10th.Percentile, 
                                S=c("X10th.Percentile"), 
                                E=10^(-fastq$X10th.Percentile/10), 
                                A=Reduce('+', 10^(-fastq$X10th.Percentile/10), accumulate = TRUE)),
                     data.frame(R=fastq$Base, 
                                Q=fastq$X90th.Percentile, 
                                S=c("X90th.Percentile"), 
                                E=10^(-fastq$X90th.Percentile/10), 
                                A=Reduce('+', 10^(-fastq$X90th.Percentile/10), accumulate = TRUE)))
  
  fastq.tmp <- data.frame(fastq.tmp, filename = basename(fileqc))
  
  return(fastq.tmp)
}


fqdata <- do.call(rbind, lapply(files, fastqStats))



fastq.tmp %>% 
  as_tibble() |>
  # mutate(group = sapply(strsplit(as.character(filename), "_"), `[`, 1)) %>%
  # mutate(filename = sapply(strsplit(as.character(filename), "_"), `[`, 1)) %>%
  ggplot(aes(color = S, group = filename)) + 
  geom_point(aes(x = R, y = Q), size=1, alpha = 3/5) + 
  labs(x="Reads position", y="Reads Quality", color = 'Stats') +
  # scale_color_brewer(palette = 'Dark2') +
  theme_bw() 
