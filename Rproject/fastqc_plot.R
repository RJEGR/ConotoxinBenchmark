

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)
library(patchwork)
library(fastqcr)

extrafont::loadfonts(device = "win")

my_custom_theme <- function(base_size = 14, legend_pos = "top", ...) {
  base_size = 14
  theme_bw(base_family = "Gillsans", base_size = base_size) +
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


dir <- "/Users/rjegr/Documents/Windows/Debian/fastqc_dir/"


fileqc <- list.files(dir, pattern = '*zip', full.names = TRUE)


plotdf <- data.frame(qc_read_collection(fileqc, modules = "Per base sequence quality", verbose = TRUE)$per_base_sequence_quality)

which_cols <- names(plotdf)[-c(1,2)]
  
plotdf.summary <- plotdf |> group_by(Base) |> 
  summarise_at(all_of(which_cols), mean) |>
  mutate(color = "Simulated (Art)")


dir <- "/Users/rjegr/Documents/Windows/Debian/fastqc_dir/2_post_limpieza/"

fileqc <- list.files(dir, pattern = '*zip', full.names = TRUE)

plotdf <- data.frame(qc_read_collection(fileqc, modules = "Per base sequence quality", verbose = TRUE)$per_base_sequence_quality)

which_cols <- names(plotdf)[-c(1,2)]

plotdf.summary <- plotdf |> group_by(Base) |> 
  summarise_at(all_of(which_cols), mean) |>
  mutate(color = "Real RNAseq") |>
  rbind(plotdf.summary)



p2 <- plotdf.summary |> 
  mutate(facet = "A) Empirical error and quality") |>
  ggplot(aes(x = Base, y = Mean,  color = color)) + 
  geom_line(aes(y = Mean)) +  #  
  geom_line(aes(y = Lower.Quartile), size = 0.5, linetype = "dashed", alpha = 0.5) +
  # geom_line(aes(y = Median), size = 0.7) +
  geom_line(aes(y = Upper.Quartile), size = 0.5, linetype = "dashed", alpha = 0.5) + 
  ylab("Phred score") + xlab("Cycle") + 
  theme_classic(base_size = 16) + 
  scale_color_manual(values =  c("#2c7fb8", "#de2d26")) +
  # theme(panel.grid = element_blank()) + 
  guides(color = guide_legend(title = NULL)) +
  ylim(c(0, NA)) + xlim(0,100) +
  facet_grid(~ facet) +
  my_custom_theme(legend.title = element_text(size = 7), 
                  legend.text = element_text(size = 7))

###
###
###
##

dir <- "C:/Users/cinai/OneDrive/Escritorio/0_simulated_data_dir/depth_stats_dir/"

dir <- "/Users/rjegr/Documents/Windows/"

df <- read_rds(paste0(dir, "depths.rds"))

df <- df %>% 
  separate(file, into = c("dataset", "coverage"), sep = "_", remove = F) |>
  mutate(coverage = factor(coverage, levels = c("10", "20", "50","100","200","500", "700", "1000")))



library(ggdist)
library(distributional)

# ggbeeswarm::geom_quasirandom(
#   size = 1, width = .33, shape = 1, color = "black", stroke = .8
# ) 
# 


g_interval <- df |>
  distinct(dataset, coverage, pct_coverage) |>
  mutate(facet = "B) Total bases covered") |>
  ggplot(aes(x = coverage, y = pct_coverage)) +
  my_custom_theme() +
  labs(y = "% Bases", x = "Number read pairs (millions)")
  


data_text <- df %>% 
  filter(pct_coverage >= 0.95) |>
  group_by(dataset, coverage) |>
  summarise(pct_coverage = min(pct_coverage), n = n()) |>
  group_by(coverage) |>
  summarise(sd = round(sd(n), digits = 0), n = round(mean(n), digits = 0), pct_coverage = min(pct_coverage)) |> 
  mutate(label = paste0("(", n, "±",sd,")"))
  

p1 <- g_interval +
  geom_text(data = data_text, aes(label = label), vjust = 10, hjust = 0.5, size= 1.5, family =  "GillSans") +
  # ggdist::stat_interval(
  #   .width = c(.25, .5, .95, 1),
  #   size = 7
  # ) +
  ggdist::stat_gradientinterval(
  # side = "left", position = "dodge", scale = 0.5,
  fill_type = "gradient",
  linewidth = 0.05, size = Inf,
  # linetype = "dashed",
  fill = "blue"
  ) +
  ggdist::stat_halfeye(scale = 0.5,
    adjust = .33, ## bandwidth
    width = .5, fill = "grey85",
    interval_colour = NA, point_colour = "black",
    shape = 1, stroke = 0.5, point_size = 0.5, point_fill = "white",
    position = position_nudge(x = .03),
    aes(thickness = after_stat(f*n))
  ) +
  # geom_line(aes(y = pct_coverage, group = 2), linewidth = 1) +
  facet_grid(~ facet)

p1

# 
# library(patchwork)

PSAVE <- p2+p1

outdir <- "/Users/rjegr/Documents/GitHub/ConotoxinBenchmark/INPUTS/"

ggsave(PSAVE, filename = 'fastqc_plot.png', path = outdir, width = 7, height = 3.5, device = png, dpi = 800)


design <- "ABC"


wrap_plots(C = p1, 
           B = p2, 
           A = p2, design = design) +
  plot_layout(heights  = c(0.65,0.15,1))


