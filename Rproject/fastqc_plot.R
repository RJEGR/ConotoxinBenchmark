

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


dir <- "//wsl.localhost/Debian/home/ricardo/fastqc_dir"


fileqc <- list.files(dir, pattern = '*zip', full.names = TRUE)


plotdf <- data.frame(qc_read_collection(fileqc, modules = "Per base sequence quality", verbose = TRUE)$per_base_sequence_quality)

which_cols <- names(plotdf)[-c(1,2)]
  
plotdf.summary <- plotdf |> group_by(Base) |> 
  summarise_at(all_of(which_cols), mean) |>
  mutate(color = "Simulated (Art)")


p2 <- plotdf.summary |> 
  mutate(facet = "A) Quality scores") |>
  ggplot(aes(x = Base, y = Mean)) + 
  # geom_tile(aes(fill = Count)) + 
  # scale_fill_gradient(low = "#F5F5F5",high = "black") +
  geom_line(aes(y = Mean, color = color)) +  #  
  geom_line(aes(y = Lower.Quartile), color = "#de2d26", size = 0.7, linetype = "dashed") +
  geom_line(aes(y = Median), color = "#de2d26", size = 0.7) + 
  geom_line(aes(y = Upper.Quartile), color = "#de2d26", size = 0.7, linetype = "dashed") + 
  ylab("Phred Score") + xlab("Cycle") + 
  theme_classic(base_size = 16) + 
  scale_color_manual(values =  c("#2c7fb8")) +
  # theme(panel.grid = element_blank()) + 
  guides(color = guide_legend(title = NULL)) +
  ylim(c(0, NA)) +
  facet_grid(~ facet) +
  my_custom_theme()

###
###
###
##

dir <- "C:/Users/cinai/OneDrive/Escritorio/0_simulated_data_dir/depth_stats_dir/"

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
  # scale_color_viridis_d(
  #   option = "mako", name = "Level:", direction = -1, 
  #   begin = .15, end = .9
  # ) +
  # guides(
  #   color = guide_legend(reverse = TRUE, title.position = "top")
  # ) + 
  my_custom_theme() +
  labs(y = "Covered bases (%)", x = "Number read pairs (millions)")
  


p1 <- g_interval +
  # ggdist::stat_interval(
  #   .width = c(.25, .5, .95, 1),
  #   size = 7
  # ) +
  # ggdist::stat_gradientinterval(
    # width = .5, 
    # point_size = 2.5,
    # fill = "skyblue"
  # ) +
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

# 
# library(patchwork)

PSAVE <- p2+p1

outdir <- "C://Users//cinai/OneDrive/Documentos/GitHub/ConotoxinBenchmark/INPUTS"

ggsave(PSAVE, filename = 'fastqc_plot.png', path = outdir, width = 7, height = 3.5, device = png, dpi = 800)


design <- "ABC"


wrap_plots(C = p1, 
           B = p2, 
           A = p2, design = design) +
  plot_layout(heights  = c(0.65,0.15,1))


