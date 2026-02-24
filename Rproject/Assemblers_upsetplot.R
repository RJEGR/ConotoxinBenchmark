
# January 2026
# Sanity check of y axis are the same between plots
# 

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

outdir <- "C://Users//cinai/OneDrive/Documentos/GitHub/ConotoxinBenchmark/INPUTS"


f <- list.files(path = outdir, pattern = "curated_nuc_conoServerDB.rds", full.names = T)

# conoServerDB <- read_rds(f) %>% dplyr::rename("hits" = "entry_id")


dir <- "C://Users//cinai/OneDrive/Escritorio/transrate_contigs_reference_dir/"

transratedf <- read_tsv( file.path(dir, "benchmark_assemblers.tsv"))


# library(ggVennDiagram)
# 
assemblers <- c("STRINGTIE","SPADES", "TRINITY", "RNABLOOM", "IDB", "MEGAHIT", "PLASS")

dat <- transratedf %>% 
  filter(Assembler %in% assemblers) %>%
  dplyr::filter(!is.na(hits)) %>%
  # dplyr::filter(reference_coverage > 0.9) %>%
  dplyr::distinct(Assembler, reference_coverage, hits) 


# dat <- transratedf %>% 
  # filter(Assembler %in% assemblers) %>%
  # filter(!is.na(hits)) %>%
  # filter(reference_coverage > 0.95) %>%
  # mutate(summarise = "< 80 % alignment") %>%
  # mutate(summarise = ifelse(reference_coverage >= 0.8, ">= 80% alignment", summarise)) %>%
  # mutate(summarise = ifelse(reference_coverage >= 0.9, ">= 90% alignment", summarise)) %>%
  # mutate(summarise = ifelse(reference_coverage >= 0.95, ">= 95% alignment", summarise)) %>%
  # mutate(summarise = ifelse(reference_coverage == 1, "100% alignment", summarise)) %>%
  # dplyr::distinct(Assembler, hits)


# Create combinations of assemblers for each hit
# 
hits_combinations <- dat |>
  distinct(Assembler, hits) |>
  group_by(hits) |>
  summarise(
    combination = paste(sort(unique(Assembler)), collapse = ","),
    n_assemblers = n_distinct(Assembler),
    .groups = "drop"
  ) |>
  count(combination, name = "count")

# Reorder by count (descending)
hits_combinations <- hits_combinations |>
  mutate(combination = fct_reorder(combination, count, .desc = TRUE))

# Get assembler counts
assembler_counts <- dat |>
  distinct(Assembler, hits) |>
  count(Assembler, name = "hits_count") |>
  arrange(hits_count) |>
  mutate(Assembler = fct_reorder(Assembler, hits_count, .desc = TRUE))
  

# Create intersection points data (assembler x combination)
points_data <- hits_combinations |>
  mutate(
    assemblers = str_split(combination, ",")
  ) |>
  unnest(assemblers) |>
  mutate(
    assemblers = factor(assemblers, levels = levels(assembler_counts$Assembler)),
    combination = factor(combination, levels = levels(hits_combinations$combination))
  ) |>
  filter(assemblers != "")

combination_lev <- rev(levels(hits_combinations$combination))

# 1. Main bar chart (hits per combination)
bar_chart <- hits_combinations |>
  mutate(label = paste0(" (", count,")")) |>
  mutate(col = ifelse( grepl(",", combination), "Intersect", "Unique")) |>
  ggplot(aes(y = combination, x = count)) +
  geom_col(width = 0.6, aes(fill = col)) +
  geom_text(aes(label = label),
            vjust = 0.5, hjust = -0.15, size= 2.5,
            color="black",
            # position=position_dodge(0.5),
            family =  "Gill Sans MT") +
  # scale_x_discrete(position = "top", drop = FALSE) + 
  scale_y_discrete(limits = combination_lev, labels = NULL) +
  scale_x_continuous(position = "top", limits = c(0, 750), breaks = c(0,300,600), labels = NULL) +
  my_custom_theme() +
  scale_fill_manual(values = c("gray90", "black")) +
  theme(
    legend.position = "none",
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x.top = element_text(angle = 0, vjust = -0.5, hjust = 1)
  ) +
  coord_cartesian(expand = FALSE) +
  labs(x = element_blank(), y = element_blank())

# 2. Intersection point chart (shows which assemblers in each combination)
# 

point_chart <- points_data |>
  mutate(col = ifelse( grepl(",", combination), "Intersect", "Unique")) |>
  mutate(facets = "Intersections") |>
  ggplot(aes(y = combination, x = assemblers)) +
  geom_line(aes(group = combination, colour = col), size = 1) +
  geom_point(aes(colour = col), size = 3.5) +
  # facet_grid(facets ~ ., switch = "y") +
  my_custom_theme() +
  ggstats::geom_stripped_cols() +
  scale_x_discrete(position = "top", drop = FALSE) + 
  scale_y_discrete(limits = combination_lev) +
  scale_color_manual(values = c("gray90", "black")) +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x.top = element_text(angle = -45, vjust = -0.5, hjust = 1)
  ) +
  labs(x = element_blank(), y = element_blank()) 




# 3. Side bar chart (total hits per assembler)

assembler_bars <- assembler_counts |>
  ggplot(aes(y = hits_count, x = Assembler)) +
  geom_col(width = 0.6, fill = 'gray90') +
  # geom_segment(aes(x = Assembler, y = 0, yend = -hits_count), size = 12) +
  geom_text(aes(label = Assembler), size = 1.5, hjust = 1.5, vjust = 0.5, angle = -90, color = "white") +
  scale_y_reverse() +
  # ggstats::geom_stripped_cols() +
  my_custom_theme() +
  theme(
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
  ) +
  coord_cartesian(expand = FALSE) +
  labs(x = element_blank(), y = "# Conotoxins")

# 4.
# 

# Create combinations of assemblers 
# 
mean_combinations <- dat |>
  group_by(hits) |>
  summarise(
    combination = paste(sort(unique(Assembler)), collapse = ","),
    n_assemblers = n_distinct(Assembler),
    reference_coverage = mean(reference_coverage),
    .groups = "drop"
  ) |>
  group_by(combination) |> rstatix::get_summary_stats(type = "mean_sd") 

mean_chart <- mean_combinations |>
  # mutate(assemblers = str_split(combination, ",")) |> unnest(assemblers) |>
  mutate(min = mean-sd, max = mean+sd) |>
  ggplot(aes(y = combination, x = mean)) +
  geom_point(size = 1.5, alpha = 0.5) +
  geom_errorbar(aes(xmin = min, xmax = max), width = 0.2) +
  ggstats::geom_stripped_rows() +
  my_custom_theme() +
  scale_y_discrete(limits = combination_lev, labels = NULL) +
  scale_x_continuous(position = "top", limits = c(0,1.15), breaks = c(0,0.5,1)) +
  theme(
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    axis.text.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
  ) +
  coord_cartesian(expand = FALSE) +
  labs(x = element_blank(), y = element_blank())
         
         
         

# Match plots to areas by name
design <- "#AC
           #B#"


PSAVE <- wrap_plots(C = bar_chart, 
                    B = assembler_bars, 
                    A = point_chart, design = design) +
  plot_layout(heights  = c(0.65,0.15,1))

PSAVE 

ggsave(PSAVE, filename = 'UPSET_FOR_PUB.png', path = outdir, width = 5, height = 7, device = png, dpi = 800)

# Match plots to areas by name
design <- "#ACD
           #B##"


PSAVE <- wrap_plots(D = mean_chart, 
                    C = bar_chart, 
                    B = assembler_bars, 
                    A = point_chart, design = design) +
  plot_layout(heights  = c(0.65,0.15,1,0.75))

ggsave(PSAVE, filename = 'UPSET_FOR_PUB_.png', path = outdir, width = 7, height = 7, device = png, dpi = 800)

