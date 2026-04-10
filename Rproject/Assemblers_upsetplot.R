
# January 2026
# Sanity check of y axis are the same between plots
# 

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)
library(patchwork)

# extrafont::loadfonts(device = "win")

my_custom_theme <- function(base_size = 14, legend_pos = "top", ...) {
  # base_size = 14
  theme_bw(base_family = "GillSans", base_size = base_size) +
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

outdir <- "/Users/rjegr/Documents/GitHub/ConotoxinBenchmark/INPUTS/"

dir <- "/Users/rjegr/Documents/Windows/Documents/ConotoxinBenchmark/INPUTS/"

f <- list.files(path = dir, pattern = "curated_nuc_conoServerDB.rds", full.names = T)

# conoServerDB <- read_rds(f) %>% dplyr::rename("hits" = "entry_id")

transratedf <- read_tsv( file.path(outdir, "benchmark_assemblers.tsv"))


dir <- "/Users/rjegr/Documents/Windows/Documents/ConotoxinBenchmark/INPUTS/BLAST_based_annotation_dir/"

f <- list.files(dir, full.names = T, pattern = ".rds")

annotation_results <- do.call(rbind, lapply(f, read_rds)) |>
  mutate(file_name = gsub("_into_Fold[0-9]+[0-9]+.[1|2].blast", "", file_name)) |> 
  dplyr::rename("contig_name" = qseqid, "Assembler" = file_name) |>
  mutate(vfold_set = sapply(strsplit(Assembler, "_"), `[`, 1) ,Assembler = sapply(strsplit(Assembler, "_"), `[`, 5))


# library(ggVennDiagram)
# 

dat <- transratedf |> ungroup() |> 
  # dplyr::filter(!is.na(hits)) %>%
  distinct(contig_name, Assembler, vfold_set, reference_coverage, hits) |>
  # Whether or not distinct to count the true number of assignments
  # distinct(Method, protein_id, Region, tab) |> 
  # Use right Join to count the number of full contigs
  left_join(annotation_results, by = c("contig_name", "Assembler", "vfold_set")) |>
  dplyr::distinct(Assembler, reference_coverage, hits, final_annotation) 


recode_to <- c("STRINGTIE","SPADES", "TRINITY","IDBA", "MEGAHIT", "RNABLOOM", "PLASS","BRIDGER", "TRANSABBYS", "BINPACKER","SOAPDENOVO" ,"CSTONE")

recode_to <- structure(c("StringTie","Spades", "Trinity", "IDBA", "MEGAHIT", "RNA-bloom", "Plass", "BRIDGER","Transabbys", "BinPacker", "SOAP-denovo", "Cstone"), names = recode_to)


dat <- dat |> dplyr::mutate(Assembler = dplyr::recode(Assembler, !!!recode_to)) 

discrete_scale <- dat %>% ungroup() %>% distinct(Assembler) %>% pull()

n <- length(discrete_scale)

scale_col <- c(ggsci::pal_startrek()(7), ggsci::pal_cosmic()(n-7))

scale_col <- structure(scale_col, names = sort(discrete_scale))

f_assembler <- c("StringTie","Spades", "Trinity", "IDBA", "MEGAHIT", "RNA-bloom", "Plass")


dat <- dat |>
  filter(Assembler %in% f_assembler) |>
  dplyr::filter(!is.na(hits)) |>
  filter(final_annotation %in% c("multi", "full"))
  # dplyr::filter(reference_coverage >= 0.95)
 

dat |> count(Assembler) 
# assembler_colors <- c(
#   BinPacker = "#CC0C00FF",
#   BRIDGER = "#5C88DAFF",
#   Cstone = "#84BD00FF",
#   IDBA = "#FFCD00FF",
#   MEGAHIT = "#7C878EFF",
#   PLASS = "#00B5E2FF",
#   `RNA-bloom` = "#00AF66FF",
#   `SOAP-denovo` = "#2E2A2BFF",
#   Spades = "#CF4E9CFF",
#   StringTie = "#8C57A2FF",
#   Transabbys = "#358DB9FF",
#   TransLiG = "#82581FFF",
#   Trinity = "#2F509EFF"
# )

# dat <- transratedf %>% 
#   filter(Assembler %in% assemblers) %>%
#   dplyr::filter(!is.na(hits)) %>%
#   dplyr::distinct(Assembler, reference_coverage, hits) 
# 


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

# hits_combinations_bar <- dat |>
#   distinct(Assembler, hits, final_annotation) |>
#   group_by(hits, final_annotation) |>
#   summarise(
#     combination = paste(sort(unique(Assembler)), collapse = ","),
#     n_assemblers = n_distinct(Assembler),
#     .groups = "drop"
#   ) |>
#   count(final_annotation, combination, name = "count") |>
#   dplyr::filter(!is.na(final_annotation))


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

bar_chart <- 
  hits_combinations |>
  mutate(label = paste0(" (", count,")")) |>
  mutate(col = ifelse( grepl(",", combination), "Intersect", "Unique")) |>
  ggplot(aes(y = combination, x = count)) +
  geom_col(width = 0.6, aes(fill = col)) +
  geom_text(aes(label = label),
            vjust = 0.5, hjust = -0.15, size= 2.5,
            color="black",
            # position=position_dodge(0.5),
            family =  "GillSans") +
  # scale_x_discrete(position = "top", drop = FALSE) + 
  scale_y_discrete(limits = combination_lev, labels = NULL) +
  scale_x_continuous(position = "top", limits = c(0, 750), labels = NULL) +
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

# bar_chart

# 2. Intersection point chart (shows which assemblers in each combination)
# 

point_chart <- points_data |>
  # dplyr::mutate(assemblers = dplyr::recode(assemblers, !!!recode_to)) |>
  mutate(col = ifelse( grepl(",", combination), "Intersect", "Unique")) |>
  mutate(facets = "Intersections") |>
  ggplot(aes(y = combination, x = assemblers)) +
  # ggstats::geom_stripped_cols() +
  geom_tile(aes(fill = assemblers), height = Inf, width = 0.9 , alpha = 0.2) +
  geom_line(aes(group = combination, colour = col), size = 1) +
  geom_point(aes(colour = col), size = 2.5) +
  # facet_grid(facets ~ ., switch = "y") +
  my_custom_theme() +
  scale_x_discrete(position = "top", drop = FALSE) + 
  scale_y_discrete(limits = combination_lev) +
  scale_color_manual(values = c("gray90", "black")) +
  scale_fill_manual(values = scales::alpha(scale_col, 0.2)) +
  guides(fill = "none", color = "none") +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x.top = element_text(angle = -45, vjust = -0.5, hjust = 1)
  ) +
  # coord_cartesian(expand = FALSE) +
  labs(x = "", y = "") 


point_chart

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

# Individual values

mean_combinations_v2 <- dat |>
  group_by(hits, Assembler) |>
  summarise(mean = mean(reference_coverage), n = n(), .groups = "drop")
  # 


# se <- function(sd, n) {sum(sd^2, na.rm = TRUE) / sqrt(n)}

mean_combinations_v2 <- dat |>
  distinct(Assembler, hits) |>
  group_by(hits) |>
  summarise(
    combination = paste(sort(unique(Assembler)), collapse = ","),
    # n_assemblers = n_distinct(Assembler),
    .groups = "drop"
  ) |>
  left_join(mean_combinations_v2) |>
  group_by(combination, Assembler) |>
  # rstatix::get_summary_stats(type = "mean_se") 
  summarise(se = sd(mean)/sum(n), n = sum(n),  mean = mean(mean),  .groups = "drop")

mean_combinations_v2|> filter(combination == "Plass,RNA-bloom,Spades")


mean_chart_v2 <- mean_combinations_v2 |>
  mutate(min = mean-se, max = mean+se) |>
  ggplot(aes(y = combination, x = mean)) +
  geom_point(aes(color = Assembler), size = 1, alpha = 0.8, shape = 3, position = position_dodge(0.8)) +
  geom_errorbar(aes(color = Assembler, xmin = min, xmax = mean), alpha = 0.8, position = position_dodge(0.8), size = 0.5, width = 0) +
  ggstats::geom_stripped_rows() +
  my_custom_theme(base_size = 10) +
  scale_y_discrete(limits = combination_lev, labels = NULL) +
  scale_x_reverse(position = "top", limits = c(0,1.15), breaks = c(0,0.5,1)) +
  scale_color_manual("", values = scale_col) +
  theme(
    axis.text.x = element_text(size = 10),
    legend.position = "none",
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    axis.text.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
  ) +
  # coord_cartesian(expand = FALSE) +
  labs(x = "Reference coverage", y = "")

mean_chart_v2

# Global

mean_combinations <- dat |>
  group_by(hits) |>
  summarise(
    combination = paste(sort(unique(Assembler)), collapse = ","),
    # n_assemblers = n_distinct(Assembler),
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
  # coord_cartesian(expand = FALSE) +
  labs(x = element_blank(), y = element_blank())
         
         
         

# Match plots to areas by name
design <- "#AC
           #B#"


PSAVE <- wrap_plots(C = bar_chart, 
                    B = assembler_bars, 
                    A = point_chart, design = design) +
  plot_layout(heights  = c(0.65,0.15,1))

PSAVE

ggsave(PSAVE, filename = 'figure_3_upset_plot_.png', path = outdir, width = 5, height = 7, device = png, dpi = 800)

# Match plots to areas by name
design <- "DAC"


PSAVE <- wrap_plots(D = mean_chart_v2, 
                    C = bar_chart, 
                    # B = assembler_bars, 
                    A = point_chart, design = design) +
  plot_layout(
    # guides = "collect", axis_titles = "collect",
              heights  = c(1), widths = c(0.25, 0.35, 0.25)) &
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"), panel.spacing = unit(0, "cm"))

# & (Apply to All): Applies a single plot modification (theme, axis scale, etc.) to all plots involved in the layout at once.
  # theme(legend.position = "top",
  #       legend.justification = c(0.5, 0),  # (x, y) coordinates - 0.5 centers horizontally
  #       legend.box = "horizontal",
  #       plot.margin = margin(t = 20)) &
  # guides(fill = "none", color = guide_legend(title = "",
  #                             title.position = "top", title.hjust = 0,
  #                             keywidth = unit(0.35, "cm"), 
  #                             keyheight = unit(0.35, "cm")))
  



# PSAVE

ggsave(PSAVE, filename = 'figure_3_upset_plot.png', path = outdir, width = 5, height = 7, device = png, dpi = 800)

