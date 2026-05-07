# January 2026
# Upset plot of assembler combinations x gene superfamily
# Sanity check: y-axis alignment across panels

rm(list = ls())
if (!is.null(dev.list())) dev.off()
options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)


# ─────────────────────────────────────────────────────────────────────────────
# CONFIGURE PATHS HERE --------------------------------------------------------
tmpdir  <- "/Users/rjegr/Documents/Windows/Documents/ConotoxinBenchmark/INPUTS/"
outdir  <- "/Users/rjegr/Documents/GitHub/ConotoxinBenchmark/INPUTS/"
blastdir <- "/Users/rjegr/Documents/Windows/Documents/ConotoxinBenchmark/INPUTS/BLAST_based_annotation_dir/"
# ─────────────────────────────────────────────────────────────────────────────


library(tidyverse)
library(patchwork)

# ── Theme ────────────────────────────────────────────────────────────────────
my_custom_theme <- function(base_size = 14, legend_pos = "top", ...) {
  # base_family <- if ("GillSans" %in% systemfonts::system_fonts()$family) "GillSans" else ""
  base_family <- "GillSans"
  theme_bw(base_family = base_family, base_size = base_size) +
    theme(
      legend.position   = legend_pos,
      strip.placement   = "outside",
      strip.background  = element_rect(fill = "gray90", color = "white"),
      strip.text        = element_text(angle = 0, size = base_size, hjust = 0),
      axis.text         = element_text(size = rel(0.7), color = "black"),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.x = element_blank(),
      ...
    )
}

# ── Load data ─────────────────────────────────────────────────────────────────
f_cono <- list.files(tmpdir, pattern = "curated_nuc_conoServerDB.rds", full.names = TRUE)

conoServerDB <- read_rds(f_cono) |>
  mutate(
    genesuperfamily = ifelse(is.na(genesuperfamily), "Other", genesuperfamily),
    genesuperfamily = gsub(" superfamily", "", genesuperfamily)
  ) |>
  dplyr::rename(gs_conoServer = genesuperfamily, hits = entry_id) |>
  distinct(hits, gs_conoServer)

transratedf <- read_tsv(file.path(outdir, "benchmark_assemblers.tsv"))

f_blast <- list.files(blastdir, full.names = TRUE, pattern = "\\.rds$")

annotation_results <- do.call(rbind, lapply(f_blast, read_rds)) |>
  mutate(
    file_name = gsub("_into_Fold[0-9]+[0-9]+.[12].blast", "", file_name)
  ) |>
  dplyr::rename(contig_name = qseqid, Assembler = file_name) |>
  mutate(
    vfold_set = sapply(strsplit(Assembler, "_"), `[`, 1),
    Assembler = sapply(strsplit(Assembler, "_"), `[`, 5)
  )

# ── Assembler name recoding ───────────────────────────────────────────────────
recode_to <- c(
  STRINGTIE = "StringTie",  SPADES    = "rnaSPAdes",
  TRINITY   = "Trinity",    IDBA      = "IDBA",
  MEGAHIT   = "MEGAHIT",    RNABLOOM  = "RNA-Bloom",
  BRIDGER   = "BRIDGER",    TRANSABBYS= "Trans-ABySS",
  BINPACKER = "BinPacker",  SOAPDENOVO= "SOAP-denovo",
  CSTONE    = "Cstone",     TRANSLIG  = "TransLiG",
  Baseline  = "Baseline",   PLASS     = "PinguiN (nuclassemble)"
)

scale_col <- c(
  BinPacker   = "#CC0C00FF", BRIDGER     = "#5C88DAFF",
  Cstone      = "#84BD00FF", IDBA        = "#FFCD00FF",
  MEGAHIT     = "#7C878EFF", PLASS       = "#00B5E2FF",
  `RNA-Bloom` = "#00AF66FF", `SOAP-denovo`= "#2E2A2BFF",
  rnaSPAdes   = "#CF4E9CFF", StringTie   = "#8C57A2FF",
  Transabbys  = "#358DB9FF", TransLiG    = "#82581FFF",
  Trinity     = "#2F509EFF"
)

# ── Build analysis dataset ────────────────────────────────────────────────────
f_assembler <- c("rnaSPAdes", "Trinity", "RNA-Bloom")

dat <- transratedf |>
  ungroup() |>
  distinct(contig_name, Assembler, vfold_set, reference_coverage, hits) |>
  left_join(annotation_results, by = c("contig_name", "Assembler", "vfold_set")) |>
  distinct(Assembler, reference_coverage, hits, final_annotation) |>
  left_join(conoServerDB, relationship = "many-to-many") |>
  mutate(Assembler = dplyr::recode(Assembler, !!!recode_to)) |>
  filter(
    Assembler %in% f_assembler,
    !is.na(hits),
    final_annotation %in% c("multi", "full")
  )

# ── Combination summaries ─────────────────────────────────────────────────────
hits_combinations <- dat |>
  distinct(Assembler, hits) |>
  group_by(hits) |>
  summarise(
    combination  = paste(sort(unique(Assembler)), collapse = ","),
    n_assemblers = n_distinct(Assembler),
    .groups      = "drop"
  ) |>
  count(combination, name = "count") |>
  mutate(combination = fct_reorder(combination, count, .desc = TRUE))

assembler_counts <- dat |>
  distinct(Assembler, hits) |>
  count(Assembler, name = "hits_count") |>
  arrange(hits_count) |>
  mutate(Assembler = fct_reorder(Assembler, hits_count, .desc = TRUE))

combination_lev <- rev(levels(hits_combinations$combination))

# Intersection matrix (which assemblers appear in each combination)
points_data <- hits_combinations |>
  mutate(assemblers = str_split(combination, ",")) |>
  unnest(assemblers) |>
  mutate(
    assemblers  = factor(assemblers, levels = levels(assembler_counts$Assembler)),
    combination = factor(combination, levels = levels(hits_combinations$combination))
  ) |>
  filter(assemblers != "")

# ── Plot 1: bar chart – hits per combination ──────────────────────────────────
bar_chart <- hits_combinations |>
  mutate(
    label = paste0(" (", count, ")"),
    col   = ifelse(grepl(",", combination), "Intersect", "Unique")
  ) |>
  ggplot(aes(y = combination, x = count)) +
  geom_col(width = 0.6, aes(fill = col)) +
  geom_text(aes(label = label), vjust = 0.5, hjust = 1.5,
            size = 1.5, color = "black") +
  scale_y_discrete(limits = combination_lev, labels = NULL) +
  scale_x_reverse(position = "top", labels = NULL) +
  scale_fill_manual(values = c("gray90", "black")) +
  my_custom_theme() +
  theme(
    legend.position      = "none",
    axis.ticks           = element_blank(),
    panel.border         = element_blank(),
    axis.text.x.top      = element_text(angle = 0, vjust = -0.5, hjust = 1)
  ) +
  coord_cartesian(expand = FALSE) +
  labs(x = NULL, y = NULL)

# ── Plot 2: intersection point chart ─────────────────────────────────────────

breaks <- hits_combinations |> mutate(label = paste0("(", count, ")")) |> 
  # mutate(label = reorder(label, count, .desc = TRUE)) |>
  pull(combination, name = label)


point_chart <- points_data |>
  mutate(col = ifelse(grepl(",", combination), "Intersect", "Unique")) |>
  ggplot(aes(y = combination, x = assemblers)) +
  geom_tile(aes(fill = assemblers), height = Inf, width = 0.9, alpha = 1) +
  geom_line(aes(group = combination, colour = col), linewidth = 0.7) +
  geom_point(aes(colour = col), size = 1.25) +
  scale_x_discrete(position = "top", drop = FALSE) +
  scale_y_discrete(limits = combination_lev, breaks = breaks) +
  scale_color_manual(values = c("gray90", "black")) +
  scale_fill_manual(values = scales::alpha(scale_col, 1)) +
  guides(fill = "none", color = "none") +
  my_custom_theme() +
  theme(
    legend.position      = "none",
    panel.border         = element_blank(),
    # axis.text.y          = element_blank(),
    axis.ticks           = element_blank(),
    axis.text.x.top      = element_text(angle = -45, vjust = -0.5, hjust = 1)
  ) +
  labs(x = "", y = "")

# ── Plot 3: total hits per assembler (side bar) ───────────────────────────────
assembler_bars <- assembler_counts |>
  ggplot(aes(y = hits_count, x = Assembler)) +
  geom_col(width = 0.6, fill = "gray90") +
  geom_text(aes(label = Assembler), size = 1.5, hjust = 1.5,
            vjust = 0.5, angle = -90, color = "white") +
  scale_y_reverse() +
  my_custom_theme() +
  theme(
    axis.ticks           = element_blank(),
    panel.border         = element_blank(),
    axis.text.x          = element_blank()
  ) +
  coord_cartesian(expand = FALSE) +
  labs(x = NULL, y = "# Conotoxins")

# ── Plot 4: gene superfamily heatmap ─────────────────────────────────────────

# Per-hit coverage + counts across assembler-combination × superfamily
mean_combinations_v2 <- dat |>
  group_by(hits, Assembler, gs_conoServer) |>
  summarise(mean = mean(reference_coverage), n = n_distinct(hits), .groups = "drop")

mean_combinations_v2 <- dat |>
  distinct(Assembler, hits, gs_conoServer) |>
  group_by(hits, gs_conoServer) |>
  summarise(
    combination = paste(sort(unique(Assembler)), collapse = ","),
    .groups = "drop"
  ) |>
  # left_join(mean_combinations_v2, by = "hits") |>
  group_by(combination, gs_conoServer) |>
  summarise(
    # se   = sd(mean) / sum(n),
    n    = n_distinct(hits),
    # mean = mean(mean),
    .groups = "drop"
  )

# Superfamily order (most abundant first)
gs_levs <- mean_combinations_v2 |>
  group_by(gs_conoServer) |>
  tally(n, sort = TRUE) |>
  pull(gs_conoServer)

# ── FIX: dynamic countfactor labels ──────────────────────────────────────────
# pretty() may return 0 as the first break; c(-1, breaks) adds a leading bin
# for zero-count cells.  We need exactly length(breaks) labels — one per
# interval in c(-1, breaks[1], breaks[2], ..., breaks[k]).
breaks <- pretty(mean_combinations_v2$n)
k      <- length(breaks)

# Build labels to match each of the k intervals created by c(-1, breaks):
#   Interval 1 : (-1,       breaks[1]] → "0"        (zero / absent)
#   Interval 2 : (breaks[1], breaks[2]] → "0-100"   (middle bins)
#   ...
#   Interval k : (breaks[k-1], breaks[k]] → ">breaks[k-1]"
labels_auto <- c(
  as.character(breaks[1]),                                    # "0"
  if (k > 2) paste0(breaks[1:(k-2)], "-", breaks[2:(k-1)]), # middle intervals
  paste0(">", breaks[k-1])                                    # open upper bin
)
# labels_auto now has exactly length(breaks) = k elements  ✓

mean_combinations_v2 <- mean_combinations_v2 |>
  mutate(
    gs_conoServer = factor(gs_conoServer, levels = gs_levs),
    countfactor   = cut(n,
                        breaks = c(-1, breaks),
                        labels = labels_auto)
  )

# Dynamic palette: one colour per factor level
n_levels   <- length(labels_auto)
fill_vals  <- rev(scales::pal_brewer(palette = "YlGnBu")(n_levels))
fill_named <- setNames(fill_vals, labels_auto)

mean_chart_v2 <- mean_combinations_v2 |>
  complete(combination, gs_conoServer) |>
  ggplot(aes(y = combination, x = gs_conoServer, fill = countfactor)) +
  geom_tile(color = "white", linewidth = 0.4, na.rm = FALSE) +
  geom_text(aes(label = scales::comma(n)), size = 1.5) +
  scale_y_discrete(limits = combination_lev, labels = NULL) +
  scale_x_discrete(position = "top") +
  scale_fill_manual(
    values   = fill_named,
    na.value = "grey90",
    breaks   = labels_auto          # legend ordered low → high
  ) +
  my_custom_theme(base_size = 10) +
  guides(fill = guide_legend(
    title          = "Number of toxins",
    title.position = "top",
    title.hjust    = 0,
    keywidth       = unit(0.35, "cm"),
    keyheight      = unit(0.35, "cm"),
    override.aes   = list(size = 5)
  )) +
  theme(
    axis.text.x      = element_text(size = 7, angle = -45, vjust = -0.5, hjust = 1),
    legend.position  = "top",
    axis.ticks       = element_blank(),
    panel.border     = element_blank(),
    axis.text.y      = element_blank()
  ) +
  labs(x = "", y = "")

# ── Plot 5: mean reference coverage per combination ───────────────────────────
# gradRGB <- encode_pattern_params_as_hex_colour(pattern_name="gradient",angle=90, 
#                                                colour1="White", colour2="#0570b0")

mean_chart <- dat |>
  group_by(hits) |>
  summarise(
    combination        = paste(sort(unique(Assembler)), collapse = ","),
    reference_coverage = mean(reference_coverage),
    .groups            = "drop"
  ) |> mutate(facet = "Sequence identiy") |>
  ggplot(aes(x = reference_coverage, y = combination, fill = after_stat(x))) + 
  # facet_grid(~facet) +
  scale_y_discrete("", limits = combination_lev, labels = NULL) +
  scale_x_continuous("", position = "top", limits = c(0, 1.15), breaks = c(0, 0.5, 1)) +
  ggridges::geom_density_ridges_gradient(
    scale = 0.8,  alpha=0.33, 
    # fill=gradRGB, 
    colour = alpha(0.1),
    jittered_points = T,
    position = ggridges::position_points_jitter(width = 0.05, height = 0),
    point_shape = '|', point_size = 0.5, point_alpha = 1) +
  # scale_fill_viridis_c(option = "D", direction = -1) +
  my_custom_theme(legend_pos = "none", base_size = 10) +
  theme(
    axis.ticks   = element_blank(),
    panel.border  = element_blank(),
    axis.text.y  = element_blank(),
    axis.text.x.top      = element_text(size = 10)
  ) 

mean_chart

mean_combinations <- dat |>
  group_by(hits) |>
  summarise(
    combination        = paste(sort(unique(Assembler)), collapse = ","),
    reference_coverage = mean(reference_coverage),
    .groups            = "drop"
  ) |>
  group_by(combination) |>
  rstatix::get_summary_stats(type = "mean_sd")

# 
# mean_chart <- mean_combinations |>
#   mutate(min = mean - sd, max = mean + sd) |>
#   ggplot(aes(y = combination, x = mean)) +
#   geom_point(size = 1.5, alpha = 0.5) +
#   geom_errorbar(aes(xmin = min, xmax = max), width = 0.2) +
#   ggstats::geom_stripped_rows() +
#   scale_y_discrete(limits = combination_lev, labels = NULL) +
#   scale_x_continuous(position = "top", limits = c(0, 1.15), breaks = c(0, 0.5, 1)) +
#   my_custom_theme() +
#   theme(
#     axis.ticks   = element_blank(),
#     panel.border = element_blank(),
#     axis.text.y  = element_blank()
#   ) +
#   labs(x = NULL, y = NULL)

# ── Assemble & save figure_3_upset_plot_.png ──────────────────────────────────
PSAVE <- wrap_plots(
  D = mean_chart_v2,
  C = mean_chart,
  A = point_chart,
  design = "CAD"
) +
  plot_layout(
    heights = 1,
    widths  = c(0.1, 0.1, 0.70)
  ) &
  theme(
    plot.margin = margin(1,0,1,0),
    # panel.spacing   = unit(0, "cm")
  )

ggsave(
  PSAVE,
  filename = "figure_3_upset_plot_.png",
  path     = outdir,
  width    = 10,
  height   = 3,
  device   = png,
  dpi      = 800
)

# message("Saved: ", file.path(outdir, "figure_3_upset_plot_.png"))
