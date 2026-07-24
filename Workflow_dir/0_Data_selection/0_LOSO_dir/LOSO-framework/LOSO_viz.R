# ====================================================================
# LOSO_viz.R
#
# Plot functions for the LOSO benchmark aggregator. All functions
# accept a tibble and return a ggplot object; none write to disk.
# The aggregator runner is responsible for ggsave() calls.
#
# Public API:
#   theme_loso(...)
#   recode_tool_factor(tbl, col = tool)
#   plot_precision_vs_sensitivity(per_sf_tbl)
#   plot_f1_by_superfamily(per_sf_tbl)
#   plot_cif_per_tool(cif_per_tool_tbl)
#   plot_three_way_stack(three_way_tbl)
#   plot_leakage_distribution(leak_tbl)
#
# ====================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(forcats)
  library(ggplot2)
})


# Display names + a stable colour palette for the assemblers reported
# in this benchmark. Tools not listed here render in grey by default.
.TOOL_RECODE <- c(
  STRINGTIE  = "StringTie",
  SPADES     = "rnaSPAdes",
  TRINITY    = "Trinity",
  IDBA       = "IDBA",
  MEGAHIT    = "MEGAHIT",
  RNABLOOM   = "RNA-Bloom",
  BRIDGER    = "Bridger",
  TRANSABBYS = "Trans-ABySS",
  BINPACKER  = "BinPacker",
  SOAPDENOVO = "SOAP-denovo",
  CSTONE     = "Cstone",
  TRANSLIG   = "TransLiG",
  Baseline   = "Baseline",
  PLASS      = "PinguiN (nuclassemble)"
)

.TOOL_COLOURS <- function() {
  pal <- c(
    "#5B1A18", "#FD6467", "#F1BB7B", "#D67236", "#7294D4",
    "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
    "#D55E00", "#CC79A7", "#999999", "#000000"
  )
  setNames(pal[seq_along(.TOOL_RECODE)], unname(.TOOL_RECODE))
}


theme_loso <- function(base_size = 10, legend_pos = "top", ...) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      legend.position    = legend_pos,
      strip.background   = ggplot2::element_rect(fill = "gray90", color = "white"),
      strip.text         = ggplot2::element_text(hjust = 0),
      axis.text          = ggplot2::element_text(color = "black"),
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      ...
    )
}


# Recode a `tool` column (uppercase keys) into pretty display labels
# and lock the factor levels to the .TOOL_RECODE order so plots stay
# consistent across panels.
recode_tool_factor <- function(tbl, col = "tool") {
  tbl[[col]] <- dplyr::recode(tbl[[col]], !!!.TOOL_RECODE)
  tbl[[col]] <- factor(tbl[[col]], levels = unname(.TOOL_RECODE))
  tbl
}


# -----------------------------------------------------------------
# plot_precision_vs_sensitivity(per_sf_tbl)
#
# Scatter of precision vs sensitivity at the manuscript threshold
# (90% reference coverage RBH), faceted by tool. Replaces the older
# proxy-based plot (mean_ref_cov on x) with the actual sensitivity
# metric TP / InputNsequences.
#
# Expects columns: tool, superfamily, mean_precision, mean_sensitivity.
# -----------------------------------------------------------------
plot_precision_vs_sensitivity <- function(per_sf_tbl) {

  d <- recode_tool_factor(per_sf_tbl)

  ggplot2::ggplot(d, ggplot2::aes(x = mean_sensitivity, y = mean_precision,
                                   color = tool)) +
    ggplot2::facet_wrap(~tool, nrow = 1) +
    ggplot2::geom_point(alpha = 0.85, size = 2.4) +
    ggplot2::scale_color_manual("", values = .TOOL_COLOURS(),
                                 na.value = "grey60") +
    ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
    ggplot2::labs(
      x     = "Mean sensitivity (TP / InputNsequences)",
      y     = "Mean precision (TP / (TP + FP))",
      title = "LOSO precision vs sensitivity at 90% ref. coverage, per held-out superfamily"
    ) +
    theme_loso(legend_pos = "none")
}


# -----------------------------------------------------------------
# plot_f1_by_superfamily(per_sf_tbl)
#
# Dot-and-line plot of F1@90% per held-out superfamily, one line per
# tool. Reveals which superfamilies are hardest to assemble OOD.
# Expects: tool, superfamily, mean_F1.
# -----------------------------------------------------------------
plot_f1_by_superfamily <- function(per_sf_tbl) {

  d <- recode_tool_factor(per_sf_tbl) %>%
    dplyr::mutate(
      superfamily = forcats::fct_reorder(superfamily, mean_F1, .desc = TRUE)
    )

  ggplot2::ggplot(d, ggplot2::aes(x = superfamily, y = mean_F1,
                                   color = tool, group = tool)) +
    ggplot2::geom_point(size = 2.5, alpha = 0.85) +
    ggplot2::geom_line(alpha = 0.45) +
    ggplot2::scale_color_manual("", values = .TOOL_COLOURS(),
                                 na.value = "grey60") +
    ggplot2::labs(
      x     = "Held-out conotoxin superfamily",
      y     = expression(F[1] ~ "at 90% ref. coverage"),
      title = "LOSO F1 by held-out superfamily, all assemblers"
    ) +
    theme_loso() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
}


# -----------------------------------------------------------------
# plot_cif_per_tool(cif_per_tool_tbl)
#
# Bar plot of median CIF per assembler. CIF == 1 (dashed reference)
# means random-CV would have given the same precision as LOSO; bars
# above 1 quantify the random-CV inflation.
# -----------------------------------------------------------------
plot_cif_per_tool <- function(cif_per_tool_tbl) {

  d <- recode_tool_factor(cif_per_tool_tbl) %>%
    dplyr::mutate(tool = forcats::fct_reorder(tool, median_CIF, .desc = TRUE))

  ggplot2::ggplot(d, ggplot2::aes(x = tool, y = median_CIF, fill = tool)) +
    ggplot2::geom_col(width = 0.65) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "grey30") +
    ggplot2::scale_fill_manual("", values = .TOOL_COLOURS(),
                                na.value = "grey60") +
    ggplot2::labs(
      x     = "Assembler",
      y     = "Median Circularity Inflation Factor",
      title = "CIF per assembler  (precision_combined / precision_OOD)",
      subtitle = "CIF = 1 means random-CV would have given the same precision as LOSO"
    ) +
    theme_loso(legend_pos = "none") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
}


# -----------------------------------------------------------------
# plot_three_way_stack(three_way_tbl)
#
# Stacked-bar plot of the per-assembly three-way classification:
# TP (test hit), FP_paralog (train-only hit), FP_spurious (neither).
# X axis = assembly, faceted by tool.
# -----------------------------------------------------------------
plot_three_way_stack <- function(three_way_tbl) {

  d <- three_way_tbl %>%
    dplyr::select(sample_id, tool, TP, FP_paralog, FP_spurious) %>%
    tidyr::pivot_longer(c(TP, FP_paralog, FP_spurious),
                        names_to = "class", values_to = "n")
  d <- recode_tool_factor(d)

  d$class <- factor(d$class, levels = c("TP", "FP_paralog", "FP_spurious"))

  ggplot2::ggplot(d, ggplot2::aes(x = sample_id, y = n, fill = class)) +
    ggplot2::facet_wrap(~tool, scales = "free", nrow = 2) +
    ggplot2::geom_col(position = "fill", width = 0.85) +
    ggplot2::scale_fill_manual(
      values = c(TP = "#0072B2", FP_paralog = "#E69F00", FP_spurious = "grey60"),
      labels = c(TP          = "TP (test SF)",
                 FP_paralog  = "FP-paralog (train-only)",
                 FP_spurious = "FP-spurious (neither)")
    ) +
    ggplot2::labs(
      x     = "Held-out fold (sample_id)",
      y     = "Proportion of contigs",
      title = "Three-way classification of LOSO contigs",
      fill  = NULL
    ) +
    theme_loso() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank())
}


# -----------------------------------------------------------------
# plot_leakage_distribution(leak_tbl)
#
# Diagnostic plot of leakage hit counts per held-out superfamily.
# A leak-free LOSO run should be near zero everywhere except where
# genuine paralog overlap is documented in the literature.
# Expects: held_out_sf, hits_pident90_qcov80.
# -----------------------------------------------------------------
plot_leakage_distribution <- function(leak_tbl) {

  d <- leak_tbl %>%
    dplyr::mutate(
      held_out_sf = forcats::fct_reorder(held_out_sf,
                                          hits_pident90_qcov80, median)
    )

  ggplot2::ggplot(d, ggplot2::aes(x = held_out_sf,
                                   y = hits_pident90_qcov80)) +
    ggplot2::geom_boxplot(fill = "grey85", outlier.size = 1) +
    ggplot2::labs(
      x     = "Held-out superfamily",
      y     = "Contigs with high-identity hit to TRAINING reference",
      title = "Leakage audit: BLAST against train_<sf>.fasta",
      subtitle = paste("pident >= 90, qcov >= 80; non-zero values usually",
                       "indicate genuine paralog overlap")
    ) +
    theme_loso() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}
