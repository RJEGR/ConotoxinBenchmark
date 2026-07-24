# ====================================================================
# LOSO_aggregate_v2.R
#
# Orchestrator script for the LOSO benchmark aggregation pipeline.
#
# Re-scores LOSO assemblies with the same F1@90%-identity full-length
# RBH metric used in manuscript section 3.2 (Behera et al., 2022),
# then layers the three-way Train + Test classification and the
# Circularity Inflation Factor (CIF) on top using the leakage audit
# (Bernett et al., 2024; Grimm et al., 2015; Roberts et al., 2017).
#
# Reads:
#   - <loso_dir>/LOSO_manifest.tsv
#   - <loso_dir>/test_*.fasta                       (for InputNsequences)
#   - <sim_dir>/LOSO_metadata.tsv
#   - <sim_dir>/transrate_contigs_dir/**/contigs.csv
#   - <sim_dir>/leakage_check/leakage_summary.tsv   (optional but
#                                                    required for CIF)
#
# Writes (to --out_dir):
#   - LOSO_per_fold_metrics_v2.tsv         (assembly x threshold; F1 etc.)
#   - LOSO_three_way_classification.tsv    (TP / FP_paralog / FP_spurious)
#   - LOSO_per_superfamily_summary_v2.tsv  (per (sf, tool) summary)
#   - LOSO_cif_per_tool.tsv                (CIF marginal per assembler)
#   - LOSO_leakage_summary_v2.tsv          (per-SF leakage statistics)
#   - LOSO_precision_vs_sensitivity.png
#   - LOSO_F1_by_superfamily.png
#   - LOSO_CIF_per_tool.png
#   - LOSO_three_way_stack.png
#   - LOSO_leakage_distribution.png
#
# Usage:
#   Rscript LOSO_aggregate_v2.R \
#     --loso_dir INPUTS/vfolds_loso_resampling_dir \
#     --sim_dir  INPUTS/vfolds_loso_resampling_dir/transrate_contigs_dir \
#     --out_dir  LOSO_summary_v2 \
#     --threshold 0.90
#
# Dependencies:
#   - LOSO_metrics_core.R   (same directory)
#   - LOSO_io.R             (same directory)
#   - LOSO_viz.R            (same directory)
#
# ====================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(ggplot2)
  library(optparse)
})

# ------------------------------------------------------------------
# Resolve the location of the module files. We assume they live in
# the same directory as this script. When sourced interactively from
# a different CWD, the user can override the modules dir explicitly.
# ------------------------------------------------------------------
.SCRIPT_DIR <- tryCatch(
  dirname(sys.frame(1)$ofile),
  error = function(e) {
    # Fallback: try the commandArgs-based path used when run via Rscript
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) dirname(sub("^--file=", "", file_arg[1])) else "."
  }
)

rm(list = ls())

if(!is.null(dev.list())) dev.off()

setwd("/Users/rjegr/Documents/GitHub/ConotoxinBenchmark/Workflow_dir/0_Data_selection/0_LOSO_dir/LOSO-framework/")

source(file.path(.SCRIPT_DIR, "LOSO_metrics_core.R"))
source(file.path(.SCRIPT_DIR, "LOSO_io.R"))
source(file.path(.SCRIPT_DIR, "LOSO_viz.R"))


# ------------------------------------------------------------------
# CLI
# ------------------------------------------------------------------
option_list <- list(
  optparse::make_option(c("--loso_dir"), type = "character", default = NULL,
              help = "Path to LOSO splits dir (contains LOSO_manifest.tsv)"),
  optparse::make_option(c("--sim_dir"),  type = "character", default = NULL,
              help = "Simulation/assembly output dir (contains transrate_contigs_dir/)"),
  optparse::make_option(c("--out_dir"),  type = "character", default = "LOSO_summary_v2",
              help = "Output directory [default: %default]"),
  optparse::make_option(c("--threshold"), type = "numeric", default = 0.90,
              help = "Reference coverage threshold for F1 (matches manuscript section 3.2) [default: %default]")
)




opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
stopifnot(!is.null(opt$loso_dir), !is.null(opt$sim_dir))


opt$loso_dir <- "/Users/rjegr/Documents/GitHub/ConotoxinBenchmark/INPUTS/vfolds_loso_resampling_dir/"
opt$sim_dir <- "/Users/rjegr/Documents/GitHub/ConotoxinBenchmark/INPUTS/vfolds_loso_resampling_dir/daf3bd43f024d30577d90973fa1544aa_loso_dir/" 
opt$trans_dir <- "/Users/rjegr/Documents/GitHub/ConotoxinBenchmark/INPUTS/vfolds_loso_resampling_dir/transrate_contigs_dir/" 
opt$threshold <- 0.90
opt$out_dir <- "LOSO_summary_v2_dir"

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------
# Step 1: Load LOSO context
# ------------------------------------------------------------------
message("\n[step 1] Loading LOSO registry and metadata...")
manifest      <- load_manifest(opt$loso_dir)
meta          <- load_metadata(opt$sim_dir)
input_n_seq   <- count_test_fasta(opt$loso_dir)
leak          <- load_leakage(opt$trans_dir)

message(sprintf("        %d folds in manifest, %d simulated samples, %d test FASTAs",
                nrow(manifest), nrow(meta), nrow(input_n_seq)))


# ------------------------------------------------------------------
# Step 2: Re-score every assembly with F1 at the manuscript threshold
# ------------------------------------------------------------------
message(sprintf("\n[step 2] Re-scoring assemblies (F1 at reference_coverage >= %.2f)...",
                opt$threshold))

contig_files <- discover_contig_files(opt$trans_dir)

# Map+filter: read each contigs.csv, score it, collect non-null rows.
per_fold_metrics <- purrr::map_dfr(contig_files, function(f) {
  d <- read_contigs(f)
  if (is.null(d)) return(NULL)
  calculate_metrics(d,
                    reference_coverage_val = opt$threshold,
                    input_n_seq_tbl        = input_n_seq)
})

if (nrow(per_fold_metrics) == 0) {
  stop("No assemblies scored successfully. Check TransRate column names.")
}

# Attach (superfamily, n_test, n_train) context
any(unique(per_fold_metrics$sample_id) %in% meta$sample_id)
any(manifest$held_out_superfamily %in% meta$superfamily)

per_fold_metrics <- join_loso_context(per_fold_metrics, manifest, meta)

readr::write_tsv(per_fold_metrics,
                  file.path(opt$out_dir, "LOSO_per_fold_metrics_v2.tsv"))

message(sprintf("        wrote %d per-assembly rows", nrow(per_fold_metrics)))


# ------------------------------------------------------------------
# Step 3: Three-way classification (TP / FP_paralog / FP_spurious)
#         + Circularity Inflation Factor.  Requires the leakage audit.
# ------------------------------------------------------------------
if (!is.null(leak)) {
  message("\n[step 3] Three-way classification + CIF...")

  three_way <- three_way_classify(per_fold_metrics, leak)
  readr::write_tsv(three_way,
                   file.path(opt$out_dir, "LOSO_three_way_classification.tsv"))
  message(sprintf("        wrote three-way table for %d assemblies",
                  nrow(three_way)))

  cif_obj <- compute_cif(three_way)
  readr::write_tsv(cif_obj$per_tool,
                   file.path(opt$out_dir, "LOSO_cif_per_tool.tsv"))
  message("        CIF per tool:")
  print(cif_obj$per_tool, n = Inf)

} else {
  three_way <- NULL
  cif_obj   <- NULL
  message("\n[step 3] SKIPPED -- leakage_summary.tsv unavailable. ",
          "Three-way classification and CIF cannot be computed.")
}


# ------------------------------------------------------------------
# Step 4: Per-superfamily summary
# ------------------------------------------------------------------
message("\n[step 4] Per-superfamily summary...")

per_sf <- per_fold_metrics %>%
  dplyr::group_by(superfamily) %>%
  dplyr::summarise(
    n_replicates      = dplyr::n(),
    mean_n_contigs    = mean(rawcontigs, na.rm = TRUE),
    mean_precision    = mean(precision,   na.rm = TRUE),
    sd_precision      = sd(precision,     na.rm = TRUE),
    mean_sensitivity  = mean(sensitivity, na.rm = TRUE),
    mean_F1           = mean(F1,          na.rm = TRUE),
    sd_F1             = sd(F1,            na.rm = TRUE),
    n_test            = dplyr::first(n_test),
    n_train           = dplyr::first(n_train),
    .groups = "drop"
  )



readr::write_tsv(per_sf,
                 file.path(opt$out_dir, "LOSO_per_superfamily_summary_v2.tsv"))

message(sprintf("        wrote %d (superfamily, tool) rows", nrow(per_sf)))


# ------------------------------------------------------------------
# Step 5: Leakage summary v2
# ------------------------------------------------------------------
if (!is.null(leak)) {
  message("\n[step 5] Per-SF leakage statistics...")
  leak_sf <- leak %>%
    dplyr::group_by(held_out_sf) %>%
    dplyr::summarise(
      n_assemblies          = dplyr::n(),
      median_high_hits      = median(hits_pident90_qcov80, na.rm = TRUE),
      mean_high_hits        = mean(hits_pident90_qcov80,   na.rm = TRUE),
      max_high_hits         = max(hits_pident90_qcov80,    na.rm = TRUE),
      pct_assemblies_clean  = mean(hits_pident90_qcov80 == 0, na.rm = TRUE),
      .groups = "drop"
    )
  readr::write_tsv(leak_sf,
                   file.path(opt$out_dir, "LOSO_leakage_summary_v2.tsv"))
  message("        leakage per held-out SF (median, mean, max, % clean):")
  print(leak_sf, n = Inf)
}


# ------------------------------------------------------------------
# Step 6: Figures
# ------------------------------------------------------------------
message("\n[step 6] Rendering figures...")

# ggplot2::ggplot(per_sf, ggplot2::aes(x = mean_sensitivity, y = mean_precision,
#                                 color = superfamily)) +
#   ggplot2::geom_point(alpha = 0.85, size = 2.4) +
#   ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
#   ggplot2::labs(
#     x     = "Mean sensitivity (TP / InputNsequences)",
#     y     = "Mean precision (TP / (TP + FP))",
#     title = "LOSO precision vs sensitivity at 90% ref. coverage, per held-out superfamily"
#   )

p_ps <- plot_precision_vs_sensitivity(per_sf)

ggplot2::ggsave(file.path(opt$out_dir, "LOSO_precision_vs_sensitivity.png"),
                p_ps, width = 13, height = 4.5, dpi = 300)

ggplot2::ggplot(per_sf, ggplot2::aes(y = superfamily, x = mean_F1)) +
  ggplot2::geom_point(size = 2.5, alpha = 0.85) +
  ggplot2::geom_line(alpha = 0.45) +
  ggplot2::labs(
    x     = "Held-out conotoxin superfamily",
    y     = expression(F[1] ~ "at 90% ref. coverage"),
    title = "LOSO F1 by held-out superfamily, all assemblers"
  ) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
  )
p_f1 <- plot_f1_by_superfamily(per_sf)

ggplot2::ggsave(file.path(opt$out_dir, "LOSO_F1_by_superfamily.png"),
                p_f1, width = 11, height = 5.5, dpi = 300)


ggplot2::ggplot(cif_obj$per_assembly, ggplot2::aes(y = superfamily, x = precision_ood)) +
  ggplot2::geom_jitter()

if (!is.null(cif_obj)) {
  p_cif <- plot_cif_per_tool(cif_obj$per_tool)
  ggplot2::ggsave(file.path(opt$out_dir, "LOSO_CIF_per_tool.png"),
                  p_cif, width = 8, height = 5, dpi = 300)

  p_stk <- plot_three_way_stack(three_way)
  ggplot2::ggsave(file.path(opt$out_dir, "LOSO_three_way_stack.png"),
                  p_stk, width = 13, height = 6, dpi = 300)
}

if (!is.null(leak)) {
  p_lk <- plot_leakage_distribution(leak)
  ggplot2::ggsave(file.path(opt$out_dir, "LOSO_leakage_distribution.png"),
                  p_lk, width = 10, height = 5, dpi = 300)
}

message(sprintf("\nLOSO aggregation v2 complete. Outputs written to: %s",
                opt$out_dir))
