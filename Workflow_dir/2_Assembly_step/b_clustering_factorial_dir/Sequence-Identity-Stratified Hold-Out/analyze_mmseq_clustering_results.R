#!/usr/bin/env Rscript
# =============================================================================
# analyze_mmseq_clustering_results.R  (Strategy 3)
#
# Sequence-Identity-Stratified Hold-Out analysis.
#
# Strategy 3 reframes the question. The ORIGINAL script plotted the number of
# clusters vs identity threshold (a descriptive measure of ConoServer
# redundancy). This version plots ASSEMBLER RECOVERY as a function of each
# held-out sequence's Nearest-Reference-Identity (NRI) -- the continuous
# divergence axis that defuses StringTie / ConoSorter circularity.
#
# Inputs:
#   --strat3_dir   directory from cluster_toxins_mmseqs.sh; must contain
#                  nearest_reference_identity.tsv (and optionally the
#                  protein NRI table + <base>_cluster_counts.csv).
#   --runlog_glob  glob for per-fold runlog.tsv files written by Assemblers.sh
#                  (fold_tag, sample_id, tool, ..., assembly_fa, contigs_csv).
#   --recovery_tsv (optional) precomputed recovery table with columns
#                  test_id, tool, recovered (0/1), precision. If absent, the
#                  script derives recovery from the TransRate contigs.csv
#                  referenced in the runlogs.
#   --out_dir      output directory (default: <strat3_dir>/strat3_analysis).
#
# Output:
#   strat3_recovery_vs_NRI.png/.pdf   headline figure (logistic fit per tool)
#   strat3_recovery_binned.tsv        recovery per NRI bin per tool
#   strat3_logistic_coefficients.tsv  fitted slope / NRI50 per tool
#   strat3_diversity_curve.png        supplementary (original cluster curve)
#
# Usage:
#   Rscript analyze_mmseq_clustering_results.R \
#       --strat3_dir curated_nuc_conoServerDB_strat3_dir \
#       --runlog_glob '*_runlog.tsv'
# =============================================================================

rm(list = ls())


options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(optparse)
})

# ----------------------------- arguments ------------------------------------

option_list <- list(
  make_option("--strat3_dir",   type = "character", default = NULL,
              help = "Directory from cluster_toxins_mmseqs.sh (Strategy 3)"),
  make_option("--rate_dir",     type = "character", default = NULL,
              help = "Directory containing TransRate contigs.csv files (recursive search)"),
  make_option("--recovery_tsv", type = "character", default = NULL,
              help = "Optional precomputed recovery table (threshold,tool,test_id,recovered,precision)"),
  make_option("--out_dir",      type = "character", default = NULL,
              help = "Output directory [default <strat3_dir>/strat3_analysis]"),
  make_option("--id_threshold", type = "double",    default = 0.90,
              help = "Min reference coverage to call a held-out seq 'recovered' [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))


opt$strat3_dir <- "~/Documents/GitHub/ConotoxinBenchmark/1_assembly_dir/sequence_identity_stratified_dir/curated_nuc_conoServerDB_strat3_dir/"

if (is.null(opt$strat3_dir)) stop("--strat3_dir is required")

if (is.null(opt$out_dir))
  opt$out_dir <- file.path(opt$strat3_dir, "strat3_analysis")

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------- theme ----------------------------------------
my_theme <- function(base_size = 12, legend_pos = "top", ...) {
  theme_bw(base_size = base_size) +
    theme(legend.position = legend_pos,
          strip.background = element_rect(fill = "gray90", color = "white"),
          strip.text = element_text(hjust = 0),
          panel.grid.minor = element_blank(), ...)
}

# ============================================================================
# 1. Load the NRI table (the Strategy 3 x-axis)
# ============================================================================
nri_path <- file.path(opt$strat3_dir, "nearest_reference_identity.tsv")

if (!file.exists(nri_path))
  stop("NRI table not found: ", nri_path,
       "\nRun cluster_toxins_mmseqs.sh first.")

nri <- read_tsv(nri_path) %>%
  # nri_nuc is a fraction 0-1; convert to percent for readability
  mutate(nri_pct = nri_nuc * 100)

cat(sprintf("Loaded NRI for %d held-out sequences across %d thresholds\n",
            n_distinct(nri$test_id), n_distinct(nri$threshold)))

# Optional protein NRI (robustness axis)
prot_path <- file.path(opt$strat3_dir, "nearest_reference_identity_protein.tsv")
nri_prot <- if (file.exists(prot_path)) {
  read_tsv(prot_path) %>% mutate(nri_prot_pct = nri_prot * 100) %>%
    select(test_id, nri_prot, nri_prot_pct, threshold)
} else NULL

# ============================================================================
# 2. Build the recovery table from TransRate contigs.csv files
#
# Directory layout assumed (two variants, both handled):
#   <rate_dir>/<tool>/<fold_tag>/contigs.csv   (two-level)
#   <rate_dir>/<tool>_<fold_tag>/contigs.csv   (one-level)
#
# tool      <- first path component below rate_dir
# threshold <- numeric value extracted from the fold_tag
#              (e.g. "fold_0.90", "0.90", "90" are all recognised)
#
# The recovery table produced has columns:
#   threshold, tool, test_id, best_cov, precision, recovered
# ============================================================================

if (is.null(opt$rate_dir))
  opt$rate_dir <- "~/Documents/GitHub/ConotoxinBenchmark/1_assembly_dir/sequence_identity_stratified_dir/transrate_contigs_dir/"

rate_dir_exp <- normalizePath(path.expand(opt$rate_dir), mustWork = FALSE)


# Column names vary across TransRate versions
REFCOV_COLS <- c("p_unique_ref", "reference_coverage", "p_refcov")
REFHIT_COLS <- c("hits", "reference", "CRBB")

# Parse tool name and identity threshold from a contigs.csv path
parse_path_meta <- function(path) {
  rel   <- sub(paste0("^", rate_dir_exp, "/?"), "",
               normalizePath(path, mustWork = FALSE))
  parts <- strsplit(rel, .Platform$file.sep)[[1]]
  parts <- parts[nzchar(parts) & parts != "contigs.csv"]

  if (length(parts) >= 2) {
    # two-level layout: parts[1] = tool, parts[2] = fold tag
    tool_raw <- parts[1]; fold_str <- parts[2]
  } else if (length(parts) == 1) {
    # one-level layout: split on the last numeric suffix, e.g. Trinity_0.90
    m        <- regmatches(parts[1], regexpr("_([0-9]+\\.?[0-9]*)$", parts[1]))
    tool_raw <- if (length(m)) sub("_[0-9]+\\.?[0-9]*$", "", parts[1]) else parts[1]
    fold_str <- if (length(m)) m else NA_character_
  } else {
    tool_raw <- NA_character_; fold_str <- NA_character_
  }

  thr_m <- regmatches(fold_str, regexpr("[0-9]+\\.?[0-9]*", fold_str))
  thr   <- if (length(thr_m)) {
    v <- as.numeric(thr_m)
    if (!is.na(v) && v > 1) v / 100 else v   # normalise percent -> fraction
  } else NA_real_

  list(tool = tool_raw, threshold = thr)
}

# Read one contigs.csv and return tidy hits (tool, threshold, ref_hit, ref_cov)
read_one_contig_file <- function(path) {
  meta <- parse_path_meta(path)
  ct   <- tryCatch(
    suppressWarnings(read_csv(path, show_col_types = FALSE)),
    error = function(e) { warning("Cannot read ", path, ": ", e$message); NULL }
  )
  if (is.null(ct) || nrow(ct) == 0) return(NULL)

  cov_col <- intersect(REFCOV_COLS, colnames(ct))[1]
  hit_col <- intersect(REFHIT_COLS, colnames(ct))[1]

  if (is.na(cov_col)) {
    warning("No reference-coverage column in: ", path,
            "\n  Available columns: ", paste(colnames(ct), collapse = ", "))
    return(NULL)
  }

  ct %>%
    transmute(
      tool      = meta$tool,
      threshold = meta$threshold,
      ref_hit   = if (!is.na(hit_col)) .data[[hit_col]] else NA_character_,
      ref_cov   = suppressWarnings(as.numeric(.data[[cov_col]]))
    ) %>%
    filter(!is.na(ref_hit), !is.na(ref_cov))
}

if (!is.null(opt$recovery_tsv) && file.exists(opt$recovery_tsv)) {

  cat("Using precomputed recovery table:", opt$recovery_tsv, "\n")
  recovery <- read_tsv(opt$recovery_tsv)

} else {

  contig_files <- list.files(rate_dir_exp, pattern = "contigs\\.csv$",
                             recursive = TRUE, full.names = TRUE)
  # fallback: look in PWD (where Assemblers.sh may write directly)
  if (length(contig_files) == 0)
    contig_files <- list.files("transrate_contigs_dir", pattern = "contigs\\.csv$",
                               recursive = TRUE, full.names = TRUE)
  if (length(contig_files) == 0)
    stop("No contigs.csv files found under: ", rate_dir_exp,
         "\nPass --rate_dir or --recovery_tsv.")

  cat(sprintf("Found %d contigs.csv file(s) under %s\n",
              length(contig_files), rate_dir_exp))

  hits <- map_dfr(contig_files, read_one_contig_file)

  if (nrow(hits) == 0)
    stop("All contigs.csv files were empty or lacked a recognised reference-coverage ",
         "column (checked: ", paste(REFCOV_COLS, collapse = ", "), ").")

  cat(sprintf("  Parsed %d contig-reference hit rows\n", nrow(hits)))

  # Diagnose any files where threshold or tool could not be parsed
  missing_meta <- hits %>% filter(is.na(tool) | is.na(threshold))
  if (nrow(missing_meta) > 0)
    warning(sprintf(
      "%d rows have NA tool/threshold — check directory naming under rate_dir.",
      nrow(missing_meta)))

  # A test sequence is 'recovered' if ≥1 assembled contig maps to it with
  # reference coverage >= id_threshold.
  recovery <- hits %>%
    mutate(test_id = str_remove(ref_hit, "\\|.*$")) %>%   # strip "|superfamily" tag
    group_by(tool, threshold, test_id) %>%
    summarise(
      best_cov  = max(ref_cov,  na.rm = TRUE),
      precision = mean(ref_cov, na.rm = TRUE),
      .groups   = "drop"
    ) %>%
    mutate(recovered = as.integer(best_cov >= opt$id_threshold))
}

cat(sprintf("Recovery table: %d rows  |  %d tool(s)  |  %d threshold(s)\n",
            nrow(recovery), n_distinct(recovery$tool),
            n_distinct(recovery$threshold)))

# ============================================================================
# 3. Join recovery to NRI -> the Strategy 3 analysis table
# ============================================================================


if (!is.null(nri_prot))
  nri_df <- nri_prot |>
  inner_join(nri %>% select(test_id, threshold, nri_nuc, nri_pct), by = c("test_id", "threshold")) else {
    nri_df <- nri %>% select(test_id, threshold, nri_nuc, nri_pct)
    
    nri_df |>
      ggplot(aes(nri_prot, nri_nuc)) + geom_point()
  }

    


analysis_tbl <- recovery %>%
  inner_join(nri_df, by = c("test_id", "threshold"))

if (nrow(analysis_tbl) == 0)
  stop("No overlap between recovery test_ids and NRI test_ids. ",
       "Check that simulation used the test_*.fasta from --strat3_dir.")

write_tsv(analysis_tbl, file.path(opt$out_dir, "strat3_analysis_table.tsv"))

# ============================================================================
# 4. Bin by NRI and summarise recovery
# ============================================================================
nri_breaks <- c(0, 40, 50, 60, 70, 80, 90, 100)

binned <- analysis_tbl %>%
  mutate(nri_bin = cut(nri_pct, breaks = nri_breaks,
                       include.lowest = TRUE, right = FALSE)) %>%
  group_by(tool, nri_bin) %>%
  summarise(n             = n(),
            n_recovered   = sum(recovered),
            sensitivity   = mean(recovered),
            se            = sqrt(sensitivity * (1 - sensitivity) / n()),
            mean_precision = mean(precision, na.rm = TRUE),
            .groups = "drop")

write_tsv(binned, file.path(opt$out_dir, "strat3_recovery_binned.tsv"))

# ============================================================================
# 5. Fit a logistic recovery curve per assembler
#    recovered ~ NRI ; report slope and NRI50 (the NRI at 50% recovery).
# ============================================================================
fit_logistic <- function(df) {
  if (n_distinct(df$recovered) < 2 || nrow(df) < 10)
    return(tibble(intercept = NA, slope = NA, nri50 = NA, n = nrow(df)))
  m <- tryCatch(glm(recovered ~ nri_pct, data = df, family = binomial()),
                error = function(e) NULL)
  if (is.null(m)) return(tibble(intercept = NA, slope = NA, nri50 = NA, n = nrow(df)))
  co <- coef(m)
  tibble(intercept = co[1], slope = co[2],
         nri50 = -co[1] / co[2],          # NRI where P(recovery) = 0.5
         n = nrow(df))
}
logit_coef <- analysis_tbl %>%
  group_by(tool) %>% group_modify(~ fit_logistic(.x)) %>% ungroup()

write_tsv(logit_coef, file.path(opt$out_dir, "strat3_logistic_coefficients.tsv"))

cat("\n=== Logistic recovery fits (NRI50 = divergence at 50% recovery) ===\n")
print(logit_coef, n = Inf)

# ============================================================================
# 6. Headline figure: recovery vs NRI
# ============================================================================
# smooth logistic prediction grid
pred_grid <- analysis_tbl %>%
  group_by(tool) %>%
  group_modify(~{
    if (n_distinct(.x$recovered) < 2 || nrow(.x) < 10)
      return(tibble(nri_pct = numeric(0), fit = numeric(0)))
    m <- glm(recovered ~ nri_pct, data = .x, family = binomial())
    g <- tibble(nri_pct = seq(0, 100, by = 1))
    g$fit <- predict(m, newdata = g, type = "response")
    g
  }) %>% ungroup()

p_main <- ggplot() +
  # binned empirical sensitivity with binomial SE
  geom_point(data = binned,
             aes(x = (as.numeric(nri_bin) * 10 + 30), y = sensitivity,
                 color = tool, size = n), alpha = 0.7) +
  geom_errorbar(data = binned,
                aes(x = (as.numeric(nri_bin) * 10 + 30),
                    ymin = pmax(0, sensitivity - se),
                    ymax = pmin(1, sensitivity + se),
                    color = tool), width = 1.5, linewidth = 0.3) +
  # logistic fit
  geom_line(data = pred_grid,
            aes(x = nri_pct, y = fit, color = tool), linewidth = 0.8) +
  scale_size_continuous(range = c(1, 4), name = "n held-out") +
  labs(x = "Nearest-reference identity of held-out toxin (%)",
       y = "Assembler sensitivity (recovery rate)",
       color = "Assembler",
       title = "Strategy 3: assembler recovery vs divergence from reference",
       subtitle = "Each held-out conotoxin scored against its closest training-set neighbour") +
  coord_cartesian(ylim = c(0, 1)) +
  my_theme(legend_pos = "right")

ggsave(file.path(opt$out_dir, "strat3_recovery_vs_NRI.png"),
       p_main, width = 9, height = 5.5, dpi = 300)
ggsave(file.path(opt$out_dir, "strat3_recovery_vs_NRI.pdf"),
       p_main, width = 9, height = 5.5)

# ============================================================================
# 7. Supplementary: original diversity-reduction curve (kept for continuity)
# ============================================================================
cc_path <- list.files(opt$strat3_dir, pattern = "_cluster_counts\\.csv$",
                       full.names = TRUE)
if (length(cc_path) >= 1) {
  cc <- read_csv(cc_path[1]) %>%
    filter(num_clusters != "ERROR") %>%
    mutate(across(c(sequence_identity, num_clusters), as.numeric))
  p_div <- ggplot(cc, aes(sequence_identity, num_clusters)) +
    geom_line() + geom_point(shape = 21, fill = "white") +
    scale_x_reverse() +
    labs(x = "Sequence identity threshold",
         y = "Number of conotoxin clusters",
         title = "Supplementary: ConoServer diversity-reduction curve") +
    my_theme()
  ggsave(file.path(opt$out_dir, "strat3_diversity_curve.png"),
         p_div, width = 5, height = 4, dpi = 300)
}

cat(sprintf("\nStrategy 3 analysis complete. Outputs in: %s\n", opt$out_dir))
