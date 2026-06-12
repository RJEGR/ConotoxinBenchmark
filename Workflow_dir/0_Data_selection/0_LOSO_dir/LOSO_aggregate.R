# ====================================================================
# LOSO_aggregate.R
#
# Aggregates per-fold metrics from a LOSO run into a per-superfamily
# performance report. Companion script to Run_LOSO_pipeline.sh.
#
# Reads:
#   - vfolds_loso_resampling_dir/LOSO_manifest.tsv     (fold registry)
#   - <sim_dir>/LOSO_metadata.tsv                      (sample <-> SF mapping)
#   - <sim_dir>/transrate_contigs_dir/**/contigs.csv   (per-assembly metrics)
#   - <sim_dir>/leakage_check/leakage_summary.tsv      (optional)
#
# Writes (to --out_dir):
#   - LOSO_per_fold_metrics.tsv
#   - LOSO_per_superfamily_summary.tsv
#   - LOSO_in_vs_out_distribution.tsv
#   - LOSO_sensitivity_by_superfamily.png
#   - LOSO_precision_vs_sensitivity.png
#
# Usage:
#   Rscript LOSO_aggregate.R \
#     --loso_dir INPUTS/vfolds_loso_resampling_dir \
#     --sim_dir  <md5>_loso_dir \
#     --out_dir  <md5>_loso_dir/LOSO_summary
#
# ====================================================================

dir <- "~/Documents/GitHub/ConotoxinBenchmark/INPUTS/vfolds_loso_resampling_dir/"

setwd(dir)

suppressPackageStartupMessages({
  library(tidyverse)
  library(optparse)
})

option_list <- list(
  make_option(c("--loso_dir"), type = "character", default = NULL,
              help = "Path to LOSO splits dir (contains LOSO_manifest.tsv)"),
  make_option(c("--sim_dir"),  type = "character", default = NULL,
              help = "Simulation output dir from Simulate_rnaseq_loso.sh"),
  make_option(c("--out_dir"),  type = "character", default = "LOSO_summary",
              help = "Output directory for the aggregated report")
)

opt <- parse_args(OptionParser(option_list = option_list))

stopifnot(!is.null(opt$loso_dir), !is.null(opt$sim_dir))

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)


theme_loso <- function(base_size = 10, legend_pos = "top", ...) {
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


# ---- Load LOSO registry + simulation metadata ----

opt$loso_dir <- dir
opt$sim_dir <- file.path(dir, "daf3bd43f024d30577d90973fa1544aa_loso_dir")

manifest <- read_tsv(file.path(opt$loso_dir, "LOSO_manifest.tsv"))
meta     <- read_tsv(file.path(opt$sim_dir,  "LOSO_metadata.tsv"))

cat(sprintf("Loaded %d LOSO folds and %d simulated datasets\n",
            nrow(manifest), nrow(meta)))

# ---- Discover per-assembly contig metrics from TransRate ----
# Assemblers.sh writes to transrate_contigs_dir/<assembly>_dir/contigs.csv
contig_files <- list.files(
  path       = "transrate_contigs_dir", #file.path(opt$sim_dir),
  pattern    = "contigs\\.csv$",
  recursive  = TRUE,
  full.names = TRUE
)
# Fallback: also look in PWD/transrate_contigs_dir (where Assemblers.sh writes)
if (length(contig_files) == 0) {
  contig_files <- list.files("transrate_contigs_dir",
                             pattern = "contigs\\.csv$",
                             recursive = TRUE, full.names = TRUE)
}

if (length(contig_files) == 0) {
  stop("No contigs.csv found. Did Assemblers.sh + TransRate run successfully?")
}
cat(sprintf("Found %d TransRate contigs.csv files\n", length(contig_files)))

# ---- Per-assembly summary ----
# TransRate's contigs.csv has per-contig rows; we collapse to per-assembly metrics.
# Key columns: contig_name, p_good (1 = good contig), CRBB (boolean reference hit)
# Sensitivity = (#test refs hit by >=1 contig) / (#test refs)
# Precision   = (#contigs with high-id ref hit) / (#total contigs)
# These map to TransRate's reference-based scoring.
summarise_assembly <- function(csv_path) {

  df <- tryCatch(suppressWarnings(read_csv(csv_path, show_col_types = FALSE)),
                 error = function(e) NULL)
  
  if (is.null(df) || nrow(df) == 0) return(NULL)

  # Identify the source assembly + LOSO sample id from the path
  asm_dir <- basename(dirname(csv_path))
  asm_bs  <- sub("_dir$", "", asm_dir)
  # asm_bs example: test_M_50x_PE_Trinity  -> sample_id = test_M_50x_PE, tool = Trinity
  parts <- strsplit(asm_bs, "_")[[1]]
  if (length(parts) < 2) return(NULL)
  tool      <- tail(parts, 1)
  sample_id <- paste(head(parts, -1), collapse = "_")

  # Per-assembly metrics
  n_contigs   <- nrow(df)
  # TransRate >= 1.0 columns may vary; guard with checks
  has_crbb    <- "CRBB" %in% colnames(df) || "has_crb" %in% colnames(df)
  crbb_col    <- if ("CRBB" %in% colnames(df)) "CRBB" else
                 if ("has_crb" %in% colnames(df)) "has_crb" else NA
  n_ref_hit   <- if (!is.na(crbb_col)) sum(as.logical(df[[crbb_col]]), na.rm = TRUE) else NA
  pct_ref_hit <- if (!is.na(n_ref_hit)) n_ref_hit / n_contigs else NA
  
  
  return(df)
  # df |> filter(length >= 200) %>%
  #   # group_by(vfold_set, Assembler) %>%
  #   calculate_metrics(reference_coverage_val = 0.9) %>%

  # Reference-coverage proxy for sensitivity
  ref_cov_col <- if ("p_unique_ref"   %in% colnames(df)) "p_unique_ref"   else
                 if ("reference_coverage"   %in% colnames(df)) "reference_coverage"   else NA
  
  mean_ref_cov <- if (!is.na(ref_cov_col)) mean(df[[ref_cov_col]], na.rm = TRUE) else NA
  
  sd_ref_cov <- if (!is.na(ref_cov_col)) sd(df[[ref_cov_col]], na.rm = TRUE) else NA

  tibble(
    sample_id     = sample_id,
    tool          = tool,
    n_contigs     = n_contigs,
    n_ref_hit     = n_ref_hit,
    precision_est = pct_ref_hit,
    mean_ref_cov  = mean_ref_cov,
    sd_ref_cov    = sd_ref_cov,
    csv_path      = csv_path
  )
}

asm_summary <- map_dfr(contig_files, summarise_assembly)


# replace with calculate_metrics()

read_csv <- function(csv_path) {
  
  df <- tryCatch(suppressWarnings(read_csv(csv_path, show_col_types = FALSE)),
                 error = function(e) NULL)
  
  if (is.null(df) || nrow(df) == 0) return(NULL) }


transratedf <- lapply(file.path(opt$sim_dir, contig_files), read_csv)


if (nrow(asm_summary) == 0) {
  stop("Could not extract any metrics from TransRate contigs.csv files. Check column names.")
}

# ---- Join to LOSO metadata ----
loso_metrics <- asm_summary %>%
  left_join(meta, by = "sample_id") %>%
  left_join(manifest |> 
              select(held_out_superfamily, n_test, n_train) |> 
              mutate(held_out_superfamily = gsub("[[:space:]]", "_", held_out_superfamily)),
            by = c("superfamily" = "held_out_superfamily"))

write_tsv(loso_metrics, file.path(opt$out_dir, "LOSO_per_fold_metrics.tsv"))

# ---- Per-superfamily summary ----
per_sf <- loso_metrics %>%
  group_by(superfamily, tool) %>%
  summarise(
    n_replicates  = n(),
    mean_n_contigs   = mean(n_contigs, na.rm = TRUE),
    mean_precision   = mean(precision_est, na.rm = TRUE),
    sd_precision     = sd(precision_est, na.rm = TRUE),
    mean_ref_cov     = mean(mean_ref_cov, na.rm = TRUE),
    n_test           = first(n_test),
    n_train          = first(n_train),
    .groups = "drop"
  )

write_tsv(per_sf, file.path(opt$out_dir, "LOSO_per_superfamily_summary.tsv"))

# ---- In-distribution vs out-of-distribution comparison ----
# "In-distribution" baseline = original 12-fold random CV results (the user's
# previous benchmark). If those aren't provided here, this section reports only
# the LOSO marginal as the OOD estimate. The user can paste the in-dist column
# from their existing accuracy table to complete this.
ood_summary <- per_sf %>%
  group_by(tool) %>%
  summarise(
    mean_precision_ood = mean(mean_precision, na.rm = TRUE),
    sd_precision_ood   = sd(mean_precision,   na.rm = TRUE),
    n_superfamilies    = n(),
    .groups = "drop"
  )

write_tsv(ood_summary, file.path(opt$out_dir, "LOSO_in_vs_out_distribution.tsv"))

# ---- Plots ----
# Per-superfamily precision (the headline LOSO figure)
p1 <- per_sf %>%
  mutate(superfamily = fct_reorder(superfamily, mean_precision, .desc = TRUE)) %>%
  ggplot(aes(x = superfamily, y = mean_precision, fill = tool)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = pmax(0, mean_precision - sd_precision),
                    ymax = mean_precision + sd_precision),
                position = position_dodge(width = 0.8), width = 0.25) +
  labs(x = "Held-out superfamily (LOSO)",
       y = "Mean precision",
       title = "LOSO precision per held-out conotoxin superfamily",
       fill  = "Assembler") +
  theme_loso() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(opt$out_dir, "LOSO_sensitivity_by_superfamily.png"),
       p1, width = 10, height = 5.5, dpi = 300)

# Precision vs sensitivity (proxy: ref coverage)

recode_to <- c("STRINGTIE","SPADES", "TRINITY","IDBA", "MEGAHIT", "RNABLOOM" ,"BRIDGER", "TRANSABBYS", "BINPACKER","SOAPDENOVO" ,"CSTONE", "TRANSLIG", "Baseline", "PLASS")

recode_to <- structure(c("StringTie","rnaSPAdes", "Trinity", "IDBA", "MEGAHIT", "RNA-Bloom","BRIDGER","Trans-ABySS", "BinPacker", "SOAP-denovo", "Cstone", "TransLiG", "Baseline", "PinguiN (nuclassemble)"), names = recode_to)

n <- length(recode_to)

scale_col <- c(ggsci::pal_startrek()(7), ggsci::pal_cosmic()(n-7))

scale_col <- structure(scale_col, names = sort(recode_to))



p2 <- per_sf %>%
  dplyr::mutate(tool = dplyr::recode_factor(tool, !!!recode_to)) |>
  ggplot(aes(x = mean_ref_cov, y = mean_precision, color = tool)) +
  # facet_grid(~tool) +
  geom_point(alpha = 0.8, size = 2.5) +
  scale_color_manual("", values = scale_col) +
  # scale_fill_manual("", values = scale_col) +
  # geom_text(aes(label = superfamily), size = 2.5, nudge_y = 0.02, show.legend = FALSE) +
  # ggrepel::geom_text_repel(aes(label = superfamily), 
  #                          max.overlaps = Inf,
  #                          # nudge_x      = -0.05,
  #                          nudge_y = 0.02,
  #                          # xlim = c(1.1, 5),
  #                          # ylim = c(0.80, 2),
  #                          # direction    = "y",
  #                          # hjust        = -0.5,
  #                          # vjust = -1, 
  #                          min.segment.length = 0,
  #                          segment.curvature = 0.1, 
  #                          segment.size = 0,
  #                          size = 3, family = "GillSans") +
  labs(x = "Mean reference coverage (sensitivity proxy)",
       y = "Precision",
       title = "LOSO precision vs sensitivity, per superfamily",
       color = "Assembler") +
  theme_loso(legend_pos = "top")

p2

ggsave(file.path(opt$out_dir, "LOSO_precision_vs_sensitivity.png"),
       p2, width = 12, height = 6, dpi = 300)

# ---- Leakage check summary (if present) ----
leak_path <- file.path(opt$sim_dir, "leakage_check", "leakage_summary.tsv")
if (file.exists(leak_path)) {
  leak <- read_tsv(leak_path, show_col_types = FALSE)
  cat("\n=== Leakage check (contigs matching the TRAINING reference) ===\n")
  cat("Under LOSO, high-identity hits to the training reference are unexpected.\n")
  cat("A leak-free run should show low values in 'hits_pident90_qcov80'.\n\n")
  leak_summary <- leak %>%
    group_by(held_out_sf) %>%
    summarise(mean_high_hits = mean(hits_pident90_qcov80, na.rm = TRUE),
              max_high_hits  = max(hits_pident90_qcov80,  na.rm = TRUE))
  print(leak_summary, n = Inf)
  write_tsv(leak_summary, file.path(opt$out_dir, "LOSO_leakage_summary.tsv"))
}

cat(sprintf("\nLOSO aggregation complete. Outputs written to: %s\n", opt$out_dir))
