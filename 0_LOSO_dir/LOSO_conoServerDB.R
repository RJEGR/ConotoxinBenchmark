# ====================================================================
# LOSO_conoServerDB.R
#
# Leave-One-Superfamily-Out (LOSO) cross-validation splitter.
# Strategy 1 of the circularity-controlled conotoxin benchmark.
#
# Adapted from rsampling_conoServerDB.R. Replaces random k-fold CV
# (which mixes superfamilies across folds and creates Grimm Type-2
# circularity) with grouped k-fold CV (rsample::group_vfold_cv) where
# v = n_distinct(genesuperfamily). Each fold holds out one entire
# superfamily.
#
# Inputs:
#   INPUTS/curated_nuc_conoServerDB.rds   (must contain: hits/entry_id,
#                                          sequence, genesuperfamily)
#
# Outputs (written to vfolds_loso_resampling_dir/):
#   - train_<superfamily>.fasta     : training reference per fold
#   - test_<superfamily>.fasta      : held-out ground truth per fold
#   - LOSO_manifest.tsv             : machine-readable fold registry
#   - LOSO_power_diagnostics.tsv    : per-fold sample-size diagnostics
#   - superfamily_distribution.tsv  : SF counts + fold eligibility flag
#
# Author: <author>
# ====================================================================

rm(list = ls())
if (!is.null(dev.list())) dev.off()

# ---- Libraries ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(rsample)
  library(Biostrings)
  library(Hmisc)
})

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

# ---- Configurable parameters ----
# Path to the curated ConoServer dataset (edit for your environment)
outdir <- "~/Documents/GitHub/ConotoxinBenchmark/INPUTS/"
# outdir <- "C:/Users/cinai/OneDrive/Documentos/GitHub/ConotoxinBenchmark/INPUTS"

MIN_SF_SIZE     <- 5      # superfamilies with fewer sequences are pooled into training
MIN_TRAIN_SIZE  <- 50     # if held-out set leaves training smaller than this, skip fold
SEED            <- 20260511
OR_TARGET       <- 2.0    # effect size for Hmisc::posamsize diagnostic
WRITE_WIDTH     <- 80     # FASTA line width

set.seed(SEED)

# ---- Setup ----
loso_dir <- file.path(outdir, "vfolds_loso_resampling_dir")
dir.create(loso_dir, recursive = TRUE, showWarnings = FALSE)

f <- file.path(outdir, "curated_nuc_conoServerDB.rds")
if (!file.exists(f)) stop("Cannot find input file: ", f)

df <- read_rds(f)

# Some prior scripts rename entry_id -> hits; accept either
if (!"hits" %in% names(df) && "entry_id" %in% names(df)) {
  df <- dplyr::rename(df, hits = entry_id)
}

required_cols <- c("hits", "sequence", "genesuperfamily")
missing_cols  <- setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# ---- Normalize superfamily labels ----
# Strip surrounding whitespace, collapse internal whitespace, replace NAs
df <- df %>%
  mutate(
    genesuperfamily = trimws(genesuperfamily),
    genesuperfamily = if_else(genesuperfamily == "" | is.na(genesuperfamily),
                              "UNCLASSIFIED", genesuperfamily)
  )

# ---- Diagnose superfamily distribution ----
sf_counts <- df %>%
  count(genesuperfamily, sort = TRUE, name = "n_sequences") %>%
  mutate(
    pct           = n_sequences / sum(n_sequences) * 100,
    fold_eligible = n_sequences >= MIN_SF_SIZE
  )

write_tsv(sf_counts, file.path(loso_dir, "superfamily_distribution.tsv"))

cat("\n=== Superfamily distribution ===\n")
print(sf_counts, n = Inf)

eligible_sfs <- sf_counts %>% filter(fold_eligible) %>% pull(genesuperfamily)

cat(sprintf("\nFold-eligible superfamilies (n >= %d): %d of %d total\n",
            MIN_SF_SIZE, length(eligible_sfs), nrow(sf_counts)))

# Sequences in too-small superfamilies are pooled into the training set of every
# fold (they are never held out as a test set, since LOSO is undefined for n<MIN).
df_eligible <- df %>% filter(genesuperfamily %in% eligible_sfs)
df_pooled   <- df %>% filter(!genesuperfamily %in% eligible_sfs)

cat(sprintf("Sequences in fold-eligible superfamilies: %d\n", nrow(df_eligible)))
cat(sprintf("Sequences pooled across all training sets: %d\n", nrow(df_pooled)))

if (length(eligible_sfs) < 2) {
  stop("Need at least 2 fold-eligible superfamilies for LOSO. Lower MIN_SF_SIZE or check input.")
}

# ---- Build LOSO folds with group_vfold_cv ----
# v = n_distinct(genesuperfamily): each fold's assessment() contains exactly one SF
n_folds <- length(eligible_sfs)

loso_folds <- group_vfold_cv(
  data  = df_eligible,
  group = genesuperfamily,
  v     = n_folds
)

cat(sprintf("\nCreated %d LOSO folds.\n", n_folds))

# ---- Helper: write a FASTA from a data frame ----
write_fasta <- function(df_in, path) {
  if (nrow(df_in) == 0) return(invisible(NULL))
  seqs <- DNAStringSet(df_in$sequence)
  # Headers: "<hits>|<superfamily>"
  names(seqs) <- paste(df_in$hits, df_in$genesuperfamily, sep = "|")
  writeXStringSet(seqs, filepath = path, format = "fasta", width = WRITE_WIDTH)
}

# ---- Per-fold export ----
manifest_rows <- vector("list", n_folds)
power_rows    <- vector("list", n_folds)

for (i in seq_len(n_folds)) {

  split_i <- loso_folds$splits[[i]]

  # rsample convention: analysis() = train, assessment() = test/held-out
  train_df <- analysis(split_i)
  test_df  <- assessment(split_i)

  held_out_sf <- unique(test_df$genesuperfamily)
  if (length(held_out_sf) != 1) {
    warning(sprintf("Fold %d: assessment set contains %d superfamilies (expected 1). Skipping.",
                    i, length(held_out_sf)))
    next
  }

  # Pool too-small superfamilies into the training side
  train_df <- bind_rows(train_df, df_pooled)

  if (nrow(train_df) < MIN_TRAIN_SIZE) {
    warning(sprintf("Fold %d (SF=%s): training set size %d < MIN_TRAIN_SIZE %d. Skipping.",
                    i, held_out_sf, nrow(train_df), MIN_TRAIN_SIZE))
    next
  }

  # Safe filename tag (strip non-alphanumeric)
  sf_safe   <- gsub("[^A-Za-z0-9]", "_", held_out_sf)
  fold_tag  <- sprintf("Fold%02d_%s", i, sf_safe)
  train_fa  <- file.path(loso_dir, sprintf("train_%s.fasta", sf_safe))
  test_fa   <- file.path(loso_dir, sprintf("test_%s.fasta",  sf_safe))

  write_fasta(train_df, train_fa)
  write_fasta(test_df,  test_fa)

  manifest_rows[[i]] <- tibble(
    fold_id              = i,
    fold_tag             = fold_tag,
    held_out_superfamily = held_out_sf,
    n_test               = nrow(test_df),
    n_train              = nrow(train_df),
    test_fasta           = normalizePath(test_fa,  mustWork = FALSE),
    train_fasta          = normalizePath(train_fa, mustWork = FALSE)
  )

  # Per-fold sample-size diagnostic
  power_rows[[i]] <- tibble(
    fold_id              = i,
    held_out_superfamily = held_out_sf,
    n_test               = nrow(test_df),
    n_train              = nrow(train_df),
    log10_n_test         = log10(nrow(test_df)),
    p_test               = nrow(test_df) / (nrow(test_df) + nrow(train_df))
  )

  cat(sprintf("Fold %2d  SF=%-25s  n_test=%4d  n_train=%5d  -> %s\n",
              i, held_out_sf, nrow(test_df), nrow(train_df), basename(test_fa)))
}

manifest_df <- bind_rows(manifest_rows)
power_df    <- bind_rows(power_rows)

# ---- Aggregate effect-size diagnostic (Hmisc::posamsize) ----
# Detectable OR using the marginal SF distribution as the ordinal reference.
ref_prop <- prop.table(table(df_eligible$genesuperfamily))

posamsize_result <- tryCatch({
  posamsize(p = as.numeric(ref_prop), odds.ratio = OR_TARGET,
            fraction = 0.5, alpha = 0.05, power = 0.9)
}, error = function(e) {
  message("posamsize failed: ", conditionMessage(e)); NULL
})

if (!is.null(posamsize_result)) {
  cat(sprintf("\nposamsize (OR=%.1f, power=0.9, alpha=0.05): required N = %s\n",
              OR_TARGET, round(posamsize_result$n)))
}

# ---- Write outputs ----
write_tsv(manifest_df, file.path(loso_dir, "LOSO_manifest.tsv"))
write_tsv(power_df,    file.path(loso_dir, "LOSO_power_diagnostics.tsv"))

cat(sprintf("\nWrote %d LOSO folds to %s\n", nrow(manifest_df), loso_dir))
cat("Done.\n")
