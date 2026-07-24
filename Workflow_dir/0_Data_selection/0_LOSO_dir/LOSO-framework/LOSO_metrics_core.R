# ====================================================================
# LOSO_metrics_core.R
#
# Pure metric functions for the LOSO benchmark. No I/O side effects;
# no globals; safe to source() from any context.
#
# Implements the §3.2 metric (F1 at 90% reference coverage, full-length
# best-reciprocal-BLAST per Behera et al. 2022) and the three-way
# Train+Test classification (TP / FP_paralog / FP_spurious) plus the
# Circularity Inflation Factor (CIF) per Bernett et al. (2024).
#
# Public API:
#   read_contigs(csv_path)
#   calculate_metrics(df, reference_coverage_val, input_n_seq_tbl)
#   three_way_classify(metrics_tbl, leak_tbl)
#   compute_cif(three_way_tbl)
#
# Author: <author>
# ====================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(tibble)
  library(stringr)
})


# -----------------------------------------------------------------
# read_contigs(csv_path)
#
# Load a single TransRate contigs.csv, tag it with (sample_id, tool)
# parsed from the parent directory name, and pre-condition the columns
# downstream functions expect.
#
# Returns:  grouped tibble (sample_id, tool) with at minimum
#           - contig_name
#           - length
#           - reference_coverage
#           - hits          (NA when no BLAST hit)
#           - CRBB / has_crb (boolean) if present
#
# Returns NULL on unreadable / empty inputs (callers must check).
# -----------------------------------------------------------------
read_contigs <- function(csv_path) {

  df <- tryCatch(
    suppressWarnings(readr::read_csv(csv_path, show_col_types = FALSE)),
    error = function(e) NULL
  )
  if (is.null(df) || nrow(df) == 0) return(NULL)

  # Parse "<sample_id>_<tool>" from the assembly directory name.
  # Strip trailing "_dir" suffix written by Assemblers.sh.
  asm_dir <- basename(dirname(csv_path))
  asm_bs  <- sub("_dir$", "", asm_dir)
  parts <- strsplit(asm_bs, "_")[[1]]
  if (length(parts) < 2) return(NULL)
  tool      <- tail(parts, 1)
  sample_id <- paste(head(parts, -1), collapse = "_")
  # Drop the read-depth/PE tag to match LOSO_metadata.tsv keys
  sample_id <- gsub("_200x_PE|_100x_PE|_50x_PE", "", sample_id)

  # Ensure required columns exist; coerce missing ones to NA so the
  # downstream calculate_metrics() never breaks on an exotic TransRate
  # version mismatch.
  for (col in c("reference_coverage", "hits", "length")) {
    if (!col %in% colnames(df)) df[[col]] <- NA
  }

  df <- df %>%
    dplyr::mutate(
      tool      = tool,
      sample_id = sample_id
    ) %>%
    dplyr::group_by(sample_id, tool)

  return(df)
}


# -----------------------------------------------------------------
# calculate_metrics(df, reference_coverage_val, input_n_seq_tbl)
#
# Compute TP / FP / TN / FN counts at one or more coverage thresholds.
# Mirrors the Behera et al. (2022) categorization used in manuscript
# section 3.2; default threshold = 0.90 (full-length 90% identity RBH).
#
# Inputs:
#   df                    grouped tibble from read_contigs()
#   reference_coverage_val numeric scalar or vector of thresholds in [0,1]
#   input_n_seq_tbl       tibble with cols (sample_id, InputNsequences) -
#                         the number of test-set reference sequences
#                         (used as the denominator of sensitivity).
#
# Returns a tibble with one row per (sample_id, tool, threshold) and
# columns: rawcontigs, TP, FP, TN, FN, sensitivity, precision, F1.
# -----------------------------------------------------------------
calculate_metrics <- function(df,
                              reference_coverage_val = 0.90,
                              input_n_seq_tbl) {

  stopifnot(
    !is.null(df),
    all(c("reference_coverage", "hits") %in% colnames(df)),
    is.numeric(reference_coverage_val),
    is.data.frame(input_n_seq_tbl),
    all(c("sample_id", "InputNsequences") %in% colnames(input_n_seq_tbl))
  )

  group_cols <- dplyr::group_vars(df)
  if (length(group_cols) == 0) group_cols <- c("sample_id", "tool")

  IS_CHIMERIC <- 0   # reference_coverage == 0 means "no BLAST evidence"

  # Behera et al. (2022) -- as used in manuscript section 2.2 -- defines
  # FP as "incorrectly assembled, including either partially correctly
  # assembled or those with no similarity with the reference". So FP
  # collapses 'partial-hit' and 'no-hit' contigs into one bucket. We keep
  # them as separate auxiliary counters (FP_partial, FP_nohit) for
  # diagnostic purposes, but the headline FP is their sum.
  #
  # Units (important): TP_contigs and FP are CONTIG counts; TP_refs and
  # FN are REFERENCE-SEQUENCE counts. Precision uses contig units;
  # sensitivity uses reference units, in line with Hoelzer & Marz (2019)
  # and Behera et al. (2022).
  result <- df %>%
    dplyr::ungroup() %>%
    tidyr::crossing(coverage_threshold = reference_coverage_val) %>%
    dplyr::mutate(
      ref_cov_clean = ifelse(is.na(reference_coverage), 0, reference_coverage)
    ) %>%
    dplyr::group_by(across(all_of(c(group_cols, "coverage_threshold")))) %>%
    dplyr::summarise(
      rawcontigs = dplyr::n(),
      # Contig-level TP: contigs hitting any reference at the threshold.
      TP         = sum(ref_cov_clean >= coverage_threshold, na.rm = TRUE),
      # Diagnostic split of FP into 'partial hit' and 'no hit at all'.
      FP_partial = sum(ref_cov_clean > IS_CHIMERIC &
                       ref_cov_clean < coverage_threshold, na.rm = TRUE),
      FP_nohit   = sum(ref_cov_clean <= IS_CHIMERIC, na.rm = TRUE),
      # Reference-level TP: distinct refs captured by at least one TP contig.
      TP_refs    = dplyr::n_distinct(hits[ref_cov_clean >= coverage_threshold]),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      FP = FP_partial + FP_nohit                # Behera FP collapsed
    ) %>%
    dplyr::left_join(input_n_seq_tbl, by = "sample_id") %>%
    dplyr::mutate(
      # FN = reference seqs not detected by any high-quality contig.
      FN = pmax(0L, InputNsequences - TP_refs),
      precision   = ifelse((TP + FP) > 0, TP / (TP + FP), NA_real_),
      sensitivity = ifelse(InputNsequences > 0, TP_refs / InputNsequences, NA_real_),
      F1 = ifelse((precision + sensitivity) > 0,
                  2 * precision * sensitivity / (precision + sensitivity),
                  NA_real_)
    ) %>%
    dplyr::rename(reference_coverage_val = coverage_threshold) %>%
    dplyr::select(
      all_of(group_cols), reference_coverage_val,
      rawcontigs, TP, FP, FP_partial, FP_nohit, TP_refs, FN,
      precision, sensitivity, F1
    )

  return(result)
}


# -----------------------------------------------------------------
# three_way_classify(metrics_tbl, leak_tbl)
#
# Augment per-assembly metrics with the three-way TP / FP_paralog /
# FP_spurious decomposition by intersecting test-reference hits with
# the BLAST-against-training leakage audit.
#
# Inputs:
#   metrics_tbl   output of calculate_metrics() (one row per assembly
#                 x threshold)
#   leak_tbl      from leakage_summary.tsv, with columns
#                   assembly, held_out_sf, total_hits, hits_pident90_qcov80
#
# Returns metrics_tbl with these added columns:
#   FP_paralog   contigs that hit train>=90/80 but did NOT hit test at
#                the current threshold. Conservatively capped at the
#                total non-TP contig count to avoid double-counting.
#   FP_spurious  contigs hitting neither test nor train (likely chimera
#                or non-conotoxin).
#   precision_ood        TP / rawcontigs                (current LOSO metric)
#   precision_combined   (TP + FP_paralog) / rawcontigs (when training-SF
#                        paralog hits are credited as biologically real)
# -----------------------------------------------------------------
three_way_classify <- function(metrics_tbl, leak_tbl) {

  stopifnot(
    is.data.frame(metrics_tbl),
    is.data.frame(leak_tbl),
    all(c("assembly", "hits_pident90_qcov80") %in% colnames(leak_tbl))
  )

  # The leakage table keys assemblies by their full assembly basename
  # (e.g. "test_A_superfamily_200x_PE_IDBA"). We reduce it to a
  # (sample_id, tool) key that matches metrics_tbl.
  leak_key <- leak_tbl %>%
    dplyr::mutate(
      asm_bs    = sub("_dir$", "", assembly),
      tool      = stringr::str_extract(asm_bs, "[^_]+$"),
      sample_id = sub("_[^_]+$", "", asm_bs),
      sample_id = gsub("_200x_PE|_100x_PE|_50x_PE", "", sample_id)
    ) %>%
    dplyr::select(held_out_sf, sample_id, tool, n_leak_high = hits_pident90_qcov80) %>%
    # dplyr::group_by(sample_id, tool) %>%
    dplyr::group_by(sample_id, held_out_sf) %>%
    dplyr::summarise(n_leak_high = sum(n_leak_high, na.rm = TRUE),
                     .groups = "drop")

  out <- metrics_tbl %>%
    # dplyr::left_join(leak_key, by = c("sample_id", "tool")) %>%
    dplyr::left_join(leak_key, by = c("sample_id")) %>%
    dplyr::mutate(
      n_leak_high = tidyr::replace_na(n_leak_high, 0L),
      # contigs that are NOT TP and could potentially be paralog hits
      n_non_tp = pmax(0L, rawcontigs - TP),
      # Conservatively cap paralog hits at the non-TP pool so we never
      # double-count a contig that's already TP against test.
      FP_paralog  = pmin(n_non_tp, n_leak_high),
      FP_spurious = pmax(0L, rawcontigs - TP - FP_paralog),
      # Two precisions: pure OOD (current LOSO definition) and combined
      # (counts paralog-grade hits to training-SFs as biologically real).
      precision_ood       = ifelse(rawcontigs > 0, TP / rawcontigs, NA_real_),
      precision_combined  = ifelse(rawcontigs > 0,
                                   (TP + FP_paralog) / rawcontigs,
                                   NA_real_)
    ) %>%
    dplyr::select(-n_non_tp)

  return(out)
}


# -----------------------------------------------------------------
# compute_cif(three_way_tbl)
#
# Circularity Inflation Factor (CIF) per (tool, threshold), summarised
# over held-out superfamilies. CIF = precision_combined / precision_ood
# at the assembly level, then averaged across SFs within a tool.
#
# CIF == 1   --> random-CV would have given the same precision as LOSO
# CIF >> 1   --> random-CV was inflated by within-family homology
# -----------------------------------------------------------------
compute_cif <- function(three_way_tbl) {

  stopifnot(
    all(c("tool", "reference_coverage_val",
          "precision_ood", "precision_combined") %in% colnames(three_way_tbl))
  )

  per_assembly <- three_way_tbl %>%
    dplyr::mutate(
      CIF_assembly = ifelse(precision_ood > 0,
                            precision_combined / precision_ood,
                            NA_real_)
    )

  per_tool <- per_assembly %>%
    # dplyr::group_by(tool, reference_coverage_val) %>%
    dplyr::group_by(reference_coverage_val, held_out_sf) %>%
    dplyr::summarise(
      n_assemblies   = dplyr::n(),
      median_CIF     = median(CIF_assembly, na.rm = TRUE),
      mean_CIF       = mean(CIF_assembly,   na.rm = TRUE),
      sd_CIF         = sd(CIF_assembly,     na.rm = TRUE),
      mean_p_ood     = mean(precision_ood,      na.rm = TRUE),
      mean_p_combined= mean(precision_combined, na.rm = TRUE),
      .groups = "drop"
    )

  list(per_assembly = per_assembly, per_tool = per_tool)
}

