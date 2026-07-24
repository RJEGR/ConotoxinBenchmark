# ====================================================================
# LOSO_io.R
#
# I/O and path-handling helpers for the LOSO benchmark aggregator.
#
# Public API:
#   discover_contig_files(sim_dir)
#   load_manifest(loso_dir)
#   load_metadata(sim_dir)
#   load_leakage(sim_dir)        -- returns NULL if not present
#   count_test_fasta(loso_dir)   -- per-fold n_test reference sequences
#   join_loso_context(metrics_tbl, manifest, meta)
#
# ====================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(Biostrings)
})


# -----------------------------------------------------------------
# discover_contig_files(sim_dir)
#
# Find all TransRate per-assembly contigs.csv files under either
# <sim_dir>/transrate_contigs_dir/ or, as a fallback, the CWD.
# Returns a character vector of absolute paths; throws if none found.
# -----------------------------------------------------------------
discover_contig_files <- function(sim_dir) {

  cand_dirs <- c(
    file.path(sim_dir, "transrate_contigs_dir"),
    file.path(sim_dir),
    "transrate_contigs_dir"
  )
  for (d in cand_dirs) {
    if (!dir.exists(d)) next
    f <- list.files(d, pattern = "contigs\\.csv$",
                    recursive = TRUE, full.names = TRUE)
    if (length(f) > 0) {
      message(sprintf("[io] Found %d contigs.csv under %s", length(f), d))
      return(f)
    }
  }
  stop("No contigs.csv found under any candidate transrate dir for ", sim_dir)
}


# -----------------------------------------------------------------
# load_manifest(loso_dir)
#
# Read LOSO_manifest.tsv emitted by LOSO_conoServerDB.R. Normalises
# the held_out_superfamily column so spaces are replaced by underscores
# (matching downstream metadata keys).
# -----------------------------------------------------------------
load_manifest <- function(loso_dir) {

  f <- file.path(loso_dir, "LOSO_manifest.tsv")
  if (!file.exists(f)) stop("LOSO_manifest.tsv not found in ", loso_dir)
  m <- readr::read_tsv(f, show_col_types = FALSE)

  m %>%
    dplyr::mutate(
      held_out_superfamily = gsub("[[:space:]]", "_", held_out_superfamily)
    )
}


# -----------------------------------------------------------------
# load_metadata(sim_dir)
#
# Read LOSO_metadata.tsv from the simulation output dir. Provides the
# (sample_id, superfamily) mapping used to attribute each assembly to
# the correct held-out fold.
# -----------------------------------------------------------------
load_metadata <- function(sim_dir) {

  f <- file.path(sim_dir, "LOSO_metadata.tsv")
  if (!file.exists(f)) stop("LOSO_metadata.tsv not found in ", sim_dir)
  readr::read_tsv(f, show_col_types = FALSE) |> mutate(sample_id = gsub("_200x_PE", "", sample_id))
}


# -----------------------------------------------------------------
# load_leakage(sim_dir)
#
# Read leakage_summary.tsv (optional). The file is produced by Step 4
# of Run_LOSO_pipeline.sh; if absent, three_way_classify() cannot run
# and the caller should fall back to OOD-only metrics.
# -----------------------------------------------------------------
load_leakage <- function(sim_dir) {

  cand <- c(
    file.path(sim_dir, "leakage_check", "leakage_summary.tsv"),
    file.path(sim_dir, "transrate_contigs_dir",
              "leakage_check", "leakage_summary.tsv")
  )
  hit <- cand[file.exists(cand)]
  if (length(hit) == 0) {
    message("[io] No leakage_summary.tsv found; CIF and three-way ",
            "classification will be skipped.")
    return(NULL)
  }

  leak <- readr::read_tsv(hit[1], show_col_types = FALSE) |>
    dplyr::mutate(
      held_out_sf = gsub("[[:space:]]", "_", held_out_sf)
    )
  message(sprintf("[io] Loaded leakage summary: %d rows", nrow(leak)))
  return(leak)
}


# -----------------------------------------------------------------
# count_test_fasta(loso_dir)
#
# For each test_<sf>.fasta in the LOSO splits directory, count the
# number of sequences with width >= 200 bp (the length filter applied
# upstream in LOSO_conoServerDB.R). Used as InputNsequences -- the
# sensitivity denominator -- in calculate_metrics().
#
# Returns a tibble: sample_id, InputNsequences.
# The sample_id is built so it matches the keys parsed by read_contigs()
# (i.e. no "_200x_PE" suffix).
# -----------------------------------------------------------------
count_test_fasta <- function(loso_dir) {

  fastas <- list.files(loso_dir,
                       pattern = "^test_.*\\.fasta$",
                       full.names = TRUE)
  if (length(fastas) == 0)
    stop("No test_*.fasta files under ", loso_dir)

  rows <- lapply(fastas, function(f) {
    dna <- Biostrings::readDNAStringSet(f)
    n   <- sum(Biostrings::width(dna) >= 200)
    tibble::tibble(
      sample_id        = sub("\\.fasta$", "", basename(f)),
      InputNsequences  = n
    )
  })
  dplyr::bind_rows(rows)
}


# -----------------------------------------------------------------
# join_loso_context(metrics_tbl, manifest, meta)
#
# Attach (superfamily, n_test, n_train) context to a per-assembly
# metrics table, joining via the simulation metadata.
# -----------------------------------------------------------------
join_loso_context <- function(metrics_tbl, manifest, meta) {

  metrics_tbl %>%
    dplyr::left_join(meta, by = "sample_id") %>%
    dplyr::left_join(
      manifest %>%
        dplyr::select(held_out_superfamily, n_test, n_train),
      by = c("superfamily" = "held_out_superfamily")
    )
}
