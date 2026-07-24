# =============================================================================
# transrate_EDA.R
# Unified Exploratory Data Analysis for transrate contigs.csv outputs.
#
# Replaces the parallel logic in Kmer.R + Subsampling.R with a single script
# that:
#   1. Walks a transrate_contigs_dir tree (or a parent of several)
#   2. Auto-detects the sweep type from directory names (kmer | sampling | both)
#   3. Computes TP / FP / FN and the four assembly metrics:
#         Accuracy = TP / (TP + FN + FP)
#         Precision = TP / (TP + FP)
#         Sensitivity = TP / (TP + FN)
#         F-score = 2*TP / (2*TP + FP + FN)
#   4. Writes a tidy TSV and produces benchmark plots.
#
# It accepts two directory conventions:
#   (a) New unified naming from benchmark.sh:
#        <factor>_<assembler>_k<KMER>_p<PROP>_dir
#   (b) Legacy naming from kmer.sh / subsampling.sh:
#        <factor>_..._<value>_dir
#
# Usage (interactive or Rscript):
#   Rscript transrate_EDA.R \
#       --input  /path/to/benchmark/root \
#       --refdir /path/to/vfolds_resampling_dir \
#       --out    /path/to/output_dir \
#       --rcov   0.95
# =============================================================================

rm(list = ls())

if(!is.null(dev.list())) dev.off()


dir <- "~/Documents/GitHub/ConotoxinBenchmark/3_kmer_dir/transrate_contigs_dir/"

setwd(dir)

suppressPackageStartupMessages({
  library(tidyverse)
  library(optparse)
})

# ----- CLI -------------------------------------------------------------------
option_list <- list(
  make_option(c("-i", "--input"),  type = "character", default = ".",
              help = "Root directory containing transrate_contigs_dir/ trees"),
  make_option(c("-r", "--refdir"), type = "character", default = NULL,
              help = "Directory with reference FASTAs (e.g. Fold01.fasta).
                      Used to count InputNsequences (length >= 200) for FN."),
  make_option(c("-o", "--out"),    type = "character", default = "EDA_out",
              help = "Output directory for tables and plots [default %default]"),
  make_option(c("-c", "--rcov"),   type = "double", default = 0.9,
              help = "reference_coverage cutoff for TP [default %default]"),
  make_option(c("-l", "--minlen"), type = "integer", default = 200,
              help = "Minimum contig length to keep [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$out, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 1. Discover and parse contigs.csv files
# =============================================================================

# Each contigs.csv lives at:
#   <input>/.../transrate_contigs_dir/<RUN_DIR>/contigs.csv
# RUN_DIR encodes the run parameters; we parse it with two regexes (new + legacy).

parse_run_dir <- function(run_dir) {
  # Strip "_dir" suffix if present
  base <- sub("_dir$", "", run_dir)
  
  # ---- (a) new unified scheme:  <factor>_<assembler>_k<KMER>_p<PROP>
  # The factor may itself contain underscores (e.g. "Fold01_200x_PE_samples"),
  # so we use a greedy `.+` for it; the `_k<val>_p<val>$` tail anchors the match.
  m <- str_match(base,
                 "^(.+)_([A-Za-z]+)_k([^_]+)_p([^_]+)$")
  if (!is.na(m[1, 1])) {
    return(tibble(
      vfold_set    = m[1, 2],
      Assembly     = m[1, 3],
      kmer         = ifelse(m[1, 4] == "NA", NA_character_, m[1, 4]),
      sampling_set = ifelse(m[1, 5] == "NA", NA_real_, as.numeric(m[1, 5]))
    ))
  }
  
  # ---- (b) legacy scheme: split on '_'; position 1 is factor, position 5 is value
  parts <- strsplit(base, "_")[[1]]
  factor_id <- parts[1]
  val <- if (length(parts) >= 5) parts[5] else NA
  is_numeric_val <- suppressWarnings(!is.na(as.numeric(val)))
  
  if (is_numeric_val && as.numeric(val) <= 1) {
    # looks like a subsampling proportion (0,1]
    tibble(vfold_set = factor_id, Assembly = NA_character_,
           kmer = NA_character_, sampling_set = as.numeric(val))
  } else {
    # treat as kmer
    tibble(vfold_set = factor_id, Assembly = NA_character_,
           kmer = as.character(val), sampling_set = NA_real_)
  }
}


read_transrate_scores <- function(file_path, root) {
  rel <- sub(paste0("^", normalizePath(root, mustWork = FALSE), "/"), "",
             normalizePath(file_path, mustWork = FALSE))
  # The directory immediately containing contigs.csv = the run dir
  run_dir <- basename(dirname(file_path))
  # The directory ABOVE transrate_contigs_dir often carries the assembler
  above <- basename(dirname(dirname(dirname(file_path))))

  meta <- parse_run_dir(run_dir)
  
  # If legacy naming didn't yield an Assembly, infer it from the parent dir
  if (is.na(meta$Assembly) && length(above) == 1 && nchar(above) > 0) {
    meta$Assembly <- sub("_dir$", "", above)
  }

  read_csv(file_path, show_col_types = FALSE) %>%
    bind_cols(meta[rep(1, nrow(.)), ]) %>%
    mutate(rel_path = rel, run_dir = run_dir)
}

message("Scanning ", opt$input, " for contigs.csv ...")

file_list <- list.files(opt$input, pattern = "^contigs\\.csv$",
                        recursive = TRUE, full.names = TRUE)

if (length(file_list) == 0) stop("No contigs.csv found under ", opt$input)
message("Found ", length(file_list), " contigs.csv files")


transratedf <- map_dfr(file_list, read_transrate_scores, root = opt$input)

message("Combined rows: ", nrow(transratedf))

# =============================================================================
# 2. Compute InputNsequences (denominator for FN)
# =============================================================================

count_reference_sequences <- function(refdir, minlen = 200) {
  if (is.null(refdir) || !dir.exists(refdir)) {
    warning("--refdir not provided or missing; FN will be set to 0.")
    return(NULL)
  }
  fastas <- list.files(refdir, pattern = "\\.fa(sta)?$", full.names = TRUE)
  if (!length(fastas)) {
    warning("No FASTA files in ", refdir)
    return(NULL)
  }
  # Use Biostrings if available, else a tiny built-in reader
  read_fa_widths <- function(f) {
    if (requireNamespace("Biostrings", quietly = TRUE)) {
      Biostrings::width(Biostrings::readDNAStringSet(f))
    } else {
      lines <- readLines(f)
      idx <- grep("^>", lines)
      starts <- idx + 1; ends <- c(idx[-1] - 1, length(lines))
      vapply(seq_along(idx),
             function(i) sum(nchar(lines[starts[i]:ends[i]])), integer(1))
    }
  }
  out <- vapply(fastas, function(f) sum(read_fa_widths(f) >= minlen), integer(1))
  names(out) <- sub("\\.fa(sta)?$", "", basename(fastas))
  out
}

opt$refdir <- outdir <- "~/Documents/GitHub/ConotoxinBenchmark/INPUTS/vfolds_resampling_dir/"

InputNsequences <- count_reference_sequences(opt$refdir, opt$minlen)

# =============================================================================
# 3. Metric computation (same definitions as Kmer.R / Subsampling.R)
# =============================================================================

calculate_metrics <- function(df, reference_coverage_val = 0.95,
                              InputN = NULL) {
  is_chimeric <- 0

  rawcontigs <- df %>% summarise(rawcontigs = n())

  TP <- df %>%
    filter(reference_coverage >= reference_coverage_val) %>%
    distinct(hits) %>%
    summarise(TP = n())

  FP <- df %>%
    mutate(reference_coverage = ifelse(is.na(hits) & is.na(reference_coverage),
                                       0, reference_coverage)) %>%
    filter(reference_coverage > is_chimeric &
           reference_coverage < reference_coverage_val) %>%
    summarise(FP = n())

  bind_cols(rawcontigs, TP, FP) %>%
    mutate(across(everything(), ~replace_na(., 0))) %>%
    mutate(FN = if (!is.null(InputN)) {
                  max(0, InputN - TP)
                } else 0) %>%
    mutate(
      Ratio       = ifelse(FP == 0, NA_real_, TP / FP),
      Accuracy    = TP / (TP + FN + FP),
      Precision   = ifelse((TP + FP) == 0, NA_real_, TP / (TP + FP)),
      Sensitivity = ifelse((TP + FN) == 0, NA_real_, TP / (TP + FN)),
      Fscore      = ifelse((2 * TP + FP + FN) == 0, NA_real_,
                           2 * TP / (2 * TP + FP + FN))
    )
}

# Detect which axes vary so we can group correctly
has_kmer    <- any(!is.na(transratedf$kmer))
has_sampset <- any(!is.na(transratedf$sampling_set))

group_vars <- c("vfold_set", "Assembly")

if (has_kmer)    group_vars <- c(group_vars, "kmer")
if (has_sampset) group_vars <- c(group_vars, "sampling_set")

message("Grouping by: ", paste(group_vars, collapse = ", "))



metricsdf <- transratedf %>%
  mutate(vfold_set =  sub("_200x_PE_samples$", "", vfold_set)) |>
  filter(length >= opt$minlen) %>%
  group_split(across(all_of(group_vars))) %>%
  map_dfr(function(g) {
    InputN <- if (!is.null(InputNsequences)) {
                InputNsequences[unique(g$vfold_set)]
              } else NULL
    key <- g %>% select(all_of(group_vars)) %>% distinct()
    bind_cols(key, calculate_metrics(g, opt$rcov, InputN))
  })

out_tsv <- file.path(opt$out, "benchmark_metrics.tsv")

write_tsv(metricsdf, out_tsv)

message("Wrote ", out_tsv)


