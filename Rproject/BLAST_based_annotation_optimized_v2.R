# ============================================================================
# BLAST Annotation for Conotoxins - Optimized & Corrected (v2)
# Based on Cahais et al., 2012 with modifications for multidomain toxins
# ============================================================================
#
# CHANGELOG (v2 fixes):
#   1. CRITICAL: Replaced case_when() with if/else if/else in
#      annotate_blast_hits(). R's case_when evaluates ALL right-hand sides
#      regardless of which condition matches, causing spurious errors when
#      overlap functions are called for non-matching categories.
#   2. Removed IRanges dependency. overlap_ratio() now uses base R arithmetic.
#      IRanges::pintersect(resolve.empty = "none") throws when intervals are
#      disjoint, which was the proximal cause of "Failed to evaluate the
#      right-hand side of formula N" errors.
#   3. Added input validation and defensive guards throughout the overlap
#      checking functions.
#   4. Expanded annotate_blast_hits() output to include n_hits,
#      n_unique_subjects, and max_coverage for downstream diagnostics.
#   5. Added verbose logging controlled by a single flag.
# ============================================================================

# library(dplyr)
# library(tidyr)
# library(ggplot2)
# library(ggsci)
# library(furrr)
# library(readr)

rm(list = ls())

gc()

Sys.setenv(R_MAX_VSIZE = "100Gb")

install_load <- function(package) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package)
    library(package, character.only = TRUE)
  }
}

packages <- c("dplyr", "tidyr", "ggplot2", "ggsci", "furrr", "readr")

sapply(packages, install_load)

if (!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

# ============================================================================
# 1. OVERLAP RATIO CALCULATION (base-R, no IRanges dependency)
# ============================================================================

#' Calculate the overlap ratio between two genomic intervals
#'
#' The ratio is defined as the intersection width divided by the width of the
#' shorter interval (Cahais et al., 2012).  This version uses only base-R
#' arithmetic, avoiding the IRanges pitfall where pintersect(...,
#' resolve.empty = "none") throws on disjoint ranges.
#'
#' @param start1,end1 Numeric.  Boundaries of the first interval (order does
#'   not matter; the function takes pmin/pmax internally).
#' @param start2,end2 Numeric.  Boundaries of the second interval.
#' @return Numeric in [0, 1].  Returns 0 on any NA / non-positive width input.

overlap_ratio <- function(start1, end1, start2, end2) {
  # --- guard: NAs --------------------------------------------------------
  if (any(is.na(c(start1, end1, start2, end2)))) return(0)

  # --- normalise so that start <= end ------------------------------------
  s1 <- pmin(start1, end1)
  e1 <- pmax(start1, end1)
  s2 <- pmin(start2, end2)
  e2 <- pmax(start2, end2)

  # --- interval widths (1-based, inclusive) -------------------------------
  w1 <- e1 - s1 + 1L
  w2 <- e2 - s2 + 1L

  if (any(c(w1, w2) <= 0)) return(0)

  # --- intersection width ------------------------------------------------
  ovl_start <- pmax(s1, s2)
  ovl_end   <- pmin(e1, e2)
  ovl_width <- pmax(0L, ovl_end - ovl_start + 1L)

  # --- ratio vs shorter interval ----------------------------------------
  shortest <- pmin(w1, w2)
  result   <- ovl_width / shortest
  result[is.na(result) | is.infinite(result)] <- 0

  return(result)
}

# ============================================================================
# 2. READ AND FILTER BLAST RESULTS
# ============================================================================

#' Read a BLAST outfmt-6 file (14-column flavour)
read_blast_outfmt6 <- function(f) {
  outfmt6_cols <- c(
    "qseqid", "sseqid", "pident", "length", "mismatches", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen"
  )

  readr::read_tsv(f, col_names = FALSE, show_col_types = FALSE) %>%
    setNames(outfmt6_cols) %>%
    mutate(db = basename(f))
}

#' Filter BLAST results on coverage and identity
filter_blast_results <- function(blast_df, min_coverage = 0.5, min_identity = 80) {
  blast_df %>%
    filter(!is.na(qlen), !is.na(slen), qlen > 0, slen > 0,
           !is.na(length), length > 0) %>%
    mutate(
      query_coverage   = length / qlen,
      subject_coverage = length / slen,
      max_coverage     = pmax(query_coverage, subject_coverage)
    ) %>%
    filter(max_coverage >= min_coverage, pident >= min_identity)
}

# ============================================================================
# 3. PRELIMINARY CATEGORY ASSIGNMENT
# ============================================================================

#' Assign each query to one of: no hit, 1->1, m->1, 1->n, m->n
assign_preliminary_category <- function(blast_df) {
  if (nrow(blast_df) == 0) {
    return(tibble(qseqid = character(), prelim_cat = character(),
                  n_hits = integer(), n_unique_subjects = integer(),
                  max_contigs_per_subject = integer()))
  }

  query_stats <- blast_df %>%
    group_by(qseqid) %>%
    summarise(
      n_hits            = n(),
      n_unique_subjects = n_distinct(sseqid),
      .groups = "drop"
    )

  subject_stats <- blast_df %>%
    group_by(sseqid) %>%
    summarise(n_contigs = n_distinct(qseqid), .groups = "drop")

  query_stats %>%
    left_join(
      blast_df %>%
        select(qseqid, sseqid) %>%
        distinct() %>%
        left_join(subject_stats, by = "sseqid") %>%
        group_by(qseqid) %>%
        summarise(max_contigs_per_subject = max(n_contigs, na.rm = TRUE),
                  .groups = "drop"),
      by = "qseqid"
    ) %>%
    mutate(
      prelim_cat = case_when(
        n_hits == 0                                  ~ "no hit",
        n_hits == 1  & max_contigs_per_subject == 1  ~ "1->1",
        n_hits >  1  & max_contigs_per_subject == 1  ~ "m->1",
        n_hits == 1  & max_contigs_per_subject >  1  ~ "1->n",
        TRUE                                         ~ "m->n"
      )
    )
}

# ============================================================================
# 4. CONOTOXIN-SPECIFIC ANNOTATION LOGIC
# ============================================================================

# For conotoxins, overlapping hits on subject sequences indicate MULTIDOMAIN
# (NOT chimera) because conotoxins often have multiple functional domains.

#' Check whether any pair of subject-coordinate intervals overlap
#'
#' Used for the **1->n** case: a single query hits multiple subjects.
#' We ask whether different subject alignments overlap on the *subject*
#' coordinate axis.
#'
#' @param hits  data.frame  Filtered BLAST rows for one query.
#' @param threshold  numeric  Minimum overlap ratio to call "overlapping".
#' @return logical

check_hit_overlaps_on_query <- function(hits, threshold = 0.5) {
  # --- guard ---
  if (is.null(hits) || nrow(hits) < 2) return(FALSE)

  hit_intervals <- hits %>%
    mutate(start = pmin(sstart, send),
           end   = pmax(sstart, send)) %>%
    filter(!is.na(start), !is.na(end), start <= end)

  n <- nrow(hit_intervals)
  if (n < 2) return(FALSE)

  for (x in seq_len(n - 1)) {
    for (y in (x + 1):n) {
      ov <- overlap_ratio(
        hit_intervals$start[x], hit_intervals$end[x],
        hit_intervals$start[y], hit_intervals$end[y]
      )
      if (!is.na(ov) && ov > threshold) return(TRUE)
    }
  }

  return(FALSE)
}

#' Check whether different contigs that map to the same subject overlap
#'
#' Used for the **m->1** case: multiple contigs hit a single subject.
#' We ask whether those contigs overlap on the *subject* coordinate axis.
#'
#' @param blast_df   data.frame  The full (filtered) BLAST table.
#' @param subject_id character   The subject sequence to inspect.
#' @param query_id   character   (Informational; not used in logic.)
#' @param threshold  numeric     Minimum overlap ratio to call "overlapping".
#' @return logical

check_contig_overlaps_on_subject <- function(blast_df, subject_id, query_id,
                                              threshold = 0.5) {
  # --- guard ---
  if (is.null(subject_id) || is.na(subject_id) || subject_id == "") {
    return(FALSE)
  }

  contigs_on_subject <- blast_df %>%
    filter(sseqid == subject_id) %>%
    select(qseqid, sstart, send) %>%
    distinct() %>%
    filter(!is.na(sstart), !is.na(send)) %>%
    mutate(start = pmin(sstart, send),
           end   = pmax(sstart, send)) %>%
    filter(start <= end)

  contig_ids <- unique(contigs_on_subject$qseqid)
  if (length(contig_ids) < 2) return(FALSE)

  n_ids <- length(contig_ids)
  for (a in seq_len(n_ids - 1)) {
    for (b in (a + 1):n_ids) {
      i1 <- contigs_on_subject %>% filter(qseqid == contig_ids[a])
      i2 <- contigs_on_subject %>% filter(qseqid == contig_ids[b])

      if (nrow(i1) == 0 || nrow(i2) == 0) next

      ov <- overlap_ratio(i1$start[1], i1$end[1],
                          i2$start[1], i2$end[1])
      if (!is.na(ov) && ov > threshold) return(TRUE)
    }
  }

  return(FALSE)
}

# ============================================================================
# 4b. CORE ANNOTATION ENGINE  (FIXED: if/else replaces case_when)
# ============================================================================

#' Classify one query sequence
#'
#' This is the per-query annotation logic extracted into its own function so
#' that it can be called safely inside a tryCatch without the case_when
#' evaluate-everything pitfall.
#'
#' @param cat_type          character  Preliminary category (1->1 etc.)
#' @param hits              data.frame Filtered BLAST rows for this query.
#' @param blast_df          data.frame The full filtered BLAST table.
#' @param q_id              character  Query sequence id.
#' @param overlap_threshold numeric    Overlap ratio threshold (default 0.5).
#' @param coverage_threshold numeric   Coverage threshold for "full" (default 0.9).
#' @return character  One of: full, fragment, allele, multi, chimera, other, no hit.

classify_query <- function(cat_type, hits, blast_df, q_id,
                            overlap_threshold  = 0.5,
                            coverage_threshold = 0.9) {

  # ---- no hit ----
  if (cat_type == "no hit") {
    return("no hit")
  }

  # ---- 1 -> 1 ----
  if (cat_type == "1->1") {
    coverage <- hits$length[1] / hits$slen[1]
    return(if (coverage >= coverage_threshold) "full" else "fragment")
  }

  # ---- m -> 1 ----
  if (cat_type == "m->1") {
    subject_id <- hits$sseqid[1]
    overlaps   <- check_contig_overlaps_on_subject(
                    blast_df, subject_id, q_id, overlap_threshold)
    return(if (overlaps) "allele" else "fragment")
  }

  # ---- 1 -> n  (conotoxin-specific) ----
  if (cat_type == "1->n") {
    overlaps <- check_hit_overlaps_on_query(hits, overlap_threshold)
    return(if (overlaps) "multi" else "chimera")
  }

  # ---- m -> n ----
  if (cat_type == "m->n") {
    n_contigs  <- n_distinct(hits$qseqid)
    n_subjects <- n_distinct(hits$sseqid)

    if (n_contigs == n_subjects) {
      coverage_vec <- hits %>%
        group_by(qseqid) %>%
        slice_max(length, n = 1, with_ties = FALSE) %>%
        ungroup() %>%
        mutate(coverage = length / slen) %>%
        pull(coverage)

      return(
        if (mean(coverage_vec >= coverage_threshold) >= 0.5) "full" else "fragment"
      )
    } else {
      return("multi")
    }
  }

  # ---- fallback ----
  return("other")
}

# ============================================================================
# 4c. ANNOTATE ALL QUERIES IN A BLAST TABLE
# ============================================================================

#' Annotate every query in a (filtered) BLAST table
#'
#' @param blast_df           data.frame  Raw BLAST outfmt6 table (will be filtered).
#' @param overlap_threshold  numeric     Default 0.5.
#' @param coverage_threshold numeric     Default 0.9.
#' @param verbose            logical     Print per-query progress.
#' @return tibble with columns: qseqid, prelim_cat, final_annotation,
#'         n_hits, n_unique_subjects, max_contigs_per_subject,
#'         best_bitscore, best_pident, best_max_coverage.

annotate_blast_hits <- function(blast_df,
                                 overlap_threshold  = 0.5,
                                 coverage_threshold = 0.9,
                                 verbose = FALSE) {

  blast_df <- filter_blast_results(blast_df)

  if (nrow(blast_df) == 0) {
    return(tibble(
      qseqid              = character(),
      prelim_cat           = character(),
      final_annotation     = character(),
      n_hits               = integer(),
      n_unique_subjects    = integer(),
      max_contigs_per_subject = integer(),
      best_bitscore        = numeric(),
      best_pident          = numeric(),
      best_max_coverage    = numeric()
    ))
  }

  prelim_cats <- assign_preliminary_category(blast_df)

  # Build annotations table with diagnostics
  annotations <- tibble(qseqid = unique(blast_df$qseqid)) %>%
    left_join(prelim_cats, by = "qseqid") %>%
    mutate(
      prelim_cat              = replace_na(prelim_cat, "no hit"),
      n_hits                  = replace_na(n_hits, 0L),
      n_unique_subjects       = replace_na(n_unique_subjects, 0L),
      max_contigs_per_subject = replace_na(max_contigs_per_subject, 0L),
      final_annotation        = NA_character_,
      best_bitscore           = NA_real_,
      best_pident             = NA_real_,
      best_max_coverage       = NA_real_
    )

  for (i in seq_along(annotations$qseqid)) {
    q_id     <- annotations$qseqid[i]
    cat_type <- annotations$prelim_cat[i]

    hits <- blast_df %>% filter(qseqid == q_id)

    # Record best-hit diagnostics
    if (nrow(hits) > 0) {
      annotations$best_bitscore[i]    <- max(hits$bitscore, na.rm = TRUE)
      annotations$best_pident[i]      <- max(hits$pident,   na.rm = TRUE)
      annotations$best_max_coverage[i] <- max(hits$max_coverage, na.rm = TRUE)
    }

    if (verbose) cat("  ", q_id, " [", cat_type, "] -> ")

    # FIX: use classify_query() which uses if/else (not case_when)
    annotation <- tryCatch(
      classify_query(cat_type, hits, blast_df, q_id,
                     overlap_threshold, coverage_threshold),
      error = function(e) {
        message("Error processing query ", q_id, ": ", e$message)
        "error"
      }
    )

    annotations$final_annotation[i] <- annotation
    if (verbose) cat(annotation, "\n")
  }

  return(annotations)
}

# ============================================================================
# 5. FILE-LEVEL PIPELINE
# ============================================================================

#' Process a list of BLAST files and annotate each
annotate_blast_files <- function(file_list, parallel = FALSE, n_workers = NULL,
                                  verbose = FALSE) {

  if (is.null(n_workers)) {
    n_workers <- min(availableCores() - 1, length(file_list))
  }

  process_file <- function(f) {
    tryCatch({
      cat("Processing: ", basename(f), "\n")
      read_blast_outfmt6(f) %>%
        annotate_blast_hits(verbose = verbose) %>%
        mutate(file_name = basename(f))
    }, error = function(e) {
      message("ERROR in ", basename(f), ": ", e$message)
      tibble(
        qseqid              = NA_character_,
        prelim_cat           = NA_character_,
        final_annotation     = NA_character_,
        n_hits               = NA_integer_,
        n_unique_subjects    = NA_integer_,
        max_contigs_per_subject = NA_integer_,
        best_bitscore        = NA_real_,
        best_pident          = NA_real_,
        best_max_coverage    = NA_real_,
        file_name            = basename(f)
      )
    })
  }

  if (parallel && require(furrr, quietly = TRUE)) {
    plan(multisession, workers = n_workers)
    on.exit(plan(sequential), add = TRUE)
    results <- future_map_dfr(file_list, process_file, .progress = TRUE)
  } else {
    results <- purrr::map_df(file_list, process_file)
  }

  return(results)
}

# ============================================================================
# 6. SUMMARY & VISUALIZATION
# ============================================================================

summarize_annotations <- function(annotation_df) {
  list(
    category_counts = annotation_df %>%
      count(final_annotation, sort = TRUE),

    by_category = annotation_df %>%
      count(prelim_cat, final_annotation),

    proportions = annotation_df %>%
      group_by(final_annotation) %>%
      summarise(n = n(), proportion = n() / nrow(annotation_df), .groups = "drop")
  )
}

plot_annotations <- function(annotation_df) {

  my_custom_theme <- function(...) {
    theme_bw(base_family = "sans", base_size = 14) +
      theme(
        legend.position    = "top",
        strip.placement    = "outside",
        strip.background   = element_rect(fill = "gray90", color = "white"),
        strip.text         = element_text(angle = 0, hjust = 0),
        axis.text          = element_text(size = rel(0.7), color = "black"),
        panel.grid         = element_blank(),
        ...
      )
  }

  annotation_df %>%
    count(final_annotation) %>%
    ggplot(aes(x = reorder(final_annotation, -n), y = n, fill = final_annotation)) +
    geom_col() +
    scale_fill_startrek() +
    labs(title = "BLAST Hit Annotations", x = "Category", y = "Count") +
    my_custom_theme() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# ============================================================================
# 7. USAGE EXAMPLE
# ============================================================================

dir <- "~/Documents/Windows/Escritorio/blast_outputs/"

pattern <- "1.blast" # this output is using contig assembled as query and reference as subject

# pattern <- "2.blast" # This output us using reference as query and contig as subject

str(file_list <- list.files(path = dir, pattern = pattern, recursive = T, full.names = TRUE))

# file_list <- file_list[!grepl("MERGEPIPE|PLASS", file_list)] # Huge sizes because many contigs in PLASS

file_list <- file_list[grepl("PLASS", file_list)] # Huge sizes because many contigs in PLASS

outdir <- "~/Documents/GitHub/ConotoxinBenchmark/INPUTS/BLAST_based_annotation_v2_dir/PLASS_dir" # PLASS_dir

dir.create(outdir, recursive = T)

annotation_results <- annotate_blast_files(file_list[2], parallel = F, verbose = T)



process_split_memory_efficient <- function(split_name) {
  cat("Processing", split_name, "...\n")
  
  files <- file_list[grepl(split_name, basename(file_list))]
  
  # Process with memory monitoring
  annotation_results <- annotate_blast_files(files, parallel = FALSE, verbose = T)
  
  filename <- paste0("BLAST_based_overlaps", split_name, ".rds")
  
  filename <- file.path(outdir, filename)
  
  readr::write_rds(annotation_results, file = filename)
  
  cat("Saved:", filename, "\n")
  
  # Cleanup
  rm(annotation_results)
  
  
  gc()
}

# Apply function without storing intermediate results

splits <-  unique(sapply(strsplit(basename(file_list), "_"), `[`, 1))


for (i in splits[c(2:3)]) {
  process_split_memory_efficient(i)
}


f <- list.files(outdir, pattern = ".rds", full.names = T)

annotation_results <- do.call(rbind, lapply(f, read_rds))


# annotation_results <- annotate_blast_files(head(file_list, 3), parallel = FALSE, verbose = TRUE)
# summary_stats <- summarize_annotations(annotation_results)
# print(summary_stats)

