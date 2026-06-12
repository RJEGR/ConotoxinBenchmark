# LOSO analysis

rm(list = ls())

if(!is.null(dev.list())) dev.off()

library(tidyverse)

dir <- "/Users/rjegr/Documents/GitHub/ConotoxinBenchmark/INPUTS/vfolds_loso_resampling_dir/transrate_contigs_dir"

contig_files <- list.files(
  path       = dir,
  pattern    = "contigs\\.csv$",
  recursive  = TRUE,
  full.names = TRUE
)


read_contigs <- function(csv_path) {
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
  sample_id <- gsub("_200x_PE", "", sample_id)
  
  
  cbind(df, tool, sample_id) |> group_by(sample_id, tool)
  
  # fix: if not reference coverage (or crbb?) omit, else return df
  crbb_col    <- if ("CRBB" %in% colnames(df)) "CRBB" else
    if ("has_crb" %in% colnames(df)) "has_crb" else NA
  
  
  }


calculate_metrics <- function(df, reference_coverage_val = c(0.70, 0.80, 0.85, 0.90, 0.95)) {
  
  is_chimeric_value = 0
  
  # Handle scalar input for backward compatibility
  if (length(reference_coverage_val) == 1) {
    reference_coverage_val <- reference_coverage_val
  }
  
  count_false <- function(df, ref_val) {
    
    count_Nsequences <- function() {
      outdir <- "/Users/rjegr/Documents/GitHub/ConotoxinBenchmark/INPUTS/vfolds_loso_resampling_dir"
      f <- list.files(path = outdir, pattern = ".fasta", full.names = T)
      
      my_func <- function(x) { 
        dna <- Biostrings::readDNAStringSet(x)
        structure(
          sum(Biostrings::width(dna) >= 200), 
          names = gsub(".fasta", "", basename(x)))
      }
      
      unlist(lapply(f, my_func))
    }
    
    InputNsequences <- count_Nsequences()
    
    data.frame(InputNsequences) %>% 
      as_tibble(rownames = "sample_id") #%>% 
    # right_join(df) %>%
    # mutate(FN = abs(InputNsequences - TP))
  }
  
  count_false_df <- count_false()
  
  # Get group columns before ungrouping
  group_cols <- dplyr::group_vars(df)
  
  # Ungroup and expand data for each reference_coverage_val
  result <- df %>%
    dplyr::ungroup() %>%
    crossing(coverage_threshold = reference_coverage_val) %>%
    dplyr::mutate(
      # Calculate TP for each threshold
      TP_temp = reference_coverage >= coverage_threshold,
      # Calculate FP for each threshold
      FP_temp = reference_coverage > is_chimeric_value & reference_coverage < coverage_threshold,
      # Handle NA values in reference_coverage
      reference_coverage_clean = ifelse(is.na(hits) & is.na(reference_coverage), 0, reference_coverage)
    ) %>%
    dplyr::group_by(across(all_of(c(group_cols, "coverage_threshold")))) %>%
    dplyr::summarise(
      # Count distinct hits for TP (Recall metric)
      TP = n_distinct(hits[reference_coverage >= coverage_threshold]),
      
      # Count contigs for FP
      FP = sum(reference_coverage_clean > is_chimeric_value & 
                 reference_coverage_clean < coverage_threshold, na.rm = TRUE),
      
      # Count contigs with no reference hits
      TN = sum(reference_coverage_clean <= is_chimeric_value, na.rm = TRUE),
      
      # Raw contig count
      rawcontigs = n(),
      
      .groups = "drop"
    ) %>%
    # Calculate FN for each threshold
    dplyr::group_by(across(all_of(group_cols))) %>%
    
    # dplyr::mutate(
    #   # FN = reference sequences not detected (InputNsequences - TP)
    #   # FN = max(TP) - TP  # Simple approximation based on max TP per group
    #   
    # ) %>%
    left_join(count_false_df) %>%
    mutate(FN = abs(InputNsequences - TP)) %>%
    dplyr::ungroup() %>%
    dplyr::rename(reference_coverage_val = coverage_threshold) %>%
    dplyr::mutate_all(~replace(., is.na(.), 0)) %>%
    dplyr::select(all_of(group_cols), reference_coverage_val, rawcontigs, TP, FP, TN, FN)
  
  return(result)
}


summarise_assembly <- function(csv_path) {
  
  cat(csv_path, "/n")
  
  df <- read_contigs(csv_path)

  df |>
    calculate_metrics()
  

  }

read_contigs(contig_files[6])
summarise_assembly(contig_files[5])

asm_summary <- map_dfr(contig_files, summarise_assembly)

df <- read_contigs(contig_files[7])

read_contigs(contig_files[7]) |>
  calculate_metrics()



###



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

