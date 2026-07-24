#' Calculate Metrics for Multiple Reference Coverage Values
#'
#' @param data A data frame with grouped data from dplyr
#' @param reference_coverage_val A numeric vector of reference coverage thresholds
#'
#' @return A data frame with calculated TP, FP, FN metrics for each coverage value
#'
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' metricsdf <- transratedf %>%
#'   group_by(vfold_set, Assembler) %>%
#'   calculate_metrics(reference_coverage_val = c(0.70, 0.80, 0.85, 0.90, 0.95)) %>%
#'   mutate(
#'     Ratio = TP / FP,
#'     Accuracy = TP / (TP + FN + FP),
#'     Precision = TP / (TP + FP),
#'     Sensitivity = TP / (TP + FN),
#'     Fscore = 2 * (TP) / (2 * (TP) + FP + FN)
#'   )
#' }
#'
#' @export
#' 
calculate_metrics <- function(df, reference_coverage_val = c(0.70, 0.80, 0.85, 0.90, 0.95)) {
  
  is_chimeric_value = 0
  
  # Handle scalar input for backward compatibility
  if (length(reference_coverage_val) == 1) {
    reference_coverage_val <- reference_coverage_val
  }
  
  count_false <- function(df, ref_val) {
    
    count_Nsequences <- function() {
      outdir <- "~/Documents/Windows/Documents/ConotoxinBenchmark/INPUTS/vfolds_resampling_dir/"
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
      as_tibble(rownames = "vfold_set") #%>% 
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


