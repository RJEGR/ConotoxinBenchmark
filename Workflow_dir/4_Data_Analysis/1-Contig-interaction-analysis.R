# OPTIMIZED CONTIG-BASED STATS ANALYSIS
# Processes contig scores, annotations, and database joins efficiently
# SCOPE: Identity patterns between blast-based annotation ~ ConoSorter + transrate_score 
# OPTIMIZATIONS:
# - Modular functions for better maintainability
# - Efficient data reading and joining
# - Configurable file paths
# - Better error handling and memory management
# - Streamlined workflow

# ============================================================================
# SETUP AND CONFIGURATION
# ============================================================================

rm(list = ls())
if(!is.null(dev.list())) dev.off()
options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

# Load required libraries
required_packages <- c("tidyverse", "data.table", "here")
for(pkg in required_packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# Configuration - Update these paths as needed
CONFIG <- list(
  conoserver_dir = "/Users/rjegr/Documents/Windows/Documents/ConotoxinBenchmark/INPUTS/",
  transrate_dir = "/Users/rjegr/Documents/Windows/Debian/transrate_contigs_dir/transrate_dir/",
  blast_dir = "~/Documents/GitHub/ConotoxinBenchmark/INPUTS/BLAST_based_annotation_v2_dir/",
  conosorter_dir = "/Users/rjegr/Documents/Windows/Documents/ConotoxinBenchmark/INPUTS/ConoSorter_dir/",
  output_dir = "~/Documents/GitHub/ConotoxinBenchmark/INPUTS/",
  protein_name = "protein_1",
  exclude_patterns = c("MERGEPIPE"),
  exclude_proteins = c("P06680")
)

# ============================================================================
# OPTIMIZED FUNCTIONS
# ============================================================================

#' Read TransRate scores efficiently
#' @param file_list Character vector of file paths
#' @param base_dir Base directory for relative path calculation
#' @return data.frame with TransRate scores
read_transrate_scores_optimized <- function(file_list, base_dir) {
  
  # Define required columns
  required_cols <- c(
    "at_skew", "contig_name", "coverage", "cpg_count", "cpg_ratio", 
    "eff_count", "eff_length", "gc_skew", "in_bridges", "length", 
    "linguistic_complexity_6", "orf_length", "p_bases_covered", 
    "p_good", "p_not_segmented", "p_seq_true", "prop_gc", 
    "score", "tpm", "hits", "reference_coverage"
  )
  
  # Read files efficiently using data.table
  cat("Reading", length(file_list), "TransRate files...\n")
  
  results <- map_dfr(file_list, function(file) {
    cat("Processing:", basename(file), "\n")
    
    # Use fread for faster reading
    tryCatch({
      dt <- fread(file, select = intersect(names(fread(file, nrows = 0)), required_cols))
      dt[, rel_path := str_remove(file, paste0("^", dirname(base_dir), "/"))]
      as_tibble(dt)
    }, error = function(e) {
      warning("Failed to read file: ", file, " - ", e$message)
      return(NULL)
    })
  })
  
  return(results)
}

#' Load and process ConoServer database
#' @param conoserver_dir Directory containing ConoServer files
#' @return Processed ConoServer data
load_conoserver_db <- function(conoserver_dir) {
  cat("Loading ConoServer database...\n")
  
  conoserver_file <- list.files(
    path = conoserver_dir, 
    pattern = "curated_nuc_conoServerDB.rds", 
    full.names = TRUE
  )
  
  if(length(conoserver_file) == 0) {
    stop("ConoServer database file not found in: ", conoserver_dir)
  }
  
  conoserver_db <- read_rds(conoserver_file[1]) %>% 
    filter(!proteinid %in% CONFIG$exclude_proteins) %>%
    mutate(
      genesuperfamily = ifelse(is.na(genesuperfamily), "Other", genesuperfamily),
      genesuperfamily = gsub(" superfamily", "", genesuperfamily)
    ) %>%
    rename(
      gs_conoServer = genesuperfamily, 
      hits = entry_id
    ) %>%
    distinct(hits, gs_conoServer, organismlatin, organismdiet, organismregion)
  
  cat("ConoServer entries loaded:", nrow(conoserver_db), "\n")
  return(conoserver_db)
}

#' Load and process BLAST annotation results
#' @param blast_dir Directory containing BLAST results
#' @return Processed BLAST annotations
load_blast_annotations <- function(blast_dir) {
  cat("Loading BLAST annotations...\n")
  
  blast_files <- list.files(blast_dir, full.names = TRUE, pattern = ".rds")
  
  if(length(blast_files) == 0) {
    stop("No BLAST annotation files found in: ", blast_dir)
  }
  
  annotation_results <- map_dfr(blast_files, function(file) {
    tryCatch({
      read_rds(file)
    }, error = function(e) {
      warning("Failed to read BLAST file: ", file, " - ", e$message)
      return(NULL)
    })
  }) %>%
    mutate(file_name = gsub("_into_Fold[0-9]+[0-9]+.[1|2].blast", "", file_name)) %>%
    rename(protein_id = qseqid, Method = file_name)
  
  cat("BLAST annotations loaded:", nrow(annotation_results), "\n")
  return(annotation_results)
}

#' Load and process ConoSorter results
#' @param conosorter_dir Directory containing ConoSorter files
#' @param protein_name Protein name identifier
#' @return Processed ConoSorter data
load_conosorter_results <- function(conosorter_dir, protein_name) {
  cat("Loading ConoSorter results...\n")
  
  conosorter_file <- file.path(conosorter_dir, paste0(protein_name, "_ConoSorter.rds"))
  
  if(!file.exists(conosorter_file)) {
    stop("ConoSorter file not found: ", conosorter_file)
  }
  
  conosorter_data <- read_rds(conosorter_file) %>%
    mutate(
      Method = gsub(".fa.transdecoder", "", Method),
      Superfamily = gsub(" ", "", Superfamily),
      protein_id = gsub(".p[0-9]+$", "", protein_id)
    ) %>%
    rename(gs_conoSorter = Superfamily)
  
  cat("ConoSorter entries loaded:", nrow(conosorter_data), "\n")
  return(conosorter_data)
}

#' Create coverage summary categories
#' @param data Data frame with reference_coverage column
#' @return Data frame with coverage summary
add_coverage_summary <- function(data) {
  data %>%
    mutate(
      coverage_summary = case_when(
        reference_coverage == 1 ~ "100% alignment",
        reference_coverage >= 0.95 ~ ">= 95% alignment",
        reference_coverage >= 0.9 ~ ">= 90% alignment", 
        reference_coverage >= 0.8 ~ ">= 80% alignment",
        TRUE ~ "< 80% alignment"
      )
    )
}

#' Calculate gene superfamily concordance
#' @param conosorter_data ConoSorter results
#' @param transrate_data TransRate data with ConoServer annotations
#' @return Concordance analysis
calculate_gs_concordance <- function(conosorter_data, transrate_data) {
  cat("Calculating gene superfamily concordance...\n")
  
  # Extract distinct gene superfamily data
  gs_df <- transrate_data %>% 
    distinct(Method, protein_id, gs_conoServer) %>% 
    drop_na()
  
  # Process ConoSorter superfamilies and calculate concordance
  gs_concordance_df <- conosorter_data %>%
    mutate(gs_conoSorter = str_split(gs_conoSorter, "_")) %>%
    unnest(gs_conoSorter) %>%
    distinct(Method, protein_id, gs_conoSorter) %>%
    right_join(gs_df, relationship = "many-to-many", by = c("protein_id", "Method")) %>%
    mutate(sf_evidence = ifelse(gs_conoSorter == gs_conoServer, "Concordance", "Ambiguous")) %>%
    distinct(Method, protein_id, sf_evidence)
  
  cat("Concordance analysis completed.\n")
  print(gs_concordance_df %>% count(sf_evidence))
  
  return(gs_concordance_df)
}

#' Main function to join all databases efficiently
#' @return Final integrated dataset
join_databases_optimized <- function() {
  
  # ==========================================
  # STEP 1: Load TransRate data
  # ==========================================
  
  transrate_files <- list.files(
    path = CONFIG$transrate_dir,
    pattern = "contigs.csv",
    recursive = TRUE,
    full.names = TRUE
  )
  
  # Filter out excluded patterns
  for(pattern in CONFIG$exclude_patterns) {
    transrate_files <- transrate_files[!grepl(pattern, transrate_files)]
  }
  
  # Read TransRate data
  transrate_data <- read_transrate_scores_optimized(transrate_files, CONFIG$transrate_dir)
  
  # Add coverage summary
  transrate_data <- add_coverage_summary(transrate_data)
  
  # ==========================================
  # STEP 2: Load reference databases
  # ==========================================
  
  conoserver_db <- load_conoserver_db(CONFIG$conoserver_dir)
  blast_annotations <- load_blast_annotations(CONFIG$blast_dir)
  conosorter_data <- load_conosorter_results(CONFIG$conosorter_dir, CONFIG$protein_name)
  
  # ==========================================
  # STEP 3: Create primary database (DB1)
  # ==========================================
  
  cat("Creating primary database...\n")
  
  DB1 <- transrate_data %>%
    mutate(rel_path = basename(dirname(rel_path))) %>%
    rename(protein_id = contig_name, Method = rel_path) %>%
    left_join(conoserver_db, relationship = "many-to-many") %>%
    left_join(blast_annotations, by = c("protein_id", "Method"))
  
  # ==========================================
  # STEP 4: Calculate concordance
  # ==========================================
  
  gs_concordance_df <- calculate_gs_concordance(conosorter_data, DB1)
  
  # ==========================================
  # STEP 5: Create final integrated database
  # ==========================================
  
  cat("Creating final integrated database...\n")
  
  final_db <- conosorter_data %>%
    ungroup() %>%
    drop_na(Region) %>%
    right_join(DB1, by = c("protein_id", "Method"), relationship = "many-to-many") %>%
    left_join(gs_concordance_df, by = c("protein_id", "Method"), relationship = "many-to-many") %>%
    mutate(
      Fold = sapply(strsplit(Method, "_"), `[`, 1),
      Method = sapply(strsplit(Method, "_"), `[`, 5)
    )
  
  cat("Final database created with", nrow(final_db), "entries.\n")
  
  return(final_db)
}

#' Generate summary statistics and plots
#' @param data Final integrated dataset
generate_summary_analysis <- function(data) {
  cat("Generating summary analysis...\n")
  
  # Summary statistics
  cat("\n=== SUMMARY STATISTICS ===\n")
  cat("Total entries:", nrow(data), "\n")
  cat("Unique proteins:", n_distinct(data$protein_id), "\n")
  cat("Unique methods:", n_distinct(data$Method), "\n")
  
  # Coverage distribution
  if("coverage_summary" %in% names(data)) {
    cat("\nCoverage distribution:\n")
    print(data %>% count(coverage_summary, sort = TRUE))
  }
  
  # Annotation distribution
  if("final_annotation" %in% names(data)) {
    cat("\nAnnotation distribution:\n")
    print(data %>% count(final_annotation, sort = TRUE))
  }
  
  # Region distribution
  if("Region" %in% names(data)) {
    cat("\nRegion distribution:\n")
    print(data %>% drop_na(Region) %>% count(Region, sort = TRUE))
  }
  
  # Superfamily evidence
  if("sf_evidence" %in% names(data)) {
    cat("\nSuperfamily evidence:\n")
    print(data %>% count(sf_evidence, sort = TRUE))
  }
}

# ============================================================================
# MAIN EXECUTION
# ============================================================================

main <- function() {
  cat("Starting optimized contig-based stats analysis...\n")
  start_time <- Sys.time()
  
  tryCatch({
    # Join all databases
    final_database <- join_databases_optimized()
    
    # Generate summary analysis
    generate_summary_analysis(final_database)
    
    # Prepare final output
    output_data <- final_database %>%
      select(-protein_id, -ID, -Conflict, Fold) %>%
      filter(!is.na(Method))  # Remove entries with missing Method
    
    # Export final dataset
    output_file <- file.path(CONFIG$output_dir, "contig-based-data-longer.csv")
    
    cat("\nExporting final dataset...\n")
    write_csv(output_data, output_file)
    
    cat("Analysis completed successfully!\n")
    cat("Output saved to:", output_file, "\n")
    cat("Final dataset contains", nrow(output_data), "entries and", ncol(output_data), "variables.\n")
    
    end_time <- Sys.time()
    cat("Total execution time:", round(difftime(end_time, start_time, units = "secs"), 2), "seconds\n")
    
    return(output_data)
    
  }, error = function(e) {
    cat("ERROR: Analysis failed - ", e$message, "\n")
    traceback()
    return(NULL)
  })
}

# ============================================================================
# EXECUTE ANALYSIS
# ============================================================================

if(interactive()) {
  cat("Running in interactive mode. Call main() to execute analysis.\n")
} else {
  # Run automatically if sourced
  result <- main()
}

# ============================================================================
# OPTIONAL: GENERATE DIAGNOSTIC PLOTS
# ============================================================================

#' Generate diagnostic plots for quality assessment
#' @param data Final integrated dataset
generate_diagnostic_plots <- function(data) {
  
  # Plot 1: P_good distribution by annotation
  if(all(c("p_good", "final_annotation") %in% names(data))) {
    p1 <- data %>%
      drop_na(final_annotation) %>%
      ggplot(aes(p_good)) + 
      geom_histogram(bins = 30, alpha = 0.7) + 
      ggforce::facet_col(~ final_annotation) +
      labs(title = "P_good Distribution by Final Annotation",
           x = "P_good Score", y = "Count") +
      theme_bw()
    
    ggsave("diagnostic_plot_1_pgood_by_annotation.png", p1, width = 12, height = 8, dpi = 300)
  }
  
  # Plot 2: Reference coverage vs P_good
  if(all(c("reference_coverage", "p_good", "final_annotation", "Method") %in% names(data))) {
    p2 <- data %>%
      drop_na(final_annotation) %>%
      ggplot(aes(reference_coverage, p_good, color = final_annotation)) +
      geom_point(alpha = 0.6) +
      facet_wrap(~ Method) +
      labs(title = "Reference Coverage vs P_good by Method",
           x = "Reference Coverage", y = "P_good Score") +
      theme_bw() +
      theme(legend.position = "bottom")
    
    ggsave("diagnostic_plot_2_coverage_vs_pgood.png", p2, width = 14, height = 10, dpi = 300)
  }
  
  # Plot 3: Linguistic complexity by annotation
  if(all(c("linguistic_complexity_6", "final_annotation") %in% names(data))) {
    p3 <- data %>%
      drop_na(final_annotation) %>%
      ggplot(aes(y = final_annotation, x = linguistic_complexity_6, fill = after_stat(x))) +
      ggridges::geom_density_ridges_gradient(alpha = 0.8) +
      scale_fill_viridis_c(option = "C") +
      labs(title = "Linguistic Complexity by Final Annotation",
           y = "Final Annotation", x = "Linguistic Complexity (6-mer)") +
      theme_bw() +
      theme(legend.position = "none")
    
    ggsave("diagnostic_plot_3_linguistic_complexity.png", p3, width = 12, height = 8, dpi = 300)
  }
  
  cat("Diagnostic plots saved to current working directory.\n")
}

result <- main()

# Uncomment to generate diagnostic plots:
# if(exists("result") && !is.null(result)) {
#   generate_diagnostic_plots(result)
# }