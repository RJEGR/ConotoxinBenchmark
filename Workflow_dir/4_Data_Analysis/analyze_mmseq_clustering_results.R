#!/usr/bin/env Rscript

# Clustering Results Analysis and Visualization Script
# Analyzes MMseqs2 clustering output for conotoxin diversity studies

#scp -r rgomez@omica.cicese.mx:/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/fernando_pub/1_assembly_dir/Folds_200x_dir/FASTA_DIR/cluster_toxins_mmseq_dir/clustering_results_20260502_013449/all_clustering_results.csv

# all_clustering_results.csv
# Load required libraries


rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)


suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(stringr)
})

# Function to list available CSV files
list_clustering_files <- function(directory = ".") {
  csv_files <- list.files(
    path = directory,
    pattern = ".*cluster_counts\\.csv$|all_clustering_results\\.csv$",
    full.names = TRUE,
    recursive = TRUE
  )
  
  if (length(csv_files) == 0) {
    stop("No clustering CSV files found in directory: ", directory)
  }
  
  cat("Found clustering results files:\n")
  for (i in seq_along(csv_files)) {
    cat(sprintf("  %d. %s\n", i, basename(csv_files[i])))
  }
  
  return(csv_files)
}

# Function to parse dataset names into components
parse_dataset_name <- function(dataset_name) {
  # Example: "Fold01_200x_PE_samples_CSTONE" -> sample="Fold01", assembler="CSTONE"
  # Adjust regex patterns based on your naming convention
  
  parsed_data <- dataset_name %>%
    as_tibble_col("dataset") %>%
    mutate(
      # Extract sample (everything before the first underscore or until _200x pattern)
      sample = str_extract(dataset, "^[^_]+") %>%
        str_replace_na("Unknown"),
      
      # Extract assembler (last component after final underscore, or specific patterns)
      assembler = case_when(
        str_detect(dataset, "CSTONE") ~ "CSTONE",
        str_detect(dataset, "Trinity") ~ "Trinity",
        str_detect(dataset, "SOAPdenovo") ~ "SOAPdenovo",
        str_detect(dataset, "Spades") ~ "Spades",
        str_detect(dataset, "MEGAHIT") ~ "MEGAHIT",
        str_detect(dataset, "Velvet") ~ "Velvet",
        # Extract last component after final underscore as fallback
        TRUE ~ str_extract(dataset, "[^_]+$") %>% str_replace_na("Unknown")
      ),
      
      # Extract additional metadata if present
      coverage = str_extract(dataset, "\\d+x") %>% str_replace("x", "") %>% as.numeric(),
      read_type = case_when(
        str_detect(dataset, "_PE_") ~ "PE",
        str_detect(dataset, "_SE_") ~ "SE",
        TRUE ~ "Unknown"
      )
    )
  
  return(parsed_data)
}

# Function to read and process clustering data
read_clustering_data <- function(file_path) {
  cat("Reading file:", basename(file_path), "\n")
  
  # Read the CSV file
  data <- read_csv(file_path, show_col_types = FALSE)
  
  # Check if it's consolidated format (has dataset column) or individual format
  if ("dataset" %in% colnames(data)) {
    cat("  -> Detected consolidated format\n")
    # Already has dataset column
    processed_data <- data
  } else {
    cat("  -> Detected individual format\n")
    # Extract dataset name from filename
    dataset_name <- basename(file_path) %>%
      str_replace("_cluster_counts\\.csv$", "")
    
    processed_data <- data %>%
      mutate(dataset = dataset_name, .before = 1)
  }
  
  # Parse dataset names
  parsed_info <- parse_dataset_name(processed_data$dataset)
  
  # Combine with original data
  result <- processed_data %>%
    left_join(parsed_info, by = "dataset") %>%
    filter(!is.na(sequence_identity), !is.na(num_clusters)) %>%
    mutate(
      sequence_identity = as.numeric(sequence_identity),
      num_clusters = as.numeric(num_clusters)
    )
  
  cat(sprintf("  -> Processed %d rows for %d datasets\n", 
              nrow(result), 
              length(unique(result$dataset))))
  
  return(result)
}

# Function to create summary statistics
create_summary_stats <- function(data) {
  summary_stats <- data %>%
    group_by(assembler, sequence_identity) %>%
    summarise(
      n_datasets = n(),
      mean_clusters = mean(num_clusters, na.rm = TRUE),
      se_clusters = sd(num_clusters, na.rm = TRUE) / sqrt(n()),
      min_clusters = min(num_clusters, na.rm = TRUE),
      max_clusters = max(num_clusters, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    # Handle cases where SE is NA (single observation)
    mutate(se_clusters = ifelse(is.na(se_clusters), 0, se_clusters))
  
  return(summary_stats)
}

# Function to create the main plot
create_clustering_plot <- function(data, summary_stats) {
  # Create the plot
  p <- ggplot(summary_stats, aes(x = sequence_identity, y = mean_clusters)) +
    
    # Add points for means
    stat_summary(data = data, aes(x = sequence_identity, y = num_clusters),
                 fun = "mean", geom = "point", shape = 1, size = 1) +
    
    # Add error bars
    stat_summary(data = data, aes(x = sequence_identity, y = num_clusters),
                 fun.data = mean_se, geom = "errorbar", width = 0.1, size = 0.2) +
    
    # Add connecting lines
    stat_summary(data = data, aes(x = sequence_identity, y = num_clusters),
                 fun.data = mean_se, geom = "line") +
    
    # Faceting by sample and assembler
    facet_grid(assembler ~ sample, scales = "free_y") +
    
    # Styling
    labs(
      title = "Sequence Clustering Analysis Across Identity Thresholds",
      subtitle = "Number of clusters vs sequence identity cutoff",
      x = "Sequence Identity Threshold",
      y = "Number of Clusters (mean ± SE)",
      caption = paste("Analysis date:", Sys.Date())
    ) +
    
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 10),
      strip.text = element_text(size = 10, face = "bold"),
      panel.grid.minor = element_blank(),
      plot.caption = element_text(size = 8, hjust = 0)
    ) +
    
    # Scale adjustments
    scale_x_continuous(breaks = seq(0.5, 1.0, 0.1)) +
    scale_y_continuous(labels = scales::comma_format())
  
  return(p)
}

# Function to create additional summary plots
create_summary_plots <- function(summary_stats) {
  # Plot 1: Heatmap of cluster numbers
  p1 <- ggplot(summary_stats, aes(x = sequence_identity, y = interaction(sample, assembler), 
                                  fill = mean_clusters)) +
    geom_tile() +
    scale_fill_viridis_c(name = "Mean\nClusters") +
    labs(
      title = "Clustering Heatmap",
      x = "Sequence Identity Threshold",
      y = "Sample × Assembler"
    ) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8))
  
  # Plot 2: Diversity reduction curves
  p2 <- ggplot(summary_stats, aes(x = sequence_identity, y = mean_clusters, 
                                  color = assembler, linetype = sample)) +
    geom_line(size = 0.8) +
    geom_point(size = 1.5) +
    labs(
      title = "Diversity Reduction Curves",
      x = "Sequence Identity Threshold",
      y = "Mean Number of Clusters",
      color = "Assembler",
      linetype = "Sample"
    ) +
    theme_bw() +
    scale_x_continuous(breaks = seq(0.5, 1.0, 0.1))
  
  return(list(heatmap = p1, curves = p2))
}

# Function to save results
save_results <- function(data, summary_stats, plots, output_dir = "clustering_analysis") {
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  
  # Save processed data
  write_csv(data, file.path(output_dir, paste0("processed_data_", timestamp, ".csv")))
  write_csv(summary_stats, file.path(output_dir, paste0("summary_stats_", timestamp, ".csv")))
  
  # Save main plot
  ggsave(
    filename = file.path(output_dir, paste0("clustering_plot_", timestamp, ".png")),
    plot = plots$main,
    width = 12, height = 8, dpi = 300
  )
  
  ggsave(
    filename = file.path(output_dir, paste0("clustering_plot_", timestamp, ".pdf")),
    plot = plots$main,
    width = 12, height = 8
  )
  
  # Save additional plots
  if (!is.null(plots$additional)) {
    ggsave(
      filename = file.path(output_dir, paste0("heatmap_", timestamp, ".png")),
      plot = plots$additional$heatmap,
      width = 10, height = 6, dpi = 300
    )
    
    ggsave(
      filename = file.path(output_dir, paste0("diversity_curves_", timestamp, ".png")),
      plot = plots$additional$curves,
      width = 10, height = 6, dpi = 300
    )
  }
  
  cat("Results saved to:", output_dir, "\n")
  return(output_dir)
}

# Main analysis function
main_analysis <- function(directory = ".", auto_process = FALSE, plot = TRUE) {
  cat("=== Conotoxin Clustering Results Analysis ===\n\n")
  
  # List available files
  csv_files <- list_clustering_files(directory)
  
  # File selection
  if (auto_process || length(csv_files) == 1) {
    selected_files <- csv_files
    cat("\nProcessing all found files...\n\n")
  } else {
    cat("\nSelect files to process (enter numbers separated by spaces, or 'all' for all files):\n")
    selection <- readline("Selection: ")
    
    if (tolower(trimws(selection)) == "all") {
      selected_files <- csv_files
    } else {
      indices <- as.numeric(strsplit(selection, "\\s+")[[1]])
      selected_files <- csv_files[indices[!is.na(indices)]]
    }
  }
  
  if (length(selected_files) == 0) {
    stop("No valid files selected.")
  }
  
  # Process all selected files
  all_data <- map_dfr(selected_files, read_clustering_data)
  
  cat("\n=== Data Summary ===\n")
  cat("Total datasets:", length(unique(all_data$dataset)), "\n")
  cat("Samples:", paste(unique(all_data$sample), collapse = ", "), "\n")
  cat("Assemblers:", paste(unique(all_data$assembler), collapse = ", "), "\n")
  cat("Identity range:", min(all_data$sequence_identity), "-", max(all_data$sequence_identity), "\n")
  cat("Cluster range:", min(all_data$num_clusters), "-", max(all_data$num_clusters), "\n\n")
  
  # Create summary statistics
  summary_stats <- create_summary_stats(all_data)
  
  # Initialize plots variable
  plots <- NULL
  
  # Create plots (conditional based on plot parameter)
  if (plot) {
    cat("Creating visualizations...\n")
    main_plot <- create_clustering_plot(all_data, summary_stats)
    additional_plots <- create_summary_plots(summary_stats)
    
    plots <- list(
      main = main_plot,
      additional = additional_plots
    )
    
    # Display main plot
    print(main_plot)
    
    # Save results with plots
    output_dir <- save_results(all_data, summary_stats, plots)
  } else {
    cat("Skipping plot generation (plot = FALSE)\n")
    
    # Save results without plots (data only)
    output_dir <- save_results(all_data, summary_stats, plots = NULL)
  }
  
  cat("\n=== Analysis Complete ===\n")
  if (plot) {
    cat("Check the plots and saved files in:", output_dir, "\n")
  } else {
    cat("Data processing complete. Files saved in:", output_dir, "\n")
    cat("To generate plots, run: main_analysis(plot = TRUE)\n")
  }
  
  return(list(
    data = all_data,
    summary = summary_stats,
    plots = plots
  ))
}

my_custom_theme <- function(base_size = 14, legend_pos = "top", ...) {
  # base_size = 14
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


dir <- "/Users/rjegr/Documents/GitHub/ConotoxinBenchmark/1_assembly_dir/cluster_toxins_mmseq_dir/"

results <- main_analysis(plot = F,directory = dir)

results$data |> count(assembler)

all_data <- results$data |> filter(!assembler %in% c("MERGEPIPE", "PLASS")) # 

all_data |> count(assembler)

all_data <- all_data |> 
  mutate(assembler = ifelse(grepl("Fold", assembler), "Baseline", assembler))
  # mutate(alpha = ifelse(grepl("Baseline", assembler), 1.0, 1))

recode_to <- c("STRINGTIE","SPADES", "TRINITY","IDBA", "MEGAHIT", "RNABLOOM" ,"BRIDGER", "TRANSABBYS", "BINPACKER","SOAPDENOVO" ,"CSTONE", "TRANSLIG", "Baseline", "PLASS")

recode_to <- structure(c("StringTie","rnaSPAdes", "Trinity", "IDBA", "MEGAHIT", "RNA-Bloom","BRIDGER","Trans-ABySS", "BinPacker", "SOAP-denovo", "Cstone", "TransLiG", "Baseline", "PinguiN (nuclassemble)"), names = recode_to)

all_data <- all_data |> dplyr::mutate(assembler = dplyr::recode_factor(assembler, !!!recode_to))

summary_stats <- create_summary_stats(all_data)  
  # mutate(alpha = ifelse(grepl("Baseline", assembler), 1.0, 1))

discrete_scale <- all_data |> ungroup() |> distinct(assembler) |> pull() |> levels()

n <- length(discrete_scale)


scale_col <- c(ggsci::pal_startrek()(7), ggsci::pal_cosmic()(n-7))

scale_col <- structure(scale_col, names = sort(discrete_scale))

text_data <- all_data |> 
  filter(sequence_identity == 1) |>
  create_summary_stats()




ggplot(summary_stats, aes(x = sequence_identity, y = mean_clusters, 
                          fill = assembler, color = assembler)) +
  # Add error bars
  stat_summary(data = all_data, aes(x = sequence_identity, y = num_clusters),
               fun.data = mean_se, geom = "errorbar", width = 0.01, size = 0.2) +
  # Add connecting lines
  stat_summary(data = all_data, aes(x = sequence_identity, y = num_clusters),
               fun.data = mean_se, geom = "line") +
  # Add points for means
  stat_summary(data = all_data, aes(x = sequence_identity, y = num_clusters),
               fun = "mean", geom = "point", shape = 21, size = 1.5, color = "white", stroke = 0.5) +
  scale_x_reverse() +
  scale_color_manual("", values = scale_col) +
  scale_fill_manual("", values = scale_col) +
  my_custom_theme(legend_pos = "none", base_size = 12) +
  ggrepel::geom_text_repel(data = text_data, aes(label = assembler), 
                           max.overlaps = Inf,
                           nudge_x      = -0.05, 
                           # xlim = c(1.1, 5),
                           # ylim = c(0.80, 2),
                           direction    = "y",
                           # hjust        = -0.5,
                           # vjust = -1, 
                           min.segment.length = 0,
                           segment.curvature = 0.1, 
                           segment.size = 0,
                           size = 1.5, family = "GillSans") +
  labs(x = "Sequence identity", y = "Number of toxin clusters (Mean±SE)") -> psave

outdir <- "~/Documents/GitHub/ConotoxinBenchmark/INPUTS/"

ggsave(psave, filename = 'Assemblers_mmseq_clustering.png', 
       path = outdir, width = 4.5, height = 4.5, dpi = 1000, device = png)

