# ============================================================================
# BLAST Annotation for Conotoxins - Optimized & Corrected
# Based on Cahis et al., 2012 with modifications for multidomain toxins
# ============================================================================

# library(dplyr)
# library(tidyr)
# library(ggplot2)
# library(ggsci)
# library(furrr)
# library(IRanges)

rm(list = ls())


Sys.setenv(R_MAX_VSIZE = "100Gb")

install_load <- function(package) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package)
    library(package, character.only = TRUE)
  }
}

# Apply the function to multiple packages
packages <- c("dplyr", "tidyr", "ggplot2", "ggsci", "furrr", "IRanges", "readr")

sapply(packages, install_load)


if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

# ============================================================================
# 1. OVERLAP RATIO CALCULATION (FIXED)
# ============================================================================

overlap_ratio <- function(start1, end1, start2, end2) {
  if(any(is.na(c(start1, end1, start2, end2)))) return(0)
  
  ir1 <- IRanges(start = pmin(start1, end1), end = pmax(start1, end1))
  ir2 <- IRanges(start = pmin(start2, end2), end = pmax(start2, end2))
  
  ovl <- width(pintersect(ir1, ir2, resolve.empty = "none"))
  ovl[is.na(ovl)] <- 0
  
  w1 <- width(ir1)
  w2 <- width(ir2)
  shortest <- pmin(w1, w2)
  
  overlap <- pmax(0, ovl)
  result <- overlap / shortest
  result[is.na(result) | is.infinite(result)] <- 0
  
  return(result)
}

# ============================================================================
# 2. READ AND FILTER BLAST RESULTS
# ============================================================================

read_blast_outfmt6 <- function(f) {
  outfmt6_cols <- c(
    "qseqid", "sseqid", "pident", "length", "mismatches", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen"
  )
  
  readr::read_tsv(f, col_names = FALSE, show_col_types = FALSE) %>%
    setNames(outfmt6_cols) %>%
    mutate(db = basename(f))
}

filter_blast_results <- function(blast_df, min_coverage = 0.5, min_identity = 80) {
  blast_df %>%
    filter(!is.na(qlen) & !is.na(slen) & qlen > 0 & slen > 0) %>%
    mutate(
      query_coverage = length / qlen,
      subject_coverage = length / slen,
      max_coverage = pmax(query_coverage, subject_coverage)
    ) %>%
    filter(max_coverage >= min_coverage & pident >= min_identity)
}

# ============================================================================
# 3. PRELIMINARY CATEGORY ASSIGNMENT
# ============================================================================

assign_preliminary_category <- function(blast_df) {
  if (nrow(blast_df) == 0) {
    return(tibble(qseqid = character(), prelim_cat = character()))
  }
  
  query_stats <- blast_df %>%
    group_by(qseqid) %>%
    summarise(
      n_hits = n(),
      n_unique_subjects = n_distinct(sseqid),
      .groups = "drop"
    )
  
  subject_stats <- blast_df %>%
    group_by(sseqid) %>%
    summarise(
      n_contigs = n_distinct(qseqid),
      .groups = "drop"
    )
  
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
        n_hits == 0 ~ "no hit",
        n_hits == 1 & max_contigs_per_subject == 1 ~ "1->1",
        n_hits > 1 & max_contigs_per_subject == 1 ~ "m->1",
        n_hits == 1 & max_contigs_per_subject > 1 ~ "1->n",
        TRUE ~ "m->n"
      )
    ) %>%
    select(qseqid, prelim_cat)
}

# ============================================================================
# 4. CONOTOXIN-SPECIFIC ANNOTATION LOGIC
# ============================================================================

# For conotoxins, overlapping hits on subject sequences indicate MULTIDOMAIN
# (NOT chimera) because conotoxins often have multiple functional domains

check_hit_overlaps_on_query <- function(hits, threshold = 0.5) {
  hit_intervals <- hits %>%
    mutate(start = pmin(sstart, send), end = pmax(sstart, send)) %>%
    filter(start <= end)
  
  if (nrow(hit_intervals) < 2) return(FALSE)
  
  for (x in seq_len(nrow(hit_intervals) - 1)) {
    for (y in (x + 1):nrow(hit_intervals)) {
      ov <- overlap_ratio(
        hit_intervals$start[x], hit_intervals$end[x],
        hit_intervals$start[y], hit_intervals$end[y]
      )
      if (!is.na(ov) && ov > threshold) return(TRUE)
    }
  }
  
  return(FALSE)
}

check_contig_overlaps_on_subject <- function(blast_df, subject_id, query_id, threshold = 0.5) {
  contigs_on_subject <- blast_df %>%
    filter(sseqid == subject_id) %>%
    select(qseqid, sstart, send) %>%
    distinct() %>%
    mutate(start = pmin(sstart, send), end = pmax(sstart, send)) %>%
    filter(start <= end)
  
  contig_ids <- unique(contigs_on_subject$qseqid)
  if (length(contig_ids) < 2) return(FALSE)
  
  for (a in seq_along(contig_ids[-length(contig_ids)])) {
    for (b in (a + 1):length(contig_ids)) {
      i1 <- contigs_on_subject %>% filter(qseqid == contig_ids[a])
      i2 <- contigs_on_subject %>% filter(qseqid == contig_ids[b])
      
      if (nrow(i1) == 0 || nrow(i2) == 0) next
      
      ov <- overlap_ratio(i1$start[1], i1$end[1], i2$start[1], i2$end[1])
      if (ov > threshold) return(TRUE)
    }
  }
  
  return(FALSE)
}

annotate_blast_hits <- function(blast_df, overlap_threshold = 0.5, coverage_threshold = 0.9) {
  
  blast_df <- filter_blast_results(blast_df)
  
  if (nrow(blast_df) == 0) {
    return(tibble(qseqid = character(), prelim_cat = character(), final_annotation = character()))
  }
  
  prelim_cats <- assign_preliminary_category(blast_df)
  
  annotations <- tibble(qseqid = unique(blast_df$qseqid)) %>%
    left_join(prelim_cats, by = "qseqid") %>%
    mutate(
      prelim_cat = replace_na(prelim_cat, "no hit"),
      final_annotation = NA_character_
    )
  
  for (i in seq_along(annotations$qseqid)) {
    q_id <- annotations$qseqid[i]
    cat_type <- annotations$prelim_cat[i]
    
    hits <- blast_df %>% filter(qseqid == q_id)
    
    annotation <- tryCatch({
      case_when(
        cat_type == "no hit" ~ 
          "no hit",
        
        cat_type == "1->1" ~ 
          {
            coverage <- hits$length[1] / hits$slen[1]
            if (coverage >= coverage_threshold) "full" else "fragment"
          },
        
        cat_type == "m->1" ~ 
          {
            subject_id <- hits$sseqid[1]
            overlaps <- check_contig_overlaps_on_subject(blast_df, subject_id, q_id, overlap_threshold)
            if (overlaps) "allele" else "fragment"
          },
        
        cat_type == "1->n" ~ 
          {
            # CONOTOXIN SPECIFIC: Overlapping hits = multidomain, not chimera
            overlaps <- check_hit_overlaps_on_query(hits, overlap_threshold)
            if (overlaps) "multi" else "chimera"
          },
        
        cat_type == "m->n" ~ 
          {
            contig_ids <- n_distinct(hits$qseqid)
            subj_ids <- n_distinct(hits$sseqid)
            
            if (contig_ids == subj_ids) {
              coverage_vec <- hits %>%
                group_by(qseqid) %>%
                slice_max(length, n = 1) %>%
                mutate(coverage = length / slen) %>%
                pull(coverage)
              
              if (mean(coverage_vec >= coverage_threshold) >= 0.5) "full" else "fragment"
            } else {
              "multi"
            }
          },
        
        TRUE ~ "other"
      )
    }, error = function(e) {
      message("Error processing query ", q_id, ": ", e$message)
      "error"
    })
    
    annotations$final_annotation[i] <- annotation
  }
  
  return(annotations)
}

# ============================================================================
# 5. PIPELINE
# ============================================================================

annotate_blast_files <- function(file_list, parallel = FALSE, n_workers = NULL) {
  
  if (is.null(n_workers)) {
    n_workers <- min(availableCores() - 1, length(file_list))
  }
  
  process_file <- function(f) {
    tryCatch({
      cat("Processing: ", basename(f), "\n")
      read_blast_outfmt6(f) %>%
        annotate_blast_hits() %>%
        mutate(file_name = basename(f))
    }, error = function(e) {
      message("ERROR in ", basename(f), ": ", e$message)
      tibble(
        qseqid = NA_character_,
        prelim_cat = NA_character_,
        final_annotation = NA_character_,
        file_name = basename(f)
      )
    })
  }
  
  if (parallel && require(furrr, quietly = TRUE)) {
    plan(multisession, workers = n_workers)
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
      summarise(n = n(), proportion = n()/nrow(annotation_df), .groups = "drop")
  )
}

#' Plot annotation results
plot_annotations <- function(annotation_df) {
  
  my_custom_theme <- function(...) {
    theme_bw(base_family = "sans", base_size = 14) +
      theme(
        legend.position = "top",
        strip.placement = "outside",
        strip.background = element_rect(fill = 'gray90', color = 'white'),
        strip.text = element_text(angle = 0, hjust = 0),
        axis.text = element_text(size = rel(0.7), color = "black"),
        panel.grid = element_blank(),
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
# 7. USAGE
# ============================================================================

# annotation_results <- annotate_blast_files(head(file_list, 3), parallel = FALSE)
# summary_stats <- summarize_annotations(annotation_results)
# print(summary_stats)

# ============================================================================
# USAGE EXAMPLE
# ============================================================================

dir <- "//wsl.localhost/Debian/home/ricardo/blast_outputs/"

pattern <- "1.blast" # this output is using contig assembled as query and reference as subject

# pattern <- "2.blast" # This output us using reference as query and contig as subject

str(file_list <- list.files(path = dir, pattern = pattern, recursive = T, full.names = TRUE))

file_list <- file_list[!grepl("MERGEPIPE", file_list)] # Huge sizes because many contigs in PLASS

file_list <- file_list[!grepl("PLASS", file_list)] # Huge sizes because many contigs in PLASS

outdir <- "C://Users//cinai/OneDrive/Documentos/GitHub/ConotoxinBenchmark/INPUTS/BLAST_based_annotation_dir/"

splits <-  unique(sapply(strsplit(basename(file_list), "_"), `[`, 1))


process_split_memory_efficient <- function(split_name) {
  cat("Processing", split_name, "...\n")
  
  files <- file_list[grepl(split_name, basename(file_list))]
  
  # Process with memory monitoring
  annotation_results <- annotate_blast_files(files, parallel = TRUE)
  
  filename <- paste0("BLAST_based_overlaps", split_name, ".rds")
  
  filename <- file.path(outdir, filename)
  
  readr::write_rds(annotation_results, file = filename)
  
  cat("Saved:", filename, "\n")
  
  # Cleanup
  rm(annotation_results)
  
  
  gc()
}

# Apply function without storing intermediate results
for (i in splits) {
  process_split_memory_efficient(i)
}


f <- list.files(outdir, full.names = T)

annotation_results <- do.call(rbind, lapply(f, read_rds))


# annotation_results <- annotate_blast_files(file_list, parallel = TRUE)

# summarize_annotations(annotation_results)

# plot_annotations(annotation_results)


DataViz <- annotation_results %>%
  filter(final_annotation != "error") %>%
  group_by(file_name) %>%
  count(final_annotation) %>%
  mutate(vfold_set = sapply(strsplit(file_name, "_"), `[`, 1)) %>%
  mutate(file_name = sapply(strsplit(file_name, "_"), `[`, 5)) %>%
  mutate(file_name = gsub(".[1|2].blast", "", file_name)) %>%
  group_by(vfold_set, file_name) %>%
  # mutate(n = n/sum(n)) %>%
  group_by(final_annotation, file_name) %>%
  rstatix::get_summary_stats(type = "mean_sd")


# annotation_results %>%
#   mutate(vfold_set = sapply(strsplit(file_name, "_"), `[`, 1)) %>%
#   mutate(file_name = sapply(strsplit(file_name, "_"), `[`, 5)) %>%
#   mutate(file_name = gsub(".[1|2].blast", "", file_name)) %>%
#   group_by(file_name) %>%
#   summarize_annotations()

# annotation_results %>%
#   mutate(vfold_set = sapply(strsplit(file_name, "_"), `[`, 1)) %>%
#   mutate(file_name = sapply(strsplit(file_name, "_"), `[`, 5)) %>%
#   mutate(file_name = gsub(".[1|2].blast", "", file_name)) %>%
#   group_by(file_name) %>%
#   plot_annotations() +
#   facet_wrap(~ file_name)



extrafont::loadfonts(device = "win")

my_custom_theme <- function(base_size = 14, legend_pos = "top", ...) {
  base_size = 14
  theme_bw(base_family = "Gill Sans MT", base_size = base_size) +
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



DataViz %>%
  ggplot(aes(y = file_name, x = n, fill = final_annotation)) +
  geom_col() +
  scale_fill_startrek() +
  my_custom_theme() +
  labs(x = "Fraction of Assemblies")


# viz 2
# 

DataViz <- annotation_results %>%
  filter(final_annotation != "error") %>%
  mutate(vfold_set = sapply(strsplit(file_name, "_"), `[`, 1)) %>%
  mutate(file_name = sapply(strsplit(file_name, "_"), `[`, 5)) %>%
  mutate(file_name = gsub(".[1|2].blast", "", file_name)) %>%
  count(file_name, final_annotation)

DataViz %>%
  mutate(x = mean) %>%
  mutate(label = paste0(file_name, " (",  round(x, 3), ")")) %>%
  filter(final_annotation%in% c("chimera", "fragment")) %>%
  mutate(xmin = x-sd, xmax = x+sd, xlab = xmax+0.1) %>%
  ggplot(aes(y = file_name, x = n, fill = final_annotation)) +
  facet_grid(~ final_annotation) +
  geom_col(aes(y = file_name, x = x), width = 0.7) + 
  # geom_text(aes(y = file_name, x = xlab, label = label), hjust = 0.1, size = 3) +
  geom_errorbar(aes(y = file_name, x = x, xmin = xmin, xmax = xmax), width = 0.15, alpha = 0.3, color = "black") + 
  scale_fill_startrek() +
  my_custom_theme() +
  labs(x = "N of Assemblies") +
  guides(fill = guide_legend(title = "", byrow = T))


