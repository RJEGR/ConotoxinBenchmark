# LOAD transrate of 67 assemblies
# proccessess contig scores per transcript and hit_group, 
# Bind transcript_id with blast-based annotation, and conoSorter values
# SCOPE: 
# Identity the patterns between blast-based annotation ~ ConoSorter + transrate_score 
# 

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

read_transrate_scores <- function(file_list) {
  
  which_cols <- c(
    "at_skew", 
    "contig_name", 
    "coverage", 
    "cpg_count",
    "cpg_ratio", 
    "eff_count", 
    "eff_length", 
    "gc_skew",
    "in_bridges", 
    "length", 
    "linguistic_complexity_6", 
    "orf_length",
    "p_bases_covered", 
    "p_good", 
    "p_not_segmented", 
    "p_seq_true",
    "prop_gc", 
    "rel_path", 
    "score", 
    "tpm"
  )
  
  cat("\nReading\n")
  cat(str_remove(file_list, paste0("^", dir, "/")))
  cat("\n")
  
  read_csv(file_list) %>%
    mutate(rel_path = str_remove(file_list, paste0("^", dirname(dir), "/"))) |>
    select(any_of(which_cols))
}

dir <- "//wsl.localhost/Debian/home/ricardo/transrate_contigs_dir/transrate_dir/"

str(file_list <- list.files(path = dir, pattern = "contigs.csv", recursive = T, full.names = TRUE))

file_list <- file_list[!grepl("MERGEPIPE", file_list)]

transratedf <- lapply(file_list, read_transrate_scores)

transratedf <- do.call(rbind, transratedf)

DB1 <- transratedf |> 
  mutate(rel_path = basename(dirname(rel_path))) %>%
  # mutate(vfold_set = sapply(strsplit(rel_path, "_"), `[`, 1)) %>%
  # mutate(Assembler = sapply(strsplit(rel_path, "_"), `[`, 5)) |>
  dplyr::rename("protein_id" = contig_name, "Method" = rel_path) |>
  select(protein_id, Method, p_bases_covered, p_good, score)

dir <- "C://Users//cinai/OneDrive/Documentos/GitHub/ConotoxinBenchmark/INPUTS/BLAST_based_annotation_dir/"

f <- list.files(dir, full.names = T, pattern = ".rds")

annotation_results <- do.call(rbind, lapply(f, read_rds)) |>
  mutate(file_name = gsub("_into_Fold[0-9]+[0-9]+.[1|2].blast", "", file_name)) |> 
  dplyr::rename("protein_id" = qseqid, "Method" = file_name)


DB1 <- DB1 |> ungroup() |> 
  left_join(annotation_results, by = c("protein_id", "Method"))

DB1 |> ggplot(aes(p_good)) + 
  geom_histogram() + 
  ggforce::facet_col(~ final_annotation) 

