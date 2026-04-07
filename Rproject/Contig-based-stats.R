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
  cat(stringr::str_remove(file_list, paste0("^", dir, "/")))
  cat("\n")
  
  read_csv(file_list) %>%
    mutate(rel_path = stringr::str_remove(file_list, paste0("^", dirname(dir), "/"))) |>
    select(any_of(which_cols))
}

dir <- "/Users/rjegr/Documents/Windows/Debian/transrate_contigs_dir/transrate_dir/"

str(file_list <- list.files(path = dir, pattern = "contigs.csv", recursive = T, full.names = TRUE))

file_list <- file_list[!grepl("MERGEPIPE", file_list)]

transratedf <- lapply(file_list, read_transrate_scores)

transratedf <- do.call(rbind, transratedf)

DB1 <- transratedf |> 
  mutate(rel_path = basename(dirname(rel_path))) %>%
  # mutate(vfold_set = sapply(strsplit(rel_path, "_"), `[`, 1)) %>%
  # mutate(Assembler = sapply(strsplit(rel_path, "_"), `[`, 5)) |>
  dplyr::rename("protein_id" = contig_name, "Method" = rel_path) |>
  select(protein_id, Method, p_bases_covered, p_good, score, linguistic_complexity_6)

dir <- "/Users/rjegr/Documents/Windows/Documents/ConotoxinBenchmark/INPUTS/BLAST_based_annotation_dir/"

f <- list.files(dir, full.names = T, pattern = ".rds")

annotation_results <- do.call(rbind, lapply(f, read_rds)) |>
  mutate(file_name = gsub("_into_Fold[0-9]+[0-9]+.[1|2].blast", "", file_name)) |> 
  dplyr::rename("protein_id" = qseqid, "Method" = file_name)


DB1 <- DB1 |> ungroup() |> 
  left_join(annotation_results, by = c("protein_id", "Method"))

DB1 |> ggplot(aes(p_good)) + 
  geom_histogram() + 
  ggforce::facet_col(~ final_annotation) 

# LOAD conoSorter protein results
# SCOPE: Identify if causal relation between contig score and conosorter annotation

# Assemblies2ConoSorterViz.R

creates_domains <- function() {
  
  library(stringr)
  
  base_domains <- c("Mature", "Pro-region", "Signal")
  
  # Generate all possible combinations of length 1, 2, and 3
  all_combos <- c()
  
  for (k in 1:length(base_domains)) {
    combos <- combn(base_domains, k, simplify = FALSE)
    formatted <- sapply(combos, function(x) {
      paste0("(", paste(x, collapse = ")_("), ")")
    })
    all_combos <- c(all_combos, formatted)
  }
  
  domains <- all_combos
  print(domains)
}

names_to <-  c(creates_domains()[creates_domains() %in% "(Mature)_(Pro-region)_(Signal)"],
               creates_domains()[!creates_domains() %in% "(Mature)_(Pro-region)_(Signal)"])

recode_to <-structure(
  c("Primary precursor", rep("Incomplete precursor", 6)),
  names = names_to
)

dir <- "/Users/rjegr/Documents/Windows/Documents/ConotoxinBenchmark/INPUTS/ConoSorter_dir/"

outName <- "protein_1" 

f <- file.path(dir, paste0(outName, "_ConoSorter.rds"))

DB2 <- read_rds(f) |> mutate(Method = gsub(".fa.transdecoder", "", Method)) 


DB1 |> distinct(Method, protein_id)
DB2 |> distinct(Method, protein_id)

DB <- DB2 |> ungroup() |> 
  mutate(protein_id = gsub(".p[0-9]+$", "", protein_id)) |>
  drop_na(Region) |>
  dplyr::mutate(Region = dplyr::recode_factor(Region, !!!recode_to)) |>
  # distinct(Method, protein_id, Region, tab)  |> 
  # Use right Join to count the number of full contigs
  right_join(DB1, by = c("protein_id", "Method")) 

DB |> 
  # drop_na() |>
  ggplot(aes(p_good, fill =Region)) + 
  geom_histogram() + 
  ggforce::facet_col(~ final_annotation, scales = "free_y") 

DB |>
  count(Region, final_annotation)


DB |>
  ggplot(aes(y = Superfamily, x = linguistic_complexity_6, fill = after_stat(x))) +
  # geom_violin() +
  # facet_grid(~ Superfamily) +
  ggridges::geom_density_ridges_gradient(
    jittered_points = F,
    position = ggridges::position_points_jitter(width = 0.05, height = 0),
    point_shape = '|', point_size = 3, point_alpha = 1, alpha = 1) +
  scale_fill_viridis_c(option = "C") +
  labs(y = "", x = "Reference coverage (TransRate)") +
  theme_bw(base_family = "GillSans", base_size = 12) + 
  theme(legend.position = "none", 
        # axis.text.y = element_blank(), 
        # axis.ticks.y = element_blank(), 
        axis.title.x = element_text(size = 7)) +
  scale_x_continuous(position = "top")

DB |>
  mutate(Method = sapply(strsplit(Method, "_"), `[`, 5)) |>
  drop_na(final_annotation, Region) |>
  ggplot(aes(linguistic_complexity_6, p_good, color = final_annotation)) +
  geom_point() +
  facet_grid(Method ~ final_annotation)
