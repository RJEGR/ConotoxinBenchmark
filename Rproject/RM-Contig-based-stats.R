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
    "tpm",
    "hits",
    "reference_coverage"
  )
  
  cat("\nReading\n")
  cat(stringr::str_remove(file_list, paste0("^", dir, "/")))
  cat("\n")
  
  read_csv(file_list) %>%
    mutate(rel_path = stringr::str_remove(file_list, paste0("^", dirname(dir), "/"))) |>
    dplyr::select(any_of(which_cols))
}


dir <- "/Users/rjegr/Documents/Windows/Documents/ConotoxinBenchmark/INPUTS/"

f <- list.files(path = dir, pattern = "curated_nuc_conoServerDB.rds", full.names = T)

conoServerDB <- read_rds(f) %>% 
  # Filter special case:
  filter(proteinid != "P06680") |>
  mutate(genesuperfamily = ifelse(is.na(genesuperfamily), "Other", genesuperfamily)) |>
  mutate(genesuperfamily = gsub(" superfamily", "", genesuperfamily)) |>
  dplyr::rename("gs_conoServer" = genesuperfamily, "hits" = "entry_id") |>
  distinct(hits, gs_conoServer, organismlatin, organismdiet, organismregion) # organismlatin, organismdiet, 

dir <- "/Users/rjegr/Documents/Windows/Debian/transrate_contigs_dir/transrate_dir/"

str(file_list <- list.files(path = dir, pattern = "contigs.csv", recursive = T, full.names = TRUE))

file_list <- file_list[!grepl("MERGEPIPE", file_list)]

transratedf <- lapply(file_list, read_transrate_scores)

# transratedf[[7]] |> distinct(rel_path)

transratedf[[7]] <- NULL

# read_csv(file_list[[7]])

transratedf <- do.call(rbind, transratedf) |> 
  mutate(summarise = "< 80 % alignment") %>%
  mutate(summarise = ifelse(reference_coverage >= 0.8, ">= 80% alignment", summarise)) %>%
  mutate(summarise = ifelse(reference_coverage >= 0.9, ">= 90% alignment", summarise)) %>%
  mutate(summarise = ifelse(reference_coverage >= 0.95, ">= 95% alignment", summarise)) %>%
  mutate(summarise = ifelse(reference_coverage == 1, "100% alignment", summarise))

DB1 <- transratedf |> 
  mutate(rel_path = basename(dirname(rel_path))) %>%
  # mutate(vfold_set = sapply(strsplit(rel_path, "_"), `[`, 1)) %>%
  # mutate(Assembler = sapply(strsplit(rel_path, "_"), `[`, 5)) |>
  dplyr::rename("protein_id" = contig_name, "Method" = rel_path) |>
  left_join(conoServerDB, relationship = "many-to-many")
  # dplyr::select(protein_id, Method, p_bases_covered, p_good, score, linguistic_complexity_6, reference_coverage, summarise)

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

dir <- "/Users/rjegr/Documents/Windows/Documents/ConotoxinBenchmark/INPUTS/ConoSorter_dir/"

outName <- "protein_1" 

f <- file.path(dir, paste0(outName, "_ConoSorter.rds"))

DB2 <- read_rds(f) |> mutate(Method = gsub(".fa.transdecoder", "", Method)) |> 
  mutate(Superfamily = gsub(" ", "", Superfamily), 
         protein_id = gsub(".p[0-9]+$", "", protein_id)) |>
  dplyr::rename("gs_conoSorter" = Superfamily) 

gs_df <- DB1 |> distinct(Method, protein_id, gs_conoServer) |> drop_na()

any(gs_df$protein_id %in% DB2$protein_id)
any(gs_df$Method %in% DB2$Method)

gs_concordance_df <- DB2 |> 
  mutate(
    gs_conoSorter = str_split(gs_conoSorter, "_")) |>
  unnest(gs_conoSorter) |>
  distinct(Method, protein_id, gs_conoSorter) |>
  right_join(gs_df, relationship = "many-to-many", by = c("protein_id", "Method")) |>
  mutate(sf_evidence = ifelse(gs_conoSorter == gs_conoServer, "Concordance", "Ambiguous")) |>
  distinct(Method, protein_id, sf_evidence)

gs_concordance_df |> count(sf_evidence)

DB1 |> distinct(Method, protein_id)
DB2 |> distinct(Method, protein_id)


any(DB1$protein_id %in% DB1$protein_id)
any(DB1$Method %in% DB2$Method)

DB <- DB2 |> ungroup() |> 
  # mutate(protein_id = gsub(".p[0-9]+$", "", protein_id)) |>
  drop_na(Region) |>
  # dplyr::mutate(Region = dplyr::recode_factor(Region, !!!recode_to)) |>
  # distinct(Method, protein_id, Region, tab)  |> 
  # Use right Join to count the number of full contigs
  right_join(DB1, by = c("protein_id", "Method"), relationship = "many-to-many") |>
  left_join(gs_concordance_df, by = c("protein_id", "Method"), relationship = "many-to-many") 

DB |> 
  drop_na(Region) |>
  ggplot(aes(p_good, fill =Region)) + 
  geom_histogram() + 
  ggforce::facet_col(~ sf_evidence, scales = "free_y") 

DB |>
  count(Region, final_annotation)


DB |>
  drop_na(final_annotation) |>
  ggplot(aes(y = final_annotation, x = linguistic_complexity_6, fill = after_stat(x))) +
  # geom_violin() +
  # facet_grid(~ Superfamily) +
  ggridges::geom_density_ridges_gradient(
    jittered_points = F,
    position = ggridges::position_points_jitter(width = 0.05, height = 0),
    point_shape = '|', point_size = 3, point_alpha = 1, alpha = 1) +
  scale_fill_viridis_c(option = "C") +
  labs(y = "", x = "linguistic_complexity_6") +
  theme_bw(base_family = "GillSans", base_size = 12) + 
  theme(legend.position = "none", 
        # axis.text.y = element_blank(), 
        # axis.ticks.y = element_blank(), 
        axis.title.x = element_text(size = 7)) +
  scale_x_continuous(position = "top")

DB |>
  mutate(Method = sapply(strsplit(Method, "_"), `[`, 5)) |>
  drop_na(final_annotation, Region) |>
  ggplot(aes(reference_coverage, p_good)) +
  # facet_grid(Method ~ final_annotation) +
  geom_point() 


# FIND BY LM what is the feature (predictor) explain well the accurate annotation (full or multi or hit > 95%)

DB <- DB |> 
  mutate(Fold = sapply(strsplit(Method, "_"), `[`, 1), 
         Method = sapply(strsplit(Method, "_"), `[`, 5))

outdir <- "~/Documents/GitHub/ConotoxinBenchmark/INPUTS/"

DB |> 
  select(-protein_id, -ID, -Conflict, Fold) |>
  # group_by(Method) |>
  # sample_n(10) |>
  write_csv(file.path(outdir, "contig-based-data-longer.csv"))
  # write_tsv(file.path(outdir, "contig-based-data.tsv"))

