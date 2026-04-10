# Count the profile of gene superfamilies and multi, and full annotations
# therefore plot dendogram of sequences per method

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

my_custom_theme <- function(base_size = 10, legend_pos = "top", ...) {
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


dir <- "/Users/rjegr/Documents/Windows/Documents/ConotoxinBenchmark/INPUTS/ConoSorter_dir/"

outdir <- "/Users/rjegr/Documents/GitHub/ConotoxinBenchmark/INPUTS/"

recode_to <- c("STRINGTIE","SPADES", "TRINITY","IDBA", "MEGAHIT", "RNABLOOM", "BRIDGER", "TRANSABBYS", "BINPACKER","SOAPDENOVO" ,"CSTONE")

recode_to <- structure(c("StringTie","Spades", "Trinity", "IDBA", "MEGAHIT", "RNA-bloom", "BRIDGER","Transabbys", "BinPacker", "SOAP-denovo", "Cstone"), names = recode_to)

recode_col <- c("full","multi", "fragment","chimera", "error")

recode_col <- structure(c("Full","Multi", "Fragment", "Chimera", "Unannotated"), names = recode_col)

outName <- "protein_1" 

f <- file.path(dir, paste0(outName, "_ConoSorter.rds"))

DB <- read_rds(f) |> mutate(Assembler = gsub(".fa.transdecoder", "", Method)) |>
  mutate(vfold_set = sapply(strsplit(Assembler, "_"), `[`, 1)) |> 
  mutate(Assembler = sapply(strsplit(Assembler, "_"), `[`, 5))

DB |> count(protein_id, sort = T)

dir <- "/Users/rjegr/Documents/Windows/Documents/ConotoxinBenchmark/INPUTS/BLAST_based_annotation_dir/"

f <- list.files(dir, full.names = T, pattern = ".rds")

annotation_results <- do.call(rbind, lapply(f, read_rds)) |>
  mutate(file_name = gsub("_into_Fold[0-9]+[0-9]+.[1|2].blast", "", file_name)) |> 
  dplyr::rename("protein_id" = qseqid, "Assembler" = file_name) |>
  mutate(vfold_set = sapply(strsplit(Assembler, "_"), `[`, 1) ,Assembler = sapply(strsplit(Assembler, "_"), `[`, 5))


dir <- "/Users/rjegr/Documents/Windows/Documents/ConotoxinBenchmark/INPUTS/"

f <- list.files(path = dir, pattern = "curated_nuc_conoServerDB.rds", full.names = T)

conoServerDB <- read_rds(f) %>% dplyr::rename("hits" = "entry_id")

transratedf <- read_tsv( file.path(outdir, "benchmark_assemblers.tsv")) |>
  dplyr::rename("protein_id" = contig_name)

DB2 <- transratedf |> ungroup() |>  
  distinct(protein_id, Assembler, vfold_set, reference_coverage, hits) |>
  left_join(annotation_results, by = c("protein_id", "Assembler", "vfold_set")) 


DB <- DB |> ungroup() |> 
  mutate(protein_id = gsub(".p[0-9]+$", "", protein_id)) |>
  left_join(DB2, by = c("protein_id", "Assembler", "vfold_set")) 

DB <- DB |>
  left_join(conoServerDB, by = c("hits")) 

# here above, we going to plot dot plot of assemblers vs top-hits-gene-superfamilies, using sequences to creates dendogram


colors_ <- c("#357EBD", "#486983", "#FEA31A", "#F1E688", "gray")


DB |>
  drop_na(hits, genesuperfamily) |>
  dplyr::mutate(final_annotation = dplyr::recode_factor(final_annotation, !!!recode_col)) |>
  ggplot(aes(x = reference_coverage, y = genesuperfamily, fill=final_annotation)) +
    # facet_wrap(~ Method) +ggplot2::stat_ecdf()
    ggridges::geom_density_ridges_gradient(
      jittered_points = T,
      position = ggridges::position_points_jitter(width = 0.05, height = 0),
      point_shape = '|', point_size = 0.5, point_alpha = 1, alpha = 0.7) +
  scale_fill_manual(values = colors_) +
  my_custom_theme()


# Dendogram

DB3 <- DB |> 
  # distinct(Assembler, hits) |>
  group_by(hits, genesuperfamily, sequence) |>
  summarise(
    combination = paste(sort(unique(Assembler)), collapse = ","),
    n_assemblers = n_distinct(Assembler),
    .groups = "drop"
  ) |> drop_na(sequence, genesuperfamily)



# How to select representative gene superfamilies sequences?

# Subset the sequence assembled by all the assemblers (19 sequencces by 11 assemblers)
# DB3

DB3 |> count(genesuperfamily, n_assemblers, name = "n_toxins") |>
  ggplot(aes(y = genesuperfamily, x = as.factor(n_assemblers), fill = n_toxins)) +
  geom_tile() + geom_text(aes(label = n_toxins))

DB3 |> filter(n_assemblers == 11)

# See conotoxin prottein fam diversity.R

sequence_df <- clusterize(conoServerDB, seq_type = "sequence")

summarisedf <- conoServerDB %>%  
  # mutate(genesuperfamily = ifelse(is.na(genesuperfamily), sequence_clus, genesuperfamily)) %>% 
  # dplyr::rename("sequence" = seq_type) %>% 
  distinct(genesuperfamily, sequence) %>%
  left_join(sequence_df)


summarise_genesuperfamily <- summarisedf %>%
  dplyr::count(genesuperfamily, sort = T) %>%
  dplyr::rename("seq_number" = "n") 

summarise_clusters <- summarisedf %>%
  distinct(genesuperfamily, sequence_clus) %>%
  dplyr::count(genesuperfamily, sort = T) %>%
  dplyr::rename("clus_number" = "n") 

summarisedf <- left_join(summarise_genesuperfamily, summarise_clusters) %>% mutate(index = clus_number / seq_number)

caption <- "index"

summarisedf %>% 
  drop_na() %>% 
  mutate(index = 1-(clus_number / seq_number)) %>%
  group_by(genesuperfamily) %>%
  ungroup() %>% arrange(desc(index)) %>% 
  mutate(genesuperfamily = factor(genesuperfamily, levels = unique(genesuperfamily))) %>%
  ggplot(aes(y = genesuperfamily, x = index, alpha = seq_number)) +
  geom_col() + labs(x = "1 -(Number of clusters/Number of sequences)",caption = caption) +
  my_custom_theme()

