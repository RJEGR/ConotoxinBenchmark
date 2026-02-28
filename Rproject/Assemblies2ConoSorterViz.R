
# De aquí puedes contar el número de inconsistencias en
#  1) aquellas cosas full o múltiple que en principio deben asignarse a conotox y 
#  2) aquellas asignaciones a conotox que son quimeras qué nunca debierom asignarse
# El analisis debe enfocarse en las 1600 conotoxinas iniciales, por alguna razon los numeros superan los 1600 en ambos, nucleotide y protein, 
# por lo que quiza haya que normalizar a ???

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

outdir <- "C://Users//cinai/OneDrive/Documentos/GitHub/ConotoxinBenchmark/INPUTS/ConoSorter_dir"

outName <- "Nucleotide_1" 

f <- file.path(outdir, paste0(outName, "_ConoSorter.rds"))

DB1 <- read_rds(f)

outName <- "protein_1" 

f <- file.path(outdir, paste0(outName, "_ConoSorter.rds"))

DB2 <- read_rds(f) |> mutate(Method = gsub(".fa.transdecoder", "", Method)) 

dir <- "C://Users//cinai/OneDrive/Documentos/GitHub/ConotoxinBenchmark/INPUTS/BLAST_based_annotation_dir/"

f <- list.files(dir, full.names = T, pattern = ".rds")

annotation_results <- do.call(rbind, lapply(f, read_rds)) |>
  mutate(file_name = gsub("_into_Fold[0-9]+[0-9]+.[1|2].blast", "", file_name)) |> 
  dplyr::rename("protein_id" = qseqid, "Method" = file_name)

# DataViz
# 
# Summarise
# 

DB1 <- DB1 |> ungroup() |> 
  distinct(Method, protein_id, Region, tab) |> 
  # Use right Join to count the number of full contigs
  right_join(annotation_results, by = c("protein_id", "Method")) |> 
  count(Method, Region, tab, final_annotation) |> mutate(Mode = "Nucleotide")
 
# 
# 

DB2 <- DB2 |> ungroup() |> 
  mutate(protein_id = gsub(".p[0-9]+$", "", protein_id)) |>
  distinct(Method, protein_id, Region, tab)  |> 
  # Use right Join to count the number of full contigs
  left_join(annotation_results, by = c("protein_id", "Method")) |> 
  count(Method, Region, tab, final_annotation) |> mutate(Mode = "Protein")



# DB2 <- DB2 |> ungroup() |> distinct(Method, protein_id, Region, tab) |> count(Method, Region, tab) |> mutate(Mode = "Protein")

DB <- rbind(DB1, DB2) |> 
  mutate(vfold_set = sapply(strsplit(Method, "_"), `[`, 1)) |> 
  mutate(Method = sapply(strsplit(Method, "_"), `[`, 5))

rm(DB1,DB2, annotation_results);gc()

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

DB <- DB |>
  filter(final_annotation != "error") |> 
  # mutate(tab = ifelse(is.na(tab), "Not assigned", "Assigned")) |>
  group_by(Method, final_annotation, Mode) %>%
  tally(n, sort = T) 

DB |> write_csv(file = paste0(outdir, "/assemblies2ConoServerSummary.csv"))

DB |>
  # group_by() |> mutate(f = )
  # if calculate mean sd
  # group_by(final_annotation, tab) |> rstatix::get_summary_stats(type = "mean_sd")
  ggplot(aes(y = Method, x = final_annotation, fill = n)) +
  geom_tile() +
  geom_text(aes(label = n), color = "white") +
  facet_grid(~ Mode) +
  # geom_col() +
  my_custom_theme()

DB |>
  filter(final_annotation %in% c("chimera")) |> 
  # group_by() |> mutate(f = )
  # if calculate mean sd
  # group_by(final_annotation, tab) |> rstatix::get_summary_stats(type = "mean_sd")
  ggplot(aes(y = Method, x = final_annotation, fill = n)) +
  geom_tile() +
  geom_text(aes(label = n), color = "white") +
  facet_grid(~ Mode) +
  # geom_col() +
  my_custom_theme()


ggsave(PSAVE, filename = 'UPSET_FOR_PUB.png', path = outdir, width = 5, height = 7, device = png, dpi = 800)
  