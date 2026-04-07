
# De aquí puedes contar el número de inconsistencias en
#  1) aquellas cosas full o múltiple que en principio deben asignarse a conotox y 
#  2) aquellas asignaciones a conotox que son quimeras qué nunca debierom asignarse
# El analisis debe enfocarse en las 1600 conotoxinas iniciales, por alguna razon los numeros superan los 1600 en ambos, nucleotide y protein, 
# por lo que quiza haya que normalizar a ???

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

# outdir <- "C://Users//cinai/OneDrive/Documentos/GitHub/ConotoxinBenchmark/INPUTS/ConoSorter_dir"

dir <- "/Users/rjegr/Documents/Windows/Documents/ConotoxinBenchmark/INPUTS/ConoSorter_dir/"

outdir <- "/Users/rjegr/Documents/GitHub/ConotoxinBenchmark/INPUTS/"

outName <- "Nucleotide_1" 

f <- file.path(dir, paste0(outName, "_ConoSorter.rds"))

DB1 <- read_rds(f)

outName <- "protein_1" 

f <- file.path(dir, paste0(outName, "_ConoSorter.rds"))

DB2 <- read_rds(f) |> mutate(Method = gsub(".fa.transdecoder", "", Method)) 

# dir <- "C://Users//cinai/OneDrive/Documentos/GitHub/ConotoxinBenchmark/INPUTS/BLAST_based_annotation_dir/"

dir <- "/Users/rjegr/Documents/Windows/Documents/ConotoxinBenchmark/INPUTS/BLAST_based_annotation_dir/"

f <- list.files(dir, full.names = T, pattern = ".rds")

annotation_results <- do.call(rbind, lapply(f, read_rds)) |>
  mutate(file_name = gsub("_into_Fold[0-9]+[0-9]+.[1|2].blast", "", file_name)) |> 
  dplyr::rename("protein_id" = qseqid, "Method" = file_name)

# DataViz
# 
# Summarise
# 

DB1 <- DB1 |> ungroup() |> 
  # Whether or not distinct to count the true number of assignments
  # distinct(Method, protein_id, Region, tab) |> 
  # Use right Join to count the number of full contigs
  right_join(annotation_results, by = c("protein_id", "Method")) |> 
  count(Method, Region, tab, final_annotation) |> mutate(Mode = "Nucleotide")
 
# 
# 

DB2 <- DB2 |> ungroup() |> 
  mutate(protein_id = gsub(".p[0-9]+$", "", protein_id)) |>
  # distinct(Method, protein_id, Region, tab)  |> 
  # Use right Join to count the number of full contigs
  left_join(annotation_results, by = c("protein_id", "Method")) |> 
  count(Method, Region, tab, final_annotation) |> mutate(Mode = "Protein")



# DB2 <- DB2 |> ungroup() |> distinct(Method, protein_id, Region, tab) |> count(Method, Region, tab) |> mutate(Mode = "Protein")

DB <- rbind(DB1, DB2) |> 
  mutate(vfold_set = sapply(strsplit(Method, "_"), `[`, 1)) |> 
  mutate(Method = sapply(strsplit(Method, "_"), `[`, 5))

.DB <- DB

rm(DB1,DB2, annotation_results);gc()

extrafont::loadfonts(device = "win")

my_custom_theme <- function(base_size = 14, legend_pos = "top", ...) {
  base_size = 14
  theme_bw(base_family = "Gillsans", base_size = base_size) +
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

# Creates all the combination of domains:
# 

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


DB <- DB |>
  filter(final_annotation != "error") |>
  drop_na(Region) |>
  # mutate(Region = ifelse(is.na(Region), "Not assigned", Region)) |>
  dplyr::mutate(Region = dplyr::recode_factor(Region, !!!recode_to)) |>
  mutate(final_annotation
         = ifelse(final_annotation %in% c("full", "multi"), "Full and Multi", final_annotation)) |>
  group_by(Region, final_annotation, Mode, Method, vfold_set) %>%   tally(n, sort = T) |>
  group_by(Region, final_annotation, Mode, Method) |> rstatix::get_summary_stats(type = "mean_se")

# 

# DB <- DB |>
#   filter(final_annotation != "error") |> 
#   mutate(Region = ifelse(is.na(Region), "Not assigned", Region)) |>
#   # mutate(tab = ifelse(is.na(tab), "Not assigned", "Assigned")) |>
#   group_by(Method, final_annotation, Mode, Region) %>%
#   # filter(vfold_set == "Fold01") |>
#   tally(n, sort = T) 

# DB |> write_csv(file = paste0(outdir, "/assemblies2ConoServerSummary.csv"))


DB |>
  mutate(label = paste0(round(mean), " ± ", round(se, digits = 2), "")) |>
  ggplot(aes(y = Method, x = Mode, fill = mean)) +
  geom_tile() +
  geom_text(aes(label = label), color = "gray60", size = 5) +
  # ggforce::facet_col(final_annotation ~ .) +
  facet_grid(final_annotation ~ Region) +
  scale_fill_continuous(palette = c("#FEE0D2", "#FC9272", "#DE2D26")) +
  my_custom_theme() -> psave

DB |>
  # filter(Region == "Primary precursor") |>
  mutate(xmin = mean-se, xmax = mean+se, xlab = xmax+0.1) %>%
  ggplot(aes(y = Method, x = mean, fill = Mode)) +
  # geom_col(, position = position_dodge()) +
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  geom_errorbar(aes(xmin = xmin, xmax = xmax), 
                position = position_dodge(width = 0.9), width = 0.15, alpha = 0.3, color = "gray20") +
  # geom_text(aes(label = n), color = "white") +
  facet_grid(Region ~ final_annotation, scales = "free_x") +
  # ggforce::facet_col(Mode ~ final_annotation, scales = "free_x") +
  ggsci::scale_fill_aaas() +
  # geom_col() +
  my_custom_theme()


  
ggsave(psave, filename = 'ConoSorter_assemblies.png', path = outdir, width = 7, height = 12, device = png, dpi = 800)
  
DB|> count(Method, Mode, Region, sort = T) 
# Calculate LnRR

lnrr <- function(x,y) {log(x) - log(y)}


# Calculate SE of LnRR using Delta method
# SE_LnRR = sqrt((SE_xT/xT)² + (SE_xC/xC)²)

se_lnrr <- function(se_x, se_y, x, y) {sqrt((se_x / x)^2 + (se_y / y)^2)}

# revisar bien esta formula:

DB |>
  group_by(Method, Region, final_annotation) |>
  pivot_wider(names_from = Mode, values_from = c(mean, se), values_fill = 0) |>
  reframe(
    se_rnrr = se_lnrr(se_Protein, se_Nucleotide, mean_Protein, mean_Nucleotide),
    lnrr = lnrr(mean_Protein, mean_Nucleotide)
         ) |> 
  drop_na() |>
  mutate(label = paste0(round(lnrr, digits = 2), " ± ", round(se_rnrr, digits = 2), "")) |> 
  ggplot(aes(y = Method, x = Region, fill = lnrr)) +
  geom_tile() +
  geom_text(aes(label = label), color = "black", size = 5) +
  # ggforce::facet_col(Method ~ .) +
  facet_grid( ~ final_annotation) +
  scale_fill_continuous(palette = c("blue", "white", "#DE2D26"), na.value = "white") +
  my_custom_theme()
  
# Is significant the differences betweeen Modes

library(rstatix)

DB |>
  group_by(Region, Method) %>%
  rstatix::shapiro_test(mean) %>%
  mutate(gauss = ifelse(p > 0.05, TRUE, FALSE)) 
 
# Priori 

DB %>% 
  # group_by(Region, Method) %>%
  rstatix::kruskal_test(mean ~ Mode) %>%
  adjust_pvalue(method = "none") %>%
  add_significance("p") -> prior.stats

# Posteriori
# 

.DB |>
  group_by(Method) %>%
  pairwise_wilcox_test(n ~ Mode,  conf.level = 0.95, ref.group = 'Protein') %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance() -> post.test

DB |>
  group_by(Method) %>%
  rstatix::pairwise_t_test(mean ~ Mode, ref.group = 'Protein') %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance() -> post.test
