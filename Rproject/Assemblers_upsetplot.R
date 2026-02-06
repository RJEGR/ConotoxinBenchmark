
# January 2026
# LOAD transratedf derived form assemblers.R
# FORMAT input to group conotoxins based on intersection between assemblers
# considered filter reference_coverage  > 0.95 or not
# for those groups, estimates mean se reference_coverage  values to second facet of contig coverage plot, preserving assembler groups
# plot upset plot in facet1 and join to facet2
# # considered to label w colors gene superfamilies or not..


rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)


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

outdir <- "C://Users//cinai/OneDrive/Documentos/GitHub/ConotoxinBenchmark/INPUTS"


f <- list.files(path = outdir, pattern = "curated_nuc_conoServerDB.rds", full.names = T)

conoServerDB <- read_rds(f) %>% dplyr::rename("hits" = "entry_id")


dir <- "C://Users//cinai/OneDrive/Escritorio/transrate_contigs_reference_dir/"


transratedf <- read_tsv( file.path(dir, "benchmark_assemblers.tsv"))



# library(ggVennDiagram)
# 
assemblers <- c("STRINGTIE","SPADES", "TRINITY", "RNABLOOM", "IDB", "MEGAHIT")

fltr_data <- transratedf %>% 
  filter(Assembler %in% assemblers) %>%
  dplyr::filter(!is.na(hits)) %>%
  dplyr::filter(reference_coverage > 0.9) %>%
  dplyr::distinct(Assembler, reference_coverage, hits) 


fltr_data <- transratedf %>% 
  filter(Assembler %in% assemblers) %>%
  filter(!is.na(hits)) %>%
  filter(reference_coverage > 0.5) %>%
  mutate(summarise = "< 80 % alignment") %>%
  mutate(summarise = ifelse(reference_coverage >= 0.8, ">= 80% alignment", summarise)) %>%
  mutate(summarise = ifelse(reference_coverage >= 0.9, ">= 90% alignment", summarise)) %>%
  mutate(summarise = ifelse(reference_coverage >= 0.95, ">= 95% alignment", summarise)) %>%
  mutate(summarise = ifelse(reference_coverage == 1, "100% alignment", summarise)) %>%
  dplyr::distinct(Assembler, summarise, reference_coverage, hits)


UPSETDF <- fltr_data %>%
  group_by(hits, summarise) %>% # omit summarise
  dplyr::distinct(Assembler) %>% 
  summarise(across(Assembler, .fns = list), n = n()) %>% arrange(desc(n))



n_pallet <- length(unique(fltr_data$summarise))

scale_col <- ggsci::pal_uchicago(alpha = 0.8)(n_pallet) 

scale_col <- structure(scale_col, names = sort(unique(fltr_data$summarise)))

scale_fill <- ggsci::pal_uchicago(alpha = 0.8)(n_pallet) 

scale_fill <- structure(scale_fill, names = sort(unique(fltr_data$summarise)))


library(ggupset)

UPSETDF %>%
  # filter(n <= 2) %>%
  # mutate(facet = "C) Intersected") %>%
  ggplot(aes(x = Assembler, fill = summarise, color = summarise)) +
  # geom_bar(position = position_dodge(width = 1), color = "black", linewidth = 0.2) +
  geom_bar(position = position_stack(), linewidth = 0.2, width = 0.5) +
  # geom_segment(aes(xend = CONTRAST_DE, yend = n, y = Inf, color = SIGN), linewidth = 1) +
  geom_text(stat='count', aes(label = after_stat(count)),color = "black",
            position = position_stack(), vjust = -0.2, family = "Gill Sans MT", size = 2.5) +
  ggupset::scale_x_upset(order_by = "freq", reverse = F) +
  ggupset::theme_combmatrix(combmatrix.panel.point.color.fill = "black",
                            combmatrix.label.make_space = FALSE,
                            combmatrix.panel.line.size = 0, base_family = "Gill Sans MT") +
  axis_combmatrix(levels = assemblers, clip = "off") +
  scale_color_manual("",values = scale_col ) +
  scale_fill_manual("",values = scale_fill) +
  my_custom_theme()

# 
# 
# 
# 
# 

paste_ <- function(x) { 
  x <- x[!is.na(x)] 
  x <- unique(sort(x))
  x <- paste(x, sep = ';', collapse = ';')
}

# UPSETDF <- fltr_data %>%
  group_by(summarise) %>% # omit summarise
  # dplyr::distinct(Assembler) %>% 
  reframe(across(Assembler, .fns = paste_), n = n(), mean_se(reference_coverage )) %>% arrange(desc(n))

# Sort by assembler groups
# 

Assembler_levs <- UPSETDF %>%
  ungroup() %>%
  dplyr::count(Assembler, sort = T)

UPSETDF %>%
  mutate(Assembler = factor(Assembler, levels = Assembler_levs$Assembler)) %>%
  ggplot() +
  geom_errorbar(aes(y = Assembler, x = y, xmin = ymin , xmax = ymax), width = 0.15, alpha = 0.3) 
  
