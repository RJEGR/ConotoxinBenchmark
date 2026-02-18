
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
  # filter(Assembler %in% assemblers) %>%
  filter(!is.na(hits)) %>%
  filter(reference_coverage > 0.5) %>%
  mutate(summarise = "< 80 % alignment") %>%
  mutate(summarise = ifelse(reference_coverage >= 0.8, ">= 80% alignment", summarise)) %>%
  mutate(summarise = ifelse(reference_coverage >= 0.9, ">= 90% alignment", summarise)) %>%
  mutate(summarise = ifelse(reference_coverage >= 0.95, ">= 95% alignment", summarise)) %>%
  mutate(summarise = ifelse(reference_coverage == 1, "100% alignment", summarise)) %>%
  dplyr::distinct(Assembler, summarise, reference_coverage, hits)


UPSETDF <- fltr_data %>%
  filter(reference_coverage > 0.95) %>%
  group_by(hits, summarise) %>% # omit summarise
  dplyr::distinct(Assembler) %>% 
  summarise(across(Assembler, .fns = list), n = n()) %>% arrange(desc(n))



n_pallet <- length(unique(fltr_data$summarise))

scale_col <- ggsci::pal_uchicago(alpha = 0.8)(n_pallet) 

scale_col <- structure(scale_col, names = sort(unique(fltr_data$summarise)))

scale_fill <- ggsci::pal_uchicago(alpha = 0.8)(n_pallet) 

scale_fill <- structure(scale_fill, names = sort(unique(fltr_data$summarise)))


library(ggupset)

# UPSETDF %>% mutate(Assembler = sapply(Assembler     , paste, collapse = "-")) %>%
#   group_by(Assembler) %>%
#   summarise(elements = list(hits), count = n())

P <- UPSETDF %>%
  # filter(n <= 2) %>%
  # mutate(facet = "C) Intersected") %>%
  ggplot(aes(x = Assembler)) +
  # geom_bar(position = position_dodge(width = 1), color = "black", linewidth = 0.2) +
  geom_bar(position = position_stack(), linewidth = 0.2, width = 0.5) +
  # geom_segment(aes(xend = CONTRAST_DE, yend = n, y = Inf, color = SIGN), linewidth = 1) +
  geom_text(stat='count', aes(label = after_stat(count)),color = "black",
            position = position_stack(), vjust = -0.2, family = "Gill Sans MT", size = 2.5) +
  ggupset::scale_x_upset(order_by = "freq", reverse = F) +
  axis_combmatrix(sep = ";", override_plotting_function = function(df){
    ggplot(df, aes(x= at, y= single_label)) +
      geom_rect(aes(fill= index %% 2 == 0), ymin=df$index-0.5,
                ymax=df$index+0.5, xmin=0, xmax=1) +
      geom_point(aes(color= observed), size = 3) +
      # geom_line(data= function(dat) dat[dat$observed, ,drop=FALSE],
      #           aes(group = labels), size= 1.2) +
      ylab("") + xlab("") +
      scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
      scale_fill_manual(values= c(`TRUE` = "white", `FALSE` = "#F7F7F7")) +
      scale_color_manual(values= c(`TRUE` = "black", `FALSE` = "#E0E0E0")) +
      # my_custom_theme()
      guides(color="none", fill="none") +
      theme(
        panel.background = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.length = unit(0, "pt"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.line = element_blank(),
        panel.border = element_blank()
      )
  })

P

# https://albert-rapp.de/posts/ggplot2-tips/26_upset_charts/26_upset_charts.html

# built_plot <- ggplot_build(P)
# 
# intersection_data <- built_plot$data[[1]]

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

# Manually edit upset 

fltr_data %>%
  group_by(hits, summarise) %>% # omit summarise
  dplyr::distinct(Assembler) %>% 
  summarise(across(Assembler, .fns = paste_), n = n()) %>% arrange(desc(n)) %>% 
  ggplot(aes(x=Assembler)) +
  axis_combmatrix(sep = ";")
  # 
fltr_data %>%
  group_by(hits, summarise) %>% # omit summarise
  dplyr::distinct(Assembler) %>% 
  summarise(across(Assembler, .fns = paste_), n = n()) %>% arrange(desc(n)) %>% 
  ggplot(aes(x=Assembler)) +
  geom_bar() +
  # axis_combmatrix(sep = ";")
  axis_combmatrix(sep = ";", override_plotting_function = function(df){
    ggplot(df, aes(x= at, y= single_label)) +
      geom_rect(aes(fill= index %% 2 == 0), ymin=df$index-0.5,
                ymax=df$index+0.5, xmin=0, xmax=1) +
      geom_point(aes(color= observed), size = 3) +
      # geom_line(data= function(dat) dat[dat$observed, ,drop=FALSE],
      #           aes(group = labels), size= 1.2) +
      ylab("") + xlab("") +
      scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
      scale_fill_manual(values= c(`TRUE` = "white", `FALSE` = "#F7F7F7")) +
      scale_color_manual(values= c(`TRUE` = "black", `FALSE` = "#E0E0E0")) +
      # my_custom_theme()
      guides(color="none", fill="none") +
      theme(
        panel.background = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.length = unit(0, "pt"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.line = element_blank(),
        panel.border = element_blank()
      )
  })

# ggplot(UPSETDF, aes(x=Assembler)) +
#   geom_bar() +



UPSETDF <- fltr_data %>%
  group_by(summarise) %>% # omit summarise
  # dplyr::distinct(Assembler) %>% 
  reframe(across(Assembler, .fns = paste_), n = n(), mean_se(reference_coverage )) %>% arrange(desc(n))

# Sort by assembler groups of intersection
# 

Assembler_levs <- UPSETDF %>%
  ungroup() %>%
  dplyr::count(Assembler, sort = T)

UPSETDF %>%
  mutate(Assembler = factor(Assembler, levels = Assembler_levs$Assembler)) %>%
  ggplot() +
  geom_errorbar(aes(y = Assembler, x = y, xmin = ymin , xmax = ymax), width = 0.15, alpha = 0.3) 

# Describe gene families founb by assembler
# 
# 

fltr_data <- transratedf %>% 
  # filter(Assembler %in% assemblers) %>%
  filter(!is.na(hits)) %>%
  filter(reference_coverage > 0.5) %>%
  mutate(summarise = "< 80 % alignment") %>%
  mutate(summarise = ifelse(reference_coverage >= 0.8, ">= 80% alignment", summarise)) %>%
  mutate(summarise = ifelse(reference_coverage >= 0.9, ">= 90% alignment", summarise)) %>%
  mutate(summarise = ifelse(reference_coverage >= 0.95, ">= 95% alignment", summarise)) %>%
  mutate(summarise = ifelse(reference_coverage == 1, "100% alignment", summarise)) %>%
  dplyr::distinct(Assembler, summarise, reference_coverage, hits)


sf_levs <- conoServerDB %>% 
  # drop_na(genesuperfamily)  %>%
  count(split_as, sort = T) 

fltr_data %>%
  filter(reference_coverage >  0.95) %>%
  # count number of isoforms per gene
  count( summarise, Assembler, hits, sort = T) %>%
  left_join(distinct(conoServerDB, hits, split_as) ) %>%
  # Count number of gene sf per isoform group
  drop_na(split_as)  %>%
  group_by(split_as, summarise, Assembler) %>% tally(n, sort = T) %>%
  mutate(split_as = factor(split_as, levels = sf_levs$split_as)) %>%
  # ggplot(aes(y = Assembler, x = n)) + geom_line()
  ggplot(aes(y = split_as, x = Assembler, fill = n, label = n)) +
  facet_grid(~ summarise) +
  geom_tile() +
  geom_text() +
  theme(axis.text.x = element_text(angle = 90, size = 7, hjust = 1))

