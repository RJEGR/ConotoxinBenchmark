
# This is the mac version
# Read transrate scores (contigs.txt) using read_transrate_scores() (to evals accuracy of assembly methods)
# Output transrateDB containing follow columns
# cd /Users/cigom/Documents/GitHub/ConotoxinBenchmark/1_assembly
# scp -r rgomez@omica.cicese.mx:/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/fernando_pub/1_assembly_dir/Folds_200x_dir/transrate_contigs_dir .

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)


# extrafont::loadfonts(device = "win")

my_custom_theme <- function(base_size = 14, legend_pos = "top", ...) {
  base_size = 14
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

read_transrate_scores <- function(file_list) {
  
  # dir <- dirname(dir)
  rel_path <- str_remove(file_list, paste0("^", dir, "/"))
  
  cat("\nReading\n")
  cat(rel_path)
  cat("\n")
  
  rel_path <- sapply(strsplit(dirname(rel_path), "/"), `[`, 9)
  vfold_set <- sapply(strsplit(rel_path, "_"), `[`, 1)
  Assembler <- sapply(strsplit(rel_path, "_"), `[`, 5)
  
  cat("\n")
  cat(vfold_set, Assembler)
  cat("\n")

  read_csv(file_list) |>
    # mutate(file_list = file_list)
    mutate(rel_path, vfold_set, Assembler)
}

recode_to <- c("STRINGTIE","SPADES", "TRINITY","IDBA", "MEGAHIT", "RNABLOOM", "BRIDGER", "TRANSABBYS", "BINPACKER","SOAPDENOVO" ,"CSTONE", "TRANSLIG")

recode_to <- structure(c("StringTie","Spades", "Trinity", "IDBA", "MEGAHIT", "RNA-bloom", "BRIDGER","Transabbys", "BinPacker", "SOAP-denovo", "Cstone", "TransLiG"), names = recode_to)


outdir <- "~/Documents/Windows/Documents/ConotoxinBenchmark/INPUTS/"

f <- list.files(path = outdir, pattern = "curated_nuc_conoServerDB.rds", full.names = T)

conoServerDB <- read_rds(f) |> dplyr::rename("hits" = "entry_id")

dir <- "~/Documents/Windows/Escritorio/transrate_contigs_reference_dir/"

str(file_list <- list.files(path = dir, pattern = "contigs.csv", recursive = T, full.names = TRUE))

file_list <- file_list[!grepl("MERGEPIPE", file_list)]

transratedf <- lapply(file_list, read_transrate_scores)

transratedf <- do.call(rbind,transratedf)

transratedf |>
  group_by(vfold_set,Assembler) |>
  dplyr::count()  

transratedf |>
  group_by(Assembler) |>
  summarise(mean = mean(length), min = min(length), max = max(length))
  # ggplot(aes(length)) +
  # facet_grid(Assembler ~., scales = "free", space = "free") +
  # geom_histogram()

transratedf |>
  dplyr::mutate(Assembler = dplyr::recode_factor(Assembler, !!!recode_to)) |>
  write_tsv(file = file.path(dir, "benchmark_assemblers.tsv"))



metricsdf <- transratedf %>% 
  filter(length >= 200) %>%
  group_by(vfold_set, Assembler) %>%
  calculate_metrics(reference_coverage_val = c(0.5,0.6,0.70, 0.80, 0.85, 0.90, 0.95, 1)) %>% 
  mutate(
    Ratio = TP / FP,
    Accuracy = TP / (TP + FN + FP),
    Precision = TP / (TP + FP),
    Sensitivity = TP / (TP + FN),
    Fscore = 2 * (TP) / (2 * (TP) + FP + FN)
  )

# View metrics across all thresholds

# metricsdf %>% 
#   group_by(Assembler, reference_coverage_val) %>%
#   summarise(across(c(TP, FP, Accuracy, Precision, Sensitivity), mean))
# 


metricsdf <- metricsdf |> dplyr::mutate(Assembler = dplyr::recode_factor(Assembler, !!!recode_to))

discrete_scale <- metricsdf |> ungroup() |> distinct(Assembler) |> pull() |> levels()

n <- length(discrete_scale)

scale_col <- c(ggsci::pal_startrek()(7), ggsci::pal_cosmic()(n-7))

scale_col <- structure(scale_col, names = sort(discrete_scale))


# PLOT


text_data <- metricsdf |> 
  filter(reference_coverage_val == 0.5) |>
  group_by(Assembler, reference_coverage_val) |>
  summarise(Accuracy = mean(Accuracy))


metricsdf %>% 
  ggplot(aes(group = Assembler, x = as.factor(reference_coverage_val), y = Accuracy, color = Assembler)) +
  # geom_text(data = text_data, aes(label = Assembler),  hjust = -0.3, color = "gray70", angle = 360, size = 3) +
  ggrepel::geom_text_repel(data = text_data, aes(label = Assembler), 
                           max.overlaps = Inf,
                           nudge_x      = -0.7, 
                           # xlim = c(1.1, 5),
                           # ylim = c(0.80, 2),
                           direction    = "y",
                           # hjust        = -0.5,
                           # vjust = -1, 
                           min.segment.length = 0,
                           segment.curvature = 0.1, 
                           segment.size = 0,
                           size = 1.5, family = "GillSans") +
  # geom_jitter(position = position_jitter(0.1), shape = 1) +
  stat_summary(fun = "mean", geom = "point", shape = 1, size =1) +
  stat_summary(fun.data=mean_se, geom="errorbar", width = 0.1, size = 0.2) +
  stat_summary(fun.data=mean_se, geom="line") +
  geom_hline(yintercept=0.5, linetype="dashed",color = "gray", size=0.5) +
  geom_vline(xintercept = "0.95", linetype="dashed", color = "gray", size=0.5) +
  scale_color_manual("", values = scale_col) +
  # scale_fill_manual("", values = scale_col) +
  my_custom_theme(legend_pos = "none") +
  labs(x = "Sequence identity threshold") -> psave

ggsave(psave, filename = 'Assemblers_Accuracy_thresholds.png', 
       path = outdir, width = 4.5, height = 4.5, dpi = 1000, device = png)


metricsdf_flt <- metricsdf |> 
  filter(reference_coverage_val == 0.95) |>
  group_by(Assembler) |> 
  # filter-out outliers
  mutate(Sensitivity = ifelse(rstatix::is_outlier(Sensitivity), NA, Sensitivity)) |>
  mutate(Precision = ifelse(rstatix::is_outlier(Precision), NA, Precision)) |>
  drop_na(Sensitivity, Precision)
  
mark_circle <- metricsdf_flt |>
  summarise(
    sens_sd = sd(Sensitivity), pre_sd = sd(Precision), n = n(),
    Sensitivity = mean(Sensitivity), Precision = mean(Precision))

line_data <- metricsdf_flt |> 
  filter(Assembler == "Trinity") |>
  summarise(Sensitivity = max(Sensitivity), Precision = max(Precision), label = "Baseline")


metricsdf |>
  filter(reference_coverage_val == 0.95) |>
  group_by(Assembler) |> 
  mutate(Sensitivity = ifelse(rstatix::is_outlier(Sensitivity), NA, Sensitivity)) |>
  mutate(Precision = ifelse(rstatix::is_outlier(Precision), NA, Precision)) |>
  drop_na(Sensitivity, Precision) |>
  mutate(facet = "A) Benchmark metrics") |>
  ggplot(aes(x = Sensitivity, y = Precision)) +  
  facet_grid(~ facet) +
  geom_segment(data = line_data, aes(xend=Sensitivity, yend=0),  
    color = "gray40", linetype = "dashed", linewidth = 0.5) +
  geom_segment(data = line_data, aes(xend=0,  yend=Precision),  
    color = "gray40", linetype = "dashed", linewidth = 0.5) +
  geom_text(data = line_data, aes(x = Sensitivity, y = 0, label = label),  
    hjust = -0.5, vjust = 1.5,color = "gray50", angle = 90, size = 3) +
  geom_point(aes(color = Assembler), alpha = 0.5) +
  scale_y_continuous(limits = c(0,1)) + scale_x_continuous(limits = c(0,1)) +
  scale_color_manual("", values = scale_col) +
  ggforce::geom_mark_circle(
    data = mark_circle,
    aes(
      Sensitivity, Precision,
      group = Assembler, 
      # label = Assembler, 
      fill= Assembler, color = Assembler),
    label.buffer = unit(0.5, 'lines'),
    # label.margin = margin(0, 0, 0, 0, "mm"),
    label.family = "GillSans",
    label.colour = "inherit",
    label.fontsize = 5,
    con.colour = "inherit",
    con.cap = 0,
    con.type = "none", #elbow
    con.size = 0.5,
    con.border = "none",
    con.linetype = 2,
    expand = unit(2.5, "mm"), #radius = expand
    alpha = 0.2) +
  my_custom_theme(legend_pos = "none") -> p

p <- p + 
  ggrepel::geom_text_repel(data = mark_circle, 
    aes(label = Assembler, color = Assembler), 
    max.overlaps = 100, family = "GillSans", size = 2.5,
      
    # nudge_x = .15,
    box.padding = 0.5,
    # nudge_y = 1,
    segment.curvature = -0.1
    # segment.ncp = 3,
    # segment.angle = 20
    ) 

p

# ggsave(p, filename = 'Assemblers_Precision_Sensitivity.png', 
#   path = outdir, width = 4.5, height = 4.5, dpi = 1000, device = png)



p2 <- metricsdf |> filter(reference_coverage_val == 0.95) |> 
  group_by(Assembler) |> 
  summarise(sd = sd(Accuracy), x = mean(Accuracy)) |>
  arrange(desc(x)) |>
  mutate(label = paste0(Assembler, " (",  round(x, 3), ")")) |>
  mutate(Assembler = factor(Assembler, levels = rev(unique(Assembler)))) |>
  mutate(facet = "B) Accuracy") |>
  mutate(xmin = x-sd, xmax = x+sd, xlab = xmax+0.1) |>
  mutate(xlab = ifelse(Assembler == "StringTie", 0.25, xlab)) |> 
  ggplot() + 
  facet_grid(~ facet) +
  geom_col(aes(y = Assembler, x = x, fill= Assembler, color = Assembler), width = 0.7) + 
  geom_text(aes(y = Assembler, x = xlab, label = label), hjust = 0.1, size = 2) +
  geom_errorbar(aes(y = Assembler, x = x, xmin = xmin, xmax = xmax), width = 0.15, alpha = 0.3, color = "gray20") +
  scale_x_continuous("",limits = c(0,1)) +
  scale_y_discrete(position = "right")+
  scale_color_manual("", values = scale_col) +
  scale_fill_manual("", values = scale_col) +
  my_custom_theme(legend_pos = "none", axis.ticks.y = element_blank(), 
    axis.text.y = element_blank(), axis.title.y = element_blank())  
  
# theme(axis.ticks = element_blank(), 
  #       axis.text.x = element_text(size = 7.5),
  #       panel.border = element_blank(),
  #       panel.grid.major.x = element_blank(),
  #       panel.grid.minor.y = element_blank()
  # )


# p2


library(patchwork)


psave <- p + inset_element(p2, left = 0.6, bottom = 0.3, right = 0.95, top = 0.75, align_to = 'full')

outdir <- "~/Documents/GitHub/ConotoxinBenchmark/INPUTS/"

ggsave(psave, filename = 'Assemblers_Precision_Sensitivity.png', 
       path = outdir, width = 5, height = 5, dpi = 1000, device = png)


#####
#####
#####
# Exit

#####
#####
#####

metricsdf |> group_by(Assembler) |> 
  select(Precision) |>
  rstatix::get_summary_stats(type = "mean_se") |>
  arrange(desc(mean)) 


metricsdf |> group_by(Assembler) |> 
  summarise(sd = sd(Ratio), x = mean(Ratio)) |>
  arrange(desc(x)) |>
  mutate(Assembler = factor(Assembler, levels = unique(Assembler))) |> #View()
  ggplot(aes(y = Assembler, x = x)) + 
  geom_col(aes(fill= Assembler, color = Assembler), width = 0.7) + 
  geom_text(aes(label = round(x, 3)), hjust = -0.15, size = 3.5) +
  scale_x_continuous("Correct/incorrect ratio (C/I)",limits = c(0,12)) +
  scale_color_manual("", values = scale_col) +
  scale_fill_manual("", values = scale_col) +
  my_custom_theme(legend_pos = "none")


# Evaluates quantiles of accuracy
# 

transratedf |>
  ggplot(aes(reference_coverage, fill = Assembler, color = Assembler)) +
  stat_ecdf(size = 5) +
  scale_color_manual("", values = scale_col) +
  scale_fill_manual("", values = scale_col) +
  my_custom_theme(legend_pos = "top")

assembler_lev <- metricsdf |> group_by(Assembler) |> 
  summarise(sd = sd(Accuracy), x = mean(Accuracy)) |>
  arrange(desc(x)) |> pull(Assembler)

p3 <- transratedf |>
  # left_join(conoServerDB,by = "hits")
  mutate(Assembler = factor(Assembler, levels = assembler_lev)) |>
  ggplot(aes(y = Assembler , x = reference_coverage, fill = Assembler, color = Assembler, alpha = after_stat(x))) +
  ggridges::geom_density_ridges_gradient(
  jittered_points = T,
  position = ggridges::position_points_jitter(width = 0.05, height = 0),
  point_shape = '|', point_size = 0.25, point_alpha = 1, alpha = 1) +
  scale_color_manual("", values = scale_col) +
  scale_fill_manual("", values = scale_col) +
  # scale_y_reverse() +
  scale_x_continuous("Contig coverage",breaks = seq(0,1, by = 0.2)) +
  my_custom_theme(legend_pos = "none")

ggsave(p3, filename = 'Assemblers_contig_coverage.png', 
       path = outdir, width = 3, height = 3.5, dpi = 1000, device = png)


transratedf |>
  drop_na(hits) |>
  distinct(Assembler, hits, reference_coverage) |>
  left_join(distinct(conoServerDB, hits, split_as),by = "hits") |>
  ggplot(aes(y = split_as  , x = reference_coverage)) +
  geom_boxplot()

# cor(transratedf$linguistic_complexity_6, transratedf$reference_coverage)

transratedf |>
  mutate(col = ifelse(is.na(hits), "Not annotated", "Annotated")) |>
  sample_frac(size = 0.15) |>
  ggplot(aes(linguistic_complexity_6, reference_coverage, color = col)) +
  geom_point()

transratedf |>
  drop_na(hits) |>
  distinct(Assembler, hits, reference_coverage) |>
  left_join(distinct(conoServerDB, hits, split_as),by = "hits") |>
  # mutate(Assembler = factor(Assembler, levels = assembler_lev)) |>
  ggplot(aes(y = split_as  , x = reference_coverage, alpha = after_stat(x))) +
  ggridges::geom_density_ridges_gradient(
    jittered_points = F,
    position = ggridges::position_points_jitter(width = 0.05, height = 0),
    point_shape = '|', point_size = 0.25, point_alpha = 1, alpha = 1) +
  # scale_y_reverse() +
  scale_x_continuous("Contig coverage",breaks = seq(0,1, by = 0.2)) +
  my_custom_theme(legend_pos = "none")


DataViz <- transratedf |> 
  drop_na() |> 
  distinct(vfold_set, hits, reference_coverage, Assembler) |>
  # left_join(conoServerDB,by = "hits") |>
  mutate(summarise = "< 80 % alignment") |>
  mutate(summarise = ifelse(reference_coverage >= 0.8, ">= 80% alignment", summarise)) |>
  mutate(summarise = ifelse(reference_coverage >= 0.9, ">= 90% alignment", summarise)) |>
  mutate(summarise = ifelse(reference_coverage >= 0.95, ">= 95% alignment", summarise)) |>
  mutate(summarise = ifelse(reference_coverage == 1, "100% alignment", summarise)) |>
  dplyr::count(Assembler, summarise)


n_pallet <- length(unique(DataViz$summarise))

scale_col <- ggsci::pal_uchicago(alpha = 0.8)(n_pallet) 

scale_col <- structure(scale_col, names = sort(unique(DataViz$summarise)))

scale_fill <- ggsci::pal_uchicago(alpha = 0.8)(n_pallet) 

scale_fill <- structure(scale_fill, names = sort(unique(DataViz$summarise)))


p4 <- DataViz  |>
  # arrange(Assembler, summarise)
  group_by(Assembler) |>
  mutate(n = n/sum(n)) |> 
  # summarise(sum(n)) |>
  group_by(Assembler,summarise ) |>
  rstatix::get_summary_stats(type = "mean_se") |>
  mutate(Assembler = factor(Assembler, levels = rev(assembler_lev))) |>
  mutate(x = mean) |>
  mutate(facet = "C) Assembler coverage") |>
  ggplot() +
  facet_grid(~ facet) +
  geom_col(aes(y = Assembler, x = x, fill= summarise), width = 0.7) + 
  xlab("Frac. of contigs") +
  # scale_color_manual("",values = scale_col ) +
  scale_fill_manual("",values = scale_fill) +
  my_custom_theme(legend.text = element_text(size = 7), legend_pos = "right")
  

ggsave(p4, filename = 'Assemblers_coverage_bar.png', 
       path = outdir, width = 5, height = 5, dpi = 1000, device = png)


DataViz  |>
  # arrange(Assembler, summarise)
  group_by(Assembler,vfold_set ) |>
  mutate(n = n/sum(n)) |> 
  # summarise(sum(n)) |>
  group_by(Assembler,summarise ) |>
  rstatix::get_summary_stats(type = "mean_se") |>
  mutate(Assembler = factor(Assembler, levels = assembler_lev)) |>
  mutate(x = mean) |>
  # mutate(facet = "B) Accuracy") |>
  mutate(xmin = x-se, xmax = x+se, xlab = xmax+0.1) |>
  ggplot() +
  facet_grid(summarise ~.) +
  geom_col(aes(y = Assembler, x = x, fill= summarise, color = summarise), width = 0.7) + 
  # geom_text(aes(y = Assembler, x = xlab, label = label), hjust = 0.1, size = 3) +
  geom_errorbar(aes(y = Assembler, x = x, xmin = xmin, xmax = xmax), width = 0.15, alpha = 0.3, color = "gray20") 

##
##

DataViz  |>
  # summarise(sum(n)) |>
  group_by(Assembler,vfold_set ) |>
  mutate(n = n/sum(n)) |> 
  group_by(Assembler,summarise ) |>
  rstatix::get_summary_stats(type = "mean_se") |>
  mutate(Assembler = factor(Assembler, levels = assembler_lev)) |>
  ggplot() +
  geom_tile(color = 'white', linewidth = 0.2, aes(y = Assembler, x= summarise, fill = mean)) +
  scale_fill_viridis_c("", option = "inferno", direction = -1)


transratedf |> 
  drop_na() |> 
  distinct(vfold_set, hits, reference_coverage, Assembler) |>
  # filter(Assembler == "TRINITY") |>
  group_by(Assembler) |>
  rstatix::get_summary_stats(type = "quantile")

library(patchwork)
# p + p2

# Inner join of conotoxins, and assemblers,


# file_out <- file.path(outdir, "transrate_assemblers.rds")

# write_rds(transratedf, file = file_out)


# According to Cahis et al., 2012
# By assembler, assess/Classify hits, to fragmented, chimeric, allelic, paralogue, and other genomic, based on overlap -----

dir <- "C://Users//cinai/OneDrive/Escritorio/blast_outputs//"

# dir <- "/Users/cigom/Documents/GitHub/ConotoxinBenchmark/1_assembly/blast_outputs/"

str(f <- list.files(path = dir, pattern = "1.blast", recursive = T, full.names = TRUE))

blast_df <- read_outfmt6(f[1])

# blast_df <- do.call(rbind, lapply(f, read_outfmt6))

# Hits were con- sidered significant when (i) the alignment length (merg- ing all high-scoring segments pairs) was at least 70% of the query sequence or at least 70% of the hit sequence, and (ii) sequence identity between query and hit was more than 70% across the aligned portion.

# Fragmented: 
# Contigs with a single significant hit shared by other contigs were called (allele)

blast_df |> 
  filter(identity >= 95) |>
  group_by(db) |>
  dplyr::count(target, sort = T) |>
  mutate()
# dplyr::rename("allele" = "n") 

# (Quimera or multi) Contigs with several significant hits, all specific to this contig, were called



summarise_df <- transratedf |> 
  filter(!is.na(hits)) |>
  filter(reference_coverage > 0.5) |>
  mutate(summarise = "< 80 % alignment") |>
  mutate(summarise = ifelse(reference_coverage >= 0.8, ">= 80% alignment", summarise)) |>
  mutate(summarise = ifelse(reference_coverage >= 0.9, ">= 90% alignment", summarise)) |>
  mutate(summarise = ifelse(reference_coverage >= 0.95, ">= 95% alignment", summarise)) |>
  mutate(summarise = ifelse(reference_coverage == 1, "100% alignment", summarise)) |>
  group_by(summarise, Assembler, vfold_set) 

## TRansrate (if read-based metrics are found)

names(transratedf) %in% c("in_bridges", "p_good", "")

# The contig will represent a single transcript, s(Cseg). This score mea- sures the probability that the coverage depth of the transcript is univariate, i.e., that it represents an assembly of a single tran- script and not a hybrid/chimeric assembly of multiple tran- scripts expressed at different expression levels. Here, the per- nucleotide coverage depth of the contig must be best modeled by a single Dirichlet distribution (described below). If the contig is better modeled by the product of two or more Dirichlet distri- butions, then this indicates that two or more contigs with differ- ent transcript abundances have been erroneously assembled together.

#  Similarly, transcripts with low s(Cseg) scores are likely to represent chimeric transcripts. Here, al- though the transcript itself may be incorrectly assembled, the com- ponent segments of the transcript may themselves be correctly assembled and of utility if separated. To help users identify and diagnose likely assembly errors affecting low scoring contigs, TransRate provides each of the separate contig scores (in addition to the overall contig score). This

transratedf |> 
  filter(!is.na(hits)) |>
  # filter(reference_coverage < 0.5) |>
  # ggplot(aes(in_bridges)) + geom_histogram()
  # ggplot(aes(p_good, score, alpha = reference_coverage)) + geom_point()
  ggplot(aes(p_not_segmented, tpm)) + geom_point()

# Correlate 
transratedf |> 
  filter(!is.na(hits)) |> 
  select_if(is.double) |>
  mutate_all(~replace(., is.na(.), 0)) |>
  cor(method = "spearman") -> M


testRes <- corrplot::cor.mtest(M, conf.level = 0.95)
# 
corrplot::corrplot(M, p.mat = testRes$p ,method = "color", type="upper", order = "hclust", insig = "label_sig")



# transratedf |> 

true_contigs_db <- transratedf |> filter(!is.na(hits)) 

true_contigs_db |> distinct(hits) |> nrow()

# all_contigs_db <- transratedf |> filter(is.na(hits)) 

sum(sort(conoServerDB$hits) %in% sort(unique(true_contigs_db$hits))) 

# Use right_join to record the reference sequence where assemblers does not support assembly

true_contigs_db <- true_contigs_db |> 
  # select(contig_name, linguistic_complexity_6, reference_coverage, hits, p_good)
  right_join(conoServerDB,by = "hits")

file_out <- file.path(outdir, "transrate_assemblers_true_contigs.rds")

write_rds(true_contigs_db, file = file_out)



