
# Read transrate scores (contigs.txt) using read_transrate_scores() (to evals accuracy of assembly methods)
# Output transrateDB containing follow columns
# cd /Users/cigom/Documents/GitHub/ConotoxinBenchmark/1_assembly
# scp -r rgomez@omica.cicese.mx:/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/fernando_pub/1_assembly_dir/Folds_200x_dir/transrate_contigs_dir .

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

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
  
  rel_path <- sapply(strsplit(dirname(rel_path), "/"), `[`, 1)
  vfold_set <- sapply(strsplit(rel_path, "_"), `[`, 1)
  Assembler <- sapply(strsplit(rel_path, "_"), `[`, 5)
  
  cat("\n")
  cat(vfold_set, Assembler)
  cat("\n")

  read_csv(file_list) %>%
    # mutate(file_list = file_list)
    mutate(rel_path, vfold_set, Assembler)
}

calculate_metrics <- function(df, reference_coverage_val = 1) {
  
  is_chimeric_value = 0
  
  calculate_false <- function(df) {
    
    # as many assemblers use width 200 to filter contigs, count number of refseq > 200
    
    count_Nsequences <- function() {
      
      # Temporal directory where vfolds_resampling_dir are found
      
      outdir <- "~/Documents/GitHub/ConotoxinBenchmark/INPUTS/vfolds_resampling_dir/"
      
      f <- list.files(path = outdir, pattern = ".fasta", full.names = T)
      
      my_func <- function(x) { 
        
        dna <- Biostrings::readDNAStringSet(x)
        
        structure(
          sum(Biostrings::width(dna) >=200), 
          names = gsub(".fasta", "", basename(x)))
        
        
      }
      
      unlist(lapply(f, my_func))
      
    }
    
    
    InputNsequences <- count_Nsequences()
    
    # False Negatives (FN): Transcripts present in the simulated data but not assembled. 
    ## TP - (N reference sequences in InputNsequences) 
    
    
    data.frame(InputNsequences) %>% 
      as_tibble(rownames = "vfold_set") %>% 
      right_join(df) %>%
      mutate(FN = abs(InputNsequences - TP))
    
  }
  
  
  Totaldf <- df %>% 
    dplyr::count() %>%
    dplyr::rename("rawcontigs" = "n")
  
  # True Positives (TP): Transcripts correctly assembled by the assembler. 
  ## TP = reference_cov >= reference_coverage_val ||reference_cov == 1 (= !is.na(hits))
  
  TP <- df %>%
    filter(reference_coverage >= reference_coverage_val) %>%
    # Important: make distinct hits, to count only unique references (Recall), not over-inflate by N contigs!!!
    distinct(hits) %>% 
    dplyr::count() %>%
    dplyr::rename("TP" = "n")
  
  # False Positives (FP): Transcripts incorrectly assembled by the assembler.
  ## FP = N contig_name where reference_cov < 1 BUT have some identity threshold (reference_coverage > 0)
  
  FP <- df %>%
    mutate(reference_coverage = ifelse(is.na(hits) & is.na(reference_coverage), 0, reference_coverage )) %>%
    filter(reference_coverage > is_chimeric_value &  reference_coverage < reference_coverage_val)  %>%
    dplyr::count() %>%
    dplyr::rename("FP" = "n")
  
  
  # Dealing with Overestimate contig number
  # No hay manera de usar estos como TN,
  # contigs where reference_coverage == 0, OR
  # Use Rawcontigs number minus TP + FP to count potential True Negative (TN)
  #     mutate(TN = rawcontigs - (TP+FP)) %>%
  
  TN <- df %>%
    mutate(reference_coverage = ifelse(is.na(hits) & is.na(reference_coverage), 0, reference_coverage )) %>%
    # mutate_all(~replace(., is.na(.), 0)) %>%
    filter(reference_coverage <= is_chimeric_value) %>%
    # filter(is.na(hits)) %>%
    # dplyr::count(vfold_set, sampling_set) %>%
    dplyr::count() %>%
    dplyr::rename("TN" = "n")
  
  
  Totaldf %>% 
    left_join(TP) %>% 
    # left_join(TN) %>% 
    left_join(FP) %>% 
    mutate_all(~replace(., is.na(.), 0)) %>%
    left_join(calculate_false(.)) %>%
    select(-InputNsequences) 
  
  
}


outdir <- "~/Documents/GitHub/ConotoxinBenchmark/INPUTS/"

f <- list.files(path = outdir, pattern = "curated_nuc_conoServerDB.rds", full.names = T)

conoServerDB <- read_rds(f) %>% dplyr::rename("hits" = "entry_id")

dir <- "/Users/cigom/Documents/GitHub/ConotoxinBenchmark/1_assembly/transrate_contigs_reference_dir/"

str(file_list <- list.files(path = dir, pattern = "contigs.csv", recursive = T, full.names = TRUE))


transratedf <- lapply(file_list, read_transrate_scores)

transratedf <- do.call(rbind,transratedf)

transratedf %>%
  group_by(Assembler) %>%
  dplyr::count()  

transratedf %>%
  group_by(Assembler) %>%
  summarise(mean = mean(length), min = min(length), max = max(length))
  # ggplot(aes(length)) +
  # facet_grid(Assembler ~., scales = "free", space = "free") +
  # geom_histogram()

# transratedf %>% 
#   group_by(vfold_set, Assembler) %>%
#   calculate_metrics(reference_coverage_val = 0.95) #%>%
  # write_tsv(file = file.path(dir, "benchmark.tsv"))

metricsdf <- transratedf %>% 
  filter(length>=200) %>%
  group_by(vfold_set, Assembler) %>%
  calculate_metrics(reference_coverage_val = 0.95) %>% 
  mutate(
    Ratio = TP/FP,
    # Tell us what percentage of positive classes were correctly identified
    Accuracy = TP / (TP + FN + FP),
    Precision = TP /(TP + FP),
    Sensitivity = TP /(TP + FN),
    Fscore = 2 * (TP) / (2 * (TP) + FP + FN), 
  ) 


cols_to <- c("Accuracy", "Precision", "Sensitivity")

recode_to <- structure(c("A) Accuracy", "B) Precision", "C) Sensitivity"), names = cols_to)

base_size <- 14

line_data <- metricsdf %>% filter(Assembler == "TRINITY") %>%
  select(-TP, -FP, -FN, -rawcontigs) %>%
  pivot_longer(cols = all_of(cols_to), values_to = "x", names_to = "facet") %>%
  dplyr::mutate(facet = dplyr::recode_factor(facet, !!!recode_to)) %>%
  drop_na() %>% group_by(Assembler, facet) %>% summarise(x = mean(x), label = "Baseline")

p1 <- metricsdf %>%
  select(-TP, -FP, -FN, -rawcontigs) %>%
  pivot_longer(cols = all_of(cols_to), values_to = "x", names_to = "facet") %>%
  dplyr::mutate(facet = dplyr::recode_factor(facet, !!!recode_to)) %>%
  drop_na() %>%
  ggplot(aes(x = x, y = Assembler)) +
  # ggforce::facet_col(~ facet) +
  facet_grid(~ facet) + 
  geom_vline(data = line_data, aes(xintercept = x), color = "gray90", linetype = "dashed") +
  geom_text(data = line_data, aes(x = x, label = label),  hjust = -0.3, color = "gray70", angle = 360, size = 3) +
  geom_jitter(position = position_jitter(0.1), shape = 1) +
  stat_summary(fun = "mean", geom = "point", aes(group = Assembler), color="blue") +
  stat_summary(fun.data=mean_sdl, geom="pointrange", aes(group = Assembler), color="blue") +
  labs(x = "Parameter", y = "", caption = "Assemblers.R") +
  # ylim(0,NA) +
  my_custom_theme()


p1

ggsave(p1, filename = 'Assemblers_.png', 
  path = outdir, width = 7, height = 5, dpi = 1000, device = png)

# Figure 2

metricsdf

discrete_scale <- metricsdf %>% ungroup() %>% distinct(Assembler) %>% pull()

n <- length(discrete_scale)

scale_col <- c(ggsci::pal_startrek()(7), ggsci::pal_cosmic()(n-7))

scale_col <- structure(scale_col, names = sort(discrete_scale))

metricsdf_flt <- metricsdf %>% 
  group_by(Assembler) %>% 
  # filter-out outliers
  mutate(Sensitivity = ifelse(rstatix::is_outlier(Sensitivity), NA, Sensitivity)) %>%
  mutate(Precision = ifelse(rstatix::is_outlier(Precision), NA, Precision)) %>%
  drop_na(Sensitivity, Precision)
  
mark_circle <- metricsdf_flt %>%
  summarise(
    sens_sd = sd(Sensitivity), pre_sd = sd(Precision), n = n(),
    Sensitivity = mean(Sensitivity), Precision = mean(Precision))

line_data <- metricsdf_flt %>% 
  filter(Assembler == "TRINITY") %>%
  summarise(Sensitivity = max(Sensitivity), Precision = max(Precision), label = "Baseline")


metricsdf %>%
  group_by(Assembler) %>% 
  mutate(Sensitivity = ifelse(rstatix::is_outlier(Sensitivity), NA, Sensitivity)) %>%
  mutate(Precision = ifelse(rstatix::is_outlier(Precision), NA, Precision)) %>%
  drop_na(Sensitivity, Precision) %>%
  mutate(facet = "A) Benchmark metrics") %>%
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

ggsave(p, filename = 'Assemblers_Precision_Sensitivity.png', 
  path = outdir, width = 4.5, height = 4.5, dpi = 1000, device = png)




# quit()


p2 <- metricsdf %>% group_by(Assembler) %>% 
  summarise(sd = sd(Accuracy), x = mean(Accuracy)) %>%
  arrange(desc(x)) %>%
  mutate(label = paste0(Assembler, " (",  round(x, 3), ")")) %>%
  mutate(Assembler = factor(Assembler, levels = rev(unique(Assembler)))) %>%
  mutate(facet = "B) Accuracy") %>%
  mutate(xmin = x-sd, xmax = x+sd, xlab = xmax+0.1) %>%
  mutate(xlab = ifelse(Assembler == "STRINGTIE", 0.25, xlab)) %>% 
  ggplot() + 
  facet_grid(~ facet) +
  geom_col(aes(y = Assembler, x = x, fill= Assembler, color = Assembler), width = 0.7) + 
  geom_text(aes(y = Assembler, x = xlab, label = label), hjust = 0.1, size = 3) +
  geom_errorbar(aes(y = Assembler, x = x, xmin = xmin, xmax = xmax), width = 0.15, alpha = 0.3, color = "gray70") +
  scale_x_continuous("",limits = c(0,1)) +
  scale_y_discrete(position = "right")+
  scale_color_manual("", values = scale_col) +
  scale_fill_manual("", values = scale_col) +
  my_custom_theme(legend_pos = "none", axis.ticks.y = element_blank(), 
    axis.text.y = element_blank(), axis.title.y = element_blank())  

p2

metricsdf %>% group_by(Assembler) %>% 
  summarise(sd = sd(Ratio), x = mean(Ratio)) %>%
  arrange(desc(x)) %>%
  mutate(Assembler = factor(Assembler, levels = unique(Assembler))) %>% #View()
  ggplot(aes(y = Assembler, x = x)) + 
  geom_col(aes(fill= Assembler, color = Assembler), width = 0.7) + 
  geom_text(aes(label = round(x, 3)), hjust = -0.15, size = 3.5) +
  scale_x_continuous("Correct/incorrect ratio (C/I)",limits = c(0,12)) +
  scale_color_manual("", values = scale_col) +
  scale_fill_manual("", values = scale_col) +
  my_custom_theme(legend_pos = "none")

library(patchwork)
# p + p2

# Inner join of conotoxins, and assemblers,

library(ggVennDiagram)

upsetdf <- transratedf %>% 
  filter(!is.na(hits)) %>%
  filter(reference_coverage > 0.95) %>%
  distinct(Assembler, hits) #%>%
# group_by(hits) %>%
# summarise(across(Assembler, .fns = list), n = n()) 

gene2ven <- split(upsetdf$hits, upsetdf$Assembler)

gene2ven <- lapply(gene2ven, unlist)

keep <- names(gene2ven) %in% c("SPADES", "TRINITY", "PLASS", "STRINGTIE", "RNABLOOM")

ggVennDiagram(gene2ven[keep],label_font = "GillSans", label_size = 5,
  relative_height = 0.5,relative_width = 0.8, force_upset = F) +
  scale_fill_gradient(low="grey90",high = "red")

ggVennDiagram(gene2ven,label_font = "GillSans", label_size = 7,
  relative_height = 1,relative_width = 2, force_upset = T) 
# scale_fill_gradient(low="grey90",high = "red",)


# file_out <- file.path(outdir, "transrate_assemblers.rds")

# write_rds(transratedf, file = file_out)


# According to Cahis et al., 2012
# By assembler, assess/Classify hits, to fragmented, chimeric, allelic, paralogue, and other genomic, based on overlap -----


dir <- "/Users/cigom/Documents/GitHub/ConotoxinBenchmark/1_assembly/blast_outputs/"

str(f <- list.files(path = dir, pattern = "1.blast", recursive = T, full.names = TRUE))

blast_df <- read_outfmt6(f[1])

# blast_df <- do.call(rbind, lapply(f, read_outfmt6))

# Hits were con- sidered significant when (i) the alignment length (merg- ing all high-scoring segments pairs) was at least 70% of the query sequence or at least 70% of the hit sequence, and (ii) sequence identity between query and hit was more than 70% across the aligned portion.

# Fragmented: 
# Contigs with a single significant hit shared by other contigs were called (allele)

blast_df %>% 
  filter(identity >= 95) %>%
  group_by(db) %>%
  dplyr::count(target, sort = T) %>%
  mutate()
# dplyr::rename("allele" = "n") 

# (Quimera or multi) Contigs with several significant hits, all specific to this contig, were called



summarise_df <- transratedf %>% 
  filter(!is.na(hits)) %>%
  filter(reference_coverage > 0.5) %>%
  mutate(summarise = "< 80 % alignment") %>%
  mutate(summarise = ifelse(reference_coverage >= 0.8, ">= 80% alignment", summarise)) %>%
  mutate(summarise = ifelse(reference_coverage >= 0.9, ">= 90% alignment", summarise)) %>%
  mutate(summarise = ifelse(reference_coverage >= 0.95, ">= 95% alignment", summarise)) %>%
  mutate(summarise = ifelse(reference_coverage == 1, "100% alignment", summarise)) %>%
  group_by(summarise, Assembler, vfold_set) 

## TRansrate (if read-based metrics are found)

names(transratedf) %in% c("in_bridges", "p_good", "")

# The contig will represent a single transcript, s(Cseg). This score mea- sures the probability that the coverage depth of the transcript is univariate, i.e., that it represents an assembly of a single tran- script and not a hybrid/chimeric assembly of multiple tran- scripts expressed at different expression levels. Here, the per- nucleotide coverage depth of the contig must be best modeled by a single Dirichlet distribution (described below). If the contig is better modeled by the product of two or more Dirichlet distri- butions, then this indicates that two or more contigs with differ- ent transcript abundances have been erroneously assembled together.

#  Similarly, transcripts with low s(Cseg) scores are likely to represent chimeric transcripts. Here, al- though the transcript itself may be incorrectly assembled, the com- ponent segments of the transcript may themselves be correctly assembled and of utility if separated. To help users identify and diagnose likely assembly errors affecting low scoring contigs, TransRate provides each of the separate contig scores (in addition to the overall contig score). This

transratedf %>% 
  filter(!is.na(hits)) %>%
  # filter(reference_coverage < 0.5) %>%
  # ggplot(aes(in_bridges)) + geom_histogram()
  # ggplot(aes(p_good, score, alpha = reference_coverage)) + geom_point()
  ggplot(aes(p_not_segmented, tpm)) + geom_point()

# Correlate 
transratedf %>% 
  filter(!is.na(hits)) %>% 
  select_if(is.double) %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  cor(method = "spearman") -> M


testRes <- corrplot::cor.mtest(M, conf.level = 0.95)
# 
corrplot::corrplot(M, p.mat = testRes$p ,method = "color", type="upper", order = "hclust", insig = "label_sig")



# transratedf %>% 

true_contigs_db <- transratedf %>% filter(!is.na(hits)) 

true_contigs_db %>% distinct(hits) %>% nrow()

# all_contigs_db <- transratedf %>% filter(is.na(hits)) 

sum(sort(conoServerDB$hits) %in% sort(unique(true_contigs_db$hits))) 

# Use right_join to record the reference sequence where assemblers does not support assembly

true_contigs_db <- true_contigs_db %>% 
  # select(contig_name, linguistic_complexity_6, reference_coverage, hits, p_good)
  right_join(conoServerDB,by = "hits")

file_out <- file.path(outdir, "transrate_assemblers_true_contigs.rds")

write_rds(true_contigs_db, file = file_out)



