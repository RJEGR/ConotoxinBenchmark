# Ricardo Gomez-Reyes
# Visualize Nx 
# Calculate nx distribution
# Include sizes after transrating


rm(list = ls())

if(!is.null(dev.list())) dev.off()

library(tidyverse, help, pos = 2, lib.loc = NULL)

outdir <- "~/Documents/GitHub/ConotoxinBenchmark/INPUTS/"

f <- list.files(path = outdir, pattern = "curated_nuc_conoServerDB.rds", full.names = T)

conoServerDB <- read_rds(f) %>% dplyr::rename("hits" = "entry_id")

dna_ref <- list.files(path = outdir, pattern = "conoServerDB.fasta", full.names = T)

dir <- "/Users/cigom/Documents/GitHub/ConotoxinBenchmark/1_assembly/FASTA_DIR/"


require(Biostrings)
require(dplyr)
require(ggplot2)

f <- list.files(dir, "fa$", full.names = T)

metrics_df <- function(f, stringSet = "DNA") {
  
  contig_Nx <- function(f, stringSet = "DNA") {
    
    if(stringSet != "DNA") {
      stringSet <- Biostrings::readAAStringSet(f)
    } else
      
      stringSet <- Biostrings::readDNAStringSet(f)
    
    contig_width <- sort(Biostrings::width(stringSet), decreasing = T)
    
    # contig_df <- data.frame(width = contig_width, Assembly = basename(f))
    
    return(contig_width)
  }
  
  cat("\nReading\n")
  cat(basename(f),"\n")
  cat("\n")
  
  # width <- sort(seq(1,100), decreasing = TRUE)
  
  width <- contig_Nx(f, stringSet = stringSet)
  
  assembly_seqs <- length(width)
  
  v <- seq(0.1,1, by = 0.1)
  
  Lx <- function(x) { sum(cumsum(width) < (sum(width) * x)) + 1 }
  
  l <- unlist( lapply(v, Lx))
  
  # add the number of sequences per Nx (field assembly_seqs)
  
  n_seqs <- sum(width > width[l[1]]) + 2
  
  # Cut the widths into chunks based on the Nx breakpoints 
  breakpts <- l # unique(l)
  
  # chunks <- cut(seq_along(width), breaks = l, labels = FALSE)
  
  chunks <- cut(rev(seq_along(width)), breaks = breakpts, labels = FALSE) 
  
  # Split the widths into chunks based on the cuts 
  chunks <- split(width, chunks) 
  
  n_seqs <- c(n_seqs,   as.vector(unlist(lapply(chunks, length))))
  
  n_frac <- n_seqs/assembly_seqs
  
  metrics_df <- data.frame(x = paste0("N", v*100), n = width[l], l = l, 
    n_seqs, n_frac, Assembly = basename(f))
  
  
  metrics_df <- as_tibble(metrics_df)
  
  return(metrics_df)
  
}

contig_Nx <- function(f, stringSet = "DNA") {
  
  if(stringSet != "DNA") {
    stringSet <- Biostrings::readAAStringSet(f)
  } else
    
    stringSet <- Biostrings::readDNAStringSet(f)
  
  contig_width <- sort(Biostrings::width(stringSet), decreasing = T)
  
  contig_df <- data.frame(width = contig_width, Assembly = basename(f))
  
  return(contig_df)
}


my_custom_theme <- function(legend_pos = "top", ...) {
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

f <- f[!grepl("RNABLOOM.fa", f)] # as rnablom resulted in some fasta files < 10 contigs, drops


df <- lapply(f, metrics_df)
df <- do.call(rbind, df)
df <- mutate(df, x = factor(x, levels = unique(df$x)))


df <- df %>% 
  mutate(vfold_set = sapply(strsplit(Assembly, "_"), `[`, 1)) %>% 
  mutate(Assembly = sapply(strsplit(Assembly, "_"), `[`, 5)) %>%
  rbind(metrics_df(dna_ref)) %>%
  mutate(Assembly = gsub(".fa$|.fasta$","", Assembly)) 


# 
# recode_to <- c(
#   "spades_hisat_superDuper.fasta",
#   "spades.fasta",
#   "Merged_polyA_hisat_SuperDuper.fasta",
#   "Merged_clusters.fasta", # MMseqs.fasta
#   "trinity_hisat_superDuper.fasta",
#   "Trinity.fasta")
# 
# recode_to <- structure(
#   c(
#     "Full-length Spades (Lace)", 
#     "Spades (S)", 
#     "Full-length T-S (Lace)", 
#     "Full-length T-S (Mmseq)",
#     "Full-length Trinity (Lace)",
#     "Trinity (T)"), 
#   names = recode_to)


# df <- mutate(df, Assembly = dplyr::recode_factor(Assembly, !!!recode_to))

# write_rds(df, file = paste0(pub_dir, "/Contig_nx.rds"))

# rnsps <- mean(contig_Nx(f[[1]]))
# trnt <- mean(contig_Nx(f[[2]]))

# df %>%  ggplot(aes(x = x, y = Assembler)) +


discrete_scale <- unique(df$Assembly)

n_pallet <- length(discrete_scale)

# scale_col <- ggsci::pal_uchicago()(n_pallet) 
scale_col <- c(ggsci::pal_startrek()(7), ggsci::pal_cosmic()(n_pallet-7))

scale_col <- structure(scale_col, names = sort(discrete_scale))

# scale_fill <- ggsci::pal_uchicago(alpha = 0.5)(n_pallet) 

scale_fill <- c(ggsci::pal_startrek(alpha = 0.5)(7), ggsci::pal_cosmic(alpha = 0.5)(n_pallet-7))

scale_fill <- structure(scale_fill, names = sort(discrete_scale))


p <- df %>% 
  ggplot(aes(x = x, y = n, color = Assembly)) + 
  # ggforce::facet_col(~ Assembly) +
  # geom_jitter(position = position_jitter(0.1), shape = 1) +
  stat_summary(fun = "mean", geom = "line", aes(group = Assembly)) +
  # stat_summary(fun = "mean", geom = "point", aes(group = Assembly)) +
  geom_vline(xintercept = "N50", linetype="dashed", alpha=0.5) +
  labs(x = "Nx", y = "Contig length", color = "Assembly method") +
  my_custom_theme(legend_pos = "none") +
  # guides(color = guide_legend(title = "", nrow = 3, ncol = 4)) +
  scale_color_manual("",values = scale_col ) 

p
# annotate("text", y = rnsps, x = "N50", angle = 90, label = "label")

summarised_df <- df %>% 
  group_by(Assembly, x) %>%
  summarise(n = mean(n)) %>%
  mutate(label = paste0(Assembly, " (", scales::comma(n), ")"))

p <- p +
  geom_point(data = filter(summarised_df, x == "N50"), shape = 21, size = 2) +
  ggrepel::geom_text_repel(data = filter(summarised_df, x == "N50"),
    aes(label= label), max.overlaps = 50,
    nudge_x      = 1,
    # nudge_y = 1,
    direction    = "y",
    hjust        = 0,
    vjust = -15,
    segment.curvature = 0.1,
    # segment.ncp = 3,
    # segment.angle = 20,
    size = 3, family = "GillSans") 



ggsave(p, filename = 'Nx-methods.png', path = outdir, width = 5, height = 5, device = png, dpi = 300)

p2 <- 
  summarised_df %>% 
  ggplot(aes(x = Assembly, y = n, group = Assembly, fill = Assembly)) +
  ggforce::facet_col(~ x, scales = "free_y") +
  ggplot2::geom_col(position = position_dodge2(reverse = T)) +
  labs(x = "Assembler", y = "Number of contigs", fill = "Assembly method") +
  scale_color_manual("",values = scale_col) +
  scale_fill_manual("",values = scale_fill) +
  # guides(fill=guide_legend(title = "", ncol = 1)) +
  # facet_wrap(~  stringSet, scales = "free_y") + 
  scale_x_discrete(position = "top") +
  my_custom_theme(legend_pos = "none")



p2


width_df <- do.call(rbind, lapply(f, contig_Nx))

width_df <- width_df %>%
  mutate(vfold_set = sapply(strsplit(Assembly, "_"), `[`, 1)) %>%
  mutate(Assembly = sapply(strsplit(Assembly, "_"), `[`, 5)) %>%
  rbind(contig_Nx(dna_ref)) %>%
  mutate(Assembly = gsub(".fa$|.fasta$","", Assembly))

vline_df <- width_df %>%
  group_by(Assembly) %>%
  summarise(mean = mean(width), min = min(width), max = max(width))


xintercept <- 105
  
p3 <- width_df %>%
  mutate(facet = ifelse(Assembly == "conoServerDB", "A) ConoServer", "B) Assemblers")) %>%
  ggplot(aes(x = width, color = Assembly, fill = Assembly)) +
  xlim(0,1000) +
  ggforce::facet_col(~ facet, scales = "free_y", space = "fixed") +
  # ggplot2::stat_ecdf(linewidth = 3) +
  geom_histogram(stat = "density", position = position_stack()) +
  geom_vline(data = vline_df, aes(xintercept = min, color = Assembly), 
    linetype="dashed", alpha=0.5) +
  my_custom_theme(legend_pos = "top") +
  scale_color_manual("",values = scale_col) +
  scale_fill_manual("",values = scale_fill)

ggsave(p3, filename = 'width_distributions.png', path = outdir, width = 7, height = 5, device = png, dpi = 300)



# p2 + guides(fill = guide_legend(title = "", nrow = 2, ncol = 4))

library(patchwork)

p2 <- p / p2


ggsave(p2, filename = 'Nx-methods-2.png', path = dir, width = 7.5, height = 10, device = png, dpi = 300)

df %>% group_by(Assembly) %>% summarise(n_seqs = sum(n_seqs), n_frac = sum(n_frac)) %>% arrange(desc(n_seqs))
