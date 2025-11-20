
# Creates a matrix of kmers per superfamily 
# Using entry_id (single conotoxin), complicates the computational cost, 
# while group conotoxins to superfamily, agglomerates conotoxins in a matrix of k dimensions
# As gene superfamily are unbias sample sizes, include fold grouping, to replicates data

# Clusterng or dim reduction

# Vizualize as MDS, or t-SNE (https://pmc.ncbi.nlm.nih.gov/articles/PMC12047656/)

rm(list = ls())

if(!is.null(dev.list())) dev.off()

library(tidyverse)

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

my_custom_theme <- function(base_size = 14, ...) {
  theme_bw(base_family = "GillSans", base_size = base_size) +
    theme(legend.position = "top",
      strip.placement = "outside",
      strip.background = element_rect(fill = 'white', color = 'white'),
      strip.text = element_text(angle = 0, size = base_size), 
      axis.text = element_text(size = rel(0.7), color = "black"),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.x = element_blank(),
      ...
    )
}

outdir <- "~/Documents/GitHub/ConotoxinBenchmark/INPUTS/"

file_out <- file.path(outdir, "curated_nuc_conoServerDB.rds")

data <- read_rds(file_out)

k <- round(log(591506, base = 4) )

window_vec <- c(k, 19, 21, 25, 33, 49, 55, 73)

kmer_merge <- function(sequences, window_vec) {
  

  # Extract k-mers from a sequence for given k
  kwindows <- function(read, k) {
    
    kmers <- c()
    for (i in 1:(nchar(read) - k + 1)) {
      kmers <- c(kmers, substr(read, i, i + k - 1))
    }
    return(kmers)
  }
  
  # Count k-mers for a list of sequences
  count_kmers <- function(seq_list, k) {
    
    # cat("\n kmers for sequence group : ", names(seq_list), "...")

    
    kmers <- lapply(seq_list, kwindows, k = k)
    
    counts <- table(unlist(kmers))
    
    return(counts)
  }
  
  merged_matrices <- lapply(window_vec, function(win) {
   
     cat("\n kmers size : ", win, "...")
    
    # For each window length, count k-mers per sample group
    count_list <- lapply(sequences, count_kmers, k = win)
    
    # Get all unique k-mers across all groups
    all_kmers <- unique(unlist(lapply(count_list, names)))
    
    # Build count matrix by populating counts per group, missing are zero
    mat <- sapply(count_list, function(cnts) {
      vals <- cnts[all_kmers]
      vals[is.na(vals)] <- 0
      vals
    })
    
    
    
    # Name rows and columns for clarity
    rownames(mat) <- all_kmers
    
    colnames(mat) <- names(sequences)
    
    return(mat)
  })
  
  names(merged_matrices) <- paste0("k", window_vec)
  
  # Filter matrix
  
  
  filter_low_abundance <- function(count_mat, threshold = 5) {
    # Compute total counts per k-mer across all samples
    total_counts <- rowSums(count_mat)
    
    # Retain only k-mers with counts >= threshold
    filtered_mat <- count_mat[total_counts >= threshold, , drop = FALSE]
    
    return(filtered_mat)
  }

  filtered_mat <- lapply(merged_matrices, function(x) {
    filter_low_abundance(x, threshold = 1)
  })
  
  # return(merged_matrices)
  return(filtered_mat)
  
}

# if splitted seqs is by entry_id, it will take a while.therefore, omit split as this

data %>% dplyr::count(genesuperfamily, sort = T) %>% filter(n > 5) %>% tail() 

DATA <- data %>% filter(genesuperfamily %in% "S superfamily")

splitted_seqs <- split(strsplit(DATA$sequence, ";") , DATA$genesuperfamily)

# splitted_seqs <- split(strsplit(data$sequence, ";") , data$entry_id)

splitted_seqs <- lapply(splitted_seqs, unlist)

kmers_combined <- kmer_merge(splitted_seqs, window_vec = 33)

colSums(kmers_combined$k33)
table(rowSums(kmers_combined$k33))

# total nkmers: 2098    
#Freq: 1    2    3    4    6    7    8 
#n kmers: 1794   63  135   90   10    1    5 

# S superfamily 2792 

splitted_seqs <- split(strsplit(DATA$sequence, ";") , DATA$entry_id)

splitted_seqs <- lapply(splitted_seqs, unlist)

kmers_combined <- kmer_merge(splitted_seqs, window_vec = 33)

colSums(kmers_combined$k33)
table(rowSums(kmers_combined$k33))

folds <- rsample::vfold_cv(data, v = 12) # strata = stratified_sampling

# folds <- rsample::bootstraps(data, times = 12)


DATA <- do.call(rbind,
  lapply(folds$splits, 
    function(s) { 
      rsample::analysis(s) %>% 
      rsample::initial_split(prop = 0.3) %>%
      rsample::analysis() %>%
      mutate(fold = labels(s)$id)
    }))


DATA <- DATA %>% 
  mutate(genesuperfamily = ifelse(grepl("Divergent", genesuperfamily), "Divergent", genesuperfamily)) %>%
  mutate(fold = paste0(fold, "__", genesuperfamily))

# data %>% dplyr::count(genesuperfamily, sort = T) 

gene_set <- c(
  "M superfamily", "O1 superfamily",
  "A superfamily", "T superfamily",
  "O2 superfamily", "I2 superfamily",
  "O3 superfamily", "J superfamily",
  "Insulin superfamily", "I1 superfamily",
  "Divergent")

DATA <- DATA %>% 
  filter(genesuperfamily %in% gene_set)

DATA %>% dplyr::count(genesuperfamily, sort = T) 

splitted_seqs <- split(strsplit(DATA$sequence, ";") , DATA$fold)

splitted_seqs <- lapply(splitted_seqs, unlist)

kmers_combined <- kmer_merge(splitted_seqs, window_vec = 33)

# threshold: minimum total count across all samples to keep a k-mer

# filtered_mat <- lapply(kmers_combined, function(x) {
#   filter_low_abundance(x, threshold = 1)
# })


# Assuming kmers_combined is a list of count matrices, one per window/k-mer size
kmers_relative_freq <- lapply(kmers_combined, function(kmers_combined) {
  # Divide each column by its sum to get relative frequencies
  sweep(kmers_combined, 2, colSums(kmers_combined), FUN = "/")
})

# colSums(kmers_combined$k10)
# colSums(filtered_mat$k10)
# colSums(kmers_relative_freq$k10)

cmd <- function(m, name) {
  
  # m[is.na(m)] <- 0

  dist_res <- dist(cor(m, method = "spearman"), method = "euclidean")
  mds_res <- cmdscale(dist_res, eig = TRUE, k = 2)
  as_tibble(mds_res$points, rownames = "var") %>% mutate(split = name)
  
  # PCA <- prcomp(m, center = F, scale. = T)
  
  # data.frame(V1 = PCA$x[,1], V2 = PCA$x[,2], split = name)
  
}

cmd_list <- lapply(names(kmers_relative_freq), 
  
  function(n) {
  
  cat("\n kmers size : ", n, "...")

  cmd(kmers_relative_freq[[n]], name = n)

})

cmd_df <- dplyr::bind_rows(cmd_list)


dataviz <- 
  cmd_df %>%
  # mutate(V2 = log(abs(V2)), V1 = log(abs(V1))) %>%
  # filter(abs(V2) <= 0.3 & abs(V1) <= 0.3) %>%
  separate(var, into = c("var", "genesuperfamily"), sep = "__") 
  # left_join(distinct(DATA, entry_id, genesuperfamily), by = "entry_id")  

mark_circle <- dataviz %>% 
  group_by(genesuperfamily) %>% summarise(V1 = mean(V1), V2 = mean(V2))

discrete_scale <- dataviz %>% distinct(genesuperfamily) %>% pull()

n <- length(discrete_scale)

scale_col <- c(ggsci::pal_startrek()(7), ggsci::pal_cosmic()(n-7))

scale_col <- structure(scale_col, names = sort(discrete_scale))


dataviz %>%
  ggplot(., aes(V1, V2)) +
  # stat_ellipse() +
  ggforce::geom_mark_circle(
    data = mark_circle,
    aes(V1, V2,
      group = genesuperfamily, label = genesuperfamily, 
      fill= genesuperfamily, color = genesuperfamily),
    label.buffer = unit(5, 'mm'),
    label.family = "GillSans",
    label.colour = "inherit",
    con.colour = "inherit",
    # concavity = 4,
    con.cap = 0,
    con.type = "elbow",
    con.size = 1,
    con.border = "none",
    con.linetype = 2,
    expand = unit(5, "mm"),
    alpha = 0.2) +
  geom_point(aes(color= genesuperfamily), alpha = 0.5, shape = 1) +
  my_custom_theme(base_size = 20) +
  scale_color_manual("Gene superfamilies",values = scale_col) +
  scale_fill_manual("Gene superfamilies",values = scale_col) +
  labs(y = "Axis 2", x = "Axis 1", caption = "Multidimensional Scaling") +
  theme( legend.position = 'none') -> p


  
ggsave(p,
  filename = 'conoServer_kmer_assessment_2.png',
  path = outdir, width = 7, height = 7.5, dpi = 1000, device = png)

# exit

s <- folds$splits[[1]]
  
DATA <- rsample::analysis(s) %>% 
  mutate(genesuperfamily = ifelse(grepl("Divergent", genesuperfamily), "Divergent", genesuperfamily)) %>%
  filter(genesuperfamily %in% gene_set) %>%
  rsample::initial_split(prop = 0.7) %>%
  rsample::analysis() %>%
  mutate(fold = labels(s)$id) 


splitted_seqs <- split(strsplit(DATA$sequence, ";") , DATA$genesuperfamily)

splitted_seqs <- lapply(splitted_seqs, unlist)

kmers_combined <- kmer_merge(splitted_seqs, window_vec = c(10, 33, 73))

kmers_relative_freq <- lapply(kmers_combined, function(kmers_combined) {
  # Divide each column by its sum to get relative frequencies
  sweep(kmers_combined, 2, colSums(kmers_combined), FUN = "/")
})

dist_res <- dist(cor(kmers_relative_freq$k10, method = "spearman"), method = "euclidean")

hclust_res <- hclust(dist_res, method='complete')

plot(hclust_res)

pca <- function(m, name) {
  
  
  PCA <- prcomp(t(m), center = F, scale. = F)
  
  tibble(var = rownames(PCA$x), split = name, V1 = PCA$x[,1], V2 = PCA$x[,2])
  
}

# pca_list <- lapply(names(kmers_relative_freq), 
#   
#   function(n) {
#     
#     cat("\n kmers size : ", n, "...")
#     
#     pca(kmers_relative_freq[[n]], name = n)
#     
#   })
# 
# pca_df <- dplyr::bind_rows(pca_list)
# 
# pca_df %>%
#   separate(var, into = c("var", "genesuperfamily"), sep = "-") %>%
#   # left_join(data, by = "entry_id") %>%
#   ggplot(., aes(V1, V2, color= genesuperfamily)) +
#   geom_point() + 
#   facet_wrap(~ split) +
#   # ggrepel::geom_text_repel(aes(label = split), max.overlaps = 100) +
#   my_custom_theme() +
#   ggsci::scale_color_atlassian()

# Evals: The pro- portion of k-mers corresponding to a unique position in the ge- nome increases with a greater k 


kmer_df <- lapply(names(kmers_combined)[c(1,3)], function(n) {
  
  cat("\n kmers size : ", n, "...")
  
  filtered_mat <- lapply(kmers_combined, function(x) {
    filter_low_abundance(x, threshold = 0)
  })
  
  
  
  OUT <- list()  # initialize list before loop
  
  for (i in 1:ncol(filtered_mat[[n]])) {
    
    var <- colnames(filtered_mat[[n]])[i]
    
    cat("\n set : ", var, "...")
    
    total_counts <- filtered_mat[[n]][, i]
    
    OUT[[i]] <- tibble(
      kmer = names(total_counts),
      count = as.numeric(total_counts),
      k = n,
      var = var
    )
  }
  
  return(OUT)
})

kmer_df <- dplyr::bind_rows(kmer_df)

# kmer_df %>% count(count)

kmer_df %>%
  # filter(count == 1) %>%
  count(k, count, var) %>%
  group_by(k, count) %>% summarise(sd = sd(n), mean = mean(n))


# e.g. k = 10 there are  49444 k-mers unique to a single region in the conotoxin dataset, for k = 19,there are 76597 unique k-mers, .... The proportion of the whole-dataset that is represented by unique k-mers remains very small for all values of k =< 19,but becomes more substantial for k => 33, and for k => 49 even dominant, which allows forestimating k-mer coverage and downstream genome modelling.

DATA %>%
  mutate(width = nchar(sequence)) %>%
  summarise(width = sum(width))

kmer_df %>%
  filter(count == 1) %>%
  count(k, count, var) %>%
  # filter(count < 10) %>%
  separate(var, into = c("var", "genesuperfamily"), sep = "-") %>%
  ggplot(aes(x = k, y = n, color = genesuperfamily, group = genesuperfamily)) + 
  # facet_grid(genesuperfamily ~., scales = "free_y") +
  geom_jitter(position = position_jitter(0.1), shape = 1, size = 0.5, alpha = 0.7) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun.data=mean_sdl, geom="pointrange", size = 0.5, alpha = 0.5) +
  labs(y = "Proportion of conotoxins represented\nby unique k-mers") +
  my_custom_theme() +
  ggsci::scale_color_atlassian()


# 
# 
# count_kmers <- function(sequences, window_vec) {
#   
#   require(Biostrings)
#   
#   sequences <- DNAStringSet(sequences)
#   
#   kmer_counts_lists <- lapply(sequences, oligonucleotideFrequency, width = window_vec)
#   
#   return(kmer_counts_lists)
# }


head(colSums(kmer_freq))

# total_counts <- kmer_counts[,1]

library(ggplot2)

# Create a data frame for plotting
kmer_df <- data.frame(
  kmer = names(total_counts),
  count = as.numeric(total_counts)
)

ggplot(kmer_df, aes(x = count)) +
  # geom_density()
geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  scale_x_log10() +
  labs(title = paste("K-mer spectrum for k =", k),
    x = "K-mer count (log scale)",
    y = "Frequency")



kmer_counts %>% 
  as_tibble(rownames = "kmer") %>%
  pivot_longer(-kmer) %>%
  ggplot(aes(x = value)) +
  facet_wrap(~ name, scales = "free_y") +
  geom_histogram(binwidth = 1, color = "black") +
  scale_x_log10() +
  labs(title = paste("K-mer spectrum for k =", k),
    x = "K-mer count (log scale)",
    y = "Frequency")

ggplot(kmer_df, aes(x = count)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  scale_x_log10() +
  labs(title = paste("K-mer spectrum for k =", k),
    x = "K-mer count (log scale)",
    y = "Frequency")




PCAdf <- PCA(t(dist_mat))

percentVar1 <- paste0(PCAdf[[2]]$percentVar[1], "% ,",PCAdf[[1]]$percentVar[1])
percentVar1 <- paste0("PC1, VarExp: ", percentVar1, "%")


percentVar2 <- paste0(PCAdf[[2]]$percentVar[2], "% ,",PCAdf[[1]]$percentVar[2])
percentVar2 <- paste0("PC2, VarExp: ", percentVar2, "%")


PCAdf <- PCAdf[[1]]


PCAdf %>%
  mutate(genesuperfamily = rownames(.)) %>% 
  # left_join(Manifest, by = "LIBRARY_ID") %>%
  ggplot(., aes(PC1, PC2)) +
  theme_bw(base_size = 14, base_family = "GillSans") +
  xlab(percentVar1) + ylab(percentVar2) +
  geom_abline(slope = 0, intercept = 0, linetype="dashed", alpha=0.5) +
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.5) +
  # ggrepel::geom_text_repel( family = "GillSans", mapping = aes(label = HPF), size = 3.5, color = "black") +
  # geom_text( family = "GillSans", mapping = aes(label = LIBRARY_ID), size = 3.5, color = "black") +
  theme( legend.position = 'top',
    legend.key.width = unit(0.2, "mm"),
    legend.key.height = unit(0.2, "mm")) -> p

p + geom_point() + 
  ggrepel::geom_text_repel(aes(label = genesuperfamily), max.overlaps = 100) +
  my_custom_theme()


