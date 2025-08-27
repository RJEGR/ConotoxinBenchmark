# 
# Calculate baseline entropy avarege per family, diversity of sequences, and cluster similarity
# Step 1, clustering and summarise number of clusters vs number of sequences in sf group
# Step 2, calculate entropy per sequence, and summarise average entropy in sf group

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)


my_custom_theme <- function(...) {
  base_size = 14
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

f <- list.files(path = outdir, pattern = "curated_nuc_conoServerDB.rds", full.names = T)

conoServerDB <- read_rds(f) 


# 1) First split, diversity of sequnces based 1.0 identity 

clusterize <- function(data, seq_type = "...") {
  
  identify_sequence_type <- function(seq_vector) {
    
    require(Biostrings)
    
    # Coerce to a single string for Biostrings functions
    seq_string <- paste(seq_vector, collapse = "")

    
    # Try to create a DNAString object
    is_dna <- tryCatch({
      DNAString(seq_string)
      TRUE
    }, error = function(e) {
      FALSE
    })
    
    # If it's not DNA, try to create an AAString object
    if (!is_dna) {
      is_aa <- tryCatch({
        AAString(seq_string)
        TRUE
      }, error = function(e) {
        FALSE
      })
      
      if (is_aa) {
        return("Amino Acid")
      } else {
        return("Not a valid DNA or Amino Acid sequence")
      }
    } else {
      return("DNA")
    }
  }
  
  seq_subject <- data %>%
    dplyr::rename("subject" = seq_type) %>%
    distinct(subject) %>% 
    pull(subject, name = subject) 
  
  is_dna <- identify_sequence_type(sample(seq_subject,1))
  
  if(is_dna == "DNA") {
    cat("\n is DNA, using DNAStringSet\n ")
    seq_subject <- Biostrings::DNAStringSet(seq_subject)
    
  } else {
    cat("\n is AA, using AAStringSet\n ")
    seq_subject <- Biostrings::AAStringSet(seq_subject)
  }
    
  # 1.1) Using decipher
  # 
  # Use a high similarity cutoff (e.g., 90% or 95%) to create strict subfamilies of very closely related proteins with similar functions.
  # 
  # Use a lower cutoff (e.g., 50% or 60%) to create broad protein families that group more distantly related homologs together.

  clusters <- DECIPHER::Clusterize(seq_subject,
    cutoff = 0.5,  
    minCoverage=0.95, 
    processors=NULL) # use all CPUs

  
  seq_name <- paste0(seq_type, "_", "clus")
  
  clustersdf <- data.frame(
    seq_type = rownames(clusters), 
    seq_name = paste0(seq_type,"_", clusters$cluster)) %>% as_tibble()
  
  names(clustersdf) <- c(seq_type, seq_name)
  
  # clustersdf %>% right_join(data, by = seq_type)
  
  return(clustersdf)
  
  
}

sequence_df <- clusterize(conoServerDB, seq_type = "sequence")

proteinsequence_df <- clusterize(conoServerDB, seq_type = "proteinsequence")

conoServerDB <- conoServerDB %>%
  left_join(sequence_df, by = "sequence") %>%
  left_join(proteinsequence_df, by = "proteinsequence")

conoServerDB %>% dplyr::count(sequence_clus, sort = T)
conoServerDB %>% dplyr::count(proteinsequence_clus, sort = T)




# Based on the number of sequence clusters and the number of members within each cluster, you can estimate parameters related to protein family size, conservation, and diversity. These metrics provide a quantitative measure of the evolutionary and functional landscape of your dataset.

# 
# Key Parameters and Their Interpretation
# Protein Family Size Distribution: By plotting a histogram of the number of members per cluster, you can visualize the distribution of protein family sizes. This helps you identify whether your dataset is dominated by a few large, highly expanded families (e.g., kinases) or by many small, singleton families.
# 

conoServerDB %>% 
  dplyr::count(proteinsequence_clus, sort = T)

# Protein Family Diversity: The total number of clusters is a direct measure of the number of distinct protein families in your dataset. A higher number of clusters relative to the total number of sequences indicates greater diversity.

seq_type <- "sequence" # proteinsequence
seq_name <- paste0(seq_type, "_", "clus")

summarise_df <- conoServerDB %>%  
  dplyr::rename("subject" = seq_type) %>% 
  distinct(genesuperfamily, subject) %>%
  dplyr::count(genesuperfamily, sort = T) %>%
  dplyr::rename_("seq_number" = "n") 

conoServerDB %>% 
  distinct_at(c("genesuperfamily", seq_name)) %>%
  dplyr::count(genesuperfamily, sort = T) %>%
  dplyr::rename("clus_number" = "n") %>%
  left_join(summarise_df, by = "genesuperfamily") %>%
  mutate(genesuperfamily = gsub("superfamily", "", genesuperfamily)) %>%
  mutate(genesuperfamily = gsub("Divergent ", "", genesuperfamily)) %>%
  ggplot(aes(clus_number, seq_number, label = genesuperfamily )) + 
  scale_x_log10() +
  scale_y_log10() +
  geom_point(shape = 1) +
  ggrepel::geom_label_repel(max.overlaps = 100) +
  my_custom_theme() +
  labs(y = "Number of conotoxins (log10)", x = "Number of clusters (log10)", caption = "Sequence_complexity.R", subtitle = seq_type) 

# Continue w/
# write a chunk of code for Simpson's Diversity Index in r, for a given total number of protein families (clusters), number of sequences in the ith protein family, and thetotal number of sequences in your dataset.

calculate_simpson_diversity <- function(family_members) {
  # Calculate the total number of sequences (N).
  # This is the sum of all members across all families.
  N <- sum(family_members)
  
  # Calculate the proportion (p_i) for each family.
  # This is the number of members in a family divided by the total number of sequences.
  proportions <- family_members / N
  
  # Calculate the sum of squared proportions.
  # This is the core component of the Simpson's index.
  sum_of_squares <- sum(proportions^2)
  
  # The Simpson's Diversity Index is typically expressed as 1 - D,
  # so a higher value indicates greater diversity.
  simpson_index <- 1 - sum_of_squares
  
  return(simpson_index)
}

summarise_df %>% 
  # group_by(genesuperfamily) %>%
  reframe(calculate_simpson_diversity(seq_number))

calculate_simpson_diversity(summarise_df$seq_number)

# A value of 

# D=1
# indicates zero diversity (a single family with all sequences).

# A value close to 

# D=0
# indicates very high diversity. For this reason, the index is often expressed as its complement, 

# 1−D
# , where a higher value signifies higher diversity.

# 
# Conservation and Expansion:
#   
#   Large Clusters: Clusters with a high number of members often correspond to protein families that are highly conserved or have undergone significant gene expansion events within the species or lineage you are studying. This suggests they play a crucial biological role.
# 
# Small Clusters/Singletons: Conversely, clusters with only one or a few members may represent lineage-specific genes or proteins that are not widely conserved.


# Correlate 
# data %>%
#   select_if(is.double) %>%
#   mutate_all(~replace(., is.na(.), 0)) %>%
#   cor(method = "spearman") -> M
# 
# 
# testRes <- corrplot::cor.mtest(M, conf.level = 0.95)
# 
# corrplot::corrplot(M, p.mat = testRes$p ,method = "color", type="upper", order = "hclust", insig = "label_sig")
