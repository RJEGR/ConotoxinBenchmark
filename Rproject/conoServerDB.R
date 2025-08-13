# Think before code
# 1. What is th problem I am trying to solve?
# 2. Understand how the data is organized and structured
# 3. Data life-cycle Query, insert, update, modify and model (prediction, stats, )

# To do
# True set data  (artificial or simulated rnaseq data)
# Split conserved nucleotide sequences by taxon (species or genus)
# For every taxon, obtain Worms taxonomy lineage for further purposes (Fernando)

# XML format
# always extract individual noe from each entry (use xml_find_first instead of xml_find_all)
# Do NOT extract all node at once from full document, or yo lose correspondence
# use extract_entry_info to split info from nucleotide and protein nodes


rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

outdir <- "~/Documents/GitHub/ConotoxinBenchmark/INPUTS/"


url <- 'https://www.conoserver.org/download/conoserver_nucleic.xml.gz'

con <- gzcon(url(url, "rb"))

file_out <- gsub(".xml.gz",".rds", basename(url))

dedup_DNAStringSet <- function(dnaset) {
  
  require(Biostrings)
  
  seq_chars <- as.character(dnaset)
  split_names <- split(names(dnaset), seq_chars)
  unique_seqs <- Biostrings::DNAStringSet(names(split_names))
  names(unique_seqs) <- sapply(split_names, paste, collapse="|")
  unique_seqs
}
extract_entry_info <- function(entry_node) {
  
  
  entry_id <- xml_text(xml_find_first(entry_node, "./id"))
  
  
  # cat(entry_id,"\n")
  
  # nucleotide_list <- lapply(entry_node, function(node) {
  #   list(
  #     organismlatin = xml_text(xml_find_first(node, "organismlatin")),
  #     organismdiet = xml_text(xml_find_first(node, "organismdiet")),
  #     entry_id = xml_text(xml_find_first(node, "id")),
  #     name = xml_text(xml_find_first(node, "name")),
  #     sequence = xml_text(xml_find_first(node, "sequence"))
  #   )
  # })
  
  xpath_vec_nuc <- c("organismlatin","organismdiet","organismregion","name","sequence")
  
  nucleotide_list <- lapply(xpath_vec_nuc, function(xpath) { xml_text(xml_find_first(entry_node, xpath)) })
  
  nucleotide_list <- as.data.frame(nucleotide_list, stringsAsFactors = FALSE)
  
  colnames(nucleotide_list) <- xpath_vec_nuc 
  
  # Find all proteins under proteinsencoded
  proteins <- xml_find_all(entry_node, "./proteinsencoded/protein")
  
  xpath_vec <- c("proteinid","proteinname","genesuperfamily","proteinsequence")
  
  if (!length(proteins) == 0) {
    
    protein_list <- lapply(xpath_vec, function(xpath) { xml_text(xml_find_first(proteins, xpath)) })
    
    # protein_list <- lapply(proteins, function(prot) {
    #   list(
    #     proteinid = xml_text(xml_find_first(prot, "./proteinid")),
    #     proteinname = xml_text(xml_find_first(prot, "./proteinname")),
    #     genesuperfamily = xml_text(xml_find_first(prot, "./genesuperfamily")),
    #     proteinsequence = xml_text(xml_find_first(prot, "./proteinsequence"))
    #   )
    # })
    
    # as for some nucleotide there are not protein related:
    
    protein_list <- as.data.frame(protein_list, stringsAsFactors = FALSE)
    
    
    
  } else {
    
    
    # Return empty data frame with the right columns
    protein_list <- as.data.frame(matrix(nrow = 1, ncol = length(xpath_vec)))
    
    
  }
  
  colnames(protein_list) <- xpath_vec
  
  # Return a list with entry id and proteins info
  
  data.frame(
    entry_id,
    nucleotide_list,
    protein_list
  )
}

extract_entry_info_ <- function(entries, ...) {
  
  #   # wrapped version from names2wormsdf_, include ProgressBar
  #   maxq <- length(entries)
  #   
  #   print(maxq)
  #   
  #   pb <- txtProgressBar(min = 0, max = maxq, style = 3, width = 50, char = "=")
  #   
  #   out <- list()
  #   
  #   for(i in 1:maxq) {
  #     
  #     
  #     # out[[i]] <- extract_entry_info(entries)
  #     
  #     setTxtProgressBar(pb, i)
  #   }
  #   
  #   close(pb) 
  #   
  #   out <- do.call(rbind, out)
  #   
  #   return(out)
  #  
}

library(xml2)

library(tidyverse)

encoding <- "UTF-8" 

head(xml_lines <- readLines(con, encoding = encoding, warn = FALSE))

doc <- read_xml(paste(xml_lines, collapse = "\n"),
  encoding = "UTF-8",as_html = T,
  options = c("RECOVER", "NOERROR", "NOBLANKS"))

close(con)

entries <- xml_find_all(doc, ".//entry")


# Apply to all entries
all_entry_data <- lapply(entries, extract_entry_info)


data <- dplyr::bind_rows(all_entry_data) %>% as_tibble() %>% mutate(sequence = str_to_upper(sequence))


recode_to <- c("Conus flavidus" = "vermivorous", 
  "Conus varius" = "vermivorous", 
  "Conus terebra" = "vermivorous", 
  "Conus sulcatus" = "molluscivorous",
  "Conus adamsonii" = "piscivorous",
  "Conus andremenezi" = "vermivorous",
  "Conus araneosus" = "molluscivorous")

data <- data %>% 
  mutate(
    organismdiet = ifelse(
      organismlatin %in% names(recode_to),
      recode_to[organismlatin],
      organismdiet
    )
  ) 

write_rds(data, file = file.path(outdir, file_out)) 

# data <- read_rds(file.path(outdir, file_out))

# EDA -----

data %>% 
  # drop_na(proteinsequence)
  count(organismlatin, sort = T) %>%
  mutate(organismlatin = factor(organismlatin, levels = unique(organismlatin))) %>%
  ggplot(aes(y = organismlatin, x = n)) + geom_col()

data %>% count(sequence, proteinsequence, sort = T) 

data %>% count(sequence, genesuperfamily, sort = T) 

data %>% ggplot(aes(nchar(sequence))) + geom_histogram() 

data %>%
  drop_na(genesuperfamily) %>%
  count(genesuperfamily, sort = T) %>%
  mutate(genesuperfamily = factor(genesuperfamily, levels = unique(genesuperfamily))) %>%
  ggplot(aes(y = genesuperfamily, x = n)) + geom_col()


# data %>% 
#   count(organismlatin, genesuperfamily, sort = T) %>%
#   drop_na(genesuperfamily) %>%
#   mutate(genesuperfamily = factor(genesuperfamily, levels = unique(genesuperfamily))) %>%
#   ggplot(aes(organismlatin, genesuperfamily, fill = n)) +
#   geom_tile(color = "white", linewidth = 0.5) +
#   theme_bw(base_family = "GillSans", base_size = 12) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# Filter and Split -----

# assume artifitial RNAseq will be 100 SE
nrow(data)

# Considere remove nucleotide without start codon (optional)
# ex.
# data %>% drop_na(proteinsequence) %>% filter(grepl("^ATG", sequence))

nrow(Nodedf <- data %>% filter(!grepl("Patent|patent", name)))

nrow(Nodedf <- Nodedf %>% drop_na(proteinsequence) %>% filter(grepl("^M", proteinsequence)))

nrow(Nodedf <- Nodedf %>% filter(!grepl("N", sequence))) # omit ambigous sequences

nrow(Nodedf <- Nodedf %>% filter(nchar(sequence) >= 100))

# nrow(Nodedf <- Nodedf %>% filter(nchar(sequence) < 1000))

nrow(Nodedf %>% distinct(sequence))

Nodedf %>% summarise(mean = mean(nchar(sequence)), sd = sd(nchar(sequence)))

Nodedf %>% ggplot(aes(nchar(sequence))) +geom_histogram()

# Computes the total sample size required to achieve a significant statistical power for ordinal outcomes.

# The total sample size from posamsize indicates how many subjects in total are needed to detect the effect size with the specified power.
# The efficiency compares the ordinal design’s efficiency to that of a continuous response design; values less than 1 indicate somewhat less efficiency.

# The power output from popower indicates the probability of correctly rejecting the null hypothesis for the given sample size.
# The standard error of the log odds ratio gives a sense of the precision of the estimated effect size under the model.


estimate_power <- function(df, strata = stratified_sampling) {
  
  require(Hmisc)
  
  # pipeline from (https://library.virginia.edu/data/articles/power-and-sample-size-calculations-ordered-categorical-data)
  

  reference_prop <- df %>% dplyr::rename("strata" = any_of(strata)) %>% drop_na(strata) %>% group_by(strata) %>% dplyr::count(strata) %>% pull(n, name = strata)
  
  
  reference_prop <- prop.table(reference_prop)
  
  # Estimate cumulative control probabilities to obtain projected cumulative feature probabilities.
  
  # cumref <- cumsum(reference)
  
  
  OR <- (0.85*(1-0.7))/(0.7*(1-0.85))
  
  
  f <- function(or, pc){
    (or * pc)/(1 - pc + or * pc)
  }
  
  experimental <- diff(c(0, sapply(cumsum(reference_prop), f, or = OR)))
  
  probs <- rbind(reference_prop, experimental)
  
  # probs
  
  marg_probs <- apply(probs, 2, mean)
  
  # marg_probs
  
  posamsize_res <- posamsize(p = marg_probs,
    odds.ratio = OR,   fraction = 0.5,
    alpha = 0.05, power = 0.9)
  
  
  
  print(posamsize_res)
  # cat("\nThe total sample size required to achieve an significant statistical power\n")

  
  # The R functions popower and posamsize (in the Hmisc package) compute power and sample size estimates for ordinal responses using the proportional odds model.
  
  n <- round(posamsize_res$n)/2
  
  popower(p = marg_probs, odds.ratio = OR, n1 = n, n2 = n, alpha = 0.05)
  
}



estimate_power(Nodedf, strata = "genesuperfamily")
estimate_power(Nodedf, strata = "organismdiet")
estimate_power(Nodedf, strata = "organismlatin")

# Based on the calculations, genesuperfamily or organismdiet is good enogth to stratified samples

# Instead of split data using bias sample size, ----
# Use V-Fold Cross-Validation method to randomly split the data into V groups of roughly equal size

set.seed(123)

library(rsample)

# Nodedf <- Nodedf %>% mutate(strata = interaction(organismdiet, genesuperfamily, drop = TRUE))


stratified_sampling <- "organismdiet"

folds <- rsample::vfold_cv(Nodedf, v = 12, strata = stratified_sampling)

# estimate and diagnose sample bias within and across your folds. 

# mean(sapply(folds$splits, function(split) nrow(analysis(split))))

# s <- folds$splits[[1]]

summarise_vfolds <- function(s, strata = stratified_sampling) {
  
  analysis(s) %>%
    dplyr::rename("strata" = any_of(strata)) %>%
    dplyr::count(strata) %>%
    mutate(split = labels(s)$id)
}

dplyr::bind_rows(lapply(folds$splits, summarise_vfolds)) %>%
  ggplot(aes(x = n, y = strata)) + geom_boxplot()

save_vfold_fasta <- function(s) {
  
  
  resampled_data <- analysis(s)
  
  seqs <- resampled_data %>%
    pull(sequence, name = entry_id) %>%
    Biostrings::DNAStringSet() %>%
    dedup_DNAStringSet()
  
  cat(length(seqs), " sequences\n")
  
  # Write to FASTA, use group name in file
  fasta_file <- file.path(outdir, paste0(labels(s)$id, ".fasta"))
  
  cat("Into ",basename(fasta_file), "\n")
  
  
  Biostrings::writeXStringSet(seqs, fasta_file)
  
  
}

sapply(folds$splits, save_vfold_fasta)

quit()

# Old version include split of data by superfamily, but it is biased sampling method

# Nodedf %>% filter(is.na(genesuperfamily)) %>% view()

Nodedf <- Nodedf %>% mutate(split_as = genesuperfamily) %>% 
  mutate(split_as = ifelse(is.na(split_as), "Conopeptides", split_as)) %>%
  mutate(split_as = str_replace_all(split_as, " ", "_"))


data %>% mutate(split_as = genesuperfamily) %>% 
  mutate(split_as = ifelse(is.na(split_as), "Conopeptides", split_as)) %>%
  mutate(split_as = str_replace_all(split_as, " ", "_")) %>%
  mutate(facet = "A) Raw") %>%
  rbind(data.frame(Nodedf, facet = "B) Filtered")) %>%
  drop_na(genesuperfamily) %>%
  count(facet, genesuperfamily, sort = T) %>%
  mutate(genesuperfamily = factor(genesuperfamily, levels = unique(genesuperfamily))) %>%
  ggplot(aes(y = genesuperfamily, x = n)) + geom_col() + facet_grid(~ facet) +
  theme_bw(base_family = "GillSans", base_size = 14) + labs(x = "n sequences")

Nodedf %>% 
  count(organismlatin, split_as, sort = T) %>%
  mutate(split_as = factor(split_as, levels = unique(split_as))) %>%
  ggplot(aes(organismlatin, split_as, fill = n)) +
  geom_tile(color = "white", linewidth = 0.5) +
  theme_bw(base_family = "GillSans", base_size = 14) +
  # scale_x_discrete(position = "top") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


Nodedf %>% filter(is.na(organismdiet)) %>% count(organismlatin, sort = T)


Nodedf %>%
  count(organismdiet, split_as, sort = T) %>%
  drop_na(split_as, organismdiet) %>%
  mutate(split_as = factor(split_as, levels = unique(split_as))) %>%
  ggplot(aes(y = organismdiet, x = split_as, fill = n)) +
  geom_tile(color = "white", linewidth = 0.5) +
  theme_bw(base_family = "GillSans", base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


file_out <- file.path(outdir, "curated_nuc_conoServerDB.rds")

write_rds(Nodedf, file = file_out)


# Stop here, and go to rsampling_conoServerDB.R to splintting stratified resampling with randomeness to capture nature variance of data without systemic bias

quit()


# Split apart sf with less than 10 sequences and omiting NA sf

well_represented <- Nodedf %>% 
  count(split_as, sort = T) %>%
  filter(n >= 10) %>%
  pull(split_as)
  
less_represented <- Nodedf %>% 
  count(split_as, sort = T) %>%
  filter(n < 10) %>%
  # drop_na() %>%
  pull(split_as)
  
  
length(c(less_represented, well_represented)) # must match 43 sf + 1 unk sf




Nodedf %>% count(split_as) %>% view()

# For each group, subset, format, and write FASTA
for (group in well_represented) {
  # Filter for the current group
  tmp <- Nodedf %>%
    filter(split_as == group) %>%
    # unite("entry_id", entry_id:name, sep = "|") %>%
    pull(sequence, name = entry_id)
  
  # Make DNAStringSet
  seqs <- Biostrings::DNAStringSet(tmp)
  
  seqs <- dedup_DNAStringSet(seqs)
  
  cat(length(seqs), "\n")
  
  # Write to FASTA, use group name in file
  fasta_file <- file.path(outdir, paste0(group, ".fasta"))
  
  Biostrings::writeXStringSet(seqs, fasta_file)
}



seqs <- Nodedf %>%
  filter(split_as %in% less_represented) %>%
  pull(sequence, name = entry_id) %>%
  Biostrings::DNAStringSet() %>%
  dedup_DNAStringSet()


# Write to FASTA, use group name in file
fasta_file <- file.path(outdir, paste0("UNDER_superfamily", ".fasta"))

Biostrings::writeXStringSet(seqs, fasta_file)

# Output the AA sequences

dedup_StringSet <- function(dnaset) {
  seq_chars <- as.character(dnaset)
  split_names <- split(names(dnaset), seq_chars)
  unique_seqs <- Biostrings::AAStringSet(names(split_names))
  names(unique_seqs) <- sapply(split_names, paste, collapse="|")
  unique_seqs
}

# For each group, subset, format, and write FASTA
for (group in well_represented) {

  # Filter for the current group
  tmp <- Nodedf %>%
    filter(split_as == group) %>%
    # unite("entry_id", entry_id:name, sep = "|") %>%
    pull(proteinsequence, name = entry_id)
  
  # Make DNAStringSet
  seqs <- Biostrings::AAStringSet(tmp)
  
  seqs <- dedup_StringSet(seqs)
  
  cat(length(seqs), "\n")
  
  # Write to FASTA, use group name in file
  fasta_file <- file.path(outdir, paste0(group, ".pep"))
  
  Biostrings::writeXStringSet(seqs, fasta_file)
}

seqs <- Nodedf %>%
  filter(split_as %in% less_represented) %>%
  pull(proteinsequence, name = entry_id) %>%
  Biostrings::AAStringSet() %>%
  dedup_StringSet()


# Write to FASTA, use group name in file
fasta_file <- file.path(outdir, paste0("UNDER_superfamily", ".pep"))
Biostrings::writeXStringSet(seqs, fasta_file)

# Run orthofinder -S diamond -op -f peptide_dir. (peptide)
# orthofinder -d -S diamond -op -f orthofinder_dir (dna)



# source https://github.com/RJEGR/Meiofauna/blob/fe5674cb246021f131d4bfcbf6d3a1744db13e07/R/DB_TAXA2WORMS_PR2.R

Nodedf %>% distinct(organismlatin) %>% pull() -> query

tax <- names2wormsdf(query, accepted = TRUE, marine = TRUE)
