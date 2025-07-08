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

# Note: the xml structure in protein is little different from nucleotide, identify node structure first


rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

outdir <- "~/Documents/GitHub/ConotoxinBenchmark/hhsuite_model_dir/"

url <- 'https://www.conoserver.org/download/conoserver_protein.xml.gz'

con <- gzcon(url(url, "rb"))


library(xml2)

library(tidyverse)

encoding <- "UTF-8" 

head(xml_lines <- readLines(con, encoding = encoding, warn = FALSE))

doc <- read_xml(paste(xml_lines, collapse = "\n"),
  encoding = "UTF-8",as_html = T,
  options = c("RECOVER", "NOERROR", "NOBLANKS"))

close(con)

entries <- xml_find_all(doc, ".//entry")

entry_node <- entries[1]

xpath_vec <- xml_name(xml2::xml_children(entry_node))

protein_list <- lapply(xpath_vec, function(xpath) { xml_text(xml_find_first(entry_node, xpath)) })


extract_entry_info <- function(entry_node) {
  
  entry_id <- xml_text(xml_find_first(entry_node, "./id"))
  
  
  cat(entry_id,"\n")
  
  xpath_vec <- xml_name(xml2::xml_children(entry_node))
  
  # xpath_vec_nuc <- c("organismlatin","organismdiet","organismregion","name","sequence")
  
  protein_list <- lapply(xpath_vec, function(xpath) { xml_text(xml_find_first(entry_node, xpath)) })
  
  protein_list <- as.data.frame(protein_list, stringsAsFactors = FALSE)
  
  colnames(protein_list) <- xpath_vec
  
  # Find all proteins under proteinsencoded
  # proteins <- xml_find_all(entry_node, "./proteinsencoded/protein")
  
  # xpath_vec <- c("proteinid","proteinname","genesuperfamily","proteinsequence")
  
  # if (!length(proteins) == 0) {
    
    # protein_list <- lapply(xpath_vec, function(xpath) { xml_text(xml_find_first(proteins, xpath)) })
    
    # as for some nucleotide there are not protein related:
    
    # nucleotide_list <- as.data.frame(protein_list, stringsAsFactors = FALSE)
    
    
    
  # } else {

    # protein_list <- as.data.frame(matrix(nrow = 1, ncol = length(xpath_vec)))
    
    
  #}
  
  # colnames(protein_list) <- xpath_vec
  
  # Return a list with entry id and proteins info
  
  data.frame(
    # entry_id,
    # nucleotide_list,
    protein_list
  )
}

dplyr::bind_rows(
  extract_entry_info(entries[10])) %>%
  as_tibble() 

# Apply to all entries
all_entry_data <- lapply(entries, extract_entry_info)

data <- dplyr::bind_rows(all_entry_data) %>% 
  as_tibble() %>% 
  mutate(sequence = str_to_upper(sequence)) 


data %>% dplyr::count(class, sort = T)

# nrow(Nodedf <- data %>% filter(class == "conotoxin"))

nrow(Nodedf <- data %>% filter(!grepl("Patent|patent", name)))

nrow(Nodedf <- Nodedf %>% filter(!grepl("unclassified", class)))

nrow(Nodedf <- Nodedf %>% drop_na(sequence) %>% filter(grepl("^M", sequence)))

nrow(Nodedf %>% distinct(sequence))

Nodedf %>% count(genesuperfamily) 

write_rds(data, )

# Write fasta using a apply to write 

Nodedf <- Nodedf %>% mutate(split_as = genesuperfamily) %>% 
  # mutate(split_as = ifelse(is.na(split_as), "Conopeptides", split_as)) %>%
  mutate(split_as = str_replace_all(split_as, " ", "_"))


dedup_StringSet <- function(dnaset) {
  seq_chars <- as.character(dnaset)
  split_names <- split(names(dnaset), seq_chars)
  unique_seqs <- Biostrings::AAStringSet(names(split_names))
  names(unique_seqs) <- sapply(split_names, paste, collapse="|")
  unique_seqs
}


# Nodedf %>% filter(is.na(genesuperfamily)) %>% count(class)

Nodedf <- Nodedf %>% 
  mutate(split_as = ifelse(is.na(split_as), 
    paste("Other",class, sep = "_"), 
    split_as))

Nodedf %>% count(class)

well_represented <- Nodedf %>% 
  count(split_as, sort = T) %>%
  pull(split_as)

for (group in well_represented) {
  # Filter for the current group
  tmp <- Nodedf %>%
    filter(split_as == group) %>%
    # unite("entry_id", entry_id:name, sep = "|") %>%
    pull(sequence, name = id)
  
  # Make DNAStringSet
  seqs <- Biostrings::AAStringSet(tmp)
  
  seqs <- dedup_StringSet(seqs)
  
  cat(length(seqs), "\n")
  
  # Write to FASTA, use group name in file
  fasta_file <- file.path(outdir, paste0(group, ".fasta"))
  
  Biostrings::writeXStringSet(seqs, fasta_file)
}


# creates a subset of nodes for californicus ====

nrow(calNodedf <- data %>% filter(grepl("californicus", organismlatin)))
nrow(calNodedf <- calNodedf %>% drop_na(genesuperfamily))
# nrow(filter(grepl("^M", sequence)))
calNodedf


calNodedf %>% count(class)

nrow(calNodedf %>% distinct(sequence))


calNodedf <- calNodedf %>% mutate(split_as = genesuperfamily) %>% 
  # mutate(split_as = ifelse(is.na(split_as), "Conopeptides", split_as)) %>%
  mutate(split_as = str_replace_all(split_as, " ", "_")) %>%
  mutate(split_as = gsub("-$","",split_as)) %>%
  mutate(cysteineframewrok = ifelse(is.na(cysteineframewrok), "cf", cysteineframewrok))

calNodedf %>% count(split_as, sort = T) 

well_represented <- calNodedf %>% 
  count(split_as, sort = T) %>%
  pull(split_as)

for (group in well_represented) {
  # Filter for the current group
  tmp <- calNodedf %>%
    filter(split_as == group) %>%
    select(sequence, id, name, class, genesuperfamily, cysteineframewrok) %>%
    unite("id", id:cysteineframewrok, sep = "|") %>%
    pull(sequence, name = id)
  
  # Make DNAStringSet
  seqs <- Biostrings::AAStringSet(tmp)
  
  seqs <- dedup_StringSet(seqs)
  
  cat(length(seqs), "\n")
  
  # Write to FASTA, use group name in file
  fasta_file <- file.path(outdir, paste0(group, "_californicus.fasta"))
  
  Biostrings::writeXStringSet(seqs, fasta_file)
}

