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


library(xml2)

library(tidyverse)

encoding <- "UTF-8" 

head(xml_lines <- readLines(con, encoding = encoding, warn = FALSE))

doc <- read_xml(paste(xml_lines, collapse = "\n"),
  encoding = "UTF-8",as_html = T,
  options = c("RECOVER", "NOERROR", "NOBLANKS"))

close(con)

entries <- xml_find_all(doc, ".//entry")

extract_entry_info <- function(entry_node) {
  
  # pb <- txtProgressBar(min = 0, max = n, style = 3) # style=3 is a nice load bar
  # 
  # results <- vector("list", n)
  # for (i in seq_len(n)) {
  #   results[[i]] <- extract_entry_info(entries[[i]])
  #   setTxtProgressBar(pb, i)
  # }
  # close(pb)
  
  # Get id (assuming one per entry)
  entry_id <- xml_text(xml_find_first(entry_node, "./id"))
  
  
  cat(entry_id,"\n")
  
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

# Apply to all entries
all_entry_data <- lapply(entries, extract_entry_info)


# Example: print data for first entry

data <- dplyr::bind_rows(all_entry_data) %>% as_tibble() %>% mutate(sequence = str_to_upper(sequence))

# EDA -----

data %>% drop_na(proteinsequence)
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

data %>% 
  count(organismlatin, genesuperfamily, sort = T) %>%
  mutate(genesuperfamily = factor(genesuperfamily, levels = unique(genesuperfamily))) %>%
  ggplot(aes(organismlatin, genesuperfamily, fill = n)) +
  geom_tile(color = "white", linewidth = 0.5) +
  theme_bw(base_family = "GillSans", base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


# Filter and Split -----

# assume artifitial RNAseq will be 100 SE
nrow(data)

# Considere remove nucleotide without start codon (optional)
# ex.
# data %>% drop_na(proteinsequence) %>% filter(grepl("^ATG", sequence))

nrow(Nodedf <- data %>% filter(!grepl("Patent|patent", name)))

nrow(Nodedf <- Nodedf %>% drop_na(proteinsequence) %>% filter(grepl("^M", proteinsequence)))

nrow(Nodedf <- Nodedf %>% filter(nchar(sequence) >= 100))

nrow(Nodedf %>% distinct(sequence))

# Write fasta using a apply to write 

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
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



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

dedup_DNAStringSet <- function(dnaset) {
  seq_chars <- as.character(dnaset)
  split_names <- split(names(dnaset), seq_chars)
  unique_seqs <- Biostrings::DNAStringSet(names(split_names))
  names(unique_seqs) <- sapply(split_names, paste, collapse="|")
  unique_seqs
}


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

# q <- c("ATGCAGACGGCCTACTGGGTGATGGTGATGATGATGGTGTGGATTGCAGCCCCTCTGTCTGAAGGTGGTAAACTGAACGATGTAATTCGGGGTTTGGTGCCAGACGACATAACCCCACAGCTCATGTTGGGAAGTCTGATTTCCCGTCGTCAATCGGAAGAGGGTGGTTCAAATGCAACCAAGAAACCCTATATTCTAAGGGCCAGCGACCAGGTTGCATCTGGGCCATAG")

# dedup_DNAStringSet(seqs)[as.character(dedup_DNAStringSet(seqs)) %in% q]


