# Think before code
# 1. What is th problem I am trying to solve?
# 2. Understand how the data is organized and structured
# 3. Data life-cycle Query, insert, update, modify and model (prediction, stats, )

# To do
# True set data  (artificial or simulated rnaseq data)
# Split conserved nucleotide sequences by taxon (species or genus)
# For every taxon, obtain Worms taxonomy lineage for further purposes (Fernando)


rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)


# Nucleotide

url <- 'https://www.conoserver.org/download/conoserver_nucleic.xml.gz'

# Open a connection to the gzipped file at the URL
con <- gzcon(url(url, "rb"))


# Install required packages if needed
# install.packages("xml2")
# Load required packages

library(xml2)

library(tidyverse)

# if as file:

outdir <- "~/Documents/GitHub/ConotoxinBenchmark/INPUTS/"

# f <- list.files(path = dir, pattern = "xml.gz", full.names = T)
# Open a gzipped connection and read the lines

# con <- gzfile(f[2], open = "rt") # 'rt' for reading text

xml_lines <- readLines(con, warn = FALSE)

close(con)

# Collapse lines and parse as XML
# (Read as html to support special character as greek letters)

doc <- read_html(paste(xml_lines, collapse = "\n")) # options = c("NOBLANKS", "NSCLEAN")


all_nodes <- xml_find_all(doc, "//*",)

# Extract all node names
node_names <- xml_name(all_nodes)

# Count occurrences of each node name
table(node_names)

# Extract fields features for nucleotide level
# must match 3073 records nrow(Nodedf %>% drop_na(id))

xpath_vec_nuc <- c(
  "proteinid", "proteinname","genesuperfamily","proteintype","proteinsequence",
  "organismlatin","organismdiet","organismregion","id","name","sequence")

result <- sapply(xpath_vec_nuc, function(xpath) { xml_text(xml_find_first(all_nodes, xpath)) })

result <- as.data.frame(result, stringsAsFactors = FALSE) %>% as_tibble()

colnames(result) <- xpath_vec_nuc # xpath_vec_prot

# Creates a fasta file and input
# Filtering criteria
# choosing only nucleotides of precursor type and proteinsequence

Nodedf <- result %>% drop_na(id) 

Nodedf %>% 
  count(organismlatin) %>% filter(n > 10) %>% view()

Nodedf <- Nodedf %>% 
  count(sequence, id, organismlatin, name, sort = T) %>%
  mutate(sequence = str_to_upper(sequence))

Nodedf %>% ggplot(aes(nchar(sequence))) + geom_histogram()

# Write fasta using a apply to write 

Nodedf <- Nodedf %>% mutate(split_as = organismlatin)

v <- Nodedf %>% count(split_as, sort = T) %>% pull(split_as)

v <- tail(v)

# For each group, subset, format, and write FASTA
for (group in v) {
  # Filter for the current group
  tmp <- Nodedf %>%
    filter(split_as == group) %>%
    unite("id", id:name, sep = "|") %>%
    pull(sequence, name = id)
  
  # Make DNAStringSet
  seqs <- Biostrings::DNAStringSet(tmp)
  
  # Write to FASTA, use group name in file
  fasta_file <- file.path(outdir, paste0("conopeptides_", group, ".fasta"))
  Biostrings::writeXStringSet(seqs, fasta_file)
}

seqs <- Nodedf %>% 
  unite("id", id:name, sep = "|") %>%
  pull(sequence, name = id)

seqs <- Biostrings::DNAStringSet(seqs)


Biostrings::writeXStringSet(seqs, file.path(pub_dir, "conopeptides.fasta"))

# 
# gene2GO <- split(strsplit(gene2GO$GOs, ";") , gene2GO$gene_id)
# gene2GO <- lapply(gene2GO, unlist)


# PRotein =====
# url <- 'https://www.conoserver.org/download/conoserver_protein.xml.gz'
# head(xml_lines[grepl("nucleicAcidId", xml_lines) ], n = 100)

# N00316N00540N02519N02520
# xml_lines[grepl("N00940|N00316|N00540|N02519|N02520", xml_lines) ]

# head(xml_lines[grepl("N00316|", xml_lines) ], n = 100)

# OR Extract features for protein level

# conoserver identifier|name|organism|protein type|toxin class|gene superfamily|cysteine framework|pharmacological
# must match 8523 records nrow(Nodedf %>% drop_na(id))

xpath_vec_prot <- c(
  "class", "cysteineframewrok","genesuperfamily","pharmacologicalfamily",
  "nucleicacid", "nucleicacidid",
  "organismlatin","organismdiet","organismregion","id","name","sequence")




