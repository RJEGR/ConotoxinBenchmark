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

outdir <- "~/Documents/GitHub/ConotoxinBenchmark/INPUTS/"


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
all_entry_data <- lapply(entries[1:10], extract_entry_info)

dplyr::bind_rows(all_entry_data) %>% as_tibble() %>% mutate(sequence = str_to_upper(sequence)) %>% view()

