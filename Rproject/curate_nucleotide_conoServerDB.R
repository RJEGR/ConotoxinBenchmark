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


# url <- 'https://www.conoserver.org/download/conoserver_structure.xml.gz'
# Open a connection to the gzipped file at the URL
# con <- gzcon(url(url, "rb"))


# Install required packages if needed
# install.packages("xml2")
# Load required packages

library(xml2)

dir <- "~/"

f <- ""

doc <- read_xml(con,encoding = "")

doc <- read_html(con)

# Read the XML content as text
xml_content <- readLines(con, warn = FALSE)

close(con)

# Parse the XML content
doc <- read_xml(paste(xml_content, collapse = "\n"))

# Extract all 'geneSuperfamily' nodes
geneSuperfamily_nodes <- xml_find_all(doc, ".//geneSuperfamily")

# Extract their text values into a vector
geneSuperfamily_vec <- xml_text(geneSuperfamily_nodes)

# Print the vector
print(geneSuperfamily_vec)