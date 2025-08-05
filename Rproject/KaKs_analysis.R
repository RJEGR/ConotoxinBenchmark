# Complexity analysis
# To determine heterogeneous domain architecture complexity of a given superfamily, information content and entropy measures were computed using XXX, and compared with the transrate complexity measure. 


rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)
library(Biostrings)
library(seqinr)
library(msa)


outdir <- "~/Documents/GitHub/ConotoxinBenchmark/INPUTS/"

f <- list.files(path = outdir, pattern = "curated_nuc_conoServerDB.rds", full.names = T)

data <- read_rds(f)

# produces an alignment of nucleic protein-coding sequences,

nucl.file <- file.path(outdir, paste0("reverse.align.nuc.tmp.fasta"))
protaln.file <- file.path(outdir, paste0("reverse.align.pep.tmp.fasta"))

data %>%
  # distinct(entry_id, sequence) %>%
  # slice_head(n = 20) %>%
  pull(sequence, name = entry_id) %>%
  Biostrings::DNAStringSet() %>%
  Biostrings::writeXStringSet(nucl.file)

msa_align <- data %>%
  # distinct(entry_id, sequence) %>%
  # slice_head(n = 50) %>%
  pull(proteinsequence, name = entry_id) %>%
  Biostrings::AAStringSet() %>%
  msa::msa(method = "ClustalW") # %>%
  # msa::msaConvert()$seq
  # write_lines(protaln.file)
  # Biostrings::writeXStringSet(protaln.file)


# creates An object of class alignment, 
msa_align <- msaConvert(msa_align)

# seqinr::write.fasta(sequences = msa_align$seq, names = msa_align$nam, file.out = protaln.file)

alignment <- seqinr::as.alignment(nb = msa_align$nb, nam = msa_align$nam, seq = msa_align$seq)
 
seqinr::write.fasta(sequences = msa_align$seq, names = msa_align$nam, file.out = protaln.file)


# 
myOutFileName <- file.path(outdir, paste0("reverse.align.fasta"))
#   
reverse.align(
  nucl.file = nucl.file, protaln.file = protaln.file,
  # align.prot = T,
  input.format = 'fasta', out.file = myOutFileName)

alignment <- seqinr::read.(myOutFileName)

# check In kaks(alignment) : sequence lengths are not a multiple of 3

res <- kaks(alignment)

# 
# devtools::install_github("kullrich/CRBHits", build_vignettes = FALSE, dependencies = TRUE)
# CRBHits::make_last()
# CRBHits::make_dagchainer()
# 
# data(hiv)
# 
# hiv_kaks.Li <- dnastring2kaks(
#   cds=hiv,
#   model="Li")
# g.Li <- plot_kaks(hiv_kaks.Li)
# g.Li