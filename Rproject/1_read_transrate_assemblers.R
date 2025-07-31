
# Read transrate scores (contigs.txt) using read_transrate_scores() (to evals accuracy of assembly methods)
# Output transrateDB containing follow columns
# 

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

read_transrate_scores <- function(path) {
  
  basedir <- gsub("_transrate_dir","",path)
  
  Superfamily <- sapply(strsplit(basedir, "_"), `[`, 1)
  Assembler <- gsub(paste0(Superfamily, 
    "_superfamily[_|.]|all_superfamilies.fixed[_|.]|Conopeptides_superfamilies[_|.]"),"",basedir)
  
  
  
  # f <- list.files(file.path(dir,path, basedir), "contigs.csv", full.names = T)
  
  f <- list.files(list.dirs(file.path(dir,path), recursive = F), "contigs.csv", full.names = T)
  
  
  read_csv(f) %>% mutate(Superfamily = Superfamily, Assembler = Assembler)
  
}

outdir <- "~/Documents/GitHub/ConotoxinBenchmark/INPUTS/"

f <- list.files(path = outdir, pattern = "curated_nuc_conoServerDB.rds", full.names = T)

conoServerDB <- read_rds(f) %>% dplyr::rename("hits" = "entry_id")

dir <- "/Users/cigom/Documents/GitHub/ConotoxinBenchmark/2_transrate_dir/2_transrate_contigs_dir"

str(subdirs <- list.files(dir, pattern = "_transrate_dir"))

# read_transrate_scores(subdirs[2])

transratedf <- lapply(subdirs, read_transrate_scores)

transratedf <- do.call(rbind,transratedf)

transratedf %>%  count(Assembler)

transratedf %>% distinct(hits) %>% arrange(desc(hits))


conoServerDB %>% distinct(hits) %>% arrange(desc(hits))
# Because some issues with all superfamily ID,
transratedf %>% filter(nchar(hits) > 6) %>% count(Superfamily)


transratedf <- transratedf %>%  
  mutate(hits = sapply(strsplit(hits, "_"), `[`, 1)) 
  # mutate(hits = strsplit(hits, "|")) %>%
  # unnest(hits) %>% 

file_out <- file.path(outdir, "transrate_assemblers.rds")

write_rds(transratedf, file = file_out)


true_contigs_db <- transratedf %>% filter(!is.na(hits)) 

true_contigs_db %>% distinct(hits) %>% nrow()

# all_contigs_db <- transratedf %>% filter(is.na(hits)) 

sum(sort(conoServerDB$hits) %in% sort(unique(true_contigs_db$hits))) 

# Use right_join to record the reference sequence where assemblers does not support assembly

true_contigs_db <- true_contigs_db %>% 
  # select(contig_name, linguistic_complexity_6, reference_coverage, hits, p_good)
  right_join(conoServerDB,by = "hits")

file_out <- file.path(outdir, "transrate_assemblers_true_contigs.rds")

write_rds(true_contigs_db, file = file_out)



