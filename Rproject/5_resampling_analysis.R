# Utilizar técnicas de muestreo probabilístico, especialmente el muestreo aleatorio, que da igualdad de oportunidades a todos los individuos de la población.

# Emplear muestreo estratificado, dividiendo la población en subgrupos homogéneos y tomando muestras proporcionales para asegurar la representación adecuada de cada grupo.


# scp -r rgomez@omica.cicese.mx:/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/fernando_pub/2_subsampling_dir/2_transrate_contigs_dir .

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

read_transrate_scores <- function(path) {
  
  boostrap <- gsub("boostrap_dir_","", path)
  
  basedir <- gsub("boostrap_dir_","", path)
  
  Superfamily <- sapply(strsplit(basedir, "_"), `[`, 1)
  
  Assembler <- gsub(paste0(Superfamily, 
    "_superfamily[_|.]|all_superfamilies.fixed[_|.]|Conopeptides_superfamilies[_|.]"),"",basedir)
  
  
  
  # f <- list.files(file.path(dir,path, basedir), "contigs.csv", full.names = T)
  
  f <- list.files(list.dirs(file.path(dir,path), recursive = F), "contigs.csv", full.names = T)
  
  
  read_csv(f) %>% mutate(Superfamily = Superfamily, Assembler = Assembler)
  
}

dir <- "/Users/cigom/Documents/GitHub/ConotoxinBenchmark/2_subsampling_dir/2_transrate_contigs_dir/"


basedir <- list.files(dir)[1]

subdirs <- list.files(file.path(dir,basedir),pattern = "_dir")


list.files(path = dir, pattern = "contigs", recursive = TRUE, full.names = TRUE)

# 
# transratedf <- lapply(subdirs, read_transrate_scores)
# 
# transratedf <- do.call(rbind,transratedf)
