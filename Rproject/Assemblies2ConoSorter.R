
# Ho: test whether or not nucleotide level inflate identification of putative conotoxins, instead of orf level
# nucleotide level
# scp -r rgomez@omica.cicese.mx:/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/fernando_pub/1_assembly_dir/Folds_200x_dir/FASTA_DIR/ConoSorter_outdir/nucleotide_dir
# # Processs .tab and count number of true and false posive annotated to putative conotoxins, 
# Split sequences annotated to the classification group of 80%, 90%, 95% and 100% identity
# 
# orf level
# scp -r rgomez@omica.cicese.mx:/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/fernando_pub/1_assembly_dir/Folds_200x_dir/FASTA_DIR/transdecoder_dir/ConoSorter_outdir/

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

# outdir <- "C://Users//cinai/OneDrive/Documentos/GitHub/ConotoxinBenchmark/INPUTS/ConoSorter_dir"

outdir <- "/Users/rjegr/Documents/Windows/Documents/ConotoxinBenchmark/INPUTS/ConoSorter_dir/"


dir.create(outdir)



which_cols <- c("Superfamily (Signal)", "Superfamily (Pro-region)", "Superfamily (Mature)")

paste_col <- function(x) { 
  x <- x[!is.na(x)] 
  x <- unique(sort(x))
  x <- paste(x, sep = '_', collapse = '_')
}

# f <- Regex_f[1]

read_regex <- function(f, Hydrophobicity_val = 60, pwidth_val = 50) {

  cat("Reading ", basename(f), "\n")
  
  DF <- read_delim(f, delim = "|", col_names = T) %>% mutate(Method = basename(f))
  
  # is_nucleotide == TRUE
  
  DF <- DF %>%
    dplyr::rename("Hydrophobicity"="% Hydrophobicity (Signal)", "ID" = "Read Name") %>%
    mutate(Hydrophobicity = gsub("%", "", Hydrophobicity), Hydrophobicity = as.double(Hydrophobicity)) %>%
    dplyr::rename("Protein_width"="# A.A", "Cys_number" = "# Cysteine(s)") %>%
    mutate(Cys_number = as.numeric(as.character(Cys_number))) %>%
    dplyr::rename("Score_sf"="Score Superfamily", "Score_class" = "Score Class") %>%
    dplyr::rename("seq"="Protein Sequence", "seq_freq" = "Sequence Frequency") %>%
    mutate(Conflict = Score_sf)
  
  #  ELSE
  
  protein_id <- DF %>% pull(ID)
  protein_id <- gsub("_[0-9]+_[0-9]+$","", protein_id)
  
  DF <- cbind(data.frame(protein_id), DF) %>% mutate(Method = gsub("_Regex.tab", "", Method)) %>% as_tibble()
  
  Conflictdf <- DF %>% filter(grepl("CONFLICT", Conflict)) %>% distinct(protein_id, Conflict)
  
  # Filter step as Borghie et al.
  
  DF <- DF %>% 
    # Omit filtering step to preserve artefacts
    # filter(Hydrophobicity > Hydrophobicity_val) %>%
    # filter(Protein_width >= pwidth_val) %>%
    # Omit conflicts
    mutate(Score_sf = gsub("[^0-9.-]", "", Score_sf)) %>%
    mutate(Score_class = gsub("[^0-9.-]", "", Score_class))
  
  
  DF1 <- DF %>% 
    # select(contains(c("protein_id","Method", "Score_sf","Cys_number","Superfamily ("), "Hydrophobicity")) %>%
    filter(Score_sf > 0) %>%
    pivot_longer(cols = all_of(which_cols), names_to = "Region", values_to = "Superfamily") %>%
    
    # Filter all transcripts not annotated by conoSorter, Score_sf > 0 == Superfamily != "-"
    
    filter(Superfamily != "-") %>%
    mutate(Region = gsub("Superfamily ", "", Region)) %>%
    group_by(Method, protein_id, ID, Cys_number, Hydrophobicity, seq_freq, Protein_width) %>%
    summarise(across(Region,.fns = paste_col), Score_sf = n(), .groups = "drop_last") %>% ungroup()
  
  
  OUT <- DF %>% 
    # select(contains(c("protein_id","Method", "Score_sf","Superfamily ("), "ID")) %>%
    filter(Score_sf > 0) %>%
    pivot_longer(cols = all_of(which_cols), names_to = "Region", values_to = "Superfamily") %>%
    filter(Superfamily != "-") %>%
    mutate(Region = gsub("Superfamily ", "", Region)) %>%
    mutate(Superfamily =  gsub("\\(.*?\\)", "", Superfamily)) %>% 
    group_by(Method, protein_id, ID) %>%
    summarise(across(Superfamily, .fns = paste_col), .groups = "drop_last") %>% 
    # arrange(desc(Score_sf)) %>%
    left_join(DF1) %>%
    mutate(tab = "Regex")
  
  OUT <- OUT %>% left_join(Conflictdf)

  return(OUT)
  
  
}

# f <- pHHM_f[1]

read_pHMM <- function(f,  Hydrophobicity_val = 60, pwidth_val = 50, eval = 0.05) {
  
  cat("Reading ", basename(f), "\n")
  
  DF <- read_delim(f, delim = "|", col_names = T) %>% mutate(Method = basename(f))
  
  # DF %>% select(contains(c("E-value", "Score","Biais"))) # how to include ??
  
  DF <- DF %>%
    dplyr::rename("Hydrophobicity"="% Hydrophobicity (Signal)", "ID" = "Read Name") %>%
    mutate(Hydrophobicity = gsub("%", "", Hydrophobicity), Hydrophobicity = as.double(Hydrophobicity)) %>%
    dplyr::rename("Protein_width"="# A.A", "Cys_number" = "# Cysteine(s)") %>%
    mutate(Cys_number = as.numeric(as.character(Cys_number))) %>%
    dplyr::rename("E_value"="E-value Superfamily", "seq_freq" = "Sequence Frequency") %>%
    mutate(Conflict = gsub("[0-9()]|.e-", "", E_value)) 
  
  protein_id <- DF %>% pull(ID)
  protein_id <- gsub("_[0-9]+_[0-9]+$","", protein_id)
  
  DF <- cbind(data.frame(protein_id), DF) %>% mutate(Method = gsub("_pHMM.tab", "", Method)) %>% as_tibble()
  
  Conflictdf <- DF %>% filter(grepl("CONFLICT", Conflict)) %>% distinct(protein_id, Conflict)
  
  # two steps clean strings
  DF <- DF %>% mutate(E_value =  gsub("!!CONFLICT!!", "", E_value))
  
  DF <- DF %>% mutate(E_value =  gsub("[()]", "", E_value))
  
  DF <- DF %>% mutate(E_value = as.numeric(E_value))
  
  
  
  # DF <- DF %>%
  #   filter(as.numeric(E_value) < eval) %>%
  #   filter(Hydrophobicity > Hydrophobicity_val) %>%
  #   filter(Protein_width >= pwidth_val)

  
  DF1 <- DF %>% 
    # select(contains(c("protein_id","Method", "Cys_number","Superfamily ("))) %>%
    pivot_longer(cols = all_of(which_cols), names_to = "Region", values_to = "Superfamily") %>%
    filter(Superfamily != "-") %>%
    mutate(Region = gsub("Superfamily ", "", Region)) %>%
    # group_by(Method, protein_id, Cys_number) %>%
    group_by(Method, protein_id, ID, Cys_number, Hydrophobicity, seq_freq, Protein_width) %>%
    summarise(across(Region,.fns = paste_col), Score_sf = n(), .groups = "drop_last") %>% ungroup()
  
  
  OUT <- DF %>% 
    # select(contains(c("protein_id","Method", "Score_sf","Superfamily ("))) %>%
    pivot_longer(cols = all_of(which_cols), names_to = "Region", values_to = "Superfamily") %>%
    filter(Superfamily != "-") %>%
    mutate(Region = gsub("Superfamily ", "", Region)) %>%
    mutate(Superfamily =  gsub("\\(.*?\\)", "", Superfamily)) %>% 
    # group_by(Method, protein_id) %>%
    group_by(Method, protein_id, ID) %>%
    summarise(across(Superfamily, .fns = paste_col), .groups = "drop_last") %>% 
    left_join(DF1) %>%
    mutate(tab = "pHMM")
  
  OUT <- OUT %>% left_join(Conflictdf)
  
  return(OUT)
  
}


# read_regex(Regex_f[1])
# 
# read_pHMM(pHHM_f[10])

# dir1 <- "//wsl.localhost/Debian/home/ricardo/ConoSorter_outdir/nucleotide_dir/"
# dir2 <- "//wsl.localhost/Debian/home/ricardo/ConoSorter_outdir/peptide_dir/"

dir1 <- "/Users/rjegr/Documents/Windows/Debian/ConoSorter_outdir/nucleotide_dir/"
dir2 <- "/Users/rjegr/Documents/Windows/Debian/ConoSorter_outdir/peptide_dir//"


Read_conoSorter <- function(dir) {

  pattern <- "pHMM.tab"
  
  pHHM_f <- list.files(dir, pattern = pattern, full.names = T) # _pHMM.tab and _Regex.tab
  
  pattern <- "Regex.tab"
  
  Regex_f <- list.files(dir, pattern = pattern, full.names = T) # _pHMM.tab and _Regex.tab
  
  Regex_f <- Regex_f[!grepl("MERGEPIPE|PLASS", Regex_f)] # Huge sizes because many contigs in PLASS
  
  pHHM_f <- pHHM_f[!grepl("MERGEPIPE|PLASS", pHHM_f)] # Huge sizes because many contigs in PLASS
  
  cat("starting w REGEX\n")
  
  Regex_df <- do.call(rbind, lapply(Regex_f, read_regex))
  
  cat("Continue w pHMM\n")
  
  pHMM_df <- do.call(rbind, lapply(pHHM_f, read_pHMM))
  
  rbind(Regex_df, pHMM_df)
  
  
  
}

View(rbind(
  read_regex(f = file.path(dir2, "Fold12_200x_PE_samples_BINPACKER.fa.transdecoder_Regex.tab")),
  read_pHMM(f = file.path(dir2, "Fold12_200x_PE_samples_BINPACKER.fa.transdecoder_pHMM.tab"))))

DB1 <- Read_conoSorter(dir1)

# DB1 |> ungroup() |> distinct(Method, protein_id, Superfamily, Region, tab)

outName <- "Nucleotide_1" 

outFile <- file.path(outdir, paste0(outName, "_ConoSorter.rds"))

DB1 |> ungroup()|> write_rds(outFile)

rm(DB1);gc()

DB2 <- Read_conoSorter(dir2)

outName <- "protein_1" 

outFile <- file.path(outdir, paste0(outName, "_ConoSorter.rds"))

DB2 |> write_rds(outFile)

rm(DB2);gc()


# Script 2
# 
# outName <- "Assemblers_2" # this is for PLASS replicates only 

# outFile <- file.path(outdir, paste0(outName, "_ConoSorter_regex_pHMM.rds"))


