
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

outdir <- "C://Users//cinai/OneDrive/Documentos/GitHub/ConotoxinBenchmark/INPUTS/ConoSorter_dir"

dir.create(outdir)



which_cols <- c("Superfamily (Signal)", "Superfamily (Pro-region)", "Superfamily (Mature)")

paste_col <- function(x) { 
  x <- x[!is.na(x)] 
  x <- unique(sort(x))
  x <- paste(x, sep = '_', collapse = '_')
}

f <- Regex_f[1]

read_regex <- function(f, Hydrophobicity_val = 60, pwidth_val = 50) {
  
  cat("Reading ", basename(f), "\n")
  
  DF <- read_delim(f, delim = "|", col_names = T) %>% mutate(Method = basename(f))
  
  # is_nucleotide == TRUE
  
  DF <- DF %>%
    dplyr::rename("Hydrophobicity"="% Hydrophobicity (Signal)", "ID" = "Read Name") %>%
    mutate(Hydrophobicity = gsub("%", "", Hydrophobicity), Hydrophobicity = as.double(Hydrophobicity)) %>%
    dplyr::rename("Protein_width"="# A.A", "Cys_number" = "# Cysteine(s)") %>%
    dplyr::rename("Score_sf"="Score Superfamily", "Score_class" = "Score Class") %>%
    dplyr::rename("seq"="Protein Sequence") %>%
    mutate(Conflict = Score_sf)
  
  #  ELSE
  
  protein_id <- DF %>% pull(ID)
  protein_id <- gsub("_[0-9]+_[0-9]+$","", protein_id)
  
  DF <- cbind(data.frame(protein_id), DF) %>% mutate(Method = gsub("_Regex.tab", "", Method)) %>% as_tibble()
  
  Conflictdf <- DF %>% filter(grepl("CONFLICT", Conflict)) %>% distinct(protein_id, Conflict)
  
  # Filter step as Borghie et al.
  
  DF <- DF %>% 
    filter(Hydrophobicity > Hydrophobicity_val) %>%
    filter(Protein_width >= pwidth_val) %>%
    # Omit conflicts
    mutate(Score_sf = gsub("[^0-9.-]", "", Score_sf)) %>%
    mutate(Score_class = gsub("[^0-9.-]", "", Score_class))
  
  
  DF1 <- DF %>% select(contains(c("protein_id","Method", "Score_sf","Cys_number","Superfamily ("))) %>%
    # filter(Score_sf > 0) %>% # Score_sf > 0 == Superfamily != "-"
    pivot_longer(cols = all_of(which_cols), names_to = "Region", values_to = "Superfamily") %>%
    filter(Superfamily != "-") %>%
    mutate(Region = gsub("Superfamily ", "", Region)) %>%
    group_by(Method, protein_id, Cys_number) %>%
    summarise(across(Region,.fns = paste_col), Score_sf = n()) %>% ungroup()
  
  
  OUT <- DF %>% 
    select(contains(c("protein_id","Method", "Score_sf","Superfamily ("))) %>%
    # filter(Score_sf > 0) %>% 
    pivot_longer(cols = all_of(which_cols), names_to = "Region", values_to = "Superfamily") %>%
    filter(Superfamily != "-") %>%
    mutate(Region = gsub("Superfamily ", "", Region)) %>%
    mutate(Superfamily =  gsub("\\(.*?\\)", "", Superfamily)) %>% 
    group_by(Method, protein_id) %>%
    summarise(across(Superfamily, .fns = paste_col)) %>% 
    # arrange(desc(Score_sf)) %>%
    left_join(DF1) %>%
    mutate(tab = "Regex")
  
  # OUT <- OUT %>% left_join(Conflictdf)
  
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
    dplyr::rename("E_value"="E-value Superfamily") %>%
    mutate(Conflict = E_value)
  
  protein_id <- DF %>% pull(ID)
  protein_id <- gsub("_[0-9]+_[0-9]+$","", protein_id)
  
  DF <- cbind(data.frame(protein_id), DF) %>% mutate(Method = gsub("_pHMM.tab", "", Method)) %>% as_tibble()
  
  Conflictdf <- DF %>% filter(grepl("CONFLICT", Conflict)) %>% distinct(protein_id, Conflict)
  
  # two steps clean strings
  DF <- DF %>% mutate(E_value =  gsub("!!CONFLICT!!", "", E_value))
  
  DF <- DF %>% mutate(E_value =  gsub("[()]", "", E_value))
  
  DF <- DF %>% mutate(E_value = as.numeric(E_value))
  
  
  
  DF <- DF %>% 
    filter(as.numeric(E_value) < eval) %>% # omit pval filtering
    filter(Hydrophobicity > Hydrophobicity_val) %>%
    filter(Protein_width >= pwidth_val)
  
  # DF %>% 
  #   select(contains(c("protein_id","Method", "Score_sf","Superfamily ("))) %>%
  #   pivot_longer(cols = all_of(which_cols), names_to = "Region", values_to = "Superfamily") %>%
  #   filter(Superfamily != "-") %>%
  #   mutate(Region = gsub("Superfamily ", "", Region)) %>%
  #   group_by(Method, protein_id) %>%
  #   summarise(across(Region,.fns = paste_col), Score_sf = n()) %>%
  #   mutate(tab = "Regex")
  # 
  
  DF1 <- DF %>% select(contains(c("protein_id","Method", "Cys_number","Superfamily ("))) %>%
    # filter(Score_sf > 0) %>% # Score_sf > 0 == Superfamily != "-"
    pivot_longer(cols = all_of(which_cols), names_to = "Region", values_to = "Superfamily") %>%
    filter(Superfamily != "-") %>%
    mutate(Region = gsub("Superfamily ", "", Region)) %>%
    group_by(Method, protein_id, Cys_number) %>%
    summarise(across(Region,.fns = paste_col), Score_sf = n()) %>% ungroup()
  
  
  OUT <- DF %>% 
    select(contains(c("protein_id","Method", "Score_sf","Superfamily ("))) %>%
    pivot_longer(cols = all_of(which_cols), names_to = "Region", values_to = "Superfamily") %>%
    filter(Superfamily != "-") %>%
    mutate(Region = gsub("Superfamily ", "", Region)) %>%
    mutate(Superfamily =  gsub("\\(.*?\\)", "", Superfamily)) %>% 
    group_by(Method, protein_id) %>%
    summarise(across(Superfamily, .fns = paste_col)) %>% 
    left_join(DF1) %>%
    mutate(tab = "pHMM")
  
  # OUT <- OUT %>% left_join(Conflictdf)
  
  return(OUT)
  
}


# read_regex(Regex_f[1])
# 
# read_pHMM(pHHM_f[10])

dir1 <- "//wsl.localhost/Debian/home/ricardo/ConoSorter_outdir/nucleotide_dir/"

dir2 <- "//wsl.localhost/Debian/home/ricardo/ConoSorter_outdir/peptide_dir/"


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


DB1 <- Read_conoSorter(dir1)


DB2 <- Read_conoSorter(dir2)

DB1 %>% write_rds(file = outFile)

# Save bkp of data, 
# summarise number of classification by 
outName <- "Assemblers_1" 

# outName <- "Assemblers_2" # this is for PLASS replicates only 

outFile <- file.path(outdir, paste0(outName, "_ConoSorter_regex_pHMM.rds"))

rbind(Regex_df, pHMM_df) %>% write_rds(file = outFile)


rbind(Regex_df, pHMM_df) %>% count(tab)

# rbind(Regex_df, pHMM_df) |> filter(protein_id %in% "BINPACKER.0.11") 

# annotation_results |> filter(qseqid %in% "BINPACKER.0.11")

DB2 <- rbind(Regex_df, pHMM_df) |> 
  mutate(vfold_set = sapply(strsplit(Method, "_"), `[`, 1)) %>%
  mutate(Method = sapply(strsplit(Method, "_"), `[`, 5)) %>%
  # mutate(file_name = gsub(".[1|2].blast", "", file_name))
  ungroup() |> distinct(protein_id, tab, Method, vfold_set) |> 
  dplyr::rename("qseqid" = protein_id)



extrafont::loadfonts(device = "win")

my_custom_theme <- function(base_size = 14, legend_pos = "top", ...) {
  base_size = 14
  theme_bw(base_family = "Gill Sans MT", base_size = base_size) +
    theme(legend.position = legend_pos,
          strip.placement = "outside", 
          strip.background = element_rect(fill = 'gray90', color = 'white'),
          strip.text = element_text(angle = 0, size = base_size, hjust = 0), 
          axis.text = element_text(size = rel(0.7), color = "black"),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          ...
    )
}



annotation_results |> 
  filter(final_annotation != "error") %>%
  # filter(final_annotation == "chimera") %>%
  mutate(vfold_set = sapply(strsplit(file_name, "_"), `[`, 1)) %>%
  mutate(Method = sapply(strsplit(file_name, "_"), `[`, 5)) %>%
  mutate(Method = gsub(".[1|2].blast", "", Method))|>
  distinct(qseqid, final_annotation, Method, vfold_set) |> 
  left_join(DB2) |> 
  group_by(final_annotation, Method, tab) %>%
  count(sort = T) |> 
  mutate(tab = ifelse(is.na(tab), "Not annotated", tab)) |>
  group_by(final_annotation, tab) %>%
  mutate(f = n/sum(n)) |>
  # rstatix::get_summary_stats(type = "mean_sd")
  ggplot(aes(y = Method, x = final_annotation, fill = f)) +
  geom_tile() +
  geom_text(aes(label = n), color = "white") +
  facet_grid(~ tab) +
  # geom_col() +
  my_custom_theme()


