
# According to Cahis et al., 2012
# By assembler (OR better by kmer), annotate hits, to fragmented, chimeric, allelic, paralogue, and other genomic, based on overlaps -----

# cd /Users/cigom/Documents/GitHub/ConotoxinBenchmark/1_assembly

# scp -r rgomez@omica.cicese.mx:/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/fernando_pub/1_assembly_dir/Folds_200x_dir/FASTA_DIR/transrate_contigs_dir/blast_outputs .

# scp -r rgomez@omica.cicese.mx:/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/fernando_pub/3_kmer_dir/multiple_analysis_dir/transrate_contigs_dir/blast_outputs

# 1. merge to best reciprocal, 
# 2. run blast annotation

# Load necessary packages
library(dplyr)
library(tidyr)
library(tidyverse)

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

my_custom_theme <- function(...) {
  base_size = 14
  theme_bw(base_family = "GillSans", base_size = base_size) +
    theme(legend.position = "top",
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

# dir <- "/Users/cigom/Documents/GitHub/ConotoxinBenchmark/1_assembly/blast_outputs/"

dir <- "/Users/cigom/Documents/GitHub/ConotoxinBenchmark/3_kmer_dir/transrate_contigs_dir/blast_outputs//" 
str(file_list_1 <- list.files(path = dir, pattern = "1.blast", recursive = T, full.names = TRUE))

str(file_list_2 <- list.files(path = dir, pattern = "2.blast", recursive = T, full.names = TRUE))

file_list <- file_list_2

file_list <- file_list[!grepl("MERGEPIPE|PLASS", file_list)] # Huge sizes because many contigs asmb

nlines <- function(f) {
  cat("\nReading\n")
  cat(basename(f))
  cat("\n")
  
  cat("\nContaining hits:\n")
  system(paste0("wc -l ", f, "| awk '{print $1}'"))
  cat("\n")
}

# lapply(file_list, nlines)

blast_annotation <- function(f) {
  
  require(tidyverse)
  
  cat("\nReading\n")
  cat(f)
  cat("\n")
  
  cat("\nContaining hits:\n")
  system(paste0("wc -l ", f, "| awk '{print $1}'"))
  cat("\n")
  
  
  read_outfmt6 <- function(f) {
    
    # seqid = transcript_id
    outfmt6.names <- c("qseqid", "sseqid", "pident", "length", "mismatches", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen")
    
    
    # qseqid (query sequence ID)
    # sseqid (subject sequence ID)
    # pident (percentage of identical matches)
    # length (alignment length)
    # mismatch (number of mismatches)
    # gapopen (number of gap openings)
    # qstart (start of the alignment in the query sequence)
    # qend (end of the alignment in the query sequence)
    # sstart (start of the alignment in the subject sequence)
    # send (end of the alignment in the subject sequence)
    # evalue (E-value)
    # bitscore (bit score)
    # qlen (query sequence length)
    # slen (subject sequence length)
    
    
    
    df <- read_tsv(f, col_names = F) 
    
    colnames(df) <- outfmt6.names
    
    # df <- df %>% select(query, target, identity, e)
    
    df <- df %>% mutate(db = basename(f))
    
    
    
    return(df)
    
    
  }
  
  blast <- read_outfmt6(f)
  
  blast <- blast %>%
    dplyr::rename("query_length" = "qlen") %>%
    dplyr::rename("subject_length" = "slen")
  
  
  blast <- blast %>%
    filter(length >= 0.5 * query_length | length >= 0.5 * subject_length) %>%
    filter(pident >= 80)
  
  
  # Function to calculate overlap ratio between two alignments on reference or query
  overlap_ratio <- function(start1, end1, start2, end2) {
    if(any(is.na(c(start1, end1, start2, end2)))) return(NA)
    overlap <- max(0, min(end1, end2) - max(start1, start2) + 1)
    shortest <- min(abs(end1 - start1) + 1, abs(end2 - start2) + 1)
    return(overlap / shortest)
  }
  

  # contig_hits <- blast %>%
  #   group_by(qseqid) %>%
  #   summarise(n_hits = n(), subjects = list(unique(sseqid)))
  
  # For each subject, get contigs
  subject_contigs <- blast %>%
    group_by(sseqid) %>%
    summarise(n_contigs = n(), queries = list(unique(qseqid)))
  
  # Assign preliminary categories (from number and uniqueness of hits):
  
  # Helper function for assignment
  assign_preliminary_category <- function(qseqid) {
    
    hits <- blast %>% filter(qseqid == !!qseqid) # !! changes TRUE to FALSE, 
    
    n_hits <- nrow(hits)
    
    unique_subj <- length(unique(hits$sseqid))
    
    # For each subject in hits, count number of contigs hitting it
    subj_contigs <- subject_contigs %>% filter(sseqid %in% hits$sseqid)
    
    n_contigs_per_subj <- subj_contigs$n_contigs
    
    if (n_hits == 0) {
      return("no hit")
    } else if (n_hits == 1 & all(n_contigs_per_subj ==1)) {
      return("1->1")
    } else if (n_hits > 1 & all(n_contigs_per_subj ==1)) {
      return("m->1")
    } else if (n_hits == 1 & any(n_contigs_per_subj > 1)) {
      return("1->n")
    } else {
      return("m->n")
    }
  }
  
  contigs <- unique(blast$qseqid)
  
  prelim_cats <- sapply(contigs, assign_preliminary_category)
  
  annotation_df <- data.frame(qseqid = contigs, prelim_cat = prelim_cats,
    stringsAsFactors = FALSE)
  
  # Function to compute overlaps between intervals of alignments for a contig or subject
  get_overlap_matrix <- function(intervals) {
    n <- nrow(intervals)
    mat <- matrix(0, n, n)
    for (i in 1:n) {
      for (j in 1:n) {
        if(i != j) {
          mat[i,j] <- overlap_ratio(intervals$start[i], intervals$end[i], intervals$start[j], intervals$end[j])
        }
      }
    }
    return(mat)
  }
  
  # Refine annotation based on overlap and alignment coverage
  
  annotation_df$final_annotation <- NA
  
  for(i in 1:nrow(annotation_df)) {
    c_id <- annotation_df$qseqid[i]
    cat <- annotation_df$prelim_cat[i]
    hits <- blast %>% filter(qseqid == c_id)
    
    if(cat == "no hit") {
      annotation_df$final_annotation[i] <- "no hit"
    } else if(cat == "1->1") {
      # full if alignment covers ≥90% of subject length, else fragment
      hit <- hits[1,]
      cover_ratio <- hit$length / hit$subject_length
      if(cover_ratio >= 0.9) {
        annotation_df$final_annotation[i] <- "full"
      } else {
        annotation_df$final_annotation[i] <- "fragment"
      }
    } else if(cat == "m->1") {
      # multiple contigs hit 1 subject
      subject_id <- hits$sseqid[1]
      # Extract all contigs hitting this subject
      contigs_hitting_subject <- blast %>% filter(sseqid == subject_id)
      # For each contig, get query alignment intervals on subject
      intervals <- contigs_hitting_subject %>%
        select(qseqid, sstart, send) %>%
        distinct() %>%
        mutate(start = pmin(sstart, send), end = pmax(sstart, send))
      
      contig_ids <- unique(intervals$qseqid)
      # Compute overlap matrix between contigs on subject interval
      ov_mat <- matrix(0, length(contig_ids), length(contig_ids))
      rownames(ov_mat) <- colnames(ov_mat) <- contig_ids
      for(a in 1:length(contig_ids)) {
        for(b in 1:length(contig_ids)) {
          if(a != b) {
            i1 <- intervals %>% filter(qseqid == contig_ids[a])
            i2 <- intervals %>% filter(qseqid == contig_ids[b])
            # Considering the union of multiple hits: here simplifying - taking first interval only
            ov_mat[a,b] <- overlap_ratio(i1$start[1], i1$end[1], i2$start[1], i2$end[1])
          }
        }
      }
      # For this contig, check if it significantly overlaps any other contig (threshold 0.5)
      overlaps <- ov_mat[c_id, -which(colnames(ov_mat) == c_id)]
      if(any(overlaps > 0.5)) {
        annotation_df$final_annotation[i] <- "allele"
      } else {
        annotation_df$final_annotation[i] <- "fragment"
      }
    } else if(cat == "1->n") {
      # one contig hits multiple subjects
      # classify hits according to overlap on subject sequences
      # Extract all hits of this contig
      hits_intervals <- hits %>%
        mutate(start = pmin(sstart, send), end = pmax(sstart, send))
      # compare overlaps between subjects on subject seq
      overlap_flag <- FALSE
      for(x in 1:(nrow(hits_intervals)-1)) {
        for(y in (x+1):nrow(hits_intervals)) {
          ov <- overlap_ratio(hits_intervals$start[x], hits_intervals$end[x], hits_intervals$start[y], hits_intervals$end[y])
          if(!is.na(ov) > 0.5) {overlap_flag <- TRUE; break}
        }
        if(overlap_flag) break
      }
      if(overlap_flag) {
        annotation_df$final_annotation[i] <- "multi"
      } else {
        annotation_df$final_annotation[i] <- "chimera"
      }
    } else if(cat == "m->n") {
      # complex cases
      # Here simplified: if number of contigs equals number of hits and hits are unique, then 'full' or 'fragment', else 'multi'
      contig_ids <- unique(hits$qseqid)
      subj_ids <- unique(hits$sseqid)
      if(length(contig_ids) == length(subj_ids)) {
        # For each contig, check coverage like 1->1
        coverage_vec <- numeric(length(contig_ids))
        for(idx in 1:length(contig_ids)) {
          ht <- hits %>% filter(qseqid == contig_ids[idx])
          # If multiple hits keep longest
          max_hit <- ht[which.max(ht$length),]
          coverage_vec[idx] <- max_hit$length / max_hit$subject_length
        }
        # classify each contig
        annotation_vec <- ifelse(coverage_vec >= 0.9, "full", "fragment")
        # Here assign majority vote to this category for all contigs
        cat_assign <- ifelse(mean(coverage_vec>=0.9) >= 0.5, "full", "fragment")
        annotation_df$final_annotation[i] <- cat_assign
      } else {
        annotation_df$final_annotation[i] <- "multi"
      }
    } else {
      annotation_df$final_annotation[i] <- "other"
    }
  }
  
  # Output annotations
  annotation_df %>% 
    dplyr::count(prelim_cat, final_annotation) %>%
    mutate(file_name = basename(f))
  
  
}

# blast_annotation(file_list_1[1])
# blast_annotation(file_list_2[1])

# blast_annotation(file_list[6])

annotation_df_2 <- lapply(file_list_2, blast_annotation)

annotation_df_2 <- do.call(rbind,annotation_df_2)

annotation_df_2 <- annotation_df_2 %>%
  as_tibble() %>%
  mutate(vfold_set = sapply(strsplit(file_name, "_"), `[`, 1)) %>%
  # mutate(Assembler = sapply(strsplit(file_name, "_"), `[`, 7)) # 5
  mutate(kmer = sapply(strsplit(file_name, "_"), `[`, 7)) %>%
  mutate(kmer = gsub(".2.blast", "", kmer))

annotation_df_1 <- lapply(file_list_1, blast_annotation)

annotation_df_1 <- do.call(rbind,annotation_df_1)

annotation_df_1 <- annotation_df_1 %>%
  as_tibble() %>%
  mutate(vfold_set = sapply(strsplit(file_name, "_"), `[`, 1)) %>%
  mutate(kmer = sapply(strsplit(file_name, "_"), `[`, 5))


annotation_df <- rbind(annotation_df_2, annotation_df_1)

annotation_df %>% count(vfold_set, kmer)

annotation_df %>% count(kmer, final_annotation)


scale_fill <- c(ggsci::pal_startrek(alpha = 0.5)(7), ggsci::pal_cosmic(alpha = 0.5)(n-7))


annotation_df %>%
  group_by(kmer) %>%
  mutate(frac = n/sum(n)) %>%
  ggplot(aes(x = frac, y = kmer, fill = final_annotation)) +
  facet_wrap(final_annotation ~ ., scales = "free", nrow = 3) +
  geom_col() +
  ggsci::scale_fill_startrek(name = "") +
  my_custom_theme()
