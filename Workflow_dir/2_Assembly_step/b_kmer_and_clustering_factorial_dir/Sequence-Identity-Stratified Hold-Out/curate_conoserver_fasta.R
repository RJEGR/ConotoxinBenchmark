# =============================================================================
# curate_conoserver_fasta.R
#
# Produces PAIRED, ID-synchronised FASTA files from the curated ConoServer
# data set:
#     curated_nuc_conoServerDB.fasta    (nucleotide precursors)
#     curated_prot_conoServerDB.fasta   (protein precursors)
#
# Both files:
#   - carry the SAME header on every record: a single "hits" id (no "|"
#     concatenation), so record N in the nuc file corresponds to record N
#     in the prot file;
#   - contain exactly the same set of records, in the same order;
#   - are directly consumable by cluster_toxins_mmseqs.sh / run_clustering_
#     batch.slurm / analyze_mmseq_clustering_results.R (Strategy 3).
#
# WHY A PAIRED DEDUP IS NEEDED
# The original dedup_DNAStringSet() collapses identical sequences and pastes
# all their entry_ids into one header with "|". Running that independently on
# nucleotide and protein sequences yields different headers and different
# record counts -> the two files no longer share an id and the workflow
# breaks. This script deduplicates ONCE, at the (nucleotide, protein) pair
# level, choosing a single representative entry_id per unique precursor.
#
# USAGE
#   Option A - run inside conoServerDB.R after Nodedf is built:
#       source("curate_conoserver_fasta.R")            # uses Nodedf, outdir
#
#   Option B - standalone, from the curated RDS:
#       Rscript curate_conoserver_fasta.R \
#           --rds   INPUTS/curated_nuc_conoServerDB.rds \
#           --outdir INPUTS
#
# REQUIRED COLUMNS in Nodedf / the RDS:
#   entry_id, sequence (nucleotide), proteinsequence, genesuperfamily
#   (organismlatin, organismdiet, organismregion optional -> metadata sidecar)
# =============================================================================

.outdir <- "~/Documents/GitHub/ConotoxinBenchmark/INPUTS/"

f <- file.path(.outdir, "curated_nuc_conoServerDB.rds")
if (!file.exists(f)) stop("Cannot find input file: ", f)

.input_df <- read_rds(f)

suppressPackageStartupMessages({
  library(tidyverse)
  library(Biostrings)
})

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

# -----------------------------------------------------------------------------
# Resolve input: either an in-memory Nodedf (Option A) or CLI args (Option B)
# -----------------------------------------------------------------------------
.args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  i <- match(flag, .args)
  if (!is.na(i) && i < length(.args)) .args[i + 1] else default
}

# if (exists("Nodedf", inherits = TRUE) && is.data.frame(get0("Nodedf"))) {
#   message("curate_conoserver_fasta.R: using in-memory 'Nodedf'")
#   .input_df <- get("Nodedf", inherits = TRUE)
#   .outdir   <- if (exists("outdir", inherits = TRUE)) get("outdir", inherits = TRUE) else "."
# } else {
#   .rds <- get_arg("--rds")
#   .outdir <- get_arg("--outdir", ".")
#   if (is.null(.rds))
#     stop("No in-memory 'Nodedf' found and no --rds supplied. ",
#          "Provide --rds <curated_nuc_conoServerDB.rds> --outdir <dir>.")
#   message("curate_conoserver_fasta.R: reading ", .rds)
#   .input_df <- readRDS(.rds)
# }

dir.create(.outdir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Validate required columns
# -----------------------------------------------------------------------------
required <- c("entry_id", "sequence", "proteinsequence")
missing  <- setdiff(required, colnames(.input_df))
if (length(missing) > 0)
  stop("Input is missing required column(s): ", paste(missing, collapse = ", "))

# genesuperfamily is needed downstream (Strategy 1 LOSO); keep if present
has_sf <- "genesuperfamily" %in% colnames(.input_df)

# -----------------------------------------------------------------------------
# 1. Clean records: drop rows with empty nuc or protein, normalise case,
#    sanitise the id so it is FASTA- and MMseqs2-safe (no whitespace).
# -----------------------------------------------------------------------------
df <- .input_df %>%
  as_tibble() %>%
  mutate(
    entry_id        = str_trim(as.character(entry_id)),
    sequence        = str_to_upper(str_trim(as.character(sequence))),
    proteinsequence = str_to_upper(str_trim(as.character(proteinsequence)))
  ) %>%
  filter(!is.na(entry_id), entry_id != "",
         !is.na(sequence), sequence != "",
         !is.na(proteinsequence), proteinsequence != "")

n_start <- nrow(df)
message(sprintf("Records after dropping empty nuc/protein: %d", n_start))

# Make the id FASTA-safe: MMseqs2 splits headers on whitespace, so the id
# must contain none. Replace any non [A-Za-z0-9._-] with "_".
df <- df %>%
  mutate(hits = str_replace_all(entry_id, "[^A-Za-z0-9._-]", "_"))

# Guard against id collisions introduced by sanitisation
dup_ids <- df %>% count(hits) %>% filter(n > 1)
if (nrow(dup_ids) > 0) {
  message(sprintf("Note: %d id(s) collided after sanitising; disambiguating with suffixes.",
                  nrow(dup_ids)))
  df <- df %>%
    group_by(hits) %>%
    mutate(hits = if (n() > 1) paste0(hits, "_", row_number()) else hits) %>%
    ungroup()
}

# -----------------------------------------------------------------------------
# 2. PAIRED dedup: collapse records identical in BOTH nucleotide and protein
#    sequence. One representative 'hits' id is kept per unique precursor.
#    (Two records with the same nuc but different protein are NOT merged.)
# -----------------------------------------------------------------------------
df_dedup <- df %>%
  arrange(hits) %>%                                  # deterministic representative
  group_by(sequence, proteinsequence) %>%
  summarise(
    hits        = dplyr::first(hits),                # representative id
    n_collapsed = n(),
    merged_ids  = paste(sort(unique(hits)), collapse = ";"),
    genesuperfamily = if (has_sf) dplyr::first(genesuperfamily) else NA_character_,
    organismlatin   = if ("organismlatin"  %in% colnames(df)) dplyr::first(organismlatin)  else NA_character_,
    organismdiet    = if ("organismdiet"   %in% colnames(df)) dplyr::first(organismdiet)   else NA_character_,
    organismregion  = if ("organismregion" %in% colnames(df)) dplyr::first(organismregion) else NA_character_,
    .groups = "drop"
  ) %>%
  arrange(hits)

n_dedup <- nrow(df_dedup)
message(sprintf("Unique (nucleotide, protein) precursors: %d  (collapsed %d duplicate record(s))",
                n_dedup, n_start - n_dedup))

# -----------------------------------------------------------------------------
# 3. Build the two StringSets with IDENTICAL headers and identical order
# -----------------------------------------------------------------------------
nuc_set  <- DNAStringSet(df_dedup$sequence)
prot_set <- AAStringSet(df_dedup$proteinsequence)
names(nuc_set)  <- df_dedup$hits
names(prot_set) <- df_dedup$hits

# Hard invariant check: same names, same order, same count
stopifnot(length(nuc_set) == length(prot_set))
stopifnot(identical(names(nuc_set), names(prot_set)))

# -----------------------------------------------------------------------------
# 4. Write outputs
# -----------------------------------------------------------------------------
nuc_fa  <- file.path(.outdir, "curated_nuc_conoServerDB.fasta")
prot_fa <- file.path(.outdir, "curated_prot_conoServerDB.fasta")
meta_tsv <- file.path(.outdir, "curated_conoServerDB_metadata.tsv")

writeXStringSet(nuc_set,  nuc_fa) # ,  width = 80
writeXStringSet(prot_set, prot_fa) # , width = 80

# Metadata sidecar keyed by 'hits' (carries superfamily for LOSO, plus the
# list of original entry_ids that were collapsed into each representative).
df_dedup %>%
  select(hits, genesuperfamily, organismlatin, organismdiet, organismregion,
         n_collapsed, merged_ids) %>%
  write_tsv(meta_tsv)

# Also refresh the curated RDS so downstream rsample/LOSO scripts stay in sync
rds_out <- file.path(.outdir, "curated_nuc_conoServerDB.rds")
df_dedup %>%
  transmute(entry_id = hits, hits, sequence, proteinsequence,
            genesuperfamily, organismlatin, organismdiet, organismregion) %>%
  saveRDS(rds_out)

# -----------------------------------------------------------------------------
# 5. Report
# -----------------------------------------------------------------------------
message("\n=== curate_conoserver_fasta.R: done ===")
message(sprintf("  nucleotide FASTA : %s  (%d records)", nuc_fa,  length(nuc_set)))
message(sprintf("  protein FASTA    : %s  (%d records)", prot_fa, length(prot_set)))
message(sprintf("  metadata sidecar : %s", meta_tsv))
message(sprintf("  refreshed RDS    : %s", rds_out))
message("\nHeader invariant: every record in both FASTAs carries one 'hits' id,")
message("no '|' concatenation, identical id set and order across the two files.")
message("\nNext (Strategy 3):")
message(sprintf("  sbatch run_clustering_batch.slurm %s %s", nuc_fa, prot_fa))
