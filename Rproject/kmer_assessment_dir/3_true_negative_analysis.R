
#  known and correctly annotated conotoxins (representing our true positive dataset) and some random sequences of venomous animals (spiders, snakes, sea urchins, jellyfish) available on UniProtKB (to represent the true negative dataset).


outdir <- "~/Documents/GitHub/ConotoxinBenchmark/INPUTS"

file_out <- file.path(outdir, "uniprotkb_taxonomy_id_33208_AND_cc_tiss_2025_02_19.fasta.gz")

# read fasta, filter conotoxins, and select random seed and sample, then perform entropy analysis
# data <- read_rds(file_out)

