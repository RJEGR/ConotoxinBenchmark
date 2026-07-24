# Smoke test for LOSO_metrics_core.R
# Builds a tiny synthetic TransRate csv + leakage table, runs the chain
# read_contigs -> calculate_metrics -> three_way_classify -> compute_cif,
# and asserts that the output columns and counts make sense.


setwd("/Users/rjegr/Documents/GitHub/ConotoxinBenchmark/Workflow_dir/0_Data_selection/0_LOSO_dir/LOSO-framework/")

source("LOSO_metrics_core.R")

tmp <- tempfile(fileext = "_dir")
dir.create(file.path(tmp, "test_A_superfamily_200x_PE_IDBA_dir"),
           recursive = TRUE)
csv_path <- file.path(tmp, "test_A_superfamily_200x_PE_IDBA_dir",
                      "contigs.csv")

# Synthetic contigs.csv: 10 contigs, mixed quality
fake <- tibble::tibble(
  contig_name        = paste0("k51_", 1:10),
  length             = c(800, 750, 900, 850, 700, 650, 720, 680, 600, 950),
  reference_coverage = c(0.95, 0.92, 0.40, 0.91, 0.10, 0.0,
                         0.0,  0.93, 0.05, 0.97),
  hits               = c("refA1","refA2","refA1","refA3","refA1",
                          NA,     NA,     "refA4","refA1","refA5")
)
readr::write_csv(fake, csv_path)

cat("\n--- read_contigs ---\n")
d <- read_contigs(csv_path)
stopifnot(!is.null(d), nrow(d) == 10)
cat("sample_id =", unique(d$sample_id), "  tool =", unique(d$tool), "\n")

cat("\n--- calculate_metrics @ 0.90 ---\n")
input_n_seq <- tibble::tibble(
  sample_id       = "test_A_superfamily",
  InputNsequences = 7   # pretend the held-out fold had 7 reference seqs
)
m <- calculate_metrics(d, reference_coverage_val = 0.90,
                      input_n_seq_tbl = input_n_seq)
print(m)
stopifnot("F1" %in% colnames(m), "precision" %in% colnames(m),
          "sensitivity" %in% colnames(m))
# Under Behera et al. 2022:
#   TP_contigs = contigs >= 0.90  -> 5
#   FP = partial (3) + no-hit (2) = 5
#   TP_refs    = distinct hits at >=0.90 -> 5  (refA1..refA5)
#   FN = 7 - 5 = 2
#   Precision = 5 / (5+5) = 0.5
#   Sensitivity = 5 / 7   = 0.7142857
stopifnot(m$TP == 5, m$FP == 5, m$TP_refs == 5, m$FN == 2)
stopifnot(abs(m$precision - 0.5) < 1e-9)
stopifnot(abs(m$sensitivity - 5/7) < 1e-9)
cat("calculate_metrics OK: TP=5  FP=5  TP_refs=5  FN=2  P=",
    round(m$precision,3), " S=", round(m$sensitivity,3),
    " F1=", round(m$F1,3), "\n")

cat("\n--- three_way_classify ---\n")
leak <- tibble::tibble(
  assembly             = "test_A_superfamily_200x_PE_IDBA",
  held_out_sf          = "A_superfamily",
  total_hits           = 6,
  hits_pident90_qcov80 = 2   # 2 contigs hit train at >=90/80
)
tw <- three_way_classify(m, leak)
print(tw)
stopifnot("FP_paralog" %in% colnames(tw))
# 5 contigs are non-TP; leakage says 2 hit train at high identity
stopifnot(tw$FP_paralog == 2)
stopifnot(tw$FP_spurious == 10 - 5 - 2)  # 3
# Now both precisions use rawcontigs as denominator -> consistent with
# calculate_metrics' precision (TP/(TP+FP) = TP/rawcontigs when
# FP = rawcontigs - TP).
stopifnot(abs(tw$precision_ood - 5/10) < 1e-9)
stopifnot(abs(tw$precision_combined - 7/10) < 1e-9)
cat("three_way OK: TP=5  FP_paralog=2  FP_spurious=3  P_ood=0.50  P_comb=0.70\n")

cat("\n--- compute_cif ---\n")
cif <- compute_cif(tw)
print(cif$per_tool)
# CIF should be 0.7 / 0.5 = 1.4
stopifnot(abs(cif$per_tool$mean_CIF - 1.4) < 1e-9)
cat("CIF OK: mean_CIF =", cif$per_tool$mean_CIF, "\n")

cat("\nALL SMOKE TESTS PASSED\n")
