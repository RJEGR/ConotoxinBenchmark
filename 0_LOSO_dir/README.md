# LOSO (Leave-One-Superfamily-Out) Pipeline

Strategy 1 of the circularity-controlled conotoxin assembly benchmark.

## What changed vs. the random 12-fold CV pipeline

| | Before (random 12-fold) | After (LOSO) |
|---|---|---|
| **rsample call** | `vfold_cv(df, v = 12)` | `group_vfold_cv(df, group = genesuperfamily, v = n_distinct(genesuperfamily))` |
| **Fold composition** | Each fold is a random ~8.3% slice — every fold contains a mix of all superfamilies | Each fold holds out one entire superfamily; the remaining superfamilies form the training set |
| **What it tests** | Within-distribution recall under sampling variance | Generalization across gene-family boundaries (Grimm Type-2 circularity) |
| **Number of folds** | 12 | n_distinct(superfamily) — typically ~15–25 after filtering for SFs with n ≥ 5 |
| **FASTA outputs** | One subsampled FASTA per fold | Two FASTAs per fold: `train_<sf>.fasta` + `test_<sf>.fasta` |

## File layout

```
loso_pipeline/
├── LOSO_conoServerDB.R          # R: splits the curated ConoServer DB into LOSO folds
├── Simulate_rnaseq_loso.sh      # bash: simulates reads from test_*.fasta only
├── Run_LOSO_pipeline.sh         # bash: top-level orchestrator (sim + assemble + leak-check)
├── LOSO_aggregate.R             # R: aggregates per-fold TransRate metrics into a per-SF report
└── README.md                    # this file
```

## Run order

```bash
# Step 1: produce the LOSO splits (run once, on the head node)
Rscript LOSO_conoServerDB.R
# Output: INPUTS/vfolds_loso_resampling_dir/{train_*.fasta, test_*.fasta,
#                                            LOSO_manifest.tsv, LOSO_power_diagnostics.tsv,
#                                            superfamily_distribution.tsv}

# Step 2: full pipeline (HPC job)
sbatch ./Run_LOSO_pipeline.sh \
    -d INPUTS/vfolds_loso_resampling_dir \
    -c run.config
# This calls:
#   2a. Simulate_rnaseq_loso.sh       (reads from test_*.fasta only)
#   2b. Assemblers.sh (your existing) (one run per fold's manifest)
#   2c. BLAST leakage check           (assemblies vs. train_*.fasta;
#                                      flags contigs that match training)

# Step 3: aggregate (head node)
Rscript LOSO_aggregate.R \
    --loso_dir INPUTS/vfolds_loso_resampling_dir \
    --sim_dir  <md5>_loso_dir \
    --out_dir  <md5>_loso_dir/LOSO_summary
```

## Key design choices

**Minimum superfamily size (MIN_SF_SIZE = 5).** Superfamilies with fewer
sequences are pooled into the training set of every fold but are never held
out, because a held-out set of n < 5 cannot yield stable per-fold metrics.
This is configurable at the top of `LOSO_conoServerDB.R`.

**The training set is not used for assembly.** It is reserved for:
(a) building BLAST databases for the leakage check (Step 2c), (b) calibrating
any classifier whose training data you control (ConoSorter itself is fixed
and cannot be retrained per fold — that limitation should be documented in
the manuscript), and (c) optionally serving as the reference for reference-
guided assemblers like StringTie (Strategy 5).

**TransRate ground truth = the held-out test FASTA.** Reads were simulated
from `test_<sf>.fasta`, so that file is the only legitimate ground truth for
scoring sensitivity/precision on this fold. The manifest written by
`Simulate_rnaseq_loso.sh` puts `test_<sf>.fasta` in column 4 to stay
compatible with the existing `Assemblers.sh`.

**Leakage check.** Under LOSO, no assembled contig should match the training
reference at high identity, because the simulated reads came from the
held-out superfamily. Any high-identity training hits indicate either
(i) genuine paralog overlap across superfamilies (a known biological feature
of conotoxins) or (ii) annotation noise in ConoServer (mis-classified
sequences). Either way, the rate is a quantitative integrity check on the
LOSO design.

## What to report in the manuscript

1. **Per-superfamily sensitivity and precision** (from
   `LOSO_per_superfamily_summary.tsv`). Show the headline plot
   `LOSO_sensitivity_by_superfamily.png` as a panel alongside Figure 1.
2. **In-distribution (your previous 12-fold random CV) vs. out-of-distribution
   (LOSO) performance per assembler**, with a paired Wilcoxon test on the
   per-fold metrics. The expected result is that LOSO precision is
   meaningfully lower than random-CV precision; this magnitude is the
   empirical answer to reviewers asking "isn't your benchmark circular?"
3. **Leakage check summary**: report the median and 95th percentile of
   `hits_pident90_qcov80` across folds. Low values (≤2–3 per assembly) are
   consistent with paralog overlap; substantially higher values would
   indicate a problem with the splits.

## Notes

- The R script normalizes `genesuperfamily` labels by stripping whitespace
  and recoding empty/NA as `"UNCLASSIFIED"`. Check
  `superfamily_distribution.tsv` after Step 1 to confirm the categories
  match your expectations.
- All scripts use the variable `SEED = 20260511` for reproducibility. Change
  in the top of `LOSO_conoServerDB.R` if you want a different random seed.
- The leakage check requires `blastn` / `makeblastdb` in PATH. Use `-S` in
  `Run_LOSO_pipeline.sh` to skip if unavailable.
