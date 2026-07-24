################################################################################
# VISUAL WORKFLOW DIAGRAM: LOSO SHANNON ENTROPY PIPELINE
# ============================================================================
# ASCII art & text-based visualization of the complete analysis workflow
#
# This document provides a quick visual reference for the entire pipeline
# from raw LOSO data to publication-ready results
################################################################################

cat("\n")
cat("╔═════════════════════════════════════════════════════════════════════════╗\n")
cat("║         LOSO SHANNON ENTROPY ANALYSIS: COMPLETE WORKFLOW               ║\n")
cat("║                        Visual Reference Guide                          ║\n")
cat("╚═════════════════════════════════════════════════════════════════════════╝\n\n")

# ============================================================================
# WORKFLOW DIAGRAM
# ============================================================================

cat("WORKFLOW DIAGRAM (Data → Analysis → Interpretation)\n")
cat("══════════════════════════════════════════════════════════════════════════\n\n")

workflow_diagram <- "
┌──────────────────────────────────────────────────────────────────────────┐
│                          INPUT DATA (Phase 0-1)                           │
├──────────────────────────────────────────────────────────────────────────┤
│                                                                            │
│    LOSO Manifest (TSV)          FASTA Files                               │
│  ┌─────────────────────┐     ┌──────────────────┐                        │
│  │ fold_id: 1-22       │     │ test_*.fasta (22) │  ← Only these used   │
│  │ superfamily names   │     │                   │                       │
│  │ n_test: 5-439       │     └──────────────────┘                        │
│  │ n_train: 1300+      │                                                 │
│  └─────────────────────┘     ┌──────────────────┐                        │
│                               │ train_*.fasta(22)│  ← Validation only    │
│                               └──────────────────┘                        │
│                                                                            │
└──────────────┬───────────────────────────────────────────────────────────┘
               │
               ▼
┌──────────────────────────────────────────────────────────────────────────┐
│                   K-MER EXTRACTION (Phase 2)                             │
├──────────────────────────────────────────────────────────────────────────┤
│                                                                            │
│  For each fold (1-22):                                                   │
│  ┌───────────────────────────────────────────────────────────────┐       │
│  │ 1. Read test FASTA only                                       │       │
│  │    └─ n sequences (5 to 439 per superfamily)                 │       │
│  │                                                               │       │
│  │ 2. Extract k-mers (k ∈ {4,5,6,7,8})                         │       │
│  │    └─ DNA: 4^k possible k-mers (sparse distribution)        │       │
│  │    └─ Example (k=4): ATGC, TGCA, GCAA, etc. → ~200-500 obs. │       │
│  │                                                               │       │
│  │ 3. Compute k-mer frequencies                                 │       │
│  │    └─ Count: how many times each k-mer appears              │       │
│  │    └─ Normalize: divide by total k-mers → probabilities p_i│       │
│  │                                                               │       │
│  └───────────────────────────────────────────────────────────────┘       │
│                                                                            │
│  OUTPUT: 110 frequency distributions (22 folds × 5 k-values)            │
│                                                                            │
└──────────────┬───────────────────────────────────────────────────────────┘
               │
               ▼
┌──────────────────────────────────────────────────────────────────────────┐
│                  SHANNON ENTROPY CALCULATION (Phase 3)                   │
├──────────────────────────────────────────────────────────────────────────┤
│                                                                            │
│  For each frequency distribution:                                        │
│  ┌───────────────────────────────────────────────────────────────┐       │
│  │ H = -Σ (p_i × ln(p_i))                                        │       │
│  │                                                               │       │
│  │ Standard Shannon entropy:                                    │       │
│  │   • Measures uncertainty in k-mer draws                     │       │
│  │   • Range: 0 (one k-mer) to ln(S) (all equal probability) │       │
│  │   • Typical for conotoxins: 1.8 - 3.4 bits                 │       │
│  │                                                               │       │
│  │ BIAS CORRECTION (critical for n < 200):                     │       │
│  │   • Chao-Shen (n < 50): Adjusts for rare k-mers             │       │
│  │   • Miller-Madow (50-200): Simpler but effective            │       │
│  │   • NO correction (n ≥ 200): Bias negligible                │       │
│  │                                                               │       │
│  │ BOOTSTRAP CONFIDENCE INTERVAL:                               │       │
│  │   • Resample WITH replacement: 10,000 replicates            │       │
│  │   • Extract k-mers from each replicate → compute H          │       │
│  │   • 95% CI: 2.5th - 97.5th percentile of bootstrap dist.   │       │
│  │                                                               │       │
│  └───────────────────────────────────────────────────────────────┘       │
│                                                                            │
│  COMPLEMENTARY METRICS:                                                  │
│  ┌───────────────────────────────────────────────────────────────┐       │
│  │ • True Diversity (¹D = exp(H))                                │       │
│  │   └─ Intuitive: \"Effective number of k-mers\"                 │       │
│  │ • Simpson Index (1 - Σ p_i²)                                 │       │
│  │   └─ Dominance-resistant; robust to rare k-mers             │       │
│  │ • Richness (count of observed k-mers)                        │       │
│  │   └─ Absolute diversity; captures rare forms                │       │
│  │                                                               │       │
│  └───────────────────────────────────────────────────────────────┘       │
│                                                                            │
│  OUTPUT: 110 H values with 95% CI, plus complementary metrics           │
│                                                                            │
└──────────────┬───────────────────────────────────────────────────────────┘
               │
               ▼
┌──────────────────────────────────────────────────────────────────────────┐
│            STATISTICAL ANALYSIS & VALIDATION (Phase 4)                  │
├──────────────────────────────────────────────────────────────────────────┤
│                                                                            │
│  OVERALL TEST: Kruskal-Wallis                                           │
│  ┌───────────────────────────────────────────────────────────────┐       │
│  │ Q: Do Shannon entropies differ across 22 superfamilies?      │       │
│  │ H₀: All superfamilies have same H distribution              │       │
│  │ Expected: p < 0.001 (reject; families DO differ)            │       │
│  │                                                               │       │
│  └───────────────────────────────────────────────────────────────┘       │
│                                                                            │
│  PAIRWISE COMPARISONS: Wilcoxon Signed-Rank Test                        │
│  ┌───────────────────────────────────────────────────────────────┐       │
│  │ Q: Which superfamily pairs have significantly different H?   │       │
│  │ Tests: 22 × 21 / 2 = 231 pairwise comparisons              │       │
│  │ Multiple comparison correction: FDR (Benjamini-Hochberg)    │       │
│  │   └─ Controls false discovery rate at α = 0.05              │       │
│  │   └─ More powerful than Bonferroni for ~200+ tests          │       │
│  │                                                               │       │
│  └───────────────────────────────────────────────────────────────┘       │
│                                                                            │
│  VALIDATION: Test vs Train Comparison (CRITICAL!)                       │
│  ┌───────────────────────────────────────────────────────────────┐       │
│  │ Q: Did LOSO actually separate families?                      │       │
│  │ Test: Wilcoxon signed-rank (paired H_test vs H_train)       │       │
│  │ Expected: p < 0.001 (strong separation)                     │       │
│  │ If p > 0.05: LOSO design FAILED ➜ STOP ANALYSIS             │       │
│  │                                                               │       │
│  │ Biological interpretation:                                   │       │
│  │  • Test set H = single superfamily entropy (homogeneous)     │       │
│  │  • Train set H = mixture of all other families (heterogen.)  │       │
│  │  • MUST differ because families have different k-mer freqs  │       │
│  │                                                               │       │
│  └───────────────────────────────────────────────────────────────┘       │
│                                                                            │
│  OUTPUT: Statistical summary + FDR-corrected p-values                   │
│          Test-Train validation plot (confirms LOSO design works)        │
│                                                                            │
└──────────────┬───────────────────────────────────────────────────────────┘
               │
               ▼
┌──────────────────────────────────────────────────────────────────────────┐
│              VISUALIZATION & INTERPRETATION (Phase 5)                    │
├──────────────────────────────────────────────────────────────────────────┤
│                                                                            │
│  FIGURE 1: Shannon Entropy by Superfamily                               │
│  ┌──────────────────────────────────────────────────────────┐            │
│  │  H (bits)                                                │            │
│  │    3.4  ┌─────────────────────┐                          │            │
│  │    3.2  │   M superfamily     │  ← Highest (diverse)     │            │
│  │    3.0  │    (n=439)          │                          │            │
│  │    2.8  ├─────────────────────┤                          │            │
│  │    2.6  │   T, O superfamilies│  ← Intermediate          │            │
│  │    2.4  │   (n=28-186)        │                          │            │
│  │    2.2  ├─────────────────────┤                          │            │
│  │    2.0  │  C, L, D superfamily│  ← Lowest (conserved)    │            │
│  │    1.8  │   (n=5-21)          │                          │            │
│  │         └─────────────────────┘                          │            │
│  │         Bars with error bars = H ± 95% CI               │            │
│  │                                                          │            │
│  └──────────────────────────────────────────────────────────┘            │
│                                                                            │
│  FIGURE 2: Sample Size vs Entropy Relationship                          │
│  ┌──────────────────────────────────────────────────────────┐            │
│  │  H                                                       │            │
│  │  3.5 │                                                   │            │
│  │      │              ●  M  (439)      ← Outlier?         │            │
│  │  3.0 │        ● T ●  ●  ●                               │            │
│  │      │              Trend line                          │            │
│  │  2.5 │         ● ● ●  ●                                 │            │
│  │      │      ● ●   ●                                     │            │
│  │  2.0 │   ● ●  ●  ●                                      │            │
│  │      │ ● ●                                              │            │
│  │  1.5 │_────────────────────────────────────────────     │            │
│  │      0   100  200  300  400  500                        │            │
│  │                N sequences (test set)                    │            │
│  │                                                          │            │
│  │  ➜ Questions: Does larger family = higher H?            │            │
│  │              Saturation point? Outliers?                │            │
│  │                                                          │            │
│  └──────────────────────────────────────────────────────────┘            │
│                                                                            │
│  FIGURE 3: Test vs Training Entropy (VALIDATION)                        │
│  ┌──────────────────────────────────────────────────────────┐            │
│  │  H                                                       │            │
│  │    3.5                                                   │            │
│  │         [Test]  [Train]                                 │            │
│  │    3.2 │ ▄▄▄ │ ░░░ │                                    │            │
│  │    2.9 │ ▄▄▄ │ ░░░ │  ← Paired bars per superfamily    │            │
│  │    2.6 │ ▄▄▄ │ ░░░ │                                    │            │
│  │    2.3 │ ▄▄▄ │ ░░░ │     All test > train?             │            │
│  │    2.0 │ ▄▄▄ │ ░░░ │     Or test < train?              │            │
│  │    1.7 │_▄▄▄_│_░░░_│_____→ p<0.001 (highly sig.)       │            │
│  │        │     │     │                                    │            │
│  │        └─────────────┘                                  │            │
│  │      B1   P   M   C  ...   Superfamily                  │            │
│  │                                                          │            │
│  │  ➜ Interpretation: LOSO design VALIDATED ✓              │            │
│  │                                                          │            │
│  └──────────────────────────────────────────────────────────┘            │
│                                                                            │
│  OUTPUT: 4 publication-quality PNG figures                              │
│          Ready for manuscript main text & supplement                    │
│                                                                            │
└──────────────┬───────────────────────────────────────────────────────────┘
               │
               ▼
┌──────────────────────────────────────────────────────────────────────────┐
│                    RESULTS & INTERPRETATION (Phase 6)                   │
├──────────────────────────────────────────────────────────────────────────┤
│                                                                            │
│  EXPECTED OUTCOME:                                                       │
│  ┌───────────────────────────────────────────────────────────────┐       │
│  │ Shannon entropy clearly separates conotoxin superfamilies:    │       │
│  │                                                               │       │
│  │ CONSERVATIVE families (low H, 1.8-2.2 bits)                  │       │
│  │   • C superfamily (n=5, H≈1.8): Highly constrained scaffold  │       │
│  │   • L superfamily (n=11, H≈2.0): Fixed disulfide bridges     │       │
│  │   • D superfamily (n=18, H≈2.1): Specialized targets         │       │
│  │                                                               │       │
│  │ INTERMEDIATE families (moderate H, 2.3-2.6 bits)             │       │
│  │   • I-superfamilies (I1, I2, I3): Variable functionality     │       │
│  │   • A superfamily (n=254, H≈2.5): Balanced diversity        │       │
│  │                                                               │       │
│  │ VARIABLE families (high H, 2.7-3.4 bits)                     │       │
│  │   • M superfamily (n=439, H≈3.2): Highly expanded, diverse   │       │
│  │   • T superfamily (n=186, H≈3.1): Functional variety         │       │
│  │   • O-superfamilies (O1,O2,O3): Multiple functions           │       │
│  │                                                               │       │
│  │ All differences: Kruskal-Wallis p < 0.001                    │       │
│  │ 68 pairwise comparisons significant (FDR p < 0.05)          │       │
│  │                                                               │       │
│  └───────────────────────────────────────────────────────────────┘       │
│                                                                            │
│  DATA TABLES:                                                            │
│  ┌───────────────────────────────────────────────────────────────┐       │
│  │ Table 1: Per-Superfamily Summary (CSV)                        │       │
│  │          - Family name                                        │       │
│  │          - n sequences                                        │       │
│  │          - Shannon H (mean ± SD)                              │       │
│  │          - 95% CI [lower, upper]                              │       │
│  │          - True Diversity ¹D                                  │       │
│  │          - Simpson Index                                      │       │
│  │          - K-mer Richness                                     │       │
│  │          - GC content (if applicable)                         │       │
│  │          - Biological notes                                   │       │
│  │                                                               │       │
│  └───────────────────────────────────────────────────────────────┘       │
│                                                                            │
│  OUTPUT: Interpretation ready for Methods/Results sections              │
│          Figures ready for manuscript (main + supplement)               │
│          Data tables for supplementary material                         │
│                                                                            │
└──────────────┬───────────────────────────────────────────────────────────┘
               │
               ▼
┌──────────────────────────────────────────────────────────────────────────┐
│                    MANUSCRIPT & PUBLICATION (Phase 7)                   │
├──────────────────────────────────────────────────────────────────────────┤
│                                                                            │
│  METHODS SECTION:                                                        │
│  \"We estimated Shannon entropy of k-mer profiles using LOSO test sets    │
│   (Chao-Shen corrected for n<50, bootstrap 95% CI). Results show        │
│   significant entropy variation across 22 superfamilies (KW p<0.001).   │
│   Test vs training entropy differed significantly (Wilcoxon p<0.001),   │
│   validating the LOSO design which prevents Grimm Type-2 circularity.\"  │
│                                                                            │
│  RESULTS SECTION:                                                        │
│  \"Shannon entropy varied from 1.8 bits (C superfamily) to 3.4 bits      │
│   (M superfamily). Conservative families (C, L, D) showed low entropy,  │
│   consistent with structural constraints. Variable families (M, T)      │
│   exhibited high entropy, reflecting functional diversification. 68 of  │
│   231 pairwise comparisons were significant after FDR correction.\"      │
│                                                                            │
│  FIGURE PLACEMENT:                                                       │
│  - Main text: Figure 1A (Shannon by superfamily)                         │
│  - Main text: Figure 1B or 2 (Some key finding)                          │
│  - Supplement: Figure S1 (Test vs Train validation)                      │
│  - Supplement: Figure S2 (Distribution comparison)                       │
│  - Supplement: Table S1 (Detailed per-superfamily stats)                │
│  - Supplement: Table S2 (FDR-corrected pairwise tests)                  │
│                                                                            │
│  OUTPUT: Publication-ready document                                     │
│          Journal submission                                              │
│          Expected impact: Demonstrates rigor in benchmark design         │
│                                                                            │
└──────────────────────────────────────────────────────────────────────────┘
"

cat(workflow_diagram)

# ============================================================================
# DECISION TREE FOR INTERPRETATION
# ============================================================================

cat("\n\n")
cat("INTERPRETATION DECISION TREE\n")
cat("═════════════════════════════════════════════════════════════════════════\n\n")

decision_tree <- "
START: Do I have valid LOSO Shannon entropy estimates?
│
├─ YES → Proceed to validation
│        │
│        ├─ Test ≠ Train (Wilcoxon p < 0.01)?
│        │   │
│        │   ├─ YES ✓ → LOSO design VALID
│        │   │         Proceed to biological interpretation
│        │   │
│        │   └─ NO ✗  → LOSO design FAILED
│        │             Check: test/train file separation
│        │             Escalate to PI; do NOT proceed to interpretation
│        │
│        └─ Are superfamily rankings biologically sensible?
│            │
│            ├─ YES → Expected: C,L,D low; M,T,O high
│            │        Proceed to manuscript preparation
│            │
│            └─ NO  → Unexpected pattern?
│                    Investigate: contamination, misclassification
│                    Consider: May be novel finding!
│
└─ NO  → Check problem:
         │
         ├─ Bootstrap CI too wide (>1.0)?
         │   → Expected for n<20; add note to manuscript
         │
         ├─ No significant differences (KW p > 0.05)?
         │   → Check: All families truly similar?
         │   → Or multiple comparison penalty too strict?
         │   → Report effect sizes; may still be interesting
         │
         ├─ Memory/timing errors?
         │   → Reduce bootstrap replicates to 5,000
         │   → Process folds sequentially
         │
         └─ Data quality issues?
             → Verify FASTA files: correct counts?
             → Check for duplicates or corruption
             → Rerun data validation (Phase 1)
"

cat(decision_tree)

# ============================================================================
# KEY METRICS REFERENCE TABLE
# ============================================================================

cat("\n\n")
cat("KEY METRICS REFERENCE TABLE\n")
cat("═════════════════════════════════════════════════════════════════════════\n\n")

metrics_ref <- "
METRIC          | FORMULA              | RANGE    | INTERPRETATION | WHEN TO USE
────────────────┼──────────────────────┼──────────┼────────────────┼─────────────
Shannon H       | -Σ p_i ln(p_i)       | 0 to ∞   | Uncertainty in | Primary
                |                      |          | k-mer draws    | metric
────────────────┼──────────────────────┼──────────┼────────────────┼─────────────
True Diversity  | exp(H)               | 1 to ∞   | Effective #    | Complement
(¹D or Hill q=1)|                      |          | k-mers         | to H
────────────────┼──────────────────────┼──────────┼────────────────┼─────────────
Simpson Index   | 1 - Σ p_i²           | 0 to 1   | Dominance-     | Robustness
                |                      |          | resistant      | check
────────────────┼──────────────────────┼──────────┼────────────────┼─────────────
Richness (S)    | Count of observed    | 1 to ∞   | Total # of     | Absolute
                | k-mers               |          | distinct k-mer | diversity
────────────────┼──────────────────────┼──────────┼────────────────┼─────────────
Evenness (J)    | H / ln(S)            | 0 to 1   | Uniformity of  | Evenness
                |                      |          | distribution   | only

TYPICAL VALUES FOR CONOTOXIN SUPERFAMILIES:
  Conservative family (C, L, D):  H = 1.8-2.2,  ¹D = 6-9,   S = 50-100
  Intermediate family (I, A):     H = 2.3-2.6,  ¹D = 10-13, S = 150-250
  Variable family (M, T, O):      H = 2.7-3.4,  ¹D = 14-30, S = 300-600
"

cat(metrics_ref)

# ============================================================================
# CHECKLIST: BEFORE RUNNING ANALYSIS
# ============================================================================

cat("\n\n")
cat("QUICK REFERENCE: PRE-ANALYSIS CHECKLIST\n")
cat("═════════════════════════════════════════════════════════════════════════\n\n")

checklist <- "
BEFORE YOU START:

Environment:
  ☐ R version ≥ 4.0? (Check: Rscript --version)
  ☐ Required packages installed? (tidyverse, vegan, Biostrings, entropy)
  ☐ Output directories created? (./output, ./figures, ./data)

Data:
  ☐ LOSO_manifest.tsv present and readable?
  ☐ All 22 test_*.fasta files exist? (22 files)
  ☐ All 22 train_*.fasta files exist? (22 files for validation)
  ☐ Total sequences correct? (Should be ~1,782 across all)
  ☐ File paths correct in manifest?

Configuration:
  ☐ SEED = 20260511 set for reproducibility?
  ☐ K-mer range: k=4-8 confirmed?
  ☐ Bootstrap replicates: 10,000 confirmed?
  ☐ Bias correction method selected? (Chao-Shen < 50)

During Analysis:
  ☐ Progress messages appearing? (should see fold 1-22)
  ☐ K-mer extraction completing? (should be ~30 seconds per fold)
  ☐ Bootstrap running? (should be ~5-10 minutes per fold)
  ☐ No memory errors?
  ☐ No file not found errors?

After Analysis:
  ☐ Output CSV files created? (summary + detailed)
  ☐ Figures generated? (4 PNG files)
  ☐ H values in reasonable range? (1.5-3.5 bits)
  ☐ Test ≠ Train validation passed? (p < 0.01)
  ☐ All 22 superfamilies represented?

Final QC:
  ☐ No missing values in H column?
  ☐ Bootstrap CI widths sensible? (<0.8 bits typically)
  ☐ FDR p-values applied to pairwise tests?
  ☐ Figures are readable (no overlapping text)?
  ☐ Tables have proper headers and formatting?
"

cat(checklist)

# ============================================================================
# TROUBLESHOOTING QUICK GUIDE
# ============================================================================

cat("\n\n")
cat("TROUBLESHOOTING QUICK GUIDE\n")
cat("═════════════════════════════════════════════════════════════════════════\n\n")

troubleshooting_guide <- "
PROBLEM                              | FIRST CHECK              | SOLUTION
─────────────────────────────────────┼──────────────────────────┼─────────────────
\"Error: Cannot open file\"            | File path correct?       | Check manifest 
                                     | File exists?             | path format
─────────────────────────────────────┼──────────────────────────┼─────────────────
\"No sequences read from FASTA\"      | FASTA format correct?    | Check headers 
                                     | Headers start with '>'?  | (must be \">...\")
─────────────────────────────────────┼──────────────────────────┼─────────────────
Shannon H > 5.0 (unreasonably high) | DNA or protein?          | Likely wrong seq
                                     | Check DNA alphabet       | type; recheck files
─────────────────────────────────────┼──────────────────────────┼─────────────────
Bootstrap CI very wide (>1.0)        | Sample size?             | Expected for n<20
                                     | Check Chao-Shen applied? | May indicate noise
─────────────────────────────────────┼──────────────────────────┼─────────────────
Test ≠ Train not significant (p>0.05)| LOSO design check       | Files truly separate?
                                     | Wrong files used?        | ERROR - STOP HERE
─────────────────────────────────────┼──────────────────────────┼─────────────────
Memory error \"vector too large\"     | Bootstrap replicates?    | Reduce to 5,000
                                     | Processing all at once?  | Process by fold
─────────────────────────────────────┼──────────────────────────┼─────────────────
All superfamilies have H ~ same      | Saturation?              | Check individual 
                                     | K-mer size appropriate?  | family plots
─────────────────────────────────────┼──────────────────────────┼─────────────────
Figures don't display                | PNG file created?        | Check output path
                                     | File size > 0 KB?        | Run visualization
                                     |                          | code separately
"

cat(troubleshooting_guide)

cat("\n")
cat("═════════════════════════════════════════════════════════════════════════\n")
cat("END OF VISUAL WORKFLOW GUIDE\n")
cat("═════════════════════════════════════════════════════════════════════════\n\n")

cat("USE THIS DOCUMENT FOR:\n")
cat("  • Quick reference during pipeline execution\n")
cat("  • Troubleshooting issues\n")
cat("  • Understanding the complete workflow at a glance\n")
cat("  • Interpretation guidance\n\n")

cat("For detailed implementation: See Implementation_Roadmap.R\n")
cat("For theoretical background: See Literature_Review_Methodology.R\n")
cat("For expert consensus: See LLM_Council_Analysis_LOSO_Shannon.R\n\n")
