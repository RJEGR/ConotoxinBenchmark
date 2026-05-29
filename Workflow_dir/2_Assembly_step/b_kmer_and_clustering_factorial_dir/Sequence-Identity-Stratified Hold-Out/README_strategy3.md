# Strategy 3 — Sequence-Identity-Stratified Hold-Out

Repurposes the MMseqs2 clustering workflow from a *descriptive* tool (counting
ConoServer redundancy) into a *held-out experiment* that measures assembler
recovery as a continuous function of how divergent a toxin is from the
reference database.

## What changed vs. the original three scripts

| File | Original behaviour | Strategy 3 behaviour |
|---|---|---|
| `cluster_toxins_mmseqs.sh` | Cluster at many thresholds, **count clusters**, write `sequence_identity,num_clusters` CSV | Cluster at thresholds, **split each cluster** into reference (representative) + one held-out test member, write `reference_<thr>.fasta` / `test_<thr>.fasta`, then `mmseqs search` to compute the **exact nearest-reference-identity (NRI)** per held-out sequence |
| `run_clustering_batch.slurm` | Loop over many input FASTAs, count clusters in each | Take **one** curated ConoServer FASTA, run the stratified splitter, optionally compute a **protein-level NRI** robustness axis, print the simulation/assembly hand-off |
| `analyze_mmseq_clustering_results.R` | Plot **cluster count vs identity threshold** (diversity-reduction curve) | Join assembler **recovery outcomes** to per-sequence NRI, bin by NRI, fit a **logistic recovery curve** per assembler, report `NRI50` (divergence at 50% recovery). Diversity curve kept as a supplementary panel |

## The key conceptual shift

The original workflow's endpoint was *how many clusters does ConoServer have*.
Strategy 3's endpoint is *how does assembler accuracy decay as a test toxin
gets further from anything in the reference set*. Clustering is no longer the
result — it is the **splitting mechanism** that produces train/test pairs with
a known, controllable identity gap.

`NRI` (nearest-reference-identity) is the continuous x-axis. A held-out toxin
with NRI = 88% has a near-twin in the reference; one with NRI = 35% is
genuinely novel. Plotting recovery against NRI shows whether an assembler
*generalises* or merely *recapitulates* sequences it already had neighbours
for — which is the direct, quantitative answer to the StringTie and ConoSorter
circularity critiques.

## Run order

```bash
# 1. Stratified hold-out split + NRI  (HPC job)
sbatch run_clustering_batch.slurm curated_nuc_conoServerDB.fasta
#   add a protein FASTA for the robustness axis:
# sbatch run_clustering_batch.slurm curated_nuc_conoServerDB.fasta curated_prot_conoServerDB.fasta
#
# Output: curated_nuc_conoServerDB_strat3_dir/
#           reference_<thr>.fasta, test_<thr>.fasta
#           nearest_reference_identity.tsv          <- the NRI x-axis
#           cluster_split_<thr>.tsv

# 2. Simulate reads from the held-out test sets only
./Simulate_rnaseq_loso.sh -d curated_nuc_conoServerDB_strat3_dir -p test_

# 3. Assemble each fold (writes <fold>_runlog.tsv)
./Assemblers.sh -s all -c run.config

# 4. Plot recovery vs NRI
Rscript analyze_mmseq_clustering_results.R \
    --strat3_dir curated_nuc_conoServerDB_strat3_dir \
    --runlog_glob '*_runlog.tsv'
#   or, if you already have a recovery table:
# Rscript analyze_mmseq_clustering_results.R \
#     --strat3_dir curated_nuc_conoServerDB_strat3_dir \
#     --recovery_tsv my_recovery.tsv
```

## Outputs of the analysis step

| File | Content |
|---|---|
| `strat3_recovery_vs_NRI.png/.pdf` | **Headline figure**: empirical sensitivity per NRI bin + logistic fit per assembler |
| `strat3_recovery_binned.tsv` | Recovery rate + binomial SE per NRI bin per tool |
| `strat3_logistic_coefficients.tsv` | Logistic slope and `NRI50` per assembler |
| `strat3_analysis_table.tsv` | Per (tool × held-out sequence) recovery joined to NRI |
| `strat3_diversity_curve.png` | Supplementary: the original cluster-count curve |

## Notes / things to verify in your environment

- **NRI is computed at the nucleotide level by default** (consistent with what
  ART simulates and what assemblers reconstruct). Supplying a protein FASTA to
  the SLURM script adds a protein-level NRI as a robustness axis — useful
  because synonymous-site saturation can depress nucleotide identity for a
  hypervariable family.
- **The "recovered" call uses a reference-coverage threshold** (`--id_threshold`,
  default 0.90) on the TransRate `contigs.csv`. TransRate column names vary by
  version; the R script tries `p_unique_ref` / `ref_coverage` / `p_refcov`. If
  your TransRate uses a different column, edit the `refcov_col` line.
- **The recovery table can be supplied directly** via `--recovery_tsv`
  (columns: `test_id, tool, recovered, precision`) if you score recovery with
  your own BLAST-back step instead of TransRate.
- The splitter is **deterministic per seed** — the held-out member chosen from
  each cluster is fixed by `-s SEED` (default 20260520), and the seed is
  perturbed per threshold so different thresholds draw independent held-out
  members.
- Singleton clusters never donate a test sequence (they would empty
  themselves); they are placed in the reference set and recorded with
  `test_id = NA` in `cluster_split_<thr>.tsv`.
