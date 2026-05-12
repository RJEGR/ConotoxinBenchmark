# Unified De Novo Assembly Benchmark

Merges `kmer.sh` + `subsampling.sh` into one configurable pipeline and replaces
the parallel `Kmer.R` / `Subsampling.R` analyses with a single EDA script.

## What changed vs. the original

| Concern | Before | After |
|---|---|---|
| K-mer sweep | `kmer.sh -s … -k 49` | `benchmark.sh -s … -m kmer -k 49` |
| Subsampling sweep | `subsampling.sh -s … -p 0.5` | `benchmark.sh -s … -m subsampling -p 0.5 -a trinity` |
| Full factorial sweep | not supported | `benchmark.sh -s … -m both -k 49 -p 0.5` |
| Choice of assembler | hard-coded per script | `-a spades\|trinity\|idba\|megahit\|rnabloom` |
| `run_reference_transrate` | duplicated in both scripts | single shared function |
| Output directory naming | inconsistent (parsed by string position) | `<factor>_<assembler>_kVAL_pVAL_dir` (regex-parseable) |
| EDA | two scripts with copy-pasted logic | one script that auto-detects sweep type |

The new directory naming is **also backward-compatible**: `transrate_EDA.R` falls
back to the legacy "split by `_`, take position 5" parser when it sees the old
layout, so existing `2_subsampling_dir/` and `3_kmer_dir/` trees still work.

## End-to-end example

```bash
# 1. K-mer sweep on all manifests, default assembler (spades)
for k in 19 21 25 33 49 55 73; do
    ./benchmark.sh -s all -m kmer -k $k -o kmer_sweep_dir
done

# 2. Subsampling sweep across multiple assemblers
for a in spades trinity idba megahit rnabloom; do
    for p in 0.1 0.2 0.3 0.5 0.7 1.0; do
        ./benchmark.sh -s all -m subsampling -p $p -a $a -o subs_sweep_dir
    done
done

# 3. Full factorial (k x proportion) for spades only
for k in 25 49 73; do
    for p in 0.3 0.5 0.7 1.0; do
        ./benchmark.sh -s all -m both -k $k -p $p -o factorial_dir
    done
done

# 4. EDA — point at whichever sweep root you want to analyze
Rscript transrate_EDA.R \
    --input  subs_sweep_dir \
    --refdir ~/Documents/GitHub/ConotoxinBenchmark/INPUTS/vfolds_resampling_dir \
    --out    subs_sweep_dir/EDA \
    --rcov   0.95 \
    --minlen 200
```

`transrate_EDA.R` writes:

- `benchmark_metrics.tsv` — tidy table with rawcontigs, TP, FP, FN, Ratio,
  Accuracy, Precision, Sensitivity, Fscore for every (vfold × assembler ×
  swept-parameter) cell. This is the unified replacement for
  `Sumbsampling_accuracy.tsv` and the per-kmer metrics tables.
- `benchmark_metrics.png` — Precision / Sensitivity / Accuracy / F-score vs.
  the swept axis (faceted, with mean ± SE).
- `precision_vs_sensitivity.png` — scatter colored by sweep value.
- `alignment_strata.png` — count of contigs in the `<80% / ≥80 / ≥90 / ≥95 /
  100%` reference-coverage strata (the plot from `Subsampling.R`).

## Manifest format (unchanged)

Whitespace-separated, no header, four columns:

```
factor_id    forward.fq[.gz]    reverse.fq[.gz]    reference.fa
```

For the "all" mode (`-s all`), drop one `*.txt` manifest per fold into the
working directory and `benchmark.sh` iterates over them.

## Notes

- `-k` accepts a comma-separated list for the SPAdes `-k` flag, e.g. `-k 21,33,55`.
  The list is preserved in the run-id with dashes (`k21-33-55`) so the EDA can
  still parse it.
- `seqkit` subsampling uses a fixed seed (`123`) so PE1/PE2 stay paired and
  results are reproducible across runs.
- The combined `-m both` mode subsamples first, then assembles with the
  requested k-mer — matching what you'd otherwise do by chaining the two
  original scripts.
