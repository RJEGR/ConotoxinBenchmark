# BLAST Annotation Script — v2 Fix Report

## Test Results Summary

| Metric | ORIGINAL | FIXED v2 | Delta |
|--------|----------|----------|-------|
| chimera | 31 | 31 | 0 |
| **error** | **3** | **0** | **-3 ✓** |
| fragment | 64 | 64 | 0 |
| full | 75 | 75 | 0 |
| multi | 333 | **336** | **+3 ✓** |

## Three Reported Error Sequences — All Fixed

| qseqid | prelim_cat | ORIGINAL | FIXED v2 | hits | subjects | bitscore | %ident | max_coverage |
|--------|-----------|----------|----------|------|----------|----------|--------|-------------|
| comp338_seq0 | m→n | **error** | **multi** | 2 | 2 | 368.0 | 98.6% | 0.972 |
| comp67_seq0 | m→n | **error** | **multi** | 32 | 32 | 420.0 | 96.6% | 0.914 |
| comp336_seq0 | m→n | **error** | **multi** | 2 | 2 | 220.0 | 100.0% | 0.580 |

All three sequences are correctly classified as `m→n` (multiple queries hitting multiple subjects where n_contigs ≠ n_subjects), which resolves to `multi`.

## Cross-tabulation: prelim_cat × final_annotation (v2)

| prelim_cat | chimera | fragment | full | multi | Total |
|-----------|---------|----------|------|-------|-------|
| 1→1 | 0 | 27 | 75 | 0 | 102 |
| 1→n | 31 | 0 | 0 | 0 | 31 |
| m→1 | 0 | 37 | 0 | 0 | 37 |
| m→n | 0 | 0 | 0 | 336 | 336 |
| **Total** | **31** | **64** | **75** | **336** | **506** |

High-confidence errors remaining in v2: **NONE**

---

## Root Cause Analysis

### Bug 1 — `case_when()` evaluates ALL right-hand sides (CRITICAL)

R's `case_when()` is a vectorised function. Even for a scalar `cat_type`, it evaluates **every** RHS formula before selecting the match:

```r
# ORIGINAL (buggy)
annotation <- case_when(
  cat_type == "no hit" ~ "no hit",                            # formula 1
  cat_type == "1->1"   ~ { ... },                             # formula 2
  cat_type == "m->1"   ~ { check_contig_overlaps_on_subject(...) },  # formula 3  ← evaluated even when cat_type ≠ "m->1"
  cat_type == "1->n"   ~ { check_hit_overlaps_on_query(...) },       # formula 4  ← evaluated even when cat_type ≠ "1->n"
  cat_type == "m->n"   ~ { ... },                             # formula 5
  TRUE ~ "other"
)
```

When `cat_type` is `"m->n"`, formulas 3 and 4 are still evaluated. If they call overlap functions that throw, the entire `case_when()` fails and `tryCatch` catches it as `"error"`.

**Error messages explained:**
- `comp338_seq0`: formula 3 (`m->1`) evaluated overlap functions that called `IRanges::pintersect()` on disjoint ranges → threw
- `comp67_seq0`: same as above
- `comp336_seq0`: formula 4 (`1->n`) triggered the same `pintersect()` failure

### Bug 2 — `IRanges::pintersect(resolve.empty = "none")` throws on disjoint ranges

The original `overlap_ratio()` function uses:

```r
ovl <- width(pintersect(ir1, ir2, resolve.empty = "none"))
```

When `ir1` and `ir2` do not overlap, `resolve.empty = "none"` causes `pintersect()` to throw an error instead of returning 0. This is documented IRanges behavior.

---

## Fixes Applied in v2

### Fix 1 — Replace `case_when()` with `if/else if/else`

Extracted annotation logic into `classify_query()`, which uses proper short-circuit `if/else if/else`:

```r
# FIXED (v2)
classify_query <- function(cat_type, hits, blast_df, q_id, ...) {
  if (cat_type == "no hit") return("no hit")
  if (cat_type == "1->1")  return(if (coverage >= threshold) "full" else "fragment")
  if (cat_type == "m->1")  return(if (overlaps) "allele" else "fragment")
  if (cat_type == "1->n")  return(if (overlaps) "multi" else "chimera")
  if (cat_type == "m->n")  { ... return("multi" or "full" or "fragment") }
  return("other")
}
```

Only the matching branch is executed. Non-matching overlap functions are never called.

### Fix 2 — Remove IRanges dependency from `overlap_ratio()`

Replaced with pure base-R arithmetic:

```r
# FIXED (v2) — no IRanges, no pintersect
overlap_ratio <- function(start1, end1, start2, end2) {
  s1 <- pmin(start1, end1);  e1 <- pmax(start1, end1)
  s2 <- pmin(start2, end2);  e2 <- pmax(start2, end2)
  ovl_start <- pmax(s1, s2)
  ovl_end   <- pmin(e1, e2)
  ovl_width <- pmax(0L, ovl_end - ovl_start + 1L)
  shortest  <- pmin(e1 - s1 + 1L, e2 - s2 + 1L)
  return(ovl_width / shortest)
}
```

Disjoint ranges simply return 0 — no exception.

### Fix 3 — Input validation and defensive guards

- `check_contig_overlaps_on_subject()`: guards against `NULL`/`NA` subject IDs, empty data frames, and `NA` coordinates
- `check_hit_overlaps_on_query()`: guards against `NULL` hits and `NA` start/end values
- `filter_blast_results()`: now also filters `NA` and non-positive `length` values
- `assign_preliminary_category()`: returns diagnostic columns (`n_hits`, `n_unique_subjects`, `max_contigs_per_subject`)

### Fix 4 — Enriched output for diagnostics

`annotate_blast_hits()` now returns additional columns: `n_hits`, `n_unique_subjects`, `max_contigs_per_subject`, `best_bitscore`, `best_pident`, `best_max_coverage`.

### Fix 5 — Performance improvement (4× faster)

Removing IRanges construction overhead reduces runtime from ~3.7 s to ~0.9 s on the 9,506-row test file.
