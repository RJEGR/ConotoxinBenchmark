#!/usr/bin/env python3
"""
Comprehensive test harness for BLAST_based_annotation_optimized.R
Faithfully replicates both the ORIGINAL (buggy) and FIXED (v2) logic,
runs them on the test data, and compares results.
"""

import pandas as pd
import numpy as np
import sys
import traceback
from collections import OrderedDict
import time

# ============================================================================
# 1. SHARED UTILITIES
# ============================================================================

def overlap_ratio(start1, end1, start2, end2):
    """Base-R style overlap ratio (v2 fix: no IRanges dependency)"""
    if pd.isna(start1) or pd.isna(end1) or pd.isna(start2) or pd.isna(end2):
        return 0
    s1, e1 = min(start1, end1), max(start1, end1)
    s2, e2 = min(start2, end2), max(start2, end2)
    w1 = e1 - s1 + 1
    w2 = e2 - s2 + 1
    if w1 <= 0 or w2 <= 0:
        return 0
    ovl_start = max(s1, s2)
    ovl_end   = min(e1, e2)
    ovl_width = max(0, ovl_end - ovl_start + 1)
    shortest = min(w1, w2)
    return ovl_width / shortest if shortest > 0 else 0


def overlap_ratio_iranges_sim(start1, end1, start2, end2):
    """Simulates IRanges::pintersect(resolve.empty='none') behavior.
    THROWS when intervals are disjoint, exactly as IRanges does."""
    if pd.isna(start1) or pd.isna(end1) or pd.isna(start2) or pd.isna(end2):
        return 0
    s1, e1 = min(start1, end1), max(start1, end1)
    s2, e2 = min(start2, end2), max(start2, end2)
    w1 = e1 - s1 + 1
    w2 = e2 - s2 + 1
    ovl_start = max(s1, s2)
    ovl_end   = min(e1, e2)
    if ovl_end < ovl_start:
        # IRanges::pintersect with resolve.empty="none" throws here!
        raise RuntimeError(
            f"IRanges::pintersect error: ranges [{s1},{e1}] and [{s2},{e2}] "
            f"do not overlap and resolve.empty='none'"
        )
    ovl_width = ovl_end - ovl_start + 1
    shortest = min(w1, w2)
    return ovl_width / shortest if shortest > 0 else 0


def filter_blast_results(blast_df, min_coverage=0.5, min_identity=80):
    df = blast_df.dropna(subset=['qlen', 'slen', 'length']).copy()
    df = df[(df['qlen'] > 0) & (df['slen'] > 0) & (df['length'] > 0)]
    df['query_coverage']   = df['length'] / df['qlen']
    df['subject_coverage'] = df['length'] / df['slen']
    df['max_coverage']     = np.maximum(df['query_coverage'], df['subject_coverage'])
    df = df[(df['max_coverage'] >= min_coverage) & (df['pident'] >= min_identity)]
    return df


def assign_preliminary_category(blast_df):
    if len(blast_df) == 0:
        return pd.DataFrame(columns=['qseqid','prelim_cat','n_hits',
                                      'n_unique_subjects','max_contigs_per_subject'])

    query_stats = blast_df.groupby('qseqid').agg(
        n_hits=('sseqid', 'count'),
        n_unique_subjects=('sseqid', 'nunique')
    ).reset_index()

    subject_stats = blast_df.groupby('sseqid')['qseqid'].nunique().reset_index()
    subject_stats.columns = ['sseqid', 'n_contigs']

    qs_pairs = blast_df[['qseqid', 'sseqid']].drop_duplicates()
    qs_merged = qs_pairs.merge(subject_stats, on='sseqid', how='left')
    max_cps = qs_merged.groupby('qseqid')['n_contigs'].max().reset_index()
    max_cps.columns = ['qseqid', 'max_contigs_per_subject']

    result = query_stats.merge(max_cps, on='qseqid', how='left')

    def assign_cat(row):
        nh  = row['n_hits']
        mcp = row['max_contigs_per_subject']
        if nh == 0:                        return "no hit"
        if nh == 1 and mcp == 1:           return "1->1"
        if nh >  1 and mcp == 1:           return "m->1"
        if nh == 1 and mcp >  1:           return "1->n"
        return "m->n"

    result['prelim_cat'] = result.apply(assign_cat, axis=1)
    return result

# ============================================================================
# 2. OVERLAP CHECKING FUNCTIONS (used by both versions)
# ============================================================================

def check_hit_overlaps_on_query(hits, threshold=0.5, use_iranges=False):
    """1->n case: check if subject-coordinate intervals overlap."""
    if hits is None or len(hits) < 2:
        return False
    intervals = []
    for _, h in hits.iterrows():
        s, e = min(h['sstart'], h['send']), max(h['sstart'], h['send'])
        if not (pd.isna(s) or pd.isna(e)) and s <= e:
            intervals.append((s, e))
    if len(intervals) < 2:
        return False
    fn = overlap_ratio_iranges_sim if use_iranges else overlap_ratio
    for x in range(len(intervals) - 1):
        for y in range(x + 1, len(intervals)):
            ov = fn(intervals[x][0], intervals[x][1],
                    intervals[y][0], intervals[y][1])
            if not pd.isna(ov) and ov > threshold:
                return True
    return False


def check_contig_overlaps_on_subject(blast_df, subject_id, query_id,
                                      threshold=0.5, use_iranges=False):
    """m->1 case: check if different contigs overlap on a shared subject."""
    if pd.isna(subject_id) or subject_id == "":
        return False
    sub = blast_df[blast_df['sseqid'] == subject_id].copy()
    sub = sub[['qseqid','sstart','send']].drop_duplicates()
    sub = sub.dropna(subset=['sstart','send'])
    sub['start'] = np.minimum(sub['sstart'], sub['send'])
    sub['end']   = np.maximum(sub['sstart'], sub['send'])
    sub = sub[sub['start'] <= sub['end']]
    contig_ids = sub['qseqid'].unique()
    if len(contig_ids) < 2:
        return False
    fn = overlap_ratio_iranges_sim if use_iranges else overlap_ratio
    n = len(contig_ids)
    for a in range(n - 1):
        for b in range(a + 1, n):
            i1 = sub[sub['qseqid'] == contig_ids[a]]
            i2 = sub[sub['qseqid'] == contig_ids[b]]
            if len(i1) == 0 or len(i2) == 0:
                continue
            ov = fn(i1.iloc[0]['start'], i1.iloc[0]['end'],
                    i2.iloc[0]['start'], i2.iloc[0]['end'])
            if not pd.isna(ov) and ov > threshold:
                return True
    return False

# ============================================================================
# 3. ORIGINAL (BUGGY) ANNOTATION ENGINE — simulates R case_when
# ============================================================================

def annotate_hits_original(blast_df, overlap_threshold=0.5,
                            coverage_threshold=0.9):
    """Faithfully replicates the ORIGINAL R code including the case_when bug.
    case_when evaluates ALL RHS blocks, so overlap functions are called even
    for non-matching categories, and they may throw via IRanges."""

    blast_df = filter_blast_results(blast_df)
    if len(blast_df) == 0:
        return pd.DataFrame(columns=['qseqid','prelim_cat','final_annotation'])

    prelim = assign_preliminary_category(blast_df)
    results = []

    for _, prow in prelim.iterrows():
        q_id     = prow['qseqid']
        cat_type = prow['prelim_cat']
        hits     = blast_df[blast_df['qseqid'] == q_id]

        try:
            # ---- SIMULATE case_when: evaluate ALL formulas ----
            formula_results = {}
            errors = {}

            # Formula 1: no hit
            formula_results[1] = "no hit"

            # Formula 2: 1->1
            try:
                coverage = hits.iloc[0]['length'] / hits.iloc[0]['slen']
                formula_results[2] = "full" if coverage >= coverage_threshold else "fragment"
            except Exception as e:
                errors[2] = str(e)

            # Formula 3: m->1  (uses IRanges-sim overlap check)
            try:
                subject_id = hits.iloc[0]['sseqid']
                overlaps = check_contig_overlaps_on_subject(
                    blast_df, subject_id, q_id, overlap_threshold,
                    use_iranges=True)
                formula_results[3] = "allele" if overlaps else "fragment"
            except Exception as e:
                errors[3] = str(e)

            # Formula 4: 1->n  (uses IRanges-sim overlap check)
            try:
                overlaps = check_hit_overlaps_on_query(
                    hits, overlap_threshold, use_iranges=True)
                formula_results[4] = "multi" if overlaps else "chimera"
            except Exception as e:
                errors[4] = str(e)

            # Formula 5: m->n
            try:
                nc = hits['qseqid'].nunique()
                ns = hits['sseqid'].nunique()
                if nc == ns:
                    cvec = []
                    for qid in hits['qseqid'].unique():
                        qh = hits[hits['qseqid'] == qid]
                        best = qh.loc[qh['length'].idxmax()]
                        cvec.append(best['length'] / best['slen'])
                    if np.mean([c >= coverage_threshold for c in cvec]) >= 0.5:
                        formula_results[5] = "full"
                    else:
                        formula_results[5] = "fragment"
                else:
                    formula_results[5] = "multi"
            except Exception as e:
                errors[5] = str(e)

            # Now: case_when picks the FIRST matching condition
            # BUT if ANY evaluated formula threw, the whole case_when throws
            formula_map = {
                "no hit": 1, "1->1": 2, "m->1": 3, "1->n": 4, "m->n": 5
            }
            target_formula = formula_map.get(cat_type, 6)

            # Check if ANY formula we evaluated threw an error
            if errors:
                # In R, case_when throws when ANY RHS formula fails
                first_err_formula = min(errors.keys())
                raise RuntimeError(
                    f"Failed to evaluate the right-hand side of formula "
                    f"{first_err_formula}.")

            annotation = formula_results.get(target_formula, "other")

        except Exception as e:
            annotation = "error"

        results.append({
            'qseqid': q_id,
            'prelim_cat': cat_type,
            'final_annotation': annotation
        })

    return pd.DataFrame(results)

# ============================================================================
# 4. FIXED (v2) ANNOTATION ENGINE — if/else, base-R overlap
# ============================================================================

def classify_query_v2(cat_type, hits, blast_df, q_id,
                       overlap_threshold=0.5, coverage_threshold=0.9):
    """Fixed classify_query using if/else (no case_when)."""

    if cat_type == "no hit":
        return "no hit"

    if cat_type == "1->1":
        coverage = hits.iloc[0]['length'] / hits.iloc[0]['slen']
        return "full" if coverage >= coverage_threshold else "fragment"

    if cat_type == "m->1":
        subject_id = hits.iloc[0]['sseqid']
        overlaps = check_contig_overlaps_on_subject(
            blast_df, subject_id, q_id, overlap_threshold, use_iranges=False)
        return "allele" if overlaps else "fragment"

    if cat_type == "1->n":
        overlaps = check_hit_overlaps_on_query(
            hits, overlap_threshold, use_iranges=False)
        return "multi" if overlaps else "chimera"

    if cat_type == "m->n":
        nc = hits['qseqid'].nunique()
        ns = hits['sseqid'].nunique()
        if nc == ns:
            cvec = []
            for qid in hits['qseqid'].unique():
                qh = hits[hits['qseqid'] == qid]
                best = qh.loc[qh['length'].idxmax()]
                cvec.append(best['length'] / best['slen'])
            return "full" if np.mean([c >= coverage_threshold for c in cvec]) >= 0.5 else "fragment"
        else:
            return "multi"

    return "other"


def annotate_hits_v2(blast_df, overlap_threshold=0.5, coverage_threshold=0.9):
    """Fixed annotation engine (v2): if/else + base-R overlap."""

    blast_df = filter_blast_results(blast_df)
    if len(blast_df) == 0:
        return pd.DataFrame(columns=['qseqid','prelim_cat','final_annotation',
                                      'best_bitscore','best_pident','best_max_coverage'])

    prelim = assign_preliminary_category(blast_df)
    results = []

    for _, prow in prelim.iterrows():
        q_id     = prow['qseqid']
        cat_type = prow['prelim_cat']
        hits     = blast_df[blast_df['qseqid'] == q_id]

        # Diagnostics
        best_bs  = hits['bitscore'].max() if len(hits) > 0 else np.nan
        best_pid = hits['pident'].max() if len(hits) > 0 else np.nan
        best_cov = hits['max_coverage'].max() if len(hits) > 0 else np.nan

        try:
            annotation = classify_query_v2(
                cat_type, hits, blast_df, q_id,
                overlap_threshold, coverage_threshold)
        except Exception as e:
            annotation = "error"

        results.append({
            'qseqid': q_id,
            'prelim_cat': cat_type,
            'final_annotation': annotation,
            'n_hits': int(prow['n_hits']),
            'n_unique_subjects': int(prow['n_unique_subjects']),
            'max_contigs_per_subject': int(prow['max_contigs_per_subject']),
            'best_bitscore': best_bs,
            'best_pident': best_pid,
            'best_max_coverage': best_cov
        })

    return pd.DataFrame(results)

# ============================================================================
# 5. TEST RUNNER
# ============================================================================

def run_tests():
    print("=" * 76)
    print("  BLAST_based_annotation_optimized.R — ORIGINAL vs FIXED (v2)")
    print("=" * 76)

    # ---- load data ----
    blast_raw = pd.read_csv(
        '/mnt/project/test_error_classes_blast_annotation.tsv', sep='\t')
    print(f"\nLoaded {len(blast_raw):,} raw BLAST rows")

    filtered = filter_blast_results(blast_raw)
    print(f"After filtering: {len(filtered):,} rows")

    prelim = assign_preliminary_category(filtered)
    print(f"Unique queries: {len(prelim):,}")
    print(f"\nPreliminary category distribution:")
    for cat, cnt in prelim['prelim_cat'].value_counts().items():
        print(f"  {cat:>6s}  {cnt:>5d}")

    # ---- run ORIGINAL ----
    print("\n" + "-" * 76)
    print("  Running ORIGINAL annotation engine (simulates case_when + IRanges)")
    print("-" * 76)
    t0 = time.time()
    orig = annotate_hits_original(blast_raw)
    t_orig = time.time() - t0
    print(f"  Completed in {t_orig:.2f}s")

    # ---- run FIXED (v2) ----
    print("\n" + "-" * 76)
    print("  Running FIXED (v2) annotation engine (if/else + base-R overlap)")
    print("-" * 76)
    t0 = time.time()
    fixed = annotate_hits_v2(blast_raw)
    t_fixed = time.time() - t0
    print(f"  Completed in {t_fixed:.2f}s")

    # ---- compare summaries ----
    print("\n" + "=" * 76)
    print("  ANNOTATION SUMMARY COMPARISON")
    print("=" * 76)

    orig_counts = orig['final_annotation'].value_counts().to_dict()
    fix_counts  = fixed['final_annotation'].value_counts().to_dict()
    all_cats = sorted(set(list(orig_counts.keys()) + list(fix_counts.keys())))

    print(f"\n  {'Category':<14s}  {'ORIGINAL':>10s}  {'FIXED v2':>10s}  {'Delta':>8s}")
    print(f"  {'-'*14}  {'-'*10}  {'-'*10}  {'-'*8}")
    for cat in all_cats:
        o = orig_counts.get(cat, 0)
        f = fix_counts.get(cat, 0)
        d = f - o
        flag = " <-- FIXED" if cat == "error" and d < 0 else ""
        print(f"  {cat:<14s}  {o:>10d}  {f:>10d}  {d:>+8d}{flag}")

    # ---- detailed diff: sequences that changed ----
    print("\n" + "=" * 76)
    print("  SEQUENCES THAT CHANGED ANNOTATION")
    print("=" * 76)

    merged = orig[['qseqid','prelim_cat','final_annotation']].rename(
        columns={'final_annotation': 'orig_annotation'}).merge(
        fixed[['qseqid','final_annotation','n_hits','n_unique_subjects',
               'best_bitscore','best_pident','best_max_coverage']].rename(
            columns={'final_annotation': 'v2_annotation'}),
        on='qseqid', how='outer')

    changed = merged[merged['orig_annotation'] != merged['v2_annotation']]
    print(f"\n  Total changed: {len(changed)} / {len(merged)} sequences\n")

    if len(changed) > 0:
        print(f"  {'qseqid':<20s}  {'prelim':>6s}  {'ORIGINAL':>10s}  {'FIXED v2':>10s}  "
              f"{'hits':>5s}  {'subjs':>5s}  {'bitscore':>9s}  {'%ident':>7s}  {'maxcov':>7s}")
        print(f"  {'-'*20}  {'-'*6}  {'-'*10}  {'-'*10}  "
              f"{'-'*5}  {'-'*5}  {'-'*9}  {'-'*7}  {'-'*7}")

        for _, r in changed.sort_values('qseqid').iterrows():
            print(f"  {r['qseqid']:<20s}  {r['prelim_cat']:>6s}  "
                  f"{r['orig_annotation']:>10s}  {r['v2_annotation']:>10s}  "
                  f"{r.get('n_hits','?'):>5}  {r.get('n_unique_subjects','?'):>5}  "
                  f"{r.get('best_bitscore',0):>9.1f}  "
                  f"{r.get('best_pident',0):>7.1f}  "
                  f"{r.get('best_max_coverage',0):>7.3f}")

    # ---- focus on the three reported sequences ----
    print("\n" + "=" * 76)
    print("  FOCUS: THREE REPORTED ERROR SEQUENCES")
    print("=" * 76)

    focus_ids = ["comp338_seq0", "comp67_seq0", "comp336_seq0"]
    for sid in focus_ids:
        orig_row = orig[orig['qseqid'] == sid]
        fix_row  = fixed[fixed['qseqid'] == sid]
        if len(orig_row) == 0 and len(fix_row) == 0:
            print(f"\n  {sid}: NOT FOUND in results (removed by filtering?)")
            continue

        o_cat  = orig_row.iloc[0]['prelim_cat'] if len(orig_row) > 0 else "?"
        o_ann  = orig_row.iloc[0]['final_annotation'] if len(orig_row) > 0 else "?"
        f_ann  = fix_row.iloc[0]['final_annotation'] if len(fix_row) > 0 else "?"
        f_hits = int(fix_row.iloc[0]['n_hits']) if len(fix_row) > 0 else 0
        f_subj = int(fix_row.iloc[0]['n_unique_subjects']) if len(fix_row) > 0 else 0
        f_bs   = fix_row.iloc[0]['best_bitscore'] if len(fix_row) > 0 else 0
        f_pid  = fix_row.iloc[0]['best_pident'] if len(fix_row) > 0 else 0
        f_cov  = fix_row.iloc[0]['best_max_coverage'] if len(fix_row) > 0 else 0

        status = "FIXED" if o_ann == "error" and f_ann != "error" else \
                 "UNCHANGED" if o_ann == f_ann else "CHANGED"

        print(f"\n  {sid}")
        print(f"    Prelim category  : {o_cat}")
        print(f"    ORIGINAL result  : {o_ann}")
        print(f"    FIXED v2 result  : {f_ann}")
        print(f"    Status           : {status}")
        print(f"    Hits / Subjects  : {f_hits} / {f_subj}")
        print(f"    Best bitscore    : {f_bs:.1f}")
        print(f"    Best % identity  : {f_pid:.1f}")
        print(f"    Best max coverage: {f_cov:.3f}")

    # ---- cross-tabulation: prelim_cat x final_annotation ----
    print("\n" + "=" * 76)
    print("  CROSS-TABULATION: prelim_cat × final_annotation (FIXED v2)")
    print("=" * 76)

    ct = pd.crosstab(fixed['prelim_cat'], fixed['final_annotation'], margins=True)
    print(f"\n{ct.to_string()}")

    # ---- high-confidence errors check ----
    print("\n" + "=" * 76)
    print("  HIGH-CONFIDENCE SEQUENCES STILL ANNOTATED AS 'error' (v2)")
    print("=" * 76)

    errors_v2 = fixed[fixed['final_annotation'] == 'error']
    if len(errors_v2) == 0:
        print("\n  NONE — all errors resolved!")
    else:
        hc_errors = errors_v2[errors_v2['best_max_coverage'] > 0.5]
        print(f"\n  Total errors remaining: {len(errors_v2)}")
        print(f"  High-confidence errors (coverage > 0.5): {len(hc_errors)}")
        if len(hc_errors) > 0:
            for _, r in hc_errors.head(10).iterrows():
                print(f"    {r['qseqid']}: {r['prelim_cat']}, "
                      f"cov={r['best_max_coverage']:.3f}, "
                      f"pid={r['best_pident']:.1f}")

    print("\n" + "=" * 76)
    print("  TEST COMPLETE")
    print("=" * 76)


if __name__ == "__main__":
    run_tests()
