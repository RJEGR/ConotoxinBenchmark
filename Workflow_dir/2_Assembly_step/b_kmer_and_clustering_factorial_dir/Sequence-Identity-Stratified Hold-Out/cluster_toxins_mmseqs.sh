#!/bin/bash
# =============================================================================
# cluster_toxins_mmseqs.sh  (Strategy 3: Sequence-Identity-Stratified Hold-Out)
#
# Repurposes MMseqs2 clustering as a SPLITTING mechanism, not an endpoint.
#
# For each identity threshold:
#   1. Cluster the input FASTA at --min-seq-id <threshold>.
#   2. For every multi-member cluster: representative -> reference set,
#      one randomly chosen non-representative member -> held-out test set.
#      Singletons are placed in the reference set (cannot donate a test seq).
#   3. Write reference_<thr>.fasta and test_<thr>.fasta.
#   4. mmseqs search the held-out test sequences against the reference set
#      to obtain the exact Nearest-Reference-Identity (NRI) per test seq.
#   5. Emit per-threshold NRI tables + a combined NRI table.
#
# The NRI table is the x-axis for Strategy 3: assembler recovery is later
# plotted as a function of NRI.
#
# Usage:
#   ./cluster_toxins_mmseqs.sh -i input.fasta [-t "1.00 0.90 0.80 0.60 0.40"]
#                              [-o outdir] [-p N] [-s SEED]
#     -i  input nucleotide FASTA (curated ConoServer precursors)   [required]
#     -t  space-separated identity thresholds                       [default below]
#     -o  output directory                                          [default: <base>_strat3_dir]
#     -p  threads                                                   [default: 24]
#     -s  random seed for held-out member selection                 [default: 20260520]
#     -h  help
#
# Requires: mmseqs, awk, python3 (only stdlib used).
# =============================================================================

set -euo pipefail

# ---- defaults ----
INPUT_FASTA=""
THRESHOLDS="1.00 0.95 0.90 0.80 0.70 0.60 0.50 0.40"
OUTDIR=""
THREADS=24
SEED=20260520
COV=0.8          # MMseqs2 coverage (-c); kept from original script
EVAL=1e-3

usage() { grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0; }

while getopts "i:t:o:p:s:h" opt; do
    case $opt in
        i) INPUT_FASTA="$OPTARG";;
        t) THRESHOLDS="$OPTARG";;
        o) OUTDIR="$OPTARG";;
        p) THREADS="$OPTARG";;
        s) SEED="$OPTARG";;
        h) usage;;
        ?) echo "Invalid option" >&2; exit 1;;
    esac
done

if [[ -z "$INPUT_FASTA" ]]; then
    echo "Error: -i input.fasta is required. Run with -h for help." >&2
    exit 1
fi
if [[ ! -f "$INPUT_FASTA" ]]; then
    echo "Error: input file not found: $INPUT_FASTA" >&2
    exit 1
fi
command -v mmseqs >/dev/null 2>&1 || { echo "Error: mmseqs not in PATH" >&2; exit 1; }
command -v python3 >/dev/null 2>&1 || { echo "Error: python3 not in PATH" >&2; exit 1; }

INPUT_FASTA="$(cd "$(dirname "$INPUT_FASTA")" && pwd)/$(basename "$INPUT_FASTA")"
BASE_NAME=$(basename "$INPUT_FASTA" | sed 's/\.[^.]*$//')
[[ -z "$OUTDIR" ]] && OUTDIR="${BASE_NAME}_strat3_dir"
mkdir -p "$OUTDIR"
OUTDIR="$(cd "$OUTDIR" && pwd)"

WORK_DIR="$OUTDIR/mmseqs_work"
mkdir -p "$WORK_DIR"

echo "============================================================"
echo "Strategy 3 stratified hold-out splitter"
echo "  input      : $INPUT_FASTA"
echo "  thresholds : $THRESHOLDS"
echo "  outdir     : $OUTDIR"
echo "  threads    : $THREADS   seed: $SEED"
echo "============================================================"

# -----------------------------------------------------------------------------
# Build the MMseqs2 database once (nucleotide; --dbtype 1)
# -----------------------------------------------------------------------------
DB="$WORK_DIR/${BASE_NAME}_db"
echo "[db] creating nucleotide database..."
mmseqs createdb "$INPUT_FASTA" "$DB" --dbtype 1 > "$WORK_DIR/createdb.log" 2>&1

# Combined NRI table header
NRI_ALL="$OUTDIR/nearest_reference_identity.tsv"
printf "threshold\ttest_id\tcluster_id\tnearest_ref_id\tnri_nuc\taln_len\tn_ref\tn_test\n" > "$NRI_ALL"

# Cluster-count summary (kept for backward compatibility with the diversity curve)
CLUSTER_CSV="$OUTDIR/${BASE_NAME}_cluster_counts.csv"
echo "sequence_identity,num_clusters,n_reference,n_test" > "$CLUSTER_CSV"

# -----------------------------------------------------------------------------
# Per-threshold processing
# -----------------------------------------------------------------------------
for thr in $THRESHOLDS; do
    thr_str=$(printf "%.2f" "$thr")
    echo
    echo "=== threshold $thr_str ==="
    tmp="$WORK_DIR/tmp_${thr_str}"
    clu="$WORK_DIR/clu_${thr_str}"
    mkdir -p "$tmp"

    # ---- cluster ----
    echo "[$thr_str] clustering..."
    mmseqs cluster "$DB" "$clu" "$tmp" \
        --min-seq-id "$thr" --threads "$THREADS" \
        --cluster-mode 0 --cov-mode 0 -c "$COV" -e "$EVAL" \
        > "$WORK_DIR/cluster_${thr_str}.log" 2>&1

    # ---- export cluster membership as TSV: rep_id <TAB> member_id ----
    tsv="$WORK_DIR/clu_${thr_str}.tsv"
    mmseqs createtsv "$DB" "$DB" "$clu" "$tsv" > /dev/null 2>&1

    num_clusters=$(cut -f1 "$tsv" | sort -u | wc -l)

    # ---- split each cluster into reference + one held-out test member ----
    # Python: deterministic per-seed selection of one non-rep member per cluster.
    ref_ids="$WORK_DIR/ref_ids_${thr_str}.txt"
    test_ids="$WORK_DIR/test_ids_${thr_str}.txt"
    map_file="$OUTDIR/cluster_split_${thr_str}.tsv"

    python3 - "$tsv" "$ref_ids" "$test_ids" "$map_file" "$SEED" "$thr_str" <<'PYEOF'
import sys, random, collections
tsv, ref_out, test_out, map_out, seed, thr = sys.argv[1:7]
random.seed(int(seed) + int(round(float(thr) * 100)))   # threshold-specific but reproducible

clusters = collections.OrderedDict()
with open(tsv) as fh:
    for line in fh:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 2:
            continue
        rep, member = parts[0], parts[1]
        clusters.setdefault(rep, []).append(member)

ref_ids, test_ids = [], []
with open(map_out, "w") as mh:
    mh.write("cluster_id\trep_id\ttest_id\tcluster_size\n")
    for rep, members in clusters.items():
        # representative always goes to reference
        non_rep = [m for m in members if m != rep]
        if non_rep:
            test = random.choice(non_rep)
            test_ids.append(test)
            # all members except the chosen test sequence -> reference
            for m in members:
                if m != test:
                    ref_ids.append(m)
            mh.write(f"{rep}\t{rep}\t{test}\t{len(members)}\n")
        else:
            # singleton cluster: representative only, no test donation
            ref_ids.extend(members)
            mh.write(f"{rep}\t{rep}\tNA\t{len(members)}\n")

with open(ref_out, "w") as fh:
    fh.write("\n".join(sorted(set(ref_ids))) + ("\n" if ref_ids else ""))
with open(test_out, "w") as fh:
    fh.write("\n".join(sorted(set(test_ids))) + ("\n" if test_ids else ""))

sys.stderr.write(f"  threshold {thr}: {len(set(ref_ids))} reference, "
                 f"{len(set(test_ids))} held-out test sequences\n")
PYEOF

    n_ref=$(wc -l < "$ref_ids" 2>/dev/null || echo 0)
    n_test=$(wc -l < "$test_ids" 2>/dev/null || echo 0)
    echo "[$thr_str] clusters=$num_clusters  reference=$n_ref  test=$n_test"
    echo "$thr_str,$num_clusters,$n_ref,$n_test" >> "$CLUSTER_CSV"

    if [[ "$n_test" -eq 0 ]]; then
        echo "[$thr_str] no multi-member clusters -> no held-out set; skipping NRI."
        rm -rf "$tmp"
        continue
    fi

    # ---- extract reference and test FASTAs by ID ----
    ref_fa="$OUTDIR/reference_${thr_str}.fasta"
    test_fa="$OUTDIR/test_${thr_str}.fasta"

    # mmseqs subdb needs a name->key lookup; simplest robust route is awk over FASTA.
    extract_fasta() {
        local idfile="$1" outfa="$2"
        awk 'BEGIN{while((getline l < ARGV[1])>0) keep[l]=1; ARGV[1]="";
                   p=0}
             /^>/{id=substr($1,2); p=(id in keep)}
             {if(p) print}' "$idfile" "$INPUT_FASTA" > "$outfa"
    }
    extract_fasta "$ref_ids"  "$ref_fa"
    extract_fasta "$test_ids" "$test_fa"

    # ---- exact NRI: search held-out test sequences vs the reference set ----
    echo "[$thr_str] computing nearest-reference-identity (NRI)..."
    ref_db="$WORK_DIR/refdb_${thr_str}"
    test_db="$WORK_DIR/testdb_${thr_str}"
    res_db="$WORK_DIR/searchres_${thr_str}"
    mmseqs createdb "$ref_fa"  "$ref_db"  --dbtype 1 > /dev/null 2>&1
    mmseqs createdb "$test_fa" "$test_db" --dbtype 1 > /dev/null 2>&1

    # search reports fident (fraction identity), alnlen; -s 7.5 = sensitive
    mmseqs search "$test_db" "$ref_db" "$res_db" "$tmp" \
        --search-type 3 -s 7.5 --threads "$THREADS" \
        > "$WORK_DIR/search_${thr_str}.log" 2>&1

    nri_tsv="$WORK_DIR/nri_raw_${thr_str}.tsv"
    mmseqs convertalis "$test_db" "$ref_db" "$res_db" "$nri_tsv" \
        --format-output "query,target,fident,alnlen,evalue,bits" \
        > /dev/null 2>&1

    # For each test seq keep the single best hit (max fident, tie-break on alnlen)
    python3 - "$nri_tsv" "$map_file" "$thr_str" "$n_ref" "$n_test" >> "$NRI_ALL" <<'PYEOF'
import sys, collections
nri_tsv, map_file, thr, n_ref, n_test = sys.argv[1:6]

best = {}
with open(nri_tsv) as fh:
    for line in fh:
        q, t, fident, alnlen, ev, bits = line.rstrip("\n").split("\t")
        fident, alnlen = float(fident), int(alnlen)
        if q not in best or (fident, alnlen) > (best[q][1], best[q][2]):
            best[q] = (t, fident, alnlen)

# map test_id -> cluster_id
clu_of = {}
with open(map_file) as fh:
    next(fh)
    for line in fh:
        cid, rep, test, size = line.rstrip("\n").split("\t")
        if test != "NA":
            clu_of[test] = cid

# emit one row per held-out test sequence (including those with no hit)
test_ids = set(clu_of)
for tid in sorted(test_ids):
    if tid in best:
        t, fident, alnlen = best[tid]
        sys.stdout.write(f"{thr}\t{tid}\t{clu_of[tid]}\t{t}\t{fident:.4f}\t{alnlen}\t{n_ref}\t{n_test}\n")
    else:
        # no detectable hit in the reference set -> NRI ~ 0 (maximally novel)
        sys.stdout.write(f"{thr}\t{tid}\t{clu_of[tid]}\tNA\t0.0000\t0\t{n_ref}\t{n_test}\n")
PYEOF

    rm -rf "$tmp"
    echo "[$thr_str] done."
done

# -----------------------------------------------------------------------------
# Wrap up
# -----------------------------------------------------------------------------
echo
echo "============================================================"
echo "Strategy 3 splitting complete."
echo "  reference/test FASTAs : $OUTDIR/reference_*.fasta , test_*.fasta"
echo "  per-threshold maps    : $OUTDIR/cluster_split_*.tsv"
echo "  NRI table (x-axis)    : $NRI_ALL"
echo "  cluster counts        : $CLUSTER_CSV"
echo "============================================================"
echo
echo "Next: simulate reads from test_*.fasta, assemble, then run"
echo "      analyze_mmseq_clustering_results.R to plot recovery vs NRI."

# keep mmseqs_work? comment the next line out to retain intermediate DBs
rm -rf "$WORK_DIR"
exit 0
