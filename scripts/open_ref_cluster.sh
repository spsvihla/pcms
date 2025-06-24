#!/usr/bin/env bash
set -euxo pipefail

# -----------------------------------------------------------------------------
# Closed-reference OTU picking with QIIME 2 + q2-vsearch
#
# Usage:
#   otu_pick.sh \
#     -r ref_seqs.fasta \
#     -q query-seqs.qza \
#     -t query_table.qza \
#     [-i 0.97] \
#     [-p your_prefix] \
#     [-o /path/to/output_dir]
#
# Examples:
#   otu_pick.sh -r gg_13_8_otus/rep_set.fa \
#               -q query-seqs.qza \
#               -t query-table.qza \
#               -i 0.99 \
#               -p myproject \
#               -o results/
# -----------------------------------------------------------------------------

# Default parameters
IDENTITY=0.97
PREFIX="cluster"
OUTDIR="."

# Print help message
usage() {
  sed -n '2,15p' "$0"
  exit 1
}

# Parse options
while getopts ":r:t:q:i:p:o:h" opt; do
  case "$opt" in
    r) REF_FASTA="$OPTARG" ;;
    t) QUERY_TABLE="$OPTARG" ;;
    q) QUERY_QZA="$OPTARG" ;;
    i) IDENTITY="$OPTARG" ;;
    p) PREFIX="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    h) usage ;;
    *)
      echo "Invalid option: -$OPTARG" >&2
      usage
      ;;
  esac
done

# Check required arguments
if [[ -z "${REF_FASTA:-}" ]] || [[ -z "${QUERY_QZA:-}" ]] || [[ -z "${QUERY_TABLE:-}" ]]; then
  echo "Error: -r <ref_seqs.fasta>, -q <query_seqs.qza> and -t <query_table.qza> are required." >&2
  usage
fi

# Ensure QIIME 2 commands are available
command -v qiime >/dev/null 2>&1 || { echo "Error: qiime not found in PATH"; exit 1; }
# The 'vsearch' binary is checked implicitly by qiime plugin; optional:
# command -v vsearch >/dev/null 2>&1 || { echo "Error: vsearch not found in PATH"; exit 1; }

# Create output directory
mkdir -p "$OUTDIR"

# Derive base names and output paths
REF_BASE=$(basename "$REF_FASTA" .fasta)
REF_QZA="$OUTDIR/${REF_BASE}.qza"

TABLE_QZA="$OUTDIR/${PREFIX}_table.qza"
SEQS_QZA="$OUTDIR/${PREFIX}_seqs.qza"
NEW_REFS_QZA="$OUTDIR/${PREFIX}_new_refs.qza"
TABLE_QZV="$OUTDIR/${PREFIX}_table.qzv"

# -----------------------------------------------------------------------------
echo "Importing reference sequences → $REF_QZA ..."
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path "$REF_FASTA" \
  --output-path "$REF_QZA" \
  --input-format DNAFASTAFormat

echo "Open-reference clustering (identity = $IDENTITY) ..."
qiime vsearch cluster-features-open-reference \
  --i-table "$QUERY_TABLE" \
  --i-sequences "$QUERY_QZA" \
  --i-reference-sequences "$REF_QZA" \
  --p-perc-identity "$IDENTITY" \
  --o-clustered-table "$TABLE_QZA" \
  --o-clustered-sequences "$SEQS_QZA" \
  --o-new-reference-sequences "$NEW_REFS_QZA"

echo "Summarizing feature table → $TABLE_QZV ..."
qiime feature-table summarize \
  --i-table "$TABLE_QZA" \
  --o-visualization "$TABLE_QZV"

echo "Done. Outputs are in: $OUTDIR"
