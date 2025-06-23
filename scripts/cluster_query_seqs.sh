#!/usr/bin/env bash
set -euxo pipefail

# -----------------------------------------------------------------------------
# Closed‐reference OTU picking with QIIME 2 + q2‐vsearch
#
# Usage:
#   otu_pick.sh \
#     -r ref_seqs.fasta \
#     -q query-seqs.qza \
#     [-i 0.97] \
#     [-p your_prefix] \
#     [-o /path/to/output_dir]
#
# Examples:
#   otu_pick.sh -r gg_13_8_otus/rep_set.fa \
#               -q query-seqs.qza \
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
  sed -n '2,13p' "$0"
  exit 1
}

# Parse options
while getopts ":r:q:i:p:o:h" opt; do
  case "$opt" in
    r) REF_FASTA="$OPTARG" ;;
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
if [[ -z "${REF_FASTA:-}" ]] || [[ -z "${QUERY_QZA:-}" ]]; then
  echo "Error: both -r <ref_seqs.fasta> and -q <query-seqs.qza> are required." >&2
  usage
fi

# Ensure QIIME 2 commands are available
command -v qiime >/dev/null 2>&1 || { echo "qiime not found in PATH"; exit 1; }
command -v qiime-vsearch >/dev/null 2>&1 || { echo "q2-vsearch plugin not available"; exit 1; }

# Create output directory
mkdir -p "$OUTDIR"

# Derive base names and output paths
REF_BASE=$(basename "$REF_FASTA" .fasta)
REF_QZA="$OUTDIR/${REF_BASE}.qza"

TABLE_QZA="$OUTDIR/${PREFIX}-table.qza"
SEQS_QZA="$OUTDIR/${PREFIX}-seqs.qza"
UNMATCHED_QZA="$OUTDIR/${PREFIX}-unmatched.qza"
TABLE_QZV="$OUTDIR/${PREFIX}-table.qzv"

# -----------------------------------------------------------------------------
echo "1) Importing reference sequences → $REF_QZA"
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path "$REF_FASTA" \
  --output-path "$REF_QZA" \
  --input-format DNAFASTAFormat

echo "2) Closed-reference clustering (identity = $IDENTITY)"
qiime vsearch cluster-features-closed-reference \
  --i-sequences "$QUERY_QZA" \
  --i-reference-sequences "$REF_QZA" \
  --p-perc-identity "$IDENTITY" \
  --o-clustered-table "$TABLE_QZA" \
  --o-clustered-sequences "$SEQS_QZA" \
  --o-unmatched-sequences "$UNMATCHED_QZA"

echo "3) Summarizing feature table → $TABLE_QZV"
qiime feature-table summarize \
  --i-table "$TABLE_QZA" \
  --o-visualization "$TABLE_QZV"

echo "Done! Outputs in: $OUTDIR"