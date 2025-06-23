#!/usr/bin/env bash
set -euxo pipefail

# -----------------------------------------------------------------------------
# Import query sequences into QIIME 2 (with optional fetch fallback)
#
# Usage:
#   import_query.sh [-f query-seqs.fasta] [-o query-seqs.qza] [-s fetch_query_seqs.sh]
#
# Examples:
#   import_query.sh
#   import_query.sh -f my-queries.fasta -o my-queries.qza
#
# If the input FASTA file is missing, the script will run fetch_query_seqs.sh
# (or a custom script via -s) to create it.
# -----------------------------------------------------------------------------

# Defaults
FASTA="query-seqs.fasta"
QZA="query-seqs.qza"
FETCH_SCRIPT="./fetch_query_seqs.sh"

# Usage help
usage() {
  sed -n '2,12p' "$0"
  exit 1
}

# Parse options
while getopts ":f:o:s:h" opt; do
  case "$opt" in
    f) FASTA="$OPTARG" ;;
    o) QZA="$OPTARG" ;;
    h) usage ;;
    *)
      echo "Invalid option: -$OPTARG" >&2
      usage
      ;;
  esac
done

# Check QIIME 2 availability
command -v qiime >/dev/null 2>&1 || { echo "Error: qiime not found in PATH"; exit 1; }

# Check for query-seqs.fasta or run fetch script
if [[ ! -f "$FASTA" ]]; then
  echo "Input FASTA $FASTA not found."
  exit 1
fi

# Import query sequences
echo "Importing $FASTA → $QZA"
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path "$FASTA" \
  --output-path "$QZA" \
  --input-format DNAFASTAFormat

echo "Query sequences imported successfully into: $QZA"
