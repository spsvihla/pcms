#!/usr/bin/env bash
set -euxo pipefail

# -----------------------------------------------------------------------------
# Import query sequences into QIIME 2 using a manifest file (multi-sample)
#
# Usage:
#   import_query_manifest.sh [-m manifest.csv] [-o query_seqs.qza]
#
# Examples:
#   import_query_manifest.sh
#   import_query_manifest.sh -m my_manifest.csv -o my_queries.qza
# -----------------------------------------------------------------------------

# Defaults
MANIFEST_FILE="manifest.csv"
QZA_FILE="query_seqs.qza"

# Usage help
usage() {
  sed -n '2,12p' "$0"
  exit 1
}

# Parse options
while getopts ":m:o:h" opt; do
  case "$opt" in
    m) MANIFEST_FILE="$OPTARG" ;;
    o) QZA_FILE="$OPTARG" ;;
    h) usage ;;
    *)
      echo "Invalid option: -$OPTARG" >&2
      usage
      ;;
  esac
done

# Check QIIME 2 availability
command -v qiime >/dev/null 2>&1 || { echo "Error: qiime not found in PATH"; exit 1; }

# Check for manifest file
if [[ ! -f "$MANIFEST_FILE" ]]; then
  echo "Error: Manifest file not found: $MANIFEST_FILE"
  exit 1
fi

# Import query sequences using manifest
echo "Importing sequences using manifest $MANIFEST_FILE → $QZA_FILE"
qiime tools import \
  --type 'SampleData[Sequences]' \
  --input-path "$MANIFEST_FILE" \
  --output-path "$QZA_FILE" \
  --input-format SingleEndFastaManifestPhred33

echo "Query sequences imported successfully into: $QZA_FILE"
