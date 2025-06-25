#!/usr/bin/env bash
set -euo pipefail

# -----------------------------------------------------------------------------
# Import a set of per-sample FASTQ sequences into QIIME 2 as SampleData[SequencesWithQuality]
# using a manifest file. The manifest CSV must have columns:
# sample-id, absolute-filepath, direction
#
# Usage:
#   ./import_query_seqs_from_manifest.sh -m manifest.csv [-o output.qza]
#
# Options:
#   -m manifest.csv       Path to manifest CSV file (required)
#   -o output.qza         Path to output QZA file (default: query_seqs.qza)
#
# This imports sequences with quality scores. FASTQ files must use Phred33 encoding.
# -----------------------------------------------------------------------------

MANIFEST_FILE=""
OUTPUT_QZA="query_seqs.qza"

usage() {
  echo "Usage: $0 -m manifest.csv [-o output.qza]"
  echo ""
  echo "Options:"
  echo "  -m manifest.csv     Required manifest file with per-sample FASTQ info"
  echo "  -o output.qza       Output path for .qza file (default: query_seqs.qza)"
  exit 1
}

# Parse command-line options
while getopts ":m:o:h" opt; do
  case "$opt" in
    m) MANIFEST_FILE="$OPTARG" ;;
    o) OUTPUT_QZA="$OPTARG" ;;
    h) usage ;;
    *) echo "Invalid option: -$OPTARG" >&2; usage ;;
  esac
done

# Check required arguments
if [[ -z "$MANIFEST_FILE" ]]; then
  echo "Error: You must provide a manifest file with -m"
  usage
fi

if [[ ! -f "$MANIFEST_FILE" ]]; then
  echo "Error: Manifest file not found: $MANIFEST_FILE"
  exit 1
fi

# Check qiime is installed
command -v qiime >/dev/null 2>&1 || {
  echo "Error: 'qiime' not found in PATH"
  exit 1
}

echo "Importing sequences from manifest: $MANIFEST_FILE"
echo "Output artifact: $OUTPUT_QZA"

# Import using QIIME 2 CLI
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path "$MANIFEST_FILE" \
  --output-path "$OUTPUT_QZA" \
  --input-format SingleEndFastqManifestPhred33

echo "Import complete. Output saved to: $OUTPUT_QZA"
