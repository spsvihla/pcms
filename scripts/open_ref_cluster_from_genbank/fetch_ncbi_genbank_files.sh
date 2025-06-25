#!/usr/bin/env bash
set -euo pipefail

# -----------------------------------------------------------------------------
# Batch fetch GenBank sequences from NCBI given an accession list file.
#
# This script reads a list of NCBI nucleotide accessions and downloads
# the corresponding GenBank sequences in batches to avoid server overload.
# Each batch is posted to NCBI's E-utilities, then sequences are fetched
# in GenBank format with parts.
#
# Usage:
#   fetch_query_seqs.sh -i accession_list.txt [-o output_fasta] [-b batch_size]
#
# Options:
#   -i    Path to a text file containing one accession per line (required).
#   -o    Output FASTA file path (default: ./query_seqs.fasta).
#   -b    Number of accessions per batch for fetching (default: 100).
#
# The script creates a temporary directory for batching and cleans up after.
# If the output file already exists, fetching is skipped.
#
# Requirements:
#   - NCBI E-utilities: epost, efetch
#   - coreutils: split
#   - realpath utility
# -----------------------------------------------------------------------------

INPUT_FILE=""
OUTPUT_FILE="query_seqs.fasta"
BATCH_SIZE=100

usage() {
  echo "Usage: $0 -i accession_list.txt [-o output_fasta] [-b batch_size]"
  exit 1
}

# Parse command line options
while getopts ":i:o:b:h" opt; do
  case $opt in
    i) INPUT_FILE="$OPTARG" ;;
    o) OUTPUT_FILE="$OPTARG" ;;
    b) BATCH_SIZE="$OPTARG" ;;
    h) usage ;;
    *) echo "Invalid option: -$OPTARG" >&2; usage ;;
  esac
done

# Validate required input file argument
if [[ -z "$INPUT_FILE" ]]; then
  echo "Error: Input accession list file (-i) is required." >&2
  usage
fi

# Check for required commands in PATH
command -v epost >/dev/null 2>&1 || { echo "Error: epost not found." >&2; exit 1; }
command -v efetch >/dev/null 2>&1 || { echo "Error: efetch not found." >&2; exit 1; }
command -v split >/dev/null 2>&1 || { echo "Error: split not found." >&2; exit 1; }
command -v realpath >/dev/null 2>&1 || { echo "Error: realpath not found." >&2; exit 1; }

# Resolve absolute paths for input and output files
INPUT_FILE="$(realpath "$INPUT_FILE")"
OUTPUT_FILE="$(realpath -m "$OUTPUT_FILE")"   # -m allows non-existing path

# Temporary directory for batching
TMP_DIR="$(dirname "$OUTPUT_FILE")/batch_tmp"

# Create output directory if needed
mkdir -p "$(dirname "$OUTPUT_FILE")"

# If output file already exists, skip fetching
if [[ -f "$OUTPUT_FILE" ]]; then
  echo "File $OUTPUT_FILE already exists. Skipping fetch."
  exit 0
fi

# Clean up old temporary directory if it exists
if [[ -d "$TMP_DIR" ]]; then
  echo "Cleaning up existing temporary directory: $TMP_DIR"
  rm -rf "$TMP_DIR"
fi
mkdir "$TMP_DIR"

echo "Splitting '$INPUT_FILE' into batches of $BATCH_SIZE lines each in $TMP_DIR..."
split -l "$BATCH_SIZE" "$INPUT_FILE" "$TMP_DIR/batch_"

count=0
total=$(ls "$TMP_DIR"/batch_* | wc -l)
echo "Starting fetch loop over $total batches..."

# Loop over each batch file and fetch sequences
for batch_file in "$TMP_DIR"/batch_*; do
  count=$((count + 1))
  echo "[$count/$total] Fetching: $(basename "$batch_file")"

  # Post accession batch, then fetch sequences in GenBank with parts format
  if ! cat "$batch_file" | epost -db nucleotide | efetch -format gbwithparts >> "$OUTPUT_FILE"; then
    echo "Warning: Failed to fetch batch: $(basename "$batch_file")" >&2
  fi

  # Sleep 1 second between requests to avoid hammering the server
  sleep 1
done

echo "Cleaning up temporary directory..."
rm -rf "$TMP_DIR"

echo "All done. Fetched sequences saved to: $OUTPUT_FILE"
