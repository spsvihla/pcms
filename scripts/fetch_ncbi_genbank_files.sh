#!/usr/bin/env bash
set -euo pipefail

# -----------------------------------------------------------------------------
# Batch fetch GenBank sequences from NCBI given an accession list file
#
# Usage:
#   fetch_query_seqs.sh -i accession_list.txt [-o output_fasta] [-b batch_size]
#
# Default output fasta: ./query_seqs.fasta
# -----------------------------------------------------------------------------

INPUT_FILE=""
OUTPUT_FILE="query_seqs.fasta"
BATCH_SIZE=100

usage() {
  echo "Usage: $0 -i accession_list.txt [-o output_fasta] [-b batch_size]"
  exit 1
}

while getopts ":i:o:b:h" opt; do
  case $opt in
    i) INPUT_FILE="$OPTARG" ;;
    o) OUTPUT_FILE="$OPTARG" ;;
    b) BATCH_SIZE="$OPTARG" ;;
    h) usage ;;
    *) echo "Invalid option: -$OPTARG" >&2; usage ;;
  esac
done

if [[ -z "$INPUT_FILE" ]]; then
  echo "Error: Input accession list file (-i) is required." >&2
  usage
fi

command -v epost >/dev/null 2>&1 || { echo "Error: epost not found." >&2; exit 1; }
command -v efetch >/dev/null 2>&1 || { echo "Error: efetch not found." >&2; exit 1; }
command -v split >/dev/null 2>&1 || { echo "Error: split not found." >&2; exit 1; }
command -v realpath >/dev/null 2>&1 || { echo "Error: realpath not found." >&2; exit 1; }

INPUT_FILE="$(realpath "$INPUT_FILE")"
OUTPUT_FILE="$(realpath -m "$OUTPUT_FILE")"   # -m allows non-existing path
TMP_DIR="$(dirname "$OUTPUT_FILE")/batch_tmp"

mkdir -p "$(dirname "$OUTPUT_FILE")"

if [[ -f "$OUTPUT_FILE" ]]; then
  echo "File $OUTPUT_FILE already exists. Skipping fetch."
  exit 0
fi

if [[ -d "$TMP_DIR" ]]; then
  echo "Cleaning up existing temporary directory: $TMP_DIR"
  rm -rf "$TMP_DIR"
fi
mkdir "$TMP_DIR"

echo "Splitting '$INPUT_FILE' into batches of $BATCH_SIZE in $TMP_DIR..."
split -l "$BATCH_SIZE" "$INPUT_FILE" "$TMP_DIR/batch_"

count=0
total=$(ls "$TMP_DIR"/batch_* | wc -l)
echo "Starting fetch loop over $total batches..."

for batch_file in "$TMP_DIR"/batch_*; do
  count=$((count + 1))
  echo "[$count/$total] Fetching: $(basename "$batch_file")"

  if ! cat "$batch_file" | epost -db nucleotide | efetch -format gbwithparts >> "$OUTPUT_FILE"; then
    echo "Warning: Failed to fetch batch: $(basename "$batch_file")" >&2
  fi

  sleep 1
done

echo "Cleaning up temporary directory..."
rm -rf "$TMP_DIR"

echo "All done. Fetched sequences saved to: $OUTPUT_FILE"
