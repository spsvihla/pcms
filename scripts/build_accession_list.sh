#!/usr/bin/env bash
set -euo pipefail

# -----------------------------------------------------------------------------
# Download accession list from NCBI given accession range or list
#
# Usage:
#   download_accessions.sh -a accession_range_or_list [-o output_dir]
#
# Example:
#   download_accessions.sh -a "JN427016:JN539989" -o /path/to/output
# -----------------------------------------------------------------------------

ACCESSION=""
OUTPUT_DIR="."

usage() {
  echo "Usage: $0 -a accession_range_or_list [-o output_dir]"
  echo ""
  echo "Options:"
  echo "  -a STRING   Accession range or comma-separated list (e.g. \"JN427016:JN539989\")"
  echo "  -o DIR      Output directory (default: current directory)"
  exit 1
}

while getopts ":a:o:h" opt; do
  case $opt in
    a) ACCESSION="$OPTARG" ;;
    o) OUTPUT_DIR="$OPTARG" ;;
    h) usage ;;
    *) echo "Invalid option: -$OPTARG" >&2; usage ;;
  esac
done

if [[ -z "$ACCESSION" ]]; then
  echo "Error: accession range/list (-a) is required." >&2
  usage
fi

command -v esearch >/dev/null 2>&1 || { echo "Error: esearch not found." >&2; exit 1; }
command -v efetch >/dev/null 2>&1 || { echo "Error: efetch not found." >&2; exit 1; }

mkdir -p "$OUTPUT_DIR"

acc_file_base=$(echo "$ACCESSION" | sed 's/[:,]/_/g')
ACCESSION_LIST="$OUTPUT_DIR/accessions_${acc_file_base}.txt"

echo "Downloading accession list for \"$ACCESSION\" into $ACCESSION_LIST ..."
esearch -db nucleotide -query "${ACCESSION}[Accession]" | efetch -format acc > "$ACCESSION_LIST"

echo "Done. Downloaded $(wc -l < "$ACCESSION_LIST") accession(s)."
