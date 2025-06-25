#!/usr/bin/env bash
set -euo pipefail

# -----------------------------------------------------------------------------
# Download accession list from NCBI given an accession range or list
#
# This script uses NCBI's EDirect tools to search nucleotide sequences
# by accession range or list and fetches the accession numbers into a text file.
#
# Usage:
#   download_accessions.sh -a accession_range_or_list [-o output_dir]
#
# Arguments:
#   -a STRING   Accession range or comma-separated list (e.g. "JN427016:JN539989")
#   -o DIR      Output directory for accession list file (default: current directory)
#
# Example:
#   download_accessions.sh -a "JN427016:JN539989" -o /path/to/output
# -----------------------------------------------------------------------------

# Initialize variables with default values
ACCESSION=""         # Accession range or list input string (required)
OUTPUT_DIR="."       # Output directory default is current directory

# Function to display usage/help information
usage() {
  echo "Usage: $0 -a accession_range_or_list [-o output_dir]"
  echo ""
  echo "Options:"
  echo "  -a STRING   Accession range or comma-separated list (e.g. \"JN427016:JN539989\")"
  echo "  -o DIR      Output directory (default: current directory)"
  exit 1
}

# Parse command-line options using getopts
while getopts ":a:o:h" opt; do
  case $opt in
    a) ACCESSION="$OPTARG" ;;  # Set accession range or list
    o) OUTPUT_DIR="$OPTARG" ;; # Set output directory
    h) usage ;;                # Show usage if -h given
    *)
      echo "Invalid option: -$OPTARG" >&2
      usage
      ;;
  esac
done

# Check that the accession argument was provided
if [[ -z "$ACCESSION" ]]; then
  echo "Error: accession range/list (-a) is required." >&2
  usage
fi

# Check required external commands are available
command -v esearch >/dev/null 2>&1 || {
  echo "Error: esearch not found in PATH. Please install NCBI EDirect." >&2
  exit 1
}
command -v efetch >/dev/null 2>&1 || {
  echo "Error: efetch not found in PATH. Please install NCBI EDirect." >&2
  exit 1
}

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Sanitize accession string for safe filename (replace colon and comma with underscore)
acc_file_base=$(echo "$ACCESSION" | sed 's/[:,]/_/g')

# Construct full path to accession list output file
ACCESSION_LIST="$OUTPUT_DIR/accessions_${acc_file_base}.txt"

echo "Downloading accession list for \"$ACCESSION\" into $ACCESSION_LIST ..."

# Use esearch to search nucleotide database for the accession(s)
# Pipe results to efetch to get accession numbers in plain text format
esearch -db nucleotide -query "${ACCESSION}[Accession]" | efetch -format acc > "$ACCESSION_LIST"

echo "Done. Downloaded $(wc -l < "$ACCESSION_LIST") accession(s)."
