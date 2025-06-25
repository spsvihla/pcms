#!/usr/bin/env bash
set -euo pipefail

# -----------------------------------------------------------------------------
# Dereplicate QIIME 2 sequences using vsearch
#
# This script takes an input QIIME 2 sequences artifact (.qza) and runs
# the vsearch dereplication plugin to produce:
#   - A feature table of dereplicated sequences
#   - A representative sequences artifact
#
# Usage:
#   dereplicate.sh -i input.qza -t table.qza -s rep_seqs.qza
#
# Arguments:
#   -i  Input sequences artifact (.qza) [required]
#   -t  Output feature table artifact (.qza) [optional, default: query_table.qza]
#   -s  Output representative sequences artifact (.qza) [optional, default: rep_seqs.qza]
# -----------------------------------------------------------------------------

# Function to display usage information
usage() {
  sed -n '2,12p' "$0"  # Print lines 2 to 12 (header + usage)
  exit 1
}

# Initialize variables with defaults or empty
INPUT_QZA=""
TABLE_QZA="query_table.qza"
SEQS_QZA="rep_seqs.qza"

# Parse command line options
while getopts ":i:t:s:h" opt; do
  case "$opt" in
    i) INPUT_QZA="$OPTARG" ;;
    t) TABLE_QZA="$OPTARG" ;;
    s) SEQS_QZA="$OPTARG" ;;
    h) usage ;;
    *)
      echo "Invalid option: -$OPTARG" >&2
      usage
      ;;
  esac
done

# Verify the required input argument is provided
if [[ -z "$INPUT_QZA" ]]; then
  echo "Error: -i input.qza is required" >&2
  usage
fi

# Check that the QIIME 2 command line tool is available
command -v qiime >/dev/null 2>&1 || {
  echo "Error: qiime not found in PATH" >&2
  exit 1
}

# Check if the input file actually exists before proceeding
if [[ ! -f "$INPUT_QZA" ]]; then
  echo "Error: Input file not found: $INPUT_QZA" >&2
  exit 1
fi

# Notify user what is going to happen
echo "Dereplicating $INPUT_QZA"
echo "Output feature table: $TABLE_QZA"
echo "Output representative sequences: $SEQS_QZA"

# Run vsearch dereplication through QIIME 2
qiime vsearch dereplicate-sequences \
  --i-sequences "$INPUT_QZA" \
  --o-dereplicated-table "$TABLE_QZA" \
  --o-dereplicated-sequences "$SEQS_QZA"

echo "Dereplication complete."
