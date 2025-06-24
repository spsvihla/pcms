#!/usr/bin/env bash
set -euxo pipefail

# -----------------------------------------------------------------------------
# Dereplicate QIIME 2 sequences using vsearch
#
# Usage:
#   dereplicate.sh -i input.qza -t table.qza -s sequences.qza
#
# Example:
#   dereplicate.sh -i query_seqs.qza -t query_table.qza -s rep_seqs.qza
# -----------------------------------------------------------------------------

# Help
usage() {
  sed -n '2,12p' "$0"
  exit 1
}

INPUT_QZA=""
TABLE_QZA="query_table.qza"
SEQS_QZA="rep_seqs.qza"

# Parse options
while getopts ":i:t:s:h" opt; do
  case "$opt" in
    i) INPUT_QZA="$OPTARG" ;;
    t) TABLE_QZA="$OPTARG" ;;
    s) SEQS_QZA="$OPTARG" ;;
    h) usage ;;
    *)
      echo "Invalid option: -$OPTARG" >&2
      usage ;;
  esac
done

# Check required argument
if [[ -z "$INPUT_QZA" ]]; then
  echo "Error: -i input.qza is required"
  usage
fi

# Check QIIME 2 availability
command -v qiime >/dev/null 2>&1 || { echo "Error: qiime not found in PATH"; exit 1; }

# Check input file exists
if [[ ! -f "$INPUT_QZA" ]]; then
  echo "Error: Input file not found: $INPUT_QZA"
  exit 1
fi

# Run dereplication
echo "Dereplicating $INPUT_QZA → $TABLE_QZA + $SEQS_QZA"
qiime vsearch dereplicate-sequences \
  --i-sequences "$INPUT_QZA" \
  --o-dereplicated-table "$TABLE_QZA" \
  --o-dereplicated-sequences "$SEQS_QZA"

echo "Dereplication complete."
