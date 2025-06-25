#!/usr/bin/env bash
set -euo pipefail

# -----------------------------------------------------------------------------
# Export a QIIME 2 feature table artifact (.qza) to TSV format for per-sample OTU counts
#
# Usage:
#   export_feature_table.sh -i feature-table.qza -o output.tsv
#
# Arguments:
#   -i FILE     Input QIIME 2 feature table artifact (.qza) [required]
#   -o FILE     Output TSV filename for per-sample OTU counts [required]
# -----------------------------------------------------------------------------

INPUT_QZA=""
OUTPUT_TSV=""

usage() {
  echo "Usage: $0 -i feature-table.qza -o output.tsv"
  echo ""
  echo "Arguments:"
  echo "  -i FILE     Input QIIME 2 feature table artifact (.qza) [required]"
  echo "  -o FILE     Output TSV filename (with path) [required]"
  exit 1
}

while getopts ":i:o:h" opt; do
  case "$opt" in
    i) INPUT_QZA="$OPTARG" ;;
    o) OUTPUT_TSV="$OPTARG" ;;
    h) usage ;;
    *) echo "Invalid option: -$OPTARG" >&2; usage ;;
  esac
done

# Validate inputs
if [[ -z "$INPUT_QZA" ]] || [[ -z "$OUTPUT_TSV" ]]; then
  echo "Error: Both -i input.qza and -o output.tsv are required." >&2
  usage
fi

if [[ ! -f "$INPUT_QZA" ]]; then
  echo "Error: Input file not found: $INPUT_QZA" >&2
  exit 1
fi

command -v qiime >/dev/null 2>&1 || { echo "Error: 'qiime' command not found in PATH." >&2; exit 1; }
command -v biom >/dev/null 2>&1 || { echo "Error: 'biom' command not found in PATH. Install with 'pip install biom-format'." >&2; exit 1; }

# Extract directory for intermediate export files
EXPORT_DIR=$(dirname "$OUTPUT_TSV")
mkdir -p "$EXPORT_DIR"

echo "Exporting feature table from $INPUT_QZA to temporary directory: $EXPORT_DIR"
qiime tools export \
  --input-path "$INPUT_QZA" \
  --output-path "$EXPORT_DIR"

BIOM_FILE="$EXPORT_DIR/feature-table.biom"
if [[ ! -f "$BIOM_FILE" ]]; then
  echo "Error: Exported BIOM file not found: $BIOM_FILE" >&2
  exit 1
fi

echo "Converting BIOM to TSV at $OUTPUT_TSV"
biom convert \
  -i "$BIOM_FILE" \
  -o "$OUTPUT_TSV" \
  --to-tsv

echo "Done. TSV file saved to: $OUTPUT_TSV"
