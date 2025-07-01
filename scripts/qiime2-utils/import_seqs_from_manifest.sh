#!/usr/bin/env bash
set -euo pipefail

# ANSI color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# --- Parse inputs --- 

MANIFEST_FILE=""                  # Manifest filename (required)
OUTPUT_QZA="query_seqs.qza"       # Output QZA filename

usage() {
  echo "Usage: $0 --i-manifest manifest.csv [--o-sequences output.qza]"
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --i-manifest) 
      [[ $# -lt 2 ]] && usage
      MANIFEST_FILE="$2" 
      shift 2
      ;;
    --o-sequences) 
      [[ $# -lt 2 ]] && usage
      OUTPUT_QZA="$2" 
      shift 2
      ;;
    -h|--help) 
      usage 
      ;;
    *)
      echo -e "${RED}Invalid option: $1${NC}" >&2
      usage
      ;;
  esac
done

if [[ -z "$MANIFEST_FILE" ]]; then
  echo -e "${RED}Error:${NC} --i-manifest argument is required." >&2
  usage
fi

# --- Perform checks ---

if [[ ! -f "$MANIFEST_FILE" ]]; then
  echo -e "${RED}Error:${NC} Manifest file not found: $MANIFEST_FILE" >&2
  exit 1
fi

if ! command -v qiime >/dev/null 2>&1; then
  echo -e "${RED}Error:${NC} 'qiime' command not found in PATH. Please install QIIME 2 and ensure it is accessible." >&2
  exit 1
fi

# --- Main ---

echo -e "${CYAN}Importing sequences from manifest: $MANIFEST_FILE${NC}"

qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path "$MANIFEST_FILE" \
  --output-path "$OUTPUT_QZA" \
  --input-format SingleEndFastqManifestPhred33

echo -e "${CYAN}Output artifact saved as:${NC} $OUTPUT_QZA"
echo -e "${GREEN}Done importing sequences from manifest.${NC}"
