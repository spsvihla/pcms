#!/usr/bin/env bash
set -euo pipefail

# ANSI color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# --- Parse inputs ---

SEQS_QZA=""                       # Input QZA filename
TABLE_QZA="query_table.qza"       # Output table QZA filename
DEREP_SEQS_QZA="rep_seqs.qza"     # Output representative sequences filename

usage() {
  echo "Usage: $0 --i-sequences [--o-table query_table.qza] [--o-sequences rep_seqs.qza]"
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --i-sequences) 
      SEQS_QZA="$2" 
      shift 2
      ;;
    --o-table) 
      TABLE_QZA="$2" 
      shift 2
      ;;
    --o-sequences) 
      DEREP_SEQS_QZA="$2" 
      shift 2
      ;;
    -h|--help) 
      usage 
      ;;
    *)
      echo -e "${RED}Invalid option: -$OPTARG${NC}" >&2
      usage
      ;;
  esac
done

if [[ -z "$SEQS_QZA" ]]; then
  echo -e "${RED}Error:${NC} --i-input_qza input.qza is required." >&2
  usage
fi

# --- Perform checks ---

if ! command -v qiime >/dev/null 2>&1; then
  echo -e "${RED}Error:${NC} 'qiime' command not found in PATH. Please install QIIME 2 and ensure it is accessible." >&2
  exit 1
fi

if [[ ! -f "$SEQS_QZA" ]]; then
  echo -e "${RED}Error:${NC} Input file not found: $SEQS_QZA" >&2
  exit 1
fi

# --- Main ---

echo -e "${CYAN}Dereplicating sequences from $SEQS_QZA${NC}"

qiime vsearch dereplicate-sequences \
  --i-sequences "$SEQS_QZA" \
  --o-dereplicated-table "$TABLE_QZA" \
  --o-dereplicated-sequences "$DEREP_SEQS_QZA"

echo -e "${CYAN}Output feature saved as:${NC} $TABLE_QZA"
echo -e "${CYAN}Output representative sequences saved as:${NC} $DEREP_SEQS_QZA"
echo -e "${GREEN}Dereplication complete.${NC}"
