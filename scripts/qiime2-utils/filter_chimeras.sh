#!/usr/bin/env bash
set -euo pipefail

# ANSI color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# --- Parse inputs ---

INPUT_SEQS_QZA=""
INPUT_TABLE_QZA=""
REF_CHIMERA_FA=""
PREFIX=""

usage() {
  echo "Usage: $0 --i-sequences clustered_seqs.qza --i-table clustered_table.qza -i-ref-chimera-fasta chimera_ref.fasta [-p prefix | --prefix prefix]"
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --i-sequences) 
      [[ $# -lt 2 ]] && usage
      INPUT_SEQS_QZA="$2" 
      shift 2
      ;;
    --i-table) 
      [[ $# -lt 2 ]] && usage
      INPUT_TABLE_QZA="$2" 
      shift 2
      ;;
    --i-ref-chimera-fasta) 
      [[ $# -lt 2 ]] && usage
      REF_CHIMERA_FA="$2" 
      shift 2
      ;;
    -p|--prefix) 
      [[ $# -lt 2 ]] && usage
      PREFIX="$2" 
      shift 2
      ;;
    --o-outdir)
      [[ $# -lt 2 ]] && usage
      OUTDIR="$2"
      shift 2
      ;;
    -h|--help) 
      usage 
      ;;
    *)
      echo -e "${RED}Invalid option: -$1${NC}" >&2
      usage
      ;;
  esac
done

if [[ -z "$INPUT_SEQS_QZA" || -z "$INPUT_TABLE_QZA" || -z "$REF_CHIMERA_FA" || -z "$PREFIX" ]]; then
  echo -e "${RED}Error:${NC} --i-sequences --i-table and --i-ref-chimera-fasta are required." >&2
  usage
fi

# --- Perform checks ---

command -v qiime >/dev/null 2>&1 || {
  echo -e "${RED}Error:${NC} 'qiime' command not found in PATH. Please install QIIME 2 and ensure it is accessible." >&2
  exit 1
}

# --- Main ---

FILTERED_SEQS="${OUTDIR}/${PREFIX}_seqs.qza"
CHIMERA_SEQS="${OUTDIR}/${PREFIX}_chimera_seqs.qza"
CHIMERA_STATS="${OUTDIR}/${PREFIX}_chimera_stats.qza"
FILTERED_TABLE="${OUTDIR}/${PREFIX}_table.qza"

CHIMERA_REF_BASE=$(basename "$REF_CHIMERA_FA")
CHIMERA_REF_BASE="${REF_CHIMERA_FA%.fasta}"
CHIMERA_REF_BASE="${REF_CHIMERA_FA%.fa}"
REF_CHIMERA_QZA="${OUTDIR}/${CHIMERA_REF_BASE}.qza"

echo -e "${CYAN}Importing chimera reference FASTA as QIIME artifact...${NC}"
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path "$REF_CHIMERA_FA" \
  --output-path "$REF_CHIMERA_QZA" \
  --input-format DNAFASTAFormat

echo -e "${CYAN}Detecting chimeras in sequences using qiime vsearch uchime-ref...${NC}"
qiime vsearch uchime-ref \
  --i-sequences "$INPUT_SEQS_QZA" \
  --i-table "$INPUT_TABLE_QZA" \
  --i-reference-sequences "$REF_CHIMERA_QZA" \
  --o-nonchimeras "$FILTERED_SEQS" \
  --o-chimeras "$CHIMERA_SEQS" \
  --o-stats "$CHIMERA_STATS"

echo -e "${CYAN}Filtering feature table to remove chimera sequences...${NC}"
qiime feature-table filter-features \
  --i-table "$INPUT_TABLE_QZA" \
  --m-metadata-file "$CHIMERA_SEQS" \
  --p-exclude-ids \
  --o-filtered-table "$FILTERED_TABLE"

echo -e "${CYAN}Non-chimeric sequences artifact:${NC} $FILTERED_SEQS"
echo -e "${CYAN}Chimeric sequences artifact:${NC} $CHIMERA_SEQS"
echo -e "${CYAN}Chimera stats artifact:${NC} $CHIMERA_STATS"
echo -e "${CYAN}Filtered feature table artifact:${NC} $FILTERED_TABLE"
echo -e "${GREEN}Done filtering chimera.${NC}"
