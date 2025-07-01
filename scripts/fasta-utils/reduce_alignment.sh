#!/usr/bin/env bash
set -euo pipefail

# ANSI color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# --- Parse inputs ---

REF_SEQS_FA=""
REF_SEQS_REDUCED_FA=""
REF_SEQS_ALIGNED_FA=""
REF_SEQS_ALIGNED_REDUCED_FA=""
TRIMMED_REF_SEQ_IDS_TXT=""
THRESH=0.05

usage() {
  echo "Usage: $0 --i-ref-sequences ref_seqs.fasta --i-ref-alignment ref_seqs_aligned.fasta [--o-ref-sequences-reduced ref_seqs_reduced.fasta] [--o-ref-alignment-reduced ref_seqs_aligned_reduced.fasta]"
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --i-ref-sequences)
      [[ $# -lt 2 ]] && usage
      REF_SEQS_FA="$2"
      shift 2
      ;;
    --i-ref-alignment) 
      [[ $# -lt 2 ]] && usage
      REF_SEQS_ALIGNED_FA="$2" 
      shift 2
      ;;
    --o-ref-sequences-reduced)
      [[ $# -lt 2 ]] && usage 
      REF_SEQS_REDUCED_FA="$2"
      shift 2
      ;;
    --o-ref-alignment-reduced)
      [[ $# -lt 2 ]] && usage
      REF_SEQS_ALIGNED_REDUCED_FA="$2"
      shift 2 
      ;;
    -t|--threshold)
      [[ $# -lt 2 ]] && usage
      THRESH="$2"
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

if [[ -z "$REF_SEQS_FA" || -z "$REF_SEQS_ALIGNED_FA" ]]; then
  echo -e "${RED}Error:${NC} --i-ref-sequences ref_seqs.fasta and --i-ref-alignment ref_seqs_aligned.fasta required." >&2
  usage
fi

if [[ -z "$REF_SEQS_ALIGNED_REDUCED_FA" ]]; then
    REF_SEQS_ALIGNED_REDUCED_FA="${REF_SEQS_ALIGNED_FA%.fasta}"
    REF_SEQS_ALIGNED_REDUCED_FA="${REF_SEQS_ALIGNED_REDUCED_FA%.fa}"
    REF_SEQS_ALIGNED_REDUCED_FA="${REF_SEQS_ALIGNED_REDUCED_FA}_reduced.fasta"
fi

if [[ -z "$REF_SEQS_REDUCED_FA" ]]; then
    REF_SEQS_REDUCED_FA="${REF_SEQS_FA%.fasta}"
    REF_SEQS_REDUCED_FA="${REF_SEQS_FA%.fa}"
    REF_SEQS_REDUCED_FA="${REF_SEQS_REDUCED_FA}_reduced.fasta"
fi

# --- Perform checks ---

command -v trimal >/dev/null 2>&1 || {
  echo -e "${RED}Error:${NC} 'trimal' command not found in PATH. Please install 'trimal' and ensure it is accessible." >&2
  exit 1
}

# --- Main ---

# Trim columns
echo -e "${CYAN}Trimming alignment sequences at a threshold of ${THRESH}${NC}"

set +e
trimal -in "$REF_SEQS_ALIGNED_FA" -out "$REF_SEQS_ALIGNED_REDUCED_FA" -gt "$THRESH"
status=$?
set -e

if [[ $status -ne 0 ]]; then
  echo -e "${RED}Error:${NC} trimal failed with exit code $status" >&2
  exit $status
fi

# Report number of trimmed columns
orig_cols=$(awk 'BEGIN{FS=""} /^[^>]/ && NF > 0 { print NF; exit }' "$REF_SEQS_ALIGNED_FA")
trimmed_cols=$(awk 'BEGIN{FS=""} /^[^>]/ && NF > 0 { print NF; exit }' "$REF_SEQS_ALIGNED_REDUCED_FA")

if ! [[ "$orig_cols" =~ ^[0-9]+$ && "$trimmed_cols" =~ ^[0-9]+$ ]]; then
  echo -e "${RED}Error:${NC} Could not determine sequence lengths. Check your input FASTA files." >&2
  echo -e "${RED}orig_cols=${orig_cols:-unset}, trimmed_cols=${trimmed_cols:-unset}${NC}" >&2
  exit 1
fi

cols_dropped=$((orig_cols - trimmed_cols))

if [[ "$cols_dropped" -eq 0 ]]; then
  echo -e "${CYAN}No columns were trimmed. Alignment remains unchanged.${NC}"
else
  echo -e "${CYAN}Trimmed ${cols_dropped} columns.${NC}"
fi

# Define output file for trimmed sequence IDs
TRIMMED_REF_SEQ_IDS_TXT="${REF_SEQS_ALIGNED_FA%.fasta}"
TRIMMED_REF_SEQ_IDS_TXT="${TRIMMED_REF_SEQ_IDS_TXT%.fa}"
TRIMMED_REF_SEQ_IDS_TXT="${TRIMMED_REF_SEQ_IDS_TXT}_trimmed_ids.txt"

# Find empty sequences (all gaps or N) after trimming, save their IDs
awk '
  BEGIN { RS=">"; ORS="" }
  NR > 1 {
    header = $1
    seq = $0
    sub(/\n.*/, "", header)    # Extract header line only
    sub(/^[^\n]*\n/, "", $0)   # Remove header line from sequence block
    gsub(/\n/, "", $0)         # Remove newlines from sequence
    total = length($0)
    gaps = gsub(/[-Nn?]/, "", $0)  # Count gaps/N/?
    if (gaps == total) {
      print header "\n"
    }
  }
' "$REF_SEQS_ALIGNED_REDUCED_FA" > "$TRIMMED_REF_SEQ_IDS_TXT"

num_empty=$(wc -l < "$TRIMMED_REF_SEQ_IDS_TXT")

if [[ "$num_empty" -eq 0 ]]; then
  echo -e "${CYAN}No empty sequences after trimming.${NC}"
  cp $REF_SEQS_FA $REF_SEQS_REDUCED_FA
else
  echo -e "${CYAN}Found ${num_empty} empty sequences after trimming.${NC}"
  
  # Prune those sequences from the unaligned FASTA
  awk -v drop_ids="$TRIMMED_REF_SEQ_IDS_TXT" '
    BEGIN {
      while ((getline line < drop_ids) > 0) {
        ids[line] = 1
      }
      close(drop_ids)
    }
    /^>/ {
      seq_id = substr($0, 2)
      keep = !(seq_id in ids)
    }
    keep { print }
  ' "$REF_SEQS_FA" > "$REF_SEQS_REDUCED_FA"
fi

echo -e "${CYAN}Reduced reference sequences stored in${NC} ${REF_SEQS_REDUCED_FA}"
echo -e "${CYAN}Reduced alignment stored in:${NC} ${REF_SEQS_ALIGNED_REDUCED_FA}"
echo -e "${CYAN}Trimmed sequence IDs stored in:${NC} ${TRIMMED_REF_SEQ_IDS_TXT}"
echo -e "${GREEN}Done trimming alignment sequences.${NC}"
