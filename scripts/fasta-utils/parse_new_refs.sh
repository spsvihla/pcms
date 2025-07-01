#!/usr/bin/env bash
set -euo pipefail

# ANSI color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# --- Parse inputs ---

REF_SEQS_FA=""
QUERY_SEQS_FA=""
NEW_REFS_FA=""
OUTDIR="."

usage() {
  echo "Usage: $0 --i-ref-sequences ref_seqs.fasta --i-query-sequences query_seqs.fasta --o-new-ref-sequences new_refs.fasta"
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --i-ref-sequences) 
      [[ $# -lt 2 ]] && usage
      REF_SEQS_FA="$2" 
      shift 2
      ;;
    --i-query-sequences) 
      [[ $# -lt 2 ]] && usage
      QUERY_SEQS_FA="$2" 
      shift 2
      ;;
    --o-new-ref-sequences)
      [[ $# -lt 2 ]] && usage
      NEW_REFS_FA="$2" 
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

if [[ -z "$REF_SEQS_FA" || -z "$QUERY_SEQS_FA" || -z "$NEW_REFS_FA" ]]; then
  echo -e "${RED}Error:${NC} --i-ref-sequences ref_seqs.fasta --i-query-sequences query_seqs.fasta and --o-new-ref-sequences new_refs.fasta are required." >&2
  usage
fi

# --- Perform checks ---

command -v seqtk >/dev/null 2>&1 || {
  echo -e "${RED}Error:${NC} 'seqtk' command not found in PATH. Please install seqtk and ensure it is accessible." >&2
  exit 1
}


# --- Main ---

REF_BASE=$(basename $REF_SEQS_FA)
REF_BASE="${REF_BASE%.fasta}"
REF_BASE="${REF_BASE%.fa}"

QUERY_BASE=$(basename $QUERY_SEQS_FA)
QUERY_BASE="${QUERY_BASE%.fasta}"
QUERY_BASE="${QUERY_BASE%.fa}"

REF_SEQ_IDS="$OUTDIR/${REF_BASE}_ids.txt"
QUERY_SEQ_IDS="$OUTDIR/${QUERY_BASE}_ids.txt"
QUERY_SEQ_IDS_UNIQUE="$OUTDIR/${QUERY_BASE}_ids_new.txt"

grep '^>' "$REF_SEQS_FA" | sed 's/^>//' | sort > "$REF_SEQ_IDS"
echo -e "${CYAN}Created $REF_SEQ_IDS with $(wc -l < "$REF_SEQ_IDS") lines${NC}"

grep '^>' "$QUERY_SEQS_FA" | sed 's/^>//' | sort > "$QUERY_SEQ_IDS"
echo -e "${CYAN}Created $QUERY_SEQ_IDS with $(wc -l < "$QUERY_SEQ_IDS") lines${NC}"

comm -23 "$QUERY_SEQ_IDS" "$REF_SEQ_IDS" > "$QUERY_SEQ_IDS_UNIQUE"
echo -e "${CYAN}Created $QUERY_SEQ_IDS_UNIQUE with $(wc -l < "$QUERY_SEQ_IDS_UNIQUE") lines${NC}"

seqtk subseq "$QUERY_SEQS_FA" "$QUERY_SEQ_IDS_UNIQUE" > "$NEW_REFS_FA"
echo -e "${CYAN}Created $NEW_REFS_FA with $(wc -l < "$NEW_REFS_FA") lines${NC}"

echo -e "${CYAN}New reference sequences located in:${NC} $NEW_REFS_FA"
echo -e "${GREEN}Done parsing new reference sequences.${NC}"