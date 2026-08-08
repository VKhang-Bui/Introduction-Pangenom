#!/usr/bin/env bash

# ==============================================================================
# run_pggb_params.sh
# Pure Bash PGGB Graph Runner (Extracts parameters directly from TXT file)
# Usage: ./run_pggb_params.sh [LOCUS] [OUT_DIR] [OVERRIDE_S] [TXT_FILE]
# Example: ./run_pggb_params.sh DRB1-3123 data/intern/small_variants/02_pggb 500
# ==============================================================================

set -euo pipefail

LOCUS="${1:-DRB1-3123}"
OUT_DIR="${2:-data/intern/small_variants/02_pggb}"
OVERRIDE_S="${3:-}"
TXT_FILE="${4:-data/intern/hla_divergence.txt}"
THREADS=8

if [[ ! -f "$TXT_FILE" ]]; then
    echo "Error: Parameter file '$TXT_FILE' does not exist."
    exit 1
fi

N_PARAM=$(awk -v locus="$LOCUS" '$1 == locus {print $2}' "$TXT_FILE")
P_PARAM=$(awk -v locus="$LOCUS" '$1 == locus {print $6}' "$TXT_FILE")
S_CALC=$(awk -v locus="$LOCUS" '$1 == locus {print $7}' "$TXT_FILE")

if [[ -z "$P_PARAM" || -z "$S_CALC" || -z "$N_PARAM" ]]; then
    echo "Error: Locus '$LOCUS' not found in '$TXT_FILE'."
    exit 1
fi

S_PARAM="${OVERRIDE_S:-$S_CALC}"

FASTA_INPUT="data/intern/partitions/hla_${LOCUS}.fasta.gz"
if [[ ! -f "$FASTA_INPUT" ]]; then
    FASTA_INPUT="data/intern/partitions/${LOCUS}.fasta.gz"
fi
if [[ ! -f "$FASTA_INPUT" ]]; then
    FASTA_INPUT="data/raw/HLA/${LOCUS}.fa"
fi

echo "================================================================="
echo "[RUN] PGGB GRAPH CONSTRUCTION"
echo "  • Locus:          $LOCUS"
echo "  • Input FASTA:    $FASTA_INPUT"
echo "  • Output Dir:     $OUT_DIR"
echo "  • Parameters:     -p ${P_PARAM}% | -s ${S_PARAM} bp | -n ${N_PARAM} samples"
echo "  • Source Table:   $TXT_FILE"
echo "================================================================="

rm -rf "$OUT_DIR"
mkdir -p "$OUT_DIR"

pggb -i "$FASTA_INPUT" \
     -p "$P_PARAM" \
     -s "$S_PARAM" \
     -n "$N_PARAM" \
     -t "$THREADS" \
     -o "$OUT_DIR"

# Dọn dẹp đồ thị trung gian, chỉ giữ lại *.final.gfa và *.final.og
rm -f "$OUT_DIR"/*.seqwish.gfa "$OUT_DIR"/*.smooth.gfa "$OUT_DIR"/*.fix.gfa
