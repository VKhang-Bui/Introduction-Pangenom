#!/usr/bin/env bash

# ==============================================================================
# main.sh
# Master Pangenome Pipeline Orchestrator
# Default Mode: Runs Stage 1 + Loops through ALL loci for Stage 2 & Stage 3
# Single-Locus Mode (-i LOCUS): Runs Stage 1 + Stage 2 & 3 for target locus only
# Output Structure: data/intern/[LOCUS]/{01_input, 02_pggb, 03_nucmer, 04_vcfeval, 05_plots}
# Usage Mode 1 (All Loci):     ./main.sh [RAW_DIR] [PREFIX] [OUT_DIR]
# Usage Mode 2 (Single Locus): ./main.sh -i DRB1-3123
# ==============================================================================

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT_DIR"

PGGB_SRC_DIR="$ROOT_DIR/src/pggb"

# Default values
RAW_DIR="data/raw/HLA"
PREFIX="hla"
OUT_DIR="data/intern"
SINGLE_LOCUS=""

# Parse flags (-i LOCUS, -r RAW_DIR, -p PREFIX, -o OUT_DIR)
while getopts "i:r:p:o:h" opt; do
  case $opt in
    i) SINGLE_LOCUS="$OPTARG" ;;
    r) RAW_DIR="$OPTARG" ;;
    p) PREFIX="$OPTARG" ;;
    o) OUT_DIR="$OPTARG" ;;
    h)
       echo "Usage (Default - All Loci): $0 [-r RAW_DIR] [-p PREFIX] [-o OUT_DIR]"
       echo "Usage (Single Locus Mode):   $0 -i <LOCUS> [-r RAW_DIR] [-p PREFIX] [-o OUT_DIR]"
       exit 0
       ;;
    *) ;;
  esac
done

# Stage 1: Partitioning & Divergence Estimation (Run ONCE)
bash "$PGGB_SRC_DIR/pangenome_partitioning_pipeline.sh" "$RAW_DIR" "$PREFIX" "$OUT_DIR"

DIVERGENCE_TXT="${OUT_DIR}/${PREFIX}_divergence.txt"

if [[ ! -f "$DIVERGENCE_TXT" ]]; then
    echo "Error: Divergence table '$DIVERGENCE_TXT' not found."
    exit 1
fi

# Single-Locus Mode (-i <LOCUS>)
if [[ -n "$SINGLE_LOCUS" ]]; then
    echo ""
    echo "================================================================="
    echo "[SINGLE LOCUS MODE] Processing Locus: $SINGLE_LOCUS"
    echo "================================================================="
    LOCUS_DIR="${OUT_DIR}/${SINGLE_LOCUS}"
    PGGB_OUT_DIR="${LOCUS_DIR}/02_pggb"
    
    bash "$PGGB_SRC_DIR/run_pggb_params.sh" -l "$SINGLE_LOCUS" -o "$PGGB_OUT_DIR" -t "$DIVERGENCE_TXT"
    bash "$PGGB_SRC_DIR/small_variants_eval.sh" -i "$SINGLE_LOCUS" -g "$PGGB_OUT_DIR" -r grch38 -o "$LOCUS_DIR"
    exit 0
fi

# Default Mode: Loop through ALL loci in divergence table
echo ""
echo "================================================================="
echo "[MULTI-LOCI LOOP] Running PGGB & Evaluation for ALL Partitioned Loci"
echo "================================================================="

tail -n +2 "$DIVERGENCE_TXT" | while read -r LOCUS REST; do
    if [[ -z "$LOCUS" ]]; then continue; fi

    LOCUS_DIR="${OUT_DIR}/${LOCUS}"
    PGGB_OUT_DIR="${LOCUS_DIR}/02_pggb"

    echo ""
    echo ">>> Processing Locus: ${LOCUS} <<<"

    bash "$PGGB_SRC_DIR/run_pggb_params.sh" -l "$LOCUS" -o "$PGGB_OUT_DIR" -t "$DIVERGENCE_TXT"
    bash "$PGGB_SRC_DIR/small_variants_eval.sh" -i "$LOCUS" -g "$PGGB_OUT_DIR" -r grch38 -o "$LOCUS_DIR"
done

echo ""
echo "================================================================="
echo "[MASTER WORKFLOW] ALL LOCI EVALUATIONS COMPLETED SUCCESSFULLY!"
echo " Assets stored under: $OUT_DIR/<LOCUS>/"
echo "================================================================="
