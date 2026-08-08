#!/usr/bin/env bash

# ==============================================================================
# run_pggb_params.sh
# Universal PGGB Graph Runner with Sleek 4-Line Terminal Status Display
# Usage Style 1 (Positional): ./run_pggb_params.sh [LOCUS] [OUT_DIR] [OVERRIDE_S] [TXT_FILE]
# Usage Style 2 (Flags):      ./run_pggb_params.sh -l DRB1-3123 -o data/intern/small_variants/02_pggb -s 500 -t data/intern/hla_divergence.txt
# ==============================================================================

set -euo pipefail

# Default values
LOCUS="DRB1-3123"
OUT_DIR="data/intern/small_variants/02_pggb"
OVERRIDE_S=""
TXT_FILE="data/intern/hla_divergence.txt"
THREADS=8

# Style 1: Positional arguments fallback (if first arg does not start with -)
if [[ $# -gt 0 && "$1" != -* ]]; then
    LOCUS="${1:-$LOCUS}"
    OUT_DIR="${2:-$OUT_DIR}"
    OVERRIDE_S="${3:-$OVERRIDE_S}"
    TXT_FILE="${4:-$TXT_FILE}"
fi

# Style 2: Named flags (-l, -o, -s, -t)
while getopts "l:o:s:t:h" opt; do
  case $opt in
    l) LOCUS="$OPTARG" ;;
    o) OUT_DIR="$OPTARG" ;;
    s) OVERRIDE_S="$OPTARG" ;;
    t) TXT_FILE="$OPTARG" ;;
    h)
       echo "Usage (Flags):      $0 -l <LOCUS> -o <OUT_DIR> [-s OVERRIDE_S] [-t TXT_FILE]"
       echo "Usage (Positional): $0 [LOCUS] [OUT_DIR] [OVERRIDE_S] [TXT_FILE]"
       exit 0
       ;;
    *) ;;
  esac
done

if [[ ! -f "$TXT_FILE" ]]; then
    echo "Error: Parameter file '$TXT_FILE' does not exist."
    exit 1
fi

# Extract parameters from TXT table using pure Bash awk
N_PARAM=$(awk -v locus="$LOCUS" '$1 == locus {print $2}' "$TXT_FILE")
P_PARAM=$(awk -v locus="$LOCUS" '$1 == locus {print $6}' "$TXT_FILE")
S_CALC=$(awk -v locus="$LOCUS" '$1 == locus {print $7}' "$TXT_FILE")

if [[ -z "$P_PARAM" || -z "$S_CALC" || -z "$N_PARAM" ]]; then
    echo "Error: Locus '$LOCUS' not found in '$TXT_FILE'."
    exit 1
fi

S_PARAM="${OVERRIDE_S:-$S_CALC}"

FASTA_INPUT=$(ls data/intern/partitions/*_"${LOCUS}".fasta.gz data/intern/partitions/*"${LOCUS}"*.fasta.gz data/raw/*/*"${LOCUS}"*.fa 2>/dev/null | head -n 1 || true)

if [[ -z "$FASTA_INPUT" ]]; then
    echo "Error: Could not find any FASTA file for locus '$LOCUS' in data/intern/partitions/ or data/raw/"
    exit 1
fi

rm -rf "$OUT_DIR"
mkdir -p "$OUT_DIR"

PGGB_LOG="$OUT_DIR/pggb_execution.log"

echo "================================================================="
echo "[RUN] PGGB GRAPH CONSTRUCTION"
echo "  • Locus:          $LOCUS"
echo "  • Input FASTA:    $FASTA_INPUT"
echo "  • Output Dir:     $OUT_DIR"
echo "  • Parameters:     -p ${P_PARAM}% | -s ${S_PARAM} bp | -n ${N_PARAM} samples"
echo "  • Execution Log:  $PGGB_LOG"
echo "================================================================="

# Start PGGB in background and redirect log
pggb -i "$FASTA_INPUT" \
     -p "$P_PARAM" \
     -s "$S_PARAM" \
     -n "$N_PARAM" \
     -t "$THREADS" \
     -o "$OUT_DIR" > "$PGGB_LOG" 2>&1 &
PGGB_PID=$!

status_wfmash="[Running... ⏳]"
status_seqwish="[Waiting...]"
status_smoothxg="[Waiting...]"
status_odgi="[Waiting...]"

print_status() {
    echo -e " 1/4  wfmash:   $status_wfmash"
    echo -e " 2/4  seqwish:  $status_seqwish"
    echo -e " 3/4  smoothxg: $status_smoothxg"
    echo -e " 4/4  odgi:     $status_odgi"
}

# Initial print
print_status

while kill -0 $PGGB_PID 2>/dev/null; do
    sleep 0.5

    if grep -q -i "seqwish" "$PGGB_LOG" 2>/dev/null; then
        status_wfmash="[Completed ✔]"
        status_seqwish="[Running... ⏳]"
    fi

    if grep -q -i "smoothxg" "$PGGB_LOG" 2>/dev/null; then
        status_seqwish="[Completed ✔]"
        status_smoothxg="[Running... ⏳]"
    fi

    if grep -q -i "odgi" "$PGGB_LOG" 2>/dev/null; then
        status_smoothxg="[Completed ✔]"
        status_odgi="[Running... ⏳]"
    fi

    # Refresh 4-line terminal block
    echo -ne "\033[4A"
    print_status
done

wait $PGGB_PID || true
PGGB_EXIT_CODE=$?

if [[ $PGGB_EXIT_CODE -eq 0 ]]; then
    status_wfmash="[Completed ✔]"
    status_seqwish="[Completed ✔]"
    status_smoothxg="[Completed ✔]"
    status_odgi="[Completed ✔]"
    echo -ne "\033[4A"
    print_status

    # Dọn dẹp đồ thị trung gian, chỉ giữ lại *.final.gfa và *.final.og
    rm -f "$OUT_DIR"/*.seqwish.gfa "$OUT_DIR"/*.smooth.gfa "$OUT_DIR"/*.fix.gfa

    echo ""
    echo "================================================================="
    echo "[SUCCESS] PGGB Graph Construction Completed!"
    echo "  • Output GFA:  $OUT_DIR/*.final.gfa"
    echo "  • Output OG:   $OUT_DIR/*.final.og"
    echo "  • Detailed Log: $PGGB_LOG"
    echo "================================================================="
else
    echo ""
    echo "================================================================="
    echo "[ERROR] PGGB Execution Failed. Inspect log: $PGGB_LOG"
    echo "================================================================="
    exit $PGGB_EXIT_CODE
fi
