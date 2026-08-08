#!/usr/bin/env bash
# ==============================================================================
# 04_run_pggb.sh: Run PGGB using parameters from divergence summary TXT
#
# Options:
#   -i <COMMUNITY_ID>   Run a single specified community (e.g. -i community_0). Default: all
#   -p <PARTITIONS_DIR> Input partition directory (default: data/intern/partitions_mash)
#   -t <THREADS>        Number of threads (default: 4)
#   -h                  Show help message
#
# Examples:
#   ./04_run_pggb.sh                             # Run all communities in partitions_mash
#   ./04_run_pggb.sh -i community_0              # Run only community_0
#   ./04_run_pggb.sh -i community_5 -t 8         # Run community_5 with 8 threads
# ==============================================================================

set -euo pipefail

TARGET_COMMUNITY="all"
PARTITIONS_DIR="data/intern/partitions_mash"
THREADS="4"

while getopts "i:p:t:h" opt; do
  case $opt in
    i) TARGET_COMMUNITY="$OPTARG" ;;
    p) PARTITIONS_DIR="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    h)
      echo "Usage: $0 [-i community_X] [-p partitions_dir] [-t threads]"
      exit 0
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

DIR_NAME=$(basename "${PARTITIONS_DIR}")
DIVERGENCE_TXT="data/intern/divergence_summary_${DIR_NAME}.txt"
OUT_DIR="data/intern/graphs_${DIR_NAME}"

if [[ ! -f "$DIVERGENCE_TXT" ]]; then
    echo "Error: Divergence summary file '$DIVERGENCE_TXT' not found."
    echo "Please run 03_estimate_divergence.sh first."
    exit 1
fi

TMP_DIR="${OUT_DIR}/tmp"
mkdir -p "${OUT_DIR}" "${TMP_DIR}"
export TMPDIR="${TMP_DIR}"

echo "[1/1] Running PGGB (Target: ${TARGET_COMMUNITY}) in ${PARTITIONS_DIR}..."

# Đọc từng dòng trong file TXT (bỏ dòng tiêu đề) và chạy pggb
tail -n +2 "${DIVERGENCE_TXT}" | while read -r COMM_ID TOTAL_SEQS NUM_SAMPLES MAX_D REC_P REC_S REC_N; do
    if [[ "${TARGET_COMMUNITY}" != "all" && "${COMM_ID}" != "${TARGET_COMMUNITY}" ]]; then
        continue
    fi

    COMM_FASTA="${PARTITIONS_DIR}/${COMM_ID}.fa.gz"
    COMM_OUT="${OUT_DIR}/${COMM_ID}"

    if [[ ! -f "${COMM_FASTA}" ]]; then
        continue
    fi

    echo ">>> Running PGGB for ${COMM_ID} (-p ${REC_P} -s ${REC_S} -n ${REC_N})..."
    
    pggb -i "${COMM_FASTA}" \
         -p "${REC_P}" \
         -s "${REC_S}" \
         -n "${REC_N}" \
         -t "${THREADS}" \
         -o "${COMM_OUT}"
done

echo "Done! PGGB graphs saved in: ${OUT_DIR}"
