#!/usr/bin/env bash
# ==============================================================================
# main.sh: Master Pipeline Script for PGGB Pangenome Analysis & Evaluation
#
# Options:
#   -i <INPUT_TARGET>   Input FASTA file or community ID (Default: data/intern/hla_all.fasta.gz)
#   -o <OUT_DIR>        Output base directory (Default: data/intern)
#   -t <THREADS>        Number of threads (Default: 4)
#   -h                  Show help message
#
# Examples:
#   ./main.sh                                # Run full pipeline with default settings
#   ./main.sh -i community_2                 # Run pipeline targeting community_2
# ==============================================================================

set -euo pipefail

# Tự động nhận diện đường dẫn Conda Environment (Git-friendly)
if [[ -n "${CONDA_PREFIX:-}" ]]; then
    export PATH="${CONDA_PREFIX}/bin/scripts:${CONDA_PREFIX}/bin:${PATH}"
fi

INPUT_TARGET="data/intern/hla_all.fasta.gz"
OUT_DIR="data/intern"
THREADS="4"

while getopts "i:o:t:h" opt; do
  case $opt in
    i) INPUT_TARGET="$OPTARG" ;;
    o) OUT_DIR="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    h)
      echo "Usage: $0 [-i input_fasta_or_community] [-o out_dir] [-t threads]"
      exit 0
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

echo "================================================================="
echo "        PGGB PANGENOME PIPELINE MASTER EXECUTION                "
echo "================================================================="

PARTITIONS_DIR="${OUT_DIR}/partitions_mash"
DIR_NAME=$(basename "${PARTITIONS_DIR}")
DIVERGENCE_TXT="${OUT_DIR}/divergence_summary_${DIR_NAME}.txt"

# Kiểm tra nếu chưa có kết quả phân vùng & ước tính phân kỳ thì mới chạy 3 bước đầu
if [[ ! -f "${DIVERGENCE_TXT}" || ! -d "${PARTITIONS_DIR}" ]]; then
    echo "[INFO] Chưa có dữ liệu phân vùng. Bắt đầu chạy 3 bước chuẩn bị..."

    # Step 1: Prep PanSN
    echo "[STEP 1] Running 01_prep_pansn.sh..."
    ./src/pggb/01_prep_pansn.sh

    # Step 2: Partition Communities
    echo "[STEP 2] Running 02_partition_communities_mash.sh..."
    ./src/pggb/02_partition_communities_mash.sh "${OUT_DIR}/hla_all.fasta.gz" "${PARTITIONS_DIR}" "${THREADS}"

    # Step 3: Estimate Divergence
    echo "[STEP 3] Running 03_estimate_divergence.sh..."
    ./src/pggb/03_estimate_divergence.sh "${PARTITIONS_DIR}"
else
    echo "[SKIP] Đã có file '${DIVERGENCE_TXT}'. Bỏ qua 3 bước đầu."
fi

# Step 4: Run PGGB Graph Building
echo "[STEP 4] Running 04_run_pggb.sh..."
if [[ "${INPUT_TARGET}" == community_* ]]; then
    ./src/pggb/04_run_pggb.sh -i "${INPUT_TARGET}" -p "${PARTITIONS_DIR}" -t "${THREADS}"
else
    ./src/pggb/04_run_pggb.sh -p "${PARTITIONS_DIR}" -t "${THREADS}"
fi

# Step 5: Call Variants
echo "[STEP 5] Running 05_call_variants.sh..."
if [[ "${INPUT_TARGET}" == community_* ]]; then
    ./src/pggb/05_call_variants.sh -i "${INPUT_TARGET}" -t "${THREADS}"
else
    ./src/pggb/05_call_variants.sh -i community_2 -t "${THREADS}"
fi

# Step 6: Evaluate & Plot Variants
echo "[STEP 6] Running 06_eval_and_plot_variants.sh..."
if [[ "${INPUT_TARGET}" == community_* ]]; then
    ./src/pggb/06_eval_and_plot_variants.sh -i "${INPUT_TARGET}" -t "${THREADS}"
else
    ./src/pggb/06_eval_and_plot_variants.sh -i community_2 -t "${THREADS}"
fi

echo "================================================================="
echo "✔ PIPELINE COMPLETED SUCCESSFULLY!"
echo "================================================================="
