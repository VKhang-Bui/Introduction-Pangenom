#!/usr/bin/env bash
# ==============================================================================
# 06_eval_and_plot_variants.sh: Variants Evaluation & Plotting with rtg vcfeval
# Directly based on official PGGB Small Variants Evaluation Tutorial.
#
# Usage: ./06_eval_and_plot_variants.sh [-i community_X] [-r ref_prefix] [-t threads]
# ==============================================================================

set -euo pipefail

# Tự động nhận diện đường dẫn Conda Environment (Git-friendly)
if [[ -n "${CONDA_PREFIX:-}" ]]; then
    export PATH="${CONDA_PREFIX}/bin/scripts:${CONDA_PREFIX}/bin:${PATH}"
fi

TARGET_COMM="community_2"
REF_PREFIX="GRCh38"
THREADS="4"

while getopts "i:r:t:h" opt; do
  case $opt in
    i) TARGET_COMM="$OPTARG" ;;
    r) REF_PREFIX="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    h)
      echo "Usage: $0 [-i community_X] [-r ref_prefix] [-t threads]"
      exit 0
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

COMM_NAME=$(basename "${TARGET_COMM}" .fa.gz)
VCF_DIR="data/intern/vcf_calls/${COMM_NAME}"
OUT_DIR="data/intern/vcf_evaluation/${COMM_NAME}"
mkdir -p "${OUT_DIR}"

REF_FASTA="${VCF_DIR}/ref_${REF_PREFIX}.fa"
PG_VCF="${VCF_DIR}/pg_variants.vcf.gz"
MUMMER_VCF="${VCF_DIR}/mummer_baseline.vcf.gz"

if [[ ! -f "$PG_VCF" || ! -f "$MUMMER_VCF" ]]; then
    echo "Lỗi: Không tìm thấy VCF trong ${VCF_DIR}. Hãy chạy 05_call_variants.sh trước."
    exit 1
fi

SDF_DIR="${OUT_DIR}/ref.sdf"
rm -rf "${SDF_DIR}"

echo "[1/3] Chuẩn bị định dạng tham chiếu SDF bằng rtg format..."
rtg format -o "${SDF_DIR}" "${REF_FASTA}"

VCFEVAL_BASE_DIR="${OUT_DIR}/vcfeval"
rm -rf "${VCFEVAL_BASE_DIR}"
mkdir -p "${VCFEVAL_BASE_DIR}"

STAT_TSV="${VCFEVAL_BASE_DIR}/statistics.tsv"
echo -e "contig\tprecision\trecall\tf1.score" > "${STAT_TSV}"

echo "[2/3] Đánh giá PGGB VCF vs MUMmer VCF bằng rtg vcfeval cho từng mẫu..."
# Lấy danh sách các mẫu query trong file VCF PGGB
SAMPLES=$(bcftools query -l "${PG_VCF}" | grep -v "^${REF_PREFIX}$")

for SAMPLE in ${SAMPLES}; do
    SAMPLE_OUT="${VCFEVAL_BASE_DIR}/${SAMPLE}"
    rm -rf "${SAMPLE_OUT}"
    
    rtg vcfeval \
        -t "${SDF_DIR}" \
        -b "${MUMMER_VCF}" \
        -c "${PG_VCF}" \
        --sample "${REF_PREFIX},${SAMPLE}" \
        --sample-ploidy 1 \
        -T "${THREADS}" \
        -o "${SAMPLE_OUT}" >/dev/null 2>&1 || true

    if [[ -f "${SAMPLE_OUT}/summary.txt" ]]; then
        LINE=$(grep "None" "${SAMPLE_OUT}/summary.txt" | tr -s ' ' | cut -f 7,8,9 -d ' ')
        if [[ -n "$LINE" ]]; then
            echo -e "${SAMPLE}\t${LINE}" | tr ' ' '\t' >> "${STAT_TSV}"
        fi
    fi
done

echo "[3/3] Vẽ biểu đồ Precision/Recall/F1 bằng R ggplot2 theo đúng tài liệu PGGB..."
R_SCRIPT="${OUT_DIR}/plot_eval.R"
cat << 'EOF' > "${R_SCRIPT}"
suppressPackageStartupMessages({
    library(ggplot2)
    library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
tsv_file <- args[1]
out_png <- args[2]

stat_df <- read.table(tsv_file, sep = '\t', header = TRUE, comment.char = '?')
if (nrow(stat_df) > 0) {
    stat_df <- pivot_longer(stat_df, precision:f1.score, names_to = "Metric")

    p <- ggplot(stat_df, aes(x = contig, y = value, fill = contig)) +
      geom_bar(stat = "identity") +
      facet_wrap(~Metric, ncol = 1) +
      theme_bw() +
      theme(
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()
      ) +
      theme(legend.position = "none")

    ggsave(out_png, plot = p, width = 7, height = 7, dpi = 300)
}
EOF

OUT_PNG="${OUT_DIR}/MHC.nucmer_vs_pggb.precision_recall_f1score.png"
Rscript "${R_SCRIPT}" "${STAT_TSV}" "${OUT_PNG}"

rm -f "${R_SCRIPT}"

echo "================================================================="
echo "[SUCCESS] Hoàn tất đánh giá rtg vcfeval cho ${COMM_NAME}!"
echo "  • Bảng điểm: ${STAT_TSV}"
echo "  • File ảnh:  ${OUT_PNG}"
echo "================================================================="
