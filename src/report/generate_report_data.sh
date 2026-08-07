#!/usr/bin/env bash
# ==============================================================================
# Script: generate_report_data.sh
# Mục đích: Thực thi các lệnh phân tích chuyên sâu cho Đồ thị Pangenome HLA-DRB1
#           và ghi nhận toàn bộ lệnh thành công & dữ liệu vào thư mục src/report/
# ==============================================================================

set -euo pipefail

BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
REPORT_DIR="${BASE_DIR}/src/report"
OUT_IMG_DIR="${BASE_DIR}/output/image"
PYTHON="/home/vkhang-bui/miniforge3/envs/pggb_env/bin/python3"
BCFTOOLS="/home/vkhang-bui/miniforge3/envs/pggb_env/bin/bcftools"

mkdir -p "${REPORT_DIR}"

COMMANDS_LOG="${REPORT_DIR}/successful_commands.log"
echo "# ==================================================================" > "${COMMANDS_LOG}"
echo "# LỊCH SỬ CÁC CÂU LỆNH PHÂN TÍCH ĐÃ THỰC THI THÀNH CÔNG" >> "${COMMANDS_LOG}"
echo "# Thời gian khởi tạo: $(date)" >> "${COMMANDS_LOG}"
echo "# ==================================================================" >> "${COMMANDS_LOG}"

log_cmd() {
    local title="$1"
    local cmd="$2"
    echo -e "\n[$(date '+%Y-%m-%d %H:%M:%S')] ${title}" >> "${COMMANDS_LOG}"
    echo "CMD: ${cmd}" >> "${COMMANDS_LOG}"
    eval "${cmd}" >> "${COMMANDS_LOG}" 2>&1
    echo "[STATUS] SUCCESS" >> "${COMMANDS_LOG}"
}

echo "[1/3] Trích xuất thống kê VCF bằng bcftools stats..."
log_cmd "Thống kê chi tiết VCF HLA-DRB1" "${BCFTOOLS} stats ${OUT_IMG_DIR}/DRB1-3123.fa.gz.e8de02e.11fba48.89cf512.smooth.final.grch38.vcf.gz > ${REPORT_DIR}/vcf_summary.stats"

echo "[2/3] Thực thi Python parse_metrics.py để trích xuất topology GFA & phân tích siêu dữ liệu ảnh..."
log_cmd "Thực thi Python Script Phân Tích Đồ Thị & Hình Ảnh" "${PYTHON} ${REPORT_DIR}/parse_metrics.py"

echo "[3/3] Xuất danh sách file và thuộc tính kích thước trong output/image/..."
log_cmd "Liệt kê danh sách file kết quả đồ thị và hình ảnh" "ls -lh ${OUT_IMG_DIR} > ${REPORT_DIR}/output_image_files_list.txt"

echo "================================================================="
echo "[HOÀN THÀNH] Toàn bộ lệnh thành công và dữ liệu đã lưu trong src/report/"
