#!/usr/bin/env bash
# ==============================================================================
# Script: hla_partitioning_pipeline.sh
# Mục đích: Tự động hóa chuẩn hóa PanSN, gộp dữ liệu HLA, DỰNG ĐỒ THỊ MẠNG LƯỚI & PHÂN CỤM LEIDEN
#           (paf2net.py & net2communities.py --plot), bóc tách FASTA và tính ma trận Mash Triangle.
# Đầu vào:  data/raw/HLA/*.fa
# Đầu ra:   data/intern/ (File gộp, sơ đồ Đồ thị Mạng lưới PNG, tệp phân đoạn, ma trận Mash, nhật ký)
# Dự án:   Introduction to Pangenomics & Practice
# ==============================================================================

set -euo pipefail

# 1. Định nghĩa các đường dẫn thư mục
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
RAW_DIR="${BASE_DIR}/data/raw/HLA"
INTERN_DIR="${BASE_DIR}/data/intern"
COMM_DIR="${INTERN_DIR}/communities"
PART_DIR="${INTERN_DIR}/partitions"
DIST_DIR="${INTERN_DIR}/distances"
LOG_DIR="${INTERN_DIR}/logs"

COMBINED_FA="${INTERN_DIR}/hla_all.fasta"
COMBINED_GZ="${INTERN_DIR}/hla_all.fasta.gz"
PAF_FILE="${COMM_DIR}/hla_all.mapping.paf"
COMMUNITY_PREFIX="${COMM_DIR}/hla_all_community"
GRAPH_PLOT_PNG="${INTERN_DIR}/hla_all_community_plot.png"
SUMMARY_LOG="${LOG_DIR}/partitioning_summary.log"

# 2. Khởi tạo thư mục đầu ra
mkdir -p "${INTERN_DIR}" "${COMM_DIR}" "${PART_DIR}" "${DIST_DIR}" "${LOG_DIR}"

echo "=================================================================" | tee "${SUMMARY_LOG}"
echo "[START] QUY TRÌNH DỰNG ĐỒ THỊ MẠNG LƯỚI & PHÂN ĐOẠN TRÌNH TỰ HLA" | tee -a "${SUMMARY_LOG}"
echo "Thời gian khởi chạy: $(date)" | tee -a "${SUMMARY_LOG}"
echo "=================================================================" | tee -a "${SUMMARY_LOG}"

# 3. Kiểm tra công cụ cần thiết
for cmd in samtools bgzip wfmash paf2net.py net2communities.py mash sed; do
    if ! command -v "$cmd" &> /dev/null; then
        echo "[ERROR] Công cụ '$cmd' chưa được cài đặt hoặc chưa active trong môi trường!" | tee -a "${SUMMARY_LOG}"
        exit 1
    fi
done

# 4. GIAI ĐOẠN 1: Chuẩn hóa tiêu đề PanSN & Gộp toàn bộ file FASTA
echo "" | tee -a "${SUMMARY_LOG}"
echo "[GIAI ĐOẠN 1] Chuẩn hóa định dạng PanSN (>sample#1#contig) & Gộp dữ liệu..." | tee -a "${SUMMARY_LOG}"

rm -f "${COMBINED_FA}" "${COMBINED_GZ}" "${COMBINED_GZ}.fai" "${COMBINED_GZ}.gzi"

sample_count=0
for fa_file in "${RAW_DIR}"/*.fa; do
    if [[ -f "${fa_file}" ]]; then
        sample_name=$(basename "${fa_file}" .fa)
        echo "  -> Đang xử lý mẫu: ${sample_name}" | tee -a "${SUMMARY_LOG}"
        sed -e "s/^>/>${sample_name}#1#/g" -e 's/ .*//g' "${fa_file}" >> "${COMBINED_FA}"
        sample_count=$((sample_count + 1))
    fi
done

echo "  -> Đã gộp tổng cộng ${sample_count} mẫu HLA vào file: ${COMBINED_FA}" | tee -a "${SUMMARY_LOG}"
bgzip -f -@ 4 "${COMBINED_FA}"
samtools faidx "${COMBINED_GZ}"
echo "[SUCCESS] File gộp hoàn chỉnh: ${COMBINED_GZ}" | tee -a "${SUMMARY_LOG}"

# 5. GIAI ĐOẠN 2: Dựng Đồ thị Mạng lưới & Phân cụm Cộng đồng Leiden (Mapping Graph & Plotting)
echo "" | tee -a "${SUMMARY_LOG}"
echo "[GIAI ĐOẠN 2] Gióng hàng All-vs-All bằng wfmash -m & Dựng Đồ thị Mạng lưới..." | tee -a "${SUMMARY_LOG}"

wfmash -m "${COMBINED_GZ}" -s 1000 -p 90 -t 8 > "${PAF_FILE}"

echo "  -> Chuyển đổi PAF sang cấu trúc Đồ thị Mạng lưới (paf2net.py)..." | tee -a "${SUMMARY_LOG}"
paf2net.py -p "${PAF_FILE}"

echo "  -> Phân cụm thuật toán Leiden & Vẽ Đồ thị Mạng lưới (net2communities.py --plot)..." | tee -a "${SUMMARY_LOG}"
net2communities.py \
    -e "${PAF_FILE}.edges.list.txt" \
    -w "${PAF_FILE}.edges.weights.txt" \
    -n "${PAF_FILE}.vertices.id2name.txt" \
    --output-prefix "${COMMUNITY_PREFIX}" \
    --accurate-detection \
    --plot | tee -a "${SUMMARY_LOG}"

# Chuyển đổi sơ đồ PDF sang ảnh PNG hiển thị trực tiếp ở data/intern/
if command -v pdftoppm &> /dev/null; then
    pdftoppm -png -r 150 "${COMMUNITY_PREFIX}.communities.pdf" "${INTERN_DIR}/hla_all_community_temp"
    mv "${INTERN_DIR}/hla_all_community_temp-1.png" "${GRAPH_PLOT_PNG}"
    echo "  [PLOT SUCCESS] Đã xuất ảnh Đồ thị Mạng lưới trực quan PNG tại: [data/intern/hla_all_community_plot.png](file://${GRAPH_PLOT_PNG})" | tee -a "${SUMMARY_LOG}"
fi

# 6. GIAI ĐOẠN 3: Phân đoạn tệp FASTA theo từng Cộng đồng / Locus
echo "" | tee -a "${SUMMARY_LOG}"
echo "[GIAI ĐOẠN 3] Bóc tách phân đoạn trình tự (Sequence Partitioning)..." | tee -a "${SUMMARY_LOG}"

mapfile -t LOCI_LIST < <(cut -f 1 -d '#' "${COMBINED_GZ}.fai" | sort | uniq)

for gene in "${LOCI_LIST[@]}"; do
    part_gz="${PART_DIR}/hla_${gene}.fasta.gz"
    mapfile -t seq_ids < <(grep "^${gene}#" "${COMBINED_GZ}.fai" | cut -f 1)
    seq_count=${#seq_ids[@]}
    
    samtools faidx "${COMBINED_GZ}" "${seq_ids[@]}" | bgzip -f -@ 4 > "${part_gz}"
    samtools faidx "${part_gz}"
    echo "  [PARTITION] Locus: ${gene} | Số trình tự: ${seq_count} -> ${part_gz}" | tee -a "${SUMMARY_LOG}"
done

# 7. GIAI ĐOẠN 4: Tính ma trận khoảng cách Mash Triangle (Divergence Estimation)
echo "" | tee -a "${SUMMARY_LOG}"
echo "[GIAI ĐOẠN 4] Tính ma trận khoảng cách di truyền Mash Triangle..." | tee -a "${SUMMARY_LOG}"

for gene in "${LOCI_LIST[@]}"; do
    part_gz="${PART_DIR}/hla_${gene}.fasta.gz"
    dist_txt="${DIST_DIR}/hla_${gene}.mash_triangle.txt"
    seq_count=$(wc -l < "${part_gz}.fai")
    
    if [[ "${seq_count}" -ge 2 ]]; then
        mash triangle "${part_gz}" > "${dist_txt}"
        echo "  [MASH TRIANGLE] Locus: ${gene} (${seq_count} seqs) -> Ma trận đã lưu tại ${dist_txt}" | tee -a "${SUMMARY_LOG}"
    fi
done

echo "" | tee -a "${SUMMARY_LOG}"
echo "=================================================================" | tee -a "${SUMMARY_LOG}"
echo "[COMPLETE] HOÀN TẤT TOÀN BỘ QUY TRÌNH DỰNG ĐỒ THỊ & PHÂN ĐOẠN HLA!" | tee -a "${SUMMARY_LOG}"
echo "  - Tệp FASTA gộp:           [data/intern/hla_all.fasta.gz](file://${COMBINED_GZ})" | tee -a "${SUMMARY_LOG}"
echo "  - Ảnh Đồ thị trực quan PNG: [data/intern/hla_all_community_plot.png](file://${GRAPH_PLOT_PNG})" | tee -a "${SUMMARY_LOG}"
echo "  - Thư mục chứa Đồ thị/Cộng đồng: [data/intern/communities](file://${COMM_DIR})" | tee -a "${SUMMARY_LOG}"
echo "  - Thư mục tệp phân đoạn:   [data/intern/partitions](file://${PART_DIR})" | tee -a "${SUMMARY_LOG}"
echo "  - Thư mục ma trận Mash:    [data/intern/distances](file://${DIST_DIR})" | tee -a "${SUMMARY_LOG}"
echo "  - File nhật ký tổng hợp:   [data/intern/logs/partitioning_summary.log](file://${SUMMARY_LOG})" | tee -a "${SUMMARY_LOG}"
echo "=================================================================" | tee -a "${SUMMARY_LOG}"
