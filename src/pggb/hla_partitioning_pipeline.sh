#!/usr/bin/env bash
# ==============================================================================
# Script: hla_partitioning_pipeline.sh
# Mục đích: Tự động hóa chuẩn hóa PanSN, gộp dữ liệu HLA, DỰNG ĐỒ THỊ MẠNG LƯỚI & PHÂN CỤM LEIDEN
#           (paf2net.py & net2communities.py --plot), bóc tách FASTA và TÍNH BỘ THAM SỐ CỜ -P & -S TỐI ƯU CHO PGGB.
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
for cmd in samtools bgzip wfmash paf2net.py net2communities.py mash sed python3; do
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
total_samples=$(ls -1 "${RAW_DIR}"/*.fa 2>/dev/null | wc -l)

for fa_file in "${RAW_DIR}"/*.fa; do
    if [[ -f "${fa_file}" ]]; then
        sample_name=$(basename "${fa_file}" .fa)
        sample_count=$((sample_count + 1))
        # Hiển thị tiến độ ghi đè trên cùng 1 dòng ở màn hình (dùng \r)
        echo -ne "  -> Tiến độ chuẩn hóa PanSN: Đã xử lý ${sample_count}/${total_samples} mẫu (${sample_name})   \r"
        sed -e "s/^>/>${sample_name}#1#/g" -e 's/ .*//g' "${fa_file}" >> "${COMBINED_FA}"
    fi
done

echo "" # Chuyển dòng sau khi kết thúc vòng lặp
bgzip -f -@ 4 "${COMBINED_FA}"
samtools faidx "${COMBINED_GZ}"
echo "  -> [SUCCESS] Đã chuẩn hóa PanSN, nén bgzip và gộp ${sample_count} mẫu HLA vào file: ${COMBINED_GZ}" | tee -a "${SUMMARY_LOG}"

# 5. GIAI ĐOẠN 2: Dựng Đồ thị Mạng lưới & Phân cụm Cộng đồng Leiden (Mapping Graph & Plotting)
echo "" | tee -a "${SUMMARY_LOG}"
echo "[GIAI ĐOẠN 2] Gióng hàng All-vs-All bằng wfmash -m & Dựng Đồ thị Mạng lưới..." | tee -a "${SUMMARY_LOG}"

wfmash -m "${COMBINED_GZ}" -s 1000 -p 90 -t 8 > "${PAF_FILE}" 2>/dev/null

echo "  -> Chuyển đổi PAF sang cấu trúc Đồ thị Mạng lưới (paf2net.py)..." | tee -a "${SUMMARY_LOG}"
paf2net.py -p "${PAF_FILE}" > /dev/null

echo "  -> Phân cụm thuật toán Leiden & Vẽ Đồ thị Mạng lưới (net2communities.py --plot)..." | tee -a "${SUMMARY_LOG}"
net2communities.py \
    -e "${PAF_FILE}.edges.list.txt" \
    -w "${PAF_FILE}.edges.weights.txt" \
    -n "${PAF_FILE}.vertices.id2name.txt" \
    --output-prefix "${COMMUNITY_PREFIX}" \
    --accurate-detection \
    --plot > /dev/null

# Chuyển đổi sơ đồ PDF sang ảnh PNG hiển thị trực tiếp ở data/intern/
if command -v pdftoppm &> /dev/null; then
    pdftoppm -png -r 150 "${COMMUNITY_PREFIX}.communities.pdf" "${INTERN_DIR}/hla_all_community_temp"
    mv "${INTERN_DIR}/hla_all_community_temp-1.png" "${GRAPH_PLOT_PNG}"
    echo "  -> [PLOT SUCCESS] Đã xuất ảnh Đồ thị Mạng lưới trực quan PNG tại: [data/intern/hla_all_community_plot.png](file://${GRAPH_PLOT_PNG})" | tee -a "${SUMMARY_LOG}"
fi

# 6. GIAI ĐOẠN 3: Phân đoạn tệp FASTA theo từng Cộng đồng / Locus
echo "" | tee -a "${SUMMARY_LOG}"
echo "[GIAI ĐOẠN 3] Bóc tách phân đoạn trình tự (Sequence Partitioning)..." | tee -a "${SUMMARY_LOG}"

mapfile -t LOCI_LIST < <(cut -f 1 -d '#' "${COMBINED_GZ}.fai" | sort | uniq)
total_loci=${#LOCI_LIST[@]}
locus_count=0

for gene in "${LOCI_LIST[@]}"; do
    part_gz="${PART_DIR}/hla_${gene}.fasta.gz"
    mapfile -t seq_ids < <(grep "^${gene}#" "${COMBINED_GZ}.fai" | cut -f 1)
    locus_count=$((locus_count + 1))
    
    # Hiển thị tiến độ ghi đè trên cùng 1 dòng ở màn hình (dùng \r)
    echo -ne "  -> Tiến độ phân đoạn trình tự: Đã bóc tách ${locus_count}/${total_loci} Locus (${gene})   \r"
    
    samtools faidx "${COMBINED_GZ}" "${seq_ids[@]}" | bgzip -f -@ 4 > "${part_gz}"
    samtools faidx "${part_gz}"
done

echo "" # Chuyển dòng sau khi kết thúc vòng lặp
echo "  -> [SUCCESS] Đã bóc tách thành công ${total_loci} Locus phân đoạn vào thư mục: ${PART_DIR}" | tee -a "${SUMMARY_LOG}"

# 7. GIAI ĐOẠN 4: Tính ma trận Mash Triangle & TÍNH BỘ CỜ -P VÀ -S TỐI ƯU CHO PGGB
echo "" | tee -a "${SUMMARY_LOG}"
echo "[GIAI ĐOẠN 4] TỰ ĐỘNG TÍNH BỘ THAM SỐ CỜ CẮT -P VÀ -S TỐI ƯU CHO LỆNH PGGB..." | tee -a "${SUMMARY_LOG}"
echo "---------------------------------------------------------------------------------------" | tee -a "${SUMMARY_LOG}"
echo -e "LOCUS\t\tMEAN_LEN\tMAX_DIV\t\tEXACT_LIMIT\tRECOMMENDED_PGGB_CMD" | tee -a "${SUMMARY_LOG}"
echo "---------------------------------------------------------------------------------------" | tee -a "${SUMMARY_LOG}"

for gene in "${LOCI_LIST[@]}"; do
    part_gz="${PART_DIR}/hla_${gene}.fasta.gz"
    part_fai="${part_gz}.fai"
    dist_txt="${DIST_DIR}/hla_${gene}.mash_triangle.txt"
    seq_count=$(wc -l < "${part_fai}")
    
    if [[ "${seq_count}" -ge 2 ]]; then
        mash triangle "${part_gz}" > "${dist_txt}" 2>/dev/null
        
        # Script Python tính toán bộ cờ -p và -s chuẩn xác
        params=$(python3 -c "
import math

# 1. Đọc phân bố độ dài từ .fai
lengths = [int(line.split()[1]) for line in open('${part_fai}')]
mean_l = sum(lengths) / len(lengths)
min_l = min(lengths)

# -s tối ưu = 25% độ dài trung bình (làm tròn hàng trăm), không vượt quá min_l
calc_s = int(round(mean_l * 0.25 / 100) * 100)
opt_s = max(200, min(calc_s, min_l))

# 2. Đọc max divergence từ mash triangle
max_d = 0.0
with open('${dist_txt}') as f:
    for line in f:
        for token in line.strip().split():
            try:
                v = float(token)
                if 0.0 < v < 0.9:
                    if v > max_d:
                        max_d = v
            except ValueError:
                pass

exact_limit = 100.0 - (max_d * 100.0)
rec_p = max(70, math.floor(exact_limit - 3.0)) if max_d > 0 else 95

print(f'{int(mean_l)}\t{max_d:.6f}\t{exact_limit:.2f}%\t-p {rec_p} -s {opt_s}')
")
        mean_len=$(echo "${params}" | cut -f 1)
        max_div=$(echo "${params}" | cut -f 2)
        exact_p=$(echo "${params}" | cut -f 3)
        pggb_flags=$(echo "${params}" | cut -f 4)

        echo -e "${gene}\t\t${mean_len} bp\t${max_div}\t${exact_p}\t\t${pggb_flags}" | tee -a "${SUMMARY_LOG}"
    else
        echo -e "${gene}\t\tN/A\t\tN/A\t\tN/A\t\t-p 95 -s 1000" | tee -a "${SUMMARY_LOG}"
    fi
done

echo "---------------------------------------------------------------------------------------" | tee -a "${SUMMARY_LOG}"

echo "" | tee -a "${SUMMARY_LOG}"
echo "=================================================================" | tee -a "${SUMMARY_LOG}"
echo "[COMPLETE] HOÀN TẤT TOÀN BỘ QUY TRÌNH DỰNG ĐỒ THỊ & TÍNH THAM SỐ PGGB!" | tee -a "${SUMMARY_LOG}"
echo "  - Tệp FASTA gộp:           [data/intern/hla_all.fasta.gz](file://${COMBINED_GZ})" | tee -a "${SUMMARY_LOG}"
echo "  - Ảnh Đồ thị trực quan PNG: [data/intern/hla_all_community_plot.png](file://${GRAPH_PLOT_PNG})" | tee -a "${SUMMARY_LOG}"
echo "  - Thư mục chứa Đồ thị/Cộng đồng: [data/intern/communities](file://${COMM_DIR})" | tee -a "${SUMMARY_LOG}"
echo "  - Thư mục tệp phân đoạn:   [data/intern/partitions](file://${PART_DIR})" | tee -a "${SUMMARY_LOG}"
echo "  - Thư mục ma trận Mash:    [data/intern/distances](file://${DIST_DIR})" | tee -a "${SUMMARY_LOG}"
echo "  - File nhật ký tổng hợp:   [data/intern/logs/partitioning_summary.log](file://${SUMMARY_LOG})" | tee -a "${SUMMARY_LOG}"
echo "=================================================================" | tee -a "${SUMMARY_LOG}"
