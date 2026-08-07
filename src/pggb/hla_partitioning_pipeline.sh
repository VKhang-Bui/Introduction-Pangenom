#!/usr/bin/env bash
# ==============================================================================
# Script: hla_partitioning_pipeline.sh
# Mục đích: Tự động hóa PanSN, Phân cụm Leiden Giai đoạn 2, Bóc tách FASTA Giai đoạn 3,
#           và TÍNH TOÁN BỘ THAM SỐ CỜ -P VÀ -S TỐI ƯU XUẤT THÀNH FILE JSON.
# Đầu vào:  data/raw/HLA/*.fa
# Đầu ra:   data/intern/ (hla_all.fasta.gz, sơ đồ PNG, partitions/, distances/, hla_loci_params.json)
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

COMBINED_FA="${INTERN_DIR}/hla_all.fasta"
COMBINED_GZ="${INTERN_DIR}/hla_all.fasta.gz"
PAF_FILE="${COMM_DIR}/hla_all.mapping.paf"
COMMUNITY_PREFIX="${COMM_DIR}/hla_all_community"
GRAPH_PLOT_PNG="${INTERN_DIR}/hla_all_community_plot.png"
JSON_PARAMS_FILE="${INTERN_DIR}/hla_loci_params.json"

# 2. Khởi tạo thư mục đầu ra
mkdir -p "${INTERN_DIR}" "${COMM_DIR}" "${PART_DIR}" "${DIST_DIR}"

echo "================================================================="
echo "[START] QUY TRÌNH TIỀN XỬ LÝ & TÍNH THAM SỐ DỰNG ĐỒ THỊ PANGENOME"
echo "================================================================="

# 3. Kiểm tra công cụ cần thiết
for cmd in samtools bgzip wfmash paf2net.py net2communities.py mash sed python3; do
    if ! command -v "$cmd" &> /dev/null; then
        echo "[ERROR] Công cụ '$cmd' chưa được cài đặt hoặc chưa active trong môi trường!"
        exit 1
    fi
done

# 4. GIAI ĐOẠN 1: Chuẩn hóa tiêu đề PanSN & Gộp toàn bộ file FASTA
rm -f "${COMBINED_FA}" "${COMBINED_GZ}" "${COMBINED_GZ}.fai" "${COMBINED_GZ}.gzi"
sample_count=0
total_samples=$(ls -1 "${RAW_DIR}"/*.fa 2>/dev/null | wc -l)

for fa_file in "${RAW_DIR}"/*.fa; do
    if [[ -f "${fa_file}" ]]; then
        sample_name=$(basename "${fa_file}" .fa)
        sample_count=$((sample_count + 1))
        echo -ne "  -> Tiến độ Giai đoạn 1 (PanSN): Đã xử lý ${sample_count}/${total_samples} mẫu (${sample_name})   \r"
        sed -e "s/^>/>${sample_name}#1#/g" -e 's/ .*//g' "${fa_file}" >> "${COMBINED_FA}"
    fi
done
echo ""
bgzip -f -@ 4 "${COMBINED_FA}"
samtools faidx "${COMBINED_GZ}"

# 5. GIAI ĐOẠN 2: Dựng Đồ thị Mạng lưới & Phân cụm Cộng đồng Leiden (Mapping Graph & Plotting)
wfmash -m "${COMBINED_GZ}" -s 1000 -p 90 -t 8 > "${PAF_FILE}" 2>/dev/null
paf2net.py -p "${PAF_FILE}" > /dev/null

comm_output=$(net2communities.py \
    -e "${PAF_FILE}.edges.list.txt" \
    -w "${PAF_FILE}.edges.weights.txt" \
    -n "${PAF_FILE}.vertices.id2name.txt" \
    --output-prefix "${COMMUNITY_PREFIX}" \
    --accurate-detection \
    --plot 2>&1)

num_communities=$(echo "${comm_output}" | grep -oP 'Detected \K[0-9]+' || echo "N/A")

if command -v pdftoppm &> /dev/null; then
    pdftoppm -png -r 150 "${COMMUNITY_PREFIX}.communities.pdf" "${INTERN_DIR}/hla_all_community_temp" > /dev/null 2>&1
    mv "${INTERN_DIR}/hla_all_community_temp-1.png" "${GRAPH_PLOT_PNG}"
fi

# 6. GIAI ĐOẠN 3: Phân đoạn tệp FASTA theo từng Cộng đồng / Locus
mapfile -t LOCI_LIST < <(cut -f 1 -d '#' "${COMBINED_GZ}.fai" | sort | uniq)
total_loci=${#LOCI_LIST[@]}
locus_count=0

for gene in "${LOCI_LIST[@]}"; do
    part_gz="${PART_DIR}/hla_${gene}.fasta.gz"
    mapfile -t seq_ids < <(grep "^${gene}#" "${COMBINED_GZ}.fai" | cut -f 1)
    locus_count=$((locus_count + 1))
    echo -ne "  -> Tiến độ Giai đoạn 3 (Partitioning): Đã bóc tách ${locus_count}/${total_loci} Locus (${gene})   \r"
    samtools faidx "${COMBINED_GZ}" "${seq_ids[@]}" | bgzip -f -@ 4 > "${part_gz}"
    samtools faidx "${part_gz}"
done
echo ""

# 7. GIAI ĐOẠN 4 & TẠO TỆP CẤU HÌNH JSON
python3 -c "
import glob, os, math, json

loci_list = [l for l in '''${LOCI_LIST[*]}'''.split() if l]
intern_dir = '${INTERN_DIR}'
part_dir = '${PART_DIR}'
dist_dir = '${DIST_DIR}'
json_out = '${JSON_PARAMS_FILE}'

params_dict = {}
p_list = []
s_list = []

for gene in loci_list:
    part_gz = os.path.join(part_dir, f'hla_{gene}.fasta.gz')
    part_fai = part_gz + '.fai'
    dist_txt = os.path.join(dist_dir, f'hla_{gene}.mash_triangle.txt')
    
    lengths = [int(line.split()[1]) for line in open(part_fai)]
    seq_count = len(lengths)
    mean_l = sum(lengths) / seq_count
    min_l = min(lengths)
    
    calc_s = int(round(mean_l * 0.25 / 100) * 100)
    opt_s = max(200, min(calc_s, min_l))
    
    max_d = 0.0
    if seq_count >= 2:
        os.system(f'mash triangle \"{part_gz}\" > \"{dist_txt}\" 2>/dev/null')
        if os.path.exists(dist_txt):
            with open(dist_txt) as f:
                for line in f:
                    for token in line.strip().split():
                        try:
                            v = float(token)
                            if 0.0 < v < 0.9 and v > max_d:
                                max_d = v
                        except ValueError:
                            pass

    exact_limit = 100.0 - (max_d * 100.0)
    rec_p = max(70, math.floor(exact_limit - 3.0)) if max_d > 0 else 95
    
    p_list.append(rec_p)
    s_list.append(opt_s)
    
    params_dict[gene] = {
        'locus': gene,
        'fasta_path': part_gz,
        'num_samples': seq_count,
        'mean_length_bp': int(mean_l),
        'max_divergence': round(max_d, 6),
        'exact_identity_limit': round(exact_limit, 2),
        'p_param': rec_p,
        's_param': opt_s,
        'vcf_ref': 'grch38'
    }

with open(json_out, 'w') as f:
    json.dump(params_dict, f, indent=2)

avg_p = sum(p_list) / len(p_list)
avg_s = sum(s_list) / len(s_list)

print(f'AVG_P={avg_p:.1f}')
print(f'AVG_S={avg_s:.0f}')
" > "${INTERN_DIR}/avg_stats.tmp"

avg_p=$(grep 'AVG_P=' "${INTERN_DIR}/avg_stats.tmp" | cut -d '=' -f 2)
avg_s=$(grep 'AVG_S=' "${INTERN_DIR}/avg_stats.tmp" | cut -d '=' -f 2)
rm -f "${INTERN_DIR}/avg_stats.tmp"

echo ""
echo "================================================================="
echo "[SUMMARY] TỔNG KẾT THÔNG TIN TIỀN XỬ LÝ DỮ LIỆU HLA:"
echo "  • Số lượng mẫu ban đầu (RAW):      ${sample_count} mẫu"
echo "  • Số lượng khu vực cộng đồng (GĐ2): ${num_communities} khu vực (Leiden Communities)"
echo "  • Số lượng file đầu ra (GĐ3):      ${total_loci} tệp phân đoạn Locus FASTA"
echo "  • Giá trị tham số trung bình:       -p ${avg_p}% | -s ${avg_s} bp"
echo "  • File cấu hình tham số JSON xuất: [data/intern/hla_loci_params.json](file://${JSON_PARAMS_FILE})"
echo "================================================================="
