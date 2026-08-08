#!/usr/bin/env bash
# ==============================================================================
# 03_estimate_divergence.sh: Estimate sequence divergence for each community
# Displays both Total_Seqs (Tổng số chuỗi) and Num_Samples (Số mẫu độc lập).
#
# Usage: ./03_estimate_divergence.sh [PARTITIONS_DIR] [OUT_TXT]
# Example: ./03_estimate_divergence.sh data/intern/partitions_mash
# -> Outputs: data/intern/divergence_summary_partitions_mash.txt
# ==============================================================================

set -euo pipefail

PARTITIONS_DIR="${1:-data/intern/partitions_mash}"

DIR_NAME=$(basename "${PARTITIONS_DIR}")
DEFAULT_OUT="data/intern/divergence_summary_${DIR_NAME}.txt"
OUT_TXT="${2:-${DEFAULT_OUT}}"

mkdir -p "$(dirname "${OUT_TXT}")"

echo "[1/1] Estimating sequence divergence for communities in ${PARTITIONS_DIR}..."

echo -e "Community\tTotal_Seqs\tNum_Samples\tMax_Divergence\tRecommended_P\tRecommended_S\tRecommended_N" > "${OUT_TXT}"

for chr_fasta in "${PARTITIONS_DIR}"/community_*.fa.gz; do
    if [[ -f "$chr_fasta" ]]; then
        comm_name=$(basename "$chr_fasta" | sed -E 's/\.fa\.gz$//')
        
        if [[ ! -f "${chr_fasta}.fai" ]]; then
            samtools faidx "$chr_fasta"
        fi

        # 1. Tổng số chuỗi contigs trong cụm
        TOTAL_SEQS=$(wc -l < "${chr_fasta}.fai")

        # 2. Số lượng mẫu độc lập (Unique Sample IDs)
        NUM_SAMPLES=$(cut -f 1 "${chr_fasta}.fai" | cut -f 1 -d '#' | sort | uniq | wc -l)
        REC_N=$(( NUM_SAMPLES > 1 ? NUM_SAMPLES - 1 : 1 ))

        # 3. Tính khoảng cách lớn nhất bằng mash triangle
        MAX_D=$(mash triangle "$chr_fasta" 2>/dev/null | sed '1d' | tr '\t' '\n' | grep -E '^[0-9]' | sort -gr | awk 'NR==1')
        MAX_D="${MAX_D:-0}"

        # 4. Tính -p và -s bằng awk đơn giản
        REC_P=$(awk -v d="$MAX_D" 'BEGIN { p = 100 - (d * 100) - 3; print (d > 0 ? (p < 70 ? 70 : int(p)) : 95) }')
        
        MEAN_LEN=$(awk '{sum += $2; count++} END {print (count > 0 ? int(sum/count) : 5000)}' "${chr_fasta}.fai")
        REC_S=$(awk -v l="$MEAN_LEN" 'BEGIN { s = int(l * 0.25); print (s < 200 ? 200 : s) }')

        echo -e "${comm_name}\t${TOTAL_SEQS}\t${NUM_SAMPLES}\t${MAX_D}\t${REC_P}\t${REC_S}\t${REC_N}" >> "${OUT_TXT}"
    fi
done

echo "Done! Divergence summary saved in: ${OUT_TXT}"
cat "${OUT_TXT}" | column -t
