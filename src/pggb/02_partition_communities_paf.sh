#!/usr/bin/env bash
# ==============================================================================
# 02_partition_communities_paf.sh: PAF-based Community Detection (wfmash + paf2net)
#
# Usage: ./02_partition_communities_paf.sh [INPUT_FASTA] [OUT_DIR] [THREADS]
# Example: ./02_partition_communities_paf.sh data/intern/hla_all.fasta.gz data/intern/partitions_paf 4
# ==============================================================================

set -euo pipefail

INPUT_FASTA="${1:-data/intern/hla_all.fasta.gz}"
OUT_DIR="${2:-data/intern/partitions_paf}"
THREADS="${3:-4}"

TMP_DIR="${OUT_DIR}/tmp_net"
mkdir -p "${OUT_DIR}" "${TMP_DIR}"

if [[ ! -f "${INPUT_FASTA}.fai" ]]; then
    samtools faidx "${INPUT_FASTA}"
fi

# 1. Tự động tính NUM_MAPPINGS (-n) = Số lượng mẫu - 1
NUM_SAMPLES=$(cut -f 1 "${INPUT_FASTA}.fai" | cut -f 1 -d '#' | sort | uniq | wc -l)
NUM_MAPPINGS=$(( NUM_SAMPLES > 1 ? NUM_SAMPLES - 1 : 1 ))

# 2. Tự động tính MAP_PCT_ID (-p) từ mash triangle: p = 100 - (max_d * 100) - 3
DIST_TXT="${TMP_DIR}/mash_dist.txt"
mash triangle -p "${THREADS}" "${INPUT_FASTA}" > "${DIST_TXT}" 2>/dev/null || true

MAX_D=$(sed '1d' "${DIST_TXT}" | tr '\t' '\n' | grep -E '^[0-9]' | sort -gr | awk 'NR==1')
MAP_PCT_ID=$(awk -v d="${MAX_D:-0}" 'BEGIN { p = 100 - (d * 100) - 3; print (p < 90 ? 90 : int(p)) }')

echo "[AUTO] Computed Parameters: -p ${MAP_PCT_ID}% | -n ${NUM_MAPPINGS} mappings | Unique Samples: ${NUM_SAMPLES}"

PAF_FILE="${TMP_DIR}/mapping.paf"
echo "[1/4] Running pairwise mapping with wfmash..."
wfmash "${INPUT_FASTA}" -p "${MAP_PCT_ID}" -n "${NUM_MAPPINGS}" -t "${THREADS}" -m > "${PAF_FILE}"

echo "[2/4] Projecting PAF mappings to network with paf2net.py..."
paf2net.py -p "${PAF_FILE}"

echo "[3/4] Detecting communities with net2communities.py..."
net2communities.py \
    -e "${PAF_FILE}.edges.list.txt" \
    -w "${PAF_FILE}.edges.weights.txt" \
    -n "${PAF_FILE}.vertices.id2name.txt" \
    --plot

mv "${PAF_FILE}.edges.weights.txt.communities.pdf" "${OUT_DIR}/communities_paf_network.pdf"

echo "[4/4] Extracting partitioned FASTA files..."
for comm_file in "${PAF_FILE}.edges.weights.txt".community.*.txt; do
    comm_id=$(basename "$comm_file" | sed -E 's/.*\.community\.([0-9]+)\.txt$/\1/')
    OUT_FASTA="${OUT_DIR}/community_${comm_id}.fa.gz"

    mapfile -t SEQ_IDS < "$comm_file"
    samtools faidx "${INPUT_FASTA}" "${SEQ_IDS[@]}" | bgzip -f -@ "${THREADS}" -c > "${OUT_FASTA}"
    samtools faidx "${OUT_FASTA}"
done

echo "[CLEANUP] Removing temporary directory ${TMP_DIR}..."
rm -rf "${TMP_DIR}"

echo "Done! PAF Partitioned FASTAs and network PDF saved in: ${OUT_DIR}"
