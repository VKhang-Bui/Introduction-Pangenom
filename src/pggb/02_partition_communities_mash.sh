#!/usr/bin/env bash
# ==============================================================================
# 02_partition_communities_mash.sh: Mash-based Community Detection (Khuyên dùng)
#
# Usage: ./02_partition_communities_mash.sh [INPUT_FASTA] [OUT_DIR] [THREADS] [SKETCH_SIZE]
# Example: ./02_partition_communities_mash.sh data/intern/hla_all.fasta.gz data/intern/partitions_mash 4 10000
# ==============================================================================

set -euo pipefail

INPUT_FASTA="${1:-data/intern/hla_all.fasta.gz}"
OUT_DIR="${2:-data/intern/partitions_mash}"
THREADS="${3:-4}"
SKETCH_SIZE="${4:-10000}"

TMP_DIR="${OUT_DIR}/tmp_net"
mkdir -p "${OUT_DIR}" "${TMP_DIR}"

if [[ ! -f "${INPUT_FASTA}.fai" ]]; then
    samtools faidx "${INPUT_FASTA}"
fi

DIST_TSV="${TMP_DIR}/distances.tsv"

echo "[1/4] Calculating pairwise Mash distances (-s ${SKETCH_SIZE})..."
mash dist "${INPUT_FASTA}" "${INPUT_FASTA}" -s "${SKETCH_SIZE}" -i -p "${THREADS}" > "${DIST_TSV}"

echo "[2/4] Projecting Mash distances to network with mash2net.py..."
mash2net.py -m "${DIST_TSV}"

echo "[3/4] Detecting communities with net2communities.py..."
net2communities.py \
    -e "${DIST_TSV}.edges.list.txt" \
    -w "${DIST_TSV}.edges.weights.txt" \
    -n "${DIST_TSV}.vertices.id2name.txt" \
    --plot

mv "${DIST_TSV}.edges.weights.txt.communities.pdf" "${OUT_DIR}/communities_mash_network.pdf"

echo "[4/4] Extracting partitioned FASTA files..."
for comm_file in "${DIST_TSV}.edges.weights.txt".community.*.txt; do
    comm_id=$(basename "$comm_file" | sed -E 's/.*\.community\.([0-9]+)\.txt$/\1/')
    OUT_FASTA="${OUT_DIR}/community_${comm_id}.fa.gz"

    mapfile -t SEQ_IDS < "$comm_file"
    samtools faidx "${INPUT_FASTA}" "${SEQ_IDS[@]}" | bgzip -f -@ "${THREADS}" -c > "${OUT_FASTA}"
    samtools faidx "${OUT_FASTA}"
done

echo "[CLEANUP] Removing temporary directory ${TMP_DIR}..."
rm -rf "${TMP_DIR}"

echo "Done! Mash Partitioned FASTAs and network PDF saved in: ${OUT_DIR}"
