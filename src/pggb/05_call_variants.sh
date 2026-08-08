#!/usr/bin/env bash
# ==============================================================================
# 05_call_variants.sh: Gọi biến thể nhỏ cho 1 Locus / Community chỉ định
#
# Options:
#   -i <COMMUNITY_ID>   Tên community hoặc thư mục (Mặc định: community_2)
#   -o <OUT_DIR>        Thư mục đầu ra VCF (Mặc định: data/intern/vcf_calls/<COMMUNITY_ID>)
#   -r <REF_PREFIX>     Tên mẫu tham chiếu (Mặc định: GRCh38)
#   -t <THREADS>        Số luồng CPU (Mặc định: 4)
#   -h                  Hiển thị hướng dẫn
#
# Examples:
#   ./05_call_variants.sh                        # Chạy mặc định cho community_2
#   ./05_call_variants.sh -i community_9         # Chạy cho community_9
# ==============================================================================

set -euo pipefail

# Tự động nhận diện đường dẫn Conda Environment (Git-friendly)
if [[ -n "${CONDA_PREFIX:-}" ]]; then
    export PATH="${CONDA_PREFIX}/bin/scripts:${CONDA_PREFIX}/bin:${PATH}"
fi

TARGET_COMM="community_2"
REF_PREFIX="GRCh38"
THREADS="4"
CUSTOM_OUT_DIR=""

while getopts "i:o:r:t:h" opt; do
  case $opt in
    i) TARGET_COMM="$OPTARG" ;;
    o) CUSTOM_OUT_DIR="$OPTARG" ;;
    r) REF_PREFIX="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    h)
      echo "Usage: $0 [-i community_X] [-o out_dir] [-r ref_prefix] [-t threads]"
      exit 0
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

COMM_NAME=$(basename "${TARGET_COMM}" .fa.gz)

GRAPH_DIR="data/intern/graphs_partitions_mash/${COMM_NAME}"
COMM_FASTA="data/intern/partitions_mash/${COMM_NAME}.fa.gz"

if [[ -n "$CUSTOM_OUT_DIR" ]]; then
    OUT_DIR="$CUSTOM_OUT_DIR"
else
    OUT_DIR="data/intern/vcf_calls/${COMM_NAME}"
fi

mkdir -p "${OUT_DIR}"

OG_FILE=$(ls "${GRAPH_DIR}"/*.final.og 2>/dev/null || ls "${GRAPH_DIR}"/*.og 2>/dev/null | head -n 1)

if [[ ! -f "$OG_FILE" ]]; then
    echo "Lỗi: Không tìm thấy file đồ thị .og tại '${GRAPH_DIR}'."
    exit 1
fi

REF_PATH=$(odgi paths -i "$OG_FILE" -L | grep "^${REF_PREFIX}" | head -n 1)
if [[ -z "$REF_PATH" ]]; then
    REF_PATH=$(odgi paths -i "$OG_FILE" -L | head -n 1)
fi


PG_VCF="${OUT_DIR}/pg_variants.vcf.gz"

odgi view -i "${OG_FILE}" -g \
    | vg deconstruct -P "${REF_PREFIX}" -a - \
    | vcfwave \
    | bcftools sort \
    | bgzip -f -@ "${THREADS}" -c > "${PG_VCF}"

tabix -f -p vcf "${PG_VCF}"
``
REF_FASTA="${OUT_DIR}/ref_${REF_PREFIX}.fa"
samtools faidx "${COMM_FASTA}" "${REF_PATH}" > "${REF_FASTA}" 2>/dev/null || samtools faidx "${COMM_FASTA}" > "${REF_FASTA}"
samtools faidx "${REF_FASTA}"

QUERY_FASTA="${OUT_DIR}/query_${COMM_NAME}.fa"
if [[ "$COMM_FASTA" == *.gz ]]; then
    zcat "${COMM_FASTA}" > "${QUERY_FASTA}"
else
    cp "${COMM_FASTA}" "${QUERY_FASTA}"
fi
samtools faidx "${QUERY_FASTA}"

MUMMER_VCF="${OUT_DIR}/mummer_baseline.vcf.gz"
MUMMER_PREFIX="${OUT_DIR}/nucmer_align"

nucmer --maxmatch "${REF_FASTA}" "${QUERY_FASTA}" -p "${MUMMER_PREFIX}"
show-snps -H -T -r -x 1 "${MUMMER_PREFIX}.delta" | cut -f 1-6,11-14 > "${MUMMER_PREFIX}.snps"

NUCMER2VCF_SCRIPT=$(which nucmer2vcf.R)
TMP_MUMMER_VCF="${MUMMER_PREFIX}.tmp.vcf"
Rscript "${NUCMER2VCF_SCRIPT}" "${MUMMER_PREFIX}.snps" "${REF_PREFIX}" "${REF_FASTA}" "MUMMER4" "${TMP_MUMMER_VCF}"

bgzip -f -@ "${THREADS}" -c "${TMP_MUMMER_VCF}" > "${MUMMER_VCF}"
tabix -f -p vcf "${MUMMER_VCF}"

rm -f "${QUERY_FASTA}" "${QUERY_FASTA}.fai" "${MUMMER_PREFIX}.delta" "${MUMMER_PREFIX}.snps" "${MUMMER_PREFIX}.ntref" "${TMP_MUMMER_VCF}"
