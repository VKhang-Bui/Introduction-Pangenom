#!/usr/bin/env bash

# ==============================================================================
# pangenome_partitioning_pipeline.sh
# Sequence Partitioning & Divergence Estimation Pipeline
# Based directly on official PGGB Divergence Estimation Tutorial standard
# Usage: ./pangenome_partitioning_pipeline.sh [RAW_DIR] [PREFIX] [OUT_DIR]
# Example: ./pangenome_partitioning_pipeline.sh data/raw/HLA hla data/intern
# ==============================================================================

set -euo pipefail

RAW_DIR="${1:-data/raw/HLA}"
PREFIX="${2:-hla}"
OUT_DIR="${3:-data/intern}"

PART_DIR="${OUT_DIR}/partitions"
DIST_DIR="${OUT_DIR}/distances"
mkdir -p "${OUT_DIR}" "${PART_DIR}" "${DIST_DIR}"

COMBINED_GZ="${OUT_DIR}/${PREFIX}_all.fasta.gz"
JSON_PARAMS="${OUT_DIR}/${PREFIX}_loci_params.json"
DIVERGENCE_TXT="${OUT_DIR}/${PREFIX}_divergence.txt"
LOG_FILE="${OUT_DIR}/${PREFIX}_divergence.log"

echo "================================================================="
echo "[START] PANGENOME SEQUENCE PARTITIONING & DIVERGENCE ESTIMATION"
echo " RAW Dir: $RAW_DIR | Prefix: $PREFIX | Output: $OUT_DIR"
echo "================================================================="


echo "[1/4] Normalizing PanSN headers and combining FASTA files..."
COMBINED_FA="${OUT_DIR}/${PREFIX}_all.fasta"
rm -f "${COMBINED_FA}" "${COMBINED_GZ}" "${COMBINED_GZ}.fai"

for f in "${RAW_DIR}"/*.fa "${RAW_DIR}"/*.fasta; do
    if [[ -f "$f" ]]; then
        sample_name=$(basename "$f" | sed -E 's/\.(fa|fasta)$//')
        sed -e "s/^>/>${sample_name}#1#/g" -e 's/ .*//g' "$f" >> "${COMBINED_FA}"
    fi
done

bgzip -f -@ 4 "${COMBINED_FA}"
samtools faidx "${COMBINED_GZ}"


echo "[2/4] Partitioning sequences by contig/locus..."
mapfile -t CHROM_LIST < <(cut -f 1 "${COMBINED_GZ}.fai" | cut -f 1 -d '#' | sort | uniq)

for CHROM in "${CHROM_LIST[@]}"; do
    CHR_FASTA="${PART_DIR}/${PREFIX}_${CHROM}.fasta.gz"
    mapfile -t SEQ_IDS < <(grep -P $"^${CHROM}#" "${COMBINED_GZ}.fai" | cut -f 1)
    samtools faidx "${COMBINED_GZ}" "${SEQ_IDS[@]}" | bgzip -f -@ 4 > "${CHR_FASTA}"
    samtools faidx "${CHR_FASTA}"
done


echo "[3/4] Estimating sequence divergence with Mash..."
echo -e "Locus\tNum_Samples\tMean_Length_bp\tMax_Divergence\tExact_Identity_Limit\tRecommended_P\tRecommended_S" > "${DIVERGENCE_TXT}"

python3 -c "
import glob, os, math, json

chrom_list = '''${CHROM_LIST[*]}'''.split()
part_dir = '${PART_DIR}'
dist_dir = '${DIST_DIR}'
prefix = '${PREFIX}'
json_out = '${JSON_PARAMS}'
div_txt = '${DIVERGENCE_TXT}'

params_dict = {}

for chrom in chrom_list:
    chr_fasta = os.path.join(part_dir, f'{prefix}_{chrom}.fasta.gz')
    chr_fai = chr_fasta + '.fai'
    dist_txt = os.path.join(dist_dir, f'{prefix}_{chrom}.mash_triangle.txt')
    
    if not os.path.exists(chr_fai):
        continue

    lengths = [int(line.split()[1]) for line in open(chr_fai)]
    seq_count = len(lengths)
    mean_l = sum(lengths) / seq_count
    min_l = min(lengths)

    # Segment length s estimation: ~25% of mean sequence length
    calc_s = int(round(mean_l * 0.25 / 100) * 100)
    opt_s = max(200, min(calc_s, min_l))

    max_d = 0.0
    if seq_count >= 2:
        os.system(f'mash triangle \"{chr_fasta}\" > \"{dist_txt}\" 2>/dev/null')
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

    params_dict[chrom] = {
        'locus': chrom,
        'fasta_path': chr_fasta,
        'num_samples': seq_count,
        'mean_length_bp': int(mean_l),
        'max_divergence': round(max_d, 6),
        'exact_identity_limit': round(exact_limit, 2),
        'p_param': rec_p,
        's_param': opt_s,
        'vcf_ref': 'grch38'
    }

    with open(div_txt, 'a') as f:
        f.write(f'{chrom}\t{seq_count}\t{int(mean_l)}\t{max_d:.6f}\t{exact_limit:.2f}\t{rec_p}\t{opt_s}\n')

with open(json_out, 'w') as f:
    json.dump(params_dict, f, indent=2)

if prefix != 'hla_loci_params':
    hla_json = os.path.join('${OUT_DIR}', 'hla_loci_params.json')
    with open(hla_json, 'w') as f:
        json.dump(params_dict, f, indent=2)
"

LOG_FILE="${OUT_DIR}/${PREFIX}_divergence.log"
EXEC_TIME=$(date '+%Y-%m-%d %H:%M:%S')

cat << EOF > "${LOG_FILE}"
=================================================================
 PANGENOME SEQUENCE PARTITIONING & DIVERGENCE ESTIMATION LOG
=================================================================
Execution Timestamp:  ${EXEC_TIME}
Pipeline Script:      pangenome_partitioning_pipeline.sh
Raw Input Directory:  ${RAW_DIR}
Dataset Prefix:       ${PREFIX}
Output Directory:     ${OUT_DIR}
Total Loci Processed: ${#CHROM_LIST[@]}

[OUTPUT ASSET LOCATIONS]
  • Partition FASTAs:  ${PART_DIR}/
  • Mash Distances:    ${DIST_DIR}/
  • Divergence Table:  ${DIVERGENCE_TXT}
  • Parameter JSON:    ${JSON_PARAMS}

[SUMMARY METRICS TABLE]
$(cat "${DIVERGENCE_TXT}" | column -t)
=================================================================
EOF

echo ""
echo "================================================================="
echo "[SUMMARY] PANGENOME PARTITIONING & DIVERGENCE ESTIMATION COMPLETE"
echo "  • Total Loci Partitioned:  ${#CHROM_LIST[@]}"
echo "  • Divergence Summary TXT:  ${DIVERGENCE_TXT}"
echo "  • Execution Log Output:    ${LOG_FILE}"
echo "  • Parameter JSON Output:   ${JSON_PARAMS}"
echo "================================================================="
cat "${DIVERGENCE_TXT}" | column -t
