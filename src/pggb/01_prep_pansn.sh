#!/usr/bin/env bash
# ==============================================================================
# 01_prep_pansn.sh: Standardize FASTA headers to true PanSN spec (sample#1#contig)
# Correctly extracts sample names (GRCh38, APD, COX, DBB, MANN, MCF, QBL, SSTO, CHM1, HuRef)
# ==============================================================================

set -euo pipefail

RAW_DIR="${1:-data/raw/HLA}"
PREFIX="$(basename "$RAW_DIR" | tr '[:upper:]' '[:lower:]')"
OUT_DIR="${3:-data/intern}"

COMBINED_FA="${OUT_DIR}/${PREFIX}_all.fasta"
COMBINED_GZ="${OUT_DIR}/${PREFIX}_all.fasta.gz"

mkdir -p "${OUT_DIR}"
rm -f "${COMBINED_FA}" "${COMBINED_GZ}" "${COMBINED_GZ}.fai"

echo "[1/2] Parsing raw FASTAs & extracting true PanSN sample names..."

python3 -c "
import os, sys, glob, gzip, re

raw_dir = '${RAW_DIR}'
out_fa_path = '${COMBINED_FA}'

with open(out_fa_path, 'w') as out_f:
    for f in sorted(glob.glob(os.path.join(raw_dir, '*'))):
        if not (f.endswith('.fa') or f.endswith('.fasta') or f.endswith('.fa.gz') or f.endswith('.fasta.gz')):
            continue
        
        gene_name = os.path.basename(f).split('/')[-1].split('-')[0]
        opener = gzip.open if f.endswith('.gz') else open
        
        with opener(f, 'rt') as in_f:
            for line in in_f:
                if line.startswith('>'):
                    raw_h = line[1:].strip()
                    
                    # 1. Check if already in PanSN spec (3+ parts separated by #)
                    if raw_h.count('#') >= 2:
                        out_f.write(f'>{raw_h}\n')
                        continue
                    
                    # 2. Extract true Sample ID from NCBI/GenBank header
                    sample = 'GRCh38'
                    if 'MHC_' in raw_h:
                        m = re.search(r'MHC_([A-Z0-9]+)_', raw_h)
                        if m:
                            sample = m.group(1)
                    elif 'CHM1' in raw_h:
                        sample = 'CHM1'
                    elif 'HuRef' in raw_h:
                        sample = 'HuRef'
                    elif 'GRCh38' in raw_h or 'Primary Assembly' in raw_h:
                        sample = 'GRCh38'
                    else:
                        sample = 'Sample_' + gene_name
                    
                    # Extract contig accession
                    accession = raw_h.split()[0].replace('|', '_').replace(':', '_')
                    clean_header = f'{sample}#1#{gene_name}_{accession}'
                    out_f.write(f'>{clean_header}\n')
                else:
                    out_f.write(line)
"

echo "[2/2] Compressing with bgzip & indexing with samtools..."
bgzip -f -@ 4 "${COMBINED_FA}"
samtools faidx "${COMBINED_GZ}"

echo "Done! Created: ${COMBINED_GZ} (Index: ${COMBINED_GZ}.fai)"
