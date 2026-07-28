#!/bin/bash
# Script xuat 2 so do 1D ODGI quan trong nhat (-M va -d)

INPUT=$1
OUT_DIR=${2:-./results}
mkdir -p $OUT_DIR

# 1. Chuyen GFA sang OG
odgi build -g $INPUT -o $OUT_DIR/graph.og

# 2. Xuat so do Haplotype Core/Accessory (-M)
odgi viz -i $OUT_DIR/graph.og -o $OUT_DIR/viz_multiqc.png -M

# 3. Xuat so do Do sau Coverage & CNVs (-d)
odgi viz -i $OUT_DIR/graph.og -o $OUT_DIR/viz_depth.png -d

echo "Da tao xong: $OUT_DIR/viz_multiqc.png va $OUT_DIR/viz_depth.png"
