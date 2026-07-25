#!/bin/bash

# Script xây dựng đồ thị Pangenome bằng PGGB và tự động xuất ảnh trực quan 1D/2D
# Cách dùng: ./run.sh <file_input.fa> <thu_muc_output> <so_haplotype>

input=$1
output=$2
n_haplo=$3

mkdir -p "$output"

export PATH=/home/vkhang-bui/miniforge3/envs/shina/bin:$PATH

pggb \
  -i "$input" \
  -o "$output" \
  -n "$n_haplo" \
  -t 8 \
  -T 1 \
  -p 70 \
  -s 3000 \
  -V "grch38#1#chr6:" || true

# Tự động xuất ảnh trực quan 1D (odgi viz) và 2D (Bandage)
gfa_file=$(ls -1 "$output"/*.smooth.gfa 2>/dev/null | head -n 1)

if [ -f "$gfa_file" ]; then
    echo "Đang tạo ảnh trực quan 1D và 2D..."
    odgi build -g "$gfa_file" -o "${output}/pangenome.og"
    odgi viz -i "${output}/pangenome.og" -o "${output}/graph_1d_viz.png" -x 1500 -y 500 -P
    Bandage image "$gfa_file" "${output}/graph_2d_bandage.png"
    echo "Đã tạo xong ảnh 1D và 2D trong thư mục $output"
fi
