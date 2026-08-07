#!/bin/bash

# Script xây dựng đồ thị Pangenome bằng PGGB và tự động xuất ảnh trực quan 1D/2D
# Cách dùng: ./run.sh <file_input.fa> <thu_muc_output> <so_haplotype>

input=$1
output=$2
n_haplo=$3

mkdir -p "$output"

export PATH=/home/vkhang-bui/miniforge3/envs/pggb_env/bin:$PATH

pggb \
  -i "$input" \
  -o "$output" \
  -n "$n_haplo" \
  -t 8 \
  # -T 1 \
  # -p 70 \
  # -s 3000 \
  -v
