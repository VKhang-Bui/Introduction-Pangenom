#!/bin/bash

# Script chuyển đổi header FASTA sang chuẩn PanSN bằng Bash thuần (dùng sed)
# Cách dùng: ./bash_native.sh <thu_muc_input> <thu_muc_output>

input=$1
output=$2

mkdir -p "$output"

for f in "$input"/*.fa; do
    [ -e "$f" ] || continue
    name=$(basename "$f" .fa)
    sed "s/>/>${name}#1#/" "$f" > "${output}/${name}_pansn.fa"
done
