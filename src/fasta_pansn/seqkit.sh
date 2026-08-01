#!/bin/bash

# Script chuyển đổi header FASTA sang chuẩn PanSN bằng SeqKit (Tương đương logic bash_native.sh)
# Cách dùng: ./seqkit.sh [thu_muc_input] [thu_muc_output]

input="${1:-/home/vkhang-bui/1.HocViec/projects/pangenom/data/raw}"
output="${2:-/home/vkhang-bui/1.HocViec/projects/pangenom/data}"

mkdir -p "$output"
outfile="$output/pansn.fasta"
> "$outfile"  # Tạo/xóa trắng file output rỗng trước khi gộp

for f in "$input"/*.fa; do
    [ -e "$f" ] || continue
    name=$(basename "$f" .fa)

    # Dùng seqkit replace để thay thế header:
    # -n : Áp dụng trên toàn bộ header (Full header name)
    # -p "^" : Khớp vị trí đầu tiên của header
    # -r "${name}#1#read_{nr}_" : Chèn prefix PanSN và chỉ số thứ tự {nr} tự động của SeqKit
    seqkit replace -p "^" -r "${name}#1#read_{nr}_" "$f" >> "$outfile"

    count=$(seqkit stats -t dna "$f" | awk 'NR==2{print $4}')
    echo "✓ Đã xử lý file '$name.fa' với $count reads/contigs."
done