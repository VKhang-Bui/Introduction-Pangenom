#!/bin/bash

# Script chuyển đổi header FASTA sang chuẩn PanSN bằng Bash thuần (dùng sed)
# Cách dùng: ./bash_native.sh <thu_muc_input> <thu_muc_output>

input="/home/vkhang-bui/1.HocViec/projects/pangenom/data/raw"
output="/home/vkhang-bui/1.HocViec/projects/pangenom/data"

mkdir -p "$output"
> "$output/pansn.fasta"  # Tạo file output rỗng
for f in "$input"/*.fa; do

    name=$(basename "$f" .fa)


    # Biến đếm thứ tự read/contig trong file (bắt đầu từ 1)
    read_idx=1
        
    # Vòng lặp 2: Đọc từng dòng trong tập tin FASTA
    while IFS= read -r line || [ -n "$line" ]; do
        # Kiểm tra nếu dòng hiện tại là Header (bắt đầu bằng ký tự '>')
        if [[ "$line" =~ ^\> ]]; then
            # Dùng sed thay thế '>' bằng '>Sample#Haplotype#read_N_' kết hợp tên gốc
            echo "$line" | sed "s/^>/>${name}#1#read_${read_idx}_/" >> "$output/pansn.fasta"   
            read_idx=$((read_idx + 1))
        else
            echo "$line" >> "$output/pansn.fasta"
        fi 
    done < "$f"
    echo "✓ Đã xử lý file '$name.fa' với $((read_idx - 1)) reads/contigs."
done
