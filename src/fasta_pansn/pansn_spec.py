#!/usr/bin/env python3
"""
Script chuyển đổi header FASTA theo đúng chuẩn PanSN (NCBI/HPRC Spec) bằng Python.
Cách dùng: python3 pansn_spec.py <thu_muc_input> <thu_muc_output>
"""
import sys
import os
import glob
import re

def sanitize_pansn_name(name):
    # Quy chuẩn PanSN: Loại bỏ ký tự đặc biệt/khoảng trắng trong tên mẫu để không làm hỏng cú pháp đồ thị GFA
    return re.sub(r'[^a-zA-Z0-9_.]', '_', name)

def process_pansn_spec(input_dir, output_dir, haplotype="1"):
    os.makedirs(output_dir, exist_ok=True)
    fasta_files = glob.glob(os.path.join(input_dir, "*.fa")) + glob.glob(os.path.join(input_dir, "*.fasta"))
    
    if not fasta_files:
        print(f"Không tìm thấy file .fa hoặc .fasta nào trong: {input_dir}")
        return

    for fasta_file in fasta_files:
        raw_name = os.path.splitext(os.path.basename(fasta_file))[0]
        sample_name = sanitize_pansn_name(raw_name)
        output_file = os.path.join(output_dir, f"{sample_name}_pansn.fa")
        
        with open(fasta_file, 'r') as infile, open(output_file, 'w') as outfile:
            for line in infile:
                if line.startswith('>'):
                    # Lấy contig ID chuẩn, lọc bỏ mô tả thừa phía sau nếu có khoảng trắng
                    contig_id = line.strip().lstrip('>').split()[0]
                    outfile.write(f">{sample_name}#{haplotype}#{contig_id}\n")
                else:
                    outfile.write(line)
        print(f"Đã chuyển đổi chuẩn PanSN Spec: {os.path.basename(fasta_file)} -> {os.path.basename(output_file)}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Cú pháp: python3 pansn_spec.py <thu_muc_input> <thu_muc_output>")
        sys.exit(1)
    
    input_dir = sys.argv[1]
    output_dir = sys.argv[2]
    process_pansn_spec(input_dir, output_dir)
