#!/usr/bin/env python3
"""
Script chuyển đổi header FASTA sang chuẩn PanSN bằng Python thuần.
Cách dùng: python3 python_pure.py <thu_muc_input> <thu_muc_output>
"""
import sys
import os
import glob

def convert_fasta_to_pansn(input_dir, output_dir, haplotype="1"):
    os.makedirs(output_dir, exist_ok=True)
    fasta_files = glob.glob(os.path.join(input_dir, "*.fa")) + glob.glob(os.path.join(input_dir, "*.fasta"))
    
    if not fasta_files:
        print(f"Không tìm thấy file .fa hoặc .fasta nào trong: {input_dir}")
        return

    for fasta_file in fasta_files:
        sample_name = os.path.splitext(os.path.basename(fasta_file))[0]
        output_file = os.path.join(output_dir, f"{sample_name}_pansn.fa")
        
        with open(fasta_file, 'r') as infile, open(output_file, 'w') as outfile:
            for line in infile:
                if line.startswith('>'):
                    contig_name = line.strip().lstrip('>')
                    outfile.write(f">{sample_name}#{haplotype}#{contig_name}\n")
                else:
                    outfile.write(line)
        print(f"Đã xử lý xong: {os.path.basename(fasta_file)} -> {os.path.basename(output_file)}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Cú pháp: python3 python_pure.py <thu_muc_input> <thu_muc_output>")
        sys.exit(1)
    
    input_dir = sys.argv[1]
    output_dir = sys.argv[2]
    convert_fasta_to_pansn(input_dir, output_dir)
