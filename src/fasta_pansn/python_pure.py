#!/usr/bin/env python3
import glob, os

input_dir = "/home/vkhang-bui/1.HocViec/projects/pangenom/data/raw"
output_dir = "/home/vkhang-bui/1.HocViec/projects/pangenom/data"

os.makedirs(output_dir, exist_ok=True)
outfile_path = os.path.join(output_dir, "pansn.fasta")
open(outfile_path, "w").close()  # Tạo file output rỗng

for f in sorted(glob.glob(f"{input_dir}/*.fa")):
    name = os.path.basename(f)[:-3]
    read_idx = 1
    
    with open(f, "r") as infile, open(outfile_path, "a") as outfile:
        for line in infile:
            if line.startswith(">"):
                outfile.write(f">{name}#1#read_{read_idx}_{line[1:]}")
                read_idx += 1
            else:
                outfile.write(line)

    print(f"✓ Đã xử lý file '{name}.fa' với {read_idx - 1} reads/contigs.")
