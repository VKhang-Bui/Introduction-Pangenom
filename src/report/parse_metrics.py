#!/usr/bin/env python3
import json
import os
import glob
import gzip
import math

REPORT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.abspath(os.path.join(REPORT_DIR, "../.."))
OUT_IMG_DIR = os.path.join(BASE_DIR, "output/image")
INTERN_DIR = os.path.join(BASE_DIR, "data/intern")

summary_results = {}

# 1. Phân tích trực tiếp GFA file (DRB1-3123.fa.gz.e8de02e.11fba48.89cf512.smooth.final.gfa.gz)
gfa_gz = os.path.join(OUT_IMG_DIR, "DRB1-3123.fa.gz.e8de02e.11fba48.89cf512.smooth.final.gfa.gz")

nodes = {} # node_id -> length
paths = {} # path_name -> list of node_ids
path_nodes = {} # path_name -> set of node_ids

if os.path.exists(gfa_gz):
    with gzip.open(gfa_gz, "rt") as f:
        for line in f:
            if line.startswith("S\t"):
                parts = line.strip().split("\t")
                node_id = parts[1]
                seq = parts[2]
                nodes[node_id] = len(seq)
            elif line.startswith("P\t"):
                parts = line.strip().split("\t")
                path_name = parts[1]
                segment_str = parts[2]
                node_list = []
                n_set = set()
                for seg in segment_str.split(","):
                    if seg:
                        nid = seg[:-1]
                        node_list.append(nid)
                        n_set.add(nid)
                paths[path_name] = node_list
                path_nodes[path_name] = n_set

    # Thống kê phân bố Node Coverage (Depth)
    node_hap_count = {nid: 0 for nid in nodes}
    for pname, nset in path_nodes.items():
        for nid in nset:
            node_hap_count[nid] += 1

    depth_dist = {}
    core_nodes = 0
    shell_nodes = 0
    unique_nodes = 0
    core_bp = 0
    shell_bp = 0
    unique_bp = 0

    total_paths = len(paths)
    for nid, count in node_hap_count.items():
        l = nodes[nid]
        depth_dist[count] = depth_dist.get(count, 0) + 1
        if count == total_paths:
            core_nodes += 1
            core_bp += l
        elif count == 1:
            unique_nodes += 1
            unique_bp += l
        else:
            shell_nodes += 1
            shell_bp += l

    total_bp = sum(nodes.values())

    # Tính ma trận Jaccard Similarity giữa các Path dựa trên tập node chung
    path_names = sorted(list(paths.keys()))
    jaccard_scores = []
    sim_matrix = {}

    for i in range(len(path_names)):
        p1 = path_names[i]
        sim_matrix[p1] = {}
        for j in range(len(path_names)):
            p2 = path_names[j]
            s1 = path_nodes[p1]
            s2 = path_nodes[p2]
            intersection = len(s1.intersection(s2))
            union = len(s1.union(s2))
            jaccard = intersection / union if union > 0 else 1.0
            sim_matrix[p1][p2] = round(jaccard, 4)
            if i < j:
                jaccard_scores.append(jaccard)

    mean_j = sum(jaccard_scores) / len(jaccard_scores) if jaccard_scores else 0.0
    min_j = min(jaccard_scores) if jaccard_scores else 0.0
    max_j = max(jaccard_scores) if jaccard_scores else 0.0

    summary_results['graph_gfa_topology'] = {
        'total_nodes': len(nodes),
        'total_graph_bp': total_bp,
        'total_paths': len(paths),
        'core_nodes_count': core_nodes,
        'core_nodes_bp': core_bp,
        'core_bp_pct': round(core_bp / total_bp * 100, 2) if total_bp > 0 else 0,
        'shell_nodes_count': shell_nodes,
        'shell_nodes_bp': shell_bp,
        'shell_bp_pct': round(shell_bp / total_bp * 100, 2) if total_bp > 0 else 0,
        'unique_nodes_count': unique_nodes,
        'unique_nodes_bp': unique_bp,
        'unique_bp_pct': round(unique_bp / total_bp * 100, 2) if total_bp > 0 else 0,
        'depth_distribution': depth_dist,
        'mean_jaccard_similarity': round(mean_j, 4),
        'min_jaccard_similarity': round(min_j, 4),
        'max_jaccard_similarity': round(max_j, 4),
    }

# 2. Phân tích thông tin tập tin hình ảnh (.png)
image_files = sorted(glob.glob(os.path.join(OUT_IMG_DIR, "*.png")) + glob.glob(os.path.join(INTERN_DIR, "*.png")))
img_metadata = {}
for img_p in image_files:
    rel_p = os.path.relpath(img_p, BASE_DIR)
    img_metadata[rel_p] = {
        'size_bytes': os.path.getsize(img_p),
        'size_kb': round(os.path.getsize(img_p) / 1024, 1)
    }

summary_results['image_metadata'] = img_metadata

# 3. Reading MultiQC Data JSON
mqc_json = os.path.join(OUT_IMG_DIR, "multiqc_data/multiqc_data.json")
if os.path.exists(mqc_json):
    with open(mqc_json) as f:
        mqc_data = json.load(f)
    summary_results['multiqc_general_stats'] = mqc_data.get('report_general_stats_data', [])

# 4. Ghi kết quả tổng hợp ra JSON và TXT trong src/report/
json_out = os.path.join(REPORT_DIR, "drb1_deep_metrics_summary.json")
with open(json_out, "w") as f:
    json.dump(summary_results, f, indent=2)

txt_out = os.path.join(REPORT_DIR, "drb1_deep_metrics_summary.txt")
with open(txt_out, "w", encoding="utf-8") as f:
    f.write("=================================================================\n")
    f.write("BÁO CÁO PHÂN TÍCH ĐỔI MỚI SỐ LIỆU & HÌNH ẢNH ĐỒ THỊ PANGENOME HLA-DRB1\n")
    f.write("=================================================================\n\n")
    
    if 'graph_gfa_topology' in summary_results:
        gt = summary_results['graph_gfa_topology']
        f.write("1. PHÂN TÍCH CHUYÊN SÂU TOPOLOGY & ĐỘ BAO PHỦ NUCLEOTIDE (GFA):\n")
        f.write(f"   • Tổng số Nodes:          {gt['total_nodes']}\n")
        f.write(f"   • Tổng dung lượng Đồ thị: {gt['total_graph_bp']} bp\n")
        f.write(f"   • Số lượng Haplotypes:   {gt['total_paths']}\n")
        f.write(f"   • Vùng Core (12/12):      {gt['core_nodes_count']} nodes | {gt['core_nodes_bp']} bp ({gt['core_bp_pct']}%)\n")
        f.write(f"   • Vùng Shell (2..11):     {gt['shell_nodes_count']} nodes | {gt['shell_nodes_bp']} bp ({gt['shell_bp_pct']}%)\n")
        f.write(f"   • Vùng Unique (1/12):     {gt['unique_nodes_count']} nodes | {gt['unique_nodes_bp']} bp ({gt['unique_bp_pct']}%)\n")
        f.write(f"   • Chỉ số Tương đồng Jaccard (Path Similarity):\n")
        f.write(f"     - Trực quan Trung bình:  {gt['mean_jaccard_similarity']}\n")
        f.write(f"     - Thấp nhất (Phân dị nhất): {gt['min_jaccard_similarity']}\n")
        f.write(f"     - Cao nhất (Tương đồng nhất): {gt['max_jaccard_similarity']}\n")
        f.write("   • Phân bố Độ bao phủ (Depth Distribution):\n")
        for depth, count in sorted(gt['depth_distribution'].items(), key=lambda x: int(x[0])):
            f.write(f"     - Haplotype Depth {depth}: {count} nodes\n")
    
    f.write("\n2. PHÂN TÍCH SIÊU DỮ LIỆU CÁC HÌNH ẢNH SƠ ĐỒ SẢN XUẤT (.PNG):\n")
    for img_name, meta in img_metadata.items():
        f.write(f"   • [{img_name}]: {meta.get('size_kb')} KB\n")

print(f"[SUCCESS] Deep analysis completed. Results written to:\n  - {json_out}\n  - {txt_out}")
