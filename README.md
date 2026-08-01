# Introduction to Pangenomics & Practice (Introduction-Pangenom)

Chào mừng đến với kho lưu trữ nghiên cứu và thực hành về **Đồ thị Hệ gen toàn diện (Pangenome Graphs)**, tập trung vào mô hình sinh học **Nấm (Fungi)** và **Người (Humans)**. Dự án được xây dựng nhằm cung cấp lý thuyết nền tảng, kịch bản dòng lệnh phân tích và quy trình trực quan hóa biểu đồ phục vụ báo cáo khoa học.

---

## 👥 Đóng Góp Của Các Thành Viên (Team Contributions)

| Thành viên | Vai trò | Đóng góp chính |
| :--- | :--- | :--- |
| **VKhang-Bui** | **Project Lead (Chủ dự án)** | Quản lý dự án, định hướng mục tiêu nghiên cứu, duyệt tài liệu và điều phối công việc nhóm. |
| **Minh** | **Tooling & Scripting** | Tổng hợp các công cụ tin sinh học (`PGGB`, `Minigraph-Cactus`, `vg`, `odgi`), viết kịch bản Bash/Python tự động hóa và vận hành pipeline. |
| **Giàu** | **Graph Visualizations & Analysis** | Tổng hợp các mô hình/dạng đồ thị Pangenome, phụ trách phân tích và vẽ biểu đồ trực quan hóa (Odgi viz, Bandage, SequenceTubeMap, R plots). |
| **Đức** | **Project Architecture & Core Doc** | Xây dựng cấu trúc dự án, tổng hợp nội dung cơ bản và đóng góp tài liệu lý thuyết nền tảng. |

---

## 🎯 Mục Tiêu Dự Án (Project Objectives)

1.  **Lý thuyết nền tảng:** Hiểu sâu bản chất Pangenome, sự khác biệt giữa hệ gen tham chiếu tuyến tính và đồ thị biến dị (Variation Graph - GFA format).
2.  **Xây dựng Đồ thị (Graph Construction):** Thực hành làm chủ hai luồng thuật toán chính:
    *   `PGGB` (Reference-free): `wfmash` $\rightarrow$ `seqwish` $\rightarrow$ `smoothxg`.
    *   `Minigraph-Cactus` (Reference-centric).
3.  **Phân tích & Gọi biến dị (Downstream Analysis):**
    *   Lập chỉ mục đồ thị (`vg autoindex`) và căn chỉnh đọc ngắn FASTQ (`vg giraffe`).
    *   Chiếu ngược tọa độ (`vg surject`) và gọi biến dị cấu trúc (`vg deconstruct`, `vg call`).
4.  **Trực quan hóa chuẩn báo cáo:** Xuất biểu đồ trực quan chẩn đoán 1D/2D phục vụ báo cáo khoa học bài bản.

---

## 📂 Cấu Trúc Thư Mục Dự Án (Repository Structure)

```text
Introduction-Pangenom/
├── huong_dan_lam_viec.md        # Nguyên tắc hoạt động & Quy định dành cho nhóm/AI
├── README.md                    # Tài liệu giới thiệu dự án & phân công nhiệm vụ
├── documents/                   # Thư mục lưu trữ tài liệu hướng dẫn & lộ trình
│   └── images/                  # Thư mục hình ảnh minh họa trực quan
├── reference/                   # Tài liệu tham khảo, tệp PDF & hướng dẫn chuẩn
│   ├── REFERENCES.md            # Tổng hợp đầy đủ nguồn tham khảo & bài báo khoa học
│   └── GFA-spec/                # Đặc tả định dạng GFAv1 & GFAv2
├── src/                         # Kịch bản dòng lệnh Bash, Python & R (Mã nguồn)
│   ├── fasta_pansn/             # Kịch bản chuẩn hóa tên chuỗi PanSN (Bash, Python, SeqKit)
│   └── pggb/                    # Kịch bản chạy pipeline PGGB
├── data/                        # Thư mục chứa dữ liệu FASTA/FASTQ thực hành (Local)
│   ├── raw/                     # Thư mục lưu dữ liệu thô ban đầu (Reads/FASTQ/FASTA)
│   └── intern/                  # Thư mục lưu kết quả xử lý trung gian
├── out1/                        # Thư mục chứa biểu đồ trực quan hóa ODGI/PGGB
└── results/                     # Thư mục kết quả phân tích & báo cáo cuối cùng
```

---

## 🚀 Lộ Trình Thực Hành (Roadmap Quickstart)

*   **Giai đoạn 1:** Nền tảng tư duy & Chuẩn hóa định dạng PanSN (`<sample>#<haplotype>#<contig>`) cho các file FASTA nấm.
*   **Giai đoạn 2:** Chạy xây dựng đồ thị Pangenome bằng `PGGB` và `Minigraph-Cactus`.
*   **Giai đoạn 3:** Trực quan hóa bằng `odgi viz`, `Bandage`, `SequenceTubeMap` và viết script R.
*   **Giai đoạn 4:** Lập chỉ mục đồ thị, ánh xạ NGS Reads và gọi biến dị cấu trúc (SVs).

👉 Xem chi tiết tại [documents/lo_trinh_hoc_pangenome.md](file:///home/vkhang-bui/1.HocViec/theory_and_resources/pangenom/documents/lo_trinh_hoc_pangenome.md).

---

## 🛠️ Hướng Dẫn Cài Đặt Môi Trường (Environment Setup)

Dưới đây là câu lệnh 1 bước để khởi tạo môi trường làm việc **`shina`** và cài đặt toàn bộ các công cụ Pangenomics & Tin sinh học được sử dụng trong dự án (`PGGB`, `VG`, `ODGI`, `Minigraph`, `Bandage`, `Samtools`, `Biopython`, `R-tidyverse`,...):

### Cách 1: Cài đặt bằng Mamba (Nhanh & Tối ưu - Khuyên dùng)
```bash
# Tạo môi trường 'shina' và cài đặt tất cả các công cụ
mamba create -n shina -c bioconda -c conda-forge pggb vg odgi minigraph bandage samtools bcftools bedtools bwa fastp seqkit gfatools pigz python biopython pandas matplotlib seaborn r-base r-tidyverse r-ggplot2 -y

# Kích hoạt môi trường
mamba activate shina
```

### Cách 2: Cài đặt bằng Conda
```bash
# Tạo môi trường 'shina' và cài đặt tất cả các công cụ
conda create -n shina -c bioconda -c conda-forge pggb vg odgi minigraph bandage samtools bcftools bedtools bwa fastp seqkit gfatools pigz python biopython pandas matplotlib seaborn r-base r-tidyverse r-ggplot2 -y

# Kích hoạt môi trường
conda activate shina
```

---

## 📖 Tài Liệu Tham Khảo (References)

> 📌 **Danh mục tổng hợp đầy đủ:** Tất cả bài báo khoa học, tài liệu kỹ thuật và liên kết được lưu trữ và cập nhật tại [reference/REFERENCES.md](file:///home/vkhang-bui/1.HocViec/projects/pangenom/reference/REFERENCES.md).

1.  **Matthews, C. A., Watson-Haigh, N. S., Burton, R. A., & Sheppard, A. E. (2024).** *A gentle introduction to pangenomics.* Briefings in Bioinformatics, 25(6), bbae588. [https://doi.org/10.1093/bib/bbae588](https://doi.org/10.1093/bib/bbae588)
2.  **Garrison, E., Guarracino, A., Heumos, S., et al. (2023).** *Building pangenome graphs.* bioRxiv, 2023.04.05.535718. [https://doi.org/10.1101/2023.04.05.535718](https://doi.org/10.1101/2023.04.05.535718)
3.  **BioinfOmics INRAE. (2024).** *Giới thiệu hướng dẫn - Cuộc thi lập trình Pangenome Hackathon.*