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
│   ├── pangenome_theory_guide.md# Cẩm nang lý thuyết & phân tích Pangenome
│   └── lo_trinh_hoc_pangenome.md # Lộ trình 4 giai đoạn làm chủ Pangenomics
├── reference/                   # Tài liệu tham khảo, tệp PDF & hướng dẫn chuẩn
│   ├── Matthews et al. - 2024 - A gentle introduction to pangenomics.pdf
│   └── GFA-spec/                # Đặc tả định dạng GFAv1 & GFAv2
├── src/                         # Kịch bản dòng lệnh Bash, Python & R (Mã nguồn)
└── data/                        # Thư mục chứa dữ liệu FASTA/FASTQ thực hành (Local)
```

---

## 🚀 Lộ Trình Thực Hành (Roadmap Quickstart)

*   **Giai đoạn 1:** Nền tảng tư duy & Chuẩn hóa định dạng PanSN (`<sample>#<haplotype>#<contig>`) cho các file FASTA nấm.
*   **Giai đoạn 2:** Chạy xây dựng đồ thị Pangenome bằng `PGGB` và `Minigraph-Cactus`.
*   **Giai đoạn 3:** Trực quan hóa bằng `odgi viz`, `Bandage`, `SequenceTubeMap` và viết script R.
*   **Giai đoạn 4:** Lập chỉ mục đồ thị, ánh xạ NGS Reads và gọi biến dị cấu trúc (SVs).

👉 Xem chi tiết tại [documents/lo_trinh_hoc_pangenome.md](file:///mnt/d/1.HocViec/theory_and_resources/pangenom/documents/lo_trinh_hoc_pangenome.md).

---

## 📖 Tài Liệu Tham Khảo (References)

1.  **Matthews, C. A., Watson-Haigh, N. S., Burton, R. A., & Sheppard, A. E. (2024).** *A gentle introduction to pangenomics.* Briefings in Bioinformatics, 25(6), bbae588. [https://doi.org/10.1093/bib/bbae588](https://doi.org/10.1093/bib/bbae588)
2.  **Garrison, E., Guarracino, A., Heumos, S., et al. (2023).** *Building pangenome graphs.* bioRxiv, 2023.04.05.535718. [https://doi.org/10.1101/2023.04.05.535718](https://doi.org/10.1101/2023.04.05.535718)
3.  **BioinfOmics INRAE. (2024).** *Giới thiệu hướng dẫn - Cuộc thi lập trình Pangenome Hackathon.*