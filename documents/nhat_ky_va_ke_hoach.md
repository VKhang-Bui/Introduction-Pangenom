# Nhật Ký Làm Việc & Kế Hoạch Hành Động (Project Log & Action Plan)

Tài liệu này ghi chép lại toàn bộ tiến độ làm việc, các cột mốc đã hoàn thành và kế hoạch triển khai tiếp theo cho dự án **Introduction-Pangenom**.

---

## 1. Nhật Ký Tiến Độ (Work Log)

### 📅 Ngày 27/07/2026: Đào Sâu Trực Quan Hóa 1D/2D ODGI & Cập Nhật Thư Mục Tài Liệu
* **Công việc đã hoàn thành:**
  1. Phân tích chi tiết toàn bộ các bức ảnh trực quan trong thư mục [out1](file:///home/vkhang-bui/1.HocViec/theory_and_resources/pangenom/out1) từ kết quả chạy thử nghiệm PGGB dữ liệu locus `DRB1-3123`.
  2. Giải thích bản chất toán học của thuật toán **PG-SGD (1D Projection)** trong ODGI.
  3. Xây dựng quy trình 4 bước chuẩn (**4-Step Cheatsheet**) để đọc và giải mã bất kỳ sơ đồ ma trận 1D `odgi viz` nào.
  4. Phân tích sinh học 6 sơ đồ ma trận 1D (`viz_multiqc`, `viz_depth`, `viz_inv`, `viz_pos`, `viz_uncalled`, `viz_O`), sơ đồ bố cục 2D (`draw_multiqc`) và hệ thống đường móc đen phía dưới (Bottom Edge Arcs).
  5. Phát hiện đặc tính di truyền quan trọng: Locus `DRB1-3123` có **Core Genome = 0%**, hoàn toàn cấu thành từ 5 phân khu biến thể cấu trúc (Accessory Blocks).
  6. Xây dựng kịch bản báo cáo chuyên nghiệp chuẩn khoa học cho Chủ nhiệm dự án.
  7. Khởi tạo và sao lưu toàn bộ hình ảnh minh họa vào thư mục [documents/images](file:///home/vkhang-bui/1.HocViec/theory_and_resources/pangenom/documents/images).
  8. Cập nhật cẩm nang tài liệu hoàn chỉnh tại [documents/odgi_1d_viz_guide.md](file:///home/vkhang-bui/1.HocViec/theory_and_resources/pangenom/documents/odgi_1d_viz_guide.md) tích hợp liên kết hình ảnh trực quan.

---

## 2. Kế Hoạch Hành Động Tiếp Theo (Action Plan)

* **Giai đoạn 3 (Trực quan hóa nâng cao & Báo cáo):**
  - [ ] Thực hành trích xuất đồ thị 2D bằng **Bandage** và đồ thị con (subgraphs) bằng **SequenceTubeMap**.
  - [ ] Thử nghiệm viết script R custom để vẽ ma trận pangenome chuẩn xuất bản bài báo khoa học.

* **Giai đoạn 4 (Phân tích hạ nguồn NGS):**
  - [ ] Thực hành lập chỉ mục đồ thị bằng `vg autoindex` (tạo `.gbz`, `.dist`, `.min`).
  - [ ] Thực hành ánh xạ đọc ngắn NGS FASTQ bằng `vg giraffe`.
  - [ ] Gọi biến thể cấu trúc (SVs) bằng `vg deconstruct` và `vg call`.

* **Ứng dụng trên Nấm (Fungi):**
  - [ ] Áp dụng quy trình chuẩn PanSN + PGGB lên tập dữ liệu FASTA bộ gen Nấm thực tế.

---

## Nguồn Tham Khảo (References)

1. **Guarracino, A., Heumos, S., Garrison, E., et al. (2022).** *ODGI: understanding pangenome graphs.* Bioinformatics, 38(13), 3319–3326. [https://doi.org/10.1093/bioinformatics/btac308](https://doi.org/10.1093/bioinformatics/btac308)
2. **Matthews, C. A., Watson-Haigh, N. S., Burton, R. A., & Sheppard, A. E. (2024).** *A gentle introduction to pangenomics.* Briefings in Bioinformatics, 25(6), bbae588. [https://doi.org/10.1093/bib/bbae588](https://doi.org/10.1093/bib/bbae588)
