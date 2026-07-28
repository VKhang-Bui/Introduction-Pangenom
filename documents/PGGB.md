# Hướng Dẫn Kỹ Thuật & Quy Trình Đầu Vào Cho PGGB (Pangenome Graph Builder)

Tài liệu này tổng hợp kiến thức nền tảng, định dạng dữ liệu đầu vào chuẩn và quy trình xử lý dữ liệu trước khi đưa vào **PGGB (Pangenome Graph Builder)**.

---

## 1. Giới Thiệu & Nhóm Phát Triển

*   **Tên công cụ:** PGGB (Pangenome Graph Builder).
*   **Nhóm phát triển:** Được phát triển bởi **TS. Erik Garrison** (Đại học Tennessee Health Science Center - UTHSC), **TS. Andrea Guarracino** (UTHSC / HPRC) cùng các cộng sự thuộc Liên minh Tham chiếu Pangenome Người (**Human Pangenome Reference Consortium - HPRC**).
*   **Mục đích chính:** Dựng đồ thị pangenome (Variation Graph - GFA) **hoàn toàn không phụ thuộc vào hệ gen tham chiếu (Reference-free)**, coi mọi haplotype đầu vào có vai trò di truyền bình đẳng và bảo tồn đầy đủ các biến dị cấu trúc (SVs).

---

## 2. Bản Chất Dữ Liệu Đầu Vào Của PGGB

Đầu vào chính thức của PGGB là **một tệp multi-FASTA** (nén `.fa.gz` hoặc `.fasta.gz`) cùng tệp chỉ mục `.fai` đi kèm.

> [!IMPORTANT]
> **Yêu cầu nguồn gốc dữ liệu FASTA:**
> Tệp FASTA đầu vào phải là kết quả **Lắp Ráp De Novo (De Novo Assembly)** từ dữ liệu đọc dài (**Long-Reads**: PacBio HiFi hoặc Oxford Nanopore). 
> PGGB **không nhận trực tiếp tệp đọc thô FASTQ**.

---

## 3. Quy Chuẩn Đặt Tên Chuỗi PanSN (PanSN Convention)

Để PGGB nhận diện chính xác mẫu và đường đi (Walks/Paths) trên đồ thị, tiêu đề (header) của từng chuỗi FASTA trong tệp gộp **bắt buộc** phải tuân theo chuẩn **PanSN (Pangenome Sequence Naming)**:

$$\text{><sample>\#<haplotype>\#<contig\_name>}$$

*   `sample`: Tên mẫu / chủng sinh vật (ví dụ: `IPO323`, `Zt09`).
*   `haplotype`: Chỉ số bộ đơn / bộ nhiễm sắc thể (ví dụ: `1` cho đơn bội, `1` hoặc `2` cho lưỡng bội).
*   `contig_name`: Tên đoạn contig hoặc nhiễm sắc thể (ví dụ: `chr01`, `ctg0001`).

*Ví dụ tiêu đề chuẩn:*
```text
>IPO323#1#chr01
ACGTACGT...
>Zt09#1#chr01
ACGTACGT...
```

---

## 4. Quy Trình Tổng Quan Xử Lý Dữ Liệu (Bioinformatics Workflow)

Dưới đây là luồng công việc chuẩn từ tệp FASTQ thô đến đồ thị PGGB:

```text
┌─────────────────────────────────────────────────────────┐
│ 1. DỮ LIỆU ĐẦU VÀO THÔ                                  │
│    Tệp FASTQ đọc dài (PacBio HiFi / Oxford Nanopore)    │
└────────────────────────────┬────────────────────────────┘
                             │
                             ▼
┌─────────────────────────────────────────────────────────┐
│ 2. LẮP RÁP DE NOVO (Assembly)                           │
│    Dùng hifiasm / Flye ➔ Xuất tệp FASTA từng mẫu riêng  │
└────────────────────────────┬────────────────────────────┘
                             │
                             ▼
┌─────────────────────────────────────────────────────────┐
│ 3. CHUẨN HÓA PANSN & GỘP TỆP                            │
│    • Đổi tên header ➔ >Sample#Hap#Contig                │
│    • Gộp tệp FASTA, nén bgzip & tạo chỉ mục .fai        │
└────────────────────────────┬────────────────────────────┘
                             │
                             ▼
┌─────────────────────────────────────────────────────────┐
│ 4. KHỞI CHẠY PGGB PIPELINE                              │
│    PGGB Engine (wfmash ➔ seqwish ➔ smoothxg)            │
│    ➔ Xuất tệp Đồ thị Pangenome GFA đầu ra (.gfa)        │
└─────────────────────────────────────────────────────────┘
```

### Các bước câu lệnh thực thi mẫu:

1.  **Lắp ráp De Novo từng mẫu (với PacBio HiFi):**
    ```bash
    hifiasm -o sampleA.asm -t 16 sampleA.fastq.gz
    gfatools gfa2fa sampleA.asm.bp.p_ctg.gfa > sampleA.fasta
    ```

2.  **Đổi tên tiêu đề FASTA theo chuẩn PanSN:**
    ```bash
    sed -i 's/>/>sampleA#1#/g' sampleA.fasta
    ```

3.  **Gộp tệp, nén bgzip và tạo chỉ mục fai:**
    ```bash
    cat *.fasta > all_samples.fasta
    bgzip -@ 8 all_samples.fasta
    samtools faidx all_samples.fasta.gz
    ```

4.  **Khởi chạy PGGB:**
    ```bash
    pggb -i all_samples.fasta.gz -p 98 -s 5000 -t 16 -o pggb_output/
    ```

---

## Tài Liệu Tham Khảo (References)

1. **Garrison, E., Guarracino, A., Heumos, S., et al. (2023).** *Building pangenome graphs.* bioRxiv, 2023.04.05.535718. [https://doi.org/10.1101/2023.04.05.535718](https://doi.org/10.1101/2023.04.05.535718)
2. **PanSN Specification (HPRC):** *Pangenome Sequence Naming standard convention.* GitHub Repository. [https://github.com/pangenome/PanSN-spec](https://github.com/pangenome/PanSN-spec)
3. **Cheng, H., Concepcion, G. T., Feng, X., Zhang, H., & Li, H. (2021).** *Haplotype-resolved de novo assembly of diploid genomes with hifiasm.* Nature Methods, 18(2), 170-175. [https://doi.org/10.1038/s41592-020-01056-5](https://doi.org/10.1038/s41592-020-01056-5)
4. **Tài liệu dự án liên quan:** [lo_trinh_hoc_pangenome.md](file:///home/vkhang-bui/1.HocViec/theory_and_resources/pangenom/documents/lo_trinh_hoc_pangenome.md) và [pangenome_theory_guide.md](file:///home/vkhang-bui/1.HocViec/theory_and_resources/pangenom/documents/pangenome_theory_guide.md).
