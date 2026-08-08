#!/usr/bin/env Rscript

# nucmer2vcf.R: Converts MUMmer show-snps (-THC) output to VCF 4.2
args <- commandArgs(trailingOnly = TRUE)
var_file     <- args[1]
query_contig <- args[2]
ref_fasta    <- args[3]
nucmer_ver   <- args[4]
output_vcf   <- args[5]

fai_df <- read.table(paste0(ref_fasta, ".fai"), stringsAsFactors = FALSE, comment.char = "")
ref_chrom  <- fai_df[1, 1]
ref_length <- fai_df[1, 2]

vcf_conn <- file(output_vcf, "w")
cat("##fileformat=VCFv4.2\n", file = vcf_conn)
cat(paste0("##source=nucmer_", nucmer_ver, "\n"), file = vcf_conn)
cat(paste0("##contig=<ID=", ref_chrom, ",length=", ref_length, ">\n"), file = vcf_conn)
cat("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n", file = vcf_conn)
cat("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t", query_contig, "\n", sep = "", file = vcf_conn)

if (file.exists(var_file) && file.info(var_file)$size > 0) {
  var_data <- read.table(var_file, header = FALSE, stringsAsFactors = FALSE, sep = "\t", quote = "", comment.char = "")
  if (nrow(var_data) > 0) {
    for (i in 1:nrow(var_data)) {
      pos <- as.numeric(var_data[i, 1])
      ref <- trimws(as.character(var_data[i, 2]))
      alt <- trimws(as.character(var_data[i, 3]))
      ref_val <- if (ref == ".") "N" else ref
      alt_val <- if (alt == ".") "N" else alt
      cat(sprintf("%s\t%d\t.\t%s\t%s\t.\tPASS\t.\tGT\t1\n", ref_chrom, pos, ref_val, alt_val), sep = "", file = vcf_conn)
    }
  }
}
close(vcf_conn)
