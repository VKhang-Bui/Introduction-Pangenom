#!/usr/bin/env Rscript

# plot_small_variants.R: Visualizes Precision, Recall, and F1-score
suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
stat_file  <- args[1]
output_png <- args[2]

if (!file.exists(stat_file) || file.info(stat_file)$size == 0) {
  cat("Warning: Statistics file is empty or missing. Skipping plot.\n")
  quit(status = 0)
}

stat_df <- read.table(stat_file, sep = '\t', header = TRUE, comment.char = '', stringsAsFactors = FALSE)
colnames(stat_df) <- gsub("\\.", "_", colnames(stat_df))

if (nrow(stat_df) == 0) {
  cat("Warning: No rows found in statistics file. Skipping plot.\n")
  quit(status = 0)
}

stat_df_long <- pivot_longer(stat_df, cols = c("precision", "recall", "f1_score"), names_to = "Metric", values_to = "Value")
stat_df_long$Value <- as.numeric(as.character(stat_df_long$Value))
stat_df_long$Value[is.na(stat_df_long$Value)] <- 0

p <- ggplot(stat_df_long, aes(x = contig, y = Value, fill = contig)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Metric, ncol = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

ggsave(output_png, plot = p, width = 8, height = 6, dpi = 300)
cat("Plot successfully saved to:", output_png, "\n")
