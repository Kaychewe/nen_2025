#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: fig6g_cpg_table.R <per_cpg_tables_dir> <out_dir>")
}

in_dir <- args[1]
out_dir <- args[2]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(in_dir, pattern = "^table_.*_perCpG_pvals\\.csv$", full.names = TRUE)
if (length(files) == 0) {
  stop("No per-CpG tables found in input directory.")
}

extract_gene <- function(path) {
  fname <- basename(path)
  str_match(fname, "^table_(.*)_perCpG_pvals\\.csv$")[, 2]
}

summaries <- lapply(files, function(f) {
  gene <- extract_gene(f)
  df <- read_csv(f, show_col_types = FALSE)

  total <- nrow(df)
  sig_NETvN <- sum(df$FDR_NETvN < 0.05, na.rm = TRUE)
  sig_NECvNET <- sum(df$FDR_NECvNET < 0.05, na.rm = TRUE)
  sig_NECvN <- sum(df$FDR_NECvN < 0.05, na.rm = TRUE)

  tibble(
    Gene = gene,
    Category = c("Significant CpGs", "Non Significant CpGs"),
    `NET - Normal` = c(
      sprintf("%d (%.0f%%)", sig_NETvN, 100 * sig_NETvN / total),
      sprintf("%d (%.0f%%)", total - sig_NETvN, 100 * (total - sig_NETvN) / total)
    ),
    `NEC - NET` = c(
      sprintf("%d (%.0f%%)", sig_NECvNET, 100 * sig_NECvNET / total),
      sprintf("%d (%.0f%%)", total - sig_NECvNET, 100 * (total - sig_NECvNET) / total)
    ),
    `NEC - Normal` = c(
      sprintf("%d (%.0f%%)", sig_NECvN, 100 * sig_NECvN / total),
      sprintf("%d (%.0f%%)", total - sig_NECvN, 100 * (total - sig_NECvN) / total)
    )
  )
})

out <- bind_rows(summaries) %>%
  arrange(Gene, desc(Category))

out_csv <- file.path(out_dir, "fig6g_cpg_significance_by_gene.csv")
out_tsv <- file.path(out_dir, "fig6g_cpg_significance_by_gene.tsv")

write_csv(out, out_csv)
write_tsv(out, out_tsv)

cat("Wrote:\n", out_csv, "\n", out_tsv, "\n")
