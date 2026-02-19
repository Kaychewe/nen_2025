#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(ggprism)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  cat("Usage: fig1b.R <NECPath.all.data.gene.h.tsv> [outdir]\n")
  quit(status = 1)
}

infile <- args[1]
outdir <- ifelse(length(args) >= 2, args[2], "figures")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# Load full NanoString dataset (36 samples × 784 genes)
raw_counts <- read_tsv(infile, show_col_types = FALSE)

# Set gene names as rownames
counts_mat <- raw_counts %>%
  column_to_rownames("Name") %>%
  mutate(across(everything(), ~ as.numeric(str_trim(.)))) %>%
  as.matrix()

stopifnot(all(is.finite(counts_mat)))

sample_info <- tibble(
  sample_id = colnames(counts_mat),
  nen_subtype = case_when(
    str_detect(sample_id, "_LG") ~ "WD",
    str_detect(sample_id, "_HG") ~ "PD",
    TRUE ~ NA_character_
  )
)

expr_log <- log2(counts_mat + 1)
expr_log_t <- t(expr_log)

pca_input <- expr_log_t %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  left_join(sample_info, by = "sample_id")

pca_res <- prcomp(
  pca_input %>% select(-sample_id, -nen_subtype),
  center = TRUE,
  scale. = TRUE
)

var_explained <- (pca_res$sdev^2) / sum(pca_res$sdev^2)
pc1_var <- round(var_explained[1] * 100, 1)
pc2_var <- round(var_explained[2] * 100, 1)

pca_plot_df <- data.frame(
  sample_id   = pca_input$sample_id,
  nen_subtype = pca_input$nen_subtype,
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2]
)

perm_text <- "PERMANOVA: R² = 0.17, F = 6.9, p = 0.001"

fig1b <- ggplot(pca_plot_df, aes(PC1, PC2, color = nen_subtype)) +
  geom_point(size = 4.5, alpha = 0.6) +
  stat_ellipse(type = "norm", linewidth = 1) +
  scale_color_manual(
    values = c("WD" = "#2C7BB6", "PD" = "#ff7f0eff"),
    name = "NEN subtype"
  ) +
  labs(
    x = paste0("PC1 (", pc1_var, "% variance)"),
    y = paste0("PC2 (", pc2_var, "% variance)")
  ) +
  annotate(
    "label",
    x = Inf, y = -Inf,
    label = perm_text,
    hjust = 1.05,
    vjust = -0.6,
    size = 4,
    fill = "white"
  ) +
  theme_prism(border = TRUE, base_rect_size = 2.5) +
  theme(
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_blank()
  )

# Save outputs
png_path <- file.path(outdir, "fig1b.png")
svg_path <- file.path(outdir, "fig1b.svg")
pdf_path <- file.path(outdir, "fig1b.pdf")
jpg_path <- file.path(outdir, "fig1b.jpg")

ggsave(plot = fig1b, filename = png_path, dpi = 1200, width = 10, height = 8)
if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(plot = fig1b, filename = svg_path, dpi = 1200, width = 10, height = 8)
} else {
  message("Package 'svglite' not available; using base svg device.")
  grDevices::svg(filename = svg_path, width = 10, height = 8)
  print(fig1b)
  grDevices::dev.off()
}
ggsave(plot = fig1b, filename = pdf_path, dpi = 1200, width = 10, height = 8)
ggsave(plot = fig1b, filename = jpg_path, dpi = 1200, width = 10, height = 8)

cat("Wrote:\n", png_path, "\n", pdf_path, "\n", jpg_path, "\n")
cat(svg_path, "\n")
