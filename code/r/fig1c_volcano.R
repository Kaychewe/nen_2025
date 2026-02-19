#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(ggprism)
  library(edgeR)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: fig1c_volcano.R <input_counts_tsv> <output_dir>")
}

input_path <- args[1]
out_dir <- args[2]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Load counts (36 samples x 784 genes)
raw_counts <- read_tsv(input_path, show_col_types = FALSE)

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

if (anyNA(sample_info$nen_subtype)) {
  stop("Found samples with unknown subtype (expected _LG or _HG in sample_id).")
}

# edgeR differential expression: PD vs WD
sample_info <- sample_info %>%
  mutate(nen_subtype = factor(nen_subtype, levels = c("WD", "PD")))

y <- DGEList(counts = counts_mat, group = sample_info$nen_subtype)

y <- calcNormFactors(y)

design <- model.matrix(~ group, data = y$samples)

y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)

qlf <- glmQLFTest(fit, coef = "groupPD")

de <- topTags(qlf, n = Inf)$table %>%
  tibble::rownames_to_column("gene") %>%
  mutate(
    FDR = p.adjust(PValue, method = "BH"),
    direction = ifelse(logFC > 0, "PD_up", "WD_up")
  )

volcano_df <- de %>%
  mutate(
    neg_log10_FDR = -log10(FDR),
    sig = case_when(
      FDR < 0.05 & logFC > 1  ~ "PD-up",
      FDR < 0.05 & logFC < -1 ~ "WD-up",
      TRUE                   ~ "NS"
    )
  )

volcano_base <- ggplot(volcano_df, aes(logFC, neg_log10_FDR)) +
  geom_point(aes(color = sig), alpha = 0.6, size = 2.5) +
  scale_color_manual(
    values = c(
      "PD-up" = "#D7191C",
      "WD-up" = "#2C7BB6",
      "NS"    = "grey70"
    ),
    name = "Differential expression"
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.5) +
  labs(
    x = "log2FC(PD/WD)",
    y = expression(-log[10]~FDR)
  ) +
  theme_prism(base_size = 16, border = TRUE, base_rect_size = 1.5) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

genes_to_label <- c(
  "FANCA", "TTK", "CDC6", "CACNB2",
  "E2F1", "RASA4", "PAK3",
  "EZH2", "CHEK1", "SFN",
  "FGF14", "MKI67", "TOP2A", "AURKA",
  "CAMK2B", "GRIA3"
)

volcano_labeled <- volcano_base +
  geom_text_repel(
    data = volcano_df %>% filter(gene %in% genes_to_label),
    aes(label = gene),
    size = 3.5,
    max.overlaps = 20,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "grey50"
  )

n_pd <- sum(volcano_df$sig == "PD-up")
n_wd <- sum(volcano_df$sig == "WD-up")

volcano_final <- volcano_labeled +
  annotate(
    "text",
    x = -3.5, y = max(volcano_df$neg_log10_FDR),
    label = paste0("WD-up: ", n_wd),
    hjust = 0, size = 4, color = "#2C7BB6"
  ) +
  annotate(
    "text",
    x = 2.5, y = max(volcano_df$neg_log10_FDR),
    label = paste0("PD-up: ", n_pd),
    hjust = 0, size = 4, color = "#D7191C"
  )

png_path <- file.path(out_dir, "fig1c.png")
pdf_path <- file.path(out_dir, "fig1c.pdf")
svg_path <- file.path(out_dir, "fig1c.svg")
jpg_path <- file.path(out_dir, "fig1c.jpg")

# Write outputs

ggsave(plot = volcano_final, filename = png_path, dpi = 1200, width = 11, height = 8)
if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(plot = volcano_final, filename = svg_path, width = 11, height = 8)
} else {
  message("Package 'svglite' not available; using base svg device.")
  grDevices::svg(filename = svg_path, width = 11, height = 8)
  print(volcano_final)
  grDevices::dev.off()
}
ggsave(plot = volcano_final, filename = pdf_path, width = 11, height = 8)
ggsave(plot = volcano_final, filename = jpg_path, dpi = 1200, width = 11, height = 8)

cat("Wrote:\n", png_path, "\n", pdf_path, "\n", jpg_path, "\n", svg_path, "\n")
