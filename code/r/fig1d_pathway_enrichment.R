#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(ggprism)
  library(ggh4x)
  library(edgeR)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})


input_path <- "../../data/NECPath.all.data.gene.h.tsv"
out_dir <- "../../figures/individuals/fig1/"

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

# Define significance thresholds
pd_genes <- de %>%
  filter(logFC > 1, FDR < 0.05) %>%
  pull(gene)

wd_genes <- de %>%
  filter(logFC < -1, FDR < 0.05) %>%
  pull(gene)

# Background = all genes measured on the NanoString panel
# (Use all genes in de table)

gene_universe <- unique(de$gene)

# Enrichment: GO BP

ego_pd <- enrichGO(
  gene          = pd_genes,
  universe      = gene_universe,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

ego_wd <- enrichGO(
  gene          = wd_genes,
  universe      = gene_universe,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

clean_go_names <- function(x) {
  x %>%
    str_replace_all("^DNA-templated DNA replication$", "DNA replication") %>%
    str_replace_all("double-strand break repair via homologous recombination",
                    "Homologous recombination repair") %>%
    str_replace_all("^cell cycle process$", "Cell cycle") %>%
    str_replace_all("^synaptic signaling$", "Synaptic signaling") %>%
    str_replace_all("^trans-synaptic signaling$", "Trans-synaptic signaling") %>%
    str_replace_all("^nervous system process$", "Nervous system process") %>%
    str_replace_all("^negative regulation of epithelial cell differentiation$",
                    "Inhibition of epithelial differentiation") %>%
    str_replace_all("^regulation of epithelial cell differentiation$",
                    "Regulation of epithelial differentiation")
}

pd_df <- ego_pd@result %>%
  as_tibble() %>%
  mutate(group = "PD-NEC")

wd_df <- ego_wd@result %>%
  as_tibble() %>%
  mutate(group = "WD-NET")

enrich_df <- bind_rows(pd_df, wd_df) %>%
  mutate(
    Description = clean_go_names(Description),
    log10FDR = -log10(p.adjust)
  )

enrich_df$group <- factor(enrich_df$group, levels = c("PD-NEC", "WD-NET"))

# Select top pathways per group
top_enrich <- enrich_df %>%
  group_by(group) %>%
  slice_max(order_by = log10FDR, n = 8) %>%
  ungroup()

fig1d <- ggplot(
  top_enrich,
  aes(
    x = log10FDR,
    y = reorder(Description, log10FDR),
    size = Count,
    color = group
  )
) +
  geom_point(alpha = 0.9) +
  facet_wrap(~ group, scales = "free_y", drop = FALSE) +
  scale_color_manual(
    values = c("PD-NEC" = "#ff7f0eff", "WD-NET" = "#2C7BB6"),
    guide = "none"
  ) +
  scale_size(range = c(2.5, 6)) +
  labs(
    x = expression(-log[10]("FDR-adjusted p value")),
    y = NULL,
    title = "Pathway enrichment of subtype-specific genes"
  ) +
  theme_prism(border = TRUE, base_rect_size = 1.2) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 9)
  )

pd_enrich <- top_enrich %>% filter(group == "PD-NEC")
wd_enrich <- top_enrich %>% filter(group == "WD-NET")

div_enrich <- bind_rows(
  pd_enrich %>% mutate(signed_log10FDR = log10FDR),
  wd_enrich %>% mutate(signed_log10FDR = -log10FDR)
)

fig1d_bar_div <- ggplot(
  div_enrich,
  aes(
    x = signed_log10FDR,
    y = reorder(Description, signed_log10FDR),
    fill = group
  )
) +
  geom_col(alpha = 0.9, width = 0.7, linewidth = 0.6, color = "black") +
  scale_fill_manual(
    values = c("PD-NEC" = "#ff7f0eff", "WD-NET" = "#2C7BB6"),
    name = NULL
  ) +
  labs(
    x = expression(-log[10]("FDR-adjusted p value")),
    y = NULL,
    title = "Pathway enrichment of subtype-specific genes"
  ) +
  theme_prism(border = TRUE, base_rect_size = 1.2) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text.y = element_text(size = 9),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.8),
    panel.border = element_rect(linewidth = 0.8, fill = NA),
    legend.position = c(0.98, 0.02),
    legend.justification = c(1, 0),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
    legend.key.size = unit(0.35, "cm")
  )

fig1d_bar_pd <- ggplot(
  pd_enrich,
  aes(
    x = log10FDR,
    y = reorder(Description, log10FDR)
  )
) +
  geom_col(alpha = 0.9, width = 0.7, fill = "#ff7f0eff") +
  scale_x_continuous(limits = c(0, 6)) +
  labs(
    x = expression(-log[10]("FDR-adjusted p value")),
    y = NULL,
    title = "PD-NEC"
  ) +
  theme_prism(border = TRUE, base_rect_size = 1.2) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.text.y = element_text(size = 9)
  )

fig1d_bar_wd <- ggplot(
  wd_enrich,
  aes(
    x = log10FDR,
    y = reorder(Description, log10FDR)
  )
) +
  geom_col(alpha = 0.9, width = 0.7, fill = "#2C7BB6") +
  scale_x_continuous(limits = c(0, 3)) +
  labs(
    x = expression(-log[10]("FDR-adjusted p value")),
    y = NULL,
    title = "WD-NET"
  ) +
  theme_prism(border = TRUE, base_rect_size = 1.2) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.text.y = element_text(size = 9)
  )

png_path <- file.path(out_dir, "fig1d.png")
pdf_path <- file.path(out_dir, "fig1d.pdf")
svg_path <- file.path(out_dir, "fig1d.svg")
jpg_path <- file.path(out_dir, "fig1d.jpg")

bar_pd_png_path <- file.path(out_dir, "fig1d_bar_pd.png")
bar_pd_pdf_path <- file.path(out_dir, "fig1d_bar_pd.pdf")
bar_pd_svg_path <- file.path(out_dir, "fig1d_bar_pd.svg")
bar_pd_jpg_path <- file.path(out_dir, "fig1d_bar_pd.jpg")

bar_wd_png_path <- file.path(out_dir, "fig1d_bar_wd.png")
bar_wd_pdf_path <- file.path(out_dir, "fig1d_bar_wd.pdf")
bar_wd_svg_path <- file.path(out_dir, "fig1d_bar_wd.svg")
bar_wd_jpg_path <- file.path(out_dir, "fig1d_bar_wd.jpg")

bar_div_png_path <- file.path(out_dir, "fig1d_bar_div.png")
bar_div_pdf_path <- file.path(out_dir, "fig1d_bar_div.pdf")
bar_div_svg_path <- file.path(out_dir, "fig1d_bar_div.svg")
bar_div_jpg_path <- file.path(out_dir, "fig1d_bar_div.jpg")

# Write outputs

ggsave(plot = fig1d, filename = png_path, dpi = 600, width = 10, height = 6)
if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(plot = fig1d, filename = svg_path, width = 10, height = 6)
} else {
  message("Package 'svglite' not available; using base svg device.")
  grDevices::svg(filename = svg_path, width = 10, height = 6)
  print(fig1d)
  grDevices::dev.off()
}
ggsave(plot = fig1d, filename = pdf_path, width = 10, height = 6)
ggsave(plot = fig1d, filename = jpg_path, dpi = 600, width = 10, height = 6)

ggsave(plot = fig1d_bar_pd, filename = bar_pd_png_path, dpi = 600, width = 8, height = 6)
if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(plot = fig1d_bar_pd, filename = bar_pd_svg_path, width = 8, height = 6)
} else {
  message("Package 'svglite' not available; using base svg device.")
  grDevices::svg(filename = bar_pd_svg_path, width = 8, height = 6)
  print(fig1d_bar_pd)
  grDevices::dev.off()
}
ggsave(plot = fig1d_bar_pd, filename = bar_pd_pdf_path, width = 8, height = 6)
ggsave(plot = fig1d_bar_pd, filename = bar_pd_jpg_path, dpi = 600, width = 8, height = 6)

cat("Wrote:\n", png_path, "\n", pdf_path, "\n", jpg_path, "\n", svg_path, "\n")
cat("Wrote:\n", bar_pd_png_path, "\n", bar_pd_pdf_path, "\n", bar_pd_jpg_path, "\n", bar_pd_svg_path, "\n")

ggsave(plot = fig1d_bar_wd, filename = bar_wd_png_path, dpi = 600, width = 8, height = 6)
if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(plot = fig1d_bar_wd, filename = bar_wd_svg_path, width = 8, height = 6)
} else {
  message("Package 'svglite' not available; using base svg device.")
  grDevices::svg(filename = bar_wd_svg_path, width = 8, height = 6)
  print(fig1d_bar_wd)
  grDevices::dev.off()
}
ggsave(plot = fig1d_bar_wd, filename = bar_wd_pdf_path, width = 8, height = 6)
ggsave(plot = fig1d_bar_wd, filename = bar_wd_jpg_path, dpi = 600, width = 8, height = 6)

cat("Wrote:\n", bar_wd_png_path, "\n", bar_wd_pdf_path, "\n", bar_wd_jpg_path, "\n", bar_wd_svg_path, "\n")

ggsave(plot = fig1d_bar_div, filename = bar_div_png_path, dpi = 600, width = 8, height = 6)
if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(plot = fig1d_bar_div, filename = bar_div_svg_path, width = 8, height = 6)
} else {
  message("Package 'svglite' not available; using base svg device.")
  grDevices::svg(filename = bar_div_svg_path, width = 8, height = 6)
  print(fig1d_bar_div)
  grDevices::dev.off()
}
ggsave(plot = fig1d_bar_div, filename = bar_div_pdf_path, width = 8, height = 6)
ggsave(plot = fig1d_bar_div, filename = bar_div_jpg_path, dpi = 600, width = 8, height = 6)

cat("Wrote:\n", bar_div_png_path, "\n", bar_div_pdf_path, "\n", bar_div_jpg_path, "\n", bar_div_svg_path, "\n")
