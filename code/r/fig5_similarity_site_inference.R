#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readxl)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: fig5_similarity_site_inference.R <expr_tsv> <meta_xlsx> <ml_split_tsv> <out_dir>")
}

expr_path <- args[1]
meta_xlsx <- args[2]
split_path <- args[3]
out_dir <- args[4]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
table_dir <- file.path(dirname(out_dir), "tables", "fig5")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Load expression (log2 counts) ----
expr <- read_tsv(expr_path, show_col_types = FALSE)
if (!"sample_id" %in% names(expr)) {
  stop("Expression file must contain a 'sample_id' column.")
}
expr_mat <- expr %>% select(-sample_id) %>% as.matrix()
rownames(expr_mat) <- expr$sample_id

# Normalize expression sample IDs for matching (strip _LG/_HG suffix)
expr_ids <- rownames(expr_mat)
expr_key <- str_replace(expr_ids, "(_LG|_HG)$", "")

# ---- Load metadata (Table S1) with header detection ----
meta_raw <- read_excel(meta_xlsx, sheet = 1, col_names = FALSE) %>% as.data.frame()
header_row <- which(apply(meta_raw, 1, function(r) any(str_detect(tolower(as.character(r)), "accession_number|sample_id|sample"))))[1]
if (is.na(header_row)) stop("Could not detect header row in metadata.")
header <- as.character(unlist(meta_raw[header_row, , drop = TRUE]))
meta <- meta_raw[(header_row + 1):nrow(meta_raw), , drop = FALSE]
header[is.na(header) | header == ""] <- paste0("X", which(is.na(header) | header == ""))
header <- make.unique(header)
names(meta) <- header
meta <- meta %>% filter(if_any(everything(), ~ !is.na(.x)))

sample_col <- names(meta)[str_detect(tolower(names(meta)), "sample|id|accession")][1]
if (is.na(sample_col)) stop("Could not detect sample ID column in metadata.")
meta <- meta %>% mutate(sample_id = as.character(.data[[sample_col]]))
meta <- meta %>% mutate(sample_key = str_replace(sample_id, "(_LG|_HG)$", ""))

# ---- Join split info ----
split_df <- read_tsv(split_path, show_col_types = FALSE) %>%
  mutate(sample_key = str_replace(sample_id, "(_LG|_HG)$", ""))
meta <- meta %>% left_join(split_df %>% select(sample_key, split), by = "sample_key")

# ---- Primary site + subtype ----
primary_col <- names(meta)[str_detect(tolower(names(meta)), "primary")][1]
if (is.na(primary_col)) stop("Could not detect primary_site column in metadata.")
meta <- meta %>% mutate(primary_site = as.character(.data[[primary_col]]))

subtype_col <- names(meta)[str_detect(tolower(names(meta)), "subtype|nen")][1]
if (is.na(subtype_col)) {
  meta <- meta %>% mutate(nen_subtype = ifelse(str_detect(sample_id, "_HG"), "PD", "WD"))
} else {
  meta <- meta %>% mutate(nen_subtype = as.character(.data[[subtype_col]]))
}

# ---- Match samples ----
idx <- match(expr_key, meta$sample_key)
keep <- which(!is.na(idx))
if (length(keep) == 0) {
  stop("No matching samples between expression matrix and metadata (after stripping _LG/_HG).")
}
expr_mat <- expr_mat[keep, , drop = FALSE]
meta <- meta[idx[keep], , drop = FALSE]
meta$expr_sample_id <- expr_ids[keep]
meta$sample_id <- meta$sample_key
rownames(expr_mat) <- meta$expr_sample_id

# ---- PCA ----
expr_scaled <- scale(expr_mat)
pca <- prcomp(expr_scaled, center = TRUE, scale. = FALSE)
pca_df <- as.data.frame(pca$x[, 1:2])
pca_df$sample_id <- rownames(pca_df)

# ---- KNN primary-site inference (train on known primary_site) ----
known <- !is.na(meta$primary_site) & meta$primary_site != ""
train_idx <- which(known)
if (length(train_idx) < 5) stop("Not enough known primary_site samples for KNN training.")

# Use first 10 PCs for KNN
pc_use <- 10
pc_mat <- as.data.frame(pca$x[, 1:pc_use])
pc_mat$sample_id <- rownames(pc_mat)

train_x <- pc_mat[train_idx, 1:pc_use, drop = FALSE]
train_y <- meta$primary_site[train_idx]
test_x <- pc_mat[, 1:pc_use, drop = FALSE]

library(class)
set.seed(332)
k <- 5
pred_site <- knn(train = train_x, test = test_x, cl = train_y, k = k)

meta <- meta %>% mutate(primary_site_inferred = pred_site)

# ---- PCA plot ----
plot_df <- pca_df %>%
  left_join(meta %>% select(expr_sample_id, primary_site_inferred, nen_subtype, split), by = c("sample_id" = "expr_sample_id")) %>%
  mutate(
    shape = ifelse(nen_subtype == "PD", 22, 21),
    fill = ifelse(split == "test", "black", "white")
  )

site_levels <- sort(unique(plot_df$primary_site_inferred))
site_colors <- setNames(grDevices::hcl.colors(length(site_levels), "Set2"), site_levels)

p_pca <- ggplot(plot_df, aes(x = PC1, y = PC2, color = primary_site_inferred)) +
  geom_point(aes(shape = shape, fill = fill), size = 3, stroke = 0.8) +
  scale_shape_identity() +
  scale_fill_identity() +
  scale_color_manual(values = site_colors) +
  labs(color = "Primary site (inferred)", x = "PC1", y = "PC2") +
  theme_classic(base_size = 12)

ggsave(file.path(out_dir, "fig5a_pca_primary_site.png"), p_pca, width = 6.0, height = 4.5, dpi = 300)
ggsave(file.path(out_dir, "fig5a_pca_primary_site.pdf"), p_pca, width = 6.0, height = 4.5)

# ---- Correlation heatmap ----
cor_mat <- cor(t(expr_mat), method = "pearson")
write_tsv(as.data.frame(cor_mat) %>% mutate(sample_id = rownames(cor_mat)) %>% relocate(sample_id),
          file.path(table_dir, "fig5c_pairwise_correlation.tsv"))

heatmap_path <- file.path(out_dir, "fig5c_pairwise_correlation_heatmap.png")
if (requireNamespace("pheatmap", quietly = TRUE)) {
  pheatmap::pheatmap(
    cor_mat,
    filename = heatmap_path,
    width = 6.5, height = 5.5
  )
} else {
  png(heatmap_path, width = 1200, height = 1000, res = 150)
  heatmap(cor_mat, symm = TRUE, Colv = NA, Rowv = NA, col = colorRampPalette(c("blue", "white", "red"))(100))
  dev.off()
}

# ---- QQ plots ----
qq_plot <- function(a, b, out_prefix) {
  xa <- expr_mat[a, ]
  xb <- expr_mat[b, ]
  df <- data.frame(
    q1 = sort(xa),
    q2 = sort(xb)
  )
  p <- ggplot(df, aes(x = q1, y = q2)) +
    geom_point(size = 0.7, alpha = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    theme_classic(base_size = 12) +
    labs(x = paste0(a, " expression rank"), y = paste0(b, " expression rank"))
  ggsave(file.path(out_dir, paste0(out_prefix, ".png")), p, width = 4.0, height = 4.0, dpi = 300)
  ggsave(file.path(out_dir, paste0(out_prefix, ".pdf")), p, width = 4.0, height = 4.0)
}

# Most/least similar pairs (excluding self)
cor_mat_off <- cor_mat
diag(cor_mat_off) <- NA
max_pair <- which(cor_mat_off == max(cor_mat_off, na.rm = TRUE), arr.ind = TRUE)[1, ]
min_pair <- which(cor_mat_off == min(cor_mat_off, na.rm = TRUE), arr.ind = TRUE)[1, ]

max_a <- rownames(cor_mat_off)[max_pair[1]]
max_b <- colnames(cor_mat_off)[max_pair[2]]
min_a <- rownames(cor_mat_off)[min_pair[1]]
min_b <- colnames(cor_mat_off)[min_pair[2]]

# Anchor sample (if present) for panel B
anchor_candidates <- c("SC14-4289", "SC14-4289_HG", "SC14-4289_LG")
anchor <- anchor_candidates[anchor_candidates %in% rownames(expr_mat)][1]
if (is.na(anchor)) anchor <- rownames(expr_mat)[1]
anchor_cor <- cor_mat_off[anchor, ]
anchor_cor[is.na(anchor_cor)] <- -Inf
nearest <- names(which.max(anchor_cor))

qq_plot(anchor, nearest, "fig5b_qq_anchor_neighbor")
qq_plot(min_a, min_b, "fig5d_qq_lowest_pair")
qq_plot(max_a, max_b, "fig5e_qq_highest_pair")

# ---- Save inferred sites table ----
write_tsv(meta %>% select(expr_sample_id, sample_id, primary_site, primary_site_inferred, nen_subtype, split),
          file.path(table_dir, "fig5a_primary_site_inference.tsv"))

cat("Wrote Fig5 A/B/C/D/E outputs to:", out_dir, "\n")
