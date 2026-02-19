#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(readxl)
  library(ComplexHeatmap)
  library(circlize)
})

# Paths
input_expr <- "../../data/NECPath.all.data.gene.h.tsv"
input_meta <- "../../supplementary/tables/Table S1-2.xlsx"
out_dir <- "../../figures/individuals/fig2/"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Load expression (genes x samples)
raw_counts <- read_tsv(input_expr, show_col_types = FALSE)
counts_mat <- raw_counts %>%
  column_to_rownames("Name") %>%
  mutate(across(everything(), ~ as.numeric(str_trim(.)))) %>%
  as.matrix()

stopifnot(all(is.finite(counts_mat)))

# Log2 transform
expr_log <- log2(counts_mat + 1)

# Z-score by gene
expr_z <- t(scale(t(expr_log)))
expr_z[is.na(expr_z)] <- 0

# Load metadata (first sheet) with robust header detection
meta_raw <- read_excel(input_meta, sheet = 1, col_names = FALSE)
meta_raw <- as.data.frame(meta_raw)

# Find a header row (look for accession_number or sample_id)
header_row <- which(apply(meta_raw, 1, function(r) any(str_detect(tolower(as.character(r)), "accession_number|sample_id|sample"))))[1]
if (is.na(header_row)) {
  stop("Could not detect header row in metadata (expected 'accession_number' or 'sample_id').")
}

header <- as.character(unlist(meta_raw[header_row, , drop = TRUE]))
meta <- meta_raw[(header_row + 1):nrow(meta_raw), , drop = FALSE]
names(meta) <- header
meta <- meta %>% filter(if_any(everything(), ~ !is.na(.x)))

# Try to detect sample ID column
sample_col <- names(meta)[str_detect(tolower(names(meta)), "sample|id|accession")][1]
if (is.na(sample_col)) {
  stop("Could not detect sample ID column in metadata. Please rename a column to include 'sample', 'id', or 'accession'.")
}

# Align samples (allow partial matching by stripping _LG/_HG suffix)
meta <- meta %>%
  mutate(sample_id = as.character(.data[[sample_col]]))

expr_ids <- colnames(expr_z)
expr_key <- str_replace(expr_ids, "(_LG|_HG)$", "")
meta_key <- meta$sample_id

idx <- match(expr_key, meta_key)
if (all(is.na(idx))) {
  stop("No matching samples between expression matrix and metadata (even after stripping _LG/_HG).")
}

keep <- which(!is.na(idx))
expr_z <- expr_z[, keep, drop = FALSE]
meta <- meta[idx[keep], , drop = FALSE]
meta$sample_id <- meta_key[idx[keep]]

# Compute sample clustering and 3 groups
sample_cor <- cor(expr_z, method = "pearson")
sample_dist <- as.dist(1 - sample_cor)
hc <- hclust(sample_dist, method = "complete")
clusters <- cutree(hc, k = 3)
cluster_labels <- paste0("group", clusters)

# Reorder by clustering
ord <- hc$order
expr_z <- expr_z[, ord, drop = FALSE]
meta <- meta[ord, , drop = FALSE]
cluster_labels <- cluster_labels[ord]

# Select top variable genes (adjust if needed)
top_n_genes <- 50
var_genes <- order(apply(expr_z, 1, var), decreasing = TRUE)
expr_top <- expr_z[var_genes[1:min(top_n_genes, length(var_genes))], , drop = FALSE]

# Helper to get column if exists
get_col <- function(df, patterns) {
  nm <- names(df)
  idx <- which(str_detect(tolower(nm), paste(patterns, collapse = "|")))
  if (length(idx) == 0) return(NULL)
  nm[idx[1]]
}

# Extract annotations if present
col_subtype <- get_col(meta, c("subtype", "nen", "wd", "pd"))
col_gender  <- get_col(meta, c("gender", "sex"))
col_race    <- get_col(meta, c("race"))
col_grade   <- get_col(meta, c("grade", "g1", "g2", "g3"))
col_stage   <- get_col(meta, c("stage"))
col_vital   <- get_col(meta, c("vital", "status", "alive", "dead"))
col_site    <- get_col(meta, c("primary", "site", "origin"))
col_age     <- get_col(meta, c("age"))

# Build top annotation list
ann_list <- list(
  Group = cluster_labels
)

if (!is.null(col_subtype)) ann_list$Subtype <- meta[[col_subtype]]
if (!is.null(col_gender))  ann_list$Gender  <- meta[[col_gender]]
if (!is.null(col_race))    ann_list$Race    <- str_trim(tolower(as.character(meta[[col_race]])))
if (!is.null(col_grade))   ann_list$Grade   <- meta[[col_grade]]
if (!is.null(col_stage)) {
  stage_raw <- str_trim(toupper(as.character(meta[[col_stage]])))
  stage_raw[stage_raw %in% c("NA", "", "N/A")] <- NA
  stage_raw <- str_replace_all(stage_raw, "\\s+", "")
  stage_raw <- str_replace(stage_raw, "^I-III$", "I-III")
  stage_raw <- str_replace(stage_raw, "^II-III$", "II-III")
  stage_raw <- str_replace(stage_raw, "^III-IV$", "III-IV")
  ann_list$Stage <- stage_raw
}
if (!is.null(col_vital))   ann_list$Vital   <- meta[[col_vital]]
if (!is.null(col_site))    ann_list$Primary_Site <- meta[[col_site]]
if (!is.null(col_age))     ann_list$Age     <- as.numeric(meta[[col_age]])

# Color maps
group_levels <- sort(unique(cluster_labels))
group_palette <- c("#ff7f0eff", "#2C7BB6", "#7f7f7f")
ann_colors <- list(
  Group = setNames(group_palette[seq_along(group_levels)], group_levels)
)

if (!is.null(col_subtype)) ann_colors$Subtype <- c("PD"="#ff7f0eff","WD"="#2C7BB6")
if (!is.null(col_gender)) ann_colors$Gender <- c("F"="#e377c2","M"="#1f77b4")
if (!is.null(col_vital)) ann_colors$Vital <- c("Alive"="#4daf4a","Dead"="#e41a1c")
if (!is.null(col_grade)) ann_colors$Grade <- c("G1"="#66c2a5","G2"="#fc8d62","G3"="#8da0cb")
if (!is.null(col_stage)) {
  stage_levels <- sort(unique(ann_list$Stage))
  stage_palette <- c(
    "I"="#a6cee3","II"="#1f78b4","III"="#b2df8a","IV"="#33a02c",
    "I-III"="#6a3d9a","II-III"="#cab2d6","III-IV"="#fb9a99"
  )
  ann_colors$Stage <- stage_palette[stage_levels]
}
if (!is.null(col_race)) {
  race_levels <- sort(unique(ann_list$Race))
  base_race_palette <- c("asian"="#1b9e77","black"="#d95f02","white"="#7570b3")
  if (all(race_levels %in% names(base_race_palette))) {
    ann_colors$Race <- base_race_palette[race_levels]
  } else {
    ann_colors$Race <- setNames(
      colorRampPalette(c("#1b9e77", "#d95f02", "#7570b3", "#66a61e"))(length(race_levels)),
      race_levels
    )
  }
}

# Primary site palette (generate automatically)
if (!is.null(col_site)) {
  sites <- sort(unique(as.character(meta[[col_site]])))
  ann_colors$Primary_Site <- setNames(colorRampPalette(c("#f1c40f","#e67e22","#9b59b6","#3498db","#2ecc71"))(length(sites)), sites)
}

# Age color ramp
if (!is.null(col_age)) {
  ann_colors$Age <- colorRamp2(range(ann_list$Age, na.rm = TRUE), c("#f7fbff", "#08306b"))
}

# Mutation tiles if columns exist
mut_genes <- c("IDH1","APC","TP53","KRAS","PIK3CA","CTNNB1","PTEN","FBXW7","BRAF","RB1")
mut_cols <- mut_genes[mut_genes %in% names(meta)]
mut_mat <- NULL
if (length(mut_cols) > 0) {
  mut_mat <- as.matrix(meta[, mut_cols, drop = FALSE])
  rownames(mut_mat) <- meta$sample_id
  mut_mat <- t(mut_mat)
}

# Top annotation
legend_params <- lapply(ann_list, function(x) {
  if (is.numeric(x)) {
    list(title_gp = grid::gpar(fontsize = 9), labels_gp = grid::gpar(fontsize = 8))
  } else {
    list(title_gp = grid::gpar(fontsize = 9), labels_gp = grid::gpar(fontsize = 8))
  }
})

top_ha <- HeatmapAnnotation(
  df = ann_list,
  col = ann_colors,
  annotation_name_side = "left",
  show_annotation_name = TRUE,
  annotation_legend_param = legend_params
)

# Bottom mutation annotation (optional)
bottom_ha <- NULL
if (!is.null(mut_mat)) {
  bottom_ha <- HeatmapAnnotation(
    Mutation = anno_heatmap(mut_mat, col = c("0" = "#f0f0f0", "1" = "#2ca25f")),
    annotation_name_side = "left",
    show_annotation_name = TRUE,
    annotation_legend_param = list(
      Mutation = list(title_gp = grid::gpar(fontsize = 9), labels_gp = grid::gpar(fontsize = 8))
    )
  )
}

# Heatmap
ht <- Heatmap(
  expr_top,
  name = "Z-score",
  top_annotation = top_ha,
  bottom_annotation = bottom_ha,
  col = colorRamp2(c(-2, 0, 2), c("#2c7bb6", "#f7f7f7", "#d7191c")),
  cluster_columns = hc,
  cluster_rows = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_title = "NET patients",
  column_names_rot = 90,
  row_names_gp = grid::gpar(fontsize = 6),
  heatmap_legend_param = list(
    at = c(-2, -1, 0, 1, 2),
    title = "Z-score",
    legend_height = grid::unit(4, "cm")
  )
)

# Output
pdf(file.path(out_dir, "fig2a_heatmap.pdf"), width = 11, height = 8)
draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legend = TRUE
)
dev.off()

png(file.path(out_dir, "fig2a_heatmap.png"), width = 2200, height = 1600, res = 200)
draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legend = TRUE
)
dev.off()

cat("Wrote fig2a_heatmap.pdf and fig2a_heatmap.png\n")
