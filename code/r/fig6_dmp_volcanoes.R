#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(tibble)
  library(ggrepel)
})

# Paths
base_dir <- "/mnt/f/research_drive/Manuscripts/src/01.NEN/analysis/methyl/01.NEN/results/tables/dmp"
out_fig <- "/mnt/f/research_drive/Manuscripts/submissions/genome_biology/supplementary/figures"
out_tab <- "/mnt/f/research_drive/Manuscripts/submissions/genome_biology/supplementary/tables"

dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(out_tab, recursive = TRUE, showWarnings = FALSE)

contrasts <- tribble(
  ~label, ~file_full, ~file_annot, ~supp_fig,
  "NEC vs NET", "dmp_NECvsNET_full.tsv.gz", "dmp_NECvsNET_annot.tsv.gz", "Supplementary Figure 4",
  "NET vs Normal", "dmp_NETvsNormal_full.tsv.gz", "dmp_NETvsNormal_annot.tsv.gz", "Supplementary Figure 5",
  "NEC vs Normal", "dmp_NECvsNormal_full.tsv.gz", "dmp_NECvsNormal_annot.tsv.gz", "Supplementary Figure 6"
)

sig_fdr <- 0.05
sig_db  <- 0.2

summaries <- list()

for (i in seq_len(nrow(contrasts))) {
  c <- contrasts[i, ]
  full_path <- file.path(base_dir, c$file_full)
  annot_path <- file.path(base_dir, c$file_annot)

  full <- read_tsv(full_path, show_col_types = FALSE)
  annot <- read_tsv(annot_path, show_col_types = FALSE)

  full <- full %>%
    mutate(
      delta_beta = ifelse(!is.na(delta_beta), delta_beta, logFC),
      neglog10p = -log10(pmax(adj.P.Val, .Machine$double.eps)),
      sig = adj.P.Val < sig_fdr & abs(delta_beta) >= sig_db,
      direction = case_when(
        sig & delta_beta > 0 ~ "Hypermethylated",
        sig & delta_beta < 0 ~ "Hypomethylated",
        TRUE ~ "Not significant"
      )
    )

  # Top promoter-proximal DMPs for labeling (top 5 hypo + top 5 hyper by FDR, then |delta_beta|)
  labels_df <- annot %>%
    mutate(
      delta_beta = ifelse(!is.na(delta_beta), delta_beta, logFC),
      label = ifelse(
        !is.na(Gene_primary) & Gene_primary != "",
        paste0(CpG, " (", Gene_primary, ")"),
        CpG
      )
    ) %>%
    filter(
      adj.P.Val < sig_fdr,
      abs(delta_beta) >= sig_db,
      is_promoter_prox == TRUE
    ) %>%
    mutate(direction = ifelse(delta_beta > 0, "hyper", "hypo")) %>%
    group_by(direction) %>%
    arrange(adj.P.Val, desc(abs(delta_beta)), .by_group = TRUE) %>%
    slice_head(n = 5) %>%
    ungroup() %>%
    select(CpG, delta_beta, adj.P.Val, label)

  # Align labels to volcano axes
  if (nrow(labels_df) > 0) {
    labels_df <- labels_df %>%
      mutate(neglog10p = -log10(pmax(adj.P.Val, .Machine$double.eps)))
  }

  # Volcano plot
  p <- ggplot(full, aes(x = delta_beta, y = neglog10p)) +
    geom_point(aes(color = direction), size = 0.8, alpha = 0.8) +
    geom_vline(xintercept = c(-sig_db, sig_db), linetype = "dashed", linewidth = 0.4, color = "black") +
    geom_hline(yintercept = -log10(sig_fdr), linetype = "dashed", linewidth = 0.4, color = "black") +
    scale_color_manual(values = c(
      "Hypermethylated" = "#d7191c",
      "Hypomethylated" = "#2c7bb6",
      "Not significant" = "#bdbdbd"
    )) +
    labs(
      title = paste0(c$label, " DMPs"),
      x = expression(Delta*beta),
      y = expression(-log[10]~FDR),
      color = NULL
    ) +
    theme_classic(base_size = 12)

  if (nrow(labels_df) > 0) {
    p <- p +
      geom_text_repel(
        data = labels_df,
        aes(x = delta_beta, y = neglog10p, label = label),
        size = 3,
        max.overlaps = Inf,
        box.padding = 0.3,
        point.padding = 0.2,
        min.segment.length = 0
      )
  }

  ggsave(file.path(out_fig, paste0(c$supp_fig, ".png")), p, width = 6.5, height = 5, dpi = 300)
  ggsave(file.path(out_fig, paste0(c$supp_fig, ".pdf")), p, width = 6.5, height = 5)

  # Summary stats
  sig_tbl <- full %>% filter(sig)
  n_total <- nrow(full)
  n_sig <- nrow(sig_tbl)
  n_hyper <- sum(sig_tbl$delta_beta > 0, na.rm = TRUE)
  n_hypo <- sum(sig_tbl$delta_beta < 0, na.rm = TRUE)

  # promoter/island counts if available
  promoter_col <- "is_promoter_prox"
  island_col <- "Relation_to_Island"
  prom_count <- NA_integer_
  island_count <- NA_integer_
  if (promoter_col %in% names(annot)) {
    prom_count <- annot %>%
      mutate(delta_beta = ifelse(!is.na(delta_beta), delta_beta, logFC)) %>%
      filter(adj.P.Val < sig_fdr, abs(delta_beta) >= sig_db, is_promoter_prox == TRUE) %>%
      nrow()
  }
  if (island_col %in% names(annot)) {
    island_count <- annot %>%
      mutate(delta_beta = ifelse(!is.na(delta_beta), delta_beta, logFC)) %>%
      filter(adj.P.Val < sig_fdr, abs(delta_beta) >= sig_db, Relation_to_Island %in% c("Island", "N_Shore", "S_Shore")) %>%
      nrow()
  }

  summaries[[i]] <- tibble(
    contrast = c$label,
    total_tests = n_total,
    dmp_sig = n_sig,
    dmp_hyper = n_hyper,
    dmp_hypo = n_hypo,
    dmp_promoter_prox = prom_count,
    dmp_island_or_shore = island_count
  )
}

summary_df <- bind_rows(summaries)
write_tsv(summary_df, file.path(out_tab, "Supplementary_Data_DMP_summary.tsv"))
write.csv(summary_df, file.path(out_tab, "Supplementary_Data_DMP_summary.csv"), row.names = FALSE)

cat("Wrote DMP volcanoes to", out_fig, "and summary tables to", out_tab, "\n")
