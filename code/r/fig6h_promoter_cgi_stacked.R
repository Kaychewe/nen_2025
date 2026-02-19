#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: fig6h_promoter_cgi_stacked.R <stacked_barplot_dir> <out_dir>")
}

in_dir <- args[1]
out_dir <- args[2]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

prop_path <- file.path(in_dir, "A_promoter_CGI_two_contrasts_proportions.tsv")
fisher_path <- file.path(in_dir, "A_promoter_CGI_two_contrasts_fisher_between.tsv")

props <- read_tsv(prop_path, show_col_types = FALSE)
fisher <- read_tsv(fisher_path, show_col_types = FALSE)

# Ensure ordering: category with NET vs Normal then NEC vs NET
props <- props %>%
  mutate(
    contrast = factor(contrast, levels = c("NET vs Normal", "NEC vs NET")),
    A_category = factor(A_category, levels = c("CGI Promoter", "CGI", "Promoter")),
    x_combo = factor(
      x_combo,
      levels = c(
        "CGI Promoter \u27c2 NET vs Normal",
        "CGI Promoter \u27c2 NEC vs NET",
        "CGI \u27c2 NET vs Normal",
        "CGI \u27c2 NEC vs NET",
        "Promoter \u27c2 NET vs Normal",
        "Promoter \u27c2 NEC vs NET"
      )
    ),
    status = factor(status, levels = c("Hypo", "Hyper"))
  )

# Compute p-value label positions (centered over each category pair)
pair_centers <- props %>%
  distinct(A_category, x_combo) %>%
  mutate(x_index = as.numeric(x_combo)) %>%
  group_by(A_category) %>%
  summarise(x = mean(x_index), .groups = "drop") %>%
  left_join(fisher %>% select(category, p_between), by = c("A_category" = "category")) %>%
  mutate(
    p_label = paste0("p=", format(p_between, scientific = TRUE, digits = 2)),
    y = 1.05
  )

cols <- c("Hypo" = "#0000FF", "Hyper" = "#FF0000")

p <- ggplot(props, aes(x = x_combo, y = prop, fill = status)) +
  geom_col(width = 0.85, color = "black", linewidth = 0.4) +
  scale_fill_manual(values = cols, name = NULL) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.12)) +
  labs(
    x = NULL,
    y = "Proportion of sig CpGs (100%)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "top",
    legend.direction = "horizontal"
  ) +
  geom_text(
    data = pair_centers,
    aes(x = x, y = y, label = p_label),
    inherit.aes = FALSE,
    size = 3
  )

png_path <- file.path(out_dir, "fig6h_promoter_cgi_two_contrasts.png")
pdf_path <- file.path(out_dir, "fig6h_promoter_cgi_two_contrasts.pdf")

ggsave(png_path, p, width = 6.5, height = 3.8, dpi = 300)
ggsave(pdf_path, p, width = 6.5, height = 3.8)

cat("Wrote:\n", png_path, "\n", pdf_path, "\n")
