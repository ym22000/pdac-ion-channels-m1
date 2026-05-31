#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tibble)
})

args <- commandArgs(trailingOnly = FALSE)
file_args <- args[grep("^--file=", args)]
script_path <- if (length(file_args) > 0) sub("^--file=", "", file_args[[1]]) else getwd()
script_dir <- normalizePath(dirname(script_path))
tables_dir <- file.path(script_dir, "tables")
fig_dir <- file.path(script_dir, "figures", "by_kcn_visual")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

edge_path <- file.path(tables_dir, "kcn_union_25_misty_network_edges_top3.tsv")
summary_path <- file.path(tables_dir, "kcn_union_25_misty_schematic_summary.tsv")

edges <- read.delim(edge_path, stringsAsFactors = FALSE)
summary_df <- read.delim(summary_path, stringsAsFactors = FALSE)

celltype_colors <- c(
  "Fibroblasts" = "#8c510a",
  "Tumor Epithelial cells" = "#b2182b",
  "TAM" = "#2166ac",
  "Normal Epithelial cells" = "#1b7837",
  "Endothelial cells" = "#762a83",
  "Monocytes" = "#4393c3",
  "PVL" = "#f4a582",
  "T-NK cells" = "#4d9221",
  "B cells" = "#92c5de",
  "Hepatocytes" = "#7f7f7f"
)

view_labels <- c(
  intra = "Intra\nsame spot",
  juxta_ct = "Juxta\nneighbor spots",
  para_ct = "Para\nbroader context"
)

make_kcn_plot <- function(kcn) {
  kcn_edges <- edges %>%
    filter(target_gene == kcn) %>%
    mutate(
      view_label = factor(view_labels[view], levels = unname(view_labels[c("intra", "juxta_ct", "para_ct")]))
    )

  if (nrow(kcn_edges) == 0) {
    return(NULL)
  }

  kcn_summary <- summary_df %>% filter(target_gene == kcn)
  dominant_view <- if (nrow(kcn_summary) > 0) kcn_summary$dominant_view[[1]] else NA_character_
  dominant_predictor <- if (nrow(kcn_summary) > 0) kcn_summary$dominant_predictor[[1]] else "NA"
  abs_fraction <- if (nrow(kcn_summary) > 0) kcn_summary$dominant_view_abs_fraction[[1]] else NA_real_

  panel_df <- kcn_edges %>%
    group_by(view_label) %>%
    arrange(importance_rank, desc(normalized_importance), .by_group = TRUE) %>%
    mutate(
      yend = c(3, 2, 1)[seq_len(n())],
      x = 0.18,
      y = 2,
      xend = 0.82,
      curve = c(0.22, 0, -0.22)[seq_len(n())],
      edge_label_x = 0.50,
      edge_label_y = (y + yend) / 2 + c(0.12, 0.00, -0.12)[seq_len(n())]
    ) %>%
    ungroup()

  kcn_node_df <- tibble(
    view_label = factor(unname(view_labels[c("intra", "juxta_ct", "para_ct")]), levels = unname(view_labels[c("intra", "juxta_ct", "para_ct")])),
    x = 0.18,
    y = 2,
    label = kcn
  )

  ct_node_df <- panel_df %>%
    transmute(
      view_label,
      x = xend,
      y = yend,
      label = Predictor_original,
      color = Predictor_original
    )

  dom_text <- if (is.na(abs_fraction)) {
    paste0("Dominant view: ", dominant_view, " | Dominant predictor: ", dominant_predictor)
  } else {
    paste0(
      "Dominant view: ", dominant_view,
      " | Dominant predictor: ", dominant_predictor,
      " | Abs fraction: ", sprintf("%.2f", abs_fraction)
    )
  }

  ggplot() +
    facet_wrap(~view_label, nrow = 1) +
    geom_curve(
      data = panel_df,
      aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend,
        color = Predictor_original,
        size = normalized_importance,
        alpha = normalized_importance,
        curvature = curve
      ),
      lineend = "round",
      show.legend = FALSE
    ) +
    geom_label(
      data = kcn_node_df,
      aes(x = x, y = y, label = label),
      fill = "#333333",
      color = "white",
      label.size = 0,
      fontface = "bold",
      size = 4.8,
      label.padding = unit(0.22, "lines")
    ) +
    geom_label(
      data = ct_node_df,
      aes(x = x, y = y, label = label, fill = color),
      color = "white",
      label.size = 0,
      fontface = "bold",
      size = 4.1,
      label.padding = unit(0.18, "lines"),
      show.legend = FALSE
    ) +
    geom_text(
      data = panel_df,
      aes(
        x = edge_label_x,
        y = edge_label_y,
        label = sprintf("norm %.2f", normalized_importance)
      ),
      size = 3.3,
      color = "#444444",
      fontface = "bold"
    ) +
    scale_fill_manual(values = celltype_colors) +
    scale_color_manual(values = celltype_colors) +
    scale_size(range = c(1.2, 4.5)) +
    scale_alpha(range = c(0.45, 0.95)) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0.5, 3.5), clip = "off") +
    labs(
      title = paste0(kcn, " spatial cell-type context from MISTy"),
      subtitle = dom_text
    ) +
    theme_void(base_size = 12) +
    theme(
      strip.text = element_text(face = "bold", size = 13, color = "#222222"),
      plot.title = element_text(face = "bold", size = 18, color = "#111111"),
      plot.subtitle = element_text(size = 11, color = "#444444"),
      plot.background = element_rect(fill = "#fcfcfa", color = NA),
      panel.background = element_rect(fill = "#fcfcfa", color = "#efefea", linewidth = 0.4),
      panel.spacing = unit(1.0, "lines"),
      plot.margin = margin(16, 18, 12, 18)
    )
}

kcn_list <- sort(unique(edges$target_gene))

for (kcn in kcn_list) {
  plot_obj <- make_kcn_plot(kcn)
  if (is.null(plot_obj)) {
    next
  }
  out_path <- file.path(fig_dir, paste0(tolower(kcn), "_misty_bipartite_network_visual.pdf"))
  ggsave(out_path, plot_obj, width = 14, height = 6.2, units = "in", device = cairo_pdf)
}

message("Saved per-KCN visual PDFs to: ", fig_dir)
