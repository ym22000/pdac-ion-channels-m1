#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

args <- commandArgs(trailingOnly = FALSE)
file_args <- args[grep("^--file=", args)]
script_path <- if (length(file_args) > 0) sub("^--file=", "", file_args[[1]]) else getwd()
script_dir <- normalizePath(dirname(script_path))
tables_dir <- file.path(script_dir, "tables")
fig_dir <- file.path(script_dir, "figures", "by_kcn")
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
  "Hepatocytes" = "#999999"
)

view_labels <- c(
  intra = "Intra (same spot)",
  juxta_ct = "Juxta (neighbor spots)",
  para_ct = "Para (broader context)"
)

all_celltypes <- c(
  "Fibroblasts",
  "Tumor Epithelial cells",
  "TAM",
  "Normal Epithelial cells",
  "Endothelial cells",
  "Monocytes",
  "PVL",
  "T-NK cells",
  "B cells",
  "Hepatocytes"
)
all_celltypes <- all_celltypes[all_celltypes %in% unique(edges$Predictor_original)]

draw_panel <- function(kcn_edges, view_name, dominant_view) {
  view_edges <- kcn_edges %>%
    filter(view == view_name) %>%
    arrange(importance_rank, desc(normalized_importance))

  if (nrow(view_edges) == 0) {
    plot.new()
    title(main = view_labels[[view_name]], font.main = 2, cex.main = 1.1)
    text(0.5, 0.5, "No edges available", cex = 1)
    return(invisible(NULL))
  }

  panel_celltypes <- unique(view_edges$Predictor_original)
  panel_celltypes <- all_celltypes[all_celltypes %in% panel_celltypes]
  max_y <- max(length(panel_celltypes), 3)

  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0.5, max_y + 1.5))

  panel_title <- view_labels[[view_name]]
  if (identical(view_name, dominant_view)) {
    panel_title <- paste0(panel_title, "  [dominant]")
  }
  title(main = panel_title, font.main = 2, cex.main = 1.1)

  kcn_y <- rev(seq_len(nrow(view_edges)))
  ct_y <- rev(match(panel_celltypes, panel_celltypes))
  names(ct_y) <- panel_celltypes

  view_edges <- view_edges %>%
    mutate(
      x = 0.15,
      xend = 0.85,
      y = kcn_y,
      yend = unname(ct_y[Predictor_original])
    )

  apply(view_edges, 1, function(row) {
    color <- celltype_colors[[row[["Predictor_original"]]]]
    segments(
      x0 = as.numeric(row[["x"]]),
      y0 = as.numeric(row[["y"]]),
      x1 = as.numeric(row[["xend"]]),
      y1 = as.numeric(row[["yend"]]),
      col = grDevices::adjustcolor(color, alpha.f = 0.22 + 0.70 * as.numeric(row[["normalized_importance"]])),
      lwd = 0.8 + 4.0 * as.numeric(row[["normalized_importance"]]),
      lend = "round"
    )
  })

  points(rep(0.15, nrow(view_edges)), kcn_y, pch = 16, cex = 1.0, col = "#4d4d4d")
  text(rep(0.11, nrow(view_edges)), kcn_y, labels = rep(unique(view_edges$target_gene), nrow(view_edges)), adj = 1, cex = 0.85)

  points(rep(0.85, length(panel_celltypes)), ct_y, pch = 16, cex = 1.2, col = celltype_colors[panel_celltypes])
  text(rep(0.89, length(panel_celltypes)), ct_y, labels = panel_celltypes, adj = 0, cex = 0.82, col = celltype_colors[panel_celltypes])

  text(0.15, max_y + 1.0, "KCN", font = 2, cex = 0.85)
  text(0.85, max_y + 1.0, "Cell types", font = 2, cex = 0.85)

  edge_labels <- paste0(
    view_edges$Predictor_original,
    " | norm=",
    sprintf("%.2f", as.numeric(view_edges$normalized_importance))
  )
  text(rep(0.50, nrow(view_edges)), kcn_y, labels = edge_labels, cex = 0.72, col = "#333333")
}

kcn_list <- sort(unique(edges$target_gene))

for (kcn in kcn_list) {
  kcn_edges <- edges %>% filter(target_gene == kcn)
  kcn_summary <- summary_df %>% filter(target_gene == kcn)
  dominant_view <- if (nrow(kcn_summary) > 0) kcn_summary$dominant_view[[1]] else NA_character_
  dominant_predictor <- if (nrow(kcn_summary) > 0) kcn_summary$dominant_predictor[[1]] else "NA"
  abs_fraction <- if (nrow(kcn_summary) > 0) kcn_summary$dominant_view_abs_fraction[[1]] else NA_real_

  out_path <- file.path(fig_dir, paste0(tolower(kcn), "_misty_bipartite_network.pdf"))
  pdf(out_path, width = 10, height = 10)
  par(mfrow = c(3, 1), mar = c(1.8, 1.2, 3.2, 1.2))

  for (view_name in c("intra", "juxta_ct", "para_ct")) {
    draw_panel(kcn_edges, view_name, dominant_view)
  }

  mtext(
    paste0(
      kcn,
      " | dominant view: ",
      dominant_view,
      " | dominant predictor: ",
      dominant_predictor,
      " | abs fraction: ",
      ifelse(is.na(abs_fraction), "NA", sprintf("%.2f", abs_fraction))
    ),
    side = 3,
    outer = TRUE,
    line = -2,
    cex = 1.0,
    font = 2
  )

  dev.off()
}

message("Saved per-KCN bipartite PDFs to: ", fig_dir)
