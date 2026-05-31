#!/usr/bin/env Rscript

###############################################################################
# Spatial marker projections for a few reference genes
#
# This script projects a small set of markers used as spatial references in the
# project, spanning epithelial, stromal, and iCAF-like territories.
###############################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  file_path <- sub(file_arg, "", args[grep(file_arg, args)])
  if (length(file_path) == 0) {
    return(normalizePath(getwd()))
  }
  normalizePath(dirname(file_path))
}

save_multipage_pdf <- function(plot_list, filename, out_dir, width = 3.4, height = 3.4) {
  pdf_path <- file.path(out_dir, filename)
  pdf(pdf_path, width = width, height = height, onefile = TRUE)
  on.exit(dev.off(), add = TRUE)
  invisible(lapply(plot_list, print))
}

write_tsv_simple <- function(x, filename, out_dir) {
  write.table(
    x,
    file = file.path(out_dir, filename),
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE,
    na = "NA"
  )
}

get_assay_matrix_compat <- function(seurat_obj, assay = "Spatial", layer_name = "data") {
  tryCatch(
    Seurat::GetAssayData(seurat_obj, assay = assay, layer = layer_name),
    error = function(e) Seurat::GetAssayData(seurat_obj, assay = assay, slot = layer_name)
  )
}

build_spatial_feature_plot <- function(seurat_obj, image_name, expr_df, feature_name, point_size = 0.65) {
  image_obj <- seurat_obj@images[[image_name]]
  image_array <- image_obj@image
  image_height <- dim(image_array)[1]
  scale_factor <- image_obj@scale.factors$hires

  coords_df <- image_obj@coordinates %>%
    rownames_to_column("spot") %>%
    left_join(expr_df, by = "spot") %>%
    mutate(
      x_plot = imagecol * scale_factor,
      y_plot = image_height - (imagerow * scale_factor)
    )

  ggplot(coords_df, aes(x = x_plot, y = y_plot)) +
    geom_point(aes(color = expression), size = point_size, alpha = 0.75) +
    scale_color_gradient(low = "grey95", high = "#8B0000", na.value = "grey80") +
    coord_fixed(expand = FALSE) +
    labs(title = image_name, color = feature_name) +
    theme_void(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
      legend.position = "none"
    )
}

script_dir <- get_script_dir()

find_project_dir <- function(start_dir) {
  current_dir <- normalizePath(start_dir)
  repeat {
    if (file.exists(file.path(current_dir, 'inputs', 'PDAC_Updated_ST.rds'))) {
      return(current_dir)
    }
    parent_dir <- dirname(current_dir)
    if (identical(parent_dir, current_dir)) {
      stop('Could not find the ST project directory.')
    }
    current_dir <- parent_dir
  }
}

project_dir <- find_project_dir(script_dir)
input_dir <- file.path(project_dir, "inputs")
fig_dir <- file.path(script_dir, "figures")
tab_dir <- file.path(script_dir, "tables")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

object_path <- file.path(input_dir, "PDAC_Updated_ST.rds")
if (!file.exists(object_path)) stop("Missing input object: ", object_path)

message("Loading ST Seurat object for marker projections...")
st_obj <- readRDS(object_path)
st_obj <- suppressWarnings(UpdateSeuratObject(st_obj))

metadata <- st_obj@meta.data %>%
  rownames_to_column("spot")

selected_images <- metadata %>%
  filter(Origin %in% c("Normal Pancreas", "Pancreas")) %>%
  distinct(orig.ident, Origin) %>%
  arrange(factor(Origin, levels = c("Normal Pancreas", "Pancreas")), orig.ident)

write_tsv_simple(selected_images, "selected_normal_and_primary_sections.tsv", tab_dir)

spatial_assay_mat <- get_assay_matrix_compat(st_obj, assay = "Spatial", layer_name = "data")

marker_map <- data.frame(
  display_name = c(
    "COL1A1", "FAP", "EPCAM", "CK20",
    "C7", "C1R", "CCL19", "CCL21", "MFAP4", "CLU", "DCN",
    "CCN2", "ACTA2", "THBS1", "COL12A1"
  ),
  gene_symbol = c(
    "COL1A1", "FAP", "EPCAM", "KRT20",
    "C7", "C1R", "CCL19", "CCL21", "MFAP4", "CLU", "DCN",
    "CCN2", "ACTA2", "THBS1", "COL12A1"
  ),
  file_stem = c(
    "col1a1", "fap", "epcam", "krt20",
    "c7", "c1r", "ccl19", "ccl21", "mfap4", "clu", "dcn",
    "ccn2", "acta2", "thbs1", "col12a1"
  ),
  stringsAsFactors = FALSE
)
marker_map$in_spatial_assay <- marker_map$gene_symbol %in% rownames(spatial_assay_mat)
write_tsv_simple(marker_map, "marker_gene_availability.tsv", tab_dir)

for (i in seq_len(nrow(marker_map))) {
  display_name <- marker_map$display_name[i]
  gene_symbol <- marker_map$gene_symbol[i]
  file_stem <- marker_map$file_stem[i]

  expr_df <- data.frame(
    spot = colnames(st_obj),
    expression = if (gene_symbol %in% rownames(spatial_assay_mat)) {
      as.numeric(spatial_assay_mat[gene_symbol, ])
    } else {
      rep(0, ncol(st_obj))
    },
    stringsAsFactors = FALSE
  )

  plot_list <- lapply(selected_images$orig.ident, function(img_name) {
    build_spatial_feature_plot(
      seurat_obj = st_obj,
      image_name = img_name,
      expr_df = expr_df,
      feature_name = display_name,
      point_size = 0.65
    )
  })

  save_multipage_pdf(
    plot_list = plot_list,
    filename = paste0(file_stem, "_spatial_expression_normal_and_primary.pdf"),
    out_dir = fig_dir,
    width = 3.4,
    height = 3.4
  )
}


message("ST PDAC marker projections complete.")


