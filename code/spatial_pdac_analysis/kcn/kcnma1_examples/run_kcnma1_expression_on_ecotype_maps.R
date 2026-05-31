#!/usr/bin/env Rscript

###############################################################################
# KCNMA1-positive spots on top of ecotype maps
#
# This small reference script overlays KCNMA1-positive spots on the ecotype
# maps. It was kept because it gives a quick visual check of where KCNMA1 falls
# without removing the original ecotype colors.
###############################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
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

save_multipage_pdf <- function(plot_list, filename, width = 4.6, height = 3.8) {
  pdf_path <- file.path(fig_dir, filename)
  pdf(pdf_path, width = width, height = height, onefile = TRUE)
  on.exit(dev.off(), add = TRUE)
  invisible(lapply(plot_list, print))
}

write_tsv_simple <- function(x, filename) {
  write.table(
    x,
    file = file.path(tab_dir, filename),
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

build_manual_color_legend <- function(color_mapping_df, title_text = "cc_ischia_10") {
  legend_df <- color_mapping_df %>%
    mutate(y = rev(seq_len(n())))

  ggplot(legend_df, aes(x = 1, y = y)) +
    geom_tile(aes(fill = label), width = 0.35, height = 0.35) +
    geom_text(aes(x = 1.35, label = label), hjust = 0, size = 3.6) +
    scale_fill_manual(values = setNames(legend_df$color, legend_df$label)) +
    guides(fill = "none") +
    coord_cartesian(xlim = c(0.7, 4.7), ylim = c(0.5, nrow(legend_df) + 0.8), clip = "off") +
    labs(title = title_text) +
    theme_void(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, size = 10),
      plot.margin = margin(10, 10, 10, 10)
    )
}

build_spatial_overlay_plot <- function(seurat_obj, image_name, metadata_df, expr_df, palette_values, point_size = 0.65) {
  image_obj <- seurat_obj@images[[image_name]]
  image_array <- image_obj@image
  image_height <- dim(image_array)[1]
  scale_factor <- image_obj@scale.factors$hires

  coords_df <- image_obj@coordinates %>%
    rownames_to_column("spot") %>%
    left_join(metadata_df %>% select(spot, cc_ischia_10), by = "spot") %>%
    left_join(expr_df %>% select(spot, kcnma1_detected), by = "spot") %>%
    mutate(
      x_plot = imagecol * scale_factor,
      y_plot = image_height - (imagerow * scale_factor),
      cc_ischia_10 = factor(as.character(cc_ischia_10), levels = names(palette_values)),
      kcnma1_detected = isTRUE(kcnma1_detected) | (!is.na(kcnma1_detected) & kcnma1_detected)
    )

  ggplot(coords_df, aes(x = x_plot, y = y_plot)) +
    geom_point(aes(color = cc_ischia_10), size = point_size, alpha = 0.95) +
    geom_point(
      data = coords_df %>% filter(kcnma1_detected),
      color = "black",
      size = point_size + 0.08,
      alpha = 0.70
    ) +
    scale_color_manual(values = palette_values, na.value = "grey80") +
    coord_fixed(expand = FALSE) +
    labs(title = image_name, color = "cc_ischia_10") +
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
paper_path <- file.path(input_dir, "spatial_transcriptomics_PDAC_paper.pdf")

if (!file.exists(object_path)) stop("Missing input object: ", object_path)
if (!file.exists(paper_path)) stop("Missing local paper PDF: ", paper_path)

message("Loading ST Seurat object for KCNMA1 ecotype overlay...")
st_obj <- readRDS(object_path)
st_obj <- suppressWarnings(UpdateSeuratObject(st_obj))

metadata <- st_obj@meta.data %>%
  rownames_to_column("spot")

selected_images <- metadata %>%
  filter(Origin %in% c("Normal Pancreas", "Pancreas")) %>%
  distinct(orig.ident, Origin) %>%
  arrange(factor(Origin, levels = c("Normal Pancreas", "Pancreas")), orig.ident)

spatial_mat <- get_assay_matrix_compat(st_obj, assay = "Spatial", layer_name = "data")
if (!("KCNMA1" %in% rownames(spatial_mat))) {
  stop("KCNMA1 is not available in the Spatial assay.")
}

kcnma1_df <- data.frame(
  spot = colnames(st_obj),
  expression = as.numeric(spatial_mat["KCNMA1", ]),
  stringsAsFactors = FALSE
) %>%
  mutate(kcnma1_detected = expression > 0) %>%
  left_join(metadata %>% select(spot, orig.ident, Origin, cc_ischia_10), by = "spot")

overlay_summary <- kcnma1_df %>%
  filter(orig.ident %in% selected_images$orig.ident) %>%
  group_by(orig.ident, Origin) %>%
  summarise(
    n_spots = n(),
    n_kcnma1_positive = sum(kcnma1_detected, na.rm = TRUE),
    pct_kcnma1_positive = round(100 * n_kcnma1_positive / n_spots, 2),
    .groups = "drop"
  )
write_tsv_simple(overlay_summary, "kcnma1_ecotype_overlay_summary_by_section.tsv")

figure13_color_mapping <- data.frame(
  color = c(
    "#3ab44b", "#bcf60b", "#f031e6", "#901eb4", "#fabdbe",
    "#e7872b", "#e6194b", "#317ec2", "#ffe119", "#11dfed"
  ),
  label = paste0("CC", seq_len(10)),
  stringsAsFactors = FALSE
)
write_tsv_simple(figure13_color_mapping, "ecotype_color_mapping.tsv")

ecotype_spatial_colors <- setNames(figure13_color_mapping$color, figure13_color_mapping$label)

plot_list <- lapply(selected_images$orig.ident, function(img_name) {
  build_spatial_overlay_plot(
    seurat_obj = st_obj,
    image_name = img_name,
    metadata_df = metadata,
    expr_df = kcnma1_df,
    palette_values = ecotype_spatial_colors,
    point_size = 0.65
  ) + build_manual_color_legend(
    color_mapping_df = figure13_color_mapping,
    title_text = "cc_ischia_10"
  ) + plot_layout(widths = c(4.2, 1.8))
})

save_multipage_pdf(
  plot_list = plot_list,
  filename = "kcnma1_expression_on_ecotype_maps.pdf",
  width = 4.4,
  height = 3.6
)


message("ST PDAC KCNMA1 ecotype overlay complete.")


