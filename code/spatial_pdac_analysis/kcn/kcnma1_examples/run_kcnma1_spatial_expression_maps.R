#!/usr/bin/env Rscript

###############################################################################
# KCNMA1 expression maps on the selected normal and primary sections
#
# This is one of the early KCNMA1-only plots kept for reference. It shows the
# same set of sections used in the main spatial figures with one shared color
# scale for KCNMA1 expression.
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

build_spatial_feature_plot <- function(
    seurat_obj,
    image_name,
    expr_df,
    feature_name,
    point_size = 0.65,
    limits = NULL,
    palette_values = NULL
) {
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
    geom_point(aes(color = expression_plot), size = point_size, alpha = 0.85) +
    scale_color_gradientn(
      colours = palette_values,
      limits = limits,
      oob = scales::squish,
      na.value = "grey80"
    ) +
    coord_fixed(expand = FALSE) +
    labs(title = image_name, color = feature_name) +
    theme_void(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 9),
      legend.text = element_text(size = 8)
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

message("Loading ST Seurat object for KCNMA1 projection...")
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
  left_join(metadata %>% select(spot, orig.ident, Origin, cc_ischia_10), by = "spot")

selected_expr <- kcnma1_df %>%
  filter(orig.ident %in% selected_images$orig.ident)

upper_limit <- as.numeric(stats::quantile(selected_expr$expression, probs = 0.99, na.rm = TRUE))
lower_limit <- 0

kcnma1_df <- kcnma1_df %>%
  mutate(expression_plot = pmin(pmax(expression, lower_limit), upper_limit))

kcnma1_sample_summary <- selected_expr %>%
  group_by(orig.ident, Origin) %>%
  summarise(
    n_spots = n(),
    n_detected = sum(expression > 0, na.rm = TRUE),
    pct_detected = round(100 * n_detected / n_spots, 2),
    mean_expression = mean(expression, na.rm = TRUE),
    median_expression = median(expression, na.rm = TRUE),
    q99_expression = stats::quantile(expression, probs = 0.99, na.rm = TRUE),
    .groups = "drop"
  )

kcnma1_ecotype_summary <- selected_expr %>%
  group_by(orig.ident, Origin, cc_ischia_10) %>%
  summarise(
    n_spots = n(),
    n_detected = sum(expression > 0, na.rm = TRUE),
    pct_detected = round(100 * n_detected / n_spots, 2),
    mean_expression = mean(expression, na.rm = TRUE),
    median_expression = median(expression, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(orig.ident, cc_ischia_10)

write_tsv_simple(selected_images, "selected_normal_and_primary_sections.tsv")
write_tsv_simple(kcnma1_sample_summary, "kcnma1_spatial_summary_by_section.tsv")
write_tsv_simple(kcnma1_ecotype_summary, "kcnma1_spatial_summary_by_section_and_ecotype.tsv")

roma_like_anchors <- c(
  "#f7f7f7",
  "#ddd9ec",
  "#c0b8db",
  "#9a8fc2",
  "#6d62a8",
  "#3f428d",
  "#17306f"
)
roma_like_palette <- grDevices::colorRampPalette(roma_like_anchors)(256)

plot_list <- lapply(selected_images$orig.ident, function(img_name) {
  build_spatial_feature_plot(
    seurat_obj = st_obj,
    image_name = img_name,
    expr_df = kcnma1_df,
    feature_name = "KCNMA1",
    point_size = 0.65,
    limits = c(lower_limit, upper_limit),
    palette_values = roma_like_palette
  )
})

save_multipage_pdf(
  plot_list = plot_list,
  filename = "kcnma1_spatial_expression_normal_and_primary.pdf",
  width = 4.4,
  height = 3.6
)


message("ST PDAC KCNMA1 projection complete.")


