#!/usr/bin/env Rscript

###############################################################################
# KCN expression maps on normal and primary sections
#
# This script makes one spatial PDF per KCN for the selected normal and primary
# sections. All spots are shown in grey and positive spots are recolored by
# expression intensity.
#
# The goal is just to keep a simple view of where each KCN is detected in the
# tissue before adding niche contours or downstream statistics.
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

save_multipage_pdf <- function(plot_list, filename, out_dir, width = 4.6, height = 3.8) {
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

build_spatial_expression_plot <- function(
    seurat_obj,
    image_name,
    gene_df,
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
    left_join(gene_df %>% select(spot, expression_plot, detected), by = "spot") %>%
    mutate(
      x_plot = imagecol * scale_factor,
      y_plot = image_height - (imagerow * scale_factor),
      detected = isTRUE(detected) | (!is.na(detected) & detected)
    )

  ggplot(coords_df, aes(x = x_plot, y = y_plot)) +
    geom_point(color = "grey80", size = point_size, alpha = 0.90) +
    geom_point(
      data = coords_df %>% filter(detected),
      aes(color = expression_plot),
      size = point_size,
      alpha = 0.95
    ) +
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
projection_fig_dir <- fig_dir

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(projection_fig_dir, recursive = TRUE, showWarnings = FALSE)

object_path <- file.path(input_dir, "PDAC_Updated_ST.rds")
paper_path <- file.path(input_dir, "spatial_transcriptomics_PDAC_paper_expression_map_normal_and_primary.pdf")
kcn_list_path <- normalizePath(file.path(project_dir, "..", "mycaf_icaf_correlations", "kcn_union_25.tsv"), mustWork = FALSE)

if (!file.exists(object_path)) stop("Missing input object: ", object_path)
if (!file.exists(paper_path)) stop("Missing local paper PDF: ", paper_path)
if (!file.exists(kcn_list_path)) stop("Missing KCN list: ", kcn_list_path)

message("Loading ST Seurat object for KCN25 expression projections...")
st_obj <- readRDS(object_path)
st_obj <- suppressWarnings(UpdateSeuratObject(st_obj))

metadata <- st_obj@meta.data %>%
  rownames_to_column("spot")

selected_images <- metadata %>%
  filter(Origin %in% c("Normal Pancreas", "Pancreas")) %>%
  distinct(orig.ident, Origin) %>%
  arrange(factor(Origin, levels = c("Normal Pancreas", "Pancreas")), orig.ident)

kcn_list <- read.delim(kcn_list_path, stringsAsFactors = FALSE) %>%
  pull(gene) %>%
  unique()

spatial_mat <- get_assay_matrix_compat(st_obj, assay = "Spatial", layer_name = "data")

gene_availability <- data.frame(
  gene = kcn_list,
  in_spatial_assay = kcn_list %in% rownames(spatial_mat),
  stringsAsFactors = FALSE
)

write_tsv_simple(gene_availability, "kcn_gene_availability.tsv", tab_dir)
write_tsv_simple(selected_images, "selected_normal_and_primary_sections.tsv", tab_dir)

sample_summaries <- list()
ecotype_summaries <- list()

expression_palette <- c("#ffe119", "#317ec2", "#e6194b")

for (gene_name in kcn_list) {
  if (!(gene_name %in% rownames(spatial_mat))) {
    next
  }

  gene_df <- data.frame(
    spot = colnames(st_obj),
    expression = as.numeric(spatial_mat[gene_name, ]),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      gene = gene_name,
      detected = expression > 0
    ) %>%
    left_join(metadata %>% select(spot, orig.ident, Origin, cc_ischia_10), by = "spot")

  gene_selected <- gene_df %>%
    filter(orig.ident %in% selected_images$orig.ident)

  upper_limit <- as.numeric(stats::quantile(gene_selected$expression, probs = 0.99, na.rm = TRUE))
  if (!is.finite(upper_limit) || upper_limit <= 0) {
    upper_limit <- max(gene_selected$expression, na.rm = TRUE)
  }
  if (!is.finite(upper_limit) || upper_limit <= 0) {
    upper_limit <- 1
  }

  gene_df <- gene_df %>%
    mutate(expression_plot = pmin(pmax(expression, 0), upper_limit))

  sample_summary <- gene_selected %>%
    group_by(gene, orig.ident, Origin) %>%
    summarise(
      n_spots = n(),
      n_detected = sum(detected, na.rm = TRUE),
      pct_detected = round(100 * n_detected / n_spots, 2),
      mean_expression = mean(expression, na.rm = TRUE),
      median_expression = median(expression, na.rm = TRUE),
      q99_expression = stats::quantile(expression, probs = 0.99, na.rm = TRUE),
      .groups = "drop"
    )

  ecotype_summary <- gene_selected %>%
    group_by(gene, Origin, cc_ischia_10) %>%
    summarise(
      n_spots = n(),
      n_detected = sum(detected, na.rm = TRUE),
      pct_detected = round(100 * n_detected / n_spots, 2),
      mean_expression = mean(expression, na.rm = TRUE),
      median_expression = median(expression, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(gene, Origin, cc_ischia_10)

  sample_summaries[[gene_name]] <- sample_summary
  ecotype_summaries[[gene_name]] <- ecotype_summary

  plot_list <- lapply(selected_images$orig.ident, function(img_name) {
    build_spatial_expression_plot(
      seurat_obj = st_obj,
      image_name = img_name,
      gene_df = gene_df,
      feature_name = gene_name,
      point_size = 0.65,
      limits = c(0, upper_limit),
      palette_values = expression_palette
    )
  })

  save_multipage_pdf(
    plot_list = plot_list,
    filename = paste0("", gene_name, "_expression_map_normal_and_primary.pdf"),
    out_dir = projection_fig_dir,
    width = 4.4,
    height = 3.6
  )
}

sample_summary_all <- bind_rows(sample_summaries)
ecotype_summary_all <- bind_rows(ecotype_summaries)

write_tsv_simple(sample_summary_all, "kcn_expression_projection_by_section.tsv", tab_dir)
write_tsv_simple(ecotype_summary_all, "kcn_expression_projection_by_ecotype.tsv", tab_dir)


message("ST PDAC KCN25 expression projections complete.")



