#!/usr/bin/env Rscript

###############################################################################
# Export normal pancreas spatial sections for Hotspot analyses
#
# This script creates the small section-level files used by the Python Hotspot
# runs. We already used the same export idea for the primary PDAC sections.
# Here the only difference is that we keep the three normal pancreas sections
# only, so the downstream Hotspot and local-correlation analyses stay
# biologicaly clean and easy to compare with the tumor-only run.
###############################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(Matrix)
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

get_assay_matrix_compat <- function(seurat_obj, assay = "Spatial", layer_name = "counts") {
  tryCatch(
    Seurat::GetAssayData(seurat_obj, assay = assay, layer = layer_name),
    error = function(e) Seurat::GetAssayData(seurat_obj, assay = assay, slot = layer_name)
  )
}

write_tsv_simple <- function(x, filename) {
  write.table(
    x,
    file = filename,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE,
    na = "NA"
  )
}

find_project_dir <- function(start_dir) {
  current_dir <- normalizePath(start_dir)
  repeat {
    if (file.exists(file.path(current_dir, "inputs", "PDAC_Updated_ST.rds"))) {
      return(current_dir)
    }
    parent_dir <- dirname(current_dir)
    if (identical(parent_dir, current_dir)) {
      stop("Could not find the ST project directory.")
    }
    current_dir <- parent_dir
  }
}

script_dir <- get_script_dir()
project_dir <- find_project_dir(script_dir)
input_dir <- file.path(project_dir, "inputs")

export_root_dir <- file.path(script_dir, "exports_normal_pancreas")
table_dir <- file.path(export_root_dir, "tables")
export_sections_dir <- file.path(export_root_dir, "exports")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(export_sections_dir, recursive = TRUE, showWarnings = FALSE)

object_path <- file.path(input_dir, "PDAC_Updated_ST.rds")
selected_images_path <- file.path(project_dir, "kcn", "projections", "all_sections", "tables", "selected_normal_and_primary_sections.tsv")
manifest_path <- file.path(table_dir, "normal_pancreas_export_manifest.tsv")

if (!file.exists(object_path)) {
  stop("Missing ST object: ", object_path)
}
if (!file.exists(selected_images_path)) {
  stop("Missing selected sections file: ", selected_images_path)
}

message("Loading spatial Seurat object...")
st_obj <- readRDS(object_path)
st_obj <- suppressWarnings(UpdateSeuratObject(st_obj))

metadata <- st_obj@meta.data %>%
  rownames_to_column("spot")

selected_images <- read.delim(selected_images_path, stringsAsFactors = FALSE) %>%
  filter(Origin == "Normal Pancreas") %>%
  distinct(orig.ident, Origin) %>%
  arrange(orig.ident)

if (nrow(selected_images) == 0) {
  stop("No normal pancreas sections were found in the selected sections file.")
}

counts_mat <- get_assay_matrix_compat(st_obj, assay = "Spatial", layer_name = "counts")
manifest_rows <- list()

for (image_name in selected_images$orig.ident) {
  if (!(image_name %in% names(st_obj@images))) {
    warning("Skipping missing image in Seurat object: ", image_name)
    next
  }

  message("Exporting normal pancreas section: ", image_name)

  image_spots <- metadata %>%
    filter(Origin == "Normal Pancreas", orig.ident == image_name) %>%
    pull(spot)

  coords_df <- st_obj@images[[image_name]]@coordinates %>%
    rownames_to_column("spot") %>%
    filter(spot %in% image_spots)

  ordered_spots <- coords_df$spot
  counts_sub <- counts_mat[, ordered_spots, drop = FALSE]

  meta_sub <- metadata %>%
    filter(spot %in% ordered_spots) %>%
    select(spot, orig.ident, Origin, cc_ischia_10, first_type) %>%
    slice(match(ordered_spots, spot))

  image_export_dir <- file.path(export_sections_dir, image_name)
  dir.create(image_export_dir, recursive = TRUE, showWarnings = FALSE)

  Matrix::writeMM(counts_sub, file.path(image_export_dir, "counts.mtx"))
  writeLines(rownames(counts_sub), file.path(image_export_dir, "genes.tsv"))
  writeLines(colnames(counts_sub), file.path(image_export_dir, "spots.tsv"))
  write_tsv_simple(coords_df, file.path(image_export_dir, "coords.tsv"))
  write_tsv_simple(meta_sub, file.path(image_export_dir, "metadata.tsv"))

  manifest_rows[[length(manifest_rows) + 1]] <- data.frame(
    orig.ident = image_name,
    n_spots = ncol(counts_sub),
    n_genes = nrow(counts_sub),
    export_dir = normalizePath(image_export_dir, winslash = "/"),
    stringsAsFactors = FALSE
  )
}

manifest_df <- bind_rows(manifest_rows)
write_tsv_simple(manifest_df, manifest_path)
write_tsv_simple(selected_images, file.path(table_dir, "selected_normal_pancreas_sections_for_hotspot.tsv"))

message("Normal pancreas exports complete.")
