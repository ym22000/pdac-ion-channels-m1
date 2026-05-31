#!/usr/bin/env Rscript

###############################################################################
# Quick overview of the ST PDAC Seurat object
#
# This script is a compact entry point for the spatial dataset. It summarizes
# the object, the section-level metadata, and the UMAP embeddings already stored
# in the curated Seurat object.
#
# It does not rebuild the atlas. It only reads the existing object and writes
# a clean overview of what is inside.
###############################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(forcats)
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
rds_dir <- file.path(script_dir, "rds")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

object_path <- file.path(input_dir, "PDAC_Updated_ST.rds")
paper_path <- file.path(input_dir, "spatial_transcriptomics_PDAC_paper.pdf")

if (!file.exists(object_path)) stop("Missing input object: ", object_path)
if (!file.exists(paper_path)) stop("Missing local paper PDF: ", paper_path)

theme_stage <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "grey95", colour = "grey70"),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.title = element_text(face = "bold")
    )
}

save_plot <- function(plot_obj, filename, width = 10, height = 8) {
  pdf_path <- file.path(fig_dir, filename)
  png_path <- sub("\\.pdf$", ".png", pdf_path, ignore.case = TRUE)
  ggsave(pdf_path, plot = plot_obj, width = width, height = height, units = "in", dpi = 300, bg = "white")
  ggsave(png_path, plot = plot_obj, width = width, height = height, units = "in", dpi = 300, bg = "white")
}

save_multipage_pdf <- function(plot_list, filename, width = 8.5, height = 8.5) {
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

extract_reduction_df <- function(seurat_obj, reduction_name = "umap", prefix = "umap") {
  emb <- as.data.frame(Embeddings(seurat_obj, reduction_name))
  out <- emb[, 1:2, drop = FALSE]
  colnames(out) <- c(paste0(prefix, "_1"), paste0(prefix, "_2"))
  rownames_to_column(out, "spot")
}

get_assay_matrix_compat <- function(seurat_obj, assay = "Spatial", layer_name = "data") {
  tryCatch(
    Seurat::GetAssayData(seurat_obj, assay = assay, layer = layer_name),
    error = function(e) Seurat::GetAssayData(seurat_obj, assay = assay, slot = layer_name)
  )
}

compress_top_labels <- function(x, top_n = 12, other_label = "Other") {
  x <- as.character(x)
  top_levels <- names(sort(table(x), decreasing = TRUE))[seq_len(min(top_n, length(unique(x))))]
  ifelse(x %in% top_levels, x, other_label)
}

build_spatial_category_plot <- function(seurat_obj, image_name, metadata_df, group_var, palette_values, point_size = 0.65, show_legend = FALSE, x_limits = NULL, y_limits = NULL) {
  image_obj <- seurat_obj@images[[image_name]]
  image_array <- image_obj@image
  image_height <- dim(image_array)[1]
  scale_factor <- image_obj@scale.factors$hires

  coords_df <- image_obj@coordinates %>%
    rownames_to_column("spot") %>%
    left_join(metadata_df %>% select(spot, all_of(group_var)), by = "spot") %>%
    mutate(
      x_plot = imagecol * scale_factor,
      y_plot = image_height - (imagerow * scale_factor)
    )

  coords_df[[group_var]] <- factor(as.character(coords_df[[group_var]]), levels = names(palette_values))

  ggplot(coords_df, aes(x = x_plot, y = y_plot)) +
    geom_point(aes(color = .data[[group_var]]), size = point_size, alpha = 0.70) +
    scale_color_manual(values = palette_values, na.value = "grey80") +
    coord_fixed(xlim = x_limits, ylim = y_limits, expand = FALSE) +
    labs(title = image_name, color = group_var) +
    theme_void(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
      legend.position = if (show_legend) "right" else "none",
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 8)
    )
}

build_spatial_feature_plot <- function(seurat_obj, image_name, expr_df, feature_name, point_size = 0.65, x_limits = NULL, y_limits = NULL) {
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
    coord_fixed(xlim = x_limits, ylim = y_limits, expand = FALSE) +
    labs(title = image_name, color = feature_name) +
    theme_void(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
      legend.position = "none"
    )
}

build_manual_color_legend <- function(color_mapping_df, title_text = "Dominant cell state", label_col = "label") {
  legend_df <- color_mapping_df %>%
    mutate(y = rev(seq_len(n())))

  ggplot(legend_df, aes(x = 1, y = y)) +
    geom_tile(aes(fill = .data[[label_col]]), width = 0.35, height = 0.35) +
    geom_text(aes(x = 1.35, label = .data[[label_col]]), hjust = 0, size = 3.6) +
    scale_fill_manual(values = setNames(legend_df$color, legend_df[[label_col]])) +
    guides(fill = "none") +
    coord_cartesian(xlim = c(0.7, 4.7), ylim = c(0.5, nrow(legend_df) + 0.8), clip = "off") +
    labs(title = title_text) +
    theme_void(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, size = 10),
      plot.margin = margin(10, 10, 10, 10)
    )
}

message("Loading ST Seurat object...")
st_obj <- readRDS(object_path)

if (!inherits(st_obj, "Seurat")) {
  stop("The input object is not a Seurat object.")
}

message("Updating ST object to the current Seurat structure for spatial plotting...")
st_obj <- suppressWarnings(UpdateSeuratObject(st_obj))

metadata <- st_obj@meta.data %>%
  rownames_to_column("spot")

umap_df <- extract_reduction_df(st_obj, reduction_name = "umap", prefix = "umap")
metadata <- metadata %>%
  left_join(umap_df, by = "spot")

assay_summary <- bind_rows(lapply(names(st_obj@assays), function(assay_name) {
  aa <- st_obj[[assay_name]]
  feats <- rownames(aa)
  notes <- dplyr::case_when(
    assay_name == "Spatial" ~ "Raw spatial transcript counts.",
    assay_name == "SCT" ~ "SCTransform-normalized gene expression.",
    assay_name == "integrated" ~ "Integrated expression space used for joint analysis.",
    grepl("^rctd", assay_name) ~ "Cell-type deconvolution assay inferred from the assay name and feature labels.",
    assay_name == "fges" ~ "Functional or pathway-like spot scores inferred from the assay name and feature labels.",
    TRUE ~ "Assay stored in the object."
  )
  data.frame(
    assay = assay_name,
    n_features = length(feats),
    n_spots = ncol(aa),
    example_features = paste(head(feats, 10), collapse = ", "),
    notes = notes,
    stringsAsFactors = FALSE
  )
}))

metadata_columns <- data.frame(
  column = colnames(metadata),
  class = vapply(metadata, function(x) paste(class(x), collapse = ", "), character(1)),
  n_unique = vapply(metadata, function(x) dplyr::n_distinct(x), numeric(1)),
  n_missing = vapply(metadata, function(x) sum(is.na(x)), numeric(1)),
  stringsAsFactors = FALSE
)

dataset_overview <- data.frame(
  metric = c(
    "n_spots",
    "n_assays",
    "n_images_sections",
    "n_slide_names",
    "default_assay",
    "available_reductions"
  ),
  value = c(
    ncol(st_obj),
    length(names(st_obj@assays)),
    length(unique(metadata$orig.ident)),
    length(unique(metadata$SlideName)),
    DefaultAssay(st_obj),
    paste(names(st_obj@reductions), collapse = ", ")
  ),
  stringsAsFactors = FALSE
)

image_summary <- metadata %>%
  group_by(orig.ident) %>%
  summarise(
    n_spots = n(),
    dominant_origin = names(sort(table(Origin), decreasing = TRUE))[1],
    dominant_slide = names(sort(table(SlideName), decreasing = TRUE))[1],
    n_slide_names = n_distinct(SlideName),
    .groups = "drop"
  ) %>%
  arrange(desc(n_spots), orig.ident)

sample_mapping <- metadata %>%
  group_by(orig.ident) %>%
  summarise(
    origin = names(sort(table(Origin), decreasing = TRUE))[1],
    slide_name = names(sort(table(SlideName), decreasing = TRUE))[1],
    sample_id2 = names(sort(table(Sample_ID2), decreasing = TRUE))[1],
    n_spots = n(),
    .groups = "drop"
  ) %>%
  mutate(
    interpretation = case_when(
      origin == "Pancreas" ~ "Primary PDAC tissue section",
      origin == "Normal Pancreas" ~ "Normal pancreas tissue section",
      origin == "Liver" ~ "Liver metastasis tissue section",
      origin == "Lymph node" ~ "Lymph node metastasis tissue section",
      TRUE ~ "Other tissue origin"
    )
  ) %>%
  arrange(origin, orig.ident)

origin_counts <- metadata %>%
  count(Origin, name = "n_spots") %>%
  mutate(pct_spots = round(100 * n_spots / sum(n_spots), 2)) %>%
  arrange(desc(n_spots))

spot_class_counts <- metadata %>%
  count(spot_class, name = "n_spots") %>%
  mutate(pct_spots = round(100 * n_spots / sum(n_spots), 2)) %>%
  arrange(desc(n_spots))

ecotype_counts <- metadata %>%
  count(CompositionCluster_CC, name = "n_spots") %>%
  mutate(pct_spots = round(100 * n_spots / sum(n_spots), 2)) %>%
  arrange(desc(n_spots))

first_type_counts <- metadata %>%
  count(first_type, name = "n_spots") %>%
  mutate(pct_spots = round(100 * n_spots / sum(n_spots), 2)) %>%
  arrange(desc(n_spots))

metadata$first_type_compact <- as.character(compress_top_labels(metadata$first_type, top_n = 12))

origin_by_tissue <- metadata %>%
  count(orig.ident, Origin, name = "n_spots") %>%
  group_by(orig.ident) %>%
  mutate(frac_spots = n_spots / sum(n_spots)) %>%
  ungroup()

ecotype_by_tissue <- metadata %>%
  count(orig.ident, CompositionCluster_CC, name = "n_spots") %>%
  group_by(orig.ident) %>%
  mutate(frac_spots = n_spots / sum(n_spots)) %>%
  ungroup()

first_type_by_tissue <- metadata %>%
  count(orig.ident, first_type_compact, name = "n_spots") %>%
  group_by(orig.ident) %>%
  mutate(frac_spots = n_spots / sum(n_spots)) %>%
  ungroup()

key_metadata <- metadata %>%
  select(
    spot, orig.ident, SlideName, Origin, spot_class,
    first_type, second_type, first_class, second_class,
    CompositionCluster_CC, cc_ischia, cc_ischia_10, cc_ischia_12, cc_ischia_16, cc_ischia_18,
    umap_1, umap_2
  )

write_tsv_simple(dataset_overview, "dataset_overview.tsv")
write_tsv_simple(assay_summary, "assay_summary.tsv")
write_tsv_simple(metadata_columns, "metadata_columns.tsv")
write_tsv_simple(image_summary, "image_summary.tsv")
write_tsv_simple(sample_mapping, "sample_mapping.tsv")
write_tsv_simple(origin_counts, "origin_counts.tsv")
write_tsv_simple(spot_class_counts, "spot_class_counts.tsv")
write_tsv_simple(ecotype_counts, "ecotype_counts.tsv")
write_tsv_simple(first_type_counts, "first_type_counts.tsv")
write_tsv_simple(origin_by_tissue, "origin_by_tissue.tsv")
write_tsv_simple(ecotype_by_tissue, "ecotype_by_tissue.tsv")
write_tsv_simple(first_type_by_tissue, "first_type_by_tissue_top12.tsv")
write_tsv_simple(key_metadata, "spot_metadata_key.tsv")

origin_colors <- c(
  "Pancreas" = "#B22222",
  "Liver" = "#1F78B4",
  "Lymph node" = "#33A02C",
  "Normal Pancreas" = "#6A3D9A"
)

ecotype_levels <- sort(unique(as.character(metadata$cc_ischia_10)))
ecotype_colors <- setNames(grDevices::hcl.colors(length(ecotype_levels), "Spectral"), ecotype_levels)

first_type_levels <- names(sort(table(metadata$first_type_compact), decreasing = TRUE))
first_type_colors <- setNames(grDevices::hcl.colors(length(first_type_levels), "Dark 3"), first_type_levels)

p_umap_origin <- ggplot(metadata, aes(x = umap_1, y = umap_2, color = Origin)) +
  geom_point(size = 0.20, alpha = 0.75) +
  scale_color_manual(values = origin_colors) +
  labs(
    title = "ST PDAC atlas UMAP",
    subtitle = "Spots colored by anatomical origin",
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Origin"
  ) +
  theme_stage(11)

save_plot(p_umap_origin, "umap_by_origin.pdf", width = 10, height = 8)

p_umap_ecotype <- ggplot(metadata, aes(x = umap_1, y = umap_2, color = cc_ischia_10)) +
  geom_point(size = 0.20, alpha = 0.80) +
  scale_color_manual(values = ecotype_colors) +
  labs(
    title = "ST PDAC atlas UMAP",
    subtitle = "Spots colored by the 10 spatial ecotypes used in the paper (cc_ischia_10)",
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Ecotype"
  ) +
  theme_stage(10)

save_plot(p_umap_ecotype, "umap_by_ecotype.pdf", width = 11, height = 8.5)

p_umap_first_type <- ggplot(metadata, aes(x = umap_1, y = umap_2, color = first_type_compact)) +
  geom_point(size = 0.20, alpha = 0.80) +
  scale_color_manual(values = first_type_colors) +
  labs(
    title = "ST PDAC atlas UMAP",
    subtitle = "Spots colored by the dominant deconvolved cell state (top 12 labels + Other)",
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Dominant state"
  ) +
  theme_stage(10)

save_plot(p_umap_first_type, "umap_by_cell_type.pdf", width = 11, height = 8.5)

p_spots_per_tissue <- image_summary %>%
  mutate(orig.ident = forcats::fct_reorder(orig.ident, n_spots)) %>%
  ggplot(aes(x = orig.ident, y = n_spots, fill = dominant_origin)) +
  geom_col(width = 0.85) +
  scale_fill_manual(values = origin_colors) +
  coord_flip() +
  labs(
    title = "Spots per spatial image",
    subtitle = "Thirty spatial images/sections stored in the object",
    x = "Image / section (orig.ident)",
    y = "Number of spots",
    fill = "Dominant origin"
  ) +
  theme_stage(11)

save_plot(p_spots_per_tissue, "spots_per_section.pdf", width = 8.5, height = 10)

p_origin_by_tissue <- origin_by_tissue %>%
  mutate(orig.ident = forcats::fct_reorder(orig.ident, frac_spots, .fun = sum)) %>%
  ggplot(aes(x = orig.ident, y = frac_spots, fill = Origin)) +
  geom_col(width = 0.85) +
  scale_fill_manual(values = origin_colors) +
  coord_flip() +
  labs(
    title = "Origin composition across the 30 spatial images",
    x = "Image / section (orig.ident)",
    y = "Fraction of spots",
    fill = "Origin"
  ) +
  theme_stage(11)

save_plot(p_origin_by_tissue, "origin_composition_by_section.pdf", width = 9, height = 10)

exploration_bundle <- list(
  dataset_overview = dataset_overview,
  assay_summary = assay_summary,
  metadata_columns = metadata_columns,
  image_summary = image_summary,
  origin_counts = origin_counts,
  spot_class_counts = spot_class_counts,
  ecotype_counts = ecotype_counts,
  first_type_counts = first_type_counts
)

saveRDS(exploration_bundle, file = file.path(rds_dir, "st_pdac_exploration_bundle.rds"))

summary_lines <- c(
  "ST PDAC exploration summary",
  paste0("Date: ", Sys.Date()),
  "",
  "Sources:",
  "- Nature Genetics 2024: https://www.nature.com/articles/s41588-024-01914-4",
  "- GitHub: https://github.com/Masood-Lab/PDAC_Mets",
  "",
  paste0("Spots: ", ncol(st_obj)),
  paste0("Assays: ", paste(names(st_obj@assays), collapse = ", ")),
  paste0("Spatial images / sections (orig.ident): ", length(unique(metadata$orig.ident))),
  paste0("SlideName values: ", length(unique(metadata$SlideName))),
  paste0("Origins: ", paste(names(sort(table(metadata$Origin), decreasing = TRUE)), collapse = ", ")),
  paste0("Ecotypes (CompositionCluster_CC): ", length(unique(metadata$CompositionCluster_CC))),
  paste0("Dominant cell-state labels (first_type): ", length(unique(metadata$first_type)))
)

writeLines(summary_lines, con = file.path(script_dir, "summary.txt"))


message("ST PDAC exploration complete.")

