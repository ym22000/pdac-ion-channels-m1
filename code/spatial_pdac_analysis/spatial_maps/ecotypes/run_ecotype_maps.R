#!/usr/bin/env Rscript

###############################################################################
# Spatial ecotype maps for the selected normal and primary sections
#
# This script plots the published `cc_ischia_10` labels back onto the tissue
# sections. These maps are the main spatial reference used later for KCN
# projections and niche interpretation.
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

save_multipage_pdf <- function(plot_list, filename, out_dir, width = 4.4, height = 3.6) {
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

build_spatial_category_plot <- function(seurat_obj, image_name, metadata_df, palette_values, point_size = 0.65) {
  image_obj <- seurat_obj@images[[image_name]]
  image_array <- image_obj@image
  image_height <- dim(image_array)[1]
  scale_factor <- image_obj@scale.factors$hires

  coords_df <- image_obj@coordinates %>%
    rownames_to_column("spot") %>%
    left_join(metadata_df %>% select(spot, cc_ischia_10), by = "spot") %>%
    mutate(
      x_plot = imagecol * scale_factor,
      y_plot = image_height - (imagerow * scale_factor)
    )

  coords_df$cc_ischia_10 <- factor(coords_df$cc_ischia_10, levels = names(palette_values))

  ggplot(coords_df, aes(x = x_plot, y = y_plot)) +
    geom_point(aes(color = cc_ischia_10), size = point_size, alpha = 0.70) +
    scale_color_manual(values = palette_values, na.value = "grey80") +
    coord_fixed(expand = FALSE) +
    labs(title = image_name, color = "Ecotype") +
    theme_void(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
      legend.position = "none"
    )
}

build_manual_color_legend <- function(color_mapping_df, title_text = "Ecotype") {
  legend_df <- color_mapping_df %>%
    mutate(y = rev(seq_len(n())))

  ggplot(legend_df, aes(x = 1, y = y)) +
    geom_tile(aes(fill = label), width = 0.35, height = 0.35) +
    geom_text(aes(x = 1.35, label = label), hjust = 0, size = 3.4) +
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

script_dir <- get_script_dir()
project_dir <- find_project_dir(script_dir)
input_dir <- file.path(project_dir, "inputs")
fig_dir <- file.path(script_dir, "figures")
tab_dir <- file.path(script_dir, "tables")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

object_path <- file.path(input_dir, "PDAC_Updated_ST.rds")
if (!file.exists(object_path)) stop("Missing input object: ", object_path)

message("Loading ST Seurat object for ecotype maps...")
st_obj <- readRDS(object_path)
st_obj <- suppressWarnings(UpdateSeuratObject(st_obj))

metadata <- st_obj@meta.data %>%
  rownames_to_column("spot") %>%
  mutate(cc_ischia_10 = as.character(cc_ischia_10))

selected_images <- metadata %>%
  filter(Origin %in% c("Normal Pancreas", "Pancreas")) %>%
  distinct(orig.ident, Origin) %>%
  arrange(factor(Origin, levels = c("Normal Pancreas", "Pancreas")), orig.ident)

write_tsv_simple(selected_images, "selected_normal_and_primary_sections.tsv", tab_dir)

ecotype_levels <- paste0("CC", seq_len(10))
ecotype_colors <- c(
  CC1 = "#3ab44b",
  CC2 = "#bcf60b",
  CC3 = "#f031e6",
  CC4 = "#901eb4",
  CC5 = "#fabdbe",
  CC6 = "#e7872b",
  CC7 = "#e6194b",
  CC8 = "#317ec2",
  CC9 = "#ffe119",
  CC10 = "#11dfed"
)

color_mapping <- data.frame(
  label = ecotype_levels,
  color = unname(ecotype_colors[ecotype_levels]),
  stringsAsFactors = FALSE
)
write_tsv_simple(color_mapping, "ecotype_color_mapping.tsv", tab_dir)

panel_list <- lapply(selected_images$orig.ident, function(image_name) {
  spatial_plot <- build_spatial_category_plot(
    seurat_obj = st_obj,
    image_name = image_name,
    metadata_df = metadata,
    palette_values = ecotype_colors,
    point_size = 0.65
  )

  spatial_plot + build_manual_color_legend(color_mapping) +
    patchwork::plot_layout(widths = c(3.5, 1.8))
})

save_multipage_pdf(
  plot_list = panel_list,
  filename = "normal_and_primary_ecotype_maps.pdf",
  out_dir = fig_dir,
  width = 4.4,
  height = 3.6
)


message("ST PDAC ecotype maps complete.")
