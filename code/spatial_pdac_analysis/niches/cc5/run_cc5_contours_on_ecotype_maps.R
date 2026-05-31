#!/usr/bin/env Rscript

###############################################################################
# Ecotype maps with CC5 contours
#
# This script keeps the ecotype colors and adds a black contour around the CC5
# spot islands. The contour follows connected CC5 clusters, not single spots.
###############################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(sf)
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

build_manual_color_legend <- function(color_mapping_df, title_text = "cc_ischia_10", label_col = "label") {
  legend_df <- color_mapping_df %>%
    mutate(y = rev(seq_len(n())))

  ggplot(legend_df, aes(x = 1, y = y)) +
    geom_tile(aes(fill = .data[[label_col]]), width = 0.35, height = 0.35) +
    geom_text(aes(x = 1.35, label = .data[[label_col]]), hjust = 0, size = 3.4) +
    scale_fill_manual(values = setNames(legend_df$color, legend_df[[label_col]])) +
    guides(fill = "none") +
    coord_cartesian(xlim = c(0.7, 4.3), ylim = c(0.5, nrow(legend_df) + 0.8), clip = "off") +
    labs(title = title_text) +
    theme_void(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, size = 10),
      plot.margin = margin(10, 10, 10, 10)
    )
}

nearest_neighbor_scale <- function(coords_mat) {
  n <- nrow(coords_mat)
  if (n <= 1) {
    return(0)
  }

  dmat <- as.matrix(stats::dist(coords_mat))
  diag(dmat) <- Inf
  median(apply(dmat, 1, min), na.rm = TRUE)
}

connected_components_from_threshold <- function(coords_mat, distance_threshold) {
  n <- nrow(coords_mat)
  if (n == 0) return(integer(0))
  if (n == 1) return(1L)

  dmat <- as.matrix(stats::dist(coords_mat))
  adj <- dmat <= distance_threshold
  diag(adj) <- TRUE

  component <- rep(NA_integer_, n)
  comp_id <- 0L

  for (i in seq_len(n)) {
    if (!is.na(component[i])) next
    comp_id <- comp_id + 1L
    queue <- i
    component[i] <- comp_id
    while (length(queue) > 0) {
      current <- queue[1]
      queue <- queue[-1]
      neigh <- which(adj[current, ] & is.na(component))
      if (length(neigh) > 0) {
        component[neigh] <- comp_id
        queue <- c(queue, neigh)
      }
    }
  }

  component
}

compute_component_hulls <- function(coords_df, target_label = "CC5") {
  cc_df <- coords_df %>%
    filter(cc_ischia_10 == target_label) %>%
    select(spot, x_plot, y_plot)

  if (nrow(cc_df) < 3) {
    return(data.frame())
  }

  coords_mat <- as.matrix(cc_df[, c("x_plot", "y_plot")])
  nn_scale <- nearest_neighbor_scale(coords_mat)
  if (!is.finite(nn_scale) || nn_scale <= 0) {
    return(data.frame())
  }

  component_id <- connected_components_from_threshold(coords_mat, distance_threshold = nn_scale * 1.6)
  cc_df$component_id <- component_id

  hull_list <- lapply(split(cc_df, cc_df$component_id), function(df_comp) {
    if (nrow(df_comp) < 3) {
      return(NULL)
    }

    pts_sf <- sf::st_as_sf(df_comp, coords = c("x_plot", "y_plot"), remove = FALSE)
    hull_geom <- tryCatch(
      {
        sf::st_concave_hull(
          sf::st_union(pts_sf),
          ratio = 0.25,
          allow_holes = FALSE
        )
      },
      error = function(e) sf::st_convex_hull(sf::st_union(pts_sf))
    )

    hull_coords <- sf::st_coordinates(hull_geom)
    if (nrow(hull_coords) == 0) {
      return(NULL)
    }

    hull_df <- data.frame(
      x_plot = hull_coords[, "X"],
      y_plot = hull_coords[, "Y"],
      component_id = df_comp$component_id[1],
      n_spots_component = nrow(df_comp),
      stringsAsFactors = FALSE
    )
    hull_df
  })

  bind_rows(hull_list)
}

build_spatial_cc5_plot <- function(seurat_obj, image_name, metadata_df, palette_values, color_mapping_df, point_size = 0.65) {
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

  coords_df$cc_ischia_10 <- factor(as.character(coords_df$cc_ischia_10), levels = names(palette_values))
  hull_df <- compute_component_hulls(coords_df, target_label = "CC5")

  spatial_plot <- ggplot(coords_df, aes(x = x_plot, y = y_plot)) +
    geom_point(aes(color = cc_ischia_10), size = point_size, alpha = 0.70) +
    scale_color_manual(values = palette_values, na.value = "grey80") +
    coord_fixed(expand = FALSE) +
    labs(title = image_name, color = "cc_ischia_10") +
    theme_void(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
      legend.position = "none"
    )

  if (nrow(hull_df) > 0) {
    spatial_plot <- spatial_plot +
      geom_path(
        data = hull_df,
        aes(x = x_plot, y = y_plot, group = component_id),
        inherit.aes = FALSE,
        color = "black",
        linewidth = 0.7,
        lineend = "round"
      )
  }

  spatial_plot +
    build_manual_color_legend(
      color_mapping_df = color_mapping_df,
      title_text = "cc_ischia_10"
    ) +
    patchwork::plot_layout(widths = c(3.6, 1.8))
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

message("Loading ST Seurat object for CC5 contour maps...")
st_obj <- readRDS(object_path)
st_obj <- suppressWarnings(UpdateSeuratObject(st_obj))

metadata <- st_obj@meta.data %>%
  rownames_to_column("spot")

selected_images <- metadata %>%
  filter(Origin %in% c("Normal Pancreas", "Pancreas")) %>%
  distinct(orig.ident, Origin) %>%
  arrange(factor(Origin, levels = c("Normal Pancreas", "Pancreas")), orig.ident)

write_tsv_simple(selected_images, "selected_normal_and_primary_sections_for_cc5_contours.tsv", tab_dir)

figure_color_mapping <- data.frame(
  label = paste0("CC", seq_len(10)),
  color = c(
    "#3ab44b", "#bcf60b", "#f031e6", "#901eb4", "#fabdbe",
    "#e7872b", "#e6194b", "#317ec2", "#ffe119", "#11dfed"
  ),
  stringsAsFactors = FALSE
)

ecotype_colors <- setNames(figure_color_mapping$color, figure_color_mapping$label)
write_tsv_simple(figure_color_mapping, "ecotype_color_mapping.tsv", tab_dir)

component_summary <- list()
plot_list <- list()

for (img_name in selected_images$orig.ident) {
  image_obj <- st_obj@images[[img_name]]
  image_height <- dim(image_obj@image)[1]
  scale_factor <- image_obj@scale.factors$hires

  coords_df <- image_obj@coordinates %>%
    rownames_to_column("spot") %>%
    left_join(metadata %>% select(spot, cc_ischia_10), by = "spot") %>%
    mutate(
      x_plot = imagecol * scale_factor,
      y_plot = image_height - (imagerow * scale_factor)
    )

  hull_df <- compute_component_hulls(coords_df, target_label = "CC5")
  if (nrow(hull_df) > 0) {
    component_summary[[img_name]] <- hull_df %>%
      distinct(component_id, n_spots_component) %>%
      mutate(image = img_name)
  } else {
    component_summary[[img_name]] <- data.frame(
      image = character(0),
      component_id = integer(0),
      n_spots_component = integer(0)
    )
  }

  plot_list[[img_name]] <- build_spatial_cc5_plot(
    seurat_obj = st_obj,
    image_name = img_name,
    metadata_df = metadata,
    palette_values = ecotype_colors,
    color_mapping_df = figure_color_mapping,
    point_size = 0.65
  )
}

write_tsv_simple(bind_rows(component_summary), "cc5_contour_components_on_ecotype_maps.tsv", tab_dir)

save_multipage_pdf(
  plot_list = plot_list,
  filename = "cc5_contours_on_ecotype_maps_normal_and_primary.pdf",
  out_dir = fig_dir,
  width = 4.6,
  height = 3.8
)


message("ST PDAC CC5 contour maps complete.")


