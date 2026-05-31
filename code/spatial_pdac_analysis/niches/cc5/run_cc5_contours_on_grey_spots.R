#!/usr/bin/env Rscript

###############################################################################
# CC5 contours on grey spot maps
#
# This script shows all spots in grey and draws only the CC5 boundaries as thin
# black contours.
###############################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(ggplot2)
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

save_multipage_pdf <- function(plot_list, filename, out_dir, width = 4.0, height = 3.2) {
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

  cc_df$component_id <- connected_components_from_threshold(
    coords_mat,
    distance_threshold = nn_scale * 1.6
  )

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

    data.frame(
      x_plot = hull_coords[, "X"],
      y_plot = hull_coords[, "Y"],
      component_id = df_comp$component_id[1],
      n_spots_component = nrow(df_comp),
      stringsAsFactors = FALSE
    )
  })

  bind_rows(hull_list)
}

build_contours_grey_spots_plot <- function(seurat_obj, image_name, metadata_df) {
  image_obj <- seurat_obj@images[[image_name]]
  image_height <- dim(image_obj@image)[1]
  scale_factor <- image_obj@scale.factors$hires

  coords_df <- image_obj@coordinates %>%
    rownames_to_column("spot") %>%
    left_join(metadata_df %>% select(spot, cc_ischia_10), by = "spot") %>%
    mutate(
      x_plot = imagecol * scale_factor,
      y_plot = image_height - (imagerow * scale_factor)
    )

  hull_df <- compute_component_hulls(coords_df, target_label = "CC5")

  plot_obj <- ggplot(coords_df, aes(x = x_plot, y = y_plot)) +
    geom_point(color = "grey78", size = 0.55, alpha = 0.85) +
    coord_fixed(expand = FALSE) +
    labs(title = image_name) +
    theme_void(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 10)
    )

  if (nrow(hull_df) > 0) {
    plot_obj <- plot_obj +
      geom_path(
        data = hull_df,
        aes(x = x_plot, y = y_plot, group = component_id),
        inherit.aes = FALSE,
        color = "black",
        linewidth = 0.35,
        lineend = "round"
      )
  }

  plot_obj
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

message("Loading ST Seurat object for CC5 contour maps with grey spots...")
st_obj <- readRDS(object_path)
st_obj <- suppressWarnings(UpdateSeuratObject(st_obj))

metadata <- st_obj@meta.data %>%
  rownames_to_column("spot")

selected_images <- metadata %>%
  filter(Origin %in% c("Normal Pancreas", "Pancreas")) %>%
  distinct(orig.ident, Origin) %>%
  arrange(factor(Origin, levels = c("Normal Pancreas", "Pancreas")), orig.ident)

write_tsv_simple(selected_images, "selected_normal_and_primary_sections_for_cc5_contours_on_grey_spots.tsv", tab_dir)

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

  plot_list[[img_name]] <- build_contours_grey_spots_plot(
    seurat_obj = st_obj,
    image_name = img_name,
    metadata_df = metadata
  )
}

write_tsv_simple(bind_rows(component_summary), "cc5_contour_components_on_grey_spots.tsv", tab_dir)

save_multipage_pdf(
  plot_list = plot_list,
  filename = "cc5_contours_on_grey_spots_normal_and_primary.pdf",
  out_dir = fig_dir,
  width = 4.0,
  height = 3.2
)


message("ST PDAC CC5 contour maps with grey spots complete.")


