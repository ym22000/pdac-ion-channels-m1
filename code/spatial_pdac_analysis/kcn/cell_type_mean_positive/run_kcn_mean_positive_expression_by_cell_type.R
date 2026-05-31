#!/usr/bin/env Rscript

###############################################################################
# Mean KCN expression in positive spots by dominant spot identity
#
# For each KCN in the 25-gene list, this script calculates the mean expression
# among positive spots only inside each `first_type` group, separately for
# normal pancreas and primary PDAC sections.
#
# This complements the ecotype summaries by using the dominant spot identity
# stored in the object metadata.
###############################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
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

save_two_page_pdf <- function(plot_normal, plot_primary, filename, out_dir, width = 8.8, height = 5.6) {
  pdf_path <- file.path(out_dir, filename)
  pdf(pdf_path, width = width, height = height, onefile = TRUE)
  on.exit(dev.off(), add = TRUE)
  print(plot_normal)
  print(plot_primary)
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

theme_barplot <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
}

make_mean_positive_barplot <- function(df, gene_name, origin_name, palette_first_type) {
  ggplot(df, aes(x = first_type, y = mean_positive_expression_plot, fill = first_type)) +
    geom_col(width = 0.82, color = "grey30", linewidth = 0.2) +
    scale_fill_manual(values = palette_first_type, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
    labs(
      title = paste0(gene_name, " in ", origin_name),
      subtitle = "Mean expression among positive spots only; cell types with no positive spots are shown as 0",
      x = "Dominant spot identity (first_type)",
      y = "Mean positive expression"
    ) +
    theme_barplot(10)
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
kcn_list_path <- normalizePath(file.path(project_dir, "..", "mycaf_icaf_correlations", "kcn_union_25.tsv"), mustWork = FALSE)

if (!file.exists(object_path)) stop("Missing input object: ", object_path)
if (!file.exists(kcn_list_path)) stop("Missing KCN list: ", kcn_list_path)

message("Loading ST Seurat object for KCN25 first_type mean-positive-expression barplots...")
st_obj <- readRDS(object_path)
st_obj <- suppressWarnings(UpdateSeuratObject(st_obj))

metadata <- st_obj@meta.data %>%
  rownames_to_column("spot") %>%
  filter(Origin %in% c("Normal Pancreas", "Pancreas")) %>%
  mutate(first_type = as.character(first_type))

kcn_list <- read.delim(kcn_list_path, stringsAsFactors = FALSE) %>%
  pull(gene) %>%
  unique()

spatial_mat <- get_assay_matrix_compat(st_obj, assay = "Spatial", layer_name = "data")
first_type_levels <- names(sort(table(metadata$first_type), decreasing = TRUE))
origin_levels <- c("Normal Pancreas", "Pancreas")
palette_first_type <- setNames(grDevices::hcl.colors(length(first_type_levels), "Dark 3"), first_type_levels)

gene_availability <- data.frame(
  gene = kcn_list,
  in_spatial_assay = kcn_list %in% rownames(spatial_mat),
  stringsAsFactors = FALSE
)
write_tsv_simple(gene_availability, "kcn_gene_availability.tsv", tab_dir)

selected_images <- metadata %>%
  distinct(orig.ident, Origin) %>%
  arrange(factor(Origin, levels = origin_levels), orig.ident)
write_tsv_simple(selected_images, "selected_normal_and_primary_sections.tsv", tab_dir)

first_type_counts <- metadata %>%
  count(Origin, first_type, name = "n_spots") %>%
  group_by(Origin) %>%
  mutate(pct_spots = 100 * n_spots / sum(n_spots)) %>%
  ungroup() %>%
  arrange(factor(Origin, levels = origin_levels), desc(n_spots))
write_tsv_simple(first_type_counts, "first_type_counts_normal_and_primary.tsv", tab_dir)

color_mapping <- data.frame(
  first_type = first_type_levels,
  color = unname(palette_first_type[first_type_levels]),
  stringsAsFactors = FALSE
)
write_tsv_simple(color_mapping, "first_type_color_mapping.tsv", tab_dir)

summary_list <- list()

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
    inner_join(metadata %>% select(spot, Origin, first_type), by = "spot")

  summary_df <- gene_df %>%
    group_by(gene, Origin, first_type) %>%
    summarise(
      n_spots = n(),
      n_positive = sum(detected, na.rm = TRUE),
      pct_detected = 100 * n_positive / n_spots,
      mean_expression_all_spots = mean(expression, na.rm = TRUE),
      mean_positive_expression = ifelse(any(detected), mean(expression[detected], na.rm = TRUE), NA_real_),
      median_positive_expression = ifelse(any(detected), median(expression[detected], na.rm = TRUE), NA_real_),
      .groups = "drop"
    ) %>%
    right_join(
      expand.grid(
        gene = gene_name,
        Origin = origin_levels,
        first_type = first_type_levels,
        stringsAsFactors = FALSE
      ),
      by = c("gene", "Origin", "first_type")
    ) %>%
    mutate(
      n_spots = ifelse(is.na(n_spots), 0, n_spots),
      n_positive = ifelse(is.na(n_positive), 0, n_positive),
      pct_detected = ifelse(is.na(pct_detected), 0, pct_detected),
      mean_expression_all_spots = ifelse(is.na(mean_expression_all_spots), 0, mean_expression_all_spots),
      mean_positive_expression_plot = ifelse(is.na(mean_positive_expression), 0, mean_positive_expression)
    ) %>%
    arrange(factor(Origin, levels = origin_levels), factor(first_type, levels = first_type_levels))

  summary_list[[gene_name]] <- summary_df

  plot_normal <- summary_df %>%
    filter(Origin == "Normal Pancreas") %>%
    mutate(first_type = factor(first_type, levels = first_type_levels)) %>%
    make_mean_positive_barplot(gene_name = gene_name, origin_name = "Normal Pancreas", palette_first_type = palette_first_type)

  plot_primary <- summary_df %>%
    filter(Origin == "Pancreas") %>%
    mutate(first_type = factor(first_type, levels = first_type_levels)) %>%
    make_mean_positive_barplot(gene_name = gene_name, origin_name = "Primary Pancreas", palette_first_type = palette_first_type)

  save_two_page_pdf(
    plot_normal = plot_normal,
    plot_primary = plot_primary,
    filename = paste0("", gene_name, "_mean_positive_expression_by_cell_type.pdf"),
    out_dir = fig_dir,
    width = 8.8,
    height = 5.6
  )
}

summary_all <- bind_rows(summary_list)
write_tsv_simple(summary_all, "summary.tsv", tab_dir)


message("ST PDAC KCN25 first_type mean-positive-expression barplots complete.")



