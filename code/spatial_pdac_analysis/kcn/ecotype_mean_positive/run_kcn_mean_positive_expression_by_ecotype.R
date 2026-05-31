#!/usr/bin/env Rscript

###############################################################################
# Mean KCN expression in positive spots by ecotype
#
# For each KCN in the 25-gene list, this script calculates the mean expression
# among positive spots only inside each spatial ecotype (`cc_ischia_10`),
# separately for normal and primary pancreas sections.
#
# This is useful next to the percent-positive metric because it shows how strong
# the signal is where the gene is actually detected.
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

save_two_page_pdf <- function(plot_normal, plot_primary, filename, out_dir, width = 6.8, height = 4.8) {
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

theme_barplot <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
}

make_mean_positive_barplot <- function(df, gene_name, origin_name, palette_ecotype) {
  ggplot(df, aes(x = cc_ischia_10, y = mean_positive_expression_plot, fill = cc_ischia_10)) +
    geom_col(width = 0.82, color = "grey30", linewidth = 0.2) +
    scale_fill_manual(values = palette_ecotype, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
    labs(
      title = paste0(gene_name, " in ", origin_name),
      subtitle = "Mean expression among positive spots only; ecotypes with no positive spots are shown as 0",
      x = "Ecotype (cc_ischia_10)",
      y = "Mean positive expression"
    ) +
    theme_barplot(11)
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

message("Loading ST Seurat object for KCN25 mean-positive-expression barplots...")
st_obj <- readRDS(object_path)
st_obj <- suppressWarnings(UpdateSeuratObject(st_obj))

metadata <- st_obj@meta.data %>%
  rownames_to_column("spot") %>%
  filter(Origin %in% c("Normal Pancreas", "Pancreas"))

kcn_list <- read.delim(kcn_list_path, stringsAsFactors = FALSE) %>%
  pull(gene) %>%
  unique()

spatial_mat <- get_assay_matrix_compat(st_obj, assay = "Spatial", layer_name = "data")
ecotype_levels <- paste0("CC", seq_len(10))
origin_levels <- c("Normal Pancreas", "Pancreas")

palette_ecotype <- c(
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
    inner_join(metadata %>% select(spot, Origin, cc_ischia_10), by = "spot")

  summary_df <- gene_df %>%
    group_by(gene, Origin, cc_ischia_10) %>%
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
        cc_ischia_10 = ecotype_levels,
        stringsAsFactors = FALSE
      ),
      by = c("gene", "Origin", "cc_ischia_10")
    ) %>%
    mutate(
      n_spots = ifelse(is.na(n_spots), 0, n_spots),
      n_positive = ifelse(is.na(n_positive), 0, n_positive),
      pct_detected = ifelse(is.na(pct_detected), 0, pct_detected),
      mean_expression_all_spots = ifelse(is.na(mean_expression_all_spots), 0, mean_expression_all_spots),
      mean_positive_expression_plot = ifelse(is.na(mean_positive_expression), 0, mean_positive_expression)
    ) %>%
    arrange(factor(Origin, levels = origin_levels), factor(cc_ischia_10, levels = ecotype_levels))

  summary_list[[gene_name]] <- summary_df

  plot_normal <- summary_df %>%
    filter(Origin == "Normal Pancreas") %>%
    mutate(cc_ischia_10 = factor(cc_ischia_10, levels = ecotype_levels)) %>%
    make_mean_positive_barplot(gene_name = gene_name, origin_name = "Normal Pancreas", palette_ecotype = palette_ecotype)

  plot_primary <- summary_df %>%
    filter(Origin == "Pancreas") %>%
    mutate(cc_ischia_10 = factor(cc_ischia_10, levels = ecotype_levels)) %>%
    make_mean_positive_barplot(gene_name = gene_name, origin_name = "Primary Pancreas", palette_ecotype = palette_ecotype)

  save_two_page_pdf(
    plot_normal = plot_normal,
    plot_primary = plot_primary,
    filename = paste0("", gene_name, "_mean_positive_expression_by_ecotype.pdf"),
    out_dir = fig_dir,
    width = 6.8,
    height = 4.8
  )
}

summary_all <- bind_rows(summary_list)
write_tsv_simple(summary_all, "kcn_mean_positive_expression_by_ecotype.tsv", tab_dir)


message("ST PDAC KCN25 mean-positive-expression barplots complete.")



