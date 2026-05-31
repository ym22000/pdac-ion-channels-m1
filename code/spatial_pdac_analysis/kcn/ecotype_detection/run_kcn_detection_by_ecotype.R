#!/usr/bin/env Rscript

###############################################################################
# Percent of KCN-positive spots by ecotype
#
# For each KCN in the 25-gene list, this script calculates the percentage of
# positive spots in each ecotype (`cc_ischia_10`), separately for normal and
# primary pancreas sections.
#
# The output is one PDF per gene plus the source tables used to make the plots.
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

save_two_page_pdf <- function(plot_normal, plot_primary, filename, out_dir, width = 6.6, height = 4.6) {
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
      strip.background = element_blank(),
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
}

make_ecotype_barplot <- function(df, gene_name, origin_name, palette_ecotype) {
  ggplot(df, aes(x = cc_ischia_10, y = pct_detected, fill = cc_ischia_10)) +
    geom_col(width = 0.82, color = "grey30", linewidth = 0.2) +
    scale_fill_manual(values = palette_ecotype, drop = FALSE) +
    scale_y_continuous(limits = c(0, max(100, max(df$pct_detected, na.rm = TRUE))), expand = expansion(mult = c(0, 0.03))) +
    labs(
      title = paste0(gene_name, " in ", origin_name),
      subtitle = "Percentage of spots with detectable expression (> 0) within each ecotype",
      x = "Ecotype (cc_ischia_10)",
      y = "Positive spots (%)"
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
barplot_fig_dir <- fig_dir

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(barplot_fig_dir, recursive = TRUE, showWarnings = FALSE)

object_path <- file.path(input_dir, "PDAC_Updated_ST.rds")
kcn_list_path <- normalizePath(file.path(project_dir, "..", "mycaf_icaf_correlations", "kcn_union_25.tsv"), mustWork = FALSE)

if (!file.exists(object_path)) stop("Missing input object: ", object_path)
if (!file.exists(kcn_list_path)) stop("Missing KCN list: ", kcn_list_path)

message("Loading ST Seurat object for KCN25 ecotype barplots...")
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
    count(gene, Origin, cc_ischia_10, wt = detected, name = "n_positive") %>%
    right_join(
      expand.grid(
        gene = gene_name,
        Origin = origin_levels,
        cc_ischia_10 = ecotype_levels,
        stringsAsFactors = FALSE
      ),
      by = c("gene", "Origin", "cc_ischia_10")
    ) %>%
    left_join(
      gene_df %>%
        count(Origin, cc_ischia_10, name = "n_spots"),
      by = c("Origin", "cc_ischia_10")
    ) %>%
    mutate(
      n_positive = ifelse(is.na(n_positive), 0, n_positive),
      n_spots = ifelse(is.na(n_spots), 0, n_spots),
      pct_detected = ifelse(n_spots > 0, 100 * n_positive / n_spots, 0)
    ) %>%
    arrange(factor(Origin, levels = origin_levels), factor(cc_ischia_10, levels = ecotype_levels))

  summary_list[[gene_name]] <- summary_df

  plot_normal <- summary_df %>%
    filter(Origin == "Normal Pancreas") %>%
    mutate(cc_ischia_10 = factor(cc_ischia_10, levels = ecotype_levels)) %>%
    make_ecotype_barplot(gene_name = gene_name, origin_name = "Normal Pancreas", palette_ecotype = palette_ecotype)

  plot_primary <- summary_df %>%
    filter(Origin == "Pancreas") %>%
    mutate(cc_ischia_10 = factor(cc_ischia_10, levels = ecotype_levels)) %>%
    make_ecotype_barplot(gene_name = gene_name, origin_name = "Primary Pancreas", palette_ecotype = palette_ecotype)

  save_two_page_pdf(
    plot_normal = plot_normal,
    plot_primary = plot_primary,
    filename = paste0("", gene_name, "_detection_by_ecotype.pdf"),
    out_dir = barplot_fig_dir,
    width = 6.6,
    height = 4.6
  )
}

summary_all <- bind_rows(summary_list)
write_tsv_simple(summary_all, "kcn_detection_by_ecotype_summary.tsv", tab_dir)


message("ST PDAC KCN25 ecotype barplots complete.")



