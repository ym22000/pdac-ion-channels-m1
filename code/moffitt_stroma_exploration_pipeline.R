###############################################################################
# Moffitt stroma exploration pipeline
#
# This script gives a first overview of the stromal-focused Moffitt atlas.
# It summarizes the metadata, QC metrics, stromal subtypes, and marker patterns
# before running more specific KCN analyses.
###############################################################################
required_packages <- c(
  "Seurat",
  "dplyr",
  "ggplot2",
  "patchwork",
  "tibble",
  "scales",
  "purrr"
)

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0) {
  stop(
    "Please install the following R packages before running the script: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(tibble)
  library(scales)
  library(purrr)
})

###############################################################################
# 1. Helpers
###############################################################################

get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE))
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

find_project_root <- function() {
  start_points <- unique(c(
    normalizePath(getwd(), winslash = "/", mustWork = TRUE),
    dirname(get_script_path())
  ))

  for (start_dir in start_points) {
    current_dir <- start_dir
    for (i in seq_len(6)) {
      dataset_candidate <- file.path(current_dir, "scripts", "dataset", "Stroma_Subset2021.rds")
      if (file.exists(dataset_candidate)) {
        return(normalizePath(current_dir, winslash = "/", mustWork = TRUE))
      }
      parent_dir <- dirname(current_dir)
      if (identical(parent_dir, current_dir)) {
        break
      }
      current_dir <- parent_dir
    }
  }

  stop("Project root not found. Expected scripts/dataset/Stroma_Subset2021.rds.")
}

save_tsv <- function(x, file_path) {
  write.table(
    x,
    file = file_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )
}

save_plot <- function(plot_object, file_path, width = 10, height = 7, dpi = 320) {
  ggsave(
    filename = file_path,
    plot = plot_object,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
}

clean_embedding_theme <- function() {
  theme_classic(base_size = 13) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 11),
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.border = element_blank()
    )
}

clean_classic_theme <- function() {
  theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.title = element_text(face = "bold"),
      panel.grid = element_blank()
    )
}

###############################################################################
# 2. Paths and colors
###############################################################################

project_root <- find_project_root()
analysis_root <- file.path(project_root, "M1_Bioinformatics_Ion_Channel", "code")
project_data_dir <- file.path(project_root, "M1_Bioinformatics_Ion_Channel", "data")
local_input_rds <- file.path(project_data_dir, "Stroma_Subset2021.rds")
fallback_input_rds <- file.path(project_root, "scripts", "dataset", "Stroma_Subset2021.rds")
input_rds <- if (file.exists(local_input_rds)) local_input_rds else fallback_input_rds

output_root <- file.path(analysis_root, "moffitt_stromal_exploration_outputs")
figure_dir <- file.path(output_root, "figures")
table_dir <- file.path(output_root, "tables")

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

state_colors <- c(
  "Activated" = "#FF686B",
  "Normal" = "#70D6FF"
)

subtype_order <- c("qPSC", "smPSC", "myCAF", "csCAF", "iCAF", "IL11.CAF", "Myocyte", "Schwann")

celltype_colors <- c(
  "qPSC" = "#70D6FF",
  "smPSC" = "#0096C7",
  "myCAF" = "#FF686B",
  "csCAF" = "#FF9770",
  "iCAF" = "#FF5400",
  "IL11.CAF" = "#E9FF70",
  "Myocyte" = "#60D394",
  "Schwann" = "#8338EC"
)

dataset_colors <- c(
  "Peng" = "#0096C7",
  "Powers" = "#EB6424",
  "Qadir" = "#8338EC"
)

feature_colors <- c("grey94", "#FF5400")

###############################################################################
# 3. Load the Seurat object
###############################################################################

if (!file.exists(input_rds)) {
  stop("Input .rds file not found: ", input_rds)
}

stroma <- readRDS(input_rds)
DefaultAssay(stroma) <- "RNA"

required_meta <- c(
  "Patient",
  "Dataset",
  "Condition",
  "CellType2",
  "nFeature_RNA",
  "nCount_RNA",
  "percent.mt",
  "Moffitt.F5_ActivatedStroma.top25",
  "Moffitt.F13_NormalStroma.top25"
)

missing_meta <- setdiff(required_meta, colnames(stroma@meta.data))
if (length(missing_meta) > 0) {
  stop("Missing metadata columns: ", paste(missing_meta, collapse = ", "))
}

if (!"umap" %in% Reductions(stroma)) {
  stop("The Seurat object does not contain a UMAP reduction.")
}

###############################################################################
# 4. Build the stromal state annotation
###############################################################################

stroma$StromaScore <- stroma$Moffitt.F5_ActivatedStroma.top25 -
  stroma$Moffitt.F13_NormalStroma.top25

stroma$StromaState <- ifelse(
  stroma$Moffitt.F5_ActivatedStroma.top25 >
    stroma$Moffitt.F13_NormalStroma.top25,
  "Activated",
  "Normal"
)

stroma$StromaState <- factor(stroma$StromaState, levels = c("Normal", "Activated"))
stroma$CellType2 <- factor(as.character(stroma$CellType2), levels = subtype_order)

meta <- stroma@meta.data %>%
  rownames_to_column("cell_id")

###############################################################################
# 5. Metadata inspection tables
###############################################################################

metadata_columns <- tibble(metadata_column = colnames(stroma@meta.data))

dataset_summary <- tibble(
  metric = c(
    "n_cells",
    "n_genes",
    "n_samples",
    "n_datasets",
    "n_celltype2_labels",
    "n_reductions",
    "pct_normal_state",
    "pct_activated_state"
  ),
  value = c(
    ncol(stroma),
    nrow(stroma),
    dplyr::n_distinct(meta$Patient),
    dplyr::n_distinct(meta$Dataset),
    dplyr::n_distinct(meta$CellType2),
    length(Reductions(stroma)),
    round(mean(meta$StromaState == "Normal") * 100, 2),
    round(mean(meta$StromaState == "Activated") * 100, 2)
  )
)

state_counts <- meta %>%
  count(StromaState, name = "n_cells") %>%
  mutate(percent_cells = round(100 * n_cells / sum(n_cells), 2))

celltype_by_state <- meta %>%
  count(CellType2, StromaState, name = "n_cells") %>%
  group_by(CellType2) %>%
  mutate(percent_within_celltype = round(100 * n_cells / sum(n_cells), 2)) %>%
  ungroup()

cells_per_sample <- meta %>%
  count(Patient, Dataset, Condition, name = "n_cells") %>%
  mutate(percent_cells = 100 * n_cells / sum(n_cells)) %>%
  arrange(desc(percent_cells)) %>%
  mutate(Patient = factor(Patient, levels = rev(Patient)))

patient_celltype_composition <- meta %>%
  count(Patient, Dataset, CellType2, name = "n_cells") %>%
  group_by(Patient, Dataset) %>%
  mutate(fraction = n_cells / sum(n_cells)) %>%
  ungroup()

qc_by_state <- meta %>%
  group_by(StromaState) %>%
  summarise(
    n_cells = n(),
    median_nFeature_RNA = median(nFeature_RNA, na.rm = TRUE),
    median_nCount_RNA = median(nCount_RNA, na.rm = TRUE),
    median_percent_mt = median(percent.mt, na.rm = TRUE),
    .groups = "drop"
  )

qc_by_sample <- meta %>%
  group_by(Patient, Dataset, Condition) %>%
  summarise(
    n_cells = n(),
    median_nFeature_RNA = median(nFeature_RNA, na.rm = TRUE),
    median_nCount_RNA = median(nCount_RNA, na.rm = TRUE),
    median_percent_mt = median(percent.mt, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_cells))

qc_by_celltype <- meta %>%
  group_by(CellType2) %>%
  summarise(
    n_cells = n(),
    median_nFeature_RNA = median(nFeature_RNA, na.rm = TRUE),
    median_nCount_RNA = median(nCount_RNA, na.rm = TRUE),
    median_percent_mt = median(percent.mt, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_cells))

save_tsv(metadata_columns, file.path(table_dir, "metadata_columns.tsv"))
save_tsv(dataset_summary, file.path(table_dir, "dataset_summary.tsv"))
save_tsv(state_counts, file.path(table_dir, "state_counts.tsv"))
save_tsv(celltype_by_state, file.path(table_dir, "celltype_by_state.tsv"))
save_tsv(cells_per_sample, file.path(table_dir, "cells_per_sample.tsv"))
save_tsv(patient_celltype_composition, file.path(table_dir, "patient_celltype_composition.tsv"))
save_tsv(qc_by_state, file.path(table_dir, "qc_by_state.tsv"))
save_tsv(qc_by_sample, file.path(table_dir, "qc_by_sample.tsv"))
save_tsv(qc_by_celltype, file.path(table_dir, "qc_by_celltype.tsv"))

summary_lines <- c(
  "Moffitt stroma exploration pipeline",
  paste("Run date:", as.character(Sys.Date())),
  paste("Input object:", input_rds),
  "",
  "Main counts",
  paste("Cells:", ncol(stroma)),
  paste("Genes:", nrow(stroma)),
  paste("Samples:", dplyr::n_distinct(meta$Patient)),
  paste("Datasets:", dplyr::n_distinct(meta$Dataset)),
  paste("Available reductions:", paste(Reductions(stroma), collapse = ", ")),
  "",
  "Stromal state rule",
  "Activated if Moffitt.F5_ActivatedStroma.top25 > Moffitt.F13_NormalStroma.top25",
  "Normal otherwise"
)
writeLines(summary_lines, con = file.path(output_root, "metadata_summary.txt"))

###############################################################################
# 6. Curated literature markers used for subtype visualization
###############################################################################

literature_marker_table <- tibble::tribble(
  ~CellType2,   ~gene,     ~literature_note,
  "qPSC",       "RGS5",    "Perivascular / quiescent stellate-like marker",
  "qPSC",       "COX4I2",  "Perivascular hypoxia-associated marker",
  "smPSC",      "ACTA2",   "Smooth-muscle-like stellate marker",
  "smPSC",      "MYH11",   "Contractile smooth-muscle-like marker",
  "myCAF",      "POSTN",   "Canonical myCAF extracellular matrix marker",
  "myCAF",      "COL11A1", "Canonical matrix-remodeling myCAF marker",
  "csCAF",      "C3",      "Complement-secreting CAF marker",
  "csCAF",      "C7",      "Complement-secreting CAF marker",
  "iCAF",       "IL6",     "Inflammatory CAF cytokine marker",
  "iCAF",       "HAS1",    "Inflammatory / secretory CAF marker",
  "IL11.CAF",   "IL11",    "IL11-high activated CAF marker",
  "IL11.CAF",   "INHBA",   "Activated IL11-associated CAF marker",
  "Myocyte",    "ACTG2",   "Muscle-like contractile marker",
  "Myocyte",    "MYLK",    "Muscle contraction marker",
  "Schwann",    "S100B",   "Canonical Schwann cell marker",
  "Schwann",    "MPZ",     "Myelinating Schwann cell marker"
) %>%
  filter(!CellType2 %in% c("IL11.CAF", "Myocyte")) %>%
  mutate(CellType2 = factor(CellType2, levels = subtype_order))

marker_genes_used <- literature_marker_table %>%
  filter(gene %in% rownames(stroma))

if (nrow(marker_genes_used) < nrow(literature_marker_table)) {
  missing_marker_genes <- setdiff(literature_marker_table$gene, marker_genes_used$gene)
  stop("Missing curated literature markers in the Seurat object: ", paste(missing_marker_genes, collapse = ", "))
}

save_tsv(marker_genes_used, file.path(table_dir, "marker_genes_used.tsv"))

###############################################################################
# 7. QC plots by stromal subtype
###############################################################################

p_qc <- VlnPlot(
  object = stroma,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "CellType2",
  cols = celltype_colors,
  pt.size = 0,
  ncol = 3
) &
  clean_classic_theme() &
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(face = "bold", angle = 45, hjust = 1),
    legend.position = "none"
  )

save_plot(
  plot_object = p_qc,
  file_path = file.path(figure_dir, "Figure_01_QC_metrics_by_subtype.png"),
  width = 16,
  height = 5.8
)

###############################################################################
# 8. UMAP plots
###############################################################################
# We reuse the UMAP coordinates already stored inside the Seurat object.
# This means we only change the coloring of cells here, not the embedding itself.

umap_df <- Embeddings(stroma, reduction = "umap") %>%
  as.data.frame() %>%
  rownames_to_column("cell_id") %>%
  left_join(
    meta %>%
      select(
        cell_id,
        CellType2,
        Moffitt.F5_ActivatedStroma.top25,
        Moffitt.F13_NormalStroma.top25
      ),
    by = "cell_id"
  )

label_positions <- umap_df %>%
  group_by(CellType2) %>%
  summarise(
    UMAP_1 = median(UMAP_1, na.rm = TRUE),
    UMAP_2 = median(UMAP_2, na.rm = TRUE),
    .groups = "drop"
  )

make_signature_umap <- function(plot_df, score_column, panel_title, high_color) {
  ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = .data[[score_column]])) +
    geom_point(size = 0.32, alpha = 0.95) +
    geom_text(
      data = label_positions,
      aes(x = UMAP_1, y = UMAP_2, label = CellType2),
      inherit.aes = FALSE,
      color = "black",
      fontface = "bold",
      size = 4.8
    ) +
    scale_color_gradient(
      low = "grey92",
      high = high_color,
      name = "Relative signature expression",
      guide = guide_colorbar(
        title.position = "top",
        direction = "horizontal",
        barwidth = grid::unit(4.2, "cm"),
        barheight = grid::unit(0.25, "cm")
      )
    ) +
    labs(title = panel_title) +
    clean_embedding_theme() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      legend.position = "bottom",
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 10, face = "bold")
    )
}

p_umap <- make_signature_umap(
  plot_df = umap_df,
  score_column = "Moffitt.F5_ActivatedStroma.top25",
  panel_title = "Activated",
  high_color = "#8C510A"
) +
  make_signature_umap(
    plot_df = umap_df,
    score_column = "Moffitt.F13_NormalStroma.top25",
    panel_title = "Normal",
    high_color = "#70D6FF"
  )

save_plot(
  plot_object = p_umap,
  file_path = file.path(figure_dir, "Figure_02_UMAP_stromal_signatures_doubleplot.png"),
  width = 12,
  height = 6.8
)

# This figure shows all CellType2 labels on the shared UMAP embedding.
p_umap_celltype <- DimPlot(
  object = stroma,
  reduction = "umap",
  group.by = "CellType2",
  cols = celltype_colors,
  pt.size = 0.35,
  label = TRUE,
  repel = TRUE,
  raster = FALSE
) +
  labs(title = "UMAP of stromal subtypes") +
  clean_embedding_theme() +
  theme(legend.position = "right")

save_plot(
  plot_object = p_umap_celltype,
  file_path = file.path(figure_dir, "Figure_03_UMAP_all_stromal_subtypes.png"),
  width = 10,
  height = 8
)

###############################################################################
# 9. Multi-panel FeaturePlot for curated literature markers
###############################################################################

featureplot_genes <- marker_genes_used %>%
  arrange(CellType2) %>%
  pull(gene)

featureplot_labels <- marker_genes_used %>%
  arrange(CellType2) %>%
  mutate(panel_title = paste0(gene, " (", CellType2, ")")) %>%
  pull(panel_title)

featureplot_list <- FeaturePlot(
  object = stroma,
  features = featureplot_genes,
  reduction = "umap",
  cols = c("grey93", "#FF5400"),
  pt.size = 0.28,
  order = TRUE,
  ncol = 4,
  combine = FALSE
)

featureplot_list <- purrr::map2(
  .x = featureplot_list,
  .y = featureplot_labels,
  .f = function(plot_object, panel_title) {
    plot_object +
      ggtitle(panel_title) +
      clean_embedding_theme() +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
        legend.position = "right"
      )
  }
)

p_feature_deg <- wrap_plots(featureplot_list, ncol = 4) +
  plot_annotation(
    title = "Literature markers by stromal subtype"
  )

save_plot(
  plot_object = p_feature_deg,
  file_path = file.path(figure_dir, "Figure_04_FeaturePlot_literature_markers_by_subtype.png"),
  width = 16,
  height = 14
)

###############################################################################
# 10. Patient composition for stromal cell subtypes
###############################################################################

# This barplot displays stacked CellType2 proportions for each patient.
p_patient_composition <- ggplot(
  patient_celltype_composition %>%
    mutate(
      Dataset = factor(Dataset, levels = c("Peng", "Powers", "Qadir")),
      Patient = factor(Patient, levels = unique(cells_per_sample$Patient)),
      CellType2 = factor(CellType2, levels = names(celltype_colors))
    ),
  aes(x = Patient, y = fraction, fill = CellType2)
) +
  geom_col(width = 0.88, color = "white", linewidth = 0.15) +
  facet_grid(~ Dataset, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = celltype_colors, drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
  labs(
    title = "Patient composition of stromal cell subtypes",
    x = "Sample",
    y = "Cell proportion",
    fill = "Cell subtype"
  ) +
  clean_classic_theme() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "bold", size = 8),
    legend.position = "top",
    panel.spacing.x = grid::unit(0.5, "lines")
  )

save_plot(
  plot_object = p_patient_composition,
  file_path = file.path(figure_dir, "Figure_05_Patient_composition_cell_subtypes.png"),
  width = 14,
  height = 7
)

###############################################################################
# 11. Session info
###############################################################################

message("Exploration finished.")
message("Outputs saved in: ", output_root)

