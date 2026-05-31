#!/usr/bin/env Rscript

###############################################################################
# USER scRNA-seq validation of Moffitt stromal states
#
# This script projects the Moffitt stromal state framework onto the local USER
# PDAC single-cell dataset.
#
# The main goal is to rebuild a fibroblast-focused view of USER, score the
# Moffitt stromal programs, and recover a cleaner four-state fibroblast panel
# (`myCAF`, `iCAF`, `csCAF`, `PSC`).
###############################################################################
suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
  library(readxl)
  library(scales)
})

options(stringsAsFactors = FALSE)
set.seed(1234)

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
input_dir <- file.path(script_dir, "inputs")
resource_dir <- file.path(script_dir, "resources")
fig_dir <- file.path(script_dir, "figures")
tab_dir <- file.path(script_dir, "tables")
rds_dir <- file.path(script_dir, "rds")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

matrix_path <- file.path(input_dir, "gene_sorted-naivedata_scp.mtx")
genes_path <- file.path(input_dir, "naivedata_scp.genes.csv")
barcodes_path <- file.path(input_dir, "naivedata_scp.barcodes.csv")
annot_path <- file.path(input_dir, "combinenaivedata-reprocessed-clean-detailed-annotations.tsv")
umap_path <- file.path(input_dir, "combinenaivedata-reprocessed-clean-detailed-UMAP.tsv")
moffitt_markers_xlsx <- file.path(resource_dir, "Moffitt_SuppData1_CellType2_Markers.xlsx")
marker_genes_of_interest <- c("CXCL14", "ACTA2", "C7", "FAP", "IL11", "CXCL12", "PDGFRA")
marker_color_high <- "#03045E"
marker_color_low <- "#FCF300"

save_plot <- function(plot_obj, filename, width = 10, height = 8) {
  pdf_path <- file.path(fig_dir, filename)
  png_path <- sub("\\.pdf$", ".png", pdf_path, ignore.case = TRUE)
  ggplot2::ggsave(pdf_path, plot = plot_obj, width = width, height = height, units = "in", dpi = 300, bg = "white")
  ggplot2::ggsave(png_path, plot = plot_obj, width = width, height = height, units = "in", dpi = 300, bg = "white")
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

is_gzip_file <- function(path) {
  con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)
  magic <- readBin(con, what = "raw", n = 2)
  identical(as.integer(magic), c(31L, 139L))
}

open_text_connection <- function(path) {
  if (is_gzip_file(path)) gzfile(path, open = "rt") else file(path, open = "rt")
}

make_palette <- function(labels) {
  base_palette <- c(
    "Tumor" = "#7400b8",
    "Fibroblast" = "#72efdd",
    "Immune" = "#ffd6a5",
    "Ductal" = "#fdffb6",
    "Endothelial" = "#bde0fe",
    "Atypical_Ductal" = "#ffc6ff",
    "Endocrine" = "#caffbf",
    "SmoothMuscle" = "#fec5bb",
    "NascentEndothelial" = "#a0c4ff",
    "Schwann" = "#d0f4de",
    "Acinar" = "#ce4257"
  )

  missing_labels <- setdiff(labels, names(base_palette))
  if (length(missing_labels) > 0) {
    extra_cols <- grDevices::hcl.colors(length(missing_labels), palette = "Pastel 1")
    names(extra_cols) <- missing_labels
    base_palette <- c(base_palette, extra_cols)
  }

  base_palette[labels]
}

get_signature_score <- function(expr_mat, genes_use) {
  genes_use <- intersect(genes_use, rownames(expr_mat))
  if (length(genes_use) == 0) {
    return(rep(NA_real_, ncol(expr_mat)))
  }
  Matrix::colMeans(expr_mat[genes_use, , drop = FALSE])
}

get_global_log_normalized_markers <- function(count_mat, genes_use) {
  genes_use <- unique(toupper(genes_use))
  lib_size <- Matrix::colSums(count_mat)
  lib_size[lib_size == 0] <- 1

  out_list <- lapply(genes_use, function(gene_i) {
    expr_i <- if (gene_i %in% rownames(count_mat)) {
      as.numeric(count_mat[gene_i, , drop = TRUE])
    } else {
      rep(0, ncol(count_mat))
    }
    data.frame(
      cell = colnames(count_mat),
      gene = gene_i,
      expression = log1p(10000 * expr_i / lib_size),
      stringsAsFactors = FALSE
    )
  })

  bind_rows(out_list)
}

get_seurat_marker_long <- function(seurat_obj, genes_use, assay = "RNA", layer_name = "data") {
  genes_use <- unique(toupper(genes_use))
  expr_mat <- get_assay_matrix(seurat_obj, assay = assay, layer_name = layer_name)

  out_list <- lapply(genes_use, function(gene_i) {
    expr_i <- if (gene_i %in% rownames(expr_mat)) {
      as.numeric(expr_mat[gene_i, , drop = TRUE])
    } else {
      rep(0, ncol(expr_mat))
    }
    data.frame(
      cell = colnames(expr_mat),
      gene = gene_i,
      expression = expr_i,
      stringsAsFactors = FALSE
    )
  })

  bind_rows(out_list)
}

summarize_marker_stats <- function(marker_df, group_var) {
  marker_df %>%
    group_by(.data[[group_var]], gene) %>%
    summarise(
      n_cells = n(),
      pct_detected = round(100 * mean(expression > 0, na.rm = TRUE), 2),
      mean_expression = mean(expression, na.rm = TRUE),
      median_expression = stats::median(expression, na.rm = TRUE),
      mean_detected_only = ifelse(any(expression > 0, na.rm = TRUE), mean(expression[expression > 0], na.rm = TRUE), NA_real_),
      .groups = "drop"
    ) %>%
    rename(group = all_of(group_var))
}

zscore_vec <- function(x) {
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) {
    return(rep(0, length(x)))
  }
  (x - mean(x, na.rm = TRUE)) / s
}

get_assay_matrix <- function(seurat_obj, assay = "RNA", layer_name = "data") {
  out <- tryCatch(
    Seurat::GetAssayData(seurat_obj, assay = assay, layer = layer_name),
    error = function(e) Seurat::GetAssayData(seurat_obj, assay = assay, slot = layer_name)
  )
  out
}

extract_umap_df <- function(seurat_obj, reduction_name = "umap", prefix = "umap") {
  emb <- as.data.frame(Embeddings(seurat_obj, reduction_name))
  if (ncol(emb) < 2) {
    stop("UMAP embedding does not contain two coordinates.")
  }
  out <- emb[, 1:2, drop = FALSE]
  colnames(out) <- c(paste0(prefix, "_1"), paste0(prefix, "_2"))
  tibble::rownames_to_column(out, "cell")
}

theme_stage <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "grey95", colour = "grey70"),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    )
}

theme_umap_clean <- function(base_size = 11) {
  theme_stage(base_size) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.key.height = grid::unit(0.45, "cm"),
      legend.key.width = grid::unit(0.45, "cm")
    )
}

message("Loading Moffitt stromal signatures...")
if (!file.exists(moffitt_markers_xlsx)) {
  stop("Moffitt marker file not found: ", moffitt_markers_xlsx)
}

moffitt_markers <- readxl::read_excel(moffitt_markers_xlsx, sheet = "Stroma") %>%
  mutate(
    cluster = as.character(cluster),
    gene = toupper(as.character(gene))
  ) %>%
  filter(
    !is.na(cluster),
    !is.na(gene),
    gene != "",
    p_val_adj < 0.05,
    avg_log2FC > 0.25,
    pct.1 > 0.10
  )

moffitt_signatures <- moffitt_markers %>%
  arrange(cluster, desc(avg_log2FC), desc(pct.1)) %>%
  group_by(cluster) %>%
  slice_head(n = 30) %>%
  summarise(
    signature_genes = list(unique(gene)),
    n_genes = length(unique(gene)),
    .groups = "drop"
  ) %>%
  rename(moffitt_label = cluster)

write_tsv_simple(
  moffitt_signatures %>% mutate(signature_genes = vapply(signature_genes, paste, collapse = ", ", character(1))),
  "moffitt_stromal_signatures_top30.tsv"
)

message("Loading USER genes...")
con_genes <- open_text_connection(genes_path)
genes <- readLines(con_genes, warn = FALSE)
close(con_genes)
genes <- trimws(genes)
genes <- genes[genes != ""]
genes <- make.unique(toupper(genes))

message("Loading USER barcodes...")
con_barcodes <- open_text_connection(barcodes_path)
barcode_lines <- readLines(con_barcodes, warn = FALSE)
close(con_barcodes)
barcode_split <- strsplit(barcode_lines, ",", fixed = TRUE)
barcode_split <- barcode_split[lengths(barcode_split) >= 2]
barcodes <- vapply(barcode_split[-1], function(x) trimws(tail(x, 1)), character(1))

message("Loading USER annotations and original UMAP coordinates...")
annot_tbl <- read.delim(annot_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
if (annot_tbl$NAME[1] == "TYPE") annot_tbl <- annot_tbl[-1, , drop = FALSE]

umap_tbl <- read.delim(umap_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
if (umap_tbl$NAME[1] == "TYPE") umap_tbl <- umap_tbl[-1, , drop = FALSE]
umap_tbl$X <- as.numeric(umap_tbl$X)
umap_tbl$Y <- as.numeric(umap_tbl$Y)

message("Loading USER sparse expression matrix...")
expr_mat <- readMM(matrix_path)
expr_mat <- as(expr_mat, "CsparseMatrix")

if (nrow(expr_mat) != length(genes)) stop("Gene count mismatch between matrix and gene file.")
if (ncol(expr_mat) != length(barcodes)) stop("Barcode count mismatch between matrix and barcode file.")

rownames(expr_mat) <- genes
colnames(expr_mat) <- barcodes

if (!identical(colnames(expr_mat), annot_tbl$NAME)) stop("Cell order mismatch between matrix and annotation table.")
if (!identical(colnames(expr_mat), umap_tbl$NAME)) stop("Cell order mismatch between matrix and UMAP table.")

user_full_meta <- annot_tbl %>%
  left_join(umap_tbl, by = "NAME")

message("Computing global QC metrics from the count matrix...")
mt_genes <- grep("^MT-", rownames(expr_mat), value = TRUE)
total_counts <- Matrix::colSums(expr_mat)
detected_features <- Matrix::colSums(expr_mat > 0)
mt_counts <- if (length(mt_genes) > 0) {
  Matrix::colSums(expr_mat[mt_genes, , drop = FALSE])
} else {
  rep(0, ncol(expr_mat))
}

user_full_meta <- user_full_meta %>%
  mutate(
    nCount_RNA = as.numeric(total_counts[NAME]),
    nFeature_RNA = as.numeric(detected_features[NAME]),
    percent.mt = 100 * as.numeric(mt_counts[NAME]) / pmax(as.numeric(total_counts[NAME]), 1)
  )

cell_subset_counts <- user_full_meta %>%
  count(cell_subsets, name = "n_cells") %>%
  mutate(pct_cells = round(100 * n_cells / sum(n_cells), 2)) %>%
  arrange(desc(n_cells))

dataset_overview <- data.frame(
  metric = c("n_total_cells", "n_total_genes", "n_fibroblast_cells", "n_unique_patients"),
  value = c(
    ncol(expr_mat),
    nrow(expr_mat),
    sum(user_full_meta$cell_subsets == "Fibroblast"),
    length(unique(user_full_meta$pid))
  ),
  stringsAsFactors = FALSE
)

message("Computing selected marker summaries on the global USER dataset...")
global_marker_long <- get_global_log_normalized_markers(expr_mat, marker_genes_of_interest)
global_marker_wide <- global_marker_long %>%
  pivot_wider(names_from = gene, values_from = expression)

user_full_meta_out <- user_full_meta %>%
  left_join(global_marker_wide, by = c("NAME" = "cell"))

global_marker_stats <- user_full_meta_out %>%
  select(NAME, cell_subsets, all_of(marker_genes_of_interest)) %>%
  pivot_longer(cols = all_of(marker_genes_of_interest), names_to = "gene", values_to = "expression") %>%
  summarize_marker_stats("cell_subsets")

write_tsv_simple(user_full_meta_out, "USER_full_metadata.tsv")
write_tsv_simple(cell_subset_counts, "USER_cell_subset_counts.tsv")
write_tsv_simple(dataset_overview, "USER_dataset_overview.tsv")
write_tsv_simple(global_marker_stats, "USER_global_marker_stats.tsv")

global_subset_levels <- unique(user_full_meta$cell_subsets)
global_subset_levels <- cell_subset_counts$cell_subsets
user_full_meta$cell_subsets <- factor(user_full_meta$cell_subsets, levels = global_subset_levels)
global_colors <- make_palette(global_subset_levels)

p_global <- ggplot(user_full_meta, aes(x = X, y = Y, color = cell_subsets)) +
  geom_point(size = 0.22, alpha = 0.8) +
  scale_color_manual(values = global_colors) +
  labs(
    title = "USER UMAP",
    color = "Cell subset"
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme_umap_clean(11)

save_plot(p_global, "Figure_01_USER_global_umap_cell_subsets.pdf", width = 11, height = 8.5)

patient_composition_df <- user_full_meta %>%
  count(pid, cell_subsets, name = "n_cells") %>%
  group_by(pid) %>%
  mutate(cell_proportion = n_cells / sum(n_cells)) %>%
  ungroup()

patient_order <- patient_composition_df %>%
  group_by(pid) %>%
  summarise(tumor_prop = sum(cell_proportion[cell_subsets == "Tumor"], na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(tumor_prop), pid) %>%
  pull(pid)

patient_composition_df$pid <- factor(patient_composition_df$pid, levels = patient_order)

p_patient_composition <- ggplot(patient_composition_df, aes(x = pid, y = cell_proportion, fill = cell_subsets)) +
  geom_col(width = 0.85, color = "white", linewidth = 0.15) +
  scale_fill_manual(values = global_colors, drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "USER cell subset composition by patient",
    x = "Patient",
    y = "Cell proportion",
    fill = "Cell subset"
  ) +
  theme_stage(11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.key.height = grid::unit(0.45, "cm"),
    legend.key.width = grid::unit(0.45, "cm")
  )

save_plot(p_patient_composition, "Figure_06_USER_patient_cell_subset_composition.pdf", width = 11.5, height = 6.8)

qc_long <- user_full_meta %>%
  select(NAME, cell_subsets, nFeature_RNA, nCount_RNA, percent.mt) %>%
  pivot_longer(
    cols = c(nFeature_RNA, nCount_RNA, percent.mt),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = factor(metric, levels = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
  )

set.seed(1234)
qc_points <- qc_long %>%
  group_by(metric, cell_subsets) %>%
  group_modify(~ dplyr::slice_sample(.x, n = min(600, nrow(.x)))) %>%
  ungroup()

p_qc <- ggplot(qc_long, aes(x = cell_subsets, y = value, fill = cell_subsets)) +
  geom_violin(scale = "width", linewidth = 0.2, color = "grey35", trim = TRUE) +
  geom_boxplot(width = 0.14, outlier.shape = NA, linewidth = 0.22, fill = "white") +
  geom_point(
    data = qc_points,
    aes(color = cell_subsets),
    size = 0.18,
    alpha = 0.18,
    stroke = 0,
    position = position_jitter(width = 0.18, height = 0)
  ) +
  facet_wrap(~ metric, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = global_colors, drop = FALSE) +
  scale_color_manual(values = global_colors, drop = FALSE) +
  labs(
    title = "USER QC metrics by cell subset",
    x = "Cell subset",
    y = NULL
  ) +
  theme_stage(10.5) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )

save_plot(p_qc, "Figure_07_USER_global_qc_by_cell_subset.pdf", width = 14, height = 6.5)

global_marker_plot_df <- user_full_meta_out %>%
  select(NAME, X, Y, all_of(marker_genes_of_interest)) %>%
  pivot_longer(cols = all_of(marker_genes_of_interest), names_to = "gene", values_to = "expression")

p_global_markers <- ggplot(global_marker_plot_df, aes(x = X, y = Y, color = expression)) +
  geom_point(size = 0.20, alpha = 0.85) +
  facet_wrap(~ gene, ncol = 4) +
  scale_color_gradient(low = marker_color_low, high = marker_color_high) +
  labs(
    title = "USER global UMAP: selected fibroblast markers",
    subtitle = "Published USER coordinates colored by log-normalized expression",
    x = "USER UMAP 1",
    y = "USER UMAP 2",
    color = "Expression"
  ) +
  theme_stage(11)

save_plot(p_global_markers, "Figure_04_USER_global_marker_umap.pdf", width = 14, height = 7.6)

message("Subsetting USER fibroblasts...")
fibro_cells <- annot_tbl$NAME[annot_tbl$cell_subsets == "Fibroblast"]
fibro_idx <- which(barcodes %in% fibro_cells)
expr_fibro <- expr_mat[, fibro_idx, drop = FALSE]

fibro_meta <- annot_tbl %>%
  filter(NAME %in% colnames(expr_fibro)) %>%
  slice(match(colnames(expr_fibro), NAME))

fibro_umap_orig <- umap_tbl %>%
  filter(NAME %in% colnames(expr_fibro)) %>%
  slice(match(colnames(expr_fibro), NAME))

if (!identical(fibro_meta$NAME, colnames(expr_fibro))) stop("Fibro metadata and expression matrix are not aligned.")
if (!identical(fibro_umap_orig$NAME, colnames(expr_fibro))) stop("Fibro UMAP table and expression matrix are not aligned.")

message("Building fibroblast Seurat object...")
fibro <- CreateSeuratObject(counts = expr_fibro, meta.data = fibro_meta)
fibro$USER_UMAP_1 <- fibro_umap_orig$X
fibro$USER_UMAP_2 <- fibro_umap_orig$Y

message("Running fibroblast preprocessing...")
fibro <- NormalizeData(fibro, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
fibro <- FindVariableFeatures(fibro, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
fibro <- ScaleData(fibro, features = VariableFeatures(fibro), verbose = FALSE)
fibro <- RunPCA(fibro, features = VariableFeatures(fibro), npcs = 30, verbose = FALSE)
fibro <- FindNeighbors(fibro, dims = 1:20, verbose = FALSE)
fibro <- FindClusters(fibro, resolution = 0.6, algorithm = 1, verbose = FALSE)
fibro <- RunUMAP(fibro, dims = 1:20, n.neighbors = 20, min.dist = 0.20, verbose = FALSE)

message("Computing fibroblast cluster markers...")
Idents(fibro) <- fibro$seurat_clusters
deg_tbl <- FindAllMarkers(
  object = fibro,
  only.pos = TRUE,
  min.pct = 0.10,
  logfc.threshold = 0.25,
  assay = "RNA",
  verbose = FALSE
) %>%
  arrange(cluster, desc(avg_log2FC), p_val_adj)

write_tsv_simple(deg_tbl, "USER_fibro_clusters_all_DEGs.tsv")

deg_top5 <- deg_tbl %>%
  group_by(cluster) %>%
  slice_head(n = 5) %>%
  summarise(top5_degs = paste(gene, collapse = ", "), .groups = "drop")

deg_top100 <- deg_tbl %>%
  group_by(cluster) %>%
  slice_head(n = 100) %>%
  ungroup() %>%
  mutate(cluster = as.character(cluster))

message("Projecting Moffitt stromal signatures...")
norm_mat <- get_assay_matrix(fibro, assay = "RNA", layer_name = "data")
score_cols <- character(0)
detected_signature_rows <- list()

for (i in seq_len(nrow(moffitt_signatures))) {
  label_i <- moffitt_signatures$moffitt_label[i]
  genes_i <- moffitt_signatures$signature_genes[[i]]
  genes_present_i <- intersect(genes_i, rownames(norm_mat))
  score_name_i <- paste0("score_", gsub("[^A-Za-z0-9]", "_", label_i))
  fibro[[score_name_i]] <- get_signature_score(norm_mat, genes_i)
  score_cols <- c(score_cols, score_name_i)
  detected_signature_rows[[i]] <- data.frame(
    moffitt_label = label_i,
    n_genes_in_signature = length(genes_i),
    n_genes_detected_in_USER = length(genes_present_i),
    detected_genes = paste(genes_present_i, collapse = ", "),
    stringsAsFactors = FALSE
  )
}

detected_signature_tbl <- bind_rows(detected_signature_rows)
write_tsv_simple(detected_signature_tbl, "USER_detected_Moffitt_signature_genes.tsv")

message("Annotating fibroblast clusters from Moffitt projection...")
cluster_score_summary <- fibro@meta.data %>%
  rownames_to_column("cell") %>%
  group_by(seurat_clusters) %>%
  summarise(across(all_of(score_cols), \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
  pivot_longer(cols = all_of(score_cols), names_to = "score_col", values_to = "mean_signature_score") %>%
  mutate(moffitt_label = gsub("^score_", "", score_col)) %>%
  mutate(moffitt_label = gsub("_", ".", moffitt_label))

cluster_overlap_summary <- bind_rows(lapply(seq_len(nrow(moffitt_signatures)), function(i) {
  label_i <- moffitt_signatures$moffitt_label[i]
  genes_i <- moffitt_signatures$signature_genes[[i]]
  deg_top100 %>%
    group_by(cluster) %>%
    summarise(
      overlap_n = sum(toupper(gene) %in% genes_i),
      overlap_genes = paste(unique(toupper(gene[toupper(gene) %in% genes_i])), collapse = ", "),
      .groups = "drop"
    ) %>%
    mutate(moffitt_label = label_i)
})) %>%
  rename(seurat_clusters = cluster)

cluster_annotation_long <- cluster_score_summary %>%
  left_join(cluster_overlap_summary, by = c("seurat_clusters", "moffitt_label")) %>%
  group_by(seurat_clusters) %>%
  mutate(
    score_rank = rank(-mean_signature_score, ties.method = "min"),
    overlap_rank = rank(-overlap_n, ties.method = "min"),
    combined_rank = score_rank + overlap_rank
  ) %>%
  ungroup()

cluster_annotation_best <- cluster_annotation_long %>%
  group_by(seurat_clusters) %>%
  arrange(combined_rank, desc(mean_signature_score), desc(overlap_n), .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  transmute(
    seurat_clusters = as.character(seurat_clusters),
    moffitt_label,
    best_mean_signature_score = mean_signature_score,
    overlap_n,
    overlap_genes
  )

cluster_annotation_second <- cluster_annotation_long %>%
  group_by(seurat_clusters) %>%
  arrange(combined_rank, desc(mean_signature_score), desc(overlap_n), .by_group = TRUE) %>%
  slice(2) %>%
  ungroup() %>%
  transmute(
    seurat_clusters = as.character(seurat_clusters),
    second_label = moffitt_label,
    second_mean_signature_score = mean_signature_score
  )

cluster_annotation_tbl <- cluster_annotation_best %>%
  left_join(cluster_annotation_second, by = "seurat_clusters") %>%
  mutate(
    score_gap = best_mean_signature_score - second_mean_signature_score,
    annotation_confidence = case_when(
      overlap_n >= 5 & score_gap >= 0.10 ~ "high",
      overlap_n >= 2 & score_gap >= 0.05 ~ "medium",
      TRUE ~ "low"
    ),
    clean_state = case_when(
      moffitt_label == "myCAF" ~ "myCAF",
      moffitt_label == "iCAF" ~ "iCAF",
      moffitt_label == "csCAF" ~ "csCAF",
      moffitt_label %in% c("qPSC", "smPSC") ~ "PSC",
      TRUE ~ "Excluded"
    ),
    # Clusters are excluded when their dominant projected identity lies outside
    # the four-state fibroblast panel retained for the USER-to-Moffitt
    # comparison, or when the cluster does not map convincingly to one of the
    # retained fibroblast states.
    cleaning_reason = case_when(
      clean_state %in% c("myCAF", "iCAF", "csCAF", "PSC") ~ "Retained for the four-state fibroblast validation panel",
      moffitt_label %in% c("Myocyte", "Schwann", "IL11.CAF") ~ "Removed because the dominant projected identity falls outside the four-state fibroblast panel",
      TRUE ~ "Removed because the cluster did not map to the retained fibroblast states"
    )
  ) %>%
  arrange(as.numeric(seurat_clusters))

cluster_to_label <- setNames(
  paste0("C", cluster_annotation_tbl$seurat_clusters, " | ", cluster_annotation_tbl$moffitt_label),
  cluster_annotation_tbl$seurat_clusters
)

cluster_to_clean_state <- setNames(cluster_annotation_tbl$clean_state, cluster_annotation_tbl$seurat_clusters)

fibro$cluster_label <- unname(cluster_to_label[as.character(fibro$seurat_clusters)])
fibro$moffitt_label <- unname(cluster_annotation_tbl$moffitt_label[match(as.character(fibro$seurat_clusters), cluster_annotation_tbl$seurat_clusters)])
fibro$clean_state <- unname(cluster_to_clean_state[as.character(fibro$seurat_clusters)])
fibro$clean_keep <- fibro$clean_state != "Excluded"

cluster_sizes <- fibro@meta.data %>%
  rownames_to_column("cell") %>%
  count(seurat_clusters, cluster_label, moffitt_label, clean_state, name = "n_cells") %>%
  mutate(
    pct_cells = round(100 * n_cells / sum(n_cells), 2),
    seurat_clusters = as.character(seurat_clusters)
  )

cluster_summary_table <- cluster_sizes %>%
  left_join(cluster_annotation_tbl, by = c("seurat_clusters", "moffitt_label", "clean_state")) %>%
  left_join(deg_top5, by = c("seurat_clusters" = "cluster")) %>%
  arrange(as.numeric(seurat_clusters))

fibro_meta_out <- fibro@meta.data %>%
  rownames_to_column("cell") %>%
  left_join(extract_umap_df(fibro, reduction_name = "umap", prefix = "fibro_umap"), by = "cell")

cleaning_rules_tbl <- cluster_annotation_tbl %>%
  select(seurat_clusters, moffitt_label, clean_state, annotation_confidence, overlap_n, score_gap, cleaning_reason)

write_tsv_simple(cluster_summary_table, "USER_fibro_cluster_annotation_summary.tsv")
write_tsv_simple(cluster_annotation_long, "USER_fibro_cluster_signature_projection.tsv")
write_tsv_simple(fibro_meta_out, "USER_fibro_metadata.tsv")
write_tsv_simple(cleaning_rules_tbl, "USER_fibro_cleaning_rules.tsv")

cluster_levels <- cluster_summary_table$cluster_label
cluster_colors <- make_palette(cluster_levels)

fibro$cluster_label <- factor(as.character(fibro$cluster_label), levels = cluster_levels)
Idents(fibro) <- fibro$cluster_label

p_fibro_full <- DimPlot(
  fibro,
  reduction = "umap",
  group.by = "cluster_label",
  label = TRUE,
  cols = cluster_colors,
  pt.size = 0.45
) +
  labs(
    title = "USER fibroblast UMAP",
    subtitle = "Fibroblast-only reclustering annotated with the dominant projected Moffitt stromal label"
  ) +
  theme_stage(11)

save_plot(p_fibro_full, "Figure_02_USER_fibro_umap_moffitt_labels.pdf", width = 11, height = 8)

message("Building cleaned fibroblast subset...")
cells_keep <- rownames(fibro@meta.data)[fibro$clean_keep]
fibro_clean <- subset(fibro, cells = cells_keep)
fibro_clean <- NormalizeData(fibro_clean, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
fibro_clean <- FindVariableFeatures(fibro_clean, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
fibro_clean <- ScaleData(fibro_clean, features = VariableFeatures(fibro_clean), verbose = FALSE)
fibro_clean <- RunPCA(fibro_clean, features = VariableFeatures(fibro_clean), npcs = 30, verbose = FALSE)
fibro_clean <- FindNeighbors(fibro_clean, dims = 1:20, verbose = FALSE)
fibro_clean <- RunUMAP(fibro_clean, dims = 1:20, n.neighbors = 20, min.dist = 0.20, verbose = FALSE)

fibro_clean$clean_state <- factor(fibro_clean$clean_state, levels = c("myCAF", "iCAF", "csCAF", "PSC"))
clean_state_colors <- c(
  myCAF = "#8c2d04",
  iCAF = "#cb181d",
  csCAF = "#238b45",
  PSC = "#2166ac"
)

p_fibro_clean <- DimPlot(
  fibro_clean,
  reduction = "umap",
  group.by = "clean_state",
  label = TRUE,
  cols = clean_state_colors,
  pt.size = 0.5
) +
  labs(
    title = "USER cleaned fibroblast UMAP",
    subtitle = "Clusters retained after Moffitt-guided cleaning and collapsed into myCAF, iCAF, csCAF and PSC"
  ) +
  theme_stage(11)

save_plot(p_fibro_clean, "Figure_03_USER_fibro_umap_clean_states.pdf", width = 10, height = 8)

clean_meta_out <- fibro_clean@meta.data %>%
  rownames_to_column("cell") %>%
  left_join(extract_umap_df(fibro_clean, reduction_name = "umap", prefix = "clean_umap"), by = "cell")

clean_marker_long <- get_seurat_marker_long(fibro_clean, marker_genes_of_interest, assay = "RNA", layer_name = "data")
clean_marker_wide <- clean_marker_long %>%
  pivot_wider(names_from = gene, values_from = expression)

clean_meta_out <- clean_meta_out %>%
  left_join(clean_marker_wide, by = "cell")

clean_state_counts <- clean_meta_out %>%
  count(clean_state, name = "n_cells") %>%
  mutate(pct_cells = round(100 * n_cells / sum(n_cells), 2))

clean_marker_stats <- clean_meta_out %>%
  select(cell, clean_state, all_of(marker_genes_of_interest)) %>%
  pivot_longer(cols = all_of(marker_genes_of_interest), names_to = "gene", values_to = "expression") %>%
  summarize_marker_stats("clean_state")

write_tsv_simple(clean_meta_out, "USER_clean_fibro_metadata.tsv")
write_tsv_simple(clean_state_counts, "USER_clean_fibro_state_counts.tsv")
write_tsv_simple(clean_marker_stats, "USER_clean_fibro_marker_stats.tsv")

clean_marker_plot_df <- clean_meta_out %>%
  select(cell, clean_umap_1, clean_umap_2, clean_state, all_of(marker_genes_of_interest)) %>%
  pivot_longer(cols = all_of(marker_genes_of_interest), names_to = "gene", values_to = "expression")

p_clean_markers <- ggplot(clean_marker_plot_df, aes(x = clean_umap_1, y = clean_umap_2, color = expression)) +
  geom_point(size = 0.30, alpha = 0.90) +
  facet_wrap(~ gene, ncol = 4) +
  scale_color_gradient(low = marker_color_low, high = marker_color_high) +
  labs(
    title = "USER cleaned fibroblast UMAP: selected markers",
    subtitle = "Moffitt-guided cleaned fibroblasts colored by normalized expression",
    x = "Fibro UMAP 1",
    y = "Fibro UMAP 2",
    color = "Expression"
  ) +
  theme_stage(11)

save_plot(p_clean_markers, "Figure_05_USER_clean_fibro_marker_umap.pdf", width = 14, height = 7.6)

saveRDS(fibro, file = file.path(rds_dir, "USER_fibro_moffitt_annotated_seurat.rds"))
saveRDS(fibro_clean, file = file.path(rds_dir, "USER_clean_fibro_moffitt_seurat.rds"))

validation_bundle <- list(
  inputs = list(
    matrix_path = basename(matrix_path),
    genes_path = basename(genes_path),
    barcodes_path = basename(barcodes_path),
    annotation_path = basename(annot_path),
    umap_path = basename(umap_path),
    moffitt_marker_xlsx = basename(moffitt_markers_xlsx)
  ),
  dataset_overview = dataset_overview,
  moffitt_signatures = moffitt_signatures,
  full_metadata = user_full_meta_out,
  fibro_cluster_annotation = cluster_summary_table,
  fibro_metadata = fibro_meta_out,
  clean_fibro_metadata = clean_meta_out,
  global_marker_stats = global_marker_stats,
  clean_marker_stats = clean_marker_stats
)
saveRDS(validation_bundle, file = file.path(rds_dir, "USER_moffitt_validation_bundle.rds"))

summary_lines <- c(
  "USER validation against Moffitt stromal states",
  paste0("Date: ", Sys.Date()),
  "",
  "Inputs used locally:",
  "- gene_sorted-naivedata_scp.mtx",
  "- naivedata_scp.genes.csv",
  "- naivedata_scp.barcodes.csv",
  "- combinenaivedata-reprocessed-clean-detailed-annotations.tsv",
  "- combinenaivedata-reprocessed-clean-detailed-UMAP.tsv",
  "- Moffitt_SuppData1_CellType2_Markers.xlsx",
  "",
  "USER resource provenance:",
  "- Preprint: https://www.biorxiv.org/content/10.1101/2020.08.25.267336v1.full",
  "- Code repository: https://github.com/karthikj89/humanpdac",
  "- Interactive visualization: Single Cell Portal",
  "",
  "Cleaning logic:",
  "- start from USER cells already annotated as Fibroblast",
  "- recluster the fibroblast compartment from scratch",
  "- project Moffitt top-30 stromal signatures",
  "- assign each USER fibro cluster to the best supported Moffitt label",
  "- keep only myCAF, iCAF, csCAF and qPSC/smPSC-like clusters",
  "- merge qPSC and smPSC into a single PSC category for the cleaned validation view",
  "",
  paste0("Total USER cells: ", ncol(expr_mat)),
  paste0("Fibroblast cells entering the Seurat workflow: ", ncol(fibro)),
  paste0("Fibroblast clusters before cleaning: ", length(unique(fibro$seurat_clusters))),
  paste0("Fibroblast cells retained after cleaning: ", ncol(fibro_clean)),
  paste0("Retained clean states: ", paste(levels(fibro_clean$clean_state), collapse = ", ")),
  paste0("Selected marker summaries: ", paste(marker_genes_of_interest, collapse = ", "))
)
writeLines(summary_lines, con = file.path(script_dir, "summary.txt"))


message("USER Moffitt validation complete.")


