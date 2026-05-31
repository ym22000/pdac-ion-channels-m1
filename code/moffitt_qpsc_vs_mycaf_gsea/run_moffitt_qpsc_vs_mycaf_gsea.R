#!/usr/bin/env Rscript

###############################################################################
# Moffitt qPSC vs myCAF pseudobulk Hallmark GSEA
#
# This script builds paired pseudobulk profiles for qPSC and myCAF samples,
# runs a patient-aware DESeq2 contrast, and performs Hallmark fgsea.
###############################################################################

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(DESeq2)
  library(msigdbr)
  library(fgsea)
  library(patchwork)
  library(stringr)
})

options(stringsAsFactors = FALSE)
set.seed(1234)

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  file_path <- sub(file_arg, "", args[grep(file_arg, args)])
  if (length(file_path) == 0) {
    return(normalizePath(getwd(), winslash = "/"))
  }
  normalizePath(dirname(file_path), winslash = "/")
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

save_plot <- function(plot_object, filename, out_dir, width = 10, height = 8, dpi = 320) {
  ggsave(
    filename = file.path(out_dir, filename),
    plot = plot_object,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
}

clean_theme <- function() {
  theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 11),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 10),
      panel.grid = element_blank()
    )
}

make_pathway_label <- function(x) {
  x %>%
    str_remove("^HALLMARK_") %>%
    str_replace_all("_", " ") %>%
    str_to_title()
}

script_dir <- get_script_dir()
code_root <- normalizePath(file.path(script_dir, ".."), winslash = "/")
project_root <- normalizePath(file.path(code_root, "..", ".."), winslash = "/")
data_root <- file.path(project_root, "M1_Bioinformatics_Ion_Channel", "data")

input_rds <- file.path(data_root, "Stroma_Subset2021_fibroblast_focus.rds")
output_root <- script_dir
figure_dir <- file.path(output_root, "figures")
table_dir <- file.path(output_root, "tables")
rds_dir <- file.path(output_root, "rds")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

subtype_colors <- c(
  "qPSC" = "#70D6FF",
  "myCAF" = "#FF686B"
)

min_cells_per_sample <- 20L

if (!file.exists(input_rds)) {
  stop("Input object not found: ", input_rds)
}

message("Loading stromal Seurat object...")
stroma <- readRDS(input_rds)
DefaultAssay(stroma) <- "RNA"

required_meta <- c("Patient", "Dataset", "CellType2")
missing_meta <- setdiff(required_meta, colnames(stroma@meta.data))
if (length(missing_meta) > 0) {
  stop("Missing metadata columns: ", paste(missing_meta, collapse = ", "))
}

meta_df <- stroma@meta.data %>%
  rownames_to_column("cell_id") %>%
  transmute(
    cell_id = cell_id,
    patient = as.character(Patient),
    dataset = as.character(Dataset),
    subtype = as.character(CellType2)
  ) %>%
  filter(subtype %in% c("qPSC", "myCAF"))

sample_info_all <- meta_df %>%
  dplyr::count(patient, dataset, subtype, name = "n_cells") %>%
  mutate(keep_sample = n_cells >= min_cells_per_sample)

paired_patients <- sample_info_all %>%
  filter(keep_sample) %>%
  dplyr::distinct(patient, subtype) %>%
  dplyr::count(patient, name = "n_subtypes") %>%
  filter(n_subtypes == 2) %>%
  pull(patient)

sample_info <- sample_info_all %>%
  filter(patient %in% paired_patients, keep_sample) %>%
  mutate(sample_id = paste(patient, subtype, sep = "__")) %>%
  arrange(patient, subtype)

if (nrow(sample_info) == 0) {
  stop("No paired qPSC/myCAF pseudobulk samples remained after filtering.")
}

meta_keep <- meta_df %>%
  inner_join(sample_info %>% select(patient, subtype, sample_id), by = c("patient", "subtype"))

message("Building paired pseudobulk matrix...")
counts_mat <- GetAssayData(stroma, assay = "RNA", layer = "counts")
counts_mat <- counts_mat[, meta_keep$cell_id, drop = FALSE]

sample_levels <- sample_info$sample_id
cell_to_sample <- sparseMatrix(
  i = seq_len(ncol(counts_mat)),
  j = match(meta_keep$sample_id, sample_levels),
  x = 1L,
  dims = c(ncol(counts_mat), length(sample_levels)),
  dimnames = list(meta_keep$cell_id, sample_levels)
)

pb_counts <- counts_mat %*% cell_to_sample
sample_info$library_size <- as.numeric(Matrix::colSums(pb_counts)[sample_info$sample_id])

write_tsv_simple(sample_info, "moffitt_qpsc_vs_mycaf_pseudobulk_samples.tsv", table_dir)

qpsc_counts <- pb_counts[, sample_info$sample_id[sample_info$subtype == "qPSC"], drop = FALSE]
mycaf_counts <- pb_counts[, sample_info$sample_id[sample_info$subtype == "myCAF"], drop = FALSE]

qpsc_patients <- sub("__qPSC$", "", colnames(qpsc_counts))
mycaf_patients <- sub("__myCAF$", "", colnames(mycaf_counts))

qpsc_counts <- qpsc_counts[, order(qpsc_patients), drop = FALSE]
mycaf_counts <- mycaf_counts[, order(mycaf_patients), drop = FALSE]
qpsc_patients <- sub("__qPSC$", "", colnames(qpsc_counts))
mycaf_patients <- sub("__myCAF$", "", colnames(mycaf_counts))

if (!identical(qpsc_patients, mycaf_patients)) {
  stop("Paired patient mismatch between qPSC and myCAF pseudobulk samples.")
}

count_matrix <- cbind(qpsc_counts, mycaf_counts)
coldata_df <- tibble(
  sample_id = colnames(count_matrix),
  patient = factor(c(qpsc_patients, mycaf_patients), levels = unique(qpsc_patients)),
  subtype = factor(
    c(rep("qPSC", ncol(qpsc_counts)), rep("myCAF", ncol(mycaf_counts))),
    levels = c("qPSC", "myCAF")
  )
) %>%
  as.data.frame(stringsAsFactors = FALSE)
rownames(coldata_df) <- coldata_df$sample_id
coldata_df$sample_id <- NULL

message("Running patient-aware DESeq2 contrast...")
dds <- DESeqDataSetFromMatrix(
  countData = round(count_matrix),
  colData = coldata_df,
  design = ~ patient + subtype
)

dds <- dds[rowSums(counts(dds)) >= 10, , drop = FALSE]
dds <- estimateSizeFactors(dds, type = "poscounts")
dds <- DESeq(dds, quiet = TRUE, sfType = "poscounts", fitType = "local")

res_tbl <- results(dds, contrast = c("subtype", "myCAF", "qPSC"), independentFiltering = TRUE) %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  transmute(
    gene = gene,
    baseMean = baseMean,
    log2FC = log2FoldChange,
    lfcSE = lfcSE,
    stat = stat,
    pvalue = pvalue,
    padj = padj
  ) %>%
  arrange(padj, desc(abs(log2FC)), gene)

rank_tbl <- res_tbl %>%
  filter(!is.na(stat)) %>%
  arrange(desc(stat))

write_tsv_simple(res_tbl, "moffitt_qpsc_vs_mycaf_deseq2_results.tsv", table_dir)
write_tsv_simple(rank_tbl, "moffitt_qpsc_vs_mycaf_ranked_stats.tsv", table_dir)

message("Loading Hallmark gene sets...")
hallmark_tbl <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  select(gs_name, gene_symbol) %>%
  distinct()

hallmark_list <- split(hallmark_tbl$gene_symbol, hallmark_tbl$gs_name)

rank_vector <- rank_tbl$stat
names(rank_vector) <- rank_tbl$gene

message("Running fgsea...")
fgsea_tbl <- fgseaMultilevel(
  pathways = hallmark_list,
  stats = rank_vector,
  eps = 0
) %>%
  as_tibble() %>%
  mutate(
    leadingEdge = vapply(leadingEdge, function(x) paste(x, collapse = ";"), character(1)),
    pathway_label = make_pathway_label(pathway),
    enriched_in = ifelse(NES >= 0, "myCAF", "qPSC"),
    direction_rank = abs(NES)
  ) %>%
  arrange(padj, desc(abs(NES)), pathway)

fgsea_sig <- fgsea_tbl %>%
  filter(!is.na(padj), padj < 0.05)

write_tsv_simple(fgsea_tbl, "moffitt_qpsc_vs_mycaf_hallmark_fgsea.tsv", table_dir)
write_tsv_simple(fgsea_sig, "moffitt_qpsc_vs_mycaf_hallmark_fgsea_significant.tsv", table_dir)

plot_tbl <- fgsea_sig %>%
  group_by(enriched_in) %>%
  slice_max(order_by = abs(NES), n = 10, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(NES) %>%
  mutate(pathway_label = factor(pathway_label, levels = pathway_label))

if (nrow(plot_tbl) == 0) {
  plot_tbl <- fgsea_tbl %>%
    group_by(enriched_in) %>%
    slice_max(order_by = abs(NES), n = 10, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(NES) %>%
    mutate(pathway_label = factor(pathway_label, levels = pathway_label))
}

summary_plot <- ggplot(plot_tbl, aes(x = NES, y = pathway_label, fill = enriched_in)) +
  geom_col(width = 0.75, color = NA) +
  scale_fill_manual(
    values = subtype_colors,
    name = "Enriched in",
    breaks = c("qPSC", "myCAF")
  ) +
  labs(
    title = "Hallmark GSEA: qPSC vs myCAF",
    x = "Normalized enrichment score"
  ) +
  clean_theme()

save_plot(
  summary_plot,
  "moffitt_qpsc_vs_mycaf_hallmark_gsea_summary.pdf",
  figure_dir,
  width = 11,
  height = 8
)

top_mycaf <- fgsea_tbl %>%
  filter(NES > 0) %>%
  slice_min(order_by = padj, n = 1, with_ties = FALSE)

top_qpsc <- fgsea_tbl %>%
  filter(NES < 0) %>%
  slice_min(order_by = padj, n = 1, with_ties = FALSE)

curve_plots <- list()

if (nrow(top_mycaf) == 1) {
  curve_plots[["myCAF"]] <- plotEnrichment(
    pathway = hallmark_list[[top_mycaf$pathway]],
    stats = rank_vector
  ) +
    labs(
      title = paste0("myCAF: ", top_mycaf$pathway_label),
      x = "Ranked genes",
      y = "Enrichment score"
    ) +
    geom_hline(yintercept = 0, linetype = "solid", color = "grey70") +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.grid = element_blank()
    )
}

if (nrow(top_qpsc) == 1) {
  curve_plots[["qPSC"]] <- plotEnrichment(
    pathway = hallmark_list[[top_qpsc$pathway]],
    stats = rank_vector
  ) +
    labs(
      title = paste0("qPSC: ", top_qpsc$pathway_label),
      x = "Ranked genes",
      y = "Enrichment score"
    ) +
    geom_hline(yintercept = 0, linetype = "solid", color = "grey70") +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.grid = element_blank()
    )
}

if (length(curve_plots) > 0) {
  curve_panel <- wrap_plots(curve_plots, ncol = length(curve_plots))
  save_plot(
    curve_panel,
    "moffitt_qpsc_vs_mycaf_hallmark_top_curves.pdf",
    figure_dir,
    width = 12,
    height = 5
  )
}

analysis_summary <- tibble(
  comparison = "qPSC_vs_myCAF",
  n_paired_patients = length(unique(qpsc_patients)),
  min_cells_per_sample = min_cells_per_sample,
  qpsc_total_cells = sum(sample_info$n_cells[sample_info$subtype == "qPSC"]),
  mycaf_total_cells = sum(sample_info$n_cells[sample_info$subtype == "myCAF"]),
  n_genes_tested = nrow(res_tbl),
  n_significant_pathways = nrow(fgsea_sig)
)

write_tsv_simple(analysis_summary, "moffitt_qpsc_vs_mycaf_analysis_summary.tsv", table_dir)

saveRDS(
  list(
    sample_info = sample_info,
    deseq_results = res_tbl,
    fgsea_results = fgsea_tbl,
    fgsea_significant = fgsea_sig,
    analysis_summary = analysis_summary
  ),
  file = file.path(rds_dir, "moffitt_qpsc_vs_mycaf_gsea_results.rds")
)

message("Moffitt qPSC vs myCAF GSEA complete.")
