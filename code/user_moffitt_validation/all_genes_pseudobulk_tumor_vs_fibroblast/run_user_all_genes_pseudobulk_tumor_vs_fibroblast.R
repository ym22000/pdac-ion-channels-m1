#!/usr/bin/env Rscript

###############################################################################
# USER Tumor vs Fibroblast pseudobulk DESeq2 analysis across all genes
# - patient-aware design: ~ patient + group
# - volcano plot with significant Tumor-up and Fibroblast-up genes
# - union KCN genes highlighted in black
###############################################################################

suppressPackageStartupMessages({
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(DESeq2)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)
set.seed(1234)

TUMOR_COLOR <- "#7400b8"
FIBRO_COLOR <- "#72efdd"
NS_COLOR <- "grey82"
PADJ_THRESHOLD <- 0.05
MIN_COUNT_SUM <- 10
LABEL_TOP_KCN <- 12

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  file_path <- sub(file_arg, "", args[grep(file_arg, args)])
  if (length(file_path) == 0) {
    return(normalizePath(getwd()))
  }
  normalizePath(dirname(file_path))
}

open_text_connection <- function(path) {
  con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)
  magic <- readBin(con, what = "raw", n = 2)
  if (identical(as.integer(magic), c(31L, 139L))) {
    return(gzfile(path, open = "rt"))
  }
  file(path, open = "rt")
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

make_kcn_label_layer <- function(plot_obj, label_df) {
  if (nrow(label_df) == 0) {
    return(plot_obj)
  }

  if (requireNamespace("ggrepel", quietly = TRUE)) {
    plot_obj +
      ggrepel::geom_text_repel(
        data = label_df,
        aes(x = log2FC, y = neg_log10_fdr, label = gene),
        size = 3.1,
        color = "black",
        max.overlaps = Inf,
        box.padding = 0.35,
        point.padding = 0.15,
        min.segment.length = 0
      )
  } else {
    plot_obj +
      geom_text(
        data = label_df,
        aes(x = log2FC, y = neg_log10_fdr, label = gene),
        size = 3.1,
        color = "black",
        vjust = -0.7,
        check_overlap = TRUE
      )
  }
}

script_dir <- get_script_dir()
user_root <- normalizePath(file.path(script_dir, ".."))
code_root <- normalizePath(file.path(user_root, ".."))
input_dir <- file.path(user_root, "inputs")
tab_dir <- file.path(script_dir, "tables")
fig_dir <- file.path(script_dir, "figures")
rds_dir <- file.path(script_dir, "rds")
kcn_list_path <- file.path(code_root, "mycaf_icaf_correlations", "kcn_union_25.tsv")

dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

matrix_path <- file.path(input_dir, "gene_sorted-naivedata_scp.mtx")
genes_path <- file.path(input_dir, "naivedata_scp.genes.csv")
barcodes_path <- file.path(input_dir, "naivedata_scp.barcodes.csv")
annot_path <- file.path(input_dir, "combinenaivedata-reprocessed-clean-detailed-annotations.tsv")

message("Loading union KCN list...")
kcn_genes <- read.delim(kcn_list_path, header = TRUE, stringsAsFactors = FALSE)$gene %>%
  unique() %>%
  toupper()

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

message("Loading USER annotations...")
annot_tbl <- read.delim(
  annot_path,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
if (annot_tbl$NAME[1] == "TYPE") {
  annot_tbl <- annot_tbl[-1, , drop = FALSE]
}

message("Loading USER sparse expression matrix...")
expr_mat <- readMM(matrix_path)
expr_mat <- as(expr_mat, "CsparseMatrix")
rownames(expr_mat) <- genes
colnames(expr_mat) <- barcodes

if (!identical(colnames(expr_mat), annot_tbl$NAME)) {
  stop("Cell order mismatch between matrix and annotation table.")
}

meta_df <- annot_tbl %>%
  transmute(
    barcode = NAME,
    cell_subset = as.character(cell_subsets),
    patient = as.character(pid)
  ) %>%
  filter(cell_subset %in% c("Tumor", "Fibroblast"))

common_patients <- meta_df %>%
  distinct(patient, cell_subset) %>%
  group_by(patient) %>%
  summarise(n_groups = dplyr::n(), .groups = "drop") %>%
  filter(n_groups == 2) %>%
  pull(patient)

meta_df <- meta_df %>% filter(patient %in% common_patients)

sample_id <- paste(meta_df$patient, meta_df$cell_subset, sep = "__")
sample_levels <- unique(sample_id)

cell_to_sample <- sparseMatrix(
  i = match(meta_df$barcode, colnames(expr_mat)),
  j = match(sample_id, sample_levels),
  x = 1L,
  dims = c(ncol(expr_mat), length(sample_levels)),
  dimnames = list(colnames(expr_mat), sample_levels)
)

message("Building Tumor/Fibroblast pseudobulk counts...")
pb_counts <- expr_mat %*% cell_to_sample

pb_sample_info <- tibble(sample_id = sample_levels) %>%
  separate(sample_id, into = c("patient", "cell_subset"), sep = "__", remove = FALSE) %>%
  mutate(
    n_cells = as.integer(Matrix::colSums(cell_to_sample)[sample_id]),
    library_size = as.numeric(Matrix::colSums(pb_counts)[sample_id])
  ) %>%
  arrange(patient, cell_subset)

write_tsv_simple(pb_sample_info, "USER_all_genes_tumor_vs_fibroblast_samples.tsv", tab_dir)

tumor_counts <- pb_counts[, pb_sample_info$sample_id[pb_sample_info$cell_subset == "Tumor"], drop = FALSE]
fibro_counts <- pb_counts[, pb_sample_info$sample_id[pb_sample_info$cell_subset == "Fibroblast"], drop = FALSE]

tumor_patients <- sub("__Tumor$", "", colnames(tumor_counts))
fibro_patients <- sub("__Fibroblast$", "", colnames(fibro_counts))

tumor_counts <- tumor_counts[, order(tumor_patients), drop = FALSE]
fibro_counts <- fibro_counts[, order(fibro_patients), drop = FALSE]
tumor_patients <- sub("__Tumor$", "", colnames(tumor_counts))
fibro_patients <- sub("__Fibroblast$", "", colnames(fibro_counts))

if (!identical(tumor_patients, fibro_patients)) {
  stop("Tumor/Fibroblast patient mismatch after ordering.")
}

count_matrix <- cbind(tumor_counts, fibro_counts)
coldata_df <- tibble(
  sample_id = colnames(count_matrix),
  patient = factor(c(tumor_patients, fibro_patients), levels = unique(tumor_patients)),
  group = factor(
    c(rep("Tumor", ncol(tumor_counts)), rep("Fibroblast", ncol(fibro_counts))),
    levels = c("Fibroblast", "Tumor")
  )
) %>%
  as.data.frame(stringsAsFactors = FALSE)
rownames(coldata_df) <- coldata_df$sample_id
coldata_df$sample_id <- NULL

message("Running DESeq2 pseudobulk Tumor vs Fibroblast on all genes...")
dds <- DESeqDataSetFromMatrix(
  countData = round(count_matrix),
  colData = coldata_df,
  design = ~ patient + group
)

dds <- dds[rowSums(counts(dds)) >= MIN_COUNT_SUM, , drop = FALSE]
dds <- estimateSizeFactors(dds, type = "poscounts")
dds <- DESeq(dds, quiet = TRUE, sfType = "poscounts", fitType = "local")

res_tbl <- results(
  dds,
  contrast = c("group", "Tumor", "Fibroblast"),
  independentFiltering = TRUE
) %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  transmute(
    gene = gene,
    baseMean = baseMean,
    log2FC = log2FoldChange,
    lfcSE = lfcSE,
    stat = stat,
    p_value = pvalue,
    FDR = padj
  ) %>%
  mutate(
    direction = case_when(
      is.na(FDR) ~ "Not tested",
      FDR < PADJ_THRESHOLD & log2FC > 0 ~ "Tumor up",
      FDR < PADJ_THRESHOLD & log2FC < 0 ~ "Fibroblast up",
      TRUE ~ "Not significant"
    ),
    is_kcn_union = gene %in% kcn_genes
  ) %>%
  arrange(FDR, desc(abs(log2FC)), gene)

kcn_tbl <- res_tbl %>%
  filter(is_kcn_union) %>%
  mutate(
    Lecture = case_when(
      is.na(log2FC) ~ "NA",
      log2FC >= 0 ~ "Tumor up",
      TRUE ~ "Fibroblast up"
    )
  ) %>%
  select(gene, log2FC, FDR, Lecture, everything())

comparison_info <- tibble(
  comparison = "Tumor_vs_Fibroblast",
  n_patients = length(unique(tumor_patients)),
  tumor_total_cells = sum(pb_sample_info$n_cells[pb_sample_info$cell_subset == "Tumor"]),
  fibroblast_total_cells = sum(pb_sample_info$n_cells[pb_sample_info$cell_subset == "Fibroblast"]),
  n_genes_tested = nrow(res_tbl),
  n_tumor_up_fdr_lt_0_05 = sum(res_tbl$direction == "Tumor up", na.rm = TRUE),
  n_fibroblast_up_fdr_lt_0_05 = sum(res_tbl$direction == "Fibroblast up", na.rm = TRUE)
)

write_tsv_simple(res_tbl, "USER_all_genes_pseudobulk_tumor_vs_fibroblast.tsv", tab_dir)
write_tsv_simple(kcn_tbl, "USER_kcn25_within_all_genes_tumor_vs_fibroblast.tsv", tab_dir)
write_tsv_simple(comparison_info, "USER_all_genes_pseudobulk_tumor_vs_fibroblast_info.tsv", tab_dir)

volcano_df <- res_tbl %>%
  mutate(
    plot_fdr = ifelse(is.na(FDR) | FDR <= 0, .Machine$double.xmin, FDR),
    neg_log10_fdr = -log10(plot_fdr)
  )

kcn_highlight_df <- volcano_df %>%
  filter(is_kcn_union, !is.na(log2FC), !is.na(FDR), FDR < PADJ_THRESHOLD) %>%
  arrange(FDR, desc(abs(log2FC)))

kcn_label_df <- kcn_highlight_df %>% slice_head(n = LABEL_TOP_KCN)

volcano_plot <- ggplot(volcano_df, aes(x = log2FC, y = neg_log10_fdr)) +
  geom_point(
    data = volcano_df %>% filter(direction == "Not significant"),
    color = NS_COLOR,
    size = 0.7,
    alpha = 0.7
  ) +
  geom_point(
    data = volcano_df %>% filter(direction == "Tumor up"),
    aes(color = direction),
    size = 0.8,
    alpha = 0.8
  ) +
  geom_point(
    data = volcano_df %>% filter(direction == "Fibroblast up"),
    aes(color = direction),
    size = 0.8,
    alpha = 0.8
  ) +
  geom_point(
    data = kcn_highlight_df,
    color = "black",
    size = 1.8,
    alpha = 1
  ) +
  geom_hline(yintercept = -log10(PADJ_THRESHOLD), linetype = "dashed", color = "grey45", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey45", linewidth = 0.35) +
  scale_color_manual(
    values = c("Tumor up" = TUMOR_COLOR, "Fibroblast up" = FIBRO_COLOR),
    breaks = c("Tumor up", "Fibroblast up"),
    labels = c("Tumor", "Fibroblast"),
    name = NULL
  ) +
  labs(
    title = "USER pseudobulk DESeq2: Tumor vs Fibroblast",
    x = "log2FC (Tumor vs Fibroblast)",
    y = expression(-log[10]("FDR"))
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.key = element_blank()
  )

volcano_plot <- make_kcn_label_layer(volcano_plot, kcn_label_df)
volcano_plot <- volcano_plot +
  guides(color = guide_legend(override.aes = list(shape = 15, size = 5, alpha = 1)))

ggsave(
  filename = file.path(fig_dir, "USER_all_genes_pseudobulk_tumor_vs_fibroblast_volcano.pdf"),
  plot = volcano_plot,
  width = 8.8,
  height = 6.2,
  useDingbats = FALSE
)

saveRDS(
  list(
    sample_info = pb_sample_info,
    comparison_info = comparison_info,
    deseq_results = res_tbl,
    kcn_results = kcn_tbl
  ),
  file = file.path(rds_dir, "USER_all_genes_pseudobulk_tumor_vs_fibroblast_results.rds")
)

message("USER all-gene pseudobulk Tumor vs Fibroblast analysis complete.")
