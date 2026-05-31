#!/usr/bin/env Rscript

###############################################################################
# USER Tumor vs Fibroblast pseudobulk DESeq2 analysis for the 25 KCN union genes
###############################################################################

suppressPackageStartupMessages({
  library(Matrix)
  library(dplyr)
  library(tibble)
  library(DESeq2)
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

script_dir <- get_script_dir()
user_root <- normalizePath(file.path(script_dir, ".."))
code_root <- normalizePath(file.path(user_root, ".."))
input_dir <- file.path(user_root, "inputs")
tab_dir <- file.path(script_dir, "tables")
rds_dir <- file.path(script_dir, "rds")
kcn_list_path <- file.path(code_root, "mycaf_icaf_correlations", "kcn_union_25.tsv")

dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

matrix_path <- file.path(input_dir, "gene_sorted-naivedata_scp.mtx")
genes_path <- file.path(input_dir, "naivedata_scp.genes.csv")
barcodes_path <- file.path(input_dir, "naivedata_scp.barcodes.csv")
annot_path <- file.path(input_dir, "combinenaivedata-reprocessed-clean-detailed-annotations.tsv")

message("Loading 25 KCN union list...")
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
annot_tbl <- read.delim(annot_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
if (annot_tbl$NAME[1] == "TYPE") annot_tbl <- annot_tbl[-1, , drop = FALSE]

message("Loading USER sparse expression matrix...")
expr_mat <- readMM(matrix_path)
expr_mat <- as(expr_mat, "CsparseMatrix")

rownames(expr_mat) <- genes
colnames(expr_mat) <- barcodes

if (!identical(colnames(expr_mat), annot_tbl$NAME)) {
  stop("Cell order mismatch between matrix and annotation table.")
}

present_kcn <- intersect(kcn_genes, rownames(expr_mat))
missing_kcn <- setdiff(kcn_genes, rownames(expr_mat))
if (length(present_kcn) == 0) {
  stop("None of the 25 KCN genes were found in USER.")
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
  summarise(n_groups = n(), .groups = "drop") %>%
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
expr_kcn <- expr_mat[present_kcn, , drop = FALSE]
pb_counts <- expr_kcn %*% cell_to_sample

pb_sample_info <- tibble(sample_id = sample_levels) %>%
  tidyr::separate(sample_id, into = c("patient", "cell_subset"), sep = "__", remove = FALSE) %>%
  mutate(
    n_cells = as.integer(Matrix::colSums(cell_to_sample)[sample_id]),
    library_size = as.numeric(Matrix::colSums(pb_counts)[sample_id])
  ) %>%
  arrange(patient, cell_subset)

write_tsv_simple(pb_sample_info, "USER_kcn25_tumor_vs_fibroblast_samples.tsv", tab_dir)

target_counts <- pb_counts[, pb_sample_info$sample_id[pb_sample_info$cell_subset == "Tumor"], drop = FALSE]
fibro_counts <- pb_counts[, pb_sample_info$sample_id[pb_sample_info$cell_subset == "Fibroblast"], drop = FALSE]

target_patients <- sub("__Tumor$", "", colnames(target_counts))
fibro_patients <- sub("__Fibroblast$", "", colnames(fibro_counts))

target_counts <- target_counts[, order(target_patients), drop = FALSE]
fibro_counts <- fibro_counts[, order(fibro_patients), drop = FALSE]
target_patients <- sub("__Tumor$", "", colnames(target_counts))
fibro_patients <- sub("__Fibroblast$", "", colnames(fibro_counts))

if (!identical(target_patients, fibro_patients)) {
  stop("Tumor/Fibroblast patient mismatch after ordering.")
}

count_matrix <- cbind(target_counts, fibro_counts)
coldata_df <- tibble(
  sample_id = colnames(count_matrix),
  patient = factor(c(target_patients, fibro_patients), levels = unique(target_patients)),
  group = factor(c(rep("Tumor", ncol(target_counts)), rep("Fibroblast", ncol(fibro_counts))), levels = c("Fibroblast", "Tumor"))
) %>%
  as.data.frame(stringsAsFactors = FALSE)
rownames(coldata_df) <- coldata_df$sample_id
coldata_df$sample_id <- NULL

message("Running DESeq2 pseudobulk Tumor vs Fibroblast...")
dds <- DESeqDataSetFromMatrix(
  countData = round(count_matrix),
  colData = coldata_df,
  design = ~ patient + group
)

dds <- dds[rowSums(counts(dds)) >= 10, , drop = FALSE]
dds <- estimateSizeFactors(dds, type = "poscounts")
dds <- DESeq(dds, quiet = TRUE, sfType = "poscounts", fitType = "local")

res_tbl <- results(dds, contrast = c("group", "Tumor", "Fibroblast"), independentFiltering = TRUE) %>%
  as.data.frame() %>%
  rownames_to_column("KCN") %>%
  transmute(
    KCN = KCN,
    log2FC = log2FoldChange,
    FDR = padj,
    Lecture = ifelse(log2FoldChange >= 0, "Tumor up", "Fibroblast up")
  ) %>%
  arrange(FDR, desc(abs(log2FC)), KCN)

res_complete <- tibble(KCN = present_kcn) %>%
  left_join(res_tbl, by = "KCN") %>%
  mutate(
    Lecture = ifelse(is.na(log2FC), "NA", ifelse(log2FC >= 0, "Tumor up", "Fibroblast up"))
  ) %>%
  arrange(FDR, desc(abs(log2FC)), KCN)

res_significant <- res_complete %>% filter(!is.na(FDR) & FDR < 0.05)

missing_rows <- tibble(
  KCN = missing_kcn,
  log2FC = NA_real_,
  FDR = NA_real_,
  Lecture = "Not detected in USER"
)

comparison_info <- tibble(
  comparison = "Tumor_vs_Fibroblast",
  n_patients = length(unique(target_patients)),
  tumor_total_cells = sum(pb_sample_info$n_cells[pb_sample_info$cell_subset == "Tumor"]),
  fibroblast_total_cells = sum(pb_sample_info$n_cells[pb_sample_info$cell_subset == "Fibroblast"])
)

write_tsv_simple(res_complete, "USER_kcn25_pseudobulk_tumor_vs_fibroblast.tsv", tab_dir)
write_tsv_simple(res_significant, "USER_kcn25_pseudobulk_tumor_vs_fibroblast_significant.tsv", tab_dir)
write_tsv_simple(missing_rows, "USER_kcn25_missing_in_user.tsv", tab_dir)
write_tsv_simple(comparison_info, "USER_kcn25_pseudobulk_tumor_vs_fibroblast_info.tsv", tab_dir)

saveRDS(
  list(
    present_kcn = present_kcn,
    missing_kcn = missing_kcn,
    sample_info = pb_sample_info,
    comparison_info = comparison_info,
    deseq_results = res_complete,
    deseq_significant = res_significant
  ),
  file = file.path(rds_dir, "USER_kcn25_pseudobulk_tumor_vs_fibroblast_results.rds")
)

message("USER KCN25 pseudobulk Tumor vs Fibroblast analysis complete.")
