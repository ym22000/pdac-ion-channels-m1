#!/usr/bin/env Rscript

###############################################################################
# USER cell-subset MAST analysis for the 25 KCN union genes
#
# Goal:
#   Identify where the 25 union KCN genes are enriched across the main USER
#   cell subsets using a simple one-vs-rest MAST design.
###############################################################################

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(MAST)
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

if (!file.exists(kcn_list_path)) {
  stop("Missing KCN list: ", kcn_list_path)
}

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

if (nrow(expr_mat) != length(genes)) stop("Gene count mismatch between matrix and gene file.")
if (ncol(expr_mat) != length(barcodes)) stop("Barcode count mismatch between matrix and barcode file.")

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

message("Building a compact Seurat object for the KCN panel...")
expr_kcn <- expr_mat[present_kcn, , drop = FALSE]
meta_df <- annot_tbl %>%
  select(NAME, cell_subsets, pid, n_genes) %>%
  as.data.frame(stringsAsFactors = FALSE)
rownames(meta_df) <- meta_df$NAME
meta_df$NAME <- NULL

user_kcn <- CreateSeuratObject(counts = expr_kcn, meta.data = meta_df)
user_kcn <- NormalizeData(user_kcn, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
Idents(user_kcn) <- user_kcn$cell_subsets

cell_subset_counts <- tibble(
  cell_subset = as.character(user_kcn$cell_subsets),
  pid = as.character(user_kcn$pid)
) %>%
  group_by(cell_subset) %>%
  summarise(n_cells = n(), n_patients = n_distinct(pid), .groups = "drop") %>%
  mutate(
    pct_cells = round(100 * n_cells / sum(n_cells), 2)
  ) %>%
  arrange(desc(n_cells))
write_tsv_simple(cell_subset_counts, "USER_kcn25_cell_subset_counts.tsv", tab_dir)

message("Running one-vs-rest MAST across USER cell subsets...")
mast_markers <- FindAllMarkers(
  object = user_kcn,
  assay = "RNA",
  only.pos = FALSE,
  features = present_kcn,
  group.by = "cell_subsets",
  test.use = "MAST",
  latent.vars = "nCount_RNA",
  logfc.threshold = 0,
  min.pct = 0,
  min.diff.pct = -Inf,
  verbose = FALSE
)

if (!"avg_log2FC" %in% colnames(mast_markers) && "avg_logFC" %in% colnames(mast_markers)) {
  mast_markers <- mast_markers %>% rename(avg_log2FC = avg_logFC)
}

mast_kcn_summary <- mast_markers %>%
  transmute(
    cell_subset = cluster,
    KCN = gene,
    log2FC = avg_log2FC,
    FDR = p_val_adj,
    Lecture = ifelse(avg_log2FC >= 0, paste0(cluster, " up"), paste0(cluster, " down"))
  ) %>%
  arrange(cell_subset, FDR, desc(abs(log2FC)), KCN)

mast_kcn_significant <- mast_kcn_summary %>%
  filter(FDR < 0.05)

missing_rows <- tibble(
  cell_subset = NA_character_,
  KCN = missing_kcn,
  log2FC = NA_real_,
  FDR = NA_real_,
  Lecture = "Not detected in USER"
) %>%
  arrange(KCN)

write_tsv_simple(mast_kcn_summary, "USER_kcn25_mast_by_cell_subset.tsv", tab_dir)
write_tsv_simple(mast_kcn_significant, "USER_kcn25_mast_by_cell_subset_significant.tsv", tab_dir)
write_tsv_simple(missing_rows, "USER_kcn25_missing_in_user.tsv", tab_dir)

saveRDS(
  list(
    present_kcn = present_kcn,
    missing_kcn = missing_kcn,
    cell_subset_counts = cell_subset_counts,
    mast_results = mast_kcn_summary,
    mast_significant = mast_kcn_significant
  ),
  file = file.path(rds_dir, "USER_kcn25_mast_cell_subsets_results.rds")
)

message("USER KCN25 MAST cell-subset analysis complete.")
