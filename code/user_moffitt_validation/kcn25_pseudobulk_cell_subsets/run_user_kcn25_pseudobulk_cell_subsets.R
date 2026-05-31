#!/usr/bin/env Rscript

###############################################################################
# USER cell-subset pseudobulk DESeq2 analysis for the 25 KCN union genes
#
# Goal:
#   Summarize the 25 union KCN genes across USER cell subsets with a simple
#   patient-aware pseudobulk design. For each cell subset, counts are compared
#   against the pooled "rest" from the same patients using ~ patient + group.
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

run_deseq_pseudobulk <- function(count_matrix, coldata_df) {
  dds <- DESeqDataSetFromMatrix(
    countData = round(count_matrix),
    colData = as.data.frame(coldata_df),
    design = ~ patient + group
  )

  dds <- dds[rowSums(counts(dds)) >= 10, , drop = FALSE]
  if (nrow(dds) == 0) {
    return(NULL)
  }

  dds <- tryCatch(
    estimateSizeFactors(dds, type = "poscounts"),
    error = function(e) dds
  )

  current_size_factors <- tryCatch(
    sizeFactors(dds),
    error = function(e) rep(NA_real_, ncol(dds))
  )

  if (any(is.na(current_size_factors))) {
    library_sizes <- colSums(counts(dds))
    positive_sizes <- library_sizes[library_sizes > 0]
    if (length(positive_sizes) == 0) {
      fallback_size_factors <- rep(1, length(library_sizes))
    } else {
      fallback_size_factors <- library_sizes / exp(mean(log(positive_sizes)))
    }
    fallback_size_factors[!is.finite(fallback_size_factors)] <- 1
    fallback_size_factors[fallback_size_factors <= 0] <- 1
    sizeFactors(dds) <- fallback_size_factors
  }

  dds <- DESeq(dds, quiet = TRUE, sfType = "poscounts", fitType = "local")
  res <- results(dds, contrast = c("group", "target", "rest"), independentFiltering = TRUE) %>%
    as.data.frame() %>%
    rownames_to_column("KCN") %>%
    transmute(
      KCN = KCN,
      log2FC = log2FoldChange,
      FDR = padj
    )

  res
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

message("Restricting to the KCN panel...")
expr_kcn <- expr_mat[present_kcn, , drop = FALSE]
meta_df <- annot_tbl %>%
  transmute(
    barcode = NAME,
    cell_subset = as.character(cell_subsets),
    patient = as.character(pid)
  )

sample_id <- paste(meta_df$patient, meta_df$cell_subset, sep = "__")
sample_levels <- unique(sample_id)
cell_to_sample <- sparseMatrix(
  i = seq_along(sample_id),
  j = match(sample_id, sample_levels),
  x = 1L,
  dims = c(length(sample_id), length(sample_levels)),
  dimnames = list(meta_df$barcode, sample_levels)
)

message("Building pseudobulk counts per patient and cell subset...")
pb_counts <- expr_kcn %*% cell_to_sample
pb_sample_info <- tibble(sample_id = sample_levels) %>%
  tidyr::separate(sample_id, into = c("patient", "cell_subset"), sep = "__", remove = FALSE) %>%
  mutate(
    n_cells = as.integer(Matrix::colSums(cell_to_sample)[sample_id]),
    library_size = as.numeric(Matrix::colSums(pb_counts)[sample_id])
  ) %>%
  arrange(cell_subset, patient)

write_tsv_simple(pb_sample_info, "USER_kcn25_pseudobulk_samples.tsv", tab_dir)

cell_subset_counts <- pb_sample_info %>%
  group_by(cell_subset) %>%
  summarise(
    n_patients = n(),
    total_cells = sum(n_cells),
    median_cells_per_patient = median(n_cells),
    .groups = "drop"
  ) %>%
  arrange(desc(total_cells))
write_tsv_simple(cell_subset_counts, "USER_kcn25_pseudobulk_cell_subset_counts.tsv", tab_dir)

cell_subsets <- unique(pb_sample_info$cell_subset)
all_results <- vector("list", length(cell_subsets))
comparison_info <- vector("list", length(cell_subsets))

message("Running one-vs-rest pseudobulk DESeq2 across USER cell subsets...")
for (i in seq_along(cell_subsets)) {
  target_subset <- cell_subsets[i]
  target_info <- pb_sample_info %>% filter(cell_subset == target_subset)
  target_patients <- unique(target_info$patient)

  rest_candidates <- pb_sample_info %>%
    filter(patient %in% target_patients, cell_subset != target_subset)

  rest_patients <- unique(rest_candidates$patient)
  common_patients <- intersect(target_patients, rest_patients)

  comparison_info[[i]] <- tibble(
    cell_subset = target_subset,
    n_patients_target = length(target_patients),
    n_patients_candidate = length(common_patients),
    n_patients_used = NA_integer_
  )

  if (length(common_patients) < 2) {
    all_results[[i]] <- tibble(
      cell_subset = target_subset,
      KCN = present_kcn,
      log2FC = NA_real_,
      FDR = NA_real_,
      Lecture = paste0(target_subset, " up")
    )
    next
  }

  target_matrix <- pb_counts[
    present_kcn,
    pb_sample_info$sample_id[pb_sample_info$patient %in% common_patients & pb_sample_info$cell_subset == target_subset],
    drop = FALSE
  ]
  target_matrix <- target_matrix[, order(colnames(target_matrix)), drop = FALSE]

  rest_matrix <- sapply(common_patients, function(patient_id) {
    rest_cols <- pb_sample_info$sample_id[pb_sample_info$patient == patient_id & pb_sample_info$cell_subset != target_subset]
    Matrix::rowSums(pb_counts[present_kcn, rest_cols, drop = FALSE])
  })
  rest_matrix <- as.matrix(rest_matrix)
  if (is.null(dim(rest_matrix))) {
    rest_matrix <- matrix(rest_matrix, nrow = length(present_kcn), dimnames = list(present_kcn, common_patients))
  }
  colnames(rest_matrix) <- paste0(common_patients, "__rest")

  target_patients_ordered <- sub("__.*$", "", colnames(target_matrix))
  rest_patients_ordered <- sub("__rest$", "", colnames(rest_matrix))
  if (!identical(sort(target_patients_ordered), sort(rest_patients_ordered))) {
    stop("Patient mismatch while building target/rest pseudobulks for ", target_subset)
  }

  target_matrix <- target_matrix[, order(target_patients_ordered), drop = FALSE]
  rest_matrix <- rest_matrix[, order(rest_patients_ordered), drop = FALSE]

  target_patients_ordered <- sub("__.*$", "", colnames(target_matrix))
  rest_patients_ordered <- sub("__rest$", "", colnames(rest_matrix))

  keep_patients <- (Matrix::colSums(target_matrix) > 0) & (colSums(rest_matrix) > 0)
  if (sum(keep_patients) < 2) {
    comparison_info[[i]]$n_patients_used <- sum(keep_patients)
    all_results[[i]] <- tibble(
      cell_subset = target_subset,
      KCN = present_kcn,
      log2FC = NA_real_,
      FDR = NA_real_,
      Lecture = paste0(target_subset, " up")
    )
    next
  }

  target_matrix <- target_matrix[, keep_patients, drop = FALSE]
  rest_matrix <- rest_matrix[, keep_patients, drop = FALSE]
  target_patients_ordered <- target_patients_ordered[keep_patients]
  rest_patients_ordered <- rest_patients_ordered[keep_patients]
  comparison_info[[i]]$n_patients_used <- sum(keep_patients)

  count_matrix <- cbind(target_matrix, rest_matrix)
  coldata_df <- tibble(
    sample_id = colnames(count_matrix),
    patient = factor(c(target_patients_ordered, rest_patients_ordered), levels = unique(c(target_patients_ordered, rest_patients_ordered))),
    group = factor(c(rep("target", ncol(target_matrix)), rep("rest", ncol(rest_matrix))), levels = c("rest", "target"))
  ) %>%
    as.data.frame(stringsAsFactors = FALSE)
  rownames(coldata_df) <- coldata_df$sample_id
  coldata_df$sample_id <- NULL

  res_tbl <- run_deseq_pseudobulk(count_matrix = count_matrix, coldata_df = coldata_df)

  if (is.null(res_tbl)) {
    all_results[[i]] <- tibble(
      cell_subset = target_subset,
      KCN = present_kcn,
      log2FC = NA_real_,
      FDR = NA_real_,
      Lecture = paste0(target_subset, " up")
    )
    next
  }

  res_tbl <- tibble(KCN = present_kcn) %>%
    left_join(res_tbl, by = "KCN") %>%
    mutate(
      cell_subset = target_subset,
      Lecture = ifelse(is.na(log2FC), paste0(target_subset, " up"), ifelse(log2FC >= 0, paste0(target_subset, " up"), paste0(target_subset, " down")))
    ) %>%
    select(cell_subset, KCN, log2FC, FDR, Lecture) %>%
    arrange(FDR, desc(abs(log2FC)), KCN)

  all_results[[i]] <- res_tbl
}

pseudobulk_summary <- bind_rows(all_results) %>%
  arrange(cell_subset, FDR, desc(abs(log2FC)), KCN)

pseudobulk_significant <- pseudobulk_summary %>%
  filter(!is.na(FDR) & FDR < 0.05)

missing_rows <- tibble(
  cell_subset = NA_character_,
  KCN = missing_kcn,
  log2FC = NA_real_,
  FDR = NA_real_,
  Lecture = "Not detected in USER"
) %>%
  arrange(KCN)

comparison_info_tbl <- bind_rows(comparison_info) %>%
  left_join(cell_subset_counts, by = "cell_subset") %>%
  arrange(desc(total_cells))

write_tsv_simple(pseudobulk_summary, "USER_kcn25_pseudobulk_by_cell_subset.tsv", tab_dir)
write_tsv_simple(pseudobulk_significant, "USER_kcn25_pseudobulk_by_cell_subset_significant.tsv", tab_dir)
write_tsv_simple(missing_rows, "USER_kcn25_missing_in_user.tsv", tab_dir)
write_tsv_simple(comparison_info_tbl, "USER_kcn25_pseudobulk_comparison_info.tsv", tab_dir)

saveRDS(
  list(
    present_kcn = present_kcn,
    missing_kcn = missing_kcn,
    pseudobulk_samples = pb_sample_info,
    cell_subset_counts = cell_subset_counts,
    comparison_info = comparison_info_tbl,
    deseq_results = pseudobulk_summary,
    deseq_significant = pseudobulk_significant
  ),
  file = file.path(rds_dir, "USER_kcn25_pseudobulk_cell_subsets_results.rds")
)

message("USER KCN25 pseudobulk cell-subset analysis complete.")
