#!/usr/bin/env Rscript

###############################################################################
# Marker ecotype-vs-rest contrasts in primary PDAC sections
#
# Section-wise Wilcoxon contrasts for selected stromal / myCAF-like markers:
# - ACTA2
# - CCN2
# - POSTN
# - TAGLN
# - COL12A1
#
# Main contrast groups:
# - CC1+CC5 vs rest
# - CC2+CC3 vs rest
###############################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

options(stringsAsFactors = FALSE)

MIN_SPOTS_PER_GROUP <- 20
LOG2FC_PSEUDOCOUNT <- 1
TARGET_GENES <- c("ACTA2", "CCN2", "POSTN", "TAGLN", "COL12A1")

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  file_path <- sub(file_arg, "", args[grep(file_arg, args)])
  if (length(file_path) == 0) {
    return(normalizePath(getwd()))
  }
  normalizePath(dirname(file_path))
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

find_project_dir <- function(start_dir) {
  current_dir <- normalizePath(start_dir)
  repeat {
    if (file.exists(file.path(current_dir, "inputs", "PDAC_Updated_ST.rds"))) {
      return(current_dir)
    }
    parent_dir <- dirname(current_dir)
    if (identical(parent_dir, current_dir)) {
      stop("Could not find the ST project directory.")
    }
    current_dir <- parent_dir
  }
}

safe_auc_from_wilcox <- function(x, y) {
  n1 <- length(x)
  n2 <- length(y)
  if (n1 == 0 || n2 == 0) {
    return(NA_real_)
  }

  combined <- c(x, y)
  ranks <- rank(combined, ties.method = "average")
  rank_sum_x <- sum(ranks[seq_len(n1)])
  u_stat <- rank_sum_x - (n1 * (n1 + 1) / 2)
  u_stat / (n1 * n2)
}

summarize_section_contrast <- function(x_in, x_out) {
  positive_in <- x_in > 0
  positive_out <- x_out > 0

  mean_in <- mean(x_in)
  mean_out <- mean(x_out)
  mean_pos_in <- if (any(positive_in)) mean(x_in[positive_in]) else NA_real_
  mean_pos_out <- if (any(positive_out)) mean(x_out[positive_out]) else NA_real_

  wt <- suppressWarnings(
    wilcox.test(x_in, x_out, exact = FALSE, correct = FALSE)
  )

  auc <- safe_auc_from_wilcox(x_in, x_out)

  list(
    n_inside = length(x_in),
    n_outside = length(x_out),
    n_inside_positive = sum(positive_in),
    n_outside_positive = sum(positive_out),
    pct_inside_positive = 100 * mean(positive_in),
    pct_outside_positive = 100 * mean(positive_out),
    mean_expression_inside = mean_in,
    mean_expression_outside = mean_out,
    mean_positive_expression_inside = mean_pos_in,
    mean_positive_expression_outside = mean_pos_out,
    delta_pct_positive = 100 * mean(positive_in) - 100 * mean(positive_out),
    log2fc_mean_expression = log2((mean_in + LOG2FC_PSEUDOCOUNT) / (mean_out + LOG2FC_PSEUDOCOUNT)),
    auc = auc,
    wilcoxon_p = wt$p.value
  )
}

script_dir <- get_script_dir()
project_dir <- find_project_dir(script_dir)
input_dir <- file.path(project_dir, "inputs")
tab_dir <- file.path(script_dir, "tables")

dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

object_path <- file.path(input_dir, "PDAC_Updated_ST.rds")

if (!file.exists(object_path)) stop("Missing input object: ", object_path)

message("Loading ST Seurat object for marker ecotype-vs-rest contrasts...")
st_obj <- readRDS(object_path)
st_obj <- suppressWarnings(UpdateSeuratObject(st_obj))

metadata <- st_obj@meta.data %>%
  rownames_to_column("spot") %>%
  filter(Origin == "Pancreas")

spatial_mat <- get_assay_matrix_compat(st_obj, assay = "Spatial", layer_name = "data")

contrast_definitions <- list(
  CC1_CC5 = c("CC1", "CC5"),
  CC2_CC3 = c("CC2", "CC3")
)

gene_availability <- data.frame(
  gene = TARGET_GENES,
  in_spatial_assay = TARGET_GENES %in% rownames(spatial_mat),
  stringsAsFactors = FALSE
)
write_tsv_simple(gene_availability, "marker_gene_availability.tsv", tab_dir)

section_contrast_availability <- bind_rows(lapply(
  names(contrast_definitions),
  function(contrast_name) {
    labels <- contrast_definitions[[contrast_name]]
    metadata %>%
      group_by(orig.ident) %>%
      summarize(
        contrast = contrast_name,
        n_inside = sum(cc_ischia_10 %in% labels, na.rm = TRUE),
        n_outside = sum(!(cc_ischia_10 %in% labels), na.rm = TRUE),
        usable = n_inside >= MIN_SPOTS_PER_GROUP & n_outside >= MIN_SPOTS_PER_GROUP,
        .groups = "drop"
      )
  }
))
write_tsv_simple(section_contrast_availability, "section_contrast_availability.tsv", tab_dir)

results_list <- list()

for (section_name in sort(unique(metadata$orig.ident))) {
  section_meta <- metadata %>% filter(orig.ident == section_name)
  if (nrow(section_meta) == 0) next

  section_spots <- section_meta$spot

  for (contrast_name in names(contrast_definitions)) {
    labels <- contrast_definitions[[contrast_name]]
    inside_mask <- section_meta$cc_ischia_10 %in% labels
    outside_mask <- !inside_mask

    if (sum(inside_mask, na.rm = TRUE) < MIN_SPOTS_PER_GROUP || sum(outside_mask, na.rm = TRUE) < MIN_SPOTS_PER_GROUP) {
      next
    }

    for (gene_name in TARGET_GENES) {
      if (!(gene_name %in% rownames(spatial_mat))) next

      expr_values <- as.numeric(spatial_mat[gene_name, section_spots])
      x_in <- expr_values[inside_mask]
      x_out <- expr_values[outside_mask]

      metrics <- summarize_section_contrast(x_in, x_out)

      results_list[[paste(section_name, contrast_name, gene_name, sep = "__")]] <- data.frame(
        orig.ident = section_name,
        contrast = contrast_name,
        contrast_labels = paste(labels, collapse = ","),
        gene = gene_name,
        n_inside = metrics$n_inside,
        n_outside = metrics$n_outside,
        n_inside_positive = metrics$n_inside_positive,
        n_outside_positive = metrics$n_outside_positive,
        pct_inside_positive = metrics$pct_inside_positive,
        pct_outside_positive = metrics$pct_outside_positive,
        delta_pct_positive = metrics$delta_pct_positive,
        mean_expression_inside = metrics$mean_expression_inside,
        mean_expression_outside = metrics$mean_expression_outside,
        mean_positive_expression_inside = metrics$mean_positive_expression_inside,
        mean_positive_expression_outside = metrics$mean_positive_expression_outside,
        log2fc_mean_expression = metrics$log2fc_mean_expression,
        auc = metrics$auc,
        wilcoxon_p = metrics$wilcoxon_p,
        stringsAsFactors = FALSE
      )
    }
  }
}

results_df <- bind_rows(results_list) %>%
  group_by(orig.ident, contrast) %>%
  mutate(
    wilcoxon_fdr = p.adjust(wilcoxon_p, method = "BH")
  ) %>%
  ungroup() %>%
  arrange(contrast, orig.ident, wilcoxon_p, desc(log2fc_mean_expression))

write_tsv_simple(results_df, "marker_ecotype_vs_rest_by_section.tsv", tab_dir)

summary_df <- results_df %>%
  group_by(gene, contrast, contrast_labels) %>%
  summarize(
    n_sections_tested = n(),
    n_sections_log2fc_positive = sum(log2fc_mean_expression > 0, na.rm = TRUE),
    n_sections_auc_gt_0_5 = sum(auc > 0.5, na.rm = TRUE),
    n_sections_fdr_lt_0_05 = sum(wilcoxon_fdr < 0.05, na.rm = TRUE),
    median_log2fc_mean_expression = median(log2fc_mean_expression, na.rm = TRUE),
    median_auc = median(auc, na.rm = TRUE),
    median_delta_pct_positive = median(delta_pct_positive, na.rm = TRUE),
    median_pct_inside_positive = median(pct_inside_positive, na.rm = TRUE),
    median_pct_outside_positive = median(pct_outside_positive, na.rm = TRUE),
    median_mean_expression_inside = median(mean_expression_inside, na.rm = TRUE),
    median_mean_expression_outside = median(mean_expression_outside, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(contrast, desc(median_log2fc_mean_expression), desc(median_auc))

write_tsv_simple(summary_df, "marker_ecotype_vs_rest_summary.tsv", tab_dir)


message("Marker ecotype-vs-rest contrasts complete.")
