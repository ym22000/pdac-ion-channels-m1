#!/usr/bin/env Rscript

###############################################################################
# Fisher tests for KCN enrichment in the main niche groups
#
# This script tests whether each KCN is enriched in CC1+CC5 or CC2+CC3 using
# primary pancreas spots only.
#
# The test is simple: for each gene and niche, build a 2x2 table based on
# inside/outside the niche and positive/negative expression, then run Fisher's
# exact test.
###############################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(tidyr)
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

safe_fisher <- function(a, b, c, d) {
  mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  res <- suppressWarnings(fisher.test(mat))
  list(
    odds_ratio = unname(res$estimate),
    p_value = res$p.value
  )
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
tab_dir <- file.path(script_dir, "tables")

dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

object_path <- file.path(input_dir, "PDAC_Updated_ST.rds")
kcn_list_path <- normalizePath(file.path(project_dir, "..", "mycaf_icaf_correlations", "kcn_union_25.tsv"), mustWork = FALSE)

if (!file.exists(object_path)) stop("Missing input object: ", object_path)
if (!file.exists(kcn_list_path)) stop("Missing KCN list: ", kcn_list_path)

message("Loading ST Seurat object for KCN niche Fisher tests...")
st_obj <- readRDS(object_path)
st_obj <- suppressWarnings(UpdateSeuratObject(st_obj))

metadata <- st_obj@meta.data %>%
  rownames_to_column("spot") %>%
  filter(Origin == "Pancreas")

spatial_mat <- get_assay_matrix_compat(st_obj, assay = "Spatial", layer_name = "data")
kcn_list <- read.delim(kcn_list_path, stringsAsFactors = FALSE) %>%
  pull(gene) %>%
  unique()

niche_definitions <- list(
  CC1_CC5 = c("CC1", "CC5"),
  CC2_CC3 = c("CC2", "CC3")
)

available_genes <- kcn_list[kcn_list %in% rownames(spatial_mat)]
write_tsv_simple(
  data.frame(gene = kcn_list, in_spatial_assay = kcn_list %in% rownames(spatial_mat), stringsAsFactors = FALSE),
  "kcn_gene_availability.tsv",
  tab_dir
)

results_list <- list()

for (gene_name in available_genes) {
  gene_df <- data.frame(
    spot = colnames(st_obj),
    expression = as.numeric(spatial_mat[gene_name, ]),
    stringsAsFactors = FALSE
  ) %>%
    inner_join(metadata %>% select(spot, cc_ischia_10, orig.ident), by = "spot") %>%
    mutate(
      gene = gene_name,
      detected = expression > 0
    )

  for (niche_name in names(niche_definitions)) {
    niche_labels <- niche_definitions[[niche_name]]
    inside <- gene_df$cc_ischia_10 %in% niche_labels
    outside <- !inside
    positive <- gene_df$detected
    negative <- !positive

    a <- sum(inside & positive, na.rm = TRUE)
    b <- sum(outside & positive, na.rm = TRUE)
    c <- sum(inside & negative, na.rm = TRUE)
    d <- sum(outside & negative, na.rm = TRUE)

    ft <- safe_fisher(a, b, c, d)

    n_inside <- a + c
    n_outside <- b + d
    n_positive <- a + b
    n_total <- n_inside + n_outside
    expected_inside_positive <- if (n_total > 0) (n_positive * n_inside / n_total) else NA_real_

    results_list[[paste(gene_name, niche_name, sep = "__")]] <- data.frame(
      gene = gene_name,
      niche = niche_name,
      niche_labels = paste(niche_labels, collapse = ","),
      n_total_spots = n_total,
      n_inside = n_inside,
      n_outside = n_outside,
      n_positive_total = n_positive,
      n_inside_positive = a,
      n_outside_positive = b,
      pct_inside_positive = ifelse(n_inside > 0, 100 * a / n_inside, NA_real_),
      pct_outside_positive = ifelse(n_outside > 0, 100 * b / n_outside, NA_real_),
      pct_of_positive_inside_niche = ifelse(n_positive > 0, 100 * a / n_positive, NA_real_),
      niche_fraction_of_all_spots = ifelse(n_total > 0, 100 * n_inside / n_total, NA_real_),
      expected_inside_positive = expected_inside_positive,
      enrichment_ratio = ifelse(!is.na(expected_inside_positive) && expected_inside_positive > 0, a / expected_inside_positive, NA_real_),
      odds_ratio = ft$odds_ratio,
      fisher_p = ft$p_value,
      stringsAsFactors = FALSE
    )
  }
}

results_df <- bind_rows(results_list) %>%
  mutate(
    fisher_fdr_global = p.adjust(fisher_p, method = "BH")
  ) %>%
  group_by(niche) %>%
  mutate(
    fisher_fdr_within_niche = p.adjust(fisher_p, method = "BH")
  ) %>%
  ungroup() %>%
  arrange(niche, fisher_p, desc(odds_ratio))

primary_preference <- results_df %>%
  select(gene, niche, odds_ratio, fisher_p, fisher_fdr_within_niche, enrichment_ratio, pct_inside_positive, pct_outside_positive) %>%
  tidyr::pivot_wider(
    names_from = niche,
    values_from = c(odds_ratio, fisher_p, fisher_fdr_within_niche, enrichment_ratio, pct_inside_positive, pct_outside_positive),
    names_sep = "__"
  ) %>%
  mutate(
    favored_niche = dplyr::case_when(
      odds_ratio__CC1_CC5 > odds_ratio__CC2_CC3 ~ "CC1_CC5",
      odds_ratio__CC2_CC3 > odds_ratio__CC1_CC5 ~ "CC2_CC3",
      TRUE ~ "tie"
    )
  )

write_tsv_simple(results_df, "kcn_niche_fisher_tests_primary.tsv", tab_dir)
write_tsv_simple(primary_preference, "kcn_preferred_niche_primary.tsv", tab_dir)


message("KCN niche Fisher tests complete.")


