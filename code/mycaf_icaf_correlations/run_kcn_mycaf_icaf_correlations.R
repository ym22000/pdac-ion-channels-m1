# KCN - myCAF/iCAF correlation analysis
#
# This script computes pseudobulk Spearman correlations between the KCN list
# and human myCAF/iCAF marker genes in the stromal single-cell dataset.
#
# The main idea is to work at the patient-level pseudobulk scale instead of at
# the single-cell scale. That keeps the analysis less noisy and reduces the
# pseudoreplication problem.
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  file_path <- sub(file_arg, "", args[grep(file_arg, args)])
  if (length(file_path) == 0) {
    return(normalizePath(getwd()))
  }
  normalizePath(dirname(file_path))
}

parse_args <- function(defaults) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0) {
    return(defaults)
  }

  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key)
    }
    key <- sub("^--", "", key)
    if (!key %in% names(defaults)) {
      stop("Unknown argument: --", key)
    }
    if (i == length(args)) {
      stop("Missing value for argument: --", key)
    }
    defaults[[key]] <- args[[i + 1]]
    i <- i + 2
  }
  defaults
}

read_gene_list <- function(path) {
  genes <- readLines(path, warn = FALSE)
  genes <- trimws(genes)
  genes <- genes[nzchar(genes)]
  unique(genes)
}

read_kcn_list <- function(path) {
  x <- read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"gene" %in% colnames(x)) {
    stop("KCN input file must contain a column named 'gene': ", path)
  }
  unique(trimws(x$gene[nzchar(trimws(x$gene))]))
}

robust_z <- function(x) {
  s <- stats::sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) {
    return(rep(0, length(x)))
  }
  as.numeric((x - mean(x, na.rm = TRUE)) / s)
}

make_signature_score <- function(expr_df, genes) {
  present <- intersect(genes, colnames(expr_df))
  if (length(present) == 0) {
    return(rep(NA_real_, nrow(expr_df)))
  }
  z_mat <- vapply(present, function(g) robust_z(expr_df[[g]]), numeric(nrow(expr_df)))
  if (is.vector(z_mat)) {
    z_mat <- matrix(z_mat, ncol = 1)
  }
  apply(z_mat, 1, stats::median, na.rm = TRUE)
}

aggregate_counts_by_sample <- function(counts, sample_ids) {
  sample_levels <- unique(sample_ids)
  model <- sparse.model.matrix(~ 0 + factor(sample_ids, levels = sample_levels))
  colnames(model) <- sample_levels
  aggregated <- counts %*% model
  aggregated
}

counts_to_logcpm <- function(count_mat) {
  lib_size <- Matrix::colSums(count_mat)
  cpm <- t(t(count_mat) / lib_size) * 1e6
  log1p(cpm)
}

safe_spearman <- function(x, y, min_n = 6L) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < min_n || length(unique(x)) < 2 || length(unique(y)) < 2) {
    return(list(rho = NA_real_, p_value = NA_real_, n = length(x)))
  }
  test <- suppressWarnings(stats::cor.test(x, y, method = "spearman", exact = FALSE))
  list(rho = unname(test$estimate), p_value = unname(test$p.value), n = length(x))
}

safe_wilcoxon <- function(df) {
  group_sizes <- table(df$kcn_group)
  if (!all(c("Low", "High") %in% names(group_sizes))) {
    return(c(statistic = NA_real_, p_value = NA_real_))
  }
  if (any(group_sizes < 3)) {
    return(c(statistic = NA_real_, p_value = NA_real_))
  }
  test <- suppressWarnings(stats::wilcox.test(score ~ kcn_group, data = df, exact = FALSE))
  c(statistic = unname(test$statistic), p_value = unname(test$p.value))
}

project_script_dir <- get_script_dir()
project_root <- normalizePath(file.path(project_script_dir, "..", ".."))

defaults <- list(
  seurat = file.path(project_root, "data", "Stroma_Subset2021_fibroblast_focus.rds"),
  kcn = file.path(project_script_dir, "kcn_union_25.tsv"),
  mycaf = file.path(project_script_dir, "genes_top30_myCAF_human.txt"),
  icaf = file.path(project_script_dir, "genes_top30_iCAF_human.txt"),
  outdir = file.path(project_script_dir, "results"),
  min_cells = "20"
)

args <- parse_args(defaults)
args$seurat <- normalizePath(args$seurat, mustWork = TRUE)
args$kcn <- normalizePath(args$kcn, mustWork = TRUE)
args$mycaf <- normalizePath(args$mycaf, mustWork = TRUE)
args$icaf <- normalizePath(args$icaf, mustWork = TRUE)
args$outdir <- normalizePath(args$outdir, mustWork = FALSE)
min_cells <- as.integer(args$min_cells)

dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)
fig_dir <- file.path(args$outdir, "figures")
tab_dir <- file.path(args$outdir, "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading inputs...")
obj <- readRDS(args$seurat)
kcn_genes <- read_kcn_list(args$kcn)
mycaf_markers <- read_gene_list(args$mycaf)
icaf_markers <- read_gene_list(args$icaf)
all_markers <- unique(c(mycaf_markers, icaf_markers))

if (!all(c("Patient", "CellType2") %in% colnames(obj@meta.data))) {
  stop("Seurat object must contain 'Patient' and 'CellType2' metadata columns.")
}

meta <- obj@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  transmute(
    cell_id = cell_id,
    Patient = as.character(Patient),
    CellType2 = as.character(CellType2),
    sample_id = paste(Patient, CellType2, sep = "__")
  )

sample_qc <- meta %>%
  count(sample_id, Patient, CellType2, name = "n_cells") %>%
  filter(n_cells >= min_cells)

if (nrow(sample_qc) < 6) {
  stop("Too few pseudobulk samples after filtering. Reduce --min_cells or inspect the object.")
}

keep_cells <- meta$cell_id[meta$sample_id %in% sample_qc$sample_id]
meta_keep <- meta %>% filter(cell_id %in% keep_cells)

rna_counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
genes_needed <- unique(c(kcn_genes, all_markers))
genes_present <- intersect(genes_needed, rownames(rna_counts))
missing_genes <- setdiff(genes_needed, genes_present)

counts_subset <- rna_counts[genes_present, meta_keep$cell_id, drop = FALSE]
aggregated_counts <- aggregate_counts_by_sample(counts_subset, meta_keep$sample_id)
logcpm <- counts_to_logcpm(aggregated_counts)

pseudobulk_expr <- as.data.frame(t(as.matrix(logcpm)))
pseudobulk_expr$sample_id <- rownames(pseudobulk_expr)
pseudobulk <- sample_qc %>%
  left_join(pseudobulk_expr, by = "sample_id")

pseudobulk$myCAF_score <- make_signature_score(pseudobulk, mycaf_markers)
pseudobulk$iCAF_score <- make_signature_score(pseudobulk, icaf_markers)

qc_table <- tibble::tibble(
  gene = genes_needed,
  category = dplyr::case_when(
    gene %in% kcn_genes ~ "KCN",
    gene %in% mycaf_markers ~ "myCAF_marker",
    gene %in% icaf_markers ~ "iCAF_marker",
    TRUE ~ "other"
  ),
  present_in_expression = gene %in% genes_present
) %>%
  arrange(category, desc(present_in_expression), gene)

message("Computing KCN-marker correlations...")
cor_results <- vector("list", length(kcn_genes) * length(all_markers))
idx <- 1L

for (kcn in kcn_genes) {
  if (!kcn %in% colnames(pseudobulk)) {
    next
  }
  for (marker in all_markers) {
    if (!marker %in% colnames(pseudobulk)) {
      next
    }
    stats_out <- safe_spearman(pseudobulk[[kcn]], pseudobulk[[marker]])
    cor_results[[idx]] <- tibble::tibble(
      kcn_gene = kcn,
      marker_gene = marker,
      marker_group = ifelse(marker %in% mycaf_markers, "myCAF", "iCAF"),
      n_samples = stats_out$n,
      spearman_rho = stats_out$rho,
      p_value = stats_out$p_value
    )
    idx <- idx + 1L
  }
}

cor_table <- dplyr::bind_rows(cor_results) %>%
  filter(!is.na(marker_gene)) %>%
  mutate(
    p_adj_bh = p.adjust(p_value, method = "BH")
  ) %>%
  arrange(marker_group, kcn_gene, desc(abs(spearman_rho)))

summary_table <- cor_table %>%
  group_by(kcn_gene, marker_group) %>%
  summarise(
    n_markers_tested = sum(!is.na(spearman_rho)),
    mean_rho = mean(spearman_rho, na.rm = TRUE),
    median_rho = median(spearman_rho, na.rm = TRUE),
    mean_abs_rho = mean(abs(spearman_rho), na.rm = TRUE),
    min_p_adj_bh = min(p_adj_bh, na.rm = TRUE),
    n_fdr_lt_0_05 = sum(p_adj_bh < 0.05, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(marker_group, desc(mean_abs_rho), kcn_gene)

message("Building pseudobulk score barplot...")
score_long <- pseudobulk %>%
  select(sample_id, Patient, CellType2, n_cells, myCAF_score, iCAF_score, all_of(intersect(kcn_genes, colnames(pseudobulk)))) %>%
  pivot_longer(cols = c(myCAF_score, iCAF_score), names_to = "signature", values_to = "score")

violin_list <- vector("list", length(kcn_genes))
comp_list <- vector("list", length(kcn_genes) * 2L)
vidx <- 1L

for (kcn in intersect(kcn_genes, colnames(pseudobulk))) {
  kcn_expr <- pseudobulk[[kcn]]
  if (all(!is.finite(kcn_expr)) || length(unique(kcn_expr[is.finite(kcn_expr)])) < 2) {
    next
  }
  med <- stats::median(kcn_expr, na.rm = TRUE)
  group <- ifelse(kcn_expr > med, "High", "Low")
  if (length(unique(group)) < 2) {
    next
  }

  df_kcn <- score_long %>%
    filter(sample_id %in% pseudobulk$sample_id) %>%
    mutate(
      kcn_gene = kcn,
      kcn_expression = rep(kcn_expr, each = 2),
      kcn_group = rep(group, each = 2)
    )

  violin_list[[vidx]] <- df_kcn
  vidx <- vidx + 1L

  for (sig in c("myCAF_score", "iCAF_score")) {
    sig_df <- df_kcn %>% filter(signature == sig)
    test_out <- safe_wilcoxon(sig_df)
    comp_list[[length(Filter(Negate(is.null), comp_list)) + 1L]] <- tibble::tibble(
      kcn_gene = kcn,
      signature = sig,
      n_samples = nrow(sig_df),
      n_low = sum(sig_df$kcn_group == "Low"),
      n_high = sum(sig_df$kcn_group == "High"),
      wilcoxon_statistic = test_out[["statistic"]],
      p_value = test_out[["p_value"]]
    )
  }
}

violin_table <- bind_rows(violin_list) %>%
  mutate(
    signature = recode(signature, myCAF_score = "myCAF score", iCAF_score = "iCAF score"),
    kcn_group = factor(kcn_group, levels = c("Low", "High"))
  )

comparison_table <- bind_rows(comp_list) %>%
  filter(!is.na(kcn_gene)) %>%
  mutate(
    p_adj_bh = p.adjust(p_value, method = "BH"),
    signature = recode(signature, myCAF_score = "myCAF score", iCAF_score = "iCAF score")
  )

annotation_table <- cor_table %>%
  filter(marker_group %in% c("myCAF", "iCAF")) %>%
  group_by(kcn_gene, marker_group) %>%
  summarise(
    mean_rho = mean(spearman_rho, na.rm = TRUE),
    min_fdr = min(p_adj_bh, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    signature = ifelse(marker_group == "myCAF", "myCAF score", "iCAF score"),
    label = paste0("mean rho=", sprintf("%.2f", mean_rho), "\nmin FDR=", formatC(min_fdr, format = "e", digits = 2))
  )

y_positions <- violin_table %>%
  group_by(kcn_gene, signature) %>%
  summarise(y_pos = max(score, na.rm = TRUE) + 0.15, .groups = "drop")

annotation_table <- annotation_table %>%
  left_join(y_positions, by = c("kcn_gene", "signature"))

bar_plot <- ggplot(violin_table, aes(x = signature, y = score, fill = kcn_group)) +
  stat_summary(
    fun = mean,
    geom = "col",
    position = position_dodge(width = 0.8),
    width = 0.72,
    alpha = 0.9,
    color = NA
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    position = position_dodge(width = 0.8),
    width = 0.18,
    linewidth = 0.35,
    color = "black"
  ) +
  geom_text(
    data = annotation_table,
    aes(x = signature, y = y_pos, label = label),
    inherit.aes = FALSE,
    size = 2.2,
    lineheight = 0.95,
    vjust = 0
  ) +
  facet_wrap(~ kcn_gene, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = c("Low" = "#b9c4cf", "High" = "#800E13")) +
  labs(
    title = "myCAF and iCAF pseudobulk signature scores across KCN-defined groups",
    subtitle = "Pseudobulk samples are defined as Patient Ã— CellType2; annotations summarize KCN-marker correlations.",
    x = NULL,
    y = "Signature score (median z-scored top30 markers)"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    axis.text.x = element_text(angle = 22, hjust = 1),
    legend.title = element_blank(),
    legend.position = "top"
  )

message("Writing moffitt_stromal_exploration_outputs...")
write.table(cor_table, file.path(tab_dir, "kcn_marker_spearman_correlations.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(summary_table, file.path(tab_dir, "kcn_marker_correlation_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(qc_table, file.path(tab_dir, "gene_filtering_qc.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(pseudobulk, file.path(tab_dir, "pseudobulk_expression_matrix.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(comparison_table, file.path(tab_dir, "kcn_signature_group_comparisons.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
ggsave(
  filename = file.path(fig_dir, "Figure_01_KCN_myCAF_iCAF_signature_barplot.png"),
  plot = bar_plot,
  width = 18,
  height = 14,
  dpi = 300
)

message("Done.")

