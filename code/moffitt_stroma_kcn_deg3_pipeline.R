###############################################################################
# Moffitt stroma KCN DEG pipeline with 3 methods
#
# This script is the main DEG workflow for the fibroblast-focused Moffitt atlas.
# It combines pseudobulk and cell-level approaches to define the KCN candidates
# linked to the stromal states of interest.
###############################################################################
required_packages <- c(
  "Seurat",
  "dplyr",
  "ggplot2",
  "patchwork",
  "tibble",
  "tidyr",
  "purrr",
  "Matrix",
  "DESeq2",
  "VennDiagram",
  "ggrepel",
  "lme4"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

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
  library(tidyr)
  library(purrr)
  library(Matrix)
  library(DESeq2)
  library(VennDiagram)
  library(ggrepel)
  library(lme4)
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

wrap_text <- function(x, width = 44) {
  paste(strwrap(x, width = width), collapse = "\n")
}

collapse_unique <- function(x) {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) {
    return(NA_character_)
  }
  paste(sort(x), collapse = "; ")
}

make_placeholder_plot <- function(title_text, body_text) {
  ggplot() +
    annotate("text", x = 0, y = 0, label = paste(title_text, body_text, sep = "\n\n"), size = 5) +
    theme_void()
}

scale_numeric_covariate <- function(x) {
  x <- log10(x + 1)
  if (length(unique(x)) <= 1) {
    return(rep(0, length(x)))
  }
  as.numeric(scale(x))
}

run_deseq2 <- function(
  count_matrix,
  coldata_df,
  coef_name,
  method_label,
  comparison_label,
  subtype_label,
  design_formula
) {
  coldata_df <- as.data.frame(coldata_df)
  coldata_df$n_cells_raw <- coldata_df$n_cells
  coldata_df$n_cells <- scale_numeric_covariate(coldata_df$n_cells)

  dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = coldata_df,
    design = design_formula
  )

  dds <- dds[rowSums(counts(dds)) >= 10, ]

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
    fallback_size_factors <- library_sizes / exp(mean(log(positive_sizes)))
    fallback_size_factors[!is.finite(fallback_size_factors)] <- 1
    sizeFactors(dds) <- fallback_size_factors
  }

  dds <- DESeq(dds, quiet = TRUE, sfType = "poscounts", fitType = "local")

  res <- results(dds, name = coef_name, independentFiltering = TRUE) %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    mutate(
      method = method_label,
      comparison = comparison_label,
      CellType2 = subtype_label,
      padj_BH = padj,
      bh_significant = !is.na(padj_BH) & padj_BH < 0.05,
      direction = case_when(
        log2FoldChange > 0 ~ "positive",
        log2FoldChange < 0 ~ "negative",
        TRUE ~ "neutral"
      )
    ) %>%
    arrange(padj_BH, desc(abs(log2FoldChange)))

  list(dds = dds, results = res)
}

capture_warnings <- function(expr) {
  warnings <- character(0)
  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = unique(warnings))
}

fit_detection_model <- function(df_gene) {
  full_fit <- capture_warnings(
    glmer(
      detected ~ is_target + log10_nCount_RNA + Dataset + (1 | Patient),
      data = df_gene,
      family = binomial(),
      control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
    )
  )

  red_fit <- capture_warnings(
    glmer(
      detected ~ log10_nCount_RNA + Dataset + (1 | Patient),
      data = df_gene,
      family = binomial(),
      control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
    )
  )

  an <- suppressWarnings(anova(red_fit$value, full_fit$value, test = "Chisq"))
  coef_tab <- summary(full_fit$value)$coefficients

  beta <- unname(coef_tab["is_targetTRUE", "Estimate"])
  se <- unname(coef_tab["is_targetTRUE", "Std. Error"])
  zval <- unname(coef_tab["is_targetTRUE", "z value"])

  list(
    beta = beta,
    se = se,
    z = zval,
    pval_lrt = an$`Pr(>Chisq)`[2],
    odds_ratio = exp(beta),
    warnings = unique(c(full_fit$warnings, red_fit$warnings))
  )
}

fit_positive_model <- function(df_gene_pos) {
  full_fit <- capture_warnings(
    lmer(
      expr_pos ~ is_target + log10_nCount_RNA + Dataset + (1 | Patient),
      data = df_gene_pos,
      REML = FALSE,
      control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
    )
  )

  red_fit <- capture_warnings(
    lmer(
      expr_pos ~ log10_nCount_RNA + Dataset + (1 | Patient),
      data = df_gene_pos,
      REML = FALSE,
      control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
    )
  )

  an <- suppressWarnings(anova(red_fit$value, full_fit$value))
  coef_tab <- summary(full_fit$value)$coefficients

  beta <- unname(coef_tab["is_targetTRUE", "Estimate"])
  se <- unname(coef_tab["is_targetTRUE", "Std. Error"])
  tval <- unname(coef_tab["is_targetTRUE", "t value"])

  list(
    beta = beta,
    se = se,
    t = tval,
    pval_lrt = an$`Pr(>Chisq)`[2],
    warnings = unique(c(full_fit$warnings, red_fit$warnings))
  )
}

pick_top_genes <- function(df, group_col, gene_col = "gene", n_per_group = 4, score_col = "padj_BH") {
  group_col <- rlang::ensym(group_col)
  gene_col <- rlang::ensym(gene_col)
  score_col <- rlang::ensym(score_col)

  df %>%
    group_by(!!group_col) %>%
    arrange(!!score_col, desc(abs(log2FoldChange %||% 0)), .by_group = TRUE) %>%
    slice_head(n = n_per_group) %>%
    ungroup() %>%
    pull(!!gene_col) %>%
    unique()
}

make_dotplot <- function(seu, features, group_by, colors_use, title_text, file_path, width = 12, height = 6) {
  if (length(features) == 0) {
    save_plot(make_placeholder_plot(title_text, "No KCN gene passed the current selection."), file_path, width, height)
    return(invisible(NULL))
  }

  p <- suppressWarnings(
    DotPlot(
      object = seu,
      features = unique(features),
      group.by = group_by,
      cols = c("grey92", "#FF5400"),
      dot.scale = 10
    ) +
      scale_size(
        range = c(5, 14),
        breaks = c(5, 25, 50, 75, 100),
        limits = c(0, 100)
      ) +
      scale_color_gradient(low = "grey92", high = "#FF5400") +
      labs(title = title_text, x = NULL, y = NULL) +
      clean_classic_theme() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(face = "bold")
      )
  )

  save_plot(p, file_path, width = width, height = height)
}

make_volcano_plot <- function(
  df,
  x_col,
  y_col,
  p_threshold = 0.05,
  lfc_threshold = 0.25,
  facet_col = NULL,
  color_map = NULL,
  label_df = NULL,
  title_text,
  x_label,
  y_label,
  file_path,
  width = 12,
  height = 8
) {
  p <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_hline(yintercept = -log10(p_threshold), linetype = "dashed", color = "grey60") +
    geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype = "dashed", color = "grey60")

  if ("plot_group" %in% colnames(df)) {
    p <- p +
      geom_point(
        data = df %>% filter(plot_group == "NS"),
        color = "grey75",
        alpha = 0.35,
        size = 1.2
      )

    if ("Significant_non_KCN" %in% df$plot_group) {
      if (!is.null(color_map) && "plot_color" %in% colnames(df)) {
        p <- p +
          geom_point(
            data = df %>% filter(plot_group == "Significant_non_KCN"),
            aes(color = plot_color),
            alpha = 0.8,
            size = 1.5
          ) +
          scale_color_manual(values = color_map)
      } else {
        p <- p +
          geom_point(
            data = df %>% filter(plot_group == "Significant_non_KCN"),
            color = "#FF5400",
            alpha = 0.8,
            size = 1.5
          )
      }
    }

    p <- p +
      geom_point(
        data = df %>% filter(plot_group == "KCN"),
        color = "black",
        alpha = 0.95,
        size = 1.7
      )
  }

  if (!is.null(label_df) && nrow(label_df) > 0) {
    p <- p +
      geom_text_repel(
        data = label_df,
        aes(label = gene),
        color = "black",
        size = 3,
        box.padding = 0.3,
        point.padding = 0.2,
        max.overlaps = 60,
        show.legend = FALSE
      )
  }

  if (!is.null(facet_col)) {
    p <- p + facet_wrap(stats::as.formula(paste("~", facet_col)), scales = "free")
  }

  p <- p +
    labs(title = title_text, x = x_label, y = y_label) +
    clean_classic_theme() +
    theme(legend.position = "none")

  save_plot(p, file_path, width = width, height = height)
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

###############################################################################
# 2. Paths, colors, and constants
###############################################################################

project_root <- find_project_root()
analysis_root <- file.path(project_root, "M1_Bioinformatics_Ion_Channel", "code")
project_data_dir <- file.path(project_root, "M1_Bioinformatics_Ion_Channel", "data")
input_rds <- file.path(project_data_dir, "Stroma_Subset2021_fibroblast_focus.rds")

output_root <- file.path(analysis_root, "moffitt_kcn_deg_3_methods")
figure_dir <- file.path(output_root, "figures")
table_dir <- file.path(output_root, "tables")
rds_dir <- file.path(output_root, "rds")

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

subtype_levels <- c("qPSC", "smPSC", "myCAF", "csCAF", "iCAF")
celltype_colors <- c(
  "qPSC" = "#70D6FF",
  "smPSC" = "#0096C7",
  "myCAF" = "#FF686B",
  "csCAF" = "#FF9770",
  "iCAF" = "#FF5400"
)
state_colors <- c(
  "Activated" = "#FF686B",
  "Normal" = "#70D6FF"
)
venn_colors <- c(
  "Test1_Activated_vs_Normal_DESeq2" = "#ACFFAC",
  "Test2_Subtype_OneVsRest_DESeq2" = "#FFAAAB",
  "Test3_Subtype_Hurdle_Mixed" = "#FDFEAF"
)

min_cells_per_state <- 20
min_cells_target <- c(qPSC = 20, smPSC = 20, myCAF = 20, csCAF = 20, iCAF = 10)
min_cells_rest <- 20
min_patients_required <- 5
padj_cutoff <- 0.05
logfc_cutoff <- 0.25
min_positive_cells_total <- 20
min_positive_patients <- 3
min_pct_overall <- 1
min_pct_subtype <- 1
min_n_subtypes <- 2

###############################################################################
# 3. Load object and metadata
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
  "nCount_RNA",
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

stroma$CellType2 <- factor(as.character(stroma$CellType2), levels = subtype_levels)
stroma$StromaState <- ifelse(
  stroma$Moffitt.F5_ActivatedStroma.top25 > stroma$Moffitt.F13_NormalStroma.top25,
  "Activated",
  "Normal"
)
stroma$StromaState <- factor(stroma$StromaState, levels = c("Normal", "Activated"))

meta <- stroma[[]] %>%
  rownames_to_column("cell_id") %>%
  filter(CellType2 %in% subtype_levels)

counts_mat <- GetAssayData(stroma, assay = "RNA", layer = "counts")
expr_mat <- GetAssayData(stroma, assay = "RNA", layer = "data")
kcn_genes <- grep("^KCN", rownames(stroma), value = TRUE)

###############################################################################
# 4. UMAP of retained stromal subtypes
###############################################################################

p_umap_subtypes <- DimPlot(
  object = stroma,
  reduction = "umap",
  group.by = "CellType2",
  cols = celltype_colors,
  pt.size = 0.55,
  label = FALSE
) +
  labs(title = "UMAP of retained stromal fibroblast subtypes") +
  clean_embedding_theme()

save_plot(
  p_umap_subtypes,
  file.path(figure_dir, "Figure_01_UMAP_retained_stromal_subtypes.png"),
  width = 9,
  height = 7
)

###############################################################################
# 5. Shared KCN prefilter for subtype-aware methods
###############################################################################

idx_list <- split(seq_len(ncol(counts_mat)), as.character(meta$CellType2))
pct_by_ct <- sapply(idx_list, function(idx) {
  Matrix::rowSums(counts_mat[kcn_genes, idx, drop = FALSE] > 0) / length(idx) * 100
})

if (is.vector(pct_by_ct)) {
  pct_by_ct <- matrix(pct_by_ct, ncol = 1)
  colnames(pct_by_ct) <- names(idx_list)
}

pct_overall <- Matrix::rowSums(counts_mat[kcn_genes, , drop = FALSE] > 0) / ncol(counts_mat) * 100
n_subtypes_ge1 <- rowSums(pct_by_ct[, colnames(pct_by_ct) %in% subtype_levels, drop = FALSE] >= min_pct_subtype)

kcn_prefilter_tbl <- data.frame(
  gene = kcn_genes,
  pct_cells_overall = round(pct_overall, 3),
  n_subtypes_ge1pct = n_subtypes_ge1,
  keep_for_subtype_methods = pct_overall >= min_pct_overall & n_subtypes_ge1 >= min_n_subtypes,
  stringsAsFactors = FALSE
)

for (ct in colnames(pct_by_ct)) {
  kcn_prefilter_tbl[[paste0("pct_", ct)]] <- round(pct_by_ct[, ct], 3)
}

filtered_kcn <- kcn_prefilter_tbl$gene[kcn_prefilter_tbl$keep_for_subtype_methods]
save_tsv(kcn_prefilter_tbl, file.path(table_dir, "KCN_prefilter_table.tsv"))

###############################################################################
# 6. Test 1 - DESeq2 pseudobulk Activated vs Normal
###############################################################################
# Pseudobulk rationale:
# - Pseudobulk aggregation is recommended for multi-sample single-cell
#   differential expression because it respects the biological replicate
#   structure and reduces pseudoreplication bias
#   (Zimmerman et al., Nat Commun 2021; Murphy and Skene, Nat Commun 2022).
# - Adding n_cells as a covariate helps control for unequal sampling depth and
#   unequal cell representation across pseudo-samples.
# - This n_cells term is a pragmatic sample-level adjustment rather than a
#   single universally mandated formula. The underlying rationale comes from
#   the pseudobulk literature showing that the number of aggregated cells
#   affects power and stability in multi-sample single-cell studies
#   (Schmid et al., Nat Commun 2021; Zimmerman et al., Nat Commun 2021;
#   Murphy and Skene, Nat Commun 2022).
#
# KCN rationale:
# - KCN genes often show moderate but consistent state-dependent effects.
# - A patient-aware pseudobulk model is therefore useful to recover these
#   signals without over-interpreting cell-level sampling noise.
#
# Multiple-testing rule:
# - BH-adjusted p-values are used with a threshold of 0.05.
# - BH controls the false discovery rate and is more appropriate than
#   Bonferroni in high-dimensional transcriptomic data, where Bonferroni would
#   be too conservative and could miss biologically relevant KCN genes.
#
# Input filtering:
# - keep pseudo-samples with at least 20 cells per stromal state
# - keep only patients represented in both Activated and Normal states
# - use ~ Patient + n_cells + condition so that the direction of log2FC is
#   directly interpretable as Activated relative to Normal
#
# Practical summary:
# - unit of analysis: one pseudo-sample per patient and stromal state
# - tested coefficient: Activated versus Normal
# - exported moffitt_stromal_exploration_outputs: full DEG table + KCN-only table + volcano + dotplot

test1_samples <- meta %>%
  dplyr::count(Patient, Dataset, StromaState, name = "n_cells") %>%
  filter(n_cells >= min_cells_per_state) %>%
  group_by(Patient) %>%
  filter(all(c("Normal", "Activated") %in% StromaState)) %>%
  ungroup() %>%
  mutate(sample_id = paste(Patient, StromaState, sep = "___"))

test1_counts_list <- list()
for (i in seq_len(nrow(test1_samples))) {
  current_cells <- meta$cell_id[
    meta$Patient == test1_samples$Patient[i] &
      meta$StromaState == test1_samples$StromaState[i]
  ]
  test1_counts_list[[test1_samples$sample_id[i]]] <- Matrix::rowSums(
    counts_mat[, current_cells, drop = FALSE]
  )
}

test1_matrix <- do.call(cbind, test1_counts_list)
storage.mode(test1_matrix) <- "integer"

test1_coldata <- test1_samples %>%
  mutate(condition = StromaState) %>%
  select(sample_id, Patient, Dataset, n_cells, condition) %>%
  column_to_rownames("sample_id")
test1_coldata$Patient <- factor(test1_coldata$Patient)
test1_coldata$condition <- factor(test1_coldata$condition, levels = c("Normal", "Activated"))

test1_out <- run_deseq2(
  count_matrix = test1_matrix,
  coldata_df = test1_coldata,
  coef_name = "condition_Activated_vs_Normal",
  method_label = "Test1_Activated_vs_Normal_DESeq2",
  comparison_label = "Activated_vs_Normal",
  subtype_label = "All_retained_subtypes",
  design_formula = ~ Patient + n_cells + condition
)

# Log2FC interpretation:
# - Positive values mean higher expression in Activated because the coefficient
#   is defined as Activated versus Normal.
# - Negative values mean higher expression in Normal.
test1_full <- test1_out$results %>%
  mutate(
    state_direction = case_when(
      log2FoldChange >= logfc_cutoff & bh_significant ~ "Activated_up",
      log2FoldChange <= -logfc_cutoff & bh_significant ~ "Normal_up",
      TRUE ~ "NS"
    )
  )

test1_kcn <- test1_full %>%
  filter(gene %in% kcn_genes) %>%
  select(gene, log2FoldChange, padj_BH, CellType2, method, comparison, everything())

save_tsv(test1_samples, file.path(table_dir, "Test1_pseudobulk_samples.tsv"))
save_tsv(test1_full, file.path(table_dir, "DEG_Test1_Activated_vs_Normal_full.tsv"))
save_tsv(test1_kcn, file.path(table_dir, "DEG_Test1_Activated_vs_Normal_KCN.tsv"))

###############################################################################
# 7. Test 2 - DESeq2 pseudobulk subtype one-vs-rest
###############################################################################
# Pseudobulk rationale:
# - The one-vs-rest comparison is still performed at the patient level rather
#   than at the cell level, which preserves the replicate structure.
# - Adding n_cells improves robustness by adjusting for different numbers of
#   cells in the target subtype and in the rest compartment.
# - As above, this is a pragmatic pseudobulk covariate used to reduce bias from
#   unbalanced aggregation rather than a single fixed canonical requirement.
#
# KCN rationale:
# - Ion channels may be subtype-restricted and only weakly shifted in global
#   analyses, so subtype-aware pseudobulk contrasts improve sensitivity.
#
# Multiple-testing rule:
# - BH-adjusted p-values are used with a threshold of 0.05.
# - BH is preferred over Bonferroni because transcriptome-wide testing is
#   high-dimensional and Bonferroni would be unnecessarily conservative.
#
# Input filtering:
# - target pseudo-samples must contain at least 20 cells for qPSC, smPSC,
#   myCAF, and csCAF, and at least 10 cells for iCAF
# - the paired rest compartment must contain at least 20 cells
# - only patients with both target and rest pseudo-samples are retained
#
# Practical summary:
# - unit of analysis: paired patient-level pseudo-samples
# - tested coefficient: target subtype versus pooled rest compartment
# - exported moffitt_stromal_exploration_outputs: full DEG table + KCN-only table + volcano + dotplot

test2_results_list <- list()
test2_sample_tables <- list()

for (current_subtype in subtype_levels) {
  contrast_table <- meta %>%
    mutate(group_label = ifelse(CellType2 == current_subtype, current_subtype, "Rest")) %>%
    dplyr::count(Patient, Dataset, group_label, name = "n_cells") %>%
    group_by(Patient) %>%
    filter(
      sum(group_label == current_subtype & n_cells >= min_cells_target[[current_subtype]]) >= 1,
      sum(group_label == "Rest" & n_cells >= min_cells_rest) >= 1
    ) %>%
    ungroup() %>%
    mutate(sample_id = paste(current_subtype, Patient, group_label, sep = "___"))

  group_support <- table(contrast_table$group_label)
  if (!all(c(current_subtype, "Rest") %in% names(group_support))) {
    next
  }
  if (any(group_support[c(current_subtype, "Rest")] < min_patients_required)) {
    next
  }

  pseudo_list <- list()
  for (i in seq_len(nrow(contrast_table))) {
    cell_selector <- if (contrast_table$group_label[i] == current_subtype) {
      meta$CellType2 == current_subtype
    } else {
      meta$CellType2 != current_subtype
    }

    current_cells <- meta$cell_id[
      meta$Patient == contrast_table$Patient[i] &
        cell_selector
    ]

    pseudo_list[[contrast_table$sample_id[i]]] <- Matrix::rowSums(
      counts_mat[, current_cells, drop = FALSE]
    )
  }

  pseudo_matrix <- do.call(cbind, pseudo_list)
  storage.mode(pseudo_matrix) <- "integer"

  coldata <- contrast_table %>%
    mutate(condition = group_label) %>%
    select(sample_id, Patient, Dataset, n_cells, condition) %>%
    column_to_rownames("sample_id")
  coldata$Patient <- factor(coldata$Patient)
  coldata$condition <- factor(coldata$condition, levels = c("Rest", current_subtype))

  de_out <- run_deseq2(
    count_matrix = pseudo_matrix,
    coldata_df = coldata,
    coef_name = paste0("condition_", current_subtype, "_vs_Rest"),
    method_label = "Test2_Subtype_OneVsRest_DESeq2",
    comparison_label = paste0(current_subtype, "_vs_rest"),
    subtype_label = current_subtype,
    design_formula = ~ Patient + n_cells + condition
  )

  test2_results_list[[current_subtype]] <- de_out$results
  test2_sample_tables[[current_subtype]] <- contrast_table
}

# Log2FC interpretation:
# - Positive values mean higher expression in the target subtype than in the
#   pooled rest compartment because the coefficient is target versus Rest.
# - Negative values mean the gene is relatively higher outside the target
#   subtype.
test2_full <- bind_rows(test2_results_list) %>%
  mutate(
    subtype_direction = case_when(
      log2FoldChange >= logfc_cutoff & bh_significant ~ "higher_in_target",
      log2FoldChange <= -logfc_cutoff & bh_significant ~ "lower_in_target",
      TRUE ~ "NS"
    )
  )

test2_kcn <- test2_full %>%
  filter(gene %in% kcn_genes) %>%
  select(gene, log2FoldChange, padj_BH, CellType2, method, comparison, everything())

save_tsv(bind_rows(test2_sample_tables), file.path(table_dir, "Test2_pseudobulk_samples.tsv"))
save_tsv(test2_full, file.path(table_dir, "DEG_Test2_Subtype_OneVsRest_full.tsv"))
save_tsv(test2_kcn, file.path(table_dir, "DEG_Test2_Subtype_OneVsRest_KCN.tsv"))

###############################################################################
# 8. Test 3 - mixed hurdle model subtype one-vs-rest
###############################################################################
# Mixed hurdle rationale:
# - The detection and positive-expression components are modeled separately.
# - This is useful when sparse genes combine dropout-driven detection changes
#   with more subtle intensity shifts among positive cells.
# - The patient random intercept preserves the nested structure of cells within
#   patients, which is critical in single-cell studies.
# - In practice, the detection component is a binomial GLMM fitted with glmer,
#   and the positive-expression component is a linear mixed model fitted with
#   lmer, both with a patient random intercept.
#
# KCN rationale:
# - KCN genes can differ through both detection frequency and expression
#   intensity, so a two-part mixed model is biologically relevant here.
#
# Multiple-testing rule:
# - BH-adjusted p-values are applied within each subtype for both model parts.
# - BH remains more appropriate than Bonferroni because the aim is to control
#   false discovery while retaining sensitivity for moderate ion-channel effects.
#
# Input filtering:
# - the mixed hurdle model is intentionally restricted to the filtered KCN panel
#   plus a compact set of non-KCN background genes selected after prevalence and
#   subtype-variability filtering
# - this broadened panel keeps the hurdle step computationally tractable while
#   still allowing a volcano-style visualization with both KCN and non-KCN
#   genes
# - genes are retained when they are detected in at least 1% of cells overall
#   and in at least 2 stromal subtypes at >=1% detection
#
# Practical summary:
# - detection model: binomial GLMM with patient random intercept
# - intensity model: linear mixed model on positive cells with patient random
#   intercept
# - exported moffitt_stromal_exploration_outputs: full hurdle table + KCN-only table + KCN dotplot +
#   detection barplot + positive-expression boxplot

# Expanded hurdle panel:
# - all filtered KCN are kept
# - non-KCN background genes are first filtered for prevalence and then ranked
#   by subtype-level expression variance
# - the top 300 non-KCN genes are added as visual context for the hurdle
#   volcano panels
idx_list_all <- split(seq_len(ncol(counts_mat)), as.character(meta$CellType2))

pct_by_ct_all <- sapply(idx_list_all, function(idx) {
  Matrix::rowSums(counts_mat[, idx, drop = FALSE] > 0) / length(idx) * 100
})

if (is.vector(pct_by_ct_all)) {
  pct_by_ct_all <- matrix(pct_by_ct_all, ncol = 1)
  colnames(pct_by_ct_all) <- names(idx_list_all)
}

gene_pct_overall <- Matrix::rowSums(counts_mat > 0) / ncol(counts_mat) * 100
gene_n_subtypes_ge5 <- rowSums(pct_by_ct_all[, colnames(pct_by_ct_all) %in% subtype_levels, drop = FALSE] >= 5)

mean_expr_by_ct <- sapply(idx_list_all, function(idx) {
  Matrix::rowMeans(expr_mat[, idx, drop = FALSE])
})

if (is.vector(mean_expr_by_ct)) {
  mean_expr_by_ct <- matrix(mean_expr_by_ct, ncol = 1)
  colnames(mean_expr_by_ct) <- names(idx_list_all)
}

hurdle_gene_panel_tbl <- tibble(
  gene = rownames(stroma),
  pct_cells_overall = as.numeric(gene_pct_overall),
  n_subtypes_ge5pct = as.numeric(gene_n_subtypes_ge5),
  subtype_mean_variance = apply(mean_expr_by_ct, 1, var),
  is_kcn = gene %in% kcn_genes
) %>%
  mutate(
    keep_background = pct_cells_overall >= 5 & n_subtypes_ge5pct >= 2
  )

background_non_kcn <- hurdle_gene_panel_tbl %>%
  filter(keep_background, !is_kcn) %>%
  arrange(desc(subtype_mean_variance), desc(pct_cells_overall)) %>%
  slice_head(n = 300)

hurdle_gene_panel <- unique(c(filtered_kcn, background_non_kcn$gene))

hurdle_gene_panel_tbl <- hurdle_gene_panel_tbl %>%
  mutate(selected_for_hurdle = gene %in% hurdle_gene_panel)

save_tsv(hurdle_gene_panel_tbl, file.path(table_dir, "Test3_hurdle_gene_panel.tsv"))

meta_hurdle <- meta %>%
  mutate(
    Dataset = factor(Dataset),
    Patient = factor(Patient),
    log10_nCount_RNA = log10(nCount_RNA + 1)
  )

contrast_design_tbl <- bind_rows(lapply(subtype_levels, function(ct) {
  meta_hurdle %>%
    mutate(group_simple = ifelse(CellType2 == ct, "target", "rest")) %>%
    dplyr::count(Patient, group_simple, name = "n_cells") %>%
    pivot_wider(names_from = group_simple, values_from = n_cells, values_fill = 0) %>%
    mutate(
      CellType2 = ct,
      min_target_required = min_cells_target[[ct]],
      keep_for_contrast = target >= min_cells_target[[ct]] & rest >= min_cells_rest
    )
})) %>%
  relocate(CellType2, Patient)

hurdle_results_list <- list()

for (ct in subtype_levels) {
  keep_patients <- contrast_design_tbl %>%
    filter(CellType2 == ct, keep_for_contrast) %>%
    pull(Patient) %>%
    as.character()

  if (length(keep_patients) < min_patients_required) {
    next
  }

  meta_ct <- meta_hurdle %>%
    filter(as.character(Patient) %in% keep_patients) %>%
    mutate(is_target = CellType2 == ct)

  for (gene in hurdle_gene_panel) {
    detected <- as.numeric(counts_mat[gene, meta_ct$cell_id] > 0)
    expr_vals <- as.numeric(expr_mat[gene, meta_ct$cell_id])

    df_gene <- meta_ct %>%
      transmute(
        cell_id,
        Patient,
        Dataset,
        is_target,
        log10_nCount_RNA,
        detected = detected,
        expr = expr_vals
      )

    mean_target_count <- mean(counts_mat[gene, meta_ct$cell_id[meta_ct$is_target]])
    mean_rest_count <- mean(counts_mat[gene, meta_ct$cell_id[!meta_ct$is_target]])
    empirical_log2fc <- log2((mean_target_count + 1) / (mean_rest_count + 1))

    det_ok <- length(unique(df_gene$Patient[df_gene$is_target])) >= min_positive_patients &&
      length(unique(df_gene$Patient[!df_gene$is_target])) >= min_positive_patients &&
      sum(df_gene$detected) >= min_positive_cells_total &&
      length(unique(df_gene$detected)) > 1

    det_res <- list(beta = NA_real_, se = NA_real_, z = NA_real_, pval_lrt = NA_real_, odds_ratio = NA_real_, warnings = "not_fitted")
    det_fit <- FALSE
    if (det_ok) {
      det_tmp <- tryCatch(fit_detection_model(df_gene), error = function(e) list(error = conditionMessage(e)))
      if (!is.null(det_tmp$error)) {
        det_res$warnings <- det_tmp$error
      } else {
        det_res <- det_tmp
        det_fit <- TRUE
      }
    }

    df_gene_pos <- df_gene %>%
      filter(detected == 1) %>%
      transmute(Patient, Dataset, is_target, log10_nCount_RNA, expr_pos = expr)

    pos_ok <- nrow(df_gene_pos) >= min_positive_cells_total &&
      sum(df_gene_pos$is_target) >= 10 &&
      sum(!df_gene_pos$is_target) >= 10 &&
      length(unique(df_gene_pos$Patient[df_gene_pos$is_target])) >= min_positive_patients &&
      length(unique(df_gene_pos$Patient[!df_gene_pos$is_target])) >= min_positive_patients &&
      n_distinct(df_gene_pos$expr_pos) > 1

    pos_res <- list(beta = NA_real_, se = NA_real_, t = NA_real_, pval_lrt = NA_real_, warnings = "not_fitted")
    pos_fit <- FALSE
    if (pos_ok) {
      pos_tmp <- tryCatch(fit_positive_model(df_gene_pos), error = function(e) list(error = conditionMessage(e)))
      if (!is.null(pos_tmp$error)) {
        pos_res$warnings <- pos_tmp$error
      } else {
        pos_res <- pos_tmp
        pos_fit <- TRUE
      }
    }

    hurdle_results_list[[paste(ct, gene, sep = "|")]] <- data.frame(
      CellType2 = ct,
      gene = gene,
      n_cells = nrow(df_gene),
      n_detected_total = sum(df_gene$detected),
      n_detected_target = sum(df_gene$detected[df_gene$is_target]),
      n_detected_rest = sum(df_gene$detected[!df_gene$is_target]),
      n_patients = dplyr::n_distinct(df_gene$Patient),
      n_patients_target = dplyr::n_distinct(df_gene$Patient[df_gene$is_target]),
      n_patients_rest = dplyr::n_distinct(df_gene$Patient[!df_gene$is_target]),
      empirical_log2FC = empirical_log2fc,
      detection_model_fit = det_fit,
      detection_beta = det_res$beta,
      detection_OR = det_res$odds_ratio,
      detection_p = det_res$pval_lrt,
      detection_warning = paste(det_res$warnings, collapse = " | "),
      positive_model_fit = pos_fit,
      positive_beta = pos_res$beta,
      positive_p = pos_res$pval_lrt,
      positive_warning = paste(pos_res$warnings, collapse = " | "),
      stringsAsFactors = FALSE
    )
  }
}

test3_full <- bind_rows(hurdle_results_list) %>%
  group_by(CellType2) %>%
  mutate(
    detection_padj = p.adjust(detection_p, method = "BH"),
    positive_padj = p.adjust(positive_p, method = "BH"),
    combined_p = case_when(
      is.finite(detection_p) & is.finite(positive_p) ~ pchisq(
        -2 * (log(detection_p) + log(positive_p)),
        df = 4,
        lower.tail = FALSE
      ),
      is.finite(detection_p) ~ detection_p,
      is.finite(positive_p) ~ positive_p,
      TRUE ~ NA_real_
    ),
    combined_padj = p.adjust(combined_p, method = "BH")
  ) %>%
  ungroup() %>%
  left_join(hurdle_gene_panel_tbl, by = "gene") %>%
  mutate(
    method = "Test3_Subtype_Hurdle_Mixed",
    comparison = paste0(CellType2, "_vs_rest"),
    padj_BH = combined_padj,
    detection_sig = !is.na(detection_padj) & detection_padj < padj_cutoff & detection_beta > 0,
    positive_sig = !is.na(positive_padj) & positive_padj < padj_cutoff & positive_beta > 0,
    support_pattern = case_when(
      detection_sig & positive_sig ~ "both_parts",
      detection_sig & !positive_sig ~ "detection_only",
      !detection_sig & positive_sig ~ "positive_only",
      TRUE ~ "none"
    ),
    combined_score = -log10(pmax(combined_padj, 1e-300)),
    effect_beta = case_when(
      !is.na(positive_beta) ~ positive_beta,
      !is.na(detection_beta) ~ detection_beta,
      TRUE ~ 0
    ),
    bh_significant = support_pattern != "none"
  ) %>%
  arrange(CellType2, desc(combined_score), gene)

test3_kcn <- test3_full %>%
  select(gene, effect_beta, detection_beta, positive_beta, padj_BH, CellType2, method, comparison, everything())

save_tsv(contrast_design_tbl, file.path(table_dir, "Test3_hurdle_eligibility.tsv"))
save_tsv(test3_full, file.path(table_dir, "DEG_Test3_Subtype_Hurdle_Mixed_full.tsv"))
save_tsv(test3_kcn, file.path(table_dir, "DEG_Test3_Subtype_Hurdle_Mixed_KCN.tsv"))

###############################################################################
# 9. Union of KCN genes across the 3 methods
###############################################################################
# Union rationale:
# - Different DEG methods capture complementary signals.
# - The union therefore increases sensitivity and robustness.
# - This is particularly important for ion channels, which may show subtle
#   but reproducible effects across complementary statistical frameworks.

test1_kcn_sig <- test1_kcn %>%
  filter(bh_significant)

test2_kcn_sig <- test2_kcn %>%
  filter(bh_significant, abs(log2FoldChange) >= logfc_cutoff)

test3_kcn_sig <- test3_kcn %>%
  filter(support_pattern != "none")

kcn_method_sets <- list(
  Test1_Activated_vs_Normal_DESeq2 = sort(unique(test1_kcn_sig$gene)),
  Test2_Subtype_OneVsRest_DESeq2 = sort(unique(test2_kcn_sig$gene)),
  Test3_Subtype_Hurdle_Mixed = sort(unique(test3_kcn_sig$gene))
)

union_kcn <- sort(unique(unlist(kcn_method_sets)))

kcn_union_table <- tibble(gene = union_kcn) %>%
  mutate(
    Test1_Activated_vs_Normal_DESeq2 = gene %in% kcn_method_sets$Test1_Activated_vs_Normal_DESeq2,
    Test2_Subtype_OneVsRest_DESeq2 = gene %in% kcn_method_sets$Test2_Subtype_OneVsRest_DESeq2,
    Test3_Subtype_Hurdle_Mixed = gene %in% kcn_method_sets$Test3_Subtype_Hurdle_Mixed
  ) %>%
  rowwise() %>%
  mutate(
    n_methods = sum(c_across(Test1_Activated_vs_Normal_DESeq2:Test3_Subtype_Hurdle_Mixed)),
    methods_present = paste(
      c("Test1_Activated_vs_Normal_DESeq2", "Test2_Subtype_OneVsRest_DESeq2", "Test3_Subtype_Hurdle_Mixed")[
        c_across(Test1_Activated_vs_Normal_DESeq2:Test3_Subtype_Hurdle_Mixed)
      ],
      collapse = "; "
    )
  ) %>%
  ungroup()

save_tsv(kcn_union_table, file.path(table_dir, "KCN_union_membership.tsv"))

###############################################################################
# 10. Venn diagram of KCN overlap
###############################################################################

venn_grob <- venn.diagram(
  x = kcn_method_sets,
  filename = NULL,
  fill = unname(venn_colors),
  alpha = 0.75,
  col = "grey25",
  lwd = 1.1,
  cex = 1.4,
  fontface = "bold",
  cat.cex = 0,
  category.names = rep("", 3),
  margin = 0.1
)

venn_legend_df <- tibble(
  method = names(venn_colors),
  color = unname(venn_colors),
  genes = vapply(kcn_method_sets, function(x) {
    if (length(x) == 0) {
      return("None")
    }
    wrap_text(paste(x, collapse = ", "), width = 42)
  }, character(1))
) %>%
  mutate(
    y = rev(seq_len(n())),
    label = paste0(method, "\n", genes)
  )

legend_plot <- ggplot(venn_legend_df, aes(x = 0, y = y)) +
  geom_point(aes(color = method), size = 6) +
  geom_text(aes(label = label), hjust = 0, nudge_x = 0.18, size = 3.2, lineheight = 0.95) +
  scale_color_manual(values = venn_colors) +
  coord_cartesian(xlim = c(0, 5.8), ylim = c(0.5, nrow(venn_legend_df) + 0.5), clip = "off") +
  theme_void() +
  theme(legend.position = "none")

p_venn <- wrap_elements(full = venn_grob) + legend_plot + plot_layout(widths = c(1.15, 1))

save_plot(
  p_venn,
  file.path(figure_dir, "Figure_02_Venn_KCN_overlap_3_methods.png"),
  width = 16,
  height = 8
)

###############################################################################
# 11. DotPlots of KCN genes
###############################################################################

test1_dot_genes <- test1_kcn %>%
  filter(bh_significant) %>%
  arrange(padj_BH, desc(abs(log2FoldChange))) %>%
  slice_head(n = 15) %>%
  pull(gene) %>%
  unique()

if (length(test1_dot_genes) == 0) {
  test1_dot_genes <- head(test1_kcn$gene, 15)
}

test2_dot_genes <- test2_kcn %>%
  filter(bh_significant, log2FoldChange > 0) %>%
  group_by(CellType2) %>%
  arrange(padj_BH, desc(log2FoldChange), .by_group = TRUE) %>%
  slice_head(n = 4) %>%
  ungroup() %>%
  pull(gene) %>%
  unique()

if (length(test2_dot_genes) == 0) {
  test2_dot_genes <- head(filtered_kcn, 15)
}

test3_dot_genes <- test3_kcn %>%
  filter(support_pattern != "none") %>%
  group_by(CellType2) %>%
  arrange(desc(combined_score), .by_group = TRUE) %>%
  slice_head(n = 4) %>%
  ungroup() %>%
  pull(gene) %>%
  unique()

if (length(test3_dot_genes) == 0) {
  test3_dot_genes <- head(filtered_kcn, 15)
}

make_dotplot(
  seu = stroma,
  features = test1_dot_genes,
  group_by = "StromaState",
  colors_use = state_colors,
  title_text = "KCN DotPlot - Activated vs Normal",
  file_path = file.path(figure_dir, "Figure_03_DotPlot_KCN_Test1_Activated_vs_Normal.png"),
  width = 12,
  height = 6
)

make_dotplot(
  seu = stroma,
  features = test2_dot_genes,
  group_by = "CellType2",
  colors_use = celltype_colors,
  title_text = "KCN DotPlot - DESeq2 subtype one-vs-rest",
  file_path = file.path(figure_dir, "Figure_04_DotPlot_KCN_Test2_Subtype_OneVsRest.png"),
  width = 14,
  height = 6
)

make_dotplot(
  seu = stroma,
  features = test3_dot_genes,
  group_by = "CellType2",
  colors_use = celltype_colors,
  title_text = "KCN DotPlot - mixed hurdle subtype one-vs-rest",
  file_path = file.path(figure_dir, "Figure_05_DotPlot_KCN_Test3_Hurdle_Mixed.png"),
  width = 14,
  height = 6
)

###############################################################################
# 12. Volcano plots
###############################################################################

test1_volcano <- test1_full %>%
  filter(!is.na(log2FoldChange), !is.na(padj_BH)) %>%
  mutate(
    neg_log10_padj = -log10(padj_BH + 1e-300),
    strong_signal = bh_significant & abs(log2FoldChange) >= logfc_cutoff,
    plot_group = case_when(
      bh_significant & gene %in% kcn_genes ~ "KCN",
      strong_signal ~ "Significant_non_KCN",
      TRUE ~ "NS"
    ),
    plot_color = case_when(
      strong_signal & log2FoldChange >= logfc_cutoff ~ "Activated",
      strong_signal & log2FoldChange <= -logfc_cutoff ~ "Normal",
      TRUE ~ NA_character_
    )
  )

test1_labels <- test1_volcano %>%
  filter(plot_group == "KCN")

make_volcano_plot(
  df = test1_volcano,
  x_col = "log2FoldChange",
  y_col = "neg_log10_padj",
  p_threshold = padj_cutoff,
  lfc_threshold = logfc_cutoff,
  color_map = state_colors,
  label_df = test1_labels,
  title_text = "Volcano plot - Activated vs Normal",
  x_label = "log2FC",
  y_label = "-log10(padj)",
  file_path = file.path(figure_dir, "Figure_06_Volcano_Test1_Activated_vs_Normal.png"),
  width = 10,
  height = 8
)

test2_volcano <- test2_full %>%
  filter(!is.na(log2FoldChange), !is.na(padj_BH)) %>%
  mutate(
    neg_log10_padj = -log10(padj_BH + 1e-300),
    strong_signal = bh_significant & abs(log2FoldChange) >= logfc_cutoff,
    plot_group = case_when(
      bh_significant & gene %in% kcn_genes ~ "KCN",
      strong_signal ~ "Significant_non_KCN",
      TRUE ~ "NS"
    ),
    plot_color = as.character(CellType2)
  )

test2_labels <- test2_volcano %>%
  filter(plot_group == "KCN")

make_volcano_plot(
  df = test2_volcano,
  x_col = "log2FoldChange",
  y_col = "neg_log10_padj",
  p_threshold = padj_cutoff,
  lfc_threshold = logfc_cutoff,
  facet_col = "CellType2",
  color_map = celltype_colors,
  label_df = test2_labels,
  title_text = "Volcano plot - DESeq2 subtype one-vs-rest",
  x_label = "log2FC",
  y_label = "-log10(padj)",
  file_path = file.path(figure_dir, "Figure_07_Volcano_Test2_Subtype_OneVsRest.png"),
  width = 14,
  height = 10
)

hurdle_top_kcn <- test3_kcn %>%
  filter(support_pattern != "none") %>%
  group_by(CellType2) %>%
  arrange(desc(combined_score), .by_group = TRUE) %>%
  slice_head(n = 2) %>%
  ungroup()

if (nrow(hurdle_top_kcn) == 0) {
  hurdle_top_kcn <- test3_kcn %>%
    group_by(CellType2) %>%
    arrange(desc(combined_score), .by_group = TRUE) %>%
    slice_head(n = 2) %>%
    ungroup()
}

save_tsv(hurdle_top_kcn, file.path(table_dir, "Test3_hurdle_top_KCN_for_plots.tsv"))

###############################################################################
# 13. Save objects and logs
###############################################################################

saveRDS(
  list(
    Test1 = test1_full,
    Test2 = test2_full,
    Test3 = test3_full,
    KCN_union = kcn_union_table
  ),
  file = file.path(rds_dir, "deg3_results_bundle.rds")
)

writeLines(
  c(
    "Moffitt stromal KCN DEG pipeline (3 methods)",
    paste("Run date:", as.character(Sys.Date())),
    paste("Input object:", input_rds),
    paste("Current union KCN genes:", if (length(union_kcn) > 0) paste(union_kcn, collapse = ", ") else "None")
  ),
  con = file.path(output_root, "README_pipeline_summary.txt")
)

message("KCN DEG pipeline finished.")
message("Outputs saved in: ", output_root)

