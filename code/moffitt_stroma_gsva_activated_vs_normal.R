###############################################################################
# Moffitt stroma GSVA pipeline: Activated vs Normal
#
# This script compares Activated and Normal stromal states with pathway-level
# scoring and then looks at how KCN-high vs KCN-low cells behave inside each
# fibroblast subtype.
###############################################################################
required_packages <- c(
  "Seurat",
  "dplyr",
  "ggplot2",
  "tibble",
  "tidyr",
  "purrr",
  "Matrix",
  "GSVA",
  "msigdbr",
  "AUCell",
  "BiocParallel"
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
  library(tibble)
  library(tidyr)
  library(purrr)
  library(Matrix)
  library(GSVA)
  library(msigdbr)
  library(AUCell)
  library(BiocParallel)
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

clean_classic_theme <- function() {
  theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.title = element_text(face = "bold"),
      panel.grid = element_blank()
    )
}

format_hallmark_name <- function(x) {
  x <- gsub("^HALLMARK_", "", x)
  gsub("_", " ", x)
}

safe_file_label <- function(x) {
  gsub("[^A-Za-z0-9]+", "_", x)
}

get_selected_threshold <- function(cells_auc_object, pathway_name) {
  threshold_obj <- AUCell_exploreThresholds(
    cells_auc_object[pathway_name, ],
    plotHist = FALSE,
    assign = TRUE,
    nCores = 1,
    verbose = FALSE
  )

  auc_values <- as.numeric(getAUC(cells_auc_object[pathway_name, ]))
  max_auc <- max(auc_values, na.rm = TRUE)

  selected_name <- names(threshold_obj[[pathway_name]]$aucThr$selected)[1]
  selected_value <- as.numeric(threshold_obj[[pathway_name]]$aucThr$selected[1])

  threshold_tbl <- as.data.frame(threshold_obj[[pathway_name]]$aucThr$thresholds)
  threshold_tbl$threshold_name <- rownames(threshold_tbl)

  valid_thresholds <- threshold_tbl %>%
    filter(nCells > 0, threshold <= max_auc) %>%
    arrange(desc(threshold))

  if (is.finite(selected_value) && selected_value <= max_auc && sum(auc_values > selected_value) > 0) {
    return(list(
      threshold = selected_value,
      threshold_name = selected_name,
      used_fallback = FALSE
    ))
  }

  if (nrow(valid_thresholds) == 0) {
    fallback_value <- stats::quantile(auc_values, probs = 0.95, na.rm = TRUE)
    return(list(
      threshold = as.numeric(fallback_value),
      threshold_name = "empirical_q95",
      used_fallback = TRUE
    ))
  }

  list(
    threshold = valid_thresholds$threshold[1],
    threshold_name = valid_thresholds$threshold_name[1],
    used_fallback = TRUE
  )
}

###############################################################################
# 2. Paths, colors, and constants
###############################################################################

project_root <- find_project_root()
analysis_root <- file.path(project_root, "M1_Bioinformatics_Ion_Channel", "code")
project_data_dir <- file.path(project_root, "M1_Bioinformatics_Ion_Channel", "data")
input_rds <- file.path(project_data_dir, "Stroma_Subset2021_fibroblast_focus.rds")

output_root <- file.path(analysis_root, "moffitt_activated_vs_normal_gsva_outputs")
figure_dir <- file.path(output_root, "figures")
table_dir <- file.path(output_root, "tables")
aucell_root <- file.path(output_root, "aucell_by_subtype")
aucell_figure_dir <- file.path(aucell_root, "figures")
aucell_table_dir <- file.path(aucell_root, "tables")
deg_table_dir <- file.path(analysis_root, "moffitt_kcn_deg_4_methods", "tables")

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(aucell_root, recursive = TRUE, showWarnings = FALSE)
dir.create(aucell_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(aucell_table_dir, recursive = TRUE, showWarnings = FALSE)

state_colors <- c(
  "Activated" = "#FF686B",
  "Normal" = "#70D6FF"
)

subtype_colors <- c(
  "qPSC" = "#70D6FF",
  "smPSC" = "#0096C7",
  "myCAF" = "#FF686B",
  "csCAF" = "#FF9770",
  "iCAF" = "#FF5400"
)
subtype_order <- names(subtype_colors)

min_cells_per_state <- 20
padj_cutoff <- 0.05
aucell_padj_cutoff <- 0.05
aucell_min_cells_per_group <- 20
aucell_min_cells_per_patient_group <- 5
aucell_min_paired_patients <- 3

###############################################################################
# 3. Load the Seurat object and build the stromal state
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
  "Moffitt.F5_ActivatedStroma.top25",
  "Moffitt.F13_NormalStroma.top25"
)

missing_meta <- setdiff(required_meta, colnames(stroma@meta.data))
if (length(missing_meta) > 0) {
  stop("Missing metadata columns: ", paste(missing_meta, collapse = ", "))
}

stroma$StromaState <- ifelse(
  stroma$Moffitt.F5_ActivatedStroma.top25 > stroma$Moffitt.F13_NormalStroma.top25,
  "Activated",
  "Normal"
)

stroma$StromaState <- factor(stroma$StromaState, levels = c("Normal", "Activated"))

meta <- stroma[[]] %>%
  rownames_to_column("cell_id")

counts_mat <- GetAssayData(stroma, assay = "RNA", layer = "counts")

###############################################################################
# 4. Build paired patient-state pseudobulks
###############################################################################
# Pseudobulk aggregation is used here because pathway scores should be compared
# at the level of biological replicates rather than treating every cell as an
# independent sample. The paired design preserves within-patient contrast
# between Normal and Activated stromal states.

pseudobulk_samples <- meta %>%
  count(Patient, Dataset, StromaState, name = "n_cells") %>%
  filter(n_cells >= min_cells_per_state) %>%
  group_by(Patient) %>%
  filter(all(c("Normal", "Activated") %in% StromaState)) %>%
  ungroup() %>%
  mutate(sample_id = paste(Patient, StromaState, sep = "___"))

pseudobulk_counts_list <- list()
for (i in seq_len(nrow(pseudobulk_samples))) {
  current_cells <- meta$cell_id[
    meta$Patient == pseudobulk_samples$Patient[i] &
      meta$StromaState == pseudobulk_samples$StromaState[i]
  ]
  pseudobulk_counts_list[[pseudobulk_samples$sample_id[i]]] <- Matrix::rowSums(
    counts_mat[, current_cells, drop = FALSE]
  )
}

pseudobulk_counts <- do.call(cbind, pseudobulk_counts_list)
storage.mode(pseudobulk_counts) <- "numeric"

library_sizes <- colSums(pseudobulk_counts)
log_cpm <- log1p(t(t(pseudobulk_counts) / library_sizes) * 1e6)

save_tsv(pseudobulk_samples, file.path(table_dir, "GSVA_pseudobulk_samples.tsv"))

###############################################################################
# 5. Build Hallmark gene sets
###############################################################################

hallmark_tbl <- msigdbr(
  species = "Homo sapiens",
  collection = "H"
) %>%
  select(gs_name, gene_symbol)

hallmark_gene_sets <- split(hallmark_tbl$gene_symbol, hallmark_tbl$gs_name) %>%
  lapply(unique)

hallmark_gene_sets <- hallmark_gene_sets[
  vapply(hallmark_gene_sets, function(gs) {
    overlap_n <- sum(gs %in% rownames(log_cpm))
    overlap_n >= 10 && overlap_n <= 500
  }, logical(1))
]

hallmark_sizes <- tibble(
  pathway = names(hallmark_gene_sets),
  n_genes_in_matrix = vapply(hallmark_gene_sets, function(gs) sum(gs %in% rownames(log_cpm)), integer(1))
)

save_tsv(hallmark_sizes, file.path(table_dir, "GSVA_hallmark_gene_set_sizes.tsv"))

###############################################################################
# 6. Run GSVA
###############################################################################

gsva_par <- gsvaParam(
  exprData = as.matrix(log_cpm),
  geneSets = hallmark_gene_sets,
  kcdf = "Gaussian"
)

gsva_scores <- gsva(gsva_par, verbose = FALSE)

gsva_scores_tbl <- as.data.frame(gsva_scores) %>%
  rownames_to_column("pathway") %>%
  pivot_longer(
    cols = -pathway,
    names_to = "sample_id",
    values_to = "gsva_score"
  ) %>%
  left_join(
    pseudobulk_samples %>% select(sample_id, Patient, Dataset, StromaState, n_cells),
    by = "sample_id"
  )

save_tsv(gsva_scores_tbl, file.path(table_dir, "GSVA_hallmark_scores_long.tsv"))

###############################################################################
# 7. Compare Activated and Normal GSVA scores
###############################################################################

gsva_stats <- gsva_scores_tbl %>%
  group_by(pathway) %>%
  group_modify(~{
    wide_df <- .x %>%
      select(Patient, StromaState, gsva_score) %>%
      distinct() %>%
      pivot_wider(names_from = StromaState, values_from = gsva_score) %>%
      filter(!is.na(Normal), !is.na(Activated))

    p_value <- wilcox.test(
      x = wide_df$Activated,
      y = wide_df$Normal,
      paired = TRUE,
      exact = FALSE
    )$p.value

    tibble(
      n_patients = nrow(wide_df),
      mean_normal = mean(wide_df$Normal),
      mean_activated = mean(wide_df$Activated),
      delta_activated_minus_normal = mean(wide_df$Activated - wide_df$Normal),
      median_delta = median(wide_df$Activated - wide_df$Normal),
      p_value = p_value
    )
  }) %>%
  ungroup() %>%
  mutate(
    padj_BH = p.adjust(p_value, method = "BH"),
    direction = ifelse(delta_activated_minus_normal >= 0, "Activated", "Normal"),
    pathway_label = format_hallmark_name(pathway),
    significant = padj_BH < padj_cutoff
  ) %>%
  arrange(desc(delta_activated_minus_normal))

gsva_stats_sig <- gsva_stats %>%
  filter(significant) %>%
  arrange(delta_activated_minus_normal)

save_tsv(gsva_stats, file.path(table_dir, "GSVA_hallmark_statistics.tsv"))
save_tsv(gsva_stats_sig, file.path(table_dir, "GSVA_hallmark_statistics_significant.tsv"))

###############################################################################
# 8. Barplot of pathway differences
###############################################################################

plot_df <- gsva_stats_sig
if (nrow(plot_df) == 0) {
  plot_df <- gsva_stats %>%
    slice_max(order_by = abs(delta_activated_minus_normal), n = 15) %>%
    arrange(delta_activated_minus_normal)
}

plot_df$pathway_label <- factor(plot_df$pathway_label, levels = plot_df$pathway_label)

p_gsva_barplot <- ggplot(
  plot_df,
  aes(x = pathway_label, y = delta_activated_minus_normal, fill = direction)
) +
  geom_col(width = 0.78) +
  coord_flip() +
  scale_fill_manual(values = state_colors) +
  labs(
    title = "GSVA pathway differences between Activated and Normal stroma",
    x = NULL,
    y = "GSVA score difference (Activated - Normal)",
    fill = NULL
  ) +
  clean_classic_theme() +
  theme(
    axis.text.y = element_text(face = "bold"),
    legend.position = "top"
  )

save_plot(
  p_gsva_barplot,
  file.path(figure_dir, "Figure_01_GSVA_barplot_activated_vs_normal.png"),
  width = 11,
  height = 8
)

###############################################################################
# 9. Prepare subtype-level KCN candidates for AUCell
###############################################################################
# AUCell is added here as a complementary single-cell functional analysis.
# It uses rank-based pathway scoring within each stromal subtype, which is
# useful when ion-channel expression varies across subsets of cells rather than
# through a uniform shift across the entire subtype.
#
# To keep the KCN panel consistent with the main differential-expression
# section of the project, the AUCell step uses the final union of KCN genes
# recovered across the four DEG frameworks. This keeps the downstream
# functional analysis aligned with the full KCN candidate space retained after
# multi-method DEG integration.

deg_test2_file <- file.path(deg_table_dir, "DEG_Test2_Subtype_OneVsRest_KCN.tsv")
deg_test3_file <- file.path(deg_table_dir, "DEG_Test3_Subtype_Hurdle_Mixed_KCN.tsv")
deg_test4_file <- file.path(deg_table_dir, "DEG_Test4_Subtype_MAST_KCN.tsv")
deg_union_file <- file.path(deg_table_dir, "KCN_union_membership.tsv")

deg_files <- c(deg_test2_file, deg_test3_file, deg_test4_file, deg_union_file)
missing_deg_files <- deg_files[!file.exists(deg_files)]
if (length(missing_deg_files) > 0) {
  stop("Missing DEG tables required for AUCell: ", paste(missing_deg_files, collapse = ", "))
}

deg_test2_kcn <- read.delim(deg_test2_file, check.names = FALSE) %>%
  filter(bh_significant, subtype_direction == "higher_in_target") %>%
  transmute(CellType2, gene, method = "Test2_DESeq2_subtype")

deg_test3_kcn <- read.delim(deg_test3_file, check.names = FALSE) %>%
  filter(bh_significant) %>%
  mutate(
    positive_support = (!is.na(detection_beta) & detection_beta > 0) |
      (!is.na(positive_beta) & positive_beta > 0)
  ) %>%
  filter(positive_support) %>%
  transmute(CellType2, gene, method = "Test3_mixed_hurdle")

deg_test4_kcn <- read.delim(deg_test4_file, check.names = FALSE) %>%
  filter(bh_significant, subtype_direction == "higher_in_target") %>%
  transmute(CellType2, gene, method = "Test4_MAST")

aucell_kcn_support_by_subtype <- bind_rows(
  deg_test2_kcn,
  deg_test3_kcn,
  deg_test4_kcn
) %>%
  distinct() %>%
  count(CellType2, gene, name = "n_deg_methods") %>%
  arrange(CellType2, desc(n_deg_methods), gene)

deg_union_kcn <- read.delim(deg_union_file, check.names = FALSE) %>%
  transmute(
    gene,
    n_union_methods = n_methods,
    union_methods_present = methods_present
  ) %>%
  arrange(desc(n_union_methods), gene)

aucell_kcn_panel <- tidyr::crossing(
  CellType2 = subtype_order,
  deg_union_kcn
) %>%
  left_join(aucell_kcn_support_by_subtype, by = c("CellType2", "gene")) %>%
  mutate(n_deg_methods = ifelse(is.na(n_deg_methods), 0L, n_deg_methods)) %>%
  arrange(CellType2, desc(n_union_methods), gene)

save_tsv(
  aucell_kcn_support_by_subtype,
  file.path(aucell_table_dir, "AUCell_KCN_support_by_subtype.tsv")
)

save_tsv(
  deg_union_kcn,
  file.path(aucell_table_dir, "AUCell_KCN_union_panel.tsv")
)

save_tsv(
  aucell_kcn_panel,
  file.path(aucell_table_dir, "AUCell_KCN_candidates_by_subtype.tsv")
)

###############################################################################
# 10. AUCell within each stromal fibroblast subtype
###############################################################################
# AUCell is computed inside each subtype separately so that pathway variation is
# interpreted within a biologically coherent stromal compartment. For each KCN,
# cells are split into KCN-high and KCN-low groups within the subtype.
#
# Statistical testing is then performed on patient-level mean AUCell scores
# using paired Wilcoxon tests. This preserves the biological replicate
# structure and avoids treating all cells as independent observations.
#
# Benjamini-Hochberg correction is used for each KCN across Hallmark pathways.
# This controls the false discovery rate while remaining less conservative than
# Bonferroni in high-dimensional transcriptomic analyses.

subtype_hallmark_sizes <- list()
aucell_split_summary <- list()
aucell_stats_all <- list()
aucell_stats_significant <- list()
aucell_cache <- list()

for (subtype in subtype_order) {
  subtype_kcn_tbl <- aucell_kcn_panel %>%
    filter(CellType2 == subtype)

  if (nrow(subtype_kcn_tbl) == 0) {
    next
  }

  subtype_obj <- subset(stroma, subset = CellType2 == subtype)
  subtype_meta <- subtype_obj[[]] %>%
    rownames_to_column("cell_id") %>%
    mutate(
      Patient = as.character(Patient),
      Dataset = as.character(Dataset),
      Condition = as.character(Condition),
      CellType2 = as.character(CellType2)
    )

  subtype_counts <- as.matrix(GetAssayData(subtype_obj, assay = "RNA", layer = "counts"))
  subtype_norm <- as.matrix(GetAssayData(subtype_obj, assay = "RNA", layer = "data"))

  subtype_hallmark_sets <- lapply(hallmark_gene_sets, intersect, y = rownames(subtype_counts))
  subtype_hallmark_sets <- subtype_hallmark_sets[
    vapply(subtype_hallmark_sets, length, integer(1)) >= 10
  ]

  subtype_hallmark_sizes[[subtype]] <- tibble(
    CellType2 = subtype,
    pathway = names(subtype_hallmark_sets),
    n_genes_in_matrix = vapply(subtype_hallmark_sets, length, integer(1))
  )

  subtype_rankings <- AUCell_buildRankings(
    exprMat = subtype_counts,
    plotStats = FALSE,
    BPPARAM = BiocParallel::SerialParam()
  )

  subtype_auc <- AUCell_calcAUC(
    geneSets = subtype_hallmark_sets,
    rankings = subtype_rankings,
    aucMaxRank = ceiling(nrow(subtype_counts) * 0.05)
  )

  subtype_auc_tbl <- as.data.frame(t(getAUC(subtype_auc))) %>%
    rownames_to_column("cell_id")

  aucell_cache[[subtype]] <- list(
    cells_auc = subtype_auc,
    auc_mat = getAUC(subtype_auc),
    counts = subtype_counts,
    norm = subtype_norm,
    meta = subtype_meta
  )

  for (i in seq_len(nrow(subtype_kcn_tbl))) {
    current_kcn <- subtype_kcn_tbl$gene[i]
    current_support <- subtype_kcn_tbl$n_deg_methods[i]

    if (!(current_kcn %in% rownames(subtype_norm))) {
      aucell_split_summary[[length(aucell_split_summary) + 1]] <- tibble(
        CellType2 = subtype,
        gene = current_kcn,
        n_deg_methods = current_support,
        median_expression = NA_real_,
        n_cells_low = NA_integer_,
        n_cells_high = NA_integer_,
        n_paired_patients = 0L,
        eligible_for_testing = FALSE,
        exclusion_reason = "gene_absent_from_expression_matrix"
      )
      next
    }

    kcn_expression <- as.numeric(subtype_norm[current_kcn, subtype_meta$cell_id])
    names(kcn_expression) <- subtype_meta$cell_id
    median_expression <- median(kcn_expression, na.rm = TRUE)

    expression_group <- if (median_expression == 0) {
      ifelse(kcn_expression > 0, "High", "Low")
    } else {
      ifelse(kcn_expression > median_expression, "High", "Low")
    }

    split_tbl <- subtype_meta %>%
      transmute(
        cell_id,
        Patient,
        Dataset,
        Condition,
        CellType2,
        gene = current_kcn,
        kcn_expression = kcn_expression[cell_id],
        group = factor(expression_group, levels = c("Low", "High"))
      )

    eligible_patients <- split_tbl %>%
      count(Patient, group, name = "n_cells_group") %>%
      pivot_wider(names_from = group, values_from = n_cells_group, values_fill = 0)

    if (!("Low" %in% colnames(eligible_patients))) {
      eligible_patients$Low <- 0L
    }
    if (!("High" %in% colnames(eligible_patients))) {
      eligible_patients$High <- 0L
    }

    eligible_patients <- eligible_patients %>%
      filter(
        Low >= aucell_min_cells_per_patient_group,
        High >= aucell_min_cells_per_patient_group
      )

    n_cells_low <- sum(split_tbl$group == "Low")
    n_cells_high <- sum(split_tbl$group == "High")
    n_paired_patients <- nrow(eligible_patients)

    is_eligible <- n_cells_low >= aucell_min_cells_per_group &&
      n_cells_high >= aucell_min_cells_per_group &&
      n_paired_patients >= aucell_min_paired_patients

    exclusion_reason <- NA_character_
    if (!is_eligible) {
      exclusion_reason <- paste(c(
        if (n_cells_low < aucell_min_cells_per_group) "too_few_low_cells" else NULL,
        if (n_cells_high < aucell_min_cells_per_group) "too_few_high_cells" else NULL,
        if (n_paired_patients < aucell_min_paired_patients) "too_few_paired_patients" else NULL
      ), collapse = ";")
    }

    aucell_split_summary[[length(aucell_split_summary) + 1]] <- tibble(
      CellType2 = subtype,
      gene = current_kcn,
      n_deg_methods = current_support,
      median_expression = median_expression,
      n_cells_low = n_cells_low,
      n_cells_high = n_cells_high,
      n_paired_patients = n_paired_patients,
      eligible_for_testing = is_eligible,
      exclusion_reason = exclusion_reason
    )

    if (!is_eligible) {
      next
    }

    patient_auc_tbl <- split_tbl %>%
      filter(Patient %in% eligible_patients$Patient) %>%
      select(cell_id, Patient, group) %>%
      left_join(subtype_auc_tbl, by = "cell_id") %>%
      pivot_longer(
        cols = -c(cell_id, Patient, group),
        names_to = "pathway",
        values_to = "auc_score"
      ) %>%
      group_by(Patient, group, pathway) %>%
      summarise(
        mean_auc = mean(auc_score),
        n_cells = n(),
        .groups = "drop"
      )

    pathway_stats_input <- patient_auc_tbl %>%
      select(Patient, group, pathway, mean_auc) %>%
      pivot_wider(names_from = group, values_from = mean_auc)

    if (!("Low" %in% colnames(pathway_stats_input))) {
      pathway_stats_input$Low <- NA_real_
    }
    if (!("High" %in% colnames(pathway_stats_input))) {
      pathway_stats_input$High <- NA_real_
    }

    pathway_stats <- pathway_stats_input %>%
      filter(!is.na(Low), !is.na(High)) %>%
      group_by(pathway) %>%
      group_modify(~{
        p_value <- wilcox.test(
          x = .x$High,
          y = .x$Low,
          paired = TRUE,
          exact = FALSE
        )$p.value

        tibble(
          n_patients = nrow(.x),
          mean_low = mean(.x$Low),
          mean_high = mean(.x$High),
          delta_high_minus_low = mean(.x$High - .x$Low),
          median_delta = median(.x$High - .x$Low),
          p_value = p_value
        )
      }) %>%
      ungroup() %>%
      mutate(
        padj_BH = p.adjust(p_value, method = "BH"),
        CellType2 = subtype,
        gene = current_kcn,
        n_deg_methods = current_support,
        pathway_label = format_hallmark_name(pathway),
        direction = ifelse(delta_high_minus_low >= 0, "KCN-high", "KCN-low"),
        significant = padj_BH < aucell_padj_cutoff
      ) %>%
      arrange(desc(abs(delta_high_minus_low)))

    aucell_stats_all[[length(aucell_stats_all) + 1]] <- pathway_stats

    significant_pathways <- pathway_stats %>%
      filter(significant) %>%
      arrange(delta_high_minus_low)

    if (nrow(significant_pathways) == 0) {
      next
    }

    aucell_stats_significant[[length(aucell_stats_significant) + 1]] <- significant_pathways
  }
}

subtype_hallmark_sizes_tbl <- bind_rows(subtype_hallmark_sizes)
aucell_split_summary_tbl <- bind_rows(aucell_split_summary)

if (length(aucell_stats_all) > 0) {
  aucell_stats_all_tbl <- bind_rows(aucell_stats_all) %>%
    arrange(CellType2, gene, padj_BH, desc(abs(delta_high_minus_low)))
} else {
  aucell_stats_all_tbl <- tibble()
}

if (length(aucell_stats_significant) > 0) {
  aucell_stats_sig_tbl <- bind_rows(aucell_stats_significant) %>%
    arrange(CellType2, gene, padj_BH, desc(abs(delta_high_minus_low)))
} else {
  aucell_stats_sig_tbl <- tibble()
}

save_tsv(
  subtype_hallmark_sizes_tbl,
  file.path(aucell_table_dir, "AUCell_hallmark_gene_set_sizes_by_subtype.tsv")
)

save_tsv(
  aucell_split_summary_tbl,
  file.path(aucell_table_dir, "AUCell_KCN_split_summary.tsv")
)

save_tsv(
  aucell_stats_all_tbl,
  file.path(aucell_table_dir, "AUCell_hallmark_statistics_all.tsv")
)

save_tsv(
  aucell_stats_sig_tbl,
  file.path(aucell_table_dir, "AUCell_hallmark_statistics_significant.tsv")
)

###############################################################################
# 11. AUCell figures integrated in the main GSVA workflow
###############################################################################

selected_examples <- tibble(
  CellType2 = c("qPSC", "qPSC", "smPSC", "smPSC"),
  KCN = c("KCNJ8", "KCNJ8", "KCNMB1", "KCNMB1"),
  pathway = c(
    "HALLMARK_ANGIOGENESIS",
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
    "HALLMARK_MYOGENESIS",
    "HALLMARK_PI3K_AKT_MTOR_SIGNALING"
  )
) %>%
  mutate(
    pathway_label = format_hallmark_name(pathway),
    example_label = paste(CellType2, KCN, pathway_label, sep = " | ")
  )

selected_examples <- selected_examples %>%
  left_join(
    aucell_stats_sig_tbl %>%
      select(CellType2, gene, pathway, delta_high_minus_low, padj_BH, direction),
    by = c("CellType2", "KCN" = "gene", "pathway")
  ) %>%
  filter(!is.na(padj_BH))

all_score_examples <- aucell_stats_sig_tbl %>%
  group_by(CellType2, gene) %>%
  arrange(padj_BH, desc(abs(delta_high_minus_low)), pathway_label, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    CellType2,
    KCN = gene,
    pathway,
    pathway_label,
    padj_BH,
    delta_high_minus_low,
    n_deg_methods,
    example_label = paste(CellType2, gene, pathway_label, sep = " | ")
  )

save_tsv(selected_examples, file.path(aucell_table_dir, "AUCell_selected_examples.tsv"))
save_tsv(all_score_examples, file.path(aucell_table_dir, "AUCell_all_score_UMAP_examples.tsv"))

detected_genes_tbl <- tibble(
  cell_id = colnames(counts_mat),
  n_detected_genes = Matrix::colSums(counts_mat > 0)
) %>%
  inner_join(
    meta %>%
      select(cell_id, CellType2) %>%
      filter(CellType2 %in% subtype_order),
    by = "cell_id"
  ) %>%
  mutate(CellType2 = factor(CellType2, levels = subtype_order))

auc_max_rank <- ceiling(nrow(counts_mat) * 0.05)

p_detection <- ggplot(
  detected_genes_tbl,
  aes(x = n_detected_genes, fill = CellType2)
) +
  geom_histogram(bins = 40, color = "white", alpha = 0.9) +
  geom_vline(xintercept = auc_max_rank, linetype = "dashed", linewidth = 0.8) +
  facet_wrap(~ CellType2, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = subtype_colors) +
  labs(
    title = "Detected genes per cell before AUCell thresholding",
    x = "Detected genes per cell",
    y = "Number of cells",
    fill = NULL
  ) +
  clean_classic_theme() +
  theme(legend.position = "none")

save_plot(
  p_detection,
  file.path(aucell_figure_dir, "Figure_01_detected_genes_per_cell.png"),
  width = 10,
  height = 7
)

threshold_summary <- list()
threshold_plot_data <- list()
score_plot_data <- list()

for (i in seq_len(nrow(selected_examples))) {
  current_subtype <- selected_examples$CellType2[i]
  current_kcn <- selected_examples$KCN[i]
  current_pathway <- selected_examples$pathway[i]
  current_cache <- aucell_cache[[current_subtype]]

  threshold_info <- get_selected_threshold(current_cache$cells_auc, current_pathway)
  threshold_value <- threshold_info$threshold
  auc_scores <- current_cache$auc_mat[current_pathway, ]
  auc_score_values <- as.numeric(auc_scores)
  auc_cell_ids <- names(auc_scores)

  threshold_summary[[i]] <- tibble(
    CellType2 = current_subtype,
    KCN = current_kcn,
    pathway = current_pathway,
    pathway_label = format_hallmark_name(current_pathway),
    threshold = threshold_value,
    threshold_rule = threshold_info$threshold_name,
    used_fallback_threshold = threshold_info$used_fallback,
    n_cells = length(auc_score_values),
    n_active_cells = sum(auc_score_values > threshold_value),
    pct_active_cells = 100 * mean(auc_score_values > threshold_value)
  )

  threshold_plot_data[[i]] <- tibble(
    CellType2 = current_subtype,
    KCN = current_kcn,
    pathway = current_pathway,
    pathway_label = format_hallmark_name(current_pathway),
    example_label = selected_examples$example_label[i],
    auc_score = auc_score_values,
    threshold = threshold_value
  )

  score_plot_data[[i]] <- tibble(
    cell_id = auc_cell_ids,
    CellType2 = current_subtype,
    KCN = current_kcn,
    pathway = current_pathway,
    pathway_label = format_hallmark_name(current_pathway),
    example_label = selected_examples$example_label[i],
    auc_score = auc_score_values,
    active = auc_score_values > threshold_value
  )
}

threshold_summary_tbl <- bind_rows(threshold_summary)
threshold_plot_tbl <- bind_rows(threshold_plot_data)
score_plot_tbl <- bind_rows(score_plot_data)

save_tsv(threshold_summary_tbl, file.path(aucell_table_dir, "AUCell_threshold_summary.tsv"))

p_thresholds <- ggplot(
  threshold_plot_tbl,
  aes(x = auc_score, fill = CellType2)
) +
  geom_histogram(bins = 35, color = "white", alpha = 0.9) +
  geom_vline(aes(xintercept = threshold), linetype = "dashed", linewidth = 0.8) +
  facet_wrap(~ example_label, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = subtype_colors) +
  labs(
    title = "Selected AUCell thresholds for representative stromal examples",
    x = "AUCell score",
    y = "Number of cells",
    fill = NULL
  ) +
  clean_classic_theme() +
  theme(legend.position = "none")

save_plot(
  p_thresholds,
  file.path(aucell_figure_dir, "Figure_02_AUCell_threshold_histograms.png"),
  width = 12,
  height = 8
)

umap_tbl <- as.data.frame(Embeddings(stroma, reduction = "umap")) %>%
  rownames_to_column("cell_id") %>%
  left_join(
    meta %>%
      select(cell_id, CellType2) %>%
      filter(CellType2 %in% subtype_order),
    by = "cell_id"
  )

all_score_plot_data <- list()
for (i in seq_len(nrow(all_score_examples))) {
  current_subtype <- all_score_examples$CellType2[i]
  current_pathway <- all_score_examples$pathway[i]
  current_cache <- aucell_cache[[current_subtype]]

  if (!(current_pathway %in% rownames(current_cache$auc_mat))) {
    next
  }

  auc_scores <- current_cache$auc_mat[current_pathway, ]
  all_score_plot_data[[length(all_score_plot_data) + 1]] <- tibble(
    cell_id = names(auc_scores),
    CellType2 = current_subtype,
    KCN = all_score_examples$KCN[i],
    pathway = current_pathway,
    pathway_label = all_score_examples$pathway_label[i],
    example_label = all_score_examples$example_label[i],
    auc_score = as.numeric(auc_scores)
  )
}

score_umap_tbl <- score_plot_tbl %>%
  left_join(umap_tbl, by = "cell_id")

all_score_umap_tbl <- bind_rows(all_score_plot_data) %>%
  left_join(umap_tbl, by = "cell_id")

p_umap_score <- ggplot() +
  geom_point(
    data = umap_tbl,
    aes(x = UMAP_1, y = UMAP_2),
    color = "grey88",
    size = 0.35,
    alpha = 0.7
  ) +
  geom_point(
    data = all_score_umap_tbl,
    aes(x = UMAP_1, y = UMAP_2, color = auc_score),
    size = 0.45,
    alpha = 0.9
  ) +
  facet_wrap(~ example_label, ncol = 4) +
  scale_color_gradient(low = "#FFF4E6", high = "#800E13") +
  labs(
    title = "UMAP colored by AUCell scores for all significant KCN contexts",
    x = "UMAP 1",
    y = "UMAP 2",
    color = "AUC"
  ) +
  clean_classic_theme()

n_score_panels <- max(1, dplyr::n_distinct(all_score_umap_tbl$example_label))

save_plot(
  p_umap_score,
  file.path(aucell_figure_dir, "Figure_03_AUCell_score_UMAP_examples.png"),
  width = 18,
  height = max(10, ceiling(n_score_panels / 4) * 3.0)
)

p_umap_binary <- ggplot() +
  geom_point(
    data = umap_tbl,
    aes(x = UMAP_1, y = UMAP_2),
    color = "grey88",
    size = 0.35,
    alpha = 0.7
  ) +
  geom_point(
    data = score_umap_tbl,
    aes(x = UMAP_1, y = UMAP_2, color = active),
    size = 0.45,
    alpha = 0.95
  ) +
  facet_wrap(~ example_label, ncol = 2) +
  scale_color_manual(values = c("TRUE" = "#800E13", "FALSE" = "#70D6FF")) +
  labs(
    title = "UMAP colored by AUCell threshold assignments",
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Active"
  ) +
  clean_classic_theme()

save_plot(
  p_umap_binary,
  file.path(aucell_figure_dir, "Figure_04_AUCell_threshold_UMAP_examples.png"),
  width = 12,
  height = 8
)

comparison_examples <- tibble(
  CellType2 = c("qPSC", "smPSC"),
  KCN = c("KCNJ8", "KCNMB1"),
  pathway = c("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_MYOGENESIS")
)

comparison_plot_tbl <- list()
for (i in seq_len(nrow(comparison_examples))) {
  current_subtype <- comparison_examples$CellType2[i]
  current_kcn <- comparison_examples$KCN[i]
  current_pathway <- comparison_examples$pathway[i]
  current_cache <- aucell_cache[[current_subtype]]
  current_cells <- colnames(current_cache$norm)
  current_genes <- intersect(hallmark_gene_sets[[current_pathway]], rownames(current_cache$norm))

  mean_expression <- colMeans(current_cache$norm[current_genes, current_cells, drop = FALSE])
  auc_scores <- score_plot_tbl %>%
    filter(CellType2 == current_subtype, pathway == current_pathway) %>%
    select(cell_id, auc_score)

  comparison_tbl <- tibble(
    cell_id = current_cells,
    mean_expression = mean_expression,
    pathway = current_pathway
  ) %>%
    left_join(auc_scores, by = "cell_id") %>%
    left_join(umap_tbl, by = "cell_id")

  comparison_plot_tbl[[length(comparison_plot_tbl) + 1]] <- bind_rows(
    comparison_tbl %>%
      transmute(
        cell_id, UMAP_1, UMAP_2,
        CellType2 = current_subtype,
        KCN = current_kcn,
        pathway_label = format_hallmark_name(current_pathway),
        metric = "AUCell score",
        value = auc_score,
        example_label = paste(current_subtype, current_kcn, sep = " | ")
      ),
    comparison_tbl %>%
      transmute(
        cell_id, UMAP_1, UMAP_2,
        CellType2 = current_subtype,
        KCN = current_kcn,
        pathway_label = format_hallmark_name(current_pathway),
        metric = "Mean expression",
        value = mean_expression,
        example_label = paste(current_subtype, current_kcn, sep = " | ")
      )
  )
}

comparison_plot_tbl <- bind_rows(comparison_plot_tbl)

p_compare <- ggplot() +
  geom_point(
    data = umap_tbl,
    aes(x = UMAP_1, y = UMAP_2),
    color = "grey88",
    size = 0.35,
    alpha = 0.7
  ) +
  geom_point(
    data = comparison_plot_tbl,
    aes(x = UMAP_1, y = UMAP_2, color = value),
    size = 0.45,
    alpha = 0.95
  ) +
  facet_grid(metric ~ example_label) +
  scale_color_gradient(low = "#FFF4E6", high = "#800E13") +
  labs(
    title = "AUCell scores versus simple pathway mean expression",
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Value"
  ) +
  clean_classic_theme()

save_plot(
  p_compare,
  file.path(aucell_figure_dir, "Figure_05_AUCell_vs_mean_expression.png"),
  width = 12,
  height = 8
)

clustered_rows <- aucell_stats_sig_tbl %>%
  transmute(
    CellType2,
    gene,
    row_label = paste(CellType2, gene, sep = " | ")
  ) %>%
  distinct()

cluster_pathways <- aucell_stats_sig_tbl %>%
  count(pathway_label, sort = TRUE) %>%
  pull(pathway_label)

cluster_matrix_tbl <- aucell_stats_all_tbl %>%
  mutate(row_label = paste(CellType2, gene, sep = " | ")) %>%
  filter(row_label %in% clustered_rows$row_label, pathway_label %in% cluster_pathways) %>%
  select(row_label, pathway_label, delta_high_minus_low) %>%
  pivot_wider(names_from = pathway_label, values_from = delta_high_minus_low, values_fill = 0)

cluster_matrix <- as.matrix(cluster_matrix_tbl[, -1, drop = FALSE])
rownames(cluster_matrix) <- cluster_matrix_tbl$row_label

row_hclust <- hclust(dist(cluster_matrix, method = "euclidean"), method = "ward.D2")
row_order <- row_hclust$labels[row_hclust$order]

dotplot_tbl <- aucell_stats_sig_tbl %>%
  mutate(
    row_label = paste(CellType2, gene, sep = " | "),
    dot_size = pmin(-log10(padj_BH), 12)
  ) %>%
  filter(row_label %in% row_order, pathway_label %in% cluster_pathways) %>%
  mutate(
    row_label = factor(row_label, levels = rev(row_order)),
    pathway_label = factor(pathway_label, levels = cluster_pathways)
  )

save_tsv(
  tibble(row_label = row_order),
  file.path(aucell_table_dir, "AUCell_clustered_KCN_row_order.tsv")
)

save_tsv(
  dotplot_tbl %>%
    select(CellType2, gene, row_label, pathway_label, delta_high_minus_low, padj_BH, direction, dot_size),
  file.path(aucell_table_dir, "AUCell_clustered_KCN_dotplot_data.tsv")
)

p_clustered_dotplot <- ggplot(
  dotplot_tbl,
  aes(x = pathway_label, y = row_label, size = dot_size, color = delta_high_minus_low)
) +
  geom_point(alpha = 0.95) +
  scale_size_continuous(
    range = c(3.5, 10.5),
    breaks = c(2, 4, 6, 8, 10),
    name = "-log10 adjusted p-value"
  ) +
  scale_color_gradient2(
    low = "#70D6FF",
    mid = "#F7F7F7",
    high = "#800E13",
    midpoint = 0,
    name = "Delta AUCell\n(High - Low)"
  ) +
  labs(
    title = "Clustered AUCell dotplot of KCN-associated pathway outputs",
    x = NULL,
    y = NULL
  ) +
  clean_classic_theme() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, face = "bold", size = 9),
    axis.text.y = element_text(face = "bold", size = 9),
    legend.position = "right"
  )

save_plot(
  p_clustered_dotplot,
  file.path(aucell_figure_dir, "Figure_06_AUCell_clustered_KCN_dotplot.png"),
  width = 18,
  height = 11
)

###############################################################################
# 12. Session info
###############################################################################

message("GSVA Activated vs Normal analysis finished.")
message("Outputs saved in: ", output_root)

