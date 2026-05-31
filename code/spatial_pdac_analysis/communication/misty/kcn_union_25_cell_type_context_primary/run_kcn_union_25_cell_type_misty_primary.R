#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(mistyR)
})

options(stringsAsFactors = FALSE)

TARGET_ASSAY <- "SCT"
TARGET_LAYER <- "data"
PREDICTOR_ASSAY <- "rctd_multi1"
PREDICTOR_LAYER <- "data"
PRIMARY_ORIGIN <- "Pancreas"
JUXTA_NEIGHBOR_THRESHOLD <- 5
PARA_L <- 10
PARA_ZOI <- 5
CV_FOLDS <- 5
SEED <- 42
TOP_NETWORK_EDGES <- 3

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  file_path <- sub(file_arg, "", args[grep(file_arg, args)])
  if (length(file_path) == 0) {
    return(normalizePath(getwd()))
  }
  normalizePath(dirname(file_path))
}

find_project_dir <- function(start_dir) {
  current_dir <- normalizePath(start_dir)
  repeat {
    if (basename(current_dir) == "spatial_pdac_analysis") {
      return(current_dir)
    }
    parent_dir <- dirname(current_dir)
    if (identical(parent_dir, current_dir)) {
      stop("Could not find spatial_pdac_analysis from script location.")
    }
    current_dir <- parent_dir
  }
}

find_repo_dir <- function(start_dir) {
  current_dir <- normalizePath(start_dir)
  repeat {
    if (file.exists(file.path(current_dir, "M1_Bioinformatics_Ion_Channel", "data", "PDAC_Updated_ST.rds"))) {
      return(current_dir)
    }
    parent_dir <- dirname(current_dir)
    if (identical(parent_dir, current_dir)) {
      stop("Could not find repository root from script location.")
    }
    current_dir <- parent_dir
  }
}

write_tsv_simple <- function(x, path) {
  write.table(
    x,
    file = path,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE,
    na = "NA"
  )
}

get_assay_matrix_compat <- function(seurat_obj, assay, layer_name = "data") {
  tryCatch(
    GetAssayData(seurat_obj, assay = assay, layer = layer_name),
    error = function(e) GetAssayData(seurat_obj, assay = assay, slot = layer_name)
  )
}

sanitize_names <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

safe_numeric <- function(x) {
  as.numeric(as.vector(x))
}

make_target_section_views <- function(
  seurat_obj,
  section_name,
  target_gene,
  target_mat,
  predictor_mat,
  predictor_lookup
) {
  coords_df <- seurat_obj@images[[section_name]]@coordinates |>
    rownames_to_column("spot") |>
    filter(spot %in% colnames(target_mat)) |>
    select(spot, row, col, imagerow, imagecol)

  spots <- coords_df$spot
  target_values <- safe_numeric(target_mat[target_gene, spots, drop = TRUE])
  predictor_values <- t(as.matrix(predictor_mat[predictor_lookup$original_name, spots, drop = FALSE]))
  colnames(predictor_values) <- predictor_lookup$sanitized_name
  predictor_values <- as.data.frame(predictor_values, check.names = FALSE, stringsAsFactors = FALSE)

  varying_predictors <- vapply(
    predictor_values,
    function(x) {
      x <- as.numeric(x)
      is.finite(stats::sd(x, na.rm = TRUE)) && stats::sd(x, na.rm = TRUE) > 0
    },
    logical(1)
  )
  predictor_values <- predictor_values[, varying_predictors, drop = FALSE]

  section_info <- data.frame(
    target_gene = target_gene,
    orig.ident = section_name,
    n_spots = length(spots),
    n_positive_spots = sum(target_values > 0, na.rm = TRUE),
    mean_target = mean(target_values, na.rm = TRUE),
    sd_target = stats::sd(target_values, na.rm = TRUE),
    n_predictors_used = ncol(predictor_values),
    predictors_removed = paste(
      predictor_lookup$original_name[!predictor_lookup$sanitized_name %in% colnames(predictor_values)],
      collapse = "; "
    ),
    stringsAsFactors = FALSE
  )

  if (!is.finite(section_info$sd_target) || section_info$sd_target == 0) {
    return(list(
      usable = FALSE,
      reason = paste0(target_gene, " has zero variance in this section."),
      section_info = section_info
    ))
  }

  if (ncol(predictor_values) == 0) {
    return(list(
      usable = FALSE,
      reason = "No varying cell-type predictor left after section-level filtering.",
      section_info = section_info
    ))
  }

  intraview_df <- data.frame(
    stringsAsFactors = FALSE,
    check.names = FALSE,
    setNames(data.frame(target_values, check.names = FALSE), target_gene),
    predictor_values
  )

  predictor_cols <- colnames(predictor_values)
  positions <- coords_df |> select(row, col)

  views <- create_initial_view(intraview_df, unique.id = paste(section_name, target_gene, sep = "__")) |>
    add_juxtaview(positions, neighbor.thr = JUXTA_NEIGHBOR_THRESHOLD, verbose = FALSE) |>
    add_paraview(positions, l = PARA_L, zoi = PARA_ZOI, family = "gaussian", verbose = FALSE) |>
    select_markers(
      paste0("juxtaview.", JUXTA_NEIGHBOR_THRESHOLD),
      all_of(predictor_cols)
    ) |>
    select_markers(
      paste0("paraview.", PARA_L),
      all_of(predictor_cols)
    ) |>
    rename_view(
      paste0("juxtaview.", JUXTA_NEIGHBOR_THRESHOLD),
      "juxtaview_cell_types",
      "juxta_ct"
    ) |>
    rename_view(
      paste0("paraview.", PARA_L),
      "paraview_cell_types",
      "para_ct"
    )

  list(
    usable = TRUE,
    views = views,
    section_info = section_info
  )
}

read_performance <- function(sample_dir, target_gene, section_name) {
  x <- read.table(file.path(sample_dir, "performance.txt"), header = TRUE, check.names = FALSE)
  tibble(
    target_gene = target_gene,
    orig.ident = section_name,
    intra.RMSE = x$`intra.RMSE`,
    intra.R2 = x$`intra.R2`,
    multi.RMSE = x$`multi.RMSE`,
    multi.R2 = x$`multi.R2`,
    p.RMSE = x$`p.RMSE`,
    p.R2 = x$`p.R2`,
    gain.RMSE = x$`multi.RMSE` - x$`intra.RMSE`,
    gain.R2 = x$`multi.R2` - x$`intra.R2`
  )
}

read_coefficients <- function(sample_dir, target_gene, section_name) {
  x <- read.table(file.path(sample_dir, "coefficients.txt"), header = TRUE, check.names = FALSE)
  tibble(
    target_gene = target_gene,
    orig.ident = section_name,
    view = c("intercept", "intra", "juxta_ct", "para_ct"),
    value = c(x$intercept, x$intra, x$juxta_ct, x$para_ct),
    p_value = c(NA_real_, x$`p.intra`, x$`p.juxta_ct`, x$`p.para_ct`)
  )
}

read_importances <- function(sample_dir, target_gene, section_name) {
  files <- list.files(sample_dir, pattern = "^importances_.*\\.txt$", full.names = TRUE)
  bind_rows(lapply(files, function(path) {
    file_name <- basename(path)
    view <- sub("^importances_[^_]+_", "", file_name)
    view <- sub("\\.txt$", "", view)
    df <- read.csv(path, stringsAsFactors = FALSE)
    tibble(
      target_gene = target_gene,
      orig.ident = section_name,
      view = view,
      Predictor = df$target,
      Importance = df$imp
    )
  }))
}

script_dir <- get_script_dir()
project_dir <- find_project_dir(script_dir)
repo_dir <- find_repo_dir(script_dir)
input_dir <- file.path(project_dir, "inputs")
union_path <- file.path(repo_dir, "M1_Bioinformatics_Ion_Channel", "code", "mycaf_icaf_correlations", "kcn_union_25.tsv")

out_dir <- script_dir
results_dir <- file.path(out_dir, "results")
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
unlink(results_dir, recursive = TRUE, force = TRUE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

object_path <- file.path(input_dir, "PDAC_Updated_ST.rds")
if (!file.exists(object_path)) {
  stop("Missing input object: ", object_path)
}
if (!file.exists(union_path)) {
  stop("Missing KCN union file: ", union_path)
}

message("Loading spatial object...")
st_obj <- readRDS(object_path)
target_mat <- get_assay_matrix_compat(st_obj, assay = TARGET_ASSAY, layer_name = TARGET_LAYER)
predictor_mat <- get_assay_matrix_compat(st_obj, assay = PREDICTOR_ASSAY, layer_name = PREDICTOR_LAYER)

target_genes <- read.delim(union_path, stringsAsFactors = FALSE)$gene
target_genes_present <- intersect(target_genes, rownames(target_mat))
target_genes_missing <- setdiff(target_genes, rownames(target_mat))

write_tsv_simple(
  data.frame(target_gene = target_genes_present, stringsAsFactors = FALSE),
  file.path(tab_dir, "kcn_union_25_targets_present_in_sct.tsv")
)
write_tsv_simple(
  data.frame(target_gene = target_genes_missing, stringsAsFactors = FALSE),
  file.path(tab_dir, "kcn_union_25_targets_missing_in_sct.tsv")
)

primary_sections <- st_obj@meta.data |>
  rownames_to_column("spot") |>
  filter(Origin == PRIMARY_ORIGIN) |>
  distinct(orig.ident) |>
  arrange(orig.ident) |>
  pull(orig.ident)

predictor_lookup <- tibble(
  original_name = rownames(predictor_mat),
  sanitized_name = sanitize_names(rownames(predictor_mat))
)
write_tsv_simple(predictor_lookup, file.path(tab_dir, "kcn_union_25_misty_predictor_lookup.tsv"))

section_info_list <- list()
status_list <- list()
result_index <- list()

for (target_gene in target_genes_present) {
  message("Processing target ", target_gene, " ...")
  for (section_name in primary_sections) {
    section_payload <- make_target_section_views(
      seurat_obj = st_obj,
      section_name = section_name,
      target_gene = target_gene,
      target_mat = target_mat,
      predictor_mat = predictor_mat,
      predictor_lookup = predictor_lookup
    )

    section_info_list[[paste(target_gene, section_name, sep = "__")]] <- section_payload$section_info

    if (!isTRUE(section_payload$usable)) {
      status_list[[paste(target_gene, section_name, sep = "__")]] <- data.frame(
        target_gene = target_gene,
        orig.ident = section_name,
        status = "skipped",
        message = section_payload$reason,
        stringsAsFactors = FALSE
      )
      next
    }

    section_results_dir <- file.path(results_dir, target_gene, section_name)
    dir.create(section_results_dir, recursive = TRUE, showWarnings = FALSE)

    message("Running MISTy for ", target_gene, " / ", section_name, " ...")
    section_result <- tryCatch(
      run_misty(
        section_payload$views,
        results.folder = section_results_dir,
        seed = SEED,
        target.subset = target_gene,
        cv.folds = CV_FOLDS,
        model.function = random_forest_model,
        num.trees = 100,
        importance = "impurity",
        num.threads = 1
      ),
      error = function(e) e
    )

    if (inherits(section_result, "error")) {
      status_list[[paste(target_gene, section_name, sep = "__")]] <- data.frame(
        target_gene = target_gene,
        orig.ident = section_name,
        status = "failed",
        message = conditionMessage(section_result),
        stringsAsFactors = FALSE
      )
      next
    }

    status_list[[paste(target_gene, section_name, sep = "__")]] <- data.frame(
      target_gene = target_gene,
      orig.ident = section_name,
      status = "completed",
      message = NA_character_,
      stringsAsFactors = FALSE
    )

    result_index[[paste(target_gene, section_name, sep = "__")]] <- data.frame(
      target_gene = target_gene,
      orig.ident = section_name,
      result_dir = section_results_dir,
      stringsAsFactors = FALSE
    )
  }
}

section_info_df <- bind_rows(section_info_list)
status_df <- bind_rows(status_list)
result_index_df <- bind_rows(result_index)

write_tsv_simple(section_info_df, file.path(tab_dir, "kcn_union_25_misty_section_info.tsv"))
write_tsv_simple(status_df, file.path(tab_dir, "kcn_union_25_misty_section_status.tsv"))
write_tsv_simple(result_index_df, file.path(tab_dir, "kcn_union_25_misty_result_index.tsv"))

if (nrow(result_index_df) == 0) {
  stop("No successful MISTy result was generated.")
}

message("Collecting raw MISTy outputs section by section...")
performance_wide_df <- bind_rows(lapply(seq_len(nrow(result_index_df)), function(i) {
  read_performance(result_index_df$result_dir[i], result_index_df$target_gene[i], result_index_df$orig.ident[i])
}))

contributions_df <- bind_rows(lapply(seq_len(nrow(result_index_df)), function(i) {
  read_coefficients(result_index_df$result_dir[i], result_index_df$target_gene[i], result_index_df$orig.ident[i])
}))

importances_df <- bind_rows(lapply(seq_len(nrow(result_index_df)), function(i) {
  read_importances(result_index_df$result_dir[i], result_index_df$target_gene[i], result_index_df$orig.ident[i])
})) |>
  mutate(
    Predictor_original = predictor_lookup$original_name[match(Predictor, predictor_lookup$sanitized_name)]
  )

improvements_df <- performance_wide_df |>
  tidyr::pivot_longer(
    cols = c(`intra.RMSE`, `intra.R2`, `multi.RMSE`, `multi.R2`, `p.RMSE`, `p.R2`, `gain.RMSE`, `gain.R2`),
    names_to = "measure",
    values_to = "value"
  )

view_fraction_df <- contributions_df |>
  filter(view != "intercept") |>
  group_by(target_gene, orig.ident) |>
  mutate(abs_fraction = abs(value) / sum(abs(value), na.rm = TRUE)) |>
  ungroup()

improvements_summary_df <- improvements_df |>
  group_by(target_gene, measure) |>
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd = stats::sd(value, na.rm = TRUE),
    cv = ifelse(abs(mean) < .Machine$double.eps, NA_real_, sd / mean),
    .groups = "drop"
  )

contributions_summary_df <- contributions_df |>
  filter(view != "intercept") |>
  group_by(target_gene, view) |>
  summarise(
    mean = mean(value, na.rm = TRUE),
    p.mean = mean(p_value, na.rm = TRUE),
    p.sd = stats::sd(p_value, na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(
    view_fraction_df |>
      group_by(target_gene, view) |>
      summarise(abs_fraction = mean(abs_fraction, na.rm = TRUE), .groups = "drop"),
    by = c("target_gene", "view")
  ) |>
  select(target_gene, view, mean, abs_fraction, p.mean, p.sd)

importances_summary_df <- importances_df |>
  group_by(target_gene, view, Predictor, Predictor_original) |>
  summarise(
    Importance = mean(Importance, na.rm = TRUE),
    nsamples = dplyr::n(),
    .groups = "drop"
  )

top_predictors_df <- importances_summary_df |>
  group_by(target_gene, view) |>
  arrange(target_gene, view, desc(Importance), .by_group = TRUE) |>
  slice_head(n = 5) |>
  ungroup()

network_edges_df <- importances_summary_df |>
  group_by(target_gene, view) |>
  mutate(
    importance_rank = dplyr::row_number(desc(Importance)),
    normalized_importance = Importance / max(Importance, na.rm = TRUE)
  ) |>
  ungroup()

network_edges_top_df <- network_edges_df |>
  filter(importance_rank <= TOP_NETWORK_EDGES)

dominant_view_df <- contributions_summary_df |>
  group_by(target_gene) |>
  arrange(desc(abs_fraction), .by_group = TRUE) |>
  slice_head(n = 1) |>
  ungroup() |>
  rename(dominant_view = view, dominant_view_abs_fraction = abs_fraction, dominant_view_mean = mean)

dominant_predictor_by_view_df <- importances_summary_df |>
  group_by(target_gene, view) |>
  arrange(desc(Importance), .by_group = TRUE) |>
  slice_head(n = 1) |>
  ungroup() |>
  rename(
    dominant_predictor = Predictor_original,
    dominant_importance = Importance,
    dominant_predictor_nsamples = nsamples
  )

schematic_summary_df <- dominant_view_df |>
  left_join(
    dominant_predictor_by_view_df |>
      select(target_gene, view, dominant_predictor, dominant_importance),
    by = c("target_gene", "dominant_view" = "view")
  ) |>
  left_join(
    improvements_summary_df |>
      filter(measure == "gain.R2") |>
      select(target_gene, gain_R2_mean = mean),
    by = "target_gene"
  )

write_tsv_simple(improvements_df, file.path(tab_dir, "kcn_union_25_misty_improvements.tsv"))
write_tsv_simple(performance_wide_df, file.path(tab_dir, "kcn_union_25_misty_performance_wide.tsv"))
write_tsv_simple(contributions_df, file.path(tab_dir, "kcn_union_25_misty_contributions.tsv"))
write_tsv_simple(view_fraction_df, file.path(tab_dir, "kcn_union_25_misty_view_fractions_by_section.tsv"))
write_tsv_simple(importances_df, file.path(tab_dir, "kcn_union_25_misty_importances.tsv"))
write_tsv_simple(improvements_summary_df, file.path(tab_dir, "kcn_union_25_misty_improvements_summary.tsv"))
write_tsv_simple(contributions_summary_df, file.path(tab_dir, "kcn_union_25_misty_contributions_summary.tsv"))
write_tsv_simple(importances_summary_df, file.path(tab_dir, "kcn_union_25_misty_importances_summary.tsv"))
write_tsv_simple(top_predictors_df, file.path(tab_dir, "kcn_union_25_misty_top_predictors_by_view.tsv"))
write_tsv_simple(network_edges_df, file.path(tab_dir, "kcn_union_25_misty_network_edges.tsv"))
write_tsv_simple(network_edges_top_df, file.path(tab_dir, "kcn_union_25_misty_network_edges_top3.tsv"))
write_tsv_simple(dominant_view_df, file.path(tab_dir, "kcn_union_25_misty_dominant_view.tsv"))
write_tsv_simple(dominant_predictor_by_view_df, file.path(tab_dir, "kcn_union_25_misty_dominant_predictor_by_view.tsv"))
write_tsv_simple(schematic_summary_df, file.path(tab_dir, "kcn_union_25_misty_schematic_summary.tsv"))

contrib_heatmap_df <- contributions_summary_df |>
  mutate(
    view = factor(view, levels = c("intra", "juxta_ct", "para_ct")),
    target_gene = factor(
      target_gene,
      levels = dominant_view_df |>
        arrange(desc(dominant_view_abs_fraction), target_gene) |>
        pull(target_gene)
    )
  )

contrib_heatmap <- ggplot(contrib_heatmap_df, aes(x = view, y = target_gene, fill = abs_fraction)) +
  geom_tile(color = "white", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.2f", abs_fraction)), size = 2.5) +
  scale_fill_viridis_c(option = "magma", end = 0.95) +
  labs(
    title = "25-KCN MISTy view fractions across primary PDAC sections",
    subtitle = "Target = KCN, predictors = rctd_multi1 cell types",
    x = NULL,
    y = NULL,
    fill = "Mean absolute\nfraction"
  ) +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold"))

ggsave(
  filename = file.path(fig_dir, "kcn_union_25_misty_view_fraction_heatmap.pdf"),
  plot = contrib_heatmap,
  width = 8,
  height = 10
)

top3_network_plot_df <- network_edges_top_df |>
  mutate(
    view = factor(view, levels = c("intra", "juxta_ct", "para_ct")),
    target_gene = factor(target_gene, levels = unique(target_gene)),
    Predictor_original = factor(Predictor_original, levels = rev(unique(Predictor_original)))
  )

network_heatmap <- ggplot(top3_network_plot_df, aes(x = Predictor_original, y = target_gene, fill = normalized_importance)) +
  geom_tile(color = "white", linewidth = 0.2) +
  facet_wrap(~ view, ncol = 1, scales = "free_x") +
  scale_fill_viridis_c(option = "plasma", end = 0.95) +
  labs(
    title = "Top cell-type neighbors per KCN and MISTy view",
    subtitle = "Top 3 cell-type edges retained per target-view",
    x = NULL,
    y = NULL,
    fill = "Normalized\nimportance"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

ggsave(
  filename = file.path(fig_dir, "kcn_union_25_misty_top3_network_heatmap.pdf"),
  plot = network_heatmap,
  width = 10,
  height = 12
)

message("25-KCN cell-type-context MISTy run complete.")
