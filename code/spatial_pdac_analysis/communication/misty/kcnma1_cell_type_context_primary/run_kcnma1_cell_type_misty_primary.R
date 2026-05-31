#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(mistyR)
})

options(stringsAsFactors = FALSE)

TARGET_GENE <- "KCNMA1"
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

make_section_views <- function(seurat_obj, section_name, target_mat, predictor_mat, predictor_lookup) {
  spots <- rownames(seurat_obj@images[[section_name]]@coordinates)
  coords_df <- seurat_obj@images[[section_name]]@coordinates |>
    rownames_to_column("spot") |>
    filter(spot %in% colnames(target_mat)) |>
    select(spot, row, col, imagerow, imagecol)

  spots <- coords_df$spot
  target_values <- safe_numeric(target_mat[TARGET_GENE, spots, drop = TRUE])
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

  intraview_df <- data.frame(
    KCNMA1 = target_values,
    predictor_values,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  # Skip sections where KCNMA1 is constant because MISTy cannot model a target
  # without variance.
  if (stats::sd(intraview_df$KCNMA1, na.rm = TRUE) == 0 || all(!is.finite(intraview_df$KCNMA1))) {
    return(list(
      usable = FALSE,
      reason = "KCNMA1 has zero variance in this section.",
      section_info = data.frame(
        orig.ident = section_name,
        n_spots = length(spots),
        n_positive_spots = sum(intraview_df$KCNMA1 > 0, na.rm = TRUE),
        mean_kcnma1 = mean(intraview_df$KCNMA1, na.rm = TRUE),
        sd_kcnma1 = stats::sd(intraview_df$KCNMA1, na.rm = TRUE),
        n_predictors_used = ncol(predictor_values),
        predictors_removed = paste(
          predictor_lookup$original_name[!predictor_lookup$sanitized_name %in% colnames(predictor_values)],
          collapse = "; "
        ),
        stringsAsFactors = FALSE
      )
    ))
  }

  if (ncol(predictor_values) == 0) {
    return(list(
      usable = FALSE,
      reason = "No varying cell-type predictor left after section-level filtering.",
      section_info = data.frame(
        orig.ident = section_name,
        n_spots = length(spots),
        n_positive_spots = sum(intraview_df$KCNMA1 > 0, na.rm = TRUE),
        mean_kcnma1 = mean(intraview_df$KCNMA1, na.rm = TRUE),
        sd_kcnma1 = stats::sd(intraview_df$KCNMA1, na.rm = TRUE),
        n_predictors_used = 0,
        predictors_removed = paste(predictor_lookup$original_name, collapse = "; "),
        stringsAsFactors = FALSE
      )
    ))
  }

  positions <- coords_df |>
    select(row, col)

  views <- create_initial_view(intraview_df, unique.id = section_name) |>
    add_juxtaview(positions, neighbor.thr = JUXTA_NEIGHBOR_THRESHOLD, verbose = FALSE) |>
    add_paraview(positions, l = PARA_L, zoi = PARA_ZOI, family = "gaussian", verbose = FALSE) |>
    select_markers(
      paste0("juxtaview.", JUXTA_NEIGHBOR_THRESHOLD),
      -all_of(TARGET_GENE)
    ) |>
    select_markers(
      paste0("paraview.", PARA_L),
      -all_of(TARGET_GENE)
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
    section_info = data.frame(
      orig.ident = section_name,
      n_spots = length(spots),
      n_positive_spots = sum(intraview_df$KCNMA1 > 0, na.rm = TRUE),
      mean_kcnma1 = mean(intraview_df$KCNMA1, na.rm = TRUE),
      sd_kcnma1 = stats::sd(intraview_df$KCNMA1, na.rm = TRUE),
      n_predictors_used = ncol(predictor_values),
      predictors_removed = paste(
        predictor_lookup$original_name[!predictor_lookup$sanitized_name %in% colnames(predictor_values)],
        collapse = "; "
      ),
      stringsAsFactors = FALSE
    )
  )
}

clean_sample_labels <- function(x) {
  basename(normalizePath(x, winslash = "/", mustWork = FALSE))
}

script_dir <- get_script_dir()
project_dir <- find_project_dir(script_dir)
input_dir <- file.path(project_dir, "inputs")
out_dir <- script_dir
results_dir <- file.path(out_dir, "results")
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
unlink(results_dir, recursive = TRUE, force = TRUE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

object_path <- file.path(input_dir, "PDAC_Updated_ST.rds")
if (!file.exists(object_path)) {
  stop("Missing input object: ", object_path)
}

message("Loading spatial object...")
st_obj <- readRDS(object_path)

target_mat <- get_assay_matrix_compat(st_obj, assay = TARGET_ASSAY, layer_name = TARGET_LAYER)
predictor_mat <- get_assay_matrix_compat(st_obj, assay = PREDICTOR_ASSAY, layer_name = PREDICTOR_LAYER)

if (!TARGET_GENE %in% rownames(target_mat)) {
  stop(TARGET_GENE, " is not available in assay ", TARGET_ASSAY, ".")
}

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

write_tsv_simple(
  predictor_lookup,
  file.path(tab_dir, "kcnma1_misty_predictor_lookup.tsv")
)

section_info_list <- list()
misty_result_dirs <- character(0)
failed_sections <- list()

for (section_name in primary_sections) {
  message("Preparing section ", section_name, " ...")
  section_payload <- make_section_views(
    seurat_obj = st_obj,
    section_name = section_name,
    target_mat = target_mat,
    predictor_mat = predictor_mat,
    predictor_lookup = predictor_lookup
  )

  section_info_list[[section_name]] <- section_payload$section_info

  if (!isTRUE(section_payload$usable)) {
    failed_sections[[section_name]] <- data.frame(
      orig.ident = section_name,
      status = "skipped",
      message = section_payload$reason,
      stringsAsFactors = FALSE
    )
    next
  }

  message("Running MISTy for ", section_name, " ...")
  section_results_dir <- file.path(results_dir, section_name)
  dir.create(section_results_dir, recursive = TRUE, showWarnings = FALSE)
  section_result_dir <- tryCatch(
    run_misty(
      section_payload$views,
      results.folder = section_results_dir,
      seed = SEED,
      target.subset = TARGET_GENE,
      cv.folds = CV_FOLDS,
      model.function = random_forest_model,
      num.trees = 100,
      importance = "impurity",
      num.threads = 1
    ),
    error = function(e) e
  )

  if (inherits(section_result_dir, "error")) {
    failed_sections[[section_name]] <- data.frame(
      orig.ident = section_name,
      status = "failed",
      message = conditionMessage(section_result_dir),
      stringsAsFactors = FALSE
    )
    next
  }

  misty_result_dirs <- c(misty_result_dirs, section_results_dir)
  failed_sections[[section_name]] <- data.frame(
    orig.ident = section_name,
    status = "completed",
    message = NA_character_,
    stringsAsFactors = FALSE
  )
}

section_info_df <- bind_rows(section_info_list)
status_df <- bind_rows(failed_sections)

write_tsv_simple(section_info_df, file.path(tab_dir, "kcnma1_misty_section_info.tsv"))
write_tsv_simple(status_df, file.path(tab_dir, "kcnma1_misty_section_status.tsv"))

if (length(misty_result_dirs) == 0) {
  stop("No section completed successfully; no MISTy result to collect.")
}

message("Collecting raw MISTy outputs section by section...")

read_performance <- function(sample_dir) {
  x <- read.table(file.path(sample_dir, "performance.txt"), header = TRUE, check.names = FALSE)
  tibble(
    sample = basename(sample_dir),
    target = x$target,
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

read_coefficients <- function(sample_dir) {
  x <- read.table(file.path(sample_dir, "coefficients.txt"), header = TRUE, check.names = FALSE)
  tibble(
    sample = basename(sample_dir),
    target = x$target,
    view = c("intercept", "intra", "juxta_ct", "para_ct"),
    value = c(x$intercept, x$intra, x$juxta_ct, x$para_ct),
    p_value = c(NA_real_, x$`p.intra`, x$`p.juxta_ct`, x$`p.para_ct`)
  )
}

read_importances <- function(sample_dir) {
  sample_name <- basename(sample_dir)
  files <- list.files(sample_dir, pattern = "^importances_.*\\.txt$", full.names = TRUE)
  bind_rows(lapply(files, function(path) {
    file_name <- basename(path)
    view <- sub("^importances_[^_]+_", "", file_name)
    view <- sub("\\.txt$", "", view)
    df <- read.csv(path, stringsAsFactors = FALSE)
    tibble(
      sample = sample_name,
      view = view,
      Predictor = df$target,
      Target = TARGET_GENE,
      Importance = df$imp
    )
  }))
}

performance_wide_df <- bind_rows(lapply(misty_result_dirs, read_performance))
contributions_df <- bind_rows(lapply(misty_result_dirs, read_coefficients))
importances_df <- bind_rows(lapply(misty_result_dirs, read_importances)) |>
  mutate(
    Predictor_original = predictor_lookup$original_name[match(Predictor, predictor_lookup$sanitized_name)]
  )

improvements_df <- performance_wide_df |>
  tidyr::pivot_longer(
    cols = c(`intra.RMSE`, `intra.R2`, `multi.RMSE`, `multi.R2`, `p.RMSE`, `p.R2`, `gain.RMSE`, `gain.R2`),
    names_to = "measure",
    values_to = "value"
  )

improvements_stats_df <- improvements_df |>
  group_by(target, measure) |>
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd = stats::sd(value, na.rm = TRUE),
    cv = ifelse(abs(mean) < .Machine$double.eps, NA_real_, sd / mean),
    .groups = "drop"
  )

contributions_stats_df <- contributions_df |>
  filter(view != "intercept") |>
  group_by(target, view) |>
  summarise(
    mean = mean(value, na.rm = TRUE),
    p.mean = mean(p_value, na.rm = TRUE),
    p.sd = stats::sd(p_value, na.rm = TRUE),
    .groups = "drop"
  )

importances_agg_df <- importances_df |>
  group_by(view, Predictor, Target, Predictor_original) |>
  summarise(
    Importance = mean(Importance, na.rm = TRUE),
    nsamples = dplyr::n(),
    .groups = "drop"
  )

view_fraction_df <- contributions_df |>
  filter(view != "intercept") |>
  group_by(sample) |>
  mutate(fraction = value / sum(value, na.rm = TRUE)) |>
  ungroup()

contributions_stats_df <- contributions_stats_df |>
  left_join(
    view_fraction_df |>
      group_by(target, view) |>
      summarise(fraction = mean(fraction, na.rm = TRUE), .groups = "drop"),
    by = c("target", "view")
  ) |>
  select(target, view, mean, fraction, p.mean, p.sd)

top_predictors_df <- importances_agg_df |>
  group_by(view) |>
  arrange(view, desc(abs(Importance))) |>
  slice_head(n = 10) |>
  ungroup()

write_tsv_simple(improvements_df, file.path(tab_dir, "kcnma1_misty_improvements.tsv"))
write_tsv_simple(contributions_df, file.path(tab_dir, "kcnma1_misty_contributions.tsv"))
write_tsv_simple(view_fraction_df, file.path(tab_dir, "kcnma1_misty_view_fractions_by_section.tsv"))
write_tsv_simple(importances_df, file.path(tab_dir, "kcnma1_misty_importances.tsv"))
write_tsv_simple(improvements_stats_df, file.path(tab_dir, "kcnma1_misty_improvements_summary.tsv"))
write_tsv_simple(contributions_stats_df, file.path(tab_dir, "kcnma1_misty_contributions_summary.tsv"))
write_tsv_simple(importances_agg_df, file.path(tab_dir, "kcnma1_misty_importances_summary.tsv"))
write_tsv_simple(top_predictors_df, file.path(tab_dir, "kcnma1_misty_top_predictors_by_view.tsv"))

contrib_plot_df <- contributions_stats_df |>
  mutate(view = factor(view, levels = view[order(fraction, decreasing = TRUE)]))

contrib_plot <- ggplot(contrib_plot_df, aes(x = view, y = fraction, fill = view)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.2f", fraction)), vjust = -0.4, size = 3.3) +
  labs(
    title = "KCNMA1 MISTy mean view contribution across primary sections",
    subtitle = "Cell-type predictors from rctd_multi1; target = SCT KCNMA1",
    x = NULL,
    y = "Mean contribution fraction"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 20, hjust = 1)
  ) +
  coord_cartesian(ylim = c(0, max(contrib_plot_df$fraction, na.rm = TRUE) * 1.15))

ggsave(
  filename = file.path(fig_dir, "kcnma1_misty_view_contributions.pdf"),
  plot = contrib_plot,
  width = 8,
  height = 5
)

importance_plot_df <- importances_agg_df |>
  filter(!is.na(Predictor_original)) |>
  mutate(
    Predictor_original = factor(
      Predictor_original,
      levels = rev(unique(Predictor_original[order(view, Importance)]))
    ),
    view = factor(view, levels = unique(view))
  )

importance_plot <- ggplot(
  importance_plot_df,
  aes(x = view, y = Predictor_original, fill = Importance)
) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b") +
  labs(
    title = "KCNMA1 MISTy predictor importance by view",
    subtitle = "Positive and negative weighted importances aggregated across primary sections",
    x = NULL,
    y = NULL,
    fill = "Importance"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 20, hjust = 1)
  )

ggsave(
  filename = file.path(fig_dir, "kcnma1_misty_predictor_importance_heatmap.pdf"),
  plot = importance_plot,
  width = 8,
  height = 6
)

message("KCNMA1 MISTy run complete.")
