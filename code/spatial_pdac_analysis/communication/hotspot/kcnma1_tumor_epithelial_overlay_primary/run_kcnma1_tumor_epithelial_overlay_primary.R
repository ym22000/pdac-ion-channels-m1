#!/usr/bin/env Rscript

###############################################################################
# KCNMA1 tumor epithelial signature overlay in primary PDAC sections
#
# Signature used here:
#   EPCAM, KRT19 (CK19), MUC1, and a PanCK-like keratin block
#   approximated with KRT7 / KRT8 / KRT18.
#
# Notes:
#   - "PanCK" is an antibody marker, not a gene symbol.
#   - We therefore use a small epithelial keratin module instead.
###############################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(BiocParallel)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(escape)
})

options(stringsAsFactors = FALSE)

TARGET_GENE <- "KCNMA1"
RUN_ESCAPE_GROUPS <- 1000
HOTSPOT_COLOR_UPPER_QUANTILE <- 0.95
PATHWAY_COLOR_LOWER_QUANTILE <- 0.02
PATHWAY_COLOR_UPPER_QUANTILE <- 0.98
SPOT_POINT_SIZE <- 1.35

PATHWAY_GENE_SETS <- list(
  Tumor_Epithelial = c("EPCAM", "KRT19", "MUC1", "KRT7", "KRT8", "KRT18")
)

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

safe_update_seurat <- function(seurat_obj) {
  seurat_obj
}

get_assay_matrix_compat <- function(seurat_obj, assay, layer_name = "data") {
  tryCatch(
    GetAssayData(seurat_obj, assay = assay, layer = layer_name),
    error = function(e) GetAssayData(seurat_obj, assay = assay, slot = layer_name)
  )
}

load_section_payload <- function(section_dir, target_gene) {
  gene_names <- read.delim(
    file.path(section_dir, "genes.tsv"),
    header = FALSE,
    stringsAsFactors = FALSE
  )[[1]]
  spot_names <- read.delim(
    file.path(section_dir, "spots.tsv"),
    header = FALSE,
    stringsAsFactors = FALSE
  )[[1]]
  coords_df <- read.delim(file.path(section_dir, "coords.tsv"), check.names = FALSE)
  if (!"spot" %in% colnames(coords_df) && ncol(coords_df) >= 1) {
    colnames(coords_df)[1] <- "spot"
  }
  coords_df <- coords_df %>%
    select(spot, imagecol, imagerow)

  counts_mat <- Matrix::readMM(file.path(section_dir, "counts.mtx"))
  gene_idx <- match(target_gene, gene_names)
  if (is.na(gene_idx)) {
    expr <- rep(0, length(spot_names))
  } else {
    expr <- as.numeric(counts_mat[gene_idx, ])
  }

  tibble(
    spot = spot_names,
    expr = expr
  ) %>%
    left_join(coords_df, by = "spot") %>%
    mutate(orig.ident = basename(section_dir))
}

make_section_map <- function(section_dirs) {
  bind_rows(lapply(section_dirs, function(section_dir) {
    tibble(
      spot = read.delim(
        file.path(section_dir, "spots.tsv"),
        header = FALSE,
        stringsAsFactors = FALSE
      )[[1]],
      orig.ident = basename(section_dir)
    )
  }))
}

score_limits_from_vector <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(c(-1, 1))
  }
  limits <- as.numeric(stats::quantile(
    x,
    probs = c(PATHWAY_COLOR_LOWER_QUANTILE, PATHWAY_COLOR_UPPER_QUANTILE),
    na.rm = TRUE
  ))
  if (!all(is.finite(limits)) || limits[1] >= limits[2]) {
    limits <- range(x, na.rm = TRUE)
  }
  if (!all(is.finite(limits)) || limits[1] >= limits[2]) {
    limits <- c(min(x, na.rm = TRUE) - 1, max(x, na.rm = TRUE) + 1)
  }
  limits
}

canonicalize_pathway_name <- function(x) {
  tolower(gsub("[^A-Za-z0-9]+", "_", x))
}

build_hotspot_panel <- function(section_df, hotspot_row) {
  vmax <- as.numeric(stats::quantile(
    section_df$expr,
    probs = HOTSPOT_COLOR_UPPER_QUANTILE,
    na.rm = TRUE
  ))
  if (!is.finite(vmax) || vmax <= 0) {
    vmax <- max(section_df$expr, na.rm = TRUE)
  }
  if (!is.finite(vmax) || vmax <= 0) {
    vmax <- 1
  }

  ggplot(section_df, aes(x = imagecol, y = imagerow)) +
    geom_point(aes(color = expr), size = SPOT_POINT_SIZE, shape = 16) +
    scale_color_gradientn(
      colors = grDevices::hcl.colors(256, palette = "Viridis"),
      limits = c(0, vmax),
      oob = scales::squish
    ) +
    coord_fixed() +
    scale_y_reverse() +
    labs(
      title = paste0(
        section_df$orig.ident[[1]], "\n",
        TARGET_GENE, ": C=", sprintf("%.3f", hotspot_row$C),
        ", Z=", sprintf("%.2f", hotspot_row$Z),
        ", FDR=", format(hotspot_row$FDR, digits = 2, scientific = TRUE)
      ),
      color = paste0(TARGET_GENE, " raw counts")
    ) +
    theme_void(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
      legend.position = "right"
    )
}

build_pathway_panel <- function(section_df, pathway_label, corr_row, score_limits) {
  ggplot(section_df, aes(x = imagecol, y = imagerow)) +
    geom_point(aes(color = score_scaled), size = SPOT_POINT_SIZE, shape = 16) +
    scale_color_gradientn(
      colors = viridisLite::magma(256),
      limits = score_limits,
      oob = scales::squish
    ) +
    coord_fixed() +
    scale_y_reverse() +
    labs(
      title = paste0(
        section_df$orig.ident[[1]], "\n",
        pathway_label, ": rho=", sprintf("%.2f", corr_row$spearman_rho),
        ", FDR=", format(corr_row$fdr, digits = 2, scientific = TRUE)
      ),
      color = paste0(pathway_label, "\nscaled score")
    ) +
    theme_void(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
      legend.position = "right"
    )
}

make_summary_page <- function(pathway_label, gene_label, gene_set_text, corr_df) {
  corr_df <- corr_df %>% arrange(desc(spearman_rho))

  top_plot <- ggplot(corr_df, aes(x = reorder(orig.ident, spearman_rho), y = spearman_rho)) +
    geom_col(fill = "#7A1FA2") +
    coord_flip() +
    labs(
      title = paste0(gene_label, " vs ", pathway_label),
      subtitle = "Spearman correlation between raw counts and scaled signature score",
      x = NULL,
      y = "Spearman rho"
    ) +
    theme_minimal(base_size = 10)

  text_df <- tibble(x = 0, y = 0, label = gene_set_text)
  bottom_plot <- ggplot(text_df, aes(x = x, y = y, label = label)) +
    geom_text(hjust = 0, vjust = 1, size = 3.2, family = "mono", lineheight = 1.1) +
    xlim(0, 1) +
    ylim(0, 1) +
    labs(
      title = "Gene set",
      subtitle = "Tumor epithelial module with PanCK approximated by keratins"
    ) +
    theme_void(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, size = 10),
      plot.subtitle = element_text(hjust = 0, size = 9)
    )

  top_plot / bottom_plot + plot_layout(heights = c(2.3, 1))
}

script_dir <- get_script_dir()
project_dir <- find_project_dir(script_dir)

input_object_path <- file.path(project_dir, "inputs", "PDAC_Updated_ST.rds")
exports_dir <- file.path(project_dir, "communication", "commot", "shared_exports", "exports")
hotspot_table_path <- file.path(
  project_dir,
  "communication", "hotspot", "all_kcns_primary", "tables", "hotspot_by_gene_and_section.tsv"
)

fig_dir <- file.path(script_dir, "figures")
tab_dir <- file.path(script_dir, "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

section_dirs <- sort(list.dirs(exports_dir, recursive = FALSE, full.names = TRUE))
section_dirs <- section_dirs[grepl("IU_PDA_T", basename(section_dirs))]
if (length(section_dirs) == 0) {
  stop("No primary PDAC exported sections found in: ", exports_dir)
}

message("Loading exported primary sections...")
section_payloads <- lapply(section_dirs, load_section_payload, target_gene = TARGET_GENE)
names(section_payloads) <- basename(section_dirs)
section_map <- make_section_map(section_dirs)

hotspot_stats <- read.delim(hotspot_table_path, check.names = FALSE) %>%
  filter(gene == TARGET_GENE, orig.ident %in% names(section_payloads)) %>%
  select(orig.ident, C, Z, FDR)
write_tsv_simple(hotspot_stats, file.path(tab_dir, "kcnma1_hotspot_stats_by_section.tsv"))

message("Loading Seurat object and keeping primary PDAC spots only...")
st_obj <- readRDS(input_object_path)
st_obj <- safe_update_seurat(st_obj)

primary_spots <- intersect(colnames(st_obj), section_map$spot)
if (length(primary_spots) == 0) {
  stop("No overlap between Seurat spots and exported primary spots.")
}

if (!"RNA" %in% names(st_obj@assays)) {
  st_obj[["RNA"]] <- st_obj[["Spatial"]]
}
input_assay <- "RNA"
DefaultAssay(st_obj) <- input_assay
input_features <- rownames(get_assay_matrix_compat(st_obj, assay = input_assay, layer_name = "data"))

gene_set_members <- bind_rows(lapply(names(PATHWAY_GENE_SETS), function(pathway_name) {
  tibble(
    pathway = pathway_name,
    gene = PATHWAY_GENE_SETS[[pathway_name]],
    in_input_assay = PATHWAY_GENE_SETS[[pathway_name]] %in% input_features
  )
}))
gene_set_summary <- gene_set_members %>%
  group_by(pathway) %>%
  summarise(
    n_genes_defined = n(),
    n_genes_in_input_assay = sum(in_input_assay),
    genes_present = paste(gene[in_input_assay], collapse = ", "),
    genes_missing = paste(gene[!in_input_assay], collapse = ", "),
    .groups = "drop"
  )
write_tsv_simple(gene_set_members, file.path(tab_dir, "pathway_gene_set_members.tsv"))
write_tsv_simple(gene_set_summary, file.path(tab_dir, "pathway_gene_set_summary.tsv"))

score_cache_path <- file.path(tab_dir, "pathway_scores_by_spot.tsv")
score_long_df <- NULL
if (file.exists(score_cache_path)) {
  score_long_df <- read.delim(score_cache_path, check.names = FALSE)
  cache_ok <- score_long_df %>%
    group_by(pathway) %>%
    summarise(non_na = sum(is.finite(score_scaled)), .groups = "drop") %>%
    summarise(
      has_all_pathways = setequal(pathway, names(PATHWAY_GENE_SETS)),
      min_non_na = min(non_na),
      .groups = "drop"
    )
  if (!isTRUE(cache_ok$has_all_pathways) || cache_ok$min_non_na <= 0) {
    score_long_df <- NULL
  } else {
    message("Using cached pathway scores from pathway_scores_by_spot.tsv ...")
  }
}

if (is.null(score_long_df)) {
  message("Computing pathway activity with escape::runEscape()...")
  st_obj <- runEscape(
    input.data = st_obj,
    gene.sets = PATHWAY_GENE_SETS,
    method = "ssGSEA",
    groups = RUN_ESCAPE_GROUPS,
    min.size = 5,
    new.assay.name = "scfea",
    BPPARAM = BiocParallel::SerialParam(progressbar = FALSE)
  )

  raw_score_mat <- as.matrix(get_assay_matrix_compat(st_obj, assay = "scfea", layer_name = "data"))
  row_name_match <- match(
    canonicalize_pathway_name(rownames(raw_score_mat)),
    canonicalize_pathway_name(names(PATHWAY_GENE_SETS))
  )
  matched_names <- names(PATHWAY_GENE_SETS)[row_name_match]
  rownames(raw_score_mat) <- ifelse(is.na(matched_names), rownames(raw_score_mat), matched_names)
  scaled_score_mat <- t(scale(t(raw_score_mat)))
  scaled_score_mat[!is.finite(scaled_score_mat)] <- 0

  raw_score_df <- as.data.frame(t(raw_score_mat)) %>%
    rownames_to_column("spot") %>%
    left_join(section_map, by = "spot") %>%
    filter(!is.na(orig.ident)) %>%
    select(spot, orig.ident, everything())
  scaled_score_df <- as.data.frame(t(scaled_score_mat)) %>%
    rownames_to_column("spot") %>%
    left_join(section_map, by = "spot") %>%
    filter(!is.na(orig.ident)) %>%
    select(spot, orig.ident, everything())

  score_long_df <- bind_rows(lapply(names(PATHWAY_GENE_SETS), function(pathway_name) {
    tibble(
      spot = scaled_score_df$spot,
      orig.ident = scaled_score_df$orig.ident,
      pathway = pathway_name,
      score_raw = raw_score_df[[pathway_name]],
      score_scaled = scaled_score_df[[pathway_name]]
    )
  }))
  write_tsv_simple(score_long_df, score_cache_path)
}

message("Computing section-wise overlap between KCNMA1 and the tumor epithelial signature...")
correlation_rows <- list()
for (pathway_name in names(PATHWAY_GENE_SETS)) {
  for (section_name in names(section_payloads)) {
    section_counts <- section_payloads[[section_name]] %>%
      select(spot, expr)
    section_scores <- score_long_df %>%
      filter(pathway == pathway_name, orig.ident == section_name) %>%
      select(spot, score_scaled)

    merged_df <- section_counts %>%
      inner_join(section_scores, by = "spot") %>%
      filter(is.finite(expr), is.finite(score_scaled))

    if (nrow(merged_df) < 3 ||
        is.na(stats::sd(merged_df$expr)) ||
        is.na(stats::sd(merged_df$score_scaled)) ||
        stats::sd(merged_df$expr) == 0 ||
        stats::sd(merged_df$score_scaled) == 0) {
      rho <- NA_real_
      pval <- NA_real_
    } else {
      corr_test <- suppressWarnings(cor.test(
        merged_df$expr,
        merged_df$score_scaled,
        method = "spearman",
        exact = FALSE
      ))
      rho <- unname(corr_test$estimate)
      pval <- corr_test$p.value
    }

    correlation_rows[[length(correlation_rows) + 1]] <- tibble(
      pathway = pathway_name,
      orig.ident = section_name,
      n_spots = nrow(merged_df),
      spearman_rho = rho,
      p_value = pval
    )
  }
}

correlation_df <- bind_rows(correlation_rows) %>%
  mutate(fdr = p.adjust(p_value, method = "BH"))
correlation_summary <- correlation_df %>%
  group_by(pathway) %>%
  summarise(
    n_sections = n(),
    n_fdr_lt_0_05 = sum(fdr < 0.05, na.rm = TRUE),
    median_rho = median(spearman_rho, na.rm = TRUE),
    max_rho = max(spearman_rho, na.rm = TRUE),
    .groups = "drop"
  )
write_tsv_simple(correlation_df, file.path(tab_dir, "kcnma1_pathway_correlations_by_section.tsv"))
write_tsv_simple(correlation_summary, file.path(tab_dir, "kcnma1_pathway_correlations_summary.tsv"))

pathway_name <- names(PATHWAY_GENE_SETS)[[1]]
pathway_scores_this <- score_long_df %>% filter(pathway == pathway_name)
pathway_limits <- score_limits_from_vector(pathway_scores_this$score_scaled)
pathway_corr_this <- correlation_df %>% filter(pathway == pathway_name)
gene_set_text <- paste(PATHWAY_GENE_SETS[[pathway_name]], collapse = ", ")

pdf_path <- file.path(fig_dir, "tumor_epithelial_kcnma1_hotspot_vs_pathway_primary.pdf")
pdf(pdf_path, width = 10.5, height = 5.3, onefile = TRUE)
print(make_summary_page(pathway_name, TARGET_GENE, gene_set_text, pathway_corr_this))

for (section_name in names(section_payloads)) {
  left_df <- section_payloads[[section_name]]
  right_df <- pathway_scores_this %>%
    filter(orig.ident == section_name) %>%
    select(spot, score_scaled) %>%
    right_join(left_df, by = "spot") %>%
    mutate(orig.ident = section_name)

  hotspot_row <- hotspot_stats %>% filter(orig.ident == section_name)
  corr_row <- pathway_corr_this %>% filter(orig.ident == section_name)
  if (nrow(hotspot_row) != 1 || nrow(corr_row) != 1) {
    next
  }

  left_plot <- build_hotspot_panel(left_df, hotspot_row)
  right_plot <- build_pathway_panel(right_df, pathway_name, corr_row, pathway_limits)
  print(left_plot + right_plot + plot_layout(ncol = 2))
}
dev.off()

message("KCNMA1 tumor epithelial overlay run complete.")
