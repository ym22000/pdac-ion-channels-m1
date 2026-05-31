#!/usr/bin/env Rscript

# SPARK-X analysis for the KCN panel in primary PDAC sections
#
# This script was added to complement the Hotspot analysis with another
# published method for spatially variable genes. Hotspot was already useful to
# show that some KCN genes form spatial patches, but it is still helpful to
# check the same idea with a second framework that is widely used in spatial
# transcriptomics.
#
# Here the goal is simple:
#   1. run SPARK-X section by section on the primary Visium exports
#   2. keep the results for the 25-gene KCN union panel
#   3. summarize which KCN genes are repeatedly significant across sections
#   4. save one PDF per KCN with a section summary and viridis expression maps
#
# We use SPARK-X instead of the older SPARK model because these Visium sections
# have more than 3,000 spots, and the official SPARK documentation recommends
# SPARK-X for larger datasets. This also keeps the code simple and scalable.
#
# Method resources
# - SPARK package:
#   https://github.com/xzhoulab/SPARK
# - Official website:
#   https://xzhoulab.github.io/SPARK/
# - SPARK-X paper:
#   Zhu, Sun and Zhou, Genome Biology 2021

suppressPackageStartupMessages({
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(tibble)
  library(SPARK)
})

find_project_dir <- function(start_dir = getwd()) {
  current <- normalizePath(start_dir, winslash = "/", mustWork = TRUE)

  repeat {
    if (basename(current) == "spatial_pdac_analysis") {
      return(current)
    }

    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not find spatial_pdac_analysis project root.")
    }
    current <- parent
  }
}

project_dir <- find_project_dir()
exports_dir <- file.path(project_dir, "communication", "commot", "shared_exports", "exports")
kcn_list_path <- file.path(dirname(project_dir), "mycaf_icaf_correlations", "kcn_union_25.tsv")
out_dir <- file.path(project_dir, "sparkx")
tables_dir <- file.path(out_dir, "tables")
figures_dir <- file.path(out_dir, "figures")
logs_dir <- file.path(out_dir, "logs")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

options(stringsAsFactors = FALSE)
theme_set(theme_minimal(base_size = 11))

num_cores <- 1
top_n_gene_labels <- 12
color_upper_quantile <- 0.95

load_section <- function(section_dir) {
  counts <- readMM(file.path(section_dir, "counts.mtx"))
  genes <- read_tsv(file.path(section_dir, "genes.tsv"), col_names = FALSE, show_col_types = FALSE)[[1]]
  spots <- read_tsv(file.path(section_dir, "spots.tsv"), col_names = FALSE, show_col_types = FALSE)[[1]]
  coords <- read_tsv(file.path(section_dir, "coords.tsv"), show_col_types = FALSE)
  metadata <- read_tsv(file.path(section_dir, "metadata.tsv"), show_col_types = FALSE)

  rownames(counts) <- genes
  colnames(counts) <- spots

  coords <- coords[match(spots, coords$spot), , drop = FALSE]
  metadata <- metadata[match(spots, metadata$spot), , drop = FALSE]

  list(
    counts = counts,
    genes = genes,
    spots = spots,
    coords = coords,
    metadata = metadata
  )
}

run_section_sparkx <- function(section_name, section_dir) {
  message("Running SPARK-X on ", section_name)
  section <- load_section(section_dir)
  positive_spot_vec <- as.numeric(Matrix::rowSums(section$counts > 0))
  names(positive_spot_vec) <- rownames(section$counts)

  loc <- as.matrix(section$coords[, c("imagecol", "imagerow")])
  sparkx_res <- sparkx(
    count_in = section$counts,
    locus_in = loc,
    numCores = num_cores,
    option = "mixture",
    verbose = FALSE
  )

  stats_tbl <- as.data.frame(sparkx_res$stats) %>%
    rownames_to_column("gene")

  single_kernel_p_tbl <- as.data.frame(sparkx_res$res_stest) %>%
    rownames_to_column("gene")

  multi_test_tbl <- as.data.frame(sparkx_res$res_mtest) %>%
    rownames_to_column("gene")

  merged <- stats_tbl %>%
    inner_join(single_kernel_p_tbl, by = "gene", suffix = c("_stat", "_pval")) %>%
    inner_join(multi_test_tbl, by = "gene") %>%
    mutate(
      orig.ident = section_name,
      n_spots = ncol(section$counts),
      n_genes_tested = nrow(section$counts),
      positive_spots = positive_spot_vec[gene],
      best_kernel = c(
        "projection", "gaus1", "gaus2", "gaus3", "gaus4", "gaus5",
        "cos1", "cos2", "cos3", "cos4", "cos5"
      )[max.col(as.matrix(select(., ends_with("_stat"))), ties.method = "first")],
      best_kernel_stat = pmax(
        projection_stat, gaus1_stat, gaus2_stat, gaus3_stat, gaus4_stat, gaus5_stat,
        cos1_stat, cos2_stat, cos3_stat, cos4_stat, cos5_stat
      )
    ) %>%
    mutate(
      best_kernel_pval = case_when(
        best_kernel == "projection" ~ projection_pval,
        best_kernel == "gaus1" ~ gaus1_pval,
        best_kernel == "gaus2" ~ gaus2_pval,
        best_kernel == "gaus3" ~ gaus3_pval,
        best_kernel == "gaus4" ~ gaus4_pval,
        best_kernel == "gaus5" ~ gaus5_pval,
        best_kernel == "cos1" ~ cos1_pval,
        best_kernel == "cos2" ~ cos2_pval,
        best_kernel == "cos3" ~ cos3_pval,
        best_kernel == "cos4" ~ cos4_pval,
        best_kernel == "cos5" ~ cos5_pval,
        TRUE ~ NA_real_
      )
    )

  list(section = section, results = merged)
}

plot_gene_pdf <- function(gene, gene_results, section_cache, pdf_path) {
  pdf(pdf_path, width = 11, height = 8.5, onefile = TRUE)
  on.exit(dev.off(), add = TRUE)

  summary_tbl <- gene_results %>%
    mutate(neg_log10_fdr = -log10(pmax(adjustedPval, 1e-300))) %>%
    arrange(desc(neg_log10_fdr))

  top_label_tbl <- summary_tbl %>%
    slice_head(n = min(top_n_gene_labels, nrow(summary_tbl)))

  p_summary <- ggplot(summary_tbl, aes(x = neg_log10_fdr, y = reorder(orig.ident, neg_log10_fdr))) +
    geom_col(fill = "#1f77b4") +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "firebrick") +
    labs(
      title = paste(gene, "SPARK-X significance across primary sections"),
      subtitle = "Bars show -log10 adjusted p-value. Dashed line marks FDR = 0.05.",
      x = expression(-log[10]("adjusted p-value")),
      y = "Section"
    ) +
    theme(panel.grid.major.y = element_blank())
  print(p_summary)

  for (section_name in gene_results$orig.ident) {
    section_obj <- section_cache[[section_name]]
    section_df <- tibble(
      spot = section_obj$spots,
      imagecol = section_obj$coords$imagecol,
      imagerow = section_obj$coords$imagerow,
      expr = as.numeric(section_obj$counts[gene, ])
    )

    stat_row <- gene_results %>% filter(orig.ident == section_name)
    upper_limit <- stats::quantile(section_df$expr, probs = color_upper_quantile, na.rm = TRUE)
    if (!is.finite(upper_limit) || upper_limit <= 0) {
      upper_limit <- max(section_df$expr, na.rm = TRUE)
    }
    if (!is.finite(upper_limit) || upper_limit <= 0) {
      upper_limit <- 1
    }

    p_map <- ggplot(section_df, aes(x = imagecol, y = -imagerow, color = expr)) +
      geom_point(size = 0.65) +
      scale_color_viridis_c(limits = c(0, upper_limit), oob = scales::squish) +
      coord_equal() +
      labs(
        title = paste(gene, "expression in", section_name),
        subtitle = paste0(
          "positive spots = ", stat_row$positive_spots,
          " | best kernel = ", stat_row$best_kernel,
          " | best stat = ", sprintf("%.2f", stat_row$best_kernel_stat),
          " | adjusted p = ", format(stat_row$adjustedPval, scientific = TRUE, digits = 3)
        ),
        x = "imagecol",
        y = "imagerow",
        color = "Counts"
      ) +
      theme_void() +
      theme(
        plot.title = element_text(face = "bold"),
        legend.position = "right"
      )
    print(p_map)
  }
}

kcn_tbl <- read_tsv(kcn_list_path, show_col_types = FALSE)
kcn_genes <- unique(kcn_tbl$gene)

section_dirs <- list.dirs(exports_dir, recursive = FALSE, full.names = TRUE)
section_names <- basename(section_dirs)

section_results <- vector("list", length(section_dirs))
names(section_results) <- section_names
section_cache <- vector("list", length(section_dirs))
names(section_cache) <- section_names

for (i in seq_along(section_dirs)) {
  out <- run_section_sparkx(section_names[i], section_dirs[i])
  section_results[[i]] <- out$results
  section_cache[[i]] <- out$section
}

all_results <- bind_rows(section_results)

availability_tbl <- expand_grid(
  orig.ident = section_names,
  gene = kcn_genes
) %>%
  left_join(
    all_results %>% select(orig.ident, gene, positive_spots, combinedPval, adjustedPval, best_kernel, best_kernel_stat, best_kernel_pval),
    by = c("orig.ident", "gene")
  ) %>%
  mutate(
    available_in_section = !is.na(combinedPval),
    positive_spots = if_else(is.na(positive_spots), 0, positive_spots)
  )

kcn_results <- availability_tbl %>%
  filter(available_in_section) %>%
  left_join(
    all_results %>% select(orig.ident, gene, n_spots, n_genes_tested),
    by = c("orig.ident", "gene")
  ) %>%
  arrange(gene, orig.ident)

summary_tbl <- kcn_results %>%
  group_by(gene) %>%
  summarise(
    n_sections_tested = n(),
    n_sections_fdr_lt_0_05 = sum(adjustedPval < 0.05, na.rm = TRUE),
    median_adjusted_pval = median(adjustedPval, na.rm = TRUE),
    min_adjusted_pval = min(adjustedPval, na.rm = TRUE),
    median_best_kernel_stat = median(best_kernel_stat, na.rm = TRUE),
    max_best_kernel_stat = max(best_kernel_stat, na.rm = TRUE),
    top_section = orig.ident[which.min(adjustedPval)],
    top_section_fdr = min(adjustedPval, na.rm = TRUE),
    top_section_best_kernel = best_kernel[which.min(adjustedPval)],
    .groups = "drop"
  ) %>%
  arrange(desc(n_sections_fdr_lt_0_05), median_adjusted_pval)

availability_summary_tbl <- availability_tbl %>%
  group_by(gene) %>%
  summarise(
    n_sections_available = sum(available_in_section),
    n_sections_missing = sum(!available_in_section),
    .groups = "drop"
  ) %>%
  arrange(desc(n_sections_available), gene)

write_tsv(kcn_results, file.path(tables_dir, "sparkx_by_gene_and_section.tsv"))
write_tsv(summary_tbl, file.path(tables_dir, "sparkx_gene_summary.tsv"))
write_tsv(availability_tbl, file.path(tables_dir, "sparkx_gene_availability_by_section.tsv"))
write_tsv(availability_summary_tbl, file.path(tables_dir, "sparkx_gene_availability_summary.tsv"))

for (gene in unique(kcn_results$gene)) {
  gene_results <- kcn_results %>% filter(gene == !!gene)
  pdf_path <- file.path(figures_dir, paste0(tolower(gene), "_sparkx_primary_sections.pdf"))
  plot_gene_pdf(gene, gene_results, section_cache, pdf_path)
}

message("SPARK-X analysis finished.")
