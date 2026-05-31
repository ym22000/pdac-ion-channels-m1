#!/usr/bin/env Rscript

###############################################################################
# Hallmark GSEA for KCN-high territories in normal pancreas sections
#
# For each KCN, this script defines a KCN-high territory inside each normal
# pancreas section, compares it to the rest of the same section, and then runs
# Hallmark GSEA on the average section-level expression difference.
#
# The analysis is kept section by section on purpose. Spots from different
# Visium slides should not be mixed into one fake tissue.
#
# The SCT assay is used here because it gives a cleaner normalized scale for
# ranking genes than the Spatial data layer in this object. The goal is to see
# which KCN-linked territory programs already exist in healthy tissue before
# comparing them with PDAC.
###############################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(Matrix)
  library(ggplot2)
  library(fgsea)
  library(msigdbr)
})

options(stringsAsFactors = FALSE)

TOP_FRACTION <- 0.10
MIN_POSITIVE_SPOTS <- 30
MIN_HIGH_SPOTS <- 20
MIN_PATHWAY_SIZE <- 10

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

select_high_spots <- function(expr_values, top_fraction = TOP_FRACTION, min_high_spots = MIN_HIGH_SPOTS) {
  positive_spots <- names(expr_values)[expr_values > 0]
  if (length(positive_spots) == 0) {
    return(character(0))
  }
  n_high <- max(min_high_spots, ceiling(length(positive_spots) * top_fraction))
  n_high <- min(n_high, length(positive_spots))
  ordered_spots <- positive_spots[order(expr_values[positive_spots], decreasing = TRUE, method = "radix")]
  ordered_spots[seq_len(n_high)]
}

plot_fgsea_bubble <- function(gsea_res, target_gene, file_path) {
  filtered <- gsea_res %>%
    filter(!is.na(NES), !is.na(padj)) %>%
    arrange(padj, desc(abs(NES)))

  positive_hits <- filtered %>% filter(NES > 0) %>% slice_head(n = 8)
  negative_hits <- filtered %>% filter(NES < 0) %>% slice_head(n = 8)
  plot_df <- bind_rows(positive_hits, negative_hits) %>%
    distinct(pathway, .keep_all = TRUE)

  if (nrow(plot_df) == 0) {
    grDevices::pdf(file_path, width = 9, height = 6)
    plot.new()
    text(0.5, 0.5, paste0("No Hallmark pathway passed the plotting filters for ", target_gene, "."))
    grDevices::dev.off()
    return(invisible(NULL))
  }

  plot_df <- plot_df %>%
    mutate(
      pathway_label = gsub("^HALLMARK_", "", pathway),
      pathway_label = gsub("_", " ", pathway_label),
      neg_log10_fdr = -log10(padj + 1e-300),
      pathway_label = factor(pathway_label, levels = rev(pathway_label[order(NES)]))
    )

  p <- ggplot(plot_df, aes(x = NES, y = pathway_label)) +
    geom_vline(xintercept = 0, linewidth = 0.4, color = "grey70") +
    geom_point(aes(size = size, color = neg_log10_fdr)) +
    scale_color_viridis_c(option = "magma", end = 0.95) +
    scale_size_continuous(range = c(3, 9)) +
    labs(
      title = paste0(target_gene, " normal-pancreas territory Hallmark GSEA"),
      subtitle = "Ranking based on mean section-level delta (high territory minus rest)",
      x = "NES",
      y = NULL,
      color = "-log10(FDR)",
      size = "Pathway size"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold")
    )

  ggplot2::ggsave(file_path, p, width = 10, height = 7)
}

plot_placeholder_page <- function(target_gene, message_text, file_path) {
  grDevices::pdf(file_path, width = 10, height = 7)
  plot.new()
  title(main = paste0(target_gene, " normal-pancreas territory Hallmark GSEA"))
  text(
    x = 0.5,
    y = 0.58,
    labels = message_text,
    cex = 1.1
  )
  grDevices::dev.off()
}

make_gene_summary_row <- function(target_gene, fgsea_res, gene_stats_df, section_summary_df) {
  pos_path <- fgsea_res %>%
    filter(NES > 0) %>%
    arrange(padj, desc(NES)) %>%
    slice_head(n = 1)

  neg_path <- fgsea_res %>%
    filter(NES < 0) %>%
    arrange(padj, NES) %>%
    slice_head(n = 1)

  pos_gene <- gene_stats_df %>%
    arrange(desc(mean_delta)) %>%
    slice_head(n = 1)

  neg_gene <- gene_stats_df %>%
    arrange(mean_delta) %>%
    slice_head(n = 1)

  data.frame(
    target_gene = target_gene,
    n_sections_used = sum(section_summary_df$used_for_gsea, na.rm = TRUE),
    n_pathways_fdr_lt_0_05 = sum(fgsea_res$padj < 0.05, na.rm = TRUE),
    n_positive_pathways_fdr_lt_0_05 = sum(fgsea_res$padj < 0.05 & fgsea_res$NES > 0, na.rm = TRUE),
    n_negative_pathways_fdr_lt_0_05 = sum(fgsea_res$padj < 0.05 & fgsea_res$NES < 0, na.rm = TRUE),
    top_positive_pathway = if (nrow(pos_path) > 0) pos_path$pathway else NA_character_,
    top_positive_nes = if (nrow(pos_path) > 0) pos_path$NES else NA_real_,
    top_positive_fdr = if (nrow(pos_path) > 0) pos_path$padj else NA_real_,
    top_negative_pathway = if (nrow(neg_path) > 0) neg_path$pathway else NA_character_,
    top_negative_nes = if (nrow(neg_path) > 0) neg_path$NES else NA_real_,
    top_negative_fdr = if (nrow(neg_path) > 0) neg_path$padj else NA_real_,
    top_positive_gene = if (nrow(pos_gene) > 0) pos_gene$gene else NA_character_,
    top_positive_gene_mean_delta = if (nrow(pos_gene) > 0) pos_gene$mean_delta else NA_real_,
    top_negative_gene = if (nrow(neg_gene) > 0) neg_gene$gene else NA_character_,
    top_negative_gene_mean_delta = if (nrow(neg_gene) > 0) neg_gene$mean_delta else NA_real_,
    stringsAsFactors = FALSE
  )
}

script_dir <- get_script_dir()
project_dir <- find_project_dir(script_dir)

input_dir <- file.path(project_dir, "inputs")
fig_dir <- file.path(script_dir, "figures")
tab_dir <- file.path(script_dir, "tables")
kcn_list_path <- normalizePath(file.path(project_dir, "..", "mycaf_icaf_correlations", "kcn_union_25.tsv"), mustWork = FALSE)

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

object_path <- file.path(input_dir, "PDAC_Updated_ST.rds")
if (!file.exists(object_path)) {
  stop("Missing input object: ", object_path)
}
if (!file.exists(kcn_list_path)) {
  stop("Missing KCN list: ", kcn_list_path)
}

message("Loading spatial object...")
st_obj <- readRDS(object_path)
st_obj <- suppressWarnings(UpdateSeuratObject(st_obj))

metadata <- st_obj@meta.data %>%
  rownames_to_column("spot") %>%
  filter(Origin == "Normal Pancreas")

target_mat <- get_assay_matrix_compat(st_obj, assay = "SCT", layer_name = "data")
kcn_genes <- read.delim(kcn_list_path, stringsAsFactors = FALSE) %>%
  pull(gene) %>%
  unique()

availability_df <- data.frame(
  target_gene = kcn_genes,
  in_sct_assay = kcn_genes %in% rownames(target_mat),
  stringsAsFactors = FALSE
)
write_tsv_simple(availability_df, "kcn_territory_gsea_gene_availability.tsv", tab_dir)

hallmark_tbl <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
hallmark_pathways_full <- split(toupper(hallmark_tbl$gene_symbol), hallmark_tbl$gs_name)

normal_sections <- sort(unique(metadata$orig.ident))
message("Found ", length(normal_sections), " normal pancreas sections.")

all_section_summaries <- list()
all_high_spots <- list()
all_gene_stats <- list()
all_fgsea <- list()
all_top_positive <- list()
all_top_negative <- list()
all_gene_summaries <- list()

for (target_gene in kcn_genes) {
  message("Processing ", target_gene, " ...")
  figure_path <- file.path(fig_dir, paste0(tolower(target_gene), "_normal_pancreas_territory_hallmark_gsea.pdf"))

  if (!(target_gene %in% rownames(target_mat))) {
    plot_placeholder_page(
      target_gene = target_gene,
      message_text = "Target gene not available in the SCT assay.\nNo territory GSEA was run for this gene.",
      file_path = figure_path
    )
    all_gene_summaries[[target_gene]] <- data.frame(
      target_gene = target_gene,
      n_sections_used = 0,
      n_pathways_fdr_lt_0_05 = NA_integer_,
      n_positive_pathways_fdr_lt_0_05 = NA_integer_,
      n_negative_pathways_fdr_lt_0_05 = NA_integer_,
      top_positive_pathway = NA_character_,
      top_positive_nes = NA_real_,
      top_positive_fdr = NA_real_,
      top_negative_pathway = NA_character_,
      top_negative_nes = NA_real_,
      top_negative_fdr = NA_real_,
      top_positive_gene = NA_character_,
      top_positive_gene_mean_delta = NA_real_,
      top_negative_gene = NA_character_,
      top_negative_gene_mean_delta = NA_real_,
      stringsAsFactors = FALSE
    )
    next
  }

  section_delta_list <- list()
  section_summary_list <- list()
  high_spot_list <- list()

  for (section_id in normal_sections) {
    section_meta <- metadata %>%
      filter(orig.ident == section_id)

    section_spots <- section_meta$spot
    section_expr <- as.numeric(target_mat[target_gene, section_spots])
    names(section_expr) <- section_spots

    positive_spots <- names(section_expr)[section_expr > 0]
    high_spots <- character(0)

    if (length(positive_spots) >= MIN_POSITIVE_SPOTS) {
      high_spots <- select_high_spots(section_expr)
    }

    rest_spots <- setdiff(section_spots, high_spots)

    section_summary_list[[section_id]] <- data.frame(
      target_gene = target_gene,
      orig.ident = section_id,
      n_spots = length(section_spots),
      n_positive = length(positive_spots),
      n_high = length(high_spots),
      min_expression_in_high = if (length(high_spots) > 0) min(section_expr[high_spots]) else NA_real_,
      mean_expression_in_high = if (length(high_spots) > 0) mean(section_expr[high_spots]) else NA_real_,
      mean_expression_in_rest = if (length(rest_spots) > 0) mean(section_expr[rest_spots]) else NA_real_,
      used_for_gsea = length(high_spots) >= MIN_HIGH_SPOTS && length(rest_spots) > 0,
      stringsAsFactors = FALSE
    )

    if (length(high_spots) == 0 || length(rest_spots) == 0) {
      next
    }

    high_spot_list[[section_id]] <- data.frame(
      target_gene = target_gene,
      orig.ident = section_id,
      spot = high_spots,
      target_expression = section_expr[high_spots],
      stringsAsFactors = FALSE
    )

    mean_high <- Matrix::rowMeans(target_mat[, high_spots, drop = FALSE])
    mean_rest <- Matrix::rowMeans(target_mat[, rest_spots, drop = FALSE])
    section_delta_list[[section_id]] <- mean_high - mean_rest
  }

  section_summary_df <- bind_rows(section_summary_list)
  all_section_summaries[[target_gene]] <- section_summary_df

  if (length(high_spot_list) > 0) {
    all_high_spots[[target_gene]] <- bind_rows(high_spot_list) %>%
      arrange(orig.ident, desc(target_expression))
  }

  if (length(section_delta_list) < 3) {
    plot_placeholder_page(
      target_gene = target_gene,
      message_text = paste0(
        "Too few normal pancreas sections passed the territory filters.\n",
        "A stable territory GSEA needs at least 3 usable sections.\n",
        "Here, n_sections_used = ", sum(section_summary_df$used_for_gsea, na.rm = TRUE), "."
      ),
      file_path = figure_path
    )
    all_gene_summaries[[target_gene]] <- data.frame(
      target_gene = target_gene,
      n_sections_used = sum(section_summary_df$used_for_gsea, na.rm = TRUE),
      n_pathways_fdr_lt_0_05 = NA_integer_,
      n_positive_pathways_fdr_lt_0_05 = NA_integer_,
      n_negative_pathways_fdr_lt_0_05 = NA_integer_,
      top_positive_pathway = NA_character_,
      top_positive_nes = NA_real_,
      top_positive_fdr = NA_real_,
      top_negative_pathway = NA_character_,
      top_negative_nes = NA_real_,
      top_negative_fdr = NA_real_,
      top_positive_gene = NA_character_,
      top_positive_gene_mean_delta = NA_real_,
      top_negative_gene = NA_character_,
      top_negative_gene_mean_delta = NA_real_,
      stringsAsFactors = FALSE
    )
    next
  }

  delta_mat <- do.call(cbind, section_delta_list)
  colnames(delta_mat) <- names(section_delta_list)

  gene_stats <- lapply(rownames(delta_mat), function(gene_id) {
    values <- as.numeric(delta_mat[gene_id, ])
    values <- values[is.finite(values)]
    n_tested <- length(values)

    mean_delta <- if (n_tested > 0) mean(values) else NA_real_
    median_delta <- if (n_tested > 0) median(values) else NA_real_
    n_sections_delta_gt_0 <- if (n_tested > 0) sum(values > 0) else NA_integer_
    n_sections_delta_lt_0 <- if (n_tested > 0) sum(values < 0) else NA_integer_

    wilcox_p <- NA_real_
    if (n_tested >= 3 && length(unique(values)) > 1) {
      wilcox_p <- suppressWarnings(wilcox.test(values, mu = 0, alternative = "two.sided", exact = FALSE)$p.value)
    }

    data.frame(
      target_gene = target_gene,
      gene = gene_id,
      n_sections_tested = n_tested,
      n_sections_delta_gt_0 = n_sections_delta_gt_0,
      n_sections_delta_lt_0 = n_sections_delta_lt_0,
      mean_delta = mean_delta,
      median_delta = median_delta,
      wilcox_p = wilcox_p,
      stringsAsFactors = FALSE
    )
  })

  gene_stats_df <- bind_rows(gene_stats) %>%
    mutate(wilcox_fdr = p.adjust(wilcox_p, method = "BH")) %>%
    arrange(desc(mean_delta))
  all_gene_stats[[target_gene]] <- gene_stats_df

  ranks <- gene_stats_df$mean_delta
  names(ranks) <- gene_stats_df$gene
  ranks <- ranks[is.finite(ranks)]
  ranks <- sort(ranks, decreasing = TRUE)

  hallmark_pathways <- lapply(hallmark_pathways_full, function(x) intersect(unique(x), names(ranks)))
  hallmark_pathways <- hallmark_pathways[sapply(hallmark_pathways, length) >= MIN_PATHWAY_SIZE]

  fgsea_res <- fgsea::fgseaMultilevel(
    pathways = hallmark_pathways,
    stats = ranks,
    eps = 0
  ) %>%
    as.data.frame() %>%
    arrange(padj, desc(abs(NES))) %>%
    mutate(target_gene = target_gene, .before = 1)

  fgsea_out <- fgsea_res
  if ("leadingEdge" %in% colnames(fgsea_out)) {
    fgsea_out$leadingEdge <- vapply(fgsea_out$leadingEdge, function(x) paste(x, collapse = ";"), character(1))
  }
  all_fgsea[[target_gene]] <- fgsea_out

  top_pos_genes <- gene_stats_df %>%
    arrange(desc(mean_delta)) %>%
    slice_head(n = 30)
  all_top_positive[[target_gene]] <- top_pos_genes

  top_neg_genes <- gene_stats_df %>%
    arrange(mean_delta) %>%
    slice_head(n = 30)
  all_top_negative[[target_gene]] <- top_neg_genes

  all_gene_summaries[[target_gene]] <- make_gene_summary_row(
    target_gene = target_gene,
    fgsea_res = fgsea_res,
    gene_stats_df = gene_stats_df,
    section_summary_df = section_summary_df
  )

  plot_fgsea_bubble(
    gsea_res = fgsea_res,
    target_gene = target_gene,
    file_path = figure_path
  )
}

if (length(all_section_summaries) > 0) {
  write_tsv_simple(bind_rows(all_section_summaries), "kcn_territory_section_summary.tsv", tab_dir)
}

if (length(all_high_spots) > 0) {
  write_tsv_simple(bind_rows(all_high_spots), "kcn_high_spots.tsv", tab_dir)
}

if (length(all_gene_stats) > 0) {
  write_tsv_simple(bind_rows(all_gene_stats), "kcn_territory_gene_statistics.tsv", tab_dir)
}

if (length(all_fgsea) > 0) {
  write_tsv_simple(bind_rows(all_fgsea), "kcn_territory_hallmark_fgsea.tsv", tab_dir)
}

if (length(all_top_positive) > 0) {
  write_tsv_simple(bind_rows(all_top_positive), "kcn_territory_top_positive_genes.tsv", tab_dir)
}

if (length(all_top_negative) > 0) {
  write_tsv_simple(bind_rows(all_top_negative), "kcn_territory_top_negative_genes.tsv", tab_dir)
}

if (length(all_gene_summaries) > 0) {
  write_tsv_simple(bind_rows(all_gene_summaries), "kcn_territory_gsea_summary.tsv", tab_dir)
}


message("All KCN territory GSEA runs complete.")
