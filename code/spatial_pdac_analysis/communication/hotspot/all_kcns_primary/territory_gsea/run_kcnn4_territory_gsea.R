#!/usr/bin/env Rscript

###############################################################################
# Hallmark GSEA for the KCNN4-high territory in primary PDAC sections
#
# Earlier maps and Hotspot results showed that KCNN4 is spatially structured in
# several primary sections. Here the goal is to describe the transcriptional
# program linked to the spots where KCNN4 is highest.
#
# The analysis is section-aware. In each primary section, KCNN4-high spots are
# defined first, then compared to the rest of the same section. This avoids
# mixing spots from different slides into one artificial tissue.
#
# Visium spots are mixed, so this should be read as a territory-level result
# and not as a pure single-cell statement.
# rather than "what is the exact intrinsic program of one pure cell type?"
#
# The comparative expression step uses the SCT data layer. In this object, the
# Spatial data layer behaves more like raw counts for many genes, which is not
# ideal for a ranked GSEA. SCT gives a normalized spot-level scale that is more
# appropriate for section-wise comparisons.
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

TARGET_GENE <- "KCNN4"
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

fmt_p <- function(x) {
  ifelse(is.na(x), "NA", formatC(x, format = "e", digits = 2))
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
  # The territory is defined among positive spots only. This keeps the focus on
  # the highest KCNN4 signal instead of letting the many zero-value spots drive
  # the threshold.
  n_high <- max(min_high_spots, ceiling(length(positive_spots) * top_fraction))
  n_high <- min(n_high, length(positive_spots))
  ordered_spots <- positive_spots[order(expr_values[positive_spots], decreasing = TRUE, method = "radix")]
  ordered_spots[seq_len(n_high)]
}

plot_fgsea_bubble <- function(gsea_res, file_path) {
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
    text(0.5, 0.5, "No Hallmark pathway passed the plotting filters.")
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
      title = "Hallmark programs in the KCNN4-high territory",
      subtitle = "Ranking based on mean section-level delta (KCNN4-high minus rest)",
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

script_dir <- get_script_dir()
project_dir <- find_project_dir(script_dir)

input_dir <- file.path(project_dir, "inputs")
fig_dir <- file.path(script_dir, "figures")
tab_dir <- file.path(script_dir, "tables")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

object_path <- file.path(input_dir, "PDAC_Updated_ST.rds")
if (!file.exists(object_path)) {
  stop("Missing input object: ", object_path)
}

message("Loading spatial object...")
st_obj <- readRDS(object_path)
st_obj <- suppressWarnings(UpdateSeuratObject(st_obj))

metadata <- st_obj@meta.data %>%
  rownames_to_column("spot") %>%
  filter(Origin == "Pancreas")

# Use SCT-normalized values for the territory comparison and for the GSEA rank.
target_mat <- get_assay_matrix_compat(st_obj, assay = "SCT", layer_name = "data")

if (!(TARGET_GENE %in% rownames(target_mat))) {
  stop("Target gene not found in SCT assay: ", TARGET_GENE)
}

primary_sections <- sort(unique(metadata$orig.ident))
message("Found ", length(primary_sections), " primary sections.")

section_delta_list <- list()
section_summary_list <- list()
high_spot_list <- list()

for (section_id in primary_sections) {
  section_meta <- metadata %>%
    filter(orig.ident == section_id)

  section_spots <- section_meta$spot
  section_expr <- as.numeric(target_mat[TARGET_GENE, section_spots])
  names(section_expr) <- section_spots

  positive_spots <- names(section_expr)[section_expr > 0]
  high_spots <- character(0)

  # Require a minimum number of positive spots before defining a territory.
  # This avoids unstable high-vs-rest comparisons in very sparse sections.
  if (length(positive_spots) >= MIN_POSITIVE_SPOTS) {
    high_spots <- select_high_spots(section_expr)
  }

  rest_spots <- setdiff(section_spots, high_spots)

  section_summary_list[[section_id]] <- data.frame(
    orig.ident = section_id,
    n_spots = length(section_spots),
    n_kcnn4_positive = length(positive_spots),
    n_kcnn4_high = length(high_spots),
    min_kcnn4_in_high = if (length(high_spots) > 0) min(section_expr[high_spots]) else NA_real_,
    mean_kcnn4_in_high = if (length(high_spots) > 0) mean(section_expr[high_spots]) else NA_real_,
    mean_kcnn4_in_rest = if (length(rest_spots) > 0) mean(section_expr[rest_spots]) else NA_real_,
    used_for_gsea = length(high_spots) >= MIN_HIGH_SPOTS && length(rest_spots) > 0,
    stringsAsFactors = FALSE
  )

  if (length(high_spots) == 0 || length(rest_spots) == 0) {
    next
  }

  high_spot_list[[section_id]] <- data.frame(
    orig.ident = section_id,
    spot = high_spots,
    kcnn4_expression = section_expr[high_spots],
    stringsAsFactors = FALSE
  )

  # Section-level delta is the core quantity used later for ranking genes.
  # Positive delta means higher average expression inside the KCNN4-rich
  # territory than in the rest of the same section.
  mean_high <- Matrix::rowMeans(target_mat[, high_spots, drop = FALSE])
  mean_rest <- Matrix::rowMeans(target_mat[, rest_spots, drop = FALSE])
  section_delta_list[[section_id]] <- mean_high - mean_rest
}

section_summary_df <- bind_rows(section_summary_list)
write_tsv_simple(section_summary_df, "kcnn4_territory_section_summary.tsv", tab_dir)

if (length(high_spot_list) > 0) {
  write_tsv_simple(
    bind_rows(high_spot_list) %>% arrange(orig.ident, desc(kcnn4_expression)),
    "kcnn4_high_spots.tsv",
    tab_dir
  )
}

if (length(section_delta_list) < 3) {
  stop("Too few sections passed the KCNN4-high definition for a stable GSEA run.")
}

delta_mat <- do.call(cbind, section_delta_list)
colnames(delta_mat) <- names(section_delta_list)

# Summarize the section-level deltas for each gene. The mean delta is used as
# the ranking statistic for fgsea, while the section counts and Wilcoxon test
# give a simple sense of consistency across slides.
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
  mutate(
    wilcox_fdr = p.adjust(wilcox_p, method = "BH")
  ) %>%
  arrange(desc(mean_delta))

write_tsv_simple(gene_stats_df, "kcnn4_territory_gene_statistics.tsv", tab_dir)

ranks <- gene_stats_df$mean_delta
names(ranks) <- gene_stats_df$gene
ranks <- ranks[is.finite(ranks)]
ranks <- sort(ranks, decreasing = TRUE)

# Hallmark gene sets are used here because they give a compact first overview
# of the main biological programs linked to the territory.
hallmark_tbl <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
hallmark_pathways <- split(toupper(hallmark_tbl$gene_symbol), hallmark_tbl$gs_name)
hallmark_pathways <- lapply(hallmark_pathways, function(x) intersect(unique(x), names(ranks)))
hallmark_pathways <- hallmark_pathways[sapply(hallmark_pathways, length) >= MIN_PATHWAY_SIZE]

message("Running fgsea on ", length(hallmark_pathways), " Hallmark pathways...")
fgsea_res <- fgsea::fgseaMultilevel(
  pathways = hallmark_pathways,
  stats = ranks,
  eps = 0
) %>%
  as.data.frame() %>%
  arrange(padj, desc(abs(NES)))

fgsea_out <- fgsea_res
if ("leadingEdge" %in% colnames(fgsea_out)) {
  fgsea_out$leadingEdge <- vapply(fgsea_out$leadingEdge, function(x) paste(x, collapse = ";"), character(1))
}

write_tsv_simple(fgsea_out, "kcnn4_territory_hallmark_fgsea.tsv", tab_dir)

top_pos_genes <- gene_stats_df %>%
  arrange(desc(mean_delta)) %>%
  slice_head(n = 30)
write_tsv_simple(top_pos_genes, "kcnn4_territory_top_positive_genes.tsv", tab_dir)

top_neg_genes <- gene_stats_df %>%
  arrange(mean_delta) %>%
  slice_head(n = 30)
write_tsv_simple(top_neg_genes, "kcnn4_territory_top_negative_genes.tsv", tab_dir)

# The bubble plot keeps the output easy to read: top positive and top negative
# Hallmark programs in the KCNN4-high territory.
plot_fgsea_bubble(
  fgsea_res,
  file.path(fig_dir, "kcnn4_territory_hallmark_gsea.pdf")
)


message("KCNN4 territory GSEA complete.")
