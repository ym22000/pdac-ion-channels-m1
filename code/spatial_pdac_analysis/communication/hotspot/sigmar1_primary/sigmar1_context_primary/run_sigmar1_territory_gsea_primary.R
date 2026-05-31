#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(Matrix)
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)

TOP_FRACTION <- 0.10
MIN_POSITIVE_SPOTS <- 30
MIN_HIGH_SPOTS <- 20
MIN_PATHWAY_SIZE <- 10
TARGET_GENE <- "SIGMAR1"

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  file_path <- sub(file_arg, "", args[grep(file_arg, args)])
  if (length(file_path) == 0) return(normalizePath(getwd()))
  normalizePath(dirname(file_path))
}

find_project_dir <- function(start_dir) {
  current_dir <- normalizePath(start_dir)
  repeat {
    if (basename(current_dir) == "spatial_pdac_analysis") return(current_dir)
    parent_dir <- dirname(current_dir)
    if (identical(parent_dir, current_dir)) stop("Could not find spatial_pdac_analysis.")
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

get_assay_matrix_compat <- function(seurat_obj, assay = "SCT", layer_name = "data") {
  tryCatch(
    Seurat::GetAssayData(seurat_obj, assay = assay, layer = layer_name),
    error = function(e) Seurat::GetAssayData(seurat_obj, assay = assay, slot = layer_name)
  )
}

select_high_spots <- function(expr_values) {
  positive_spots <- names(expr_values)[expr_values > 0]
  if (length(positive_spots) < MIN_POSITIVE_SPOTS) return(character(0))
  n_high <- max(MIN_HIGH_SPOTS, ceiling(length(positive_spots) * TOP_FRACTION))
  n_high <- min(n_high, length(positive_spots))
  ordered_spots <- positive_spots[order(expr_values[positive_spots], decreasing = TRUE, method = "radix")]
  ordered_spots[seq_len(n_high)]
}

make_summary_row <- function(fgsea_res, gene_stats_df, section_summary_df) {
  pos_path <- fgsea_res %>% filter(NES > 0) %>% arrange(padj, desc(NES)) %>% slice_head(n = 1)
  neg_path <- fgsea_res %>% filter(NES < 0) %>% arrange(padj, NES) %>% slice_head(n = 1)
  pos_gene <- gene_stats_df %>% arrange(desc(mean_delta)) %>% slice_head(n = 1)
  neg_gene <- gene_stats_df %>% arrange(mean_delta) %>% slice_head(n = 1)
  data.frame(
    target_gene = TARGET_GENE,
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
    text(0.5, 0.5, paste0("No Hallmark pathway passed the plotting filters for ", TARGET_GENE, "."))
    grDevices::dev.off()
    return(invisible(NULL))
  }

  plot_df <- plot_df %>%
    mutate(
      pathway_label = gsub("^HALLMARK_", "", pathway),
      pathway_label = gsub("_", " ", pathway_label),
      neg_log10_fdr = -log10(padj + 1e-300),
      pathway_label = factor(pathway_label, levels = rev(pathway_label[order(NES)])),
      size = size
    )

  p <- ggplot(plot_df, aes(x = NES, y = pathway_label)) +
    geom_vline(xintercept = 0, linewidth = 0.4, color = "grey70") +
    geom_point(aes(size = size, color = neg_log10_fdr)) +
    scale_color_viridis_c(option = "magma", end = 0.95) +
    scale_size_continuous(range = c(3, 9)) +
    labs(
      title = paste0(TARGET_GENE, " territory Hallmark GSEA"),
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

plot_placeholder_page <- function(message_text, file_path) {
  grDevices::pdf(file_path, width = 10, height = 7)
  plot.new()
  title(main = paste0(TARGET_GENE, " territory Hallmark GSEA"))
  text(x = 0.5, y = 0.58, labels = message_text, cex = 1.1)
  grDevices::dev.off()
}

script_dir <- get_script_dir()
project_dir <- find_project_dir(script_dir)
input_dir <- file.path(project_dir, "inputs")
out_dir <- file.path(script_dir, "territory_gsea")
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

object_path <- file.path(input_dir, "PDAC_Updated_ST.rds")
if (!file.exists(object_path)) stop("Missing input object: ", object_path)

message("Loading spatial object for SIGMAR1 territory GSEA...")
st_obj <- readRDS(object_path)
st_obj <- suppressWarnings(UpdateSeuratObject(st_obj))

metadata <- st_obj@meta.data %>%
  rownames_to_column("spot") %>%
  filter(Origin == "Pancreas")

target_mat <- get_assay_matrix_compat(st_obj, assay = "SCT", layer_name = "data")
availability_df <- data.frame(
  target_gene = TARGET_GENE,
  in_sct = TARGET_GENE %in% rownames(target_mat),
  stringsAsFactors = FALSE
)
write_tsv_simple(availability_df, "sigmar1_territory_gsea_availability.tsv", tab_dir)

hallmark_tbl <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
hallmark_pathways_full <- split(toupper(hallmark_tbl$gene_symbol), hallmark_tbl$gs_name)
primary_sections <- sort(unique(metadata$orig.ident))

figure_path <- file.path(fig_dir, "sigmar1_territory_hallmark_gsea.pdf")

if (!(TARGET_GENE %in% rownames(target_mat))) {
  plot_placeholder_page("SIGMAR1 is absent from the SCT assay matrix.", figure_path)
  quit(save = "no", status = 0)
}

section_rows <- list()
gene_deltas <- list()

for (section_name in primary_sections) {
  section_spots <- metadata %>% filter(orig.ident == section_name) %>% pull(spot)
  gene_expr <- as.numeric(target_mat[TARGET_GENE, section_spots, drop = TRUE])
  names(gene_expr) <- section_spots
  high_spots <- select_high_spots(gene_expr)

  section_rows[[section_name]] <- data.frame(
    target_gene = TARGET_GENE,
    orig.ident = section_name,
    n_spots = length(section_spots),
    n_positive_spots = sum(gene_expr > 0),
    n_high_spots = length(high_spots),
    used_for_gsea = length(high_spots) >= MIN_HIGH_SPOTS,
    stringsAsFactors = FALSE
  )

  if (length(high_spots) < MIN_HIGH_SPOTS) next
  rest_spots <- setdiff(section_spots, high_spots)
  if (length(rest_spots) == 0) next

  high_mean <- Matrix::rowMeans(target_mat[, high_spots, drop = FALSE])
  rest_mean <- Matrix::rowMeans(target_mat[, rest_spots, drop = FALSE])
  delta <- as.numeric(high_mean - rest_mean)
  gene_deltas[[section_name]] <- data.frame(
    target_gene = TARGET_GENE,
    orig.ident = section_name,
    gene = rownames(target_mat),
    delta = delta,
    stringsAsFactors = FALSE
  )
}

section_summary_df <- bind_rows(section_rows)
write_tsv_simple(section_summary_df, "sigmar1_territory_section_summary.tsv", tab_dir)

if (length(gene_deltas) == 0) {
  plot_placeholder_page(
    "Too few primary sections passed the territory filters.\nNo territory GSEA could be run for SIGMAR1.",
    figure_path
  )
  summary_df <- data.frame(
    target_gene = TARGET_GENE,
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
  write_tsv_simple(summary_df, "sigmar1_territory_gsea_summary.tsv", tab_dir)
  quit(save = "no", status = 0)
}

gene_df <- bind_rows(gene_deltas)
gene_stats_df <- gene_df %>%
  group_by(target_gene, gene) %>%
  summarise(
    n_sections = n(),
    mean_delta = mean(delta, na.rm = TRUE),
    median_delta = median(delta, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_delta))
write_tsv_simple(gene_stats_df, "sigmar1_territory_gene_statistics.tsv", tab_dir)

ranks <- gene_stats_df$mean_delta
names(ranks) <- toupper(gene_stats_df$gene)
ranks <- sort(ranks, decreasing = TRUE)
present_pathways <- hallmark_pathways_full[
  vapply(hallmark_pathways_full, function(gs) sum(gs %in% names(ranks)) >= MIN_PATHWAY_SIZE, logical(1))
]
if (length(present_pathways) == 0) {
  plot_placeholder_page(
    "No Hallmark pathway reached the minimum gene overlap for SIGMAR1 territory GSEA.",
    figure_path
  )
  summary_df <- data.frame(
    target_gene = TARGET_GENE,
    n_sections_used = sum(section_summary_df$used_for_gsea, na.rm = TRUE),
    n_pathways_fdr_lt_0_05 = 0,
    n_positive_pathways_fdr_lt_0_05 = 0,
    n_negative_pathways_fdr_lt_0_05 = 0,
    top_positive_pathway = NA_character_,
    top_positive_nes = NA_real_,
    top_positive_fdr = NA_real_,
    top_negative_pathway = NA_character_,
    top_negative_nes = NA_real_,
    top_negative_fdr = NA_real_,
    top_positive_gene = gene_stats_df$gene[[1]],
    top_positive_gene_mean_delta = gene_stats_df$mean_delta[[1]],
    top_negative_gene = gene_stats_df$gene[[nrow(gene_stats_df)]],
    top_negative_gene_mean_delta = gene_stats_df$mean_delta[[nrow(gene_stats_df)]],
    stringsAsFactors = FALSE
  )
  write_tsv_simple(summary_df, "sigmar1_territory_gsea_summary.tsv", tab_dir)
  quit(save = "no", status = 0)
}
fgsea_res <- suppressWarnings(
  fgsea::fgsea(pathways = present_pathways, stats = ranks, minSize = MIN_PATHWAY_SIZE)
) %>%
  as.data.frame() %>%
  arrange(padj, desc(abs(NES))) %>%
  mutate(target_gene = TARGET_GENE)
if ("leadingEdge" %in% colnames(fgsea_res)) {
  fgsea_res$leadingEdge <- vapply(
    fgsea_res$leadingEdge,
    function(x) paste(x, collapse = ";"),
    character(1)
  )
}
write_tsv_simple(fgsea_res, "sigmar1_territory_hallmark_fgsea.tsv", tab_dir)

summary_df <- make_summary_row(
  fgsea_res = fgsea_res,
  gene_stats_df = gene_stats_df,
  section_summary_df = section_summary_df
)
write_tsv_simple(summary_df, "sigmar1_territory_gsea_summary.tsv", tab_dir)

plot_fgsea_bubble(fgsea_res, figure_path)

message("SIGMAR1 territory GSEA complete.")
