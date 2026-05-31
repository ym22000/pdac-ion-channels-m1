###############################################################################
# Moffitt stroma KCN DEG pipeline update: Test 4 with MAST
#
# This script updates the previous DEG workflow by adding a cell-level MAST test
# for subtype one-vs-rest comparisons.
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
  "VennDiagram",
  "MAST"
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
  library(VennDiagram)
  library(MAST)
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

wrap_text <- function(x, width = 44) {
  paste(strwrap(x, width = width), collapse = "\n")
}

make_placeholder_plot <- function(title_text, body_text) {
  ggplot() +
    annotate("text", x = 0, y = 0, label = paste(title_text, body_text, sep = "\n\n"), size = 5) +
    theme_void()
}

save_placeholder_if_empty <- function(df, title_text, file_path, width, height) {
  if (nrow(df) > 0) {
    return(FALSE)
  }
  save_plot(make_placeholder_plot(title_text, "No eligible KCN gene passed the current filter."), file_path, width, height)
  TRUE
}

make_dotplot <- function(seu, features, group_by, title_text, file_path, width = 12, height = 6) {
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

read_deg_table <- function(file_path) {
  read.delim(file_path, check.names = FALSE, stringsAsFactors = FALSE)
}

###############################################################################
# 2. Paths, colors, and constants
###############################################################################

project_root <- find_project_root()
analysis_root <- file.path(project_root, "M1_Bioinformatics_Ion_Channel", "code")
project_data_dir <- file.path(project_root, "M1_Bioinformatics_Ion_Channel", "data")
input_rds <- file.path(project_data_dir, "Stroma_Subset2021_fibroblast_focus.rds")

previous_output_root <- file.path(analysis_root, "moffitt_kcn_deg_3_methods")
previous_figure_dir <- file.path(previous_output_root, "figures")
previous_table_dir <- file.path(previous_output_root, "tables")

output_root <- file.path(analysis_root, "moffitt_kcn_deg_4_methods")
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
venn_colors <- c(
  "Test1_Activated_vs_Normal_DESeq2" = "#ACFFAC",
  "Test2_Subtype_OneVsRest_DESeq2" = "#FFAAAB",
  "Test3_Subtype_Hurdle_Mixed" = "#FDFEAF",
  "Test4_Subtype_MAST" = "#ADAAFF"
)

padj_cutoff <- 0.05
logfc_cutoff <- 0.25

###############################################################################
# 3. Input checks and base object
###############################################################################

if (!file.exists(input_rds)) {
  stop("Input .rds file not found: ", input_rds)
}

required_previous_files <- c(
  file.path(previous_table_dir, "DEG_Test1_Activated_vs_Normal_full.tsv"),
  file.path(previous_table_dir, "DEG_Test1_Activated_vs_Normal_KCN.tsv"),
  file.path(previous_table_dir, "DEG_Test2_Subtype_OneVsRest_full.tsv"),
  file.path(previous_table_dir, "DEG_Test2_Subtype_OneVsRest_KCN.tsv"),
  file.path(previous_table_dir, "DEG_Test3_Subtype_Hurdle_Mixed_full.tsv"),
  file.path(previous_table_dir, "DEG_Test3_Subtype_Hurdle_Mixed_KCN.tsv"),
  file.path(previous_table_dir, "Test3_hurdle_gene_panel.tsv")
)

missing_previous_files <- required_previous_files[!file.exists(required_previous_files)]
if (length(missing_previous_files) > 0) {
  stop(
    "The DEG3 output folder is incomplete. Missing files:\n",
    paste(missing_previous_files, collapse = "\n")
  )
}

stroma <- readRDS(input_rds)
DefaultAssay(stroma) <- "RNA"
stroma <- subset(stroma, subset = CellType2 %in% subtype_levels)
stroma$CellType2 <- factor(as.character(stroma$CellType2), levels = subtype_levels)
Idents(stroma) <- stroma$CellType2

meta <- stroma[[]] %>%
  rownames_to_column("cell_id") %>%
  filter(CellType2 %in% subtype_levels)

counts_mat <- GetAssayData(stroma, assay = "RNA", layer = "counts")
expr_mat <- GetAssayData(stroma, assay = "RNA", layer = "data")
kcn_genes <- grep("^KCN", rownames(stroma), value = TRUE)

###############################################################################
# 4. Bring forward the existing DEG1-3 moffitt_stromal_exploration_outputs
###############################################################################

previous_figures_to_copy <- list.files(previous_figure_dir, pattern = "\\.png$", full.names = TRUE)
previous_figures_to_copy <- previous_figures_to_copy[!grepl("Figure_02_Venn_KCN_overlap_3_methods\\.png$", previous_figures_to_copy)]
previous_figures_to_copy <- previous_figures_to_copy[!grepl("Figure_09_Hurdle_positive_expression_boxplot\\.png$", previous_figures_to_copy)]
invisible(file.copy(previous_figures_to_copy, figure_dir, overwrite = TRUE))

previous_tables_to_copy <- list.files(previous_table_dir, pattern = "\\.tsv$", full.names = TRUE)
previous_tables_to_copy <- previous_tables_to_copy[!grepl("KCN_union_membership\\.tsv$", previous_tables_to_copy)]
invisible(file.copy(previous_tables_to_copy, table_dir, overwrite = TRUE))

test1_full <- read_deg_table(file.path(previous_table_dir, "DEG_Test1_Activated_vs_Normal_full.tsv"))
test1_kcn <- read_deg_table(file.path(previous_table_dir, "DEG_Test1_Activated_vs_Normal_KCN.tsv")) %>%
  filter(grepl("^KCN", gene))
test2_full <- read_deg_table(file.path(previous_table_dir, "DEG_Test2_Subtype_OneVsRest_full.tsv"))
test2_kcn <- read_deg_table(file.path(previous_table_dir, "DEG_Test2_Subtype_OneVsRest_KCN.tsv")) %>%
  filter(grepl("^KCN", gene))
test3_full <- read_deg_table(file.path(previous_table_dir, "DEG_Test3_Subtype_Hurdle_Mixed_full.tsv"))
test3_kcn <- read_deg_table(file.path(previous_table_dir, "DEG_Test3_Subtype_Hurdle_Mixed_KCN.tsv")) %>%
  filter(grepl("^KCN", gene))

hurdle_gene_panel_tbl <- read_deg_table(file.path(previous_table_dir, "Test3_hurdle_gene_panel.tsv"))
mast_gene_panel <- hurdle_gene_panel_tbl %>%
  filter(selected_for_hurdle) %>%
  pull(gene) %>%
  unique()

save_tsv(hurdle_gene_panel_tbl, file.path(table_dir, "Test4_mast_gene_panel.tsv"))

###############################################################################
# 5. Test 4 - MAST subtype one-vs-rest
###############################################################################
# MAST rationale:
# - MAST uses a hurdle framework tailored to sparse single-cell expression
#   profiles (Finak et al., Genome Biology, 2015).
# - This is useful for KCN genes because ion channels can change through a
#   combination of detection frequency and positive-cell expression shifts.
# - Unlike the mixed hurdle model, this MAST implementation does not include a
#   patient random effect. It therefore complements, rather than replaces, the
#   patient-aware methods already used in Tests 1 to 3.
#
# Technical choice:
# - The current update runs MAST only once because Tests 1 to 3 were already
#   computed in the previous step of the project.
# - The tested feature panel matches the broadened hurdle panel used in Test 3:
#   all filtered KCN genes plus a compact set of non-KCN background genes.
#   This keeps the cell-level analysis tractable while avoiding a KCN-only
#   multiple-testing universe.
#
# Multiple-testing rule:
# - BH-adjusted p-values are computed from the raw MAST p-values with a
#   threshold of 0.05.
# - BH is preferred over Bonferroni because the goal is to control FDR while
#   remaining sensitive to moderate KCN signals.
# - Seurat also returns p_val_adj, but that value is not used here because the
#   project requires BH-adjusted p-values across the tested panel.
#
# Contrast interpretation:
# - Positive avg_log2FC means higher expression in the target subtype than in
#   the pooled rest compartment.
# - Negative avg_log2FC means relatively higher expression outside the target
#   subtype.

test4_results_list <- list()

for (current_subtype in subtype_levels) {
  markers_ct <- FindMarkers(
    object = stroma,
    ident.1 = current_subtype,
    ident.2 = NULL,
    test.use = "MAST",
    latent.vars = "nCount_RNA",
    features = mast_gene_panel,
    logfc.threshold = 0,
    min.pct = 0,
    verbose = FALSE
  )

  if (nrow(markers_ct) == 0) {
    next
  }

  markers_ct <- markers_ct %>%
    rownames_to_column("gene") %>%
    mutate(
      CellType2 = current_subtype,
      method = "Test4_Subtype_MAST",
      comparison = paste0(current_subtype, "_vs_rest"),
      padj_BH = p.adjust(p_val, method = "BH"),
      bh_significant = !is.na(padj_BH) & padj_BH < padj_cutoff,
      subtype_direction = case_when(
        avg_log2FC >= logfc_cutoff & bh_significant ~ "higher_in_target",
        avg_log2FC <= -logfc_cutoff & bh_significant ~ "lower_in_target",
        TRUE ~ "NS"
      ),
      is_kcn = gene %in% kcn_genes
    ) %>%
    arrange(padj_BH, desc(abs(avg_log2FC)))

  test4_results_list[[current_subtype]] <- markers_ct
}

test4_full <- bind_rows(test4_results_list)
test4_kcn <- test4_full %>%
  filter(is_kcn) %>%
  filter(grepl("^KCN", gene)) %>%
  select(gene, avg_log2FC, padj_BH, CellType2, method, comparison, everything())

save_tsv(test4_full, file.path(table_dir, "DEG_Test4_Subtype_MAST_full.tsv"))
save_tsv(test4_kcn, file.path(table_dir, "DEG_Test4_Subtype_MAST_KCN.tsv"))

###############################################################################
# 6. Update the 4-method KCN union table
###############################################################################
# Union rationale:
# - Different DEG frameworks capture complementary aspects of stromal biology.
# - The union therefore maximizes sensitivity and reduces method-specific bias.
# - This is particularly important for ion channels, whose effects can remain
#   moderate but consistent across complementary models.

test1_kcn_sig <- test1_kcn %>%
  filter(bh_significant)

test2_kcn_sig <- test2_kcn %>%
  filter(bh_significant, abs(log2FoldChange) >= logfc_cutoff)

test3_kcn_sig <- test3_kcn %>%
  filter(support_pattern != "none")

test4_kcn_sig <- test4_kcn %>%
  filter(bh_significant, abs(avg_log2FC) >= logfc_cutoff)

kcn_method_sets <- list(
  Test1_Activated_vs_Normal_DESeq2 = sort(unique(test1_kcn_sig$gene)),
  Test2_Subtype_OneVsRest_DESeq2 = sort(unique(test2_kcn_sig$gene)),
  Test3_Subtype_Hurdle_Mixed = sort(unique(test3_kcn_sig$gene)),
  Test4_Subtype_MAST = sort(unique(test4_kcn_sig$gene))
)

union_kcn <- sort(unique(unlist(kcn_method_sets)))

kcn_union_table <- tibble(gene = union_kcn) %>%
  mutate(
    Test1_Activated_vs_Normal_DESeq2 = gene %in% kcn_method_sets$Test1_Activated_vs_Normal_DESeq2,
    Test2_Subtype_OneVsRest_DESeq2 = gene %in% kcn_method_sets$Test2_Subtype_OneVsRest_DESeq2,
    Test3_Subtype_Hurdle_Mixed = gene %in% kcn_method_sets$Test3_Subtype_Hurdle_Mixed,
    Test4_Subtype_MAST = gene %in% kcn_method_sets$Test4_Subtype_MAST
  ) %>%
  rowwise() %>%
  mutate(
    n_methods = sum(c_across(Test1_Activated_vs_Normal_DESeq2:Test4_Subtype_MAST)),
    methods_present = paste(
      c(
        "Test1_Activated_vs_Normal_DESeq2",
        "Test2_Subtype_OneVsRest_DESeq2",
        "Test3_Subtype_Hurdle_Mixed",
        "Test4_Subtype_MAST"
      )[c_across(Test1_Activated_vs_Normal_DESeq2:Test4_Subtype_MAST)],
      collapse = "; "
    )
  ) %>%
  ungroup()

save_tsv(kcn_union_table, file.path(table_dir, "KCN_union_membership.tsv"))

###############################################################################
# 7. Updated 4-set Venn diagram
###############################################################################

venn_grob <- venn.diagram(
  x = kcn_method_sets,
  filename = NULL,
  fill = unname(venn_colors),
  alpha = 0.75,
  col = "grey25",
  lwd = 1.1,
  cex = 1.2,
  fontface = "bold",
  cat.cex = 0,
  category.names = rep("", 4),
  margin = 0.08
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
  geom_text(aes(label = label), hjust = 0, nudge_x = 0.18, size = 3.1, lineheight = 0.95) +
  scale_color_manual(values = venn_colors) +
  coord_cartesian(xlim = c(0, 5.9), ylim = c(0.5, nrow(venn_legend_df) + 0.5), clip = "off") +
  theme_void() +
  theme(legend.position = "none")

p_venn <- wrap_elements(full = venn_grob) + legend_plot + plot_layout(widths = c(1.15, 1))

save_plot(
  p_venn,
  file.path(figure_dir, "Figure_02_Venn_KCN_overlap_4_methods.png"),
  width = 16,
  height = 8
)

###############################################################################
# 8. MAST KCN DotPlot
###############################################################################

test4_dot_genes <- test4_kcn %>%
  filter(bh_significant, avg_log2FC > 0) %>%
  group_by(CellType2) %>%
  arrange(padj_BH, desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 4) %>%
  ungroup() %>%
  pull(gene) %>%
  unique()

if (length(test4_dot_genes) == 0) {
  test4_dot_genes <- head(filtered_kcn <- sort(unique(test4_kcn$gene)), 15)
}

make_dotplot(
  seu = stroma,
  features = test4_dot_genes,
  group_by = "CellType2",
  title_text = "KCN DotPlot - MAST subtype one-vs-rest",
  file_path = file.path(figure_dir, "Figure_10_DotPlot_KCN_Test4_MAST.png"),
  width = 14,
  height = 6
)

###############################################################################
# 9. Save bundle and logs
###############################################################################

saveRDS(
  list(
    Test1 = test1_full,
    Test2 = test2_full,
    Test3 = test3_full,
    Test4 = test4_full,
    KCN_union = kcn_union_table
  ),
  file = file.path(rds_dir, "deg4_results_bundle.rds")
)

writeLines(
  c(
    "Moffitt stromal KCN DEG pipeline update (4 methods)",
    paste("Run date:", as.character(Sys.Date())),
    paste("Input object:", input_rds),
    paste("Current union KCN genes:", if (length(union_kcn) > 0) paste(union_kcn, collapse = ", ") else "None")
  ),
  con = file.path(output_root, "README_pipeline_summary.txt")
)

message("MAST update finished.")
message("Outputs saved in: ", output_root)

