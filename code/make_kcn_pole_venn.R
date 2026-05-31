###############################################################################
# KCN pole Venn diagram
#
# This script builds the Venn diagram used to summarize the normal-like and
# activated-like KCN poles from the different DEG and support tests.
###############################################################################
required_packages <- c("dplyr", "ggplot2", "patchwork", "tibble", "VennDiagram")

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
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(tibble)
  library(VennDiagram)
})

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

wrap_text <- function(x, width = 44) {
  paste(strwrap(x, width = width), collapse = "\n")
}

project_root <- find_project_root()
output_root <- file.path(project_root, "M1_Bioinformatics_Ion_Channel", "code", "moffitt_kcn_deg_4_methods")
figure_dir <- file.path(output_root, "figures")
table_dir <- file.path(output_root, "tables")

t1 <- read.delim(file.path(table_dir, "DEG_Test1_Activated_vs_Normal_KCN.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
t2 <- read.delim(file.path(table_dir, "DEG_Test2_Subtype_OneVsRest_KCN.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
t3 <- read.delim(file.path(table_dir, "DEG_Test3_Subtype_Hurdle_Mixed_KCN.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
t4 <- read.delim(file.path(table_dir, "DEG_Test4_Subtype_MAST_KCN.tsv"), check.names = FALSE, stringsAsFactors = FALSE)

normal_subtypes <- c("qPSC", "smPSC")
activated_subtypes <- c("myCAF", "csCAF", "iCAF")

normal_evidence <- bind_rows(
  t1 %>%
    filter(state_direction == "Normal_up") %>%
    transmute(gene, source = "Test1_Normal", effect = state_direction),
  t2 %>%
    filter(CellType2 %in% normal_subtypes, subtype_direction == "higher_in_target") %>%
    transmute(gene, source = paste0("Test2_", CellType2), effect = subtype_direction),
  t3 %>%
    filter(CellType2 %in% normal_subtypes, support_pattern != "none") %>%
    transmute(gene, source = paste0("Test3_", CellType2), effect = support_pattern),
  t4 %>%
    filter(CellType2 %in% normal_subtypes, subtype_direction == "higher_in_target") %>%
    transmute(gene, source = paste0("Test4_", CellType2), effect = subtype_direction)
) %>%
  filter(grepl("^KCN", gene))

activated_evidence <- bind_rows(
  t1 %>%
    filter(state_direction == "Activated_up") %>%
    transmute(gene, source = "Test1_Activated", effect = state_direction),
  t2 %>%
    filter(CellType2 %in% activated_subtypes, subtype_direction == "higher_in_target") %>%
    transmute(gene, source = paste0("Test2_", CellType2), effect = subtype_direction),
  t3 %>%
    filter(CellType2 %in% activated_subtypes, support_pattern != "none") %>%
    transmute(gene, source = paste0("Test3_", CellType2), effect = support_pattern),
  t4 %>%
    filter(CellType2 %in% activated_subtypes, subtype_direction == "higher_in_target") %>%
    transmute(gene, source = paste0("Test4_", CellType2), effect = subtype_direction)
) %>%
  filter(grepl("^KCN", gene))

normal_genes <- sort(unique(normal_evidence$gene))
activated_genes <- sort(unique(activated_evidence$gene))
intersection_genes <- intersect(normal_genes, activated_genes)

pole_membership <- tibble(
  gene = sort(unique(c(normal_genes, activated_genes)))
) %>%
  mutate(
    Normal_pole = gene %in% normal_genes,
    Activated_pole = gene %in% activated_genes,
    membership = case_when(
      Normal_pole & Activated_pole ~ "Both",
      Normal_pole ~ "Normal_only",
      Activated_pole ~ "Activated_only",
      TRUE ~ "None"
    )
  ) %>%
  left_join(
    normal_evidence %>%
      group_by(gene) %>%
      summarise(normal_sources = paste(sort(unique(source)), collapse = "; "), .groups = "drop"),
    by = "gene"
  ) %>%
  left_join(
    activated_evidence %>%
      group_by(gene) %>%
      summarise(activated_sources = paste(sort(unique(source)), collapse = "; "), .groups = "drop"),
    by = "gene"
  )

save_tsv(pole_membership, file.path(table_dir, "KCN_normal_vs_activated_poles_membership.tsv"))

venn_grob <- venn.diagram(
  x = list(
    Normal_pole = normal_genes,
    Activated_pole = activated_genes
  ),
  filename = NULL,
  fill = c("#70D6FF", "#FF686B"),
  alpha = 0.78,
  col = "grey25",
  lwd = 1.2,
  cex = 1.5,
  fontface = "bold",
  cat.cex = 1.3,
  cat.fontface = "bold",
  category.names = c("Normal pole", "Activated pole"),
  cat.pos = c(-20, 20),
  cat.dist = c(0.03, 0.03),
  margin = 0.08
)

legend_df <- tibble(
  block = c("Normal pole", "Activated pole", "Intersection"),
  color = c("#70D6FF", "#FF686B", "#7A7A7A"),
  genes = c(
    if (length(normal_genes) == 0) "None" else wrap_text(paste(normal_genes, collapse = ", "), 40),
    if (length(activated_genes) == 0) "None" else wrap_text(paste(activated_genes, collapse = ", "), 40),
    if (length(intersection_genes) == 0) "None" else wrap_text(paste(intersection_genes, collapse = ", "), 40)
  )
) %>%
  mutate(
    y = rev(seq_len(n())),
    label = paste0(block, "\n", genes)
  )

legend_plot <- ggplot(legend_df, aes(x = 0, y = y)) +
  geom_point(aes(color = block), size = 6) +
  geom_text(aes(label = label), hjust = 0, nudge_x = 0.18, size = 3.4, lineheight = 0.95) +
  scale_color_manual(values = stats::setNames(legend_df$color, legend_df$block)) +
  coord_cartesian(xlim = c(0, 6.1), ylim = c(0.5, nrow(legend_df) + 0.5), clip = "off") +
  theme_void() +
  theme(legend.position = "none")

p_venn <- wrap_elements(full = venn_grob) + legend_plot + plot_layout(widths = c(1.05, 1))

save_plot(
  p_venn,
  file.path(figure_dir, "Figure_13_Venn_KCN_normal_vs_activated_poles.png"),
  width = 15,
  height = 8
)

message("Pole Venn saved in: ", figure_dir)

