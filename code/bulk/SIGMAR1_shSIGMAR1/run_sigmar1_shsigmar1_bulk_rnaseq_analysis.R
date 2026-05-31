###############################################################################
# Bulk RNA-seq analysis for the SigmaR1 / shSigmaR1 dataset
#
# This script works on the local CTRL vs shSigmaR1 bulk RNA-seq data. The goal
# is to summarize the perturbation with QC plots, differential expression,
# CAF-oriented signatures, and pathway-level results.
#
# Here the count-aware DESeq2 objects are available locally, so the differential
# expression step is run from those inputs.
###############################################################################
suppressPackageStartupMessages({
  library(DESeq2)
  library(SummarizedExperiment)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(apeglm)
  library(ggplot2)
  library(dplyr)
  library(fgsea)
  library(msigdbr)
  library(pheatmap)
  library(viridis)
  library(viridisLite)
})

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  file_path <- sub(file_arg, "", args[grep(file_arg, args)])
  if (length(file_path) == 0) {
    return(normalizePath(getwd()))
  }
  normalizePath(dirname(file_path))
}

script_dir <- get_script_dir()
input_dir <- file.path(script_dir, "inputs")
resource_dir <- file.path(script_dir, "resources")

# Main inputs: local DESeq2 objects for the SigmaR1 perturbation dataset.
dds_rdata <- file.path(input_dir, "dds.RData")
vsd_rdata <- file.path(input_dir, "vsd.RData")

# External marker resources used to reconstruct CAF-oriented signatures.
top30_icaf_path <- file.path(resource_dir, "genes_top30_iCAF_human.txt")
top30_mycaf_path <- file.path(resource_dir, "genes_top30_myCAF_human.txt")
top30_cscaf_path <- file.path(resource_dir, "genes_top30_csCAF_human.txt")
apcaf_markers_path <- file.path(resource_dir, "apCAFsmarkers.txt")

output_dir <- script_dir
fig_dir <- file.path(output_dir, "figures")
tab_dir <- file.path(output_dir, "tables")
rds_dir <- file.path(output_dir, "rds")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

message0 <- function(...) cat(..., "\n")

ctrl_color <- "#2166ac"
sh_color <- "#b2182b"
neutral_color <- "grey70"
diverging_cols <- grDevices::colorRampPalette(c(ctrl_color, "white", sh_color))(101)
seismic_cols <- grDevices::hcl.colors(101, "PuOr")
heat_blues <- grDevices::colorRampPalette(rev(c("#08306B", "#2171B5", "#6BAED6", "#C6DBEF", "#F7FBFF")))(100)

signature_colors <- c(
  SigmaR1 = "#2b2b2b",
  myCAF = "#8c2d04",
  iCAF = "#cb181d",
  csCAF = "#238b45",
  apCAF = "#6a51a3",
  ifCAF = "#dd3497"
)
sig_colors <- c("DE_FDR<0.05" = "#252525", "W_FDR<0.10" = "#525252", "p<0.05_only" = "#969696", "NS" = "#f0f0f0")

read_gene_set <- function(path) {
  if (!file.exists(path)) return(character(0))
  unique(toupper(trimws(readLines(path, warn = FALSE))))
}

safe_num <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

zscore_vec <- function(x) {
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

heatmap_scale_spec <- function(mat, n = 101) {
  vals <- as.numeric(mat)
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0) {
    return(list(colors = heat_blues, breaks = seq(0, 1, length.out = length(heat_blues) + 1), palette_name = "Blues"))
  }
  rng <- range(vals)
  if (rng[1] < 0 && rng[2] > 0) {
    lim <- max(abs(rng))
    return(list(colors = seismic_cols, breaks = seq(-lim, lim, length.out = length(seismic_cols) + 1), palette_name = "seismic"))
  }
  lower <- rng[1]
  upper <- rng[2]
  if (!is.finite(lower) || !is.finite(upper) || identical(lower, upper)) {
    lower <- 0
    upper <- max(1, upper, na.rm = TRUE)
  }
  if (upper <= lower) upper <- lower + 1e-6
  list(colors = heat_blues, breaks = seq(lower, upper, length.out = length(heat_blues) + 1), palette_name = "Blues")
}

save_portrait_block_heatmap <- function(display_mat,
                                        annotation_col,
                                        annotation_colors,
                                        gaps_row,
                                        block_labels,
                                        block_lengths,
                                        file_path,
                                        width = 7.5,
                                        height = 11.5,
                                        title_text = "CAF signature heatmap") {
  pal <- heatmap_scale_spec(display_mat)
  starts <- c(1, cumsum(block_lengths)[-length(block_lengths)] + 1)
  ends <- cumsum(block_lengths)
  centers_top <- (((starts + ends) / 2) - 0.5) / nrow(display_mat)
  centers_y <- 1 - centers_top

  draw_heatmap <- function() {
    ph <- pheatmap::pheatmap(
      display_mat,
      color = pal$colors,
      breaks = pal$breaks,
      cluster_rows = FALSE,
      cluster_cols = TRUE,
      clustering_method = "ward.D2",
      scale = "none",
      border_color = NA,
      annotation_col = annotation_col,
      annotation_colors = annotation_colors,
      annotation_legend = FALSE,
      gaps_row = gaps_row,
      treeheight_row = 0,
      show_rownames = TRUE,
      show_colnames = TRUE,
      angle_col = 90,
      cellwidth = 22,
      cellheight = 8,
      fontsize = 9,
      fontsize_row = 8.5,
      fontsize_col = 8,
      main = "",
      silent = TRUE
    )
    gt <- ph[["gtable"]]
    drop_idx <- which(gt$layout$name %in% c("annotation_legend", "col_annotation_names"))
    if (length(drop_idx) > 0) {
      gt$grobs <- gt$grobs[-drop_idx]
      gt$layout <- gt$layout[-drop_idx, , drop = FALSE]
    }
    matrix_idx <- which(gt$layout$name == "matrix")[1]
    matrix_t <- gt$layout$t[matrix_idx]
    matrix_b <- gt$layout$b[matrix_idx]
    matrix_l <- gt$layout$l[matrix_idx]
    matrix_r <- gt$layout$r[matrix_idx]
    gt <- gtable::gtable_add_rows(gt, heights = grid::unit(0.65, "cm"), pos = 0)
    gt <- gtable::gtable_add_cols(gt, widths = grid::unit(0.9, "cm"), pos = matrix_l - 1)
    gt <- gtable::gtable_add_grob(
      gt,
      grobs = grid::textGrob(title_text, gp = grid::gpar(fontsize = 14, fontface = "bold")),
      t = 1, l = matrix_l + 1, r = matrix_r + 1, clip = "off", name = "portrait_title"
    )
    label_grobs <- lapply(seq_along(block_labels), function(i) {
      grid::textGrob(
        block_labels[i],
        x = grid::unit(0.5, "npc"),
        y = grid::unit(centers_y[i], "npc"),
        rot = 90,
        gp = grid::gpar(fontsize = 10, fontface = "bold")
      )
    })
    gt <- gtable::gtable_add_grob(
      gt,
      grobs = label_grobs,
      t = matrix_t + 1, b = matrix_b + 1, l = matrix_l, r = matrix_l,
      clip = "off", name = paste0("block_label_", seq_along(block_labels))
    )
    grid::grid.newpage()
    grid::grid.draw(gt)
  }

  if (grepl("\\.pdf$", file_path, ignore.case = TRUE)) {
    grDevices::pdf(file_path, width = width, height = height, useDingbats = FALSE)
    draw_heatmap()
    grDevices::dev.off()
  } else {
    grDevices::png(file_path, width = width, height = height, units = "in", res = 300)
    draw_heatmap()
    grDevices::dev.off()
  }
}

collapse_matrix_by_mean <- function(mat, groups) {
  rs <- rowsum(mat, group = groups)
  denom <- as.numeric(table(groups)[rownames(rs)])
  rs / denom
}

collapse_res_table <- function(df) {
  ord <- order(
    is.na(df$padj),
    ifelse(is.na(df$padj), Inf, df$padj),
    -abs(df$log2FoldChange),
    decreasing = FALSE
  )
  df <- df[ord, , drop = FALSE]
  df[!duplicated(df$gene_symbol), , drop = FALSE]
}

fmt_p <- function(x) {
  if (!is.finite(x)) return("NA")
  if (x < 1e-4) return(format(x, scientific = TRUE, digits = 2))
  sprintf("%.4f", x)
}

sig_star <- function(p = NA_real_, padj = NA_real_) {
  if (is.finite(padj) && padj < 0.05) return("**")
  if (is.finite(p) && p < 0.05) return("*")
  ""
}

condition_from_sample <- function(x) {
  ifelse(grepl("^sh11", x, ignore.case = TRUE), "shSigmaR1", "CTRL")
}

replicate_from_sample <- function(x) {
  sub(".*?([0-9]+)$", "\\1", x)
}

sample_summary_stats <- function(mat) {
  data.frame(
    sample = colnames(mat),
    mean = apply(mat, 2, mean, na.rm = TRUE),
    median = apply(mat, 2, stats::median, na.rm = TRUE),
    sd = apply(mat, 2, stats::sd, na.rm = TRUE),
    iqr = apply(mat, 2, stats::IQR, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

score_signature_z <- function(mat, genes) {
  genes <- intersect(toupper(genes), rownames(mat))
  if (length(genes) == 0) return(rep(NA_real_, ncol(mat)))
  zmat <- t(apply(mat[genes, , drop = FALSE], 1, zscore_vec))
  colMeans(zmat, na.rm = TRUE)
}

score_signature_raw <- function(mat, genes) {
  genes <- intersect(toupper(genes), rownames(mat))
  if (length(genes) == 0) return(rep(NA_real_, ncol(mat)))
  colMeans(mat[genes, , drop = FALSE], na.rm = TRUE)
}

rank_biserial <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (length(x) == 0 || length(y) == 0) return(NA_real_)
  r <- rank(c(x, y))
  rx <- sum(r[seq_along(x)])
  u <- rx - length(x) * (length(x) + 1) / 2
  2 * u / (length(x) * length(y)) - 1
}

safe_wilcox_p <- function(x, group) {
  ok <- is.finite(x) & !is.na(group)
  x <- x[ok]
  group <- droplevels(as.factor(group[ok]))
  if (length(unique(group)) != 2) return(NA_real_)
  if (any(table(group) < 2)) return(NA_real_)
  stats::wilcox.test(x ~ group, exact = FALSE)$p.value
}

safe_spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 4) {
    return(c(rho = NA_real_, pvalue = NA_real_, n = length(x)))
  }
  if (stats::sd(x) == 0 || stats::sd(y) == 0) {
    return(c(rho = NA_real_, pvalue = NA_real_, n = length(x)))
  }
  ct <- suppressWarnings(stats::cor.test(x, y, method = "spearman"))
  c(rho = unname(ct$estimate), pvalue = ct$p.value, n = length(x))
}

sig_label <- function(de_padj, w_p, w_padj) {
  if (is.finite(de_padj) && de_padj < 0.05) return("DE_FDR<0.05")
  if (is.finite(w_padj) && w_padj < 0.10) return("W_FDR<0.10")
  if (is.finite(w_p) && w_p < 0.05) return("p<0.05_only")
  "NS"
}

cluster_order <- function(mat) {
  if (nrow(mat) <= 2) return(rownames(mat))
  rownames(mat)[stats::hclust(stats::dist(mat), method = "ward.D2")$order]
}

draw_violin <- function(values, at, fill_col, border_col = NA, width = 0.35) {
  values <- values[is.finite(values)]
  if (length(values) == 0) return(invisible(NULL))
  if (length(unique(values)) == 1) {
    graphics::rect(at - width / 2, values[1] - 0.01, at + width / 2, values[1] + 0.01,
                   col = grDevices::adjustcolor(fill_col, alpha.f = 0.45), border = border_col)
    return(invisible(NULL))
  }
  dens <- stats::density(values, na.rm = TRUE, bw = "nrd0")
  dens$y <- dens$y / max(dens$y) * width
  graphics::polygon(
    c(at - dens$y, rev(at + dens$y)),
    c(dens$x, rev(dens$x)),
    col = grDevices::adjustcolor(fill_col, alpha.f = 0.45),
    border = border_col
  )
}

plot_density_lines <- function(mat, phenotype, file_path) {
  pdf(file_path, width = 9, height = 6)
  xlim <- range(mat, na.rm = TRUE)
  ylim <- c(0, 1.05 * max(sapply(seq_len(ncol(mat)), function(i) max(stats::density(mat[, i], na.rm = TRUE)$y))))
  plot(NA, xlim = xlim, ylim = ylim, xlab = "log2(count + 1)", ylab = "Density",
       main = "Per-sample density curves")
  for (i in seq_len(ncol(mat))) {
    col_i <- if (phenotype$condition[i] == "CTRL") ctrl_color else sh_color
    lines(stats::density(mat[, i], na.rm = TRUE), col = grDevices::adjustcolor(col_i, alpha.f = 0.65), lwd = 2)
  }
  legend("topright", legend = c("CTRL", "shSigmaR1"), col = c(ctrl_color, sh_color), lwd = 2, bty = "n")
  dev.off()
}

plot_qc_boxplots <- function(mat, phenotype, file_path) {
  pdf(file_path, width = 9, height = 6)
  cols <- ifelse(phenotype$condition == "CTRL", ctrl_color, sh_color)
  boxplot(mat, las = 2, col = cols, outline = FALSE,
          ylab = "log2(count + 1)", main = "Per-sample expression distributions")
  legend("topright", legend = c("CTRL", "shSigmaR1"), fill = c(ctrl_color, sh_color), bty = "n")
  dev.off()
}

plot_qc_distributions <- function(expr_raw, expr_log, file_path) {
  pdf(file_path, width = 10, height = 5)
  par(mfrow = c(1, 2))
  hist(as.vector(expr_raw), breaks = 200, main = "Normalized counts", xlab = "Counts",
       col = "steelblue", border = "white")
  hist(as.vector(expr_log), breaks = 200, main = "log2(count + 1)", xlab = "Expression",
       col = "darkorange", border = "white")
  par(mfrow = c(1, 1))
  dev.off()
}

plot_sample_correlation_heatmap <- function(mat, file_path) {
  corr <- stats::cor(mat, method = "pearson")
  pdf(file_path, width = 7, height = 6)
  heatmap(corr, Rowv = NA, Colv = NA, scale = "none", col = heat_blues,
          margins = c(8, 8), main = "Sample-sample correlation")
  dev.off()
  corr
}

plot_sample_distance_heatmap <- function(mat, file_path) {
  d <- stats::dist(t(mat))
  dm <- as.matrix(d)
  pdf(file_path, width = 7, height = 6)
  heatmap(dm, Rowv = NA, Colv = NA, scale = "none", col = heat_blues,
          margins = c(8, 8), main = "Euclidean sample distance")
  dev.off()
  dm
}

plot_pca <- function(mat, phenotype, file_path) {
  pca <- stats::prcomp(t(mat), scale. = TRUE)
  percent <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)
  pdf(file_path, width = 7, height = 6)
  plot(pca$x[, 1], pca$x[, 2],
       pch = ifelse(phenotype$condition == "CTRL", 21, 24),
       bg = ifelse(phenotype$condition == "CTRL", ctrl_color, sh_color),
       col = "black", cex = 1.8,
       xlab = paste0("PC1 (", percent[1], "%)"),
       ylab = paste0("PC2 (", percent[2], "%)"),
       main = "PCA - bulk SigmaR1 / shSigmaR1")
  text(pca$x[, 1], pca$x[, 2], labels = phenotype$sample, pos = 3, cex = 0.8)
  legend("topright", legend = c("CTRL", "shSigmaR1"), pt.bg = c(ctrl_color, sh_color),
         pch = c(21, 24), bty = "n")
  dev.off()
  data.frame(sample = rownames(pca$x), PC1 = pca$x[, 1], PC2 = pca$x[, 2],
             PC3 = pca$x[, 3], condition = phenotype$condition,
             stringsAsFactors = FALSE,
             variance_PC1 = percent[1], variance_PC2 = percent[2], variance_PC3 = percent[3])
}

plot_pca_variance <- function(mat, file_path) {
  pca <- stats::prcomp(t(mat), scale. = TRUE)
  var_exp <- pca$sdev^2 / sum(pca$sdev^2)
  pdf(file_path, width = 6, height = 5)
  barplot(var_exp[1:min(10, length(var_exp))], names.arg = paste0("PC", seq_len(min(10, length(var_exp)))),
          col = "steelblue", ylab = "Variance explained", main = "PCA variance")
  dev.off()
  data.frame(PC = paste0("PC", seq_along(var_exp)), variance = var_exp, stringsAsFactors = FALSE)
}

plot_volcano <- function(res_df, file_path) {
  plot_df <- res_df[is.finite(res_df$padj) & is.finite(res_df$log2FoldChange), , drop = FALSE]
  plot_df$y <- -log10(pmax(plot_df$padj, 1e-300))
  plot_df$class <- "NS"
  plot_df$class[plot_df$padj < 0.05 & plot_df$log2FoldChange > 1] <- "Up in shSigmaR1"
  plot_df$class[plot_df$padj < 0.05 & plot_df$log2FoldChange < -1] <- "Down in shSigmaR1"
  cols <- c("Down in shSigmaR1" = ctrl_color, "NS" = neutral_color, "Up in shSigmaR1" = sh_color)
  labs_df <- plot_df[order(plot_df$padj, -abs(plot_df$log2FoldChange)), , drop = FALSE]
  labs_df <- unique(rbind(
    labs_df[seq_len(min(12, nrow(labs_df))), ],
    plot_df[plot_df$gene_symbol %in% c("SIGMAR1", "ACTA2", "IL6", "C3", "CD74"), , drop = FALSE]
  ))
  pdf(file_path, width = 8, height = 7)
  plot(plot_df$log2FoldChange, plot_df$y,
       pch = 16, cex = 0.75, col = cols[plot_df$class],
       xlab = "log2 Fold Change (shSigmaR1 / CTRL)",
       ylab = "-log10(FDR)",
       main = "Volcano plot - shSigmaR1 vs CTRL")
  abline(v = c(-1, 1), lty = 2, col = "grey50")
  abline(h = -log10(0.05), lty = 2, col = "grey50")
  if (nrow(labs_df) > 0) {
    text(labs_df$log2FoldChange, labs_df$y, labels = labs_df$gene_symbol, pos = 3, cex = 0.7)
  }
  legend("topright", legend = names(cols), col = cols, pch = 16, bty = "n")
  dev.off()
}

plot_ma_like <- function(res_df, file_path) {
  pdf(file_path, width = 7, height = 6)
  cols <- ifelse(res_df$padj < 0.05 & res_df$log2FoldChange > 0, sh_color,
                 ifelse(res_df$padj < 0.05 & res_df$log2FoldChange < 0, ctrl_color, neutral_color))
  plot(res_df$mean_expression, res_df$log2FoldChange,
       pch = 16, cex = 0.7, col = cols,
       xlab = "Mean log2 expression", ylab = "log2 Fold Change (shSigmaR1 / CTRL)",
       main = "MA-like plot")
  abline(h = 0, lty = 2, col = "grey40")
  legend("topright", legend = c("Down in shSigmaR1", "NS", "Up in shSigmaR1"),
         col = c(ctrl_color, neutral_color, sh_color), pch = 16, bty = "n")
  dev.off()
}

plot_matrix_heatmap <- function(mat, file_path, main, row_side_cols = NULL, col_side_cols = NULL,
                                margins = c(8, 8), row_cex = 0.7, col_cex = 0.9) {
  pdf(file_path, width = 8, height = max(6, 0.16 * nrow(mat)))
  heatmap(mat, Rowv = NA, Colv = NA, scale = "none", col = diverging_cols,
          RowSideColors = row_side_cols, ColSideColors = col_side_cols,
          margins = margins, cexRow = row_cex, cexCol = col_cex, main = main)
  dev.off()
}

plot_marker_violin_panels <- function(expr_mat, phenotype, genes, marker_stats, file_path) {
  genes <- intersect(genes, rownames(expr_mat))
  if (length(genes) == 0) return(invisible(NULL))
  n_panels <- length(genes)
  ncol <- 4
  nrow <- ceiling(n_panels / ncol)
  pdf(file_path, width = 14, height = max(8, 2.8 * nrow))
  par(mfrow = c(nrow, ncol), mar = c(4, 4, 3, 1))
  for (g in genes) {
    vals_ctrl <- expr_mat[g, phenotype$condition == "CTRL"]
    vals_sh <- expr_mat[g, phenotype$condition == "shSigmaR1"]
    y_all <- c(vals_ctrl, vals_sh)
    ylim <- range(y_all, na.rm = TRUE)
    ylim[2] <- ylim[2] + 0.22 * diff(ylim + c(0, 1e-6))
    plot(NA, xlim = c(0.5, 2.5), ylim = ylim, xaxt = "n", xlab = "", ylab = "log2(count + 1)",
         main = g, cex.main = 1)
    axis(1, at = c(1, 2), labels = c("CTRL", "shSigmaR1"))
    draw_violin(vals_ctrl, at = 1, fill_col = ctrl_color)
    draw_violin(vals_sh, at = 2, fill_col = sh_color)
    boxplot(list(vals_ctrl, vals_sh), add = TRUE, at = c(1, 2), xaxt = "n", yaxt = "n",
            outline = FALSE, boxwex = 0.12, border = "black", col = NA)
    stripchart(list(vals_ctrl, vals_sh), add = TRUE, at = c(1, 2), vertical = TRUE,
               method = "jitter", jitter = 0.08, pch = 21,
               bg = c(rep(ctrl_color, length(vals_ctrl)), rep(sh_color, length(vals_sh))),
               col = "black", cex = 0.9)
    st <- marker_stats[marker_stats$gene_symbol == g, , drop = FALSE]
    if (nrow(st) == 1) {
      title(main = paste0(g, sig_star(st$wilcox_p[1], st$wilcox_padj[1])), line = 0.2, cex.main = 1)
      txt <- c(
        paste0("LFC=", round(st$log2FoldChange[1], 2)),
        paste0("DE FDR=", fmt_p(st$padj[1])),
        paste0("W p/FDR=", fmt_p(st$wilcox_p[1]), "/", fmt_p(st$wilcox_padj[1]))
      )
      usr <- par("usr")
      text(usr[1] + 0.04 * diff(usr[1:2]), usr[4] - 0.03 * diff(usr[3:4]),
           labels = paste(txt, collapse = "\n"), adj = c(0, 1), cex = 0.72)
    }
  }
  dev.off()
}

plot_bulk_dotplot <- function(dot_df, file_path, main) {
  genes <- unique(dot_df$gene)
  conds <- c("CTRL", "shSigmaR1")
  x_pos <- match(dot_df$condition, conds)
  y_pos <- match(dot_df$gene, rev(genes))
  color_values <- diverging_cols[pmin(101, pmax(1, round((dot_df$mean_z + 3) / 6 * 100) + 1))]
  cex_values <- 0.9 + 2.5 * dot_df$frac_high
  pdf(file_path, width = 7, height = max(6, 0.24 * length(genes)))
  plot(NA, xlim = c(0.5, 2.5), ylim = c(0.5, length(genes) + 0.5), xaxt = "n", yaxt = "n",
       xlab = "", ylab = "", main = main)
  axis(1, at = 1:2, labels = conds)
  axis(2, at = seq_along(genes), labels = rev(genes), las = 2, cex.axis = 0.75)
  abline(v = 1:2, col = "grey90", lty = 3)
  abline(h = seq_along(genes), col = "grey95", lty = 3)
  points(x_pos, y_pos, pch = 21, bg = color_values, col = "black", cex = cex_values)
  legend("topright", legend = c("High z-score", "Low z-score"), pt.bg = c(sh_color, ctrl_color),
         pch = 21, bty = "n", cex = 0.8)
  text(2.45, 1, labels = "Point size = fraction of replicates\nabove gene-wise median", adj = c(1, 0), cex = 0.7)
  dev.off()
}

plot_signature_score_panels <- function(score_mat, phenotype, stats_df, file_path) {
  sigs <- colnames(score_mat)
  pdf(file_path, width = 11, height = 6)
  par(mfrow = c(2, ceiling(length(sigs) / 2)), mar = c(4, 4, 3, 1))
  for (sig in sigs) {
    x_ctrl <- score_mat[phenotype$condition == "CTRL", sig]
    x_sh <- score_mat[phenotype$condition == "shSigmaR1", sig]
    y_all <- c(x_ctrl, x_sh)
    ylim <- range(y_all, na.rm = TRUE)
    ylim[2] <- ylim[2] + 0.2 * diff(ylim + c(0, 1e-6))
    plot(NA, xlim = c(0.5, 2.5), ylim = ylim, xaxt = "n", xlab = "", ylab = "Signature z-score",
         main = sig)
    axis(1, at = c(1, 2), labels = c("CTRL", "shSigmaR1"))
    draw_violin(x_ctrl, 1, ctrl_color)
    draw_violin(x_sh, 2, sh_color)
    boxplot(list(x_ctrl, x_sh), add = TRUE, at = c(1, 2), xaxt = "n", yaxt = "n",
            outline = FALSE, boxwex = 0.12, border = "black", col = NA)
    stripchart(list(x_ctrl, x_sh), add = TRUE, at = c(1, 2), vertical = TRUE,
               method = "jitter", jitter = 0.08, pch = 21,
               bg = c(rep(ctrl_color, length(x_ctrl)), rep(sh_color, length(x_sh))),
               col = "black", cex = 0.9)
    st <- stats_df[stats_df$signature == sig, , drop = FALSE]
    if (nrow(st) == 1) {
      title(main = paste0(sig, sig_star(st$wilcox_p[1], st$wilcox_padj[1])), line = 0.2, cex.main = 1)
      txt <- c(
        paste0("Delta=", round(st$mean_shSigmaR1[1] - st$mean_CTRL[1], 2)),
        paste0("W p/FDR=", fmt_p(st$wilcox_p[1]), "/", fmt_p(st$wilcox_padj[1]))
      )
      usr <- par("usr")
      text(usr[1] + 0.04 * diff(usr[1:2]), usr[4] - 0.04 * diff(usr[3:4]),
           labels = paste(txt, collapse = "\n"), adj = c(0, 1), cex = 0.75)
    }
  }
  dev.off()
}

plot_condition_matrixplot <- function(score_condition_mat, file_path, main) {
  zlim <- max(abs(score_condition_mat), na.rm = TRUE)
  nr <- nrow(score_condition_mat)
  nc <- ncol(score_condition_mat)
  pdf(file_path, width = 5.5, height = max(4, 0.45 * nr))
  par(mar = c(5, 10, 3, 2))
  image(
    x = seq_len(nc),
    y = seq_len(nr),
    z = t(score_condition_mat[nr:1, , drop = FALSE]),
    col = diverging_cols,
    zlim = c(-zlim, zlim),
    axes = FALSE,
    xlab = "", ylab = "", main = main
  )
  axis(1, at = seq_len(nc), labels = colnames(score_condition_mat))
  axis(2, at = seq_len(nr), labels = rev(rownames(score_condition_mat)), las = 2, cex.axis = 0.85)
  box()
  dev.off()
}

plot_gene_contribution_barplot <- function(contrib_df, file_path) {
  if (nrow(contrib_df) == 0) return(invisible(NULL))
  contrib_df <- contrib_df[order(contrib_df$signature, contrib_df$delta), , drop = FALSE]
  contrib_df$gene_facet <- paste(contrib_df$signature, contrib_df$gene, sep = "___")
  contrib_df$gene_facet <- factor(contrib_df$gene_facet, levels = contrib_df$gene_facet)
  contrib_df$gene_label <- contrib_df$gene
  p <- ggplot(contrib_df, aes(x = gene_facet, y = delta, fill = delta)) +
    geom_col(width = 0.78, colour = "black", linewidth = 0.15) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey45") +
    coord_flip() +
    facet_wrap(~ signature, scales = "free_y", ncol = 2) +
    scale_x_discrete(labels = contrib_df$gene_label) +
    scale_fill_gradient2(low = ctrl_color, mid = "white", high = sh_color, midpoint = 0, name = "Delta\n(shSigmaR1 - CTRL)") +
    labs(
      title = "Top marker contributions by CAF signature",
      subtitle = "Bars show the mean log2(count + 1) difference between shSigmaR1 and CTRL",
      x = NULL,
      y = "Mean expression difference"
    ) +
    theme_stage() +
    theme(axis.text.y = element_text(size = 7))
  save_plot(p, file_path, 13.5, max(8.5, 2 + 0.22 * nrow(contrib_df)))
}

plot_sigmar1_marker_summary_matrix <- function(summary_mat, file_path, main) {
  nr <- nrow(summary_mat)
  nc <- ncol(summary_mat)
  zlim <- max(abs(summary_mat), na.rm = TRUE)
  pdf(file_path, width = 6.5, height = max(6, 0.24 * nr))
  par(mar = c(5, 12, 3, 2))
  image(seq_len(nc), seq_len(nr), t(summary_mat[nr:1, , drop = FALSE]),
        col = diverging_cols, zlim = c(-zlim, zlim), axes = FALSE, xlab = "", ylab = "", main = main)
  axis(1, at = seq_len(nc), labels = colnames(summary_mat))
  axis(2, at = seq_len(nr), labels = rev(rownames(summary_mat)), las = 2, cex.axis = 0.75)
  box()
  dev.off()
}

plot_key_scatter_panels <- function(expr_mat, phenotype, pairs_df, file_path) {
  genes <- pairs_df$gene_symbol
  genes <- genes[genes %in% rownames(expr_mat)]
  if (length(genes) == 0) return(invisible(NULL))
  ncol <- 3
  nrow <- ceiling(length(genes) / ncol)
  x <- expr_mat["SIGMAR1", phenotype$sample]
  pdf(file_path, width = 12, height = max(6, 3.5 * nrow))
  par(mfrow = c(nrow, ncol), mar = c(4, 4, 3, 1))
  for (g in genes) {
    y <- expr_mat[g, phenotype$sample]
    cols <- ifelse(phenotype$condition == "CTRL", ctrl_color, sh_color)
    plot(x, y, pch = 21, bg = cols, col = "black", cex = 1.2,
         xlab = "SIGMAR1 log2(count + 1)", ylab = paste0(g, " log2(count + 1)"),
         main = g)
    text(x, y, labels = phenotype$sample, pos = 3, cex = 0.7)
    fit <- stats::lm(y ~ x)
    abline(fit, col = "grey30", lwd = 2)
    st <- pairs_df[pairs_df$gene_symbol == g, , drop = FALSE]
    if (nrow(st) == 1) {
      title(main = paste0(g, sig_star(st$rho_pvalue[1], st$rho_padj[1])), line = 0.2, cex.main = 1)
      usr <- par("usr")
      text(usr[1] + 0.04 * diff(usr[1:2]), usr[4] - 0.04 * diff(usr[3:4]),
           labels = paste0("rho=", round(st$rho_all[1], 2), "\np/FDR=", fmt_p(st$rho_pvalue[1]), "/", fmt_p(st$rho_padj[1])),
           adj = c(0, 1), cex = 0.72)
    }
  }
  dev.off()
}

plot_radial_sigmar1_network <- function(edge_df, node_categories, file_path) {
  if (nrow(edge_df) == 0) return(invisible(NULL))
  node_df <- data.frame(gene = c("SIGMAR1", edge_df$gene_symbol), stringsAsFactors = FALSE)
  node_df <- node_df[!duplicated(node_df$gene), , drop = FALSE]
  outer_nodes <- node_df$gene[node_df$gene != "SIGMAR1"]
  angles <- seq(0, 2 * pi, length.out = length(outer_nodes) + 1)[-1]
  coords <- data.frame(
    gene = c("SIGMAR1", outer_nodes),
    x = c(0, cos(angles)),
    y = c(0, sin(angles)),
    stringsAsFactors = FALSE
  )
  coords$category <- node_categories$category[match(coords$gene, node_categories$gene_symbol)]
  coords$category[is.na(coords$category)] <- "SigmaR1"
  coords$col <- signature_colors[coords$category]
  pdf(file_path, width = 8, height = 8)
  plot(NA, xlim = c(-1.35, 1.35), ylim = c(-1.35, 1.35), xaxt = "n", yaxt = "n",
       xlab = "", ylab = "", bty = "n",
       main = "SigmaR1-marker correlation network\n(overall Spearman across 6 samples; exploratory)")
  for (i in seq_len(nrow(edge_df))) {
    g <- edge_df$gene_symbol[i]
    rho <- edge_df$rho_all[i]
    xy <- coords[coords$gene == g, , drop = FALSE]
    graphics::segments(0, 0, xy$x, xy$y,
                       col = if (rho >= 0) sh_color else ctrl_color,
                       lwd = 1 + 4 * abs(rho))
  }
  points(coords$x, coords$y, pch = 21,
         bg = coords$col, col = "black",
         cex = ifelse(coords$gene == "SIGMAR1", 3.2, 2.2))
  text(coords$x, coords$y, labels = coords$gene, pos = 3, cex = 0.8)
  legend("topright", legend = names(signature_colors), pt.bg = signature_colors,
         pch = 21, bty = "n", cex = 0.8)
  legend("bottomleft", legend = c("Positive rho", "Negative rho"),
         col = c(sh_color, ctrl_color), lwd = 3, bty = "n")
  dev.off()
}

compute_preranked_gsea <- function(ranks, pathways, nperm = 2000, seed = 42) {
  set.seed(seed)
  ranks <- sort(ranks, decreasing = TRUE)
  gene_names <- names(ranks)
  abs_stats <- abs(ranks)

  calc_es <- function(hit_idx) {
    hit_vec <- rep(FALSE, length(ranks))
    hit_vec[hit_idx] <- TRUE
    nh <- sum(hit_vec)
    if (nh < 2 || nh >= length(ranks)) {
      return(list(ES = NA_real_, running = rep(NA_real_, length(ranks)), leading_edge = character(0)))
    }
    hit_weights <- abs_stats
    hit_weights[!hit_vec] <- 0
    Phit <- cumsum(hit_weights / sum(hit_weights))
    Pmiss <- cumsum((!hit_vec) / sum(!hit_vec))
    running <- Phit - Pmiss
    max_pos <- max(running)
    min_neg <- min(running)
    if (abs(max_pos) >= abs(min_neg)) {
      peak_idx <- which.max(running)
      es <- max_pos
      leading_edge <- gene_names[seq_len(peak_idx)][hit_vec[seq_len(peak_idx)]]
    } else {
      peak_idx <- which.min(running)
      es <- min_neg
      leading_edge <- gene_names[seq_len(peak_idx)][hit_vec[seq_len(peak_idx)]]
    }
    list(ES = es, running = running, leading_edge = leading_edge)
  }

  out <- vector("list", length(pathways))
  for (i in seq_along(pathways)) {
    pw <- names(pathways)[i]
    genes <- intersect(unique(toupper(pathways[[i]])), gene_names)
    if (length(genes) < 5) {
      out[[i]] <- data.frame(
        pathway = pw, size = length(genes), ES = NA_real_, NES = NA_real_,
        pvalue = NA_real_, padj = NA_real_, direction = NA_character_,
        leading_edge = "", stringsAsFactors = FALSE
      )
      next
    }
    hit_idx <- which(gene_names %in% genes)
    obs <- calc_es(hit_idx)
    null_es <- numeric(nperm)
    for (b in seq_len(nperm)) {
      perm_idx <- sort(sample.int(length(ranks), length(hit_idx), replace = FALSE))
      null_es[b] <- calc_es(perm_idx)$ES
    }
    if (obs$ES >= 0) {
      pos_null <- null_es[null_es >= 0]
      pvalue <- if (length(pos_null) == 0) NA_real_ else mean(pos_null >= obs$ES)
      nes <- if (length(pos_null) == 0) NA_real_ else obs$ES / mean(pos_null)
      direction <- "Enriched in shSigmaR1"
    } else {
      neg_null <- abs(null_es[null_es < 0])
      pvalue <- if (length(neg_null) == 0) NA_real_ else mean(neg_null >= abs(obs$ES))
      nes <- if (length(neg_null) == 0) NA_real_ else obs$ES / mean(neg_null)
      direction <- "Enriched in CTRL"
    }
    out[[i]] <- data.frame(
      pathway = pw,
      size = length(genes),
      ES = obs$ES,
      NES = nes,
      pvalue = pvalue,
      direction = direction,
      leading_edge = paste(obs$leading_edge, collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
  res <- do.call(rbind, out)
  res$padj <- stats::p.adjust(res$pvalue, method = "BH")
  res[order(res$padj, -abs(res$NES)), , drop = FALSE]
}

plot_gsea_bubble <- function(gsea_res, file_path) {
  df <- gsea_res[is.finite(gsea_res$NES), , drop = FALSE]
  if (nrow(df) == 0) return(invisible(NULL))
  df <- df[order(df$NES), , drop = FALSE]
  y_pos <- seq_len(nrow(df))
  cols <- ifelse(df$NES > 0, sh_color, ctrl_color)
  sizes <- 0.8 + 3 * pmin(1, -log10(pmax(df$padj, 1e-8)) / 5)
  pdf(file_path, width = 8, height = max(5, 0.35 * nrow(df)))
  plot(df$NES, y_pos, pch = 21, bg = cols, cex = sizes, col = "black",
       yaxt = "n", xlab = "NES", ylab = "", main = "Preranked enrichment of CAF programs")
  axis(2, at = y_pos, labels = df$pathway, las = 2, cex.axis = 0.8)
  abline(v = 0, lty = 2, col = "grey50")
  legend("bottomright", legend = c("Enriched in CTRL", "Enriched in shSigmaR1"),
         pt.bg = c(ctrl_color, sh_color), pch = 21, bty = "n")
  text(max(df$NES, na.rm = TRUE), min(y_pos), labels = "Point size = -log10(FDR)", pos = 2, cex = 0.75)
  dev.off()
}

plot_enrichment_curve <- function(ranks, genes, title_text, file_path) {
  genes <- intersect(names(ranks), toupper(genes))
  if (length(genes) < 5) return(invisible(NULL))
  ranks <- sort(ranks, decreasing = TRUE)
  gene_names <- names(ranks)
  hit_idx <- which(gene_names %in% genes)
  abs_stats <- abs(ranks)
  hit_vec <- seq_along(ranks) %in% hit_idx
  hit_weights <- abs_stats
  hit_weights[!hit_vec] <- 0
  Phit <- cumsum(hit_weights / sum(hit_weights))
  Pmiss <- cumsum((!hit_vec) / sum(!hit_vec))
  running <- Phit - Pmiss
  pdf(file_path, width = 7, height = 5)
  plot(running, type = "l", lwd = 2, col = if (sum(ranks[hit_idx]) >= 0) sh_color else ctrl_color,
       xlab = "Ranked genes", ylab = "Running enrichment score", main = title_text)
  abline(h = 0, lty = 2, col = "grey50")
  rug(hit_idx, col = "black")
  dev.off()
}

###############################################################################
# Plotting helper functions
###############################################################################

theme_stage <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      legend.position = "right",
      legend.box = "vertical",
      strip.background = element_rect(fill = "grey95", colour = "grey70"),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold")
    )
}

save_plot <- function(plot_obj, file_path, width, height) {
  ggplot2::ggsave(filename = file_path, plot = plot_obj, width = width, height = height, units = "in")
  png_path <- sub("\\.pdf$", ".png", file_path)
  if (!identical(png_path, file_path)) {
    ggplot2::ggsave(filename = png_path, plot = plot_obj, width = width, height = height, units = "in", dpi = 300)
  }
}

plot_gene_contribution_barplot <- function(contrib_df, file_path) {
  if (nrow(contrib_df) == 0) return(invisible(NULL))
  contrib_df <- contrib_df[order(contrib_df$signature, contrib_df$delta), , drop = FALSE]
  contrib_df$gene_facet <- paste(contrib_df$signature, contrib_df$gene, sep = "___")
  contrib_df$gene_facet <- factor(contrib_df$gene_facet, levels = contrib_df$gene_facet)
  contrib_df$gene_label <- contrib_df$gene
  p <- ggplot(contrib_df, aes(x = gene_facet, y = delta, fill = delta)) +
    geom_col(width = 0.78, colour = "black", linewidth = 0.15) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey45") +
    coord_flip() +
    facet_wrap(~ signature, scales = "free_y", ncol = 2) +
    scale_x_discrete(labels = contrib_df$gene_label) +
    scale_fill_gradient2(low = ctrl_color, mid = "white", high = sh_color, midpoint = 0, name = "Delta\n(shSigmaR1 - CTRL)") +
    labs(
      title = "Top marker contributions by CAF signature",
      subtitle = "Bars show the mean log2(count + 1) difference between shSigmaR1 and CTRL",
      x = NULL,
      y = "Mean expression difference"
    ) +
    theme_stage() +
    theme(axis.text.y = element_text(size = 7))
  save_plot(p, file_path, 13.5, max(8.5, 2 + 0.22 * nrow(contrib_df)))
}

###############################################################################
# 1. Load inputs and define sample metadata
###############################################################################

message0("=== SigmaR1 / shSigmaR1 bulk RNA-seq analysis ===")
if (!file.exists(dds_rdata)) stop("dds.RData not found: ", dds_rdata)
if (!file.exists(vsd_rdata)) stop("vsd.RData not found: ", vsd_rdata)

load(dds_rdata)
load(vsd_rdata)

if (!exists("dds")) stop("Object 'dds' not found after loading dds.RData")
if (!exists("vsd3")) stop("Object 'vsd3' not found after loading vsd.RData")

# Sample metadata are recovered directly from the DESeq2 colData. This is the
# authoritative source for condition labels and replicate annotations.
coldata_df <- as.data.frame(SummarizedExperiment::colData(dds))
sample_cols <- rownames(coldata_df)
ctrl_cols <- sample_cols[tolower(coldata_df$condition) == "ctrl"]
sh_cols <- sample_cols[tolower(coldata_df$condition) == "sh11"]

if (length(sh_cols) == 0 || length(ctrl_cols) == 0) stop("Unable to identify CTRL/sh11 sample columns in dds colData.")

phenotype <- data.frame(
  sample = sample_cols,
  condition = factor(ifelse(sample_cols %in% ctrl_cols, "CTRL", "shSigmaR1"), levels = c("CTRL", "shSigmaR1")),
  replicate = coldata_df$replicat[match(sample_cols, rownames(coldata_df))],
  stringsAsFactors = FALSE
)

###############################################################################
# 2. Build the gene-level expression matrix and DE summary table
###############################################################################

# Map ENSEMBL identifiers to HGNC symbols and collapse duplicated symbols so
# that downstream interpretation operates at the gene level.
ensg_ids_vsd <- sub("\\..*", "", rownames(vsd3))
gene_symbols_vsd <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = ensg_ids_vsd,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)
SummarizedExperiment::rowData(vsd3)$gene_symbol <- gene_symbols_vsd

ensg_ids_dds <- sub("\\..*", "", rownames(dds))
gene_symbols_dds <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = ensg_ids_dds,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)
SummarizedExperiment::rowData(dds)$gene_symbol <- gene_symbols_dds

expr_raw <- assay(vsd3)
has_symbol <- !is.na(SummarizedExperiment::rowData(vsd3)$gene_symbol) & SummarizedExperiment::rowData(vsd3)$gene_symbol != ""
expr_raw <- expr_raw[has_symbol, phenotype$sample, drop = FALSE]
rownames(expr_raw) <- toupper(SummarizedExperiment::rowData(vsd3)$gene_symbol[has_symbol])
expr_raw <- collapse_matrix_by_mean(expr_raw, rownames(expr_raw))
expr_log <- expr_raw

# Reuse the DESeq2 design already embedded in `dds`, then extract classical
# statistics plus shrink log2 fold-changes for more stable ranking.
if (!"results" %in% names(metadata(dds))) {
  dds <- DESeq(dds)
}
res <- DESeq2::results(dds, contrast = c("condition", "sh11", "ctrl"))
res_shrink <- DESeq2::lfcShrink(dds, coef = "condition_sh11_vs_ctrl", type = "apeglm")

res_df <- data.frame(
  ID = rownames(res),
  gene_symbol = toupper(gene_symbols_dds[match(sub("\\..*", "", rownames(res)), ensg_ids_dds)]),
  log2FoldChange = safe_num(res_shrink$log2FoldChange),
  pvalue = safe_num(res$pvalue),
  padj = safe_num(res$padj),
  stringsAsFactors = FALSE
)
res_df <- res_df[!is.na(res_df$gene_symbol) & res_df$gene_symbol != "", , drop = FALSE]
res_df <- collapse_res_table(res_df)
res_df <- res_df[res_df$gene_symbol %in% rownames(expr_log), , drop = FALSE]
res_df$mean_CTRL <- rowMeans(expr_log[res_df$gene_symbol, ctrl_cols, drop = FALSE], na.rm = TRUE)
res_df$mean_shSigmaR1 <- rowMeans(expr_log[res_df$gene_symbol, sh_cols, drop = FALSE], na.rm = TRUE)
res_df$mean_expression <- rowMeans(expr_log[res_df$gene_symbol, sample_cols, drop = FALSE], na.rm = TRUE)
res_df$delta_sh_minus_ctrl <- res_df$mean_shSigmaR1 - res_df$mean_CTRL

###############################################################################
# 3. Define CAF marker panels and signature gene sets
###############################################################################

iCAF_top30 <- read_gene_set(top30_icaf_path)
myCAF_top30 <- read_gene_set(top30_mycaf_path)
csCAF_top30 <- read_gene_set(top30_cscaf_path)
apCAF_local_extended <- read_gene_set(apcaf_markers_path)

# Signature provenance:
# - myCAF and iCAF concepts originate from Ã–hlund et al. (J Exp Med, 2017),
#   and were further refined in human PDAC single-cell analyses by Elyada et al.
#   (Cancer Discovery, 2019).
# - The local top30 files are curated human marker lists used previously in the
#   project to define these PDAC CAF states.
# - The apCAF component is anchored in MHC-II-like antigen-presentation genes
#   described by Elyada et al., complemented here by the local project resource.
# - The csCAF list comes from the local project resource and is used here as a
#   stromal signature reference.
ifCAF_manual <- c("RSAD2", "IFIT1", "IFIT2", "IFIT3", "ISG15", "MX1", "OAS1", "OAS2", "STAT1", "STAT2",
                  "IRF7", "IRF9", "DDX58", "IFI44", "IFI44L", "CXCL10", "CXCL11")
apCAF_elyada_mhcii <- c(
  "CD74", "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1",
  "HLA-DQA1", "HLA-DQB1", "HLA-DMA", "HLA-DMB", "CIITA"
)
apCAF_chen2025_extended <- c(
  "SLPI", "MSLN", "CLU", "ENPP2", "UPK3B", "KRT19", "KRT8", "KRT18",
  "KRT7", "SPINT2", "CXADR", "SLC9A3R1", "CAV1", "PTGIS", "PDPN", "DCN", "COL1A1"
)
apCAF_core <- unique(c(apCAF_elyada_mhcii, apCAF_chen2025_extended, apCAF_local_extended))

signature_list <- list(
  myCAF_top30 = myCAF_top30,
  iCAF_top30 = iCAF_top30,
  csCAF_top30 = csCAF_top30,
  apCAF_core = apCAF_core,
  ifCAF_manual = ifCAF_manual
)
signature_list <- lapply(signature_list, function(x) intersect(unique(toupper(x)), rownames(expr_log)))

signature_reference <- data.frame(
  signature = c("myCAF_top30", "iCAF_top30", "csCAF_top30", "apCAF_core", "ifCAF_manual"),
  source = c(
    "Ohlund et al. 2017 / Elyada et al. 2019",
    "Ohlund et al. 2017 / Elyada et al. 2019",
    "Complement-secreting CAFs in PDAC single-cell multiomics",
    "Elyada et al. 2019 MHC-II apCAF anchors + Chen X, Huang H et al. 2025 Cancer Cell apCAF niche markers",
    "Curated stage interferon / ifCAF-like panel"
  ),
  interpretation_note = c(
    "Contractile / myofibroblastic axis",
    "Inflammatory fibroblast axis",
    "Complement / secretory fibroblast axis",
    "Bulk readout interpreted cautiously because only a subset of apCAF genes is detectable",
    "Work signature assembled to track interferon-like fibro-inflammatory signal"
  ),
  n_defined = c(
    length(unique(myCAF_top30)),
    length(unique(iCAF_top30)),
    length(unique(csCAF_top30)),
    length(unique(apCAF_core)),
    length(unique(ifCAF_manual))
  ),
  n_present_bulk = c(
    length(signature_list$myCAF_top30),
    length(signature_list$iCAF_top30),
    length(signature_list$csCAF_top30),
    length(signature_list$apCAF_core),
    length(signature_list$ifCAF_manual)
  ),
  stringsAsFactors = FALSE
)

core_marker_panel <- c(
  "SIGMAR1",
  "ACTA2", "TAGLN", "COL1A1", "POSTN",
  "IL6", "CXCL8", "CXCL12", "CCL2",
  "C3", "CFD", "OGN", "GPX3",
  "CD74", "HLA-DRA", "HLA-DPA1", "HLA-DPB1", "SLPI", "MSLN", "CLU",
  "RSAD2"
)
core_marker_panel <- intersect(core_marker_panel, rownames(expr_log))

marker_core_map <- data.frame(
  gene_symbol = core_marker_panel,
  category = ifelse(core_marker_panel %in% c("ACTA2", "TAGLN", "COL1A1", "POSTN"), "myCAF",
             ifelse(core_marker_panel %in% c("IL6", "CXCL8", "CXCL12", "CCL2"), "iCAF",
             ifelse(core_marker_panel %in% c("C3", "CFD", "OGN", "GPX3"), "csCAF",
             ifelse(core_marker_panel %in% c("CD74", "HLA-DRA", "HLA-DPA1", "HLA-DPB1", "SLPI", "MSLN", "CLU"), "apCAF",
             ifelse(core_marker_panel %in% c("RSAD2"), "ifCAF", "SigmaR1"))))),
  stringsAsFactors = FALSE
)
marker_signature_map <- dplyr::bind_rows(
  data.frame(gene_symbol = signature_list$myCAF_top30, category = "myCAF", stringsAsFactors = FALSE),
  data.frame(gene_symbol = signature_list$iCAF_top30, category = "iCAF", stringsAsFactors = FALSE),
  data.frame(gene_symbol = signature_list$csCAF_top30, category = "csCAF", stringsAsFactors = FALSE),
  data.frame(gene_symbol = signature_list$apCAF_core, category = "apCAF", stringsAsFactors = FALSE),
  data.frame(gene_symbol = signature_list$ifCAF_manual, category = "ifCAF", stringsAsFactors = FALSE)
)
marker_category_map <- dplyr::bind_rows(
  data.frame(gene_symbol = "SIGMAR1", category = "SigmaR1", stringsAsFactors = FALSE),
  marker_core_map,
  marker_signature_map
)
marker_category_map <- marker_category_map[!duplicated(marker_category_map$gene_symbol), , drop = FALSE]

###############################################################################
# 4. Sample-level QC moffitt_stromal_exploration_outputs
###############################################################################

plot_qc_distributions(expr_raw, expr_log, file.path(fig_dir, "QC_distributions.pdf"))
plot_qc_boxplots(expr_log, phenotype, file.path(fig_dir, "QC_boxplot_samples.pdf"))
plot_density_lines(expr_log, phenotype, file.path(fig_dir, "QC_density_samples.pdf"))
corr_mat <- stats::cor(expr_log[, phenotype$sample, drop = FALSE], method = "pearson", use = "pairwise.complete.obs")
pca_df <- plot_pca(expr_log, phenotype, file.path(fig_dir, "PCA_samples.pdf"))
pca_var_df <- plot_pca_variance(expr_log, file.path(fig_dir, "PCA_variance.pdf"))

###############################################################################
# 5. Summarize differential expression
###############################################################################

plot_volcano(res_df, file.path(fig_dir, "Volcano_plot.pdf"))
plot_ma_like(res_df, file.path(fig_dir, "MA_like_plot.pdf"))

top_up <- res_df[is.finite(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange > 0, , drop = FALSE]
top_up <- top_up[order(top_up$padj, -top_up$log2FoldChange), , drop = FALSE]
top_down <- res_df[is.finite(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange < 0, , drop = FALSE]
top_down <- top_down[order(top_down$padj, top_down$log2FoldChange), , drop = FALSE]
deg_heat_genes <- unique(c(head(top_down$gene_symbol, 25), head(top_up$gene_symbol, 25)))
deg_heat_genes <- deg_heat_genes[deg_heat_genes %in% rownames(expr_log)]
deg_heat_mat <- t(apply(expr_log[deg_heat_genes, sample_cols, drop = FALSE], 1, zscore_vec))
deg_row_cols <- ifelse(res_df$log2FoldChange[match(rownames(deg_heat_mat), res_df$gene_symbol)] > 0, sh_color, ctrl_color)
deg_col_cols <- ifelse(phenotype$condition == "CTRL", ctrl_color, sh_color)

###############################################################################
# 6. Compute marker-level statistics
###############################################################################

all_marker_genes <- unique(c(core_marker_panel, unlist(signature_list)))
all_marker_genes <- intersect(all_marker_genes, rownames(expr_log))
sigmar1_expr <- if ("SIGMAR1" %in% rownames(expr_log)) expr_log["SIGMAR1", phenotype$sample] else rep(NA_real_, nrow(phenotype))

# At the marker level, the goal is to quantify both differential behavior
# between conditions and monotonic association with the perturbed SigmaR1 gene.
marker_stats <- do.call(rbind, lapply(all_marker_genes, function(g) {
  expr_vec <- expr_log[g, phenotype$sample]
  ctrl_vals <- expr_vec[phenotype$condition == "CTRL"]
  sh_vals <- expr_vec[phenotype$condition == "shSigmaR1"]
  wt <- safe_wilcox_p(expr_vec, phenotype$condition)
  cor_all <- safe_spearman(sigmar1_expr, expr_vec)
  cor_ctrl <- safe_spearman(sigmar1_expr[phenotype$condition == "CTRL"], ctrl_vals)
  cor_sh <- safe_spearman(sigmar1_expr[phenotype$condition == "shSigmaR1"], sh_vals)
  data.frame(
    gene_symbol = g,
    category = marker_category_map$category[match(g, marker_category_map$gene_symbol)],
    mean_CTRL = mean(ctrl_vals, na.rm = TRUE),
    mean_shSigmaR1 = mean(sh_vals, na.rm = TRUE),
    median_CTRL = stats::median(ctrl_vals, na.rm = TRUE),
    median_shSigmaR1 = stats::median(sh_vals, na.rm = TRUE),
    delta_sh_minus_ctrl = mean(sh_vals, na.rm = TRUE) - mean(ctrl_vals, na.rm = TRUE),
    wilcox_p = wt,
    rank_biserial = rank_biserial(sh_vals, ctrl_vals),
    log2FoldChange = res_df$log2FoldChange[match(g, res_df$gene_symbol)],
    pvalue = res_df$pvalue[match(g, res_df$gene_symbol)],
    padj = res_df$padj[match(g, res_df$gene_symbol)],
    rho_all = cor_all["rho"],
    rho_pvalue = cor_all["pvalue"],
    rho_n = cor_all["n"],
    rho_CTRL = cor_ctrl["rho"],
    rho_CTRL_pvalue = cor_ctrl["pvalue"],
    rho_shSigmaR1 = cor_sh["rho"],
    rho_shSigmaR1_pvalue = cor_sh["pvalue"],
    stringsAsFactors = FALSE
  )
}))

marker_stats$wilcox_padj <- stats::p.adjust(marker_stats$wilcox_p, method = "BH")
marker_stats$rho_padj <- stats::p.adjust(marker_stats$rho_pvalue, method = "BH")
marker_stats <- marker_stats[order(marker_stats$category, marker_stats$padj, -abs(marker_stats$rho_all)), , drop = FALSE]

###############################################################################
# 7. Generate marker-level visualizations retained in the cleaned output
###############################################################################

plot_marker_violin_panels(expr_log, phenotype, core_marker_panel, marker_stats,
                          file.path(fig_dir, "CAF_marker_violin_panels.pdf"))

dot_genes <- unique(c(
  head(signature_list$myCAF_top30, 8),
  head(signature_list$iCAF_top30, 8),
  head(signature_list$csCAF_top30, 8),
  head(signature_list$apCAF_core, 8),
  head(signature_list$ifCAF_manual, 8)
))
dot_genes <- intersect(dot_genes, rownames(expr_log))
dot_df <- do.call(rbind, lapply(dot_genes, function(g) {
  gene_vals <- expr_log[g, phenotype$sample]
  gene_z <- zscore_vec(gene_vals)
  gene_med <- stats::median(gene_vals, na.rm = TRUE)
  do.call(rbind, lapply(levels(phenotype$condition), function(cond) {
    idx <- phenotype$condition == cond
    data.frame(
      gene = g,
      condition = cond,
      mean_z = mean(gene_z[idx], na.rm = TRUE),
      frac_high = mean(gene_vals[idx] > gene_med, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
}))

heatmap_genes <- unique(c(
  head(signature_list$myCAF_top30, 10),
  head(signature_list$iCAF_top30, 10),
  head(signature_list$csCAF_top30, 10),
  head(signature_list$apCAF_core, 10),
  head(signature_list$ifCAF_manual, 10)
))
heatmap_genes <- intersect(heatmap_genes, rownames(expr_log))
heatmap_order_df <- marker_stats[match(heatmap_genes, marker_stats$gene_symbol), , drop = FALSE]
heatmap_order_df <- heatmap_order_df[order(heatmap_order_df$category, heatmap_order_df$delta_sh_minus_ctrl), , drop = FALSE]
heatmap_genes <- heatmap_order_df$gene_symbol
heatmap_mat <- t(apply(expr_log[heatmap_genes, sample_cols, drop = FALSE], 1, zscore_vec))
row_side_cols <- signature_colors[heatmap_order_df$category]
col_side_cols <- ifelse(phenotype$condition == "CTRL", ctrl_color, sh_color)

###############################################################################
# 8. Score CAF-related signatures across samples
###############################################################################

score_mat <- cbind(
  myCAF_top30 = score_signature_z(expr_log, signature_list$myCAF_top30),
  iCAF_top30 = score_signature_z(expr_log, signature_list$iCAF_top30),
  csCAF_top30 = score_signature_z(expr_log, signature_list$csCAF_top30),
  apCAF_core = score_signature_z(expr_log, signature_list$apCAF_core),
  ifCAF_manual = score_signature_z(expr_log, signature_list$ifCAF_manual)
)
rownames(score_mat) <- phenotype$sample

signature_stats <- do.call(rbind, lapply(colnames(score_mat), function(sig) {
  x <- score_mat[, sig]
  # Wilcoxon testing is appropriate here because the number of replicates is
  # limited and we do not want to rely on strong normality assumptions.
  data.frame(
    signature = sig,
    mean_CTRL = mean(x[phenotype$condition == "CTRL"], na.rm = TRUE),
    mean_shSigmaR1 = mean(x[phenotype$condition == "shSigmaR1"], na.rm = TRUE),
    median_CTRL = stats::median(x[phenotype$condition == "CTRL"], na.rm = TRUE),
    median_shSigmaR1 = stats::median(x[phenotype$condition == "shSigmaR1"], na.rm = TRUE),
    wilcox_p = safe_wilcox_p(x, phenotype$condition),
    stringsAsFactors = FALSE
  )
}))
signature_stats$wilcox_padj <- stats::p.adjust(signature_stats$wilcox_p, method = "BH")
plot_signature_score_panels(score_mat, phenotype, signature_stats,
                            file.path(fig_dir, "CAF_signature_score_panels.pdf"))

score_condition_mat <- rbind(
  CTRL = colMeans(score_mat[phenotype$condition == "CTRL", , drop = FALSE], na.rm = TRUE),
  shSigmaR1 = colMeans(score_mat[phenotype$condition == "shSigmaR1", , drop = FALSE], na.rm = TRUE)
)
contrib_signatures <- list(
  myCAF = signature_list$myCAF_top30,
  iCAF = signature_list$iCAF_top30,
  csCAF = signature_list$csCAF_top30,
  apCAF = signature_list$apCAF_core,
  ifCAF = signature_list$ifCAF_manual
)
contrib_genes <- unique(unlist(contrib_signatures))
contrib_genes <- contrib_genes[contrib_genes %in% rownames(expr_log)]
contrib_mat <- expr_log[contrib_genes, phenotype$sample, drop = FALSE]
ctrl_mean <- rowMeans(contrib_mat[, phenotype$condition == "CTRL", drop = FALSE], na.rm = TRUE)
sh_mean <- rowMeans(contrib_mat[, phenotype$condition == "shSigmaR1", drop = FALSE], na.rm = TRUE)
contrib_df <- data.frame(
  gene = names(ctrl_mean),
  mean_CTRL = ctrl_mean,
  mean_shSigmaR1 = sh_mean,
  delta = sh_mean - ctrl_mean,
  signature = NA_character_,
  stringsAsFactors = FALSE
)
for (sig in names(contrib_signatures)) {
  contrib_df$signature[contrib_df$gene %in% contrib_signatures[[sig]]] <- sig
}
contrib_df <- contrib_df[!is.na(contrib_df$signature), , drop = FALSE]
contrib_df$signature <- factor(contrib_df$signature, levels = c("myCAF", "iCAF", "csCAF", "apCAF", "ifCAF"))
contrib_df <- contrib_df[order(contrib_df$signature, contrib_df$delta), , drop = FALSE]
plot_gene_contribution_barplot(contrib_df, file.path(fig_dir, "CAF_gene_contribution_barplot.pdf"))

###############################################################################
# 9. Summarize SigmaR1-centered marker correlations
###############################################################################

sigmar1_corr_focus <- marker_stats[marker_stats$gene_symbol != "SIGMAR1", , drop = FALSE]
sigmar1_corr_focus <- sigmar1_corr_focus[order(sigmar1_corr_focus$rho_padj, -abs(sigmar1_corr_focus$rho_all)), , drop = FALSE]

summary_matrix_genes <- unique(c(
  head(sigmar1_corr_focus$gene_symbol[sigmar1_corr_focus$rho_all > 0], 8),
  head(sigmar1_corr_focus$gene_symbol[sigmar1_corr_focus$rho_all < 0], 8)
))
summary_matrix_genes <- summary_matrix_genes[summary_matrix_genes %in% sigmar1_corr_focus$gene_symbol]
summary_matrix <- cbind(
  rho_all = sigmar1_corr_focus$rho_all[match(summary_matrix_genes, sigmar1_corr_focus$gene_symbol)],
  delta_sh_minus_ctrl = sigmar1_corr_focus$delta_sh_minus_ctrl[match(summary_matrix_genes, sigmar1_corr_focus$gene_symbol)],
  log2FoldChange = sigmar1_corr_focus$log2FoldChange[match(summary_matrix_genes, sigmar1_corr_focus$gene_symbol)]
)
rownames(summary_matrix) <- summary_matrix_genes

scatter_genes <- unique(c(
  head(sigmar1_corr_focus$gene_symbol[sigmar1_corr_focus$category == "myCAF"], 2),
  head(sigmar1_corr_focus$gene_symbol[sigmar1_corr_focus$category == "iCAF"], 1),
  head(sigmar1_corr_focus$gene_symbol[sigmar1_corr_focus$category == "csCAF"], 1),
  head(sigmar1_corr_focus$gene_symbol[sigmar1_corr_focus$category == "apCAF"], 1),
  head(sigmar1_corr_focus$gene_symbol[sigmar1_corr_focus$category == "ifCAF"], 1)
))
scatter_df <- sigmar1_corr_focus[sigmar1_corr_focus$gene_symbol %in% scatter_genes, , drop = FALSE]

###############################################################################
# 10. Figures and Hallmark enrichment analysis
###############################################################################

expr_long <- data.frame(
  gene_symbol = rep(rownames(expr_log), times = ncol(expr_log)),
  sample = rep(colnames(expr_log), each = nrow(expr_log)),
  expression = as.vector(expr_log),
  stringsAsFactors = FALSE
)
expr_long$condition <- phenotype$condition[match(expr_long$sample, phenotype$sample)]

sample_summary_df <- sample_summary_stats(expr_log)
sample_summary_df$condition <- phenotype$condition[match(sample_summary_df$sample, phenotype$sample)]

annotation_col <- data.frame(Condition = phenotype$condition, row.names = phenotype$sample, stringsAsFactors = FALSE)

plot_qc_boxplots(expr_log, phenotype, file.path(fig_dir, "QC_boxplot_samples.pdf"))
plot_density_lines(expr_log, phenotype, file.path(fig_dir, "QC_density_samples.pdf"))

pca_df <- plot_pca(expr_log[, phenotype$sample, drop = FALSE], phenotype, file.path(fig_dir, "PCA_samples.pdf"))
tmp_pca <- stats::prcomp(t(expr_log[, phenotype$sample, drop = FALSE]), scale. = TRUE)
pca_var_df <- data.frame(
  PC = paste0("PC", seq_along(tmp_pca$sdev)),
  variance = tmp_pca$sdev^2 / sum(tmp_pca$sdev^2),
  stringsAsFactors = FALSE
)
plot_pca_variance(expr_log[, phenotype$sample, drop = FALSE], file.path(fig_dir, "PCA_variance.pdf"))

plot_marker_violin_panels(expr_log, phenotype, core_marker_panel, marker_stats,
                             file.path(fig_dir, "CAF_marker_violin_panels.pdf"))
plot_signature_score_panels(score_mat, phenotype, signature_stats,
                               file.path(fig_dir, "CAF_signature_score_panels.pdf"))

signature_blocks <- list(
  SigmaR1 = intersect(c("SIGMAR1"), rownames(expr_log)),
  myCAF = signature_list$myCAF_top30,
  iCAF = signature_list$iCAF_top30,
  csCAF = signature_list$csCAF_top30,
  apCAF = signature_list$apCAF_core,
  ifCAF = signature_list$ifCAF_manual
)

row_orders <- lapply(names(signature_blocks), function(cat_name) {
  genes <- unique(signature_blocks[[cat_name]])
  genes <- genes[genes %in% rownames(expr_log)]
  if (length(genes) == 0) return(character(0))
  submat <- t(apply(expr_log[genes, phenotype$sample, drop = FALSE], 1, zscore_vec))
  cluster_order(submat)
})
names(row_orders) <- names(signature_blocks)
row_order_all <- unlist(row_orders, use.names = FALSE)

row_meta <- data.frame(
  gene_symbol = row_order_all,
  Category = rep(names(row_orders), times = vapply(row_orders, length, integer(1))),
  stringsAsFactors = FALSE
)
row_meta$DE_padj <- marker_stats$padj[match(row_meta$gene_symbol, marker_stats$gene_symbol)]
row_meta$W_p <- marker_stats$wilcox_p[match(row_meta$gene_symbol, marker_stats$gene_symbol)]
row_meta$W_padj <- marker_stats$wilcox_padj[match(row_meta$gene_symbol, marker_stats$gene_symbol)]
row_meta$Significance <- mapply(sig_label, row_meta$DE_padj, row_meta$W_p, row_meta$W_padj, USE.NAMES = FALSE)
row_meta$label <- make.unique(row_meta$gene_symbol, sep = " | ")

annotation_col_block <- data.frame(
  Condition = factor(phenotype$condition, levels = c("CTRL", "shSigmaR1")),
  row.names = phenotype$sample,
  stringsAsFactors = FALSE
)
ann_colors_block <- list(
  Condition = c(CTRL = ctrl_color, shSigmaR1 = sh_color)
)

if (length(row_order_all) > 0) {
  all_heatmap_mat <- t(apply(expr_log[row_order_all, phenotype$sample, drop = FALSE], 1, zscore_vec))
  rownames(all_heatmap_mat) <- row_meta$label
  annotation_row_all <- data.frame(
    Category = factor(row_meta$Category, levels = names(signature_blocks)),
    Significance = factor(row_meta$Significance, levels = names(sig_colors)),
    row.names = row_meta$label,
    stringsAsFactors = FALSE
  )
  gaps_row <- cumsum(vapply(row_orders, length, integer(1)))
  gaps_row <- gaps_row[gaps_row < nrow(all_heatmap_mat)]
}

combined_blocks <- c("SigmaR1", "myCAF", "iCAF", "csCAF", "ifCAF")
row_order_combined <- unlist(row_orders[combined_blocks], use.names = FALSE)
if (length(row_order_combined) > 0) {
  combined_mat <- t(apply(expr_log[row_order_combined, phenotype$sample, drop = FALSE], 1, zscore_vec))
  combined_meta <- row_meta[row_meta$gene_symbol %in% row_order_combined, , drop = FALSE]
  combined_meta <- combined_meta[match(row_order_combined, combined_meta$gene_symbol), , drop = FALSE]
  rownames(combined_mat) <- combined_meta$label
  block_lengths <- vapply(row_orders[combined_blocks], length, integer(1))
  keep_blocks <- block_lengths > 0
  block_lengths <- block_lengths[keep_blocks]
  block_labels <- names(block_lengths)
  gaps_row_combined <- cumsum(block_lengths)
  gaps_row_combined <- gaps_row_combined[gaps_row_combined < nrow(combined_mat)]
  annotation_col_combined <- data.frame(
    Condition = factor(phenotype$condition, levels = c("CTRL", "shSigmaR1")),
    row.names = phenotype$sample,
    stringsAsFactors = FALSE
  )
  save_portrait_block_heatmap(
    display_mat = combined_mat,
    annotation_col = annotation_col_combined,
    annotation_colors = ann_colors_block,
    gaps_row = gaps_row_combined,
    block_labels = block_labels,
    block_lengths = block_lengths,
    file_path = file.path(fig_dir, "CAF_signature_heatmap.pdf"),
    width = 7.5,
    height = 11.5,
    title_text = "CAF signature heatmap"
  )
  save_portrait_block_heatmap(
    display_mat = combined_mat,
    annotation_col = annotation_col_combined,
    annotation_colors = ann_colors_block,
    gaps_row = gaps_row_combined,
    block_labels = block_labels,
    block_lengths = block_lengths,
    file_path = file.path(fig_dir, "CAF_signature_heatmap.png"),
    width = 7.5,
    height = 11.5,
    title_text = "CAF signature heatmap"
  )
}

if (length(summary_matrix_genes) >= 3) {
  kcnma1_pca_genes <- unique(c("SIGMAR1", head(summary_matrix_genes, 12)))
  kcnma1_pca_genes <- kcnma1_pca_genes[kcnma1_pca_genes %in% rownames(expr_log)]
  pca_kcnma1_df <- plot_pca(
    expr_log[kcnma1_pca_genes, phenotype$sample, drop = FALSE],
    phenotype,
    file.path(fig_dir, "PCA_SIGMAR1_module.pdf")
  )
} else {
  pca_kcnma1_df <- NULL
}

ranks <- res_df$log2FoldChange
names(ranks) <- res_df$gene_symbol
ranks <- tapply(ranks, names(ranks), mean)
ranks <- sort(unlist(ranks), decreasing = TRUE)

# FGSEA complements single-gene DE by highlighting whether coordinated
# Hallmark programs shift after SigmaR1 perturbation.
hallmark_tbl <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
hallmark_pathways <- split(toupper(hallmark_tbl$gene_symbol), hallmark_tbl$gs_name)
hallmark_pathways <- lapply(hallmark_pathways, function(x) intersect(unique(x), names(ranks)))
hallmark_pathways <- hallmark_pathways[sapply(hallmark_pathways, length) >= 10]
gsea_res <- fgsea::fgseaMultilevel(pathways = hallmark_pathways, stats = ranks, eps = 0) %>%
  as.data.frame() %>%
  arrange(padj, desc(abs(NES)))
plot_gsea_bubble(gsea_res, file.path(fig_dir, "GSEA_Hallmark_bubbleplot.pdf"))

###############################################################################
# 11. Export tables and serialized objects
###############################################################################

expr_summary <- data.frame(
  metric = c("n_input_rows", "n_gene_symbols", "n_collapsed_genes", "n_samples", "n_ctrl", "n_shSigmaR1",
             "n_deg_padj_lt_0_05", "n_up_shSigmaR1", "n_down_shSigmaR1"),
  value = c(
    nrow(dds),
    length(unique(res_df$gene_symbol)),
    nrow(expr_log),
    ncol(expr_log),
    length(ctrl_cols),
    length(sh_cols),
    sum(res_df$padj < 0.05, na.rm = TRUE),
    sum(res_df$padj < 0.05 & res_df$log2FoldChange > 0, na.rm = TRUE),
    sum(res_df$padj < 0.05 & res_df$log2FoldChange < 0, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)

sample_summary_df <- sample_summary_stats(expr_log)
sample_summary_df$condition <- phenotype$condition[match(sample_summary_df$sample, phenotype$sample)]
sample_summary_df$SIGMAR1 <- if ("SIGMAR1" %in% rownames(expr_log)) as.numeric(expr_log["SIGMAR1", sample_summary_df$sample]) else NA_real_

phenotype_out <- phenotype
for (sig in colnames(score_mat)) phenotype_out[[sig]] <- score_mat[, sig]

marker_presence <- do.call(rbind, lapply(names(signature_list), function(sig) {
  data.frame(
    signature = sig,
    n_defined = signature_reference$n_defined[match(sig, signature_reference$signature)],
    n_present = length(signature_list[[sig]]),
    present_genes = paste(signature_list[[sig]], collapse = ";"),
    stringsAsFactors = FALSE
  )
}))
row_export <- if (exists("row_meta")) {
  row_meta[, c("gene_symbol", "Category", "Significance", "DE_padj", "W_p", "W_padj", "label")]
} else {
  data.frame()
}

utils::write.csv(expr_summary, file.path(tab_dir, "expression_summary.csv"), row.names = FALSE, na = "")
utils::write.csv(phenotype_out, file.path(tab_dir, "sample_metadata_and_scores.csv"), row.names = FALSE, na = "")
utils::write.csv(sample_summary_df, file.path(tab_dir, "sample_qc_summary.csv"), row.names = FALSE, na = "")
utils::write.csv(data.frame(sample1 = rep(colnames(corr_mat), each = ncol(corr_mat)),
                            sample2 = rep(colnames(corr_mat), times = ncol(corr_mat)),
                            correlation = as.vector(corr_mat),
                            stringsAsFactors = FALSE),
                 file.path(tab_dir, "sample_correlation_matrix_long.csv"), row.names = FALSE, na = "")
utils::write.csv(pca_df, file.path(tab_dir, "PCA_coordinates.csv"), row.names = FALSE, na = "")
utils::write.csv(pca_var_df, file.path(tab_dir, "PCA_variance.csv"), row.names = FALSE, na = "")
utils::write.csv(res_df[order(res_df$padj, -abs(res_df$log2FoldChange)), ],
                 file.path(tab_dir, "DESeq2_results_gene_level.csv"), row.names = FALSE, na = "")
utils::write.csv(marker_stats, file.path(tab_dir, "CAF_marker_statistics.csv"), row.names = FALSE, na = "")
utils::write.csv(signature_stats, file.path(tab_dir, "CAF_signature_statistics.csv"), row.names = FALSE, na = "")
utils::write.csv(marker_presence, file.path(tab_dir, "CAF_signature_presence.csv"), row.names = FALSE, na = "")
utils::write.csv(signature_reference, file.path(tab_dir, "CAF_signature_reference.csv"), row.names = FALSE, na = "")
utils::write.csv(row_export, file.path(tab_dir, "CAF_signature_heatmap_row_annotations.csv"), row.names = FALSE, na = "")
utils::write.csv(contrib_df, file.path(tab_dir, "CAF_gene_contribution.csv"), row.names = FALSE, na = "")
utils::write.csv(sigmar1_corr_focus, file.path(tab_dir, "SIGMAR1_marker_correlations.csv"), row.names = FALSE, na = "")
gsea_res_out <- gsea_res
if ("leadingEdge" %in% colnames(gsea_res_out)) {
  gsea_res_out$leadingEdge <- vapply(gsea_res_out$leadingEdge, function(x) paste(x, collapse = ";"), character(1))
}
utils::write.csv(gsea_res_out, file.path(tab_dir, "Hallmark_fgsea_results.csv"), row.names = FALSE, na = "")
if (exists("pca_kcnma1_df") && !is.null(pca_kcnma1_df)) {
  utils::write.csv(pca_kcnma1_df, file.path(tab_dir, "PCA_SIGMAR1_module_coordinates.csv"), row.names = FALSE, na = "")
}

saveRDS(
  list(
    expr_log = expr_log,
    phenotype = phenotype_out,
    res_df = res_df,
    marker_stats = marker_stats,
    signature_stats = signature_stats,
    gsea_res = gsea_res_out,
    pca_sigmar1 = if (exists("pca_kcnma1_df")) pca_kcnma1_df else NULL
  ),
  file.path(rds_dir, "sigmar1_shsigmar1_bulk_analysis_bundle.rds")
)


###############################################################################
# 12. Write the README for the output folder
###############################################################################

best_pos <- marker_stats[marker_stats$gene_symbol != "SIGMAR1" & marker_stats$rho_all > 0, , drop = FALSE]
best_neg <- marker_stats[marker_stats$gene_symbol != "SIGMAR1" & marker_stats$rho_all < 0, , drop = FALSE]
best_pos <- best_pos[order(best_pos$rho_padj, -best_pos$rho_all), , drop = FALSE]
best_neg <- best_neg[order(best_neg$rho_padj, best_neg$rho_all), , drop = FALSE]
best_sig <- signature_stats[order(signature_stats$wilcox_padj, -abs(signature_stats$mean_shSigmaR1 - signature_stats$mean_CTRL)), , drop = FALSE]
best_gsea <- gsea_res[order(gsea_res$padj, -abs(gsea_res$NES)), , drop = FALSE]

readme_lines <- c(
  "# SigmaR1 / shSigmaR1 bulk RNA-seq analysis",
  "",
  "## Summary",
  "This folder contains a cleaned and reproducible downstream analysis of the SigmaR1 perturbation bulk RNA-seq dataset, centered on essential QC, CAF-oriented transcriptional programs, and the combined CAF heatmap without apCAF.",
  "",
  "## Key points",
  paste0("- collapsed genes: ", nrow(expr_log)),
  paste0("- samples: ", ncol(expr_log), " (", length(ctrl_cols), " CTRL, ", length(sh_cols), " shSigmaR1)"),
  paste0("- DEG `padj < 0.05`: ", sum(res_df$padj < 0.05, na.rm = TRUE)),
  paste0("- top Wilcoxon signature: `", best_sig$signature[1], "` (p=", fmt_p(best_sig$wilcox_p[1]), ", FDR=", fmt_p(best_sig$wilcox_padj[1]), ")"),
  paste0("- strongest positive `SIGMAR1` correlation: `", best_pos$gene_symbol[1], "` (rho=", round(best_pos$rho_all[1], 3), ", FDR=", fmt_p(best_pos$rho_padj[1]), ")"),
  paste0("- strongest negative `SIGMAR1` correlation: `", best_neg$gene_symbol[1], "` (rho=", round(best_neg$rho_all[1], 3), ", FDR=", fmt_p(best_neg$rho_padj[1]), ")"),
  paste0("- top Hallmark pathway: `", best_gsea$pathway[1], "` (NES=", round(best_gsea$NES[1], 3), ", FDR=", fmt_p(best_gsea$padj[1]), ")"),
  "",
  "## Main figures",
  "- `figures/PCA_samples.pdf`",
  "- `figures/PCA_SIGMAR1_module.pdf`",
  "- `figures/Volcano_plot.pdf`",
  "- `figures/CAF_marker_violin_panels.pdf`",
  "- `figures/CAF_gene_contribution_barplot.pdf`",
  "- `figures/CAF_signature_heatmap.pdf`",
  "- `figures/CAF_signature_score_panels.pdf`",
  "- `figures/GSEA_Hallmark_bubbleplot.pdf`",
  "",
  "## Main tables",
  "- `tables/expression_summary.csv`",
  "- `tables/DESeq2_results_gene_level.csv`",
  "- `tables/CAF_marker_statistics.csv`",
  "- `tables/CAF_signature_statistics.csv`",
  "- `tables/CAF_signature_reference.csv`",
  "- `tables/CAF_gene_contribution.csv`",
  "- `tables/SIGMAR1_marker_correlations.csv`",
  "- `tables/Hallmark_fgsea_results.csv`",
  "- `rds/sigmar1_shsigmar1_bulk_analysis_bundle.rds`",
  "",
  "Last updated: 2026-05-04"
)

writeLines(readme_lines, file.path(output_dir, "README.md"))
message0("SigmaR1 bulk RNA-seq analysis complete.")

