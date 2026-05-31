#!/usr/bin/env Rscript

###############################################################################
# Multicohort Kaplan-Meier plots for the 25 KCN genes
#
# For each gene:
# - split samples within each cohort at the cohort-specific median
# - draw a 2x2 panel with one Kaplan-Meier curve per cohort
# - report sample counts, events, median OS, and log-rank p-value
###############################################################################

suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
})

options(stringsAsFactors = FALSE)

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  file_path <- sub(file_arg, "", args[grep(file_arg, args)])
  if (length(file_path) == 0) {
    return(normalizePath(getwd(), winslash = "/"))
  }
  normalizePath(dirname(file_path), winslash = "/")
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
  if (!is.finite(x)) {
    return("NA")
  }
  if (x < 1e-3) {
    return(formatC(x, format = "e", digits = 2))
  }
  sprintf("%.3f", x)
}

extract_group_medians <- function(fit_km) {
  km_tab <- summary(fit_km)$table
  out <- c(median_high = NA_real_, median_low = NA_real_)
  if (is.matrix(km_tab)) {
    if ("group=High" %in% rownames(km_tab)) out["median_high"] <- km_tab["group=High", "median"]
    if ("group=Low" %in% rownames(km_tab)) out["median_low"] <- km_tab["group=Low", "median"]
  }
  out
}

compute_logrank <- function(df) {
  fit <- tryCatch(
    survdiff(Surv(os_months, event) ~ group, data = df),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    return(list(chisq = NA_real_, pvalue = NA_real_))
  }

  pvalue <- 1 - pchisq(fit$chisq, df = length(fit$n) - 1)
  list(chisq = unname(fit$chisq), pvalue = unname(pvalue))
}

make_panel_title <- function(display_name, n_samples) {
  paste0(display_name, " (n=", n_samples, ")")
}

script_dir <- get_script_dir()
bundle_path <- file.path(script_dir, "integrated_survival_bundle", "fused_kcn_survival_multicohort_bundle.rds")
kcn_path <- file.path(script_dir, "..", "moffitt_kcn_deg_4_methods", "tables", "KCN_union_membership.tsv")
out_dir <- file.path(script_dir, "survival_kcn25_km_within_cohort_median")
tab_dir <- file.path(out_dir, "tables")
fig_dir <- file.path(out_dir, "figures")

dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

bundle <- readRDS(bundle_path)
kcn_table <- read.delim(kcn_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
kcn_genes <- unique(kcn_table$gene)

cohort_order <- bundle$cohort_order
display_names <- vapply(
  cohort_order,
  function(x) {
    current <- bundle$original_bundles[[x]]
    if (is.null(current) || is.null(current$display_name)) {
      return(x)
    }
    current$display_name
  },
  character(1)
)
names(display_names) <- cohort_order

ph <- bundle$phenotype_harmonized %>%
  mutate(
    cohort = as.character(cohort),
    sample_id = as.character(sample_id),
    os_months = as.numeric(os_months),
    event = as.numeric(event)
  )

km_summary_rows <- list()
km_status_rows <- list()
row_idx <- 1L
status_idx <- 1L

high_col <- "#800E13"
low_col <- "#2D6A9F"

for (gene in kcn_genes) {
  pdf(
    file = file.path(fig_dir, paste0(gene, "_km_within_cohort_median.pdf")),
    width = 11,
    height = 8.5,
    bg = "white"
  )

  par(mfrow = c(2, 2), mar = c(4.2, 4.2, 3.2, 1.2), oma = c(0, 0, 2, 0))

  for (cohort_name in cohort_order) {
    expr_mat <- bundle$z_common_by_cohort[[cohort_name]]

    if (is.null(expr_mat) || !gene %in% rownames(expr_mat)) {
      plot.new()
      title(main = make_panel_title(display_names[[cohort_name]], 0))
      text(0.5, 0.5, "Gene not available", cex = 1)

      km_summary_rows[[row_idx]] <- data.frame(
        gene = gene,
        cohort = cohort_name,
        display_name = display_names[[cohort_name]],
        cutoff_type = "within_cohort_median",
        cutoff_value = NA_real_,
        n_total = 0L,
        n_high = NA_integer_,
        n_low = NA_integer_,
        events_high = NA_integer_,
        events_low = NA_integer_,
        median_os_high = NA_real_,
        median_os_low = NA_real_,
        logrank_p = NA_real_,
        logrank_chisq = NA_real_,
        status = "gene_not_available",
        stringsAsFactors = FALSE
      )
      row_idx <- row_idx + 1L
      next
    }

    ph_cohort <- ph %>% filter(cohort == cohort_name)
    ph_cohort <- ph_cohort[match(colnames(expr_mat), ph_cohort$sample_id), , drop = FALSE]

    if (!identical(colnames(expr_mat), ph_cohort$sample_id)) {
      stop("Sample order mismatch for cohort: ", cohort_name)
    }

    df <- ph_cohort %>%
      transmute(
        cohort = cohort,
        sample_id = sample_id,
        os_months = os_months,
        event = event,
        expr = as.numeric(expr_mat[gene, ])
      ) %>%
      filter(is.finite(expr), is.finite(os_months), os_months > 0, is.finite(event))

    n_total <- nrow(df)
    panel_title <- make_panel_title(display_names[[cohort_name]], n_total)

    if (n_total < 4 || stats::sd(df$expr, na.rm = TRUE) == 0) {
      plot.new()
      title(main = panel_title)
      text(0.5, 0.5, "Insufficient expression variation", cex = 1)

      km_summary_rows[[row_idx]] <- data.frame(
        gene = gene,
        cohort = cohort_name,
        display_name = display_names[[cohort_name]],
        cutoff_type = "within_cohort_median",
        cutoff_value = NA_real_,
        n_total = n_total,
        n_high = NA_integer_,
        n_low = NA_integer_,
        events_high = NA_integer_,
        events_low = NA_integer_,
        median_os_high = NA_real_,
        median_os_low = NA_real_,
        logrank_p = NA_real_,
        logrank_chisq = NA_real_,
        status = "insufficient_variation",
        stringsAsFactors = FALSE
      )
      row_idx <- row_idx + 1L
      next
    }

    cutoff_value <- stats::median(df$expr, na.rm = TRUE)
    df$group <- ifelse(df$expr >= cutoff_value, "High", "Low")
    df$group <- factor(df$group, levels = c("Low", "High"))

    n_low <- sum(df$group == "Low")
    n_high <- sum(df$group == "High")
    events_low <- sum(df$event[df$group == "Low"] == 1, na.rm = TRUE)
    events_high <- sum(df$event[df$group == "High"] == 1, na.rm = TRUE)

    if (n_low == 0 || n_high == 0) {
      plot.new()
      title(main = panel_title)
      text(0.5, 0.5, "Median split not informative", cex = 1)

      km_summary_rows[[row_idx]] <- data.frame(
        gene = gene,
        cohort = cohort_name,
        display_name = display_names[[cohort_name]],
        cutoff_type = "within_cohort_median",
        cutoff_value = cutoff_value,
        n_total = n_total,
        n_high = n_high,
        n_low = n_low,
        events_high = events_high,
        events_low = events_low,
        median_os_high = NA_real_,
        median_os_low = NA_real_,
        logrank_p = NA_real_,
        logrank_chisq = NA_real_,
        status = "uninformative_split",
        stringsAsFactors = FALSE
      )
      row_idx <- row_idx + 1L
      next
    }

    fit_km <- survfit(Surv(os_months, event) ~ group, data = df)
    logrank <- compute_logrank(df)
    medians <- extract_group_medians(fit_km)

    plot(
      fit_km,
      col = c(low_col, high_col),
      lwd = 2,
      mark.time = TRUE,
      xlab = "Months",
      ylab = "Overall survival probability",
      main = panel_title,
      cex.main = 1.02,
      cex.lab = 0.9,
      cex.axis = 0.82
    )

    legend(
      "topright",
      legend = c(
        paste0("Low (n=", n_low, ", events=", events_low, ")"),
        paste0("High (n=", n_high, ", events=", events_high, ")")
      ),
      col = c(low_col, high_col),
      lwd = 2,
      bty = "n",
      cex = 0.75
    )

    legend(
      "bottomleft",
      legend = c(
        paste0("Median cutoff z = ", sprintf("%.3f", cutoff_value)),
        paste0("Log-rank p = ", fmt_p(logrank$pvalue)),
        paste0("Median OS High = ", ifelse(is.finite(medians["median_high"]), sprintf("%.1f", medians["median_high"]), "NA")),
        paste0("Median OS Low = ", ifelse(is.finite(medians["median_low"]), sprintf("%.1f", medians["median_low"]), "NA"))
      ),
      bty = "n",
      cex = 0.72
    )

    km_summary_rows[[row_idx]] <- data.frame(
      gene = gene,
      cohort = cohort_name,
      display_name = display_names[[cohort_name]],
      cutoff_type = "within_cohort_median",
      cutoff_value = cutoff_value,
      n_total = n_total,
      n_high = n_high,
      n_low = n_low,
      events_high = events_high,
      events_low = events_low,
      median_os_high = unname(medians["median_high"]),
      median_os_low = unname(medians["median_low"]),
      logrank_p = logrank$pvalue,
      logrank_chisq = logrank$chisq,
      status = "ok",
      stringsAsFactors = FALSE
    )
    row_idx <- row_idx + 1L
  }

  mtext(
    text = paste0(gene, " | Kaplan-Meier with within-cohort median split"),
    side = 3,
    outer = TRUE,
    line = 0.5,
    cex = 1.25,
    font = 2
  )

  dev.off()

  gene_ok <- do.call(rbind, km_summary_rows)[do.call(rbind, km_summary_rows)$gene == gene, , drop = FALSE]
  gene_ok <- gene_ok[gene_ok$status == "ok", , drop = FALSE]

  if (nrow(gene_ok) > 0) {
    km_status_rows[[status_idx]] <- data.frame(
      gene = gene,
      best_cohort = gene_ok$cohort[which.min(gene_ok$logrank_p)],
      best_display_name = gene_ok$display_name[which.min(gene_ok$logrank_p)],
      best_logrank_p = min(gene_ok$logrank_p, na.rm = TRUE),
      n_ok_cohorts = nrow(gene_ok),
      stringsAsFactors = FALSE
    )
  } else {
    km_status_rows[[status_idx]] <- data.frame(
      gene = gene,
      best_cohort = NA_character_,
      best_display_name = NA_character_,
      best_logrank_p = NA_real_,
      n_ok_cohorts = 0L,
      stringsAsFactors = FALSE
    )
  }
  status_idx <- status_idx + 1L
}

km_summary_table <- do.call(rbind, km_summary_rows)
km_gene_summary <- do.call(rbind, km_status_rows) %>%
  mutate(best_logrank_p_adj_bh = p.adjust(best_logrank_p, method = "BH")) %>%
  arrange(best_logrank_p, gene)

write_tsv_simple(km_summary_table, "kcn25_multicohort_km_within_cohort_median.tsv", tab_dir)
write_tsv_simple(km_gene_summary, "kcn25_multicohort_km_within_cohort_median_gene_summary.tsv", tab_dir)
write_tsv_simple(kcn_table, "kcn25_input_from_venn.tsv", tab_dir)

