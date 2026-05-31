#!/usr/bin/env Rscript

###############################################################################
# Two-stage multicohort survival meta-analysis for the 25 KCN genes
#
# Stage 1: fit one adjusted Cox model per cohort
# Stage 2: pool cohort-specific log(HR) estimates with a random-effects model
#          using REML and Knapp-Hartung inference
###############################################################################

suppressPackageStartupMessages({
  library(survival)
  library(metafor)
  library(ggplot2)
  library(dplyr)
  library(tibble)
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

clean_theme <- function() {
  theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.title = element_text(face = "bold"),
      panel.grid = element_blank()
    )
}

safe_cox_extract <- function(fit, term = "gene_value") {
  sm <- summary(fit)
  coefs <- sm$coefficients
  conf <- sm$conf.int

  if (!term %in% rownames(coefs)) {
    return(tibble(
      beta = NA_real_,
      hr = NA_real_,
      hr_low = NA_real_,
      hr_high = NA_real_,
      se = NA_real_,
      z = NA_real_,
      p_value = NA_real_,
      concordance = sm$concordance[1]
    ))
  }

  tibble(
    beta = unname(coefs[term, "coef"]),
    hr = unname(conf[term, "exp(coef)"]),
    hr_low = unname(conf[term, "lower .95"]),
    hr_high = unname(conf[term, "upper .95"]),
    se = unname(coefs[term, "se(coef)"]),
    z = unname(coefs[term, "z"]),
    p_value = unname(coefs[term, "Pr(>|z|)"]),
    concordance = sm$concordance[1]
  )
}

build_cox_formula <- function(df) {
  covariates <- "gene_value"

  if (dplyr::n_distinct(df$stage_group) > 1) {
    covariates <- c(covariates, "stage_group")
  }
  if (dplyr::n_distinct(df$grade_group) > 1) {
    covariates <- c(covariates, "grade_group")
  }

  list(
    formula = as.formula(
      paste("Surv(os_months, event) ~", paste(covariates, collapse = " + "))
    ),
    covariates = paste(covariates, collapse = " + ")
  )
}

run_single_cohort_cox <- function(df, gene) {
  formula_info <- build_cox_formula(df)

  fit <- tryCatch(
    coxph(formula_info$formula, data = df, ties = "efron"),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    return(tibble(
      gene = gene,
      cohort = unique(df$cohort),
      n_samples = nrow(df),
      n_events = sum(df$event == 1, na.rm = TRUE),
      covariates_used = formula_info$covariates,
      beta = NA_real_,
      hr = NA_real_,
      hr_low = NA_real_,
      hr_high = NA_real_,
      se = NA_real_,
      z = NA_real_,
      p_value = NA_real_,
      concordance = NA_real_,
      fit_status = paste("fit_failed:", fit$message)
    ))
  }

  safe_cox_extract(fit) %>%
    mutate(
      gene = gene,
      cohort = unique(df$cohort),
      n_samples = nrow(df),
      n_events = sum(df$event == 1, na.rm = TRUE),
      covariates_used = formula_info$covariates,
      fit_status = "ok",
      .before = 1
    )
}

run_random_effects_meta <- function(gene_df) {
  meta_df <- gene_df %>%
    filter(is.finite(beta), is.finite(se), se > 0)

  if (nrow(meta_df) < 2) {
    return(tibble(
      gene = unique(gene_df$gene),
      n_cohorts = nrow(meta_df),
      beta_pooled = NA_real_,
      se_pooled = NA_real_,
      hr = NA_real_,
      hr_low = NA_real_,
      hr_high = NA_real_,
      p_value = NA_real_,
      tau2 = NA_real_,
      I2 = NA_real_,
      H2 = NA_real_,
      QE = NA_real_,
      QEp = NA_real_,
      pred_hr_low = NA_real_,
      pred_hr_high = NA_real_,
      model_status = "fewer_than_two_cohorts"
    ))
  }

  fit <- tryCatch(
    rma.uni(
      yi = meta_df$beta,
      sei = meta_df$se,
      method = "REML",
      test = "knha",
      slab = meta_df$cohort
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    return(tibble(
      gene = unique(gene_df$gene),
      n_cohorts = nrow(meta_df),
      beta_pooled = NA_real_,
      se_pooled = NA_real_,
      hr = NA_real_,
      hr_low = NA_real_,
      hr_high = NA_real_,
      p_value = NA_real_,
      tau2 = NA_real_,
      I2 = NA_real_,
      H2 = NA_real_,
      QE = NA_real_,
      QEp = NA_real_,
      pred_hr_low = NA_real_,
      pred_hr_high = NA_real_,
      model_status = paste("meta_failed:", fit$message)
    ))
  }

  pred <- tryCatch(
    predict(fit),
    error = function(e) NULL
  )

  tibble(
    gene = unique(gene_df$gene),
    n_cohorts = nrow(meta_df),
    beta_pooled = as.numeric(fit$b[1, 1]),
    se_pooled = fit$se,
    hr = exp(as.numeric(fit$b[1, 1])),
    hr_low = exp(fit$ci.lb),
    hr_high = exp(fit$ci.ub),
    p_value = fit$pval,
    tau2 = fit$tau2,
    I2 = fit$I2,
    H2 = fit$H2,
    QE = fit$QE,
    QEp = fit$QEp,
    pred_hr_low = if (is.null(pred)) NA_real_ else exp(pred$pi.lb),
    pred_hr_high = if (is.null(pred)) NA_real_ else exp(pred$pi.ub),
    model_status = "ok"
  )
}

make_gene_forest_plot <- function(gene, gene_df, meta_row, out_dir) {
  meta_df <- gene_df %>%
    filter(gene == !!gene, is.finite(beta), is.finite(se), se > 0)

  if (nrow(meta_df) < 2 || !is.finite(meta_row$beta_pooled)) {
    return(invisible(NULL))
  }

  fit <- rma.uni(
    yi = meta_df$beta,
    sei = meta_df$se,
    method = "REML",
    test = "knha",
    slab = meta_df$cohort
  )

  pdf(
    file = file.path(out_dir, paste0(gene, "_two_stage_meta_forest.pdf")),
    width = 8.5,
    height = 5.5,
    bg = "white"
  )

  forest(
    fit,
    atransf = exp,
    refline = 1,
    addpred = TRUE,
    xlab = "Hazard ratio per 1 SD increase",
    slab = meta_df$cohort,
    mlab = "Random-effects pooled HR",
    header = c("Cohort", "HR [95% CI]"),
    main = paste0(gene, ": two-stage multicohort meta-analysis")
  )

  text(
    x = par("usr")[1],
    y = -1.5,
    pos = 4,
    cex = 0.9,
    labels = paste0(
      "I2 = ", sprintf("%.1f", meta_row$I2), "%; tau2 = ", sprintf("%.4f", meta_row$tau2),
      "; pooled p = ", formatC(meta_row$p_value, format = "e", digits = 2)
    )
  )

  dev.off()
}

script_dir <- get_script_dir()
bundle_path <- file.path(script_dir, "integrated_survival_bundle", "fused_kcn_survival_multicohort_bundle.rds")
kcn_path <- file.path(script_dir, "..", "moffitt_kcn_deg_4_methods", "tables", "KCN_union_membership.tsv")
out_dir <- file.path(script_dir, "survival_kcn25_two_stage_meta")
tab_dir <- file.path(out_dir, "tables")
fig_dir <- file.path(out_dir, "figures")
gene_fig_dir <- file.path(fig_dir, "by_gene")

dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(gene_fig_dir, recursive = TRUE, showWarnings = FALSE)

bundle <- readRDS(bundle_path)
kcn_table <- read.delim(kcn_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
kcn_genes <- unique(kcn_table$gene)

cohort_order <- bundle$cohort_order
ph <- bundle$phenotype_harmonized %>%
  mutate(
    cohort = as.character(cohort),
    sample_id = as.character(sample_id),
    os_months = as.numeric(os_months),
    event = as.numeric(event),
    stage_group = factor(as.character(stage_group)),
    grade_group = factor(as.character(grade_group))
  )

stage1_results <- vector("list", length(kcn_genes) * length(cohort_order))
idx <- 1L

for (gene in kcn_genes) {
  for (cohort_name in cohort_order) {
    expr_mat <- bundle$z_common_by_cohort[[cohort_name]]
    if (is.null(expr_mat) || !gene %in% rownames(expr_mat)) {
      stage1_results[[idx]] <- tibble(
        gene = gene,
        cohort = cohort_name,
        n_samples = NA_integer_,
        n_events = NA_integer_,
        covariates_used = NA_character_,
        beta = NA_real_,
        hr = NA_real_,
        hr_low = NA_real_,
        hr_high = NA_real_,
        se = NA_real_,
        z = NA_real_,
        p_value = NA_real_,
        concordance = NA_real_,
        fit_status = "gene_not_available"
      )
      idx <- idx + 1L
      next
    }

    ph_cohort <- ph %>%
      filter(cohort == cohort_name)

    ph_cohort <- ph_cohort[match(colnames(expr_mat), ph_cohort$sample_id), , drop = FALSE]

    if (!identical(colnames(expr_mat), ph_cohort$sample_id)) {
      stop("Sample order mismatch for cohort: ", cohort_name)
    }

    gene_df <- ph_cohort %>%
      transmute(
        cohort = cohort,
        sample_id = sample_id,
        os_months = os_months,
        event = event,
        stage_group = droplevels(stage_group),
        grade_group = droplevels(grade_group),
        gene_value = as.numeric(expr_mat[gene, ])
      )

    stage1_results[[idx]] <- run_single_cohort_cox(gene_df, gene)
    idx <- idx + 1L
  }
}

stage1_table <- bind_rows(stage1_results) %>%
  arrange(gene, cohort)

meta_table <- stage1_table %>%
  group_by(gene) %>%
  group_modify(~ run_random_effects_meta(.x)) %>%
  ungroup() %>%
  mutate(
    p_adj_bh = p.adjust(p_value, method = "BH"),
    effect_direction = case_when(
      hr > 1 ~ "Adverse_when_high",
      hr < 1 ~ "Protective_when_high",
      TRUE ~ "Neutral"
    ),
    significant_fdr = p_adj_bh < 0.05
  ) %>%
  arrange(p_adj_bh, p_value, desc(abs(log(hr))))

significant_meta <- meta_table %>%
  filter(is.finite(p_adj_bh), p_adj_bh < 0.05)

leave_one_out_list <- vector("list", length(kcn_genes))

for (i in seq_along(kcn_genes)) {
  gene <- kcn_genes[[i]]
  gene_df <- stage1_table %>%
    filter(gene == !!gene, is.finite(beta), is.finite(se), se > 0)

  if (nrow(gene_df) < 3) {
    leave_one_out_list[[i]] <- tibble(
      gene = gene,
      omitted_cohort = NA_character_,
      hr = NA_real_,
      hr_low = NA_real_,
      hr_high = NA_real_,
      p_value = NA_real_,
      tau2 = NA_real_,
      I2 = NA_real_,
      status = "not_enough_cohorts"
    )
    next
  }

  loo_rows <- vector("list", nrow(gene_df))

  for (j in seq_len(nrow(gene_df))) {
    loo_df <- gene_df[-j, , drop = FALSE]
    fit <- tryCatch(
      rma.uni(
        yi = loo_df$beta,
        sei = loo_df$se,
        method = "REML",
        test = "knha",
        slab = loo_df$cohort
      ),
      error = function(e) e
    )

    if (inherits(fit, "error")) {
      loo_rows[[j]] <- tibble(
        gene = gene,
        omitted_cohort = gene_df$cohort[j],
        hr = NA_real_,
        hr_low = NA_real_,
        hr_high = NA_real_,
        p_value = NA_real_,
        tau2 = NA_real_,
        I2 = NA_real_,
        status = paste("fit_failed:", fit$message)
      )
    } else {
      loo_rows[[j]] <- tibble(
        gene = gene,
        omitted_cohort = gene_df$cohort[j],
        hr = exp(as.numeric(fit$b[1, 1])),
        hr_low = exp(fit$ci.lb),
        hr_high = exp(fit$ci.ub),
        p_value = fit$pval,
        tau2 = fit$tau2,
        I2 = fit$I2,
        status = "ok"
      )
    }
  }

  leave_one_out_list[[i]] <- bind_rows(loo_rows)
}

leave_one_out_table <- bind_rows(leave_one_out_list) %>%
  arrange(gene, omitted_cohort)

summary_plot_table <- meta_table %>%
  filter(is.finite(hr), is.finite(hr_low), is.finite(hr_high)) %>%
  mutate(
    gene = factor(gene, levels = rev(gene))
  )

summary_plot <- ggplot(summary_plot_table, aes(y = gene, x = hr, color = effect_direction)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey55", linewidth = 0.4) +
  geom_segment(
    aes(x = pred_hr_low, xend = pred_hr_high, yend = gene),
    color = "grey80",
    linewidth = 1.1,
    na.rm = TRUE
  ) +
  geom_segment(
    aes(x = hr_low, xend = hr_high, yend = gene),
    linewidth = 0.7
  ) +
  geom_point(aes(size = -log10(p_value)), alpha = 0.95) +
  scale_x_log10() +
  scale_color_manual(
    values = c(
      "Adverse_when_high" = "#800E13",
      "Protective_when_high" = "#2D6A9F",
      "Neutral" = "grey55"
    ),
    name = "Direction"
  ) +
  scale_size_continuous(range = c(2.5, 6), name = "-log10(p)") +
  labs(
    title = "Two-stage multicohort survival meta-analysis of the 25 KCN genes",
    subtitle = "Per-cohort adjusted Cox models pooled with REML random effects; grey bars show prediction intervals",
    x = "Hazard ratio per 1 SD increase in within-cohort expression",
    y = NULL
  ) +
  clean_theme()

ggsave(
  filename = file.path(fig_dir, "Figure_01_KCN25_two_stage_meta_summary.pdf"),
  plot = summary_plot,
  width = 11,
  height = 8.5,
  dpi = 320,
  bg = "white"
)

ggsave(
  filename = file.path(fig_dir, "Figure_01_KCN25_two_stage_meta_summary.png"),
  plot = summary_plot,
  width = 11,
  height = 8.5,
  dpi = 320,
  bg = "white"
)

for (i in seq_len(nrow(meta_table))) {
  make_gene_forest_plot(
    gene = meta_table$gene[i],
    gene_df = stage1_table,
    meta_row = meta_table[i, , drop = FALSE],
    out_dir = gene_fig_dir
  )
}

cohort_summary_table <- ph %>%
  group_by(cohort) %>%
  summarise(
    n_samples = n(),
    n_events = sum(event == 1, na.rm = TRUE),
    median_os_months = median(os_months, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(match(cohort, cohort_order))

analysis_summary <- tibble(
  n_cohorts = length(cohort_order),
  n_samples = nrow(ph),
  n_events = sum(ph$event == 1, na.rm = TRUE),
  n_kcn_genes = length(kcn_genes),
  n_significant_fdr = sum(meta_table$significant_fdr, na.rm = TRUE)
)

write_tsv_simple(stage1_table, "kcn25_multicohort_two_stage_stage1_cox.tsv", tab_dir)
write_tsv_simple(meta_table, "kcn25_multicohort_two_stage_meta.tsv", tab_dir)
write_tsv_simple(significant_meta, "kcn25_multicohort_two_stage_meta_significant.tsv", tab_dir)
write_tsv_simple(leave_one_out_table, "kcn25_multicohort_two_stage_leave_one_out.tsv", tab_dir)
write_tsv_simple(cohort_summary_table, "kcn25_multicohort_two_stage_cohort_summary.tsv", tab_dir)
write_tsv_simple(analysis_summary, "kcn25_multicohort_two_stage_analysis_summary.tsv", tab_dir)
write_tsv_simple(kcn_table, "kcn25_input_from_venn.tsv", tab_dir)

saveRDS(
  list(
    stage1_table = stage1_table,
    meta_table = meta_table,
    significant_meta = significant_meta,
    leave_one_out_table = leave_one_out_table,
    cohort_summary = cohort_summary_table,
    analysis_summary = analysis_summary
  ),
  file = file.path(out_dir, "kcn25_multicohort_two_stage_meta_results.rds")
)

