# Fused Multicohort Survival Bundle Builder
#
# Purpose:
#   Build a single patient-level survival bundle from four public PDAC bulk
#   cohorts while preserving cohort-specific information and avoiding naive
#   cross-platform intensity merging.
#
# Why this approach:
#   - The input cohorts mix microarray and RNA-seq platforms.
#   - Directly concatenating raw expression values would create strong platform
#     artifacts and would not be statistically defensible.
#   - A robust strategy is to keep the original cohort bundles intact, extract
#     the common KCN gene space, restrict the survival dataset to eligible tumor
#     samples, and standardize each gene *within each cohort* before pooling.
#   - This preserves biological ranking within cohort while making effect sizes
#     comparable across cohorts.
#
# Output bundle:
#   The resulting RDS stores:
#   - the original cohort bundles unchanged
#   - a harmonized phenotype table for survival-ready tumor samples
#   - a pooled z-scored expression matrix on the common 93 KCN genes
#   - cohort-specific raw matrices on the same common genes
#   - QC and exclusion tables
#
# Usage:
#   Rscript build_fused_survival_bundle.R
#
# Inputs expected in the same folder:
#   - PH_GSE183795_allKCN_bundle.rds
#   - Hussain_GSE62452_bundle.rds
#   - Zhang_GSE28735_bundle.rds
#   - Bailey_QCMG_UQ_2016_allKCN_bundle.rds
#
# Main outputs:
#   - integrated_survival_bundle/fused_kcn_survival_multicohort_bundle.rds
#   - integrated_survival_bundle/tables/fused_survival_cohort_summary.tsv
#   - integrated_survival_bundle/tables/fused_survival_gene_availability.tsv
#   - integrated_survival_bundle/tables/fused_survival_sample_inclusion.tsv

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
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

safe_scale_vector <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) {
    return(rep(0, length(x)))
  }
  as.numeric((x - mean(x, na.rm = TRUE)) / s)
}

harmonize_stage <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "N/A", "<NA>")] <- NA
  x <- toupper(gsub("[[:space:]]+", "", x))
  x <- gsub("^>", "", x)
  dplyr::case_when(
    is.na(x) ~ NA_character_,
    grepl("^IV", x) ~ "IV",
    grepl("^III", x) ~ "III",
    grepl("^II", x) ~ "II",
    grepl("^I", x) ~ "I",
    TRUE ~ NA_character_
  )
}

harmonize_grade <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "N/A", "<NA>")] <- NA
  digit <- gsub(".*([1-4]).*", "\\1", x)
  digit[!grepl("[1-4]", x)] <- NA
  digit
}

harmonize_margin <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "N/A", "<NA>")] <- NA
  x <- toupper(x)
  dplyr::case_when(
    grepl("^R0$", x) ~ "R0",
    grepl("^R1$", x) ~ "R1",
    TRUE ~ NA_character_
  )
}

load_bundle <- function(path) {
  x <- readRDS(path)
  required <- c("cohort", "display_name", "source_label", "source_url", "expr_gene", "phenotype")
  missing <- setdiff(required, names(x))
  if (length(missing) > 0) {
    stop("Bundle missing required fields: ", paste(missing, collapse = ", "), " in ", path)
  }
  x
}

extract_survival_ready_samples <- function(bundle) {
  ph <- bundle$phenotype
  ph$sample_id <- as.character(ph$sample_id)
  ph$sample_group <- as.character(ph$sample_group)
  ph$survival_months <- as.numeric(ph$survival_months)
  ph$survival_status <- as.numeric(ph$survival_status)

  ph$cohort <- bundle$cohort
  ph$display_name <- bundle$display_name
  ph$source_label <- bundle$source_label
  ph$source_url <- bundle$source_url

  ph$is_tumor <- !is.na(ph$sample_group) & ph$sample_group == "Tumor"
  ph$has_survival <- !is.na(ph$survival_months) & ph$survival_months > 0 & !is.na(ph$survival_status)
  ph$eligible_survival <- ph$is_tumor & ph$has_survival

  ph$stage_harmonized <- harmonize_stage(ph$stage_geo)
  ph$grade_harmonized <- harmonize_grade(ph$grade_geo)
  ph$margin_harmonized <- harmonize_margin(ph$resection_margin_geo)
  ph$stage_group <- ifelse(is.na(ph$stage_harmonized), "Unknown", ph$stage_harmonized)
  ph$grade_group <- dplyr::case_when(
    is.na(ph$grade_harmonized) ~ "Unknown",
    ph$grade_harmonized %in% c("1", "2") ~ "Low_1_2",
    ph$grade_harmonized %in% c("3", "4") ~ "High_3_4",
    TRUE ~ "Unknown"
  )
  ph$margin_group <- ifelse(is.na(ph$margin_harmonized), "Unknown", ph$margin_harmonized)

  ph$eligibility_reason <- dplyr::case_when(
    ph$eligible_survival ~ "included",
    !ph$is_tumor ~ "excluded_non_tumor",
    !ph$has_survival ~ "excluded_missing_survival",
    TRUE ~ "excluded_other"
  )

  ph
}

script_dir <- get_script_dir()
out_dir <- file.path(script_dir, "integrated_survival_bundle")
tab_dir <- file.path(out_dir, "tables")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

bundle_paths <- c(
  file.path(script_dir, "PH_GSE183795_allKCN_bundle.rds"),
  file.path(script_dir, "Hussain_GSE62452_bundle.rds"),
  file.path(script_dir, "Zhang_GSE28735_bundle.rds"),
  file.path(script_dir, "Bailey_QCMG_UQ_2016_allKCN_bundle.rds")
)

bundles <- lapply(bundle_paths, load_bundle)
names(bundles) <- vapply(bundles, function(x) x$cohort, character(1))

# Common gene space across all cohorts. This avoids platform-specific missingness
# in the pooled survival matrix while retaining original raw matrices separately.
common_genes <- Reduce(intersect, lapply(bundles, function(x) rownames(x$expr_gene)))
common_genes <- sort(common_genes)

phenotype_all <- bind_rows(lapply(bundles, extract_survival_ready_samples))
phenotype_kept <- phenotype_all %>%
  filter(eligible_survival) %>%
  mutate(
    integrated_sample_id = paste(cohort, sample_id, sep = "__"),
    event = survival_status,
    os_months = survival_months
  )

cohort_summary <- phenotype_all %>%
  group_by(cohort, display_name, source_label, source_url) %>%
  summarise(
    total_samples = n(),
    tumor_samples = sum(is_tumor, na.rm = TRUE),
    survival_eligible_samples = sum(eligible_survival, na.rm = TRUE),
    events = sum(survival_status[eligible_survival] == 1, na.rm = TRUE),
    median_os_months = stats::median(survival_months[eligible_survival], na.rm = TRUE),
    stage_missing_fraction = mean(is.na(stage_harmonized[eligible_survival])),
    grade_missing_fraction = mean(is.na(grade_harmonized[eligible_survival])),
    margin_missing_fraction = mean(is.na(margin_harmonized[eligible_survival])),
    .groups = "drop"
  ) %>%
  arrange(cohort)

raw_common_by_cohort <- list()
z_common_by_cohort <- list()
gene_availability <- list()

for (nm in names(bundles)) {
  bundle <- bundles[[nm]]
  ph <- phenotype_all %>% filter(cohort == nm)
  keep_ids <- ph %>% filter(eligible_survival) %>% pull(sample_id)

  expr <- bundle$expr_gene[common_genes, keep_ids, drop = FALSE]
  raw_common_by_cohort[[nm]] <- expr

  z_expr <- t(apply(expr, 1, safe_scale_vector))
  colnames(z_expr) <- keep_ids
  rownames(z_expr) <- rownames(expr)
  z_common_by_cohort[[nm]] <- z_expr

  gene_availability[[nm]] <- tibble::tibble(
    cohort = nm,
    gene = common_genes,
    present = TRUE,
    variance_raw = apply(expr, 1, stats::sd, na.rm = TRUE),
    n_non_missing = apply(expr, 1, function(v) sum(!is.na(v)))
  )
}

expr_common_z <- do.call(cbind, z_common_by_cohort)
expr_common_raw_tumor <- do.call(cbind, raw_common_by_cohort)

phenotype_kept <- phenotype_kept %>%
  filter(sample_id %in% colnames(expr_common_z)) %>%
  arrange(match(integrated_sample_id, paste(cohort, colnames(expr_common_z), sep = "__")))

colnames(expr_common_z) <- paste(phenotype_kept$cohort, phenotype_kept$sample_id, sep = "__")
colnames(expr_common_raw_tumor) <- paste(phenotype_kept$cohort, phenotype_kept$sample_id, sep = "__")

gene_availability_table <- bind_rows(gene_availability) %>%
  arrange(gene, cohort)

sample_inclusion_table <- phenotype_all %>%
  transmute(
    cohort,
    display_name,
    sample_id,
    sample_group,
    survival_months,
    survival_status,
    eligible_survival,
    eligibility_reason
  ) %>%
  arrange(cohort, sample_id)

fused_bundle <- list(
  bundle_version = "2026-05-03",
  integration_strategy = paste(
    "Tumor-only survival-ready samples pooled across cohorts.",
    "Expression is standardized within cohort on the common KCN gene set",
    "to preserve within-cohort ranking and reduce cross-platform scale bias.",
    "Original cohort bundles are retained unchanged to avoid information loss."
  ),
  cohort_order = names(bundles),
  common_genes = common_genes,
  original_bundles = bundles,
  phenotype_harmonized = phenotype_kept,
  expr_common_z = expr_common_z,
  expr_common_raw_tumor = expr_common_raw_tumor,
  raw_common_by_cohort = raw_common_by_cohort,
  z_common_by_cohort = z_common_by_cohort,
  cohort_summary = cohort_summary,
  gene_availability = gene_availability_table,
  sample_inclusion = sample_inclusion_table
)

saveRDS(fused_bundle, file.path(out_dir, "fused_kcn_survival_multicohort_bundle.rds"))
write.table(cohort_summary, file.path(tab_dir, "fused_survival_cohort_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(gene_availability_table, file.path(tab_dir, "fused_survival_gene_availability.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(sample_inclusion_table, file.path(tab_dir, "fused_survival_sample_inclusion.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
