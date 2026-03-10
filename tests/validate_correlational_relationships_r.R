#!/usr/bin/env Rscript

# Validation harness for `analyze_correlational_relationships.R`.
#
# Purpose
# - Verify subject-ID normalization + exclusion handling
# - Verify pairwise-complete-case filtering (no imputation)
# - Verify ROI averaging over available channels only
# - Verify Pearson / family-wise BH-FDR outputs on deterministic synthetic data
# - Verify per-pair figure artifact generation

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(jsonlite)
})

fail <- function(msg) {
  writeLines(paste0("[FAIL] ", msg))
  quit(status = 1)
}

assert_true <- function(cond, msg) {
  if (!isTRUE(cond)) fail(msg)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

bh_qvalues_manual <- function(p_values) {
  # Manual Benjamini-Hochberg q-values.
  # Reference: Benjamini & Hochberg (1995), see `CITATIONS.md`.
  p <- as.numeric(p_values)
  q <- rep(NA_real_, length(p))
  finite <- is.finite(p)
  if (!any(finite)) return(q)

  p_f <- p[finite]
  o <- order(p_f)
  ranked <- p_f[o]
  m <- length(ranked)
  q_ranked <- rep(NA_real_, m)
  prev <- 1.0
  for (i in seq(m, 1)) {
    val <- (m / i) * ranked[[i]]
    prev <- min(prev, val)
    q_ranked[[i]] <- prev
  }
  q_ranked <- pmin(pmax(q_ranked, 0.0), 1.0)
  q_f <- rep(NA_real_, m)
  q_f[o] <- q_ranked
  q[finite] <- q_f
  q
}

make_subject_ids <- function(n_subjects) {
  ids_int <- seq_len(n_subjects)
  list(
    combined_subject = sprintf("%04d", ids_int),
    exclusion_subject = sprintf("sub_%04d", ids_int),
    ids_int = ids_int
  )
}

cond_codes <- function() c("01", "02", "03", "04")

make_predictor_frame <- function(n_subjects) {
  ids <- make_subject_ids(n_subjects)
  idx <- ids$ids_int
  roi_signal <- seq(5, 14, length.out = n_subjects)

  tibble::tibble(
    subject_id = ids$combined_subject,
    age = 20 + idx,
    sfv_daily_duration = 40 + idx * 2,
    sfv_frequency = roi_signal,
    phq_total = 30 - idx,
    gad_total = 10 + idx * 0.2,
    asrs_total = 2 * idx,
    yang_pu_total = 3 + idx * 0.4,
    yang_mot_total = 5, # constant predictor: must be skipped, not analyzed
    diff_short_form_education = idx * 0.5,
    diff_short_form_entertainment = c(1, 2, 3, 4, 5, rep(NA_real_, n_subjects - 5)),
    diff_long_form_education = 50 + idx,
    diff_long_form_entertainment = 100 - idx,
    lf_education_engagement = 1 + idx * 0.1,
    lf_entertainment_engagement = 2 + idx * 0.2,
    sf_education_engagement = 3 + idx * 0.3,
    sf_entertainment_engagement = 4 + idx * 0.4
  )
}

moderate_signal_hbo_cond02 <- function(n_subjects) {
  vals <- c(1.5, 2.9, 2.1, 4.7, 4.0, 6.4, 5.9, 8.3, 9.1, 8.8)
  if (n_subjects > length(vals)) stop("moderate_signal_hbo_cond02 needs more seed values.")
  vals[seq_len(n_subjects)]
}

make_beta_matrix <- function(n_subjects) {
  idx <- seq_len(n_subjects)
  roi_signal <- seq(5, 14, length.out = n_subjects)
  asrs_signal <- 2 * idx

  out <- tibble::tibble(.rows = n_subjects)

  add_beta_col <- function(name, values) {
    out[[name]] <<- values
  }

  # Target channel S04_D02: both chromophores.
  for (cond in cond_codes()) {
    cond_num <- as.numeric(cond)
    add_beta_col(paste0("S04_D02_Cond", cond, "_HbO"), 100 + idx + cond_num)
    add_beta_col(paste0("S04_D02_Cond", cond, "_HbR"), 200 + idx + cond_num)
  }
  out[["S04_D02_Cond02_HbO"]] <- moderate_signal_hbo_cond02(n_subjects)
  # Force one neural missing value for pairwise-complete-case testing.
  out[["S04_D02_Cond01_HbR"]][[4]] <- NA_real_

  # L_DLPFC (HbO): ROI mean should equal sfv_frequency for Cond01 even when one channel is missing.
  out[["S01_D01_Cond01_HbO"]] <- roi_signal
  out[["S02_D01_Cond01_HbO"]] <- roi_signal
  out[["S02_D01_Cond01_HbO"]][[3]] <- NA_real_ # partial ROI pruning; mean should use S01_D01 only
  for (cond in c("02", "03", "04")) {
    cond_num <- as.numeric(cond)
    out[[paste0("S01_D01_Cond", cond, "_HbO")]] <- roi_signal + cond_num
    out[[paste0("S02_D01_Cond", cond, "_HbO")]] <- roi_signal + cond_num
  }

  # R_DLPFC (HbR): Cond02 mean should be perfectly negatively correlated with phq_total.
  out[["S07_D05_Cond01_HbR"]] <- idx * 0.1
  out[["S07_D07_Cond01_HbR"]] <- idx * 0.1 + 1
  out[["S07_D05_Cond02_HbR"]] <- -make_predictor_frame(n_subjects)$phq_total - 1
  out[["S07_D07_Cond02_HbR"]] <- -make_predictor_frame(n_subjects)$phq_total + 1
  out[["S07_D05_Cond03_HbR"]] <- idx * 0.2
  out[["S07_D07_Cond03_HbR"]] <- idx * 0.2 + 1
  out[["S07_D05_Cond04_HbR"]] <- idx * 0.3
  out[["S07_D07_Cond04_HbR"]] <- idx * 0.3 + 1

  # M_DMPFC (HbO): single-channel ROI.
  out[["S04_D04_Cond02_HbO"]] <- rep(42, n_subjects)
  for (cond in cond_codes()) {
    if (cond != "02") {
      cond_num <- as.numeric(cond)
      out[[paste0("S04_D04_Cond", cond, "_HbO")]] <- 10 + idx * 0.5 + cond_num
    }
  }

  # L_DMPFC (HbO): two-channel ROI.
  for (cond in cond_codes()) {
    cond_num <- as.numeric(cond)
    out[[paste0("S03_D02_Cond", cond, "_HbO")]] <- 20 + idx + cond_num
    out[[paste0("S03_D04_Cond", cond, "_HbO")]] <- 22 + idx + cond_num
  }
  # Entire ROI missing for one subject-condition should propagate as missing only for that pair.
  out[["S03_D02_Cond03_HbO"]][[5]] <- NA_real_
  out[["S03_D04_Cond03_HbO"]][[5]] <- NA_real_

  # Strengthen one exact correlation at the channel level.
  out[["S04_D02_Cond01_HbR"]] <- asrs_signal
  out[["S04_D02_Cond01_HbR"]][[4]] <- NA_real_

  out
}

build_input <- function(n_subjects) {
  predictors <- make_predictor_frame(n_subjects)
  betas <- make_beta_matrix(n_subjects)
  bind_cols(predictors, betas)
}

analysis_input_after_exclusion <- function(input_df) {
  input_df %>% filter(.data$subject_id != "0002")
}

compute_expected_pair_stats <- function(input_df, predictor_col, neural_col) {
  df <- analysis_input_after_exclusion(input_df) %>%
    transmute(
      predictor_value = .data[[predictor_col]],
      neural_value = .data[[neural_col]]
    ) %>%
    filter(!is.na(.data$predictor_value), !is.na(.data$neural_value))

  cor_fit <- stats::cor.test(df$predictor_value, df$neural_value, method = "pearson")
  lm_fit <- stats::lm(neural_value ~ predictor_value, data = df)
  coef_fit <- coef(lm_fit)
  list(
    n_complete = nrow(df),
    pearson_r = unname(cor_fit$estimate),
    p_unc = cor_fit$p.value,
    ci95_low = cor_fit$conf.int[[1]],
    ci95_high = cor_fit$conf.int[[2]],
    slope = unname(coef_fit[["predictor_value"]]),
    intercept = unname(coef_fit[["(Intercept)"]])
  )
}

write_roi_json <- function(path) {
  roi_obj <- list(
    R_DLPFC = c("S07_D05", "S07_D07"),
    L_DLPFC = c("S01_D01", "S02_D01"),
    M_DMPFC = c("S04_D04"),
    L_DMPFC = c("S03_D02", "S03_D04")
  )
  jsonlite::write_json(roi_obj, path = path, auto_unbox = TRUE, pretty = TRUE)
}

write_roi_overlap_json <- function(path) {
  roi_obj <- list(
    R_DLPFC = c("S07_D05", "S07_D07"),
    L_DLPFC = c("S01_D01", "S02_D01"),
    M_DMPFC = c("S04_D04"),
    L_DMPFC = c("S03_D02", "S04_D04")
  )
  jsonlite::write_json(roi_obj, path = path, auto_unbox = TRUE, pretty = TRUE)
}

write_roi_missing_channel_json <- function(path) {
  roi_obj <- list(
    R_DLPFC = c("S07_D05", "S07_D07"),
    L_DLPFC = c("S01_D01", "S99_D99"),
    M_DMPFC = c("S04_D04"),
    L_DMPFC = c("S03_D02", "S03_D04")
  )
  jsonlite::write_json(roi_obj, path = path, auto_unbox = TRUE, pretty = TRUE)
}

write_analysis_plan_json <- function(path) {
  plan_obj <- list(
    version = 1,
    description = "Synthetic validation fixture for analyze_correlational_relationships.R",
    predictors = c(
      "age",
      "sfv_daily_duration",
      "sfv_frequency",
      "phq_total",
      "gad_total",
      "asrs_total",
      "yang_pu_total",
      "yang_mot_total",
      "diff_short_form_education",
      "diff_short_form_entertainment",
      "diff_long_form_education",
      "diff_long_form_entertainment",
      "lf_education_engagement",
      "lf_entertainment_engagement",
      "sf_education_engagement",
      "sf_entertainment_engagement"
    ),
    multiple_testing = list(
      adjust_method = "BH",
      family_grouping = c("neural_level", "neural_name", "chrom")
    ),
    figures = list(
      policy = "significant_only"
    )
  )
  jsonlite::write_json(plan_obj, path = path, auto_unbox = TRUE, pretty = TRUE)
}

run_analysis <- function(input_csv, roi_json, analysis_plan_json, exclude_json, out_csv, out_fig_dir, alpha, min_subjects) {
  script <- "analyze_correlational_relationships.R"
  args <- c(
    "--input_csv", input_csv,
    "--roi_json", roi_json,
    "--analysis_plan_json", analysis_plan_json,
    "--exclude_subjects_json", exclude_json,
    "--out_csv", out_csv,
    "--out_fig_dir", out_fig_dir,
    "--alpha", as.character(alpha),
    "--min_subjects", as.character(min_subjects)
  )
  status <- suppressWarnings(system2("Rscript", c(script, args), stdout = TRUE, stderr = TRUE))
  attr(status, "status") %||% 0
}

main <- function() {
  tmp <- file.path(tempdir(), "correlational_relationships_validation")
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)

  input_csv <- file.path(tmp, "merged.csv")
  roi_json <- file.path(tmp, "roi_definition.json")
  analysis_plan_json <- file.path(tmp, "correlational_analysis_plan.json")
  exclude_json <- file.path(tmp, "excluded_subjects.json")
  out_csv <- file.path(tmp, "pairwise_correlations_r.csv")
  out_fig_dir <- file.path(tmp, "figures")

  n_subjects <- 10
  input <- build_input(n_subjects)

  # Missing predictor: should reduce only affected pair counts.
  input$age[[1]] <- NA_real_

  write_csv(input, input_csv)
  write_roi_json(roi_json)
  write_analysis_plan_json(analysis_plan_json)
  writeLines('["sub_0002"]', exclude_json)

  status <- run_analysis(
    input_csv = input_csv,
    roi_json = roi_json,
    analysis_plan_json = analysis_plan_json,
    exclude_json = exclude_json,
    out_csv = out_csv,
    out_fig_dir = out_fig_dir,
    alpha = 0.05,
    min_subjects = 6
  )

  assert_true(status == 0, "analysis script did not exit cleanly")
  assert_true(file.exists(out_csv), "missing correlation output CSV")
  assert_true(dir.exists(out_fig_dir), "missing figures output directory")

  out <- read_csv(out_csv, show_col_types = FALSE)

  # 1) Output schema + expected row count.
  required_cols <- c(
    "neural_level", "neural_name", "chrom", "condition_code", "condition_label",
    "predictor", "analysis_status", "skip_reason", "n_complete",
    "pearson_r", "r_squared", "p_unc", "ci95_low", "ci95_high",
    "slope", "intercept", "plot_file", "family_id", "family_adjust_method",
    "family_n_tested", "p_fdr"
  )
  assert_true(all(required_cols %in% names(out)), "output CSV missing required columns")
  assert_true(nrow(out) == 16 * 24, "unexpected number of tested predictor x neural-target pairs")

  # 2) Pairwise-complete-case with missing predictor only.
  row_age <- out %>%
    filter(
      neural_level == "channel",
      neural_name == "S04_D02",
      chrom == "HbO",
      condition_code == "Cond01",
      predictor == "age"
    )
  assert_true(nrow(row_age) == 1, "missing age/channel row")
  assert_true(row_age$n_complete[[1]] == 8, "age missingness should exclude exactly one remaining subject for this pair")
  expected_age <- compute_expected_pair_stats(input, predictor_col = "age", neural_col = "S04_D02_Cond01_HbO")
  assert_true(abs(row_age$pearson_r[[1]] - expected_age$pearson_r) < 1e-12, "age/channel r mismatch against manual calculation")
  assert_true(abs(row_age$p_unc[[1]] - expected_age$p_unc) < 1e-12, "age/channel p mismatch against manual calculation")
  assert_true(abs(row_age$ci95_low[[1]] - expected_age$ci95_low) < 1e-12, "age/channel CI low mismatch")
  assert_true(abs(row_age$ci95_high[[1]] - expected_age$ci95_high) < 1e-12, "age/channel CI high mismatch")
  assert_true(abs(row_age$slope[[1]] - expected_age$slope) < 1e-12, "age/channel slope mismatch")
  assert_true(abs(row_age$intercept[[1]] - expected_age$intercept) < 1e-12, "age/channel intercept mismatch")

  # 3) Pairwise-complete-case with missing neural value only.
  row_asrs <- out %>%
    filter(
      neural_level == "channel",
      neural_name == "S04_D02",
      chrom == "HbR",
      condition_code == "Cond01",
      predictor == "asrs_total"
    )
  assert_true(nrow(row_asrs) == 1, "missing asrs/channel row")
  assert_true(row_asrs$n_complete[[1]] == 8, "neural missingness should exclude exactly one remaining subject for this pair")
  assert_true(row_asrs$analysis_status[[1]] == "tested", "strong channel correlation should be tested")
  assert_true(abs(row_asrs$pearson_r[[1]] - 1.0) < 1e-12, "expected perfect positive channel correlation")
  assert_true(!grepl("0002", paste(row_asrs$plot_file[[1]], collapse = "")), "excluded subject should not appear in plot path metadata")

  # 3b) Non-perfect correlation should match manual Pearson/CI/LM values exactly.
  row_moderate <- out %>%
    filter(
      neural_level == "channel",
      neural_name == "S04_D02",
      chrom == "HbO",
      condition_code == "Cond02",
      predictor == "sfv_frequency"
    )
  assert_true(nrow(row_moderate) == 1, "missing moderate-correlation row")
  expected_moderate <- compute_expected_pair_stats(input, predictor_col = "sfv_frequency", neural_col = "S04_D02_Cond02_HbO")
  assert_true(row_moderate$n_complete[[1]] == expected_moderate$n_complete, "moderate-correlation n mismatch")
  assert_true(abs(row_moderate$pearson_r[[1]] - expected_moderate$pearson_r) < 1e-12, "moderate-correlation r mismatch")
  assert_true(abs(row_moderate$p_unc[[1]] - expected_moderate$p_unc) < 1e-12, "moderate-correlation p mismatch")
  assert_true(abs(row_moderate$ci95_low[[1]] - expected_moderate$ci95_low) < 1e-12, "moderate-correlation CI low mismatch")
  assert_true(abs(row_moderate$ci95_high[[1]] - expected_moderate$ci95_high) < 1e-12, "moderate-correlation CI high mismatch")
  assert_true(abs(row_moderate$slope[[1]] - expected_moderate$slope) < 1e-12, "moderate-correlation slope mismatch")
  assert_true(abs(row_moderate$intercept[[1]] - expected_moderate$intercept) < 1e-12, "moderate-correlation intercept mismatch")

  # 4) ROI averaging over available channels only.
  row_roi <- out %>%
    filter(
      neural_level == "roi",
      neural_name == "L_DLPFC",
      chrom == "HbO",
      condition_code == "Cond01",
      predictor == "sfv_frequency"
    )
  assert_true(nrow(row_roi) == 1, "missing L_DLPFC ROI row")
  assert_true(row_roi$n_complete[[1]] == 9, "partial ROI pruning should not drop the subject when one channel remains")
  assert_true(row_roi$analysis_status[[1]] == "tested", "ROI pair should be tested")
  assert_true(abs(row_roi$pearson_r[[1]] - 1.0) < 1e-12, "ROI averaging should recover the expected perfect correlation")

  # 4b) Entire ROI missing for one subject-condition should exclude only that subject for that pair.
  row_roi_all_missing <- out %>%
    filter(
      neural_level == "roi",
      neural_name == "L_DMPFC",
      chrom == "HbO",
      condition_code == "Cond03",
      predictor == "sfv_frequency"
    )
  assert_true(nrow(row_roi_all_missing) == 1, "missing L_DMPFC all-missing ROI row")
  assert_true(row_roi_all_missing$n_complete[[1]] == 8, "all-missing ROI condition should exclude exactly one non-excluded subject")
  assert_true(row_roi_all_missing$analysis_status[[1]] == "tested", "all-missing ROI case should still be tested when enough subjects remain")

  # 5) Negative correlation check on another ROI.
  row_neg <- out %>%
    filter(
      neural_level == "roi",
      neural_name == "R_DLPFC",
      chrom == "HbR",
      condition_code == "Cond02",
      predictor == "phq_total"
    )
  assert_true(nrow(row_neg) == 1, "missing negative ROI row")
  assert_true(row_neg$analysis_status[[1]] == "tested", "negative ROI pair should be tested")
  assert_true(abs(row_neg$pearson_r[[1]] + 1.0) < 1e-12, "expected perfect negative ROI correlation")

  # 6) Constant predictor should be skipped, not analyzed.
  row_const <- out %>%
    filter(
      neural_level == "channel",
      neural_name == "S04_D02",
      chrom == "HbO",
      condition_code == "Cond01",
      predictor == "yang_mot_total"
    )
  assert_true(nrow(row_const) == 1, "missing constant predictor row")
  assert_true(row_const$analysis_status[[1]] == "skipped_constant_input", "constant predictor should be skipped explicitly")
  assert_true(row_const$skip_reason[[1]] == "predictor_has_zero_variance", "constant predictor should record an explicit skip reason")
  assert_true(is.na(row_const$p_fdr[[1]]), "skipped constant predictor should not receive BH-adjusted p")
  assert_true(is.na(row_const$plot_file[[1]]), "skipped constant predictor should not create a plot")

  # 7) Constant neural target should be skipped explicitly.
  row_const_neural <- out %>%
    filter(
      neural_level == "roi",
      neural_name == "M_DMPFC",
      chrom == "HbO",
      condition_code == "Cond02",
      predictor == "sfv_frequency"
    )
  assert_true(nrow(row_const_neural) == 1, "missing constant-neural row")
  assert_true(row_const_neural$analysis_status[[1]] == "skipped_constant_input", "constant neural target should be skipped explicitly")
  assert_true(row_const_neural$skip_reason[[1]] == "neural_value_has_zero_variance", "constant neural target should record the correct skip reason")

  # 8) Too-few-subjects predictor should be skipped by min_subjects.
  row_min <- out %>%
    filter(
      neural_level == "channel",
      neural_name == "S04_D02",
      chrom == "HbO",
      condition_code == "Cond01",
      predictor == "diff_short_form_entertainment"
    )
  assert_true(nrow(row_min) == 1, "missing min-subjects row")
  assert_true(row_min$analysis_status[[1]] == "skipped_min_subjects", "expected min-subjects skip for sparse predictor")

  # 9) BH-FDR must match a manual implementation within each configured family.
  analyzed <- out %>% filter(.data$analysis_status == "tested")
  assert_true(all(analyzed$family_adjust_method == "BH"), "unexpected adjustment method in analyzed rows")
  family_counts <- analyzed %>% count(.data$family_id, name = "n_family")
  joined_counts <- analyzed %>% select("family_id", "family_n_tested") %>% distinct()
  assert_true(
    all(joined_counts$family_n_tested == family_counts$n_family[match(joined_counts$family_id, family_counts$family_id)]),
    "family_n_tested does not match observed tested rows per family"
  )
  family_summaries <- analyzed %>%
    group_by(.data$family_id) %>%
    summarize(max_abs_diff = max(abs(bh_qvalues_manual(.data$p_unc) - .data$p_fdr), na.rm = TRUE), .groups = "drop")
  assert_true(
    all(family_summaries$max_abs_diff < 1e-12),
    "family-wise BH-FDR values do not match manual implementation"
  )

  # 10) Figures should be created only for FDR-significant analyzed rows under the configured policy.
  plotted <- analyzed %>% filter(!is.na(.data$plot_file))
  significant <- analyzed %>% filter(is.finite(.data$p_fdr), .data$p_fdr < 0.05)
  non_significant <- analyzed %>% filter(!is.finite(.data$p_fdr) | .data$p_fdr >= 0.05)
  assert_true(nrow(significant) > 0, "synthetic fixture should yield at least one significant analyzed row")
  assert_true(all(!is.na(significant$plot_file)), "all significant analyzed rows should have a plot file path")
  assert_true(all(file.exists(significant$plot_file)), "one or more significant analyzed plot files do not exist")
  assert_true(all(is.na(non_significant$plot_file)), "non-significant analyzed rows should not create plots under significant_only policy")

  png_files <- list.files(out_fig_dir, pattern = "\\.png$", full.names = TRUE)
  assert_true(length(png_files) == nrow(plotted), "expected one PNG per plotted analyzed row")

  # 11) Duplicate subject IDs must fail hard.
  input_dup <- bind_rows(input, input[1, , drop = FALSE])
  input_dup_csv <- file.path(tmp, "merged_dup.csv")
  write_csv(input_dup, input_dup_csv)
  status_dup <- run_analysis(
    input_csv = input_dup_csv,
    roi_json = roi_json,
    analysis_plan_json = analysis_plan_json,
    exclude_json = exclude_json,
    out_csv = file.path(tmp, "dup_out.csv"),
    out_fig_dir = file.path(tmp, "dup_figs"),
    alpha = 0.05,
    min_subjects = 6
  )
  assert_true(status_dup != 0, "expected failure on duplicate subject IDs")

  # 12) Non-numeric required predictor must fail hard.
  input_bad <- input
  input_bad$phq_total <- as.character(input_bad$phq_total)
  input_bad$phq_total[[3]] <- "not_numeric"
  input_bad_csv <- file.path(tmp, "merged_bad_predictor.csv")
  write_csv(input_bad, input_bad_csv)
  status_bad <- run_analysis(
    input_csv = input_bad_csv,
    roi_json = roi_json,
    analysis_plan_json = analysis_plan_json,
    exclude_json = exclude_json,
    out_csv = file.path(tmp, "bad_out.csv"),
    out_fig_dir = file.path(tmp, "bad_figs"),
    alpha = 0.05,
    min_subjects = 6
  )
  assert_true(status_bad != 0, "expected failure on non-numeric predictor values")

  # 13) Non-numeric required beta must fail hard.
  input_bad_beta <- input
  input_bad_beta$S04_D02_Cond04_HbR <- as.character(input_bad_beta$S04_D02_Cond04_HbR)
  input_bad_beta$S04_D02_Cond04_HbR[[6]] <- "bad_beta"
  input_bad_beta_csv <- file.path(tmp, "merged_bad_beta.csv")
  write_csv(input_bad_beta, input_bad_beta_csv)
  status_bad_beta <- run_analysis(
    input_csv = input_bad_beta_csv,
    roi_json = roi_json,
    analysis_plan_json = analysis_plan_json,
    exclude_json = exclude_json,
    out_csv = file.path(tmp, "bad_beta_out.csv"),
    out_fig_dir = file.path(tmp, "bad_beta_figs"),
    alpha = 0.05,
    min_subjects = 6
  )
  assert_true(status_bad_beta != 0, "expected failure on non-numeric beta values")

  # 14) Missing required S04_D02 target beta column must fail hard.
  input_missing_target <- input %>% select(-S04_D02_Cond04_HbR)
  input_missing_target_csv <- file.path(tmp, "merged_missing_target.csv")
  write_csv(input_missing_target, input_missing_target_csv)
  status_missing_target <- run_analysis(
    input_csv = input_missing_target_csv,
    roi_json = roi_json,
    analysis_plan_json = analysis_plan_json,
    exclude_json = exclude_json,
    out_csv = file.path(tmp, "missing_target_out.csv"),
    out_fig_dir = file.path(tmp, "missing_target_figs"),
    alpha = 0.05,
    min_subjects = 6
  )
  assert_true(status_missing_target != 0, "expected failure on missing required S04_D02 beta column")

  # 15) Malformed ROI JSON must fail hard.
  roi_bad_json <- file.path(tmp, "roi_bad.json")
  writeLines("{ bad json", roi_bad_json)
  status_roi_bad <- run_analysis(
    input_csv = input_csv,
    roi_json = roi_bad_json,
    analysis_plan_json = analysis_plan_json,
    exclude_json = exclude_json,
    out_csv = file.path(tmp, "roi_bad_out.csv"),
    out_fig_dir = file.path(tmp, "roi_bad_figs"),
    alpha = 0.05,
    min_subjects = 6
  )
  assert_true(status_roi_bad != 0, "expected failure on malformed ROI JSON")

  # 16) Overlapping ROI channel assignments must fail hard.
  roi_overlap_json <- file.path(tmp, "roi_overlap.json")
  write_roi_overlap_json(roi_overlap_json)
  status_roi_overlap <- run_analysis(
    input_csv = input_csv,
    roi_json = roi_overlap_json,
    analysis_plan_json = analysis_plan_json,
    exclude_json = exclude_json,
    out_csv = file.path(tmp, "roi_overlap_out.csv"),
    out_fig_dir = file.path(tmp, "roi_overlap_figs"),
    alpha = 0.05,
    min_subjects = 6
  )
  assert_true(status_roi_overlap != 0, "expected failure on overlapping ROI channel assignments")

  # 17) ROI definition with a missing channel must fail hard.
  roi_missing_channel_json <- file.path(tmp, "roi_missing_channel.json")
  write_roi_missing_channel_json(roi_missing_channel_json)
  status_roi_missing_channel <- run_analysis(
    input_csv = input_csv,
    roi_json = roi_missing_channel_json,
    analysis_plan_json = analysis_plan_json,
    exclude_json = exclude_json,
    out_csv = file.path(tmp, "roi_missing_channel_out.csv"),
    out_fig_dir = file.path(tmp, "roi_missing_channel_figs"),
    alpha = 0.05,
    min_subjects = 6
  )
  assert_true(status_roi_missing_channel != 0, "expected failure on ROI definition containing absent channels")

  writeLines(paste0("[OK] Correlation pipeline validation passed. Outputs in: ", tmp))
}

main()
