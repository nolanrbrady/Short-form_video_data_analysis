#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
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

expected_variables <- c(
  "sf_education_engagement",
  "sf_entertainment_engagement",
  "lf_entertainment_engagement",
  "lf_education_engagement",
  "diff_short_form_education",
  "diff_short_form_entertainment",
  "diff_long_form_education",
  "diff_long_form_entertainment",
  "age",
  "education_years",
  "sfv_frequency",
  "sfv_daily_duration",
  "asrs_total",
  "yang_pu_total",
  "yang_mot_total",
  "phq_total",
  "gad_total"
)

write_plan_json <- function(path) {
  plan_obj <- list(
    version = 1,
    description = "Validation plan for the behavior pairwise correlation matrix.",
    variables = expected_variables,
    figures = list(
      lower_triangle = list(
        filename_stem = "behavior_pairwise_correlation_lower_triangle"
      )
    )
  )
  jsonlite::write_json(plan_obj, path = path, auto_unbox = TRUE, pretty = TRUE)
}

write_variable_figure_names_json <- function(path, variables = expected_variables) {
  label_obj <- as.list(stats::setNames(
    paste("Configured figure label for", variables),
    variables
  ))
  jsonlite::write_json(label_obj, path = path, auto_unbox = TRUE, pretty = TRUE)
}

make_input <- function() {
  n <- 10
  age <- seq_len(n)
  tibble::tibble(
    subject_id = sprintf("sub_%04d", seq_len(n)),
    sf_education_engagement = age,
    sf_entertainment_engagement = 20 - age,
    lf_entertainment_engagement = c(3, 6, 2, 8, 5, 7, 4, 9, 10, 1),
    lf_education_engagement = c(2, 5, 8, 3, 6, 9, 4, 7, 10, 1),
    diff_short_form_education = age + 0.5,
    diff_short_form_entertainment = c(1, 3, 2, 4, 6, 5, 7, 8, 10, 9),
    diff_long_form_education = c(10, 8, 9, 7, 5, 6, 4, 2, 3, 1),
    diff_long_form_entertainment = c(1, 4, 2, 5, 3, 7, 6, 8, 9, 10),
    age = age,
    education_years = age * 2,
    sfv_frequency = c(NA, NA, NA, 0, 1, 2, 3, 0, 1, 2),
    sfv_daily_duration = c(0, 1, 2, NA, NA, NA, 3, 2, 1, 0),
    asrs_total = c(10, 12, 11, 14, 13, 16, 15, 18, 17, 20),
    yang_pu_total = c(21, 24, 28, 30, 35, 37, 42, 46, 49, 52),
    yang_mot_total = c(22, 26, 25, 30, 31, 35, 34, 39, 41, 43),
    phq_total = age * 2,
    gad_total = rev(age),
    pd_status = c(0, 0, 1, 0, 1, 0, 1, 1, 0, 1)
  )
}

run_script <- function(input_csv, plan_json, exclude_json, out_dir, variable_figure_names_json = NULL) {
  args <- c(
    "analyze_behavior_pairwise_correlations.R",
    "--input_csv", input_csv,
    "--analysis_plan_json", plan_json,
    "--exclude_subjects_json", exclude_json,
    "--out_dir", out_dir
  )
  if (!is.null(variable_figure_names_json)) {
    args <- c(args, "--variable_figure_names_json", variable_figure_names_json)
  }
  output <- suppressWarnings(system2(
    "Rscript",
    args,
    stdout = TRUE,
    stderr = TRUE
  ))
  list(
    status = attr(output, "status") %||% 0,
    stdout = output
  )
}

find_pair <- function(df, var_a, var_b) {
  df %>%
    filter(
      (.data$var_x == var_a & .data$var_y == var_b) |
        (.data$var_x == var_b & .data$var_y == var_a)
    )
}

normalize_subject_id_reference <- function(x, column_name) {
  x_chr <- as.character(x)
  match_pos <- regexpr("[0-9]+", x_chr)
  match_len <- attr(match_pos, "match.length")
  extracted <- rep(NA_character_, length(x_chr))
  ok <- !is.na(match_pos) & match_pos > 0
  extracted[ok] <- substring(
    x_chr[ok],
    first = match_pos[ok],
    last = match_pos[ok] + match_len[ok] - 1
  )
  if (any(is.na(extracted))) {
    bad <- unique(x_chr[is.na(extracted)])
    bad <- bad[!is.na(bad)]
    fail(paste0(
      "Reference audit could not normalize IDs from ", column_name,
      ". Examples: ", paste(head(bad, 10), collapse = ", ")
    ))
  }
  as.integer(extracted)
}

coerce_numeric_reference <- function(df, cols) {
  out <- df
  allowed_missing_tokens <- c("", "NA", "NaN", "NAN")
  for (col_name in cols) {
    raw <- out[[col_name]]
    raw_chr <- trimws(as.character(raw))
    missing_token <- is.na(raw) | raw_chr %in% allowed_missing_tokens
    suppressWarnings(num <- as.numeric(raw_chr))
    num[is.nan(num)] <- NA_real_
    bad <- !missing_token & is.na(num)
    if (any(bad)) {
      fail(paste0(
        "Reference audit found non-numeric values in ", col_name,
        ". Examples: ", paste(head(unique(raw_chr[bad]), 10), collapse = ", ")
      ))
    }
    out[[col_name]] <- num
  }
  out
}

load_plan_reference <- function(plan_json) {
  plan <- jsonlite::fromJSON(plan_json, simplifyVector = FALSE)
  variables <- as.character(unlist(plan$variables, use.names = FALSE))
  assert_true(
    identical(variables, expected_variables),
    "real-data plan variable order does not match the validated behavior-pairwise contract"
  )
  assert_true(
    !("recruitment_order_proxy" %in% variables),
    "real-data plan must exclude recruitment_order_proxy"
  )
  list(variables = variables)
}

load_reference_input <- function(input_csv, exclude_json, plan) {
  df <- read_csv(input_csv, show_col_types = FALSE)
  required_input_cols <- unique(c(plan$variables, "subject_id"))
  missing <- setdiff(required_input_cols, names(df))
  assert_true(
    length(missing) == 0,
    paste("reference audit input is missing required columns:", paste(missing, collapse = ", "))
  )

  df$subject_id <- normalize_subject_id_reference(df$subject_id, "subject_id")
  dup_ids <- sort(unique(df$subject_id[duplicated(df$subject_id)]))
  assert_true(
    length(dup_ids) == 0,
    paste("reference audit found duplicate normalized subject IDs:", paste(dup_ids, collapse = ", "))
  )
  df <- coerce_numeric_reference(df, plan$variables)

  excluded_payload <- jsonlite::fromJSON(exclude_json, simplifyVector = TRUE)
  excluded_ids <- normalize_subject_id_reference(excluded_payload, "excluded_subjects_json")
  df <- df[!(df$subject_id %in% excluded_ids), , drop = FALSE]
  df
}

compute_reference_pair <- function(df, var_x, var_y, alpha, min_subjects) {
  complete <- df[!is.na(df[[var_x]]) & !is.na(df[[var_y]]), c(var_x, var_y), drop = FALSE]
  n_complete <- nrow(complete)
  if (n_complete < min_subjects) {
    return(tibble::tibble(
      var_x = var_x,
      var_y = var_y,
      analysis_status = "skipped_min_subjects",
      skip_reason = paste0("n_complete<", min_subjects),
      n_complete = n_complete,
      pearson_r = NA_real_,
      p_unc = NA_real_,
      ci95_low = NA_real_,
      ci95_high = NA_real_
    ))
  }
  if (length(unique(complete[[var_x]])) < 2) {
    return(tibble::tibble(
      var_x = var_x,
      var_y = var_y,
      analysis_status = "skipped_constant_input",
      skip_reason = "var_x_has_zero_variance",
      n_complete = n_complete,
      pearson_r = NA_real_,
      p_unc = NA_real_,
      ci95_low = NA_real_,
      ci95_high = NA_real_
    ))
  }
  if (length(unique(complete[[var_y]])) < 2) {
    return(tibble::tibble(
      var_x = var_x,
      var_y = var_y,
      analysis_status = "skipped_constant_input",
      skip_reason = "var_y_has_zero_variance",
      n_complete = n_complete,
      pearson_r = NA_real_,
      p_unc = NA_real_,
      ci95_low = NA_real_,
      ci95_high = NA_real_
    ))
  }

  cor_fit <- stats::cor.test(
    x = complete[[var_x]],
    y = complete[[var_y]],
    method = "pearson",
    alternative = "two.sided",
    conf.level = 1 - alpha
  )
  tibble::tibble(
    var_x = var_x,
    var_y = var_y,
    analysis_status = "tested",
    skip_reason = NA_character_,
    n_complete = n_complete,
    pearson_r = unname(cor_fit$estimate),
    p_unc = cor_fit$p.value,
    ci95_low = cor_fit$conf.int[[1]],
    ci95_high = cor_fit$conf.int[[2]]
  )
}

compute_reference_results <- function(df, variables, alpha = 0.05, min_subjects = 6L) {
  pairs <- utils::combn(variables, 2)
  rows <- vector("list", ncol(pairs))
  for (i in seq_len(ncol(pairs))) {
    rows[[i]] <- compute_reference_pair(df, pairs[1, i], pairs[2, i], alpha, min_subjects)
  }
  reference <- bind_rows(rows)
  tested <- reference$analysis_status == "tested" & is.finite(reference$p_unc)
  reference$p_fdr <- NA_real_
  reference$significant_fdr <- NA
  reference$p_fdr[tested] <- p.adjust(reference$p_unc[tested], method = "BH")
  reference$significant_fdr[tested] <- reference$p_fdr[tested] < alpha
  reference
}

assert_character_columns_equal <- function(actual, expected, col_name) {
  actual_chr <- as.character(actual[[col_name]])
  expected_chr <- as.character(expected[[col_name]])
  both_na <- is.na(actual_chr) & is.na(expected_chr)
  same <- both_na | actual_chr == expected_chr
  assert_true(
    all(same),
    paste0("column ", col_name, " differs from independent reference for at least one pair")
  )
}

assert_numeric_columns_close <- function(actual, expected, col_name, tolerance = 1e-12) {
  actual_num <- actual[[col_name]]
  expected_num <- expected[[col_name]]
  both_na <- is.na(actual_num) & is.na(expected_num)
  close <- both_na | abs(actual_num - expected_num) <= tolerance
  assert_true(
    all(close),
    paste0("column ", col_name, " differs from independent reference beyond tolerance ", tolerance)
  )
}

assert_logical_columns_equal <- function(actual, expected, col_name) {
  actual_logical <- actual[[col_name]]
  expected_logical <- expected[[col_name]]
  both_na <- is.na(actual_logical) & is.na(expected_logical)
  same <- both_na | actual_logical == expected_logical
  assert_true(
    all(same),
    paste0("column ", col_name, " differs from independent reference for at least one pair")
  )
}

audit_real_dataset_against_independent_reference <- function(tmp) {
  # This audit is intentionally independent of `analyze_behavior_pairwise_correlations.R`:
  # it reads the public input/config files, reapplies exclusions, recomputes
  # every Pearson test and BH-FDR q-value, and
  # compares the exported CSVs pair-by-pair. The scientific risk is silent
  # drift in the real publication dataset that a small synthetic fixture cannot catch.
  real_out_dir <- file.path(tmp, "real_output")
  real_run <- run_script(
    input_csv = "data/tabular/generated_data/homer3_betas_plus_combined_sfv_data_inner_join.csv",
    plan_json = "data/config/behavior_pairwise_correlation_plan.json",
    exclude_json = "data/config/excluded_subjects.json",
    out_dir = real_out_dir
  )
  assert_true(
    real_run$status == 0,
    paste("real-data analysis script did not exit cleanly:", paste(real_run$stdout, collapse = "\n"))
  )

  real_out_csv <- file.path(real_out_dir, "behavior_pairwise_correlations_r.csv")
  real_fdr_csv <- file.path(real_out_dir, "behavior_pairwise_correlations_fdr_r.csv")
  assert_true(file.exists(real_out_csv), "real-data audit missing all-attempted output CSV")
  assert_true(file.exists(real_fdr_csv), "real-data audit missing all-tested FDR output CSV")

  plan <- load_plan_reference("data/config/behavior_pairwise_correlation_plan.json")
  reference_input <- load_reference_input(
    "data/tabular/generated_data/homer3_betas_plus_combined_sfv_data_inner_join.csv",
    "data/config/excluded_subjects.json",
    plan
  )
  reference <- compute_reference_results(reference_input, plan$variables)

  actual <- read_csv(real_out_csv, show_col_types = FALSE)
  actual_fdr <- read_csv(real_fdr_csv, show_col_types = FALSE)
  actual$key <- paste(actual$var_x, actual$var_y, sep = "||")
  actual_fdr$key <- paste(actual_fdr$var_x, actual_fdr$var_y, sep = "||")
  reference$key <- paste(reference$var_x, reference$var_y, sep = "||")

  assert_true(
    identical(sort(actual$key), sort(reference$key)),
    "real-data all-attempted CSV pair set differs from independent reference"
  )
  actual_sorted <- actual[match(reference$key, actual$key), , drop = FALSE]
  reference_sorted <- reference

  for (col_name in c("var_x", "var_y", "analysis_status", "skip_reason")) {
    assert_character_columns_equal(actual_sorted, reference_sorted, col_name)
  }
  for (col_name in c("n_complete", "pearson_r", "p_unc", "p_fdr", "ci95_low", "ci95_high")) {
    assert_numeric_columns_close(actual_sorted, reference_sorted, col_name)
  }
  assert_logical_columns_equal(actual_sorted, reference_sorted, "significant_fdr")

  tested_reference <- reference_sorted %>% filter(.data$analysis_status == "tested")
  assert_true(
    identical(sort(actual_fdr$key), sort(tested_reference$key)),
    "real-data FDR CSV should contain exactly every independently tested pair"
  )
  actual_fdr_sorted <- actual_fdr[match(tested_reference$key, actual_fdr$key), , drop = FALSE]
  for (col_name in c("var_x", "var_y", "analysis_status", "skip_reason")) {
    assert_character_columns_equal(actual_fdr_sorted, tested_reference, col_name)
  }
  for (col_name in c("n_complete", "pearson_r", "p_unc", "p_fdr", "ci95_low", "ci95_high")) {
    assert_numeric_columns_close(actual_fdr_sorted, tested_reference, col_name)
  }
  assert_logical_columns_equal(actual_fdr_sorted, tested_reference, "significant_fdr")
}

main <- function() {
  tmp <- file.path(tempdir(), "validate_behavior_pairwise_correlations")
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)

  input_dir <- file.path(tmp, "input")
  out_dir <- file.path(tmp, "output")
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  input_csv <- file.path(input_dir, "merged.csv")
  plan_json <- file.path(input_dir, "plan.json")
  exclude_json <- file.path(input_dir, "excluded.json")
  label_json <- file.path(input_dir, "variable_figure_names.json")
  incomplete_label_json <- file.path(input_dir, "variable_figure_names_incomplete.json")
  recruitment_order_plan_json <- file.path(input_dir, "plan_with_recruitment_order.json")
  out_csv <- file.path(out_dir, "behavior_pairwise_correlations_r.csv")
  out_fdr_csv <- file.path(out_dir, "behavior_pairwise_correlations_fdr_r.csv")
  out_sig_csv <- file.path(out_dir, "behavior_pairwise_correlations_significant_r.csv")
  out_fig_dir <- file.path(out_dir, "figures")
  out_matrix_png <- file.path(out_fig_dir, "behavior_pairwise_correlation_lower_triangle.png")
  out_matrix_pdf <- file.path(out_fig_dir, "behavior_pairwise_correlation_lower_triangle.pdf")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)
  writeLines("stale", file.path(out_dir, "stale.csv"))
  writeLines("stale", file.path(out_fig_dir, "stale.txt"))

  input_df <- make_input()
  write_csv(input_df, input_csv)
  write_plan_json(plan_json)
  recruitment_order_plan <- list(
    version = 1,
    description = "Invalid validation plan that attempts to restore recruitment order.",
    variables = c(expected_variables, "recruitment_order_proxy"),
    figures = list(
      lower_triangle = list(
        filename_stem = "behavior_pairwise_correlation_lower_triangle"
      )
    )
  )
  jsonlite::write_json(
    recruitment_order_plan,
    path = recruitment_order_plan_json,
    auto_unbox = TRUE,
    pretty = TRUE
  )
  write_variable_figure_names_json(label_json)
  write_variable_figure_names_json(incomplete_label_json, setdiff(expected_variables, "gad_total"))
  writeLines("[]", exclude_json)

  incomplete_label_run <- run_script(input_csv, plan_json, exclude_json, file.path(tmp, "incomplete_label_output"), incomplete_label_json)
  assert_true(incomplete_label_run$status != 0, "analysis script should fail when the figure-label config omits a plotted variable")
  assert_true(
    any(grepl("missing labels.*gad_total|gad_total.*missing labels", incomplete_label_run$stdout, ignore.case = TRUE)),
    "missing figure-label config failure should identify the omitted variable"
  )

  recruitment_order_run <- run_script(
    input_csv,
    recruitment_order_plan_json,
    exclude_json,
    file.path(tmp, "recruitment_order_output"),
    label_json
  )
  assert_true(
    recruitment_order_run$status != 0,
    "analysis script should fail when a plan includes recruitment_order_proxy"
  )
  assert_true(
    any(grepl("recruitment order is excluded", recruitment_order_run$stdout, ignore.case = TRUE)),
    "recruitment-order plan failure should state that recruitment order is excluded"
  )

  run <- run_script(input_csv, plan_json, exclude_json, out_dir, label_json)
  assert_true(run$status == 0, paste("analysis script did not exit cleanly:", paste(run$stdout, collapse = "\n")))
  assert_true(file.exists(out_csv), "missing all-attempted output CSV")
  assert_true(file.exists(out_fdr_csv), "missing all-tested FDR output CSV")
  assert_true(!file.exists(out_sig_csv), "legacy uncorrected significant-only CSV should not be written")
  assert_true(file.exists(out_matrix_png), "missing lower-triangle PNG")
  assert_true(file.exists(out_matrix_pdf), "missing lower-triangle PDF")
  assert_true(file.info(out_matrix_png)$size > 0, "lower-triangle PNG is empty")
  assert_true(file.info(out_matrix_pdf)$size > 0, "lower-triangle PDF is empty")
  assert_true(!file.exists(file.path(out_dir, "stale.csv")), "stale output CSV was not removed before rerun")
  assert_true(!file.exists(file.path(out_fig_dir, "stale.txt")), "stale figure artifact was not removed before rerun")

  out <- read_csv(out_csv, show_col_types = FALSE)
  fdr <- read_csv(out_fdr_csv, show_col_types = FALSE)

  required_cols <- c(
    "var_x", "var_y", "analysis_status", "skip_reason", "n_complete",
    "pearson_r", "p_unc", "p_fdr", "significant_fdr",
    "ci95_low", "ci95_high"
  )
  assert_true(all(required_cols %in% names(out)), "all-attempted CSV missing required columns")
  assert_true(all(required_cols %in% names(fdr)), "FDR CSV missing required columns")
  removed_regression_cols <- c("r_squared", "slope", "intercept", "plot_file")
  assert_true(
    !any(removed_regression_cols %in% names(out)),
    "old per-pair regression/plot columns should be removed"
  )
  assert_true(!any(out$var_x == "pd_status" | out$var_y == "pd_status"), "pd_status should be absent")
  assert_true(
    !any(out$var_x == "recruitment_order_proxy" | out$var_y == "recruitment_order_proxy"),
    "recruitment_order_proxy should be absent"
  )

  observed_variables <- unique(c(out$var_x, out$var_y))
  assert_true(
    identical(sort(observed_variables), sort(expected_variables)),
    paste("unexpected variable set:", paste(sort(observed_variables), collapse = ", "))
  )
  assert_true(nrow(out) == choose(length(expected_variables), 2), "unexpected number of pairwise rows")
  assert_true(all(fdr$analysis_status == "tested"), "FDR CSV should contain only tested rows")
  assert_true(nrow(fdr) == sum(out$analysis_status == "tested"), "FDR CSV should contain every tested pair")

  underpowered_pair <- find_pair(out, "sfv_frequency", "sfv_daily_duration")
  assert_true(nrow(underpowered_pair) == 1, "missing sfv_frequency/sfv_daily_duration row")
  assert_true(underpowered_pair$analysis_status[[1]] == "skipped_min_subjects", "underpowered pair should be skipped")
  assert_true(underpowered_pair$skip_reason[[1]] == "n_complete<6", "underpowered pair should record the min-subject skip reason")
  assert_true(underpowered_pair$n_complete[[1]] == 4, "underpowered pair should report its pairwise complete-case count")
  assert_true(is.na(underpowered_pair$p_fdr[[1]]), "underpowered skipped pair should not receive a BH-FDR q-value")
  assert_true(is.na(underpowered_pair$significant_fdr[[1]]), "underpowered skipped pair should not receive an FDR significance flag")
  assert_true(nrow(find_pair(fdr, "sfv_frequency", "sfv_daily_duration")) == 0, "FDR CSV should exclude skipped underpowered pairs")

  age_phq <- find_pair(out, "age", "phq_total")
  assert_true(nrow(age_phq) == 1, "missing age/phq_total row")
  assert_true(age_phq$analysis_status[[1]] == "tested", "age/phq_total should be tested")
  assert_true(abs(age_phq$pearson_r[[1]] - 1.0) < 1e-12, "expected perfect positive Pearson correlation for age/phq_total")

  age_education_years <- find_pair(out, "age", "education_years")
  assert_true(nrow(age_education_years) == 1, "missing age/education_years row")
  assert_true(age_education_years$analysis_status[[1]] == "tested", "age/education_years should be tested")
  assert_true(
    abs(age_education_years$pearson_r[[1]] - 1.0) < 1e-12,
    "expected perfect positive Pearson correlation for age/education_years in the synthetic proxy fixture"
  )

  age_gad <- find_pair(out, "age", "gad_total")
  assert_true(nrow(age_gad) == 1, "missing age/gad_total row")
  assert_true(abs(age_gad$pearson_r[[1]] + 1.0) < 1e-12, "expected perfect negative Pearson correlation for age/gad_total")

  age_frequency <- find_pair(out, "age", "sfv_frequency")
  complete_frequency <- input_df %>% filter(!is.na(.data$age), !is.na(.data$sfv_frequency))
  expected_frequency_r <- stats::cor(complete_frequency$age, complete_frequency$sfv_frequency, method = "pearson")
  assert_true(age_frequency$n_complete[[1]] == nrow(complete_frequency), "pairwise complete cases for age/sfv_frequency are incorrect")
  assert_true(
    abs(age_frequency$pearson_r[[1]] - expected_frequency_r) < 1e-12,
    "sfv_frequency should be treated as numeric Pearson input in this workflow"
  )

  tested <- out %>% filter(.data$analysis_status == "tested", is.finite(.data$p_unc))
  expected_q <- p.adjust(tested$p_unc, method = "BH")
  assert_true(
    max(abs(tested$p_fdr - expected_q), na.rm = TRUE) < 1e-12,
    "p_fdr does not match global Benjamini-Hochberg correction"
  )
  assert_true(
    identical(tested$significant_fdr, tested$p_fdr < 0.05),
    "significant_fdr should be TRUE exactly when p_fdr < alpha"
  )

  audit_real_dataset_against_independent_reference(tmp)

  writeLines(paste0("[OK] Behavior pairwise correlation validation passed. Outputs in: ", tmp))
}

main()
