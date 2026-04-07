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

write_plan_json <- function(path) {
  plan_obj <- list(
    version = 1,
    variables = c(
      "age",
      "phq_total",
      "gad_total",
      "pd_status",
      "sfv_frequency",
      "yang_pu_total"
    ),
    figures = list(
      policy = "significant_only"
    )
  )
  jsonlite::write_json(plan_obj, path = path, auto_unbox = TRUE, pretty = TRUE)
}

make_input <- function() {
  tibble::tibble(
    subject_id = sprintf("%04d", 1:8),
    age = c(1, 2, 3, 4, 5, 6, 7, 8),
    phq_total = c(2, 4, 6, 8, 10, 12, 14, 16),
    gad_total = c(8, 7, 6, 5, 4, 3, 2, 1),
    pd_status = c(0, 0, 0, 0, 1, 1, 1, 1),
    sfv_frequency = c(1, 2, 3, 4, 5, 6, 7, NA),
    yang_pu_total = rep(5, 8)
  )
}

run_script <- function(input_csv, plan_json, exclude_json, out_dir) {
  output <- suppressWarnings(system2(
    "Rscript",
    c(
      "analyze_behavior_pairwise_correlations.R",
      "--input_csv", input_csv,
      "--analysis_plan_json", plan_json,
      "--exclude_subjects_json", exclude_json,
      "--out_dir", out_dir
    ),
    stdout = TRUE,
    stderr = TRUE
  ))
  list(
    status = attr(output, "status") %||% 0,
    stdout = output
  )
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
  out_csv <- file.path(out_dir, "behavior_pairwise_correlations_r.csv")
  out_sig_csv <- file.path(out_dir, "behavior_pairwise_correlations_significant_r.csv")
  out_fig_dir <- file.path(out_dir, "figures")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)
  writeLines("stale", file.path(out_dir, "stale.csv"))
  writeLines("stale", file.path(out_fig_dir, "stale.txt"))

  write_csv(make_input(), input_csv)
  write_plan_json(plan_json)
  writeLines("[]", exclude_json)

  run <- run_script(input_csv, plan_json, exclude_json, out_dir)
  assert_true(run$status == 0, paste("analysis script did not exit cleanly:", paste(run$stdout, collapse = "\n")))
  assert_true(file.exists(out_csv), "missing output CSV")
  assert_true(file.exists(out_sig_csv), "missing significant-only output CSV")
  assert_true(!file.exists(file.path(out_dir, "stale.csv")), "stale output CSV was not removed before rerun")
  assert_true(!file.exists(file.path(out_fig_dir, "stale.txt")), "stale figure artifact was not removed before rerun")

  out <- read_csv(out_csv, show_col_types = FALSE)
  sig <- read_csv(out_sig_csv, show_col_types = FALSE)

  required_cols <- c(
    "var_x", "var_y", "analysis_status", "skip_reason", "n_complete",
    "pearson_r", "r_squared", "p_unc", "ci95_low", "ci95_high",
    "slope", "intercept", "plot_file"
  )
  assert_true(all(required_cols %in% names(out)), "output CSV missing required columns")
  assert_true(nrow(out) == choose(6, 2), "unexpected number of pairwise rows")

  age_phq <- out %>% filter(.data$var_x == "age", .data$var_y == "phq_total")
  assert_true(nrow(age_phq) == 1, "missing age/phq_total row")
  assert_true(age_phq$analysis_status[[1]] == "tested", "age/phq_total should be tested")
  assert_true(abs(age_phq$pearson_r[[1]] - 1.0) < 1e-12, "expected perfect positive correlation for age/phq_total")

  age_gad <- out %>% filter(.data$var_x == "age", .data$var_y == "gad_total")
  assert_true(nrow(age_gad) == 1, "missing age/gad_total row")
  assert_true(abs(age_gad$pearson_r[[1]] + 1.0) < 1e-12, "expected perfect negative correlation for age/gad_total")

  age_freq <- out %>% filter(.data$var_x == "age", .data$var_y == "sfv_frequency")
  assert_true(age_freq$n_complete[[1]] == 7, "pairwise complete-case handling for age/sfv_frequency is incorrect")

  const_row <- out %>% filter(.data$var_x == "age", .data$var_y == "yang_pu_total")
  assert_true(const_row$analysis_status[[1]] == "skipped_constant_input", "constant variable row should be skipped")
  assert_true(const_row$skip_reason[[1]] == "var_y_has_zero_variance", "unexpected skip reason for constant variable row")

  binary_row <- out %>% filter(.data$var_x == "age", .data$var_y == "pd_status")
  assert_true(binary_row$analysis_status[[1]] == "tested", "binary/continuous pair should still be tested")
  assert_true(is.finite(binary_row$pearson_r[[1]]), "binary/continuous Pearson result should be finite")

  assert_true(all(sig$analysis_status == "tested"), "significant CSV should only contain tested rows")
  assert_true(all(sig$p_unc < 0.05), "significant CSV should only contain uncorrected-significant rows")
  assert_true(any(!is.na(sig$plot_file)), "significant rows should produce plot files")
  assert_true(all(file.exists(sig$plot_file[!is.na(sig$plot_file)])), "expected significant plot files to exist")

  writeLines(paste0("[OK] Behavior pairwise correlation validation passed. Outputs in: ", tmp))
}

main()
