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

write_roi_json <- function(path) {
  roi_obj <- list(
    R_DLPFC = c("S07_D05", "S07_D07"),
    L_DLPFC = c("S01_D01", "S02_D01"),
    M_DMPFC = c("S04_D04"),
    L_DMPFC = c("S03_D02", "S03_D04")
  )
  jsonlite::write_json(roi_obj, path = path, auto_unbox = TRUE, pretty = TRUE)
}

write_plan_json <- function(path) {
  plan_obj <- list(
    predictors = c(
      "lf_education_engagement",
      "lf_entertainment_engagement",
      "sf_education_engagement",
      "sf_entertainment_engagement",
      "long_form_engagement",
      "short_form_engagement",
      "education_engagement",
      "entertainment_engagement",
      "diff_short_form_education",
      "diff_short_form_entertainment",
      "diff_long_form_education",
      "diff_long_form_entertainment"
    ),
    multiple_testing = list(
      adjust_method = "BH",
      family_grouping = c("neural_level", "neural_name", "chrom", "predictor")
    ),
    figures = list(
      policy = "significant_only"
    )
  )
  jsonlite::write_json(plan_obj, path = path, auto_unbox = TRUE, pretty = TRUE)
}

make_input <- function() {
  n <- 8
  long_signal <- c(1, 3, 2, 8, 5, 7, 4, 6)
  short_signal <- c(6, 2, 7, 1, 8, 3, 5, 4)
  edu_signal <- c(4, 8, 1, 6, 2, 7, 3, 5)
  ent_signal <- c(5, 1, 6, 3, 7, 2, 8, 4)
  df <- tibble::tibble(
    subject_id = sprintf("%04d", seq_len(n)),
    age = c(22, 31, 27, 35, 24, 33, 29, 26),
    lf_education_engagement = edu_signal + c(0, 1, 0, 2, 1, 0, 2, 1),
    lf_entertainment_engagement = ent_signal + c(2, 0, 1, 0, 2, 1, 0, 1),
    sf_education_engagement = short_signal + c(1, 0, 2, 1, 0, 2, 1, 0),
    sf_entertainment_engagement = rev(ent_signal) + c(0, 2, 1, 0, 1, 0, 2, 1),
    long_form_engagement = long_signal,
    short_form_engagement = short_signal,
    education_engagement = edu_signal,
    entertainment_engagement = ent_signal,
    diff_short_form_education = c(3, 9, 2, 7, 1, 8, 4, 6),
    diff_short_form_entertainment = c(8, 1, 6, 4, 7, 2, 5, 3),
    diff_long_form_education = c(2, 5, 1, 8, 4, 7, 3, 6),
    diff_long_form_entertainment = c(7, 3, 8, 2, 6, 1, 5, 4)
  )

  add_col <- function(name, values) {
    df[[name]] <<- values
  }

  conds <- sprintf("%02d", 1:4)

  for (cond in conds) {
    add_col(paste0("S04_D02_Cond", cond, "_HbO"), df$long_form_engagement * 2 + as.integer(cond))
    add_col(paste0("S04_D02_Cond", cond, "_HbR"), rev(df$short_form_engagement) + as.integer(cond))
  }

  for (channel in c("S01_D01", "S02_D01")) {
    for (cond in conds) {
      add_col(paste0(channel, "_Cond", cond, "_HbO"), df$lf_education_engagement + as.integer(cond))
    }
  }
  for (cond in conds) {
    add_col(paste0("S04_D04_Cond", cond, "_HbO"), rep(10 + as.integer(cond), n))
  }
  for (channel in c("S03_D02", "S03_D04")) {
    for (cond in conds) {
      add_col(paste0(channel, "_Cond", cond, "_HbO"), df$diff_long_form_education + as.integer(cond))
    }
  }
  for (channel in c("S07_D05", "S07_D07")) {
    for (cond in conds) {
      add_col(paste0(channel, "_Cond", cond, "_HbR"), df$sf_entertainment_engagement + as.integer(cond))
    }
  }

  df
}

run_script <- function(input_csv, roi_json, plan_json, exclude_json, out_csv, out_fig_dir) {
  output <- suppressWarnings(system2(
    "Rscript",
    c(
      "analyze_correlational_relationships.R",
      "--input_csv", input_csv,
      "--roi_json", roi_json,
      "--analysis_plan_json", plan_json,
      "--exclude_subjects_json", exclude_json,
      "--out_csv", out_csv,
      "--out_fig_dir", out_fig_dir
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
  tmp <- file.path(tempdir(), "validate_correlational_relationships_old")
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)

  input_dir <- file.path(tmp, "input")
  out_dir <- file.path(tmp, "output")
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  input_csv <- file.path(input_dir, "merged.csv")
  roi_json <- file.path(input_dir, "roi.json")
  plan_json <- file.path(input_dir, "plan.json")
  exclude_json <- file.path(input_dir, "excluded.json")
  out_csv <- file.path(out_dir, "pairwise_correlations_r.csv")
  out_fig_dir <- file.path(out_dir, "figures")
  dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
  dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)
  writeLines("stale", file.path(dirname(out_csv), "stale.csv"))
  writeLines("stale", file.path(out_fig_dir, "stale.txt"))

  write_csv(make_input(), input_csv)
  write_roi_json(roi_json)
  write_plan_json(plan_json)
  writeLines("[]", exclude_json)

  run <- run_script(input_csv, roi_json, plan_json, exclude_json, out_csv, out_fig_dir)
  assert_true(run$status == 0, "analysis script did not exit cleanly")
  assert_true(file.exists(out_csv), "missing output CSV")
  assert_true(!file.exists(file.path(dirname(out_csv), "stale.csv")), "stale output CSV was not removed before rerun")
  assert_true(!file.exists(file.path(out_fig_dir, "stale.txt")), "stale figure artifact was not removed before rerun")

  out <- read_csv(out_csv, show_col_types = FALSE)
  required_cols <- c(
    "neural_level", "neural_name", "chrom", "condition", "condition_code", "condition_label",
    "predictor", "analysis_status", "skip_reason", "n_complete", "pearson_r", "r_squared",
    "p_unc", "ci95_low", "ci95_high", "slope", "intercept", "plot_file", "family_id",
    "family_adjust_method", "family_n_tested", "p_fdr"
  )
  assert_true(all(required_cols %in% names(out)), "output CSV missing required columns")
  assert_true(any(out$analysis_status == "tested"), "synthetic fixture should yield tested rows")
  assert_true(all(out$family_adjust_method == "BH"), "expected BH family adjustment")

  expected_rows <- 12 * (2 * 4 + 4 * 4)
  assert_true(nrow(out) == expected_rows, "unexpected number of predictor x target x condition rows")

  row_check <- out %>%
    filter(
      neural_level == "channel",
      neural_name == "S04_D02",
      chrom == "HbO",
      condition == "Cond01",
      predictor == "long_form_engagement"
    )
  assert_true(nrow(row_check) == 1, "missing known positive channel row")
  assert_true(abs(row_check$pearson_r[[1]] - 1.0) < 1e-12, "expected perfect positive Pearson correlation")
  assert_true(row_check$n_complete[[1]] == 8, "expected all synthetic subjects to be complete")
  assert_true(row_check$analysis_status[[1]] == "tested", "known positive row should be tested")

  family_rows <- out %>%
    filter(
      neural_level == "channel",
      neural_name == "S04_D02",
      chrom == "HbO",
      predictor == "long_form_engagement"
    )
  assert_true(nrow(family_rows) == 4, "expected one four-condition family for each predictor x neural target")
  assert_true(all(family_rows$family_n_tested == 4), "expected BH family size of four conditions")

  writeLines(paste0("[OK] Correlation pipeline validation passed. Outputs in: ", tmp))
}

main()
