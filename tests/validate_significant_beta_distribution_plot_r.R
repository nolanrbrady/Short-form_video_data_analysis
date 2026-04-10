#!/usr/bin/env Rscript

# Validation harness for plot_significant_beta_value_distribution.R.
#
# Purpose
# - Verify that significant-hit discovery uses the tidy LMM outputs with p_fdr < alpha.
# - Verify subject exclusions and channel/ROI complete-case filtering match the
#   inferential sample logic.
# - Verify ROI aggregation uses the mean across available non-missing member
#   channels.
# - Verify main-effect plots use subject-level marginal means across the
#   orthogonal factor.
# - Verify interaction plots retain the 4 raw conditions rather than collapsing.
# - Verify literal zero beta placeholders fail hard.

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

script_path_from_args <- function() {
  args_all <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args_all, value = TRUE)
  if (length(file_arg) == 0) {
    stop("Could not resolve current script path from commandArgs().")
  }
  normalizePath(sub("^--file=", "", file_arg[[1]]))
}

ROOT <- dirname(dirname(script_path_from_args()))
SCRIPT_PATH <- file.path(ROOT, "plot_significant_beta_value_distribution.R")
RSCRIPT_BIN <- "/Library/Frameworks/R.framework/Resources/bin/Rscript"

plot_env <- new.env(parent = globalenv())
sys.source(SCRIPT_PATH, envir = plot_env)

write_toy_merged_csv <- function(path) {
  df <- tibble::tibble(
    subject_id = c("0001", "0002", "0003", "0004"),
    age = c(20, 21, 22, 23),
    S07_D07_Cond01_HbR = c(1.0, 2.0, 3.0, 9.0),
    S07_D07_Cond02_HbR = c(3.0, 4.0, 5.0, 9.0),
    S07_D07_Cond03_HbR = c(5.0, 6.0, 7.0, 9.0),
    S07_D07_Cond04_HbR = c(7.0, 8.0, NA_real_, 9.0),
    S04_D02_Cond01_HbR = c(-5.0, -4.0, -3.0, -9.0),
    S04_D02_Cond02_HbR = c(-7.0, -6.0, -4.0, -9.0),
    S04_D02_Cond03_HbR = c(-1.0, -2.0, -2.0, -9.0),
    S04_D02_Cond04_HbR = c(-3.0, -3.5, NA_real_, -9.0),
    S03_D02_Cond01_HbO = c(1.0, 2.0, 3.0, 9.0),
    S03_D02_Cond02_HbO = c(2.0, NA_real_, 4.0, 9.0),
    S03_D02_Cond03_HbO = c(3.0, 4.0, 5.0, 9.0),
    S03_D02_Cond04_HbO = c(4.0, 5.0, 6.0, 9.0),
    S03_D04_Cond01_HbO = c(2.0, 4.0, 5.0, 9.0),
    S03_D04_Cond02_HbO = c(3.0, 4.0, 6.0, 9.0),
    S03_D04_Cond03_HbO = c(4.0, 6.0, 7.0, 9.0),
    S03_D04_Cond04_HbO = c(5.0, 7.0, 8.0, 9.0)
  )
  write_csv(df, path, na = "")
}

write_toy_roi_json <- function(path) {
  payload <- list(
    L_DMPFC = c("S03_D02", "S03_D04"),
    R_DLPFC = c("S07_D07")
  )
  writeLines(toJSON(payload, auto_unbox = TRUE), path)
}

write_toy_exclusions <- function(path) {
  writeLines('["sub_0004"]', path)
}

write_toy_channel_results <- function(path) {
  df <- tibble::tibble(
    channel = c("S04_D02", "S07_D07", "S01_D01"),
    chrom = c("HbR", "HbR", "HbO"),
    effect = c("interaction", "format", "content"),
    n_subjects = c(2L, 2L, 4L),
    n_obs = c(8L, 8L, 16L),
    converged = c(TRUE, TRUE, TRUE),
    singular_fit = c(FALSE, FALSE, FALSE),
    estimate = c(-0.00001, -0.000005, 0.000001),
    ci95_low = c(-0.00002, -0.000009, -0.000002),
    ci95_high = c(-0.000003, -0.000001, 0.000004),
    p_unc = c(0.001, 0.002, 0.2),
    p_fdr = c(0.02, 0.04, 0.2)
  )
  write_csv(df, path)
}

write_toy_roi_results <- function(path) {
  df <- tibble::tibble(
    roi = c("L_DMPFC", "L_DMPFC"),
    chrom = c("HbO", "HbR"),
    effect = c("format", "interaction"),
    n_subjects = c(3L, 3L),
    n_obs = c(12L, 12L),
    converged = c(TRUE, TRUE),
    singular_fit = c(FALSE, FALSE),
    estimate = c(0.00001, -0.000002),
    ci95_low = c(0.000002, -0.000008),
    ci95_high = c(0.000018, 0.000003),
    p_unc = c(0.01, 0.3),
    p_fdr = c(0.03, 0.3)
  )
  write_csv(df, path)
}

test_significant_hit_discovery <- function(channel_results, roi_results) {
  sig_channel <- plot_env$load_significant_hits(channel_results, "channel", 0.05)
  sig_roi <- plot_env$load_significant_hits(roi_results, "roi", 0.05)

  assert_true(nrow(sig_channel) == 2, "Expected exactly 2 significant channel hits.")
  assert_true(
    identical(sig_channel$channel, c("S04_D02", "S07_D07")),
    "Channel hits should be sorted by p_fdr and include only FDR-significant rows."
  )
  assert_true(nrow(sig_roi) == 1, "Expected exactly 1 significant ROI hit.")
  assert_true(sig_roi$roi[[1]] == "L_DMPFC", "ROI significant hit should be L_DMPFC.")
}

test_zero_placeholder_failure <- function(merged_csv) {
  zero_csv <- tempfile(fileext = ".csv")
  df <- read_csv(merged_csv, show_col_types = FALSE)
  df$S07_D07_Cond01_HbR[[1]] <- 0
  write_csv(df, zero_csv, na = "")

  ok <- FALSE
  msg <- NULL
  tryCatch(
    {
      plot_env$load_merged_input(zero_csv)
    },
    error = function(e) {
      ok <<- TRUE
      msg <<- conditionMessage(e)
    }
  )
  assert_true(ok, "Expected literal beta-value zeros to fail hard.")
  assert_true(grepl("literal 0 values", msg, fixed = TRUE), "Zero-placeholder failure message should mention literal 0 values.")
}

test_run_plotting_outputs <- function(
  merged_csv,
  roi_json,
  exclusions_json,
  channel_results,
  roi_results,
  out_dir
) {
  outputs <- plot_env$run_plotting(
    input_csv = merged_csv,
    roi_json = roi_json,
    exclude_subjects_json = exclusions_json,
    channel_results_tidy_csv = channel_results,
    roi_results_tidy_csv = roi_results,
    alpha = 0.05,
    out_dir = out_dir
  )

  expected_files <- file.path(
    out_dir,
    c(
      "channel_S04_D02_HbR_interaction_beta_distribution.png",
      "channel_S07_D07_HbR_format_beta_distribution.png",
      "roi_L_DMPFC_HbO_format_beta_distribution.png"
    )
  )

  assert_true(length(outputs$figure_paths) == 3, "Expected 3 output figures.")
  assert_true(all(file.exists(expected_files)), "Expected named PNG outputs were not created.")
  assert_true(file.exists(outputs$audit_csv_path), "Expected audit CSV was not created.")

  audit <- outputs$audit_df
  expected_plot_keys <- tibble::tribble(
    ~analysis_level, ~unit_id, ~chrom, ~effect, ~plot_mode,
    "channel", "S04_D02", "HbR", "interaction", "interaction_conditions",
    "channel", "S07_D07", "HbR", "format", "main_effect_marginal",
    "roi", "L_DMPFC", "HbO", "format", "main_effect_marginal"
  )
  actual_plot_keys <- audit %>%
    distinct(analysis_level, unit_id, chrom, effect, plot_mode) %>%
    arrange(analysis_level, unit_id, chrom, effect, plot_mode)
  assert_true(
    identical(actual_plot_keys, arrange(expected_plot_keys, analysis_level, unit_id, chrom, effect, plot_mode)),
    "Audit CSV should contain exactly the expected plotted analysis-level/unit/chrom/effect combinations."
  )

  channel_format <- audit %>%
    filter(
      analysis_level == "channel",
      unit_id == "S07_D07",
      chrom == "HbR",
      effect == "format"
    )
  assert_true(
    identical(sort(unique(channel_format$subject_id)), c(1L, 2L)),
    "Channelwise main-effect plot should keep only complete-case, non-excluded subjects 1 and 2."
  )
  assert_true(nrow(channel_format) == 8, "Channelwise format audit rows should contain 4 raw rows per retained subject.")
  sub2_short <- unique(channel_format$beta_plot_value[channel_format$subject_id == 2 & channel_format$format == "Short"])
  sub2_long <- unique(channel_format$beta_plot_value[channel_format$subject_id == 2 & channel_format$format == "Long"])
  assert_true(length(sub2_short) == 1 && abs(sub2_short[[1]] - 3.0) < 1e-9, "Subject 2 short marginal mean should equal 3.0.")
  assert_true(length(sub2_long) == 1 && abs(sub2_long[[1]] - 7.0) < 1e-9, "Subject 2 long marginal mean should equal 7.0.")

  channel_interaction <- audit %>%
    filter(
      analysis_level == "channel",
      unit_id == "S04_D02",
      chrom == "HbR",
      effect == "interaction"
    )
  assert_true(
    identical(sort(unique(channel_interaction$condition)), c("LF_Edu", "LF_Ent", "SF_Edu", "SF_Ent")),
    "Interaction audit rows should retain all four conditions."
  )
  assert_true(
    max(abs(channel_interaction$beta_plot_value - channel_interaction$beta_raw)) < 1e-12,
    "Interaction beta_plot_value should equal beta_raw for all rows."
  )

  roi_format <- audit %>%
    filter(
      analysis_level == "roi",
      unit_id == "L_DMPFC",
      chrom == "HbO",
      effect == "format"
    )
  assert_true(
    identical(sort(unique(roi_format$subject_id)), c(1L, 2L, 3L)),
    "ROI plot should retain subjects 1-3 after exclusions."
  )
  assert_true(nrow(roi_format) == 12, "ROI format audit rows should contain 4 raw rows per retained subject.")
  roi_sub2_short <- unique(roi_format$beta_plot_value[roi_format$subject_id == 2 & roi_format$format == "Short"])
  roi_sub2_long <- unique(roi_format$beta_plot_value[roi_format$subject_id == 2 & roi_format$format == "Long"])
  assert_true(length(roi_sub2_short) == 1 && abs(roi_sub2_short[[1]] - 3.5) < 1e-9, "Subject 2 ROI short marginal mean should equal 3.5.")
  assert_true(length(roi_sub2_long) == 1 && abs(roi_sub2_long[[1]] - 5.5) < 1e-9, "Subject 2 ROI long marginal mean should equal 5.5.")
}

test_cli_execution <- function(
  merged_csv,
  roi_json,
  exclusions_json,
  channel_results,
  roi_results,
  out_dir
) {
  audit_csv <- file.path(out_dir, "plotted_beta_values.csv")
  args <- c(
    shQuote(SCRIPT_PATH),
    "--input_csv", shQuote(merged_csv),
    "--roi_json", shQuote(roi_json),
    "--exclude_subjects_json", shQuote(exclusions_json),
    "--channel_results_tidy_csv", shQuote(channel_results),
    "--roi_results_tidy_csv", shQuote(roi_results),
    "--alpha", "0.05",
    "--out_dir", shQuote(out_dir)
  )
  output <- suppressWarnings(system2(RSCRIPT_BIN, args, stdout = TRUE, stderr = TRUE))
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L
  assert_true(status == 0L, paste("CLI execution failed:", paste(output, collapse = "\n")))
  assert_true(file.exists(audit_csv), "CLI execution should write the audit CSV.")
}

main <- function() {
  tmp <- file.path(tempdir(), "significant_beta_distribution_plot_validation")
  dir.create(tmp, showWarnings = FALSE, recursive = TRUE)

  merged_csv <- file.path(tmp, "merged.csv")
  roi_json <- file.path(tmp, "roi_definition.json")
  exclusions_json <- file.path(tmp, "excluded_subjects.json")
  channel_results <- file.path(tmp, "channel_results.csv")
  roi_results <- file.path(tmp, "roi_results.csv")
  out_dir <- file.path(tmp, "out")
  cli_out_dir <- file.path(tmp, "out_cli")

  write_toy_merged_csv(merged_csv)
  write_toy_roi_json(roi_json)
  write_toy_exclusions(exclusions_json)
  write_toy_channel_results(channel_results)
  write_toy_roi_results(roi_results)

  test_significant_hit_discovery(channel_results, roi_results)
  test_zero_placeholder_failure(merged_csv)
  test_run_plotting_outputs(
    merged_csv = merged_csv,
    roi_json = roi_json,
    exclusions_json = exclusions_json,
    channel_results = channel_results,
    roi_results = roi_results,
    out_dir = out_dir
  )
  test_cli_execution(
    merged_csv = merged_csv,
    roi_json = roi_json,
    exclusions_json = exclusions_json,
    channel_results = channel_results,
    roi_results = roi_results,
    out_dir = cli_out_dir
  )

  writeLines("[PASS] validate_significant_beta_distribution_plot_r")
}

main()
