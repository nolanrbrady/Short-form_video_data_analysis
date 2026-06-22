#!/usr/bin/env Rscript

# Validation harness for plot_behavior_score_distributions.R.
#
# Purpose
# - Verify shared subject exclusions are applied before plotting.
# - Verify engagement and recall/retention condition mappings.
# - Verify zero scores remain valid observations.
# - Verify non-numeric values and missing required columns fail hard.
# - Verify complete-case filtering is domain-specific and does not impute.
# - Verify plotted subject sets and condition mappings match the engagement and
#   retention LMM preprocessing functions.
# - Verify content marginal contrasts match LMM content estimates in controlled
#   no-age-effect/no-noise synthetic data where equality is expected.
# - Verify retention length marginal contrasts match LMM length estimates in the same
#   controlled synthetic setting.
# - Verify PNG figure and audit outputs are created.
# - Verify the figure source uses violin + jittered raw points with mean +/- 1 SD overlays.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(jsonlite)
  library(lmerTest)
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
SCRIPT_PATH <- file.path(ROOT, "plot_behavior_score_distributions.R")
RSCRIPT_BIN <- file.path(R.home("bin"), "Rscript")

plot_env <- new.env(parent = globalenv())
sys.source(SCRIPT_PATH, envir = plot_env)

# Source analysis scripts without executing their CLI `main()` calls. The
# validation harness needs the production LMM preprocessing/model functions for
# direct comparison, but it must not rewrite production result CSVs as a side
# effect of sourcing.
source_without_terminal_main <- function(path, envir) {
  lines <- readLines(path, warn = FALSE)
  main_calls <- which(trimws(lines) == "main()")
  if (length(main_calls) == 0) {
    stop(paste0("Expected terminal main() call in script: ", path))
  }
  terminal_main <- tail(main_calls, 1)
  lines <- lines[-terminal_main]
  eval(parse(text = paste(lines, collapse = "\n")), envir = envir)
}

engagement_lmm_env <- new.env(parent = globalenv())
source_without_terminal_main(file.path(ROOT, "analyze_engagement_format_content_lmm.R"), engagement_lmm_env)

retention_lmm_env <- new.env(parent = globalenv())
source_without_terminal_main(file.path(ROOT, "analyze_retention_format_content_lmm.R"), retention_lmm_env)

# Create a compact final-merged-data fixture with:
# - one excluded subject,
# - one engagement-incomplete subject,
# - one retention-incomplete subject,
# - valid zero values in both domains.
# This single fixture pressure-tests exclusion, complete-case filtering, and
# zero-as-observed behavior without needing production data.
write_toy_input <- function(path) {
  df <- tibble::tibble(
    subject_id = c("0001", "0002", "0003", "0004", "0005"),
    age = c(18, 19, 20, 21, 22),
    sf_education_engagement = c(1.0, 0.0, 3.0, NA_real_, 9.0),
    sf_entertainment_engagement = c(2.0, 1.0, 4.0, 4.0, 9.0),
    lf_education_engagement = c(3.0, 2.0, 5.0, 5.0, 9.0),
    lf_entertainment_engagement = c(4.0, 3.0, 6.0, 6.0, 9.0),
    diff_short_form_education = c(0.10, 0.00, NA_real_, 0.40, 0.90),
    diff_short_form_entertainment = c(0.20, 0.10, 0.30, 0.50, 0.90),
    diff_long_form_education = c(0.30, 0.20, 0.40, 0.60, 0.90),
    diff_long_form_entertainment = c(0.40, 0.30, 0.50, 0.70, 0.90),
    S01_D01_Cond01_HbO = c(1, 2, 3, 4, 5)
  )
  write_csv(df, path, na = "")
}

# Write an exclusion manifest in the same Homer-style string format used by the
# real repo config. The shared exclusion helper normalizes this to integer ID 5.
write_toy_exclusions <- function(path) {
  writeLines('["sub_0005"]', path)
}

# Test that the plotting script still uses the intended publication geometry:
# violin density envelopes, raw points, and mean +/- 1 SD overlays. This guards
# against accidental replacement with summary-only bars or connected subject
# lines, which would weaken distribution transparency.
test_plot_geometry_source <- function() {
  script_src <- paste(readLines(SCRIPT_PATH, warn = FALSE), collapse = "\n")
  assert_true(grepl("geom_violin", script_src, fixed = TRUE), "Behavior plot should use violin geometry.")
  assert_true(grepl("geom_point", script_src, fixed = TRUE), "Behavior plot should show raw jittered subject points.")
  assert_true(grepl("mean_sdl_1", script_src, fixed = TRUE), "Behavior plot should define a mean +/- 1 SD helper.")
  assert_true(grepl('geom = "errorbar"', script_src, fixed = TRUE), "Behavior plot should show +/- 1 SD intervals.")
  assert_true(grepl('geom = "point"', script_src, fixed = TRUE), "Behavior plot should show mean markers.")
}

# End-to-end plotting smoke/semantic test on the toy fixture. It verifies output
# artifacts, exclusion behavior, domain-specific complete-case membership,
# condition ordering, zero preservation, exact content marginal mean arithmetic,
# and summary CSV row counts.
test_run_plotting_outputs <- function(input_csv, exclusions_json, out_dir) {
  outputs <- plot_env$run_plotting(
    input_csv = input_csv,
    exclude_subjects_json = exclusions_json,
    out_dir = out_dir,
    width = 6,
    height = 4,
    dpi = 120
  )

  expected_figures <- file.path(
    out_dir,
    c(
      "engagement_score_distribution.png",
      "retention_score_distribution.png",
      "engagement_content_marginal_score_distribution.png",
      "retention_content_marginal_score_distribution.png",
      "retention_length_marginal_score_distribution.png"
    )
  )
  assert_true(all(file.exists(expected_figures)), "Expected PNG behavior figures were not created.")
  assert_true(file.exists(outputs$audit_csv_path), "Expected plotted behavior audit CSV was not created.")
  assert_true(file.exists(outputs$summary_csv_path), "Expected behavior summary CSV was not created.")

  audit <- outputs$audit_df
  assert_true(!any(audit$subject_id == 5L), "Excluded subject 5 should not appear in plotted audit rows.")

  engagement <- audit %>% filter(.data$plot_type == "raw_condition", .data$domain == "engagement")
  retention <- audit %>% filter(.data$plot_type == "raw_condition", .data$domain == "retention")
  assert_true(
    identical(sort(unique(engagement$subject_id)), c(1L, 2L, 3L)),
    "Engagement plot should keep complete-case non-excluded subjects 1, 2, and 3."
  )
  assert_true(
    identical(sort(unique(retention$subject_id)), c(1L, 2L, 4L)),
    "Retention plot should keep complete-case non-excluded subjects 1, 2, and 4."
  )
  assert_true(nrow(engagement) == 12, "Engagement audit should contain 4 rows per retained subject.")
  assert_true(nrow(retention) == 12, "Retention audit should contain 4 rows per retained subject.")

  observed_conditions <- unique(engagement$condition)
  assert_true(
    identical(observed_conditions[1:4], c("SF_Edu", "SF_Ent", "LF_Edu", "LF_Ent")),
    "Condition order should match the requested four-condition display order."
  )
  assert_true(
    any(engagement$subject_id == 2L & engagement$condition == "SF_Edu" & engagement$score == 0),
    "Engagement zero values should remain valid plotted scores."
  )
  assert_true(
    any(retention$subject_id == 2L & retention$condition == "SF_Edu" & retention$score == 0),
    "Retention zero values should remain valid plotted scores."
  )

  engagement_content <- audit %>% filter(.data$plot_type == "content_marginal", .data$domain == "engagement")
  retention_content <- audit %>% filter(.data$plot_type == "content_marginal", .data$domain == "retention")
  assert_true(nrow(engagement_content) == 6, "Engagement content marginal audit should contain 2 rows per retained subject.")
  assert_true(nrow(retention_content) == 6, "Retention content marginal audit should contain 2 rows per retained subject.")
  assert_true(
    all(engagement_content$n_source_conditions == 2L) && all(retention_content$n_source_conditions == 2L),
    "Content marginal rows should be computed from exactly two source conditions."
  )
  sub1_edu_engagement <- engagement_content$score[
    engagement_content$subject_id == 1L & engagement_content$content_level == "Education"
  ]
  sub1_ent_engagement <- engagement_content$score[
    engagement_content$subject_id == 1L & engagement_content$content_level == "Entertainment"
  ]
  sub2_edu_retention <- retention_content$score[
    retention_content$subject_id == 2L & retention_content$content_level == "Education"
  ]
  sub2_ent_retention <- retention_content$score[
    retention_content$subject_id == 2L & retention_content$content_level == "Entertainment"
  ]
  assert_true(length(sub1_edu_engagement) == 1 && abs(sub1_edu_engagement[[1]] - 2.0) < 1e-9, "Subject 1 engagement Education marginal mean should equal mean(SF_Edu, LF_Edu).")
  assert_true(length(sub1_ent_engagement) == 1 && abs(sub1_ent_engagement[[1]] - 3.0) < 1e-9, "Subject 1 engagement Entertainment marginal mean should equal mean(SF_Ent, LF_Ent).")
  assert_true(length(sub2_edu_retention) == 1 && abs(sub2_edu_retention[[1]] - 0.10) < 1e-9, "Subject 2 retention Education marginal mean should equal mean(SF_Edu, LF_Edu).")
  assert_true(length(sub2_ent_retention) == 1 && abs(sub2_ent_retention[[1]] - 0.20) < 1e-9, "Subject 2 retention Entertainment marginal mean should equal mean(SF_Ent, LF_Ent).")

  retention_length <- audit %>% filter(.data$plot_type == "length_marginal", .data$domain == "retention")
  assert_true(nrow(retention_length) == 6, "Retention length marginal audit should contain 2 rows per retained subject.")
  assert_true(all(retention_length$n_source_conditions == 2L), "Length marginal rows should be computed from exactly two source conditions.")
  sub1_short_retention <- retention_length$score[
    retention_length$subject_id == 1L & retention_length$length_level == "Short"
  ]
  sub1_long_retention <- retention_length$score[
    retention_length$subject_id == 1L & retention_length$length_level == "Long"
  ]
  sub2_short_retention <- retention_length$score[
    retention_length$subject_id == 2L & retention_length$length_level == "Short"
  ]
  sub2_long_retention <- retention_length$score[
    retention_length$subject_id == 2L & retention_length$length_level == "Long"
  ]
  assert_true(length(sub1_short_retention) == 1 && abs(sub1_short_retention[[1]] - 0.15) < 1e-9, "Subject 1 retention Short marginal mean should equal mean(SF_Edu, SF_Ent).")
  assert_true(length(sub1_long_retention) == 1 && abs(sub1_long_retention[[1]] - 0.35) < 1e-9, "Subject 1 retention Long marginal mean should equal mean(LF_Edu, LF_Ent).")
  assert_true(length(sub2_short_retention) == 1 && abs(sub2_short_retention[[1]] - 0.05) < 1e-9, "Subject 2 retention Short marginal mean should equal mean(SF_Edu, SF_Ent).")
  assert_true(length(sub2_long_retention) == 1 && abs(sub2_long_retention[[1]] - 0.25) < 1e-9, "Subject 2 retention Long marginal mean should equal mean(LF_Edu, LF_Ent).")

  summary_df <- outputs$summary_df
  assert_true(nrow(summary_df) == 14, "Summary CSV should contain 8 raw-condition rows, 4 content-marginal rows, and 2 retention length-marginal rows.")
  assert_true(all(summary_df$n_subjects == 3L), "Each synthetic summary condition should have 3 retained subjects.")
}

# Compare the plotting script's preprocessing to the actual engagement and
# retention LMM preprocessing functions. This is the strongest guard against
# "confabulated" plot data: the raw plotted rows must exactly match the rows that
# the LMM scripts would model after exclusions, reshaping, and complete-case
# filtering.
test_lmm_preprocessing_agreement <- function(input_csv, exclusions_json) {
  plot_df <- plot_env$load_behavior_input(input_csv, exclusions_json)
  plot_long <- plot_env$reshape_behavior_scores(plot_df)
  plot_engagement <- plot_env$complete_case_domain_scores(plot_long, "engagement")
  plot_retention <- plot_env$complete_case_domain_scores(plot_long, "retention")

  lmm_engagement_input <- engagement_lmm_env$load_engagement_input(input_csv, exclusions_json)
  lmm_engagement_long <- engagement_lmm_env$reshape_to_long(lmm_engagement_input)
  lmm_engagement_cc <- engagement_lmm_env$complete_case_subjects(lmm_engagement_long)

  lmm_retention_input <- retention_lmm_env$load_retention_input(input_csv, exclusions_json)
  lmm_retention_long <- retention_lmm_env$reshape_to_long(lmm_retention_input)
  lmm_retention_cc <- retention_lmm_env$complete_case_subjects(lmm_retention_long)

  assert_true(
    identical(sort(unique(plot_engagement$subject_id)), sort(unique(lmm_engagement_cc$subject_id))),
    "Plot engagement complete-case subject IDs should match engagement LMM complete-case subject IDs."
  )
  assert_true(
    identical(sort(unique(plot_retention$subject_id)), sort(unique(lmm_retention_cc$subject_id))),
    "Plot retention complete-case subject IDs should match retention LMM complete-case subject IDs."
  )

  plot_eng_map <- plot_engagement %>%
    transmute(
      subject_id = .data$subject_id,
      condition = as.character(.data$condition),
      score_col = .data$score_col,
      score = .data$score
    ) %>%
    arrange(.data$subject_id, .data$condition)
  lmm_eng_map <- lmm_engagement_cc %>%
    transmute(
      subject_id = .data$subject_id,
      condition = as.character(.data$condition),
      score_col = .data$eng_col,
      score = .data$engagement
    ) %>%
    arrange(.data$subject_id, .data$condition)
  assert_true(
    identical(plot_eng_map, lmm_eng_map),
    "Plot engagement raw rows should match engagement LMM reshaped rows exactly."
  )

  plot_ret_map <- plot_retention %>%
    transmute(
      subject_id = .data$subject_id,
      condition = as.character(.data$condition),
      score_col = .data$score_col,
      score = .data$score
    ) %>%
    arrange(.data$subject_id, .data$condition)
  lmm_ret_map <- lmm_retention_cc %>%
    transmute(
      subject_id = .data$subject_id,
      condition = as.character(.data$condition),
      score_col = .data$diff_col,
      score = .data$retention_diff
    ) %>%
    arrange(.data$subject_id, .data$condition)
  assert_true(
    identical(plot_ret_map, lmm_ret_map),
    "Plot retention raw rows should match retention LMM reshaped rows exactly."
  )

  assert_true(
    all(lmm_engagement_cc$content_c[lmm_engagement_cc$condition %in% c("SF_Edu", "LF_Edu")] == 0.5) &&
      all(lmm_engagement_cc$content_c[lmm_engagement_cc$condition %in% c("SF_Ent", "LF_Ent")] == -0.5),
    "Engagement LMM content coding should be +0.5 for Education and -0.5 for Entertainment."
  )
  assert_true(
    all(lmm_engagement_cc$length_c[lmm_engagement_cc$condition %in% c("SF_Edu", "SF_Ent")] == -0.5) &&
      all(lmm_engagement_cc$length_c[lmm_engagement_cc$condition %in% c("LF_Edu", "LF_Ent")] == 0.5),
    "Engagement LMM length coding should be -0.5 for Short and +0.5 for Long."
  )
  assert_true(
    all(lmm_retention_cc$content_c[lmm_retention_cc$condition %in% c("SF_Edu", "LF_Edu")] == 0.5) &&
      all(lmm_retention_cc$content_c[lmm_retention_cc$condition %in% c("SF_Ent", "LF_Ent")] == -0.5),
    "Retention LMM content coding should be +0.5 for Education and -0.5 for Entertainment."
  )
  assert_true(
    all(lmm_retention_cc$length_c[lmm_retention_cc$condition %in% c("SF_Edu", "SF_Ent")] == -0.5) &&
      all(lmm_retention_cc$length_c[lmm_retention_cc$condition %in% c("LF_Edu", "LF_Ent")] == 0.5),
    "Retention LMM length coding should be -0.5 for Short and +0.5 for Long."
  )
}

# Build a controlled final-merged-data fixture where the true content effect is
# constant within subject and independent of length. The fixture still includes
# non-content variation so the mixed model is not a degenerate all-equal table.
# In this special case, the raw plotted Education-minus-Entertainment contrast
# should equal the LMM `content_c` estimate exactly.
write_controlled_lmm_alignment_input <- function(path) {
  subject_id <- sprintf("%04d", seq_len(10))
  age <- seq(18, by = 1, length.out = length(subject_id))
  subject_shift <- c(-0.35, 0.22, -0.08, 0.41, -0.27, 0.10, -0.44, 0.31, -0.16, 0.26)
  short_deviation <- c(-0.20, 0.05, 0.18, -0.12, 0.09, -0.05, 0.16, -0.18, 0.11, -0.04)
  long_deviation <- c(0.14, -0.10, 0.07, 0.20, -0.16, 0.04, -0.09, 0.12, -0.06, -0.15)

  df <- tibble::tibble(
    subject_id = subject_id,
    age = age,
    sf_education_engagement = 3.0 + subject_shift + short_deviation + 0.5,
    sf_entertainment_engagement = 3.0 + subject_shift + short_deviation - 0.5,
    lf_education_engagement = 3.0 + subject_shift + long_deviation + 0.5,
    lf_entertainment_engagement = 3.0 + subject_shift + long_deviation - 0.5,
    diff_short_form_education = 0.40 + subject_shift * 0.05 + short_deviation * 0.1 + 0.15,
    diff_short_form_entertainment = 0.40 + subject_shift * 0.05 + short_deviation * 0.1 - 0.15,
    diff_long_form_education = 0.40 + subject_shift * 0.05 + long_deviation * 0.1 + 0.15,
    diff_long_form_entertainment = 0.40 + subject_shift * 0.05 + long_deviation * 0.1 - 0.15,
    S01_D01_Cond01_HbO = seq_along(subject_id)
  )
  write_csv(df, path, na = "")
}

# Compute the mean subject-level Education minus Entertainment contrast from the
# plot audit rows. This intentionally uses the audit output rather than
# recomputing from source columns, so the test validates the values that would be
# inspected downstream.
content_contrast_from_plot_audit <- function(audit_df, domain_name) {
  wide <- audit_df %>%
    filter(.data$plot_type == "content_marginal", .data$domain == domain_name) %>%
    select("subject_id", "content_level", "score") %>%
    tidyr::pivot_wider(names_from = "content_level", values_from = "score")
  mean(wide$Education - wide$Entertainment)
}

length_contrast_from_plot_audit <- function(audit_df, domain_name) {
  wide <- audit_df %>%
    filter(.data$plot_type == "length_marginal", .data$domain == domain_name) %>%
    select("subject_id", "length_level", "score") %>%
    tidyr::pivot_wider(names_from = "length_level", values_from = "score")
  mean(wide$Long - wide$Short)
}

# Fit the same LMM functions used by the inferential scripts on the controlled
# fixture and compare their `content_c` coefficient to the plotted marginal
# contrast. This does not claim real-data plot points are age-adjusted; it proves
# alignment in a case where raw marginal and model-estimated content effects
# should be identical.
test_content_marginal_contrast_matches_lmm_when_expected <- function(input_csv, exclusions_json) {
  outputs <- plot_env$run_plotting(
    input_csv = input_csv,
    exclude_subjects_json = exclusions_json,
    out_dir = tempfile("behavior_lmm_alignment_out_"),
    width = 6,
    height = 4,
    dpi = 120
  )

  engagement_input <- engagement_lmm_env$load_engagement_input(input_csv, exclusions_json)
  engagement_cc <- engagement_lmm_env$complete_case_subjects(engagement_lmm_env$reshape_to_long(engagement_input))
  engagement_model <- engagement_lmm_env$fit_factorial_lmm(engagement_cc)
  engagement_content_est <- as.numeric(summary(engagement_model)$coefficients["content_c", "Estimate"])
  engagement_plot_contrast <- content_contrast_from_plot_audit(outputs$audit_df, "engagement")

  retention_input <- retention_lmm_env$load_retention_input(input_csv, exclusions_json)
  retention_cc <- retention_lmm_env$complete_case_subjects(retention_lmm_env$reshape_to_long(retention_input))
  retention_model <- retention_lmm_env$fit_factorial_lmm(retention_cc)
  retention_content_est <- as.numeric(summary(retention_model)$coefficients["content_c", "Estimate"])
  retention_plot_contrast <- content_contrast_from_plot_audit(outputs$audit_df, "retention")

  assert_true(
    abs(engagement_plot_contrast - engagement_content_est) < 1e-9,
    "In controlled no-age-effect engagement data, plotted Education-Entertainment contrast should equal LMM content estimate."
  )
  assert_true(
    abs(retention_plot_contrast - retention_content_est) < 1e-9,
    "In controlled no-age-effect retention data, plotted Education-Entertainment contrast should equal LMM content estimate."
  )
}

test_length_marginal_contrast_matches_lmm_when_expected <- function(input_csv, exclusions_json) {
  outputs <- plot_env$run_plotting(
    input_csv = input_csv,
    exclude_subjects_json = exclusions_json,
    out_dir = tempfile("behavior_length_lmm_alignment_out_"),
    width = 6,
    height = 4,
    dpi = 120,
    plot_types = c("length_marginal", "raw_condition")
  )

  retention_input <- retention_lmm_env$load_retention_input(input_csv, exclusions_json)
  retention_cc <- retention_lmm_env$complete_case_subjects(retention_lmm_env$reshape_to_long(retention_input))
  retention_model <- retention_lmm_env$fit_factorial_lmm(retention_cc)
  retention_length_est <- as.numeric(summary(retention_model)$coefficients["length_c", "Estimate"])
  retention_plot_contrast <- length_contrast_from_plot_audit(outputs$audit_df, "retention")

  assert_true(
    abs(retention_plot_contrast - retention_length_est) < 1e-9,
    "In controlled no-age-effect retention data, plotted Long-Short contrast should equal LMM length estimate."
  )
}

# Verify fail-fast numeric coercion. A non-numeric score token must not become an
# NA and then silently change complete-case membership or plotted values.
test_non_numeric_failure <- function(input_csv, exclusions_json) {
  bad_csv <- tempfile(fileext = ".csv")
  df <- read_csv(input_csv, show_col_types = FALSE, col_types = cols(.default = col_character()))
  df$sf_education_engagement[[1]] <- "bad"
  write_csv(df, bad_csv, na = "")

  ok <- FALSE
  msg <- NULL
  tryCatch(
    {
      plot_env$load_behavior_input(bad_csv, exclusions_json)
    },
    error = function(e) {
      ok <<- TRUE
      msg <<- conditionMessage(e)
    }
  )
  assert_true(ok, "Expected non-numeric engagement value to fail hard.")
  assert_true(grepl("non-numeric values", msg, fixed = TRUE), "Failure should mention non-numeric values.")
}

# Verify fail-fast schema validation. Dropping one required behavioral column
# would make the four-condition and content-marginal plots scientifically
# invalid, so the script must stop before plotting.
test_missing_column_failure <- function(input_csv, exclusions_json) {
  bad_csv <- tempfile(fileext = ".csv")
  df <- read_csv(input_csv, show_col_types = FALSE)
  df$diff_long_form_entertainment <- NULL
  write_csv(df, bad_csv, na = "")

  ok <- FALSE
  msg <- NULL
  tryCatch(
    {
      plot_env$load_behavior_input(bad_csv, exclusions_json)
    },
    error = function(e) {
      ok <<- TRUE
      msg <<- conditionMessage(e)
    }
  )
  assert_true(ok, "Expected missing required retention column to fail hard.")
  assert_true(grepl("Missing required columns", msg, fixed = TRUE), "Failure should mention missing required columns.")
}

# Exercise the CLI path, not only in-memory function calls. This catches argument
# parsing, default output writing, and script entry-point regressions.
test_cli_execution <- function(input_csv, exclusions_json, out_dir) {
  args <- c(
    shQuote(SCRIPT_PATH),
    "--input_csv", shQuote(input_csv),
    "--exclude_subjects_json", shQuote(exclusions_json),
    "--out_dir", shQuote(out_dir),
    "--width", "6",
    "--height", "4",
    "--dpi", "120"
  )
  output <- suppressWarnings(system2(RSCRIPT_BIN, args, stdout = TRUE, stderr = TRUE))
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L
  assert_true(status == 0L, paste("CLI execution failed:", paste(output, collapse = "\n")))
  assert_true(file.exists(file.path(out_dir, "plotted_behavior_scores.csv")), "CLI execution should write audit CSV.")
}

# Assemble synthetic fixtures and run the checks in an order that catches cheap
# source-level mistakes first, then semantic/data-integrity failures, then CLI
# regressions.
main <- function() {
  tmp <- file.path(tempdir(), "behavior_score_distribution_plot_validation")
  dir.create(tmp, showWarnings = FALSE, recursive = TRUE)

  input_csv <- file.path(tmp, "homer3_betas_plus_combined_sfv_data_inner_join.csv")
  exclusions_json <- file.path(tmp, "excluded_subjects.json")
  out_dir <- file.path(tmp, "out")
  cli_out_dir <- file.path(tmp, "out_cli")
  controlled_csv <- file.path(tmp, "controlled_lmm_alignment.csv")
  exclude_none_json <- file.path(tmp, "excluded_none.json")

  write_toy_input(input_csv)
  write_toy_exclusions(exclusions_json)
  write_controlled_lmm_alignment_input(controlled_csv)
  writeLines("[]", exclude_none_json)

  test_plot_geometry_source()
  test_run_plotting_outputs(input_csv, exclusions_json, out_dir)
  test_lmm_preprocessing_agreement(input_csv, exclusions_json)
  test_content_marginal_contrast_matches_lmm_when_expected(controlled_csv, exclude_none_json)
  test_length_marginal_contrast_matches_lmm_when_expected(controlled_csv, exclude_none_json)
  test_non_numeric_failure(input_csv, exclusions_json)
  test_missing_column_failure(input_csv, exclusions_json)
  test_cli_execution(input_csv, exclusions_json, cli_out_dir)

  writeLines("[PASS] validate_behavior_score_distribution_plot_r")
}

main()
