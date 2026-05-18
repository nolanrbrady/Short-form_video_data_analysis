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

bh_qvalues_manual <- function(p_values) {
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

write_roi_json <- function(path) {
  roi_obj <- list(R_DLPFC = c("S10_D01", "S10_D02"))
  jsonlite::write_json(roi_obj, path = path, auto_unbox = TRUE, pretty = TRUE)
}

write_target_csvs <- function(channel_path, roi_path) {
  channel_df <- tibble::tibble(
    channel = c("S01_D01", "S02_D01", "S09_D09"),
    chrom = c("HbO", "HbR", "HbO"),
    effect = c("format", "interaction", "content"),
    estimate = c(1.2e-5, -8.0e-6, 1.0e-6),
    p_fdr = c(0.01, 0.02, 0.30)
  )
  roi_df <- tibble::tibble(
    roi = c("R_DLPFC", "R_DLPFC"),
    chrom = c("HbO", "HbR"),
    effect = c("content", "interaction"),
    estimate = c(6.5e-6, 1.0e-6),
    p_fdr = c(0.03, 0.40)
  )
  write_csv(channel_df, channel_path)
  write_csv(roi_df, roi_path)
}

build_input <- function() {
  n <- 8
  subject_id <- sprintf("sub_%04d", seq_len(n))
  base <- c(40, 42, 44, 46, 48, 50, 52, 54)

  short_signal <- c(1, 2, 3, 4, 5, 6, 7, 8)
  long_signal <- c(8, 7, 6, 5, 4, 3, 2, 1)
  edu_signal <- c(4, 1, 6, 2, 7, 3, 8, 5)

  df <- tibble::tibble(
    subject_id = subject_id,
    sf_education_engagement = base + short_signal,
    sf_entertainment_engagement = base + short_signal,
    lf_entertainment_engagement = base + long_signal,
    lf_education_engagement = base + long_signal,
    diff_short_form_education = base + 30 + short_signal,
    diff_short_form_entertainment = base + 30 + short_signal,
    diff_long_form_entertainment = base + 30 + long_signal,
    diff_long_form_education = base + 30 + long_signal
  )

  add_channel <- function(channel, chrom, sfe, sft, lft, lfe) {
    df[[paste0(channel, "_Cond01_", chrom)]] <<- sfe
    df[[paste0(channel, "_Cond02_", chrom)]] <<- sft
    df[[paste0(channel, "_Cond03_", chrom)]] <<- lft
    df[[paste0(channel, "_Cond04_", chrom)]] <<- lfe
  }

  add_channel("S01_D01", "HbO", base + short_signal, base + short_signal, base + long_signal, base + long_signal)
  add_channel("S02_D01", "HbR", 200 - (base + short_signal), 202 - (base + short_signal), 300 - (base + long_signal), 302 - (base + long_signal))
  add_channel("S10_D01", "HbO", base + edu_signal, base + 100, base + 100, base + edu_signal)
  add_channel("S10_D02", "HbO", base + edu_signal, base + 100, base + 100, base + edu_signal)
  add_channel("S09_D09", "HbO", base + 1, base + 2, base + 3, base + 4)

  df[8, "S01_D01_Cond03_HbO"] <- 0
  df[4, "S10_D02_Cond02_HbO"] <- 0
  df[5, "S10_D01_Cond03_HbO"] <- NA_real_
  df
}

run_script <- function(args) {
  output <- suppressWarnings(system2("Rscript", c("analyze_pooled_mean_correlations.R", args), stdout = TRUE, stderr = TRUE))
  list(status = attr(output, "status") %||% 0, stdout = output)
}

main <- function() {
  tmp <- file.path(tempdir(), "pooled_mean_correlation_validation")
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)

  input_csv <- file.path(tmp, "merged.csv")
  roi_json <- file.path(tmp, "roi.json")
  channel_csv <- file.path(tmp, "channel_tidy.csv")
  roi_csv <- file.path(tmp, "roi_tidy.csv")
  exclude_json <- file.path(tmp, "excluded.json")
  out_dir <- file.path(tmp, "out")

  write_csv(build_input(), input_csv, na = "NA")
  write_roi_json(roi_json)
  write_target_csvs(channel_csv, roi_csv)
  writeLines("[]", exclude_json)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  writeLines("stale", file.path(out_dir, "stale.txt"))

  run <- run_script(c(
    "--input_csv", input_csv,
    "--roi_json", roi_json,
    "--channel_results_csv", channel_csv,
    "--roi_results_csv", roi_csv,
    "--exclude_subjects_json", exclude_json,
    "--out_dir", out_dir,
    "--figure_policy", "significant_only",
    "--min_subjects", "3"
  ))
  assert_true(run$status == 0, "pooled-mean correlation script did not exit cleanly")
  assert_true(!file.exists(file.path(out_dir, "stale.txt")), "output directory was not rebuilt at run start")

  targets <- read_csv(file.path(out_dir, "selected_pooled_mean_targets_r.csv"), show_col_types = FALSE)
  pairs <- read_csv(file.path(out_dir, "subject_level_pooled_mean_pairs_r.csv"), show_col_types = FALSE)
  results <- read_csv(file.path(out_dir, "pooled_mean_correlations_r.csv"), show_col_types = FALSE)

  assert_true(nrow(targets) == 2, "expected exactly two selected pooled-mean main-effect targets")
  assert_true(setequal(targets$selected_effect, c("format", "content")), "pooled targets should retain only significant main effects")
  assert_true(identical(unique(results$association_method), "pearson"), "pooled-mean script should be Pearson-only")
  assert_true(setequal(unique(results$pool_name), c("short", "long", "education", "entertainment")), "unexpected pool names")

  gated_out <- results %>%
    filter((selected_effect == "format" & pool_name %in% c("education", "entertainment")) |
             (selected_effect == "content" & pool_name %in% c("short", "long")))
  assert_true(nrow(gated_out) == 0, "pooled rows should be gated to the matching main-effect pool family")

  short_pair <- pairs %>%
    filter(analysis_level == "channel", unit_id == "S01_D01", chrom == "HbO", behavior_domain == "engagement", pool_name == "short", subject_id == 1)
  assert_true(nrow(short_pair) == 1, "missing expected short pooled row")
  assert_true(abs(short_pair$neural_value[[1]] - 41) < 1e-12, "short pooled neural mean reconstructed incorrectly")
  assert_true(abs(short_pair$behavior_value[[1]] - 41) < 1e-12, "short pooled behavior mean reconstructed incorrectly")

  dropped_long <- pairs %>%
    filter(analysis_level == "channel", unit_id == "S01_D01", chrom == "HbO", behavior_domain == "engagement", pool_name == "long", subject_id == 8)
  assert_true(nrow(dropped_long) == 1, "missing expected dropped long pooled row")
  assert_true(is.na(dropped_long$neural_value[[1]]), "zero placeholder should invalidate the long pooled neural mean")

  roi_edu <- pairs %>%
    filter(analysis_level == "roi", unit_id == "R_DLPFC", chrom == "HbO", behavior_domain == "engagement", pool_name == "education", subject_id == 4)
  assert_true(nrow(roi_edu) == 1, "missing expected ROI education row")
  assert_true(abs(roi_edu$neural_value[[1]] - 48) < 1e-12, "ROI education mean should use available member channels when one entertainment cell is pruned")
  assert_true(roi_edu$roi_member_count[[1]] == 2, "ROI member count metadata should reflect both ROI channels")

  known_positive <- results %>%
    filter(analysis_level == "channel", unit_id == "S01_D01", chrom == "HbO", behavior_domain == "engagement", pool_name == "short")
  assert_true(nrow(known_positive) == 1, "missing known positive pooled Pearson row")
  assert_true(known_positive$analysis_status[[1]] == "tested", "known positive pooled row should be tested")
  assert_true(known_positive$n_complete[[1]] == 8, "known positive pooled short row should retain all subjects")
  assert_true(abs(known_positive$association_estimate[[1]] - 1.0) < 1e-12, "expected perfect positive pooled Pearson correlation")

  known_long <- results %>%
    filter(analysis_level == "channel", unit_id == "S02_D01", chrom == "HbR", behavior_domain == "engagement", pool_name == "long")
  assert_true(nrow(known_long) == 0, "interaction-gated channel target should not produce pooled long rows")

  known_long_positive <- results %>%
    filter(analysis_level == "channel", unit_id == "S01_D01", chrom == "HbO", behavior_domain == "engagement", pool_name == "long")
  assert_true(nrow(known_long_positive) == 1, "missing known gated long pooled Pearson row")
  assert_true(known_long_positive$n_complete[[1]] == 7, "known gated long row should drop exactly one subject due to the zero placeholder")
  assert_true(abs(known_long_positive$association_estimate[[1]] - 1.0) < 1e-12, "expected perfect positive gated long Pearson correlation")

  low_signal <- results %>%
    filter(analysis_level == "channel", unit_id == "S09_D09", chrom == "HbO", behavior_domain == "engagement", pool_name == "education")
  assert_true(nrow(low_signal) == 0, "non-significant pooled targets should not be retained")

  interaction_target <- results %>%
    filter(analysis_level == "channel", unit_id == "S02_D01", chrom == "HbR")
  assert_true(nrow(interaction_target) == 0, "interaction-only targets should be excluded from pooled main-effect follow-up")

  family_summary <- results %>%
    filter(.data$analysis_status == "tested") %>%
    group_by(.data$family_id) %>%
    summarize(max_abs = max(abs(bh_qvalues_manual(.data$p_unc) - .data$p_fdr), na.rm = TRUE), .groups = "drop")
  assert_true(all(family_summary$max_abs < 1e-12), "manual BH check failed for pooled means")

  plotted <- results %>%
    filter(.data$analysis_status == "tested", is.finite(.data$p_unc), .data$p_unc < 0.05)
  if (nrow(plotted) > 0) {
    assert_true(all(!is.na(plotted$plot_file)), "significant pooled-mean rows should record figure paths")
    assert_true(all(file.exists(plotted$plot_file)), "significant pooled-mean figure files were not created")
  }

  writeLines(paste0("[OK] Pooled-mean correlation validation passed. Outputs in: ", tmp))
}

main()
