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

method_specific_out_csv <- function(out_csv, association_method) {
  dir_name <- dirname(out_csv)
  file_name <- basename(out_csv)
  stem <- sub("\\.csv$", "", file_name, ignore.case = TRUE)
  ext <- sub("^.*(\\.csv)$", "\\1", file_name, ignore.case = TRUE)
  if (!grepl("\\.csv$", file_name, ignore.case = TRUE)) {
    ext <- ".csv"
  }
  file.path(dir_name, paste0(stem, "_", association_method, ext))
}

run_script <- function(script, args) {
  output <- suppressWarnings(system2("Rscript", c(script, args), stdout = TRUE, stderr = TRUE))
  list(
    status = attr(output, "status") %||% 0,
    stdout = output
  )
}

normalize_channel_id <- function(x) {
  x_chr <- toupper(trimws(as.character(x)))
  parts <- regexec("^S(\\d+)_D(\\d+)$", x_chr)
  matches <- regmatches(x_chr, parts)
  out <- vapply(matches, function(m) {
    if (length(m) != 3) {
      stop("Failed to normalize ROI channel IDs for forced-signal validation.")
    }
    paste0("S", sprintf("%02d", as.integer(m[[2]])), "_D", sprintf("%02d", as.integer(m[[3]])))
  }, character(1))
  out
}

# Load the production functions without executing main(), so controlled fixtures
# can check ROI eligibility and numerical aggregation independently of real data.
load_analysis_functions <- function() {
  env <- new.env(parent = globalenv())
  expressions <- parse("analyze_correlational_relationships_roi_means.R")
  for (expression in expressions) {
    if (!identical(expression, quote(main()))) eval(expression, envir = env)
  }
  env
}

assert_error_contains <- function(expr, expected) {
  message <- tryCatch({ force(expr); NULL }, error = function(e) conditionMessage(e))
  assert_true(!is.null(message) && grepl(expected, message, fixed = TRUE),
              paste0("expected failure containing: ", expected))
}

test_roi_configuration_preserves_planned_comparisons <- function(env) {
  # A changed ROI definition must remove undefined targets without broadening
  # chromophore or ROI hypotheses. Check both removal and restoration: a planned
  # target is eligible exactly when defined; an unrelated new ROI is not eligible.
  roi_map <- tibble::tribble(
    ~roi, ~channel,
    "R_DLPFC", "S01_D01", "R_DLPFC", "S02_D01",
    "L_DLPFC", "S03_D01", "R_VMPFC", "S04_D01"
  )
  active <- env$active_roi_specs(roi_map)
  assert_true(setequal(paste(active$neural_name, active$chrom),
                      c("R_DLPFC HbR", "L_DLPFC HbO")),
              "configuration changed the planned ROI/chromophore comparisons")
  restored <- bind_rows(roi_map, tibble(roi = "M_DMPFC", channel = "S05_D01"))
  assert_true("M_DMPFC" %in% env$active_roi_specs(restored)$neural_name,
              "a defined planned target should remain eligible")
  assert_error_contains(env$active_roi_specs(filter(roi_map, roi == "R_VMPFC")),
                        "none of the planned ROI/chromophore targets")
  required <- env$collect_required_target_beta_columns(roi_map)
  expected <- c(paste0("S01_D01_Cond", sprintf("%02d", 1:4), "_HbR"),
                paste0("S02_D01_Cond", sprintf("%02d", 1:4), "_HbR"),
                paste0("S03_D01_Cond", sprintf("%02d", 1:4), "_HbO"))
  assert_true(setequal(required, expected),
              "required beta columns must follow eligible targets and configured channels exactly")
}

test_roi_channel_updates_and_missingness_preserve_means <- function(env) {
  # Hand-calculated channel values distinguish correct JSON membership from a
  # stale channel map or wrong chromophore. Pruned channels remain missing: one
  # available member supplies its value, all-missing cells remain NA, and a pooled
  # format needs both condition cells. These are the existing aggregation rules.
  roi_map <- tibble::tribble(~roi, ~channel,
    "R_DLPFC", "S01_D01", "R_DLPFC", "S02_D01", "L_DLPFC", "S03_D01")
  long <- tidyr::crossing(subject_id = 1L, cond = sprintf("%02d", 1:4),
                          channel = sprintf("S%02d_D01", 1:4), chrom = c("HbO", "HbR")) %>%
    left_join(env$CONDITION_MAP, by = "cond") %>%
    mutate(beta = case_when(
      channel == "S01_D01" & chrom == "HbR" ~ 2,
      channel == "S02_D01" & chrom == "HbR" ~ 6,
      channel == "S03_D01" & chrom == "HbO" ~ 10,
      TRUE ~ 100
    ), beta = ifelse(chrom == "HbR" &
      ((channel == "S02_D01" & cond == "02") |
       (channel %in% c("S01_D01", "S02_D01") & cond == "03")), NA_real_, beta))
  actual <- env$build_roi_condition_targets(long, roi_map)
  right <- actual %>% filter(neural_name == "R_DLPFC") %>% arrange(cond)
  assert_true(isTRUE(all.equal(right$neural_value, c(4, 2, NA_real_, 4))),
              "condition means or pruned-channel handling changed")
  assert_true(all(filter(actual, neural_name == "L_DLPFC")$neural_value == 10),
              "wrong chromophore or channel entered the left ROI")
  pooled <- env$build_neural_pooled_values(actual) %>% filter(neural_name == "R_DLPFC")
  assert_true(filter(pooled, format_pool == "short")$neural_value == 3,
              "short pooled mean should equal mean(4,2)")
  assert_true(is.na(filter(pooled, format_pool == "long")$neural_value),
              "one missing condition must invalidate its pooled format")
  updated <- roi_map %>% mutate(channel = ifelse(roi == "L_DLPFC", "S04_D01", channel))
  updated_values <- env$build_roi_condition_targets(long, updated)
  assert_true(all(filter(updated_values, neural_name == "L_DLPFC")$neural_value == 100),
              "ROI mean did not follow updated JSON channel membership")
  assert_error_contains(env$build_roi_condition_targets(filter(long, channel != "S03_D01"), roi_map),
                        "channels absent from merged input")
}

test_missing_retained_beta_column_still_fails <- function(env, tmp) {
  # Removing an ROI from the configuration is allowed; losing an input column
  # needed by a retained ROI is not. Otherwise its channel denominator could
  # change silently and alter neural correlations.
  roi_map <- env$load_roi_definition("data/config/roi_definition.json")
  required <- env$collect_required_target_beta_columns(roi_map)
  input <- read_csv("data/tabular/generated_data/homer3_betas_plus_combined_sfv_data_inner_join.csv", show_col_types = FALSE)
  path <- file.path(tmp, "missing_retained_beta.csv")
  write_csv(select(input, -all_of(required[[1]])), path)
  plan <- env$load_analysis_plan("data/config/correlational_analysis_plan_roi_means.json")
  assert_error_contains(env$load_merged_input(path, "data/config/excluded_subjects.json", plan, roi_map),
                        required[[1]])
}

test_no_eligible_targets_preserves_existing_outputs <- function(tmp) {
  # An unrelated-only configuration must fail explicitly, not generate an empty
  # successful analysis or erase earlier results before the mistake is reported.
  # The deliberately unrelated ROI makes eligibility empty without malformed JSON.
  roi_path <- file.path(tmp, "unrelated_only_roi.json")
  writeLines('{"UNPLANNED_ROI": ["S01_D01"]}', roi_path)
  output_root <- file.path(tmp, "no_eligible_targets")
  dir.create(output_root)
  sentinel <- file.path(output_root, "existing_results.txt")
  writeLines("preserve", sentinel)
  run <- run_script("analyze_correlational_relationships_roi_means.R", c(
    "--roi_json", roi_path,
    "--out_csv", file.path(output_root, "results.csv"),
    "--out_fig_dir", file.path(output_root, "figures")
  ))
  assert_true(run$status != 0, "configuration with no eligible targets must fail")
  assert_true(any(grepl("none of the planned ROI/chromophore targets", run$stdout, fixed = TRUE)),
              "failure must identify the empty eligible target set")
  assert_true(file.exists(sentinel) && identical(readLines(sentinel), "preserve"),
              "invalid target configuration cleared existing outputs")
}

build_forced_signal_input <- function(source_csv, roi_json, out_csv) {
  df <- read_csv(source_csv, show_col_types = FALSE)
  roi_obj <- jsonlite::fromJSON(roi_json, simplifyVector = FALSE)
  l_dlpfc_channels <- normalize_channel_id(unlist(roi_obj$L_DLPFC, use.names = FALSE))
  short_cols <- unlist(lapply(l_dlpfc_channels, function(channel) {
    paste0(channel, c("_Cond01_HbO", "_Cond02_HbO"))
  }), use.names = FALSE)
  missing_cols <- setdiff(short_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(paste0("Forced-signal validation is missing expected beta columns: ", paste(missing_cols, collapse = ", ")))
  }
  cond01_cols <- paste0(l_dlpfc_channels, "_Cond01_HbO")
  cond02_cols <- paste0(l_dlpfc_channels, "_Cond02_HbO")
  cond01_mat <- as.matrix(df[, cond01_cols, drop = FALSE])
  cond02_mat <- as.matrix(df[, cond02_cols, drop = FALSE])
  cond01_mat[cond01_mat == 0] <- NA_real_
  cond02_mat[cond02_mat == 0] <- NA_real_

  cond01_roi <- rowMeans(cond01_mat, na.rm = TRUE)
  cond02_roi <- rowMeans(cond02_mat, na.rm = TRUE)
  cond01_roi[!is.finite(cond01_roi)] <- NA_real_
  cond02_roi[!is.finite(cond02_roi)] <- NA_real_
  forced_signal <- ifelse(
    is.finite(cond01_roi) & is.finite(cond02_roi),
    rowMeans(cbind(cond01_roi, cond02_roi)),
    NA_real_
  )
  finite_n <- sum(is.finite(forced_signal))
  if (finite_n < 10) {
    stop("Forced-signal validation did not retain enough finite pooled ROI short means.")
  }
  df$sf_education_engagement <- forced_signal
  df$sf_entertainment_engagement <- forced_signal
  write_csv(df, out_csv, na = "NA")
}

main <- function() {
  # End-to-end real and forced-positive runs below preserve Pearson inference,
  # participant exclusions, figure policy, and BH families while the target set
  # follows the configuration. Manual BH checks guard the changed family size.
  tmp <- file.path(tempdir(), "correlational_relationships_roi_means_validation")
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  env <- load_analysis_functions()
  test_roi_configuration_preserves_planned_comparisons(env)
  test_roi_channel_updates_and_missingness_preserve_means(env)
  test_missing_retained_beta_column_still_fails(env, tmp)
  test_no_eligible_targets_preserves_existing_outputs(tmp)
  configured_names <- names(jsonlite::fromJSON("data/config/roi_definition.json"))
  planned <- c(R_DLPFC = "HbR", L_DLPFC = "HbO", M_DMPFC = "HbO", L_DMPFC = "HbO")
  eligible <- planned[names(planned) %in% configured_names]
  expected_targets <- length(eligible)
  assert_true(expected_targets > 0, "real-data validation needs eligible planned targets")

  out_dir <- file.path(tmp, "observed")
  out_csv <- file.path(out_dir, "pairwise_correlations_r.csv")
  out_fig_dir <- file.path(out_dir, "figures")
  dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
  dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)
  writeLines("stale", file.path(dirname(out_csv), "stale.csv"))
  writeLines("stale", file.path(out_fig_dir, "stale.txt"))

  run <- run_script(
    "analyze_correlational_relationships_roi_means.R",
    c(
      "--input_csv", "data/tabular/generated_data/homer3_betas_plus_combined_sfv_data_inner_join.csv",
      "--roi_json", "data/config/roi_definition.json",
      "--analysis_plan_json", "data/config/correlational_analysis_plan_roi_means.json",
      "--exclude_subjects_json", "data/config/excluded_subjects.json",
      "--out_csv", out_csv,
      "--out_fig_dir", out_fig_dir
    )
  )

  assert_true(run$status == 0, paste("standalone ROI-mean script did not exit cleanly:", paste(run$stdout, collapse = "\n")))
  assert_true(any(grepl("[scope] active ROI/chromophore targets:", run$stdout, fixed = TRUE)),
              "script must report the active ROI/chromophore scope")
  for (retired in setdiff(names(planned), configured_names)) {
    assert_true(any(grepl(paste0(retired, " (", planned[[retired]], ")"),
                         run$stdout[grepl("skipped planned targets", run$stdout)], fixed = TRUE)),
                paste0("undefined target was not disclosed: ", retired))
  }
  assert_true(any(grepl("\\[out\\] results CSV:", run$stdout)), "script did not report its output path")
  assert_true(!file.exists(file.path(dirname(out_csv), "stale.csv")), "stale ROI-mean output CSV was not removed before rerun")
  assert_true(!file.exists(file.path(out_fig_dir, "stale.txt")), "stale ROI-mean figure artifact was not removed before rerun")

  pearson_csv <- method_specific_out_csv(out_csv, "pearson")
  spearman_csv <- method_specific_out_csv(out_csv, "spearman")
  assert_true(file.exists(out_csv), "missing combined ROI-mean CSV")
  assert_true(file.exists(pearson_csv), "missing ROI-mean Pearson CSV")
  assert_true(!file.exists(spearman_csv), "ROI-mean Spearman CSV should not be written")

  results <- read_csv(out_csv, show_col_types = FALSE)
  pearson <- read_csv(pearson_csv, show_col_types = FALSE)

  assert_true(nrow(results) == 4 * expected_targets, "unexpected number of standalone ROI-mean rows")
  assert_true(nrow(pearson) == 4 * expected_targets, "unexpected number of standalone Pearson rows")
  assert_true(setequal(paste(results$neural_name, results$chrom), paste(names(eligible), eligible)),
              "output contains an undefined or unplanned ROI/chromophore comparison")
  assert_true(all(results$behavior_run_type == "pooled_format"), "standalone output should contain only pooled behavior rows")
  assert_true(all(results$neural_level == "roi"), "standalone output should contain only ROI rows")
  assert_true(all(results$analysis_tier == "primary"), "standalone output should contain only primary-tier rows")
  assert_true(setequal(unique(results$behavior_run), c("engagement", "retention")), "unexpected behavior runs in standalone output")
  assert_true(setequal(unique(results$format_pool), c("long", "short")), "unexpected format pools in standalone output")
  assert_true(identical(unique(results$association_method), "pearson"), "standalone ROI-mean output should be Pearson-only")

  family_counts <- results %>% count(.data$family_id, name = "n_family")
  assert_true(nrow(family_counts) == 4, "unexpected number of multiple-testing families")
  assert_true(all(family_counts$n_family == expected_targets),
              "each standalone family should contain exactly the eligible planned ROI targets")

  joined_counts <- results %>% select("family_id", "family_n_tested") %>% distinct()
  assert_true(
    all(joined_counts$family_n_tested == family_counts$n_family[match(joined_counts$family_id, family_counts$family_id)]),
    "family_n_tested does not match tested rows per family"
  )

  family_summaries <- results %>%
    group_by(.data$family_id) %>%
    summarize(max_abs_diff = max(abs(bh_qvalues_manual(.data$p_unc) - .data$p_fdr), na.rm = TRUE), .groups = "drop")
  assert_true(all(family_summaries$max_abs_diff < 1e-12), "standalone BH-FDR values do not match manual implementation")

  assert_true(all(pearson$association_method == "pearson"), "Pearson CSV should contain only Pearson rows")

  sig_rows <- results %>% filter(.data$analysis_status == "tested", is.finite(.data$p_unc), .data$p_unc < 0.05)
  if (nrow(sig_rows) > 0) {
    assert_true(all(!is.na(sig_rows$plot_file)), "significant rows should record figure paths")
    assert_true(all(file.exists(sig_rows$plot_file)), "significant figure files were not created")
  }

  nonsig_rows <- results %>% filter(.data$analysis_status == "tested", is.finite(.data$p_unc), .data$p_unc >= 0.05)
  if (nrow(nonsig_rows) > 0) {
    assert_true(all(is.na(nonsig_rows$plot_file)), "non-significant rows should not receive figures under the default policy")
  }

  forced_input_csv <- file.path(tmp, "forced_signal_input.csv")
  build_forced_signal_input(
    source_csv = "data/tabular/generated_data/homer3_betas_plus_combined_sfv_data_inner_join.csv",
    roi_json = "data/config/roi_definition.json",
    out_csv = forced_input_csv
  )
  forced_dir <- file.path(tmp, "forced")
  forced_out_csv <- file.path(forced_dir, "pairwise_correlations_r.csv")
  forced_out_fig_dir <- file.path(forced_dir, "figures")
  forced_run <- run_script(
    "analyze_correlational_relationships_roi_means.R",
    c(
      "--input_csv", forced_input_csv,
      "--roi_json", "data/config/roi_definition.json",
      "--analysis_plan_json", "data/config/correlational_analysis_plan_roi_means.json",
      "--exclude_subjects_json", "data/config/excluded_subjects.json",
      "--out_csv", forced_out_csv,
      "--out_fig_dir", forced_out_fig_dir
    )
  )
  assert_true(forced_run$status == 0, "forced-signal ROI-mean run did not exit cleanly")
  forced <- read_csv(forced_out_csv, show_col_types = FALSE)
  forced_sig <- forced %>%
    filter(
      .data$behavior_run == "engagement",
      .data$format_pool == "short",
      .data$neural_name == "L_DLPFC",
      .data$chrom == "HbO"
    )
  assert_true(nrow(forced_sig) == 1, "forced-signal validation did not find the expected Pearson target row")
  assert_true(forced_sig$association_method[[1]] == "pearson", "forced-signal validation should be Pearson-only")
  assert_true(all(is.finite(forced_sig$association_estimate) & forced_sig$association_estimate > 0.99), "forced-signal correlations were weaker than expected")
  assert_true(all(is.finite(forced_sig$p_unc) & forced_sig$p_unc < 0.05), "forced-signal rows were not uncorrected-significant")
  assert_true(all(!is.na(forced_sig$plot_file)), "forced-signal significant rows did not record plot paths")
  assert_true(all(file.exists(forced_sig$plot_file)), "forced-signal figure files were not created")

  writeLines(paste0("[OK] ROI-mean correlation validation passed. Outputs in: ", tmp))
}

main()
