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
  short_mat <- as.matrix(df[, short_cols, drop = FALSE])
  short_mat[short_mat == 0] <- NA_real_
  forced_signal <- rowMeans(short_mat, na.rm = TRUE)
  finite_n <- sum(is.finite(forced_signal))
  if (finite_n < 10) {
    stop("Forced-signal validation did not retain enough finite pooled ROI short means.")
  }
  df$sf_education_engagement <- forced_signal
  df$sf_entertainment_engagement <- forced_signal
  write_csv(df, out_csv, na = "NA")
}

main <- function() {
  tmp <- file.path(tempdir(), "correlational_relationships_roi_means_validation")
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)

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
      "--input_csv", "data/tabular/homer3_betas_plus_combined_sfv_data_inner_join.csv",
      "--roi_json", "data/config/roi_definition.json",
      "--analysis_plan_json", "data/config/correlational_analysis_plan_roi_means.json",
      "--exclude_subjects_json", "data/config/excluded_subjects.json",
      "--out_csv", out_csv,
      "--out_fig_dir", out_fig_dir
    )
  )

  assert_true(run$status == 0, "standalone ROI-mean script did not exit cleanly")
  assert_true(any(grepl("\\[out\\] results CSV:", run$stdout)), "script did not report its output path")
  assert_true(!file.exists(file.path(dirname(out_csv), "stale.csv")), "stale ROI-mean output CSV was not removed before rerun")
  assert_true(!file.exists(file.path(out_fig_dir, "stale.txt")), "stale ROI-mean figure artifact was not removed before rerun")

  pearson_csv <- method_specific_out_csv(out_csv, "pearson")
  spearman_csv <- method_specific_out_csv(out_csv, "spearman")
  assert_true(file.exists(out_csv), "missing combined ROI-mean CSV")
  assert_true(file.exists(pearson_csv), "missing ROI-mean Pearson CSV")
  assert_true(file.exists(spearman_csv), "missing ROI-mean Spearman CSV")

  results <- read_csv(out_csv, show_col_types = FALSE)
  pearson <- read_csv(pearson_csv, show_col_types = FALSE)
  spearman <- read_csv(spearman_csv, show_col_types = FALSE)

  assert_true(nrow(results) == 32, "unexpected number of standalone ROI-mean rows")
  assert_true(nrow(pearson) == 16, "unexpected number of standalone Pearson rows")
  assert_true(nrow(spearman) == 16, "unexpected number of standalone Spearman rows")
  assert_true(all(results$behavior_run_type == "pooled_format"), "standalone output should contain only pooled behavior rows")
  assert_true(all(results$neural_level == "roi"), "standalone output should contain only ROI rows")
  assert_true(all(results$analysis_tier == "primary"), "standalone output should contain only primary-tier rows")
  assert_true(setequal(unique(results$behavior_run), c("engagement", "retention")), "unexpected behavior runs in standalone output")
  assert_true(setequal(unique(results$format_pool), c("long", "short")), "unexpected format pools in standalone output")
  assert_true(setequal(unique(results$association_method), c("pearson", "spearman")), "unexpected association methods in standalone output")

  family_counts <- results %>% count(.data$family_id, name = "n_family")
  assert_true(nrow(family_counts) == 8, "unexpected number of multiple-testing families")
  assert_true(all(family_counts$n_family == 4), "each standalone family should contain four ROI targets")

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
  assert_true(all(spearman$association_method == "spearman"), "Spearman CSV should contain only Spearman rows")

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
    source_csv = "data/tabular/homer3_betas_plus_combined_sfv_data_inner_join.csv",
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
  assert_true(nrow(forced_sig) == 2, "forced-signal validation did not find the expected Pearson and Spearman target rows")
  assert_true(all(is.finite(forced_sig$association_estimate) & forced_sig$association_estimate > 0.99), "forced-signal correlations were weaker than expected")
  assert_true(all(is.finite(forced_sig$p_unc) & forced_sig$p_unc < 0.05), "forced-signal rows were not uncorrected-significant")
  assert_true(all(!is.na(forced_sig$plot_file)), "forced-signal significant rows did not record plot paths")
  assert_true(all(file.exists(forced_sig$plot_file)), "forced-signal figure files were not created")

  writeLines(paste0("[OK] ROI-mean correlation validation passed. Outputs in: ", tmp))
}

main()
