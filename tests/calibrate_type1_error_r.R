#!/usr/bin/env Rscript

# Monte Carlo type-I error calibration for all R inferential pipelines.
#
# Purpose
# - Run repeated null-effect synthetic datasets through:
#   - analyze_format_content_lmm_channelwise.R
#   - analyze_format_content_lmm_roi.R
#   - analyze_retention_format_content_lmm.R
#   - analyze_engagement_format_content_lmm.R
# - Estimate empirical false-positive rates from adjusted p-values.
# - Fail when observed rates exceed a configured upper bound.
#
# Methodological rationale
# - Simulation-based calibration of statistical procedures is a standard way to
#   evaluate operating characteristics (e.g., type-I error) under controlled
#   data-generating mechanisms (Burton et al., 2006; see CITATIONS.md).
# - Multiple-testing procedures are those used in the pipelines themselves:
#   BH-FDR (Benjamini & Hochberg, 1995) and Holm correction (Holm, 1979).

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

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  defaults <- list(
    n_reps = 60L,
    alpha = 0.05,
    type1_upper_bound = 0.12,
    n_subjects = 32L,
    noise_sd_channel = 0.20,
    noise_sd_retention = 0.10,
    noise_sd_engagement = 0.12,
    min_subjects = 6L,
    seed = 20260214L,
    out_csv = "data/results/type1_error_calibration_summary_r.csv"
  )
  if (length(args) == 0) return(defaults)
  if (length(args) %% 2 != 0) stop("Expected --key value argument pairs.")
  parsed <- defaults
  for (i in seq(1, length(args), by = 2)) {
    key <- args[[i]]
    val <- args[[i + 1]]
    if (!startsWith(key, "--")) stop(paste0("Invalid argument: ", key))
    parsed[[substring(key, 3)]] <- val
  }
  parsed$n_reps <- as.integer(parsed$n_reps)
  parsed$alpha <- as.numeric(parsed$alpha)
  parsed$type1_upper_bound <- as.numeric(parsed$type1_upper_bound)
  parsed$n_subjects <- as.integer(parsed$n_subjects)
  parsed$noise_sd_channel <- as.numeric(parsed$noise_sd_channel)
  parsed$noise_sd_retention <- as.numeric(parsed$noise_sd_retention)
  parsed$noise_sd_engagement <- as.numeric(parsed$noise_sd_engagement)
  parsed$min_subjects <- as.integer(parsed$min_subjects)
  parsed$seed <- as.integer(parsed$seed)
  parsed
}

make_subject_ids <- function(n_subjects) {
  ids_int <- seq_len(n_subjects)
  list(
    homer_subject = sprintf("sub_%04d", ids_int),
    combined_subject = sprintf("%04d", ids_int)
  )
}

make_age_years <- function(n_subjects) {
  seq(18, by = 1, length.out = n_subjects)
}

cond_grid <- function() {
  tibble::tibble(
    cond = c("01", "02", "03", "04"),
    condition = c("SF_Edu", "SF_Ent", "LF_Ent", "LF_Edu"),
    format_c = c(-0.5, -0.5, +0.5, +0.5),
    content_c = c(+0.5, -0.5, -0.5, +0.5)
  )
}

generate_homer3_wide_null <- function(n_subjects, channels, chroms, noise_sd, seed, age_years, age_beta = 0.10) {
  set.seed(seed)
  ids <- make_subject_ids(n_subjects)
  grid <- cond_grid()

  # Shared subject random intercepts to mimic repeated-measures structure.
  subj_re <- rnorm(n_subjects, mean = 0, sd = 0.6)
  age_centered <- age_years - mean(age_years)

  out <- tibble::tibble(Subject = ids$homer_subject)
  for (ch in channels) {
    for (chrom in chroms) {
      b0 <- rnorm(1, mean = 0, sd = 0.3)
      for (i in seq_len(nrow(grid))) {
        col_name <- paste0(ch, "_Cond", grid$cond[[i]], "_", chrom)
        out[[col_name]] <- b0 + age_beta * age_centered + subj_re + rnorm(n_subjects, mean = 0, sd = noise_sd)
      }
    }
  }
  out
}

generate_combined <- function(n_subjects, age_years = make_age_years(n_subjects)) {
  ids <- make_subject_ids(n_subjects)
  tibble::tibble(
    subject_id = ids$combined_subject,
    age = age_years,
    pd_status = rep(0, n_subjects)
  )
}

build_merged_input <- function(homer, combined) {
  homer_norm <- homer %>%
    mutate(subject_id_norm = as.integer(str_extract(.data$Subject, "\\d+"))) %>%
    select(-Subject)
  combined_norm <- combined %>%
    mutate(subject_id_norm = as.integer(str_extract(.data$subject_id, "\\d+")))

  merged <- combined_norm %>%
    inner_join(homer_norm, by = "subject_id_norm", relationship = "one-to-one") %>%
    select(-subject_id_norm)
  merged
}

generate_retention_null <- function(n_subjects, noise_sd, seed, age_years = make_age_years(n_subjects), age_beta = 0.10) {
  set.seed(seed)
  sid <- sprintf("%04d", seq_len(n_subjects))
  subj_re <- rnorm(n_subjects, mean = 0, sd = 0.35)
  b0 <- 0.30
  age_centered <- age_years - mean(age_years)

  make_cell <- function() {
    b0 + age_beta * age_centered + subj_re + rnorm(n_subjects, mean = 0, sd = noise_sd)
  }

  tibble::tibble(
    subject_id = sid,
    age = age_years,
    diff_short_form_education = make_cell(),
    diff_short_form_entertainment = make_cell(),
    diff_long_form_education = make_cell(),
    diff_long_form_entertainment = make_cell()
  )
}

generate_engagement_null <- function(n_subjects, noise_sd, seed, age_years = make_age_years(n_subjects), age_beta = 0.10) {
  set.seed(seed)
  sid <- sprintf("%04d", seq_len(n_subjects))
  subj_re <- rnorm(n_subjects, mean = 0, sd = 0.50)
  b0 <- 3.0
  age_centered <- age_years - mean(age_years)

  make_cell <- function() {
    b0 + age_beta * age_centered + subj_re + rnorm(n_subjects, mean = 0, sd = noise_sd)
  }

  tibble::tibble(
    subject_id = sid,
    age = age_years,
    sf_education_engagement = make_cell(),
    sf_entertainment_engagement = make_cell(),
    lf_education_engagement = make_cell(),
    lf_entertainment_engagement = make_cell()
  )
}

run_analysis <- function(script, args, env = NULL) {
  output <- suppressWarnings(system2("Rscript", c(script, args), stdout = TRUE, stderr = TRUE, env = env))
  list(status = attr(output, "status") %||% 0, output = output)
}

assert_output_file <- function(path, run_result, label, rep_idx) {
  if (!file.exists(path)) {
    tail_msg <- paste(utils::tail(run_result$output, 20), collapse = " | ")
    fail(
      paste0(
        label, " did not produce expected output file at replicate ", rep_idx, ": ", path,
        ". status=", run_result$status, ". tail(output)=", tail_msg
      )
    )
  }
}

append_effect_rows <- function(acc, pipeline_name, rep_idx, effect_name, p_fdr_values, alpha) {
  acc[[length(acc) + 1]] <- tibble::tibble(
    pipeline = pipeline_name,
    replicate = rep_idx,
    effect = effect_name,
    n_tests = length(p_fdr_values),
    n_reject = sum(p_fdr_values < alpha, na.rm = TRUE)
  )
  acc
}

main <- function() {
  args <- parse_args()
  if (args$n_reps < 10) stop("n_reps must be >= 10 for a useful calibration run.")
  if (args$n_subjects < args$min_subjects) stop("n_subjects must be >= min_subjects.")

  tmp <- file.path(tempdir(), "type1_error_calibration_r")
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)

  channels <- c("S01_D01", "S01_D02", "S02_D01", "S02_D02")
  chroms <- c("HbO", "HbR")
  roi_map <- list(
    VMPFC = c("S01_D01", "S01_D02"),
    DLPFC = c("S02_D01", "S02_D02")
  )
  roi_json <- file.path(tmp, "roi_definition.json")
  write_json(roi_map, roi_json, auto_unbox = TRUE, pretty = TRUE)

  exclude_none <- file.path(tmp, "excluded_none.json")
  writeLines("[]", exclude_none)

  effect_rows <- list()

  for (rep_idx in seq_len(args$n_reps)) {
    if (rep_idx %% 10 == 0) {
      cat("[progress] replicate ", rep_idx, "/", args$n_reps, "\n", sep = "")
    }

    seed_rep <- args$seed + rep_idx
    age_years <- make_age_years(args$n_subjects)

    homer <- generate_homer3_wide_null(
      n_subjects = args$n_subjects,
      channels = channels,
      chroms = chroms,
      noise_sd = args$noise_sd_channel,
      seed = seed_rep,
      age_years = age_years
    )
    combined <- generate_combined(args$n_subjects, age_years = age_years)
    merged <- build_merged_input(homer, combined)
    merged_csv <- file.path(tmp, sprintf("merged_rep_%03d.csv", rep_idx))
    write_csv(merged, merged_csv)

    # Channelwise null run.
    c_main <- file.path(tmp, sprintf("channel_main_%03d.csv", rep_idx))
    c_tidy <- file.path(tmp, sprintf("channel_tidy_%03d.csv", rep_idx))
    c_posthoc <- file.path(tmp, sprintf("channel_posthoc_%03d.csv", rep_idx))
    run_c <- run_analysis(
      script = "analyze_format_content_lmm_channelwise.R",
      args = c(
        "--input_csv", merged_csv,
        "--exclude_subjects_json", exclude_none,
        "--out_main_csv", c_main,
        "--out_main_tidy_csv", c_tidy,
        "--out_posthoc_csv", c_posthoc,
        "--alpha", as.character(args$alpha),
        "--min_subjects", as.character(args$min_subjects)
      ),
      env = c("FILTER_MAIN_EFFECTS_TO_FDR_SIG_ONLY=FALSE")
    )
    assert_true(run_c$status == 0, paste0("channelwise run failed at replicate ", rep_idx))
    assert_output_file(c_tidy, run_c, "channelwise", rep_idx)
    channel_tidy <- read_csv(c_tidy, show_col_types = FALSE)
    for (eff in c("format", "content", "interaction")) {
      vals <- channel_tidy %>% filter(effect == eff) %>% pull(p_fdr)
      effect_rows <- append_effect_rows(effect_rows, "channelwise", rep_idx, eff, vals, args$alpha)
    }

    # ROI null run.
    r_main <- file.path(tmp, sprintf("roi_main_%03d.csv", rep_idx))
    r_tidy <- file.path(tmp, sprintf("roi_tidy_%03d.csv", rep_idx))
    r_posthoc <- file.path(tmp, sprintf("roi_posthoc_%03d.csv", rep_idx))
    run_r <- run_analysis(
      script = "analyze_format_content_lmm_roi.R",
      args = c(
        "--input_csv", merged_csv,
        "--roi_json", roi_json,
        "--exclude_subjects_json", exclude_none,
        "--out_main_csv", r_main,
        "--out_main_tidy_csv", r_tidy,
        "--out_posthoc_csv", r_posthoc,
        "--alpha", as.character(args$alpha),
        "--min_subjects", as.character(args$min_subjects)
      ),
      env = c("FILTER_MAIN_EFFECTS_TO_FDR_SIG_ONLY=FALSE")
    )
    assert_true(run_r$status == 0, paste0("ROI run failed at replicate ", rep_idx))
    assert_output_file(r_tidy, run_r, "ROI", rep_idx)
    roi_tidy <- read_csv(r_tidy, show_col_types = FALSE)
    for (eff in c("format", "content", "interaction")) {
      vals <- roi_tidy %>% filter(effect == eff) %>% pull(p_fdr)
      effect_rows <- append_effect_rows(effect_rows, "roi", rep_idx, eff, vals, args$alpha)
    }

    # Retention null run.
    retention_df <- generate_retention_null(args$n_subjects, args$noise_sd_retention, seed_rep + 100000L, age_years = age_years)
    retention_csv <- file.path(tmp, sprintf("retention_%03d.csv", rep_idx))
    ret_main <- file.path(tmp, sprintf("ret_main_%03d.csv", rep_idx))
    ret_posthoc <- file.path(tmp, sprintf("ret_posthoc_%03d.csv", rep_idx))
    write_csv(retention_df, retention_csv)
    run_ret <- run_analysis(
      script = "analyze_retention_format_content_lmm.R",
      args = c(
        "--input_csv", retention_csv,
        "--exclude_subjects_json", exclude_none,
        "--out_main_csv", ret_main,
        "--out_posthoc_csv", ret_posthoc,
        "--alpha", as.character(args$alpha),
        "--min_subjects", as.character(args$min_subjects)
      )
    )
    assert_true(run_ret$status == 0, paste0("retention run failed at replicate ", rep_idx))
    assert_output_file(ret_main, run_ret, "retention", rep_idx)
    ret_df <- read_csv(ret_main, show_col_types = FALSE)
    for (eff in c("length", "content", "interaction")) {
      vals <- ret_df %>% filter(effect == eff) %>% pull(p_fdr)
      effect_rows <- append_effect_rows(effect_rows, "retention", rep_idx, eff, vals, args$alpha)
    }

    # Engagement null run.
    engagement_df <- generate_engagement_null(args$n_subjects, args$noise_sd_engagement, seed_rep + 200000L, age_years = age_years)
    engagement_csv <- file.path(tmp, sprintf("engagement_%03d.csv", rep_idx))
    eng_main <- file.path(tmp, sprintf("eng_main_%03d.csv", rep_idx))
    eng_posthoc <- file.path(tmp, sprintf("eng_posthoc_%03d.csv", rep_idx))
    write_csv(engagement_df, engagement_csv)
    run_eng <- run_analysis(
      script = "analyze_engagement_format_content_lmm.R",
      args = c(
        "--input_csv", engagement_csv,
        "--exclude_subjects_json", exclude_none,
        "--out_main_csv", eng_main,
        "--out_posthoc_csv", eng_posthoc,
        "--alpha", as.character(args$alpha),
        "--min_subjects", as.character(args$min_subjects)
      )
    )
    assert_true(run_eng$status == 0, paste0("engagement run failed at replicate ", rep_idx))
    assert_output_file(eng_main, run_eng, "engagement", rep_idx)
    eng_df <- read_csv(eng_main, show_col_types = FALSE)
    for (eff in c("length", "content", "interaction")) {
      vals <- eng_df %>% filter(effect == eff) %>% pull(p_fdr)
      effect_rows <- append_effect_rows(effect_rows, "engagement", rep_idx, eff, vals, args$alpha)
    }
  }

  all_effects <- bind_rows(effect_rows)
  summary_df <- all_effects %>%
    group_by(pipeline, effect) %>%
    summarize(
      n_replicates = n(),
      total_tests = sum(.data$n_tests),
      total_reject = sum(.data$n_reject),
      empirical_type1_rate = .data$total_reject / .data$total_tests,
      .groups = "drop"
    ) %>%
    arrange(.data$pipeline, .data$effect)

  dir.create(dirname(args$out_csv), recursive = TRUE, showWarnings = FALSE)
  write_csv(summary_df, args$out_csv)

  cat("[write] calibration summary:", args$out_csv, "\n")
  print(summary_df, n = nrow(summary_df))

  viol <- summary_df %>% filter(.data$empirical_type1_rate > args$type1_upper_bound)
  if (nrow(viol) > 0) {
    fail(
      paste0(
        "Type-I calibration threshold exceeded (upper_bound=",
        args$type1_upper_bound,
        "). Offending rows: ",
        paste0(viol$pipeline, "/", viol$effect, "=", signif(viol$empirical_type1_rate, 4), collapse = "; ")
      )
    )
  }

  cat(
    "[PASS] empirical type-I rates are <= ", args$type1_upper_bound,
    " for all pipelines/effects (alpha=", args$alpha, ", n_reps=", args$n_reps, ").\n",
    sep = ""
  )
}

main()
