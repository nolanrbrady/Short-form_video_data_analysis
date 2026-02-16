#!/usr/bin/env Rscript

# Monte Carlo type-II error (power) calibration for all R inferential pipelines.
#
# Purpose
# - Run repeated non-null synthetic datasets through:
#   - analyze_format_content_lmm_channelwise.R
#   - analyze_format_content_lmm_roi.R
#   - analyze_retention_format_content_lmm.R
#   - analyze_engagement_format_content_lmm.R
# - Estimate empirical detection power from adjusted p-values.
# - Fail when observed type-II error (1 - power) exceeds a configured upper bound.
#
# Methodological rationale
# - Simulation-based calibration is a standard approach for evaluating operating
#   characteristics (power / type-II behavior) of statistical procedures under
#   controlled data-generating mechanisms (Burton et al., 2006; see CITATIONS.md).
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
    type2_upper_bound = 0.40,
    n_subjects = 32L,
    noise_sd_channel = 0.20,
    noise_sd_retention = 0.10,
    noise_sd_engagement = 0.12,
    min_subjects = 6L,
    seed = 20260214L,
    out_csv = "data/results/type2_error_calibration_summary_r.csv",
    # Non-null effect sizes used for calibration scenarios.
    effect_channel_format = 0.45,
    effect_channel_content = -0.35,
    effect_channel_interaction = 0.40,
    effect_retention_length = 0.45,
    effect_retention_content = -0.30,
    effect_retention_interaction = 0.50,
    effect_engagement_length = -0.40,
    effect_engagement_content = 0.35,
    effect_engagement_interaction = 0.45
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
  int_keys <- c("n_reps", "n_subjects", "min_subjects", "seed")
  num_keys <- c(
    "alpha", "type2_upper_bound",
    "noise_sd_channel", "noise_sd_retention", "noise_sd_engagement",
    "effect_channel_format", "effect_channel_content", "effect_channel_interaction",
    "effect_retention_length", "effect_retention_content", "effect_retention_interaction",
    "effect_engagement_length", "effect_engagement_content", "effect_engagement_interaction"
  )
  for (k in int_keys) parsed[[k]] <- as.integer(parsed[[k]])
  for (k in num_keys) parsed[[k]] <- as.numeric(parsed[[k]])
  parsed
}

make_subject_ids <- function(n_subjects) {
  ids_int <- seq_len(n_subjects)
  list(
    homer_subject = sprintf("sub_%04d", ids_int),
    combined_subject = sprintf("%04d", ids_int)
  )
}

cond_grid <- function() {
  tibble::tibble(
    cond = c("01", "02", "03", "04"),
    condition = c("SF_Edu", "SF_Ent", "LF_Ent", "LF_Edu"),
    format_c = c(-0.5, -0.5, +0.5, +0.5),
    content_c = c(+0.5, -0.5, -0.5, +0.5)
  )
}

generate_homer3_wide_alt <- function(
  n_subjects,
  channels,
  chroms,
  bF,
  bC,
  bI,
  noise_sd,
  seed
) {
  set.seed(seed)
  ids <- make_subject_ids(n_subjects)
  grid <- cond_grid()
  subj_re <- rnorm(n_subjects, mean = 0, sd = 0.6)

  out <- tibble::tibble(Subject = ids$homer_subject)
  for (ch in channels) {
    for (chrom in chroms) {
      b0 <- rnorm(1, mean = 0, sd = 0.3)
      bF_j <- bF + rnorm(1, mean = 0, sd = 0.05)
      bC_j <- bC + rnorm(1, mean = 0, sd = 0.05)
      bI_j <- bI + rnorm(1, mean = 0, sd = 0.05)
      for (i in seq_len(nrow(grid))) {
        mu <- b0 +
          bF_j * grid$format_c[[i]] +
          bC_j * grid$content_c[[i]] +
          bI_j * (grid$format_c[[i]] * grid$content_c[[i]])
        col_name <- paste0(ch, "_Cond", grid$cond[[i]], "_", chrom)
        out[[col_name]] <- mu + subj_re + rnorm(n_subjects, mean = 0, sd = noise_sd)
      }
    }
  }
  out
}

generate_combined <- function(n_subjects) {
  ids <- make_subject_ids(n_subjects)
  tibble::tibble(
    subject_id = ids$combined_subject,
    age = rep(20, n_subjects),
    pd_status = rep(0, n_subjects)
  )
}

build_merged_input <- function(homer, combined) {
  homer_norm <- homer %>%
    mutate(subject_id_norm = as.integer(str_extract(.data$Subject, "\\d+"))) %>%
    select(-Subject)
  combined_norm <- combined %>%
    mutate(subject_id_norm = as.integer(str_extract(.data$subject_id, "\\d+")))

  combined_norm %>%
    inner_join(homer_norm, by = "subject_id_norm", relationship = "one-to-one") %>%
    select(-subject_id_norm)
}

generate_retention_alt <- function(n_subjects, bL, bC, bI, noise_sd, seed) {
  set.seed(seed)
  sid <- sprintf("%04d", seq_len(n_subjects))
  subj_re <- rnorm(n_subjects, mean = 0, sd = 0.35)
  b0 <- 0.30
  grid <- cond_grid()
  make_cell <- function(format_c, content_c) {
    mu <- b0 + bL * format_c + bC * content_c + bI * (format_c * content_c)
    mu + subj_re + rnorm(n_subjects, mean = 0, sd = noise_sd)
  }
  tibble::tibble(
    subject_id = sid,
    diff_short_form_education = make_cell(grid$format_c[[1]], grid$content_c[[1]]),
    diff_short_form_entertainment = make_cell(grid$format_c[[2]], grid$content_c[[2]]),
    diff_long_form_entertainment = make_cell(grid$format_c[[3]], grid$content_c[[3]]),
    diff_long_form_education = make_cell(grid$format_c[[4]], grid$content_c[[4]])
  )
}

generate_engagement_alt <- function(n_subjects, bL, bC, bI, noise_sd, seed) {
  set.seed(seed)
  sid <- sprintf("%04d", seq_len(n_subjects))
  subj_re <- rnorm(n_subjects, mean = 0, sd = 0.50)
  b0 <- 3.0
  grid <- cond_grid()
  make_cell <- function(format_c, content_c) {
    mu <- b0 + bL * format_c + bC * content_c + bI * (format_c * content_c)
    mu + subj_re + rnorm(n_subjects, mean = 0, sd = noise_sd)
  }
  tibble::tibble(
    subject_id = sid,
    sf_education_engagement = make_cell(grid$format_c[[1]], grid$content_c[[1]]),
    sf_entertainment_engagement = make_cell(grid$format_c[[2]], grid$content_c[[2]]),
    lf_entertainment_engagement = make_cell(grid$format_c[[3]], grid$content_c[[3]]),
    lf_education_engagement = make_cell(grid$format_c[[4]], grid$content_c[[4]])
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
    n_true_alt_tests = length(p_fdr_values),
    n_detected = sum(p_fdr_values < alpha, na.rm = TRUE)
  )
  acc
}

main <- function() {
  args <- parse_args()
  if (args$n_reps < 10) stop("n_reps must be >= 10 for a useful calibration run.")
  if (args$n_subjects < args$min_subjects) stop("n_subjects must be >= min_subjects.")

  tmp <- file.path(tempdir(), "type2_error_calibration_r")
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

    homer <- generate_homer3_wide_alt(
      n_subjects = args$n_subjects,
      channels = channels,
      chroms = chroms,
      bF = args$effect_channel_format,
      bC = args$effect_channel_content,
      bI = args$effect_channel_interaction,
      noise_sd = args$noise_sd_channel,
      seed = seed_rep
    )
    combined <- generate_combined(args$n_subjects)
    merged <- build_merged_input(homer, combined)
    merged_csv <- file.path(tmp, sprintf("merged_rep_%03d.csv", rep_idx))
    write_csv(merged, merged_csv)

    # Channelwise alternative run.
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

    # ROI alternative run.
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

    # Retention alternative run.
    retention_df <- generate_retention_alt(
      n_subjects = args$n_subjects,
      bL = args$effect_retention_length,
      bC = args$effect_retention_content,
      bI = args$effect_retention_interaction,
      noise_sd = args$noise_sd_retention,
      seed = seed_rep + 100000L
    )
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

    # Engagement alternative run.
    engagement_df <- generate_engagement_alt(
      n_subjects = args$n_subjects,
      bL = args$effect_engagement_length,
      bC = args$effect_engagement_content,
      bI = args$effect_engagement_interaction,
      noise_sd = args$noise_sd_engagement,
      seed = seed_rep + 200000L
    )
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
      total_true_alt_tests = sum(.data$n_true_alt_tests),
      total_detected = sum(.data$n_detected),
      empirical_power = .data$total_detected / .data$total_true_alt_tests,
      empirical_type2_error = 1.0 - .data$empirical_power,
      .groups = "drop"
    ) %>%
    arrange(.data$pipeline, .data$effect)

  dir.create(dirname(args$out_csv), recursive = TRUE, showWarnings = FALSE)
  write_csv(summary_df, args$out_csv)

  cat("[write] calibration summary:", args$out_csv, "\n")
  print(summary_df, n = nrow(summary_df))

  viol <- summary_df %>% filter(.data$empirical_type2_error > args$type2_upper_bound)
  if (nrow(viol) > 0) {
    fail(
      paste0(
        "Type-II calibration threshold exceeded (upper_bound=",
        args$type2_upper_bound,
        "). Offending rows: ",
        paste0(viol$pipeline, "/", viol$effect, "=", signif(viol$empirical_type2_error, 4), collapse = "; ")
      )
    )
  }

  cat(
    "[PASS] empirical type-II error is <= ", args$type2_upper_bound,
    " for all pipelines/effects (alpha=", args$alpha, ", n_reps=", args$n_reps, ").\n",
    sep = ""
  )
}

main()
