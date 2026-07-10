#!/usr/bin/env Rscript

# Real-data mixed-model diagnostics for the neural primary LMMs.
#
# This validator refits the publication-facing channelwise and ROI neural
# models in memory, compares their fixed-effect statistics against exported
# result tables, and checks non-degenerate residual/fitted-value diagnostics.
# It intentionally writes no files and never mutates `data/results`.
#
# Scientific rationale:
# Bates et al. (2015) emphasize inspecting convergence and model diagnostics
# for mixed-effects models. Kenward & Roger (1997) and Halekoh & Højsgaard
# (2014) support the denominator-df/F-test approximation used here. See
# CITATIONS.md for complete citation details.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(jsonlite)
  library(lme4)
  library(lmerTest)
})

script_path <- function() {
  args_all <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args_all, value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE))
  }
  normalizePath("tests/validate_real_model_diagnostics_r.R", mustWork = TRUE)
}

ROOT <- normalizePath(file.path(dirname(script_path()), ".."), mustWork = TRUE)

source(file.path(ROOT, "r_subject_exclusions.R"), local = TRUE)
source(file.path(ROOT, "r_lmm_convergence_helpers.R"), local = TRUE)

INPUT_CSV <- file.path(ROOT, "data/tabular/generated_data/homer3_betas_plus_combined_sfv_data_inner_join.csv")
EXCLUSIONS_JSON <- file.path(ROOT, "data/config/excluded_subjects.json")
ROI_JSON <- file.path(ROOT, "data/config/roi_definition.json")
CHANNEL_TIDY_CSV <- file.path(ROOT, "data/results/format_content_lmm_main_effects_tidy_r.csv")
ROI_TIDY_CSV <- file.path(ROOT, "data/results/format_content_lmm_roi_main_effects_tidy_r.csv")
NEURAL_LMM_RESPONSE_SCALE <- 1e6
TOL <- 1e-10

normalize_subject_id <- function(x, column_name) {
  x_chr <- as.character(x)
  extracted <- str_extract(x_chr, "\\d+")
  if (any(is.na(extracted))) {
    bad <- head(unique(x_chr[is.na(extracted)]), 10)
    stop("Failed to parse numeric IDs from column '", column_name, "'. Examples: ", paste(bad, collapse = ", "))
  }
  as.integer(extracted)
}

normalize_channel_id <- function(x, context_label) {
  x_chr <- toupper(str_trim(as.character(x)))
  m <- str_match(x_chr, "^S(\\d+)_D(\\d+)$")
  if (any(is.na(m[, 1]))) {
    bad <- head(unique(x_chr[is.na(m[, 1])]), 10)
    stop("Failed to parse channel IDs in ", context_label, ". Examples: ", paste(bad, collapse = ", "))
  }
  paste0("S", sprintf("%02d", as.integer(m[, 2])), "_D", sprintf("%02d", as.integer(m[, 3])))
}

load_input <- function() {
  df <- read_csv(INPUT_CSV, show_col_types = FALSE)
  required <- c("subject_id", "age", "education_years")
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop("Merged input is missing required columns: ", paste(missing, collapse = ", "))
  }
  df <- df %>%
    mutate(
      subject_id = normalize_subject_id(.data$subject_id, "subject_id"),
      age = suppressWarnings(as.numeric(.data$age)),
      education_years = suppressWarnings(as.numeric(.data$education_years))
    )
  if (any(is.na(df$age)) || any(is.na(df$education_years))) {
    stop("Age or education_years contains missing/non-numeric values after coercion.")
  }
  apply_subject_exclusions(df, "subject_id", EXCLUSIONS_JSON, "real_model_diagnostics")$data
}

reshape_channel_long <- function(df) {
  beta_cols <- names(df)[str_detect(names(df), "^S\\d+_D\\d+_Cond\\d{2}_(HbO|HbR)$")]
  if (length(beta_cols) == 0) stop("No Homer beta columns found in merged input.")

  df %>%
    select(subject_id, age, education_years, all_of(beta_cols)) %>%
    pivot_longer(all_of(beta_cols), names_to = "beta_col", values_to = "beta") %>%
    extract(
      "beta_col",
      into = c("channel", "cond", "chrom"),
      regex = "^(S\\d+_D\\d+)_Cond(\\d{2})_(HbO|HbR)$",
      remove = TRUE
    ) %>%
    mutate(
      channel = normalize_channel_id(.data$channel, "beta column names"),
      beta = suppressWarnings(as.numeric(.data$beta)),
      condition = case_when(
        .data$cond == "01" ~ "SF_Edu",
        .data$cond == "02" ~ "SF_Ent",
        .data$cond == "03" ~ "LF_Ent",
        .data$cond == "04" ~ "LF_Edu",
        TRUE ~ NA_character_
      ),
      format_c = case_when(.data$cond %in% c("01", "02") ~ -0.5, .data$cond %in% c("03", "04") ~ 0.5),
      content_c = case_when(.data$cond %in% c("02", "03") ~ -0.5, .data$cond %in% c("01", "04") ~ 0.5)
    ) %>%
    filter(!is.na(.data$condition), !is.na(.data$chrom), !is.na(.data$channel))
}

load_roi_map <- function() {
  roi_obj <- jsonlite::fromJSON(ROI_JSON, simplifyVector = FALSE)
  rows <- list()
  for (roi_name in names(roi_obj)) {
    channels <- unlist(roi_obj[[roi_name]], use.names = FALSE)
    rows[[length(rows) + 1]] <- tibble(
      roi = roi_name,
      channel = normalize_channel_id(channels, paste0("roi_definition[", roi_name, "]"))
    )
  }
  bind_rows(rows)
}

aggregate_roi_long <- function(channel_long, roi_map) {
  channel_long %>%
    inner_join(roi_map, by = "channel", relationship = "many-to-one") %>%
    group_by(subject_id, age, education_years, roi, chrom, condition, format_c, content_c) %>%
    summarize(beta = if (all(is.na(beta))) NA_real_ else mean(beta, na.rm = TRUE), .groups = "drop")
}

complete_case_subjects <- function(sub) {
  non_missing <- sub %>% filter(!is.na(.data$beta))
  keep_ids <- non_missing %>%
    group_by(subject_id) %>%
    summarize(n_cond = n_distinct(.data$condition), .groups = "drop") %>%
    filter(.data$n_cond == 4) %>%
    pull(.data$subject_id)
  non_missing %>% filter(.data$subject_id %in% keep_ids)
}

fit_factorial <- function(sub_complete) {
  sub_complete <- sub_complete %>% mutate(beta = .data$beta * NEURAL_LMM_RESPONSE_SCALE)
  lmerTest::lmer(
    beta ~ format_c * content_c + age + education_years + (1 | subject_id),
    data = sub_complete,
    REML = TRUE
  )
}

back_transform <- function(x) {
  as.numeric(x) / NEURAL_LMM_RESPONSE_SCALE
}

extract_effects <- function(model, anova_kr) {
  fixed <- lme4::fixef(model)
  bind_rows(lapply(
    c(format = "format_c", content = "content_c", interaction = "format_c:content_c"),
    function(term) {
      est <- as.numeric(fixed[[term]])
      f_val <- as.numeric(anova_kr[term, "F value"])
      df <- as.numeric(anova_kr[term, "DenDF"])
      se <- if (f_val > 0) {
        abs(est) / sqrt(f_val)
      } else {
        as.numeric(summary(model, ddf = "Kenward-Roger")$coefficients[term, "Std. Error"])
      }
      t_val <- sign(est) * sqrt(f_val)
      ci_half <- qt(0.975, df = df) * se
      tibble(
        effect = names(which(c(format = "format_c", content = "content_c", interaction = "format_c:content_c") == term)),
        estimate = back_transform(est),
        se = back_transform(se),
        df = df,
        t = t_val,
        p_unc = as.numeric(anova_kr[term, "Pr(>F)"]),
        ci95_low = back_transform(est - ci_half),
        ci95_high = back_transform(est + ci_half)
      )
    }
  ))
}

assert_close <- function(actual, expected, label, tolerance = TOL) {
  if (!isTRUE(all.equal(actual, expected, tolerance = tolerance, scale = 1))) {
    stop(label, " mismatch: actual=", actual, ", expected=", expected)
  }
}

diagnose_model <- function(model) {
  residuals <- resid(model)
  fitted_values <- fitted(model)
  if (!all(is.finite(residuals)) || !all(is.finite(fitted_values))) {
    stop("Model produced non-finite residuals or fitted values.")
  }
  residual_sd <- sd(residuals)
  fitted_sd <- sd(fitted_values)
  if (!is.finite(residual_sd) || residual_sd <= 0) stop("Model residuals are degenerate.")
  if (!is.finite(fitted_sd) || fitted_sd <= 0) stop("Model fitted values are degenerate.")
  shapiro_p <- if (length(residuals) >= 3 && length(residuals) <= 5000) shapiro.test(residuals)$p.value else NA_real_
  tibble(
    residual_sd = residual_sd,
    fitted_sd = fitted_sd,
    max_abs_standardized_residual = max(abs(as.numeric(scale(residuals)))),
    shapiro_p = shapiro_p
  )
}

audit_family <- function(long_df, exported_tidy, unit_col, label) {
  units <- sort(unique(long_df[[unit_col]]))
  chroms <- c("HbO", "HbR")
  diagnostics <- list()
  n_checked <- 0L

  for (chrom_name in chroms) {
    for (unit_name in units) {
      sub <- long_df %>% filter(.data$chrom == chrom_name, .data[[unit_col]] == unit_name)
      sub_cc <- complete_case_subjects(sub)
      if (length(unique(sub_cc$subject_id)) < 6) next

      fit_result <- capture_lmm_fit(function() fit_factorial(sub_cc))
      model <- fit_result$model
      singular_fit <- lme4::isSingular(model, tol = 1e-4)
      anova_kr <- anova(model, ddf = "Kenward-Roger", type = 3)
      effects <- extract_effects(model, anova_kr)

      exported_rows <- exported_tidy %>%
        filter(.data$chrom == chrom_name, .data[[unit_col]] == unit_name)
      if (nrow(exported_rows) != 3) {
        stop(label, " exported table has ", nrow(exported_rows), " rows for ", unit_name, " ", chrom_name, "; expected 3.")
      }

      for (i in seq_len(nrow(effects))) {
        eff <- effects[i, ]
        exported <- exported_rows %>% filter(.data$effect == eff$effect)
        if (nrow(exported) != 1) {
          stop(label, " exported table missing unique effect row for ", unit_name, " ", chrom_name, " ", eff$effect)
        }
        assert_close(exported$n_subjects, length(unique(sub_cc$subject_id)), paste(label, unit_name, chrom_name, eff$effect, "n_subjects"))
        assert_close(exported$n_obs, nrow(sub_cc), paste(label, unit_name, chrom_name, eff$effect, "n_obs"))
        if (!identical(as.logical(exported$converged), fit_result$converged)) {
          stop(label, " convergence flag mismatch for ", unit_name, " ", chrom_name)
        }
        if (!identical(as.logical(exported$singular_fit), singular_fit)) {
          stop(label, " singular-fit flag mismatch for ", unit_name, " ", chrom_name)
        }
        for (metric in c("estimate", "se", "df", "t", "p_unc", "ci95_low", "ci95_high")) {
          assert_close(exported[[metric]], eff[[metric]], paste(label, unit_name, chrom_name, eff$effect, metric))
        }
      }

      diagnostics[[length(diagnostics) + 1]] <- diagnose_model(model) %>%
        mutate(
          analysis = label,
          unit = unit_name,
          chrom = chrom_name,
          n_subjects = length(unique(sub_cc$subject_id)),
          n_obs = nrow(sub_cc),
          converged = fit_result$converged,
          singular_fit = singular_fit
        )
      n_checked <- n_checked + 1L
    }
  }

  if (length(diagnostics) == 0) stop("No ", label, " real-data models were audited.")
  diag_df <- bind_rows(diagnostics)
  if (any(!diag_df$converged)) stop(label, " contains non-converged real-data models.")
  list(n_models = n_checked, diagnostics = diag_df)
}

main <- function() {
  if (!requireNamespace("pbkrtest", quietly = TRUE)) {
    stop("Kenward-Roger diagnostics require the pbkrtest package.")
  }

  input <- load_input()
  channel_long <- reshape_channel_long(input)
  roi_long <- aggregate_roi_long(channel_long, load_roi_map())
  channel_exported <- read_csv(CHANNEL_TIDY_CSV, show_col_types = FALSE)
  roi_exported <- read_csv(ROI_TIDY_CSV, show_col_types = FALSE)

  channel_audit <- audit_family(channel_long, channel_exported, "channel", "channelwise")
  roi_audit <- audit_family(roi_long, roi_exported, "roi", "roi")
  diagnostics <- bind_rows(channel_audit$diagnostics, roi_audit$diagnostics)

  cat("[diagnostics] channelwise models audited:", channel_audit$n_models, "\n")
  cat("[diagnostics] ROI models audited:", roi_audit$n_models, "\n")
  cat("[diagnostics] max absolute standardized residual:", max(diagnostics$max_abs_standardized_residual), "\n")
  cat("[diagnostics] minimum residual SD:", min(diagnostics$residual_sd), "\n")
  cat("[diagnostics] singular fits:", sum(diagnostics$singular_fit), "\n")
  cat("[diagnostics] minimum Shapiro-Wilk residual p-value:", min(diagnostics$shapiro_p, na.rm = TRUE), "\n")
  cat("[PASS] validate_real_model_diagnostics_r\n")
}

main()
