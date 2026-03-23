#!/usr/bin/env Rscript

# Validation harness for engagement Length x Content LMM (R).
#
# Purpose
# - Verify `analyze_engagement_format_content_lmm.R` implements:
#   - Required-column and subject_id integrity checks (fail hard)
#   - 2x2 condition mapping and +/-0.5 effect coding
#   - age-adjusted omnibus model (`+ age`) with fail-hard covariate validation
#   - Complete-case rule (subject must have all 4 engagement values)
#   - Engagement zero handling as valid data (not missing)
#   - Holm correction across the 3 planned omnibus effects
#   - Post-hoc gating based on interaction adjusted p-value threshold
#
# This script generates synthetic datasets in /tmp and runs the analysis end-to-end.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
})

fail <- function(msg) {
  writeLines(paste0("[FAIL] ", msg))
  quit(status = 1)
}

assert_true <- function(cond, msg) {
  if (!isTRUE(cond)) fail(msg)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

source("r_emmeans_posthoc_helpers.R", local = TRUE)

make_age_years <- function(n_subjects) {
  seq(18, by = 1, length.out = n_subjects)
}

holm_adjust_manual <- function(p_values) {
  # Manual Holm adjusted p-values.
  # Reference: Holm (1979), see `CITATIONS.md`.
  p <- as.numeric(p_values)
  out <- rep(NA_real_, length(p))
  finite <- is.finite(p)
  if (!any(finite)) return(out)

  p_f <- p[finite]
  m <- length(p_f)
  o <- order(p_f)
  ranked <- p_f[o]
  adj_ranked <- rep(NA_real_, m)
  prev <- 0.0
  for (i in seq_len(m)) {
    val <- (m - i + 1) * ranked[[i]]
    prev <- max(prev, val)
    adj_ranked[[i]] <- min(prev, 1.0)
  }
  adj <- rep(NA_real_, m)
  adj[o] <- adj_ranked
  out[finite] <- adj
  out
}

generate_engagement <- function(n_subjects, b0, bL, bC, bI, b_age = 0.10, noise_sd = 0.0, seed = 123, age_years = make_age_years(n_subjects)) {
  set.seed(seed)
  sid <- sprintf("%04d", seq_len(n_subjects))
  subj_re <- seq(-0.50, 0.50, length.out = n_subjects)
  age_centered <- age_years - mean(age_years)

  make_cell <- function(length_c, content_c) {
    mu <- b0 + bL * length_c + bC * content_c + bI * (length_c * content_c) + b_age * age_centered
    mu + subj_re + rnorm(n_subjects, mean = 0, sd = noise_sd)
  }

  tibble::tibble(
    subject_id = sid,
    age = age_years,
    sf_education_engagement = make_cell(-0.5, +0.5),
    sf_entertainment_engagement = make_cell(-0.5, -0.5),
    lf_education_engagement = make_cell(+0.5, +0.5),
    lf_entertainment_engagement = make_cell(+0.5, -0.5)
  )
}

compute_effects_from_cell_means <- function(df) {
  m_sfe <- mean(df$sf_education_engagement)
  m_sfeent <- mean(df$sf_entertainment_engagement)
  m_lfe <- mean(df$lf_education_engagement)
  m_lfent <- mean(df$lf_entertainment_engagement)

  length_eff <- mean(c(m_lfe, m_lfent)) - mean(c(m_sfe, m_sfeent))
  content_eff <- mean(c(m_sfe, m_lfe)) - mean(c(m_sfeent, m_lfent))
  interaction_eff <- (m_lfe - m_lfent) - (m_sfe - m_sfeent)
  list(length = length_eff, content = content_eff, interaction = interaction_eff)
}

run_analysis <- function(input_csv, out_main, out_posthoc, alpha = 0.05, min_subjects = 6, exclude_json = NULL) {
  script <- "analyze_engagement_format_content_lmm.R"
  args <- c(
    "--input_csv", input_csv,
    "--out_main_csv", out_main,
    "--out_posthoc_csv", out_posthoc,
    "--alpha", as.character(alpha),
    "--min_subjects", as.character(min_subjects)
  )
  if (!is.null(exclude_json)) {
    args <- c(args, "--exclude_subjects_json", exclude_json)
  }
  output <- suppressWarnings(system2("Rscript", c(script, args), stdout = TRUE, stderr = TRUE))
  list(status = attr(output, "status") %||% 0, output = output)
}

main <- function() {
  tmp <- file.path(tempdir(), "engagement_pipeline_validation")
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  exclude_none <- file.path(tmp, "excluded_none.json")
  exclude_one <- file.path(tmp, "excluded_one.json")
  exclude_absent <- file.path(tmp, "excluded_absent.json")
  writeLines("[]", exclude_none)
  writeLines('["sub_0001"]', exclude_one)
  writeLines('["sub_9999"]', exclude_absent)

  # 0) Posthoc formatting helper must preserve t-vs-z semantics.
  helper_t <- standardize_emmeans_pairwise_output(
    tibble::tibble(
      contrast = "SF_Edu - SF_Ent",
      estimate = 1.25,
      SE = 0.40,
      df = 18,
      t.ratio = 3.125,
      p.value = 0.006
    ),
    context_label = "engagement helper t-branch"
  )
  assert_true(helper_t$stat_type[[1]] == "t", "expected t.ratio helper branch to record stat_type=t")
  assert_true(helper_t$t[[1]] == "3.125", "t.ratio helper branch should preserve the reported t statistic")

  helper_z <- standardize_emmeans_pairwise_output(
    tibble::tibble(
      contrast = "LF_Edu - LF_Ent",
      estimate = 0.75,
      SE = 0.20,
      z.ratio = 3.75,
      p.value = 0.001
    ),
    context_label = "engagement helper z-branch"
  )
  assert_true(helper_z$stat_type[[1]] == "z", "expected z.ratio helper branch to record stat_type=z")
  assert_true(helper_z$t[[1]] == "t can't be calculated", "z.ratio helper branch should not relabel z as t")
  assert_true(is.na(helper_z$df[[1]]), "z.ratio helper branch should fill missing df with NA")

  # 1) Deterministic analytic ground truth (no noise): exact algebra checks.
  det <- generate_engagement(
    n_subjects = 20,
    b0 = 3.0,
    bL = -0.20,
    bC = 0.35,
    bI = 0.80,
    b_age = 0.12,
    noise_sd = 0.0,
    seed = 1
  )
  eff_det <- compute_effects_from_cell_means(det)
  assert_true(abs(eff_det$length - (-0.20)) < 1e-12, "analytic length effect mismatch")
  assert_true(abs(eff_det$content - 0.35) < 1e-12, "analytic content effect mismatch")
  assert_true(abs(eff_det$interaction - 0.80) < 1e-12, "analytic interaction effect mismatch")

  # 2) End-to-end run with stochastic noise and known generating parameters.
  n_subjects <- 24
  truth <- list(length = -0.25, content = 0.30, interaction = 0.75)
  age_years <- make_age_years(n_subjects)
  df <- generate_engagement(
    n_subjects = n_subjects,
    b0 = 3.0,
    bL = truth$length,
    bC = truth$content,
    bI = truth$interaction,
    b_age = 0.14,
    noise_sd = 0.08,
    seed = 42,
    age_years = age_years
  )
  input_csv <- file.path(tmp, "engagement.csv")
  out_main <- file.path(tmp, "main.csv")
  out_posthoc <- file.path(tmp, "posthoc.csv")
  write_csv(df, input_csv)

  run <- run_analysis(
    input_csv,
    out_main,
    out_posthoc,
    alpha = 0.05,
    min_subjects = 6,
    exclude_json = exclude_none
  )
  assert_true(run$status == 0, "analysis script did not exit cleanly on primary synthetic data")
  assert_true(file.exists(out_main), "missing main effects output CSV")
  assert_true(file.exists(out_posthoc), "missing posthoc output CSV")

  main_df <- read_csv(out_main, show_col_types = FALSE)
  posthoc_df <- read_csv(out_posthoc, show_col_types = FALSE)

  required_main <- c(
    "effect", "n_subjects", "n_obs", "singular_fit", "estimate", "se", "df",
    "t", "p_unc", "p_fdr", "ci95_low", "ci95_high"
  )
  assert_true(all(required_main %in% names(main_df)), "main output missing required columns")
  assert_true(nrow(main_df) == 3, "expected exactly 3 omnibus effect rows")
  assert_true(setequal(main_df$effect, c("length", "content", "interaction")), "unexpected effect labels")
  assert_true(all(main_df$n_subjects == n_subjects), "expected full subject count in main output")
  assert_true(all(main_df$n_obs == 4 * n_subjects), "expected n_obs == 4 * n_subjects")

  tol <- 0.20
  est_length <- main_df %>% filter(effect == "length") %>% pull(estimate)
  est_content <- main_df %>% filter(effect == "content") %>% pull(estimate)
  est_inter <- main_df %>% filter(effect == "interaction") %>% pull(estimate)
  assert_true(abs(est_length[[1]] - truth$length) < tol, "length estimate not close to generating truth")
  assert_true(abs(est_content[[1]] - truth$content) < tol, "content estimate not close to generating truth")
  assert_true(abs(est_inter[[1]] - truth$interaction) < tol, "interaction estimate not close to generating truth")

  # 2a) Direct reference-model agreement for the age-adjusted omnibus fit.
  ref_long <- df %>%
    select(subject_id, age, sf_education_engagement, sf_entertainment_engagement, lf_education_engagement, lf_entertainment_engagement) %>%
    pivot_longer(
      cols = c(sf_education_engagement, sf_entertainment_engagement, lf_education_engagement, lf_entertainment_engagement),
      names_to = "eng_col",
      values_to = "engagement"
    ) %>%
    mutate(
      condition = case_when(
        eng_col == "sf_education_engagement" ~ "SF_Edu",
        eng_col == "sf_entertainment_engagement" ~ "SF_Ent",
        eng_col == "lf_entertainment_engagement" ~ "LF_Ent",
        eng_col == "lf_education_engagement" ~ "LF_Edu",
        TRUE ~ NA_character_
      ),
      length_c = case_when(
        condition %in% c("SF_Edu", "SF_Ent") ~ -0.5,
        condition %in% c("LF_Ent", "LF_Edu") ~ +0.5,
        TRUE ~ NA_real_
      ),
      content_c = case_when(
        condition %in% c("SF_Ent", "LF_Ent") ~ -0.5,
        condition %in% c("SF_Edu", "LF_Edu") ~ +0.5,
        TRUE ~ NA_real_
      )
    )
  ref_model <- lmerTest::lmer(engagement ~ length_c * content_c + age + (1 | subject_id), data = ref_long, REML = TRUE)
  ref_coefs <- summary(ref_model)$coefficients
  for (term_map in list(c("length", "length_c"), c("content", "content_c"), c("interaction", "length_c:content_c"))) {
    out_row <- main_df %>% filter(effect == term_map[[1]])
    assert_true(nrow(out_row) == 1, paste0("missing output row for ", term_map[[1]], " reference fit"))
    assert_true(abs(out_row$estimate[[1]] - ref_coefs[term_map[[2]], "Estimate"]) < 1e-10, paste0("estimate mismatch for ", term_map[[1]], " reference fit"))
    assert_true(abs(out_row$se[[1]] - ref_coefs[term_map[[2]], "Std. Error"]) < 1e-10, paste0("SE mismatch for ", term_map[[1]], " reference fit"))
    assert_true(abs(out_row$p_unc[[1]] - ref_coefs[term_map[[2]], "Pr(>|t|)"]) < 1e-10, paste0("p_unc mismatch for ", term_map[[1]], " reference fit"))
  }

  # 2b) Explicit inferential TP check on strong-effect synthetic data.
  p_inter <- main_df %>% filter(effect == "interaction") %>% pull(p_fdr)
  assert_true(length(p_inter) == 1, "missing interaction p_fdr row in TP run")
  assert_true(p_inter[[1]] < 0.05, "expected interaction to be significant in TP synthetic dataset")

  # 2c) Central exclusion list: excluding one present subject should reduce analyzed N by exactly 1.
  out_main_excl <- file.path(tmp, "main_excluded_one.csv")
  out_posthoc_excl <- file.path(tmp, "posthoc_excluded_one.csv")
  run_excl <- run_analysis(
    input_csv,
    out_main_excl,
    out_posthoc_excl,
    alpha = 0.05,
    min_subjects = 6,
    exclude_json = exclude_one
  )
  assert_true(run_excl$status == 0, "analysis failed with one-subject exclusion list")
  main_excl <- read_csv(out_main_excl, show_col_types = FALSE)
  assert_true(all(main_excl$n_subjects == (n_subjects - 1)), "one excluded subject should reduce n_subjects by 1")

  # 2d) Missing exclusion IDs should warn/continue and preserve N.
  out_main_abs <- file.path(tmp, "main_excluded_absent.csv")
  out_posthoc_abs <- file.path(tmp, "posthoc_excluded_absent.csv")
  run_abs <- run_analysis(
    input_csv,
    out_main_abs,
    out_posthoc_abs,
    alpha = 0.05,
    min_subjects = 6,
    exclude_json = exclude_absent
  )
  assert_true(run_abs$status == 0, "analysis failed when exclusion ID was absent")
  main_abs <- read_csv(out_main_abs, show_col_types = FALSE)
  assert_true(all(main_abs$n_subjects == n_subjects), "absent exclusion ID should not change n_subjects")

  # 2e) Explicit inferential TN check with a null-effect synthetic dataset.
  df_null <- generate_engagement(
    n_subjects = 30,
    b0 = 3.0,
    bL = 0.00,
    bC = 0.00,
    bI = 0.00,
    b_age = 0.14,
    noise_sd = 0.12,
    seed = 888
  )
  input_null <- file.path(tmp, "engagement_null.csv")
  out_main_null <- file.path(tmp, "main_null.csv")
  out_posthoc_null <- file.path(tmp, "posthoc_null.csv")
  write_csv(df_null, input_null)
  run_null <- run_analysis(
    input_null,
    out_main_null,
    out_posthoc_null,
    alpha = 0.05,
    min_subjects = 6,
    exclude_json = exclude_none
  )
  assert_true(run_null$status == 0, "analysis failed on TN null synthetic data")
  main_null <- read_csv(out_main_null, show_col_types = FALSE)
  posthoc_null <- read_csv(out_posthoc_null, show_col_types = FALSE)
  assert_true(all(main_null$p_fdr >= 0.05), "expected all omnibus effects to be non-significant in TN null dataset")
  assert_true(nrow(posthoc_null) == 0, "posthoc should be empty for TN null dataset")

  # 3) Holm correctness across exactly the three omnibus effects.
  p_adj_manual <- holm_adjust_manual(main_df$p_unc)
  max_abs <- max(abs(p_adj_manual - main_df$p_fdr), na.rm = TRUE)
  assert_true(max_abs < 1e-10, "Holm adjusted p-values mismatch in omnibus output")

  # 4) Post-hoc gating positive case: strong interaction should yield 6 pairwise rows.
  assert_true(nrow(posthoc_df) == 6, "expected 6 posthoc contrasts when interaction gate passes")
  required_post <- c("condition_a", "condition_b", "stat_type", "mean_diff", "se", "df", "t", "p_unc")
  assert_true(all(required_post %in% names(posthoc_df)), "posthoc output missing required columns")

  # 5) mean_diff sign convention check for SF_Edu - SF_Ent.
  expected_diff <- mean(df$sf_education_engagement - df$sf_entertainment_engagement)
  row_sfe <- posthoc_df %>% filter(condition_a == "SF_Edu", condition_b == "SF_Ent")
  assert_true(nrow(row_sfe) == 1, "missing SF_Edu - SF_Ent contrast")
  assert_true(abs(row_sfe$mean_diff[[1]] - expected_diff) < 0.20, "posthoc mean_diff sign/magnitude mismatch")

  # 6) Post-hoc gating negative case via strict alpha: should force no posthoc rows.
  out_main_alpha0 <- file.path(tmp, "main_alpha0.csv")
  out_posthoc_alpha0 <- file.path(tmp, "posthoc_alpha0.csv")
  run_alpha0 <- run_analysis(
    input_csv,
    out_main_alpha0,
    out_posthoc_alpha0,
    alpha = 0.0,
    min_subjects = 6,
    exclude_json = exclude_none
  )
  assert_true(run_alpha0$status == 0, "analysis failed in alpha=0 gating test")
  posthoc_alpha0 <- read_csv(out_posthoc_alpha0, show_col_types = FALSE)
  assert_true(nrow(posthoc_alpha0) == 0, "posthoc should be empty when alpha=0")

  # 7) Zero as valid engagement value: should not trigger complete-case drop.
  df_zero <- df
  df_zero$lf_entertainment_engagement[[1]] <- 0
  input_zero <- file.path(tmp, "engagement_zero_valid.csv")
  out_main_zero <- file.path(tmp, "main_zero.csv")
  out_posthoc_zero <- file.path(tmp, "posthoc_zero.csv")
  write_csv(df_zero, input_zero)
  run_zero <- run_analysis(
    input_zero,
    out_main_zero,
    out_posthoc_zero,
    alpha = 0.05,
    min_subjects = 6,
    exclude_json = exclude_none
  )
  assert_true(run_zero$status == 0, "analysis failed with valid zero engagement value")
  main_zero <- read_csv(out_main_zero, show_col_types = FALSE)
  assert_true(all(main_zero$n_subjects == n_subjects), "valid zero engagement value should not reduce subject count")

  # 8) Complete-case drop on NA in one condition.
  df_na <- df
  df_na$lf_entertainment_engagement[[1]] <- NA_real_
  input_na <- file.path(tmp, "engagement_na.csv")
  out_main_na <- file.path(tmp, "main_na.csv")
  out_posthoc_na <- file.path(tmp, "posthoc_na.csv")
  write_csv(df_na, input_na)
  run_na <- run_analysis(
    input_na,
    out_main_na,
    out_posthoc_na,
    alpha = 0.05,
    min_subjects = 6,
    exclude_json = exclude_none
  )
  assert_true(run_na$status == 0, "analysis failed for NA complete-case test")
  main_na <- read_csv(out_main_na, show_col_types = FALSE)
  assert_true(all(main_na$n_subjects == (n_subjects - 1)), "NA should drop one subject via complete-case rule")
  assert_true(all(main_na$n_obs == 4 * (n_subjects - 1)), "n_obs should reflect complete-case drop")

  # 9) Fail hard on duplicate subject IDs.
  df_dup <- bind_rows(df, df %>% slice(1))
  input_dup <- file.path(tmp, "engagement_dup.csv")
  write_csv(df_dup, input_dup)
  run_dup <- run_analysis(
    input_dup,
    file.path(tmp, "main_dup.csv"),
    file.path(tmp, "posthoc_dup.csv"),
    exclude_json = exclude_none
  )
  assert_true(run_dup$status != 0, "expected failure on duplicate subject_id")

  # 10) Fail hard on missing required column.
  df_missing <- df %>% select(-lf_education_engagement)
  input_missing <- file.path(tmp, "engagement_missing_col.csv")
  write_csv(df_missing, input_missing)
  run_missing <- run_analysis(
    input_missing,
    file.path(tmp, "main_missing.csv"),
    file.path(tmp, "posthoc_missing.csv"),
    exclude_json = exclude_none
  )
  assert_true(run_missing$status != 0, "expected failure on missing required column")

  # 11) Age is required and must be complete/numeric after exclusions.
  df_missing_age <- df %>% select(-age)
  input_missing_age <- file.path(tmp, "engagement_missing_age.csv")
  write_csv(df_missing_age, input_missing_age)
  run_missing_age <- run_analysis(
    input_missing_age,
    file.path(tmp, "main_missing_age.csv"),
    file.path(tmp, "posthoc_missing_age.csv"),
    exclude_json = exclude_none
  )
  assert_true(run_missing_age$status != 0, "expected failure on missing age column")

  df_age_na <- df
  df_age_na$age[[3]] <- NA_real_
  input_age_na <- file.path(tmp, "engagement_age_na.csv")
  write_csv(df_age_na, input_age_na)
  run_age_na <- run_analysis(
    input_age_na,
    file.path(tmp, "main_age_na.csv"),
    file.path(tmp, "posthoc_age_na.csv"),
    exclude_json = exclude_none
  )
  assert_true(run_age_na$status != 0, "expected failure on missing age values")

  df_age_bad <- df
  df_age_bad$age <- as.character(df_age_bad$age)
  df_age_bad$age[[2]] <- "not_numeric"
  input_age_bad <- file.path(tmp, "engagement_age_bad.csv")
  write_csv(df_age_bad, input_age_bad)
  run_age_bad <- run_analysis(
    input_age_bad,
    file.path(tmp, "main_age_bad.csv"),
    file.path(tmp, "posthoc_age_bad.csv"),
    exclude_json = exclude_none
  )
  assert_true(run_age_bad$status != 0, "expected failure on non-numeric age values")

  # 12) Fail hard on non-numeric engagement entry.
  df_bad <- df
  df_bad$sf_education_engagement <- as.character(df_bad$sf_education_engagement)
  df_bad$sf_education_engagement[[2]] <- "not_numeric"
  input_bad <- file.path(tmp, "engagement_non_numeric.csv")
  write_csv(df_bad, input_bad)
  run_bad <- run_analysis(
    input_bad,
    file.path(tmp, "main_bad.csv"),
    file.path(tmp, "posthoc_bad.csv"),
    exclude_json = exclude_none
  )
  assert_true(run_bad$status != 0, "expected failure on non-numeric engagement value")

  writeLines(paste0("[OK] Engagement pipeline validation passed. Outputs in: ", tmp))
}

main()
