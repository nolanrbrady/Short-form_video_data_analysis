#!/usr/bin/env Rscript

# Retention within-subject inference for Length x Content effects.
#
# Input
#   - data/results/homer3_betas_plus_combined_sfv_data_inner_join.csv
#     (must include `subject_id` plus retention columns:
#      diff_short_form_education, diff_short_form_entertainment,
#      diff_long_form_education, diff_long_form_entertainment)
#
# Design
#   - Within-subjects 2x2 factorial:
#       Length:  Short vs Long
#       Content: Education vs Entertainment
#
# Model
#   - LMM: retention_diff ~ length_c * content_c + (1|subject_id)
#     with numeric sum coding:
#       length_c  = -0.5 (Short), +0.5 (Long)
#       content_c = -0.5 (Entertainment), +0.5 (Education)
#
# Missingness policy
#   - Complete-case within subject for retention outcomes:
#     retain only subjects with all 4 retention condition values present.
#   - IMPORTANT: retention zeros are treated as valid values (not missing).
#
# Multiple testing correction
#   - Holm correction across the 3 planned omnibus effects:
#       Length, Content, Length x Content interaction.
#
# Post-hoc (only if interaction adjusted p < alpha)
#   - All pairwise comparisons among 4 condition combinations (6 contrasts),
#     uncorrected (adjust="none").
#
# Citations (see CITATIONS.md)
#   - Holm (1979): sequentially rejective multiple testing correction.
#   - Laird & Ware (1982): mixed-effects models.
#   - Bates et al. (2015): lme4 implementation.
#   - Kuznetsova et al. (2017): lmerTest / Satterthwaite df.
#   - Lenth (2016); Searle et al. (1980): EMMs and pairwise contrasts.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(lme4)
  library(lmerTest)
  library(emmeans)
})

# Shared subject-exclusion helpers (single exclusion manifest across analyses).
source_exclusion_helpers <- function() {
  args_all <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args_all, value = TRUE)
  script_dir <- if (length(file_arg) > 0) {
    dirname(normalizePath(sub("^--file=", "", file_arg[[1]])))
  } else {
    getwd()
  }
  candidates <- c(
    file.path(getwd(), "r_subject_exclusions.R"),
    file.path(script_dir, "r_subject_exclusions.R")
  )
  helper_path <- candidates[file.exists(candidates)][1]
  if (is.na(helper_path)) {
    stop("Could not locate r_subject_exclusions.R. Run from repo root or place helper beside the script.")
  }
  source(helper_path, local = parent.frame())
}
source_exclusion_helpers()

REQUIRED_DIFF_COLS <- c(
  "diff_short_form_education",
  "diff_short_form_entertainment",
  "diff_long_form_education",
  "diff_long_form_entertainment"
)

# Parse `--key value` CLI arguments with strict typing for numeric knobs.
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  defaults <- list(
    input_csv = "data/tabular/homer3_betas_plus_combined_sfv_data_inner_join.csv",
    exclude_subjects_json = "data/tabular/excluded_subjects.json",
    out_main_csv = "data/results/retention_format_content_lmm_main_effects_r.csv",
    out_posthoc_csv = "data/results/retention_format_content_lmm_posthoc_pairwise_r.csv",
    alpha = 0.05,
    min_subjects = 6L,
    pbkrtest_limit = NA_real_,
    lmerTest_limit = NA_real_
  )
  if (length(args) == 0) return(defaults)
  if (length(args) %% 2 != 0) stop("Expected --key value argument pairs.")
  parsed <- defaults
  for (i in seq(1, length(args), by = 2)) {
    key <- args[[i]]
    val <- args[[i + 1]]
    if (!startsWith(key, "--")) stop(paste0("Invalid argument: ", key))
    key <- substring(key, 3)
    parsed[[key]] <- val
  }
  parsed$alpha <- as.numeric(parsed$alpha)
  parsed$min_subjects <- as.integer(parsed$min_subjects)
  parsed$pbkrtest_limit <- as.numeric(parsed$pbkrtest_limit)
  parsed$lmerTest_limit <- as.numeric(parsed$lmerTest_limit)
  parsed
}

# Normalize subject identifiers to integer IDs by extracting digit sequences.
# This keeps joins/uniqueness checks robust to formats like `0017` vs `sub_0017`.
normalize_subject_id <- function(x, column_name) {
  x_chr <- as.character(x)
  extracted <- stringr::str_extract(x_chr, "\\d+")
  if (any(is.na(extracted))) {
    bad <- unique(x_chr[is.na(extracted)])
    bad <- bad[!is.na(bad)]
    bad <- head(bad, 10)
    stop(
      paste0(
        "Failed to parse numeric IDs from column '", column_name, "'. ",
        "Examples of unparseable values: ", paste(bad, collapse = ", ")
      )
    )
  }
  as.integer(extracted)
}

# Fail fast if expected schema columns are missing.
assert_required_columns <- function(df, required_cols, file_label) {
  missing <- setdiff(required_cols, names(df))
  if (length(missing) > 0) {
    stop(
      paste0(
        "Missing required columns in ", file_label, ": ",
        paste(sort(missing), collapse = ", ")
      )
    )
  }
}

# Coerce specified columns to numeric and fail hard on non-numeric tokens.
# This avoids silent coercion that could alter inference.
coerce_numeric_strict <- function(df, cols) {
  out <- df
  for (col_name in cols) {
    raw <- out[[col_name]]
    raw_chr <- as.character(raw)
    suppressWarnings(num <- as.numeric(raw_chr))
    bad <- is.na(num) & !is.na(raw) & nzchar(trimws(raw_chr))
    if (any(bad)) {
      examples <- unique(raw_chr[bad])
      examples <- head(examples, 10)
      stop(
        paste0(
          "Column '", col_name, "' contains non-numeric values. ",
          "Examples: ", paste(examples, collapse = ", ")
        )
      )
    }
    out[[col_name]] <- num
  }
  out
}

# Load and validate the merged retention input table.
# Enforces one-row-per-subject integrity before any modeling.
load_retention_input <- function(input_csv, exclude_subjects_json) {
  df <- read_csv(input_csv, show_col_types = FALSE)
  assert_required_columns(df, c("subject_id", REQUIRED_DIFF_COLS), input_csv)

  df <- df %>%
    mutate(subject_id = normalize_subject_id(.data$subject_id, "subject_id"))
  df <- coerce_numeric_strict(df, REQUIRED_DIFF_COLS)

  dup <- df %>%
    count(subject_id, name = "n_rows") %>%
    filter(.data$n_rows > 1) %>%
    arrange(desc(.data$n_rows), .data$subject_id)
  if (nrow(dup) > 0) {
    examples <- dup %>% head(10)
    stop(
      paste0(
        "Duplicate subject_id values detected after normalization. ",
        "Expected exactly one row per subject. ",
        "Example duplicates (subject_id, n_rows): ",
        paste0(examples$subject_id, "=", examples$n_rows, collapse = ", "),
        ". Fix the upstream CSV before running retention analysis."
      )
    )
  }
  excluded <- apply_subject_exclusions(
    df = df,
    subject_col = "subject_id",
    exclude_json_path = exclude_subjects_json,
    context_label = "retention"
  )
  excluded$data
}

# Convert wide `diff_*` retention columns into long repeated-measures format.
# Mapping is explicit to avoid ambiguous category parsing.
reshape_to_long <- function(df) {
  df %>%
    select(subject_id, all_of(REQUIRED_DIFF_COLS)) %>%
    pivot_longer(cols = all_of(REQUIRED_DIFF_COLS), names_to = "diff_col", values_to = "retention_diff") %>%
    mutate(
      condition = case_when(
        diff_col == "diff_short_form_education" ~ "SF_Edu",
        diff_col == "diff_short_form_entertainment" ~ "SF_Ent",
        diff_col == "diff_long_form_entertainment" ~ "LF_Ent",
        diff_col == "diff_long_form_education" ~ "LF_Edu",
        TRUE ~ NA_character_
      ),
      length = case_when(
        diff_col %in% c("diff_short_form_education", "diff_short_form_entertainment") ~ "Short",
        diff_col %in% c("diff_long_form_education", "diff_long_form_entertainment") ~ "Long",
        TRUE ~ NA_character_
      ),
      content = case_when(
        diff_col %in% c("diff_short_form_education", "diff_long_form_education") ~ "Education",
        diff_col %in% c("diff_short_form_entertainment", "diff_long_form_entertainment") ~ "Entertainment",
        TRUE ~ NA_character_
      ),
      length_c = case_when(
        .data$length == "Short" ~ -0.5,
        .data$length == "Long" ~ +0.5,
        TRUE ~ NA_real_
      ),
      content_c = case_when(
        .data$content == "Entertainment" ~ -0.5,
        .data$content == "Education" ~ +0.5,
        TRUE ~ NA_real_
      )
    ) %>%
    filter(!is.na(.data$condition), !is.na(.data$length_c), !is.na(.data$content_c))
}

# Complete-case rule: a subject contributes only if all four conditions are present.
complete_case_subjects <- function(df_long) {
  non_missing <- df_long %>% filter(!is.na(.data$retention_diff))
  keep_ids <- non_missing %>%
    group_by(subject_id) %>%
    summarize(n_cond = n_distinct(condition), .groups = "drop") %>%
    filter(.data$n_cond == 4) %>%
    pull(subject_id)
  non_missing %>% filter(.data$subject_id %in% keep_ids)
}

fit_factorial_lmm <- function(df_cc) {
  # Random-intercept mixed model for repeated measures.
  # References: Laird & Ware (1982); Bates et al. (2015), see `CITATIONS.md`.
  # Inference: Satterthwaite-approximated df via lmerTest (Kuznetsova et al., 2017).
  lmerTest::lmer(
    retention_diff ~ length_c * content_c + (1 | subject_id),
    data = df_cc,
    REML = TRUE
  )
}

# Condition-coded model is only used to generate post-hoc pairwise contrasts.
fit_condition_lmm <- function(df_cc) {
  df_cc <- df_cc %>%
    mutate(condition = factor(.data$condition, levels = c("SF_Edu", "SF_Ent", "LF_Ent", "LF_Edu")))
  lmerTest::lmer(retention_diff ~ condition + (1 | subject_id), data = df_cc, REML = TRUE)
}

# Extract fixed-effect summary statistics and Wald CI for a named model term.
extract_fixed <- function(model, term) {
  coefs <- summary(model)$coefficients
  if (!(term %in% rownames(coefs))) stop(paste0("Expected term '", term, "' in model coefficients."))
  est <- as.numeric(coefs[term, "Estimate"])
  se <- as.numeric(coefs[term, "Std. Error"])
  df <- as.numeric(coefs[term, "df"])
  t <- as.numeric(coefs[term, "t value"])
  p <- as.numeric(coefs[term, "Pr(>|t|)"])
  ci <- suppressMessages(confint(model, parm = term, method = "Wald"))
  ci_low <- as.numeric(ci[1])
  ci_high <- as.numeric(ci[2])
  list(estimate = est, se = se, df = df, stat = t, p_unc = p, ci95_low = ci_low, ci95_high = ci_high)
}

main <- function() {
  args <- parse_args()

  # Optional overrides for emmeans df safeguards (large-N behavior).
  if (!is.na(args$pbkrtest_limit)) {
    emmeans::emm_options(pbkrtest.limit = args$pbkrtest_limit)
  }
  if (!is.na(args$lmerTest_limit)) {
    emmeans::emm_options(lmerTest.limit = args$lmerTest_limit)
  }

  # 1) Load + validate input, then reshape to repeated-measures long format.
  df <- load_retention_input(args$input_csv, args$exclude_subjects_json)
  df_long <- reshape_to_long(df)
  df_cc <- complete_case_subjects(df_long)

  n_subjects_raw <- dplyr::n_distinct(df$subject_id)
  n_subjects_cc <- dplyr::n_distinct(df_cc$subject_id)
  n_obs_cc <- nrow(df_cc)
  cat("[data] input subjects:", n_subjects_raw, "\n")
  cat("[data] complete-case subjects:", n_subjects_cc, "\n")
  cat("[data] complete-case observations:", n_obs_cc, "\n")

  if (n_subjects_cc < args$min_subjects) {
    stop(
      paste0(
        "Insufficient complete-case subjects after filtering: n_subjects=",
        n_subjects_cc,
        " < min_subjects=",
        args$min_subjects
      )
    )
  }

  # 2) Fit the primary factorial model and extract the three planned omnibus effects.
  model <- fit_factorial_lmm(df_cc)
  singular_fit <- lme4::isSingular(model, tol = 1e-4)
  eff_length <- extract_fixed(model, "length_c")
  eff_content <- extract_fixed(model, "content_c")
  eff_inter <- extract_fixed(model, "length_c:content_c")

  main_df <- tibble::tibble(
    effect = c("length", "content", "interaction"),
    estimate = c(eff_length$estimate, eff_content$estimate, eff_inter$estimate),
    se = c(eff_length$se, eff_content$se, eff_inter$se),
    df = c(eff_length$df, eff_content$df, eff_inter$df),
    t = c(eff_length$stat, eff_content$stat, eff_inter$stat),
    p_unc = c(eff_length$p_unc, eff_content$p_unc, eff_inter$p_unc),
    ci95_low = c(eff_length$ci95_low, eff_content$ci95_low, eff_inter$ci95_low),
    ci95_high = c(eff_length$ci95_high, eff_content$ci95_high, eff_inter$ci95_high)
  ) %>%
    mutate(
      # One Holm-adjusted family over the three pre-planned omnibus effects.
      # Column name `p_fdr` is retained for downstream compatibility.
      p_fdr = p.adjust(.data$p_unc, method = "holm"),
      n_subjects = n_subjects_cc,
      n_obs = n_obs_cc,
      singular_fit = singular_fit
    ) %>%
    select(
      effect, n_subjects, n_obs, singular_fit,
      estimate, se, df, t, p_unc, p_fdr, ci95_low, ci95_high
  )

  # 3) Gate post-hoc contrasts on Holm-significant interaction only.
  p_inter_q <- main_df %>%
    filter(.data$effect == "interaction") %>%
    pull(.data$p_fdr)
  do_posthoc <- is.finite(p_inter_q) && (p_inter_q < args$alpha)
  cat(
    "[results] interaction adjusted p=",
    signif(p_inter_q, 6),
    " (alpha=",
    args$alpha,
    "); posthoc=",
    ifelse(do_posthoc, "yes", "no"),
    "\n",
    sep = ""
  )

  posthoc_df <- tibble::tibble(
    condition_a = character(),
    condition_b = character(),
    stat_type = character(),
    mean_diff = numeric(),
    se = numeric(),
    df = numeric(),
    t = numeric(),
    p_unc = numeric()
  )

  if (do_posthoc) {
    model_cond <- fit_condition_lmm(df_cc)
    em <- emmeans::emmeans(model_cond, ~ condition)
    # EMM contrasts for condition-level comparisons.
    # Reference: Lenth (2016); Searle et al. (1980), see `CITATIONS.md`.
    pw <- emmeans::contrast(em, method = "pairwise", adjust = "none") %>% as.data.frame()
    stat_type <- if ("t.ratio" %in% names(pw)) "t" else if ("z.ratio" %in% names(pw)) "z" else NA_character_
    if (is.na(stat_type)) stop("Expected 't.ratio' or 'z.ratio' in emmeans pairwise output.")
    if (stat_type == "z") {
      pw[["t.ratio"]] <- pw[["z.ratio"]]
      cat("[warn] emmeans returned z.ratio (asymptotic). Recording as t.\n")
    }
    if (!("df" %in% names(pw))) pw[["df"]] <- NA_real_

    # `mean_diff` is defined as `condition_a - condition_b`.
    posthoc_df <- pw %>%
      mutate(
        stat_type = stat_type,
        condition_a = stringr::str_trim(stringr::str_split_fixed(.data$contrast, " - ", 2)[, 1]),
        condition_b = stringr::str_trim(stringr::str_split_fixed(.data$contrast, " - ", 2)[, 2])
      ) %>%
      select(condition_a, condition_b, stat_type, estimate, SE, df, t.ratio, p.value)
    names(posthoc_df) <- c("condition_a", "condition_b", "stat_type", "mean_diff", "se", "df", "t", "p_unc")
    posthoc_df <- posthoc_df %>% arrange(.data$condition_a, .data$condition_b)
  }

  # 4) Write outputs (main effects always; post-hoc may be empty if gate fails).
  dir.create(dirname(args$out_main_csv), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(args$out_posthoc_csv), recursive = TRUE, showWarnings = FALSE)
  write_csv(main_df, args$out_main_csv)
  write_csv(posthoc_df, args$out_posthoc_csv)

  cat("[write] main effects:", args$out_main_csv, "\n")
  cat("[write] posthoc:     ", args$out_posthoc_csv, ifelse(nrow(posthoc_df) == 0, " (empty)", ""), "\n", sep = "")
}

main()
