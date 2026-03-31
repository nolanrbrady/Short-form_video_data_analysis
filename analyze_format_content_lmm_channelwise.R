#!/usr/bin/env Rscript

# Channelwise within-subject inference for Format x Content effects on prefrontal activation.
#
# Inputs
#   - data/tabular/generated_data/homer3_betas_plus_combined_sfv_data_inner_join.csv
#     (one-row-per-subject merged file with `subject_id` plus Homer beta columns)
#
# Design
#   - Within-subjects 2x2 factorial:
#       Format:  Short vs Long
#       Content: Education vs Entertainment
#   - Condition mapping:
#       Cond01 = Short-Form Education
#       Cond02 = Short-Form Entertainment
#       Cond03 = Long-Form Entertainment
#       Cond04 = Long-Form Education
#
# Model (per channel x chromophore)
#   - Omnibus LMM: beta ~ format_c * content_c + age + (1|subject_id)
#     with numeric sum coding:
#       format_c  = -0.5 (Short), +0.5 (Long)
#       content_c = -0.5 (Entertainment), +0.5 (Education)
#
# Missingness / pruned channels (repo policy)
#   - Derived FIR-to-AUC beta tables encode pruned channels as NA (do not impute).
#   - Default: complete-case within channel/chrom (subject must have all 4 conditions present).
#   - Age is a required omnibus covariate and must be complete after subject exclusions.
#
# Multiple testing correction
#   - BH-FDR separately per chromophore and per effect across channels.
#
# Post-hoc (only if interaction q < alpha for that channel/chrom)
#   - All pairwise comparisons among the 4 condition combinations (6 contrasts),
#     uncorrected (adjust="none").
#
# Citations (see CITATIONS.md)
#   - Benjamini & Hochberg (1995): BH-FDR.
#   - Laird & Ware (1982): mixed-effects models.
#   - Kuznetsova et al. (2017): lmerTest integration with mixed-model inference tooling.
#   - Kenward & Roger (1997); Halekoh & Højsgaard (2014): Kenward-Roger df/p-value approximation.
#   - Lenth (2021): emmeans for contrasts.

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

source_lmm_convergence_helpers <- function() {
  args_all <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args_all, value = TRUE)
  script_dir <- if (length(file_arg) > 0) {
    dirname(normalizePath(sub("^--file=", "", file_arg[[1]])))
  } else {
    getwd()
  }
  candidates <- c(
    file.path(getwd(), "r_lmm_convergence_helpers.R"),
    file.path(script_dir, "r_lmm_convergence_helpers.R")
  )
  helper_path <- candidates[file.exists(candidates)][1]
  if (is.na(helper_path)) {
    stop("Could not locate r_lmm_convergence_helpers.R. Run from repo root or place helper beside the script.")
  }
  source(helper_path, local = parent.frame())
}
source_lmm_convergence_helpers()

# Fixed response-unit scaling used only inside the neural mixed-model fits.
#
# Rationale:
# the neural beta/AUC values in this project are on a very small numerical scale
# (roughly 1e-5), which can trigger optimizer/Hessian conditioning warnings in
# otherwise interpretable fits. Multiplying the Gaussian response by a fixed
# global constant is only a unit change, so the model is fit in numerically
# friendlier units and the reported estimates/CIs are back-transformed to the
# original beta space before writing outputs.
#
# Bates et al. (2015) motivate checking and addressing convergence/conditioning
# issues in mixed models rather than treating every returned fit as trustworthy.
NEURAL_LMM_RESPONSE_SCALE <- 1e6

# Output filtering toggle (main-effects CSVs only)
#
# If TRUE, the script writes only rows with FDR-significant p-values (p_fdr < threshold)
# to the main-effects CSV outputs. This is useful for quickly reviewing "hits", but should
# be used with care because it suppresses non-significant results in the saved tables.
FILTER_MAIN_EFFECTS_TO_FDR_SIG_ONLY <- FALSE
FILTER_MAIN_EFFECTS_FDR_THRESHOLD <- 0.05

# Optional environment-variable override (useful for automated validation runs).
# Accepts: "TRUE"/"FALSE" (case-insensitive). When unset/empty, the constants above apply.
env_filter <- Sys.getenv("FILTER_MAIN_EFFECTS_TO_FDR_SIG_ONLY", unset = "")
if (nzchar(env_filter)) {
  FILTER_MAIN_EFFECTS_TO_FDR_SIG_ONLY <- toupper(env_filter) %in% c("TRUE", "T", "1", "YES", "Y")
}

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  defaults <- list(
    input_csv = "data/tabular/generated_data/homer3_betas_plus_combined_sfv_data_inner_join.csv",
    exclude_subjects_json = "data/config/excluded_subjects.json",
    out_main_csv = "data/results/format_content_lmm_main_effects_r.csv",
    out_main_tidy_csv = "data/results/format_content_lmm_main_effects_tidy_r.csv",
    out_posthoc_csv = "data/results/format_content_lmm_posthoc_pairwise_r.csv",
    alpha = 0.05,
    min_subjects = 6,
    # Optional: enable/raise df adjustment limits for emmeans/lmerTest in large datasets.
    # See emmeans messages about pbkrtest.limit and lmerTest.limit. Defaults keep safeguards on.
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

normalize_subject_id <- function(x, column_name) {
  # Normalize to an integer subject ID by extracting digits.
  #
  # Examples accepted:
  #   - "sub_0001" -> 1
  #   - "0017" -> 17
  #
  # This function is intentionally strict: if any row cannot be mapped to a numeric ID,
  # we fail fast rather than performing a partial/unsafe merge.
  x_chr <- as.character(x)
  extracted <- str_extract(x_chr, "\\d+")
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

assert_covariate_complete <- function(df, covariate_col, context_label) {
  missing_rows <- which(is.na(df[[covariate_col]]))
  if (length(missing_rows) > 0) {
    subject_examples <- unique(df$subject_id[missing_rows])
    subject_examples <- head(subject_examples, 10)
    stop(
      paste0(
        "Required covariate '", covariate_col, "' contains missing values after subject exclusions in ",
        context_label, ". Example subject_id values: ",
        paste(subject_examples, collapse = ", "),
        ". Resolve the covariate upstream before running this analysis."
      )
    )
  }
}

load_merged_input <- function(input_csv) {
  # Preferred path for channelwise analysis:
  # consume pre-merged wide table so covariates and beta columns are already co-located.
  df <- read_csv(input_csv, show_col_types = FALSE)
  assert_required_columns(df, c("subject_id", "age"), input_csv)
  beta_cols <- names(df)[str_detect(names(df), "^S\\d+_D\\d+_Cond\\d{2}_(HbO|HbR)$")]
  if (length(beta_cols) == 0) {
    stop(
      paste0(
        "No beta columns matched pattern like 'S01_D01_Cond01_HbO' in merged input: ",
        input_csv
      )
    )
  }

  df <- df %>%
    mutate(subject_id = normalize_subject_id(.data$subject_id, "subject_id"))
  df <- coerce_numeric_strict(df, "age")

  dup <- df %>%
    count(subject_id, name = "n_rows") %>%
    filter(.data$n_rows > 1) %>%
    arrange(desc(.data$n_rows), .data$subject_id)
  if (nrow(dup) > 0) {
    examples <- dup %>% head(10)
    stop(
      paste0(
        "Duplicate subject_id values detected after normalization in merged input. ",
        "Expected exactly one row per subject. ",
        "Example duplicates (subject_id, n_rows): ",
        paste0(examples$subject_id, "=", examples$n_rows, collapse = ", "),
        ". Fix the upstream CSV before running the channelwise analysis."
      )
    )
  }

  df
}

reshape_to_long <- function(df_merged) {
  # Convert wide beta columns to long format for per-(channel × chrom) modeling.
  #
  # Each wide beta column is one (channel × condition × chrom) measurement; long format creates one row per:
  #   subject × channel × chrom × condition
  #
  # Derived FIR-to-AUC beta tables preserve only NA as pruned-channel missingness.
  beta_cols <- names(df_merged)[str_detect(names(df_merged), "^S\\d+_D\\d+_Cond\\d{2}_(HbO|HbR)$")]
  if (length(beta_cols) == 0) stop("No beta columns matched pattern like 'S01_D01_Cond01_HbO'.")

  df_merged %>%
    select(subject_id, age, all_of(beta_cols)) %>%
    pivot_longer(cols = all_of(beta_cols), names_to = "beta_col", values_to = "beta") %>%
    extract(
      col = "beta_col",
      into = c("channel", "cond", "chrom"),
      regex = "^(S\\d+_D\\d+)_Cond(\\d{2})_(HbO|HbR)$",
      remove = TRUE
    ) %>%
    {
      # Fail hard on any unexpected condition codes to avoid silently dropping columns.
      bad_conds <- setdiff(unique(.$cond), c("01", "02", "03", "04"))
      if (length(bad_conds) > 0) {
        stop(paste0("Unexpected condition codes in beta columns: ", paste(sort(bad_conds), collapse = ", ")))
      }
      .
    } %>%
    mutate(beta = suppressWarnings(as.numeric(beta))) %>%
    mutate(
      condition = case_when(
        cond == "01" ~ "SF_Edu",
        cond == "02" ~ "SF_Ent",
        cond == "03" ~ "LF_Ent",
        cond == "04" ~ "LF_Edu",
        TRUE ~ NA_character_
      ),
      format = case_when(
        cond %in% c("01", "02") ~ "Short",
        cond %in% c("03", "04") ~ "Long",
        TRUE ~ NA_character_
      ),
      content = case_when(
        cond %in% c("01", "04") ~ "Education",
        cond %in% c("02", "03") ~ "Entertainment",
        TRUE ~ NA_character_
      ),
      format_c = case_when(cond %in% c("01", "02") ~ -0.5, cond %in% c("03", "04") ~ +0.5, TRUE ~ NA_real_),
      content_c = case_when(cond %in% c("02", "03") ~ -0.5, cond %in% c("01", "04") ~ +0.5, TRUE ~ NA_real_)
    ) %>%
    filter(!is.na(condition), !is.na(chrom), !is.na(channel))
}

complete_case_subjects <- function(sub) {
  # Spec complete-case rule (within channel × chrom):
  # keep only subjects with all 4 conditions present (non-missing beta).
  non_missing <- sub %>% filter(!is.na(beta))
  keep_ids <- non_missing %>%
    group_by(subject_id) %>%
    summarize(n_cond = n_distinct(condition), .groups = "drop") %>%
    filter(n_cond == 4) %>%
    pull(subject_id)
  non_missing %>% filter(subject_id %in% keep_ids)
}

fit_factorial_lmm <- function(sub_complete) {
  # Per-channel × chrom LMM with random intercept for subject (within-subject design).
  # References: Laird & Ware (1982); Bates et al. (2015), see `CITATIONS.md`.
  # Age enters additively as a subject-level omnibus covariate.
  # Inference: fixed-effect tests are extracted with Kenward-Roger denominator df
  # via lmerTest + pbkrtest (Kenward & Roger, 1997; Halekoh & Højsgaard, 2014).
  sub_complete <- sub_complete %>% mutate(beta = .data$beta * NEURAL_LMM_RESPONSE_SCALE)
  lmerTest::lmer(beta ~ format_c * content_c + age + (1 | subject_id), data = sub_complete, REML = TRUE)
}

fit_condition_lmm <- function(sub_complete) {
  # Condition-coded model used only for post-hoc estimated marginal means / pairwise contrasts.
  sub_complete <- sub_complete %>%
    mutate(
      beta = .data$beta * NEURAL_LMM_RESPONSE_SCALE,
      condition = factor(condition, levels = c("SF_Edu", "SF_Ent", "LF_Ent", "LF_Edu"))
    )
  lmerTest::lmer(beta ~ condition + (1 | subject_id), data = sub_complete, REML = TRUE)
}

back_transform_neural_metric <- function(x) {
  as.numeric(x) / NEURAL_LMM_RESPONSE_SCALE
}

extract_fixed <- function(model, term, anova_kr) {
  # Extract fixed-effect results using:
  # - Estimate/SE from model coefficient table
  # - Kenward-Roger DenDF and p-value from Type-III ANOVA
  # For 1-df terms, report t as sign(beta) * sqrt(F).
  # - Wald 95% CI
  if (!(term %in% rownames(anova_kr))) {
    stop(paste0("Expected term '", term, "' in Kenward-Roger ANOVA table."))
  }

  coefs <- summary(model)$coefficients
  if (!(term %in% rownames(coefs))) stop(paste0("Expected term '", term, "' in model coefficients."))
  est <- as.numeric(coefs[term, "Estimate"])
  se <- as.numeric(coefs[term, "Std. Error"])
  num_df <- as.numeric(anova_kr[term, "NumDF"])
  if (!is.finite(num_df) || abs(num_df - 1) > 1e-8) {
    stop(paste0("Expected NumDF=1 for term '", term, "' but got ", num_df, "."))
  }
  df <- as.numeric(anova_kr[term, "DenDF"])
  f_val <- as.numeric(anova_kr[term, "F value"])
  if (!is.finite(f_val) || f_val < 0) {
    stop(paste0("Invalid Kenward-Roger F statistic for term '", term, "': ", f_val))
  }
  t <- sign(est) * sqrt(f_val)
  p <- as.numeric(anova_kr[term, "Pr(>F)"])

  ci <- suppressMessages(confint(model, parm = term, method = "Wald"))
  ci_low <- back_transform_neural_metric(ci[1])
  ci_high <- back_transform_neural_metric(ci[2])
  list(
    estimate = back_transform_neural_metric(est),
    se = back_transform_neural_metric(se),
    df = df,
    stat = t,
    p_unc = p,
    ci95_low = ci_low,
    ci95_high = ci_high
  )
}

record_warning_example <- function(examples, label, messages, max_examples = 10) {
  if (length(examples) >= max_examples) {
    return(examples)
  }

  msg <- if (length(messages) > 0) {
    paste(messages, collapse = " | ")
  } else {
    "non-convergence warning captured"
  }
  c(examples, paste0(label, " [", msg, "]"))
}

main <- function() {
  args <- parse_args()
  if (!requireNamespace("pbkrtest", quietly = TRUE)) {
    stop(
      paste0(
        "Kenward-Roger inference requires the 'pbkrtest' package, but it is not installed. ",
        "Install it before running this script."
      )
    )
  }

  # Optional: override emmeans' safeguards for df calculations in large datasets.
  # This can increase computation/memory usage substantially.
  if (!is.na(args$pbkrtest_limit)) {
    emmeans::emm_options(pbkrtest.limit = args$pbkrtest_limit)
  }
  if (!is.na(args$lmerTest_limit)) {
    emmeans::emm_options(lmerTest.limit = args$lmerTest_limit)
  }
  emmeans::emm_options(lmer.df = "kenward-roger")

  df_merged <- load_merged_input(args$input_csv)
  cat("[data] merged input file:", args$input_csv, "\n")
  excluded <- apply_subject_exclusions(
    df = df_merged,
    subject_col = "subject_id",
    exclude_json_path = args$exclude_subjects_json,
    context_label = "format_content_channelwise"
  )
  df_merged <- excluded$data
  assert_covariate_complete(df_merged, "age", "format_content_channelwise")
  df_long <- reshape_to_long(df_merged)

  channels <- sort(unique(df_long$channel))
  chroms <- c("HbO", "HbR")

  cat("[data] merged subjects:", length(unique(df_merged$subject_id)), "\n")
  cat("[data] channels detected:", length(channels), "\n")
  cat("[data] chromophores detected:", paste(sort(unique(df_long$chrom)), collapse = ", "), "\n")
  cat("[model] omnibus covariate adjustment: age\n")
  cat("[model] internal neural response scaling: beta * ", NEURAL_LMM_RESPONSE_SCALE, " for fitting; outputs are back-transformed to original beta units\n", sep = "")

  main_rows <- list()
  gated_count <- 0
  gated_examples <- c()
  nonconverged_count <- 0
  nonconverged_examples <- c()

  for (chrom_name in chroms) {
    for (channel_name in channels) {
      # Important: avoid name collisions with df_long$chrom / df_long$channel.
      # Use distinct loop variable names to ensure filtering truly subsets per channel × chrom.
      sub <- df_long %>% filter(.data$chrom == chrom_name, .data$channel == channel_name)
      sub_cc <- complete_case_subjects(sub)
      n_subj <- length(unique(sub_cc$subject_id))
      if (n_subj < args$min_subjects) {
        gated_count <- gated_count + 1
        if (length(gated_examples) < 10) {
          gated_examples <- c(
            gated_examples,
            paste0(chrom_name, " ", channel_name, " (n_subjects=", n_subj, " < min_subjects=", args$min_subjects, ")")
          )
        }
        next
      }
      n_obs <- nrow(sub_cc)

      fit_result <- capture_lmm_fit(function() fit_factorial_lmm(sub_cc))
      model <- fit_result$model
      if (!fit_result$converged) {
        nonconverged_count <- nonconverged_count + 1
        nonconverged_examples <- record_warning_example(
          examples = nonconverged_examples,
          label = paste0(chrom_name, " ", channel_name),
          messages = fit_result$nonconvergence_messages
        )
      }
      anova_kr <- anova(model, ddf = "Kenward-Roger", type = 3)
      singular_fit <- lme4::isSingular(model, tol = 1e-4)
      fmt <- extract_fixed(model, "format_c", anova_kr)
      cnt <- extract_fixed(model, "content_c", anova_kr)
      inter <- extract_fixed(model, "format_c:content_c", anova_kr)

      main_rows[[length(main_rows) + 1]] <- tibble(
        channel = channel_name,
        chrom = chrom_name,
        n_subjects = n_subj,
        n_obs = n_obs,
        converged = fit_result$converged,
        singular_fit = singular_fit,
        estimate_format = fmt$estimate,
        se_format = fmt$se,
        df_format = fmt$df,
        t_format = fmt$stat,
        p_format_unc = fmt$p_unc,
        ci95_format_low = fmt$ci95_low,
        ci95_format_high = fmt$ci95_high,
        estimate_content = cnt$estimate,
        se_content = cnt$se,
        df_content = cnt$df,
        t_content = cnt$stat,
        p_content_unc = cnt$p_unc,
        ci95_content_low = cnt$ci95_low,
        ci95_content_high = cnt$ci95_high,
        estimate_interaction = inter$estimate,
        se_interaction = inter$se,
        df_interaction = inter$df,
        t_interaction = inter$stat,
        p_interaction_unc = inter$p_unc,
        ci95_interaction_low = inter$ci95_low,
        ci95_interaction_high = inter$ci95_high
      )
    }
  }

  if (length(main_rows) == 0) {
    msg <- "No models were fit (min_subjects too high or missing data)."
    if (gated_count > 0) {
      msg <- paste0(
        msg,
        " All ", gated_count, " channel/chrom pairs evaluated were skipped due to min_subjects gating (min_subjects=",
        args$min_subjects, ")."
      )
      if (length(gated_examples) > 0) {
        msg <- paste0(msg, " Examples: ", paste(gated_examples, collapse = "; "))
      }
    }
    stop(msg)
  }

  main_df <- bind_rows(main_rows) %>%
    group_by(chrom) %>%
    mutate(
      # Multiple testing correction:
      # BH-FDR separately per chromophore and per effect, across channels.
      # Reference: Benjamini & Hochberg (1995), see `CITATIONS.md`.
      p_format_fdr = p.adjust(p_format_unc, method = "BH"),
      p_content_fdr = p.adjust(p_content_unc, method = "BH"),
      p_interaction_fdr = p.adjust(p_interaction_unc, method = "BH")
    ) %>%
    ungroup() %>%
    arrange(chrom, channel)

  if (gated_count > 0) {
    cat(
      "[warn] skipped ", gated_count,
      " channel/chrom models due to min_subjects gating (min_subjects=",
      args$min_subjects, "); showing up to 10 examples: ",
      paste(gated_examples, collapse = "; "),
      "\n",
      sep = ""
    )
  }
  if (nonconverged_count > 0) {
    cat(
      "[warn] ", nonconverged_count,
      " channel/chrom models emitted mixed-model non-convergence warnings; rows were retained with converged=FALSE. Examples: ",
      paste(nonconverged_examples, collapse = "; "),
      "\n",
      sep = ""
    )
  }

  sig_inter <- main_df %>% filter(.data$p_interaction_fdr < args$alpha)
  cat("[results] interaction FDR-significant (q<", args$alpha, "): ", nrow(sig_inter), " channel/chrom pairs\n", sep = "")

  posthoc_rows <- list()
  posthoc_nonconverged_count <- 0
  posthoc_nonconverged_examples <- c()
  if (nrow(sig_inter) > 0) {
    for (i in seq_len(nrow(sig_inter))) {
      chrom_name <- sig_inter$chrom[[i]]
      channel_name <- sig_inter$channel[[i]]
      sub <- df_long %>% filter(.data$chrom == chrom_name, .data$channel == channel_name)
      sub_cc <- complete_case_subjects(sub) %>% mutate(condition = factor(condition, levels = c("SF_Edu", "SF_Ent", "LF_Ent", "LF_Edu")))
      fit_result <- capture_lmm_fit(function() fit_condition_lmm(sub_cc))
      model_cond <- fit_result$model
      if (!fit_result$converged) {
        posthoc_nonconverged_count <- posthoc_nonconverged_count + 1
        posthoc_nonconverged_examples <- record_warning_example(
          examples = posthoc_nonconverged_examples,
          label = paste0(chrom_name, " ", channel_name),
          messages = fit_result$nonconvergence_messages
        )
      }
      em <- emmeans::emmeans(model_cond, ~ condition, lmer.df = "kenward-roger")
      # Use contrast(..., method="pairwise") instead of emmeans::pairs(...) for compatibility
      # across emmeans versions (pairs() is an S3 method and may not be exported from the namespace).
      #
      # Reference for EMMs/contrasts: Lenth (2016); Searle et al. (1980), see `CITATIONS.md`.
      pw <- emmeans::contrast(em, method = "pairwise", adjust = "none") %>% as.data.frame()
      pw <- pw %>%
        mutate(
          estimate = back_transform_neural_metric(.data$estimate),
          SE = back_transform_neural_metric(.data$SE),
          converged = fit_result$converged
        )
      # Depending on df availability, emmeans may report either t.ratio (finite df) or z.ratio (asymptotic).
      # Do not substitute z as t: if only z is available, write a clear sentinel string in `t`.
      stat_type <- if ("t.ratio" %in% names(pw)) "t" else if ("z.ratio" %in% names(pw)) "z" else NA_character_
      if (is.na(stat_type)) stop("Expected 't.ratio' or 'z.ratio' in emmeans pairwise contrast output.")
      if (stat_type == "z") {
        pw[["t_display"]] <- "t can't be calculated"
        cat(
          "[warn] emmeans returned z.ratio (asymptotic). Writing \"t can't be calculated\" in posthoc t column for ",
          chrom_name, " ", channel_name, ".\n",
          sep = ""
        )
      } else {
        pw[["t_display"]] <- as.character(pw[["t.ratio"]])
      }
      if (!("df" %in% names(pw))) pw[["df"]] <- NA_real_
      # pw columns include: contrast, estimate, SE, df, p.value and either t.ratio or z.ratio.
      pw <- pw %>%
        mutate(
          channel = channel_name,
          chrom = chrom_name,
          converged = .data$converged,
          stat_type = stat_type,
          condition_a = str_trim(str_split_fixed(.data$contrast, " - ", 2)[, 1]),
          condition_b = str_trim(str_split_fixed(.data$contrast, " - ", 2)[, 2])
        ) %>%
        select(channel, chrom, converged, condition_a, condition_b, stat_type, estimate, SE, df, t_display, p.value)

      names(pw) <- c("channel", "chrom", "converged", "condition_a", "condition_b", "stat_type", "mean_diff", "se", "df", "t", "p_unc")
      posthoc_rows[[length(posthoc_rows) + 1]] <- pw
    }
  }
  if (posthoc_nonconverged_count > 0) {
    cat(
      "[warn] ", posthoc_nonconverged_count,
      " posthoc channel/chrom fits emitted mixed-model non-convergence warnings; rows were retained with converged=FALSE. Examples: ",
      paste(posthoc_nonconverged_examples, collapse = "; "),
      "\n",
      sep = ""
    )
  }

  dir.create(dirname(args$out_main_csv), showWarnings = FALSE, recursive = TRUE)
  dir.create(dirname(args$out_main_tidy_csv), showWarnings = FALSE, recursive = TRUE)
  dir.create(dirname(args$out_posthoc_csv), showWarnings = FALSE, recursive = TRUE)
  main_df_out <- main_df
  if (isTRUE(FILTER_MAIN_EFFECTS_TO_FDR_SIG_ONLY)) {
    main_df_out <- main_df %>%
      filter(
        .data$p_format_fdr < FILTER_MAIN_EFFECTS_FDR_THRESHOLD |
          .data$p_content_fdr < FILTER_MAIN_EFFECTS_FDR_THRESHOLD |
          .data$p_interaction_fdr < FILTER_MAIN_EFFECTS_FDR_THRESHOLD
      )
    cat(
      "[write] filtering main effects (wide) to rows with any p_*_fdr < ",
      FILTER_MAIN_EFFECTS_FDR_THRESHOLD,
      "\n",
      sep = ""
    )
  }
  write_csv(main_df_out, args$out_main_csv)
  cat("[write] main effects:", args$out_main_csv, "\n")

  # Spec-compliant tidy main-effects output: one row per (channel x chrom x effect)
  main_tidy <- bind_rows(
    main_df %>%
      transmute(
        channel = .data$channel,
        chrom = .data$chrom,
        effect = "format",
        n_subjects = .data$n_subjects,
        n_obs = .data$n_obs,
        converged = .data$converged,
        singular_fit = .data$singular_fit,
        estimate = .data$estimate_format,
        ci95_low = .data$ci95_format_low,
        ci95_high = .data$ci95_format_high,
        p_unc = .data$p_format_unc,
        p_fdr = .data$p_format_fdr
      ),
    main_df %>%
      transmute(
        channel = .data$channel,
        chrom = .data$chrom,
        effect = "content",
        n_subjects = .data$n_subjects,
        n_obs = .data$n_obs,
        converged = .data$converged,
        singular_fit = .data$singular_fit,
        estimate = .data$estimate_content,
        ci95_low = .data$ci95_content_low,
        ci95_high = .data$ci95_content_high,
        p_unc = .data$p_content_unc,
        p_fdr = .data$p_content_fdr
      ),
    main_df %>%
      transmute(
        channel = .data$channel,
        chrom = .data$chrom,
        effect = "interaction",
        n_subjects = .data$n_subjects,
        n_obs = .data$n_obs,
        converged = .data$converged,
        singular_fit = .data$singular_fit,
        estimate = .data$estimate_interaction,
        ci95_low = .data$ci95_interaction_low,
        ci95_high = .data$ci95_interaction_high,
        p_unc = .data$p_interaction_unc,
        p_fdr = .data$p_interaction_fdr
      )
  ) %>%
    arrange(.data$p_unc, .data$chrom, .data$effect, .data$channel)

  main_tidy_out <- main_tidy
  if (isTRUE(FILTER_MAIN_EFFECTS_TO_FDR_SIG_ONLY)) {
    main_tidy_out <- main_tidy %>% filter(.data$p_fdr < FILTER_MAIN_EFFECTS_FDR_THRESHOLD)
    cat(
      "[write] filtering main effects (tidy/spec) to rows with p_fdr < ",
      FILTER_MAIN_EFFECTS_FDR_THRESHOLD,
      "\n",
      sep = ""
    )
  }
  write_csv(main_tidy_out, args$out_main_tidy_csv)
  cat("[write] main effects (tidy/spec):", args$out_main_tidy_csv, "\n")

  if (length(posthoc_rows) > 0) {
    posthoc_df <- bind_rows(posthoc_rows) %>% arrange(chrom, channel, condition_a, condition_b)
    write_csv(posthoc_df, args$out_posthoc_csv)
    cat("[write] posthoc:     ", args$out_posthoc_csv, "\n")
  } else {
    empty <- tibble(
      channel = character(),
      chrom = character(),
      converged = logical(),
      condition_a = character(),
      condition_b = character(),
      stat_type = character(),
      mean_diff = numeric(),
      se = numeric(),
      df = numeric(),
      t = character(),
      p_unc = numeric()
    )
    write_csv(empty, args$out_posthoc_csv)
    cat("[write] posthoc:     ", args$out_posthoc_csv, " (empty; no significant interactions)\n")
  }
}

main()
