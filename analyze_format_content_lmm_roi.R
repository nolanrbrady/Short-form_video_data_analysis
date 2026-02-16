#!/usr/bin/env Rscript

# ROI-wise within-subject inference for Format x Content effects on prefrontal activation.
#
# Inputs
#   - data/tabular/homer3_betas_plus_combined_sfv_data_inner_join.csv
#     (one-row-per-subject merged file with `subject_id` plus Homer beta columns)
#   - data/config/roi_definition.json
#     (ROI name -> array of Homer channel IDs, e.g., "S01_D01")
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
# ROI beta construction
#   - Per subject x ROI x chrom x condition, beta is the arithmetic mean across
#     available (non-missing) channels in that ROI.
#   - This follows a standard ROI signal-summary approach in neuroimaging
#     (Poldrack, 2007; see CITATIONS.md).
#
# Model (per ROI x chromophore)
#   - LMM: beta ~ format_c * content_c + (1|subject_id)
#     with numeric sum coding:
#       format_c  = -0.5 (Short), +0.5 (Long)
#       content_c = -0.5 (Entertainment), +0.5 (Education)
#
# Missingness / pruned channels (repo policy)
#   - Treat BOTH 0 and NA in beta columns as pruned/missing (do not impute).
#   - Complete-case within ROI/chrom (subject must have all 4 conditions present).
#
# Multiple testing correction
#   - BH-FDR separately per chromophore and per effect across ROIs.
#
# Post-hoc (only if interaction q < alpha for that ROI/chrom)
#   - All pairwise comparisons among the 4 condition combinations (6 contrasts),
#     uncorrected (adjust="none").
#
# Citations (see CITATIONS.md)
#   - Poldrack (2007): ROI signal extraction framework.
#   - Benjamini & Hochberg (1995): BH-FDR.
#   - Laird & Ware (1982): mixed-effects models.
#   - Kuznetsova et al. (2017): lmerTest / Satterthwaite df.
#   - Lenth (2016): emmeans for contrasts.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(jsonlite)
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
    input_csv = "data/tabular/homer3_betas_plus_combined_sfv_data_inner_join.csv",
    roi_json = "data/config/roi_definition.json",
    exclude_subjects_json = "data/config/excluded_subjects.json",
    out_main_csv = "data/results/format_content_lmm_roi_main_effects_r.csv",
    out_main_tidy_csv = "data/results/format_content_lmm_roi_main_effects_tidy_r.csv",
    out_posthoc_csv = "data/results/format_content_lmm_roi_posthoc_pairwise_r.csv",
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

normalize_channel_id <- function(x, context_label) {
  # Normalize Homer channel IDs into canonical S##_## format for reliable matching.
  x_chr <- toupper(str_trim(as.character(x)))
  m <- str_match(x_chr, "^S(\\d+)_D(\\d+)$")
  bad <- which(is.na(m[, 1]))
  if (length(bad) > 0) {
    bad_vals <- unique(x_chr[bad])
    bad_vals <- bad_vals[!is.na(bad_vals)]
    bad_vals <- head(bad_vals, 10)
    stop(
      paste0(
        "Failed to parse channel IDs in ", context_label, ". ",
        "Expected format like S01_D01. Examples: ",
        paste(bad_vals, collapse = ", ")
      )
    )
  }
  paste0("S", sprintf("%02d", as.integer(m[, 2])), "_D", sprintf("%02d", as.integer(m[, 3])))
}

load_roi_definition <- function(roi_json_path) {
  if (!file.exists(roi_json_path)) {
    stop(paste0("ROI definition file not found: ", roi_json_path))
  }

  roi_obj <- tryCatch(
    jsonlite::fromJSON(roi_json_path, simplifyVector = FALSE),
    error = function(e) {
      stop(
        paste0(
          "Failed to parse ROI JSON at ", roi_json_path,
          ". Use strict JSON syntax (double quotes, comma-separated arrays). Error: ",
          conditionMessage(e)
        )
      )
    }
  )

  if (!is.list(roi_obj) || length(roi_obj) == 0 || is.null(names(roi_obj)) || any(!nzchar(names(roi_obj)))) {
    stop("ROI definition must be a non-empty JSON object mapping ROI names to channel arrays.")
  }

  roi_rows <- list()
  for (roi_name in names(roi_obj)) {
    channels_raw <- roi_obj[[roi_name]]
    channels <- if (is.list(channels_raw) && !is.data.frame(channels_raw)) {
      unlist(channels_raw, use.names = FALSE)
    } else {
      channels_raw
    }
    if (!is.character(channels) || length(channels) == 0) {
      stop(paste0("ROI '", roi_name, "' must map to a non-empty array of channel IDs."))
    }
    channels_norm <- normalize_channel_id(channels, paste0("roi_definition[", roi_name, "]"))
    if (anyDuplicated(channels_norm) > 0) {
      dup <- unique(channels_norm[duplicated(channels_norm)])
      stop(
        paste0(
          "ROI '", roi_name, "' contains duplicate channel IDs after normalization: ",
          paste(dup, collapse = ", ")
        )
      )
    }
    roi_rows[[length(roi_rows) + 1]] <- tibble(roi = roi_name, channel = channels_norm)
  }

  roi_map <- bind_rows(roi_rows)

  overlap <- roi_map %>%
    count(channel, name = "n_rois") %>%
    filter(.data$n_rois > 1)

  if (nrow(overlap) > 0) {
    offenders <- overlap %>% head(10)
    stop(
      paste0(
        "Each channel must belong to exactly one ROI. Overlaps detected: ",
        paste0(offenders$channel, "=", offenders$n_rois, collapse = ", ")
      )
    )
  }

  roi_map %>% arrange(roi, channel)
}

load_merged_input <- function(input_csv) {
  # Preferred path for ROI analysis:
  # consume pre-merged wide table so covariates and beta columns are already co-located.
  df <- read_csv(input_csv, show_col_types = FALSE)
  if (!("subject_id" %in% names(df))) {
    stop(paste0("Expected column 'subject_id' in merged input: ", input_csv))
  }
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
        ". Fix the upstream CSV before running ROI analysis."
      )
    )
  }

  df
}

reshape_to_long <- function(df_merged) {
  # Convert wide beta columns to long format for per-(ROI × chrom) modeling.
  #
  # Pruned-channel policy: treat beta == 0 as missing (do not impute).
  beta_cols <- names(df_merged)[str_detect(names(df_merged), "^S\\d+_D\\d+_Cond\\d{2}_(HbO|HbR)$")]
  if (length(beta_cols) == 0) stop("No beta columns matched pattern like 'S01_D01_Cond01_HbO'.")

  df_merged %>%
    select(subject_id, all_of(beta_cols)) %>%
    pivot_longer(cols = all_of(beta_cols), names_to = "beta_col", values_to = "beta") %>%
    extract(
      col = "beta_col",
      into = c("channel", "cond", "chrom"),
      regex = "^(S\\d+_D\\d+)_Cond(\\d{2})_(HbO|HbR)$",
      remove = TRUE
    ) %>%
    {
      bad_conds <- setdiff(unique(.$cond), c("01", "02", "03", "04"))
      if (length(bad_conds) > 0) {
        stop(paste0("Unexpected condition codes in beta columns: ", paste(sort(bad_conds), collapse = ", ")))
      }
      .
    } %>%
    mutate(
      channel = normalize_channel_id(.data$channel, "beta column names"),
      beta = suppressWarnings(as.numeric(.data$beta)),
      beta = na_if(.data$beta, 0)
    ) %>%
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

validate_roi_channels <- function(df_long, roi_map) {
  available_channels <- sort(unique(df_long$channel))
  roi_channels <- sort(unique(roi_map$channel))

  missing_in_data <- setdiff(roi_channels, available_channels)
  if (length(missing_in_data) > 0) {
    stop(
      paste0(
        "ROI definition contains channels not found in beta columns: ",
        paste(head(missing_in_data, 20), collapse = ", "),
        if (length(missing_in_data) > 20) " ..." else ""
      )
    )
  }

  unmapped <- setdiff(available_channels, roi_channels)
  if (length(unmapped) > 0) {
    cat(
      "[warn] ", length(unmapped),
      " channels in the input are not mapped to any ROI and will be excluded from ROI analysis.\n",
      sep = ""
    )
  }
}

aggregate_to_roi <- function(df_long, roi_map) {
  # Collapse channel-level betas to ROI-level summaries for inference.
  #
  # Signal summary choice (mean across available channels) is a standard ROI
  # extraction strategy in neuroimaging; see Poldrack (2007) in CITATIONS.md.
  df_long %>%
    inner_join(roi_map, by = "channel", relationship = "many-to-one") %>%
    group_by(subject_id, roi, chrom, condition, format, content, format_c, content_c) %>%
    summarize(
      n_channels_in_roi = n_distinct(channel),
      n_channels_nonmissing = sum(!is.na(beta)),
      beta = if (all(is.na(beta))) NA_real_ else mean(beta, na.rm = TRUE),
      .groups = "drop"
    )
}

complete_case_subjects <- function(sub) {
  # Spec complete-case rule (within ROI × chrom):
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
  # Per-ROI × chrom LMM with random intercept for subject (within-subject design).
  # References: Laird & Ware (1982); Bates et al. (2015), see `CITATIONS.md`.
  # Inference: lmerTest provides approximate p-values via Satterthwaite df (Kuznetsova et al., 2017).
  lmerTest::lmer(beta ~ format_c * content_c + (1 | subject_id), data = sub_complete, REML = TRUE)
}

fit_condition_lmm <- function(sub_complete) {
  # Condition-coded model used only for post-hoc estimated marginal means / pairwise contrasts.
  sub_complete <- sub_complete %>% mutate(condition = factor(condition, levels = c("SF_Edu", "SF_Ent", "LF_Ent", "LF_Edu")))
  lmerTest::lmer(beta ~ condition + (1 | subject_id), data = sub_complete, REML = TRUE)
}

extract_fixed <- function(model, term) {
  # Extract fixed-effect results from lmerTest summary:
  # - Estimate/SE/t/df/p from lmerTest (Satterthwaite df; see `CITATIONS.md`)
  # - Wald 95% CI
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

  if (!is.na(args$pbkrtest_limit)) {
    emmeans::emm_options(pbkrtest.limit = args$pbkrtest_limit)
  }
  if (!is.na(args$lmerTest_limit)) {
    emmeans::emm_options(lmerTest.limit = args$lmerTest_limit)
  }

  roi_map <- load_roi_definition(args$roi_json)

  df_merged <- load_merged_input(args$input_csv)
  cat("[data] merged input file:", args$input_csv, "\n")
  cat("[data] ROI definition:", args$roi_json, "\n")

  excluded <- apply_subject_exclusions(
    df = df_merged,
    subject_col = "subject_id",
    exclude_json_path = args$exclude_subjects_json,
    context_label = "format_content_roi"
  )
  df_merged <- excluded$data

  df_long <- reshape_to_long(df_merged)
  validate_roi_channels(df_long, roi_map)
  df_roi <- aggregate_to_roi(df_long, roi_map)

  rois <- sort(unique(roi_map$roi))
  chroms <- c("HbO", "HbR")

  cat("[data] merged subjects:", length(unique(df_merged$subject_id)), "\n")
  cat("[data] channels mapped to ROIs:", length(unique(roi_map$channel)), "\n")
  cat("[data] ROIs detected:", length(rois), "\n")
  cat("[data] chromophores detected:", paste(sort(unique(df_roi$chrom)), collapse = ", "), "\n")

  main_rows <- list()
  gated_count <- 0
  gated_examples <- c()

  for (chrom_name in chroms) {
    for (roi_name in rois) {
      sub <- df_roi %>% filter(.data$chrom == chrom_name, .data$roi == roi_name)
      sub_cc <- complete_case_subjects(sub)
      n_subj <- length(unique(sub_cc$subject_id))
      if (n_subj < args$min_subjects) {
        gated_count <- gated_count + 1
        if (length(gated_examples) < 10) {
          gated_examples <- c(
            gated_examples,
            paste0(chrom_name, " ", roi_name, " (n_subjects=", n_subj, " < min_subjects=", args$min_subjects, ")")
          )
        }
        next
      }
      n_obs <- nrow(sub_cc)

      model <- fit_factorial_lmm(sub_cc)
      singular_fit <- lme4::isSingular(model, tol = 1e-4)
      fmt <- extract_fixed(model, "format_c")
      cnt <- extract_fixed(model, "content_c")
      inter <- extract_fixed(model, "format_c:content_c")

      main_rows[[length(main_rows) + 1]] <- tibble(
        roi = roi_name,
        chrom = chrom_name,
        n_subjects = n_subj,
        n_obs = n_obs,
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
    msg <- "No ROI models were fit (min_subjects too high or missing data)."
    if (gated_count > 0) {
      msg <- paste0(
        msg,
        " All ", gated_count, " ROI/chrom pairs evaluated were skipped due to min_subjects gating (min_subjects=",
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
      # BH-FDR separately per chromophore and per effect, across ROIs.
      # Reference: Benjamini & Hochberg (1995), see `CITATIONS.md`.
      p_format_fdr = p.adjust(p_format_unc, method = "BH"),
      p_content_fdr = p.adjust(p_content_unc, method = "BH"),
      p_interaction_fdr = p.adjust(p_interaction_unc, method = "BH")
    ) %>%
    ungroup() %>%
    arrange(chrom, roi)

  if (gated_count > 0) {
    cat(
      "[warn] skipped ", gated_count,
      " ROI/chrom models due to min_subjects gating (min_subjects=",
      args$min_subjects, "); showing up to 10 examples: ",
      paste(gated_examples, collapse = "; "),
      "\n",
      sep = ""
    )
  }

  sig_inter <- main_df %>% filter(.data$p_interaction_fdr < args$alpha)
  cat("[results] interaction FDR-significant (q<", args$alpha, "): ", nrow(sig_inter), " ROI/chrom pairs\n", sep = "")

  posthoc_rows <- list()
  if (nrow(sig_inter) > 0) {
    for (i in seq_len(nrow(sig_inter))) {
      chrom_name <- sig_inter$chrom[[i]]
      roi_name <- sig_inter$roi[[i]]
      sub <- df_roi %>% filter(.data$chrom == chrom_name, .data$roi == roi_name)
      sub_cc <- complete_case_subjects(sub) %>% mutate(condition = factor(condition, levels = c("SF_Edu", "SF_Ent", "LF_Ent", "LF_Edu")))
      model_cond <- fit_condition_lmm(sub_cc)
      em <- emmeans::emmeans(model_cond, ~ condition)
      # Use contrast(..., method="pairwise") for compatibility across emmeans versions.
      # Reference for EMMs/contrasts: Lenth (2016); Searle et al. (1980), see `CITATIONS.md`.
      pw <- emmeans::contrast(em, method = "pairwise", adjust = "none") %>% as.data.frame()
      stat_type <- if ("t.ratio" %in% names(pw)) "t" else if ("z.ratio" %in% names(pw)) "z" else NA_character_
      if (is.na(stat_type)) stop("Expected 't.ratio' or 'z.ratio' in emmeans pairwise contrast output.")
      if (stat_type == "z") {
        pw[["t.ratio"]] <- pw[["z.ratio"]]
        cat(
          "[warn] emmeans returned z.ratio (asymptotic). Recording it as t.ratio in posthoc output for ",
          chrom_name, " ", roi_name, ".\n",
          sep = ""
        )
      }
      if (!("df" %in% names(pw))) pw[["df"]] <- NA_real_

      pw <- pw %>%
        mutate(
          roi = roi_name,
          chrom = chrom_name,
          stat_type = stat_type,
          condition_a = str_trim(str_split_fixed(.data$contrast, " - ", 2)[, 1]),
          condition_b = str_trim(str_split_fixed(.data$contrast, " - ", 2)[, 2])
        ) %>%
        select(roi, chrom, condition_a, condition_b, stat_type, estimate, SE, df, t.ratio, p.value)

      names(pw) <- c("roi", "chrom", "condition_a", "condition_b", "stat_type", "mean_diff", "se", "df", "t", "p_unc")
      posthoc_rows[[length(posthoc_rows) + 1]] <- pw
    }
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

  # Spec-compliant tidy main-effects output: one row per (ROI x chrom x effect)
  main_tidy <- bind_rows(
    main_df %>%
      transmute(
        roi = .data$roi,
        chrom = .data$chrom,
        effect = "format",
        n_subjects = .data$n_subjects,
        n_obs = .data$n_obs,
        singular_fit = .data$singular_fit,
        estimate = .data$estimate_format,
        ci95_low = .data$ci95_format_low,
        ci95_high = .data$ci95_format_high,
        p_unc = .data$p_format_unc,
        p_fdr = .data$p_format_fdr
      ),
    main_df %>%
      transmute(
        roi = .data$roi,
        chrom = .data$chrom,
        effect = "content",
        n_subjects = .data$n_subjects,
        n_obs = .data$n_obs,
        singular_fit = .data$singular_fit,
        estimate = .data$estimate_content,
        ci95_low = .data$ci95_content_low,
        ci95_high = .data$ci95_content_high,
        p_unc = .data$p_content_unc,
        p_fdr = .data$p_content_fdr
      ),
    main_df %>%
      transmute(
        roi = .data$roi,
        chrom = .data$chrom,
        effect = "interaction",
        n_subjects = .data$n_subjects,
        n_obs = .data$n_obs,
        singular_fit = .data$singular_fit,
        estimate = .data$estimate_interaction,
        ci95_low = .data$ci95_interaction_low,
        ci95_high = .data$ci95_interaction_high,
        p_unc = .data$p_interaction_unc,
        p_fdr = .data$p_interaction_fdr
      )
  ) %>%
    arrange(chrom, effect, roi)

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
    posthoc_df <- bind_rows(posthoc_rows) %>% arrange(chrom, roi, condition_a, condition_b)
    write_csv(posthoc_df, args$out_posthoc_csv)
    cat("[write] posthoc:     ", args$out_posthoc_csv, "\n")
  } else {
    empty <- tibble(
      roi = character(),
      chrom = character(),
      condition_a = character(),
      condition_b = character(),
      stat_type = character(),
      mean_diff = numeric(),
      se = numeric(),
      df = numeric(),
      t = numeric(),
      p_unc = numeric()
    )
    write_csv(empty, args$out_posthoc_csv)
    cat("[write] posthoc:     ", args$out_posthoc_csv, " (empty; no significant interactions)\n")
  }
}

main()
