#!/usr/bin/env Rscript

# Exploratory pooled-mean neural-behavior correlations for the SFV study.
#
# Scope
#   - Select the currently FDR-significant channel and ROI hits from the tidy
#     LMM outputs.
#   - Convert each selected neural target into subject-level pooled means for:
#       * short
#       * long
#       * education
#       * entertainment
#   - Convert engagement and retention into subject-level pooled means for the
#     same four pools.
#   - Correlate matched pools only (e.g., neural short vs recall short).
#
# Missingness policy
#   - Channel beta values of `0` and `NA` are treated as pruned/missing for
#     this workflow.
#   - ROI scores are built from arithmetic means across available non-missing
#     member channels.
#   - A pooled mean requires both constituent condition cells for that subject.
#   - No imputation is allowed.
#
# Interpretation
#   - This workflow is exploratory because it reuses the same dataset for
#     target selection and follow-up association testing.
#
# Citations (see CITATIONS.md)
#   - Kriegeskorte et al. (2009): selective-inference caution.
#   - Searle et al. (1980): equal-weight marginal-mean logic behind the pooled
#     main-effect means.
#   - Poldrack (2007): ROI averaging across pre-specified channels.
#   - Pearson (1896): Pearson product-moment correlation.
#   - Fisher (1921): Fisher-z confidence intervals for Pearson r.
#   - Benjamini & Hochberg (1995): BH-FDR.
#   - Bender & Lange (2001): multiplicity family definition for exploratory
#     association analyses.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(jsonlite)
})

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

CONDITION_MAP <- tibble::tribble(
  ~cond, ~condition_label, ~length, ~content,
  "01", "SF_Edu", "Short", "Education",
  "02", "SF_Ent", "Short", "Entertainment",
  "03", "LF_Ent", "Long", "Entertainment",
  "04", "LF_Edu", "Long", "Education"
)

BEHAVIOR_COLUMN_MAP <- tibble::tribble(
  ~behavior_domain, ~source_col, ~condition_label,
  "engagement", "sf_education_engagement", "SF_Edu",
  "engagement", "sf_entertainment_engagement", "SF_Ent",
  "engagement", "lf_entertainment_engagement", "LF_Ent",
  "engagement", "lf_education_engagement", "LF_Edu",
  "retention", "diff_short_form_education", "SF_Edu",
  "retention", "diff_short_form_entertainment", "SF_Ent",
  "retention", "diff_long_form_entertainment", "LF_Ent",
  "retention", "diff_long_form_education", "LF_Edu"
)

`%||%` <- function(a, b) if (!is.null(a)) a else b

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  defaults <- list(
    input_csv = "data/tabular/generated_data/homer3_betas_plus_combined_sfv_data_inner_join.csv",
    roi_json = "data/config/roi_definition.json",
    channel_results_csv = "data/results/format_content_lmm_main_effects_tidy_r.csv",
    roi_results_csv = "data/results/format_content_lmm_roi_main_effects_tidy_r.csv",
    exclude_subjects_json = "data/config/excluded_subjects.json",
    out_dir = "data/results/pooled_mean_correlations",
    alpha = 0.05,
    target_alpha = 0.05,
    min_subjects = 6L,
    figure_policy = "significant_only"
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
  parsed$alpha <- as.numeric(parsed$alpha)
  parsed$target_alpha <- as.numeric(parsed$target_alpha)
  parsed$min_subjects <- as.integer(parsed$min_subjects)
  parsed$figure_policy <- as.character(parsed$figure_policy)
  if (!(parsed$figure_policy %in% c("significant_only", "all_tested"))) {
    stop("figure_policy must be 'significant_only' or 'all_tested'.")
  }
  parsed
}

normalize_subject_id <- function(x, column_name) {
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
  x_chr <- toupper(str_trim(as.character(x)))
  m <- str_match(x_chr, "^S(\\d+)_D(\\d+)$")
  bad <- which(is.na(m[, 1]))
  if (length(bad) > 0) {
    bad_vals <- unique(x_chr[bad])
    bad_vals <- bad_vals[!is.na(bad_vals)]
    bad_vals <- head(bad_vals, 10)
    stop(
      paste0(
        "Failed to parse channel IDs in ", context_label,
        ". Expected format like S01_D01. Examples: ",
        paste(bad_vals, collapse = ", ")
      )
    )
  }
  paste0("S", sprintf("%02d", as.integer(m[, 2])), "_D", sprintf("%02d", as.integer(m[, 3])))
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
  allowed_missing_tokens <- c("", "NA", "NaN", "NAN")
  for (col_name in cols) {
    raw <- out[[col_name]]
    raw_chr <- trimws(as.character(raw))
    missing_token <- is.na(raw) | raw_chr %in% allowed_missing_tokens
    suppressWarnings(num <- as.numeric(raw_chr))
    num[is.nan(num)] <- NA_real_
    bad <- !missing_token & is.na(num)
    if (any(bad)) {
      examples <- unique(raw_chr[bad])
      examples <- head(examples, 10)
      stop(
        paste0(
          "Column '", col_name, "' contains non-numeric values. Examples: ",
          paste(examples, collapse = ", ")
        )
      )
    }
    out[[col_name]] <- num
  }
  out
}

clear_output_root <- function(out_dir) {
  if (!dir.exists(out_dir)) {
    return(FALSE)
  }
  deleted <- unlink(out_dir, recursive = TRUE, force = TRUE)
  if (deleted != 0 || dir.exists(out_dir)) {
    stop("Failed to clear previous pooled-mean correlation outputs.")
  }
  TRUE
}

derive_output_paths <- function(out_dir) {
  list(
    out_dir = out_dir,
    results_csv = file.path(out_dir, "pooled_mean_correlations_r.csv"),
    selected_targets_csv = file.path(out_dir, "selected_pooled_mean_targets_r.csv"),
    subject_pairs_csv = file.path(out_dir, "subject_level_pooled_mean_pairs_r.csv"),
    out_fig_dir = file.path(out_dir, "figures")
  )
}

load_roi_definition <- function(roi_json_path) {
  if (!file.exists(roi_json_path)) {
    stop(paste0("ROI definition file not found: ", roi_json_path))
  }
  roi_obj <- tryCatch(
    jsonlite::fromJSON(roi_json_path, simplifyVector = FALSE),
    error = function(e) {
      stop(paste0("Failed to parse ROI JSON at ", roi_json_path, ". Error: ", conditionMessage(e)))
    }
  )
  if (!is.list(roi_obj) || length(roi_obj) == 0 || is.null(names(roi_obj)) || any(!nzchar(names(roi_obj)))) {
    stop("ROI definition must be a non-empty JSON object mapping ROI names to channel arrays.")
  }

  roi_rows <- list()
  for (roi_name in names(roi_obj)) {
    channels_raw <- roi_obj[[roi_name]]
    channels <- if (is.list(channels_raw) && !is.data.frame(channels_raw)) unlist(channels_raw, use.names = FALSE) else channels_raw
    if (!is.character(channels) || length(channels) == 0) {
      stop(paste0("ROI '", roi_name, "' must map to a non-empty array of channel IDs."))
    }
    channels_norm <- normalize_channel_id(channels, paste0("roi_definition[", roi_name, "]"))
    if (anyDuplicated(channels_norm) > 0) {
      dup <- unique(channels_norm[duplicated(channels_norm)])
      stop(paste0("ROI '", roi_name, "' contains duplicate channel IDs after normalization: ", paste(dup, collapse = ", ")))
    }
    roi_rows[[length(roi_rows) + 1]] <- tibble::tibble(unit_id = roi_name, channel = channels_norm)
  }

  roi_map <- bind_rows(roi_rows)
  overlap <- roi_map %>% count(channel, name = "n_rois") %>% filter(.data$n_rois > 1)
  if (nrow(overlap) > 0) {
    offenders <- overlap %>% head(10)
    stop(paste0("Each channel must belong to exactly one ROI. Overlaps detected: ", paste0(offenders$channel, "=", offenders$n_rois, collapse = ", ")))
  }
  roi_map %>% arrange(.data$unit_id, .data$channel)
}

load_significant_targets <- function(path, analysis_level, unit_col, target_alpha) {
  df <- read_csv(path, show_col_types = FALSE)
  assert_required_columns(df, c(unit_col, "chrom", "effect", "estimate", "p_fdr"), path)
  unit_values <- if (analysis_level == "channel") normalize_channel_id(df[[unit_col]], path) else as.character(df[[unit_col]])
  out <- df %>%
    transmute(
      analysis_level = analysis_level,
      unit_id = unit_values,
      chrom = as.character(.data$chrom),
      selected_effect = as.character(.data$effect),
      selection_estimate = as.numeric(.data$estimate),
      selection_p_fdr = as.numeric(.data$p_fdr),
      selection_source = path
    ) %>%
    filter(.data$selected_effect %in% c("format", "content", "interaction")) %>%
    filter(is.finite(.data$selection_p_fdr) & .data$selection_p_fdr < target_alpha)
  if (nrow(out) == 0) return(out)
  out %>%
    group_by(.data$analysis_level, .data$unit_id, .data$chrom) %>%
    arrange(.data$selection_p_fdr, .by_group = TRUE) %>%
    slice(1) %>%
    ungroup()
}

load_merged_input <- function(input_csv, exclude_subjects_json) {
  required_behavior_cols <- BEHAVIOR_COLUMN_MAP$source_col
  df <- read_csv(input_csv, show_col_types = FALSE)
  assert_required_columns(df, c("subject_id", required_behavior_cols), input_csv)

  beta_cols <- names(df)[str_detect(names(df), "^S\\d+_D\\d+_Cond\\d{2}_(HbO|HbR)$")]
  if (length(beta_cols) == 0) {
    stop(paste0("No beta columns matched pattern like 'S01_D01_Cond01_HbO' in merged input: ", input_csv))
  }

  df <- df %>% mutate(subject_id = normalize_subject_id(.data$subject_id, "subject_id"))
  dup <- df %>% count(.data$subject_id, name = "n_rows") %>% filter(.data$n_rows > 1)
  if (nrow(dup) > 0) {
    offenders <- dup %>% head(10)
    stop(paste0("Duplicate subject_id values detected after normalization in merged input. Examples: ", paste0(offenders$subject_id, "=", offenders$n_rows, collapse = ", ")))
  }

  df <- coerce_numeric_strict(df, unique(c(required_behavior_cols, beta_cols)))
  excluded <- apply_subject_exclusions(
    df = df,
    subject_col = "subject_id",
    exclude_json_path = exclude_subjects_json,
    context_label = "pooled_mean_correlations"
  )
  list(data = excluded$data, beta_cols = beta_cols)
}

reshape_beta_long <- function(df_merged, beta_cols) {
  df_merged %>%
    select(subject_id, all_of(beta_cols)) %>%
    pivot_longer(cols = all_of(beta_cols), names_to = "beta_col", values_to = "beta") %>%
    extract(
      col = "beta_col",
      into = c("channel", "cond", "chrom"),
      regex = "^(S\\d+_D\\d+)_Cond(\\d{2})_(HbO|HbR)$",
      remove = TRUE
    ) %>%
    mutate(
      channel = normalize_channel_id(.data$channel, "beta columns"),
      beta = as.numeric(.data$beta),
      beta = ifelse(is.na(.data$beta) | .data$beta == 0, NA_real_, .data$beta)
    ) %>%
    left_join(CONDITION_MAP, by = "cond")
}

build_behavior_condition_rows <- function(df_merged) {
  bind_rows(lapply(seq_len(nrow(BEHAVIOR_COLUMN_MAP)), function(i) {
    spec <- BEHAVIOR_COLUMN_MAP[i, ]
    tibble::tibble(
      subject_id = df_merged$subject_id,
      behavior_domain = spec$behavior_domain[[1]],
      condition_label = spec$condition_label[[1]],
      behavior_value = df_merged[[spec$source_col[[1]]]]
    )
  })) %>%
    left_join(CONDITION_MAP %>% select(condition_label, length, content), by = "condition_label")
}

compute_pooled_means <- function(condition_rows, group_cols, value_col) {
  bind_rows(
    condition_rows %>%
      filter(.data$length == "Short") %>%
      group_by(across(all_of(group_cols))) %>%
      summarize(
        pool_name = "short",
        n_nonmissing = sum(!is.na(.data[[value_col]])),
        pooled_value = if (n_nonmissing == 2) mean(.data[[value_col]], na.rm = TRUE) else NA_real_,
        .groups = "drop"
      ),
    condition_rows %>%
      filter(.data$length == "Long") %>%
      group_by(across(all_of(group_cols))) %>%
      summarize(
        pool_name = "long",
        n_nonmissing = sum(!is.na(.data[[value_col]])),
        pooled_value = if (n_nonmissing == 2) mean(.data[[value_col]], na.rm = TRUE) else NA_real_,
        .groups = "drop"
      ),
    condition_rows %>%
      filter(.data$content == "Education") %>%
      group_by(across(all_of(group_cols))) %>%
      summarize(
        pool_name = "education",
        n_nonmissing = sum(!is.na(.data[[value_col]])),
        pooled_value = if (n_nonmissing == 2) mean(.data[[value_col]], na.rm = TRUE) else NA_real_,
        .groups = "drop"
      ),
    condition_rows %>%
      filter(.data$content == "Entertainment") %>%
      group_by(across(all_of(group_cols))) %>%
      summarize(
        pool_name = "entertainment",
        n_nonmissing = sum(!is.na(.data[[value_col]])),
        pooled_value = if (n_nonmissing == 2) mean(.data[[value_col]], na.rm = TRUE) else NA_real_,
        .groups = "drop"
      )
  )
}

build_channel_neural_means <- function(beta_long, selected_targets) {
  selected_channel_targets <- selected_targets %>% filter(.data$analysis_level == "channel") %>% distinct(.data$unit_id, .data$chrom)
  if (nrow(selected_channel_targets) == 0) {
    return(tibble::tibble(subject_id = integer(), analysis_level = character(), unit_id = character(), chrom = character(), pool_name = character(), neural_n_nonmissing = integer(), neural_value = double()))
  }
  beta_long %>%
    inner_join(selected_channel_targets, by = c("channel" = "unit_id", "chrom" = "chrom")) %>%
    transmute(subject_id, analysis_level = "channel", unit_id = .data$channel, chrom, condition_label, length, content, neural_value = .data$beta) %>%
    compute_pooled_means(group_cols = c("subject_id", "analysis_level", "unit_id", "chrom"), value_col = "neural_value") %>%
    rename(neural_n_nonmissing = n_nonmissing, neural_value = pooled_value)
}

build_roi_neural_means <- function(beta_long, roi_map, selected_targets) {
  selected_roi_targets <- selected_targets %>% filter(.data$analysis_level == "roi") %>% distinct(.data$unit_id, .data$chrom)
  if (nrow(selected_roi_targets) == 0) {
    return(tibble::tibble(subject_id = integer(), analysis_level = character(), unit_id = character(), chrom = character(), pool_name = character(), neural_n_nonmissing = integer(), neural_value = double()))
  }
  missing_rois <- setdiff(unique(selected_roi_targets$unit_id), unique(roi_map$unit_id))
  if (length(missing_rois) > 0) {
    stop(paste0("Selected ROI targets are missing from the ROI definition JSON: ", paste(sort(missing_rois), collapse = ", ")))
  }
  beta_long %>%
    inner_join(roi_map, by = "channel") %>%
    inner_join(selected_roi_targets, by = c("unit_id", "chrom")) %>%
    group_by(.data$subject_id, .data$unit_id, .data$chrom, .data$condition_label, .data$length, .data$content) %>%
    summarize(neural_value = if (all(is.na(.data$beta))) NA_real_ else mean(.data$beta, na.rm = TRUE), .groups = "drop") %>%
    mutate(analysis_level = "roi") %>%
    compute_pooled_means(group_cols = c("subject_id", "analysis_level", "unit_id", "chrom"), value_col = "neural_value") %>%
    rename(neural_n_nonmissing = n_nonmissing, neural_value = pooled_value)
}

build_behavior_means <- function(df_merged) {
  build_behavior_condition_rows(df_merged) %>%
    compute_pooled_means(group_cols = c("subject_id", "behavior_domain"), value_col = "behavior_value") %>%
    rename(behavior_n_nonmissing = n_nonmissing, behavior_value = pooled_value)
}

build_subject_level_pairs <- function(neural_means, behavior_means, selected_targets) {
  neural_means %>%
    inner_join(selected_targets, by = c("analysis_level", "unit_id", "chrom")) %>%
    inner_join(behavior_means, by = c("subject_id", "pool_name"), relationship = "many-to-many") %>%
    mutate(complete_pair = is.finite(.data$neural_value) & is.finite(.data$behavior_value)) %>%
    select(subject_id, analysis_level, unit_id, chrom, pool_name, behavior_domain, selected_effect, selection_estimate, selection_p_fdr, selection_source, neural_n_nonmissing, behavior_n_nonmissing, neural_value, behavior_value, complete_pair)
}

compute_pairwise_correlation <- function(sub_complete, alpha, min_subjects) {
  n_complete <- nrow(sub_complete)
  if (n_complete < min_subjects) {
    return(list(analysis_status = "skipped_min_subjects", skip_reason = paste0("n_complete<", min_subjects), n_complete = n_complete, association_estimate = NA_real_, r_squared = NA_real_, p_unc = NA_real_, ci95_low = NA_real_, ci95_high = NA_real_, slope = NA_real_, intercept = NA_real_))
  }
  if (dplyr::n_distinct(sub_complete$behavior_value) < 2) {
    return(list(analysis_status = "skipped_constant_input", skip_reason = "behavior_value_has_zero_variance", n_complete = n_complete, association_estimate = NA_real_, r_squared = NA_real_, p_unc = NA_real_, ci95_low = NA_real_, ci95_high = NA_real_, slope = NA_real_, intercept = NA_real_))
  }
  if (dplyr::n_distinct(sub_complete$neural_value) < 2) {
    return(list(analysis_status = "skipped_constant_input", skip_reason = "neural_value_has_zero_variance", n_complete = n_complete, association_estimate = NA_real_, r_squared = NA_real_, p_unc = NA_real_, ci95_low = NA_real_, ci95_high = NA_real_, slope = NA_real_, intercept = NA_real_))
  }

  cor_fit <- suppressWarnings(stats::cor.test(x = sub_complete$behavior_value, y = sub_complete$neural_value, method = "pearson", alternative = "two.sided", conf.level = 1 - alpha))
  lm_fit <- stats::lm(neural_value ~ behavior_value, data = sub_complete)
  lm_coef <- stats::coef(lm_fit)
  lm_summary <- summary(lm_fit)
  list(
    analysis_status = "tested",
    skip_reason = NA_character_,
    n_complete = n_complete,
    association_estimate = unname(cor_fit$estimate),
    r_squared = unname(lm_summary$r.squared),
    p_unc = cor_fit$p.value,
    ci95_low = if (!is.null(cor_fit$conf.int)) cor_fit$conf.int[[1]] else NA_real_,
    ci95_high = if (!is.null(cor_fit$conf.int)) cor_fit$conf.int[[2]] else NA_real_,
    slope = unname(lm_coef[["behavior_value"]]),
    intercept = unname(lm_coef[["(Intercept)"]])
  )
}

sanitize_slug <- function(x) {
  x %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("^_+|_+$", "") %>%
    str_to_lower()
}

plot_pairwise_correlation <- function(sub_complete, row, out_fig_dir) {
  annotation_lines <- c(
    paste0("n = ", row$n_complete),
    paste0("r = ", formatC(row$association_estimate, digits = 3, format = "f")),
    paste0("p = ", format.pval(row$p_unc, digits = 3, eps = 1e-4)),
    paste0("q = ", format.pval(row$p_fdr, digits = 3, eps = 1e-4))
  )
  file_stub <- sanitize_slug(paste(row$analysis_level, row$unit_id, row$chrom, row$behavior_domain, row$pool_name, sep = "_"))
  file_path <- file.path(out_fig_dir, paste0(file_stub, ".png"))
  p <- ggplot(sub_complete, aes(x = .data$behavior_value, y = .data$neural_value)) +
    geom_point(size = 2.2, alpha = 0.85, color = "#1b4d3e") +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#c04b2c", fill = "#f1c9b8") +
    labs(
      title = paste0(row$unit_id, " ", row$chrom, " ", row$pool_name, " mean vs ", row$behavior_domain),
      subtitle = paste0("Target: ", row$analysis_level, " | Selected effect: ", row$selected_effect),
      x = paste0(row$behavior_domain, " ", row$pool_name, " mean"),
      y = paste0("Neural ", row$pool_name, " mean")
    ) +
    theme_minimal(base_size = 11)
  x_anchor <- min(sub_complete$behavior_value, na.rm = TRUE)
  y_anchor <- max(sub_complete$neural_value, na.rm = TRUE)
  p <- p + annotate("label", x = x_anchor, y = y_anchor, label = paste(annotation_lines, collapse = "\n"), hjust = 0, vjust = 1, linewidth = 0.25, size = 3.2)
  ggsave(file_path, p, width = 7.2, height = 5.2, dpi = 160)
  file_path
}

main <- function() {
  args <- parse_args()
  outputs <- derive_output_paths(args$out_dir)

  selected_targets <- bind_rows(
    load_significant_targets(args$channel_results_csv, "channel", "channel", args$target_alpha),
    load_significant_targets(args$roi_results_csv, "roi", "roi", args$target_alpha)
  ) %>%
    arrange(.data$analysis_level, .data$unit_id, .data$chrom)
  if (nrow(selected_targets) == 0) {
    stop("No FDR-significant channel or ROI targets were found for the requested target_alpha.")
  }

  clear_output_root(outputs$out_dir)
  dir.create(outputs$out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(outputs$out_fig_dir, recursive = TRUE, showWarnings = FALSE)

  merged <- load_merged_input(args$input_csv, args$exclude_subjects_json)
  beta_long <- reshape_beta_long(merged$data, merged$beta_cols)
  roi_map <- load_roi_definition(args$roi_json)

  neural_means <- bind_rows(
    build_channel_neural_means(beta_long, selected_targets),
    build_roi_neural_means(beta_long, roi_map, selected_targets)
  )
  behavior_means <- build_behavior_means(merged$data)

  subject_pairs <- build_subject_level_pairs(neural_means, behavior_means, selected_targets) %>%
    arrange(.data$analysis_level, .data$unit_id, .data$chrom, .data$behavior_domain, .data$pool_name, .data$subject_id)

  result_groups <- subject_pairs %>%
    distinct(.data$analysis_level, .data$unit_id, .data$chrom, .data$pool_name, .data$behavior_domain, .data$selected_effect, .data$selection_estimate, .data$selection_p_fdr, .data$selection_source) %>%
    mutate(association_method = "pearson")

  results <- bind_rows(lapply(seq_len(nrow(result_groups)), function(i) {
    meta <- result_groups[i, ]
    sub_complete <- subject_pairs %>%
      filter(
        .data$analysis_level == meta$analysis_level[[1]],
        .data$unit_id == meta$unit_id[[1]],
        .data$chrom == meta$chrom[[1]],
        .data$pool_name == meta$pool_name[[1]],
        .data$behavior_domain == meta$behavior_domain[[1]]
      ) %>%
      filter(is.finite(.data$neural_value), is.finite(.data$behavior_value)) %>%
      select(subject_id, neural_value, behavior_value)
    stats <- compute_pairwise_correlation(sub_complete, alpha = args$alpha, min_subjects = args$min_subjects)
    tibble::tibble(
      analysis_level = meta$analysis_level[[1]],
      unit_id = meta$unit_id[[1]],
      chrom = meta$chrom[[1]],
      pool_name = meta$pool_name[[1]],
      behavior_domain = meta$behavior_domain[[1]],
      selected_effect = meta$selected_effect[[1]],
      selection_estimate = meta$selection_estimate[[1]],
      selection_p_fdr = meta$selection_p_fdr[[1]],
      selection_source = meta$selection_source[[1]],
      association_method = "pearson",
      analysis_status = stats$analysis_status,
      skip_reason = stats$skip_reason,
      n_complete = stats$n_complete,
      association_estimate = stats$association_estimate,
      r_squared = stats$r_squared,
      p_unc = stats$p_unc,
      ci95_low = stats$ci95_low,
      ci95_high = stats$ci95_high,
      slope = stats$slope,
      intercept = stats$intercept,
      plot_file = NA_character_
    )
  }))

  results <- results %>%
    mutate(
      family_id = paste0("behavior_domain=", .data$behavior_domain, " | pool_name=", .data$pool_name),
      family_adjust_method = "BH",
      family_n_tested = NA_integer_,
      p_fdr = NA_real_
    )

  tested_rows <- which(results$analysis_status == "tested" & is.finite(results$p_unc))
  if (length(tested_rows) > 0) {
    for (family_id in unique(results$family_id[tested_rows])) {
      idx <- which(results$family_id == family_id & results$analysis_status == "tested" & is.finite(results$p_unc))
      results$family_n_tested[idx] <- length(idx)
      results$p_fdr[idx] <- p.adjust(results$p_unc[idx], method = "BH")
    }
  }

  plot_idx <- if (args$figure_policy == "all_tested") {
    which(results$analysis_status == "tested")
  } else {
    which(results$analysis_status == "tested" & is.finite(results$p_fdr) & results$p_fdr < args$alpha)
  }
  if (length(plot_idx) > 0) {
    for (idx in plot_idx) {
      row <- results[idx, , drop = FALSE]
      sub_complete <- subject_pairs %>%
        filter(
          .data$analysis_level == row$analysis_level[[1]],
          .data$unit_id == row$unit_id[[1]],
          .data$chrom == row$chrom[[1]],
          .data$pool_name == row$pool_name[[1]],
          .data$behavior_domain == row$behavior_domain[[1]]
        ) %>%
        filter(is.finite(.data$neural_value), is.finite(.data$behavior_value)) %>%
        select(subject_id, neural_value, behavior_value)
      results$plot_file[[idx]] <- plot_pairwise_correlation(sub_complete, row, outputs$out_fig_dir)
    }
  }

  selected_targets_out <- selected_targets %>% arrange(.data$selection_p_fdr, .data$analysis_level, .data$unit_id, .data$chrom)
  results_out <- results %>% arrange(is.na(.data$p_fdr), .data$p_fdr, .data$p_unc, .data$behavior_domain, .data$pool_name, .data$analysis_level, .data$unit_id, .data$chrom)

  write_csv(selected_targets_out, outputs$selected_targets_csv, na = "NA")
  write_csv(subject_pairs, outputs$subject_pairs_csv, na = "NA")
  write_csv(results_out, outputs$results_csv, na = "NA")

  cat("[summary] selected targets:", nrow(selected_targets_out), "\n")
  cat("[summary] tested rows:", sum(results_out$analysis_status == "tested"), "\n")
  cat("[summary] plotted rows:", length(plot_idx), "\n")
  cat("[out] selected targets:", outputs$selected_targets_csv, "\n")
  cat("[out] subject-level pairs:", outputs$subject_pairs_csv, "\n")
  cat("[out] results CSV:", outputs$results_csv, "\n")
  cat("[out] figures dir:", outputs$out_fig_dir, "\n")
}

main()
