#!/usr/bin/env Rscript

# Exploratory post-hoc ROI-mean x behavioral-mean correlations for the SFV study.
#
# Scope
#   - Restrict the correlation workflow to pooled-format behavioral rows only:
#       * engagement_long, engagement_short
#       * retention_long, retention_short
#   - Restrict the neural side to pooled ROI means only:
#       * R_DLPFC (HbR), L_DLPFC (HbO), M_DMPFC (HbO), L_DMPFC (HbO)
#
# Rationale
#   - ROI-level summary follows standard pre-specified ROI signal extraction
#     practice (Poldrack, 2007).
#   - Pearson is the primary metric, with Fisher-z confidence intervals
#     (Pearson, 1896; Fisher, 1921).
#   - Spearman is retained as a sensitivity metric for bounded outcomes
#     (Spearman, 1904).
#   - BH-FDR is applied within declared inferential families rather than across
#     unrelated questions (Benjamini & Hochberg, 1995; Bender & Lange, 2001).

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(jsonlite)
  library(ggplot2)
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
  ~cond, ~condition_label, ~format,
  "01", "SF_Edu", "Short",
  "02", "SF_Ent", "Short",
  "03", "LF_Ent", "Long",
  "04", "LF_Edu", "Long"
)

TARGET_ROI_SPECS <- tibble::tribble(
  ~neural_name, ~chrom,
  "R_DLPFC", "HbR",
  "L_DLPFC", "HbO",
  "M_DMPFC", "HbO",
  "L_DMPFC", "HbO"
)

`%||%` <- function(a, b) if (!is.null(a)) a else b

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  defaults <- list(
    input_csv = "data/tabular/homer3_betas_plus_combined_sfv_data_inner_join.csv",
    roi_json = "data/config/roi_definition.json",
    analysis_plan_json = "data/config/correlational_analysis_plan_roi_means.json",
    exclude_subjects_json = "data/config/excluded_subjects.json",
    out_csv = "data/results/correlational_relationships_roi_means/pairwise_correlations_r.csv",
    out_fig_dir = "data/results/correlational_relationships_roi_means/figures",
    alpha = 0.05,
    min_subjects = 6L
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
  parsed$min_subjects <- as.integer(parsed$min_subjects)
  parsed
}

normalize_json_string_array <- function(x, field_name) {
  values <- if (is.list(x) && !is.data.frame(x)) unlist(x, use.names = FALSE) else x
  if (!is.atomic(values) || length(values) == 0) {
    stop(paste0("Correlation analysis plan must define a non-empty '", field_name, "' array."))
  }
  values <- as.character(values)
  if (any(is.na(values) | trimws(values) == "")) {
    stop(paste0("Correlation analysis plan contains empty values in '", field_name, "'."))
  }
  values
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
        "Failed to parse numeric IDs from column '", column_name, "'. Examples: ",
        paste(bad, collapse = ", ")
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
        "Failed to parse channel IDs in ", context_label, ". Expected format like S01_D01. Examples: ",
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

as_single_string <- function(x, field_name) {
  value <- as.character(x %||% NA_character_)
  if (length(value) != 1 || is.na(value) || !nzchar(trimws(value))) {
    stop(paste0("Correlation analysis plan must define a non-empty string for '", field_name, "'."))
  }
  trimws(value)
}

as_single_flag <- function(x, field_name) {
  if (is.null(x) || length(x) != 1 || is.na(x)) {
    stop(paste0("Correlation analysis plan must define a boolean '", field_name, "'."))
  }
  as.logical(x)
}

validate_unique_names <- function(values, field_name) {
  if (anyDuplicated(values) > 0) {
    dup <- unique(values[duplicated(values)])
    stop(
      paste0(
        "Correlation analysis plan contains duplicate values in '", field_name, "': ",
        paste(dup, collapse = ", ")
      )
    )
  }
}

load_analysis_plan <- function(plan_json_path) {
  if (!file.exists(plan_json_path)) {
    stop(paste0("Correlation analysis plan file not found: ", plan_json_path))
  }
  plan_obj <- tryCatch(
    jsonlite::fromJSON(plan_json_path, simplifyVector = FALSE),
    error = function(e) {
      stop(
        paste0(
          "Failed to parse correlation analysis plan JSON at ", plan_json_path,
          ". Error: ", conditionMessage(e)
        )
      )
    }
  )
  if (!is.list(plan_obj)) {
    stop("Correlation analysis plan must be a JSON object.")
  }

  version <- as.numeric(plan_obj$version %||% NA_real_)
  if (!is.finite(version) || version != 6) {
    stop("Correlation analysis plan must define version = 6.")
  }

  behavior_runs_raw <- plan_obj$behavior_runs
  if (!is.list(behavior_runs_raw) || length(behavior_runs_raw) == 0) {
    stop("Correlation analysis plan must define a non-empty 'behavior_runs' array.")
  }
  behavior_runs <- lapply(seq_along(behavior_runs_raw), function(i) {
    run_obj <- behavior_runs_raw[[i]]
    if (!is.list(run_obj)) {
      stop("Each behavior_runs entry must be a JSON object.")
    }
    name <- as_single_string(run_obj$name, paste0("behavior_runs[", i, "].name"))
    run_type <- as_single_string(run_obj$run_type, paste0("behavior_runs[", i, "].run_type"))
    if (run_type != "pooled_format") {
      return(NULL)
    }
    short_columns <- normalize_json_string_array(
      run_obj$short_columns,
      paste0("behavior_runs[", i, "].short_columns")
    )
    long_columns <- normalize_json_string_array(
      run_obj$long_columns,
      paste0("behavior_runs[", i, "].long_columns")
    )
    if (length(short_columns) != 2 || length(long_columns) != 2) {
      stop(
        paste0(
          "Behavior run '", name,
          "' must define exactly two short_columns and two long_columns for the 2x2 design."
        )
      )
    }
    list(
      name = name,
      short_columns = short_columns,
      long_columns = long_columns,
      analysis_tier = "primary"
    )
  })
  behavior_runs <- Filter(Negate(is.null), behavior_runs)
  if (length(behavior_runs) == 0) {
    stop("Standalone ROI-mean script requires at least one pooled_format behavior run.")
  }
  behavior_names <- vapply(behavior_runs, `[[`, character(1), "name")
  validate_unique_names(behavior_names, "behavior_runs.name")

  association_methods <- normalize_json_string_array(plan_obj$association_methods, "association_methods")
  supported_methods <- c("pearson", "spearman")
  bad_methods <- setdiff(association_methods, supported_methods)
  if (length(bad_methods) > 0) {
    stop(
      paste0(
        "Unsupported association_methods values: ",
        paste(sort(bad_methods), collapse = ", "),
        ". Supported values: ", paste(supported_methods, collapse = ", ")
      )
    )
  }
  validate_unique_names(association_methods, "association_methods")

  mt <- plan_obj$multiple_testing
  if (!is.list(mt)) {
    stop("Correlation analysis plan must define a 'multiple_testing' object.")
  }
  adjust_method <- as_single_string(mt$adjust_method, "multiple_testing.adjust_method")
  family_grouping <- normalize_json_string_array(mt$family_grouping, "multiple_testing.family_grouping")
  validate_unique_names(family_grouping, "multiple_testing.family_grouping")

  figures_obj <- plan_obj$figures
  if (!is.list(figures_obj)) {
    stop("Correlation analysis plan must define a 'figures' object.")
  }
  figure_policy <- as_single_string(figures_obj$policy, "figures.policy")
  plot_methods <- normalize_json_string_array(
    figures_obj$plot_association_methods,
    "figures.plot_association_methods"
  )
  validate_unique_names(plot_methods, "figures.plot_association_methods")

  beta_missingness <- plan_obj$beta_missingness
  if (!is.list(beta_missingness)) {
    stop("Correlation analysis plan must define a 'beta_missingness' object.")
  }

  list(
    version = version,
    behavior_runs = behavior_runs,
    association_methods = association_methods,
    multiple_testing = list(
      adjust_method = adjust_method,
      family_grouping = family_grouping
    ),
    figures = list(
      policy = figure_policy,
      plot_association_methods = plot_methods
    ),
    beta_missingness = list(
      zero_is_missing = as_single_flag(beta_missingness$zero_is_missing, "beta_missingness.zero_is_missing"),
      audit_exact_zero = as_single_flag(beta_missingness$audit_exact_zero, "beta_missingness.audit_exact_zero")
    )
  )
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
          ". Error: ", conditionMessage(e)
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
    roi_rows[[length(roi_rows) + 1]] <- tibble::tibble(roi = roi_name, channel = channels_norm)
  }
  roi_map <- bind_rows(roi_rows)
  overlap <- roi_map %>% count(channel, name = "n_rois") %>% filter(.data$n_rois > 1)
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

collect_required_behavior_columns <- function(df_names, analysis_plan) {
  required_cols <- unique(unlist(lapply(analysis_plan$behavior_runs, function(run) {
    c(run$short_columns, run$long_columns)
  })))
  missing_cols <- setdiff(required_cols, df_names)
  if (length(missing_cols) > 0) {
    stop(
      paste0(
        "Merged input is missing required behavioral source columns: ",
        paste(sort(missing_cols), collapse = ", ")
      )
    )
  }
  required_cols
}

target_roi_channels <- function(roi_map) {
  roi_map %>% semi_join(TARGET_ROI_SPECS, by = c("roi" = "neural_name"))
}

collect_required_target_beta_columns <- function(roi_map) {
  roi_target_map <- target_roi_channels(roi_map) %>%
    inner_join(TARGET_ROI_SPECS, by = c("roi" = "neural_name"))
  unique(unlist(mapply(
    function(channel, chrom) {
      paste0(channel, "_Cond", CONDITION_MAP$cond, "_", chrom)
    },
    channel = roi_target_map$channel,
    chrom = roi_target_map$chrom,
    SIMPLIFY = FALSE
  ), use.names = FALSE))
}

load_merged_input <- function(input_csv, exclude_subjects_json, analysis_plan, roi_map) {
  df <- read_csv(input_csv, show_col_types = FALSE)
  if (!("subject_id" %in% names(df))) {
    stop(paste0("Expected column 'subject_id' in merged input: ", input_csv))
  }
  required_behavior_cols <- collect_required_behavior_columns(names(df), analysis_plan)
  required_target_beta_cols <- collect_required_target_beta_columns(roi_map)
  assert_required_columns(df, c(required_behavior_cols, required_target_beta_cols), input_csv)

  beta_cols <- names(df)[str_detect(names(df), "^S\\d+_D\\d+_Cond\\d{2}_(HbO|HbR)$")]
  if (length(beta_cols) == 0) {
    stop(paste0("No beta columns matched the expected Homer naming pattern in ", input_csv))
  }

  df <- df %>% mutate(subject_id = normalize_subject_id(.data$subject_id, "subject_id"))
  dup <- df %>% count(subject_id, name = "n_rows") %>% filter(.data$n_rows > 1)
  if (nrow(dup) > 0) {
    stop("Duplicate subject_id values detected after normalization in merged input.")
  }

  numeric_cols <- unique(c(required_behavior_cols, beta_cols))
  df <- coerce_numeric_strict(df, numeric_cols)
  excluded <- apply_subject_exclusions(
    df = df,
    subject_col = "subject_id",
    exclude_json_path = exclude_subjects_json,
    context_label = "roi_mean_correlation"
  )

  list(
    data = excluded$data,
    behavior_source_cols = required_behavior_cols,
    beta_cols = beta_cols
  )
}

audit_target_beta_missingness <- function(df_merged, roi_map, analysis_plan) {
  target_cols <- collect_required_target_beta_columns(roi_map)
  beta_values <- as.numeric(unlist(df_merged[, target_cols, drop = FALSE], use.names = FALSE))
  list(
    analyzed_beta_columns = length(target_cols),
    analyzed_beta_cells = length(beta_values),
    exact_zero_cells = if (analysis_plan$beta_missingness$audit_exact_zero) sum(beta_values == 0, na.rm = TRUE) else NA_integer_,
    explicit_na_cells = sum(is.na(beta_values))
  )
}

reshape_beta_long <- function(df_merged, beta_cols, analysis_plan) {
  df_merged %>%
    select(subject_id, all_of(beta_cols)) %>%
    pivot_longer(cols = all_of(beta_cols), names_to = "beta_col", values_to = "beta_raw") %>%
    extract(
      col = "beta_col",
      into = c("channel", "cond", "chrom"),
      regex = "^(S\\d+_D\\d+)_Cond(\\d{2})_(HbO|HbR)$",
      remove = TRUE
    ) %>%
    mutate(channel = normalize_channel_id(.data$channel, "beta columns")) %>%
    left_join(CONDITION_MAP, by = "cond") %>%
    mutate(
      beta = if (analysis_plan$beta_missingness$zero_is_missing) {
        dplyr::if_else(.data$beta_raw == 0, NA_real_, .data$beta_raw)
      } else {
        .data$beta_raw
      }
    ) %>%
    filter(!is.na(.data$condition_label))
}

build_roi_condition_targets <- function(beta_long, roi_map) {
  roi_map_target <- target_roi_channels(roi_map)
  missing_roi <- setdiff(TARGET_ROI_SPECS$neural_name, roi_map_target$roi)
  if (length(missing_roi) > 0) {
    stop(
      paste0(
        "Requested ROI(s) missing from ROI definition JSON: ",
        paste(sort(unique(missing_roi)), collapse = ", ")
      )
    )
  }
  available_channels <- unique(beta_long$channel)
  absent_channels <- setdiff(roi_map_target$channel, available_channels)
  if (length(absent_channels) > 0) {
    stop(
      paste0(
        "ROI definition references channels absent from merged input: ",
        paste(sort(unique(absent_channels)), collapse = ", ")
      )
    )
  }

  beta_long %>%
    inner_join(roi_map_target, by = "channel") %>%
    rename(neural_name = roi) %>%
    inner_join(TARGET_ROI_SPECS, by = c("neural_name", "chrom")) %>%
    group_by(.data$subject_id, .data$neural_name, .data$chrom, .data$cond, .data$condition_label, .data$format) %>%
    summarize(
      neural_value = if (all(is.na(.data$beta))) NA_real_ else mean(.data$beta, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    transmute(
      subject_id,
      neural_level = "roi",
      neural_name,
      chrom,
      cond,
      condition_label,
      format,
      neural_value
    )
}

build_neural_pooled_values <- function(condition_targets) {
  condition_targets %>%
    group_by(.data$subject_id, .data$neural_level, .data$neural_name, .data$chrom) %>%
    summarize(
      short_mean = if (sum(.data$format == "Short" & !is.na(.data$neural_value)) == 2) {
        mean(.data$neural_value[.data$format == "Short"], na.rm = TRUE)
      } else {
        NA_real_
      },
      long_mean = if (sum(.data$format == "Long" & !is.na(.data$neural_value)) == 2) {
        mean(.data$neural_value[.data$format == "Long"], na.rm = TRUE)
      } else {
        NA_real_
      },
      n_nonmissing_short_conditions = sum(.data$format == "Short" & !is.na(.data$neural_value)),
      n_nonmissing_long_conditions = sum(.data$format == "Long" & !is.na(.data$neural_value)),
      .groups = "drop"
    ) %>%
    tidyr::crossing(
      tibble::tibble(
        format_pool = c("long", "short"),
        neural_condition_label = c("Long", "Short")
      )
    ) %>%
    transmute(
      subject_id,
      neural_level,
      neural_name,
      chrom,
      format_pool,
      neural_condition_label,
      neural_effect = "pooled_format_mean",
      neural_value = ifelse(.data$format_pool == "long", .data$long_mean, .data$short_mean),
      neural_short_mean = .data$short_mean,
      neural_long_mean = .data$long_mean,
      n_nonmissing_conditions = ifelse(.data$format_pool == "long", .data$n_nonmissing_long_conditions, .data$n_nonmissing_short_conditions),
      neural_run_type = "pooled_format"
    )
}

build_behavior_analysis_rows <- function(df_merged, analysis_plan) {
  bind_rows(lapply(analysis_plan$behavior_runs, function(run) {
    short_mat <- as.matrix(df_merged[, run$short_columns, drop = FALSE])
    long_mat <- as.matrix(df_merged[, run$long_columns, drop = FALSE])
    short_complete <- rowSums(is.na(short_mat)) == 0
    long_complete <- rowSums(is.na(long_mat)) == 0
    short_mean <- ifelse(short_complete, rowMeans(short_mat), NA_real_)
    long_mean <- ifelse(long_complete, rowMeans(long_mat), NA_real_)

    base_rows <- tibble::tibble(
      subject_id = df_merged$subject_id,
      behavior_run = run$name,
      behavior_run_type = "pooled_format",
      behavior_effect = "pooled_format_mean",
      behavior_condition_code = NA_character_,
      behavior_condition_label = NA_character_,
      behavior_short_mean = short_mean,
      behavior_long_mean = long_mean,
      analysis_tier = run$analysis_tier
    )

    bind_rows(
      base_rows %>%
        transmute(
          subject_id,
          behavior_run,
          format_pool = "long",
          behavior_run_type,
          behavior_effect,
          behavior_condition_code,
          behavior_condition_label,
          behavior_short_mean,
          behavior_long_mean,
          behavior_value = behavior_long_mean,
          n_nonmissing_conditions = rowSums(!is.na(long_mat)),
          analysis_tier
        ),
      base_rows %>%
        transmute(
          subject_id,
          behavior_run,
          format_pool = "short",
          behavior_run_type,
          behavior_effect,
          behavior_condition_code,
          behavior_condition_label,
          behavior_short_mean,
          behavior_long_mean,
          behavior_value = behavior_short_mean,
          n_nonmissing_conditions = rowSums(!is.na(short_mat)),
          analysis_tier
        )
    )
  }))
}

make_pair_key <- function(behavior_run, format_pool, association_method, neural_name, chrom) {
  paste(behavior_run, format_pool, association_method, neural_name, chrom, sep = "||")
}

make_family_id <- function(df, family_grouping) {
  missing_grouping <- setdiff(family_grouping, names(df))
  if (length(missing_grouping) > 0) {
    stop(
      paste0(
        "Configured family grouping columns are missing from the correlation results: ",
        paste(missing_grouping, collapse = ", ")
      )
    )
  }
  apply(df[, family_grouping, drop = FALSE], 1, function(row) {
    paste(family_grouping, row, sep = "=", collapse = " | ")
  })
}

sanitize_slug <- function(x) {
  x %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("^_+|_+$", "") %>%
    str_to_lower()
}

clear_output_root <- function(out_csv, out_fig_dir) {
  output_root <- dirname(out_csv)
  output_root_norm <- normalizePath(output_root, winslash = "/", mustWork = FALSE)
  out_fig_norm <- normalizePath(out_fig_dir, winslash = "/", mustWork = FALSE)
  if (!(identical(out_fig_norm, output_root_norm) || startsWith(out_fig_norm, paste0(output_root_norm, "/")))) {
    stop("For safety, out_fig_dir must be inside dirname(out_csv) for ROI-mean output cleanup.")
  }
  if (!dir.exists(output_root)) {
    return(FALSE)
  }
  deleted <- unlink(output_root, recursive = TRUE, force = TRUE)
  if (deleted != 0 || dir.exists(output_root)) {
    stop("Failed to clear the previous ROI-mean output directory.")
  }
  TRUE
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

compute_association <- function(sub_complete, alpha, min_subjects, association_method) {
  n_complete <- nrow(sub_complete)
  if (n_complete < min_subjects) {
    return(list(
      status = "skipped_min_subjects",
      skip_reason = paste0("n_complete<", min_subjects),
      n_complete = n_complete,
      association_estimate = NA_real_,
      r_squared = NA_real_,
      p_unc = NA_real_,
      ci95_low = NA_real_,
      ci95_high = NA_real_,
      slope = NA_real_,
      intercept = NA_real_
    ))
  }
  if (dplyr::n_distinct(sub_complete$behavior_value) < 2) {
    return(list(
      status = "skipped_constant_input",
      skip_reason = "behavior_value_has_zero_variance",
      n_complete = n_complete,
      association_estimate = NA_real_,
      r_squared = NA_real_,
      p_unc = NA_real_,
      ci95_low = NA_real_,
      ci95_high = NA_real_,
      slope = NA_real_,
      intercept = NA_real_
    ))
  }
  if (dplyr::n_distinct(sub_complete$neural_value) < 2) {
    return(list(
      status = "skipped_constant_input",
      skip_reason = "neural_value_has_zero_variance",
      n_complete = n_complete,
      association_estimate = NA_real_,
      r_squared = NA_real_,
      p_unc = NA_real_,
      ci95_low = NA_real_,
      ci95_high = NA_real_,
      slope = NA_real_,
      intercept = NA_real_
    ))
  }

  if (association_method == "pearson") {
    cor_fit <- suppressWarnings(stats::cor.test(
      x = sub_complete$behavior_value,
      y = sub_complete$neural_value,
      method = "pearson",
      alternative = "two.sided",
      conf.level = 1 - alpha
    ))
    lm_fit <- stats::lm(neural_value ~ behavior_value, data = sub_complete)
    lm_coef <- stats::coef(lm_fit)
    lm_summary <- summary(lm_fit)
    return(list(
      status = "tested",
      skip_reason = NA_character_,
      n_complete = n_complete,
      association_estimate = unname(cor_fit$estimate),
      r_squared = unname(lm_summary$r.squared),
      p_unc = cor_fit$p.value,
      ci95_low = if (!is.null(cor_fit$conf.int)) cor_fit$conf.int[[1]] else NA_real_,
      ci95_high = if (!is.null(cor_fit$conf.int)) cor_fit$conf.int[[2]] else NA_real_,
      slope = unname(lm_coef[["behavior_value"]]),
      intercept = unname(lm_coef[["(Intercept)"]])
    ))
  }

  if (association_method == "spearman") {
    cor_fit <- suppressWarnings(stats::cor.test(
      x = sub_complete$behavior_value,
      y = sub_complete$neural_value,
      method = "spearman",
      exact = FALSE,
      alternative = "two.sided"
    ))
    return(list(
      status = "tested",
      skip_reason = NA_character_,
      n_complete = n_complete,
      association_estimate = unname(cor_fit$estimate),
      r_squared = NA_real_,
      p_unc = cor_fit$p.value,
      ci95_low = NA_real_,
      ci95_high = NA_real_,
      slope = NA_real_,
      intercept = NA_real_
    ))
  }

  stop(paste0("Unhandled association method: ", association_method))
}

plot_association <- function(sub_complete, row, out_fig_dir) {
  estimate_label <- if (row$association_method[[1]] == "spearman") "rho" else "r"
  annotation_lines <- c(
    paste0("n = ", row$n_complete),
    paste0(estimate_label, " = ", formatC(row$association_estimate, digits = 3, format = "f")),
    paste0("p = ", format.pval(row$p_unc, digits = 3, eps = 1e-4)),
    paste0("q = ", format.pval(row$p_fdr, digits = 3, eps = 1e-4))
  )
  file_stub <- sanitize_slug(paste(
    row$behavior_run,
    row$format_pool,
    row$association_method,
    row$neural_name,
    row$chrom,
    sep = "_"
  ))
  file_path <- file.path(out_fig_dir, paste0(file_stub, ".png"))
  x_anchor <- min(sub_complete$behavior_value, na.rm = TRUE)
  y_anchor <- max(sub_complete$neural_value, na.rm = TRUE)

  p <- ggplot(sub_complete, aes(x = .data$behavior_value, y = .data$neural_value)) +
    geom_point(size = 2.2, alpha = 0.85, color = "#1b4d3e")

  if (row$association_method[[1]] == "pearson") {
    p <- p + geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#c04b2c", fill = "#f1c9b8")
  } else {
    # Use a local smoother for Spearman plots so the visual summary matches the
    # monotonic, rank-based intent more closely than a forced linear fit.
    p <- p + geom_smooth(
      method = "loess",
      formula = y ~ x,
      se = TRUE,
      color = "#c04b2c",
      fill = "#f1c9b8",
      linetype = "dashed"
    )
  }

  p <- p +
    annotate(
      "label",
      x = x_anchor,
      y = y_anchor,
      label = paste(annotation_lines, collapse = "\n"),
      hjust = 0,
      vjust = 1,
      linewidth = 0.25,
      size = 3.2
    ) +
    labs(
      title = paste0(row$neural_name, " ", row$chrom, " ", row$format_pool, " mean vs ", row$behavior_run),
      subtitle = paste0(
        "Target: ROI | Metric: ", row$association_method, " | Tier: ", row$analysis_tier,
        if (row$association_method[[1]] == "spearman") " | LOESS smoother shown for rank-based trend visualization" else ""
      ),
      x = paste0(row$behavior_run, " ", row$format_pool, " mean"),
      y = paste0("Neural ", row$format_pool, " mean")
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )

  suppressMessages(ggplot2::ggsave(filename = file_path, plot = p, width = 7, height = 5, dpi = 300))
  file_path
}

main <- function() {
  args <- parse_args()
  cleared_output_root <- clear_output_root(args$out_csv, args$out_fig_dir)
  dir.create(dirname(args$out_csv), recursive = TRUE, showWarnings = FALSE)
  dir.create(args$out_fig_dir, recursive = TRUE, showWarnings = FALSE)

  roi_map <- load_roi_definition(args$roi_json)
  analysis_plan <- load_analysis_plan(args$analysis_plan_json)
  loaded <- load_merged_input(args$input_csv, args$exclude_subjects_json, analysis_plan, roi_map)
  df_merged <- loaded$data
  beta_cols <- loaded$beta_cols

  audit <- audit_target_beta_missingness(df_merged, roi_map, analysis_plan)
  beta_long <- reshape_beta_long(df_merged, beta_cols, analysis_plan)
  roi_targets <- build_roi_condition_targets(beta_long, roi_map)
  neural_analysis_targets <- build_neural_pooled_values(roi_targets) %>% filter(.data$neural_level == "roi")
  behavior_analysis_rows <- build_behavior_analysis_rows(df_merged, analysis_plan)

  pair_grid <- behavior_analysis_rows %>%
    inner_join(
      neural_analysis_targets,
      by = c("subject_id", "format_pool"),
      relationship = "many-to-many"
    )

  split_pairs <- pair_grid %>%
    tidyr::crossing(
      tibble::tibble(
        association_method = analysis_plan$association_methods,
        association_method_tier = ifelse(
          analysis_plan$association_methods == "pearson",
          "primary_metric",
          "sensitivity_metric"
        )
      )
    ) %>%
    group_by(
      .data$behavior_run,
      .data$format_pool,
      .data$behavior_run_type,
      .data$behavior_effect,
      .data$behavior_condition_code,
      .data$behavior_condition_label,
      .data$association_method,
      .data$association_method_tier,
      .data$analysis_tier,
      .data$neural_level,
      .data$neural_name,
      .data$chrom,
      .data$neural_effect,
      .data$neural_condition_label
    ) %>%
    group_split()

  plot_data_map <- list()
  result_rows <- vector("list", length(split_pairs))
  for (i in seq_along(split_pairs)) {
    sub <- split_pairs[[i]]
    meta <- sub[1, c(
      "behavior_run", "format_pool", "behavior_run_type", "behavior_effect", "behavior_condition_code",
      "behavior_condition_label", "association_method", "association_method_tier", "analysis_tier",
      "neural_level", "neural_name", "chrom", "neural_effect", "neural_condition_label"
    )]
    sub_complete <- sub %>%
      filter(!is.na(.data$behavior_value), !is.na(.data$neural_value)) %>%
      select(subject_id, behavior_value, neural_value)
    pair_key <- make_pair_key(
      behavior_run = meta$behavior_run[[1]],
      format_pool = meta$format_pool[[1]],
      association_method = meta$association_method[[1]],
      neural_name = meta$neural_name[[1]],
      chrom = meta$chrom[[1]]
    )
    stats_row <- compute_association(
      sub_complete = sub_complete,
      alpha = args$alpha,
      min_subjects = args$min_subjects,
      association_method = meta$association_method[[1]]
    )
    plot_data_map[[pair_key]] <- sub_complete
    result_rows[[i]] <- tibble::tibble(
      behavior_run = meta$behavior_run[[1]],
      format_pool = meta$format_pool[[1]],
      behavior_run_type = meta$behavior_run_type[[1]],
      behavior_effect = meta$behavior_effect[[1]],
      behavior_condition_code = meta$behavior_condition_code[[1]],
      behavior_condition_label = meta$behavior_condition_label[[1]],
      association_method = meta$association_method[[1]],
      association_method_tier = meta$association_method_tier[[1]],
      analysis_tier = meta$analysis_tier[[1]],
      neural_level = meta$neural_level[[1]],
      neural_name = meta$neural_name[[1]],
      chrom = meta$chrom[[1]],
      neural_effect = meta$neural_effect[[1]],
      neural_condition_label = meta$neural_condition_label[[1]],
      analysis_status = stats_row$status,
      skip_reason = stats_row$skip_reason,
      n_complete = stats_row$n_complete,
      association_estimate = stats_row$association_estimate,
      r_squared = stats_row$r_squared,
      p_unc = stats_row$p_unc,
      ci95_low = stats_row$ci95_low,
      ci95_high = stats_row$ci95_high,
      slope = stats_row$slope,
      intercept = stats_row$intercept,
      plot_file = NA_character_
    )
  }

  results <- bind_rows(result_rows)
  results$family_id <- make_family_id(results, analysis_plan$multiple_testing$family_grouping)
  results$family_adjust_method <- analysis_plan$multiple_testing$adjust_method
  results$family_n_tested <- NA_integer_
  results$p_fdr <- NA_real_

  tested_idx <- which(results$analysis_status == "tested" & is.finite(results$p_unc))
  if (length(tested_idx) > 0) {
    tested_family_ids <- unique(results$family_id[tested_idx])
    for (family_id in tested_family_ids) {
      family_idx <- which(results$family_id == family_id & results$analysis_status == "tested" & is.finite(results$p_unc))
      results$family_n_tested[family_idx] <- length(family_idx)
      results$p_fdr[family_idx] <- p.adjust(results$p_unc[family_idx], method = analysis_plan$multiple_testing$adjust_method)
    }
  }

  plot_idx <- if (analysis_plan$figures$policy == "all_tested") {
    which(results$analysis_status == "tested" & results$association_method %in% analysis_plan$figures$plot_association_methods)
  } else {
    which(
      results$analysis_status == "tested" &
        is.finite(results$p_unc) &
        results$p_unc < args$alpha &
        results$association_method %in% analysis_plan$figures$plot_association_methods
    )
  }

  for (idx in plot_idx) {
    row <- results[idx, , drop = FALSE]
    pair_key <- make_pair_key(
      behavior_run = row$behavior_run[[1]],
      format_pool = row$format_pool[[1]],
      association_method = row$association_method[[1]],
      neural_name = row$neural_name[[1]],
      chrom = row$chrom[[1]]
    )
    results$plot_file[[idx]] <- plot_association(
      sub_complete = plot_data_map[[pair_key]],
      row = row,
      out_fig_dir = args$out_fig_dir
    )
  }

  results <- results %>%
    arrange(
      dplyr::if_else(is.finite(.data$association_estimate), 0L, 1L),
      desc(.data$association_estimate),
      .data$association_method,
      .data$p_unc,
      .data$behavior_run,
      .data$format_pool,
      .data$neural_name,
      .data$chrom
    )

  write_csv(results, args$out_csv, na = "NA")
  method_out_csvs <- vapply(
    analysis_plan$association_methods,
    function(method) {
      out_path <- method_specific_out_csv(args$out_csv, method)
      write_csv(results %>% filter(.data$association_method == method), out_path, na = "NA")
      out_path
    },
    character(1)
  )

  cat("[data] merged input file:", args$input_csv, "\n")
  cat("[data] ROI definition:", args$roi_json, "\n")
  cat("[data] correlation analysis plan:", args$analysis_plan_json, "\n")
  cat("[data] subjects after exclusions:", length(unique(df_merged$subject_id)), "\n")
  cat("[audit] analyzed target beta columns:", audit$analyzed_beta_columns, "\n")
  cat("[audit] analyzed target beta cells:", audit$analyzed_beta_cells, "\n")
  cat("[audit] exact zero beta cells treated as missing:", audit$exact_zero_cells, "\n")
  cat("[audit] explicit NA beta cells:", audit$explicit_na_cells, "\n")
  cat("[data] behavior runs analyzed:", length(unique(results$behavior_run)), "\n")
  cat("[data] ROI targets analyzed:", results %>% distinct(.data$neural_name, .data$chrom) %>% nrow(), "\n")
  cat("[data] multiple-testing families:", length(unique(results$family_id)), "\n")
  cat("[data] figure policy:", analysis_plan$figures$policy, "\n")
  cat("[out] cleared previous output directory:", cleared_output_root, "\n")
  cat("[out] results CSV:", args$out_csv, "\n")
  for (method in names(method_out_csvs)) {
    cat("[out] ", method, " CSV: ", method_out_csvs[[method]], "\n", sep = "")
  }
  cat("[out] figures:", args$out_fig_dir, "\n")
  cat("[summary] tested rows:", sum(results$analysis_status == "tested"), "\n")
  cat("[summary] plotted rows:", length(plot_idx), "\n")
  cat("[summary] skipped rows:", sum(results$analysis_status != "tested"), "\n")
}

main()
