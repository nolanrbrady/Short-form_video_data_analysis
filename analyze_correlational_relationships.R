#!/usr/bin/env Rscript

# Pairwise Pearson correlations between selected sociodemographic variables and
# selected raw-condition neural summaries from the merged SFV dataset.
#
# Inputs
#   - data/tabular/generated_data/homer3_betas_plus_combined_sfv_data_inner_join.csv
#   - data/config/roi_definition.json
#   - data/config/correlational_analysis_plan.json
#
# Neural targets
#   - Channel S04_D02: HbO and HbR, Cond01-Cond04
#   - ROI means (raw condition betas):
#       R_DLPFC (HbR), L_DLPFC (HbO), M_DMPFC (HbO), L_DMPFC (HbO)
#
# Predictor set
#   - Loaded explicitly from `data/config/correlational_analysis_plan.json`
#   - This keeps the tested predictor family stable across merged-schema changes.
#
# Missingness policy
#   - Pairwise complete cases only: for a given predictor x neural target test,
#     exclude any subject missing either value for that pair.
#   - Pruned channels remain missing; do not impute.
#
# Multiple testing correction
#   - BH-FDR within each configured neural-target x chromophore family.
#   - By default, one family pools all condition x predictor pairs for a given
#     `(neural_level, neural_name, chrom)` combination.
#
# Citations (see CITATIONS.md)
#   - Pearson (1896): product-moment correlation.
#   - Fisher (1921): Fisher z transformation for confidence intervals on r.
#   - Benjamini & Hochberg (1995): BH-FDR.
#   - Poldrack (2007): ROI signal extraction by averaging across pre-specified channels.

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
  ~cond, ~condition_label,
  "01", "SF_Edu",
  "02", "SF_Ent",
  "03", "LF_Ent",
  "04", "LF_Edu"
)

TARGET_CHANNEL <- "S04_D02"
TARGET_CHANNEL_CHROMS <- c("HbO", "HbR")
TARGET_ROI_SPECS <- tibble::tribble(
  ~neural_name, ~chrom,
  "R_DLPFC", "HbR",
  "L_DLPFC", "HbO",
  "M_DMPFC", "HbO",
  "L_DMPFC", "HbO"
)

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  defaults <- list(
    input_csv = "data/tabular/generated_data/homer3_betas_plus_combined_sfv_data_inner_join.csv",
    roi_json = "data/config/roi_definition.json",
    analysis_plan_json = "data/config/correlational_analysis_plan.json",
    exclude_subjects_json = "data/config/excluded_subjects.json",
    out_csv = "data/results/correlational_relationships/pairwise_correlations_r.csv",
    out_fig_dir = "data/results/correlational_relationships/figures",
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

`%||%` <- function(a, b) if (!is.null(a)) a else b

normalize_json_string_array <- function(x, field_name) {
  values <- if (is.list(x) && !is.data.frame(x)) {
    unlist(x, use.names = FALSE)
  } else {
    x
  }
  if (!is.atomic(values) || length(values) == 0) {
    stop(paste0("Correlation analysis plan must define a non-empty '", field_name, "' array."))
  }
  values <- as.character(values)
  if (any(is.na(values) | trimws(values) == "")) {
    stop(paste0("Correlation analysis plan contains empty values in '", field_name, "'."))
  }
  values
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
          ". Use strict JSON syntax. Error: ", conditionMessage(e)
        )
      )
    }
  )

  if (!is.list(plan_obj)) {
    stop("Correlation analysis plan must be a JSON object.")
  }

  predictors <- plan_obj$predictors
  if (is.null(predictors)) {
    stop("Correlation analysis plan must define a non-empty 'predictors' array.")
  }
  predictors <- normalize_json_string_array(predictors, "predictors")
  if (anyDuplicated(predictors) > 0) {
    dup <- unique(predictors[duplicated(predictors)])
    stop(
      paste0(
        "Correlation analysis plan contains duplicate predictors: ",
        paste(dup, collapse = ", ")
      )
    )
  }

  mt <- plan_obj$multiple_testing
  if (!is.list(mt)) {
    stop("Correlation analysis plan must define a 'multiple_testing' object.")
  }
  adjust_method <- as.character(mt$adjust_method %||% NA_character_)
  family_grouping <- mt$family_grouping
  if (!nzchar(adjust_method)) {
    stop("Correlation analysis plan must define multiple_testing.adjust_method.")
  }
  supported_methods <- c("BH", "holm", "bonferroni", "hochberg", "hommel", "BY", "fdr", "none")
  if (!(adjust_method %in% supported_methods)) {
    stop(
      paste0(
        "Unsupported multiple_testing.adjust_method '", adjust_method,
        "'. Supported methods: ", paste(supported_methods, collapse = ", ")
      )
    )
  }
  if (is.null(family_grouping)) {
    stop("Correlation analysis plan must define a non-empty multiple_testing.family_grouping array.")
  }
  family_grouping <- normalize_json_string_array(family_grouping, "multiple_testing.family_grouping")
  if (anyDuplicated(family_grouping) > 0) {
    dup <- unique(family_grouping[duplicated(family_grouping)])
    stop(
      paste0(
        "Correlation analysis plan contains duplicate family grouping keys: ",
        paste(dup, collapse = ", ")
      )
    )
  }

  list(
    predictors = predictors,
    multiple_testing = list(
      adjust_method = adjust_method,
      family_grouping = family_grouping
    )
  )
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
        "Failed to parse channel IDs in ", context_label, ". ",
        "Expected format like S01_D01. Examples: ",
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
          "Column '", col_name, "' contains non-numeric values. ",
          "Examples: ", paste(examples, collapse = ", ")
        )
      )
    }
    out[[col_name]] <- num
  }
  out
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
    roi_rows[[length(roi_rows) + 1]] <- tibble::tibble(roi = roi_name, channel = channels_norm)
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

collect_planned_predictor_columns <- function(df_names, analysis_plan) {
  planned_predictors <- analysis_plan$predictors
  missing_predictors <- setdiff(planned_predictors, df_names)
  if (length(missing_predictors) > 0) {
    stop(
      paste0(
        "Merged input is missing required predictors: ",
        paste(sort(missing_predictors), collapse = ", ")
      )
    )
  }

  discovered_dynamic <- sort(df_names[str_detect(df_names, "^diff_|_engagement$")])
  unexpected_dynamic <- setdiff(discovered_dynamic, planned_predictors)
  if (length(unexpected_dynamic) > 0) {
    stop(
      paste0(
        "Merged input contains derived predictors that are not declared in the correlation analysis plan: ",
        paste(unexpected_dynamic, collapse = ", "),
        ". Update data/config/correlational_analysis_plan.json before re-running."
      )
    )
  }

  planned_predictors
}

load_merged_input <- function(input_csv, exclude_subjects_json, analysis_plan) {
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

  predictor_cols <- collect_planned_predictor_columns(names(df), analysis_plan)
  required_channel_cols <- as.vector(outer(
    paste0(TARGET_CHANNEL, "_Cond", CONDITION_MAP$cond),
    TARGET_CHANNEL_CHROMS,
    paste,
    sep = "_"
  ))
  assert_required_columns(df, c(predictor_cols, required_channel_cols), input_csv)

  df <- df %>%
    mutate(subject_id = normalize_subject_id(.data$subject_id, "subject_id"))

  dup <- df %>%
    count(subject_id, name = "n_rows") %>%
    filter(.data$n_rows > 1) %>%
    arrange(desc(n_rows), subject_id)
  if (nrow(dup) > 0) {
    examples <- dup %>% head(10)
    stop(
      paste0(
        "Duplicate subject_id values detected after normalization in merged input. ",
        "Expected exactly one row per subject. ",
        "Example duplicates (subject_id, n_rows): ",
        paste0(examples$subject_id, "=", examples$n_rows, collapse = ", "),
        ". Fix the upstream CSV before running the correlation analysis."
      )
    )
  }

  numeric_cols <- unique(c(predictor_cols, beta_cols))
  df <- coerce_numeric_strict(df, numeric_cols)

  excluded <- apply_subject_exclusions(
    df = df,
    subject_col = "subject_id",
    exclude_json_path = exclude_subjects_json,
    context_label = "correlation"
  )

  list(data = excluded$data, predictor_cols = predictor_cols, beta_cols = beta_cols)
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
    mutate(channel = normalize_channel_id(.data$channel, "beta columns")) %>%
    left_join(CONDITION_MAP, by = "cond") %>%
    {
      bad_conds <- setdiff(unique(.$cond), CONDITION_MAP$cond)
      if (length(bad_conds) > 0) {
        stop(paste0("Unexpected condition codes in beta columns: ", paste(sort(bad_conds), collapse = ", ")))
      }
      .
    }
}

build_channel_targets <- function(beta_long) {
  out <- beta_long %>%
    filter(.data$channel == TARGET_CHANNEL, .data$chrom %in% TARGET_CHANNEL_CHROMS) %>%
    transmute(
      subject_id,
      neural_level = "channel",
      neural_name = .data$channel,
      chrom,
      condition = paste0("Cond", .data$cond),
      condition_label,
      neural_value = .data$beta
    )

  expected <- length(TARGET_CHANNEL_CHROMS) * nrow(CONDITION_MAP)
  seen <- out %>% distinct(neural_name, chrom, condition) %>% nrow()
  if (seen != expected) {
    stop(
      paste0(
        "Expected ", expected, " requested S04_D02 channel targets but found ", seen, "."
      )
    )
  }
  out
}

build_roi_targets <- function(beta_long, roi_map) {
  roi_map_target <- roi_map %>%
    semi_join(TARGET_ROI_SPECS, by = c("roi" = "neural_name"))

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

  out <- beta_long %>%
    inner_join(roi_map_target, by = "channel") %>%
    rename(neural_name = roi) %>%
    inner_join(TARGET_ROI_SPECS, by = c("neural_name", "chrom")) %>%
    group_by(.data$subject_id, .data$neural_name, .data$chrom, .data$cond, .data$condition_label) %>%
    summarize(
      # ROI averaging across available channels follows a standard ROI summary
      # strategy in neuroimaging; see Poldrack (2007) in `CITATIONS.md`.
      neural_value = if (all(is.na(.data$beta))) NA_real_ else mean(.data$beta, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    transmute(
      subject_id,
      neural_level = "roi",
      neural_name,
      chrom,
      condition = paste0("Cond", .data$cond),
      condition_label,
      neural_value
    )

  expected <- nrow(TARGET_ROI_SPECS) * nrow(CONDITION_MAP)
  seen <- out %>% distinct(neural_name, chrom, condition) %>% nrow()
  if (seen != expected) {
    stop(
      paste0(
        "Expected ", expected, " requested ROI targets but found ", seen, "."
      )
    )
  }
  out
}

build_predictor_long <- function(df_merged, predictor_cols) {
  df_merged %>%
    select(subject_id, all_of(predictor_cols)) %>%
    pivot_longer(cols = all_of(predictor_cols), names_to = "predictor", values_to = "predictor_value")
}

make_pair_key <- function(neural_level, neural_name, chrom, condition, predictor) {
  paste(neural_level, neural_name, chrom, condition, predictor, sep = "||")
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

  apply(
    df[, family_grouping, drop = FALSE],
    1,
    function(row) paste(family_grouping, row, sep = "=", collapse = " | ")
  )
}

sanitize_slug <- function(x) {
  x %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("^_+|_+$", "") %>%
    str_to_lower()
}

compute_pairwise_correlation <- function(sub_complete, alpha, min_subjects) {
  n_complete <- nrow(sub_complete)
  if (n_complete < min_subjects) {
    return(list(
      status = "skipped_min_subjects",
      skip_reason = paste0("n_complete<", min_subjects),
      n_complete = n_complete,
      pearson_r = NA_real_,
      r_squared = NA_real_,
      p_unc = NA_real_,
      ci95_low = NA_real_,
      ci95_high = NA_real_,
      slope = NA_real_,
      intercept = NA_real_
    ))
  }

  if (dplyr::n_distinct(sub_complete$predictor_value) < 2) {
    return(list(
      status = "skipped_constant_input",
      skip_reason = "predictor_has_zero_variance",
      n_complete = n_complete,
      pearson_r = NA_real_,
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
      pearson_r = NA_real_,
      r_squared = NA_real_,
      p_unc = NA_real_,
      ci95_low = NA_real_,
      ci95_high = NA_real_,
      slope = NA_real_,
      intercept = NA_real_
    ))
  }

  # Pearson product-moment correlation and its Fisher-z confidence interval.
  # References: Pearson (1896); Fisher (1921), see `CITATIONS.md`.
  cor_fit <- suppressWarnings(stats::cor.test(
    x = sub_complete$predictor_value,
    y = sub_complete$neural_value,
    method = "pearson",
    alternative = "two.sided",
    conf.level = 1 - alpha
  ))

  lm_fit <- stats::lm(neural_value ~ predictor_value, data = sub_complete)
  lm_coef <- stats::coef(lm_fit)

  list(
    status = "tested",
    skip_reason = NA_character_,
    n_complete = n_complete,
    pearson_r = unname(cor_fit$estimate),
    r_squared = unname(cor_fit$estimate)^2,
    p_unc = cor_fit$p.value,
    ci95_low = if (!is.null(cor_fit$conf.int)) cor_fit$conf.int[[1]] else NA_real_,
    ci95_high = if (!is.null(cor_fit$conf.int)) cor_fit$conf.int[[2]] else NA_real_,
    slope = unname(lm_coef[["predictor_value"]]),
    intercept = unname(lm_coef[["(Intercept)"]])
  )
}

plot_pairwise_correlation <- function(sub_complete, row, out_fig_dir) {
  annotation_lines <- c(
    paste0("n = ", row$n_complete),
    paste0("r = ", formatC(row$pearson_r, digits = 3, format = "f")),
    paste0("p = ", format.pval(row$p_unc, digits = 3, eps = 1e-4)),
    paste0("q = ", format.pval(row$p_fdr, digits = 3, eps = 1e-4))
  )
  annotation <- paste(annotation_lines, collapse = "\n")

  file_stub <- sanitize_slug(paste(
    row$neural_level,
    row$neural_name,
    row$chrom,
    row$condition,
    row$predictor,
    sep = "_"
  ))
  file_path <- file.path(out_fig_dir, paste0(file_stub, ".png"))

  plot_title <- paste0(
    row$neural_name, " ", row$chrom, " ", row$condition,
    " vs ", row$predictor
  )
  subtitle <- paste0(
    "Condition: ", row$condition_label,
    " | Target: ", row$neural_level
  )

  x_anchor <- min(sub_complete$predictor_value, na.rm = TRUE)
  y_anchor <- max(sub_complete$neural_value, na.rm = TRUE)

  p <- ggplot(sub_complete, aes(x = .data$predictor_value, y = .data$neural_value)) +
    geom_point(size = 2.2, alpha = 0.85, color = "#1b4d3e") +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#c04b2c", fill = "#f1c9b8") +
    annotate(
      "label",
      x = x_anchor,
      y = y_anchor,
      label = annotation,
      hjust = 0,
      vjust = 1,
      linewidth = 0.25,
      size = 3.2
    ) +
    labs(
      title = plot_title,
      subtitle = subtitle,
      x = row$predictor,
      y = "Beta / ROI mean beta"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )

  suppressMessages(
    ggplot2::ggsave(filename = file_path, plot = p, width = 7, height = 5, dpi = 300)
  )
  file_path
}

main <- function() {
  args <- parse_args()
  dir.create(dirname(args$out_csv), recursive = TRUE, showWarnings = FALSE)
  dir.create(args$out_fig_dir, recursive = TRUE, showWarnings = FALSE)

  analysis_plan <- load_analysis_plan(args$analysis_plan_json)
  loaded <- load_merged_input(args$input_csv, args$exclude_subjects_json, analysis_plan)
  df_merged <- loaded$data
  predictor_cols <- loaded$predictor_cols
  beta_cols <- loaded$beta_cols
  roi_map <- load_roi_definition(args$roi_json)

  beta_long <- reshape_beta_long(df_merged, beta_cols)
  channel_targets <- build_channel_targets(beta_long)
  roi_targets <- build_roi_targets(beta_long, roi_map)
  neural_targets <- bind_rows(channel_targets, roi_targets)

  predictor_long <- build_predictor_long(df_merged, predictor_cols)
  pair_grid <- predictor_long %>%
    inner_join(neural_targets, by = "subject_id", relationship = "many-to-many")

  split_pairs <- pair_grid %>%
    group_by(.data$neural_level, .data$neural_name, .data$chrom, .data$condition, .data$condition_label, .data$predictor) %>%
    group_split()

  plot_data_map <- list()
  result_rows <- vector("list", length(split_pairs))

  for (i in seq_along(split_pairs)) {
    sub <- split_pairs[[i]]
    meta <- sub[1, c("neural_level", "neural_name", "chrom", "condition", "condition_label", "predictor")]
    sub_complete <- sub %>%
      filter(!is.na(.data$predictor_value), !is.na(.data$neural_value)) %>%
      select(subject_id, predictor_value, neural_value)

    pair_key <- make_pair_key(
      neural_level = meta$neural_level[[1]],
      neural_name = meta$neural_name[[1]],
      chrom = meta$chrom[[1]],
      condition = meta$condition[[1]],
      predictor = meta$predictor[[1]]
    )
    stats_row <- compute_pairwise_correlation(sub_complete, alpha = args$alpha, min_subjects = args$min_subjects)
    plot_data_map[[pair_key]] <- sub_complete

    result_rows[[i]] <- tibble::tibble(
      neural_level = meta$neural_level[[1]],
      neural_name = meta$neural_name[[1]],
      chrom = meta$chrom[[1]],
      condition = meta$condition[[1]],
      condition_code = meta$condition[[1]],
      condition_label = meta$condition_label[[1]],
      predictor = meta$predictor[[1]],
      analysis_status = stats_row$status,
      skip_reason = stats_row$skip_reason,
      n_complete = stats_row$n_complete,
      pearson_r = stats_row$pearson_r,
      r_squared = stats_row$r_squared,
      p_unc = stats_row$p_unc,
      ci95_low = stats_row$ci95_low,
      ci95_high = stats_row$ci95_high,
      slope = stats_row$slope,
      intercept = stats_row$intercept,
      plot_file = NA_character_
    )
  }

  results <- bind_rows(result_rows) %>%
    arrange(neural_level, neural_name, chrom, condition, predictor)
  results$family_id <- make_family_id(results, analysis_plan$multiple_testing$family_grouping)
  results$family_adjust_method <- analysis_plan$multiple_testing$adjust_method
  results$family_n_tested <- NA_integer_

  tested_idx <- which(results$analysis_status == "tested" & is.finite(results$p_unc))
  results$p_fdr <- NA_real_
  if (length(tested_idx) > 0) {
    tested_family_ids <- unique(results$family_id[tested_idx])
    for (family_id in tested_family_ids) {
      family_idx <- which(results$family_id == family_id & results$analysis_status == "tested" & is.finite(results$p_unc))
      results$family_n_tested[family_idx] <- length(family_idx)
      results$p_fdr[family_idx] <- p.adjust(
        results$p_unc[family_idx],
        method = analysis_plan$multiple_testing$adjust_method
      )
    }
  }

  for (idx in tested_idx) {
    row <- results[idx, , drop = FALSE]
    pair_key <- make_pair_key(
      neural_level = row$neural_level[[1]],
      neural_name = row$neural_name[[1]],
      chrom = row$chrom[[1]],
      condition = row$condition[[1]],
      predictor = row$predictor[[1]]
    )
    results$plot_file[[idx]] <- plot_pairwise_correlation(
      sub_complete = plot_data_map[[pair_key]],
      row = row,
      out_fig_dir = args$out_fig_dir
    )
  }

  write_csv(results, args$out_csv, na = "NA")

  target_count <- neural_targets %>%
    distinct(neural_level, neural_name, chrom, condition) %>%
    nrow()

  cat("[data] merged input file:", args$input_csv, "\n")
  cat("[data] ROI definition:", args$roi_json, "\n")
  cat("[data] correlation analysis plan:", args$analysis_plan_json, "\n")
  cat("[data] subjects after exclusions:", length(unique(df_merged$subject_id)), "\n")
  cat("[data] predictors analyzed:", length(predictor_cols), "\n")
  cat("[data] neural targets analyzed:", target_count, "\n")
  cat("[data] multiple-testing families:", length(unique(results$family_id)), "\n")
  cat("[out] results CSV:", args$out_csv, "\n")
  cat("[out] figures:", args$out_fig_dir, "\n")
  cat("[summary] tested pairs:", sum(results$analysis_status == "tested"), "\n")
  cat("[summary] skipped pairs:", sum(results$analysis_status != "tested"), "\n")
}

main()
