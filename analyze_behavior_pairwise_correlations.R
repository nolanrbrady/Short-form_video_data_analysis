#!/usr/bin/env Rscript

# Exploratory pairwise behavioral correlations for the SFV study.
#
# Scope
#   - Use the same merged input table consumed by
#     `analyze_correlational_relationships.R`.
#   - Restrict the analysis to an explicit behavior-only variable list declared
#     in `data/config/behavior_pairwise_correlation_plan.json`.
#   - Compute one Pearson correlation per unique unordered behavioral pair.
#
# Rationale
#   - Pearson (1896): product-moment correlation for continuous pairwise
#     behavioral association screening.
#   - Fisher (1921): confidence intervals for Pearson r via Fisher-z.
#   - Bonett (2020): binary-vs-continuous Pearson correlation with 0/1 coding
#     is on the point-biserial effect-size scale, so `pd_status` can remain in
#     the same exploratory table without a separate estimator.
#   - Kriegeskorte et al. (2009): the resulting associations are exploratory
#     and should not be overinterpreted as confirmatory evidence.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
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

`%||%` <- function(a, b) if (!is.null(a)) a else b

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  defaults <- list(
    input_csv = "data/tabular/generated_data/homer3_betas_plus_combined_sfv_data_inner_join.csv",
    analysis_plan_json = "data/config/behavior_pairwise_correlation_plan.json",
    exclude_subjects_json = "data/config/excluded_subjects.json",
    out_dir = "data/results/behavior_pairwise_correlations",
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
    stop(paste0("Behavior pairwise correlation plan must define a non-empty '", field_name, "' array."))
  }
  values <- as.character(values)
  if (any(is.na(values) | trimws(values) == "")) {
    stop(paste0("Behavior pairwise correlation plan contains empty values in '", field_name, "'."))
  }
  values
}

load_analysis_plan <- function(plan_json_path) {
  if (!file.exists(plan_json_path)) {
    stop(paste0("Behavior pairwise correlation plan file not found: ", plan_json_path))
  }

  plan_obj <- tryCatch(
    jsonlite::fromJSON(plan_json_path, simplifyVector = FALSE),
    error = function(e) {
      stop(
        paste0(
          "Failed to parse behavior pairwise correlation plan JSON at ", plan_json_path,
          ". Error: ", conditionMessage(e)
        )
      )
    }
  )

  if (!is.list(plan_obj)) {
    stop("Behavior pairwise correlation plan must be a JSON object.")
  }

  version <- as.integer(plan_obj$version %||% NA_integer_)
  if (!is.finite(version) || version != 1L) {
    stop("Behavior pairwise correlation plan must define version = 1.")
  }

  variables <- normalize_json_string_array(plan_obj$variables, "variables")
  if (anyDuplicated(variables) > 0) {
    dup <- unique(variables[duplicated(variables)])
    stop(
      paste0(
        "Behavior pairwise correlation plan contains duplicate variables: ",
        paste(dup, collapse = ", ")
      )
    )
  }

  figures_obj <- plan_obj$figures
  if (!is.list(figures_obj)) {
    stop("Behavior pairwise correlation plan must define a 'figures' object.")
  }
  figure_policy <- as.character(figures_obj$policy %||% NA_character_)
  if (!nzchar(figure_policy) || !(figure_policy %in% c("significant_only", "all_tested"))) {
    stop("Behavior pairwise correlation plan figures.policy must be 'significant_only' or 'all_tested'.")
  }

  list(
    version = version,
    variables = variables,
    figures = list(policy = figure_policy)
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
        "Failed to parse numeric IDs from column '", column_name, "'. Examples: ",
        paste(bad, collapse = ", ")
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
    stop("Failed to clear the previous behavioral pairwise output directory.")
  }
  TRUE
}

derive_output_paths <- function(out_dir) {
  list(
    out_dir = out_dir,
    out_csv = file.path(out_dir, "behavior_pairwise_correlations_r.csv"),
    out_sig_csv = file.path(out_dir, "behavior_pairwise_correlations_significant_r.csv"),
    out_fig_dir = file.path(out_dir, "figures")
  )
}

load_behavior_input <- function(input_csv, exclude_subjects_json, analysis_plan) {
  df <- read_csv(input_csv, show_col_types = FALSE)
  if (!("subject_id" %in% names(df))) {
    stop(paste0("Expected column 'subject_id' in merged input: ", input_csv))
  }

  assert_required_columns(df, analysis_plan$variables, input_csv)

  df <- df %>%
    mutate(subject_id = normalize_subject_id(.data$subject_id, "subject_id"))

  dup <- df %>%
    count(subject_id, name = "n_rows") %>%
    filter(.data$n_rows > 1)
  if (nrow(dup) > 0) {
    offenders <- dup %>% head(10)
    stop(
      paste0(
        "Duplicate subject_id values detected after normalization in merged input. Examples: ",
        paste0(offenders$subject_id, "=", offenders$n_rows, collapse = ", ")
      )
    )
  }

  df <- coerce_numeric_strict(df, analysis_plan$variables)

  excluded <- apply_subject_exclusions(
    df = df,
    subject_col = "subject_id",
    exclude_json_path = exclude_subjects_json,
    context_label = "behavior_pairwise"
  )

  excluded$data
}

is_binary_zero_one <- function(x) {
  vals <- sort(unique(x[!is.na(x)]))
  length(vals) == 2 && identical(as.numeric(vals), c(0, 1))
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

  if (dplyr::n_distinct(sub_complete$var_x_value) < 2) {
    return(list(
      status = "skipped_constant_input",
      skip_reason = "var_x_has_zero_variance",
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

  if (dplyr::n_distinct(sub_complete$var_y_value) < 2) {
    return(list(
      status = "skipped_constant_input",
      skip_reason = "var_y_has_zero_variance",
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

  cor_fit <- suppressWarnings(stats::cor.test(
    x = sub_complete$var_x_value,
    y = sub_complete$var_y_value,
    method = "pearson",
    alternative = "two.sided",
    conf.level = 1 - alpha
  ))

  lm_fit <- stats::lm(var_y_value ~ var_x_value, data = sub_complete)
  lm_coef <- stats::coef(lm_fit)
  lm_summary <- summary(lm_fit)

  list(
    status = "tested",
    skip_reason = NA_character_,
    n_complete = n_complete,
    pearson_r = unname(cor_fit$estimate),
    r_squared = unname(lm_summary$r.squared),
    p_unc = cor_fit$p.value,
    ci95_low = if (!is.null(cor_fit$conf.int)) cor_fit$conf.int[[1]] else NA_real_,
    ci95_high = if (!is.null(cor_fit$conf.int)) cor_fit$conf.int[[2]] else NA_real_,
    slope = unname(lm_coef[["var_x_value"]]),
    intercept = unname(lm_coef[["(Intercept)"]])
  )
}

plot_pairwise_correlation <- function(sub_complete, row, out_fig_dir) {
  x_binary <- is_binary_zero_one(sub_complete$var_x_value)
  y_binary <- is_binary_zero_one(sub_complete$var_y_value)

  annotation_lines <- c(
    paste0("n = ", row$n_complete),
    paste0("r = ", formatC(row$pearson_r, digits = 3, format = "f")),
    paste0("p = ", format.pval(row$p_unc, digits = 3, eps = 1e-4))
  )
  annotation <- paste(annotation_lines, collapse = "\n")

  file_stub <- sanitize_slug(paste(row$var_x, "vs", row$var_y, sep = "_"))
  file_path <- file.path(out_fig_dir, paste0(file_stub, ".png"))

  point_layer <- if (x_binary || y_binary) {
    geom_jitter(
      width = if (x_binary) 0.06 else 0,
      height = if (y_binary) 0.06 else 0,
      size = 2.2,
      alpha = 0.85,
      color = "#1b4d3e"
    )
  } else {
    geom_point(size = 2.2, alpha = 0.85, color = "#1b4d3e")
  }

  subtitle_bits <- c("Pearson exploratory association")
  if (x_binary || y_binary) {
    subtitle_bits <- c(subtitle_bits, "binary axis jittered for visibility")
  }

  x_anchor <- min(sub_complete$var_x_value, na.rm = TRUE)
  y_anchor <- max(sub_complete$var_y_value, na.rm = TRUE)

  p <- ggplot(sub_complete, aes(x = .data$var_x_value, y = .data$var_y_value)) +
    point_layer +
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
      title = paste0(row$var_x, " vs ", row$var_y),
      subtitle = paste(subtitle_bits, collapse = " | "),
      x = row$var_x,
      y = row$var_y
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
  outputs <- derive_output_paths(args$out_dir)
  cleared_output_root <- clear_output_root(outputs$out_dir)
  dir.create(outputs$out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(outputs$out_fig_dir, recursive = TRUE, showWarnings = FALSE)

  analysis_plan <- load_analysis_plan(args$analysis_plan_json)
  df <- load_behavior_input(args$input_csv, args$exclude_subjects_json, analysis_plan)

  pair_matrix <- utils::combn(analysis_plan$variables, 2)
  result_rows <- vector("list", ncol(pair_matrix))
  plot_data_map <- list()

  for (i in seq_len(ncol(pair_matrix))) {
    var_x <- pair_matrix[1, i]
    var_y <- pair_matrix[2, i]
    sub_complete <- df %>%
      transmute(
        subject_id,
        var_x_value = .data[[var_x]],
        var_y_value = .data[[var_y]]
      ) %>%
      filter(!is.na(.data$var_x_value), !is.na(.data$var_y_value))

    stats_row <- compute_pairwise_correlation(sub_complete, alpha = args$alpha, min_subjects = args$min_subjects)
    pair_key <- paste(var_x, var_y, sep = "||")
    plot_data_map[[pair_key]] <- sub_complete

    result_rows[[i]] <- tibble::tibble(
      var_x = var_x,
      var_y = var_y,
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
    mutate(abs_pearson_r = abs(.data$pearson_r)) %>%
    arrange(is.na(.data$p_unc), .data$p_unc, desc(.data$abs_pearson_r), .data$var_x, .data$var_y) %>%
    select(-abs_pearson_r)

  plot_idx <- if (analysis_plan$figures$policy == "all_tested") {
    which(results$analysis_status == "tested")
  } else {
    which(results$analysis_status == "tested" & is.finite(results$p_unc) & results$p_unc < args$alpha)
  }

  for (idx in plot_idx) {
    row <- results[idx, , drop = FALSE]
    pair_key <- paste(row$var_x[[1]], row$var_y[[1]], sep = "||")
    results$plot_file[[idx]] <- plot_pairwise_correlation(
      sub_complete = plot_data_map[[pair_key]],
      row = row,
      out_fig_dir = outputs$out_fig_dir
    )
  }

  significant_results <- results %>%
    filter(.data$analysis_status == "tested", is.finite(.data$p_unc), .data$p_unc < args$alpha)

  write_csv(results, outputs$out_csv, na = "NA")
  write_csv(significant_results, outputs$out_sig_csv, na = "NA")

  cat("[data] merged input file:", args$input_csv, "\n")
  cat("[data] behavior pairwise plan:", args$analysis_plan_json, "\n")
  cat("[data] subjects after exclusions:", length(unique(df$subject_id)), "\n")
  cat("[data] behavioral variables analyzed:", length(analysis_plan$variables), "\n")
  cat("[data] cleared output root:", if (cleared_output_root) "yes" else "no_existing_dir", "\n")
  cat("[out] results CSV:", outputs$out_csv, "\n")
  cat("[out] significant CSV:", outputs$out_sig_csv, "\n")
  cat("[out] figures:", outputs$out_fig_dir, "\n")
  cat("[summary] tested pairs:", sum(results$analysis_status == "tested"), "\n")
  cat("[summary] significant uncorrected pairs:", nrow(significant_results), "\n")
  cat("[summary] plotted pairs:", length(plot_idx), "\n")
  cat("[summary] skipped pairs:", sum(results$analysis_status != "tested"), "\n")
}

main()
