#!/usr/bin/env Rscript

# Exploratory pairwise behavioral correlations for the SFV study.
#
# Scope
#   - Use the merged input table consumed by the primary LMMs and the selected
#     `analyze_pooled_mean_correlations.R` follow-up.
#   - Restrict the analysis to an explicit behavior-only variable list declared
#     in `data/config/behavior_pairwise_correlation_plan.json`.
#   - Compute one Pearson correlation per unique unordered behavioral pair.
#   - Exclude recruitment order and subject-ID-derived proxies from the
#     behavioral correlation family.
#
# Rationale
#   - Pearson (1896): product-moment correlation for continuous pairwise
#     behavioral association screening.
#   - Fisher (1921): confidence intervals for Pearson r via Fisher-z.
#   - Benjamini & Hochberg (1995): global FDR correction across the full
#     behavioral pairwise screening family.
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
    variable_figure_names_json = "data/config/variable_figure_names.json",
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

load_variable_figure_labels <- function(label_json_path, variables) {
  if (!file.exists(label_json_path)) {
    stop(paste0("Variable figure-name config file not found: ", label_json_path))
  }

  label_obj <- tryCatch(
    jsonlite::fromJSON(label_json_path, simplifyVector = FALSE),
    error = function(e) {
      stop(
        paste0(
          "Failed to parse variable figure-name JSON at ", label_json_path,
          ". Error: ", conditionMessage(e)
        )
      )
    }
  )
  if (!is.list(label_obj) || is.data.frame(label_obj) || length(label_obj) == 0) {
    stop("Variable figure-name config must be a non-empty JSON object mapping variable names to labels.")
  }
  if (is.null(names(label_obj)) || any(!nzchar(names(label_obj)))) {
    stop("Variable figure-name config must use variable names as object keys.")
  }

  labels <- vapply(names(label_obj), function(name) {
    value <- label_obj[[name]]
    if (!is.atomic(value) || length(value) != 1) {
      stop(paste0("Variable figure-name label for '", name, "' must be a single string."))
    }
    as.character(value)
  }, character(1))

  if (any(is.na(labels) | trimws(labels) == "")) {
    bad <- names(labels)[is.na(labels) | trimws(labels) == ""]
    stop(
      paste0(
        "Variable figure-name config contains empty labels for: ",
        paste(head(bad, 10), collapse = ", ")
      )
    )
  }

  missing_labels <- setdiff(variables, names(labels))
  if (length(missing_labels) > 0) {
    stop(
      paste0(
        "Variable figure-name config is missing labels for: ",
        paste(missing_labels, collapse = ", ")
      )
    )
  }

  labels
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

  if ("recruitment_order_proxy" %in% variables) {
    stop(
      paste0(
        "Behavior pairwise correlation plan must not include ",
        "'recruitment_order_proxy'; recruitment order is excluded from this analysis."
      )
    )
  }

  figures_obj <- plan_obj$figures %||% list()
  if (!is.list(figures_obj)) {
    stop("Behavior pairwise correlation plan must define a 'figures' object.")
  }
  lower_triangle_obj <- figures_obj$lower_triangle %||% list()
  if (!is.list(lower_triangle_obj)) {
    stop("Behavior pairwise correlation plan figures.lower_triangle must be an object when provided.")
  }
  filename_stem <- as.character(lower_triangle_obj$filename_stem %||% "behavior_pairwise_correlation_lower_triangle")
  if (!nzchar(filename_stem) || grepl("[/\\\\]", filename_stem)) {
    stop("Behavior pairwise correlation plan figures.lower_triangle.filename_stem must be a non-empty file stem.")
  }

  list(
    version = version,
    variables = variables,
    figures = list(
      lower_triangle = list(filename_stem = filename_stem)
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
    out_fdr_csv = file.path(out_dir, "behavior_pairwise_correlations_fdr_r.csv"),
    out_fig_dir = file.path(out_dir, "figures")
  )
}

add_matrix_output_paths <- function(outputs, filename_stem) {
  outputs$matrix_png <- file.path(outputs$out_fig_dir, paste0(filename_stem, ".png"))
  outputs$matrix_pdf <- file.path(outputs$out_fig_dir, paste0(filename_stem, ".pdf"))
  outputs
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

compute_pairwise_correlation <- function(sub_complete, alpha, min_subjects) {
  n_complete <- nrow(sub_complete)
  if (n_complete < min_subjects) {
    return(list(
      status = "skipped_min_subjects",
      skip_reason = paste0("n_complete<", min_subjects),
      n_complete = n_complete,
      pearson_r = NA_real_,
      p_unc = NA_real_,
      ci95_low = NA_real_,
      ci95_high = NA_real_
    ))
  }

  if (dplyr::n_distinct(sub_complete$var_x_value) < 2) {
    return(list(
      status = "skipped_constant_input",
      skip_reason = "var_x_has_zero_variance",
      n_complete = n_complete,
      pearson_r = NA_real_,
      p_unc = NA_real_,
      ci95_low = NA_real_,
      ci95_high = NA_real_
    ))
  }

  if (dplyr::n_distinct(sub_complete$var_y_value) < 2) {
    return(list(
      status = "skipped_constant_input",
      skip_reason = "var_y_has_zero_variance",
      n_complete = n_complete,
      pearson_r = NA_real_,
      p_unc = NA_real_,
      ci95_low = NA_real_,
      ci95_high = NA_real_
    ))
  }

  cor_fit <- suppressWarnings(stats::cor.test(
    x = sub_complete$var_x_value,
    y = sub_complete$var_y_value,
    method = "pearson",
    alternative = "two.sided",
    conf.level = 1 - alpha
  ))

  list(
    status = "tested",
    skip_reason = NA_character_,
    n_complete = n_complete,
    pearson_r = unname(cor_fit$estimate),
    p_unc = cor_fit$p.value,
    ci95_low = if (!is.null(cor_fit$conf.int)) cor_fit$conf.int[[1]] else NA_real_,
    ci95_high = if (!is.null(cor_fit$conf.int)) cor_fit$conf.int[[2]] else NA_real_
  )
}

strip_leading_zero <- function(x) {
  x %>%
    str_replace("^0\\.", ".") %>%
    str_replace("^-0\\.", "-.")
}

format_r_value <- function(x) {
  # Round first so that tiny negatives (e.g. -0.001) do not render as "-.00".
  rounded <- round(x, 2)
  if (isTRUE(rounded == 0)) {
    rounded <- 0
  }
  strip_leading_zero(formatC(rounded, digits = 2, format = "f"))
}

# Significance markers keyed to the global Benjamini-Hochberg FDR q-value so the
# figure stays readable while exact q-values remain available in the CSV outputs.
significance_stars <- function(q) {
  if (is.na(q)) {
    return("")
  }
  if (q < 0.001) {
    return("***")
  }
  if (q < 0.01) {
    return("**")
  }
  if (q < 0.05) {
    return("*")
  }
  ""
}

display_variable_label <- function(x, variable_labels) {
  if (!(x %in% names(variable_labels))) {
    stop(paste0("No figure label loaded for variable: ", x))
  }
  str_wrap(unname(variable_labels[[x]]), width = 14)
}

format_cell_label <- function(r, q) {
  if (is.na(r)) {
    return("")
  }
  paste0(format_r_value(r), significance_stars(q))
}

plot_lower_triangle_correlation_matrix <- function(results, variables, outputs, variable_labels) {
  variable_index <- seq_along(variables)
  names(variable_index) <- variables
  diagonal_labels <- vapply(variables, display_variable_label, character(1), variable_labels = variable_labels)
  n_variables <- length(variables)

  # Diagonal tiles carry the wrapped variable names so the matrix is
  # self-documenting; this removes the cramped rotated axis labels entirely.
  diagonal_df <- tibble::tibble(
    idx = variable_index,
    label = diagonal_labels
  )

  plot_df <- results %>%
    filter(.data$analysis_status == "tested") %>%
    mutate(
      var_x_index = unname(variable_index[.data$var_x]),
      var_y_index = unname(variable_index[.data$var_y]),
      cell_label = vapply(
        seq_len(n()),
        function(i) format_cell_label(.data$pearson_r[[i]], .data$p_fdr[[i]]),
        character(1)
      ),
      label_face = if_else(.data$significant_fdr %in% TRUE, "bold", "plain"),
      # White text stays legible on the saturated (strong-correlation) tiles,
      # near-black text on the pale mid-range tiles.
      text_color = if_else(abs(.data$pearson_r) > 0.5, "white", "grey15")
    )

  p <- ggplot(plot_df, aes(x = .data$var_x_index, y = .data$var_y_index)) +
    geom_tile(aes(fill = .data$pearson_r), color = "white", linewidth = 0.6) +
    # Neutral diagonal tiles hosting the variable names.
    geom_tile(
      data = diagonal_df,
      aes(x = .data$idx, y = .data$idx),
      fill = "grey93",
      color = "white",
      linewidth = 0.6,
      inherit.aes = FALSE
    ) +
    geom_text(
      aes(label = .data$cell_label, fontface = .data$label_face, color = .data$text_color),
      size = 2.85
    ) +
    geom_text(
      data = diagonal_df,
      aes(x = .data$idx, y = .data$idx, label = .data$label),
      inherit.aes = FALSE,
      size = 2.5,
      lineheight = 0.9,
      fontface = "bold",
      color = "grey20"
    ) +
    scale_color_identity() +
    scale_fill_gradient2(
      low = "#2166ac",
      mid = "#f7f7f7",
      high = "#b2182b",
      midpoint = 0,
      limits = c(-1, 1),
      breaks = c(-1, -0.5, 0, 0.5, 1),
      name = expression("Pearson " * italic(r)),
      guide = guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = grid::unit(15, "lines"),
        barheight = grid::unit(1.2, "lines"),
        ticks.colour = "grey30",
        frame.colour = "grey30"
      )
    ) +
    scale_x_continuous(limits = c(0.5, n_variables + 0.5), expand = c(0, 0)) +
    scale_y_reverse(limits = c(n_variables + 0.5, 0.5), expand = c(0, 0)) +
    coord_fixed(clip = "off") +
    labs(
      title = "Pairwise Behavioral Correlations",
      subtitle = bquote("Pearson " * italic(r) * "; asterisks denote Benjamini-Hochberg FDR-adjusted significance (" *
        "*" * italic(q) * " < .05, ** " * italic(q) * " < .01, *** " * italic(q) * " < .001)."),
      x = NULL,
      y = NULL
    ) +
    # Park the legend inside the otherwise-empty upper-right triangle.
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(face = "bold", size = 17, margin = margin(b = 3)),
      plot.subtitle = element_text(size = 10, margin = margin(b = 6), color = "grey25"),
      legend.position = c(0.8, 0.8),
      legend.direction = "horizontal",
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 11),
      plot.margin = margin(14, 14, 12, 14)
    )

  suppressMessages({
    ggplot2::ggsave(filename = outputs$matrix_png, plot = p, width = 12, height = 12.5, dpi = 300, bg = "white")
    ggplot2::ggsave(filename = outputs$matrix_pdf, plot = p, width = 12, height = 12.5, bg = "white")
  })
}

main <- function() {
  args <- parse_args()
  analysis_plan <- load_analysis_plan(args$analysis_plan_json)
  variable_labels <- load_variable_figure_labels(args$variable_figure_names_json, analysis_plan$variables)
  outputs <- derive_output_paths(args$out_dir)
  outputs <- add_matrix_output_paths(outputs, analysis_plan$figures$lower_triangle$filename_stem)
  cleared_output_root <- clear_output_root(outputs$out_dir)
  dir.create(outputs$out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(outputs$out_fig_dir, recursive = TRUE, showWarnings = FALSE)

  df <- load_behavior_input(args$input_csv, args$exclude_subjects_json, analysis_plan)

  pair_matrix <- utils::combn(analysis_plan$variables, 2)
  result_rows <- vector("list", ncol(pair_matrix))

  for (i in seq_len(ncol(pair_matrix))) {
    var_x <- pair_matrix[1, i]
    var_y <- pair_matrix[2, i]
    sub_complete <- df %>%
      transmute(
        var_x_value = .data[[var_x]],
        var_y_value = .data[[var_y]]
      ) %>%
      filter(!is.na(.data$var_x_value), !is.na(.data$var_y_value))

    stats_row <- compute_pairwise_correlation(sub_complete, alpha = args$alpha, min_subjects = args$min_subjects)

    result_rows[[i]] <- tibble::tibble(
      var_x = var_x,
      var_y = var_y,
      analysis_status = stats_row$status,
      skip_reason = stats_row$skip_reason,
      n_complete = stats_row$n_complete,
      pearson_r = stats_row$pearson_r,
      p_unc = stats_row$p_unc,
      p_fdr = NA_real_,
      significant_fdr = NA,
      ci95_low = stats_row$ci95_low,
      ci95_high = stats_row$ci95_high
    )
  }

  results <- bind_rows(result_rows)
  tested_idx <- results$analysis_status == "tested" & is.finite(results$p_unc)
  if (any(tested_idx)) {
    # Benjamini & Hochberg (1995; see CITATIONS.md): one global FDR family for
    # this exploratory behavioral-pairwise screen.
    results$p_fdr[tested_idx] <- stats::p.adjust(results$p_unc[tested_idx], method = "BH")
    results$significant_fdr[tested_idx] <- results$p_fdr[tested_idx] < args$alpha
  }

  results <- results %>%
    mutate(abs_pearson_r = abs(.data$pearson_r)) %>%
    arrange(
      is.na(.data$p_fdr),
      .data$p_fdr,
      is.na(.data$p_unc),
      .data$p_unc,
      desc(.data$abs_pearson_r),
      .data$var_x,
      .data$var_y
    ) %>%
    select(-abs_pearson_r)

  fdr_results <- results %>%
    filter(.data$analysis_status == "tested")

  plot_lower_triangle_correlation_matrix(results, analysis_plan$variables, outputs, variable_labels)

  write_csv(results, outputs$out_csv, na = "NA")
  write_csv(fdr_results, outputs$out_fdr_csv, na = "NA")

  cat("[data] merged input file:", args$input_csv, "\n")
  cat("[data] behavior pairwise plan:", args$analysis_plan_json, "\n")
  cat("[data] variable figure-name config:", args$variable_figure_names_json, "\n")
  cat("[data] subjects after exclusions:", length(unique(df$subject_id)), "\n")
  cat("[data] behavioral variables analyzed:", length(analysis_plan$variables), "\n")
  cat("[data] cleared output root:", if (cleared_output_root) "yes" else "no_existing_dir", "\n")
  cat("[out] results CSV:", outputs$out_csv, "\n")
  cat("[out] FDR CSV:", outputs$out_fdr_csv, "\n")
  cat("[out] lower triangle PNG:", outputs$matrix_png, "\n")
  cat("[out] lower triangle PDF:", outputs$matrix_pdf, "\n")
  cat("[summary] tested pairs:", sum(results$analysis_status == "tested"), "\n")
  cat("[summary] significant FDR pairs:", sum(results$significant_fdr %in% TRUE), "\n")
  cat("[summary] skipped pairs:", sum(results$analysis_status != "tested"), "\n")
}

main()
