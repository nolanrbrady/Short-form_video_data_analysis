#!/usr/bin/env Rscript

# Plot publication-ready subject-level behavioral score distributions for the
# four study conditions and planned descriptive marginals for the score factors.
#
# Inputs
#   - data/tabular/generated_data/homer3_betas_plus_combined_sfv_data_inner_join.csv
#   - data/config/excluded_subjects.json
#
# Behavior
#   - Applies the shared subject-exclusion manifest so plotted subjects align
#     with inferential analyses (Sandve et al., 2013; see CITATIONS.md).
#   - Plots engagement scores and retention/recall improvement scores
#     (`post - pre`) as separate PNG figures across the four conditions.
#   - Adds a retention length-marginal figure by averaging short- and
#     long-form conditions within subject.
#   - Uses violin density envelopes with raw jittered subject points and
#     mean +/- 1 SD overlays instead of summary-only bars so readers can inspect
#     spread, sample size, and unusual values (Hintze & Nelson, 1998;
#     Weissgerber et al., 2015; see CITATIONS.md).
#
# Missingness policy
#   - Engagement and retention zeros are valid observed values.
#   - Missing values are not imputed. Each behavioral domain uses its own
#     complete-case subject set: a subject must have all four condition values
#     present for that domain to be plotted.
#
# Citations (see CITATIONS.md)
#   - Sandve et al. (2013): centralized, auditable exclusion manifest.
#   - Hintze & Nelson (1998): violin plots for distribution visualization.
#   - Weissgerber et al. (2015): raw-data visualization over summary-only plots.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
})

# ---------------------------------------------------------------------------
# Study constants and explicit column mapping
# ---------------------------------------------------------------------------
#
# The analysis scripts use short labels (`SF_Edu`, `SF_Ent`, `LF_Edu`,
# `LF_Ent`) internally. This script keeps those labels as the computational
# source of truth and defines display labels separately so changing plot text
# cannot silently change condition logic.
CONDITION_LEVELS <- c("SF_Edu", "SF_Ent", "LF_Edu", "LF_Ent")
ENGAGEMENT_CONDITION_LEVELS <- c("SF_Edu", "LF_Edu", "SF_Ent", "LF_Ent")
CONTENT_LEVELS <- c("Education", "Entertainment")
LENGTH_LEVELS <- c("Short", "Long")
CONDITION_LABELS <- c(
  SF_Edu = "Short-Form\nEducation",
  SF_Ent = "Short-Form\nEntertainment",
  LF_Edu = "Long-Form\nEducation",
  LF_Ent = "Long-Form\nEntertainment"
)
CONTENT_LABELS <- c(
  Education = "Education",
  Entertainment = "Entertainment"
)
LENGTH_LABELS <- c(
  Short = "Short",
  Long = "Long"
)
LENGTH_PLOTTED_DOMAINS <- c("retention")
PLOT_PALETTE <- c(
  SF_Edu = "#1b4d3e",
  SF_Ent = "#3b7a57",
  LF_Edu = "#8f2d56",
  LF_Ent = "#c04b2c",
  Short = "#1b4d3e",
  Long = "#c04b2c",
  Education = "#1b4d3e",
  Entertainment = "#c04b2c"
)

# Explicitly map every plotted behavioral column to its domain and condition.
# This avoids parsing condition labels from column names, which would be more
# compact but scientifically fragile: a typo or upstream naming change should
# fail through required-column validation rather than be guessed.
BEHAVIOR_COL_MAP <- tibble::tribble(
  ~domain, ~domain_label, ~score_col, ~condition, ~condition_label, ~score_label,
  "engagement", "Engagement", "sf_education_engagement", "SF_Edu", "Short-Form Education", "Engagement score",
  "engagement", "Engagement", "sf_entertainment_engagement", "SF_Ent", "Short-Form Entertainment", "Engagement score",
  "engagement", "Engagement", "lf_education_engagement", "LF_Edu", "Long-Form Education", "Engagement score",
  "engagement", "Engagement", "lf_entertainment_engagement", "LF_Ent", "Long-Form Entertainment", "Engagement score",
  "retention", "Recall / Retention", "diff_short_form_education", "SF_Edu", "Short-Form Education", "Recall improvement (post - pre)",
  "retention", "Recall / Retention", "diff_short_form_entertainment", "SF_Ent", "Short-Form Entertainment", "Recall improvement (post - pre)",
  "retention", "Recall / Retention", "diff_long_form_education", "LF_Edu", "Long-Form Education", "Recall improvement (post - pre)",
  "retention", "Recall / Retention", "diff_long_form_entertainment", "LF_Ent", "Long-Form Entertainment", "Recall improvement (post - pre)"
)

# Locate and source the shared participant-exclusion helper from either the
# working directory or the directory containing this script. This makes CLI
# execution robust from the repo root and from test harnesses that sys.source()
# the script.
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

# Parse simple `--key value` CLI arguments. The defaults intentionally mirror
# the final merged Homer3 + SFV dataset used by the inferential behavioral
# scripts, so a bare `Rscript plot_behavior_score_distributions.R` uses the
# publication analysis input.
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  defaults <- list(
    input_csv = "data/tabular/generated_data/homer3_betas_plus_combined_sfv_data_inner_join.csv",
    exclude_subjects_json = "data/config/excluded_subjects.json",
    out_dir = "data/results/behavior_score_distribution",
    width = 7.0,
    height = 5.0,
    dpi = 300L
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
  parsed$width <- as.numeric(parsed$width)
  parsed$height <- as.numeric(parsed$height)
  parsed$dpi <- as.integer(parsed$dpi)
  parsed
}

# Fail fast when the input is not the expected final behavioral dataset. Silent
# dropping of missing columns would be a direct threat to figure validity.
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

# Convert behavioral score columns to numeric while refusing non-numeric tokens.
# Base R/readr coercion can turn bad strings into NA with only a warning; here
# that is treated as a hard data-integrity error because it would alter plotted
# values or complete-case membership.
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

# Enforce the one-row-per-subject invariant before reshaping. Duplicated subject
# rows would overweight participants in the distribution plots and break parity
# with the LMM scripts.
assert_one_row_per_subject <- function(df, context_label) {
  dup <- df %>%
    count(subject_id, name = "n_rows") %>%
    filter(.data$n_rows > 1) %>%
    arrange(desc(.data$n_rows), .data$subject_id)
  if (nrow(dup) > 0) {
    examples <- dup %>% head(10)
    stop(
      paste0(
        "Duplicate subject_id values detected after normalization in ", context_label, ". ",
        "Expected exactly one row per subject. Example duplicates (subject_id, n_rows): ",
        paste0(examples$subject_id, "=", examples$n_rows, collapse = ", "),
        ". Fix the upstream CSV before plotting."
      )
    )
  }
}

# Read the final merged dataset, normalize subject IDs, validate behavioral
# columns, and apply the central exclusion manifest before any reshaping or
# plotting. This order matches the inferential analysis scripts: exclusions are
# participant-level decisions, not plot-level decisions.
load_behavior_input <- function(input_csv, exclude_subjects_json) {
  df <- read_csv(input_csv, show_col_types = FALSE)
  behavior_cols <- unique(BEHAVIOR_COL_MAP$score_col)
  assert_required_columns(df, c("subject_id", behavior_cols), input_csv)

  df <- df %>%
    mutate(subject_id = normalize_subject_id_exclusion(.data$subject_id, "subject_id"))
  df <- coerce_numeric_strict(df, behavior_cols)
  assert_one_row_per_subject(df, "behavior plotting input")

  excluded <- apply_subject_exclusions(
    df = df,
    subject_col = "subject_id",
    exclude_json_path = exclude_subjects_json,
    context_label = "plot_behavior_score_distributions"
  )
  excluded$data
}

# Convert the wide behavioral table to a tidy subject x domain x condition table.
# The join to BEHAVIOR_COL_MAP is the only place where raw column names become
# study design labels, making mapping errors easy to inspect and test.
reshape_behavior_scores <- function(df) {
  behavior_cols <- unique(BEHAVIOR_COL_MAP$score_col)

  df %>%
    select(subject_id, all_of(behavior_cols)) %>%
    pivot_longer(cols = all_of(behavior_cols), names_to = "score_col", values_to = "score") %>%
    left_join(BEHAVIOR_COL_MAP, by = "score_col", relationship = "many-to-one") %>%
    mutate(
      condition = factor(.data$condition, levels = CONDITION_LEVELS),
      condition_display = factor(
        CONDITION_LABELS[as.character(.data$condition)],
        levels = CONDITION_LABELS[CONDITION_LEVELS]
      )
    ) %>%
    arrange(.data$domain, .data$subject_id, .data$condition)
}

# Apply the same complete-case rule used by the behavioral LMMs, separately for
# engagement and retention. A subject can appear in one domain plot but not the
# other if only one domain has all four condition values.
complete_case_domain_scores <- function(df_long, domain_name) {
  domain_df <- df_long %>%
    filter(.data$domain == domain_name)
  non_missing <- domain_df %>% filter(!is.na(.data$score))

  keep_ids <- non_missing %>%
    group_by(subject_id) %>%
    summarize(n_conditions = n_distinct(condition), .groups = "drop") %>%
    filter(.data$n_conditions == length(CONDITION_LEVELS)) %>%
    pull(subject_id)

  dropped_ids <- setdiff(unique(domain_df$subject_id), keep_ids)
  if (length(dropped_ids) > 0) {
    cat(
      "[plot] ", domain_name,
      ": excluded ", length(dropped_ids),
      " subject(s) for incomplete four-condition behavioral data: ",
      paste(sort(dropped_ids), collapse = ", "),
      "\n",
      sep = ""
    )
  }

  non_missing %>%
    filter(.data$subject_id %in% keep_ids) %>%
    arrange(.data$subject_id, .data$condition)
}

# Collapse the four condition labels to the content factor used by the LMMs.
# Education is the positive sum-coded level (`content_c = +0.5`) in the
# inferential scripts; Entertainment is the negative level (`content_c = -0.5`).
condition_content_level <- function(condition) {
  case_when(
    as.character(condition) %in% c("SF_Edu", "LF_Edu") ~ "Education",
    as.character(condition) %in% c("SF_Ent", "LF_Ent") ~ "Entertainment",
    TRUE ~ NA_character_
  )
}

# Collapse the four condition labels to the length factor used by the LMMs.
# Short is the negative sum-coded level (`length_c = -0.5`) and Long is the
# positive level (`length_c = +0.5`) in the inferential scripts.
condition_length_level <- function(condition) {
  case_when(
    as.character(condition) %in% c("SF_Edu", "SF_Ent") ~ "Short",
    as.character(condition) %in% c("LF_Edu", "LF_Ent") ~ "Long",
    TRUE ~ NA_character_
  )
}

# Build subject-level content marginal means for the content main-effect plots.
# This is deliberately done after four-condition complete-case filtering, so an
# Education or Entertainment value is never computed from a partial subject.
build_content_marginal_scores <- function(domain_df) {
  # Content main-effect visualization mirrors the LMM's content factor:
  # Education and Entertainment are averaged across Short and Long within each
  # complete-case subject before plotting (Searle et al., 1980; see CITATIONS.md).
  out <- domain_df %>%
    mutate(content_level = condition_content_level(.data$condition)) %>%
    group_by(.data$domain, .data$domain_label, .data$subject_id, .data$content_level, .data$score_label) %>%
    summarize(
      score = mean(.data$score),
      n_source_conditions = n_distinct(.data$condition),
      source_score_cols = paste(sort(unique(.data$score_col)), collapse = ";"),
      .groups = "drop"
    )

  bad <- out %>% filter(is.na(.data$content_level) | .data$n_source_conditions != 2)
  if (nrow(bad) > 0) {
    stop(
      paste0(
        "Content marginal means require exactly two source conditions per subject/content level. ",
        "Unexpected rows detected for subject_id values: ",
        paste(head(unique(bad$subject_id), 10), collapse = ", ")
      )
    )
  }

  out %>%
    mutate(
      condition = NA_character_,
      condition_label = NA_character_,
      content_level = factor(.data$content_level, levels = CONTENT_LEVELS),
      content_label = CONTENT_LABELS[as.character(.data$content_level)],
      score_col = "content_marginal_mean"
    ) %>%
    arrange(.data$domain, .data$subject_id, .data$content_level)
}

# Build subject-level length marginal means for retention length main-effect plots.
# Each point is the within-subject mean across content levels for Short or Long
# (Searle et al., 1980; see CITATIONS.md).
build_length_marginal_scores <- function(domain_df) {
  out <- domain_df %>%
    mutate(length_level = condition_length_level(.data$condition)) %>%
    group_by(.data$domain, .data$domain_label, .data$subject_id, .data$length_level, .data$score_label) %>%
    summarize(
      score = mean(.data$score),
      n_source_conditions = n_distinct(.data$condition),
      source_score_cols = paste(sort(unique(.data$score_col)), collapse = ";"),
      .groups = "drop"
    )

  bad <- out %>% filter(is.na(.data$length_level) | .data$n_source_conditions != 2)
  if (nrow(bad) > 0) {
    stop(
      paste0(
        "Length marginal means require exactly two source conditions per subject/length. ",
        "Unexpected rows detected for subject_id values: ",
        paste(head(unique(bad$subject_id), 10), collapse = ", ")
      )
    )
  }

  out %>%
    mutate(
      condition = NA_character_,
      condition_label = NA_character_,
      length_level = factor(.data$length_level, levels = LENGTH_LEVELS),
      length_label = LENGTH_LABELS[as.character(.data$length_level)],
      content_level = NA_character_,
      content_label = NA_character_,
      score_col = "length_marginal_mean"
    ) %>%
    arrange(.data$domain, .data$subject_id, .data$length_level)
}

# Helper used by stat_summary() to show mean +/- 1 SD on top of raw data. The
# summary is descriptive only; inferential quantities remain in the LMM result
# tables.
mean_sdl_1 <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(data.frame(y = NA_real_, ymin = NA_real_, ymax = NA_real_))
  }
  mean_x <- mean(x)
  sd_x <- if (length(x) > 1) stats::sd(x) else 0
  data.frame(y = mean_x, ymin = mean_x - sd_x, ymax = mean_x + sd_x)
}

# Keep output file names stable and shell-friendly even if labels contain spaces
# or punctuation.
sanitize_filename_component <- function(x) {
  gsub("[^A-Za-z0-9_]+", "_", as.character(x))
}

# Internal generic writer for violin/jitter + mean/SD distribution figures.
# Plot-level differences are limited to x-axis grouping and labels.
write_behavior_distribution_plot <- function(
  plot_df,
  domain_name,
  out_dir,
  x_level_col,
  x_display_col,
  x_levels,
  x_labels,
  filename_suffix,
  plot_title,
  plot_subtitle,
  width,
  height,
  dpi
) {
  if (nrow(plot_df) == 0) {
    stop(paste0("No plotted rows remained for domain '", domain_name, "'."))
  }

  if (!all(c(x_level_col, x_display_col) %in% names(plot_df))) {
    stop(paste0("Missing required plotting columns for domain '", domain_name, "'."))
  }

  domain_label <- unique(plot_df$domain_label)
  score_label <- unique(plot_df$score_label)
  if (length(domain_label) != 1 || length(score_label) != 1) {
    stop(paste0("Domain metadata is inconsistent for domain '", domain_name, "'."))
  }

  n_subjects <- n_distinct(plot_df$subject_id)

  # Re-apply factor levels immediately before plotting so any downstream joins
  # or CSV round trips cannot reorder the x-axis alphabetically.
  render_df <- plot_df %>%
    mutate(
      x_level = factor(as.character(.data[[x_level_col]]), levels = x_levels),
      x_display = factor(.data[[x_display_col]], levels = x_labels[x_levels])
    )

  # Violin + jitter + mean/SD overlays follows the existing beta distribution
  # plot style and avoids summary-only bars for publication figures.
  p <- ggplot(render_df, aes(x = .data$x_display, y = .data$score)) +
    geom_violin(
      aes(fill = .data$x_level),
      trim = FALSE,
      alpha = 0.24,
      color = "#595959",
      linewidth = 0.35,
      na.rm = TRUE
    ) +
    geom_point(
      aes(color = .data$x_level),
      position = position_jitter(width = 0.08, height = 0, seed = 1),
      size = 2.1,
      alpha = 0.88,
      na.rm = TRUE
    ) +
    stat_summary(
      fun.data = mean_sdl_1,
      geom = "errorbar",
      width = 0.16,
      linewidth = 0.6,
      color = "#222222",
      na.rm = TRUE
    ) +
    stat_summary(
      fun = mean,
      geom = "point",
      shape = 23,
      size = 3.25,
      stroke = 0.8,
      fill = "white",
      color = "#222222",
      na.rm = TRUE
    ) +
    scale_fill_manual(values = PLOT_PALETTE[x_levels], guide = "none") +
    scale_color_manual(values = PLOT_PALETTE[x_levels], guide = "none") +
    labs(
      title = plot_title,
      subtitle = paste0(plot_subtitle, " | n=", n_subjects),
      x = NULL,
      y = score_label
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(color = "#4b4b4b", size = 10.5),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(face = "bold", color = "#222222"),
      axis.title.y = element_text(face = "bold"),
      plot.margin = margin(12, 14, 12, 14)
    )

  base_name <- paste0(sanitize_filename_component(domain_name), "_", filename_suffix)
  png_path <- file.path(out_dir, paste0(base_name, ".png"))
  suppressMessages(
    ggplot2::ggsave(filename = png_path, plot = p, width = width, height = height, dpi = dpi)
  )
  png_path
}

# Write a four-condition raw-score distribution plot for one behavioral domain.
# These figures answer: "What did subjects score in each of the four study
# conditions after exclusions and complete-case filtering?"
write_behavior_plot <- function(domain_df, domain_name, out_dir, width, height, dpi) {
  # For engagement, place short- and long-form entertainment together to align
  # with a content-major comparison in the same figure. Retention keeps the
  # historical order to preserve existing reporting convention.
  condition_levels <- if (unique(domain_df$domain) == "engagement") {
    ENGAGEMENT_CONDITION_LEVELS
  } else {
    CONDITION_LEVELS
  }

  write_behavior_distribution_plot(
    plot_df = domain_df,
    domain_name = domain_name,
    out_dir = out_dir,
    x_level_col = "condition",
    x_display_col = "condition_display",
    x_levels = condition_levels,
    x_labels = CONDITION_LABELS,
    filename_suffix = "score_distribution",
    plot_title = paste0(unique(domain_df$domain_label), " scores by video condition"),
    plot_subtitle = "Excluded-subject manifest applied | complete-case within domain",
    width = width,
    height = height,
    dpi = dpi
  )
}

# Summarize every plotted data layer for auditability. The `plot_type` column is
# critical: raw-condition, content-marginal, and length-marginal rows represent
# different estimands and should never be pooled accidentally.
summarize_behavior_scores <- function(audit_df) {
  audit_df %>%
    group_by(
      .data$plot_type,
      .data$domain,
      .data$domain_label,
      .data$condition,
      .data$condition_label,
      .data$length_level,
      .data$length_label,
      .data$content_level,
      .data$content_label,
      .data$score_label
    ) %>%
    summarize(
      n_subjects = n_distinct(.data$subject_id),
      n_obs = n(),
      mean = mean(.data$score),
      sd = if (n() > 1) stats::sd(.data$score) else 0,
      min = min(.data$score),
      max = max(.data$score),
      .groups = "drop"
    ) %>%
    arrange(
      .data$plot_type,
      .data$domain,
      .data$condition,
      .data$content_level,
      .data$length_level
    )
}

# Main programmatic entry point used by both the CLI and validation harness.
# It returns the audit data frames in memory and writes the same rows to CSV so
# tests and manuscript reviewers can inspect exactly what was plotted.
run_plotting <- function(
  input_csv,
  exclude_subjects_json,
  out_dir,
  width = 7.0,
  height = 5.0,
  dpi = 300L,
  plot_types = c("raw_condition", "content_marginal", "length_marginal")
) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  supported_plot_types <- c("raw_condition", "content_marginal", "length_marginal")
  requested_plot_types <- match.arg(plot_types, choices = supported_plot_types, several.ok = TRUE)

  df <- load_behavior_input(input_csv, exclude_subjects_json)
  df_long <- reshape_behavior_scores(df)

  # Apply complete-case filtering once per domain, then reuse those exact rows
  # for both raw-condition plots and content-marginal plots. This guarantees the
  # two views are derived from the same subject sets.
  domain_dfs <- lapply(unique(BEHAVIOR_COL_MAP$domain), function(domain_name) {
    complete_case_domain_scores(df_long, domain_name)
  })

  # Audit raw rows before marginalization. These rows preserve original score
  # columns so the plotted points can be traced back to the final merged dataset.
  raw_audit_df <- bind_rows(domain_dfs) %>%
    transmute(
      plot_type = "raw_condition",
      domain = .data$domain,
      domain_label = .data$domain_label,
      subject_id = .data$subject_id,
      condition = as.character(.data$condition),
      condition_label = .data$condition_label,
      content_level = condition_content_level(.data$condition),
      content_label = CONTENT_LABELS[condition_content_level(.data$condition)],
      score_col = .data$score_col,
      source_score_cols = .data$score_col,
      n_source_conditions = 1L,
      length_level = NA_character_,
      length_label = NA_character_,
      score_label = .data$score_label,
      score = .data$score,
      included_in_plot = TRUE
    ) %>%
    arrange(.data$domain, .data$subject_id, factor(.data$condition, levels = CONDITION_LEVELS))

  # Audit content marginal rows separately. `source_score_cols` records the two
  # raw condition columns used for each marginal point.
  content_audit_df <- bind_rows(lapply(domain_dfs, build_content_marginal_scores)) %>%
    transmute(
      plot_type = "content_marginal",
      domain = .data$domain,
      domain_label = .data$domain_label,
      subject_id = .data$subject_id,
      condition = .data$condition,
      condition_label = .data$condition_label,
      content_level = as.character(.data$content_level),
      content_label = .data$content_label,
      length_level = NA_character_,
      length_label = NA_character_,
      score_col = .data$score_col,
      source_score_cols = .data$source_score_cols,
      n_source_conditions = .data$n_source_conditions,
      score_label = .data$score_label,
      score = .data$score,
      included_in_plot = TRUE
    ) %>%
    arrange(.data$domain, .data$subject_id, factor(.data$content_level, levels = CONTENT_LEVELS))

  # Retention-specific length marginals are used for the recall interpretation of
  # short-vs-long differences (Searle et al., 1980; see CITATIONS.md).
  length_audit_df <- bind_rows(lapply(domain_dfs, build_length_marginal_scores)) %>%
    filter(.data$domain %in% LENGTH_PLOTTED_DOMAINS) %>%
    transmute(
      plot_type = "length_marginal",
      domain = .data$domain,
      domain_label = .data$domain_label,
      subject_id = .data$subject_id,
      condition = .data$condition,
      condition_label = .data$condition_label,
      length_level = as.character(.data$length_level),
      length_label = as.character(.data$length_label),
      content_level = as.character(.data$content_level),
      content_label = .data$content_label,
      score_col = .data$score_col,
      source_score_cols = .data$source_score_cols,
      n_source_conditions = .data$n_source_conditions,
      score_label = .data$score_label,
      score = .data$score,
      included_in_plot = TRUE
    ) %>%
    arrange(.data$domain, .data$subject_id, factor(.data$length_level, levels = LENGTH_LEVELS))

  # Use explicit factor ordering in the audit CSV so visual x-axis order and
  # machine-readable row order tell the same story.
  audit_df <- bind_rows(raw_audit_df, content_audit_df, length_audit_df) %>%
    arrange(
      .data$plot_type,
      .data$domain,
      .data$subject_id,
      factor(.data$length_level, levels = LENGTH_LEVELS),
      factor(.data$condition, levels = CONDITION_LEVELS),
      factor(.data$content_level, levels = CONTENT_LEVELS)
    )

  audit_csv_path <- file.path(out_dir, "plotted_behavior_scores.csv")
  summary_csv_path <- file.path(out_dir, "behavior_score_distribution_summary.csv")
  write_csv(audit_df, audit_csv_path, na = "")
  summary_df <- summarize_behavior_scores(audit_df)
  write_csv(summary_df, summary_csv_path, na = "")

  figure_paths <- list()
  for (domain_name in unique(BEHAVIOR_COL_MAP$domain)) {
    if ("raw_condition" %in% requested_plot_types) {
      # Four-condition descriptive figure.
      domain_df <- audit_df %>%
        filter(.data$plot_type == "raw_condition", .data$domain == domain_name) %>%
        mutate(
          condition = factor(.data$condition, levels = CONDITION_LEVELS),
          condition_display = factor(
            CONDITION_LABELS[as.character(.data$condition)],
            levels = CONDITION_LABELS[CONDITION_LEVELS]
          )
        )
      figure_paths[[domain_name]] <- write_behavior_plot(
        domain_df = domain_df,
        domain_name = domain_name,
        out_dir = out_dir,
        width = width,
        height = height,
        dpi = dpi
      )
    }

    if ("content_marginal" %in% requested_plot_types) {
      # Content marginal display: each point is a subject-level mean across video
      # length. These are raw descriptive values, not age-adjusted model EMMs.
      content_df <- audit_df %>%
        filter(.data$plot_type == "content_marginal", .data$domain == domain_name)
      figure_paths[[paste0(domain_name, "_content_marginal")]] <- write_behavior_distribution_plot(
        plot_df = content_df,
        domain_name = domain_name,
        out_dir = out_dir,
        x_level_col = "content_level",
        x_display_col = "content_label",
        x_levels = CONTENT_LEVELS,
        x_labels = CONTENT_LABELS,
        filename_suffix = "content_marginal_score_distribution",
        plot_title = paste0(unique(content_df$domain_label), " content marginal scores"),
        plot_subtitle = "Subject-level content marginal means averaged across video length",
        width = width,
        height = height,
        dpi = dpi
      )
    }

    if ("length_marginal" %in% requested_plot_types && domain_name %in% LENGTH_PLOTTED_DOMAINS) {
      # Retention length marginal display: each point is a subject-level mean
      # across content. These are raw descriptive values, not age-adjusted model
      # EMMs.
      length_df <- audit_df %>%
        filter(.data$plot_type == "length_marginal", .data$domain == domain_name) %>%
        mutate(
          length_level = factor(.data$length_level, levels = LENGTH_LEVELS),
          length_label = factor(
            as.character(.data$length_level),
            levels = LENGTH_LEVELS,
            labels = as.character(LENGTH_LABELS[LENGTH_LEVELS])
          )
        )
      figure_paths[[paste0(domain_name, "_length_marginal")]] <- write_behavior_distribution_plot(
        plot_df = length_df,
        domain_name = domain_name,
        out_dir = out_dir,
        x_level_col = "length_level",
        x_display_col = "length_label",
        x_levels = LENGTH_LEVELS,
        x_labels = LENGTH_LABELS,
        filename_suffix = "length_marginal_score_distribution",
        plot_title = paste0(unique(length_df$domain_label), " length marginal scores"),
        plot_subtitle = "Subject-level length marginal means averaged across content",
        width = width,
        height = height,
        dpi = dpi
      )
    }
  }

  list(
    figure_paths = figure_paths,
    audit_df = audit_df,
    summary_df = summary_df,
    audit_csv_path = audit_csv_path,
    summary_csv_path = summary_csv_path
  )
}

# CLI wrapper: parse arguments, run the pipeline, and print the written artifact
# locations. The heavy lifting stays in run_plotting() so tests can exercise the
# same code path without shelling out.
main <- function() {
  args <- parse_args()
  outputs <- run_plotting(
    input_csv = args$input_csv,
    exclude_subjects_json = args$exclude_subjects_json,
    out_dir = args$out_dir,
    width = args$width,
    height = args$height,
    dpi = args$dpi
  )

  cat("[plot] input CSV:", args$input_csv, "\n")
  cat("[plot] output directory:", args$out_dir, "\n")
  cat("[plot] audit CSV:", outputs$audit_csv_path, "\n")
  cat("[plot] summary CSV:", outputs$summary_csv_path, "\n")
  cat("[plot] figure files written:", sum(vapply(outputs$figure_paths, length, integer(1))), "\n")
}

if (sys.nframe() == 0) {
  main()
}
