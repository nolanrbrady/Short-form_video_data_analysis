#!/usr/bin/env Rscript

# Plot simple subject-level beta-value distributions for FDR-significant
# Format/Content/Interaction hits from the channelwise and ROI LMM outputs.
#
# Inputs
#   - data/tabular/generated_data/homer3_betas_plus_combined_sfv_data_inner_join.csv
#   - data/config/roi_definition.json
#   - data/config/excluded_subjects.json
#   - data/results/format_content_lmm_main_effects_tidy_r.csv
#   - data/results/format_content_lmm_roi_main_effects_tidy_r.csv
#
# Behavior
#   - Applies the shared subject-exclusion manifest so plotted subjects align
#     with the inferential analyses (Sandve et al., 2013; see CITATIONS.md).
#   - Reconstructs channel-level and ROI-level beta rows from the merged wide
#     beta table used by the LMM scripts.
#   - For main effects, plots subject-level marginal means across the orthogonal
#     factor so the displayed points map to the tested main-effect contrast
#     (Searle et al., 1980; see CITATIONS.md).
#   - For interaction effects, plots the four raw condition distributions rather
#     than collapsing across the interacting factor.
#   - Shows raw point distributions directly instead of summary-only bars so
#     readers can inspect spread, sample size, and unusual values
#     (Weissgerber et al., 2015; see CITATIONS.md).
#
# Missingness / pruned channels
#   - Pruned channels must remain explicit missing values.
#   - If the merged input contains literal beta-value zeros, the script fails
#     hard because this project treats zero placeholders as pruned/missing
#     stand-ins rather than true zero activation for downstream neural analyses
#     (Yucel et al., 2021; see CITATIONS.md).
#
# Citations (see CITATIONS.md)
#   - Yucel et al. (2021): transparent fNIRS reporting and explicit missingness handling.
#   - Poldrack (2007): ROI signal extraction by averaging pre-specified channels.
#   - Searle et al. (1980): equal-weight marginal means for main-effect displays.
#   - Weissgerber et al. (2015): raw-data visualization over summary-only plots.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(jsonlite)
  library(ggplot2)
})

CONDITION_LEVELS <- c("SF_Edu", "SF_Ent", "LF_Ent", "LF_Edu")
FORMAT_LEVELS <- c("Short", "Long")
CONTENT_LEVELS <- c("Education", "Entertainment")
SUPPORTED_EFFECTS <- c("format", "content", "interaction")
PLOT_PALETTE <- c(
  Short = "#1b4d3e",
  Long = "#c04b2c",
  Education = "#1b4d3e",
  Entertainment = "#c04b2c",
  SF_Edu = "#1b4d3e",
  SF_Ent = "#3b7a57",
  LF_Ent = "#c04b2c",
  LF_Edu = "#8f2d56"
)

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

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  defaults <- list(
    input_csv = "data/tabular/generated_data/homer3_betas_plus_combined_sfv_data_inner_join.csv",
    roi_json = "data/config/roi_definition.json",
    exclude_subjects_json = "data/config/excluded_subjects.json",
    channel_results_tidy_csv = "data/results/format_content_lmm_main_effects_tidy_r.csv",
    roi_results_tidy_csv = "data/results/format_content_lmm_roi_main_effects_tidy_r.csv",
    alpha = 0.05,
    out_dir = "data/results/beta_value_distribution"
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
          ". Use strict JSON syntax. Error: ",
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

  roi_map %>% arrange(.data$roi, .data$channel)
}

assert_no_zero_placeholders <- function(df, beta_cols, input_csv) {
  zero_mask <- df[beta_cols] == 0
  zero_mask[is.na(zero_mask)] <- FALSE
  if (!any(zero_mask)) return(invisible(NULL))

  idx <- which(as.matrix(zero_mask), arr.ind = TRUE)
  examples <- apply(head(idx, 10), 1, function(row_idx) {
    subj <- df$subject_id[[row_idx[[1]]]]
    col_name <- beta_cols[[row_idx[[2]]]]
    paste0("subject_id=", subj, " column=", col_name)
  })
  stop(
    paste0(
      "Merged beta input contains literal 0 values in beta columns. ",
      "This plotting script treats 0 placeholders as invalid stand-ins for pruned channels and refuses to guess whether they are true zeros. ",
      "Resolve the upstream export so pruned channels are NA before plotting. ",
      "Examples: ", paste(examples, collapse = "; "),
      ". Input: ", input_csv
    )
  )
}

load_merged_input <- function(input_csv) {
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
  df <- coerce_numeric_strict(df, c("age", beta_cols))
  assert_no_zero_placeholders(df, beta_cols, input_csv)

  dup <- df %>%
    count(subject_id, name = "n_rows") %>%
    filter(.data$n_rows > 1) %>%
    arrange(desc(.data$n_rows), .data$subject_id)
  if (nrow(dup) > 0) {
    examples <- dup %>% head(10)
    stop(
      paste0(
        "Duplicate subject_id values detected after normalization in merged input. ",
        "Expected exactly one row per subject. Example duplicates (subject_id, n_rows): ",
        paste0(examples$subject_id, "=", examples$n_rows, collapse = ", "),
        ". Fix the upstream CSV before plotting."
      )
    )
  }

  df
}

reshape_to_long <- function(df_merged) {
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
      bad_conds <- setdiff(unique(.$cond), c("01", "02", "03", "04"))
      if (length(bad_conds) > 0) {
        stop(paste0("Unexpected condition codes in beta columns: ", paste(sort(bad_conds), collapse = ", ")))
      }
      .
    } %>%
    mutate(
      channel = normalize_channel_id(.data$channel, "beta column names"),
      beta = as.numeric(.data$beta),
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
      )
    ) %>%
    filter(!is.na(.data$condition), !is.na(.data$chrom), !is.na(.data$channel))
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
}

aggregate_to_roi <- function(df_long, roi_map) {
  # ROI beta is the arithmetic mean across available member channels, matching
  # the ROI inferential script (Poldrack, 2007; see CITATIONS.md).
  df_long %>%
    inner_join(roi_map, by = "channel", relationship = "many-to-one") %>%
    group_by(subject_id, age, roi, chrom, condition, format, content) %>%
    summarize(
      n_channels_in_roi = n_distinct(channel),
      n_channels_nonmissing = sum(!is.na(beta)),
      beta = if (all(is.na(beta))) NA_real_ else mean(beta, na.rm = TRUE),
      .groups = "drop"
    )
}

complete_case_subjects <- function(sub) {
  non_missing <- sub %>% filter(!is.na(.data$beta))
  keep_ids <- non_missing %>%
    group_by(subject_id) %>%
    summarize(n_cond = n_distinct(condition), .groups = "drop") %>%
    filter(.data$n_cond == 4) %>%
    pull(subject_id)
  non_missing %>% filter(.data$subject_id %in% keep_ids)
}

empty_audit_df <- function() {
  tibble::tibble(
    analysis_level = character(),
    unit_id = character(),
    chrom = character(),
    effect = character(),
    plot_mode = character(),
    subject_id = integer(),
    condition = character(),
    format = character(),
    content = character(),
    beta_raw = numeric(),
    beta_plot_value = numeric(),
    included_in_plot = logical()
  )
}

load_significant_hits <- function(results_csv, unit_col, alpha) {
  df <- read_csv(results_csv, show_col_types = FALSE)
  assert_required_columns(df, c(unit_col, "chrom", "effect", "p_fdr"), results_csv)
  bad_effects <- setdiff(unique(df$effect), SUPPORTED_EFFECTS)
  if (length(bad_effects) > 0) {
    stop(
      paste0(
        "Unsupported effect labels in ", results_csv, ": ",
        paste(sort(bad_effects), collapse = ", ")
      )
    )
  }

  df %>%
    filter(is.finite(.data$p_fdr), .data$p_fdr < alpha) %>%
    arrange(.data$p_fdr, .data[[unit_col]], .data$chrom, .data$effect)
}

sanitize_filename_component <- function(x) {
  gsub("[^A-Za-z0-9_]+", "_", as.character(x))
}

effect_group_column <- function(effect) {
  if (identical(effect, "format")) return("format")
  if (identical(effect, "content")) return("content")
  stop(paste0("Main-effect plotting does not support effect '", effect, "'."))
}

effect_group_levels <- function(effect) {
  if (identical(effect, "format")) return(FORMAT_LEVELS)
  if (identical(effect, "content")) return(CONTENT_LEVELS)
  if (identical(effect, "interaction")) return(CONDITION_LEVELS)
  stop(paste0("Unsupported effect: ", effect))
}

build_main_effect_plot_bundle <- function(sub_complete, analysis_level, unit_id, chrom, effect) {
  group_col <- effect_group_column(effect)
  plot_points <- sub_complete %>%
    mutate(plot_level = .data[[group_col]]) %>%
    group_by(subject_id, plot_level) %>%
    summarize(beta_plot_value = mean(beta), .groups = "drop") %>%
    mutate(plot_level = factor(.data$plot_level, levels = effect_group_levels(effect)))

  audit_df <- sub_complete %>%
    mutate(plot_level = .data[[group_col]]) %>%
    left_join(plot_points, by = c("subject_id", "plot_level")) %>%
    transmute(
      analysis_level = analysis_level,
      unit_id = unit_id,
      chrom = chrom,
      effect = effect,
      plot_mode = "main_effect_marginal",
      subject_id = .data$subject_id,
      condition = .data$condition,
      format = .data$format,
      content = .data$content,
      beta_raw = .data$beta,
      beta_plot_value = .data$beta_plot_value,
      included_in_plot = TRUE
    )

  list(
    point_df = plot_points,
    audit_df = audit_df,
    plot_mode = "main_effect_marginal"
  )
}

build_interaction_plot_bundle <- function(sub_complete, analysis_level, unit_id, chrom, effect) {
  plot_points <- sub_complete %>%
    transmute(
      subject_id = .data$subject_id,
      plot_level = factor(.data$condition, levels = CONDITION_LEVELS),
      beta_plot_value = .data$beta
    )

  audit_df <- sub_complete %>%
    transmute(
      analysis_level = analysis_level,
      unit_id = unit_id,
      chrom = chrom,
      effect = effect,
      plot_mode = "interaction_conditions",
      subject_id = .data$subject_id,
      condition = .data$condition,
      format = .data$format,
      content = .data$content,
      beta_raw = .data$beta,
      beta_plot_value = .data$beta,
      included_in_plot = TRUE
    )

  list(
    point_df = plot_points,
    audit_df = audit_df,
    plot_mode = "interaction_conditions"
  )
}

build_plot_bundle <- function(sub_complete, analysis_level, unit_id, chrom, effect) {
  if (identical(effect, "interaction")) {
    return(build_interaction_plot_bundle(sub_complete, analysis_level, unit_id, chrom, effect))
  }
  build_main_effect_plot_bundle(sub_complete, analysis_level, unit_id, chrom, effect)
}

plot_subtitle <- function(effect, plot_mode, n_subjects) {
  if (identical(plot_mode, "interaction_conditions")) {
    return(paste0("FDR-significant interaction | raw condition betas | n=", n_subjects))
  }
  paste0("FDR-significant ", effect, " main effect | subject-level marginal means | n=", n_subjects)
}

write_distribution_plot <- function(point_df, analysis_level, unit_id, chrom, effect, out_dir) {
  n_subjects <- n_distinct(point_df$subject_id)
  plot_mode <- if (identical(effect, "interaction")) "interaction_conditions" else "main_effect_marginal"
  x_levels <- effect_group_levels(effect)
  x_label <- if (identical(effect, "interaction")) "Condition" else str_to_title(effect)
  title_prefix <- if (identical(analysis_level, "channel")) "Channel" else "ROI"
  title_text <- paste(title_prefix, unit_id, chrom, "beta distribution")
  subtitle_text <- plot_subtitle(effect, plot_mode, n_subjects)
  width_in <- if (identical(effect, "interaction")) 7.5 else 6.5

  plot_df <- point_df %>%
    mutate(plot_level = factor(.data$plot_level, levels = x_levels))

  p <- ggplot(plot_df, aes(x = .data$plot_level, y = .data$beta_plot_value)) +
    geom_line(
      aes(group = .data$subject_id),
      linewidth = 0.35,
      alpha = 0.35,
      color = "#8c8c8c"
    ) +
    geom_point(
      aes(color = .data$plot_level),
      position = position_jitter(width = 0.08, height = 0, seed = 1),
      size = 2.1,
      alpha = 0.88
    ) +
    scale_color_manual(values = PLOT_PALETTE[x_levels], guide = "none") +
    labs(
      title = title_text,
      subtitle = subtitle_text,
      x = x_label,
      y = "Beta value"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(face = "bold")
    )

  file_name <- paste0(
    analysis_level, "_",
    sanitize_filename_component(unit_id), "_",
    sanitize_filename_component(chrom), "_",
    sanitize_filename_component(effect), "_beta_distribution.png"
  )
  out_path <- file.path(out_dir, file_name)
  suppressMessages(
    ggplot2::ggsave(filename = out_path, plot = p, width = width_in, height = 5, dpi = 300)
  )
  out_path
}

plot_significant_hit <- function(df_source, analysis_level, unit_col, hit_row, out_dir) {
  unit_id <- hit_row[[unit_col]][[1]]
  chrom_name <- hit_row[["chrom"]][[1]]
  effect_name <- hit_row[["effect"]][[1]]

  sub <- df_source %>%
    filter(.data[[unit_col]] == .env$unit_id, .data$chrom == .env$chrom_name)
  sub_complete <- complete_case_subjects(sub)
  n_subjects <- n_distinct(sub_complete$subject_id)
  if (n_subjects == 0) {
    stop(
      paste0(
        "No complete-case subjects remained for plotted hit: ",
        analysis_level, " ", unit_id, " ", chrom_name, " ", effect_name,
        ". This indicates a mismatch between the result table and the plotting input."
      )
    )
  }

  bundle <- build_plot_bundle(
    sub_complete = sub_complete,
    analysis_level = analysis_level,
    unit_id = unit_id,
    chrom = chrom_name,
    effect = effect_name
  )
  figure_path <- write_distribution_plot(
    point_df = bundle$point_df,
    analysis_level = analysis_level,
    unit_id = unit_id,
    chrom = chrom_name,
    effect = effect_name,
    out_dir = out_dir
  )

  list(
    figure_path = figure_path,
    audit_df = bundle$audit_df
  )
}

run_plotting <- function(
  input_csv,
  roi_json,
  exclude_subjects_json,
  channel_results_tidy_csv,
  roi_results_tidy_csv,
  alpha,
  out_dir
) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  df_merged <- load_merged_input(input_csv)
  excluded <- apply_subject_exclusions(
    df = df_merged,
    subject_col = "subject_id",
    exclude_json_path = exclude_subjects_json,
    context_label = "plot_significant_beta_value_distribution"
  )
  df_merged <- excluded$data
  df_long <- reshape_to_long(df_merged)

  roi_map <- load_roi_definition(roi_json)
  validate_roi_channels(df_long, roi_map)
  df_roi <- aggregate_to_roi(df_long, roi_map)

  sig_channel_hits <- load_significant_hits(channel_results_tidy_csv, "channel", alpha)
  sig_roi_hits <- load_significant_hits(roi_results_tidy_csv, "roi", alpha)

  audit_rows <- list()
  figure_paths <- character()

  if (nrow(sig_channel_hits) > 0) {
    for (i in seq_len(nrow(sig_channel_hits))) {
      hit <- sig_channel_hits[i, , drop = FALSE]
      plotted <- plot_significant_hit(
        df_source = df_long,
        analysis_level = "channel",
        unit_col = "channel",
        hit_row = hit,
        out_dir = out_dir
      )
      figure_paths <- c(figure_paths, plotted$figure_path)
      audit_rows[[length(audit_rows) + 1]] <- plotted$audit_df
    }
  }

  if (nrow(sig_roi_hits) > 0) {
    for (i in seq_len(nrow(sig_roi_hits))) {
      hit <- sig_roi_hits[i, , drop = FALSE]
      plotted <- plot_significant_hit(
        df_source = df_roi,
        analysis_level = "roi",
        unit_col = "roi",
        hit_row = hit,
        out_dir = out_dir
      )
      figure_paths <- c(figure_paths, plotted$figure_path)
      audit_rows[[length(audit_rows) + 1]] <- plotted$audit_df
    }
  }

  audit_df <- if (length(audit_rows) > 0) {
    bind_rows(audit_rows) %>%
      arrange(.data$analysis_level, .data$unit_id, .data$chrom, .data$effect, .data$subject_id, .data$condition)
  } else {
    empty_audit_df()
  }

  audit_csv_path <- file.path(out_dir, "plotted_beta_values.csv")
  write_csv(audit_df, audit_csv_path, na = "")

  list(
    figure_paths = figure_paths,
    audit_df = audit_df,
    audit_csv_path = audit_csv_path,
    significant_channel_hits = sig_channel_hits,
    significant_roi_hits = sig_roi_hits
  )
}

main <- function() {
  args <- parse_args()
  outputs <- run_plotting(
    input_csv = args$input_csv,
    roi_json = args$roi_json,
    exclude_subjects_json = args$exclude_subjects_json,
    channel_results_tidy_csv = args$channel_results_tidy_csv,
    roi_results_tidy_csv = args$roi_results_tidy_csv,
    alpha = args$alpha,
    out_dir = args$out_dir
  )

  cat("[plot] output directory:", args$out_dir, "\n")
  cat("[plot] channelwise significant hits:", nrow(outputs$significant_channel_hits), "\n")
  cat("[plot] ROI significant hits:", nrow(outputs$significant_roi_hits), "\n")
  cat("[plot] figures written:", length(outputs$figure_paths), "\n")
  cat("[plot] audit CSV:", outputs$audit_csv_path, "\n")
}

if (sys.nframe() == 0) {
  main()
}
