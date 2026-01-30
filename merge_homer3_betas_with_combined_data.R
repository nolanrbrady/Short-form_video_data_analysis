#!/usr/bin/env Rscript

# Merge Homer3 GLM betas (wide) with tabular combined dataset (demographics/behavior).
#
# This script performs an INNER JOIN between:
#   - data/tabular/homer3_glm_betas_wide.csv (ID column: Subject, e.g. sub_0001)
#   - data/tabular/combined_sfv_data.csv     (ID column: subject_id, sometimes zero-padded)
#
# Key behavior:
#   - Extracts the numeric subject id from both ID columns (handles padding like 0017 vs 17,
#     and Homer-style IDs like sub_0017).
#   - Produces a single merged row per subject containing all columns from both inputs.
#   - Does NOT impute missingness. Downstream modeling must explicitly treat 0/NaN as
#     pruned/missing channels per repo policy.
#
# Scientific integrity note:
#   - homer3_glm_betas_wide.csv can contain both 0 and NaN as stand-ins for pruned channels.
#     Do not interpret these as true zero activation; handle as missing downstream.
#
# Usage:
#   Rscript merge_homer3_betas_with_combined_data.R \
#     --homer_csv data/tabular/homer3_glm_betas_wide.csv \
#     --combined_csv data/tabular/combined_sfv_data.csv \
#     --out_csv data/results/homer3_betas_plus_combined_sfv_data_inner_join.csv

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  defaults <- list(
    homer_csv = "data/tabular/homer3_glm_betas_wide.csv",
    combined_csv = "data/tabular/combined_sfv_data.csv",
    out_csv = "data/results/homer3_betas_plus_combined_sfv_data_inner_join.csv"
  )
  if (length(args) == 0) return(defaults)

  # Minimal arg parser: --key value pairs
  if (length(args) %% 2 != 0) stop("Expected --key value argument pairs.")
  parsed <- defaults
  for (i in seq(1, length(args), by = 2)) {
    key <- args[[i]]
    val <- args[[i + 1]]
    if (!startsWith(key, "--")) stop(paste0("Invalid argument: ", key))
    key <- substring(key, 3)
    parsed[[key]] <- val
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

coerce_beta_columns_to_numeric <- function(homer) {
  beta_cols <- names(homer)[str_detect(names(homer), "^S\\d+_D\\d+_Cond\\d{2}_(HbO|HbR)$")]
  if (length(beta_cols) == 0) stop("No beta columns matched pattern like 'S01_D01_Cond01_HbO'.")

  allowed_missing <- c("", "NA", "NaN", "nan", "NAN", "NULL", "null")
  for (col in beta_cols) {
    raw <- as.character(homer[[col]])
    raw_trim <- str_trim(raw)
    num <- suppressWarnings(as.numeric(raw_trim))
    bad <- !is.na(raw_trim) & is.na(num) & !(raw_trim %in% allowed_missing)
    if (any(bad)) {
      examples <- unique(raw_trim[bad])
      examples <- examples[!is.na(examples)]
      examples <- head(examples, 10)
      stop(
        paste0(
          "Non-numeric values found in Homer3 beta column '", col, "'. ",
          "Examples: ", paste(examples, collapse = ", "),
          ". Fix the upstream CSV export before merging."
        )
      )
    }
    homer[[col]] <- num
  }
  homer
}

main <- function() {
  args <- parse_args()

  homer <- read_csv(args$homer_csv, show_col_types = FALSE)
  combined <- read_csv(args$combined_csv, show_col_types = FALSE)

  if (!("Subject" %in% names(homer))) stop(paste0("Expected column 'Subject' in ", args$homer_csv))
  if (!("subject_id" %in% names(combined))) stop(paste0("Expected column 'subject_id' in ", args$combined_csv))

  homer <- coerce_beta_columns_to_numeric(homer)

  homer <- homer %>%
    mutate(subject_id = normalize_subject_id(.data$Subject, "Subject")) %>%
    rename(homer_subject = Subject)

  combined <- combined %>%
    mutate(subject_id = normalize_subject_id(.data$subject_id, "subject_id"))

  # Fail hard if either input contains duplicate subject IDs. The merge spec assumes one row per subject.
  for (nm in c("combined", "homer")) {
    df <- if (nm == "combined") combined else homer
    dup <- df %>%
      count(subject_id, name = "n_rows") %>%
      filter(.data$n_rows > 1) %>%
      arrange(desc(.data$n_rows), .data$subject_id)
    if (nrow(dup) > 0) {
      examples <- dup %>% head(10)
      stop(
        paste0(
          "Duplicate subject_id values detected in ", nm, " dataset after normalization. ",
          "Expected exactly one row per subject for an inner one-to-one join. ",
          "Example duplicates (subject_id, n_rows): ",
          paste0(examples$subject_id, "=", examples$n_rows, collapse = ", "),
          ". Fix the upstream CSV before running the merge."
        )
      )
    }
  }

  merged <- combined %>%
    inner_join(homer, by = "subject_id", relationship = "one-to-one")

  combined_ids <- unique(combined$subject_id)
  homer_ids <- unique(homer$subject_id)
  merged_ids <- unique(merged$subject_id)

  cat("[merge] combined subjects:", length(combined_ids), "\n")
  cat("[merge] homer subjects:   ", length(homer_ids), "\n")
  cat("[merge] inner-joined:     ", length(merged_ids), "\n")
  cat("[merge] combined-only (dropped):", length(setdiff(combined_ids, merged_ids)), "\n")
  cat("[merge] homer-only (dropped):   ", length(setdiff(homer_ids, merged_ids)), "\n")

  out_dir <- dirname(args$out_csv)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  write_csv(merged, args$out_csv)
  cat("[merge] wrote:", args$out_csv, "\n")
}

main()
