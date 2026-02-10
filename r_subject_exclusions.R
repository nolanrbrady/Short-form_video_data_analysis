#!/usr/bin/env Rscript

# Shared participant-exclusion helpers for inferential analysis scripts.
#
# Reproducibility guardrail:
# Maintain one exclusion manifest consumed by all downstream analyses so that
# subject removal rules do not diverge across endpoints (Sandve et al., 2013;
# see CITATIONS.md).

suppressPackageStartupMessages({
  library(jsonlite)
})

normalize_subject_id_exclusion <- function(x, column_name) {
  x_chr <- as.character(x)
  extracted <- stringr::str_extract(x_chr, "\\d+")
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

read_excluded_subject_ids_json <- function(path) {
  if (!file.exists(path)) {
    stop(
      paste0(
        "Exclusion file not found: ", path,
        ". Provide --exclude_subjects_json with a valid JSON file."
      )
    )
  }

  payload <- tryCatch(
    jsonlite::fromJSON(path, simplifyVector = TRUE),
    error = function(e) {
      stop(paste0("Failed to parse exclusion JSON at ", path, ": ", e$message))
    }
  )

  if (is.null(payload)) {
    payload <- character()
  }
  if (is.list(payload) && length(payload) == 0) {
    payload <- character()
  }
  if (!is.atomic(payload) || is.list(payload)) {
    stop(
      paste0(
        "Exclusion file must be a top-level JSON array of subject IDs. File: ",
        path
      )
    )
  }

  ids_raw <- as.character(payload)
  if (length(ids_raw) == 0) {
    return(list(raw_ids = character(), norm_ids = integer()))
  }
  if (any(is.na(ids_raw) | trimws(ids_raw) == "")) {
    stop(paste0("Exclusion file contains empty or NA subject IDs: ", path))
  }

  ids_norm <- normalize_subject_id_exclusion(ids_raw, paste0("exclude_subjects_json(", path, ")"))
  dup_norm <- sort(unique(ids_norm[duplicated(ids_norm)]))
  if (length(dup_norm) > 0) {
    stop(
      paste0(
        "Duplicate exclusion IDs detected after normalization in ",
        path,
        ". Duplicate normalized IDs: ",
        paste(dup_norm, collapse = ", ")
      )
    )
  }

  list(raw_ids = ids_raw, norm_ids = ids_norm)
}

apply_subject_exclusions <- function(df, subject_col, exclude_json_path, context_label = "analysis") {
  if (!(subject_col %in% names(df))) {
    stop(paste0("Column '", subject_col, "' not found in data passed to apply_subject_exclusions()."))
  }

  out <- df
  out[[subject_col]] <- normalize_subject_id_exclusion(out[[subject_col]], subject_col)

  excl <- read_excluded_subject_ids_json(exclude_json_path)
  n_before <- length(unique(out[[subject_col]]))

  excl_unique <- sort(unique(excl$norm_ids))
  present_ids <- sort(unique(out[[subject_col]]))
  matched_ids <- intersect(present_ids, excl_unique)
  missing_ids <- setdiff(excl_unique, present_ids)

  keep <- !(out[[subject_col]] %in% excl_unique)
  out <- out[keep, , drop = FALSE]
  n_after <- length(unique(out[[subject_col]]))
  n_removed <- n_before - n_after

  missing_raw <- unique(excl$raw_ids[excl$norm_ids %in% missing_ids])

  cat(
    "[exclude] ", context_label,
    ": listed=", length(excl_unique),
    ", matched=", length(matched_ids),
    ", removed=", n_removed,
    ", remaining=", n_after,
    "\n",
    sep = ""
  )
  if (length(missing_ids) > 0) {
    cat(
      "[warn] ", context_label,
      ": exclusion IDs not present in dataset: ",
      paste(missing_raw, collapse = ", "),
      "\n",
      sep = ""
    )
  }

  list(
    data = out,
    report = list(
      listed = length(excl_unique),
      matched = length(matched_ids),
      removed = n_removed,
      remaining = n_after,
      missing_norm = missing_ids,
      missing_raw = missing_raw
    )
  )
}
