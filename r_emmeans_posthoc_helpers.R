standardize_emmeans_pairwise_output <- function(pw, context_label = NULL) {
  stat_type <- if ("t.ratio" %in% names(pw)) "t" else if ("z.ratio" %in% names(pw)) "z" else NA_character_
  if (is.na(stat_type)) {
    stop("Expected 't.ratio' or 'z.ratio' in emmeans pairwise output.")
  }

  if (stat_type == "z") {
    pw[["t_display"]] <- "t can't be calculated"
    msg <- paste0(
      "[warn] emmeans returned z.ratio (asymptotic). Writing \"t can't be calculated\" in posthoc t column",
      if (!is.null(context_label) && nzchar(context_label)) paste0(" for ", context_label) else "",
      ".\n"
    )
    cat(msg)
  } else {
    pw[["t_display"]] <- as.character(pw[["t.ratio"]])
  }
  if (!("df" %in% names(pw))) {
    pw[["df"]] <- NA_real_
  }

  pw %>%
    mutate(
      stat_type = stat_type,
      condition_a = stringr::str_trim(stringr::str_split_fixed(.data$contrast, " - ", 2)[, 1]),
      condition_b = stringr::str_trim(stringr::str_split_fixed(.data$contrast, " - ", 2)[, 2])
    ) %>%
    select(condition_a, condition_b, stat_type, estimate, SE, df, t_display, p.value) %>%
    rename(mean_diff = estimate, se = SE, t = t_display, p_unc = p.value)
}
