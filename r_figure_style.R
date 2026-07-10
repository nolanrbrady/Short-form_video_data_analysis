#!/usr/bin/env Rscript

# Shared publication figure styling for the SFV study plotting scripts.
#
# Sourced by:
#   - plot_behavior_score_distributions.R
#   - plot_significant_beta_value_distribution.R
#   - analyze_correlational_relationships_roi_means.R
#
# Purpose
#   - Keep behavioral, neural, and correlation figures on a single visual
#     system so they read as one coherent set in the manuscript.
#   - Centralize the theme, the shared descriptive-overlay caption, and the
#     figure display-name loader so a styling change is a one-file edit.
#
# This file defines functions/objects only; it performs no side effects on
# source() beyond those definitions.

suppressPackageStartupMessages({
  library(jsonlite)
  library(ggplot2)
})

# Shared publication theme applied to every figure so the panels read as one
# coherent set in the manuscript.
theme_sfv_pub <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.22), color = "grey10", margin = margin(b = 2)),
      plot.subtitle = element_text(size = rel(0.9), color = "grey35", margin = margin(b = 10)),
      plot.caption = element_text(size = rel(0.72), color = "grey40", hjust = 0, lineheight = 1.1, margin = margin(t = 12)),
      axis.title.x = element_text(face = "bold", size = rel(0.95), color = "grey20", margin = margin(t = 8)),
      axis.title.y = element_text(face = "bold", size = rel(0.95), color = "grey20", margin = margin(r = 8)),
      axis.text = element_text(color = "grey30"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
      plot.margin = margin(14, 18, 12, 14)
    )
}

# Load publication display labels (ROI names, condition/format/content labels,
# behavior-run names, format pools) from the shared figure-name config so no raw
# analysis codes reach the figures. `required_sections` are the top-level keys a
# given caller depends on; validation fails fast if any are absent.
load_figure_display_names <- function(path, required_sections) {
  if (!file.exists(path)) {
    stop(paste0("Figure display-name config not found: ", path))
  }
  obj <- tryCatch(
    jsonlite::fromJSON(path, simplifyVector = TRUE),
    error = function(e) {
      stop(paste0("Failed to parse figure display-name JSON at ", path, ". Error: ", conditionMessage(e)))
    }
  )
  missing <- setdiff(required_sections, names(obj))
  if (length(missing) > 0) {
    stop(paste0("Figure display-name config is missing sections: ", paste(missing, collapse = ", ")))
  }
  obj
}

# Map a p-value to conventional significance asterisks; non-significant returns
# "n.s." so a labeled bracket is still informative.
significance_stars <- function(p) {
  if (is.na(p)) {
    return("n.s.")
  }
  if (p < 0.001) {
    return("***")
  }
  if (p < 0.01) {
    return("**")
  }
  if (p < 0.05) {
    return("*")
  }
  "n.s."
}

# Build a horizontal significance bracket (connecting bar with two down-ticks)
# and a centered label above it, returned as a list of ggplot layers in data
# coordinates. `y` is the bar height, `tick` the drop of the bracket ends.
significance_bracket <- function(x1, x2, y, label, tick, color = "grey20") {
  star_only <- !identical(label, "n.s.")
  list(
    annotate("segment", x = x1, xend = x2, y = y, yend = y, color = color, linewidth = 0.5),
    annotate("segment", x = x1, xend = x1, y = y, yend = y - tick, color = color, linewidth = 0.5),
    annotate("segment", x = x2, xend = x2, y = y, yend = y - tick, color = color, linewidth = 0.5),
    annotate(
      "text",
      x = (x1 + x2) / 2,
      y = y,
      label = label,
      vjust = if (star_only) 0.2 else -0.3,
      size = if (star_only) 5.2 else 3.6,
      fontface = if (star_only) "bold" else "plain",
      color = color
    )
  )
}
