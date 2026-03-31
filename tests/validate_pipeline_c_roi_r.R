#!/usr/bin/env Rscript

# Validation harness for Pipeline C ROI (R): Homer3 betas + Format×Content ROI-wise LMM.
#
# Purpose
# - Verify that `analyze_format_content_lmm_roi.R` implements:
#   - strict ROI JSON parsing and channel-map validation
#   - subject ID normalization (`sub_0001` aligns with `0001`/`1`)
#   - condition mapping + ±0.5 effect coding
#   - age-adjusted omnibus model (`+ age`) with fail-hard covariate validation
#   - pruned-channel policy: treat NA as missing
#   - ROI beta aggregation as mean across available channels
#   - complete-case within ROI×chrom (all 4 conditions required)
#   - BH-FDR families per chrom × effect across ROIs
#   - post-hoc gating only for interaction-FDR-significant ROI×chrom pairs

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(jsonlite)
})

fail <- function(msg) {
  writeLines(paste0("[FAIL] ", msg))
  quit(status = 1)
}

assert_true <- function(cond, msg) {
  if (!isTRUE(cond)) fail(msg)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

bh_qvalues_manual <- function(p_values) {
  p <- as.numeric(p_values)
  q <- rep(NA_real_, length(p))
  finite <- is.finite(p)
  if (!any(finite)) return(q)

  p_f <- p[finite]
  o <- order(p_f)
  ranked <- p_f[o]
  m <- length(ranked)
  q_ranked <- rep(NA_real_, m)
  prev <- 1.0
  for (i in seq(m, 1)) {
    rank <- i
    val <- (m / rank) * ranked[[i]]
    prev <- min(prev, val)
    q_ranked[[i]] <- prev
  }
  q_ranked <- pmin(pmax(q_ranked, 0.0), 1.0)
  q_f <- rep(NA_real_, m)
  q_f[o] <- q_ranked
  q[finite] <- q_f
  q
}

make_subject_ids <- function(n_subjects) {
  ids_int <- seq_len(n_subjects)
  list(
    homer_subject = sprintf("sub_%04d", ids_int),
    combined_subject = sprintf("%04d", ids_int),
    ids_int = ids_int
  )
}

make_age_years <- function(n_subjects) {
  seq(18, by = 1, length.out = n_subjects)
}

cond_grid <- function() {
  tibble::tibble(
    cond = c("01", "02", "03", "04"),
    condition = c("SF_Edu", "SF_Ent", "LF_Ent", "LF_Edu"),
    format_c = c(-0.5, -0.5, +0.5, +0.5),
    content_c = c(+0.5, -0.5, -0.5, +0.5)
  )
}

expected_mean <- function(b0, bF, bC, bI, format_c, content_c) {
  b0 + bF * format_c + bC * content_c + bI * (format_c * content_c)
}

generate_homer3_wide <- function(
  n_subjects,
  channels,
  chroms,
  true_effects,
  noise_sd,
  age_years,
  age_beta = 0.12
) {
  ids <- make_subject_ids(n_subjects)
  grid <- cond_grid()

  set.seed(123)
  subj_re <- rnorm(n_subjects, mean = 0, sd = 0.6)
  names(subj_re) <- ids$homer_subject
  age_centered <- age_years - mean(age_years)

  out <- tibble::tibble(Subject = ids$homer_subject)

  for (ch in channels) {
    for (chrom in chroms) {
      pars <- true_effects[[ch]][[chrom]]
      for (i in seq_len(nrow(grid))) {
        cond <- grid$cond[[i]]
        mu <- expected_mean(
          pars$b0,
          pars$bF,
          pars$bC,
          pars$bI,
          grid$format_c[[i]],
          grid$content_c[[i]]
        )
        vals <- mu + age_beta * age_centered + subj_re + rnorm(n_subjects, mean = 0, sd = noise_sd)
        col <- paste0(ch, "_Cond", cond, "_", chrom)
        out[[col]] <- vals
      }
    }
  }
  out
}

generate_combined <- function(n_subjects, age_years = make_age_years(n_subjects)) {
  ids <- make_subject_ids(n_subjects)
  tibble::tibble(
    subject_id = ids$combined_subject,
    age = age_years,
    pd_status = rep(0, n_subjects)
  )
}

build_reference_roi_long <- function(merged, roi_channels, chrom_name, roi_name) {
  conds <- c("01", "02", "03", "04")
  cond_labels <- c("SF_Edu", "SF_Ent", "LF_Ent", "LF_Edu")

  rows <- lapply(seq_along(conds), function(i) {
    cond <- conds[[i]]
    cond_label <- cond_labels[[i]]
    beta_cols <- paste0(roi_channels, "_Cond", cond, "_", chrom_name)
    beta_mat <- as.matrix(merged[, beta_cols, drop = FALSE])
    beta_mean <- apply(beta_mat, 1, function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE))
    tibble::tibble(
      subject_id = as.integer(str_extract(merged$subject_id, "\\d+")),
      age = merged$age,
      roi = roi_name,
      chrom = chrom_name,
      condition = cond_label,
      beta = beta_mean
    )
  })

  bind_rows(rows) %>%
    mutate(
      format_c = case_when(
        condition %in% c("SF_Edu", "SF_Ent") ~ -0.5,
        condition %in% c("LF_Ent", "LF_Edu") ~ +0.5,
        TRUE ~ NA_real_
      ),
      content_c = case_when(
        condition %in% c("SF_Ent", "LF_Ent") ~ -0.5,
        condition %in% c("SF_Edu", "LF_Edu") ~ +0.5,
        TRUE ~ NA_real_
      )
    ) %>%
    filter(!is.na(.data$beta))
}

build_merged_input <- function(homer, combined) {
  homer_norm <- homer %>%
    mutate(subject_id_norm = as.integer(str_extract(.data$Subject, "\\d+"))) %>%
    select(-Subject)
  combined_norm <- combined %>%
    mutate(subject_id_norm = as.integer(str_extract(.data$subject_id, "\\d+")))

  merged <- combined_norm %>%
    inner_join(homer_norm, by = "subject_id_norm", relationship = "one-to-one") %>%
    select(-subject_id_norm)

  dup <- merged %>%
    count(subject_id, name = "n_rows") %>%
    filter(n_rows > 1)
  if (nrow(dup) > 0) {
    stop("build_merged_input produced duplicate subject_id rows unexpectedly.")
  }
  merged
}

run_analysis <- function(
  input_csv,
  roi_json,
  out_main,
  out_main_tidy,
  out_posthoc,
  alpha,
  min_subjects,
  exclude_json = NULL
) {
  script <- "analyze_format_content_lmm_roi.R"
  args <- c(
    "--input_csv", input_csv,
    "--roi_json", roi_json,
    "--out_main_csv", out_main,
    "--out_main_tidy_csv", out_main_tidy,
    "--out_posthoc_csv", out_posthoc,
    "--alpha", as.character(alpha),
    "--min_subjects", as.character(min_subjects)
  )
  if (!is.null(exclude_json)) {
    args <- c(args, "--exclude_subjects_json", exclude_json)
  }
  env <- c("FILTER_MAIN_EFFECTS_TO_FDR_SIG_ONLY=FALSE")
  status <- suppressWarnings(system2("Rscript", c(script, args), stdout = TRUE, stderr = TRUE, env = env))
  attr(status, "status") %||% 0
}

main <- function() {
  tmp <- file.path(tempdir(), "pipeline_c_roi_validation")
  dir.create(tmp, showWarnings = FALSE, recursive = TRUE)

  input_csv <- file.path(tmp, "homer3_betas_plus_combined_sfv_data_inner_join.csv")
  roi_json <- file.path(tmp, "roi_definition.json")
  roi_json_invalid <- file.path(tmp, "roi_definition_invalid.json")
  roi_json_missing_channel <- file.path(tmp, "roi_definition_missing_channel.json")
  roi_json_overlap <- file.path(tmp, "roi_definition_overlap.json")

  out_main <- file.path(tmp, "main.csv")
  out_main_tidy <- file.path(tmp, "main_tidy.csv")
  out_posthoc <- file.path(tmp, "posthoc.csv")

  exclude_none <- file.path(tmp, "excluded_none.json")
  writeLines("[]", exclude_none)

  n_subjects <- 12
  channels <- c("S01_D01", "S01_D02", "S02_D01", "S02_D02")
  chroms <- c("HbO", "HbR")
  age_years <- make_age_years(n_subjects)

  # ROI 1 channels have strong interaction; ROI 2 channels have null interaction.
  true_effects <- list(
    S01_D01 = list(
      HbO = list(b0 = 0.2, bF = 0.8, bC = -0.3, bI = 2.4),
      HbR = list(b0 = -0.1, bF = -0.5, bC = 0.2, bI = -2.2)
    ),
    S01_D02 = list(
      HbO = list(b0 = 0.3, bF = 0.6, bC = -0.2, bI = 2.0),
      HbR = list(b0 = -0.2, bF = -0.3, bC = 0.1, bI = -1.8)
    ),
    S02_D01 = list(
      HbO = list(b0 = 0.0, bF = 0.1, bC = 0.0, bI = 0.0),
      HbR = list(b0 = 0.0, bF = -0.1, bC = 0.0, bI = 0.0)
    ),
    S02_D02 = list(
      HbO = list(b0 = 0.0, bF = 0.0, bC = 0.1, bI = 0.0),
      HbR = list(b0 = 0.0, bF = 0.0, bC = -0.1, bI = 0.0)
    )
  )

  homer <- generate_homer3_wide(
    n_subjects = n_subjects,
    channels = channels,
    chroms = chroms,
    true_effects = true_effects,
    noise_sd = 0.15,
    age_years = age_years
  )

  # Pruning injections for ROI mean/complete-case checks:
  # - subject 1 has both ROI1 channels pruned for Cond03 HbO -> ROI1/HbO complete-case drop.
  # - subject 2 has one ROI1 channel pruned for Cond02 HbO only -> should remain included.
  homer[["S01_D01_Cond03_HbO"]][[1]] <- NA_real_
  homer[["S01_D02_Cond03_HbO"]][[1]] <- NA_real_
  homer[["S01_D01_Cond02_HbO"]][[2]] <- NA_real_
  # Inject one literal zero beta in a non-pruned ROI channel; this should remain valid data.
  homer[["S02_D01_Cond01_HbO"]][[3]] <- 0

  combined <- generate_combined(n_subjects, age_years = age_years)
  merged <- build_merged_input(homer, combined)
  write_csv(merged, input_csv)

  roi_map <- list(
    VMPFC = c("S01_D01", "S01_D02"),
    DLPFC = c("S02_D01", "S02_D02")
  )
  write_json(roi_map, roi_json, auto_unbox = TRUE, pretty = TRUE)

  status <- run_analysis(
    input_csv = input_csv,
    roi_json = roi_json,
    out_main = out_main,
    out_main_tidy = out_main_tidy,
    out_posthoc = out_posthoc,
    alpha = 0.05,
    min_subjects = 6,
    exclude_json = exclude_none
  )

  assert_true(status == 0, "ROI analysis script did not exit cleanly")
  assert_true(file.exists(out_main), "missing main output CSV")
  assert_true(file.exists(out_main_tidy), "missing tidy main output CSV")
  assert_true(file.exists(out_posthoc), "missing posthoc output CSV")

  main_tidy <- read_csv(out_main_tidy, show_col_types = FALSE)
  posthoc <- read_csv(out_posthoc, show_col_types = FALSE)

  required_main_cols <- c(
    "roi",
    "chrom",
    "effect",
    "n_subjects",
    "n_obs",
    "converged",
    "singular_fit",
    "estimate",
    "ci95_low",
    "ci95_high",
    "p_unc",
    "p_fdr"
  )
  assert_true(all(required_main_cols %in% names(main_tidy)), "tidy main output missing required columns")
  assert_true(all(main_tidy$converged), "expected all synthetic ROI fits to converge cleanly")
  assert_true(!is.unsorted(main_tidy$p_unc), "expected tidy main output to be sorted by ascending p_unc")

  # Exactly 2 ROIs x 2 chrom x 3 effects
  assert_true(nrow(main_tidy) == 12, "unexpected number of tidy rows for 2 ROIs x 2 chroms x 3 effects")

  # ROI complete-case check (subject 1 dropped only for VMPFC HbO)
  row_vmpfc_hbo_int <- main_tidy %>% filter(roi == "VMPFC", chrom == "HbO", effect == "interaction")
  assert_true(nrow(row_vmpfc_hbo_int) == 1, "expected one row for VMPFC HbO interaction")
  assert_true(
    row_vmpfc_hbo_int$n_subjects[[1]] == (n_subjects - 1),
    "expected VMPFC HbO to drop exactly one subject due to all-ROI-channel pruning"
  )

  # Partial pruning should not drop an additional subject when at least one ROI channel remains.
  assert_true(
    row_vmpfc_hbo_int$n_subjects[[1]] == (n_subjects - 1),
    "partial ROI-channel pruning appears to be treated incorrectly (all-channels-required behavior detected)"
  )
  row_dlpfc_hbo_int <- main_tidy %>% filter(roi == "DLPFC", chrom == "HbO", effect == "interaction")
  assert_true(nrow(row_dlpfc_hbo_int) == 1, "expected one row for DLPFC HbO interaction")
  assert_true(
    row_dlpfc_hbo_int$n_subjects[[1]] == n_subjects,
    "literal zero ROI beta appears to have been treated as missing"
  )

  # Coefficient recovery for ROI means (tolerance relaxed for random effects + noise).
  tol <- 0.30
  expected_roi <- list(
    VMPFC = list(
      HbO = list(bF = mean(c(0.8, 0.6)), bC = mean(c(-0.3, -0.2)), bI = mean(c(2.4, 2.0))),
      HbR = list(bF = mean(c(-0.5, -0.3)), bC = mean(c(0.2, 0.1)), bI = mean(c(-2.2, -1.8)))
    ),
    DLPFC = list(
      HbO = list(bF = mean(c(0.1, 0.0)), bC = mean(c(0.0, 0.1)), bI = mean(c(0.0, 0.0))),
      HbR = list(bF = mean(c(-0.1, 0.0)), bC = mean(c(0.0, -0.1)), bI = mean(c(0.0, 0.0)))
    )
  )

  for (roi_name in names(expected_roi)) {
    for (chrom_name in names(expected_roi[[roi_name]])) {
      pars <- expected_roi[[roi_name]][[chrom_name]]
      est_fmt <- main_tidy %>% filter(roi == roi_name, chrom == chrom_name, effect == "format") %>% pull(estimate)
      est_cnt <- main_tidy %>% filter(roi == roi_name, chrom == chrom_name, effect == "content") %>% pull(estimate)
      est_int <- main_tidy %>% filter(roi == roi_name, chrom == chrom_name, effect == "interaction") %>% pull(estimate)
      assert_true(length(est_fmt) == 1, paste0("missing format estimate for ", roi_name, " ", chrom_name))
      assert_true(length(est_cnt) == 1, paste0("missing content estimate for ", roi_name, " ", chrom_name))
      assert_true(length(est_int) == 1, paste0("missing interaction estimate for ", roi_name, " ", chrom_name))
      assert_true(abs(est_fmt[[1]] - pars$bF) < tol, paste0("format estimate off for ", roi_name, " ", chrom_name))
      assert_true(abs(est_cnt[[1]] - pars$bC) < tol, paste0("content estimate off for ", roi_name, " ", chrom_name))
      assert_true(abs(est_int[[1]] - pars$bI) < tol, paste0("interaction estimate off for ", roi_name, " ", chrom_name))
    }
  }

  # Direct reference-model agreement for one representative ROI/chrom with age adjustment.
  ref_long <- build_reference_roi_long(
    merged,
    roi_channels = c("S01_D01", "S01_D02"),
    chrom_name = "HbO",
    roi_name = "VMPFC"
  ) %>%
    group_by(subject_id) %>%
    filter(n_distinct(condition) == 4) %>%
    ungroup()
  ref_model <- lmerTest::lmer(beta ~ format_c * content_c + age + (1 | subject_id), data = ref_long, REML = TRUE)
  ref_anova <- anova(ref_model, ddf = "Kenward-Roger", type = 3)
  ref_coefs <- summary(ref_model)$coefficients
  for (term_map in list(c("format", "format_c"), c("content", "content_c"), c("interaction", "format_c:content_c"))) {
    out_row <- main_tidy %>% filter(roi == "VMPFC", chrom == "HbO", effect == term_map[[1]])
    assert_true(nrow(out_row) == 1, paste0("missing output row for reference term ", term_map[[1]]))
    assert_true(abs(out_row$estimate[[1]] - ref_coefs[term_map[[2]], "Estimate"]) < 1e-10, paste0("estimate mismatch for ", term_map[[1]], " reference fit"))
    assert_true(abs(out_row$p_unc[[1]] - ref_anova[term_map[[2]], "Pr(>F)"]) < 1e-7, paste0("p_unc mismatch for ", term_map[[1]], " reference fit"))
  }

  # Explicit inferential TP/TN checks:
  # - True positive: VMPFC interaction should be FDR-significant.
  # - True negative: DLPFC interaction should remain non-significant.
  for (chrom_name in unique(main_tidy$chrom)) {
    p_sig <- main_tidy %>%
      filter(roi == "VMPFC", chrom == chrom_name, effect == "interaction") %>%
      pull(p_fdr)
    p_null <- main_tidy %>%
      filter(roi == "DLPFC", chrom == chrom_name, effect == "interaction") %>%
      pull(p_fdr)
    assert_true(length(p_sig) == 1, paste0("missing TP p_fdr row for VMPFC ", chrom_name))
    assert_true(length(p_null) == 1, paste0("missing TN p_fdr row for DLPFC ", chrom_name))
    assert_true(p_sig[[1]] < 0.05, paste0("expected TP interaction significance for VMPFC ", chrom_name))
    assert_true(p_null[[1]] >= 0.05, paste0("expected TN non-significance for DLPFC ", chrom_name))
  }

  # Post-hoc gating: only the interaction-positive ROI should appear.
  assert_true(nrow(posthoc) > 0, "expected non-empty posthoc output")
  assert_true(all(posthoc$roi %in% c("VMPFC")), "posthoc contains ROIs that should not pass interaction gate")
  assert_true("converged" %in% names(posthoc), "ROI posthoc output missing converged column")
  assert_true(all(posthoc$converged), "expected all synthetic ROI posthoc fits to converge cleanly")

  # BH-FDR checks per chrom × effect across ROIs.
  for (chrom_name in unique(main_tidy$chrom)) {
    for (eff in c("format", "content", "interaction")) {
      sub <- main_tidy %>% filter(chrom == chrom_name, effect == eff) %>% arrange(roi)
      q_manual <- bh_qvalues_manual(sub$p_unc)
      max_abs_diff <- max(abs(sub$p_fdr - q_manual), na.rm = TRUE)
      if (!is.finite(max_abs_diff)) max_abs_diff <- 0
      assert_true(max_abs_diff < 1e-10, paste0("BH-FDR mismatch for ", chrom_name, " / ", eff))
    }
  }

  # Age is required and must be complete/numeric after exclusions.
  merged_missing_age_col <- merged %>% select(-age)
  merged_missing_age_col_csv <- file.path(tmp, "merged_missing_age_col.csv")
  write_csv(merged_missing_age_col, merged_missing_age_col_csv)
  status_missing_age_col <- run_analysis(
    input_csv = merged_missing_age_col_csv,
    roi_json = roi_json,
    out_main = file.path(tmp, "main_missing_age_col.csv"),
    out_main_tidy = file.path(tmp, "main_tidy_missing_age_col.csv"),
    out_posthoc = file.path(tmp, "posthoc_missing_age_col.csv"),
    alpha = 0.05,
    min_subjects = 6,
    exclude_json = exclude_none
  )
  assert_true(status_missing_age_col != 0, "missing age column should fail fast")

  merged_age_na <- merged
  merged_age_na$age[[3]] <- NA_real_
  merged_age_na_csv <- file.path(tmp, "merged_age_na.csv")
  write_csv(merged_age_na, merged_age_na_csv)
  status_age_na <- run_analysis(
    input_csv = merged_age_na_csv,
    roi_json = roi_json,
    out_main = file.path(tmp, "main_age_na.csv"),
    out_main_tidy = file.path(tmp, "main_tidy_age_na.csv"),
    out_posthoc = file.path(tmp, "posthoc_age_na.csv"),
    alpha = 0.05,
    min_subjects = 6,
    exclude_json = exclude_none
  )
  assert_true(status_age_na != 0, "missing age values should fail fast")

  merged_age_bad <- merged
  merged_age_bad$age <- as.character(merged_age_bad$age)
  merged_age_bad$age[[4]] <- "not_numeric"
  merged_age_bad_csv <- file.path(tmp, "merged_age_bad.csv")
  write_csv(merged_age_bad, merged_age_bad_csv)
  status_age_bad <- run_analysis(
    input_csv = merged_age_bad_csv,
    roi_json = roi_json,
    out_main = file.path(tmp, "main_age_bad.csv"),
    out_main_tidy = file.path(tmp, "main_tidy_age_bad.csv"),
    out_posthoc = file.path(tmp, "posthoc_age_bad.csv"),
    alpha = 0.05,
    min_subjects = 6,
    exclude_json = exclude_none
  )
  assert_true(status_age_bad != 0, "non-numeric age values should fail fast")

  # Strict JSON parsing should fail on invalid JSON syntax.
  writeLines("{'VMPFC': ['S01_D01', 'S01_D02']}", roi_json_invalid)
  status_bad_json <- run_analysis(
    input_csv = input_csv,
    roi_json = roi_json_invalid,
    out_main = file.path(tmp, "main_bad_json.csv"),
    out_main_tidy = file.path(tmp, "main_tidy_bad_json.csv"),
    out_posthoc = file.path(tmp, "posthoc_bad_json.csv"),
    alpha = 0.05,
    min_subjects = 6,
    exclude_json = exclude_none
  )
  assert_true(status_bad_json != 0, "invalid ROI JSON should fail fast")

  # ROI channels missing from data should fail hard.
  roi_missing <- list(VMPFC = c("S01_D01"), DLPFC = c("S99_D99"))
  write_json(roi_missing, roi_json_missing_channel, auto_unbox = TRUE, pretty = TRUE)
  status_missing_channel <- run_analysis(
    input_csv = input_csv,
    roi_json = roi_json_missing_channel,
    out_main = file.path(tmp, "main_missing_channel.csv"),
    out_main_tidy = file.path(tmp, "main_tidy_missing_channel.csv"),
    out_posthoc = file.path(tmp, "posthoc_missing_channel.csv"),
    alpha = 0.05,
    min_subjects = 6,
    exclude_json = exclude_none
  )
  assert_true(status_missing_channel != 0, "ROI map with missing channels should fail fast")

  # Overlapping ROI assignments should fail hard.
  roi_overlap <- list(VMPFC = c("S01_D01", "S01_D02"), DLPFC = c("S01_D02", "S02_D01"))
  write_json(roi_overlap, roi_json_overlap, auto_unbox = TRUE, pretty = TRUE)
  status_overlap <- run_analysis(
    input_csv = input_csv,
    roi_json = roi_json_overlap,
    out_main = file.path(tmp, "main_overlap.csv"),
    out_main_tidy = file.path(tmp, "main_tidy_overlap.csv"),
    out_posthoc = file.path(tmp, "posthoc_overlap.csv"),
    alpha = 0.05,
    min_subjects = 6,
    exclude_json = exclude_none
  )
  assert_true(status_overlap != 0, "overlapping ROI channel assignments should fail fast")

  writeLines("[PASS] ROI Pipeline C validation checks passed")
}

main()
