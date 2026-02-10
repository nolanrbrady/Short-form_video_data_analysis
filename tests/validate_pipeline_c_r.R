#!/usr/bin/env Rscript

# Validation harness for Pipeline C (R): Homer3 betas + Format×Content channelwise LMM.
#
# Purpose
# - Verify that `analyze_format_content_lmm_channelwise.R` adheres to `ANALYSIS_SPEC.md`:
#   - Subject ID normalization (`sub_0001` aligns with `0001`/`1`)
#   - Condition mapping + ±0.5 effect coding
#   - Pruned-channel policy: treat 0 as missing
#   - Complete-case within channel×chrom (all 4 conditions required)
#   - BH-FDR families per chrom × effect across channels
#   - Post-hoc gating: only interaction-FDR-significant channel×chrom pairs
#
# This script generates synthetic inputs under /tmp and runs the analysis script end-to-end.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

fail <- function(msg) {
  writeLines(paste0("[FAIL] ", msg))
  quit(status = 1)
}

assert_true <- function(cond, msg) {
  if (!isTRUE(cond)) fail(msg)
}

bh_qvalues_manual <- function(p_values) {
  # Manual Benjamini–Hochberg (BH) q-values for a numeric vector.
  # Reference: Benjamini & Hochberg (1995), see `CITATIONS.md`.
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
  # named list: channel -> chrom -> list(b0,bF,bC,bI)
  true_effects,
  noise_sd
) {
  ids <- make_subject_ids(n_subjects)
  grid <- cond_grid()

  # Subject random intercepts (shared across all channels/chroms for realism)
  set.seed(123)
  subj_re <- rnorm(n_subjects, mean = 0, sd = 0.6)
  names(subj_re) <- ids$homer_subject

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
        # generate per-subject observations
        vals <- mu + subj_re + rnorm(n_subjects, mean = 0, sd = noise_sd)
        col <- paste0(ch, "_Cond", cond, "_", chrom)
        out[[col]] <- vals
      }
    }
  }
  out
}

generate_combined <- function(n_subjects) {
  ids <- make_subject_ids(n_subjects)
  tibble::tibble(
    subject_id = ids$combined_subject,
    age = rep(20, n_subjects),
    pd_status = rep(0, n_subjects)
  )
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
  out_main,
  out_main_tidy,
  out_posthoc,
  alpha,
  min_subjects,
  exclude_json = NULL
) {
  script <- "analyze_format_content_lmm_channelwise.R"
  args <- c(
    "--input_csv", input_csv,
    "--out_main_csv", out_main,
    "--out_main_tidy_csv", out_main_tidy,
    "--out_posthoc_csv", out_posthoc,
    "--alpha", as.character(alpha),
    "--min_subjects", as.character(min_subjects)
  )
  if (!is.null(exclude_json)) {
    args <- c(args, "--exclude_subjects_json", exclude_json)
  }
  # Force full output tables for validation (do not allow filtering to drop non-significant rows).
  env <- c("FILTER_MAIN_EFFECTS_TO_FDR_SIG_ONLY=FALSE")
  status <- suppressWarnings(system2("Rscript", c(script, args), stdout = TRUE, stderr = TRUE, env = env))
  attr(status, "status") %||% 0
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

main <- function() {
  tmp <- file.path(tempdir(), "pipeline_c_validation")
  dir.create(tmp, showWarnings = FALSE, recursive = TRUE)

  input_csv <- file.path(tmp, "homer3_betas_plus_combined_sfv_data_inner_join.csv")
  out_main <- file.path(tmp, "main.csv")
  out_main_tidy <- file.path(tmp, "main_tidy.csv")
  out_posthoc <- file.path(tmp, "posthoc.csv")
  exclude_none <- file.path(tmp, "excluded_none.json")
  exclude_one <- file.path(tmp, "excluded_one.json")
  exclude_absent <- file.path(tmp, "excluded_absent.json")
  writeLines("[]", exclude_none)
  writeLines('["sub_0002"]', exclude_one)
  writeLines('["sub_9999"]', exclude_absent)

  n_subjects <- 12
  channels <- c("S01_D01", "S02_D01")
  chroms <- c("HbO", "HbR")

  # Channel 1: strong interaction; Channel 2: no interaction.
  true_effects <- list(
    S01_D01 = list(
      HbO = list(b0 = 0.2, bF = 0.7, bC = -0.3, bI = 2.5),
      HbR = list(b0 = -0.1, bF = -0.4, bC = 0.2, bI = -2.2)
    ),
    S02_D01 = list(
      HbO = list(b0 = 0.0, bF = 0.1, bC = 0.0, bI = 0.0),
      HbR = list(b0 = 0.0, bF = -0.1, bC = 0.0, bI = 0.0)
    )
  )

  homer <- generate_homer3_wide(
    n_subjects = n_subjects,
    channels = channels,
    chroms = chroms,
    true_effects = true_effects,
    noise_sd = 0.15
  )
  combined <- generate_combined(n_subjects)

  # Inject a pruned channel (0) for one subject in one condition: should force complete-case drop
  # for that channel×chrom only.
  homer[["S01_D01_Cond03_HbO"]][[1]] <- 0

  merged <- build_merged_input(homer, combined)
  write_csv(merged, input_csv)

  status <- run_analysis(
    input_csv = input_csv,
    out_main = out_main,
    out_main_tidy = out_main_tidy,
    out_posthoc = out_posthoc,
    alpha = 0.05,
    min_subjects = 6,
    exclude_json = exclude_none
  )
  assert_true(status == 0, "analysis script did not exit cleanly")
  assert_true(file.exists(out_main), "missing main output CSV")
  assert_true(file.exists(out_main_tidy), "missing tidy main output CSV")
  assert_true(file.exists(out_posthoc), "missing posthoc output CSV")

  main_tidy <- read_csv(out_main_tidy, show_col_types = FALSE)
  posthoc <- read_csv(out_posthoc, show_col_types = FALSE)

  # 1) Tidy table has required columns
  required_main_cols <- c(
    "channel",
    "chrom",
    "effect",
    "n_subjects",
    "n_obs",
    "singular_fit",
    "estimate",
    "ci95_low",
    "ci95_high",
    "p_unc",
    "p_fdr"
  )
  assert_true(all(required_main_cols %in% names(main_tidy)), "tidy main output missing required columns")

  # 1b) Complete-case invariant: in this pipeline, each subject contributes 4 observations per channel×chrom.
  assert_true(all(main_tidy$n_obs == 4 * main_tidy$n_subjects), "expected n_obs == 4 * n_subjects for all tidy rows")

  # 2) Pruned channel policy + complete-case: HbO S01_D01 should drop exactly 1 subject (n_subjects = 11)
  row <- main_tidy %>% filter(channel == "S01_D01", chrom == "HbO", effect == "interaction")
  assert_true(nrow(row) == 1, "expected exactly one tidy row for S01_D01 HbO interaction")
  assert_true(row$n_subjects[[1]] == (n_subjects - 1), "expected HbO S01_D01 to drop one subject due to pruned beta==0")

  # 3) Coefficient recovery for strong-effect channel (within tolerance)
  tol <- 0.25
  for (chrom_name in chroms) {
    pars <- true_effects[["S01_D01"]][[chrom_name]]
    est_fmt <- main_tidy %>% filter(channel == "S01_D01", chrom == chrom_name, effect == "format") %>% pull(estimate)
    est_cnt <- main_tidy %>% filter(channel == "S01_D01", chrom == chrom_name, effect == "content") %>% pull(estimate)
    est_int <- main_tidy %>% filter(channel == "S01_D01", chrom == chrom_name, effect == "interaction") %>% pull(estimate)
    assert_true(length(est_fmt) == 1, paste0("missing format estimate for S01_D01 ", chrom_name))
    assert_true(length(est_cnt) == 1, paste0("missing content estimate for S01_D01 ", chrom_name))
    assert_true(length(est_int) == 1, paste0("missing interaction estimate for S01_D01 ", chrom_name))

    assert_true(abs(est_fmt[[1]] - pars$bF) < tol, paste0("format estimate off for S01_D01 ", chrom_name))
    assert_true(abs(est_cnt[[1]] - pars$bC) < tol, paste0("content estimate off for S01_D01 ", chrom_name))
    assert_true(abs(est_int[[1]] - pars$bI) < tol, paste0("interaction estimate off for S01_D01 ", chrom_name))
  }

  # 3b) Channel subsetting check: S01_D01 and S02_D01 should not produce identical estimates.
  # This catches accidental fitting on the full dataset repeatedly.
  fmt_s01 <- main_tidy %>% filter(channel == "S01_D01", chrom == "HbO", effect == "format") %>% pull(estimate)
  fmt_s02 <- main_tidy %>% filter(channel == "S02_D01", chrom == "HbO", effect == "format") %>% pull(estimate)
  assert_true(abs(fmt_s01[[1]] - fmt_s02[[1]]) > 0.2, "channelwise subsetting appears broken (format estimates too similar)")

  # 4) Post-hoc gating: should only include the channel with strong interaction (S01_D01), not S02_D01
  assert_true(all(posthoc$channel %in% c("S01_D01")), "posthoc contains channels that should not have passed interaction gate")

  # 5) Post-hoc contrast count: 6 contrasts per chrom for one channel
  expected_posthoc_rows <- 6 * length(chroms) * 1
  assert_true(nrow(posthoc) == expected_posthoc_rows, "unexpected number of posthoc rows")

  # 5b) BH-FDR implementation check: p_fdr must match a manual BH implementation within each (chrom × effect) family.
  # Reference: Benjamini & Hochberg (1995), see `CITATIONS.md`.
  for (chrom_name in chroms) {
    for (effect_name in c("format", "content", "interaction")) {
      sub <- main_tidy %>% filter(chrom == chrom_name, effect == effect_name) %>% arrange(channel)
      q_manual <- bh_qvalues_manual(sub$p_unc)
      max_abs <- max(abs(q_manual - sub$p_fdr), na.rm = TRUE)
      assert_true(max_abs < 1e-10, paste0("BH q-values mismatch for ", chrom_name, " ", effect_name))
    }
  }

  # 6) Validate mean_diff sign convention for one contrast:
  # mean_diff = condition_a - condition_b, so for SF_Edu vs SF_Ent it should match
  # the average within-subject (SF_Edu - SF_Ent) difference.
  # Recompute from the synthetic inputs for HbO and compare within tolerance.
  # Build long form from homer only for S01_D01 HbO.
  grid <- cond_grid()
  ids <- make_subject_ids(n_subjects)$homer_subject
  long <- tibble::tibble(subject = ids) %>%
    transmute(
      subject_id = as.integer(str_extract(subject, "\\d+")),
      SF_Edu = homer[["S01_D01_Cond01_HbO"]],
      SF_Ent = homer[["S01_D01_Cond02_HbO"]],
      LF_Ent = homer[["S01_D01_Cond03_HbO"]],
      LF_Edu = homer[["S01_D01_Cond04_HbO"]]
    ) %>%
    pivot_longer(cols = c("SF_Edu", "SF_Ent", "LF_Ent", "LF_Edu"), names_to = "condition", values_to = "beta") %>%
    mutate(beta = na_if(beta, 0)) %>%
    pivot_wider(names_from = condition, values_from = beta)
  # Complete cases only (matches analysis)
  long_cc <- long %>% drop_na(SF_Edu, SF_Ent, LF_Ent, LF_Edu)
  avg_diff <- mean(long_cc$SF_Edu - long_cc$SF_Ent)
  ph_row <- posthoc %>%
    filter(chrom == "HbO", condition_a == "SF_Edu", condition_b == "SF_Ent")
  assert_true(nrow(ph_row) == 1, "missing SF_Edu vs SF_Ent posthoc row for HbO")
  assert_true(abs(ph_row$mean_diff[[1]] - avg_diff) < 0.35, "posthoc mean_diff sign/magnitude appears inconsistent")

  # 6b) Central exclusion list: one present subject should reduce modeled subject counts by one.
  status_excl <- run_analysis(
    input_csv = input_csv,
    out_main = file.path(tmp, "main_excluded.csv"),
    out_main_tidy = file.path(tmp, "main_excluded_tidy.csv"),
    out_posthoc = file.path(tmp, "posthoc_excluded.csv"),
    alpha = 0.05,
    min_subjects = 6,
    exclude_json = exclude_one
  )
  assert_true(status_excl == 0, "analysis failed with one-subject exclusion list")
  main_excl <- read_csv(file.path(tmp, "main_excluded_tidy.csv"), show_col_types = FALSE)
  row_excl <- main_excl %>% filter(channel == "S01_D01", chrom == "HbO", effect == "interaction")
  assert_true(row_excl$n_subjects[[1]] == (n_subjects - 2), "excluded subject plus one pruned-case drop should yield n_subjects=10")
  row_excl_unaffected <- main_excl %>% filter(channel == "S02_D01", chrom == "HbO", effect == "interaction")
  assert_true(row_excl_unaffected$n_subjects[[1]] == (n_subjects - 1), "exclusion should reduce unaffected channel n_subjects by one")

  # 6c) Missing exclusion IDs should warn/continue and preserve results.
  status_abs <- run_analysis(
    input_csv = input_csv,
    out_main = file.path(tmp, "main_absent.csv"),
    out_main_tidy = file.path(tmp, "main_absent_tidy.csv"),
    out_posthoc = file.path(tmp, "posthoc_absent.csv"),
    alpha = 0.05,
    min_subjects = 6,
    exclude_json = exclude_absent
  )
  assert_true(status_abs == 0, "analysis failed when exclusion ID was absent")
  main_abs <- read_csv(file.path(tmp, "main_absent_tidy.csv"), show_col_types = FALSE)
  assert_true(
    all(main_abs$n_subjects == main_tidy$n_subjects),
    "absent exclusion ID should not change per-row n_subjects"
  )

  # 7) Duplicate subject_id should fail hard (data integrity)
  merged_dup <- bind_rows(merged, merged %>% slice(1))
  merged_dup_csv <- file.path(tmp, "merged_DUP.csv")
  write_csv(merged_dup, merged_dup_csv)
  status2 <- run_analysis(
    input_csv = merged_dup_csv,
    out_main = file.path(tmp, "main2.csv"),
    out_main_tidy = file.path(tmp, "main2_tidy.csv"),
    out_posthoc = file.path(tmp, "posthoc2.csv"),
    alpha = 0.05,
    min_subjects = 6,
    exclude_json = exclude_none
  )
  assert_true(status2 != 0, "expected analysis to fail hard on duplicate subject_id in merged input")

  # 8) Unexpected condition codes should fail hard (avoid silent dropping)
  merged_bad <- merged
  # Add a single extra beta column with an invalid condition code (Cond05) that matches the beta naming pattern.
  merged_bad[["S01_D01_Cond05_HbO"]] <- 0.1
  merged_bad_csv <- file.path(tmp, "merged_BAD_COND.csv")
  write_csv(merged_bad, merged_bad_csv)
  status3 <- run_analysis(
    input_csv = merged_bad_csv,
    out_main = file.path(tmp, "main3.csv"),
    out_main_tidy = file.path(tmp, "main3_tidy.csv"),
    out_posthoc = file.path(tmp, "posthoc3.csv"),
    alpha = 0.05,
    min_subjects = 6,
    exclude_json = exclude_none
  )
  assert_true(status3 != 0, "expected analysis to fail hard on unexpected condition codes")

  writeLines(paste0("[OK] Pipeline C (R) validation passed. Outputs in: ", tmp))
}

main()
