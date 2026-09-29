#!/usr/bin/env Rscript

# Scientific/data-integrity validation for exhaustive equal-item sensitivity.
# Run from repo root. Fixtures and outputs are confined to a temporary directory.
# Scoring expectations are computed independently below; primary inference has
# additional positive/null and posthoc gate tests in validate_retention_pipeline_r.R.
# Holm (1979), https://doi.org/10.2307/4615733; Bates et al. (2015),
# https://doi.org/10.18637/jss.v067.i01. See CITATIONS.md.
source("analyze_retention_sensitivity.R")
scratch <- tempfile("validate-retention-sensitivity-")
dir.create(scratch)

# Each assertion names the scientific invariant so failures can be diagnosed
# without treating a successful model fit as evidence of correct input handling.
check <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}
near <- function(x, y, tol = 1e-12) {
  isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = tol, scale = 1))
}
expect_error <- function(expr, pattern) {
  message <- tryCatch({ force(expr); NULL }, error = function(e) conditionMessage(e))
  check(!is.null(message) && grepl(pattern, message), paste("Expected error:", pattern, "got:", message))
}

args <- parse_sensitivity_args(character())
args$out_dir <- file.path(scratch, "results")
prepared <- prepare_sensitivity(args)
invalid <- jsonlite::fromJSON(args$invalid_questions_json)$invalid_questions
pre <- read_sensitivity_audit(args$audit_pre_csv, invalid)
post <- read_sensitivity_audit(args$audit_post_csv, invalid)
raw_pre <- readr::read_csv(args$audit_pre_csv, col_types = readr::cols(.default = "c"), show_col_types = FALSE)
raw_post <- readr::read_csv(args$audit_post_csv, col_types = readr::cols(.default = "c"), show_col_types = FALSE)

# Check all subsets against participant/item loops, independent of the scoring
# implementation's grouped means. Unchanged outcomes/covariates must be identical.
test_complete_enumeration <- function() {
  check(identical(prepared$questions, c("Q3", "Q4", "Q11", "Q12", "Q19", "Q20", "Q27")),
        "Unexpected valid target items in the study fixture.")
  check(length(prepared$versions) == 8L, "Must include baseline and all seven omissions.")
  baseline <- prepared$versions$baseline
  unchanged <- setdiff(names(baseline), "diff_short_form_entertainment")
  for (question in prepared$questions) {
    v <- prepared$versions[[paste0("omit_", question)]]
    check(identical(v[unchanged], baseline[unchanged]), "Omission changed covariates, cohort, or other conditions.")
    expected <- vapply(v$subject_id, function(id) {
      retained <- setdiff(prepared$questions, question)
      delta <- vapply(retained, function(q) {
        post$score[post$subject_id == id & post$question_id == q] -
          pre$score[pre$subject_id == id & pre$question_id == q]
      }, numeric(1))
      sum(delta) / 6
    }, numeric(1))
    check(near(v$diff_short_form_entertainment, expected), "Six-item paired scoring disagrees with manual calculation.")
  }
  check(near(rowMeans(sapply(prepared$versions[-1], function(x) x$diff_short_form_entertainment)),
             baseline$diff_short_form_entertainment), "Equal subset-inclusion identity failed.")
}

# Known item changes establish denominator/sign behavior, including a valid zero:
# one post-only correct item gives +1/7 at baseline, 0 when dropped, +1/6 otherwise;
# reversing the phase yields the negative change. No missing-value imputation.
test_known_item_changes <- function() {
  a <- pre
  b <- post
  a$score[!a$excluded] <- 0
  b$score[!b$excluded] <- 0
  id <- a$subject_id[[1]]
  b$score[b$subject_id == id & b$question_id == "Q3"] <- 1
  check(near(score_sensitivity_audits(a, b, id)[, 2], 1 / 7), "Baseline denominator/sign incorrect.")
  check(near(score_sensitivity_audits(a, b, id, "Q3")[, 2], 0), "Dropped improvement must become a valid zero.")
  check(near(score_sensitivity_audits(a, b, id, "Q4")[, 2], 1 / 6), "Retained improvement must use six-item denominator.")
  b$score[b$subject_id == id & b$question_id == "Q3"] <- 0
  a$score[a$subject_id == id & a$question_id == "Q3"] <- 1
  check(near(score_sensitivity_audits(a, b, id, "Q4")[, 2], -1 / 6), "Pre-only correctness must yield negative change.")
}

# Mutations simulate common corruption/staleness risks. Every case must fail
# before fitting rather than silently adjusting the item denominator or cohort.
test_audit_failures <- function() {
  validate <- function(x) {
    path <- file.path(scratch, "corrupt_audit.csv")
    readr::write_csv(x, path)
    read_sensitivity_audit(path, invalid)
  }
  expect_error(validate(dplyr::bind_rows(raw_pre, raw_pre[1, ])), "Duplicate")
  expect_error(validate(raw_pre[-1, ]), "coverage")
  i <- which(raw_pre$excluded_from_retention_score == "False")[[1]]
  j <- which(raw_pre$excluded_from_retention_score == "True")[[1]]
  bad <- raw_pre; bad$score[i] <- NA
  expect_error(validate(bad), "binary valid scores")
  bad <- raw_pre; bad$score[i] <- "0.5"
  expect_error(validate(bad), "binary valid scores")
  bad <- raw_pre; bad$score[j] <- "0"
  expect_error(validate(bad), "missing excluded scores")
  bad <- raw_pre; bad$excluded_from_retention_score[j] <- "False"
  expect_error(validate(bad), "manifest")
  bad <- raw_pre; bad$key_answer[i] <- "changed key"
  expect_error(validate(bad), "Inconsistent item")
  bad <- raw_pre; bad$method[j] <- "exact_match"
  expect_error(validate(bad), "scoring method")
  expect_error(validate(raw_pre[raw_pre$question_id != "Q3", ]), "valid counts")
  bad <- raw_post; bad$key_answer[bad$question_id == "Q3"] <- "changed for everyone"
  path <- file.path(scratch, "mismatched_post.csv")
  readr::write_csv(bad, path)
  opts <- args; opts$audit_post_csv <- path
  expect_error(prepare_sensitivity(opts), "answer keys differ")
  opts <- args
  path <- file.path(scratch, "missing_participant.csv")
  readr::write_csv(raw_post[raw_post$Q34 != raw_post$Q34[[1]], ], path)
  opts$audit_post_csv <- path
  expect_error(prepare_sensitivity(opts), "participant coverage differs")
}

# A stale merged outcome must fail; a genuine missing merged condition preserves
# the primary complete-case rule across every subset, never restoring that subject.
test_merged_consistency <- function() {
  df <- readr::read_csv(args$input_csv, show_col_types = FALSE)
  id <- prepared$versions$baseline$subject_id[[1]]
  i <- which(retention_primary$normalize_subject_id(df$subject_id, "subject_id") == id)
  path <- file.path(scratch, "merged.csv")
  opts <- args; opts$input_csv <- path
  bad <- df; bad$diff_short_form_entertainment[i] <- bad$diff_short_form_entertainment[i] + 0.01
  readr::write_csv(bad, path)
  expect_error(prepare_sensitivity(opts), "do not reproduce")
  df$diff_short_form_entertainment[i] <- NA
  readr::write_csv(df, path)
  reduced <- prepare_sensitivity(opts)
  check(all(vapply(reduced$versions, function(x) !id %in% x$subject_id &&
        nrow(x) == prepared$n_complete - 1L, logical(1))), "Complete-case cohort not fixed.")
}

# Inject optimizer/diagnostic failures into a usable fitted model to verify that
# numerical failures cannot be reported as successful sensitivity fits. A singular
# boundary notification alone remains reportable, per the primary analysis policy.
test_model_diagnostics <- function() {
  data <- retention_primary$reshape_to_long(prepared$versions$baseline)
  model <- retention_primary$fit_factorial_lmm(data)
  check_sensitivity_fit(model)
  bad <- model; bad@optinfo$conv$opt <- 1L
  expect_error(check_sensitivity_fit(bad), "did not converge")
  bad <- model; bad@optinfo$conv$lme4$messages <- "Model failed to converge"
  expect_error(check_sensitivity_fit(bad), "did not converge")
  bad <- model; bad@optinfo$conv$lme4$messages <- "Unrecognized numerical diagnostic"
  expect_error(check_sensitivity_fit(bad), "did not converge")
  boundary <- model
  boundary@optinfo$conv$lme4$messages <- "boundary (singular) fit: see help('isSingular')"
  check_sensitivity_fit(boundary)
  captured <- capture_lmm_fit(function() { warning("Model failed to converge"); model })
  check(!captured$converged, "Convergence warnings must invalidate a run.")
  data$education_years <- data$age
  rank_deficient <- suppressMessages(retention_primary$fit_factorial_lmm(data))
  expect_error(check_sensitivity_fit(rank_deficient), "Rank-deficient")
}

# With complete balanced within-subject conditions and additive subject covariates,
# factorial coefficients equal mean within-subject contrasts. This independent
# reference detects coding/sign mistakes; manual Holm guards the three-test family.
test_end_to_end <- function() {
  result <- run_retention_sensitivity(args)
  check(nrow(result$main) == 24 && nrow(result$summary) == 3, "Incomplete fit reporting.")
  for (scenario in names(prepared$versions)) {
    d <- prepared$versions[[scenario]]
    y <- as.matrix(d[retention_primary$REQUIRED_DIFF_COLS])
    expected <- c(mean((y[, 3] + y[, 4] - y[, 1] - y[, 2]) / 2),
                  mean((y[, 1] + y[, 3] - y[, 2] - y[, 4]) / 2),
                  mean(y[, 3] - y[, 4] - y[, 1] + y[, 2]))
    r <- result$main[result$main$scenario == scenario, ]
    check(near(r$estimate, expected, 1e-10), "Independent factorial contrasts disagree.")
    o <- order(r$p_unc)
    expected_p <- pmin(1, cummax(r$p_unc[o] * c(3, 2, 1)))
    check(near(r$p_holm[o], expected_p), "Holm must use exactly three tests per fit.")
    check(near(r$ci95_low, r$estimate - qnorm(.975) * r$se), "Expected unadjusted Wald CIs.")
    check(all(r$n_subjects == prepared$n_complete & r$n_obs == 4 * prepared$n_complete), "Cohort changed.")
    posthoc <- readr::read_csv(file.path(args$out_dir, "scenarios", paste0(scenario, "_posthoc.csv")), show_col_types = FALSE)
    check(nrow(posthoc) == ifelse(r$p_holm[r$effect == "interaction"] < args$alpha, 6, 0), "Posthoc gate failed.")
  }
  # Regression against saved primary inference with explicit platform tolerances.
  saved <- readr::read_csv("data/results/retention_format_content_lmm_main_effects_r.csv", show_col_types = FALSE)
  b <- result$main[result$main$scenario == "baseline", ]
  saved <- saved[match(b$effect, saved$effect), ]
  for (col in c("estimate", "se", "df", "t", "p_unc", "ci95_low", "ci95_high", "eta2_p")) {
    tol <- if (col == "df") 1e-3 else if (col == "estimate") 1e-10 else 1e-7
    check(near(b[[col]], saved[[col]], tol), paste("Saved primary regression:", col))
  }
  check(near(b$p_holm, saved$p_fdr, 1e-7), "Saved adjusted p-values changed.")
  expect_error(run_retention_sensitivity(args), "new or empty")
  expect_error(parse_sensitivity_args(c("--min_subjects", "6.5")), "integer")
  expect_error(parse_sensitivity_args(c("--alpha", "0")), "between")
  expect_error(parse_sensitivity_args(c("--misspelled", "x")), "Unknown option")

  # Raise alpha only in this test to exercise the reporting path when all eight
  # interaction gates open. This is branch validation, not a scientific result.
  gate_args <- args
  gate_args$alpha <- .9999
  gate_args$out_dir <- file.path(scratch, "posthoc_branch_only")
  run_retention_sensitivity(gate_args)
  pooled <- readr::read_csv(file.path(gate_args$out_dir, "posthoc.csv"), show_col_types = FALSE)
  check(nrow(pooled) == 48L && all(table(pooled$scenario) == 6L),
        "Combined posthoc report must retain six contrasts for each gated scenario.")
}

test_complete_enumeration()
test_known_item_changes()
test_audit_failures()
test_merged_consistency()
test_model_diagnostics()
test_end_to_end()
cat("[PASS] Exhaustive sensitivity scoring, audit integrity, cohort, model and reporting checks.\n")
cat("[artifacts] ", scratch, "\n", sep = "")
