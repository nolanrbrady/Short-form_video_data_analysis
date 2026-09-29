#!/usr/bin/env Rscript

# Exhaustive equal-question-count sensitivity for the existing retention LMM.
# Run from the repository root, after recall scoring and the merge pipeline.
# The scored audits, not raw response text, are the inputs: this does not regrade.
#
# All seven ways of retaining six Short-form Entertainment items are evaluated;
# each omission is shared across participants and pre/post. The other conditions
# already contain six items. These overlapping subsets are not independent
# replications, bootstrap samples, or permutations of a null hypothesis.
# Item-count balance cannot remove stimulus-selection confounding or establish
# generalization over stimuli: Judd, Westfall & Kenny (2012),
# https://doi.org/10.1037/a0028347. For the distinct exchangeability requirements
# of permutation inference, see Winkler et al. (2014),
# https://doi.org/10.1016/j.neuroimage.2014.01.060. See CITATIONS.md.
#
# Reuse the primary REML/Satterthwaite model, unadjusted Wald CIs, three-effect
# Holm family, and interaction-gated posthoc contrasts without redefining them:
# Bates et al. (2015), https://doi.org/10.18637/jss.v067.i01;
# Kuznetsova et al. (2017), https://doi.org/10.18637/jss.v082.i13;
# Holm (1979), https://doi.org/10.2307/4615733.

retention_primary <- new.env(parent = globalenv())
sys.source("analyze_retention_format_content_lmm.R", envir = retention_primary)
source("r_lmm_convergence_helpers.R")

SENSITIVITY_CONDITIONS <- c("Short-form Education", "Short-form Entertainment",
                          "Long-form Education", "Long-form Entertainment")
SENSITIVITY_COUNTS <- c(6L, 7L, 6L, 6L)
SENSITIVITY_TARGET <- "Short-form Entertainment"

# Reject unknown options and ambiguous thresholds rather than silently changing
# the analysis. Item-count design is deliberately fixed; changing it needs review.
parse_sensitivity_args <- function(argv = commandArgs(trailingOnly = TRUE)) {
  args <- list(
    input_csv = "data/tabular/generated_data/homer3_betas_plus_combined_sfv_data_inner_join.csv",
    audit_pre_csv = "demographic/recall_assessment_audit_pre.csv",
    audit_post_csv = "demographic/recall_assessment_audit_post.csv",
    invalid_questions_json = "data/config/recall_invalid_questions.json",
    exclude_subjects_json = "data/config/excluded_subjects.json",
    out_dir = "data/results/retention_sensitivity",
    alpha = 0.05, min_subjects = 6L
  )
  if (length(argv) %% 2 != 0) stop("Expected --key value pairs.")
  for (i in seq.int(1L, length(argv) + 1L, by = 2L)) {
    if (i > length(argv)) break
    key <- sub("^--", "", argv[[i]])
    if (!startsWith(argv[[i]], "--") || !key %in% names(args)) stop("Unknown option: ", argv[[i]])
    args[[key]] <- argv[[i + 1L]]
  }
  args$alpha <- suppressWarnings(as.numeric(args$alpha))
  n <- suppressWarnings(as.numeric(args$min_subjects))
  if (!is.finite(args$alpha) || args$alpha <= 0 || args$alpha >= 1) stop("alpha must be between 0 and 1.")
  if (!is.finite(n) || n < 6 || n != floor(n)) stop("min_subjects must be an integer >= 6.")
  args$min_subjects <- as.integer(n)
  args
}

# Validate every participant, including excluded participants, before subsetting.
# Explicit binary scores preserve upstream blank-response=incorrect scoring;
# missing valid scores are errors, never silently omitted from denominators.
read_sensitivity_audit <- function(path, invalid) {
  audit <- readr::read_csv(path, col_types = readr::cols(.default = "c"),
                          na = c("", "NA"), show_col_types = FALSE)
  if (nrow(readr::problems(audit))) stop("Malformed audit CSV: ", path)
  retention_primary$assert_required_columns(audit, c("Q34", "question_id", "condition",
    "key_answer", "score", "method", "excluded_from_retention_score"), path)
  if (!nrow(audit)) stop("Empty audit: ", path)
  audit$subject_id <- retention_primary$normalize_subject_id(audit$Q34, "Q34")
  if (anyNA(audit$question_id) || anyNA(audit$key_answer) ||
      any(!audit$condition %in% SENSITIVITY_CONDITIONS)) stop("Missing/unknown audit item, key, or condition.")
  if (anyDuplicated(audit[c("subject_id", "question_id")])) stop("Duplicate participant/question rows.")
  flags <- tolower(audit$excluded_from_retention_score)
  if (anyNA(flags) || any(!flags %in% c("true", "false"))) stop("Invalid exclusion flags.")
  audit$excluded <- flags == "true"
  manifest_match <- match(audit$question_id, invalid$question_id)
  expected <- !is.na(manifest_match)
  if (any(audit$excluded != expected)) stop("Audit exclusion flags disagree with invalid-question manifest.")
  if (any(audit$condition[expected] != invalid$condition[manifest_match[expected]])) {
    stop("Invalid-question condition disagrees with manifest.")
  }
  raw_scores <- audit$score
  audit$score <- suppressWarnings(as.numeric(raw_scores))
  if (any(!is.na(raw_scores[audit$excluded])) ||
      any(!is.finite(audit$score[!audit$excluded])) ||
      any(!audit$score[!audit$excluded] %in% c(0, 1))) stop("Expected missing excluded scores and binary valid scores.")
  if (anyNA(audit$method) || any((audit$method == "excluded_invalid_question") != audit$excluded)) {
    stop("Audit scoring method disagrees with exclusion status.")
  }
  # Each question must have one key/condition and the same participant coverage.
  dictionary <- unique(audit[c("question_id", "condition", "key_answer", "excluded")])
  if (anyDuplicated(dictionary$question_id)) stop("Inconsistent item keys or conditions across participants.")
  ids <- unique(audit$subject_id)
  if (any(table(audit$question_id) != length(ids))) stop("Incomplete participant/item coverage in audit.")
  counts <- table(factor(dictionary$condition[!dictionary$excluded], levels = SENSITIVITY_CONDITIONS))
  originals <- table(factor(dictionary$condition, levels = SENSITIVITY_CONDITIONS))
  if (!identical(as.integer(counts), SENSITIVITY_COUNTS) || any(originals != 8L)) {
    stop("Expected 8 original items per condition and valid counts 6,7,6,6; review sensitivity design.")
  }
  audit
}

# Reconstruct post-minus-pre using a common item set within each phase. Long-form
# Qualtrics IDs differ by phase; only the resampled condition requires paired IDs.
score_sensitivity_audits <- function(pre, post, ids, omitted = character()) {
  phase_means <- function(audit) {
    matrix(vapply(SENSITIVITY_CONDITIONS, function(condition) {
      keep <- !audit$excluded & audit$condition == condition &
        !(audit$condition == SENSITIVITY_TARGET & audit$question_id %in% omitted)
      values <- tapply(audit$score[keep], audit$subject_id[keep], mean)
      as.numeric(values[as.character(ids)])
    }, numeric(length(ids))), nrow = length(ids), ncol = 4L)
  }
  out <- phase_means(post) - phase_means(pre)
  if (any(!is.finite(out))) stop("Missing participant or nonfinite reconstructed score.")
  colnames(out) <- retention_primary$REQUIRED_DIFF_COLS
  out
}

# Fail on stale audits/merged scores before any fitting. Match the primary
# exclusions and complete-case cohort once and hold that cohort fixed in all fits.
prepare_sensitivity <- function(args) {
  invalid <- jsonlite::fromJSON(args$invalid_questions_json)$invalid_questions
  if (!is.data.frame(invalid) || !all(c("question_id", "condition") %in% names(invalid)) ||
      anyNA(invalid$question_id) || anyDuplicated(invalid$question_id) ||
      any(!invalid$condition %in% SENSITIVITY_CONDITIONS)) stop("Invalid question-exclusion manifest.")
  pre <- read_sensitivity_audit(args$audit_pre_csv, invalid)
  post <- read_sensitivity_audit(args$audit_post_csv, invalid)
  if (!setequal(pre$subject_id, post$subject_id)) stop("Pre/post participant coverage differs.")
  dictionary <- function(audit) {
    x <- unique(audit[!audit$excluded & audit$condition == SENSITIVITY_TARGET,
                      c("question_id", "key_answer")])
    x[order(x$question_id), ]
  }
  if (!identical(dictionary(pre), dictionary(post))) stop("Resampled pre/post item IDs or answer keys differ.")
  questions <- dictionary(pre)$question_id
  if (any(!grepl("^Q[0-9]+$", questions))) stop("Unexpected resampled question IDs.")
  questions <- questions[order(as.numeric(sub("^Q", "", questions)))]
  df <- retention_primary$load_retention_input(args$input_csv, args$exclude_subjects_json)
  if (any(!df$subject_id %in% pre$subject_id)) stop("Merged participant missing from audits.")
  baseline <- score_sensitivity_audits(pre, post, df$subject_id)
  observed <- as.matrix(df[retention_primary$REQUIRED_DIFF_COLS])
  if (any(is.infinite(observed)) || any(abs(baseline - observed) > 1e-12, na.rm = TRUE)) {
    stop("Audit scores do not reproduce merged retention scores (tolerance 1e-12).")
  }
  cc_ids <- unique(retention_primary$complete_case_subjects(retention_primary$reshape_to_long(df))$subject_id)
  cohort <- df[df$subject_id %in% cc_ids, c("subject_id", "age", "education_years",
                                          retention_primary$REQUIRED_DIFF_COLS)]
  if (nrow(cohort) < args$min_subjects) stop("Insufficient complete-case subjects.")
  if (any(!is.finite(as.matrix(cohort[c("age", "education_years")])))) stop("Nonfinite covariates.")
  versions <- list(baseline = cohort)
  for (question in questions) {
    version <- cohort
    scores <- score_sensitivity_audits(pre, post, cohort$subject_id, question)
    # Keep the other three outcomes and all covariates byte-for-byte unchanged.
    version$diff_short_form_entertainment <- scores[, "diff_short_form_entertainment"]
    versions[[paste0("omit_", question)]] <- version
  }
  # This identity follows from equal item inclusion across all seven subsets; it
  # is an arithmetic integrity check, not evidence against stimulus-related bias.
  subset_mean <- rowMeans(vapply(versions[-1], function(x) x$diff_short_form_entertainment,
                                 numeric(nrow(cohort))))
  if (any(abs(subset_mean - cohort$diff_short_form_entertainment) > 1e-12)) stop("Subset mean identity failed.")
  list(versions = versions, questions = questions,
       n_after_exclusions = nrow(df), n_complete = nrow(cohort))
}

# Require usable fits. Singularity is explicitly reported, not silently discarded.
# Numerical convergence and design rank are checked in addition to warnings:
# Bates et al. (2015), https://doi.org/10.18637/jss.v067.i01.
check_sensitivity_fit <- function(model) {
  if (is.null(model)) return(invisible(NULL))
  conv <- model@optinfo$conv
  messages <- unlist(conv$lme4$messages)
  # lme4 also stores singular-boundary notifications here; retain and report them.
  unexpected_messages <- messages[!grepl("boundary.*singular", messages, ignore.case = TRUE)]
  if (any(unlist(conv$opt) != 0) ||
      length(unexpected_messages) > 0L ||
      any(vapply(messages, is_lmm_nonconvergence_warning, logical(1)))) stop("Model did not converge.")
  if (length(attr(lme4::getME(model, "X"), "col.dropped"))) stop("Rank-deficient model design.")
  invisible(NULL)
}

# Summarize the full range without selecting a favorable deletion or interpreting
# the fraction significant as a probability. Each scenario keeps its own family
# of three Holm tests, exactly as in the primary analysis (Holm, 1979).
summarize_sensitivity <- function(results, alpha) {
  dplyr::bind_rows(lapply(unique(results$effect), function(effect) {
    b <- results[results$effect == effect & results$scenario == "baseline", ]
    s <- results[results$effect == effect & results$scenario != "baseline", ]
    tibble::tibble(effect = effect, baseline_estimate = b$estimate,
      baseline_p_holm = b$p_holm, subset_estimate_min = min(s$estimate),
      subset_estimate_max = max(s$estimate), subset_p_holm_min = min(s$p_holm),
      subset_p_holm_max = max(s$p_holm), n_subsets = nrow(s),
      n_significant = sum(s$p_holm < alpha),
      n_same_direction = sum(sign(s$estimate) == sign(b$estimate)),
      n_same_significance = sum((s$p_holm < alpha) == (b$p_holm < alpha)),
      n_singular = sum(s$singular_fit), alpha = alpha)
  }))
}

# Stage the entire run before publishing results to a new/empty output directory.
# Never overwrite primary inputs, saved primary results, or an older sensitivity run.
run_retention_sensitivity <- function(args = parse_sensitivity_args()) {
  if (file.exists(args$out_dir) && (!dir.exists(args$out_dir) ||
      length(list.files(args$out_dir, all.files = TRUE, no.. = TRUE)))) {
    stop("out_dir must be new or empty; choose a new directory to preserve prior results.")
  }
  prepared <- prepare_sensitivity(args)
  staging <- tempfile("retention-sensitivity-")
  dir.create(staging)
  on.exit(unlink(staging, recursive = TRUE), add = TRUE)
  dir.create(file.path(staging, "scenarios"))
  results <- posthocs <- scores <- diagnostics <- list()
  for (scenario in names(prepared$versions)) {
    cat("[sensitivity] ", scenario, "\n", sep = "")
    prefix <- file.path(staging, "scenarios", scenario)
    input <- paste0(prefix, "_input.csv")
    readr::write_csv(prepared$versions[[scenario]], input)
    opts <- list(input_csv = input, exclude_subjects_json = args$exclude_subjects_json,
      out_main_csv = paste0(prefix, "_main.csv"), out_posthoc_csv = paste0(prefix, "_posthoc.csv"),
      alpha = args$alpha, min_subjects = args$min_subjects,
      pbkrtest_limit = NA_real_, lmerTest_limit = NA_real_)
    fit <- capture_lmm_fit(function() retention_primary$main(opts))
    if (!fit$converged) stop("Nonconvergence in ", scenario, ": ", paste(fit$nonconvergence_messages, collapse = "; "))
    value <- fit$model
    check_sensitivity_fit(value$model)
    check_sensitivity_fit(value$posthoc_model)
    if (any(!is.finite(as.matrix(value$main[c("estimate", "se", "df", "t", "p_unc", "p_fdr",
                                             "ci95_low", "ci95_high", "eta2_p")])))) stop("Nonfinite inference in ", scenario)
    if (any(value$main$n_subjects != prepared$n_complete) ||
        any(value$main$n_obs != 4 * prepared$n_complete)) stop("Model cohort changed.")
    results[[scenario]] <- dplyr::mutate(value$main, scenario = scenario, .before = 1)
    posthocs[[scenario]] <- dplyr::mutate(value$posthoc, scenario = rep(scenario, nrow(value$posthoc)), .before = 1)
    scores[[scenario]] <- dplyr::mutate(prepared$versions[[scenario]], scenario = scenario, .before = 1)
    diagnostics[[scenario]] <- tibble::tibble(scenario = scenario, converged = TRUE,
      singular_fit = lme4::isSingular(value$model, tol = 1e-4),
      posthoc_run = !is.null(value$posthoc_model),
      posthoc_singular_fit = if (is.null(value$posthoc_model)) NA else lme4::isSingular(value$posthoc_model, tol = 1e-4),
      warnings = paste(fit$warning_messages, collapse = " | "))
  }
  combined <- dplyr::rename(dplyr::bind_rows(results), p_holm = p_fdr)
  summary <- summarize_sensitivity(combined, args$alpha)
  readr::write_csv(combined, file.path(staging, "main_effects.csv"))
  readr::write_csv(summary, file.path(staging, "summary.csv"))
  readr::write_csv(dplyr::bind_rows(posthocs), file.path(staging, "posthoc.csv"))
  readr::write_csv(dplyr::bind_rows(scores), file.path(staging, "subject_scores.csv"))
  readr::write_csv(dplyr::bind_rows(diagnostics), file.path(staging, "diagnostics.csv"))
  plan <- tibble::tibble(scenario = names(prepared$versions),
    omitted_question = c(NA_character_, prepared$questions),
    retained_short_entertainment = vapply(c(NA_character_, prepared$questions), function(q) {
      paste(prepared$questions[is.na(q) | prepared$questions != q], collapse = ";")
    }, character(1)), n_short_education = 6L, n_short_entertainment = c(7L, rep(6L, 7)),
    n_long_education = 6L, n_long_entertainment = 6L)
  readr::write_csv(plan, file.path(staging, "omission_plan.csv"))
  provenance_paths <- c(unlist(args[c("input_csv", "audit_pre_csv", "audit_post_csv",
    "invalid_questions_json", "exclude_subjects_json")]),
    "analyze_retention_sensitivity.R", "analyze_retention_format_content_lmm.R",
    "r_subject_exclusions.R", "r_emmeans_posthoc_helpers.R", "r_lmm_convergence_helpers.R")
  jsonlite::write_json(list(arguments = args, n_after_exclusions = prepared$n_after_exclusions,
    n_complete = prepared$n_complete, input_and_code_md5 = as.list(tools::md5sum(provenance_paths)),
    interpretation = paste("Exhaustive common-item subset sensitivity, not a permutation test.",
      "No multiplicity correction across overlapping sensitivity scenarios; Holm within each three-effect family.",
      "Wald CIs are unadjusted. Equal counts do not remove missing-stimulus confounding.",
      "Scored audits are reused, not regraded.")),
    file.path(staging, "metadata.json"), pretty = TRUE, auto_unbox = TRUE)
  writeLines(capture.output(sessionInfo()), file.path(staging, "session_info.txt"))
  dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)
  if (!all(file.copy(list.files(staging, full.names = TRUE), args$out_dir, recursive = TRUE))) stop("Failed to publish sensitivity outputs.")
  cat("[write] sensitivity results: ", args$out_dir, "\n", sep = "")
  print(summary)
  invisible(list(main = combined, summary = summary, prepared = prepared))
}

if (sys.nframe() == 0) {
  run_retention_sensitivity()
}
