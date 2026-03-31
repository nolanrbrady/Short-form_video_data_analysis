#!/usr/bin/env Rscript

source("r_lmm_convergence_helpers.R")

assert_true <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

main <- function() {
  positive_examples <- c(
    "Model failed to converge with 1 negative eigenvalue: -7.0e+01",
    "Model failed to converge with max|grad| = 0.013 (tol = 0.002, component 1)",
    "unable to evaluate scaled gradient",
    "Model is nearly unidentifiable: very large eigenvalue"
  )
  for (msg in positive_examples) {
    assert_true(is_lmm_nonconvergence_warning(msg), paste0("expected non-convergence classification for: ", msg))
  }

  negative_examples <- c(
    "boundary (singular) fit: see help('isSingular')",
    "some unrelated warning",
    "NaNs produced"
  )
  for (msg in negative_examples) {
    assert_true(!is_lmm_nonconvergence_warning(msg), paste0("unexpected non-convergence classification for: ", msg))
  }

  fit_result <- capture_lmm_fit(function() {
    warning("Model failed to converge with 1 negative eigenvalue: -7.0e+01")
    warning("Model failed to converge with 1 negative eigenvalue: -7.0e+01")
    warning("boundary (singular) fit: see help('isSingular')")
    structure(list(dummy = TRUE), class = "dummy_fit")
  })

  assert_true(identical(class(fit_result$model), "dummy_fit"), "expected captured model to be returned")
  assert_true(!fit_result$converged, "expected converged=FALSE when a non-convergence warning is captured")
  assert_true(length(fit_result$warning_messages) == 2, "expected duplicate warnings to be deduplicated")
  assert_true(length(fit_result$nonconvergence_messages) == 1, "expected one deduplicated non-convergence warning")

  clean_result <- capture_lmm_fit(function() structure(list(dummy = TRUE), class = "dummy_fit"))
  assert_true(clean_result$converged, "expected converged=TRUE for clean fit")
  assert_true(length(clean_result$warning_messages) == 0, "expected no warnings for clean fit")

  cat("[PASS] LMM convergence helper validation passed\n")
}

main()
