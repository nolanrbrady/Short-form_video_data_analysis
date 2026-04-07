# Shared mixed-model convergence helpers for the neural LMM scripts.
#
# Bates et al. (2015) emphasize checking optimizer/convergence diagnostics for
# mixed models rather than assuming that a returned fit is numerically clean.

is_lmm_nonconvergence_warning <- function(message) {
  msg <- tolower(as.character(message))
  patterns <- c(
    "failed to converge",
    "degenerate .* hessian",
    "unable to evaluate scaled gradient",
    "unable to evaluate hessian",
    "max\\|grad\\|",
    "negative eigenvalue",
    "very large eigenvalue",
    "nearly unidentifiable",
    "non-positive-definite"
  )
  any(vapply(patterns, function(pattern) grepl(pattern, msg, perl = TRUE), logical(1)))
}

capture_lmm_fit <- function(fit_fn) {
  warning_messages <- character()
  nonconvergence_messages <- character()

  model <- withCallingHandlers(
    fit_fn(),
    warning = function(w) {
      msg <- conditionMessage(w)
      warning_messages <<- c(warning_messages, msg)
      if (is_lmm_nonconvergence_warning(msg)) {
        nonconvergence_messages <<- c(nonconvergence_messages, msg)
      }
      invokeRestart("muffleWarning")
    }
  )

  list(
    model = model,
    converged = length(nonconvergence_messages) == 0,
    warning_messages = unique(warning_messages),
    nonconvergence_messages = unique(nonconvergence_messages)
  )
}
