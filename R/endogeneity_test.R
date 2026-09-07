## endogeneity_test(): are the regressors endogenous at all?
##
## Hou, Ramalho and Roseta-Palma (2025), Energy Economics 151:108922, Eq. (25).
## They compare five estimators and note that each of the endogeneity-correcting
## ones carries a "correction parameter" that is zero under exogeneity, so a
## Wald test on it is a Durbin-Wu-Hausman-style test for whether the correction
## was needed. See notes/code_history/endogeneity_test.md.
##
## In this package the correction parameter is `rho`, the correlation between
## the noise v and the reduced-form errors xi:
##
##   H0: rho = 0   (regressors exogenous, plain sfm() is consistent)
##   H1: rho != 0  (they are not, and ivsfm() is doing real work)
##
## The null is INTERIOR -- rho lives in the open unit ball -- so this is an
## ordinary chi-square, not the chi-bar-square mixture inefficiency_test()
## needs for a null on the boundary.

endogeneity_test <- function(object, level = 0.05) {
  if (!inherits(object, "sfareg")) {
    stop("`object` must be a fit from ivsfm().", call. = FALSE)
  }
  mn <- object$model_name
  if (is.null(object$rho) || is.null(object$vcov_rho)) {
    stop("endogeneity_test() needs a fit from ivsfm(model_name = \"IVLIML\") ",
      "or \"IVCF\"", if (!is.null(mn)) paste0("; this is a ", mn, " fit") else "",
      ". \"C2SLS\" is a moment estimator with no correction parameter to test, ",
      "and the other entry points do not model endogeneity at all.",
      call. = FALSE
    )
  }
  rho <- as.numeric(object$rho)
  V <- as.matrix(object$vcov_rho)
  m <- length(rho)
  nm <- colnames(V)
  if (is.null(nm)) nm <- paste0("rho_", seq_len(m))

  if (!all(is.finite(V))) {
    stop("endogeneity_test(): this fit has no usable covariance for rho ",
      "(the Hessian was not available or was not invertible), so the Wald ",
      "statistic cannot be formed. Refit with optHessian = TRUE.",
      call. = FALSE
    )
  }

  ## Wald, Eq. (25): W = rho' Var(rho)^{-1} rho, asymptotically chi2(m) with m
  ## the number of endogenous regressors. solve() rather than a pseudo-inverse:
  ## a singular Var(rho) means the correction is not identified, and inverting
  ## it anyway would return a statistic that looks fine and means nothing.
  Vi <- tryCatch(solve(V), error = function(e) NULL)
  if (is.null(Vi)) {
    stop("endogeneity_test(): the estimated covariance of rho is singular, so ",
      "the correction parameters are not jointly identified. Check the ",
      "instruments.",
      call. = FALSE
    )
  }
  W <- as.numeric(t(rho) %*% Vi %*% rho)
  p <- stats::pchisq(W, df = m, lower.tail = FALSE)

  se <- sqrt(diag(V))
  out <- list(
    statistic = c("W" = W),
    parameter = c("df" = m),
    p.value = p,
    method = paste0(
      "Wald test for endogeneity in a stochastic frontier\n",
      "H0: the noise is uncorrelated with the reduced-form errors (rho = 0)"
    ),
    data.name = deparse(object$call$formula),
    alternative = "at least one regressor is endogenous",
    rho = stats::setNames(rho, nm),
    se = stats::setNames(se, nm),
    z = stats::setNames(rho / se, nm),
    model_name = mn,
    endogenous = object$endogenous,
    nobs = object$nobs,
    level = level,
    reject = isTRUE(p < level)
  )
  class(out) <- c("sfa_endog_test", "htest")
  out
}

print.sfa_endog_test <- function(x, ...) {
  cat("\n", x$method, "\n\n", sep = "")
  cat("data:  ", x$data.name, "  (", x$model_name, ", n = ", x$nobs, ")\n", sep = "")
  cat(sprintf("W = %.4f, df = %d, p-value = %.4g\n", x$statistic, x$parameter, x$p.value))
  cat("alternative hypothesis: ", x$alternative, "\n\n", sep = "")
  tab <- cbind(rho = x$rho, `std. error` = x$se, z = x$z)
  print(round(tab, 4))
  cat("\n")
  if (x$reject) {
    cat("  Reject exogeneity at the ", 100 * x$level, "% level: the correction\n",
      "  is doing work, and an uncorrected sfm() fit is inconsistent here.\n",
      sep = ""
    )
  } else {
    cat("  Exogeneity is not rejected at the ", 100 * x$level, "% level. That is\n",
      "  not evidence FOR exogeneity -- with weak instruments rho is estimated\n",
      "  imprecisely and this test has little power.\n",
      sep = ""
    )
  }
  cat("\n")
  invisible(x)
}
