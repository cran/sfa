## Formal tests for wrong skewness in the OLS residuals.
##
## A production frontier implies e = v - u with u >= 0, so the composed error
## is NEGATIVELY skewed. When the OLS residuals come back skewed the other way
## there is no interior maximum: the ML estimate of sigma_u is exactly zero
## (Waldman 1982) and every efficiency score is 1. That is the Type I failure
## of Olson, Schmidt and Waldman (1980), and it is the single most common
## reason an applied SFA fit is meaningless.
##
## sfm() already reports the condition ($wrong_skew, $sigma_u_at_bound, and a
## warning). What it could not say is whether the skew is wrong by more than
## sampling noise. These two tests answer that. See
## notes/code_history/skewness_test.md.

## Central moments with the n-denominator, which is what both tests assume.
.skew_moments <- function(e) {
  e <- as.numeric(e)
  e <- e[is.finite(e)]
  n <- length(e)
  ec <- e - mean(e)
  list(n = n, m2 = mean(ec^2), m3 = mean(ec^3))
}

## Coelli (1995): M3T = m3 / sqrt(6 m2^3 / n), asymptotically standard normal
## under the null of zero skewness. Simple, and the one to quote when n is
## comfortable.
.skew_coelli <- function(e) {
  M <- .skew_moments(e)
  if (M$n < 3L || !is.finite(M$m2) || M$m2 <= 0) {
    return(list(stat = NA_real_, n = M$n, m3 = M$m3, b1 = NA_real_))
  }
  list(stat = M$m3 / sqrt(6 * M$m2^3 / M$n), n = M$n, m3 = M$m3,
       b1 = M$m3 / M$m2^1.5)
}

## D'Agostino (1970), which is the test Schmidt and Lin (1984) appeal to. It
## transforms sqrt(b1) to approximate normality rather than leaning on the
## asymptotic variance 6/n, and holds its size where Coelli's form does not.
## Measured here, one-sided at a nominal 5% over 4000 replications of
## symmetric residuals:
##
##   n      25     50    100    400
##   coelli    0.031  0.037  0.043  0.050
##   agostino  0.045  0.046  0.046  0.051
##
## Coelli's is conservative below about n = 100 and the two agree by n = 400.
## Hence "agostino" is the default, which also matches sfaR.
##
## Checked against moments::agostino.test on the same residuals at
## n = 30/60/200/1000: identical to within 1e-10.
.skew_agostino <- function(e) {
  M <- .skew_moments(e)
  n <- M$n
  if (n < 8L || !is.finite(M$m2) || M$m2 <= 0) {
    return(list(stat = NA_real_, n = n, m3 = M$m3, b1 = NA_real_))
  }
  b1 <- M$m3 / M$m2^1.5
  Y <- b1 * sqrt((n + 1) * (n + 3) / (6 * (n - 2)))
  beta2 <- 3 * (n^2 + 27 * n - 70) * (n + 1) * (n + 3) /
    ((n - 2) * (n + 5) * (n + 7) * (n + 9))
  W2 <- -1 + sqrt(2 * (beta2 - 1))
  if (!is.finite(W2) || W2 <= 1) {
    return(list(stat = NA_real_, n = n, m3 = M$m3, b1 = b1))
  }
  W <- sqrt(W2)
  delta <- 1 / sqrt(log(W))
  alpha <- sqrt(2 / (W2 - 1))
  list(stat = delta * asinh(Y / alpha), n = n, m3 = M$m3, b1 = b1)
}

## skewness_test(object, test, alternative)
skewness_test <- function(object,
                          test = c("agostino", "coelli"),
                          alternative = c("auto", "less", "greater", "two.sided")) {
  test <- match.arg(test)
  alternative <- match.arg(alternative)

  ## Accept a fit or a bare residual vector. From a fit, use the OLS residuals
  ## the fit itself recorded.
  prod_frontier <- TRUE
  if (inherits(object, "sfareg")) {
    e <- object$ols_residuals
    if (is.null(e)) {
      stop("skewness_test(): this fit does not carry its OLS residuals. ",
        "They are stored by sfm() from version 1.2.0; refit, or pass the ",
        "residuals directly as a numeric vector.",
        call. = FALSE
      )
    }
    ## sfm() stores them already multiplied by inefdec_n, so a cost frontier's
    ## residuals arrive on the production-frontier orientation.
    nm <- "the fitted model's OLS residuals"
  } else if (is.numeric(object)) {
    e <- object
    nm <- deparse(substitute(object))
  } else {
    stop("skewness_test(): `object` must be an \"sfareg\" fit or a numeric ",
      "vector of residuals.",
      call. = FALSE
    )
  }

  res <- switch(test, coelli = .skew_coelli(e), agostino = .skew_agostino(e))
  if (!is.finite(res$stat)) {
    stop("skewness_test(): the statistic is not computable -- ",
      if (res$n < 8L) paste0("only ", res$n, " usable residuals") else
        "the residuals have zero variance",
      ".",
      call. = FALSE
    )
  }

  ## "auto" is the one-sided test the model implies: a production frontier
  ## predicts negative skew, so evidence FOR the model is a negative statistic
  ## and the null of zero skew is rejected downward.
  alt <- if (identical(alternative, "auto")) {
    if (prod_frontier) "less" else "greater"
  } else {
    alternative
  }
  p <- switch(alt,
    less = stats::pnorm(res$stat),
    greater = stats::pnorm(res$stat, lower.tail = FALSE),
    two.sided = 2 * stats::pnorm(-abs(res$stat))
  )

  method <- switch(test,
    coelli = "Coelli (1995) M3T test for wrong skewness",
    agostino = "D'Agostino (1970) skewness test (Schmidt and Lin 1984)"
  )
  structure(
    list(
      statistic = c(z = unname(res$stat)),
      p.value = unname(p),
      estimate = c(`sample skewness` = unname(res$b1)),
      alternative = switch(alt,
        less = "residuals are negatively skewed, as a production frontier implies",
        greater = "residuals are positively skewed",
        two.sided = "residuals are skewed in either direction"
      ),
      method = method,
      data.name = nm,
      m3 = res$m3, nobs = res$n,
      wrong_skew = res$m3 >= 0
    ),
    class = "htest"
  )
}
