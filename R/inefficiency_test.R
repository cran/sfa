## Coelli (1995), Journal of Productivity Analysis 6:247-268.
##
## Section 3's tests of H0: no technical inefficiency, and the appendix's
## variance for the COLS estimator. See notes/code_history/inefficiency_test.md.
##
## The point of the paper for us is that the two tests everybody actually runs
## -- the Wald ratio and the ordinary LR test -- have the WRONG SIZE, because
## H0 puts gamma on the boundary of the parameter space. Coelli's Monte Carlo
## shows the one-sided LR test and the third-moment test are the ones with
## correct size, and that the one-sided LR has the better power of the two.
## skewness_test(test = "coelli") already provides the third moment test; this
## file adds the likelihood-based ones so the comparison can be made in-package.

## Battese-Corra gamma = sigma_u^2 / (sigma_u^2 + sigma_v^2), which is what
## Coelli tests, from the (lambda, sigma) pair sfm() reports. lambda =
## sigma_u/sigma_v, so gamma = lambda^2/(1 + lambda^2) and the delta-method
## derivative is 2*lambda/(1 + lambda^2)^2.
.coelli_gamma <- function(object) {
  cf <- object$out[, "par"]
  if (!"lambda" %in% names(cf)) {
    return(NULL)
  }
  lam <- unname(cf[["lambda"]])
  g <- lam^2 / (1 + lam^2)
  se <- NA_real_
  V <- tryCatch(stats::vcov(object), error = function(e) NULL)
  if (!is.null(V) && "lambda" %in% rownames(V)) {
    d <- 2 * lam / (1 + lam^2)^2
    v <- V["lambda", "lambda"]
    if (is.finite(v) && v >= 0) se <- sqrt(d^2 * v)
  }
  list(gamma = g, se = se, lambda = lam)
}

## Gaussian log-likelihood of the RESTRICTED model, which under H0: gamma = 0
## is just the OLS fit. Computed from the residuals sfm() already stored, so
## this never refits anything.
.coelli_ll0 <- function(e) {
  n <- length(e)
  s2 <- sum((e - mean(e))^2) / n
  -0.5 * n * (log(2 * pi) + log(s2) + 1)
}

## p-value for a LR statistic whose null distribution is the Gourieroux, Holly
## and Monfort (1982) chi-bar-square mixture. `q` is the number of restrictions
## that sit ON the boundary.
##   q = 1 (gamma = 0, half-normal):        0.5 chi2_0 + 0.5 chi2_1
##   q = 2 (gamma = 0 and mu = 0, trunc.):  0.25 chi2_0 + 0.5 chi2_1 + 0.25 chi2_2
## chi2_0 is a point mass at zero, so it contributes nothing to an upper-tail
## probability at any positive statistic.
.chi_bar_p <- function(stat, q) {
  if (!is.finite(stat) || stat <= 0) {
    return(1)
  }
  if (q == 1L) {
    0.5 * stats::pchisq(stat, df = 1, lower.tail = FALSE)
  } else {
    0.5 * stats::pchisq(stat, df = 1, lower.tail = FALSE) +
      0.25 * stats::pchisq(stat, df = 2, lower.tail = FALSE)
  }
}

#' Test whether there is any technical inefficiency
#'
#' @param object An `"sfareg"` fit from [sfm()].
#' @param test Which statistics to compute.
#' @param level Significance level used for the `reject` column.
#' @return A data frame, one row per test.
#' @export
inefficiency_test <- function(object,
                              test = c("lr_1sided", "lr", "wald", "m3t"),
                              level = 0.05) {
  if (!inherits(object, "sfareg")) {
    stop("`object` must be an \"sfareg\" fit.", call. = FALSE)
  }
  test <- match.arg(test, several.ok = TRUE)
  mn <- object$model_name
  if (!identical(mn, "NHN") && !identical(mn, "NTN")) {
    stop("inefficiency_test(): Coelli (1995) is derived for the half-normal ",
      "model, and the chi-bar-square mixture below is stated for half-normal ",
      "(one boundary restriction) and truncated normal (two). This fit is ",
      "model_name = ", dQuote(mn), ".",
      call. = FALSE
    )
  }
  q <- if (identical(mn, "NTN")) 2L else 1L

  e <- object$ols_residuals
  if (is.null(e)) {
    stop("inefficiency_test(): this fit did not store its OLS residuals.", call. = FALSE)
  }
  n <- length(e)
  ll1 <- tryCatch(as.numeric(stats::logLik(object)), error = function(err) NA_real_)
  ll0 <- .coelli_ll0(e)
  LR <- 2 * (ll1 - ll0)

  rows <- list()

  if ("lr_1sided" %in% test) {
    rows[[length(rows) + 1L]] <- data.frame(
      test = "LR (one-sided)", statistic = LR,
      null = if (q == 1L) "0.5 chi2(0) + 0.5 chi2(1)" else "0.25/0.5/0.25 chi2(0,1,2)",
      p.value = .chi_bar_p(LR, q), stringsAsFactors = FALSE
    )
  }
  if ("lr" %in% test) {
    rows[[length(rows) + 1L]] <- data.frame(
      test = "LR (naive, two-sided)", statistic = LR,
      null = paste0("chi2(", q, ")"),
      p.value = stats::pchisq(LR, df = q, lower.tail = FALSE), stringsAsFactors = FALSE
    )
  }
  if ("wald" %in% test) {
    g <- .coelli_gamma(object)
    W <- if (is.null(g) || !is.finite(g$se) || g$se <= 0) NA_real_ else g$gamma / g$se
    rows[[length(rows) + 1L]] <- data.frame(
      test = "Wald on gamma (one-sided)", statistic = W, null = "N(0,1)",
      p.value = if (is.na(W)) NA_real_ else stats::pnorm(W, lower.tail = FALSE),
      stringsAsFactors = FALSE
    )
  }
  if ("m3t" %in% test) {
    m2 <- mean((e - mean(e))^2)
    m3 <- mean((e - mean(e))^3)
    M3T <- m3 / sqrt(6 * m2^3 / n)
    ## Production: inefficiency makes the composed error NEGATIVELY skewed, so
    ## the alternative is one-sided in that direction.
    rows[[length(rows) + 1L]] <- data.frame(
      test = "M3T (third moment)", statistic = M3T, null = "N(0,1)",
      p.value = stats::pnorm(M3T, lower.tail = TRUE), stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, rows)
  out$reject <- !is.na(out$p.value) & out$p.value < level
  attr(out, "logLik_H1") <- ll1
  attr(out, "logLik_H0") <- ll0
  attr(out, "level") <- level
  rownames(out) <- NULL
  out
}
