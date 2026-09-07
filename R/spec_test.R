## spec_test(): is this PAIR of distributions defensible? Gap L2, from
## Papadopoulos and Parmeter (2023), Economics Letters 233:111390, resting on
## Papadopoulos and Parmeter (2021), EJOR 293:990-1001.
## See notes/code_history/spec_test.md and notes/L2_spec_test_design.md.
##
## `skewness_test()` already asks whether the residuals are skewed the way
## inefficiency makes them. That is a test of whether there IS inefficiency,
## not of whether the assumed distributions are the right ones -- and the
## package offers fifteen cross-sectional distributions.
##
## This test asks the second question, and asks it BEFORE anything is fitted:
## it needs only the OLS residuals, so it costs nothing and cannot be
## contaminated by a frontier optimizer that failed.

## Skewness and excess kurtosis of the two components. gamma1(v) is 0 for any
## symmetric noise, so only its excess kurtosis matters.
.SPEC_NOISE <- c(uniform = -6 / 5, normal = 0, logistic = 6 / 5, laplace = 3)
.SPEC_INEFF <- list(
  halfnormal = c(g1 = 0.9953, g2 = 0.8692),
  exponential = c(g1 = 2, g2 = 6),
  genexponential = c(g1 = 1.61, g2 = 4.08)
)

## The restriction being tested.
##
## With R = SNR^2/(1 + SNR^2), the share of composed-error variance coming from
## u, PP2021's equations (5) and (10) become
##
##   |gamma1(eps)| = gamma1(u) R^{3/2}                    -- identifies R
##   gamma2(eps)   = gamma2(v)(1-R)^2 + gamma2(u) R^2     -- the restriction
##
## so the skewness pins R down and the kurtosis then has nothing left to fit.
## Gamma is the gap between the two, and is 0 when the assumed pair is right.
##
## Note |gamma1|^{2/3} = |m3|^{2/3}/m2, which is what keeps the derivatives
## below short.
.spec_gamma <- function(m2, m3, m4, g1u, g2u, g2v) {
  R <- (abs(m3)^(2 / 3) / m2) * g1u^(-2 / 3)
  g2 <- m4 / m2^2 - 3
  list(R = R, g2 = g2, Gamma = g2 - g2v * (1 - R)^2 - g2u * R^2)
}

## Analytic Jacobian of Gamma with respect to (m2, m3, m4). Checked against a
## numerical one in the tests, because a hand-differentiated expression that is
## never compared to anything is a guess.
.spec_jacobian <- function(m2, m3, m4, g1u, g2u, g2v) {
  R <- (abs(m3)^(2 / 3) / m2) * g1u^(-2 / 3)
  dR_dm2 <- -R / m2
  dR_dm3 <- (2 / 3) * R / m3
  dG_dR <- 2 * g2v * (1 - R) - 2 * g2u * R
  c(
    m2 = -2 * m4 / m2^3 + dG_dR * dR_dm2,
    m3 = dG_dR * dR_dm3,
    m4 = 1 / m2^2
  )
}

## Asymptotic covariance of the central moments (m2, m3, m4). These are the
## standard expressions, and are the ones PP2023 reprints. PP2021's own
## equation (13) is deliberately NOT used: PP2023 records that its H_k terms
## are missing squares, and there is no reason to transcribe an equation whose
## authors have published a correction to it.
.spec_moment_vcov <- function(m) {
  m2 <- m[2]; m3 <- m[3]; m4 <- m[4]; m5 <- m[5]
  m6 <- m[6]; m7 <- m[7]; m8 <- m[8]
  S <- matrix(0, 3, 3, dimnames = list(c("m2", "m3", "m4"), c("m2", "m3", "m4")))
  S[1, 1] <- m4 - m2^2
  S[2, 2] <- m6 - m3^2 - 6 * m2 * m4 + 9 * m2^3
  S[3, 3] <- m8 - m4^2 - 8 * m3 * m5 + 16 * m2 * m3^2
  S[1, 2] <- S[2, 1] <- m5 - 4 * m2 * m3
  S[1, 3] <- S[3, 1] <- m6 - m2 * m4 - 4 * m3^2
  S[2, 3] <- S[3, 2] <- m7 - 3 * m2 * m5 - 5 * m3 * m4 + 12 * m2^2 * m3
  S
}


## Draws from each component with a GIVEN standard deviation. The scale
## parameter is recovered numerically from a unit-scale sample rather than
## derived per distribution: three closed forms hand-derived is three chances
## to be wrong, and this is exact to the precision of the reference draw.
.SPEC_SD1 <- local({
  set.seed(1L)
  N <- 2e6
  c(
    halfnormal = stats::sd(abs(stats::rnorm(N))),
    exponential = stats::sd(stats::rexp(N)),
    genexponential = stats::sd(-log(1 - sqrt(stats::runif(N))))
  )
})

.spec_rnoise <- function(n, family, sd) {
  switch(family,
    normal = stats::rnorm(n, 0, sd),
    logistic = stats::rlogis(n, 0, sd * sqrt(3) / pi),
    laplace = {
      b <- sd / sqrt(2)
      ifelse(stats::runif(n) < 0.5, 1, -1) * stats::rexp(n, 1 / b)
    },
    uniform = stats::runif(n, -sd * sqrt(3), sd * sqrt(3))
  )
}

.spec_rineff <- function(n, family, sd) {
  sc <- sd / .SPEC_SD1[[family]]
  switch(family,
    halfnormal = abs(stats::rnorm(n, 0, sc)),
    exponential = stats::rexp(n, 1 / sc),
    genexponential = -sc * log(1 - sqrt(stats::runif(n)))
  )
}

spec_test <- function(object,
                      noise = c("normal", "logistic", "laplace", "uniform"),
                      inefficiency = c("halfnormal", "exponential",
                                       "genexponential"),
                      null = c("bootstrap", "asymptotic"),
                      B = 999, level = 0.05, seed = NULL) {
  noise <- match.arg(noise)
  inefficiency <- match.arg(inefficiency)
  null <- match.arg(null)
  if (length(B) != 1L || !is.numeric(B) || !is.finite(B) || B < 19) {
    stop("`B` must be a single number >= 19.", call. = FALSE)
  }
  B <- as.integer(B)

  ## OLS residuals. From a fit, the ones recorded at fit time -- NOT the
  ## composed residuals of the fitted model, which already embody the
  ## distributional assumption under test and would make the test circular.
  if (inherits(object, "sfareg")) {
    e <- object$epsilon_hat
    if (is.null(e)) e <- object$ols_residuals
    if (is.null(e)) {
      stop("spec_test(): this fit does not retain its OLS residuals. Pass the ",
        "residuals directly, or use a fit from sfm().",
        call. = FALSE
      )
    }
    nm <- deparse(object$call$formula)
  } else if (is.numeric(object)) {
    e <- object
    nm <- deparse(substitute(object))
  } else {
    stop("`object` must be an \"sfareg\" fit or a numeric vector of OLS ",
      "residuals.",
      call. = FALSE
    )
  }
  e <- as.numeric(e)
  e <- e[is.finite(e)]
  n <- length(e)
  if (n < 50) {
    stop("spec_test(): ", n, " usable residuals is too few. The statistic ",
      "needs moments up to the 8th for its variance.",
      call. = FALSE
    )
  }

  ec <- e - mean(e)
  m <- vapply(1:8, function(k) mean(ec^k), numeric(1))
  m2 <- m[2]; m3 <- m[3]; m4 <- m[4]
  g1u <- .SPEC_INEFF[[inefficiency]][["g1"]]
  g2u <- .SPEC_INEFF[[inefficiency]][["g2"]]
  g2v <- .SPEC_NOISE[[noise]]

  gg <- .spec_gamma(m2, m3, m4, g1u, g2u, g2v)
  g1 <- m3 / m2^(3 / 2)

  ## R > 1 is Type II failure: the residuals are more skewed than this
  ## inefficiency distribution can be, so the pair is refuted outright and a
  ## p-value would be beside the point.
  type2 <- gg$R > 1

  J <- .spec_jacobian(m2, m3, m4, g1u, g2u, g2v)
  S <- .spec_moment_vcov(m)
  vG <- as.numeric(t(J) %*% S %*% J) / n
  stat <- if (is.finite(vG) && vG > 0) gg$Gamma / sqrt(vG) else NA_real_

  boot <- NULL
  if (identical(null, "asymptotic") || isTRUE(type2)) {
    ## The delta-method p-value. Measurably unreliable when the inefficiency
    ## term is heavy-tailed: the variance needs moments up to the 8th, and for
    ## an exponential u at n = 500-2000 that estimate is far too noisy. See the
    ## size table in ?spec_test. Under Type II failure the pair is refuted by
    ## the skewness alone and there is nothing to bootstrap under.
    p <- if (is.finite(stat)) 2 * stats::pnorm(-abs(stat)) else NA_real_
  } else {
    ## Parametric bootstrap under the ASSUMED pair. The skewness identifies the
    ## variance share R and m2 the total variance, so the two component scales
    ## are pinned down and the null can be simulated directly. Cheap, because
    ## nothing is fitted -- each replicate is n draws and eight moments.
    if (!is.null(seed)) {
      .st <- .rng_snapshot()
      on.exit(.rng_restore(.st), add = TRUE)
      set.seed(seed)
    }
    Rc <- min(max(gg$R, 1e-6), 1 - 1e-6)
    s_u <- sqrt(Rc * m2)
    s_v <- sqrt((1 - Rc) * m2)
    bs <- rep(NA_real_, B)
    for (b in seq_len(B)) {
      eb <- .spec_rnoise(n, noise, s_v) - .spec_rineff(n, inefficiency, s_u)
      eb <- eb - mean(eb)
      mb <- vapply(1:8, function(k) mean(eb^k), numeric(1))
      gb <- .spec_gamma(mb[2], mb[3], mb[4], g1u, g2u, g2v)
      Jb <- .spec_jacobian(mb[2], mb[3], mb[4], g1u, g2u, g2v)
      vb <- as.numeric(t(Jb) %*% .spec_moment_vcov(mb) %*% Jb) / n
      bs[b] <- if (is.finite(vb) && vb > 0) gb$Gamma / sqrt(vb) else NA_real_
    }
    boot <- bs[is.finite(bs)]
    if (length(boot) < 19L) {
      stop("spec_test(): only ", length(boot), " of ", B, " bootstrap ",
        "replicates were usable.",
        call. = FALSE
      )
    }
    ## Two-sided, on the studentised statistic so the bootstrap corrects the
    ## variance estimate rather than merely relocating it.
    p <- (1 + sum(abs(boot) >= abs(stat))) / (length(boot) + 1)
  }

  out <- list(
    statistic = c("z" = stat),
    p.value = p,
    method = paste0(
      "Papadopoulos-Parmeter specification test for the composed error\n",
      "pair: ", noise, " noise / ", inefficiency, " inefficiency"
    ),
    data.name = nm,
    alternative = "the assumed pair of distributions is not consistent with the residuals",
    Gamma = gg$Gamma, se = sqrt(vG),
    R = gg$R, skewness = g1, excess_kurtosis = gg$g2,
    noise = noise, inefficiency = inefficiency,
    g1u = g1u, g2u = g2u, g2v = g2v,
    nobs = n, level = level,
    null = if (isTRUE(type2)) "asymptotic" else null,
    B = if (identical(null, "bootstrap") && !isTRUE(type2)) B else NA_integer_,
    boot = boot,
    type2_failure = type2,
    reject = isTRUE(p < level)
  )
  class(out) <- c("sfa_spec_test", "htest")
  out
}

## Run every pair at once. The paper's own advice: "the researcher can run all
## 12 and use the obtained p-values as comparative evidence". They are not
## sequential tests, so no multiplicity correction is implied by the procedure
## itself -- but reading twelve p-values and picking the largest is a selection,
## and the print method says so.
spec_test_all <- function(object, level = 0.05) {
  rows <- list()
  for (v in names(.SPEC_NOISE)) {
    for (u in names(.SPEC_INEFF)) {
      tt <- tryCatch(spec_test(object, noise = v, inefficiency = u,
        level = level
      ), error = function(e) NULL)
      if (is.null(tt)) next
      rows[[length(rows) + 1]] <- data.frame(
        noise = v, inefficiency = u,
        Gamma = tt$Gamma, z = unname(tt$statistic), p.value = tt$p.value,
        R = tt$R, type2 = tt$type2_failure,
        stringsAsFactors = FALSE
      )
    }
  }
  r <- do.call(rbind, rows)
  r <- r[order(-r$p.value), ]
  rownames(r) <- NULL
  attr(r, "skewness") <- if (nrow(r)) NULL else NULL
  class(r) <- c("sfa_spec_grid", "data.frame")
  r
}

print.sfa_spec_test <- function(x, ...) {
  cat("\n", x$method, "\n\n", sep = "")
  cat("data:  ", x$data.name, "\n", sep = "")
  cat(sprintf("z = %.4f, p-value = %.4g\n", x$statistic, x$p.value))
  cat("alternative hypothesis: ", x$alternative, "\n\n", sep = "")
  cat(sprintf("  sample skewness            : %+.4f  (this u allows |g1| <= %.4f)\n",
    x$skewness, x$g1u
  ))
  cat(sprintf("  sample excess kurtosis     : %+.4f\n", x$excess_kurtosis))
  cat(sprintf("  implied variance share R   : %.4f\n", x$R))
  cat(sprintf("  Gamma (0 if pair is right) : %+.5f  (se %.5f)\n", x$Gamma, x$se))
  cat(sprintf("  observations               : %d\n\n", x$nobs))
  if (isTRUE(x$type2_failure)) {
    cat("  TYPE II FAILURE: the residuals are MORE skewed than a ",
      x$inefficiency, "\n  inefficiency term can be, so this pair is refuted ",
      "outright.\n  R > 1 means the skewness alone rules it out; the p-value ",
      "is\n  beside the point.\n\n",
      sep = ""
    )
  } else if (x$reject) {
    cat("  Reject this pair at the ", 100 * x$level, "% level.\n", sep = "")
  } else {
    cat("  This pair is not rejected at the ", 100 * x$level, "% level.\n",
      sep = ""
    )
  }
  cat("\n")
  invisible(x)
}

print.sfa_spec_grid <- function(x, ...) {
  cat("\nPapadopoulos-Parmeter specification test, all pairs\n")
  cat("Sorted by p-value; the pair at the top is the least contradicted.\n\n")
  d <- as.data.frame(x)
  d$Gamma <- round(d$Gamma, 4)
  d$z <- round(d$z, 3)
  d$p.value <- signif(d$p.value, 3)
  d$R <- round(d$R, 3)
  print(d, row.names = FALSE)
  cat("\nThese are NOT sequential tests, so one does not depend on another.\n")
  cat("But choosing the largest p-value from twelve IS a selection: treat the\n")
  cat("ranking as comparative evidence, not as a significance level.\n")
  if (any(d$type2)) {
    cat("\nRows with type2 = TRUE are refuted by the skewness alone.\n")
  }
  cat("\n")
  invisible(x)
}
