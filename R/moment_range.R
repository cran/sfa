## moment_range(): what skewness and excess kurtosis can this pair PRODUCE?
## Gap L19, from Papadopoulos and Parmeter (2021), EJOR 293:990-1001, sections
## 2.1-2.3. See notes/code_history/moment_range.md.
##
## spec_test() (gap L2) is the formal test from the same authors. This is the
## cruder thing that comes before it and that PP2021 derive first: the composed
## error of a given (v, u) pair can only reach certain skewness and excess
## kurtosis values at all, so residual moments outside that set refute the pair
## without any test statistic. PP2021 call that a Type II failure, to sit
## alongside the Type I failure (wrong skew) of Waldman (1982).
##
## Their own warning, section 2.3.4, is the reason this is a diagnostic and not
## a test: the check "may carry a high Type I error" in small samples. The
## print method repeats it and points at spec_test().

## Excess kurtosis of the symmetric component; its skewness is 0 by symmetry.
## Student t is included with a degrees-of-freedom argument because the package
## fits it (model_name = "tHN"); gamma2 = 6/(df - 4) is finite only for df > 4.
.MR_NOISE <- c(uniform = -6 / 5, normal = 0, logistic = 6 / 5, laplace = 3)

## d^j/dz^j log(2 * Phi(z)), PP2021 equation (7). zeta1 is the Mills ratio,
## formed in logs so that z well into the left tail does not become 0/0.
.mr_zeta <- function(z) {
  z1 <- exp(stats::dnorm(z, log = TRUE) - stats::pnorm(z, log.p = TRUE))
  z2 <- -z * z1 - z1^2
  z3 <- -z1 - z * z2 - 2 * z1 * z2
  z4 <- -2 * z2 - z * z3 - 2 * z2^2 - 2 * z1 * z3
  list(z1 = z1, z2 = z2, z3 = z3, z4 = z4)
}

## Skewness and excess kurtosis of a truncated normal u ~ TN(mu, sigma_u)
## truncated below at 0, as a function of z = mu/sigma_u. z = 0 is the
## half-normal and returns (0.9953, 0.8692); z -> -Inf tends to (2, 6).
##
## Accurate to ~1e-5 for |z| <= 10, checked against direct integration. Further
## into the left tail zeta4 is a difference of terms of order z*zeta3 against a
## result of order 6/z^4, and precision goes: at z = -30 this gives 5.925,
## quadrature gives 5.948 and PP2021's own table says 5.909. Nothing here
## depends on that region -- the supremum 6 is used analytically, and
## .MR_TN_G2MIN searches only z in [0, 6].
.mr_tn_moments <- function(z) {
  ze <- .mr_zeta(z)
  d <- 1 + ze$z2
  list(g1 = ze$z3 / d^(3 / 2), g2 = ze$z4 / d^2)
}

## The attainable set of (gamma1(u), gamma2(u)) for each inefficiency family.
## `g1max`/`g2max` are suprema, not maxima: they are approached as the shape
## parameter runs to its limit and are not attained inside the family.
##
## Half-normal uses the closed forms sqrt(2)(4-pi)/(pi-2)^{3/2} and
## 8(pi-3)/(pi-2)^2 rather than the rounded 0.9953/0.8692 that spec_test()
## carries. Generalized exponential has no clean closed form and is computed
## once from its own quantile function.
##
## Gamma is the interesting row: gamma1(u) = 2/sqrt(k) and gamma2(u) = 6/k are
## unbounded as k -> 0, so a gamma inefficiency term can NEVER suffer a Type II
## failure. PP2021 note the price -- k < 1 makes the density log-convex, which
## costs the composed error its log-concavity.
.mr_ineff_set <- function(family) {
  switch(family,
    halfnormal = {
      g1 <- sqrt(2) * (4 - pi) / (pi - 2)^(3 / 2)
      g2 <- 8 * (pi - 3) / (pi - 2)^2
      list(g1min = g1, g1max = g1, g2min = g2, g2max = g2, fixed = TRUE)
    },
    exponential = list(g1min = 2, g1max = 2, g2min = 6, g2max = 6, fixed = TRUE),
    genexponential = {
      g <- .MR_GENEXP
      list(g1min = g[["g1"]], g1max = g[["g1"]],
           g2min = g[["g2"]], g2max = g[["g2"]], fixed = TRUE)
    },
    truncatednormal = list(
      g1min = 0, g1max = 2,
      ## The one number here that is not a limit: gamma2(u) dips BELOW zero
      ## near z = 2, so the truncated normal is the only u in this list whose
      ## composed error can be platykurtic under normal noise.
      g2min = .MR_TN_G2MIN, g2max = 6, fixed = FALSE
    ),
    gamma = list(g1min = 0, g1max = Inf, g2min = 0, g2max = Inf, fixed = FALSE),
    stop("unknown inefficiency family: ", family, call. = FALSE)
  )
}

## Generalized exponential with shape 2, the one spec_test() draws as
## -log(1 - sqrt(U)). Its cumulant generating function is
## K(t) = log Gamma(3) + log Gamma(1-t) - log Gamma(3-t), so every cumulant is
## a difference of polygammas and nothing needs simulating:
##   kappa_r = (-1)^r [psi^(r-1)(1) - psi^(r-1)(3)].
## These come out at 1.6099 and 4.0800, the rounded pair spec_test() carries.
.MR_GENEXP <- local({
  k2 <- psigamma(1, 1) - psigamma(3, 1)
  k3 <- -psigamma(1, 2) + psigamma(3, 2)
  k4 <- psigamma(1, 3) - psigamma(3, 3)
  c(g1 = k3 / k2^(3 / 2), g2 = k4 / k2^2)
})

.MR_TN_G2MIN <- local({
  f <- function(z) .mr_tn_moments(z)$g2
  o <- stats::optimize(f, interval = c(0, 6))
  o$objective
})

## Range of gamma2(eps) = gamma2(v)(1-R)^2 + gamma2(u) R^2 over R in (0,1) and
## over the family's gamma2(u), PP2021 equation (10). Linear in gamma2(u), so
## its extremes sit at the endpoints; quadratic in R, so the candidates are the
## two endpoints plus the stationary point where that falls inside (0,1).
.mr_g2_range <- function(g2v, g2u_lo, g2u_hi) {
  cand <- numeric(0)
  for (g2u in c(g2u_lo, g2u_hi)) {
    if (!is.finite(g2u)) {
      cand <- c(cand, Inf)
      next
    }
    q <- function(R) g2v * (1 - R)^2 + g2u * R^2
    Rs <- c(0, 1)
    den <- g2v + g2u
    if (is.finite(den) && abs(den) > 1e-12) {
      Rstar <- g2v / den
      if (Rstar > 0 && Rstar < 1) Rs <- c(Rs, Rstar)
    }
    cand <- c(cand, q(Rs))
  }
  c(min(cand), max(cand))
}

## Sample skewness and excess kurtosis of whatever residuals we were given.
.mr_residuals <- function(object, nm) {
  if (inherits(object, "sfareg")) {
    e <- object$epsilon_hat
    if (is.null(e)) e <- object$ols_residuals
    if (is.null(e)) {
      stop("moment_range(): this fit does not retain its OLS residuals. Pass ",
        "the residuals directly, or use a fit from sfm().",
        call. = FALSE
      )
    }
    list(e = as.numeric(e), name = deparse(object$call$formula),
         production = .sfa_inefdec(object))
  } else if (is.numeric(object)) {
    list(e = as.numeric(object), name = nm, production = NA)
  } else {
    stop("`object` must be an \"sfareg\" fit or a numeric vector of OLS ",
      "residuals.",
      call. = FALSE
    )
  }
}

moment_range <- function(object = NULL,
                         noise = c("normal", "logistic", "laplace", "uniform", "t"),
                         inefficiency = c("halfnormal", "exponential",
                                          "genexponential", "truncatednormal",
                                          "gamma"),
                         inefdec = NULL, df = 8) {
  noise <- match.arg(noise, several.ok = TRUE)
  inefficiency <- match.arg(inefficiency, several.ok = TRUE)
  if (length(df) != 1L || !is.numeric(df) || !is.finite(df) || df <= 4) {
    stop("`df` must be a single number greater than 4: the Student t has no ",
      "finite kurtosis at or below 4 degrees of freedom.",
      call. = FALSE
    )
  }

  ## Sign convention. Production (eps = v - u) puts the skewness of the
  ## composed error on the NEGATIVE side; cost (eps = v + u) on the positive.
  g1 <- g2 <- NA_real_
  n <- 0L
  nm <- deparse(substitute(object))
  production <- TRUE
  if (!is.null(object)) {
    r <- .mr_residuals(object, nm)
    nm <- r$name
    if (!is.na(r$production)) production <- r$production
    e <- r$e[is.finite(r$e)]
    n <- length(e)
    if (n < 10L) {
      stop("moment_range(): ", n, " usable residuals is too few.", call. = FALSE)
    }
    ec <- e - mean(e)
    m2 <- mean(ec^2); m3 <- mean(ec^3); m4 <- mean(ec^4)
    g1 <- m3 / m2^(3 / 2)
    g2 <- m4 / m2^2 - 3
  }
  if (!is.null(inefdec)) {
    production <- if (is.character(inefdec)) {
      !grepl("cost", inefdec, ignore.case = TRUE)
    } else {
      isTRUE(inefdec)
    }
  }

  rows <- list()
  for (v in noise) {
    g2v <- if (identical(v, "t")) 6 / (df - 4) else .MR_NOISE[[v]]
    for (u in inefficiency) {
      su <- .mr_ineff_set(u)
      ## |gamma1(eps)| = gamma1(u) R^{3/2} with R in (0,1), so the composed
      ## skewness spans zero up to the family's own supremum, and the noise
      ## distribution does not enter at all. PP2021 equation (5).
      g1lo <- if (production) -su$g1max else 0
      g1hi <- if (production) 0 else su$g1max
      g2r <- .mr_g2_range(g2v, su$g2min, su$g2max)

      skew_ok <- kurt_ok <- NA
      type1 <- NA
      if (n > 0L) {
        skew_ok <- g1 >= g1lo && g1 <= g1hi
        kurt_ok <- g2 >= g2r[1] && g2 <= g2r[2]
        type1 <- if (production) g1 > 0 else g1 < 0
      }
      rows[[length(rows) + 1L]] <- data.frame(
        noise = v, inefficiency = u,
        g1_lo = g1lo, g1_hi = g1hi,
        g2_lo = g2r[1], g2_hi = g2r[2],
        skew_ok = skew_ok, kurt_ok = kurt_ok,
        verdict = if (n == 0L) NA_character_ else {
          bad <- character(0)
          if (isTRUE(type1)) {
            "type I (wrong skew)"
          } else {
            if (!isTRUE(skew_ok)) bad <- c(bad, "skew")
            if (!isTRUE(kurt_ok)) bad <- c(bad, "kurt")
            if (!length(bad)) "ok" else paste0("type II (", paste(bad, collapse = "+"), ")")
          }
        },
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  attr(out, "skewness") <- g1
  attr(out, "excess_kurtosis") <- g2
  attr(out, "nobs") <- n
  attr(out, "data.name") <- nm
  attr(out, "production") <- production
  attr(out, "df") <- df
  class(out) <- c("sfa_moment_range", "data.frame")
  out
}

print.sfa_moment_range <- function(x, ...) {
  g1 <- attr(x, "skewness"); g2 <- attr(x, "excess_kurtosis")
  n <- attr(x, "nobs")
  cat("\nAttainable skewness and excess kurtosis of the composed error\n",
    "Papadopoulos and Parmeter (2021), sections 2.1-2.2\n\n",
    sep = ""
  )
  cat("orientation: ", if (isTRUE(attr(x, "production"))) {
    "production (eps = v - u, skewness <= 0)"
  } else {
    "cost (eps = v + u, skewness >= 0)"
  }, "\n", sep = "")
  if (n > 0L) {
    cat("data:  ", attr(x, "data.name"), "  (n = ", n, ")\n", sep = "")
    cat(sprintf("  sample skewness        : %+.4f\n", g1))
    cat(sprintf("  sample excess kurtosis : %+.4f\n", g2))
  }
  cat("\n")
  d <- as.data.frame(x)
  d$g1_lo <- sprintf("%+.3f", d$g1_lo)
  d$g1_hi <- sprintf("%+.3f", d$g1_hi)
  d$g2_lo <- sprintf("%+.3f", d$g2_lo)
  d$g2_hi <- sprintf("%+.3f", d$g2_hi)
  print(d, row.names = FALSE)
  if (n > 0L) {
    nok <- sum(x$verdict == "ok", na.rm = TRUE)
    cat("\n  ", nok, " of ", nrow(x), " pairs can produce these two moments.\n",
      sep = ""
    )
    cat("  This is a RANGE CHECK, not a test: PP2021 section 2.3 shows it\n",
      "  over-rejects in small samples (up to 26% at n = 200). Use\n",
      "  spec_test() for the formal statistic.\n",
      sep = ""
    )
  }
  cat("\n")
  invisible(x)
}
