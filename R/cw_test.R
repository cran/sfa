## Chen and Wang (2012), Econometric Reviews 31(6):625-653.
## Centered-residuals moment estimator and specification test.
## See notes/code_history/cw_test.md.
##
## The selling point over the goodness-of-fit tests in gof_test() is that this
## one never touches the composed-error PDF or CDF. It works from the
## ARITHMETIC MOMENTS and the CHARACTERISTIC FUNCTION of the two components
## separately, which are in closed form for far more distributional pairs than
## their convolution is -- the normal-gamma model being the standard example of
## a pair with no closed-form composed density.
##
## Two things make it work:
##
##   * CENTERED residuals. The intercept of a frontier is not identified
##     separately from E[u], so the LS intercept estimates alpha - E[u], not
##     alpha. Centering the residuals cancels the inconsistent intercept
##     entirely, which is why the test is valid off an ordinary LS fit and does
##     not need maximum likelihood.
##   * A variance that corrects for BOTH nuisances: estimating theta, and
##     centering. The centering correction is the -E[d psi / d eps] * eps_c
##     term in xi below; leaving it out understates the variance and oversizes
##     the test.

## Theoretical quantities for the centered composed error, all from the same
## inefficiency quadrature the composed CDF uses.
##
##   E[cos(tau eps_c)] =  E[cos(tau v)] * hc(tau)
##   E[sin(tau eps_c)] =  E[cos(tau v)] * hs(tau)
##   hc(tau) = E[cos(tau u)] cos(tau mu1) + E[sin(tau u)] sin(tau mu1)
##   hs(tau) = E[cos(tau u)] sin(tau mu1) - E[sin(tau u)] cos(tau mu1)
##
## which are Chen and Wang's (39) and (35). E[sin(tau v)] = 0 is what the
## normality (strictly, symmetry) of the noise buys.
.cw_theory <- function(model_name, par, taus, kmax = 3L, n_nodes = 256L) {
  sp <- .composed_u_spec(model_name, par)
  if (!identical(sp$noise$type, "normal") || !is.null(sp$mix_df)) {
    stop("cw_test(): the sine/cosine restriction uses E[sin(tau v)] = 0, ",
      "which needs symmetric noise with a closed-form characteristic ",
      "function. model_name = ", dQuote(model_name), " has Student-t noise.",
      call. = FALSE
    )
  }
  sv <- sp$noise$sigma_v
  nd <- .composed_u_nodes(sp, n_nodes)
  u <- nd$u
  w <- nd$w

  mu <- vapply(0:kmax, function(k) sum(w * u^k), numeric(1)) ## E[u^k], k = 0..kmax
  mu1 <- mu[2]
  ## Central moments of u, then of eps_c = v - (u - mu1).
  uc <- vapply(0:kmax, function(s) {
    sum(vapply(0:s, function(l) (-1)^l * choose(s, l) * mu[s - l + 1] * mu1^l, numeric(1)))
  }, numeric(1))
  ## E[v^m] for v ~ N(0, sv^2): 0 for odd m, sv^m (m-1)!! for even m.
  nv <- function(m) {
    if (m == 0L) {
      return(1)
    }
    if (m %% 2L == 1L) {
      return(0)
    }
    sv^m * prod(seq.int(1L, m - 1L, by = 2L))
  }
  eps_m <- vapply(2:kmax, function(k) {
    sum(vapply(0:k, function(s) (-1)^s * choose(k, s) * nv(k - s) * uc[s + 1], numeric(1)))
  }, numeric(1))

  Ecv <- exp(-sv^2 * taus^2 / 2) ## E[cos(tau v)]
  Ecu <- vapply(taus, function(t) sum(w * cos(t * u)), numeric(1))
  Esu <- vapply(taus, function(t) sum(w * sin(t * u)), numeric(1))
  hc <- Ecu * cos(taus * mu1) + Esu * sin(taus * mu1)
  hs <- Ecu * sin(taus * mu1) - Esu * cos(taus * mu1)

  list(
    mu1 = mu1, moments = eps_m, Ecos = Ecv * hc, Esin = Ecv * hs,
    sigma_v = sv
  )
}

## The two shape parameters, in the order (sigma_v, second), for the models
## this is implemented for. Everything else about the model comes from the
## registry, so this is only the naming.
.cw_par_names <- function(model_name) {
  switch(model_name,
    "NHN" = c("lambda", "sigma"),
    "NE" = , "NR" = , "NGE" = c("sigv", "sigu"),
    "NU" = c("sigv", "theta"),
    stop("cw_test(): implemented for the two-parameter normal-noise models ",
      "NHN, NE, NR, NU and NGE. Got model_name = ", dQuote(model_name), ".",
      call. = FALSE
    )
  )
}

## Build a named par vector of the model's own shape parameters from the
## working vector (sigma_v, second_scale), so .composed_u_spec() can read it.
.cw_build_par <- function(model_name, th) {
  sv <- th[1]
  s2 <- th[2]
  switch(model_name,
    "NHN" = c(lambda = s2 / sv, sigma = sqrt(sv^2 + s2^2)),
    "NU" = c(sigv = sv, theta = s2),
    c(sigv = sv, sigu = s2)
  )
}

## Solve (1/n) sum psi_1 = 0: match the sample central moments of order
## 2..(p+1) to their theoretical values. For NHN and NE this has the closed
## form Chen and Wang give in (32) and it IS the COLS inversion, so that is
## used as the starting value and the solve only has to polish it.
.cw_estimate <- function(m_hat, model_name, start, n_nodes = 256L) {
  obj <- function(l) {
    th <- exp(l)
    tt <- tryCatch(.cw_theory(model_name, .cw_build_par(model_name, th),
      taus = 1, kmax = length(m_hat) + 1L, n_nodes = n_nodes
    ), error = function(e) NULL)
    if (is.null(tt)) {
      return(1e10)
    }
    d <- (tt$moments - m_hat) / pmax(abs(m_hat), 1e-8)
    sum(d^2)
  }
  op <- stats::optim(log(pmax(start, 1e-6)), obj,
    method = "Nelder-Mead",
    control = list(reltol = 1e-12, maxit = 2000)
  )
  list(theta = exp(op$par), value = op$value, convergence = op$convergence)
}

#' Chen and Wang (2012) centered-residuals moment estimator and test
#'
#' @param x An `"sfareg"` fit, or a numeric vector of regression residuals.
#' @param model_name The distributional pair being tested. Required when `x`
#'   is a numeric vector; taken from the fit otherwise.
#' @param tau Frequencies at which the characteristic function is compared.
#'   Chen and Wang suggest values around 1, and a SINGLE frequency: combining
#'   several oversizes the test in finite samples (see Details).
#' @param type `"cosine"` (their recommendation) or `"sine"`.
#' @param data Data the model was fitted to, when `x` is a fit.
#' @param n_nodes Quadrature nodes for the theoretical moments.
#' @return An object of class `"cw_test"` with the statistic, its p-value and
#'   the moment estimates.
#' @details
#' The test never evaluates the composed-error density or distribution
#' function. It compares the empirical characteristic function of the centred
#' residuals with its theoretical value, which follows from the characteristic
#' functions of the two components separately -- available in closed form for
#' many pairs whose convolution is not.
#'
#' Centring the residuals is what makes the test valid off an ordinary least
#' squares fit: a frontier's intercept is not identified separately from
#' `E[u]`, so the LS intercept estimates `alpha - E[u]`, and centring cancels
#' it. The variance corrects for two nuisances, estimating the scale parameters
#' and centring; omitting the second understates it and oversizes the test.
#'
#' \bold{Use one frequency, not several.} Measured here over 500 replications
#' at a nominal 5\%, half-normal data:
#'
#' \tabular{lrrrr}{
#' test \tab n \tab tau=1 \tab tau=1.5 \tab tau=2 \cr
#' cosine \tab 500 \tab 0.070 \tab 0.054 \tab 0.050 \cr
#' cosine \tab 2000 \tab 0.046 \tab 0.044 \tab 0.050 \cr
#' sine \tab 500 \tab 0.092 \tab 0.064 \tab 0.042 \cr
#' sine \tab 2000 \tab 0.068 \tab 0.062 \tab 0.064
#' }
#'
#' Against exponential inefficiency the cosine test's power at tau = 1 is 0.231
#' at n = 500 and 0.952 at n = 2000, against 0.108 and 0.706 for the sine test.
#' Both of Chen and Wang's conclusions reproduce: the cosine test is the less
#' sensitive to the choice of tau and the more powerful. Passing several
#' frequencies at once gives rejection rates of 0.09 to 0.15 at a nominal 5\%,
#' because the cosine moments at nearby frequencies are close to collinear and
#' the covariance matrix is then near-singular; the authors report the same.
#' @export
cw_test <- function(x, model_name = NULL, tau = 1,
                    type = c("cosine", "sine"), data = NULL, n_nodes = 256L) {
  type <- match.arg(type)
  if (inherits(x, "sfareg")) {
    if (is.null(model_name)) model_name <- x$model_name
    e <- .sfa_eps_hat(x, data, parent.frame())
  } else {
    if (!is.numeric(x)) {
      stop("`x` must be an \"sfareg\" fit or a numeric residual vector.", call. = FALSE)
    }
    if (is.null(model_name)) {
      stop("cw_test(): `model_name` is required when `x` is a residual vector.",
        call. = FALSE
      )
    }
    e <- as.numeric(x)
  }
  .cw_par_names(model_name) ## validates the model early
  if (!length(tau) || !all(is.finite(tau)) || any(tau <= 0)) {
    stop("cw_test(): `tau` must be positive and finite.", call. = FALSE)
  }
  ec <- e - mean(e) ## the centred residual: this is the whole point
  n <- length(ec)
  q <- length(tau)

  m_hat <- c(mean(ec^2), mean(ec^3))
  ## The moment inversion needs a negatively skewed residual for a production
  ## frontier: sigma_u solves sigma_u^3 = m3/k3 with k3 < 0. With m3 >= 0 there
  ## is no admissible solution, and continuing would report a test statistic
  ## computed at a meaningless parameter. This is the Type I failure of Olson,
  ## Schmidt and Waldman (1980).
  if (!is.finite(m_hat[2]) || m_hat[2] >= 0) {
    stop("cw_test(): the third central moment of the residuals is ",
      format(m_hat[2], digits = 3), " >= 0, so the residuals are skewed the ",
      "WRONG WAY for a production frontier and the moment equations have no ",
      "admissible solution. See skewness_test() for whether the sign is wrong ",
      "by more than sampling noise.",
      call. = FALSE
    )
  }
  ## Start from the closed-form moment inversion, which for NHN and NE is
  ## exactly Chen and Wang's (32).
  k3 <- switch(model_name,
    "NHN" = (1 - 4 / pi) * sqrt(2 / pi),
    "NE" = -2,
    NA_real_
  )
  s_start <- if (is.finite(k3) && m_hat[2] < 0) {
    su <- (m_hat[2] / k3)^(1 / 3)
    c(sqrt(max(m_hat[1] - su^2 * (1 - 2 / pi) * (model_name == "NHN") -
      su^2 * (model_name == "NE"), 1e-4)), su)
  } else {
    c(stats::sd(ec) / sqrt(2), stats::sd(ec) / sqrt(2))
  }
  fitm <- .cw_estimate(m_hat, model_name, s_start, n_nodes)
  th <- fitm$theta
  par <- .cw_build_par(model_name, th)

  th_fun <- function(t) .cw_theory(model_name, .cw_build_par(model_name, t),
    taus = tau, kmax = 3L, n_nodes = n_nodes
  )
  TT <- th_fun(th)

  ## psi_1 and psi_2 at the estimate.
  psi1 <- cbind(ec^2 - TT$moments[1], ec^3 - TT$moments[2])
  Mth <- if (type == "cosine") TT$Ecos else TT$Esin
  psi2 <- outer(ec, tau, function(a, b) if (type == "cosine") cos(a * b) else sin(a * b))
  psi2 <- sweep(psi2, 2, Mth, "-")
  Mbar <- colMeans(psi2)

  ## The centering correction: xi = psi - E[d psi / d eps_c] * eps_c.
  ## For psi_1 that is (0, 3*mu2); for the cosine it is -tau*E[sin(tau eps_c)],
  ## for the sine +tau*E[cos(tau eps_c)]. Dropping it oversizes the test.
  xi1 <- psi1
  xi1[, 2] <- xi1[, 2] - 3 * TT$moments[1] * ec
  dpsi2 <- if (type == "cosine") -tau * TT$Esin else tau * TT$Ecos
  xi2 <- psi2 - outer(ec, dpsi2, "*")

  ## K_j = E[d psi_j / d theta']; the sign is carried by differentiating the
  ## theoretical term, which enters psi with a minus.
  K1 <- -numDeriv::jacobian(function(t) th_fun(t)$moments, th)
  K2 <- -numDeriv::jacobian(function(t) {
    z <- th_fun(t)
    if (type == "cosine") z$Ecos else z$Esin
  }, th)
  if (q == 1L) K2 <- matrix(K2, nrow = 1L)

  A <- tryCatch(K2 %*% solve(K1), error = function(err) NULL)
  if (is.null(A)) {
    stop("cw_test(): the moment Jacobian K1 is singular at the estimate.",
      call. = FALSE
    )
  }
  Z <- xi2 - xi1 %*% t(A)
  Om <- crossprod(Z) / n
  Oi <- tryCatch(solve(Om), error = function(err) NULL)
  if (is.null(Oi)) {
    stop("cw_test(): the covariance matrix is singular; try fewer or better ",
      "separated `tau` values.",
      call. = FALSE
    )
  }
  D <- as.numeric(n * t(Mbar) %*% Oi %*% Mbar)

  sigmas <- c(sigma_v = unname(th[1]), sigma_u = unname(th[2]))
  out <- list(
    statistic = D, df = q, p.value = stats::pchisq(D, df = q, lower.tail = FALSE),
    type = type, tau = tau, model_name = model_name,
    estimate = sigmas, moment_fit = fitm, nobs = n, Mbar = Mbar
  )
  class(out) <- "cw_test"
  out
}

#' @export
print.cw_test <- function(x, ...) {
  cat("Chen and Wang (2012) ", x$type, " test\n", sep = "")
  cat("  H0: composed error is ", x$model_name, "\n", sep = "")
  cat("  tau = ", paste(format(x$tau), collapse = ", "),
    "   n = ", x$nobs, "\n",
    sep = ""
  )
  cat("  moment estimates: sigma_v = ", format(x$estimate[["sigma_v"]], digits = 4),
    ", sigma_u = ", format(x$estimate[["sigma_u"]], digits = 4), "\n",
    sep = ""
  )
  cat("  D = ", format(x$statistic, digits = 5), " on ", x$df,
    " df,  p = ", format.pval(x$p.value, digits = 4), "\n",
    sep = ""
  )
  invisible(x)
}
