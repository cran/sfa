## Kim and Schmidt (2008), Journal of Econometrics 144:409-427.
## Valid tests of whether technical inefficiency depends on firm characteristics.
## See notes/code_history/uhet_test.md.
##
## The two-step procedure everyone actually runs is: fit the frontier, take the
## JLMS predictor u-hat = E[u | eps], regress it on z, and test the slopes.
## That test is INVALID as usually performed, and not by a little: u-hat is a
## GENERATED dependent variable, so the first-step estimation error in psi
## enters the second-step variance. The naive test omits it.
##
## Kim and Schmidt show the omitted term is E[z grad_psi b] * r_i, so the naive
## test is valid only when G = E[z grad_psi b] = 0 -- which holds when z is
## independent of x, and generally fails when z and x are correlated. That is
## exactly the case applied work is usually in.
##
## There is no distribution-free correction: the form of G and of the score
## depend on the assumed distribution of u. Implemented here for the two cases
## the paper works through, NHN and NE.

## b(psi) = E[u | eps] - E[u], as an explicit function of the parameter vector,
## so its Jacobian in psi can be taken. Demeaning by the MODEL's E[u] rather
## than by the sample mean is what makes b a clean function of psi.
.ks_b_of_psi <- function(psi, Y, X, model_name, inefdec) {
  k <- ncol(X)
  p <- length(psi)
  beta <- psi[(p - k + 1L):p]
  eps <- as.numeric(Y - X %*% beta)
  if (!isTRUE(inefdec)) eps <- -eps

  if (identical(model_name, "NHN")) {
    lam <- psi[1]
    sg <- psi[2]
    s2 <- sg^2
    su <- sg * lam / sqrt(1 + lam^2)
    sv <- sg / sqrt(1 + lam^2)
    mu_star <- -eps * su^2 / s2
    sig_star <- su * sv / sg
    Eu <- su * sqrt(2 / pi)
  } else {
    sv <- psi[1]
    su <- psi[2]
    mu_star <- -eps - sv^2 / su
    sig_star <- rep_len(sv, length(eps))
    Eu <- su
  }
  .jlms_u(mu_star, sig_star) - Eu
}

#' Test whether inefficiency depends on firm characteristics
#'
#' @param object An `"sfareg"` fit from [sfm()], fitted WITHOUT `uhet` -- this
#'   is the first step, estimated under the null.
#' @param z A one-sided formula (e.g. `~ z1 + z2`) or a numeric matrix of the
#'   characteristics whose effect is being tested.
#' @param data The data the model was fitted to.
#' @return A data frame with the corrected and naive statistics.
#' @export
uhet_test <- function(object, z, data = NULL) {
  if (!inherits(object, "sfareg")) {
    stop("`object` must be an \"sfareg\" fit.", call. = FALSE)
  }
  mn <- object$model_name
  if (!mn %in% c("NHN", "NE")) {
    stop("uhet_test(): the correction is not distribution-free -- it depends ",
      "on the assumed distribution of u through both the JLMS predictor and ",
      "the score. Kim and Schmidt work through the half-normal and ",
      "exponential cases, and those are what is implemented. Got model_name = ",
      dQuote(mn), ".",
      call. = FALSE
    )
  }
  ## At the wrong-skew boundary sigma_u -> 0 there is no estimated
  ## inefficiency to explain: b-hat = E[u|eps] - E[u] collapses to zero for
  ## every firm, so "does inefficiency depend on z" has no content. Worse, the
  ## information matrix is SINGULAR there (Waldman 1982), so the r_i in the
  ## correction blow up and the variance with them. Measured over 400
  ## replications this inflated the mean corrected standard error to 70 times
  ## the true sampling spread while its median stayed near 1.1 -- a handful of
  ## boundary fits, not a systematic error.
  if (isTRUE(object$wrong_skew) || isTRUE(object$sigma_u_at_bound)) {
    stop("uhet_test(): this fit sits on the wrong-skew boundary, where ",
      "sigma_u is zero, every E[u | eps] is identical and the information ",
      "matrix is singular (Waldman 1982). There is no estimated inefficiency ",
      "for `z` to explain. See skewness_test().",
      call. = FALSE
    )
  }
  if (is.null(data)) {
    data <- tryCatch(eval(object$call$data, envir = parent.frame()),
      error = function(e) NULL
    )
  }
  if (is.null(data)) {
    stop("uhet_test() needs the data the model was fitted to; pass `data =`.",
      call. = FALSE
    )
  }
  data <- as.data.frame(data)

  ## Step 1 quantities.
  psi <- object$out[, "par"]
  fx <- stats::formula(Formula::Formula(object$formula), lhs = 1, rhs = 1)
  mf <- stats::model.frame(fx, data = data)
  Y <- as.numeric(stats::model.response(mf))
  X <- stats::model.matrix(fx, data = data)
  bn <- intersect(colnames(X), names(psi))
  if (length(bn) != ncol(X)) {
    stop("uhet_test(): could not match the design to the fitted coefficients.",
      call. = FALSE
    )
  }
  X <- X[, bn, drop = FALSE]
  inefdec <- .sfa_inefdec(object)
  N <- length(Y)
  ## `data` has to be the data the model was FITTED to, not merely a data frame
  ## with the right columns: the correction combines per-observation second-step
  ## quantities with the fit's own per-observation scores, and a mismatch would
  ## otherwise surface as "non-conformable arrays" further down.
  n_fit <- tryCatch(as.integer(stats::nobs(object)), error = function(e) NA_integer_)
  if (is.finite(n_fit) && n_fit != N) {
    stop("uhet_test(): `data` has ", N, " usable rows but the fit used ",
      n_fit, ". Pass the data the model was fitted to.",
      call. = FALSE
    )
  }

  ## The characteristics, with an intercept: b is demeaned in the population
  ## but not in the sample, and Kim and Schmidt's own motivation is an F test
  ## on the coefficients OTHER than the intercept.
  Z <- if (inherits(z, "formula")) {
    stats::model.matrix(z, data = data)
  } else {
    zz <- as.matrix(z)
    if (is.null(colnames(zz))) colnames(zz) <- paste0("z", seq_len(ncol(zz)))
    cbind("(Intercept)" = 1, zz)
  }
  if (nrow(Z) != N) {
    stop("uhet_test(): `z` has ", nrow(Z), " rows but the fit has ", N, ".",
      call. = FALSE
    )
  }
  if (!"(Intercept)" %in% colnames(Z)) Z <- cbind("(Intercept)" = 1, Z)
  keep <- setdiff(colnames(Z), "(Intercept)")
  q <- length(keep)
  if (q < 1L) {
    stop("uhet_test(): `z` must contain at least one non-constant variable.",
      call. = FALSE
    )
  }

  bhat <- .ks_b_of_psi(psi, Y, X, mn, inefdec)
  B <- crossprod(Z) / N
  Bi <- tryCatch(solve(B), error = function(e) NULL)
  if (is.null(Bi)) {
    stop("uhet_test(): `z` is collinear.", call. = FALSE)
  }
  gam <- as.numeric(Bi %*% (crossprod(Z, bhat) / N))
  names(gam) <- colnames(Z)

  ## The naive variance treats bhat as if it were observed.
  res <- as.numeric(bhat - Z %*% gam)
  A_naive <- crossprod(Z * res) / N

  ## The correction. G = E[z grad_psi b'], r_i = Iinv s_i with I the OPG
  ## information, and the relevant score is the FIRST-step one -- which is
  ## exactly what estfun() builds.
  S <- tryCatch(estfun.sfareg(object), error = function(e) NULL)
  if (is.null(S)) {
    stop("uhet_test(): the corrected test needs the first-step score, which ",
      "needs the fitted likelihood. Refit with keep_objective = TRUE.",
      call. = FALSE
    )
  }
  S <- S[, names(psi), drop = FALSE]
  Io <- crossprod(S) / N
  ## Near the boundary the OPG information is not merely singular but close
  ## enough to it that solve() succeeds and returns nonsense, so test the
  ## conditioning rather than waiting for an error.
  rc <- tryCatch(rcond(Io), error = function(e) 0)
  if (!is.finite(rc) || rc < 1e-10) {
    stop("uhet_test(): the outer-product information matrix is numerically ",
      "singular (reciprocal condition number ", format(rc, digits = 3),
      "), so the first-step correction cannot be computed. This usually means ",
      "the fit is at or near a boundary; see sfa_diagnostics().",
      call. = FALSE
    )
  }
  Ii <- solve(Io)
  R <- S %*% t(Ii) ## r_i', one row per observation

  Jb <- numDeriv::jacobian(function(t) .ks_b_of_psi(t, Y, X, mn, inefdec), psi)
  G <- crossprod(Z, Jb) / N ## k x p
  U <- Z * res + R %*% t(G)
  A_corr <- crossprod(U) / N

  vc <- function(A) (Bi %*% A %*% Bi) / N
  stat <- function(V) {
    Vk <- V[keep, keep, drop = FALSE]
    Vi <- tryCatch(solve(Vk), error = function(e) NULL)
    if (is.null(Vi)) {
      return(NA_real_)
    }
    as.numeric(t(gam[keep]) %*% Vi %*% gam[keep])
  }
  Vc <- vc(A_corr)
  Vn <- vc(A_naive)
  dimnames(Vc) <- dimnames(Vn) <- list(colnames(Z), colnames(Z))
  Dc <- stat(Vc)
  Dn <- stat(Vn)

  out <- data.frame(
    test = c("two-step (corrected)", "two-step (naive)"),
    statistic = c(Dc, Dn), df = q,
    p.value = stats::pchisq(c(Dc, Dn), df = q, lower.tail = FALSE),
    stringsAsFactors = FALSE
  )
  attr(out, "gamma") <- gam[keep]
  attr(out, "se_corrected") <- sqrt(pmax(diag(Vc)[keep], 0))
  attr(out, "se_naive") <- sqrt(pmax(diag(Vn)[keep], 0))
  rownames(out) <- NULL
  out
}
