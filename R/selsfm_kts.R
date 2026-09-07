## Kumbhakar, Tsionas and Sipilainen (2009), J Prod Anal 31:151-161: joint
## estimation of TECHNOLOGY CHOICE and technical efficiency.
## See notes/code_history/selsfm.md.
##
## Different model from Greene (2010), which selsfm() otherwise fits, and the
## difference is the point. Greene has ONE frontier, observed only for the
## selected, with the selection correlated with the NOISE v. Here both regimes
## are observed, each with its OWN frontier and its own two scales, and the
## choice depends on the INEFFICIENCY itself -- u plays a "dual role", lowering
## output through the frontier and shifting the technology decision through the
## choice equation. That is why neither two-step order works: the choice
## equation cannot be a probit because u is unobserved, and fitting the
## frontiers first ignores the endogeneity of the choice.
##
##   y_i | I_i, u_i ~ N(x_i'beta_{I_i} - u_i, sigma_{v,I_i}^2)
##   P(I_i = 1 | u_i) = Phi(z_i'gamma + delta u_i)
##   u_i | I_i       ~ N+(0, sigma_{u,I_i}^2)
##
## delta is the parameter the paper exists for: delta > 0 means the less
## efficient are MORE likely to choose technology 1.

## Half-normal-weighted integral over u >= 0, mapped to [0,1] by the
## half-normal's own CDF: with w = 2 Phi(u/sigma) - 1 the weight becomes
## uniform, so int_0^inf g(u) f+(u|sigma) du = int_0^1 g(sigma qnorm((1+w)/2)) dw
## and Gauss-Legendre on the unit interval integrates it. Same device copsfm()
## uses, and it needs no truncation point.
.kts_nodes <- function(n_nodes) {
  gl <- .gauss_legendre_01(as.integer(n_nodes))
  list(z = stats::qnorm((1 + gl$nodes) / 2), w = gl$weights)
}

## The joint density of (y_i, I_i), the paper's equation (14).
##
## A(sigma^2)   = int_0^inf Phi(z'gamma + delta u) f+(u|sigma) du          (10)
## phi_i        = A(sigma_u0) / (1 + A(sigma_u0) - A(sigma_u1))            (9)
## J(sigma^2)   = int_0^inf f_N(y | x'beta_I - u, sigma_vI) Phi(.)^I
##                          (1-Phi(.))^(1-I) f+(u|sigma) du                (13)
## p(y_i, I_i)  = phi_i J(sigma_u1) + (1 - phi_i) J(sigma_u0)              (14)
##
## phi_i does not depend on u, which is what lets (11) split into (12): the
## marginal density of u is a MIXTURE of the two regimes' half-normals, so the
## observation's density is the same mixture of the two inner integrals.
.kts_loglik <- function(theta, y, X0, X1, Z, I, nd, S = 1) {
  k <- ncol(X0)
  m <- ncol(Z)
  b0 <- theta[1:k]
  b1 <- theta[(k + 1):(2 * k)]
  sv0 <- exp(theta[2 * k + 1]); sv1 <- exp(theta[2 * k + 2])
  su0 <- exp(theta[2 * k + 3]); su1 <- exp(theta[2 * k + 4])
  gam <- theta[(2 * k + 5):(2 * k + 4 + m)]
  del <- theta[2 * k + 5 + m]

  a <- as.numeric(Z %*% gam)                       # n
  ## u values on the quadrature grid, one column per node, per regime scale.
  u0 <- outer(rep(1, length(y)), su0 * nd$z)       # n x q
  u1 <- outer(rep(1, length(y)), su1 * nd$z)

  ## log Phi(a + delta u) and its complement, at both sets of nodes.
  lF <- function(uu) {
    arg <- pmin(pmax(a + del * uu, .SFA_CONSTANTS$CLIP_Z1_LOWER),
      .SFA_CONSTANTS$CLIP_Z1_UPPER)
    list(hi = stats::pnorm(arg, log.p = TRUE),
      lo = stats::pnorm(arg, log.p = TRUE, lower.tail = FALSE))
  }
  F0 <- lF(u0); F1 <- lF(u1)

  ## A(sigma^2): the choice probability marginalised over inefficiency.
  A0 <- as.numeric(exp(F0$hi) %*% nd$w)
  A1 <- as.numeric(exp(F1$hi) %*% nd$w)
  den <- 1 + A0 - A1
  den[abs(den) < .SFA_CONSTANTS$MIN_POSITIVE] <- .SFA_CONSTANTS$MIN_POSITIVE
  phi <- pmin(pmax(A0 / den, 0), 1)

  ## The frontier residual uses the regime's OWN coefficients and scale.
  mu <- ifelse(I == 1, as.numeric(X1 %*% b1), as.numeric(X0 %*% b0))
  sv <- ifelse(I == 1, sv1, sv0)
  ## NOT ifelse(I == 1, FF$hi, FF$lo). ifelse() returns the shape of its TEST,
  ## so a length-n test against these n x q matrices silently collapses to the
  ## FIRST COLUMN -- the choice probability would then be evaluated at one
  ## quadrature node instead of varying with u, which removes the coupling
  ## between u and the technology decision that this model is entirely about.
  ## The likelihood still integrated to 1 with that bug, because the selection
  ## factor became a constant in u and the two regimes still summed to one, so
  ## the properness check could not catch it. What caught it was comparing the
  ## fitted parameters against the truth on FRESH data.
  sel_rows <- function(FF) {
    m <- FF$hi
    if (any(I != 1)) m[I != 1, ] <- FF$lo[I != 1, , drop = FALSE]
    m
  }
  Jof <- function(uu, FF) {
    e <- (y - (mu - S * uu)) / sv
    lg <- -0.5 * e^2 - log(sv) - .SFA_CONSTANTS$LOG_SQRT_2PI + sel_rows(FF)
    as.numeric(exp(lg) %*% nd$w)
  }
  J0 <- Jof(u0, F0)
  J1 <- Jof(u1, F1)

  p <- phi * J1 + (1 - phi) * J0
  p <- pmax(p, .SFA_CONSTANTS$MIN_POSITIVE)
  log(p)
}

## E[u | y_i, I_i], by the same quadrature: the posterior mean of u under the
## mixture density the likelihood integrates, so the predictor and the
## likelihood cannot drift apart.
.kts_jlms <- function(theta, y, X0, X1, Z, I, nd, S = 1) {
  k <- ncol(X0); m <- ncol(Z)
  b0 <- theta[1:k]; b1 <- theta[(k + 1):(2 * k)]
  sv0 <- exp(theta[2 * k + 1]); sv1 <- exp(theta[2 * k + 2])
  su0 <- exp(theta[2 * k + 3]); su1 <- exp(theta[2 * k + 4])
  gam <- theta[(2 * k + 5):(2 * k + 4 + m)]; del <- theta[2 * k + 5 + m]
  a <- as.numeric(Z %*% gam)
  mu <- ifelse(I == 1, as.numeric(X1 %*% b1), as.numeric(X0 %*% b0))
  sv <- ifelse(I == 1, sv1, sv0)
  num <- 0; den <- 0
  for (j in seq_along(nd$z)) {
    for (reg in 1:2) {
      su <- if (reg == 1) su0 else su1
      uu <- su * nd$z[j]
      arg <- pmin(pmax(a + del * uu, .SFA_CONSTANTS$CLIP_Z1_LOWER),
        .SFA_CONSTANTS$CLIP_Z1_UPPER)
      lf <- ifelse(I == 1, stats::pnorm(arg, log.p = TRUE),
        stats::pnorm(arg, log.p = TRUE, lower.tail = FALSE))
      e <- (y - (mu - S * uu)) / sv
      wgt <- exp(-0.5 * e^2 - log(sv) + lf) * nd$w[j]
      num <- num + wgt * uu
      den <- den + wgt
    }
  }
  num / pmax(den, .SFA_CONSTANTS$MIN_POSITIVE)
}

## Fit the KTS model. Both regimes contribute, so unlike Greene's branch this
## uses every row: the likelihood is a joint density of (y, I), not a density
## for the selected subsample.
.selsfm_kts_fit <- function(selection, frontier, data, n_nodes, inefdec,
                            maxit.bobyqa, maxit.psoptim, maxit.optim,
                            start_val, PSopt, optHessian, Method, verbose,
                            rand.psoptim, call) {
  Start.Time <- Sys.time()
  S <- if (isTRUE(inefdec)) 1 else -1

  mf_s <- stats::model.frame(selection, data = data, na.action = stats::na.pass)
  mf_f <- stats::model.frame(frontier, data = data, na.action = stats::na.pass)
  Z <- stats::model.matrix(stats::terms(mf_s), mf_s)
  X <- stats::model.matrix(stats::terms(mf_f), mf_f)
  I <- as.numeric(stats::model.response(mf_s))
  y <- as.numeric(stats::model.response(mf_f))

  ok <- stats::complete.cases(cbind(y, I, X, Z))
  if (!all(ok)) {
    y <- y[ok]; I <- I[ok]; X <- X[ok, , drop = FALSE]; Z <- Z[ok, , drop = FALSE]
  }
  if (!all(I %in% c(0, 1))) {
    stop("selsfm(model_name = \"kts\"): the technology indicator must be binary.",
      call. = FALSE)
  }
  ## The Greene branch tolerates a missing y for the unselected because it never
  ## uses them. Here both regimes enter the likelihood, so a missing outcome is
  ## a dropped observation and the user should know which model they asked for.
  if (min(sum(I), sum(1 - I)) < ncol(X) + 2L) {
    stop("selsfm(model_name = \"kts\"): each technology needs more observations ",
      "than it has frontier parameters; got ", sum(1 - I), " and ", sum(I),
      " for ", ncol(X), " coefficients. This model fits a SEPARATE frontier to ",
      "each regime, unlike Greene's, which fits one.", call. = FALSE)
  }

  nd <- .kts_nodes(n_nodes)
  k <- ncol(X); m <- ncol(Z)

  ## Starting values: OLS within each regime for the two frontiers, the OLS
  ## residual sd split between the two components the way start_cs() does, and
  ## a probit for gamma with delta at zero -- delta = 0 is the model without the
  ## dual role, so the search starts from the hypothesis it exists to test.
  st_reg <- function(sel) {
    fit <- stats::lm.fit(X[sel, , drop = FALSE], y[sel])
    list(b = fit$coefficients, s = stats::sd(fit$residuals))
  }
  r0 <- st_reg(I == 0); r1 <- st_reg(I == 1)
  pb <- suppressWarnings(stats::glm.fit(Z, I,
    family = stats::binomial(link = "probit")))
  start_v <- c(r0$b, r1$b,
    log(max(0.6 * r0$s, 1e-3)), log(max(0.6 * r1$s, 1e-3)),
    log(max(0.8 * r0$s, 1e-3)), log(max(0.8 * r1$s, 1e-3)),
    pb$coefficients, 0)
  if (isTRUE(is.numeric(start_val))) start_v <- start_val
  start_v[!is.finite(start_v)] <- 0

  like.fn <- function(th) {
    ll <- tryCatch(.kts_loglik(th, y, X, X, Z, I, nd, S), error = function(e) NULL)
    if (is.null(ll) || any(!is.finite(ll))) return(1e12)
    -sum(ll)
  }

  lower <- c(rep(-Inf, 2 * k), rep(log(1e-6), 4), rep(-Inf, m + 1))
  Opt.Bobyqa <- opt.bobyqa(fn = like.fn, start_v = start_v,
    lower.bobyqa = lower, maxit.bobyqa = maxit.bobyqa, bob.TF = TRUE,
    verbose = verbose)
  start_v <- Opt.Bobyqa$start_v
  L1 <- lower.start(start_v, "NHN", differ = 1)
  Opt.Ps <- opt.psoptim(fn = like.fn, start_v, lower.psoptim = lower,
    rand.psoptim = rand.psoptim, upper.psoptim = Inf,
    maxit.psoptim = maxit.psoptim, psopt.TF = PSopt, rand.order = FALSE,
    verbose = verbose)
  start_v <- Opt.Ps$start_v
  Opt.Optim <- opt.optim(fn = like.fn, start_v = start_v, lower.optim = lower,
    upper.optim = Inf, maxit.optim = maxit.optim, opt.TF = optHessian,
    method = Method, optHessian = TRUE, verbose = verbose)
  opt <- Opt.Optim$opt
  th <- Opt.Optim$start_v

  nm <- c(paste0("t0.", colnames(X)), paste0("t1.", colnames(X)),
    "sigma_v0", "sigma_v1", "sigma_u0", "sigma_u1",
    paste0("choice.", colnames(Z)), "delta")
  se <- tryCatch(sqrt(pmax(diag(solve(opt$hessian)), 0)),
    error = function(e) rep(NA_real_, length(th)))
  ## The four scales are estimated as logs and reported on the natural scale,
  ## so their standard errors move with them: d(exp b)/db = exp(b).
  iss <- (2 * k + 1):(2 * k + 4)
  par <- th
  par[iss] <- exp(th[iss])
  se[iss] <- se[iss] * par[iss]

  out <- rbind(par, se, par / se)
  rownames(out) <- c("par", "st_err", "t-val")
  colnames(out) <- nm

  jlms <- .kts_jlms(th, y, X, X, Z, I, nd, S)
  results <- list(
    t(out), c(opt), end.time(Start.Time), start_v, "KTS", frontier, selection,
    jlms, exp(-jlms), like.fn, n_nodes, length(y), sum(I), S,
    out["par", ], out["st_err", ], out["t-val", ], call
  )
  class(results) <- "sfareg"
  names(results) <- c("out", "opt", "total_time", "start_v", "model_name",
    "formula", "selection_formula", "jlms", "efficiency", "objective",
    "n_nodes", "nobs", "n_selected", "S",
    "coefficients", "std.errors", "t.values", "call")
  results
}
