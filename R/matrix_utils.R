## ------------------------------------------------------------
## Helper: column-wise products of a matrix (replaces matrixStats::colProds)
##
## Inlined to drop the matrixStats dependency (only colProds was ever used
## from it, in the GTRE/TRE simulated-ML likelihoods). Folds the rows with
## element-wise multiplication down the columns -- this is byte-identical to
## matrixStats::colProds(method = "direct") (the default), i.e. the same
## left-to-right product order, verified with identical(); it is NOT the
## expSumLog method (which would differ at ~1e-15 and change results in the
## 15th digit). For the small row counts here (Ti = number of time periods,
## typically <=10-20) this pure-R fold is actually faster than the C call it
## replaces, since it avoids matrixStats' per-call setup overhead.
## ------------------------------------------------------------
.col_prods <- function(x) {
  x <- as.matrix(x)
  cp <- x[1, ]
  if (nrow(x) > 1) {
    for (r in 2:nrow(x)) cp <- cp * x[r, ]
  }
  cp
}


## Helper: column-wise (or vector) demeaning (replaces Jmisc::demean)
.demean <- function(x) {
  if (is.vector(x)) {
    return(x - mean(x))
  }
  apply(x, 2, function(col) col - mean(col))
}


## ------------------------------------------------------------
## Helper: validate the arguments common to every exported model fitter
##
## These functions previously accepted anything and failed deep inside the
## optimizer with messages that named none of the user's arguments -- a
## data.frame with no rows, for instance, surfaced as a subscript error from
## inside the likelihood. Fail early, name the argument, and say what was
## expected. `call. = FALSE` throughout: the internal call is noise to a user
## who typed sfm().
##
## Deliberately permissive about the data class. plm::pdata.frame, tibbles and
## data.tables all inherit from data.frame, and psfm() accepts a pdata.frame by
## design, so the test is inherits(data, "data.frame") rather than an identity
## check on class().
## ------------------------------------------------------------
.validate_sfa_call <- function(formula, data, fn, n_min = 1L,
                               maxit = list(), flags = list()) {
  if (missing(formula) || !inherits(formula, "formula")) {
    stop(fn, "(): `formula` must be a formula, e.g. y ~ x1 + x2. ",
      if (missing(formula)) {
        "It was not supplied."
      } else {
        paste0(
          "Got an object of class ",
          paste(class(formula), collapse = "/"), "."
        )
      },
      call. = FALSE
    )
  }
  if (missing(data) || is.null(data)) {
    stop(fn, "(): `data` must be supplied.", call. = FALSE)
  }
  if (!inherits(data, "data.frame")) {
    stop(fn, "(): `data` must be a data.frame (a tibble, data.table or ",
      "plm::pdata.frame is fine). Got an object of class ",
      paste(class(data), collapse = "/"), ".",
      call. = FALSE
    )
  }
  if (nrow(data) < n_min) {
    stop(fn, "(): `data` has ", nrow(data), " row",
      if (nrow(data) == 1L) "" else "s",
      "; at least ", n_min, " is required to fit a model.",
      call. = FALSE
    )
  }

  for (nm in names(maxit)) {
    v <- maxit[[nm]]
    if (!is.numeric(v) || length(v) != 1L || is.na(v) || v < 1) {
      stop(fn, "(): `", nm, "` must be a single positive number of iterations. ",
        "Got ", paste(utils::capture.output(utils::str(v)), collapse = " "), ".",
        call. = FALSE
      )
    }
  }
  for (nm in names(flags)) {
    v <- flags[[nm]]
    if (!is.logical(v) || length(v) != 1L || is.na(v)) {
      stop(fn, "(): `", nm, "` must be TRUE or FALSE. Got ",
        paste(utils::capture.output(utils::str(v)), collapse = " "), ".",
        call. = FALSE
      )
    }
  }
  invisible(TRUE)
}


## ------------------------------------------------------------
## Helpers: leave the caller's random number stream as we found it
##
## Several functions here accept a seed argument (`rand`, `rand.gtre`,
## `rand.psoptim`, `seed_offset`) and call set.seed() so their own draws are
## reproducible. Doing that overwrites .Random.seed in the global environment,
## which silently changes every random number the USER draws afterwards --
## their own seed no longer controls their session once they have fitted a
## model. Demonstrated before this was added: with set.seed(42) fixed by the
## caller, rnorm(3) returned 1.37096 -0.56470 0.36313 on its own but
## 0.99032 0.11227 1.14963 after an intervening data_gen_cs(rand = 7).
##
## Writing to the global environment at all is discouraged, so restoring it is
## the minimum obligation. Snapshot before seeding, restore on exit:
##
##     if (!is.null(rand)) {
##       .rng_state <- .rng_snapshot()
##       on.exit(.rng_restore(.rng_state), add = TRUE)
##       set.seed(rand)
##     }
##
## Only when a seed was actually supplied. With no seed the function should
## draw from, and advance, the caller's stream like any other R code.
##
## This does NOT change any result: set.seed(rand) still precedes the draws, so
## the values produced are identical. Only the state left behind differs.
## ------------------------------------------------------------
.rng_snapshot <- function() {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  } else {
    NULL
  }
}

.rng_restore <- function(state) {
  if (is.null(state)) {
    ## There was no stream before we seeded; leave none behind.
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  } else {
    assign(".Random.seed", state, envir = globalenv())
  }
  invisible(NULL)
}


## ------------------------------------------------------------
## Helper: physicists' Gauss-Hermite quadrature nodes/weights
##
## Golub-Welsch: nodes are eigenvalues of the (n x n) tridiagonal Jacobi
## matrix with zero diagonal and off-diagonal sqrt(i/2), weights are
## sqrt(pi) times the squared first component of each normalized
## eigenvector. Base R only (eigen()); no dependency on statmod/gaussquad.
## Used by .log_mvn_cdf_rank1() below.
## ------------------------------------------------------------
.gauss_hermite_nodes <- function(n = 64L) {
  n <- as.integer(n)
  if (n < 5L) stop("n must be at least 5.", call. = FALSE)
  i <- seq_len(n - 1L)
  J <- matrix(0, n, n)
  off <- sqrt(i / 2)
  J[cbind(i, i + 1L)] <- off
  J[cbind(i + 1L, i)] <- off
  eg <- eigen(J, symmetric = TRUE)
  nodes <- eg$values
  weights <- sqrt(pi) * (eg$vectors[1L, ]^2)
  ord <- order(nodes)
  list(nodes = nodes[ord], weights = weights[ord])
}


## ------------------------------------------------------------
## Helper: Gauss-Legendre quadrature nodes/weights on (0, 1)
##
## Same Golub-Welsch construction as above, for the Legendre weight: zero
## diagonal and off-diagonal i/sqrt(4i^2 - 1) on (-1, 1), then shifted and
## scaled to (0, 1) with weights halved so they sum to 1.
##
## Used by sfm()'s THT efficiency predictor, which integrates on the
## probability scale: for any g, E[g(u)] = int_0^1 g(Q(p)) dp with Q the
## truncated-t quantile function. Applying a fixed rule there turns one
## integrate() call per observation into a single matrix operation.
## ------------------------------------------------------------
.gauss_legendre_01 <- function(n = 64L) {
  n <- as.integer(n)
  if (n < 5L) stop("n must be at least 5.", call. = FALSE)
  i <- seq_len(n - 1L)
  J <- matrix(0, n, n)
  off <- i / sqrt(4 * i^2 - 1)
  J[cbind(i, i + 1L)] <- off
  J[cbind(i + 1L, i)] <- off
  eg <- eigen(J, symmetric = TRUE)
  ord <- order(eg$values)
  list(
    nodes = (eg$values[ord] + 1) / 2,
    weights = eg$vectors[1L, ord]^2
  )
}


## ------------------------------------------------------------
## Helper: Student-t--half-normal (tHN) composed-error density
##
## e = v - u with v ~ sigma_v * t_nu and u ~ |N(0, sigma_u^2)|, INDEPENDENT.
## Unlike THT (Tancredi 2002), where a single scale mixture is shared by both
## components and the composed error is a closed-form skew-t, here the two
## components have different tail behaviour and there is no closed form:
##
##     f(e) = int_0^Inf f_v(e + u) f_u(u) du
##
## Evaluated by fixed Gauss-Legendre quadrature after u = sigma_u * s, so the
## range s in [0, 8] is eight half-normal standard deviations whatever sigma_u
## is, and the half-normal factor in s units is sqrt(2/pi)*exp(-s^2/2).
##
## NODE COUNT IS NOT FIXED, and this matters. The t factor
## dt((e + sigma_u*s)/sigma_v) is a peak in s of width sigma_v/sigma_u = 1/lambda,
## so a rule that resolves it at lambda = 3 does not resolve it at lambda = 30.
## Measured against adaptive integrate() with the panel boundaries placed at the
## peak, a fixed 96-node rule carries relative error 4e-2 at lambda = 20 and
## 6e-1 at lambda = 62 -- and the model really does visit that region (a
## normal-like optimum near lambda = 63-72 is documented on real data).
##
## The total mass check does NOT catch this: the error redistributes across e
## and integrates away, so int f(e) de still reads 1.000 while f is 40% wrong
## pointwise. Do not "verify" this density by integrating it alone.
##
## Empirically ~21 nodes per unit of lambda holds the relative error near 1e-6,
## so the count scales with lambda and is capped. Worst case over
## lambda in [0.05, 150] and nu in [2.05, 200] is 4e-5, reached only in the
## degenerate lambda > 100 corner; it is ~1e-6 or better for lambda <= 70.
## ------------------------------------------------------------
.thn_m_for <- function(lambda) {
  ## The optimizer WILL visit sigma_v -> 0, so lambda arrives as Inf or NaN.
  ## Clamp before the arithmetic: 64*ceiling(24*lambda/64) overflows the
  ## integer range for large lambda and as.integer() then returns NA, which
  ## surfaces far away as "missing value where TRUE/FALSE needed" inside the
  ## quadrature.
  lambda <- suppressWarnings(as.numeric(lambda))
  if (!is.finite(lambda) || lambda < 1) lambda <- 1
  lambda <- min(lambda, 200)
  max(96L, min(4096L, as.integer(64 * ceiling(24 * lambda / 64))))
}

.thn_nodes <- local({
  cache <- list()
  function(m = 96L, s_max = 8) {
    key <- paste(m, s_max, sep = "_")
    if (is.null(cache[[key]])) {
      gl <- .gauss_legendre_01(m)
      ## .gauss_legendre_01() weights sum to 1 on (0,1); rescale to [0, s_max].
      cache[[key]] <<- list(s = s_max * gl$nodes, w = s_max * gl$weights)
    }
    cache[[key]]
  }
})

## Returns the quadrature pieces shared by the density and the efficiency
## predictor: the observation-by-node matrix of f_v values and the node weights
## already multiplied by the half-normal factor.
.thn_parts <- function(e, sigma_v, sigma_u, nu) {
  nd <- .thn_nodes(.thn_m_for(sigma_u / sigma_v))
  s <- nd$s
  ## The t density is written out rather than calling dt() on the whole
  ## observations-by-nodes matrix. Its normalising constant is a scalar, so
  ## only exp(-(nu+1)/2 * log1p(x^2/nu)) has to touch every element; dt()
  ## redoes the lgamma work per call. Identical to 1e-12 and 2.7x faster at
  ## m = 128, n = 1000, which matters because this is the whole cost of the
  ## model -- it is evaluated thousands of times per fit.
  x <- outer(as.numeric(e), sigma_u * s, "+") / sigma_v
  kc <- exp(lgamma((nu + 1) / 2) - lgamma(nu / 2)) / (sqrt(nu * pi) * sigma_v)
  list(
    fv = kc * exp(-((nu + 1) / 2) * log1p(x * x / nu)),
    wu = nd$w * sqrt(2 / pi) * exp(-s^2 / 2),
    s = s
  )
}

## Largest lambda = sigma_u/sigma_v the quadrature is validated at. At the
## node-count rule above this is where m hits its 4096 cap while still holding
## ~24 nodes per unit of lambda; past it the relative error grows without
## bound and the density can no longer be trusted.
##
## Refusing the region is not cosmetic. Profiling a single n = 1000 fit showed
## the optimizer driving sigma_v to zero and evaluating at lambda up to 2e16:
## only 3.8% of calls, but they pinned the node count at its cap and consumed
## most of a 259-second fit. They are also exactly the calls whose answer is
## meaningless. sigma_v -> 0 with sigma_u fixed is the deterministic-frontier
## degeneracy (e = -u, so any observation with a positive residual has density
## zero), not a competing optimum, so returning a large finite penalty both
## bounds the cost and keeps the search in the region where the density is
## real. Note this bounds sigma_v away from 0 RELATIVE to sigma_u; it does not
## bound sigma_u away from 0, which is the collapse the model is supposed to be
## able to report.
.THN_LAMBDA_MAX <- 170

.log_d_thn <- function(e, sigma_v, sigma_u, nu) {
  n <- length(e)
  lam <- suppressWarnings(sigma_u / sigma_v)
  if (!is.finite(lam) || lam > .THN_LAMBDA_MAX) {
    return(rep(-1e6, n))
  }
  p <- .thn_parts(e, sigma_v, sigma_u, nu)
  log(pmax(as.numeric(p$fv %*% p$wu), .Machine$double.xmin))
}

## E[exp(-u)|e] and E[u|e] by the same node set -- numerator and denominator
## share it, so the Bayes rule is exact up to the quadrature error above.
.thn_eff <- function(e, sigma_v, sigma_u, nu) {
  p <- .thn_parts(e, sigma_v, sigma_u, nu)
  den <- pmax(as.numeric(p$fv %*% p$wu), .Machine$double.xmin)
  u <- sigma_u * p$s
  list(
    exp_u_hat = pmin(pmax(as.numeric(p$fv %*% (p$wu * exp(-u))) / den, 0), 1),
    u_hat = pmax(as.numeric(p$fv %*% (p$wu * u)) / den, 0)
  )
}


## ------------------------------------------------------------
## Helper: log MVN CDF for the rank-one/equicorrelated covariance
## Sigma = I_T + c*11' (replaces mnormt::pmnorm() for this special case)
##
## For Sigma of this form, X = Z + sqrt(c)*W*1 with Z ~ N(0,I_T), W ~ N(0,1)
## independent, so Pr(X <= upper) = E_W[ prod_t Phi(upper_t - sqrt(c)*W) ] --
## an exact 1-dimensional integral, here evaluated by Gauss-Hermite
## quadrature (default 64 nodes) rather than mnormt::pmnorm()'s general
## T-dimensional (stochastic, Genz-algorithm) integration. This same
## covariance structure arises in psfm()'s TFE model (see like.tfe() in
## psfm.R). Matches mnormt::pmnorm() to within pmnorm's own default
## tolerance (~1e-3 on the log scale) while being ~30x faster and
## deterministic (no Monte Carlo noise).
## ------------------------------------------------------------
.log_mvn_cdf_rank1 <- function(upper, c, gh) {
  upper <- as.numeric(upper)
  if (length(upper) == 0L) {
    return(0)
  }
  if (c < 1e-14) {
    return(sum(pnorm(upper, log.p = TRUE)))
  }

  sc <- sqrt(c)
  s_nodes <- sqrt(2) * gh$nodes
  log_terms <- numeric(length(s_nodes))
  for (r in seq_along(s_nodes)) {
    log_terms[r] <- log(gh$weights[r]) + sum(pnorm(upper - sc * s_nodes[r], log.p = TRUE))
  }
  m <- max(log_terms)
  (m + log(sum(exp(log_terms - m)))) - 0.5 * log(pi)
}


## ------------------------------------------------------------
## Helper: log density of the within-transformed disturbance vector
## for psfm()'s TFE model (replaces mnormt::dmnorm() for this special case)
##
## Closed form for the density of an (T-1)-vector eps_star with covariance
## sigma2*(I_m - (1/T)*11'_m), m = T-1, using det = 1/T and inverse =
## I_m + 11'_m (both by the matrix determinant lemma / Sherman-Morrison,
## exact given the specific 1/T normalizer used in like.tfe()'s E_t1).
## Matches mnormt::dmnorm() to floating-point precision.
## ------------------------------------------------------------
.log_within_mvn_density <- function(eps_star, sigma2) {
  m <- length(eps_star)
  q <- sum(eps_star^2) + sum(eps_star)^2
  -0.5 * m * log(2 * pi) - 0.5 * m * log(sigma2) + 0.5 * log(m + 1) - q / (2 * sigma2)
}


## ------------------------------------------------------------
## Helper: cross-sectional intercept-only normal-half-normal (Aigner,
## Lovell & Schmidt 1977) MLE, with Hessian-based standard errors.
##
## Used by psfm()'s GTRE_SEQ1 branch to decompose composite residuals
## (epsilon_hat, alpha_hat) into their normal/half-normal variance
## components. sfm()'s own NHN likelihood can't be reused directly here
## because sfm()'s formula-parsing path does not currently support an
## intercept-only (y~1) formula, so this reimplements the same closed-form
## density directly -- the T_i=1 special case of the error-components-
## frontier likelihood used for PL80/BC92, parameterized the same way,
## directly over (sigma_v, sigma_u, intercept) rather than
## (lambda=sigma_u/sigma_v, sigma).
## ------------------------------------------------------------
.fit_nhn_intercept <- function(y, inefdec_n = 1, maxit.bobyqa = 200, maxit.psoptim = 20,
                               maxit.optim = 100, PSopt = FALSE, verbose = FALSE) {
  like_fn <- function(x) {
    sigma_v <- x[1]
    sigma_u <- x[2]
    intercept <- x[3]
    if (!is.finite(sigma_v) || !is.finite(sigma_u) || sigma_v <= 0 || sigma_u <= 0) {
      return(1e12)
    }
    eps <- inefdec_n * (y - intercept)
    sigmaSq <- sigma_v * sigma_v + sigma_u * sigma_u
    z <- -sigma_u * eps / (sigma_v * sqrt(sigmaSq))
    ll <- log(2) - 0.5 * log(2 * pi) - 0.5 * log(sigmaSq) - (eps * eps) / (2 * sigmaSq) + pnorm(z, log.p = TRUE)
    -sum(ll[is.finite(ll)])
  }
  start_v <- c(sd(y) / sqrt(2), sd(y) / sqrt(2), mean(y))
  lower_v <- c(.SFA_CONSTANTS$MIN_POSITIVE, .SFA_CONSTANTS$MIN_POSITIVE, -Inf)

  Opt.Bobyqa <- opt.bobyqa(
    fn = like_fn, start_v = start_v, lower.bobyqa = lower_v,
    maxit.bobyqa = maxit.bobyqa, bob.TF = TRUE, verbose = verbose
  )
  start_v <- Opt.Bobyqa$start_v

  differ <- 2
  Opt.Psoptim <- opt.psoptim(
    fn = like_fn, start_v, lower.psoptim = c(lower_v[1:2], start_v[3] - differ),
    upper.psoptim = c(start_v[1:2] + differ, start_v[3] + differ),
    maxit.psoptim = maxit.psoptim, psopt.TF = PSopt, verbose = verbose
  )
  start_v <- Opt.Psoptim$start_v

  Opt.Optim <- opt.optim(
    fn = like_fn, start_v = start_v, lower.optim = lower_v, upper.optim = rep(Inf, 3),
    maxit.optim = maxit.optim, opt.TF = TRUE, method = "L-BFGS-B", optHessian = TRUE, verbose = verbose
  )
  opt <- Opt.Optim$opt
  list(par = opt$par, hessian = opt$hessian, value = opt$value)
}


## ------------------------------------------------------------
## Helper: time-decay pattern B_it for psfm()'s error-components-frontier
## models (PL80/BC92/K1990/K1990modified), y_it = x_it'beta + v_it - B_it*u_i.
##
## Centralizes the one thing that differs across these four models -- the
## rest of the closed-form log-likelihood in psfm.R's "PL80"/"BC92"/"K1990"/
## "K1990modified" branch is identical for all of them (it was already
## written in terms of a general B_it before K1990/K1990modified existed).
##
##   PL80:          B_it = 1                                  (no decay)
##   BC92:          B_it = exp(-eta*(t - Tref))                (Battese &
##                  Coelli 1992; Tref = the last period in the WHOLE panel,
##                  not per-firm -- see psfm.R's PL80/BC92 branch comment
##                  for how that convention was pinned down against
##                  frontier::sfa())
##   K1990:         B_it = 1 / (1 + exp(b*t + c*t^2))          (Kumbhakar
##                  1990; matches npsf::sf()'s "K1990" pattern)
##   K1990modified: B_it = 1 + d*(t - T_i) + e*(t - T_i)^2     (Kumbhakar
##                  1990's modified pattern; T_i is each FIRM's OWN last
##                  observed period -- unlike BC92's whole-panel Tref, this
##                  is simply max(t_i) since t_i already holds only this
##                  firm's own observed periods, so no separate whole-panel
##                  reference is needed or used here)
##
## `t_i` is the vector of time periods observed for one firm; `decay_par` is
## the model's extra parameter(s) beyond (sigma_v, sigma_u, beta) -- numeric(0)
## for PL80, length 1 for BC92, length 2 for K1990/K1990modified. `Tref` is
## only used by BC92.
## ------------------------------------------------------------
.build_Bit <- function(model_name, t_i, Tref = NA_real_, decay_par = numeric(0)) {
  switch(model_name,
    "PL80" = rep(1, length(t_i)),
    "BC92" = exp(-decay_par[1] * (t_i - Tref)),
    "K1990" = 1 / (1 + exp(decay_par[1] * t_i + decay_par[2] * t_i^2)),
    "K1990modified" = {
      Ti_own <- max(t_i)
      1 + decay_par[1] * (t_i - Ti_own) + decay_par[2] * (t_i - Ti_own)^2
    },
    stop("Unknown decay model_name: ", model_name, call. = FALSE)
  )
}


## ------------------------------------------------------------
## Helper: pre-estimation collinearity diagnostic for panel models.
##
## psfm()'s random-effects models (GTRE, TRE, and the _Z/SEQ variants) get
## their starting values from plm(formula, data, effect = "individual",
## model = "random"). plm's Swamy-Arora error-components step inverts the
## crossproduct of the BETWEEN-individual design matrix. That matrix can be
## rank deficient even when the pooled design the user actually specified is
## full rank, and when it is, the failure surfaces as
##
##   Error in solve.default(crossprod(ZBeta)) :
##     Lapack routine dgesv: system is exactly singular: U[18,18] = 0
##
## which is impossible to interpret without knowing psfm() calls plm()
## internally. The classic trigger is time dummies: averaging a year
## indicator within a firm gives (number of times that firm is observed in
## that year)/T_i, which is CONSTANT across firms in a balanced panel and
## near-constant in many unbalanced ones, so the year columns collapse onto
## the intercept in the between dimension while remaining perfectly
## estimable pooled.
##
## Note this check is worth running even when plm does NOT error. In the
## reference case here the between matrix has rank 3 of 8 with a condition
## number around 1e49; depending on the exact pattern LAPACK may return a
## garbage inverse instead of failing, which yields meaningless starting
## values rather than a visible error -- the worse of the two outcomes.
##
## Returns pooled and between rank, the offending column names (via the QR
## pivot ordering), and the between condition number.
## ------------------------------------------------------------
.check_collinearity <- function(formula, data, individual, tol = 1e-8) {
  vars_needed <- intersect(unique(c(all.vars(formula), individual)), names(data))
  df_est <- data[stats::complete.cases(as.data.frame(data)[, vars_needed, drop = FALSE]), ,
    drop = FALSE
  ]

  X <- tryCatch(stats::model.matrix(formula, data = df_est), error = function(e) NULL)
  if (is.null(X) || !ncol(X)) {
    return(NULL)
  }

  qx <- qr(X, tol = tol)
  pooled_drop <- if (qx$rank < ncol(X)) colnames(X)[qx$pivot[(qx$rank + 1L):ncol(X)]] else character(0)

  ## Between-individual means of each column -- the design plm's
  ## error-components step actually inverts.
  idv <- as.character(df_est[[individual]])
  if (length(idv) != nrow(X)) {
    return(NULL)
  }
  Z <- rowsum(as.data.frame(X), group = idv, reorder = FALSE)
  Z <- as.matrix(Z / as.numeric(table(idv)[rownames(Z)]))

  qz <- qr(Z, tol = tol)
  between_drop <- if (qz$rank < ncol(Z)) colnames(Z)[qz$pivot[(qz$rank + 1L):ncol(Z)]] else character(0)

  list(
    pooled_rank = qx$rank, pooled_cols = ncol(X), pooled_drop = pooled_drop,
    between_rank = qz$rank, between_cols = ncol(Z), between_drop = between_drop,
    between_condition_number = tryCatch(kappa(Z), error = function(e) NA_real_)
  )
}

## Human-readable rendering of the above, used by psfm()'s error/warning.
.collinearity_message <- function(chk, action) {
  paste0(
    "Collinearity detected in the random-effects starting-value regression ",
    "used to initialize this model.\n",
    "  The requested frontier design matrix has ", chk$pooled_cols,
    " columns (rank ", chk$pooled_rank, ").\n",
    "  The between-individual starting-value matrix has ", chk$between_cols,
    " columns but rank ", chk$between_rank,
    if (is.finite(chk$between_condition_number)) {
      paste0(" (condition number ", format(chk$between_condition_number, digits = 3), ")")
    } else {
      ""
    },
    ".\n",
    "  Linearly dependent in the starting-value step: ",
    paste(chk$between_drop, collapse = ", "), "\n",
    switch(action,
      "error" = paste0(
        "  Estimation stopped (collinear_action = \"error\"). Options: simplify the\n",
        "  time controls (e.g. broader period groups or a time trend), supply\n",
        "  start_val manually, or use collinear_action = \"start_only\" to keep the\n",
        "  requested formula and drop these columns from initialization only."
      ),
      "warn_drop" = paste0(
        "  These columns were DROPPED FROM THE MODEL (collinear_action = \"warn_drop\")."
      ),
      "start_only" = paste0(
        "  The likelihood is still estimated using the requested formula; these\n",
        "  columns were dropped only from the starting-value regression, and their\n",
        "  starting values were taken from a pooled OLS fit."
      )
    )
  )
}


## Map offending MODEL-MATRIX column names back to the FORMULA TERMS that
## produced them. A factor term like `yr` expands to columns yr2, yr3, ...;
## dropping the term requires the term label, not the column name. A term is
## dropped only if ALL of its columns are implicated, so a single aliased
## level never silently removes an entire factor.
.terms_for_columns <- function(formula, data, cols) {
  if (!length(cols)) {
    return(character(0))
  }
  X <- tryCatch(stats::model.matrix(formula, data = data), error = function(e) NULL)
  if (is.null(X)) {
    return(character(0))
  }
  asg <- attr(X, "assign")
  tl <- attr(stats::terms(formula), "term.labels")
  if (is.null(asg) || !length(tl)) {
    return(character(0))
  }
  hit <- match(cols, colnames(X))
  hit <- hit[!is.na(hit)]
  if (!length(hit)) {
    return(character(0))
  }
  out <- character(0)
  for (k in unique(asg[hit])) {
    if (k == 0L) next ## intercept, never dropped here
    cols_k <- which(asg == k)
    if (all(cols_k %in% hit)) out <- c(out, tl[k])
  }
  out
}

## Re-expand a reduced starting-value coefficient vector to the full set of
## requested regressors, filling anything the reduced regression could not
## identify with its pooled OLS estimate (falling back to 0 if even that is
## unavailable). Keeps the starting vector conformable with the likelihood's
## parameter layout, which is what the optimizer scaffold expects.
.expand_start_beta <- function(beta_raw, x_vars_vec, formula, data) {
  if (is.null(x_vars_vec) || !length(x_vars_vec)) {
    return(beta_raw)
  }
  out <- stats::setNames(rep(NA_real_, length(x_vars_vec)), x_vars_vec)
  common <- intersect(names(beta_raw), x_vars_vec)
  out[common] <- beta_raw[common]
  if (anyNA(out)) {
    pooled <- tryCatch(stats::coef(stats::lm(formula, data = as.data.frame(data))),
      error = function(e) NULL
    )
    miss <- names(out)[is.na(out)]
    if (!is.null(pooled)) {
      got <- intersect(miss, names(pooled))
      out[got] <- pooled[got]
    }
    out[is.na(out)] <- 0
  }
  out
}


## ------------------------------------------------------------
## Helper: coerce user data into the panel form psfm() needs.
##
## psfm()'s downstream machinery (data_proc(), start_panel()'s plm() calls,
## the per-firm splits inside every panel likelihood) assumes the panel index
## is already attached to `data`. Rather than requiring callers to build a
## plm::pdata.frame themselves -- and failing with an opaque "empty model"
## error when they don't -- this attaches the index from the `individual`
## and `time` arguments psfm() already takes.
##
## Also normalizes the index column TYPE. plm stores the index as factors,
## and a numeric time column silently becomes a factor whose levels sort
## lexicographically ("10" < "2"), which reorders periods for any panel with
## more than nine of them. Converting through character preserves the numeric
## ordering.
## ------------------------------------------------------------
.as_panel_data <- function(data, individual, time = NULL) {
  if (missing(individual) || is.null(individual)) {
    stop("psfm() needs `individual`: the name of the column identifying each ",
      "cross-sectional unit (firm, region, ...).",
      call. = FALSE
    )
  }

  ## Already indexed: leave it alone, so existing pdata.frame callers keep
  ## their exact current behaviour.
  if (inherits(data, "pdata.frame")) {
    return(data)
  }

  data <- as.data.frame(data)

  if (!(is.character(individual) && length(individual) == 1L)) {
    stop("`individual` must be a single column name given as a character string.",
      call. = FALSE
    )
  }
  if (!(individual %in% names(data))) {
    stop("`individual` column '", individual, "' was not found in `data`. ",
      "Available columns: ", paste(names(data), collapse = ", "),
      call. = FALSE
    )
  }

  has_time <- !is.null(time) && !identical(time, FALSE)
  if (has_time) {
    if (!(is.character(time) && length(time) == 1L)) {
      stop("`time` must be a single column name given as a character string.",
        call. = FALSE
      )
    }
    if (!(time %in% names(data))) {
      stop("`time` column '", time, "' was not found in `data`. ",
        "Available columns: ", paste(names(data), collapse = ", "),
        call. = FALSE
      )
    }
  }

  ## Sort numerically where the column is numeric, so that period 10 follows
  ## period 9 rather than period 1.
  .norm_index <- function(v) {
    if (is.factor(v)) v <- as.character(v)
    if (is.numeric(v)) {
      return(factor(v, levels = sort(unique(v))))
    }
    num <- suppressWarnings(as.numeric(v))
    if (!anyNA(num)) {
      return(factor(v, levels = as.character(sort(unique(num)))))
    }
    factor(v)
  }

  data[[individual]] <- .norm_index(data[[individual]])
  if (has_time) data[[time]] <- .norm_index(data[[time]])

  idx <- if (has_time) c(individual, time) else individual
  out <- tryCatch(plm::pdata.frame(data, index = idx),
    error = function(e) {
      stop("Could not construct a panel index from `individual`='",
        individual,
        if (has_time) paste0("' and `time`='", time) else "",
        "': ", conditionMessage(e),
        ". Check for duplicate individual-time pairs.",
        call. = FALSE
      )
    }
  )
  out
}


## ------------------------------------------------------------
## Helper: analytic gradient of the normal/half-normal negative summed
## log-likelihood, in sfm()'s own (lambda, sigma, beta) parameterization.
##
## With eps_i = s*(y_i - x_i'beta) (s = +1 production, -1 cost),
##   a_i = eps_i/sigma,  b_i = -lambda*a_i,  m_i = phi(b_i)/Phi(b_i),
##   l_i = log2 - log(sigma) + log phi(a_i) + log Phi(b_i)
## gives
##   dl/dlambda  = -a_i m_i
##   dl/dsigma   = (a_i^2 - 1 + lambda a_i m_i)/sigma
##   dl/dbeta_j  = (s x_ij/sigma)(a_i + lambda m_i)
## and this returns the gradient of the NEGATIVE SUM, matching the sign
## convention every optimizer in opts.R expects.
##
## The inverse Mills ratio is formed as exp(dnorm(log) - pnorm(log.p)) rather
## than as a raw ratio: for a very efficient unit b_i is far into the left
## tail where Phi(b_i) underflows to 0 and the direct ratio is 0/0.
##
## Verified against numDeriv::grad at multiple parameter points to ~3e-8
## (numDeriv's own accuracy). Supplying this to nlminb roughly halves the
## runtime again relative to a numerically differentiated gradient.
## ------------------------------------------------------------
.grad_nhn <- function(x, Y, data_i_vars, inefdec_n = 1) {
  lam <- x[1]
  sig <- x[2]
  beta <- x[-c(1, 2)]
  if (!is.finite(lam) || !is.finite(sig) || lam <= 0 || sig <= 0) {
    return(rep(0, length(x)))
  }
  Xm <- as.matrix(data_i_vars)
  eps <- as.numeric(inefdec_n * (Y - Xm %*% beta))
  a <- eps / sig
  b <- -lam * a
  m <- exp(dnorm(b, log = TRUE) - pnorm(b, log.p = TRUE))
  g <- c(
    sum(a * m),
    -sum((a^2 - 1 + lam * a * m) / sig),
    -colSums((inefdec_n * Xm / sig) * (a + lam * m))
  )
  g[!is.finite(g)] <- 0
  g
}


## ------------------------------------------------------------
## Helper: conditional inefficiency predictors.
##
## For every composed-error model in this package the posterior of the
## one-sided component given the composed residual is a normal truncated at
## zero, u | e ~ N+(mu_star, sigma_star^2). Only (mu_star, sigma_star) differ
## by model:
##
##   NHN  mu_star = -sigma_u^2 e / sigma^2          sigma_star = sigma_u sigma_v / sigma
##   NTN  mu_star = (sigma_v^2 mu - sigma_u^2 e)/sigma^2   (same sigma_star)
##   NE   mu_star = -e - sigma_v^2/sigma_u          sigma_star = sigma_v
##
## Given those, the two standard predictors are model-free, so they live here
## once instead of being re-derived per branch:
##
##   .te_battese_coelli()  E[exp(-u)|e]  (Battese & Coelli 1988) -- the
##                         technical efficiency score reported as exp_u_hat
##   .jlms_u()             E[u|e]        (Jondrow et al. 1982)
##
## Both are evaluated in logs (pnorm(log.p = TRUE)) rather than by flooring a
## ratio of raw normal CDFs. In the far tail -- which is exactly where a very
## efficient or very inefficient unit sits -- Phi(z) underflows to 0 and the
## naive ratio becomes 0/0; the log form stays finite and accurate.
## ------------------------------------------------------------
## log(Phi(b) - Phi(a)) for b > a, computed without catastrophic cancellation.
##
## Needed by the doubly-truncated (uniform-u) predictors. Evaluating
## pnorm(b) - pnorm(a) directly is fine in the middle of the distribution but
## fails when BOTH limits sit far out in the same tail: each CDF rounds to the
## same value and the difference underflows to zero, so a ratio built on it
## blows up. Working in whichever tail is better conditioned and subtracting
## in log space keeps full precision -- without this, ~15% of uniform-model
## efficiency scores came back pinned at the 1.0 clamp instead of their true
## value.
.log_pnorm_diff <- function(a, b) {
  use_upper <- (a + b) > 0
  la <- ifelse(use_upper, pnorm(-b, log.p = TRUE), pnorm(a, log.p = TRUE))
  lb <- ifelse(use_upper, pnorm(-a, log.p = TRUE), pnorm(b, log.p = TRUE))
  lb + log(-expm1(pmin(la - lb, -.Machine$double.eps)))
}

.te_battese_coelli <- function(mu_star, sigma_star) {
  sigma_star <- pmax(sigma_star, .SFA_CONSTANTS$MIN_POSITIVE)
  z <- mu_star / sigma_star
  out <- exp(-mu_star + 0.5 * sigma_star^2 +
    pnorm(z - sigma_star, log.p = TRUE) - pnorm(z, log.p = TRUE))
  pmin(pmax(out, 0), 1)
}

.jlms_u <- function(mu_star, sigma_star) {
  sigma_star <- pmax(sigma_star, .SFA_CONSTANTS$MIN_POSITIVE)
  z <- mu_star / sigma_star
  pmax(mu_star + sigma_star * exp(dnorm(z, log = TRUE) - pnorm(z, log.p = TRUE)), 0)
}


## ------------------------------------------------------------
## Helpers: Greene (2005) true fixed effects stochastic frontier
##
## The TFE model is the ordinary normal/half-normal composed-error frontier
## with one intercept per individual estimated jointly with everything else:
##
##   y_it = alpha_i + x_it' beta + v_it - u_it,   v ~ N(0, sigma_v^2),
##                                                u ~ N+(0, sigma_u^2)
##
##   ln L = sum_it [ ln 2 - ln sigma + ln phi(e_it/sigma)
##                   + ln Phi(-lambda e_it/sigma) ],   e_it = y_it - alpha_i - x_it'beta
##
## psfm() estimates it as a PROFILE likelihood in (lambda, sigma, beta), with
## the N firm effects concentrated out at every evaluation. Two reasons to
## profile rather than hand all N + K + 2 parameters to the optimizer:
##
##   1. The parameter vector, the `out` matrix layout and psfm_bootstrap()'s
##      row-layout assumptions all stay identical to "TFE_WMLE" (lambda|gamma,
##      sig, beta), so nothing downstream needs a special case.
##   2. optim()'s numerical Hessian over 2 + K parameters is usable; over
##      2 + K + N it is neither affordable nor well conditioned. For the
##      concentrated parameters the profile Hessian is the right one anyway --
##      it equals the inverse of the corresponding block of the inverse of the
##      full Hessian.
##
## NOTE this is the estimator subject to the incidental-parameters problem
## (sigma_u is biased upward at small T); that is a property of the model, not
## of this implementation, and is precisely why Chen, Schmidt and Wang's
## within MLE ("TFE_WMLE") also exists in this package.
## ------------------------------------------------------------

## Group sums by an integer group id in 1..ngroups, independent of the order in
## which the groups happen to appear in the data. rowsum(reorder = TRUE) would
## sort the group labels as CHARACTER ("10" before "2"), which silently
## misassigns firm effects; reorder = FALSE returns first-appearance order,
## which is only 1..ngroups if the data happen to be sorted. Reading the
## returned rownames back as integers is correct under either arrangement.
.gsum <- function(x, gid, ngroups) {
  s <- rowsum(as.numeric(x), gid, reorder = FALSE)
  out <- numeric(ngroups)
  out[as.integer(rownames(s))] <- s[, 1L]
  out
}

## Concentrated firm effects for the TFE profile likelihood.
##
## `r` is the production-oriented residual inefdec_n*(y - x'beta), so the
## composed error is e_it = r_it - alpha_i. Holding (lambda, sigma) fixed the
## objective is separable across individuals and STRICTLY CONCAVE in alpha_i
## (phi and Phi are both log-concave and e_it is affine in alpha_i), so each
## alpha_i is the unique root of a strictly decreasing score:
##
##   dL_i/dalpha = sum_t [ e_it/sigma^2 + (lambda/sigma) m_it ]
##   d2L_i/dalpha^2 = sum_t [ -1/sigma^2 - (lambda^2/sigma^2) m_it (a_it + m_it) ]
##
## with a_it = -lambda e_it/sigma and m_it = phi(a_it)/Phi(a_it).
##
## Solved by safeguarded Newton rather than a nested optimize() call: the
## profile likelihood is differentiated numerically by optim() to get standard
## errors, and a nested optimizer that stops at its own (loose, and
## parameter-dependent) tolerance injects exactly the kind of small
## non-smoothness that makes those finite differences meaningless. Newton on a
## concave score converges quadratically to machine precision, so the profile
## is smooth to the accuracy the outer Hessian needs.
##
## The groupwise mean of `r` is always a valid lower bracket: there the first
## score term sums to zero by construction and the second is strictly positive.
.tfe_alpha_profile <- function(r, gid, ngroups, lambda, sigma,
                               tol = 1e-11, maxit = 100L) {
  if (!is.finite(lambda) || !is.finite(sigma) || lambda <= 0 || sigma <= 0) {
    return(NULL)
  }

  s2 <- sigma * sigma
  ls <- lambda / sigma

  score_curv <- function(alpha) {
    e <- r - alpha[gid]
    a <- -ls * e
    m <- exp(dnorm(a, log = TRUE) - pnorm(a, log.p = TRUE))
    m[!is.finite(m)] <- 0
    ## m(a + m) is the truncated-normal variance factor, in [0, 1] on paper.
    ## Far out in the upper tail the two terms nearly cancel and it can come
    ## back slightly negative or above 1 in floating point; clamping keeps the
    ## Newton step a descent direction rather than letting rounding flip its
    ## sign.
    w <- m * (a + m)
    w[!is.finite(w)] <- 0
    w <- pmin(pmax(w, 0), 1)
    list(
      g = .gsum(e / s2 + ls * m, gid, ngroups),
      h = .gsum(-1 / s2 - (ls * ls) * w, gid, ngroups)
    )
  }

  n_i <- .gsum(rep(1, length(r)), gid, ngroups)
  lo <- .gsum(r, gid, ngroups) / pmax(n_i, 1) ## score > 0 here, by construction
  hi <- lo
  need <- rep(TRUE, ngroups) ## still to be bracketed from above
  step <- max(4 * sigma * (1 + lambda), 1e-6)
  for (k in seq_len(80L)) {
    hi[need] <- hi[need] + step
    need <- score_curv(hi)$g >= 0
    if (!any(need)) break
    step <- step * 2
  }

  alpha <- 0.5 * (lo + hi)
  for (it in seq_len(maxit)) {
    sc <- score_curv(alpha)
    pos <- sc$g > 0
    lo[pos] <- pmax(lo[pos], alpha[pos])
    hi[!pos] <- pmin(hi[!pos], alpha[!pos])

    stp <- ifelse(is.finite(sc$h) & sc$h < 0, -sc$g / sc$h, 0)
    new <- alpha + stp
    bad <- !is.finite(new) | new <= lo | new >= hi
    new[bad] <- 0.5 * (lo[bad] + hi[bad])

    done <- max(abs(new - alpha)) < tol * (1 + max(abs(alpha)))
    alpha <- new
    if (done) break
  }
  if (!all(is.finite(alpha))) {
    return(NULL)
  }
  alpha
}


## ------------------------------------------------------------
## Helper: construct variance-design matrices
##
## If vars is empty, return an intercept-only design so that
## the same code path nests the homoskedastic special case.
## ------------------------------------------------------------
.make_var_design <- function(data, vars, rows = NULL, int_name = "int") {
  if (length(vars) == 0) {
    n_rows <- if (is.null(rows)) nrow(data) else length(rows)
    out <- matrix(1, nrow = n_rows, ncol = 1)
    colnames(out) <- int_name
    return(out)
  }

  out <- as.matrix(data[, vars, drop = FALSE])

  if (!is.null(rows)) {
    out <- out[rows, , drop = FALSE]
  }

  out
}


## ------------------------------------------------------------
## Helper: numerical symmetry enforcement
##
## Averaging M with t(M) makes the *values* symmetric, but not necessarily
## dimnames(M) -- R's `+` keeps the left operand's dimnames, so if M carries
## inconsistent/spurious dimnames (e.g. list(NULL, c("","",...)), a known
## cbind()/solve() artifact when combining an unnamed vector with an
## unnamed matrix -- see the GTRE Lambda_i computation this guards), the
## result still does too. tmvtnorm::checkSymmetricPositiveDefinite() calls
## isSymmetric() with its default check.attributes=TRUE, which fails on
## attributes alone even when every number is symmetric to machine
## precision, raising a spurious "sigma must be a symmetric matrix" error.
## Since a genuinely-value-asymmetric matrix would already fail this same
## isSymmetric() check regardless of dimnames, it's always safe to drop
## dimnames here when they're the only thing making it look asymmetric --
## no value is ever altered, and a real breakdown still isSymmetric()==
## FALSE afterwards.
## ------------------------------------------------------------
.safe_symmetrize <- function(M) {
  S <- 0.5 * (M + t(M))
  if (!is.null(dimnames(S)) &&
    !isTRUE(isSymmetric(S, tol = sqrt(.Machine$double.eps))) &&
    isTRUE(isSymmetric(S, tol = sqrt(.Machine$double.eps), check.attributes = FALSE))) {
    dimnames(S) <- NULL
  }
  S
}


## ------------------------------------------------------------
## Helper: safe inverse for covariance-type matrices
##
## Strategy:
##   1. symmetrize
##   2. try Cholesky
##   3. progressively ridge if needed
##   4. fall back to solve()
## ------------------------------------------------------------
.safe_inverse <- function(M,
                          base_ridge = 1e-10,
                          ridge_mult = 10,
                          max_tries = 8,
                          name = "matrix") {
  M <- .safe_symmetrize(M)
  p <- nrow(M)
  I_p <- diag(p)

  ## First try: Cholesky path
  for (k in 0:max_tries) {
    ridge <- if (k == 0) 0 else base_ridge * ridge_mult^(k - 1)
    M_try <- if (ridge == 0) M else M + ridge * I_p

    chol_obj <- tryCatch(chol(M_try), error = function(e) NULL)

    if (!is.null(chol_obj)) {
      M_inv <- chol2inv(chol_obj)
      return(list(
        value   = .safe_symmetrize(M_inv),
        ridge   = ridge,
        success = TRUE,
        method  = if (ridge == 0) "chol" else "chol_ridge"
      ))
    }
  }

  ## Second try: direct solve fallback
  for (k in 0:max_tries) {
    ridge <- if (k == 0) 0 else base_ridge * ridge_mult^(k - 1)
    M_try <- if (ridge == 0) M else M + ridge * I_p

    sol <- tryCatch(solve(M_try), error = function(e) NULL)

    if (!is.null(sol)) {
      return(list(
        value   = .safe_symmetrize(sol),
        ridge   = ridge,
        success = TRUE,
        method  = if (ridge == 0) "solve" else "solve_ridge"
      ))
    }
  }

  stop(sprintf("Unable to invert %s even after ridging.", name), call. = FALSE)
}


## ------------------------------------------------------------
## Helper: compute Lambda_i and ARR_i robustly
##
## Lambda_i = [VEE_i^{-1} + A_i' SIG_i^{-1} A_i]^{-1}
## ARR_i    = Lambda_i A_i' SIG_i^{-1}
## ------------------------------------------------------------
.safe_linear_combo <- function(VEE_i, A_i, SIG_i,
                               base_ridge = 1e-10,
                               ridge_mult = 10,
                               max_tries = 8,
                               name = "posterior system") {
  invVEE <- .safe_inverse(VEE_i,
    base_ridge = base_ridge,
    ridge_mult = ridge_mult,
    max_tries  = max_tries,
    name       = "VEE_i"
  )

  invSIG <- .safe_inverse(SIG_i,
    base_ridge = base_ridge,
    ridge_mult = ridge_mult,
    max_tries  = max_tries,
    name       = "SIG_i"
  )

  K <- .safe_symmetrize(invVEE$value + t(A_i) %*% invSIG$value %*% A_i)

  invK <- .safe_inverse(K,
    base_ridge = base_ridge,
    ridge_mult = ridge_mult,
    max_tries  = max_tries,
    name       = name
  )

  ARR <- invK$value %*% t(A_i) %*% invSIG$value

  list(
    LAM = .safe_symmetrize(invK$value),
    ARR = ARR,
    invVEE_method = invVEE$method,
    invSIG_method = invSIG$method,
    invK_method = invK$method,
    ridge_VEE = invVEE$ridge,
    ridge_SIG = invSIG$ridge,
    ridge_K = invK$ridge
  )
}


## Function for extending user formulas
.format_formula <- function(input_val) {
  # Convert formula object to a single string
  input_string <- paste(format(input_val), collapse = "")

  # Split by pipes
  parts <- trimws(strsplit(input_string, "\\|")[[1]])

  # Replace empty parts with "1"
  parts[parts == ""] <- "1"

  # Pad with "1" until there are 3 components
  while (length(parts) < 3) {
    parts <- c(parts, "1")
  }

  # Return as a formula
  return(as.formula(paste(parts, collapse = " | ")))
}


## ------------------------------------------------------------
## Single, robust entry point for parsing sfa's pipe-delimited formula
## syntax (y ~ x1 + x2 | z1 + z2 | zp1).
##
## Returns:
##   n_parts    : number of RHS parts (1, 2, or 3)
##   y_var      : character, name of the response variable
##   formula_x  : two-sided formula for the main regression part (y ~ x...)
##   formula_z  : two-sided formula for the first pipe segment (y ~ z...),
##                NULL if n_parts < 2
##   formula_zp : two-sided formula for the second pipe segment (y ~ zp...),
##                NULL if n_parts < 3
##
## Deliberately does NOT return variable-name vectors for the z/zp parts --
## data_proc() derives those downstream from `colnames(data_conform(...))`
## instead (already robust to interactions/transforms/etc, since it comes
## from an actual model.matrix rather than text parsing), so duplicating
## that here would just be a second source of truth.
## ------------------------------------------------------------
.parse_pipe_formula <- function(formula) {
  f <- Formula::Formula(formula)
  n_parts <- length(f)[2]

  formula_x <- formula(f, lhs = 1, rhs = 1)
  formula_z <- if (n_parts >= 2) formula(f, lhs = 1, rhs = 2) else NULL
  formula_zp <- if (n_parts >= 3) formula(f, lhs = 1, rhs = 3) else NULL

  y_var <- all.vars(formula_x)[1]

  list(
    n_parts    = n_parts,
    y_var      = y_var,
    formula_x  = formula_x,
    formula_z  = formula_z,
    formula_zp = formula_zp
  )
}


## function to handle gtre lower bounds for sigmas
.generate_sfa_bounds <- function(input_formula, prep, inf_sub = -Inf) {
  # 1. Split formula into parts
  f_char <- paste(format(input_formula), collapse = "")
  parts <- trimws(strsplit(f_char, "\\|")[[1]])
  while (length(parts) < 3) parts <- c(parts, "1")
  parts[parts == ""] <- "1"

  # 2. Start with Variance Components (2 params)
  lower_bounds <- rep(.SFA_CONSTANTS$MIN_POSITIVE, 2)

  # 3. Beta Section (prep$n_x_vars params)
  lower_bounds <- c(lower_bounds, rep(inf_sub, prep$n_x_vars))

  # 4. Delta Section
  # If "1", it's 1 param. If variables, it's n_z_eff params.
  if (parts[2] == "1") {
    # lower_bounds <- c(lower_bounds, .SFA_CONSTANTS$MIN_POSITIVE)
    lower_bounds <- c(lower_bounds, inf_sub)
  } else {
    lower_bounds <- c(lower_bounds, rep(inf_sub, prep$n_z_vars))
  }

  # 5. Delta_p Section
  # If "1", it's 1 param. If variables, it's n_zp_eff params.
  if (parts[3] == "1") {
    # lower_bounds <- c(lower_bounds, .SFA_CONSTANTS$MIN_POSITIVE)
    lower_bounds <- c(lower_bounds, inf_sub)
  } else {
    lower_bounds <- c(lower_bounds, rep(inf_sub, prep$n_zp_vars))
  }

  return(lower_bounds)
}


.check_model_formula_pipes <- function(formula, model_name) {
  # 1. Map model names to their maximum ALLOWED RHS parts (pipes + 1)
  # 1 pipe = 2 parts, 2 pipes = 3 parts, 0 pipes = 1 part
  #
  # NHN/NE/GTRE/TRE (the non-Z models) are capped at 1 part -- i.e. no pipes
  # -- by design, to keep a sharp separation between the plain and Z
  # (heteroskedastic-inefficiency) model variants: a user who wants pipe
  # syntax must name the "_Z" model explicitly (e.g. "GTRE_Z") rather than
  # relying on implicit upgrading from a pipe-bearing formula passed to the
  # plain model name. data_proc() (data.processing.R) still contains the
  # underlying auto-upgrade logic (model_name=="NHN" -> "NHN_Z", etc.), but
  # since this function validates model_name before data_proc() ever runs,
  # that path is unreachable from sfm()/psfm() and is effectively dead code
  # kept in place rather than removed.
  max_parts_map <- c(
    # --- 0 Pipes Allowed (1 Part) ---
    "TFE"        = 1,     "TFE_WMLE"   = 1,   "GTRE_SEQ1"     = 1,
    "GTRE_SEQ2"  = 1,     "SSFE"       = 1,   "PL80"          = 1,
    "BC92"       = 1,     "K1990"      = 1,   "K1990modified" = 1,
    "PL80_MVTN"  = 1,
    "NR"         = 1,     "THT"        = 1,   "NTN"           = 1,
    "NG"         = 1,     "ZISF"       = 1,   "NNAK"          = 1,
    "NHN"        = 1,     "NE"         = 1,   "NU"            = 1,
    "NGE"        = 1,     "NLN"        = 1,   "NW"            = 1,
    "tHN"        = 1,
    "TRE"        = 1,     "GTRE"       = 1,   "GTRE_FML"      = 1,

    # --- 1 Pipe Allowed (2 Parts) ---
    "TRE_Z"      = 2,     "FD"         = 2,   "ZISF_Z"        = 2,
    "NHN_Z"      = 2,     "NE_Z"       = 2,

    # --- 2 Pipes Allowed (3 Parts) ---
    "GTRE_Z"     = 3,     "TTNE"       = 3,   "TTHN"          = 3,
    "TTNLS"      = 3
  )

  # 2. Check if the provided model name is valid
  if (!(model_name %in% names(max_parts_map))) {
    stop(paste("Unknown model_name:", model_name), call. = FALSE)
  }

  # 3. Parse the formula and count RHS parts separated by '|'
  f <- Formula::Formula(formula)
  rhs_parts <- length(f)[2] # Index 2 extracts the RHS parts vector length

  # 4. Get the max allowed parts for this specific model
  allowed_parts <- max_parts_map[model_name]

  # 5. Enforce the limit
  if (rhs_parts > allowed_parts) {
    formula_str <- paste(deparse(formula), collapse = " ")
    max_pipes <- allowed_parts - 1

    stop(paste0(
      "Invalid formula structure for model '", model_name, "'!\n",
      "Formula provided: ", formula_str, "\n",
      "The '", model_name, "' model allows a maximum of ", max_pipes, " pipe separator(s) (", allowed_parts, " parts).\n",
      "Found ", (rhs_parts - 1), " pipe separator(s) (", rhs_parts, " parts) instead."
    ), call. = FALSE)
  }

  return(invisible(TRUE))
}


## ------------------------------------------------------------
## Helper: corrected ordinary least squares (COLS)
##
## The moment estimator of Olson, Schmidt and Waldman (1980). OLS is consistent
## for the SLOPES of a composed-error frontier whatever the one-sided
## distribution -- only the intercept is biased, by E[u], because the composed
## error has a non-zero mean. COLS therefore takes the OLS slopes as they are,
## inverts the second and third central moments of the OLS residuals for the
## scale parameters, and shifts the intercept up by the implied E[u].
##
## It is a genuine alternative to ML rather than a starting-value device: no
## optimizer runs, the answer is closed form and deterministic, and it is the
## conventional robustness check against a maximum-likelihood fit that may have
## found a local optimum. It is less efficient than ML when the distributional
## assumption is right.
##
## The moment inversions, all from the third central moment m3 of eps = v - u:
##
##   NHN  u ~ N+(0, su^2)      m3 = su^3 sqrt(2/pi) (1 - 4/pi)
##                             m2 = sv^2 + su^2 (1 - 2/pi),  E[u] = su sqrt(2/pi)
##   NE   u ~ Exp(mean su)     m3 = -2 su^3
##                             m2 = sv^2 + su^2,             E[u] = su
##   NG   u ~ Gamma(mu, su)    uses the cumulants k3 = -2 mu su^3,
##                             k4 = 6 mu su^4, so su = -k4/(3 k3) and
##                             mu = -k3/(2 su^3);            E[u] = mu su
##
## THE WRONG-SKEW CASE IS THE POINT OF FAILURE, and it is reported rather than
## hidden. A production frontier implies m3 < 0. When a sample comes out with
## m3 >= 0 there is no real cube root that keeps su > 0, the estimator has no
## solution, and the honest answer is su = 0 -- no evidence of inefficiency in
## these data -- with sv^2 taking the whole residual variance. Olson, Schmidt
## and Waldman call this "Type I failure"; it is common in small samples and is
## precisely the diagnostic COLS is useful for. Returns $wrong_skew = TRUE so
## the caller can warn.
.cols_fit <- function(Y, X, model_name, intercept_col = 1L) {
  Y <- as.numeric(Y)
  X <- as.matrix(X)
  fit <- stats::lm.fit(X, Y)
  b <- fit$coefficients
  e <- as.numeric(fit$residuals)
  n <- length(e)
  k <- sum(!is.na(b))
  e <- e - mean(e)

  m2 <- mean(e^2)
  m3 <- mean(e^3)
  m4 <- mean(e^4)

  wrong <- !is.finite(m3) || m3 >= 0
  pars <- switch(model_name,
    "NHN" = {
      su <- if (wrong) 0 else (m3 / (sqrt(2 / pi) * (1 - 4 / pi)))^(1 / 3)
      list(
        sigma_v = sqrt(max(m2 - su^2 * (1 - 2 / pi), .Machine$double.eps)),
        sigma_u = su, extra = NULL, eu = su * sqrt(2 / pi)
      )
    },
    "NE" = {
      su <- if (wrong) 0 else (-m3 / 2)^(1 / 3)
      list(
        sigma_v = sqrt(max(m2 - su^2, .Machine$double.eps)),
        sigma_u = su, extra = NULL, eu = su
      )
    },
    "NG" = {
      k4 <- m4 - 3 * m2^2
      bad <- wrong || !is.finite(k4) || k4 <= 0
      su <- if (bad) 0 else -k4 / (3 * m3)
      sh <- if (bad || su <= 0) 1 else -m3 / (2 * su^3)
      list(
        sigma_v = sqrt(max(m2 - sh * su^2, .Machine$double.eps)),
        sigma_u = su, extra = c(mu = sh), eu = sh * su
      )
    },
    stop("COLS is implemented for model_name \"NHN\", \"NE\" and \"NG\" only. ",
      "The moment inversion is distribution-specific and no closed form is ",
      "available for \"", model_name, "\".",
      call. = FALSE
    )
  )

  ## The only coefficient COLS corrects. OLS slopes are already consistent.
  b_cols <- b
  if (!is.na(intercept_col) && intercept_col >= 1 && intercept_col <= length(b)) {
    b_cols[intercept_col] <- b[intercept_col] + pars$eu
  }

  ## Slope standard errors are the ordinary OLS ones and are valid as such.
  ## The intercept's is NOT -- it inherits the sampling error of a third-moment
  ## estimate, which OLS knows nothing about -- and neither sigma has one at
  ## all in closed form here, so both are returned as NA rather than as a
  ## number that would be read as inference. Use cols_boot for those.
  s2 <- sum(e^2) / (n - k)
  R <- qr.R(fit$qr)
  V <- tryCatch(chol2inv(R) * s2, error = function(err) matrix(NA_real_, k, k))
  se_b <- sqrt(pmax(diag(V), 0))
  if (!is.na(intercept_col) && intercept_col >= 1 && intercept_col <= length(se_b)) {
    se_b[intercept_col] <- NA_real_
  }

  list(
    beta = b_cols, se_beta = se_b, sigma_v = pars$sigma_v,
    sigma_u = pars$sigma_u, extra = pars$extra, eu = pars$eu,
    residuals = as.numeric(Y - X %*% b_cols), wrong_skew = wrong,
    moments = c(m2 = m2, m3 = m3, m4 = m4)
  )
}


## ------------------------------------------------------------
## log(exp(a) + exp(b)), elementwise, without leaving the log scale.
##
## Used by zsfm()'s two-component mixture. Forming prob*exp(f1) + (1-prob)*exp(f2)
## and then taking the log underflows both terms to 0 whenever an observation is
## unlikely under BOTH regimes, which is exactly where the optimizer needs
## gradient information. The old code papered over that with log(f + 1e-10),
## which does not just guard the log -- it floors the per-observation
## contribution at -23.03 and makes the objective FLAT across the whole region
## beyond it, so every badly-fitting parameter vector scores the same.
.log_add2 <- function(a, b) {
  m <- pmax(a, b)
  out <- m + log(exp(a - m) + exp(b - m))
  ## m = -Inf means both terms are -Inf; the arithmetic above gives NaN there.
  out[!is.finite(m)] <- m[!is.finite(m)]
  out
}

## ------------------------------------------------------------
## Helper: starting values for the normal-Rayleigh model
##
## sfm("NR") fits eps = v - u with v ~ N(0, sigma_v^2) and u ~ Rayleigh. Its
## likelihood is exact: it integrates to 1 and reproduces a numerical
## normal-Rayleigh convolution to 1.2e-8. What failed was the START. start_cs()
## hard-codes sigma_u = sigma_v = 0.1, and from there NR reached a point WORSE
## than the true parameter vector in 9 of 14 replications at n = 4000, in two
## distinct modes: a hard collapse to sigma_u = 1e-7 with sigma_v inflated to
## absorb the spread (a ~50 log-likelihood deficit), and a partial stall at
## sigma_u ~ 0.85 with the intercept ~0.4 too low (a 6-25 deficit). Seeding at
## the truth fixed all 14, which is what makes this a starting-value problem
## and not an identification one. The moment start below also fixes all 14, and
## lands on the same optimum the truth-seeded run finds.
##
## The moment equations. sfm()'s NR likelihood carries sigma_u on the
## second-raw-moment convention E[u^2] = sigma_u^2 (the same one NHN uses), so
## the Rayleigh scale is theta = sigma_u/sqrt(2) and
##   Var(u) = (1 - pi/4) sigma_u^2,   E[u] = sigma_u sqrt(pi)/2
## The Rayleigh skewness is a CONSTANT, g1 = 2 sqrt(pi) (pi-3)/(4-pi)^(3/2)
## = 0.6311, so the third central moment of eps identifies Var(u) outright:
##   m3 = -g1 Var(u)^(3/2)  =>  Var(u) = (-m3/g1)^(2/3)
## and then sigma_v^2 = m2 - Var(u). That the skewness is fixed rather than a
## free function of the parameters is exactly why NR and NHN are different
## families -- the half-normal's constant is 0.9953 -- so no reparameterization
## can carry one likelihood onto the other.
##
## Returns the full start vector (sigma_v, sigma_u, betas), or NULL when the
## residuals are wrongly skewed and the moment equations have no admissible
## solution; the caller then falls back to the old flat start.
.nr_start <- function(epsilon_hat, beta_0_st, beta_hat) {
  e <- as.numeric(epsilon_hat)
  if (!length(e) || any(!is.finite(e))) {
    return(NULL)
  }
  e <- e - mean(e)
  m2 <- mean(e^2)
  m3 <- mean(e^3)
  if (!is.finite(m2) || !is.finite(m3) || m2 <= 0 || m3 >= 0) {
    return(NULL)
  }

  g1 <- 2 * sqrt(pi) * (pi - 3) / (4 - pi)^1.5 ## 0.6311, Rayleigh skewness
  vu <- (-m3 / g1)^(2 / 3) ## Var(u)
  su <- sqrt(vu / (1 - pi / 4))
  ## A wrong-skew-free sample can still put all of the spread in u. Keep a
  ## floor under sigma_v rather than starting at zero: the likelihood divides
  ## eps by sigma_v, so a zero start is not merely a poor guess, it makes the
  ## first objective evaluation non-finite.
  sv <- sqrt(max(m2 - vu, 0.01 * m2))
  if (!is.finite(su) || !is.finite(sv) || su <= 0) {
    return(NULL)
  }

  ## OLS puts the intercept E[u] below the frontier; shift it back up.
  b0 <- if (is.na(beta_0_st)) NA else unname(beta_0_st) + su * sqrt(pi) / 2
  out <- if (is.na(beta_0_st)) {
    unname(c(sv, su, beta_hat))
  } else {
    unname(c(sv, su, b0, beta_hat))
  }
  if (any(!is.finite(out))) {
    return(NULL)
  }
  out
}

## ------------------------------------------------------------
## Helper: starting values for the normal-gamma model
##
## sfm("NG") fits eps = v - u with v ~ N(0, sigma_v^2) and u ~ Gamma(shape =
## mu, scale = sigma_u). Its likelihood is exact -- checked against numerical
## convolution of the two densities to 7e-14 across eps from +2 to -10, both
## branches of .log_pcf() -- and its profile in mu is smooth and unimodal at
## the truth. The model nevertheless failed every convergence test, because it
## was STARTED in the wrong place: start_cs() hard-codes sigma_u = sigma_v = 0.1
## and start_v_ng pinned the shape at 1, so the search began from E[u] = 0.1
## against a true 1, with the intercept at the raw OLS value (itself E[u] below
## the frontier). From there the optimizer drives sigma_u to its lower bound --
## the "no inefficiency" corner -- inflates sigma_v to absorb the spread, and
## stops 710 log-likelihood units WORSE than the true parameter vector.
##
## The difficulty is the (mu, sigma_u) ridge: the data pin E[u] = mu*sigma_u
## far better than they pin either factor, which is the weak identification
## Ritter and Simar (1997) document for this model. So the candidates below
## walk ALONG that ridge -- E[u] held at a moment estimate while mu varies --
## rather than across it.
##
## Two moment estimators supply the anchor:
##   * gamma cumulants. Gamma(P, th) has k_r = P (r-1)! th^r, and the normal
##     contributes only to k2, so for eps = v - u:
##         k2 = sigma_v^2 + P th^2,  k3 = -2 P th^3,  k4 = 6 P th^4
##     giving th = -k4/(3 k3), P = -k3/(2 th^3), sigma_v^2 = k2 - P th^2.
##     Exact in principle, but k4 is a high-variance statistic: at n = 1000 it
##     returned mu = 4.1 against a true 2. Used as one candidate, not as truth.
##   * the half-normal third-moment estimator, which needs only k3 and is far
##     steadier. It gets E[u] roughly right whatever the true shape is, which
##     is all the anchor has to do.
##
## Returns a LIST of candidate parameter vectors (sigma_v, sigma_u, mu, betas).
## The caller evaluates the likelihood at each and keeps the best.
.ng_start_candidates <- function(epsilon_hat, beta_0_st, beta_hat,
                                 mu_grid = c(0.5, 1, 2, 4, 8)) {
  e <- as.numeric(epsilon_hat)
  e <- e - mean(e)
  n <- length(e)
  if (!n || any(!is.finite(e))) {
    return(list())
  }

  k2 <- mean(e^2)
  k3 <- mean(e^3)
  k4 <- mean(e^4) - 3 * k2^2
  sd_e <- sqrt(max(k2, .Machine$double.eps))

  ## Anchor for E[u]. Half-normal third moment first (steady), gamma cumulants
  ## as a cross-check. Both need negative skew; a wrong-skew sample carries no
  ## inefficiency signal, so fall back to a fraction of the residual spread.
  kk <- sqrt(2 / pi) * (1 - 4 / pi) ## < 0
  eu_hn <- if (is.finite(k3) && k3 < 0) sqrt(2 / pi) * (k3 / kk)^(1 / 3) else NA_real_
  eu <- if (is.finite(eu_hn) && eu_hn > 0) eu_hn else 0.5 * sd_e
  eu <- min(max(eu, 1e-3), 10 * sd_e)

  mk <- function(sv, su, m) {
    if (!is.finite(sv) || !is.finite(su) || !is.finite(m) ||
      sv <= 0 || su <= 0 || m <= 0) {
      return(NULL)
    }
    b0 <- if (is.na(beta_0_st)) NULL else unname(beta_0_st) + m * su
    unname(c(sv, su, m, b0, beta_hat))
  }

  cands <- list()

  ## The gamma-cumulant solution, when the moments admit one.
  if (is.finite(k3) && k3 < 0 && is.finite(k4) && k4 > 0) {
    th <- -k4 / (3 * k3)
    P <- -k3 / (2 * th^3)
    sv <- sqrt(max(k2 - P * th^2, (0.05 * sd_e)^2))
    cands <- c(cands, list(mk(sv, th, P)))
  }

  ## Along the ridge: E[u] fixed at the anchor, shape swept over the grid, and
  ## sigma_v taking whatever variance the gamma component does not explain.
  for (m in mu_grid) {
    su <- eu / m
    sv <- sqrt(max(k2 - m * su^2, (0.05 * sd_e)^2))
    cands <- c(cands, list(mk(sv, su, m)))
  }

  Filter(function(z) !is.null(z) && all(is.finite(z)), cands)
}


## ------------------------------------------------------------
## Helper: the GTRE two-step (moment) decomposition
##
## Second and third central moments of the within residuals and of the
## individual effects, inverted for the normal/half-normal variance split at
## each level. This is what psfm(model_name = "GTRE_SEQ2") reports, and it is
## the "two step" / "pseudo maximum likelihood" estimator of Colombi (2010,
## sec. 3) and Colombi, Martini and Vittadini (2011, sec. 3.2).
##
## Both papers recommend it specifically as the STARTING POINT for the
## closed-skew-normal FIML search -- Colombi (2010) calls it "more stable with
## respect to bad starting values" and notes the two-step estimates "are good
## starting point for the search of the maximum likelihood estimates". That is
## not a stylistic preference: the FIML surface carries a boundary optimum at
## sigma_h = 0, where the model collapses to TRE and the intercept absorbs the
## missing E[h] = sigma_h*sqrt(2/pi). A fit that falls in has a LOWER
## likelihood than the truth, so it is an optimizer failure rather than an
## identification problem, and a better start avoids it.
##
## Returns the same (gamma, sigmaSq) parameterization GTRE_SEQ2 reports:
## gamma is the one-sided share of the variance at that level and sigmaSq the
## total, so sigma_u = sqrt(gamma_uv*sigmaSq_uv), sigma_v = sqrt((1 -
## gamma_uv)*sigmaSq_uv), and likewise sigma_h/sigma_r from the _hr pair.
##
## PRECONDITION: both inputs must be MEAN ZERO. The estimator uses raw second
## and third moments as if they were central ones, which holds for regression
## residuals and for ranef() by construction but not for an arbitrary vector.
## Feeding it uncentred data silently inflates m2 and biases m3, and gamma
## saturates at its min(1, .) cap.
.gtre_two_step <- function(epsilon_hat, alpha_hat, beta_0_st) {
  ## pi - 4 < 0, so k < 0 and k*m3 >= 0 given the min(0, .) truncation below;
  ## the fractional powers stay real.
  k <- sqrt(pi / 2) * (pi / (pi - 4))
  lvl <- function(z) {
    z <- as.numeric(z)
    m2 <- mean(z^2)
    m3 <- min(0, mean(z^3)) ## wrong-skew draws carry no inefficiency signal
    list(
      gamma = min(1, 1 / (m2 * (k * m3)^(-2 / 3) + (2 / pi))),
      sigmaSq = m2 + (2 / pi) * (k * m3)^(2 / 3),
      shift = sqrt(2 / pi) * (k * m3)^(1 / 3)
    )
  }
  e <- lvl(epsilon_hat)
  a <- lvl(alpha_hat)
  list(
    gamma_uv = e$gamma, sigmaSq_uv = e$sigmaSq,
    gamma_hr = a$gamma, sigmaSq_hr = a$sigmaSq,
    beta_0 = beta_0_st + e$shift + a$shift
  )
}


## ------------------------------------------------------------
## Helpers: GTRE by full information maximum likelihood (GTRE_FML)
##
## The four-component model
##   y_it = b0 + x_it'beta + r_i - h_i + v_it - u_it
## with r_i ~ N(0, sigma_r^2), v_it ~ N(0, sigma_v^2) two-sided and
## h_i ~ N+(0, sigma_h^2), u_it ~ N+(0, sigma_u^2) one-sided, has a firm-level
## density in CLOSED FORM as a closed-skew-normal (CSN). Writing
##   A = -[1_T, I_T]   (T x (T+1)),
##   V = diag(sigma_h^2, sigma_u^2, ..., sigma_u^2),
##   Sigma = sigma_v^2 I_T + sigma_r^2 1_T 1_T',
## the T-vector y_i is CSN with
##   scale   Sig    = Sigma + A V A'
##   skew    R      = Lambda A' Sigma^{-1}
##   spread  Lambda = V - V A' Sig^{-1} A V
## and density
##   f(y_i) = c^{-1} * phi_T(y_i; mu_i, Sig) * Phi_{T+1}(R(y_i - mu_i); 0, Lambda).
##
## This is FULL information: unlike psfm()'s simulated-ML "GTRE", nothing is
## integrated out by Monte Carlo, so the likelihood is deterministic. The cost
## is one (T+1)-dimensional normal CDF per firm per evaluation, which is
## expensive and, for T beyond roughly 10, the binding constraint.
##
## The normalizing constant c = Phi_{T+1}(0; 0, Lambda + R Sig R') is exactly
## 2^{-(T+1)}: that covariance comes out DIAGONAL for this A/V/Sigma structure,
## and a zero-mean orthant probability at the origin under independence is
## (1/2)^q. Verified numerically at T = 3, 5 and 8 to eight decimals. Taking it
## as a constant rather than evaluating another CDF per firm is the single
## biggest saving available here.
##
## Ported from `base code/FML.R`, with mnormt (already imported) in place of
## mvtnorm, the per-firm loop reduced to the CDF only, and the multivariate
## normal density evaluated for all firms at once.
## ------------------------------------------------------------

## Build the CSN pieces that do not vary across firms. Returns NULL if the
## parameters put Sigma or Sig outside the positive-definite region.
##
## The key structural fact, and the reason this estimator is usable at all:
## Sig = Sigma + A V A' collapses to (sigma_v^2 + sigma_u^2) I + (sigma_r^2 +
## sigma_h^2) 11', i.e. EQUICORRELATED. Propagating that through
## Lambda = V - V A' Sig^{-1} A V leaves an arrow-plus-exchangeable matrix,
## which is exactly DIAGONAL PLUS RANK ONE:
##
##   Lambda = D + w w',   w = (q/sqrt(e), sqrt(e), ..., sqrt(e)),
##                        D = diag(p - q^2/e, d - e, ..., d - e)
##
## writing p = Lambda[1,1], q = Lambda[1,2], d = Lambda[2,2], e = Lambda[2,3].
## Verified to 1e-16 for T = 3, 5, 8, 10 across a range of parameter values.
##
## That matters because a (T+1)-dimensional normal CDF with a rank-one
## covariance is a ONE-dimensional Gaussian integral -- the same reduction
## .log_mvn_cdf_rank1() already makes for the TFE model. The alternative is
## mnormt::pmnorm()'s Genz algorithm, which is both slow (it dominated the
## whole fit) and STOCHASTIC, so it would leave the log-likelihood noisy and
## the numerical Hessian unreliable. This route is deterministic.
.csn_gtre_parts <- function(sig_r, sig_v, sig_h, sig_u, BigT) {
  if (any(!is.finite(c(sig_r, sig_v, sig_h, sig_u))) ||
    sig_v <= 0 || sig_r < 0 || sig_h <= 0 || sig_u <= 0) {
    return(NULL)
  }

  a <- sig_v^2 + sig_u^2 ## Sig = a I + b 11'
  b <- sig_r^2 + sig_h^2
  if (a <= 0 || a + BigT * b <= 0) {
    return(NULL)
  }

  ## Sig^{-1} = (1/a) I - cc 11'   (Sherman-Morrison on the equicorrelated form)
  cc <- b / (a * (a + BigT * b))

  ## Lambda blocks, from V - (AV)' Sig^{-1} (AV) with AV = -[sig_h^2 1, sig_u^2 I]
  s2h <- sig_h^2
  s2u <- sig_u^2
  p <- s2h - s2h^2 * (BigT / a - cc * BigT^2) ## Lambda[1,1]
  q <- -s2h * s2u * (1 / a - cc * BigT) ## Lambda[1, j+1]
  ## The T-block of Lambda is  (d + e) on the diagonal and e off it, so the
  ## rank-one part contributes e everywhere and what is left on the diagonal is
  ## d itself. Subtracting e a second time here (i.e. writing d - e) drives D
  ## negative at perfectly ordinary parameter values -- at sigma = (0.2, 0.3,
  ## 0.4, 1.0) with T = 5 it gives -0.005 -- which trips the guard below and
  ## makes the likelihood return its penalty exactly where the truth lies.
  d <- s2u - s2u^2 / a ## Lambda[2,2] - e
  e <- s2u^2 * cc ## off-diagonal of the T-block
  if (!is.finite(e) || e <= 0) {
    return(NULL)
  }

  w <- c(q / sqrt(e), rep(sqrt(e), BigT))
  D <- c(p - q^2 / e, rep(d, BigT))
  if (any(!is.finite(D)) || any(D <= 0)) {
    return(NULL)
  }

  ## R = Lambda A' Sigma^{-1}, still needed to form the CDF argument.
  OneT <- matrix(1, BigT, 1)
  EyeT <- diag(BigT)
  A <- -cbind(OneT, EyeT)
  V <- diag(c(s2h, rep(s2u, BigT)))
  Sigma <- sig_v^2 * EyeT + sig_r^2 * tcrossprod(OneT)
  Sig <- a * EyeT + b * tcrossprod(OneT)
  Sigma_i <- tryCatch(solve(Sigma), error = function(err) NULL)
  if (is.null(Sigma_i)) {
    return(NULL)
  }
  ## Guard this the same way as solve(Sigma) directly above. It was a bare
  ## solve(), so an optimizer probing a parameter vector where Sig is singular
  ## got a Lapack EXCEPTION instead of the NULL that every caller here already
  ## treats as "bad parameters, return the penalty". That surfaced as
  ## psfm(model_name = "GTRE_FML") dying mid-fit with
  ## "Lapack routine dgesv: system is exactly singular" under tight optimizer
  ## caps, rather than simply stepping away from the bad point.
  Sig_i <- tryCatch(solve(Sig), error = function(err) NULL)
  if (is.null(Sig_i)) {
    return(NULL)
  }
  Lambda <- V - V %*% t(A) %*% Sig_i %*% A %*% V
  R <- Lambda %*% t(A) %*% Sigma_i

  list(
    Sig = .safe_symmetrize(Sig), R = R, w = w, D = D,
    log_c = -(BigT + 1) * log(2)
  )
}

## log Phi_q(z; 0, D + w w') for every firm at once, by the rank-one reduction:
## with X = sqrt(D) xi + w eta, Pr(X <= z) = E_eta[ prod_j Phi((z_j - w_j eta)/sqrt(D_j)) ].
.log_csn_cdf_rank1 <- function(Z, w, D, gh) {
  sd_j <- sqrt(D)
  nodes <- sqrt(2) * gh$nodes
  lw <- log(gh$weights) - 0.5 * log(pi)
  acc <- matrix(NA_real_, nrow(Z), length(nodes))
  for (r in seq_along(nodes)) {
    ## sweep the node shift across columns, then sum log Phi along each row
    acc[, r] <- lw[r] + rowSums(pnorm(sweep(Z, 2, w * nodes[r], "-") /
      rep(sd_j, each = nrow(Z)), log.p = TRUE))
  }
  m <- apply(acc, 1, max)
  m + log(rowSums(exp(acc - m)))
}

## Negative summed log-likelihood. `gid` is an integer group id in 1..ngroups
## whose rows are already ordered within firm; the panel must be balanced.
.csn_gtre_loglik <- function(par, Y, X, gid, ngroups, BigT, gh) {
  k <- ncol(X)
  b0 <- par[1]
  beta <- par[2:(k + 1)]
  sig_r <- par[k + 2]
  sig_v <- par[k + 3]
  sig_h <- par[k + 4]
  sig_u <- par[k + 5]

  P <- .csn_gtre_parts(sig_r, sig_v, sig_h, sig_u, BigT)
  if (is.null(P)) {
    return(1e12)
  }

  eps <- as.numeric(Y - b0 - X %*% beta)
  E <- matrix(eps[order(gid)], nrow = ngroups, ncol = BigT, byrow = TRUE)

  ## Density term: same covariance for every firm, so evaluate all at once on
  ## the centred residuals rather than once per firm.
  ld <- tryCatch(mnormt::dmnorm(E, mean = rep(0, BigT), varcov = P$Sig, log = TRUE),
    error = function(e) NULL
  )
  if (is.null(ld) || any(!is.finite(ld))) {
    return(1e12)
  }

  ## CDF term. Rank-one reduction, evaluated for all firms at once.
  Z <- E %*% t(P$R) ## ngroups x (T+1)
  lp <- tryCatch(.log_csn_cdf_rank1(Z, P$w, P$D, gh), error = function(e) NULL)
  if (is.null(lp) || any(!is.finite(lp))) {
    return(1e12)
  }

  ll <- sum(ld + lp - P$log_c)
  if (!is.finite(ll)) {
    return(1e12)
  }
  -ll
}


## ------------------------------------------------------------
## Helper: phased starting values for the expensive panel likelihoods
##
## The simulated-ML and CSN panel models spend essentially all their time in
## likelihood evaluations, so the way to make them fast is to need fewer of
## them. Three observations drive this:
##
##   1. The frontier coefficients are nearly free. A random-effects or within
##      regression already puts beta within a few percent of the MLE, and beta
##      is the part of the parameter vector the composed-error structure cares
##      least about. Optimizing over it in the early stages just multiplies the
##      dimension of a very expensive search.
##   2. The variance/efficiency parameters are where the likelihood actually
##      has curvature, where the bad local optima live, and where a poor start
##      costs the most iterations.
##   3. Evaluating the likelihood on a coarse GRID over those few parameters is
##      cheap relative to letting a derivative-free optimizer discover the same
##      region, because the grid needs no trial steps and cannot wander.
##
## So: hold beta fixed, grid-search the efficiency block, polish it with a
## low-dimensional nlminb, and only then hand the full vector to the ordinary
## optimizer stack. The full maximization still happens -- this only supplies
## it with a better starting point.
##
## MEASURED RESULT: this is OFF by default, because it costs more time than it
## saves. Standard panel DGP, N = 200, T = 10, three seeds, median seconds and
## rmse-to-truth:
##
##   model      baseline      phased, beta_0 FIXED   phased, beta_0 FREE
##   GTRE       46.1s 0.157     55.6s 0.172            58.7s 0.122
##   GTRE_FML    5.8s 0.043      2.8s 0.535             8.8s 0.038
##   TRE         8.3s 0.121     12.0s 0.090            14.2s 0.120
##
## THE INTERCEPT MUST BE FREE. The middle column is what happens when it is not:
## GTRE_FML came back twice as fast and twelve times less accurate. In a
## composed-error model the intercept and the mean of the one-sided component
## are confounded -- E[y - x'beta] = beta_0 - E[u] -- so every candidate value
## of the variance block implies a DIFFERENT intercept. Comparing candidates at
## a shared, stale beta_0 evaluates each at the wrong location, and the phased
## point then lands somewhere the downstream optimizer reads as converged, so it
## stops early at a non-optimum. Profiling beta_0 at every grid point (int_idx)
## removes the pathology completely: 0.535 -> 0.038.
##
## But even done correctly it is not a SPEEDUP for these models -- every one is
## slower than simply using better starting values. They are not start-limited:
## the staged optimizer already reaches the same optimum, so the grid and the
## profiling are added cost rather than saved cost. GTRE does reach its best
## accuracy this way (0.122 against 0.130), at +21% time, which is the only
## case where the trade might be worth taking.
##
## Kept and gated behind options(sfa.phased_start = TRUE). The idea is sound
## for likelihoods that genuinely ARE start-limited; these are not.
##
##   fn       the FULL negative log-likelihood, over the complete parameter vector
##   start_v  current full starting vector
##   idx      indices of the efficiency/variance block to search over
##   lower    lower bounds for those indices (recycled if length 1)
##   grid     multiplicative factors applied to start_v[idx], expanded over the
##            block; the identity factor 1 is always included so the incoming
##            start is never discarded
##   max_pts  cap on grid points evaluated, so cost stays bounded as the block grows
## ------------------------------------------------------------
.phased_start <- function(fn, start_v, idx, lower = .SFA_CONSTANTS$MIN_POSITIVE,
                          grid = c(0.5, 1, 2), int_idx = NULL, max_pts = 81L,
                          polish = TRUE, maxit = 60L, verbose = FALSE) {
  if (!length(idx)) {
    return(start_v)
  }
  base <- start_v

  ## Profile the intercept at a candidate point. The intercept and the mean of
  ## the one-sided component are confounded -- E[y - x'beta] = beta_0 - E[u] --
  ## so ANY change to the variance block implies a different intercept. Holding
  ## it fixed evaluates every off-centre grid point at the wrong location, which
  ## makes good candidates look bad and is why the first version of this helper
  ## made things worse rather than better. The likelihood is smooth and
  ## single-peaked in the intercept, so a 1-D optimize() is reliable and cheap.
  prof <- function(cand) {
    if (is.null(int_idx)) {
      return(list(par = cand, val = tryCatch(fn(cand), error = function(e) Inf)))
    }
    b0 <- cand[int_idx]
    w <- max(abs(b0), 1) * 2
    o <- tryCatch(
      stats::optimize(
        function(v) {
          cc <- cand
          cc[int_idx] <- v
          vv <- tryCatch(fn(cc), error = function(e) Inf)
          if (is.finite(vv)) vv else 1e12
        },
        interval = c(b0 - w, b0 + w), tol = 1e-4
      ),
      error = function(e) NULL
    )
    if (is.null(o)) {
      return(list(par = cand, val = tryCatch(fn(cand), error = function(e) Inf)))
    }
    cand[int_idx] <- o$minimum
    list(par = cand, val = o$objective)
  }

  p0 <- prof(base)
  best <- p0$par
  fbest <- if (is.finite(p0$val)) p0$val else Inf

  ## --- stage 1: COORDINATE-WISE grid over the efficiency block -------------
  ## One parameter varied at a time rather than a full factorial. With the
  ## intercept profiled at every point the per-point cost is ~10 evaluations,
  ## so a 3^4 factorial would cost more than the optimizer it is meant to save.
  ## Coordinate-wise finds the same scale corrections at a fraction of that.
  g <- setdiff(unique(grid), 1)
  for (j in idx) {
    for (fac in g) {
      cand <- best
      cand[j] <- max(best[j] * fac, lower)
      pj <- prof(cand)
      if (is.finite(pj$val) && pj$val < fbest) {
        fbest <- pj$val
        best <- pj$par
      }
    }
  }
  if (verbose) cat(sprintf("  phased: grid -> %.4f\n", fbest))

  ## --- stage 2: optimize the efficiency block AND the intercept together ---
  if (polish) {
    blk <- c(idx, int_idx)
    sub_fn <- function(p) {
      cand <- best
      cand[blk] <- p
      v <- tryCatch(fn(cand), error = function(e) Inf)
      if (is.finite(v)) v else 1e12
    }
    lo <- c(rep_len(lower, length(idx)), rep(-Inf, length(int_idx)))
    op <- tryCatch(
      stats::nlminb(
        start = best[blk], objective = sub_fn, lower = lo,
        control = list(iter.max = maxit, eval.max = 2L * maxit)
      ),
      error = function(e) NULL
    )
    if (!is.null(op) && is.finite(op$objective) && op$objective < fbest) {
      best[blk] <- op$par
      fbest <- op$objective
      if (verbose) cat(sprintf("  phased: block polish -> %.4f\n", fbest))
    }
  }
  best
}


## ------------------------------------------------------------
## Helper: per-firm simulated log-density for psfm()'s GTRE/TRE likelihood,
## evaluated for EVERY firm at once.
##
## The original likelihood looped lapply(1:N, ...) and called dnorm()/pnorm()
## separately on each firm's Ti x R matrix, so a single objective evaluation
## made N pairs of vectorized calls on small matrices. R's per-call overhead
## then dominates: for N = 200, T = 10 that is 400 calls on 10 x R blocks
## instead of 2 calls on one (N*T) x R block.
##
## Two further changes come with it:
##
##   * The per-firm product over t is taken in LOGS (rowsum of log terms)
##     rather than by multiplying Ti densities and flooring at double.xmin.
##     Same quantity, but it cannot underflow for long panels, and it removes
##     the pmax() floor entirely -- after the existing z-clipping the smallest
##     attainable pnorm() is ~5.7e-300, comfortably above double.xmin, so the
##     floor never bound in practice anyway.
##   * The average over draws is a log-sum-exp rather than mean-then-log.
##
## `E` is the (n x R) matrix of composed errors, rows in firm order, `gid` the
## firm index (non-decreasing), and the return value is the length-N vector of
## per-firm log simulated densities.
## ------------------------------------------------------------
.gtre_sim_logdens <- function(E, lambda, sigma, gid, ngroups) {
  z1 <- E / sigma
  z2 <- -E * lambda / sigma
  z1[z1 > .SFA_CONSTANTS$CLIP_Z1_UPPER] <- .SFA_CONSTANTS$CLIP_Z1_UPPER
  z1[z1 < .SFA_CONSTANTS$CLIP_Z1_LOWER] <- .SFA_CONSTANTS$CLIP_Z1_LOWER
  z2[z2 > .SFA_CONSTANTS$CLIP_Z2_UPPER] <- .SFA_CONSTANTS$CLIP_Z2_UPPER
  z2[z2 < .SFA_CONSTANTS$CLIP_Z2_LOWER] <- .SFA_CONSTANTS$CLIP_Z2_LOWER

  lt <- log(2 / sigma) + dnorm(z1, log = TRUE) + pnorm(z2, log.p = TRUE)
  lp <- rowsum(lt, gid, reorder = FALSE) ## N x R, log product over t
  m <- do.call(pmax, as.data.frame(lp)) ## row maxima, no apply()
  m + log(rowMeans(exp(lp - m)))
}


## ------------------------------------------------------------
## Helper: log parabolic cylinder function log D_nu(z), for nu < 0
##
## Used by sfm()'s NG (normal-gamma) and NNAK (normal-Nakagami) likelihoods and
## their efficiency predictors, which are the only models in the package that
## need it. Both call it with nu = -mu or -2*mu for a positive shape mu, so
## a = -nu > 0 throughout.
##
## THE PROBLEM THIS REPLACES. Both models previously used the series form
##
##   D_nu(z) = 2^(nu/2) e^(-z^2/4) sqrt(pi) [ 1F1(-nu/2; 1/2; z^2/2)/Gamma((1-nu)/2)
##                                  - sqrt(2) z 1F1((1-nu)/2; 3/2; z^2/2)/Gamma(-nu/2) ]
##
## which is exact but, for z > 0, computes an O(exp(-z^2/4)) answer as the
## difference of two O(exp(z^2/2)) terms. The cancellation is total by z ~ 8.
## Checked against numerical integration of the defining integral:
##
##   z        series form      truth
##   6        -12.659922      -12.659941
##   8               NaN      -20.203426
##  10         -6.613448      -29.634184     <- wrong by 23 log units
##  12               NaN      -40.990162
##
## Since z = eps/sigma_v + sigma_v/sigma_u for NG, any observation more than a
## few residual standard deviations out lands in that region, so the likelihood
## was returning NaN or nonsense on a subset of the sample at most parameter
## values. That is the leading explanation for NG and NNAK failing their
## convergence tests while every other closed-form model passed, and for NG's
## shape parameter never moving off its starting value.
##
## THE FIX. For z > 0 use the confluent hypergeometric function of the SECOND
## kind, which IS the decaying solution, so nothing cancels:
##
##   D_{-a}(z) = 2^(-a/2) exp(-z^2/4) U(a/2, 1/2, z^2/2)
##
## verified against numerical integration to all printed digits for z in
## [1, 12]. For z <= 0 the two series terms ADD rather than cancel (-sqrt(2) z
## is positive there), so the original form is well conditioned and is kept;
## its argument is clipped only to stop 1F1 overflowing for very negative z.
## ------------------------------------------------------------
.log_pcf <- function(nu, z) {
  z <- as.numeric(z)
  a <- -nu
  out <- rep(NA_real_, length(z))
  ok <- is.finite(z)
  if (!any(ok)) {
    return(out)
  }

  hi <- ok & z > 0.5
  lo <- ok & !(z > 0.5)

  if (any(hi)) {
    u <- tryCatch(gsl::hyperg_U(a / 2, 0.5, z[hi]^2 / 2), error = function(e) NULL)
    out[hi] <- if (is.null(u)) {
      NA_real_
    } else {
      (nu / 2) * log(2) - z[hi]^2 / 4 + log(pmax(u, .Machine$double.xmin))
    }
  }
  if (any(lo)) {
    zz <- z[lo]
    ## Clip the 1F1 argument: exp(z^2/2) overflows past z ~ 37, and the series
    ## branch only ever sees z <= 0.5 where the result is finite anyway.
    q <- pmin(zz^2 / 2, .SFA_CONSTANTS$EXP_CLIP_UPPER)
    br <- gsl::hyperg_1F1(-nu / 2, 0.5, q) / gamma((1 - nu) / 2) -
      sqrt(2) * zz * gsl::hyperg_1F1((1 - nu) / 2, 1.5, q) / gamma(-nu / 2)
    out[lo] <- (nu / 2) * log(2) + 0.5 * log(pi) - zz^2 / 4 +
      log(pmax(br, .Machine$double.xmin))
  }
  out
}


## ---------------------------------------------------------------------------
## Case-insensitive matching for model / estimator names.
##
## match.arg() is case SENSITIVE, so psfm(model_name = "gtre") errored with a
## list of valid names that visibly contained what the user had just typed.
## None of the five entry points has a case collision among its choices --
## sfm()'s "THT" and "tHN" differ in more than case -- so folding case is
## unambiguous and the canonical spelling is always recoverable.
##
## Exact (case-insensitive) matches win over partial ones, which matters here:
## "GTRE" is a prefix of "GTRE_Z", "GTRE_FML", "GTRE_SEQ1" and "GTRE_SEQ2", and
## must resolve to itself rather than being called ambiguous.
## ---------------------------------------------------------------------------
.match_model_name <- function(x, choices, arg = "model_name") {
  if (missing(x) || is.null(x) || !length(x)) {
    return(choices[1])
  }
  ## The unevaluated default is the whole choice vector; match.arg() treats
  ## that as "take the first", and so must this.
  if (length(x) > 1L) {
    if (identical(as.character(x), as.character(choices))) {
      return(choices[1])
    }
    stop("'", arg, "' must be a single value, not ", length(x), ".", call. = FALSE)
  }
  x <- as.character(x)
  if (is.na(x)) {
    stop("'", arg, "' must not be NA.", call. = FALSE)
  }

  lo <- tolower(choices)
  hit <- match(tolower(x), lo)
  if (!is.na(hit)) {
    return(choices[hit])
  }
  ## Preserve match.arg()'s partial matching, folded to lower case.
  p <- pmatch(tolower(x), lo)
  if (!is.na(p)) {
    return(choices[p])
  }

  ## Rank suggestions by edit distance and offer at most two. A wider net
  ## (agrep at max.distance 0.35) returned seven candidates for one typo,
  ## which is noise rather than help.
  dd <- utils::adist(tolower(x), lo)[1, ]
  near <- choices[order(dd)][seq_len(min(2L, length(choices)))]
  near <- near[sort(dd)[seq_along(near)] <= max(2, ceiling(nchar(x) / 2))]
  stop("'", arg, "' = \"", x, "\" is not a recognized choice.",
    if (length(near)) paste0(" Did you mean ", paste0("\"", near, "\"", collapse = " or "), "?") else "",
    "\n  Valid choices (case does not matter): ", paste(choices, collapse = ", "),
    call. = FALSE
  )
}
