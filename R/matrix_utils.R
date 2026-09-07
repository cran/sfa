## Helper: column-wise products of a matrix (replaces matrixStats::colProds)
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


## Helper: validate the arguments common to every exported model fitter
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


## Helpers: leave the caller's random number stream as we found it
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


## Helper: physicists' Gauss-Hermite quadrature nodes/weights
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


## Helper: Gauss-Legendre quadrature nodes/weights on (0, 1)
## Gauss-Legendre nodes and weights on (0, 1), MEMOISED on n.
##
## The rule depends on nothing but n, and building it is an eigendecomposition
## of an n x n matrix -- which numDeriv, differencing a quadrature-based moment
## function, would otherwise repeat hundreds of times per call. Caching also
## keeps eigen() out of forked workers: on macOS the threaded BLAS is not
## fork-safe and mclapply() segfaulted inside eigen() before this was added.
.GL_CACHE <- new.env(parent = emptyenv())

.gauss_legendre_01 <- function(n = 64L) {
  n <- as.integer(n)
  if (n < 5L) stop("n must be at least 5.", call. = FALSE)
  key <- as.character(n)
  hit <- .GL_CACHE[[key]]
  if (!is.null(hit)) {
    return(hit)
  }
  i <- seq_len(n - 1L)
  J <- matrix(0, n, n)
  off <- i / sqrt(4 * i^2 - 1)
  J[cbind(i, i + 1L)] <- off
  J[cbind(i + 1L, i)] <- off
  eg <- eigen(J, symmetric = TRUE)
  ord <- order(eg$values)
  out <- list(
    nodes = (eg$values[ord] + 1) / 2,
    weights = eg$vectors[1L, ord]^2
  )
  assign(key, out, envir = .GL_CACHE)
  out
}


## Helper: Student-t--half-normal (tHN) composed-error density
.thn_m_for <- function(lambda) {
  ## The optimizer WILL visit sigma_v -> 0, so lambda arrives as Inf or NaN.
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
## predictor.
.thn_parts <- function(e, sigma_v, sigma_u, nu) {
  nd <- .thn_nodes(.thn_m_for(sigma_u / sigma_v))
  s <- nd$s
  ## The t density is written out rather than calling dt() on the whole
  ## observations-by-nodes matrix.
  x <- outer(as.numeric(e), sigma_u * s, "+") / sigma_v
  kc <- exp(lgamma((nu + 1) / 2) - lgamma(nu / 2)) / (sqrt(nu * pi) * sigma_v)
  list(
    fv = kc * exp(-((nu + 1) / 2) * log1p(x * x / nu)),
    wu = nd$w * sqrt(2 / pi) * exp(-s^2 / 2),
    s = s
  )
}

## Largest lambda = sigma_u/sigma_v the quadrature is validated at.
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


## Helper: log MVN CDF for the rank-one/equicorrelated covariance Sigma = I_T
## + c*11' (replaces mnormt::pmnorm() for this special case)
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


## Helper: log density of the within-transformed disturbance vector for
## psfm()'s TFE model (replaces mnormt::dmnorm() for this special case)
.log_within_mvn_density <- function(eps_star, sigma2) {
  m <- length(eps_star)
  q <- sum(eps_star^2) + sum(eps_star)^2
  -0.5 * m * log(2 * pi) - 0.5 * m * log(sigma2) + 0.5 * log(m + 1) - q / (2 * sigma2)
}


## Helper: cross-sectional intercept-only normal-half-normal (Aigner, Lovell &
## Schmidt 1977) MLE, with Hessian-based standard errors.
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


## Helper: time-decay pattern B_it for psfm()'s error-components-frontier
## models (PL80/BC92/K1990/K1990modified), y_it = x_it'beta + v_it - B_it*u_i.
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


## Helper: pre-estimation collinearity diagnostic for panel models.
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
      "warn_drop_partial" = paste0(
        "  collinear_action = \"warn_drop\" could not remove these: they are SOME\n",
        "  (not all) of the columns of a single factor term, and a formula can only\n",
        "  drop the term whole, which would discard identified columns too.\n",
        "  They were therefore kept in the likelihood and dropped from the\n",
        "  starting-value regression only. To remove them from the model, replace\n",
        "  the factor with explicit dummies for the periods you want, or use\n",
        "  broader period groups / a time trend."
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
## produced them.
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
## requested regressors.
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


## Helper: coerce user data into the panel form psfm() needs.
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


## Helper: analytic gradient of the normal/half-normal negative summed
## log-likelihood, in sfm()'s own (lambda, sigma, beta) parameterization.
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


## Helper: conditional inefficiency predictors.
.log_pnorm_diff <- function(a, b) {
  use_upper <- (a + b) > 0
  la <- ifelse(use_upper, pnorm(-b, log.p = TRUE), pnorm(a, log.p = TRUE))
  lb <- ifelse(use_upper, pnorm(-a, log.p = TRUE), pnorm(b, log.p = TRUE))
  lb + log(-expm1(pmin(la - lb, -.Machine$double.eps)))
}

## Draws for simulated maximum likelihood.
.sml_draws <- function(n_units, n_draws, dim = 1L,
                       sim_type = c("halton", "sobol", "torus", "uniform"),
                       antithetics = FALSE,
                       burn = .SFA_CONSTANTS$HALTON_DISCARD,
                       scrambling = 0L, prime = NULL, seed = NULL,
                       clamp = 1e-6) {
  sim_type <- match.arg(sim_type)
  ## Checked before coercion: as.integer(NULL) is integer(0) and
  ## as.integer(NA) is NA.
  for (nm in c("n_units", "n_draws", "dim")) {
    v <- get(nm)
    if (is.null(v) || length(v) != 1L || !is.numeric(v) || !is.finite(v)) {
      stop("`", nm, "` must be a single positive number.", call. = FALSE)
    }
  }
  n_units <- as.integer(n_units)
  n_draws <- as.integer(n_draws)
  dim <- as.integer(dim)
  if (n_units < 1L || n_draws < 1L || dim < 1L) {
    stop("`n_units`, `n_draws` and `dim` must all be at least 1.", call. = FALSE)
  }

  ## Antithetics halve the number of INDEPENDENT draws needed: take n_base and
  ## mirror them.
  n_base <- if (antithetics) as.integer(ceiling(n_draws / 2)) else n_draws
  n_need <- n_units * n_base

  if (!is.null(seed)) {
    .st <- .rng_snapshot()
    on.exit(.rng_restore(.st), add = TRUE)
    set.seed(seed)
  }

  raw <- switch(sim_type,
    "halton"  = randtoolbox::halton(n_need + burn, dim = dim, start = 1, normal = FALSE),
    "sobol"   = randtoolbox::sobol(n_need + burn, dim = dim, scrambling = scrambling,
                                   seed = if (is.null(seed)) 4711L else as.integer(seed),
                                   normal = FALSE, init = TRUE),
    "torus"   = if (is.null(prime)) {
                  randtoolbox::torus(n_need + burn, dim = dim, normal = FALSE)
                } else {
                  randtoolbox::torus(n_need + burn, dim = dim, prime = prime, normal = FALSE)
                },
    "uniform" = matrix(stats::runif((n_need + burn) * dim), ncol = dim)
  )
  raw <- matrix(as.numeric(raw), ncol = dim)
  if (burn > 0L && nrow(raw) > burn) raw <- raw[-seq_len(burn), , drop = FALSE]
  raw <- raw[seq_len(n_need), , drop = FALSE]

  ## RANDOMIZE BY SHIFTING, for the deterministic sequences.
  ##
  ## halton() and torus() ignore set.seed() entirely -- they are fixed lattices
  ## -- so before this, passing `seed` with the default sim_type did precisely
  ## nothing while looking as though it controlled the draws. Worse, it left
  ## the cross-sectional simulated-ML models with UNrandomized draws, so the
  ## asymptotic theory for simulation-based estimators (which assumes random
  ## draws) was formally unjustified for them, even though the panel path had
  ## been randomized since 1.2.0 by .gtre_halton_draws().
  ##
  ## A shift mod 1 rather than a permutation: it preserves the lattice
  ## structure the sequence exists for, which is Train's (9.3.4) recommended
  ## randomization and exactly what the panel path already does. With seed
  ## NULL nothing shifts, so the default draws are unchanged.
  ## No snapshot/restore here: the `seed` block above already took one and
  ## registered its on.exit. A second pair would run AFTER that one and restore
  ## the post-set.seed(seed) state instead of the caller's, quietly leaking the
  ## seed into the caller's stream -- which is what a first attempt at this did.
  if (!is.null(seed) && sim_type %in% c("halton", "torus")) {
    set.seed(seed)
    sh <- stats::runif(dim)
    for (d in seq_len(dim)) raw[, d] <- (raw[, d] + sh[d]) %% 1
  }

  lapply(seq_len(dim), function(d) {
    ## byrow = TRUE is essential, not cosmetic: filling column-major hands
    ## unit i the stride-n subsequence of a van der Corput sequence.
    M <- matrix(raw[, d], nrow = n_units, ncol = n_base, byrow = TRUE)
    if (antithetics) {
      M <- cbind(M, 1 - M)[, seq_len(n_draws), drop = FALSE]
    }
    ## Clamped away from 0 and 1: the callers apply qnorm() or -log1p(-u), and
    ## an exact endpoint sends those to +/-Inf and kills the optimizer.
    pmin(pmax(M, clamp), 1 - clamp)
  })
}


.te_battese_coelli <- function(mu_star, sigma_star) {
  sigma_star <- pmax(sigma_star, .SFA_CONSTANTS$MIN_POSITIVE)
  z <- mu_star / sigma_star
  out <- exp(-mu_star + 0.5 * sigma_star^2 +
    pnorm(z - sigma_star, log.p = TRUE) - pnorm(z, log.p = TRUE))
  pmin(pmax(out, 0), 1)
}


## G2/I3: seed a hard model from a simpler one already fitted.
##
## Matching is BY PARAMETER NAME, not by position. Position would be actively
## dangerous here: `NHN` reports (lambda, sigma) where `NE` reports
## (sigv, sigu), so a positional copy would slide a variance ratio into a
## standard deviation and produce a confidently wrong start. By name, the two
## documented idioms both work -- NE -> NG shares sigv/sigu, NHN -> NTN shares
## lambda/sigma -- and a pair that shares nothing but the betas carries only
## the betas, which is the honest answer rather than a silent misalignment.
##
## Anything not matched keeps the start the package computed for it.
.start_from <- function(start_v, out, from, model_name, quiet = FALSE) {
  if (!inherits(from, "sfareg")) {
    stop("`start_from` must be a fitted \"sfareg\" object.", call. = FALSE)
  }
  tgt <- colnames(out)
  src <- from$out
  if (is.null(tgt) || length(tgt) != length(start_v) ||
      is.null(src) || is.null(rownames(src)) || !("par" %in% colnames(src))) {
    warning("start_from: could not read parameter names from one of the two ",
      "fits, so the supplied fit was ignored and the default starting values ",
      "were used.", call. = FALSE)
    return(start_v)
  }

  vals <- stats::setNames(as.numeric(src[, "par"]), rownames(src))
  vals <- vals[is.finite(vals)]
  hit <- intersect(tgt, names(vals))

  if (!length(hit)) {
    warning("start_from: the fitted ", from$model_name, " model shares no ",
      "parameter NAMES with ", model_name, " (", paste(tgt, collapse = ", "),
      "), so nothing could be carried over and the default starting values ",
      "were used.", call. = FALSE)
    return(start_v)
  }
  start_v[match(hit, tgt)] <- vals[hit]
  if (!quiet) {
    message("start_from: carried ", length(hit), " of ", length(tgt),
      " starting values from the ", from$model_name, " fit (",
      paste(hit, collapse = ", "), ").")
  }
  start_v
}

.jlms_u <- function(mu_star, sigma_star) {
  sigma_star <- pmax(sigma_star, .SFA_CONSTANTS$MIN_POSITIVE)
  z <- mu_star / sigma_star
  pmax(mu_star + sigma_star * exp(dnorm(z, log = TRUE) - pnorm(z, log.p = TRUE)), 0)
}


## Helper: the composed density of normal noise and a one-sided inefficiency,
## by multiple importance sampling. Half the draws come from the noise and half
## from the inefficiency, each weighted by the balance heuristic, so the
## estimate stays accurate whichever of the two densities is the narrower.
.sml_mis <- function(eps, sigma_v, FiMat, ldens, qdens) {
  m <- ncol(FiMat) %/% 2L
  e_v <- as.numeric(eps)
  n <- length(e_v)
  ## Kept in logs: p_hi underflows for a very efficient unit.
  lp_hi <- pnorm(e_v / sigma_v, lower.tail = FALSE, log.p = TRUE)
  Tm <- qnorm(lp_hi + log1p(-FiMat[, seq_len(m), drop = FALSE]),
    lower.tail = FALSE, log.p = TRUE
  )
  u <- cbind(sigma_v * Tm - e_v, qdens(FiMat[, m + seq_len(m), drop = FALSE]))
  ## A non-finite draw is given zero weight rather than aborting, so the
  ## efficiency predictor can always be formed; `ok` lets the likelihood steer
  ## the optimizer away from the region instead.
  ok_u <- all(is.finite(u))

  ## Balance-heuristic weight for equal counts from the two proposals. With
  ## a = log of the noise kernel and b = log p_hi + log f_u, the algebra
  ## collapses to a + b - logsumexp(a, b) -- one exp and one log1p per draw.
  zz <- (e_v + u) / sigma_v
  a <- -0.5 * zz * zz - .SFA_CONSTANTS$LOG_SQRT_2PI - log(sigma_v)
  b <- lp_hi + ldens(u)
  w_log <- pmin(a, b) - log1p(exp(-abs(a - b)))
  w_log[!is.finite(w_log)] <- -Inf

  ## Row log-sum-exp gives the density; the same weights give any posterior
  ## mean, so likelihood and efficiency predictor cannot drift apart.
  mx <- w_log[cbind(seq_len(n), max.col(w_log, ties.method = "first"))]
  ok <- is.finite(mx)
  w <- exp(w_log - ifelse(ok, mx, 0))
  w[!is.finite(w)] <- 0
  ldens_hat <- rep(-Inf, n)
  if (any(ok)) {
    ldens_hat[ok] <- mx[ok] + log(rowSums(w[ok, , drop = FALSE])) - log(m)
  }
  list(u = u, w = w, ldens = ldens_hat, ok = ok_u)
}

## Helper: a posterior mean over the draws .sml_mis() already weighted
.sml_mis_mean <- function(mis, gu) {
  den <- pmax(rowSums(mis$w), .Machine$double.xmin)
  rowSums(mis$w * gu) / den
}

## Helper: Horrace and Schmidt (1996) intervals for u_i
.horrace_schmidt_ci <- function(mu_star, sigma_star, level = 0.95) {
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
    level <= 0 || level >= 1) {
    stop("`level` must be a single number strictly between 0 and 1.", call. = FALSE)
  }
  sigma_star <- pmax(sigma_star, .SFA_CONSTANTS$MIN_POSITIVE)
  n <- max(length(mu_star), length(sigma_star))
  mu_star <- rep_len(mu_star, n)
  sigma_star <- rep_len(sigma_star, n)

  alpha <- 1 - level
  log_A <- pnorm(mu_star / sigma_star, log.p = TRUE)

  lower <- mu_star + sigma_star *
    qnorm(log(1 - alpha / 2) + log_A, lower.tail = FALSE, log.p = TRUE)
  upper <- mu_star + sigma_star *
    qnorm(log(alpha / 2) + log_A, lower.tail = FALSE, log.p = TRUE)

  ## u is non-negative by construction; only rounding can push the lower bound
  ## below zero, and the ordering can only invert the same way.
  lower <- pmax(lower, 0)
  upper <- pmax(upper, lower)
  list(lower = lower, upper = upper)
}


## Helpers: Greene (2005) true fixed effects stochastic frontier

## Group sums by an integer group id in 1..ngroups, independent of the order
## in which the groups happen to appear in the data.
.gsum <- function(x, gid, ngroups) {
  s <- rowsum(as.numeric(x), gid, reorder = FALSE)
  out <- numeric(ngroups)
  out[as.integer(rownames(s))] <- s[, 1L]
  out
}

## Concentrated firm effects for the TFE profile likelihood.
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


## Helper: construct variance-design matrices
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


## Helper: numerical symmetry enforcement
.safe_symmetrize <- function(M) {
  S <- 0.5 * (M + t(M))
  if (!is.null(dimnames(S)) &&
    !isTRUE(isSymmetric(S, tol = sqrt(.Machine$double.eps))) &&
    isTRUE(isSymmetric(S, tol = sqrt(.Machine$double.eps), check.attributes = FALSE))) {
    dimnames(S) <- NULL
  }
  S
}


## Helper: safe inverse for covariance-type matrices
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


## Helper: compute Lambda_i and ARR_i robustly
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


## Single, robust entry point for parsing sfa's pipe-delimited formula syntax
## (y ~ x1 + x2 | z1 + z2 | zp1).
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
  ## 1. Map model names to their maximum ALLOWED RHS parts (pipes + 1) 1 pipe
  ## = 2 parts, 2 pipes = 3 parts, 0 pipes = 1 part
  max_parts_map <- c(
    # --- 0 Pipes Allowed (1 Part) ---
    "TFE"        = 1,     "TFE_WMLE"   = 1,   "GTRE_SEQ1"     = 1,
    "GTRE_SEQ2"  = 1,     "SSFE"       = 1,   "PL80"          = 1,
    "BC92"       = 1,     "K1990"      = 1,   "K1990modified" = 1,
    "PL80_MVTN"  = 1,
    "NR"         = 1,     "THT"        = 1,   "NTN"           = 1,
    "NG"         = 1,     "ZISF"       = 1,   "NNAK"          = 1,
    "LCM"        = 1,    "LCM_CN"     = 1,
    "NHN"        = 1,     "NE"         = 1,   "NU"            = 1,
    "NGE"        = 1,     "NLN"        = 1,   "NW"            = 1,
    "tHN"        = 1,     "TSL"        = 1,
    "TRE"        = 1,     "GTRE"       = 1,   "GTRE_FML"      = 1,

    "CSS"        = 1,     "LS"         = 1,
    "SSRE"       = 1,     "SSCRE"      = 1,   "KSS"           = 1,

    # --- 1 Pipe Allowed (2 Parts) ---
    "TRE_Z"      = 2,     "FD"         = 2,   "ZISF_Z"        = 2,
    "NHN_Z"      = 2,     "NE_Z"       = 2,   "LCM_Z"         = 2,

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


## Helper: corrected ordinary least squares (COLS)
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


## log Phi(z) + z^2/2, without ever forming either term.
##
## Every exponentially tilted Gaussian in this package (NE, NGE, and the
## two-tier likelihoods) needs log Phi(z) added to a tilt that equals z^2/2 up
## to terms in eps^2/sigma_v^2.  Computed separately the two diverge with
## opposite signs as the one-sided scale goes to zero, and their sum is a
## catastrophic cancellation: at z = -5.9e9 the terms are -/+1.7e19, where
## consecutive doubles are 2048 apart, so the sum returns rounding noise --
## positive noise, which an optimiser will happily maximise.
##
## For z -> -Inf,
##   log Phi(z) = -z^2/2 - log(2 pi)/2 - log(-z) + log1p(-1/z^2 + 3/z^4 - ...)
## so the z^2/2 cancels analytically and nothing large is ever formed.  The
## direct form is exact for moderate z and is kept there; the switch at -1e3
## is conservative (the naive sum only starts losing digits near |z| ~ 1e6).
.log_phi_tilt <- function(z, z_switch = -1e3) {
  out <- numeric(length(z))
  ok <- z > z_switch
  if (any(ok)) {
    out[ok] <- stats::pnorm(z[ok], log.p = TRUE) + z[ok]^2 / 2
  }
  if (any(!ok)) {
    zz <- z[!ok]
    z2 <- zz^2
    out[!ok] <- -0.5 * log(2 * pi) - log(-zz) +
      log1p(-1 / z2 + 3 / z2^2 - 15 / z2^3)
  }
  out
}


## log(exp(a) + exp(b)), elementwise, without leaving the log scale.
.log_add2 <- function(a, b) {
  m <- pmax(a, b)
  out <- m + log(exp(a - m) + exp(b - m))
  ## m = -Inf means both terms are -Inf; the arithmetic above gives NaN there.
  out[!is.finite(m)] <- m[!is.finite(m)]
  out
}

## log(sum_j exp(M[i, j])) for each row of a matrix -- the J-class
## generalisation of .log_add2, used to mix the latent-class components and to
## normalise the multinomial-logit class probabilities. Kept on the log scale
## throughout: a class contribution can easily underflow to zero in absolute
## terms while still carrying the posterior weight for that observation.
.log_row_sum_exp <- function(M) {
  M <- as.matrix(M)
  mx <- do.call(pmax, c(lapply(seq_len(ncol(M)), function(j) M[, j]), na.rm = TRUE))
  ## An all -Inf row is a genuine zero density; the shift below would make it
  ## NaN, so carry the -Inf through untouched.
  bad <- !is.finite(mx)
  mx[bad] <- 0
  out <- mx + log(rowSums(exp(M - mx)))
  out[bad] <- -Inf
  out
}

## Helper: starting values for the normal-Rayleigh model
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
  ## A wrong-skew-free sample can still put all of the spread in u.
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

## Helper: starting values for the normal-gamma model
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

  ## Anchor for E[u].
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


## Helper: the GTRE two-step (moment) decomposition
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


## Helpers: GTRE by full information maximum likelihood (GTRE_FML)

## Build the CSN pieces that do not vary across firms. Returns NULL if the
## parameters put Sigma or Sig outside the positive-definite region.
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
  ## The T-block of Lambda is (d + e) on the diagonal and e off it.
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
  ## Guard this the same way as solve(Sigma) directly above.
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


## Helper: phased starting values for the expensive panel likelihoods
.phased_start <- function(fn, start_v, idx, lower = .SFA_CONSTANTS$MIN_POSITIVE,
                          grid = c(0.5, 1, 2), int_idx = NULL, max_pts = 81L,
                          polish = TRUE, maxit = 60L, verbose = FALSE) {
  if (!length(idx)) {
    return(start_v)
  }
  base <- start_v

  ## Profile the intercept at a candidate point.
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
  ## One parameter varied at a time rather than a full factorial.
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


## Helper: per-firm simulated log-density for psfm()'s GTRE/TRE likelihood,
## evaluated for EVERY firm at once.
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



## log D_nu(z) for nu < 0 by the integral representation
##
##   D_nu(z) = e^{-z^2/4} / Gamma(-nu) * int_0^inf t^{-nu-1} e^{-t^2/2 - z t} dt
##
## which holds for every nu < 0 and needs no special functions at all. Used
## where gsl::hyperg_U throws or underflows, i.e. exactly where the shape
## parameter is large.
##
## The integrand is evaluated in LOGS with its own maximum factored out, so
## nothing overflows however large the shape gets: g(t) = (a-1)log t - t^2/2 -
## z t peaks at the positive root of t^2 + z t - (a-1) = 0, and integrating
## exp(g(t) - g(t*)) keeps every evaluated term at or below 1.
.log_pcf_integral <- function(nu, z) {
  z <- as.numeric(z)
  nu <- rep_len(as.numeric(nu), length(z))
  a_v <- -nu
  vapply(seq_along(z), function(i) {
    zz <- z[i]
    a <- a_v[i]
    if (!is.finite(zz) || !is.finite(a) || a <= 0) return(NA_real_)
    g <- function(t) (a - 1) * log(t) - t^2 / 2 - zz * t
    ## g'(t) = (a-1)/t - t - z. For a > 1 that has one positive root and the
    ## integrand has an interior peak worth factoring out. For a <= 1 it has
    ## none -- g is decreasing on (0, Inf) and diverges as t -> 0 -- so there is
    ## no finite maximum to factor, and none is needed: a small shape keeps the
    ## integrand at moderate values and the plain integral is safe.
    hi_t <- 0
    if (a > 1) {
      tstar <- (-zz + sqrt(zz^2 + 4 * (a - 1))) / 2
      if (!is.finite(tstar) || tstar <= 0) tstar <- .Machine$double.eps
      gmax <- g(tstar)
      ## The integrand is log-concave, so a generous window either side of the
      ## peak captures it to well past double precision.
      w <- 40 * max(sqrt(1 / (a + zz^2 + 1)), 1) + 12
      lo_t <- max(tstar - w, .Machine$double.eps)
      hi_t <- tstar + w
    } else {
      gmax <- 0
      lo_t <- 0
      hi_t <- 40 + 8 * max(zz, 0)
    }
    val <- tryCatch(
      stats::integrate(function(t) exp(g(t) - gmax), lo_t, hi_t,
        rel.tol = 1e-10, stop.on.error = FALSE
      )$value,
      error = function(e) NA_real_
    )
    if (!is.finite(val) || val <= 0) return(NA_real_)
    -zz^2 / 4 - lgamma(a) + gmax + log(val)
  }, numeric(1))
}

## Helper: log parabolic cylinder function log D_nu(z), for nu < 0
.log_pcf <- function(nu, z) {
  z <- as.numeric(z)
  ## nu may be a VECTOR, one order per observation: that is what lets the
  ## gamma/Nakagami SHAPE depend on covariates (gap G5). A scalar recycles, so
  ## every existing caller is unaffected.
  nu <- rep_len(as.numeric(nu), length(z))
  a <- -nu
  out <- rep(NA_real_, length(z))
  ok <- is.finite(z) & is.finite(nu) & nu < 0
  if (!any(ok)) {
    return(out)
  }

  hi <- ok & z > 0.5
  lo <- ok & !(z > 0.5)

  if (any(hi)) {
    ## gsl::hyperg_U THROWS for a/2 >~ 16 at small arguments, and the call is
    ## vectorized, so a single bad element used to null the whole vector and
    ## take the entire fit down with it -- which is why NNAK "failed outright"
    ## once its shape wandered above about 8. Two changes: fall back per
    ## ELEMENT rather than all-or-nothing, and have something to fall back TO.
    u <- tryCatch(gsl::hyperg_U(a[hi] / 2, 0.5, z[hi]^2 / 2), error = function(e) NULL)
    v <- if (is.null(u)) rep(NA_real_, sum(hi)) else as.numeric(u)
    val <- (nu[hi] / 2) * log(2) - z[hi]^2 / 4 + log(pmax(v, .Machine$double.xmin))
    ## Underflow of U to exactly 0 is a silent precision failure, not an error:
    ## log(xmin) is finite and wrong. Treat it as a miss too.
    bad <- !is.finite(val) | !is.finite(v) | v <= 0
    if (any(bad)) val[bad] <- .log_pcf_integral(nu[hi][bad], z[hi][bad])
    out[hi] <- val
  }
  if (any(lo)) {
    zz <- z[lo]
    nl <- nu[lo]
    ## Clip the 1F1 argument: exp(z^2/2) overflows past z ~ 37, and the series
    ## branch only ever sees z <= 0.5 where the result is finite anyway.
    q <- pmin(zz^2 / 2, .SFA_CONSTANTS$EXP_CLIP_UPPER)
    br <- gsl::hyperg_1F1(-nl / 2, 0.5, q) / gamma((1 - nl) / 2) -
      sqrt(2) * zz * gsl::hyperg_1F1((1 - nl) / 2, 1.5, q) / gamma(-nl / 2)
    ## NOTE: the floor below is deliberately LEFT IN PLACE. The series is a
    ## difference and can go non-positive where the two terms nearly cancel,
    ## and log(pmax(br, xmin)) is then finite and wrong -- so replacing it with
    ## the integral representation looks like a strict improvement. It is not.
    ## `NG` evaluates on this branch for 99.8% of observations at its own
    ## optimum, and the floor acts as a barrier keeping the optimizer out of
    ## the sigma_v -> 0 corner where the composed likelihood is unbounded.
    ## Removing it was measured: NG moved from (sigv 0.185, x1 0.560) to
    ## (sigv 0.000, x1 0.775) against a true x1 of 0.5, reaching a HIGHER
    ## likelihood (-438.02 vs -440.71) at a degenerate point. The barrier is
    ## accidental but load-bearing; the reported NNAK failure was on the `hi`
    ## branch and is fixed there. Revisit only alongside a proper boundary
    ## penalty or constraint on sigma_v.
    out[lo] <- (nl / 2) * log(2) + 0.5 * log(pi) - zz^2 / 4 +
      log(pmax(br, .Machine$double.xmin))
  }
  out
}


## Case-insensitive matching for model / estimator names.
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

  ## Rank suggestions by edit distance and offer at most two.
  dd <- utils::adist(tolower(x), lo)[1, ]
  near <- choices[order(dd)][seq_len(min(2L, length(choices)))]
  near <- near[sort(dd)[seq_along(near)] <= max(2, ceiling(nchar(x) / 2))]
  stop("'", arg, "' = \"", x, "\" is not a recognized choice.",
    if (length(near)) paste0(" Did you mean ", paste0("\"", near, "\"", collapse = " or "), "?") else "",
    "\n  Valid choices (case does not matter): ", paste(choices, collapse = ", "),
    call. = FALSE
  )
}

## Helpers: the inefficiency density and quantile for the simulated-ML models
.nsml_ldens <- function(model_name, sigma_u, shape) {
  if (model_name == "NLN") {
    function(u) dlnorm(u, meanlog = shape, sdlog = sigma_u, log = TRUE)
  } else {
    function(u) dweibull(u, shape = shape, scale = sigma_u, log = TRUE)
  }
}

.nsml_qdens <- function(model_name, sigma_u, shape) {
  if (model_name == "NLN") {
    function(p) qlnorm(p, meanlog = shape, sdlog = sigma_u)
  } else {
    function(p) qweibull(p, shape = shape, scale = sigma_u)
  }
}

## Helper: is a fitted one-sided scale sitting on the zero boundary because the
## residuals are skewed the wrong way?
##
## sigma_u = 0 under wrong skew is the CORRECT maximum likelihood estimate, not
## a numerical failure -- Waldman (1982) showed the likelihood has a stationary
## point there, and Olson, Schmidt and Waldman (1980) call it the Type I
## failure. It is common in small samples: at lambda = 0.75 with N = 100, 13.2%
## of samples come out wrongly skewed and 17.7% of THOSE put sigma_u on the
## boundary, while not one correctly skewed sample does. So the right response
## is to say so, not to bound sigma_u away from zero -- a bound would corrupt
## the very samples where the boundary is the answer.
##
## `scale_hat` is the fitted one-sided scale (or lambda), `ref` the residual
## scale it is judged against.
## `rel_tol` is how small the one-sided scale has to be, relative to the
## residual SD, to count as collapsed. 1e-2 rather than the 1e-3 this used at
## first: on the NE design the test uses, the collapsed sample sits at
## 8.6e-4 and the five non-collapsed ones at 0.47-0.71. 1e-3 therefore cleared
## the boundary case by a factor of 1.16 while clearing the others by 470 --
## knife-edge on the side that matters, and the optimizer's stopping point
## moves by more than 16% between BLAS implementations. It passed here and
## failed R CMD check on CI's macOS. 1e-2 sits in the middle of a gap spanning
## nearly three orders of magnitude: 12x clear of the boundary case, 47x clear
## of the nearest interior one.
##
## Loosening this is low-risk because it is only ever the SECOND condition.
## The warning needs wrong skew (m3 >= 0) as well, and that is a deterministic
## function of the data.
.wrong_skew_boundary <- function(resid, scale_hat, ref, model_name,
                                 rel_tol = 1e-2) {
  e <- as.numeric(resid)
  e <- e[is.finite(e)]
  if (length(e) < 4L || !is.finite(scale_hat) || !is.finite(ref) || ref <= 0) {
    return(list(wrong_skew = NA, at_bound = NA, m3 = NA_real_))
  }
  m3 <- mean((e - mean(e))^3)
  list(wrong_skew = m3 >= 0,
       at_bound = scale_hat / ref < rel_tol,
       m3 = m3)
}

## Helper: the warning text for the case above
.warn_wrong_skew_boundary <- function(ws, model_name, scale_name = "sigma_u") {
  warning("sfm(model_name = \"", model_name, "\"): ", scale_name,
    " has collapsed to the boundary and the OLS residuals are skewed the WRONG ",
    "way (third central moment ", signif(ws$m3, 3), " >= 0). A production ",
    "frontier implies negative skew, so there is no interior maximum and the ",
    "boundary IS the maximum likelihood estimate -- this is the Type I failure ",
    "of Olson, Schmidt and Waldman (1980), not a numerical problem. Read it as ",
    "no evidence of inefficiency in these data rather than as an estimate of ",
    "zero, and treat the efficiency scores as uninformative. Inspect ",
    "$wrong_skew and $sigma_u_at_bound.",
    call. = FALSE
  )
}

## Helper: the variance-determinant block of a psfm() _Z fit, in the shape
## marginal_effects() expects. psfm() places z'delta on the VARIANCE, so the
## link is "var"; sfm()'s default is "sd". The delta coefficients are the
## trailing n_z entries of the parameter vector, named "(Intercept u)" and then
## the z variables.
.psfm_z_spec <- function(data, z_vars, par, family = "halfnormal",
                        anchor = "(Intercept u)") {
  if (!length(z_vars)) {
    return(NULL)
  }
  ## data_proc() returns z_vars WITH "(Intercept)" in it, which is not a column
  ## of `data`; build the intercept column rather than looking it up.
  dd <- data.frame(data, check.names = FALSE)
  has_int <- "(Intercept)" %in% z_vars
  zc <- setdiff(z_vars, "(Intercept)")
  if (!all(zc %in% names(dd))) {
    return(NULL)
  }
  Z <- as.matrix(dd[, zc, drop = FALSE])
  if (has_int) {
    Z <- cbind(`(Intercept)` = 1, Z)
  }
  if (!ncol(Z)) {
    return(NULL)
  }
  k <- ncol(Z)
  ## Locate the block BY NAME, not by position. TRE_Z ends with the sigma_u
  ## block, but GTRE_Z carries two -- "(Intercept u)", z..., "(Intercept h)",
  ## zp... -- so taking the trailing k entries would silently return the
  ## sigma_h coefficients for GTRE_Z.
  nm <- names(par)
  start <- if (!is.null(nm)) match(anchor, nm) else NA_integer_
  if (is.na(start) || start + k - 1L > length(par)) {
    return(NULL)
  }
  delta <- unname(par[seq.int(start, start + k - 1L)])
  if (!all(is.finite(delta))) {
    return(NULL)
  }
  list(Z = Z, delta = delta, link = "var", family = family)
}

## Halton draws for the panel simulated-ML models, one INDEPENDENT block per
## firm. See notes/code_history/matrix_utils.md.
##
## Two departures from what psfm() did before, both fixing defects recorded as
## gaps J1 and J2.
##
## J1: every firm used to be integrated against the SAME R draws, because a
## single R x 2 block was built once and recycled. Train (2002, p. 228) is
## explicit that Halton's advantage comes from its coverage *and the negative
## correlation it induces across observations* -- and that second half only
## exists when observations get different blocks of the sequence. Sharing one
## block also means each firm's simulation error is the same realization
## rather than an independent one, so the errors cannot average out as N
## grows. sfm() already blocks by observation and says why; this makes psfm()
## consistent with it.
##
## J2: the old code drew 9999 random permutations of dimension 1 and kept
## whichever correlated least with dimension 2. A 2-d Halton sequence is worth
## having because the PAIRS cover the unit square evenly; permuting one column
## preserves each margin, destroys the pairing, and leaves joint coverage no
## better than random. Measured at R = 150 with primes 2 and 3 and 1000
## discarded -- psfm()'s exact configuration -- joint discrepancy went from
## 0.0217 to 0.0550 against purely random pairing's 0.1400, to remove a
## correlation of 0.008 that was never a problem. The randomization here is
## instead a uniform shift modulo 1 (Tuffin 1996; Train 9.3.4), which moves the
## lattice without disturbing its structure: same measurement, 0.0250.
.gtre_halton_draws <- function(N, R, rand.gtre = NULL,
                               burn = .SFA_CONSTANTS$HALTON_DISCARD,
                               clamp = 1e-6) {
  n_pt <- as.integer(N) * as.integer(R)
  U <- randtoolbox::halton(n_pt + burn, 2, start = 1, normal = FALSE)
  U <- U[-seq_len(burn), 1:2, drop = FALSE]

  ## Randomize by SHIFTING, not by permuting: a shift preserves the lattice
  ## structure that the sequence exists for.
  if (!is.null(rand.gtre)) {
    .rng_state <- .rng_snapshot()
    on.exit(.rng_restore(.rng_state), add = TRUE)
    set.seed(rand.gtre)
    sh <- stats::runif(2L)
    U[, 1] <- (U[, 1] + sh[1]) %% 1
    U[, 2] <- (U[, 2] + sh[2]) %% 1
  }
  ## qnorm() of an exact endpoint is infinite, and a shift can land on one.
  U[] <- pmin(pmax(U, clamp), 1 - clamp)

  Z <- cbind(stats::qnorm(U[, 1]), sqrt(2) * pracma::erfinv(U[, 2]))
  ## Firm ii takes rows ((ii-1)R+1):(ii R) -- a contiguous block, as sfm() does.
  lapply(seq_len(N), function(ii) Z[((ii - 1L) * R + 1L):(ii * R), , drop = FALSE])
}

## A persistent scale sitting on the zero boundary in a panel fit.
##
## The panel analogue of the cross-sectional wrong-skew case. There, negative
## skew is absent and sigma_u = 0 is the CORRECT maximum likelihood estimate
## (Waldman 1982), which sfm() reports rather than suppresses. Here, the
## likelihood can prefer to merge the two persistent components -- setting
## sigma_h to zero and letting sigma_r absorb the persistent variation -- and
## that too can be the genuine maximum.
##
## Measured on a GTRE replication whose realized sigma_h scale was 0.3695, so
## persistent inefficiency was certainly present: the boundary solution beat
## the TRUE parameter vector by 3.66 log units, a likelihood ratio of about 39.
## It is not an optimizer failure and must not be reported as one.
##
## The tell is that sigma_r comes back inflated whenever this happens: 0.30-0.34
## against a true 0.20 in every collapsed replication observed.
.panel_scale_at_bound <- function(scale_hat, ref, rel_tol = 1e-2) {
  if (!is.finite(scale_hat) || !is.finite(ref) || ref <= 0) {
    return(NA)
  }
  scale_hat / ref < rel_tol
}

## Set $sigh_at_bound / $sigr_at_bound on a fitted GTRE object and warn if
## either persistent scale has collapsed.
##
## Factored out because psfm() reaches this point by two different routes with
## two different parameter layouts: estimator = "sml" reports the
## lambda/sigma reparameterization with sigh and sigr alongside, while
## estimator = "fiml" -- the DEFAULT -- reports the four raw scales. The check
## used to live only in the "sml" branch, so the default estimator carried no
## boundary flag at all. `ref` is the transient scale sqrt(su^2 + sv^2) under
## both layouts, which is what makes the 1% rule mean the same thing in each.
.report_panel_boundary <- function(results, model_name, sig_h, sig_r, ref) {
  results$sigh_at_bound <- .panel_scale_at_bound(sig_h, ref)
  results$sigr_at_bound <- .panel_scale_at_bound(sig_r, ref)
  if (isTRUE(results$sigh_at_bound)) {
    .warn_panel_boundary(model_name, "sigh", "sigr", sig_r)
  } else if (isTRUE(results$sigr_at_bound)) {
    .warn_panel_boundary(model_name, "sigr", "sigh", sig_h)
  }
  results
}

.warn_panel_boundary <- function(model_name, scale_name, absorber = NULL,
                                 absorber_value = NA_real_) {
  msg <- paste0(
    "psfm(model_name = \"", model_name, "\"): ", scale_name, " has collapsed ",
    "to the zero boundary. For a panel with two persistent components this is ",
    "often the CORRECT maximum likelihood estimate rather than a numerical ",
    "failure -- the likelihood can genuinely prefer to merge them -- so it is ",
    "reported rather than suppressed, in the same spirit as sfm()'s ",
    "wrong-skew boundary (Waldman 1982 is the cross-sectional analogue)."
  )
  if (!is.null(absorber) && is.finite(absorber_value)) {
    msg <- paste0(msg, " Note that ", absorber, " = ", signif(absorber_value, 4),
      ", which typically comes back inflated in this situation: the persistent ",
      "variation is not lost but reclassified. Read the two persistent scales ",
      "together rather than separately, and treat the split between them as ",
      "unidentified in this sample."
    )
  }
  ## The intercept is not a bystander here. E[y] = beta0 - E[h] - E[u], so
  ## moving mass between sigma_h and sigma_r moves E[h] = sigma_h*sqrt(2/pi)
  ## and the fitted intercept absorbs the difference -- measured at -0.30
  ## against a predicted -0.319 when sigma_h collapses from a true 0.40. A user
  ## reading a level off coef() needs to be told that, not just the scales.
  msg <- paste0(msg, " The frontier intercept has shifted with it: it carries ",
    "E[h] = sigh*sqrt(2/pi), so a collapsed persistent scale displaces the ",
    "estimated intercept by roughly that amount. Treat levels from coef() ",
    "with the same caution as the split itself."
  )
  warning(msg, call. = FALSE)
}


## Build the reduced random-effects starting-value regression at COLUMN
## granularity. Dropping whole formula TERMS cannot express "keep 13 of a
## factor's 24 dummies", which is exactly what factor(year) in an unbalanced
## panel needs, so the reduced design is assembled from the model matrix and
## handed to plm() under placeholder names. Returns NULL if it cannot be built,
## in which case the caller keeps the unreduced formula.
.re_start_design <- function(formula_x, data, drop_cols) {
  df <- as.data.frame(data)
  idx <- tryCatch(as.data.frame(plm::index(data)), error = function(e) NULL)
  if (is.null(idx) || nrow(idx) != nrow(df) || ncol(idx) < 2L) {
    return(NULL)
  }
  yname <- all.vars(formula_x)[1]
  if (!yname %in% names(df)) {
    return(NULL)
  }
  vars_needed <- intersect(all.vars(formula_x), names(df))
  ok <- stats::complete.cases(df[, vars_needed, drop = FALSE])
  if (!any(ok)) {
    return(NULL)
  }
  X <- tryCatch(stats::model.matrix(formula_x, data = df[ok, , drop = FALSE]),
    error = function(e) NULL
  )
  if (is.null(X) || !ncol(X) || nrow(X) != sum(ok)) {
    return(NULL)
  }
  keep <- setdiff(colnames(X), c("(Intercept)", drop_cols))
  if (!length(keep)) {
    return(NULL)
  }
  safe <- paste0(".v", seq_along(keep))
  nd <- data.frame(
    .id = idx[ok, 1], .time = idx[ok, 2], .y = df[[yname]][ok],
    X[, keep, drop = FALSE], check.names = FALSE
  )
  names(nd) <- c(".id", ".time", ".y", safe)
  pd <- tryCatch(plm::pdata.frame(nd, index = c(".id", ".time")),
    error = function(e) NULL
  )
  if (is.null(pd)) {
    return(NULL)
  }
  list(
    data = pd,
    formula = stats::reformulate(safe, response = ".y"),
    map = stats::setNames(keep, safe)
  )
}

## Rename a plm coefficient/SE vector fitted on .re_start_design()'s
## placeholder columns back to the original model-matrix column names.
.remap_start_names <- function(v, map) {
  if (is.null(v) || !length(v)) {
    return(v)
  }
  nm <- names(v)
  hit <- !is.na(match(nm, names(map)))
  nm[hit] <- unname(map[nm[hit]])
  stats::setNames(v, nm)
}

## Align a named SE vector to the requested coefficient names, leaving NA for
## columns the starting-value regression could not identify. The point
## estimates for those come from a pooled OLS fit (.expand_start_beta), which
## carries no comparable standard error, so NA is the honest entry.
.expand_start_se <- function(se_named, target) {
  out <- stats::setNames(rep(NA_real_, length(target)), target)
  if (!is.null(se_named) && length(se_named) && !is.null(names(se_named))) {
    common <- intersect(names(se_named), target)
    out[common] <- se_named[common]
  }
  out
}


## Delta-method standard errors for the GTRE two-step (moment) decomposition.
## Extracted from psfm()'s GTRE_SEQ2 branch so the algebra is unit-testable.
## `n_eps`/`n_alp` are the counts the two moment sets are averaged over.
##
## Note the maintained assumption: the moment variances below are the iid ones,
## so they describe sampling variability around the estimator's PROBABILITY
## LIMIT, and for fixed T that plim is not the truth (see the GTRE_SEQ note in
## ?psfm). They are also optimistic for the epsilon side, whose n residuals are
## correlated within firm.
.gtre_two_step_se <- function(eps_hat, alp_hat, n_eps, n_alp, beta_se1 = NA_real_) {
  k <- sqrt(pi / 2) * (pi / (pi - 4))
  na5 <- stats::setNames(rep(NA_real_, 5), c(
    "gamma_uv", "sigmaSq_uv", "gamma_hr", "sigmaSq_hr", "beta_0"
  ))

  side <- function(z, nn) {
    z <- as.numeric(z)
    m2 <- mean(z^2)
    m3 <- min(0, mean(z^3))
    ## m3 == 0 is the wrong-skew boundary: the inversion is not differentiable
    ## there and every derivative below diverges, so report NA rather than Inf.
    if (!is.finite(m2) || !is.finite(m3) || m3 >= 0 || nn <= 0) {
      return(NULL)
    }
    g <- min(1, 1 / (m2 * (k * m3)^(-2 / 3) + (2 / pi)))
    ssq <- m2 + (2 / pi) * (k * m3)^(2 / 3)

    ## Central moments of the implied normal-half-normal composite, shared with
    ## the COLS standard errors (.cols_se_nhn) so there is one derivation.
    .m <- .nhn_central_moments(g, ssq)
    mu2 <- .m[["mu2"]]; mu3 <- .m[["mu3"]]; mu4 <- .m[["mu4"]]
    mu5 <- .m[["mu5"]]; mu6 <- .m[["mu6"]]

    var_m2 <- (mu4 - mu2^2) / nn
    var_m3 <- (mu6 - mu3^2 - 6 * mu2 * mu4 + 9 * mu2^3) / nn
    cov_23 <- (mu5 - 4 * mu2 * mu3) / nn

    ## gamma = 1/(m2 A + 2/pi), A = (k m3)^(-2/3); sigmaSq = m2 + (2/pi)(k m3)^(2/3)
    A <- (k * m3)^(-2 / 3)
    gg <- 1 / (m2 * A + (2 / pi))
    d_g_m2 <- -A * gg^2
    d_g_m3 <- (2 / 3) * m2 * k * (k * m3)^(-5 / 3) * gg^2
    d_s_m2 <- 1
    d_s_m3 <- sqrt(2 / pi) * (pi / (pi - 4)) * (2 / 3) * (k * m3)^(-1 / 3)
    d_b_m3 <- (pi / (pi - 4)) * (1 / 3) * (k * m3)^(-2 / 3)

    dm <- function(a, b) a^2 * var_m2 + b^2 * var_m3 + 2 * a * b * cov_23
    list(
      gamma = sqrt(max(0, dm(d_g_m2, d_g_m3))),
      sigmaSq = sqrt(max(0, dm(d_s_m2, d_s_m3))),
      var_b = d_b_m3^2 * var_m3
    )
  }

  e <- side(eps_hat, n_eps)
  a <- side(alp_hat, n_alp)
  if (is.null(e) || is.null(a)) {
    return(na5)
  }
  vb <- if (is.finite(beta_se1)) beta_se1^2 else NA_real_
  c(
    gamma_uv = e$gamma, sigmaSq_uv = e$sigmaSq,
    gamma_hr = a$gamma, sigmaSq_hr = a$sigmaSq,
    beta_0 = sqrt(vb + e$var_b + a$var_b)
  )
}


## GTRE_SEQ1/GTRE_SEQ2 are LARGE-T estimators. Both take plm's random-effects
## decomposition and treat the two generated series as if they were draws from
## the latent composite errors, which they are not at finite T:
##
##   alpha_hat is the shrunken BLUP  c * (alpha_i + ebar_i),
##       c = T s2a / (T s2a + s2e)
##   eps_hat   is the quasi-demeaned residual  (1-th) alpha_i + e_it - th ebar_i,
##       th = 1 - s2e^0.5 / (T s2a + s2e)^0.5
##
## so the second stage inverts the moments of the PREDICTORS. Measured against
## the latent draws this attenuates the third central moment by
## (1 - th/T)^3 + (T-1)(-th/T)^3 on the epsilon side (0.81 at T = 10) and
## scales the alpha side by c^2 and c^3. Neither factor goes to 1 as N grows
## with T held fixed, so the estimators converge to the wrong constants.
.gtre_seq_finiteT_warning <- function(N, n, model_name) {
  Tbar <- if (is.finite(N) && N > 0) n / N else NA_real_
  warning(
    model_name, " is a LARGE-T estimator and is NOT consistent for fixed T.\n",
    "  Its two stages treat plm's random-effects output as if it were the\n",
    "  latent composite errors, but alpha_hat is a SHRUNKEN predictor and\n",
    "  eps_hat is a QUASI-DEMEANED residual. Their second and third moments\n",
    "  are attenuated by factors that depend on T and do not vanish as N grows,\n",
    "  so the variance components are biased toward zero",
    if (is.finite(Tbar)) paste0(" (mean T = ", format(Tbar, digits = 3), ")") else "",
    ".\n",
    "  sigmaSq_hr is affected most. Increasing N tightens these estimates\n",
    "  around the WRONG values; only larger T removes the bias.\n",
    "  For a consistent alternative use model_name = \"GTRE\" (estimator =\n",
    "  \"fiml\" on a balanced panel, \"sml\" otherwise), which estimates all\n",
    "  four variance components jointly by maximum likelihood.",
    call. = FALSE
  )
}


## Central moments 2..6 of the normal-half-normal composite error, in the
## (gamma, sigmaSq) parameterization. Factored out of .gtre_two_step_se() so
## the COLS standard errors below use the SAME expressions rather than a second
## hand-derivation of the same algebra.
.nhn_central_moments <- function(g, ssq) {
  mu2 <- ssq * ((1 - g) + g * ((pi - 2) / pi))
  mu3 <- ssq^(3 / 2) * (sqrt(2 / pi) * (1 - (4 / pi)) * g^(3 / 2))
  mu4 <- ssq^2 * (3 * (1 - g)^2 + ((6 * (pi - 2) * g * (1 - g)) / pi) +
    g^2 * (3 - (4 / pi) - (12 / pi^2)))
  mu5 <- ssq^(5 / 2) * g^(3 / 2) * sqrt(2 / pi) *
    (10 * (1 - (4 / pi)) * (1 - g) + (7 - (20 / pi) - (16 / pi^2)) * g)
  mu6 <- ssq^3 * (15 * (1 - g)^3 + (45 * (pi - 2) * (1 - g)^2 * g / pi) +
    15 * (3 - (4 / pi) - (12 / pi^2)) * (1 - g) * g^2 +
    (15 - (6 / pi) - (100 / pi^2) - (40 / pi^3)) * g^3)
  c(mu2 = mu2, mu3 = mu3, mu4 = mu4, mu5 = mu5, mu6 = mu6)
}

## Coelli (1995, Appendix 1) standard errors for COLS, expressed in the
## (sigma_v, sigma_u) parameterization sfm(estimator = "cols") reports rather
## than the (gamma, sigma^2) one the paper writes them in.
##
## Most applications of COLS quote the OLS standard errors for the variance
## parameters. Those are not appropriate: sigma_u and sigma_v are non-linear
## functions of the residual moments m2 and m3, not regression coefficients.
## The delta method on (m2, m3) is, with all of these divided by n,
##
##   V(m2) = mu4 - mu2^2,  V(m3) = mu6 - mu3^2 - 6 mu2 mu4 + 9 mu2^3,
##   C(m2, m3) = mu5 - 4 mu2 mu3,
##
## and the higher central moments evaluated at the estimates themselves.
## HALF-NORMAL ONLY: the paper derives it for that case and the moment
## inversion is distribution-specific, so NE and NG keep NA and the bootstrap.
.cols_se_nhn <- function(sigma_u, sigma_v, n) {
  na3 <- c(sigma_v = NA_real_, sigma_u = NA_real_, eu = NA_real_)
  if (!is.finite(sigma_u) || !is.finite(sigma_v) ||
    sigma_u <= 0 || sigma_v <= 0 || !is.finite(n) || n < 5) {
    return(na3)
  }
  s2u <- sigma_u^2
  ssq <- s2u + sigma_v^2
  g <- s2u / ssq
  ## The shared formula ALREADY returns the production-frontier signs: its
  ## (1 - 4/pi) factor is negative, so mu3 and mu5 come back negative. Flipping
  ## them here would be wrong twice over.
  m <- .nhn_central_moments(g, ssq)
  mu2 <- m[["mu2"]]; mu3 <- m[["mu3"]]; mu4 <- m[["mu4"]]
  mu5 <- m[["mu5"]]; mu6 <- m[["mu6"]]

  Vm2 <- (mu4 - mu2^2) / n
  Vm3 <- (mu6 - mu3^2 - 6 * mu2 * mu4 + 9 * mu2^3) / n
  Cm <- (mu5 - 4 * mu2 * mu3) / n

  ## sigma_u = (c m3)^(1/3) with c = 1/(sqrt(2/pi)(1 - 4/pi)); it depends on
  ## m3 alone, which is why the wrong-skew boundary is a m3 condition.
  ## cst is NEGATIVE, because 1 - 4/pi is, so cst^(1/3) is NaN in R even though
  ## (cst * m3)^(1/3) is perfectly real. Write the derivative through sigma_u
  ## instead: sigma_u^3 = cst * m3, so d sigma_u / d m3 = cst / (3 sigma_u^2),
  ## which never takes a fractional power of a negative number.
  cst <- 1 / (sqrt(2 / pi) * (1 - 4 / pi))
  du_m3 <- cst / (3 * s2u)
  ## sigma_v = sqrt(m2 - (1 - 2/pi) sigma_u^2)
  cc <- 1 - 2 / pi
  dv_m2 <- 1 / (2 * sigma_v)
  dv_m3 <- -cc * sigma_u * du_m3 / sigma_v

  var_u <- du_m3^2 * Vm3
  var_v <- dv_m2^2 * Vm2 + dv_m3^2 * Vm3 + 2 * dv_m2 * dv_m3 * Cm
  var_eu <- (2 / pi) * var_u ## E[u] = sigma_u sqrt(2/pi) shifts the intercept
  c(
    sigma_v = sqrt(max(0, var_v)), sigma_u = sqrt(max(0, var_u)),
    eu = sqrt(max(0, var_eu))
  )
}
