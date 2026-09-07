## sfma(): average over inefficiency distributions instead of picking one.
## Gap L8, from Parmeter, Wan and Zhang (2019), J Prod Anal 51:91-103.
## See notes/code_history/sfma.md.
##
## TIC() and vuong() answer "which distribution?". This answers the harder
## question of what to do when the honest reply is "the data does not say".
## Selecting a model and then reporting its efficiencies as though the choice
## had been free understates the uncertainty: the standard errors condition on
## a decision that was itself estimated.
##
## The package is unusually well placed for this. The whole point of putting
## ten inefficiency distributions behind one interface is that fitting all of
## them is a loop, and the averaging then costs one small quadratic program.

## Euclidean projection onto the unit simplex (Duchi et al. 2008). Sorting
## once and finding the threshold is exact -- not an iterative approximation --
## which matters because the weights are the output, and a weight that should
## be exactly zero should BE exactly zero.
.proj_simplex <- function(v) {
  u <- sort(v, decreasing = TRUE)
  cs <- cumsum(u) - 1
  idx <- seq_along(u)
  rho <- max(which(u - cs / idx > 0))
  pmax(v - cs[rho] / rho, 0)
}

## min over the simplex of ||A w - r||^2 + c' w, by projected gradient.
##
## Convex, so a fixed step of 1/L with L the Lipschitz constant of the
## gradient converges without a line search. Written out rather than taken
## from a QP package because the only ones available are not in Imports, and
## fifteen lines that can be checked against a brute-force grid is a better
## trade than a new dependency.
.simplex_qp <- function(A, r, cvec, maxit = 5000, tol = 1e-12) {
  S <- ncol(A)
  AtA <- crossprod(A)
  Atr <- as.numeric(crossprod(A, r))
  L <- 2 * max(abs(eigen(AtA, symmetric = TRUE, only.values = TRUE)$values))
  if (!is.finite(L) || L <= 0) L <- 1
  w <- rep(1 / S, S)
  for (it in seq_len(maxit)) {
    g <- 2 * (AtA %*% w - Atr) + cvec
    w_new <- .proj_simplex(as.numeric(w - g / L))
    if (max(abs(w_new - w)) < tol) {
      w <- w_new
      break
    }
    w <- w_new
  }
  as.numeric(w)
}

## Per-observation inefficiency prediction, E[u | eps], for one fit.
.sfma_rho <- function(fit) {
  v <- fit$u_hat
  if (is.null(v)) v <- fit$jlms
  if (is.null(v)) {
    post <- fit$u_posterior
    if (!is.null(post) && !is.null(post$mu_star)) {
      v <- .jlms_u(post$mu_star, post$sigma_star)
    }
  }
  if (is.null(v)) {
    e <- fit$exp_u_hat
    if (!is.null(e)) v <- -log(pmax(as.numeric(e), .SFA_CONSTANTS$MIN_POSITIVE))
  }
  if (is.null(v)) return(NULL)
  as.numeric(v)
}

sfma <- function(formula, data,
                 models = c("NHN", "NE", "NTN", "NR", "NGE"),
                 weights = c("sfma", "tic", "aic", "bic", "equal"),
                 inefdec = TRUE, quiet = FALSE, ...) {
  weights <- match.arg(weights)
  if (!is.character(models) || length(models) < 2L) {
    stop("`models` must name at least two model_names to average over.",
      call. = FALSE
    )
  }
  models <- unique(models)

  fits <- list()
  for (m in models) {
    f <- tryCatch(
      suppressWarnings(sfm(formula,
        model_name = m, data = data,
        inefdec = inefdec, keep_objective = TRUE, ...
      )),
      error = function(e) NULL
    )
    if (is.null(f)) {
      if (!quiet) message("sfma: ", m, " failed to fit and is dropped")
      next
    }
    if (is.null(.sfma_rho(f))) {
      if (!quiet) message("sfma: ", m, " reports no E[u | e] and is dropped")
      next
    }
    fits[[m]] <- f
  }
  if (length(fits) < 2L) {
    stop("sfma(): fewer than two candidate models fitted successfully, so ",
      "there is nothing to average.",
      call. = FALSE
    )
  }

  R <- vapply(fits, .sfma_rho, numeric(length(.sfma_rho(fits[[1L]]))))
  if (!is.matrix(R)) R <- matrix(R, ncol = length(fits))
  colnames(R) <- names(fits)
  n <- nrow(R)
  k <- vapply(fits, function(f) length(f$coefficients), numeric(1))
  ll <- vapply(fits, function(f) {
    v <- tryCatch(as.numeric(stats::logLik(f)), error = function(e) NA_real_)
    if (is.null(v)) NA_real_ else v
  }, numeric(1))

  ## Information-criterion weights, in Buckland et al.'s (1997) smoothed form:
  ## w_s proportional to exp(-IC_s / 2), normalised. TIC is the interesting one
  ## here, because unlike AIC it does not assume the candidate is correct --
  ## which is the whole premise of averaging in the first place.
  ic_weights <- function(ic) {
    d <- ic - min(ic, na.rm = TRUE)
    w <- exp(-d / 2)
    w[!is.finite(w)] <- 0
    w / sum(w)
  }

  crit <- NULL
  if (weights == "equal") {
    w <- rep(1 / length(fits), length(fits))
  } else if (weights == "aic") {
    w <- ic_weights(-2 * ll + 2 * k)
  } else if (weights == "bic") {
    w <- ic_weights(-2 * ll + log(n) * k)
  } else if (weights == "tic") {
    tic <- vapply(fits, function(f) {
      tryCatch(TIC(f), error = function(e) NA_real_)
    }, numeric(1))
    if (all(!is.finite(tic))) {
      stop("sfma(weights = \"tic\"): no candidate produced a usable Takeuchi ",
        "penalty. Refit with optHessian = TRUE.",
        call. = FALSE
      )
    }
    tic[!is.finite(tic)] <- max(tic, na.rm = TRUE) + 1e3
    w <- ic_weights(tic)
  } else {
    ## Parmeter, Wan and Zhang's criterion (their eq. 13):
    ##
    ##   C(w) = || rho(w) - rho_full ||^2  +  n^{1/2} log(n) k'w
    ##
    ## The unobservable target is replaced by the prediction from the FULL
    ## model, exactly as Mallows' Cp replaces an unknown sigma^2 with the full
    ## model's estimate. Without the penalty all the weight goes to the largest
    ## candidate by construction.
    ##
    ## "Full" is the candidate with the most parameters, which is the nested
    ## reading of their setup -- their Theorem 1 is about models that NEST the
    ## truth.
    ##
    ## That requires a unique largest candidate, and refusing without one is
    ## not pedantry. If every candidate has the same k the penalty
    ## n^{1/2} log(n) k'w is the SAME for all of them, so it drops out as a
    ## constant and the criterion collapses to "be as close as possible to
    ## rho_full" -- which is minimised by putting weight 1 on whichever model
    ## happened to be chosen as full. Measured on exponential data with
    ## candidates NHN/NE/NR (all k = 5) that put 0.879 on NHN, the wrong
    ## distribution, while TIC gave NE 0.934. A silent tie-break would have
    ## produced a confident wrong answer.
    if (sum(k == max(k)) > 1L) {
      stop("sfma(weights = \"sfma\"): the criterion of Parmeter, Wan and ",
        "Zhang needs a single most general candidate to stand in for the ",
        "unobservable target, and ", sum(k == max(k)), " of these models tie ",
        "at ", max(k), " parameters. With equal dimensions its penalty is ",
        "constant and the criterion degenerates to putting all the weight on ",
        "whichever model is nominated. Use weights = \"tic\" (or \"aic\"), ",
        "which compare the candidates directly, or include a genuinely more ",
        "general model such as \"NTN\".",
        call. = FALSE
      )
    }
    full <- names(which.max(k))
    r <- R[, full]
    cvec <- sqrt(n) * log(n) * as.numeric(k)
    w <- .simplex_qp(R, r, cvec)
    crit <- list(
      full_model = full,
      objective = sum((R %*% w - r)^2) + sum(cvec * w),
      penalty = sum(cvec * w)
    )
  }
  names(w) <- names(fits)

  u_avg <- as.numeric(R %*% w)
  out <- list(
    weights = w,
    u_hat = u_avg,
    efficiency = exp(-u_avg),
    u_by_model = R,
    fits = fits,
    models = names(fits),
    k = k, logLik = ll, nobs = n,
    weight_type = weights,
    criterion = crit,
    formula = formula,
    call = match.call()
  )
  class(out) <- "sfma"
  out
}

print.sfma <- function(x, ...) {
  cat("\nStochastic frontier model averaging over inefficiency distributions\n")
  cat("Parmeter, Wan & Zhang (2019), J Prod Anal 51:91-103\n\n")
  cat("  weighting : ", x$weight_type, "\n", sep = "")
  cat("  candidates: ", length(x$models), " of ",
    length(x$models), " fitted\n", sep = ""
  )
  cat("  obs       : ", x$nobs, "\n\n", sep = "")
  d <- data.frame(
    model = x$models,
    k = as.integer(x$k),
    logLik = round(x$logLik, 2),
    weight = round(x$weights, 4)
  )
  d <- d[order(-d$weight), ]
  print(d, row.names = FALSE)
  if (!is.null(x$criterion)) {
    cat("\n  full model used as the target: ", x$criterion$full_model, "\n",
      sep = ""
    )
  }
  cat(sprintf("\n  mean averaged efficiency: %.4f  (range %.4f to %.4f)\n",
    mean(x$efficiency), min(x$efficiency), max(x$efficiency)
  ))
  ## The spread across candidates is the reason to average at all: if the
  ## models disagree, a single-model standard error is understating things.
  sp <- apply(exp(-x$u_by_model), 1, function(z) max(z) - min(z))
  cat(sprintf("  per-firm spread across candidates: median %.4f, max %.4f\n",
    stats::median(sp), max(sp)
  ))
  cat("\n")
  invisible(x)
}
