## TIC() and vuong(): choosing among model_names when the candidates are not
## nested and none of them need be correct. Gap L1, from Lai and Huang (2010),
## J Prod Anal 34:3-13. See notes/code_history/model_selection.md.
##
## The package ships fifty-odd model_names and, until now, no principled way to
## choose between them. AIC() and BIC() are already available but both assume
## the fitted model is correctly specified -- which is precisely the assumption
## in doubt when someone is deciding between NHN, NE, NG and NTN. Worse, most
## of those pairs are NOT nested, so the standard likelihood ratio test has no
## chi-square limit either.
##
## Two things are implemented, and they answer different questions:
##
##   TIC()   -- Takeuchi's criterion. Like AIC, but the penalty is estimated
##              from the data rather than assumed equal to the parameter count.
##              Lower is better. Comparable across non-nested models.
##   vuong() -- a pairwise TEST, with a p-value, for whether two competing
##              specifications are distinguishable at all.
##
## A criterion always names a winner; a test is allowed to say "these two
## models are not distinguishable on this data", which is frequently the honest
## answer and is not available from AIC alone.

## tr[H(theta) I(theta)^-1], the Takeuchi penalty.
##
## I is the sample Fisher information -- the averaged negative Hessian of the
## log-likelihood -- and H the averaged outer product of the per-observation
## scores. Under CORRECT specification the information matrix equality gives
## H = I, the trace collapses to the number of parameters, and TIC equals AIC.
## Under misspecification the two disagree and the penalty moves; that gap is
## the whole content of the criterion, and `TIC(fit, detail = TRUE)` returns it
## so the equality can be inspected rather than trusted.
##
## Both matrices are used UNNORMALIZED: I = I_sum/n and H = H_sum/n, so the
## 1/n cancels in H I^-1 and dividing twice would be a wasted step, not a
## correction.
.tic_penalty <- function(object) {
  if (is.null(object$opt) || is.null(object$opt$hessian)) {
    stop("TIC(): this fit stores no Hessian, so the Takeuchi penalty cannot ",
      "be formed. psfm()'s moment-based and FE estimators have none by ",
      "construction; a maximum-likelihood fit may have lost it to a failed ",
      "third optimizer stage (convergence code 99), in which case refit with ",
      "different starting values.",
      call. = FALSE
    )
  }
  G <- estfun.sfareg(object)
  I_sum <- .safe_symmetrize(object$opt$hessian)
  if (any(!is.finite(I_sum))) {
    stop("TIC(): the stored Hessian is not finite, so the Takeuchi penalty ",
      "is undefined for this fit.",
      call. = FALSE
    )
  }
  H_sum <- crossprod(G)
  Iinv <- .safe_inverse(I_sum, name = "Hessian")$value
  sum(diag(H_sum %*% Iinv))
}

TIC <- function(object, detail = FALSE) {
  if (!inherits(object, "sfareg")) {
    stop("`object` must be an \"sfareg\" fit.", call. = FALSE)
  }
  ll <- as.numeric(stats::logLik(object))
  if (!is.finite(ll)) {
    stop("TIC(): this fit has no log-likelihood (psfm()'s GTRE_SEQ1, ",
      "GTRE_SEQ2 and SSFE are not maximum likelihood), so no ",
      "likelihood-based criterion applies to it.",
      call. = FALSE
    )
  }
  pen <- .tic_penalty(object)
  val <- -2 * ll + 2 * pen
  if (!isTRUE(detail)) {
    return(val)
  }
  ## `ratio` is the diagnostic worth reading: 1 means the information matrix
  ## equality holds and TIC has nothing to add over AIC. Far from 1 means the
  ## distributional assumption is doing real damage.
  df <- length(object$coefficients)
  list(
    TIC = val,
    AIC = -2 * ll + 2 * df,
    logLik = ll,
    df = df,
    penalty = pen,
    ratio = pen / df
  )
}

## Per-observation log-likelihood, which both routines below are built on.
.ll_i <- function(object, what = "vuong()") {
  if (is.null(object$objective)) {
    stop(what, ": this fit does not retain its likelihood, so per-observation ",
      "log-likelihoods are unavailable. Refit with keep_objective = TRUE.",
      call. = FALSE
    )
  }
  par <- object$opt$par
  if (is.null(par)) {
    stop(what, ": this fit has no parameter vector.", call. = FALSE)
  }
  v <- tryCatch(object$objective(par, per_obs = TRUE), error = function(e) NULL)
  if (is.null(v)) {
    stop(what, ": the stored likelihood does not accept `per_obs`. It was ",
      "retained by a version older than 1.2.0; refit.",
      call. = FALSE
    )
  }
  as.numeric(v)
}

## Vuong's (1989) non-nested likelihood ratio test, in the form Lai and Huang
## (2010) give as their equations (13) and (14).
##
##   m_i = log f_i(theta_f) - log g_i(theta_g)
##   T   = n^-1/2 * sum(m_i) / sd(m)          -> N(0, 1) under H0
##
## H0 is that the two models are equally close to the truth in the
## Kullback-Leibler sense. NEITHER model has to be correctly specified, which
## is the entire reason to prefer this over a standard LR test here.
##
## The variance in (14) is the ML one, dividing by n rather than n-1. Kept as
## published rather than "corrected" to the unbiased form: the statistic's
## asymptotic distribution is stated for this estimator, and at the sample
## sizes this package sees the difference is far below the decision margin
## anyway.
vuong <- function(object1, object2,
                  correction = c("none", "aic", "tic"),
                  level = 0.05) {
  correction <- match.arg(correction)
  for (o in list(object1, object2)) {
    if (!inherits(o, "sfareg")) {
      stop("`object1` and `object2` must both be \"sfareg\" fits.",
        call. = FALSE
      )
    }
  }
  if (length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) {
    stop("`level` must be a single number in (0, 1).", call. = FALSE)
  }

  m1 <- .ll_i(object1)
  m2 <- .ll_i(object2)
  if (length(m1) != length(m2)) {
    stop("vuong(): the two fits use different numbers of observations (",
      length(m1), " and ", length(m2), "). Vuong's test compares two models ",
      "on the SAME data, observation by observation; a difference here ",
      "usually means one fit dropped rows to missingness that the other kept.",
      call. = FALSE
    )
  }
  ## Same length is necessary but not sufficient -- two fits on different
  ## datasets of equal size would pass it and produce a meaningless number.
  d1 <- if (is.null(object1$call$data)) "" else deparse(object1$call$data)
  d2 <- if (is.null(object2$call$data)) "" else deparse(object2$call$data)
  if (nzchar(d1) && nzchar(d2) && !identical(d1, d2)) {
    warning("vuong(): the two fits name different `data` arguments (",
      d1, " and ", d2, "). The test is only meaningful when both models are ",
      "fitted to the same observations.",
      call. = FALSE
    )
  }

  n <- length(m1)
  m <- m1 - m2
  if (any(!is.finite(m))) {
    stop("vuong(): ", sum(!is.finite(m)), " of ", n, " per-observation ",
      "log-likelihood differences are not finite.",
      call. = FALSE
    )
  }
  lr <- sum(m)
  ## Equation (14): the ML variance, dividing by n.
  omega <- sqrt(mean(m^2) - mean(m)^2)
  if (!is.finite(omega) || omega <= 0) {
    stop("vuong(): the two models give identical per-observation ",
      "log-likelihoods, so the test statistic is 0/0. This happens when the ",
      "same model is passed twice, or when two model_names are ",
      "reparameterizations of one another.",
      call. = FALSE
    )
  }

  ## Lai and Huang note that the plain statistic penalizes nothing, so it can
  ## drift toward whichever model has more parameters. The Akaike and Takeuchi
  ## bias corrections address that -- but they say in the same paragraph that
  ## "further investigation of the asymptotic distribution of the modified
  ## statistic is needed". So the corrections are offered and the p-value that
  ## comes with them is flagged, rather than being reported as though the
  ## normal limit had been established for it.
  pen1 <- pen2 <- 0
  if (correction == "aic") {
    pen1 <- length(object1$coefficients)
    pen2 <- length(object2$coefficients)
  } else if (correction == "tic") {
    pen1 <- .tic_penalty(object1)
    pen2 <- .tic_penalty(object2)
  }
  lr_adj <- lr - (pen1 - pen2)
  stat <- lr_adj / (sqrt(n) * omega)
  crit <- stats::qnorm(1 - level / 2)
  p <- 2 * stats::pnorm(-abs(stat))

  nm1 <- if (is.null(object1$model_name)) "model 1" else object1$model_name
  nm2 <- if (is.null(object2$model_name)) "model 2" else object2$model_name
  decision <- if (stat > crit) {
    nm1
  } else if (stat < -crit) {
    nm2
  } else {
    "neither"
  }

  out <- list(
    statistic = stat, p.value = p, n = n,
    lr = lr, lr_adjusted = lr_adj, omega = omega,
    correction = correction, penalty = c(pen1, pen2),
    critical = crit, level = level,
    models = c(nm1, nm2), favoured = decision,
    exact_null = correction == "none"
  )
  class(out) <- "sfa_vuong"
  out
}

print.sfa_vuong <- function(x, ...) {
  cat("\nVuong test for non-nested stochastic frontier models\n")
  cat("Lai and Huang (2010), J Prod Anal 34:3-13, eq. (13)\n\n")
  cat("  model 1:      ", x$models[1], "\n", sep = "")
  cat("  model 2:      ", x$models[2], "\n", sep = "")
  cat("  observations: ", x$n, "\n", sep = "")
  cat("  correction:   ", x$correction, "\n\n", sep = "")
  cat(sprintf("  log-likelihood difference : %.4f\n", x$lr))
  if (!identical(x$correction, "none")) {
    cat(sprintf("  after bias correction     : %.4f  (penalties %.2f, %.2f)\n",
      x$lr_adjusted, x$penalty[1], x$penalty[2]
    ))
  }
  cat(sprintf("  z                         : %.4f\n", x$statistic))
  cat(sprintf("  p-value (two-sided)       : %.4g\n", x$p.value))
  cat(sprintf("  critical value at %.0f%%     : +/- %.4f\n\n",
    100 * x$level, x$critical
  ))
  if (identical(x$favoured, "neither")) {
    cat("  Neither model is favoured: the data does not distinguish them.\n")
  } else {
    cat("  Favoured: ", x$favoured, "\n", sep = "")
  }
  if (!x$exact_null) {
    cat("\n  NOTE: the bias-corrected statistic shifts the mean of the test.\n")
    cat("  Lai and Huang derive the N(0, 1) limit for the UNCORRECTED\n")
    cat("  statistic only, and state that the corrected one needs further\n")
    cat("  investigation. Read this p-value as indicative, not exact.\n")
  }
  cat("\n")
  invisible(x)
}
