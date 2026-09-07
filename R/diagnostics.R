## Optimizer diagnostics.

## Hessian of the NEGATIVE log-likelihood, as optim() returns it.
.sfa_hess_diag <- function(H, pnames) {
  if (is.null(H) || !is.matrix(H) || any(!is.finite(H))) {
    return(NULL)
  }
  H <- (H + t(H)) / 2 ## symmetrize: optim's numeric Hessian is not exactly so
  ev <- eigen(H, symmetric = TRUE)
  vals <- ev$values
  amin <- min(abs(vals))
  cond <- if (amin > 0) max(abs(vals)) / amin else Inf
  ## Loadings of the flattest direction: which parameters move together along
  ## the least-curved axis of the likelihood.
  flat <- ev$vectors[, which.min(abs(vals))]
  names(flat) <- pnames
  list(
    eigenvalues = stats::setNames(vals, paste0("e", seq_along(vals))),
    condition = cond,
    pos_def = all(vals > 0),
    flat_direction = flat[order(-abs(flat))]
  )
}

## Implied parameter correlations. Reuses vcov.sfareg so that a fit whose
## Hessian was not requested still reports whatever it can.
.sfa_corr_diag <- function(object, pnames) {
  V <- tryCatch(suppressWarnings(stats::vcov(object)), error = function(e) NULL)
  if (is.null(V) || nrow(V) < 2) {
    return(NULL)
  }
  ## A near-singular Hessian -- exactly the condition this diagnostic exists
  ## to report.
  v <- diag(V)
  keep <- is.finite(v) & v > 0 &
    vapply(seq_len(nrow(V)), function(i) all(is.finite(V[i, ])), logical(1))
  if (sum(keep) < 2L) {
    return(NULL)
  }
  dropped <- pnames[!keep]
  V <- V[keep, keep, drop = FALSE]
  pn <- pnames[keep]
  d <- sqrt(diag(V))
  R <- V / outer(d, d)
  dimnames(R) <- list(pn, pn)
  off <- R
  diag(off) <- 0
  idx <- which(abs(off) == max(abs(off)), arr.ind = TRUE)[1, ]
  list(
    correlation = R,
    worst_pair = c(pn[idx[1]], pn[idx[2]]),
    worst_value = off[idx[1], idx[2]],
    dropped = dropped
  )
}

sfa_diagnostics <- function(object, ...) {
  if (!inherits(object, "sfareg")) {
    stop("sfa_diagnostics() expects an object of class \"sfareg\".", call. = FALSE)
  }
  pnames <- names(object$coefficients)
  if (is.null(pnames)) pnames <- paste0("par", seq_along(object$coefficients))
  opt <- object$opt

  ## ---- convergence ----------------------------------------------------
  code <- if (!is.null(opt$convergence)) opt$convergence else NA_integer_
  conv <- list(
    code = code,
    ## optim()'s documented codes.
    meaning = switch(as.character(code),
      "0" = "converged",
      "1" = "ITERATION LIMIT REACHED -- not a converged optimum",
      "10" = "Nelder-Mead simplex degenerated",
      "51" = "L-BFGS-B warning",
      "52" = "L-BFGS-B ERROR",
      "99" = "final optim() polish stage did not produce a usable result -- preceding stage's estimates reported",
      if (is.na(code)) "no optimizer output (moment-based or FE estimator)" else "unknown code"
    ),
    message = if (!is.null(opt$message)) opt$message else NA_character_,
    counts = opt$counts,
    logLik = if (is.numeric(opt$value)) -opt$value else NA_real_
  )

  hess <- .sfa_hess_diag(opt$hessian, pnames)
  corr <- .sfa_corr_diag(object, pnames)

  ## ---- gradient, only if the objective was retained --------------------
  grad <- NULL
  if (is.function(object$objective) && !is.null(opt$par)) {
    g <- tryCatch(numDeriv::grad(object$objective, opt$par), error = function(e) NULL)
    if (!is.null(g) && all(is.finite(g))) {
      names(g) <- pnames
      ## Scale-free version: a gradient of 1e-3 means very different things
      ## for a parameter of size 0.001 and one of size 1000.
      rel <- abs(g) * pmax(abs(opt$par), 1) / max(abs(conv$logLik), 1)
      grad <- list(gradient = g, max_abs = max(abs(g)),
                   relative = stats::setNames(rel, pnames), max_rel = max(rel))
    }
  }

  ## ---- actionable flags ----------------------------------------------- The
  ## convergence code on its own is NOT diagnostic.
  flags <- character(0)
  hess_ok <- !is.null(hess) && hess$pos_def
  grad_ok <- !is.null(grad) && grad$max_rel <= 1e-3

  ## Code 1 is the iteration limit: curvature at the stopping point cannot
  ## show whether the search had actually finished.
  linesearch <- !is.na(code) && code %in% c(51L, 52L)
  benign <- linesearch && grad_ok && hess_ok
  ## "benign" is a positive claim, so it requires positive evidence.
  unverified <- !is.na(code) && code != 0 && !benign && is.null(grad) && linesearch
  conv$benign_nonzero <- benign
  conv$unverified_nonzero <- unverified

  if (!is.na(code) && code != 0) {
    if (benign) {
      flags <- c(flags, sprintf(paste(
        "Optimizer returned code %s (%s), but the gradient and Hessian both say",
        "this is an optimum. The staged minimizer had already converged and the",
        "final stage could not improve on it; treat the code as noise here."
      ), code, conv$meaning))
    } else if (unverified) {
      flags <- c(flags, sprintf(paste(
        "Optimizer returned code %s (%s). This is often harmless -- the final",
        "stage failing to improve on an already-converged point -- but that",
        "cannot be confirmed without the gradient. Refit with",
        "keep_objective = TRUE to settle it."
      ), code, conv$meaning))
    } else {
      flags <- c(flags, sprintf(
        "Optimizer returned code %s (%s), and the %s says the fit is not at an optimum.",
        code, conv$meaning,
        if (!grad_ok && !is.null(grad)) "gradient" else if (!hess_ok) "Hessian" else "iteration limit"
      ))
    }
  }
  if (!is.null(hess)) {
    if (!hess$pos_def) {
      flags <- c(flags, paste(
        "Hessian is NOT positive definite, so this point is not a minimum of",
        "the negative log-likelihood in every direction. Standard errors from",
        "it are not trustworthy."
      ))
    }
    if (is.finite(hess$condition) && hess$condition > 1e8) {
      flags <- c(flags, sprintf(
        paste("Hessian condition number is %.3g. The likelihood is nearly flat",
              "along at least one direction; the parameters loading on it are",
              "weakly identified (see $hessian$flat_direction)."),
        hess$condition
      ))
    } else if (!is.finite(hess$condition)) {
      flags <- c(flags, "Hessian is singular: at least one direction has zero curvature.")
    }
  } else {
    flags <- c(flags, "No usable Hessian was stored (fit with optHessian = TRUE to get one).")
  }
  if (!is.null(corr) && abs(corr$worst_value) > 0.95) {
    flags <- c(flags, sprintf(
      "%s and %s are correlated at %.3f. The data largely identify a combination of them rather than each separately.",
      corr$worst_pair[1], corr$worst_pair[2], corr$worst_value
    ))
  }
  if (!is.null(grad) && grad$max_rel > 1e-3) {
    flags <- c(flags, sprintf(
      "Largest relative gradient component is %.3g, which is large for a converged optimum.",
      grad$max_rel
    ))
  }
  if (is.null(object$objective)) {
    flags <- c(flags, "Objective not retained; refit with keep_objective = TRUE for the likelihood slice and gradient.")
  }

  out <- list(
    model_name = object$model_name, call = object$call, pnames = pnames,
    estimates = object$coefficients, convergence = conv, hessian = hess,
    correlation = corr, gradient = grad, flags = flags,
    has_objective = is.function(object$objective)
  )
  class(out) <- "sfadiag"
  out
}

print.sfadiag <- function(x, ...) {
  cat("--- sfa optimizer diagnostics ---\n")
  cat("Model:      ", x$model_name, "\n", sep = "")
  cat("Parameters: ", length(x$pnames), " (", paste(x$pnames, collapse = ", "), ")\n", sep = "")

  cat("\nConvergence\n")
  cat("  code       : ", x$convergence$code, "  -- ", x$convergence$meaning, "\n", sep = "")
  if (!is.na(x$convergence$message)) cat("  message    : ", x$convergence$message, "\n", sep = "")
  if (!is.null(x$convergence$counts)) {
    cat("  evaluations: ",
      paste(names(x$convergence$counts), x$convergence$counts, sep = "=", collapse = "  "),
      "\n", sep = "")
  }
  if (is.finite(x$convergence$logLik)) cat("  logLik     : ", signif(x$convergence$logLik, 8), "\n", sep = "")

  if (!is.null(x$hessian)) {
    cat("\nHessian (of the negative log-likelihood)\n")
    cat("  eigenvalues     : ", paste(signif(x$hessian$eigenvalues, 4), collapse = "  "), "\n", sep = "")
    cat("  condition number: ", format(x$hessian$condition, digits = 5), "\n", sep = "")
    cat("  positive definite: ", x$hessian$pos_def, "\n", sep = "")
    fd <- x$hessian$flat_direction
    keep <- fd[abs(fd) > 0.1]
    if (length(keep)) {
      cat("  flattest direction: ",
        paste(sprintf("%s %+.2f", names(keep), keep), collapse = ", "), "\n", sep = "")
    }
  }

  if (!is.null(x$correlation)) {
    cat("\nParameter correlation\n")
    cat("  largest |corr|: ", sprintf("%s vs %s = %+.3f", x$correlation$worst_pair[1],
                                      x$correlation$worst_pair[2], x$correlation$worst_value),
      "\n", sep = "")
    if (length(x$correlation$dropped)) {
      cat("  excluded (no usable variance): ",
        paste(x$correlation$dropped, collapse = ", "), "\n", sep = "")
    }
  }

  if (!is.null(x$gradient)) {
    cat("\nGradient at the reported optimum\n")
    cat("  max |component| : ", format(x$gradient$max_abs, digits = 4), "\n", sep = "")
    cat("  max relative    : ", format(x$gradient$max_rel, digits = 4), "\n", sep = "")
  }

  if (length(x$flags)) {
    cat("\nFlags\n")
    for (f in x$flags) cat("  * ", paste(strwrap(f, width = 76, exdent = 4), collapse = "\n"), "\n", sep = "")
  } else {
    cat("\nNo flags raised.\n")
  }
  invisible(x)
}


## Diagnostic plots.
plot.sfareg <- function(x, which = 1:4, n_grid = 41, span = 0.25, ...) {
  d <- sfa_diagnostics(x)
  avail <- integer(0)
  if (!is.null(d$hessian)) avail <- c(avail, 1L)
  if (!is.null(d$correlation)) avail <- c(avail, 2L)
  if (d$has_objective) avail <- c(avail, 3L, 4L)
  which <- intersect(which, avail)
  if (!length(which)) {
    stop("Nothing to plot: no usable Hessian and no retained objective. ",
      "Refit with optHessian = TRUE and/or keep_objective = TRUE.",
      call. = FALSE
    )
  }

  np <- length(d$pnames)
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  npanel <- length(which) + if (3L %in% which) np - 1L else 0L
  graphics::par(mfrow = grDevices::n2mfrow(npanel), mar = c(4, 4, 2.5, 1))

  if (1L %in% which) {
    ev <- abs(d$hessian$eigenvalues)
    graphics::plot(seq_along(ev), ev,
      log = "y", type = "b", pch = 19, xlab = "eigenvalue (ordered)",
      ylab = "|eigenvalue|", main = "Hessian spectrum"
    )
    graphics::mtext(sprintf("condition number %.3g", d$hessian$condition), side = 3, cex = 0.7)
  }

  if (2L %in% which) {
    ## Dimension from the matrix, not from d$pnames: .sfa_corr_diag() drops
    ## parameters with no usable variance, so R can be smaller than the fit.
    R <- d$correlation$correlation
    nr <- nrow(R)
    rn <- rownames(R)
    graphics::image(seq_len(nr), seq_len(nr), abs(R[, nr:1, drop = FALSE]),
      zlim = c(0, 1), axes = FALSE, xlab = "", ylab = "",
      main = "|parameter correlation|",
      col = grDevices::hcl.colors(20, "Blues", rev = TRUE)
    )
    graphics::axis(1, at = seq_len(nr), labels = rn, las = 2, cex.axis = 0.7)
    graphics::axis(2, at = seq_len(nr), labels = rev(rn), las = 2, cex.axis = 0.7)
    graphics::box()
  }

  if (3L %in% which) {
    par_hat <- x$opt$par
    for (j in seq_len(np)) {
      lo <- par_hat[j] - span * max(abs(par_hat[j]), 1)
      hi <- par_hat[j] + span * max(abs(par_hat[j]), 1)
      grid <- seq(lo, hi, length.out = n_grid)
      ll <- vapply(grid, function(v) {
        p <- par_hat
        p[j] <- v
        val <- tryCatch(x$objective(p), error = function(e) NA_real_)
        if (is.null(val) || !is.finite(val)) NA_real_ else -val
      }, numeric(1))
      graphics::plot(grid, ll,
        type = "l", xlab = d$pnames[j], ylab = "log-likelihood",
        main = paste("slice:", d$pnames[j])
      )
      graphics::abline(v = par_hat[j], lty = 2)
      graphics::points(par_hat[j], -x$opt$value, pch = 19, col = "red")
    }
  }

  if (4L %in% which) {
    g <- d$gradient$gradient
    graphics::barplot(g,
      names.arg = d$pnames, las = 2, cex.names = 0.7,
      main = "gradient at the optimum", ylab = "d(-logLik)/d(par)"
    )
    graphics::abline(h = 0)
  }
  invisible(d)
}
