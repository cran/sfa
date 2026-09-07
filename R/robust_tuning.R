## Choosing the robustness tuning parameter: data-driven (H-score) and fixed
## (weight-matching calibration). Gap L17, from Bernstein, Parmeter and Wright
## (2026). See notes/code_history/robust_tuning.md.
##
## Two silent failure modes are guarded here, both of which otherwise return a
## converged fit and a finite criterion value:
##
##   - Grid loss. Evaluated naively the H-score is non-finite wherever a fitted
##     density underflows, which discards candidates NON-randomly -- the ones
##     nearest maximum likelihood survive -- and so understates the robustness
##     the data call for. hscore() evaluates in log space instead.
##   - Warm-start capture. Walking the grid upward and warm-starting each fit
##     from the previous one is efficient but unguarded: once one fit falls into
##     the degenerate sigma_u -> 0 basin, every later fit inherits it. With
##     multistart = TRUE each candidate is fitted from both the warm start and
##     the original MLE fit, keeping the better attained objective.
##
## The starts are passed as "sfareg" OBJECTS, not parameter vectors, because
## that is what sfm(start_from = ) takes -- it matches by parameter name (see
## .start_from() in matrix_utils.R), so a numeric vector is rejected outright.

hscore_select <- function(object,
                          method = c("mlqe", "psi", "mdpd"),
                          range = c(0, 0.60),
                          coarse = 0.01,
                          fine = 0.001,
                          half_width = 0.02,
                          multistart = TRUE,
                          verbose = FALSE) {
  method <- match.arg(method)
  stopifnot(length(range) == 2L, range[2] > range[1], range[1] >= 0)

  refit <- .robust_refit_fn(object, method)
  ref_start <- .robust_ref_start(object)

  score_at <- function(cc, warm) {
    fits <- list(refit(cc, warm))
    if (isTRUE(multistart) && !identical(warm, ref_start))
      fits[[2]] <- refit(cc, ref_start)
    ok <- vapply(fits, function(f) isTRUE(f$ok), logical(1))
    if (!any(ok)) return(list(fit = NULL, hs = NA_real_))
    fits <- fits[ok]
    best <- fits[[which.max(vapply(fits, function(f) f$objective, numeric(1)))]]
    list(fit = best,
         hs = hscore(best$residuals, best$sigma_v, best$sigma_u, cc))
  }

  grid <- seq(range[1], range[2], by = coarse)
  path <- data.frame(c = grid, hscore = NA_real_, stage = "coarse",
                     stringsAsFactors = FALSE)
  warm <- ref_start
  best <- list(c = NA_real_, hs = Inf, fit = NULL)
  for (i in seq_along(grid)) {
    r <- score_at(grid[i], warm)
    path$hscore[i] <- r$hs
    if (!is.null(r$fit)) warm <- r$fit$fit
    if (is.finite(r$hs) && r$hs < best$hs)
      best <- list(c = grid[i], hs = r$hs, fit = r$fit)
    if (verbose) message(sprintf("  c = %.3f   H = %s", grid[i],
                                 format(r$hs, digits = 6)))
  }
  if (!is.finite(best$hs))
    stop("The H-score was not finite at any candidate. ",
         "Check the fitted scale parameters.", call. = FALSE)

  lo <- max(range[1], best$c - half_width)
  hi <- min(range[2], best$c + half_width)
  fgrid <- setdiff(round(seq(lo, hi, by = fine), 10), round(grid, 10))
  if (length(fgrid)) {
    fpath <- data.frame(c = fgrid, hscore = NA_real_, stage = "fine",
                        stringsAsFactors = FALSE)
    for (i in seq_along(fgrid)) {
      r <- score_at(fgrid[i], best$fit$fit)
      fpath$hscore[i] <- r$hs
      if (is.finite(r$hs) && r$hs < best$hs)
        best <- list(c = fgrid[i], hs = r$hs, fit = r$fit)
    }
    path <- rbind(path, fpath)
  }
  path <- path[order(path$c), ]
  rownames(path) <- NULL

  out <- list(c = best$c, hscore = best$hs, path = path, fit = best$fit$fit,
              method = method, range = range,
              n_finite = c(finite = sum(is.finite(path$hscore)),
                           tried = nrow(path)))
  class(out) <- "sfa_hscore"
  out
}

## Report how much of the grid was scorable: a search that lost most of its
## candidates has not chosen among the ones it was asked to.
print.sfa_hscore <- function(x, ...) {
  cat("Hyvarinen-score tuning selection\n")
  cat("  criterion : ", toupper(x$method), "\n", sep = "")
  cat("  search    : [", x$range[1], ", ", x$range[2], "]\n", sep = "")
  cat("  scorable  : ", x$n_finite[["finite"]], " of ",
      x$n_finite[["tried"]], " candidates\n", sep = "")
  cat("  selected c: ", format(x$c, digits = 4),
      if (isTRUE(all.equal(x$c, 0))) "   (maximum likelihood endpoint)" else "",
      "\n", sep = "")
  cat("  H(c)      : ", format(x$hscore, digits = 6), "\n", sep = "")
  invisible(x)
}

calibrate_c <- function(sigma_v, sigma_u,
                        method = c("mlqe", "psi", "mdpd", "all"),
                        target = 0.10, k = 3, range = c(1e-4, 0.6)) {
  method <- match.arg(method)
  if (!is.finite(sigma_v) || sigma_v <= 0) stop("sigma_v must be positive.")
  if (!is.finite(sigma_u) || sigma_u <= 0) stop("sigma_u must be positive.")
  ms <- if (method == "all") c("mlqe", "psi", "mdpd") else method
  res <- lapply(ms, function(m) .calibrate_one(m, sigma_v, sigma_u, target, k, range))
  out <- stats::setNames(vapply(res, as.numeric, numeric(1)), ms)
  attr(out, "roots") <- stats::setNames(lapply(res, attr, "roots"), ms)
  out
}

## Per-observation data term at the mean-objective scale, reparameterised in
## (intercept shift, log sigma_v, log sigma_u) so it can be differentiated in
## all three directions. The integral correction is NOT divided by n: these are
## mean-scale objectives, so it sits at the same order as the data term. Scaling
## it down erases the very asymmetry across criteria the calibration exists to
## capture, and collapses the answer to a common c.
.per_obs_term <- function(method, e_base, b, log_sv, log_su, c) {
  sigma_v <- exp(log_sv); sigma_u <- exp(log_su)
  e <- e_base - b
  f <- .dens_nhn(e, sigma_v, sigma_u)
  switch(method,
    mlqe = (f^c - 1) / c,
    psi  = f^c / c - .nhn_power_integral(sigma_v, sigma_u, 1 + c) / (1 + c),
    mdpd = ((1 + c) / c) * f^c - .nhn_power_integral(sigma_v, sigma_u, 1 + c)
  )
}

.score_grad <- function(method, e_base, log_sv, log_su, c, h = 1e-4) {
  f <- function(b, lv, lu) .per_obs_term(method, e_base, b, lv, lu, c)
  c((f( h, log_sv, log_su) - f(-h, log_sv, log_su)) / (2 * h),
    (f(0, log_sv + h, log_su) - f(0, log_sv - h, log_su)) / (2 * h),
    (f(0, log_sv, log_su + h) - f(0, log_sv, log_su - h)) / (2 * h))
}

.weight_ratio <- function(method, sigma_v, sigma_u, c, k) {
  sigma <- sqrt(sigma_v^2 + sigma_u^2)
  lv <- log(sigma_v); lu <- log(sigma_u)
  g0 <- .score_grad(method, 0, lv, lu, c)
  gk <- .score_grad(method, k * sigma, lv, lu, c)
  sqrt(sum(gk^2)) / sqrt(sum(g0^2))
}

.calibrate_one <- function(method, sigma_v, sigma_u, target, k, range) {
  obj <- function(c) .weight_ratio(method, sigma_v, sigma_u, c, k) - target
  ## The ratio is NOT monotone in c for the Fisher-consistency-corrected
  ## criteria: the integral correction pulls it back up again, so Psi and MDPD
  ## cross the target twice, at sigma_v = 0.3, sigma_u = 0.6 near c = 0.239 and
  ## again near c = 0.539. uniroot() on the whole interval is therefore refused
  ## for a reason -- the endpoints share a sign -- and the earlier version fell
  ## through to picking the grid point of smallest |obj|, which is accurate only
  ## to the grid step (it returned 0.5397 for a root at 0.53854, and matched the
  ## target to 1.5e-4 rather than to solver precision).
  ##
  ## Scan for sign changes instead and solve each bracket exactly. The LARGEST
  ## root is returned: it is the calibration reported in the paper, and it is
  ## the one on the far side of the dip, so the target down-weighting holds
  ## there rather than being passed through on the way. All roots are attached
  ## so the non-monotonicity is visible rather than silently resolved.
  g <- seq(range[1], range[2], length.out = 400)
  v <- vapply(g, obj, numeric(1))
  fin <- is.finite(v)
  cuts <- which(fin[-length(fin)] & fin[-1] & diff(sign(v)) != 0)
  if (!length(cuts)) {
    ## No crossing anywhere in range: report the closest approach, and say so.
    if (!any(fin)) {
      stop("calibrate_c(): the influence ratio could not be evaluated anywhere ",
        "in [", range[1], ", ", range[2], "].", call. = FALSE)
    }
    best <- g[fin][which.min(abs(v[fin]))]
    warning("calibrate_c(): the influence ratio never reaches the target of ",
      target, " in [", range[1], ", ", range[2], "]; returning the closest ",
      "approach (ratio ", format(min(abs(v[fin])) + target, digits = 4), ").",
      call. = FALSE)
    return(structure(best, roots = numeric(0)))
  }
  roots <- vapply(cuts, function(i) {
    stats::uniroot(obj, interval = c(g[i], g[i + 1]), tol = 1e-10)$root
  }, numeric(1))
  structure(max(roots), roots = roots)
}
