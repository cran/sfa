## simulation_se(): how much of a fitted parameter's uncertainty is SIMULATION
## noise rather than sampling noise. Gap J4.
## See notes/code_history/simulation_se.md.
##
## Simulated maximum likelihood replaces an integral with an average over draws,
## so a reported standard error from the Hessian answers "how much would this
## estimate move in a new SAMPLE" and says nothing about "how much would it move
## with the same sample and different DRAWS". Randomizing the sequence makes the
## second question answerable: refit K times with independent randomizations and
## take the standard deviation across refits.
##
## This is only possible because the draws are randomized by SHIFTING. An
## unrandomized Halton lattice gives the same answer every time, so the
## diagnostic would report exactly zero and mean nothing by it.

## Which argument randomizes the draws, and whether this model simulates at all.
## NULL means the fit is closed-form, in which case there is no simulation
## variance to measure and saying so beats reporting a column of zeros.
.sim_seed_arg <- function(object) {
  fn <- tryCatch(as.character(object$call[[1L]]), error = function(e) "")
  fn <- fn[length(fn)]
  m <- object$model_name
  switch(fn,
    sfm = if (m %in% c("NLN", "NW")) "sim_seed" else NULL,
    psfm = if (m %in% c("GTRE", "TRE", "TRE_Z", "GTRE_Z")) "rand.gtre" else NULL,
    selsfm = "seed",
    NULL
  )
}

simulation_se <- function(object, K = 10, seeds = NULL, quiet = FALSE,
                          envir = parent.frame()) {
  if (!inherits(object, "sfareg")) {
    stop("`object` must be an \"sfareg\" fit.", call. = FALSE)
  }
  if (length(K) != 1L || !is.numeric(K) || !is.finite(K) || K < 2) {
    stop("`K` must be a single number >= 2: a standard deviation across ",
      "refits needs at least two refits.",
      call. = FALSE
    )
  }
  K <- as.integer(K)
  arg <- .sim_seed_arg(object)
  if (is.null(arg)) {
    stop("simulation_se(): model_name \"",
      if (is.null(object$model_name)) "unknown" else object$model_name,
      "\" is not estimated by simulation, so it has no simulation variance. ",
      "This diagnostic applies to sfm(\"NLN\"/\"NW\"), psfm(\"GTRE\"/\"TRE\") ",
      "and selsfm().",
      call. = FALSE
    )
  }
  if (is.null(object$call)) {
    stop("simulation_se(): this fit did not retain its call, so it cannot be ",
      "refitted.",
      call. = FALSE
    )
  }
  if (is.null(seeds)) seeds <- seq_len(K) * 7919L else K <- length(seeds)
  if (length(seeds) != K) {
    stop("`seeds` must have length `K`.", call. = FALSE)
  }

  base <- object$coefficients
  nm <- names(base)
  M <- matrix(NA_real_, K, length(base), dimnames = list(NULL, nm))

  for (k in seq_len(K)) {
    cl <- object$call
    cl[[arg]] <- seeds[k]
    fit <- tryCatch(suppressWarnings(eval(cl, envir = envir)),
      error = function(e) NULL
    )
    if (is.null(fit) || is.null(fit$coefficients)) next
    if (length(fit$coefficients) == length(base)) M[k, ] <- fit$coefficients
    if (!quiet) message("simulation_se: refit ", k, " of ", K)
  }

  ok <- stats::complete.cases(M)
  if (sum(ok) < 2L) {
    stop("simulation_se(): only ", sum(ok), " of ", K, " refits succeeded, ",
      "which is not enough to form a standard deviation. If the data cannot ",
      "be recovered from the stored call, pass `envir`.",
      call. = FALSE
    )
  }
  M <- M[ok, , drop = FALSE]

  sim <- apply(M, 2L, stats::sd)
  samp <- object$std.errors
  if (is.null(samp) || length(samp) != length(base)) samp <- rep(NA_real_, length(base))

  out <- data.frame(
    parameter = nm,
    estimate = as.numeric(base),
    sampling_se = as.numeric(samp),
    simulation_se = as.numeric(sim),
    ## The number a reader actually wants: what share of the reported standard
    ## error is an artefact of having simulated the integral rather than a
    ## statement about the population. Large values mean "raise Nsim", not
    ## "collect more data".
    ratio = as.numeric(sim / samp),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  attr(out, "n_refits") <- nrow(M)
  attr(out, "seed_arg") <- arg
  out
}
