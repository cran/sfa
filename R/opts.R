## Convergence code for "the final optim() polish stage did not run or did not
## produce a usable result". optim() itself returns 0, 1, 10, 51 or 52, so 99
## cannot collide with one of its own codes.
.SFA_OPTIM_SKIPPED <- 99L

## opt.nlminb() -- primary optimization stage.
opt.nlminb <- function(fn, start_v, lower.nlminb, upper.nlminb = Inf,
                       gr = NULL, maxit.nlminb = 500, nlminb.TF = TRUE,
                       verbose = FALSE) {
  start_feval <- fn(start_v)
  nlm1 <- NULL
  if (isTRUE(nlminb.TF)) {
    nlm1 <- tryCatch(
      stats::nlminb(
        start = start_v, objective = fn, gradient = gr,
        lower = lower.nlminb, upper = upper.nlminb,
        control = list(
          iter.max = maxit.nlminb,
          eval.max = 2L * maxit.nlminb,
          trace = if (verbose) 1L else 0L
        )
      ),
      error = function(e) NULL
    )
    ## Only accept the new point if it actually improved the objective.
    if (!is.null(nlm1) && is.finite(nlm1$objective) && isTRUE(start_feval > nlm1$objective)) {
      start_v <- nlm1$par
      start_feval <- nlm1$objective
    }
  }
  results <- list(start_v, start_feval, nlm1)
  names(results) <- c("start_v", "start_feval", "nlm1")
  return(results)
}

opt.bobyqa <- function(fn, start_v, lower.bobyqa, upper.bobyqa = Inf, maxit.bobyqa, bob.TF, rhobeg = NA, rhoend = NA, verbose = verbose) {
  start_feval <- fn(start_v)
  bob1 <- NULL
  if (isTRUE(bob.TF == TRUE)) {
    bob1 <- bobyqa(
      par = start_v,
      fn = fn,
      lower = lower.bobyqa,
      upper = upper.bobyqa,
      control = list(
        iprint = if (verbose) 2 else 0,
        maxfun = maxit.bobyqa,
        rhobeg = rhobeg,
        rhoend = rhoend
      )
    )

    if (isTRUE(start_feval > bob1$fval)) {
      start_v <- bob1$par
      start_feval <- fn(start_v)
    }
  }

  results <- list(start_v, start_feval, bob1)
  names(results) <- c("start_v", "start_feval", "bob1")
  return(results)
}

## Stage 3 is a refinement pass over a point stages 1-2 already produced, so a
## failure here should cost the polish, not the fit. See
## notes/code_history/opts.md for the two failures this guards.
##
## NOT here, deliberately: a guard that falls back to the stage-2 point whenever
## optim() ends at a higher objective than it started from. It looks obviously
## right -- stages 1 and 2 refuse to move backwards, so why not stage 3 -- and it
## is wrong, because "the point it started from" is not always trustworthy. On
## sfma()'s own test data the NR likelihood has a spike: the stage-2 point
## evaluates to -3.6e17, optim() correctly escapes it to 541.8, and such a guard
## drags the fit back onto the spike and hands sfma() a log-likelihood of
## +3.6e17. L-BFGS-B is a descent method; when it ends above its own start, that
## is evidence about the START, not about the search. See
## test-ttsfm-numerics.R, which pins this.
opt.optim <- function(fn, start_v, lower.optim, upper.optim, maxit.optim, opt.TF, method, optHessian, trace, verbose = verbose) {
  start_feval <- fn(start_v)
  opt <- NULL
  if (isTRUE(opt.TF == TRUE)) {
    np <- length(start_v)

    run_optim <- function(use_hessian) {
      tryCatch(
        optim(
          par = start_v,
          fn = fn,
          lower = lower.optim,
          upper = upper.optim,
          hessian = use_hessian,
          method = method,
          control = list(
            maxit = maxit.optim,
            REPORT = base::ceiling(maxit.optim / 10),
            trace = if (verbose) {
              1
            } else {
              0
            }
          )
        ),
        error = function(e) NULL
      )
    }

    opt <- run_optim(optHessian)

    ## optim() threw, or stopped at a point the objective cannot evaluate.
    if (!is.null(opt) && (!is.numeric(opt$value) || !is.finite(opt$value))) {
      opt <- NULL
    }

    ## The point is usable but the numerical Hessian is not. Keep the point and
    ## hand downstream an all-NA Hessian, which every SE path already degrades
    ## to NA on.
    if (!is.null(opt) && isTRUE(optHessian) &&
      (is.null(opt$hessian) || !all(is.finite(opt$hessian)))) {
      opt$hessian <- matrix(NA_real_, np, np)
      opt$convergence <- .SFA_OPTIM_SKIPPED
      opt$message <- "final optim() stage returned a non-finite Hessian; standard errors are unavailable"
    }

    ## Nothing usable came back: return an optim-shaped result at the stage-2
    ## point so callers that index opt$par / opt$hessian keep working.
    if (is.null(opt)) {
      opt <- list(
        par = start_v,
        value = start_feval,
        counts = c(`function` = NA_integer_, gradient = NA_integer_),
        convergence = .SFA_OPTIM_SKIPPED,
        message = "final optim() stage skipped: non-finite objective or optim() error; the stage-2 result is returned",
        hessian = if (isTRUE(optHessian)) matrix(NA_real_, np, np) else NULL
      )
    }

    ## convergence is visible in print()/summary(), but coef() and friends do
    ## not go through those, so say it once out loud as well.
    if (identical(opt$convergence, .SFA_OPTIM_SKIPPED)) {
      warning(
        "The final optim() stage did not produce a usable result (", opt$message,
        "). Estimates come from the preceding stage and standard errors are NA. ",
        "Try optHessian = FALSE, vcov(type = \"bhhh\"), or different starting values.",
        call. = FALSE
      )
    }

    if (isTRUE(start_feval > opt$value)) {
      start_v <- opt$par
      start_feval <- fn(start_v)
    }
  }

  results <- list(start_v, start_feval, opt)
  names(results) <- c("start_v", "start_feval", "opt")
  return(results)
}

opt.psoptim <- function(fn, start_v, lower.psoptim, upper.psoptim = NA, maxit.psoptim = maxit.psoptim,
                        psopt.TF, rand.order = TRUE, verbose = verbose, rand.psoptim = rand.psoptim) {
  start_feval <- fn(start_v)
  opt00 <- NULL
  report <- base::ifelse(verbose, base::ceiling(maxit.psoptim / 10), 0)
  if (isTRUE(psopt.TF == TRUE)) {
    if (!is.null(rand.psoptim)) {
      .rng_state <- .rng_snapshot()
      on.exit(.rng_restore(.rng_state), add = TRUE)
      set.seed(rand.psoptim)
    }

    opt00 <- psoptim(
      par = start_v,
      fn = fn,
      lower = lower.psoptim,
      upper = upper.psoptim,
      control = list(
        trace = if (verbose) {
          1
        } else {
          0
        },
        REPORT = report,
        trace.stats = if (verbose) {
          TRUE
        } else {
          FALSE
        },
        maxit = maxit.psoptim,
        rand.order = rand.order
      )
    )

    if (isTRUE(start_feval > opt00$value)) {
      start_v <- opt00$par
      start_feval <- fn(start_v)
    }
  }

  results <- list(start_v, start_feval, opt00)
  names(results) <- c("start_v", "start_feval", "opt00")
  return(results)
}
