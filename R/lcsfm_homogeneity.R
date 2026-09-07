## lcsfm_homogeneity(): how many latent classes are there? Gap L3, from Stead,
## Wheat and Greene (2023), J Prod Anal 60:37-48.
## See notes/code_history/lcsfm_homogeneity.md and notes/L3_mlrt_design.md.
##
## `lcsfm()` fits J classes and, until now, gave no way to test whether J of
## them are needed -- which is the first question anyone asks of a latent class
## fit. The ordinary likelihood ratio test does not answer it, for three
## reasons that stack:
##
##   1. Under homogeneity the second class's parameters are UNIDENTIFIED: if
##      p = 0 then theta_2 can be anything at all.
##   2. The mixing proportion sits on a BOUNDARY of the parameter space.
##   3. The classes are INTERCHANGEABLE -- swapping theta_1 with theta_2 and p
##      with 1 - p leaves the likelihood unchanged.
##
## So `2 * (logLik_J - logLik_1)` compared against a chi-square on the
## parameter-count difference is wrong twice over, and wrong in the
## anti-conservative direction: it rejects homogeneity too often, which means
## reporting classes that are not there.

## The fix is a PENALISED likelihood. Chen et al. (2001) add
## c * log(J^J * prod_j p_j), which diverges to -Inf as any class probability
## approaches 0 or 1 and so keeps every class occupied; `lcsfm(penalty_c = )`
## maximises that modified likelihood. The penalty is exactly 0 at equal
## probabilities, so the null model is recovered without a handicap.
##
## WHY THE PUBLISHED chi^2_{0:1} NULL IS NOT THE DEFAULT HERE.
##
## Stead, Wheat and Greene report the modified statistic is asymptotically
## chi^2_{0:1} -- a 50:50 mixture of a point mass at zero and chi^2_1. That
## result is stated for the case where a SCALAR parameter differs between
## classes: Zhu and Zhang (2004) require theta_1 and theta_2 scalar, and the
## paper's own application restricts every parameter except the noise scale to
## be common across classes.
##
## `lcsfm()`'s "LCM" does not do that. It lets EVERY parameter vary by class --
## sigma_v, sigma_u and the whole frontier vector, per class -- so the model
## here is not the one the theorem covers, and the limit does not carry over.
##
## This was measured rather than assumed. 200 replications of a one-class DGP
## (n = 400, y = 1 + x + v - |w|), fitted with n_class = 2:
##
##   nominal   10%   ->  actual rejection  82.0%
##   nominal    5%   ->  actual rejection  63.5%
##   nominal    1%   ->  actual rejection  40.5%
##
## and the sampling distribution matches chi^2_5 (KS 0.122) -- 5 being exactly
## the parameter-count difference -- not chi^2_{0:1}, whose median is 0 against
## an observed 4.21. In other words the penalty cannot regularise a problem
## with nine free class-specific parameters the way it regularises one.
##
## So the default null is a PARAMETRIC BOOTSTRAP, which is valid whatever the
## parameter structure: simulate from the fitted one-class model, refit both
## models on each simulated sample, and read the p-value off the resulting
## distribution of the same statistic. `null = "chisq01"` is still available
## for the restricted case the theory covers, and warns.


## Null distribution of the statistic when the theorem does apply: a 50:50
## mixture of a point mass at zero and chi^2_1, so the p-value is HALF the
## one-degree-of-freedom tail and does not depend on how many parameters the
## two models differ by.
.chi2_01_p <- function(x) {
  ifelse(x <= 0, 1, 0.5 * stats::pchisq(x, df = 1, lower.tail = FALSE))
}

## sigma_u and sigma_v from an sfm("NHN") fit, which reports the
## lambda/sigma reparameterization rather than the two scales.
.nhn_scales <- function(fit) {
  p <- fit$out[, "par"]
  lam <- unname(p[["lambda"]])
  sig <- unname(p[["sigma"]])
  su <- lam * sig / sqrt(1 + lam^2)
  sv <- sig / sqrt(1 + lam^2)
  list(sigma_u = su, sigma_v = sv,
       beta = p[setdiff(names(p), c("lambda", "sigma"))])
}

lcsfm_homogeneity <- function(object, null = c("bootstrap", "chisq01"),
                              B = 199, c = 1, level = 0.05,
                              seed = NULL, quiet = FALSE,
                              envir = parent.frame()) {
  null <- match.arg(null)
  if (!inherits(object, "sfareg")) {
    stop("`object` must be an \"sfareg\" fit from lcsfm().", call. = FALSE)
  }
  if (!object$model_name %in% c("LCM", "LCM_CN")) {
    stop("lcsfm_homogeneity(): this test is for model_name \"LCM\" or ",
      "\"LCM_CN\". \"LCM_Z\" drives the class probabilities with ",
      "covariates, which is outside the published result. This fit is ",
      "model_name \"",
      if (is.null(object$model_name)) "unknown" else object$model_name, "\".",
      call. = FALSE
    )
  }
  ## "LCM_CN" IS the restricted specification the chi^2_{0:1} limit is stated
  ## for -- one scalar parameter (the noise scale) differing between
  ## components, everything else common. "LCM" is not.
  .restricted <- identical(object$model_name, "LCM_CN")
  if (length(c) != 1L || !is.numeric(c) || !is.finite(c) || c <= 0) {
    stop("`c` must be a single positive number. Chen et al. (2001) report ",
      "their simulations were insensitive to it; Stead, Wheat and Greene ",
      "(2023) use 1 and 5.",
      call. = FALSE
    )
  }
  if (length(B) != 1L || !is.numeric(B) || !is.finite(B) || B < 19) {
    stop("`B` must be a single number >= 19: the smallest bootstrap p-value ",
      "is 1/(B + 1), so B = 19 is the least that can reach 0.05.",
      call. = FALSE
    )
  }
  B <- as.integer(B)
  J <- object$n_class
  if (is.null(J) || J < 2) {
    stop("lcsfm_homogeneity(): this fit has ", if (is.null(J)) "no" else J,
      " classes, so there is no homogeneity to test against.",
      call. = FALSE
    )
  }
  if (is.null(object$call)) {
    stop("lcsfm_homogeneity(): this fit did not retain its call, so the ",
      "models it must be compared against cannot be built.",
      call. = FALSE
    )
  }

  frm <- object$formula
  if (is.null(frm)) frm <- object$call$formula
  dat <- tryCatch(eval(object$call$data, envir = envir), error = function(e) NULL)
  if (!is.data.frame(dat)) {
    stop("lcsfm_homogeneity(): the data behind this fit could not be ",
      "recovered from its call. Pass `envir`.",
      call. = FALSE
    )
  }
  inefdec <- if (is.null(object$call$inefdec)) TRUE else object$call$inefdec

  ## One evaluation of the statistic on a given data frame: the penalised
  ## J-class maximum against the one-class maximum. Written once so that the
  ## observed value and every bootstrap replicate are produced by identical
  ## code -- a bootstrap whose replicates are computed differently from the
  ## observed statistic tests nothing.
  .mn <- object$model_name
  .stat_on <- function(d) {
    fJ <- tryCatch(suppressWarnings(
      lcsfm(frm, model_name = .mn, data = d, n_class = J,
        inefdec = inefdec, penalty_c = c)
    ), error = function(e) NULL)
    f1 <- tryCatch(suppressWarnings(
      sfm(frm, model_name = "NHN", data = d, inefdec = inefdec)
    ), error = function(e) NULL)
    if (is.null(fJ) || is.null(f1)) return(list(stat = NA_real_))
    list(stat = 2 * (-fJ$opt$value - (-f1$opt$value)),
         fJ = fJ, f1 = f1)
  }

  obs <- .stat_on(dat)
  if (!is.finite(obs$stat)) {
    stop("lcsfm_homogeneity(): the observed statistic could not be computed ",
      "because one of the two fits failed on the original data.",
      call. = FALSE
    )
  }
  stat <- obs$stat
  fit_1 <- obs$f1
  fit_J <- obs$fJ

  boot <- NULL
  n_ok <- NA_integer_
  if (identical(null, "chisq01") && .restricted) {
    ## "LCM_CN" IS the specification the limit is stated for, and restricting
    ## to it takes the measured size from 63.5% to 8.0% at a nominal 5% -- so
    ## no warning. It is still mildly LIBERAL in finite samples (150
    ## replications, n = 400: 13.3% at 10%, 8.0% at 5%, 3.3% at 1%), partly
    ## because the two-component optimizer lands below the one-component fit
    ## on about half of null samples, piling 65% of the mass at exactly zero
    ## against the 50% the mixture predicts. Reported through `size_note` so a
    ## reader is not left thinking the p-value is exact.
    p <- .chi2_01_p(stat)
  } else if (identical(null, "chisq01")) {
    warning("lcsfm_homogeneity(null = \"chisq01\"): the chi^2_{0:1} limit is ",
      "stated for a SCALAR parameter differing between classes, and ",
      "lcsfm()'s \"LCM\" lets every parameter vary by class. Measured on a ",
      "one-class DGP at n = 400, this null rejects 63.5% of the time at a ",
      "nominal 5%. Use null = \"bootstrap\" unless you have restricted the ",
      "model to the case the theory covers.",
      call. = FALSE
    )
    p <- .chi2_01_p(stat)
  } else {
    ## Parametric bootstrap under the fitted ONE-class model, which is the
    ## null. Valid whatever the parameter structure, because it never appeals
    ## to an asymptotic distribution: it simulates the null the model actually
    ## implies and re-runs the whole procedure.
    if (!is.null(seed)) {
      .st <- .rng_snapshot()
      on.exit(.rng_restore(.st), add = TRUE)
      set.seed(seed)
    }
    sc <- .nhn_scales(fit_1)
    mf <- stats::model.frame(stats::as.formula(frm), data = dat)
    X <- stats::model.matrix(stats::as.formula(frm), data = mf)
    S <- if (isTRUE(inefdec)) 1 else -1
    n <- nrow(X)
    yname <- all.vars(stats::as.formula(frm))[1L]

    bs <- rep(NA_real_, B)
    for (b in seq_len(B)) {
      db <- dat
      db[[yname]] <- as.numeric(X %*% sc$beta) +
        stats::rnorm(n, 0, sc$sigma_v) -
        S * abs(stats::rnorm(n, 0, sc$sigma_u))
      bs[b] <- .stat_on(db)$stat
      if (!quiet && b %% 25 == 0) {
        message("lcsfm_homogeneity: bootstrap ", b, " of ", B)
      }
    }
    boot <- bs[is.finite(bs)]
    n_ok <- length(boot)
    if (n_ok < 19L) {
      stop("lcsfm_homogeneity(): only ", n_ok, " of ", B, " bootstrap ",
        "replicates produced a usable statistic, which is too few for a ",
        "p-value.",
        call. = FALSE
      )
    }
    ## The +1 in both places is the standard finite-sample correction: it keeps
    ## the p-value away from an exact zero, which a bootstrap can never justify.
    p <- (1 + sum(boot >= stat)) / (n_ok + 1)
  }

  out <- list(
    statistic = c("MLR" = stat),
    p.value = p,
    parameter = c("classes" = J),
    method = paste0(
      "Modified likelihood ratio test for homogeneity in a latent class\n",
      "stochastic frontier (Stead, Wheat & Greene 2023); null by ",
      if (identical(null, "bootstrap")) {
        paste0("parametric bootstrap (B = ", B, ")")
      } else {
        "chi^2_{0:1}"
      }
    ),
    data.name = paste0(deparse(frm), collapse = ""),
    alternative = paste0(J, " latent classes rather than one"),
    null = null, B = if (identical(null, "bootstrap")) B else NA_integer_,
    boot_ok = n_ok, boot = boot,
    penalty_c = c,
    logLik_penalised_J = -fit_J$opt$value,
    logLik_unpenalised_J = fit_J$logLik_unpenalised,
    penalty_at_max = fit_J$penalty,
    logLik_1 = -fit_1$opt$value,
    n_class = J,
    class_prob = fit_J$class_prob,
    level = level,
    reject = isTRUE(p < level),
    provisional = identical(null, "chisq01") && !.restricted,
    restricted = .restricted,
    size_note = if (identical(null, "chisq01") && .restricted) {
      paste0("measured size at n = 400: 13.3% at nominal 10%, 8.0% at 5%, ",
        "3.3% at 1% -- mildly liberal")
    } else {
      NULL
    }
  )
  class(out) <- c("sfa_mlrt", "htest")
  out
}

print.sfa_mlrt <- function(x, ...) {
  cat("\n", x$method, "\n\n", sep = "")
  cat("data:  ", x$data.name, "\n", sep = "")
  cat(sprintf("MLR = %.4f, p-value = %.4g\n", x$statistic, x$p.value))
  cat("alternative hypothesis: ", x$alternative, "\n\n", sep = "")
  cat(sprintf("  penalised log-likelihood (%d classes) : %.4f\n",
    x$n_class, x$logLik_penalised_J
  ))
  if (!is.null(x$penalty_at_max)) {
    cat(sprintf("    of which penalty (c = %g)           : %.4f\n",
      x$penalty_c, x$penalty_at_max
    ))
  }
  cat(sprintf("  log-likelihood (1 class, sfm \"NHN\")   : %.4f\n", x$logLik_1))
  if (!is.null(x$class_prob)) {
    cat("  fitted class probabilities            : ",
      paste(sprintf("%.3f", x$class_prob), collapse = ", "), "\n",
      sep = ""
    )
  }
  if (identical(x$null, "bootstrap")) {
    cat(sprintf("  usable bootstrap replicates           : %d of %d\n",
      x$boot_ok, x$B
    ))
  }
  cat("\n")
  if (x$statistic < 0) {
    cat("  WARNING: the statistic is NEGATIVE, so the ", x$n_class,
      "-class fit landed\n  BELOW the one-class fit. That is an optimizer ",
      "failure, not evidence\n  for homogeneity. Refit before reading the ",
      "p-value.\n\n",
      sep = ""
    )
  } else if (x$reject) {
    cat("  Reject homogeneity at the ", 100 * x$level,
      "% level: more than one class.\n", sep = ""
    )
  } else {
    cat("  Do not reject homogeneity at the ", 100 * x$level,
      "% level.\n", sep = ""
    )
  }
  if (!is.null(x$size_note)) {
    cat("\n  NOTE: \"LCM_CN\" is the restricted case the chi^2_{0:1} limit is\n")
    cat("  stated for, so it applies here -- but it is not exact in finite\n")
    cat("  samples. ", x$size_note, ".\n", sep = "")
    cat("  Use null = \"bootstrap\" if the p-value is close to your cutoff.\n")
  }
  if (isTRUE(x$provisional)) {
    cat("\n  NOTE: the chi^2_{0:1} limit assumes a SCALAR parameter differs\n")
    cat("  between classes. \"LCM\" lets every parameter vary, and on a\n")
    cat("  one-class DGP this null rejected 63.5% of the time at a nominal\n")
    cat("  5%. Prefer null = \"bootstrap\".\n")
  }
  cat("\n")
  invisible(x)
}
