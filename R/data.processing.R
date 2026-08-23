# Package-level constants
.SFA_CONSTANTS <- list(
  CLIP_Z1_UPPER  =  37,
  CLIP_Z1_LOWER  = -37,
  CLIP_Z2_UPPER  =   8,
  CLIP_Z2_LOWER  = -37,
  HALTON_DISCARD =  1000,
  MIN_POSITIVE   = .Machine$double.eps,
  MAX_VALUE      = .Machine$double.xmax,
  # Safe upper bound for arguments passed to exp() -- exp(700) ~= 1.01e304,
  # comfortably under .Machine$double.xmax (~1.80e308), so clipping an
  # exponent argument here prevents exp() from overflowing to Inf while
  # still preserving its relative ordering for any argument that wouldn't
  # have overflowed anyway. Added to fix a real "non-finite value supplied
  # by optim" error reported from ttsfm(model_name="TTNE") -- see ttsfm.R.
  EXP_CLIP_UPPER =  700
  # HALTON_PRIMES = c(2, 3)
)

## Second data cleaning for missing values
data_proc2 <- function(data, data_x, fancy_vars, fancy_vars_z, data_z, y_var,
                       x_vars_vec, halton_num, individual, N, model_name, rand.gtre) {
  data <- data[rownames(data_x), ]

  if (isTRUE(length(fancy_vars) > 0)) {
    data_inter <- cbind(data, data_x[, fancy_vars])
    colnames(data_inter) <- c(colnames(data), fancy_vars)
    data <- data_inter
  } else {}

  if (isTRUE(length(fancy_vars_z) > 0)) {
    data_inter <- cbind(data, data_z[, fancy_vars_z])
    colnames(data_inter) <- c(colnames(data), fancy_vars_z)
    data <- data_inter
  } else {}

  Y <- as.matrix(subset(data, select = y_var))
  data_i_vars <- as.matrix(data.frame(subset(data, select = x_vars_vec)))

  if (model_name %in% c("GTRE", "TRE")) {
    if (isTRUE(is.numeric(halton_num))) {
      R <- halton_num
    } else {
      R <- ceiling(sqrt(nrow(data))) + 100
    } ## Integral reps

    R_H <- randtoolbox::halton(R + 1000, 2, start = 1, normal = FALSE)[-c(1:1000), c(1:2)]
    R_H <- cbind(qnorm(R_H[, 1]), sqrt(2) * pracma::erfinv(R_H[, 2])) ## using inverse error function for R_H2

    if (!is.null(rand.gtre)) {
      .rng_state <- .rng_snapshot()
      on.exit(.rng_restore(.rng_state), add = TRUE)
      set.seed(rand.gtre)
    }

    mat <- matrix(0, nrow = R, ncol = 9999)
    for (v in 1:9999) {
      mat[, v] <- sample(R_H[, 1])
    }

    cor_vec <- abs(cor(mat, R_H[, 2]))

    R_H <- cbind(mat[, which.min(cor_vec)], R_H[, 2])
    rm(cor_vec, v, mat)

    # print(paste( "Primes 2 and 3 are in use, with 1,000 discards.  Correlation between R and H draws is:", round(cor(R_H)[1,2],10), sep = "" ))

    indiv <- noquote(as.vector(unique(data[, c(individual)])))
    t <- rep(0, N)
    data_i <- Y <- eps <- data_i_vars <- R_h1 <- R_h2 <- as.list(rep(0, N))

    for (ii in seq_len(N)) {
      data_i[[ii]] <- data[which(data[, c(individual)] == indiv[ii]), ]
      t[ii] <- nrow(data_i[[ii]])
      R_h1[[ii]] <- t(matrix(rep(R_H[, 1], t[[ii]]), R, t[[ii]]))
      R_h2[[ii]] <- abs(t(matrix(rep(R_H[, 2], t[[ii]]), R, t[[ii]])))
      Y[[ii]] <- matrix(rep(data_i[[ii]][, y_var], R), t[[ii]], R)
      data_i_vars[[ii]] <- as.matrix(data_i[[ii]][, c(x_vars_vec), drop = FALSE])
    }

    result_gtre_tre <- list(R, R_H, indiv, t, data_i, eps, R_h1, R_h2, data_x)
    names(result_gtre_tre) <- c("R", "R_H", "indiv", "t", "data_i", "eps", "R_h1", "R_h2", "data_x")
  }

  result <- list(data, Y, data_i_vars)
  names(result) <- c("data", "Y", "data_i_vars")

  if (model_name %in% c("GTRE", "TRE")) {
    result <- c(result, result_gtre_tre)
  }

  return(result)
}

## Inital Data cleaning and Processing
data_proc <- function(formula, data, model_name, individual = NULL, inefdec) {
  if (model_name == "GTRE_Z") {
    formula <- .format_formula(formula)
  }

  data_z <- formula_z <- z_vars <- intercept_z <- z_vars_vec <- z_z_vec <- n_z_vars <- formula_zp <- data_zp <- zp_vars <-
    intercept_zp <- zp_vars_vec <- n_zp_vars <- zp_zp_vec <- fancy_vars_zp <- NULL
  n_z_vars <- NA ## maybe make condition here making this default to null
  data_orig <- data

  ## Formula parsing via .parse_pipe_formula() (matrix_utils.R), which uses
  ## Formula::Formula() -- see that function's header comment for why this
  ## replaced a manual strsplit(as.character(formula),"|") + paste()
  ## reconstruction that used to live here (it silently misparsed long
  ## formulas). `form_parts` is kept as a synthetic placeholder purely so
  ## the `length(unlist(form_parts)) > 3` / `> 4` part-count checks used by
  ## this function's callers (sfm.R, zsfm.R, ttsfm.R, psfm.R) and further
  ## down in this function keep working unchanged -- nothing anywhere reads
  ## its actual content, only its unlisted length (1 RHS part -> length 3,
  ## 2 parts -> length 4, 3 parts -> length 5, matching the old object's
  ## shape).
  parsed_f <- .parse_pipe_formula(formula)
  formula_x <- parsed_f$formula_x
  y_var <- parsed_f$y_var
  form_parts <- list("~", y_var, rep(NA_character_, parsed_f$n_parts))

  if (model_name %in% c("TTNE", "TTHN", "TTNLS")) {
    formula0 <- .format_formula(formula)
    formula <- formula0
    parsed_f <- .parse_pipe_formula(formula)
    formula_x <- parsed_f$formula_x
    y_var <- parsed_f$y_var
    form_parts <- list("~", y_var, rep(NA_character_, parsed_f$n_parts))
  }

  if (parsed_f$n_parts >= 2) {
    formula_z <- parsed_f$formula_z
    data_z <- data_conform(formula = formula_z, data = data)
    if (model_name == "NHN") {
      model_name <- "NHN_Z"
    }
    if (model_name == "NE") {
      model_name <- "NE_Z"
    }
    if (model_name == "GTRE") {
      model_name <- "GTRE_Z"
      formula0 <- .format_formula(formula) ## redefine formula for GTRE_Z
      formula <- formula0
      parsed_f <- .parse_pipe_formula(formula)
      formula_x <- parsed_f$formula_x
      y_var <- parsed_f$y_var
      formula_z <- parsed_f$formula_z
      data_z <- data_conform(formula = formula_z, data = data)
      form_parts <- list("~", y_var, rep(NA_character_, parsed_f$n_parts))
    }
    if (model_name == "TRE") {
      model_name <- "TRE_Z"
    }
    if (model_name %in% c("THT", "NTN", "CHC", "NU", "tHN")) {
      return(c("Currently building this functionality"))
    } # else{z_vars <- NULL}
  }

  if (parsed_f$n_parts >= 3) {
    formula_zp <- parsed_f$formula_zp
    data_zp <- data_conform(formula = formula_zp, data = data)
  }


  data_x <- data_conform(formula = formula_x, data = data)

  intercept <- attr(terms(formula_x), "intercept")
  inefdec_n <- if (isTRUE(inefdec)) {
    1
  } else {
    -1
  }
  inefdec_TF <- if (isTRUE(inefdec)) {
    TRUE
  } else {
    FALSE
  }
  x_vars_vec <- if (model_name %in% c("TFE", "TFE_WMLE", "FD", "SSFE") & intercept == 1) {
    colnames(data_x)[-c(1)]
  } else {
    colnames(data_x)
  }
  n_x_vars <- length(x_vars_vec)
  x_vars <- x_vars_vec
  x_x_vec <- rep(0, length = n_x_vars)
  fancy_vars <- setdiff(colnames(data_x), colnames(data))
  fancy_vars_z <- NULL

  N <- NULL
  if (is.null(individual)) {} else {
    N <- length(unique(data[, c(individual)]))
  }

  if (length(unlist(form_parts)) > 3) {
    intercept_z <- attr(terms(formula_z), "intercept")
    z_vars_vec <- if (model_name %in% c("TFE", "TFE_WMLE", "FD") & intercept_z == 1) {
      colnames(data_z)[-c(1)]
    } else {
      colnames(data_z)
    }
    n_z_vars <- length(z_vars_vec)
    z_vars <- z_vars_vec
    z_z_vec <- rep(0, length = n_z_vars)
    fancy_vars_z <- setdiff(colnames(data_z), colnames(data))
  }

  if (length(unlist(form_parts)) > 4) {
    intercept_zp <- attr(terms(formula_zp), "intercept")
    zp_vars_vec <- if (model_name %in% c("TFE", "TFE_WMLE", "FD") & intercept_zp == 1) {
      colnames(data_zp)[-c(1)]
    } else {
      colnames(data_zp)
    }
    n_zp_vars <- length(zp_vars_vec)
    zp_vars <- zp_vars_vec
    zp_zp_vec <- rep(0, length = n_zp_vars)
    fancy_vars_zp <- setdiff(colnames(data_zp), colnames(data))
  }


  # 1. Identify rows that are complete across all variables in the original data
  # Find all column names present in your generated dataframes
  all_cols <- unique(c(
    colnames(data_x),
    if (!is.null(data_z)) colnames(data_z),
    if (!is.null(data_zp)) colnames(data_zp)
  ))

  # Filter out intercepts or calculated columns not found in the original data
  valid_orig_cols <- intersect(all_cols, colnames(data))

  # Get the logical vector of complete cases based on the original data rows
  is_complete <- complete.cases(data[, valid_orig_cols, drop = FALSE])

  # 2. Re-run data_conform using ONLY the complete rows from the original data
  data_x <- data_conform(formula = formula_x, data = data[is_complete, , drop = FALSE])

  if (!is.null(data_z)) {
    data_z <- data_conform(formula = formula_z, data = data[is_complete, , drop = FALSE])
  }
  if (!is.null(data_zp)) {
    data_zp <- data_conform(formula = formula_zp, data = data[is_complete, , drop = FALSE])
  }

  result <- list(
    data_orig, form_parts, formula_x, y_var, model_name, data_x, intercept, inefdec_n, inefdec_TF,
    x_vars_vec, n_x_vars, x_vars, x_x_vec, fancy_vars, fancy_vars_z, n_z_vars, N,
    formula_z, z_vars, data_z, intercept_z, z_vars_vec, z_z_vec, n_z_vars, formula_zp,
    data_zp, zp_vars, intercept_zp, zp_vars_vec, n_zp_vars, zp_zp_vec, fancy_vars_zp, formula
  )

  names(result) <- c(
    "data_orig", "form_parts", "formula_x", "y_var", "model_name", "data_x", "intercept", "inefdec_n",
    "inefdec_TF", "x_vars_vec", "n_x_vars", "x_vars", "x_x_vec", "fancy_vars", "fancy_vars_z", "n_z_vars", "N",
    "formula_z", "z_vars", "data_z", "intercept_z", "z_vars_vec", "z_z_vec", "n_z_vars", "formula_zp",
    "data_zp", "zp_vars", "intercept_zp", "zp_vars_vec", "n_zp_vars", "zp_zp_vec", "fancy_vars_zp", "formula"
  )
  return(result)
}

## lm code for data
data_conform <- function(formula, data, na.action, method = "qr",
                         model = TRUE, x = FALSE, y = FALSE, qr = TRUE, singular.ok = TRUE,
                         contrasts = NULL, offset, ...) {
  ret.x <- x
  ret.y <- y
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c(
    "formula", "data", "subset", "weights", "na.action",
    "offset"
  ), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  if (method == "model.frame") {
    return(mf)
  } else if (method != "qr") {
    warning(gettextf("method = '%s' is not supported. Using 'qr'", method), domain = NA, call. = FALSE)
  }
  mt <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  x <- model.matrix(mt, mf, contrasts)
  return(x)
}

## Basic PCS code for regression on intercept
pcs_c <- function(Y, inefdec = TRUE, Method = formals(psfm)$Method) {
  sigma_u <- runif(n = 1, min = 0.04, max = 1) ## Random starting values
  sigma_v <- runif(n = 1, min = 0.04, max = 1)
  beta_0 <- runif(n = 1, min = 0.04, max = 1)

  lambda <- sigma_u / sigma_v
  sigma <- sqrt(sigma_u^2 + sigma_v^2)

  start_v <- unname(c(lambda, sigma, beta_0))
  inefdec_n <- if (isTRUE(inefdec)) {
    1
  } else {
    -1
  }

  fn <- function(x) {
    eps <- inefdec_n * (Y - x[3])
    like <- sum(log((2 / x[2]) *
      pmax(dnorm(eps / x[2]), eps * 0 + .Machine$double.eps) *
      pmax(pnorm(-eps * x[1] / x[2]), eps * 0 + .Machine$double.eps)))
    return(-like[is.finite(like)])
  }
  opt <- optim(
    par = start_v,
    fn = fn,
    lower = c(rep(0.01, 2), -Inf),
    control = list(maxit = 300, REPORT = 1, trace = 0, pgtol = 0),
    hessian = FALSE,
    method = Method
  )
  return(list(opt, Y))
}


summary.sfareg <- function(object, ...) {
  # Perform specific summary calculations and formatting
  cat("--- SFA Regression Model Summary ---\n")
  cat("Formula:", deparse(object$formula), "\n")
  cat("Total time:", object$total_time, "\n")
  cat("Model Output:\n")
  print(object$out)
  if (is.numeric(object$opt$value)) {
    cat("log likelihood:", -object$opt$value, "\n")
  }
  .sfa_report_convergence(object)
  ## return original object
  invisible(object)
}

print.sfareg <- function(x, ...) {
  # Perform specific summary calculations and formatting
  cat("--- SFA Regression Model Summary ---\n")
  cat("Formula:", deparse(x$formula), "\n")
  cat("Total time:", x$total_time, "\n")
  cat("Model Output:\n")

  print(x$out)
  if (is.numeric(x$opt$value)) {
    cat("log likelihood:", -x$opt$value, "\n")
  }
  .sfa_report_convergence(x)
}

## One line on how the optimizer finished, shown by both print() and
## summary(). Previously neither said anything at all: a fit that stopped on
## the iteration cap looked exactly like a converged one.
.sfa_report_convergence <- function(x) {
  cc <- x$opt$convergence
  if (is.null(cc) || !length(cc) || is.na(cc)) {
    return(invisible(NULL))
  }
  if (identical(as.integer(cc), 0L)) {
    cat("convergence: 0 (converged)\n")
  } else {
    cat("convergence: ", cc, " -- ",
      switch(as.character(cc),
        "1" = "ITERATION LIMIT REACHED; this is not a converged optimum",
        "10" = "Nelder-Mead simplex degenerated",
        "51" = "L-BFGS-B warning",
        "52" = "L-BFGS-B error",
        "see ?optim for this code"
      ), "\n",
      sep = ""
    )
    if (!is.null(x$opt$message)) cat("  optimizer message: ", x$opt$message, "\n", sep = "")
    ## A non-zero code frequently means only that the final polish stage could
    ## not improve on an already-converged point; the gradient and Hessian are
    ## what settle it, and sfa_diagnostics() combines all three.
    cat("  a non-zero code does not by itself mean the fit failed --\n")
    cat("  run sfa_diagnostics() on this fit to see the gradient and Hessian.\n")
  }
  invisible(NULL)
}
