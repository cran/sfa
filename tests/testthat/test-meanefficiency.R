## meanefficiency(): model-implied E[exp(-U)] and supra-percentile means (G6).

me_data <- function() {
  data_gen_cs(N = 400, rand = 1, sig_u = 1, sig_v = 0.3, cons = 0.5,
    beta1 = 0.5, beta2 = 0.5, a = 1, mu = 0.5
  )
}
me_cols <- c(NHN = "y_pcs", NE = "y_pcs_e", NTN = "y_pcs_tn", NR = "y_pcs_r",
             NG = "y_pcs_g", NU = "y_pcs_u", NGE = "y_pcs_ge",
             NLN = "y_pcs_ln", NW = "y_pcs_wb", NNAK = "y_pcs_nak")

test_that("every closed form agrees with integrating the same density", {
  skip_on_cran()
  d <- me_data()
  ## This is the reason the numeric engine is kept even where a closed form
  ## exists. An algebraic slip in a transform cannot survive being compared
  ## against the integral of the density it was derived from.
  for (m in names(me_cols)) {
    f <- suppressWarnings(
      sfm(stats::as.formula(paste0(me_cols[[m]], "~ x1 + x2")),
        model_name = m, data = d
      )
    )
    a <- meanefficiency(f, use_closed_form = TRUE)$mean_efficiency
    b <- meanefficiency(f, use_closed_form = FALSE)$mean_efficiency
    expect_equal(a, b, tolerance = 1e-6, info = m)
    expect_true(a > 0 && a < 1, info = m)
  }
})

test_that("the density and the quantile function describe the same distribution", {
  skip_on_cran()
  ## A shared parameter-extraction error would make the closed form and the
  ## integral agree with each other and both be wrong. Drawing through the
  ## QUANTILE function and averaging exp(-u) is independent of the density, so
  ## it catches that.
  d <- me_data()
  set.seed(99)
  U <- runif(2e5)
  for (m in names(me_cols)) {
    f <- suppressWarnings(
      sfm(stats::as.formula(paste0(me_cols[[m]], "~ x1 + x2")),
        model_name = m, data = d
      )
    )
    r <- meanefficiency(f)
    mc <- mean(exp(-.u_spec(f)$qf(U)))
    expect_equal(r$mean_efficiency, mc, tolerance = 0.01, info = m)
  }
})

test_that("supra-percentile means decrease in p and exceed the overall mean", {
  skip_on_cran()
  d <- me_data()
  f <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = d)
  r <- meanefficiency(f, p = c(0.05, 0.10, 0.25, 0.50, 0.90))

  ## Conditioning on the most efficient p: a smaller p is a more select group,
  ## so mean efficiency must fall as p rises, and every one of them must sit
  ## above the unconditional mean.
  expect_true(all(diff(r$supra$mean_efficiency) < 0))
  expect_true(all(r$supra$mean_efficiency > r$mean_efficiency))
  expect_true(all(r$supra$mean_efficiency <= 1))
  ## p = 1 is the whole distribution, so it IS the unconditional mean.
  expect_equal(meanefficiency(f, p = 1)$supra$mean_efficiency,
    r$mean_efficiency,
    tolerance = 1e-8
  )
})

test_that("it reports what it did, and refuses what it cannot do", {
  skip_on_cran()
  d <- me_data()
  f <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = d)
  r <- meanefficiency(f)
  expect_identical(r$distribution, "half-normal")
  expect_identical(r$method, "closed form")
  expect_identical(meanefficiency(f, use_closed_form = FALSE)$method,
    "numerical integration"
  )

  ## Lognormal has no closed form, and says so rather than inventing one.
  g <- suppressWarnings(sfm(y_pcs_ln ~ x1 + x2, model_name = "NLN", data = d))
  expect_identical(meanefficiency(g)$method, "numerical integration")

  expect_error(meanefficiency(list()), "must be an \"sfareg\" fit")
  expect_error(meanefficiency(f, p = 0), "must be numbers")
  expect_error(meanefficiency(f, p = 1.5), "must be numbers")

  ## A model whose u-distribution is not implemented is named, not guessed at.
  h <- sfm(y_pcs_t ~ x1 + x2, model_name = "THT", data = d)
  expect_error(meanefficiency(h), "no distribution of u is implemented")
})
