## marginal_effects() on panel fits.
##
## psfm()'s _Z models put z'delta on the VARIANCE, sigma_u = sqrt(exp(z'delta)),
## which is the opposite of sfm()'s default -- see sfm()'s z_link and gap C1.
## The block also has to be found by NAME: GTRE_Z carries two determinant
## blocks, "(Intercept u)", z... then "(Intercept h)", zp..., so taking the
## trailing coefficients would hand back the sigma_h coefficients when the user
## asked about sigma_u.

pdat <- function(N = 120, t = 10, seed = 1) {
  as.data.frame(data_gen_p(t = t, N = N, rand = seed, sig_u = 1, sig_v = 0.3,
                           sig_r = 0.2, sig_h = 0.4, cons = 0.5,
                           beta1 = 0.5, beta2 = 0.5, eta = 0.1))
}

test_that("TRE_Z carries a variance-link z_spec and reports marginal effects", {
  skip_on_cran()
  d <- pdat()
  f <- suppressWarnings(psfm(y_tre_z ~ x1 + x2 | z_gtre, model_name = "TRE_Z",
                             data = d, individual = "name"))
  expect_false(is.null(f$z_spec))
  expect_identical(f$z_spec$link, "var")
  ## The delta stored must be the sigma_u block of the fitted vector.
  expect_equal(f$z_spec$delta,
               unname(coef(f)[c("(Intercept u)", "z_gtre")]), tolerance = 1e-12)
  me <- marginal_effects(f)
  expect_true(all(c("sigma_u", "E_u", "Var_u", "dE_u.dz_gtre", "dVar_u.dz_gtre") %in% names(me)))
  expect_equal(nrow(me), nrow(d))
  expect_true(all(is.finite(me$sigma_u)) && all(me$sigma_u > 0))
  ## The DGP has a positive delta on z, so E[u] must rise with z.
  expect_gt(marginal_effects(f, average = TRUE)[["dE_u.dz_gtre"]], 0)
})

test_that("GTRE_Z takes the sigma_u block, not the sigma_h block", {
  skip_on_cran()
  ## This is the regression that positional slicing would cause: the two
  ## blocks are adjacent and the sigma_h one is LAST.
  d <- pdat(N = 80, t = 8)
  f <- suppressWarnings(psfm(y_gtre_zz ~ x1 + x2 | z_gtre | zp_gtre,
                             model_name = "GTRE_Z", data = d,
                             individual = "name", halton_num = 40))
  cf <- coef(f)
  expect_true(all(c("(Intercept u)", "z_gtre", "(Intercept h)", "zp_gtre") %in% names(cf)))
  expect_false(is.null(f$z_spec))
  expect_equal(f$z_spec$delta, unname(cf[c("(Intercept u)", "z_gtre")]), tolerance = 1e-12)
  ## And explicitly NOT the h block.
  expect_false(isTRUE(all.equal(f$z_spec$delta, unname(cf[c("(Intercept h)", "zp_gtre")]))))
})

test_that("a panel model with no z still refuses politely", {
  skip_on_cran()
  d <- pdat(N = 60, t = 8)
  f <- suppressWarnings(psfm(y_tre ~ x1 + x2, model_name = "TRE", data = d,
                             individual = "name"))
  expect_true(is.null(f$z_spec))
  expect_error(marginal_effects(f), "homoskedastic")
})

test_that("the psfm z_spec helper is defensive about its inputs", {
  d <- data.frame(a = 1:5, b = rnorm(5))
  ## No z variables at all.
  expect_null(sfa:::.psfm_z_spec(d, character(0), c(1, 2)))
  ## A z variable that is not a column.
  expect_null(sfa:::.psfm_z_spec(d, c("(Intercept)", "nope"), c(1, 2)))
  ## Anchor absent from the parameter names.
  expect_null(sfa:::.psfm_z_spec(d, c("(Intercept)", "b"), c(x = 1, y = 2)))
  ## A well-formed call finds the block by name.
  par <- c(sigv = 0.3, `(Intercept u)` = 0.4, b = 0.6)
  z <- sfa:::.psfm_z_spec(d, c("(Intercept)", "b"), par)
  expect_equal(z$delta, c(0.4, 0.6))
  expect_identical(z$link, "var")
})
