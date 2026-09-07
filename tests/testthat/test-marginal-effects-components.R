## The persistent block of psfm("GTRE_Z").
##
## GTRE_Z separates persistent inefficiency h_i from transient u_it, each with
## its own determinant block. Reporting effects on only one of them would be
## half the model, and confusing the two would be worse: the sigma_h
## coefficients are the TRAILING ones, so a positional implementation returns
## them under a sigma_u label.

test_that("GTRE_Z exposes both blocks and keeps them apart", {
  skip_on_cran()
  d <- as.data.frame(data_gen_p(t = 8, N = 80, rand = 1, sig_u = 1, sig_v = 0.3,
                                sig_r = 0.2, sig_h = 0.4, cons = 0.5,
                                beta1 = 0.5, beta2 = 0.5, eta = 0.1))
  f <- suppressWarnings(psfm(y_gtre_zz ~ x1 + x2 | z_gtre | zp_gtre,
                             model_name = "GTRE_Z", data = d,
                             individual = "name", halton_num = 40))
  cf <- coef(f)
  expect_equal(f$z_spec$delta,   unname(cf[c("(Intercept u)", "z_gtre")]),  tolerance = 1e-12)
  expect_equal(f$z_spec_h$delta, unname(cf[c("(Intercept h)", "zp_gtre")]), tolerance = 1e-12)
  ## The two blocks are different, so a mix-up would be visible.
  expect_false(isTRUE(all.equal(f$z_spec$delta, f$z_spec_h$delta)))

  mu <- marginal_effects(f)
  mh <- marginal_effects(f, component = "h")
  ## Column names carry the component, so a table cannot be misread once it is
  ## separated from the call that produced it.
  expect_true(all(c("sigma_u", "E_u", "Var_u", "dE_u.dz_gtre") %in% names(mu)))
  expect_true(all(c("sigma_h", "E_h", "Var_h", "dE_h.dzp_gtre") %in% names(mh)))
  expect_identical(attr(mu, "component"), "u")
  expect_identical(attr(mh, "component"), "h")
  ## Both are real, finite, positive scales.
  expect_true(all(is.finite(mu$sigma_u)) && all(mu$sigma_u > 0))
  expect_true(all(is.finite(mh$sigma_h)) && all(mh$sigma_h > 0))
})

test_that("component = \"h\" is refused where there is no persistent block", {
  skip_on_cran()
  d <- as.data.frame(data_gen_p(t = 10, N = 100, rand = 1, sig_u = 1, sig_v = 0.3,
                                sig_r = 0.2, sig_h = 0.4, cons = 0.5,
                                beta1 = 0.5, beta2 = 0.5, eta = 0.1))
  f <- suppressWarnings(psfm(y_tre_z ~ x1 + x2 | z_gtre, model_name = "TRE_Z",
                             data = d, individual = "name"))
  expect_null(f$z_spec_h)
  expect_error(marginal_effects(f, component = "h"), "GTRE_Z")
  ## The default is unchanged.
  expect_identical(attr(marginal_effects(f), "component"), "u")
})

test_that("the u-component output is unchanged for cross-sectional fits", {
  skip_on_cran()
  d <- cs_small(N = 400)
  d$z <- runif(nrow(d), -1, 1)
  f <- suppressWarnings(sfm(y_pcs ~ x1 + x2 | z, model_name = "NHN_Z", data = d))
  me <- marginal_effects(f)
  expect_true(all(c("sigma_u", "E_u", "Var_u", "dE_u.dz", "dVar_u.dz") %in% names(me)))
  expect_error(marginal_effects(f, component = "h"), "GTRE_Z")
})
