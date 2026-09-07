## Marginal effects of the variance determinants on E[u] and Var[u].
##
## The analytic expressions are checked against a central finite difference of
## the moment itself, which is the only way to catch an algebra slip that keeps
## the right shape. The variance-link branch is exercised through a synthetic
## z_spec, because no sfm() model uses that link -- psfm()'s do, and the factor
## of two between the two conventions is the whole reason this function exists.

test_that("marginal effects match a finite difference of the moments", {
  skip_on_cran()
  d <- cs_small(N = 300)
  for (mn in c("NHN_Z", "NE_Z")) {
    fit <- sfm(y_pcs_z ~ x1 + x2 | z, data = d, model_name = mn)
    me <- marginal_effects(fit)
    zs <- fit$z_spec
    expect_equal(nrow(me), nrow(d), info = mn)

    ## Rebuild sigma_u(z) independently and difference the moment directly.
    dl <- zs$delta
    s_of <- function(zval) exp(dl[[1]] * 1 + dl[[2]] * zval)
    mom_e <- if (mn == "NHN_Z") {
      function(zval) s_of(zval) * sqrt(2 / pi)
    } else {
      function(zval) s_of(zval)
    }
    mom_v <- if (mn == "NHN_Z") {
      function(zval) s_of(zval)^2 * (1 - 2 / pi)
    } else {
      function(zval) s_of(zval)^2
    }
    h <- 1e-6
    for (i in c(1L, 7L, 50L)) {
      z0 <- zs$Z[i, 2]
      expect_equal(me$dE_u.dz[i], (mom_e(z0 + h) - mom_e(z0 - h)) / (2 * h),
                   tolerance = 1e-5, info = paste(mn, i))
      expect_equal(me$dVar_u.dz[i], (mom_v(z0 + h) - mom_v(z0 - h)) / (2 * h),
                   tolerance = 1e-5, info = paste(mn, i))
    }
  }
})

test_that("the reported moments are consistent with the reported sigma_u", {
  skip_on_cran()
  fit <- sfm(y_pcs_z ~ x1 + x2 | z, data = cs_small(N = 200), model_name = "NHN_Z")
  me <- marginal_effects(fit)
  expect_equal(me$E_u, me$sigma_u * sqrt(2 / pi), tolerance = 1e-12)
  expect_equal(me$Var_u, me$sigma_u^2 * (1 - 2 / pi), tolerance = 1e-12)
  ## Under the SD link the effect on the variance is exactly twice delta*Var.
  expect_equal(me$dVar_u.dz / me$Var_u, rep(2 * fit$z_spec$delta[[2]], nrow(me)),
               tolerance = 1e-12)
})

test_that("average = TRUE agrees with the column means", {
  skip_on_cran()
  fit <- sfm(y_pcs_z ~ x1 + x2 | z, data = cs_small(N = 200), model_name = "NHN_Z")
  me <- marginal_effects(fit)
  avg <- marginal_effects(fit, average = TRUE)
  expect_equal(avg, attr(me, "average"))
  expect_equal(unname(avg[["dE_u.dz"]]), mean(me$dE_u.dz), tolerance = 1e-12)
  expect_equal(unname(avg[["dVar_u.dz"]]), mean(me$dVar_u.dz), tolerance = 1e-12)
})

test_that("the variance link halves the effect on E[u] relative to the SD link", {
  ## Same delta, same z, only the convention differs. This is entry C1 of the
  ## gap list expressed as a test: sfm()'s _Z models are on the SD scale and
  ## psfm()'s on the variance scale, and the marginal effects differ by two.
  Z <- cbind("(Intercept)" = rep(1, 5), z = c(-1, -0.5, 0, 0.5, 1))
  dl <- c("(Intercept)" = 0.9, z = 0.6)
  mk <- function(link) {
    structure(list(model_name = "fake",
                   z_spec = list(Z = Z, delta = dl, link = link,
                                 family = "halfnormal")),
              class = "sfareg")
  }
  sd_me <- marginal_effects(mk("sd"))
  var_me <- marginal_effects(mk("variance"))

  ## sigma_u itself differs: exp(eta) vs sqrt(exp(eta)).
  eta <- as.numeric(Z %*% dl)
  expect_equal(sd_me$sigma_u, exp(eta), tolerance = 1e-12)
  expect_equal(var_me$sigma_u, sqrt(exp(eta)), tolerance = 1e-12)

  ## Holding the moment fixed, the derivative factor is halved.
  expect_equal(sd_me$dE_u.dz / sd_me$E_u, rep(dl[["z"]], 5), tolerance = 1e-12)
  expect_equal(var_me$dE_u.dz / var_me$E_u, rep(dl[["z"]] / 2, 5), tolerance = 1e-12)
  expect_equal(sd_me$dVar_u.dz / sd_me$Var_u, rep(2 * dl[["z"]], 5), tolerance = 1e-12)
  expect_equal(var_me$dVar_u.dz / var_me$Var_u, rep(dl[["z"]], 5), tolerance = 1e-12)
})

test_that("intercept-only and homoskedastic fits are refused by name", {
  skip_on_cran()
  fit <- sfm(y_pcs ~ x1 + x2, data = cs_small(N = 150), model_name = "NHN")
  expect_null(fit$z_spec)
  expect_error(marginal_effects(fit), "NHN")
  expect_error(marginal_effects(fit), "homoskedastic")
  expect_error(marginal_effects(structure(list(), class = "lm")), "sfareg")

  ## A z design with no varying column has nothing to differentiate.
  const <- structure(
    list(model_name = "fake",
         z_spec = list(Z = cbind("(Intercept)" = rep(1, 4)),
                       delta = c("(Intercept)" = 0.5),
                       link = "sd", family = "halfnormal")),
    class = "sfareg")
  expect_error(marginal_effects(const), "constant")
})
