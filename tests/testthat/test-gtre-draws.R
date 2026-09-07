## Per-firm Halton blocks (gap J1) and the two GTRE likelihood paths.

test_that("each firm gets its own block of draws", {
  D <- sfa:::.gtre_halton_draws(N = 5, R = 8, rand.gtre = 7L)
  expect_length(D, 5L)
  expect_true(all(vapply(D, function(m) identical(dim(m), c(8L, 2L)), logical(1))))
  ## The whole point: distinct blocks. Before this they were all the same
  ## object, which is what made the vectorized likelihood's shortcut valid.
  for (i in 2:5) expect_false(isTRUE(all.equal(D[[1]], D[[i]])))
  expect_true(all(is.finite(unlist(D))))
})

test_that("the randomization is a shift, and is reproducible", {
  a <- sfa:::.gtre_halton_draws(N = 3, R = 16, rand.gtre = 7L)
  b <- sfa:::.gtre_halton_draws(N = 3, R = 16, rand.gtre = 7L)
  expect_equal(a, b)
  ## A different seed must actually move them.
  cc <- sfa:::.gtre_halton_draws(N = 3, R = 16, rand.gtre = 99L)
  expect_false(isTRUE(all.equal(a, cc)))
  ## And the unrandomized sequence is deterministic.
  expect_equal(sfa:::.gtre_halton_draws(N = 2, R = 8),
               sfa:::.gtre_halton_draws(N = 2, R = 8))
})

test_that("the vectorized and loop likelihoods agree", {
  skip_on_cran()
  ## sfa.gtre_vectorized hoists the draws out of the per-firm loop, which was
  ## only valid while every firm shared one block. When firms got independent
  ## blocks the shortcut silently evaluated a DIFFERENT likelihood -- using
  ## firm 1's draws for everyone -- and nothing tested the two against each
  ## other. This is that test.
  d <- data_gen_p(t = 6, N = 40, rand = 7, sig_u = 1, sig_v = 0.3, sig_r = 0.2,
                  sig_h = 0.4, cons = 0.5, beta1 = 0.5, beta2 = 0.5)
  pd <- plm::pdata.frame(d, index = c("name", "year"))
  on.exit(options(sfa.gtre_vectorized = FALSE), add = TRUE)

  for (m in c("GTRE", "TRE")) {
    f <- if (identical(m, "GTRE")) y_gtre ~ x1 + x2 else y_tre ~ x1 + x2
    est <- if (identical(m, "GTRE")) "sml" else "fiml"
    options(sfa.gtre_vectorized = FALSE)
    a <- suppressWarnings(psfm(f, m, pd, individual = "name", halton_num = 30,
                               rand.gtre = 7L, estimator = est))
    options(sfa.gtre_vectorized = TRUE)
    b <- suppressWarnings(psfm(f, m, pd, individual = "name", halton_num = 30,
                               rand.gtre = 7L, estimator = est))
    expect_equal(as.numeric(logLik(a)), as.numeric(logLik(b)),
      tolerance = 1e-6, info = m)
    expect_equal(unname(coef(a)), unname(coef(b)), tolerance = 1e-3, info = m)
  }
})
