## Shared fixtures. Kept small on purpose: these run inside R CMD check, so
## every second here is a second on every CRAN platform. Anything that needs a
## statistically meaningful sample size generates its own data behind
## skip_on_cran().

cs_small <- function(N = 300, rand = 42)
  data_gen_cs(N = N, rand = rand, sig_u = 1, sig_v = 0.3,
              cons = 0.5, beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.5)

panel_small <- function(t = 5, N = 40, rand = 42)
  as.data.frame(data_gen_p(t = t, N = N, rand = rand, sig_u = 1, sig_v = 0.3,
                           sig_r = 0.2, sig_h = 0.4, cons = 0.5,
                           beta1 = 0.5, beta2 = 0.5))

## psfm(model_name = "TFE") deliberately warns that the name changed meaning in
## 1.1.3. That warning is the point of a dedicated test elsewhere; everywhere
## else it is noise.
fit_tfe_quietly <- function(expr)
  withCallingHandlers(expr, warning = function(w)
    if (grepl("now fits Greene", conditionMessage(w))) invokeRestart("muffleWarning"))

## Reparameterization used throughout: models report lambda = sigma_u/sigma_v
## and sigma = sqrt(sigma_u^2 + sigma_v^2) rather than the raw sigmas.
sig_u_from <- function(lambda, sigma) (lambda * sigma) / sqrt(1 + lambda^2)
sig_v_from <- function(lambda, sigma) sig_u_from(lambda, sigma) / lambda
