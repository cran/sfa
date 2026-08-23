## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup1-------------------------------------------------------------------
library(sfa)

## ----setup2-------------------------------------------------------------------
data_trial <- data_gen_p(t=6, N=70, rand = 16, sig_u = 0.3, sig_v = 0.1, sig_r = 0.1, sig_h = 0.3, cons = 0.5, beta1 = 0.5, beta2 = 0.5)

p.gtre_sml   <- psfm(formula      = y_gtre ~ x1 + x2,
                     model_name   = "GTRE",
                     estimator    = "sml",
                     data         = data_trial,
                     individual   = "name",
                     PSopt        = TRUE,
                     optHessian   = TRUE,
                     halton_num   = 50,
                     rand.gtre    = 1,
                     rand.psoptim = 1,
                     maxit.bobyqa = 150,
                     maxit.psoptim= 10,
                     maxit.optim  = 10)
summary(p.gtre_sml)
mean(p.gtre_sml$U)
mean(p.gtre_sml$H)

## ----setup3-------------------------------------------------------------------
plot(density(p.gtre_sml$U),main="Density of Transient TE")
plot(density(p.gtre_sml$H),main="Density of Persistent TE")

## ----setup4-------------------------------------------------------------------
total_te <- rep(p.gtre_sml$H, each=6) * p.gtre_sml$U
plot(density(total_te),main="Density of Total TE")

## ----setup5-------------------------------------------------------------------
coef(p.gtre_sml)                 # named vector of point estimates
vcov(p.gtre_sml)                 # variance-covariance matrix (from the Hessian)
logLik(p.gtre_sml)               # log-likelihood, with df/nobs attributes set
AIC(p.gtre_sml); BIC(p.gtre_sml) # available "for free" once logLik() works

## ----setup6-------------------------------------------------------------------
## data_trial already holds every column data_gen_p() produces, including the
## y_gtre_zz/z_gtre/zp_gtre trio this model needs -- no need to simulate again.
p.gtre_z <- psfm(formula      = y_gtre_zz ~ x1 + x2 | z_gtre | zp_gtre,
                 model_name   = "GTRE_Z",
                 data         = data_trial,
                 individual   = "name",
                 PSopt        = TRUE,
                 optHessian   = TRUE,
                 halton_num   = 50,
                 rand.gtre    = 1,
                 rand.psoptim = 1,
                 maxit.bobyqa = 150,
                 maxit.psoptim= 10,
                 maxit.optim  = 10)
summary(p.gtre_z)

## ----setup7-------------------------------------------------------------------
data_trial_tre <- data_gen_p(t=5, N=30, rand=16, sig_u=0.3, sig_v=0.1, sig_r=0.1, sig_h=0.3,
                             cons=0.5, beta1=0.5, beta2=0.5)

p.tre <- psfm(formula = y_tre ~ x1 + x2, model_name = "TRE",
              data = data_trial_tre, individual = "name",
              halton_num = 50, rand.gtre = 1, maxit.bobyqa = 300)

set.seed(1)
boot <- psfm_bootstrap(p.tre,
                       numCores      = 2,
                       BOOT          = 5,    # a real analysis should use far more, e.g. 199-999
                       individual    = "name",
                       inefdec       = TRUE,
                       maxit.bobyqa  = 150,
                       maxit.psoptim = 30)
boot$se                    # bootstrap standard errors, one per parameter in coef(p.tre)
boot$model$out             # a copy of p.tre$out with bootstrap SEs/t-values written in

