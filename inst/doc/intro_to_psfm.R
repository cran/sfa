## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup 1------------------------------------------------------------------
library(sfa)

## ----setup 2------------------------------------------------------------------
data_trial <- data_gen_p(t=10,N=100,rand = 16, sig_u = 0.3,sig_v = 0.1, sig_r = 0.1, sig_h = 0.3, cons = 0.5, beta1 = 0.5, beta2 = 0.5)

p.gtre_sml   <- psfm(formula      = y_gtre ~ x1 + x2,       
                     model_name   = "GTRE",
                     data         = data_trial,           
                     individual   = "name",
                     PSopt        = TRUE,
                     optHessian   = TRUE,
                     maxit.bobyqa = 150,
                     maxit.psoptim= 10,
                     maxit.optim  = 10)
summary(p.gtre_sml)
mean(p.gtre_sml$U)
mean(p.gtre_sml$H)

## ----setup 3------------------------------------------------------------------
plot(density(p.gtre_sml$U),main="Density of Transient TE")
plot(density(p.gtre_sml$H),main="Density of Persistent TE")

## ----setup 4------------------------------------------------------------------
total_te <- rep(p.gtre_sml$H, each=10) * p.gtre_sml$U
plot(density(total_te),main="Density of Total TE")

