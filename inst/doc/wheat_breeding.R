## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
options(rmarkdown.html_vignette.check_title = FALSE)

## ----setup--------------------------------------------------------------------
library(makemyprior)

## -----------------------------------------------------------------------------
wheat_data_scaled <- wheat_data
wheat_data_scaled$Q_a <- scale_precmat(wheat_data$Q_a)
wheat_data_scaled$Q_d <- scale_precmat(wheat_data$Q_d)
wheat_data_scaled$Q_x <- scale_precmat(wheat_data$Q_x)

formula <- y ~
  mc(a, model = "generic0", Cmatrix = Q_a, constr = TRUE) +
  mc(d, model = "generic0", Cmatrix = Q_d, constr = TRUE) +
  mc(x, model = "generic0", Cmatrix = Q_x, constr = TRUE)


## -----------------------------------------------------------------------------
prior1 <- make_prior(formula, wheat_data_scaled, prior = list(
  tree = "s1 = (d, x); s2 = (a, s1); s3 = (s2, eps)",
  w = list(s1 = list(prior = "pcM", param = c(0.67, 0.8)),
           s2 = list(prior = "pcM", param = c(0.85, 0.8)),
           s3 = list(prior = "pc0", param = 0.25))))

prior1


## ---- fig.width = 5, fig.height = 3-------------------------------------------
plot_prior(prior1) # or plot(prior1)

## ---- fig.width = 5, fig.height = 3-------------------------------------------
plot_tree_structure(prior1)

## ---- eval = FALSE------------------------------------------------------------
#  posterior1 <- inference_stan(prior1, iter = 15000, warmup = 5000,
#                              chains = 1, seed = 1)
#  
#  plot_posterior_stan(posterior1, param = "prior", prior = TRUE)
#  

## ---- eval = FALSE------------------------------------------------------------
#  posterior1_inla <- inference_inla(prior1)
#  plot_posterior_stdev(posterior1_inla)
#  

## -----------------------------------------------------------------------------
prior2 <- make_prior(formula, wheat_data_scaled, prior = list(
  tree = "s1 = (d, x); s2 = (a, s1); (eps)",
  w = list(s1 = list(prior = "dirichlet"),
           s2 = list(prior = "pc1", param = c(0.8))),
  V = list(s2 = list(prior = "pc", param = c(3, 0.05)),
           eps = list(prior = "pc", param = c(3, 0.05)))))

prior2


## ---- fig.width = 5, fig.height = 3-------------------------------------------
plot_prior(prior2)

## ---- fig.width = 5, fig.height = 3-------------------------------------------
plot_tree_structure(prior2)

## ---- eval = FALSE------------------------------------------------------------
#  posterior2 <- inference_stan(prior2, iter = 15000, warmup = 5000,
#                              chains = 1, seed = 1)
#  
#  plot_posterior_stan(posterior2, param = "prior", prior = TRUE)
#  

## -----------------------------------------------------------------------------
sessionInfo()

