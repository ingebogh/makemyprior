## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
options(rmarkdown.html_vignette.check_title = FALSE)

## ----setup--------------------------------------------------------------------
library(makemyprior)

## -----------------------------------------------------------------------------
formula <- y ~ lin + mc(row) + mc(col) + mc(iid, constr = TRUE) +
  mc(rw2, model = "rw2", constr = TRUE, lin_constr = TRUE)

## -----------------------------------------------------------------------------
prior1 <- make_prior(
  formula, latin_data,
  prior = list(tree = "s1 = (iid, rw2); s2 = (row, col, s1); s3 = (s2, eps)",
               w = list(s1 = list(prior = "pc1", param = 0.75),
                        s2 = list(prior = "dirichlet"),
                        s3 = list(prior = "pc0", param = 0.25))))

prior1


## ---- fig.width = 6, fig.height = 3-------------------------------------------
plot_prior(prior1) # or plot(prior)

## ---- fig.width = 3, fig.height = 3-------------------------------------------
plot_tree_structure(prior1)

## ---- eval = FALSE------------------------------------------------------------
#  posterior1 <- inference_stan(prior1, iter = 15000, warmup = 5000,
#                              seed = 1, init = "0", chains = 1)
#  plot_posterior_stan(posterior1, param = "prior", prior = TRUE)
#  

## ---- eval = FALSE------------------------------------------------------------
#  formula_inla <- y ~ lin + mc(row) + mc(col) + mc(iid, constr = TRUE) +
#    mc(rw2, model = "rw2", constr = TRUE, extraconstr = list(A = matrix(1:9, 1, 9), e = matrix(0, 1, 1)))
#  prior1_inla <- make_prior(
#    formula_inla, latin_data,
#    prior = list(tree = "s1 = (iid, rw2); s2 = (row, col, s1); s3 = (s2, eps)",
#                 w = list(s1 = list(prior = "pc1", param = 0.75),
#                          s2 = list(prior = "dirichlet"),
#                          s3 = list(prior = "pc0", param = 0.25))))
#  
#  posterior1_inla <- inference_inla(prior1_inla)
#  plot_posterior_stdev(posterior1_inla)
#  

## -----------------------------------------------------------------------------
prior2 <- make_prior(
  formula, latin_data,
  prior = list(tree = "s1 = (iid, rw2); s2 = (row, col, s1); (eps)",
               w = list(s1 = list(prior = "pc1", param = 0.75),
                        s2 = list(prior = "dirichlet")),
               V = list(s2 = list(prior = "pc", param = c(sqrt(0.2), 0.05)),
                        eps = list(prior = "pc", param = c(sqrt(0.2), 0.05)))))

prior2


## ---- fig.width = 6, fig.height = 3-------------------------------------------
plot_prior(prior2) # or plot(prior2)

## ---- fig.width = 3, fig.height = 3-------------------------------------------
plot_tree_structure(prior2)

## ---- eval = FALSE------------------------------------------------------------
#  posterior2 <- inference_stan(prior2, iter = 15000, warmup = 5000,
#                              seed = 1, init = "0", chains = 1)
#  plot_posterior_stan(posterior2)
#  

## -----------------------------------------------------------------------------
prior3 <- make_prior(
  formula, latin_data,
  prior = list(tree = "(row); (col); (iid); (rw2); (eps)",
               V = list(
                 row = list(prior = "pc", param = c(sqrt(0.1), 0.05)),
                 col = list(prior = "pc", param = c(sqrt(0.1), 0.05)),
                 iid = list(prior = "pc", param = c(sqrt(0.1), 0.05)),
                 rw2 = list(prior = "pc", param = c(sqrt(0.1), 0.05)),
                 eps = list(prior = "pc", param = c(sqrt(0.1), 0.05))
               )))

prior3


## ---- fig.width = 6, fig.height = 3-------------------------------------------
plot_prior(prior3) # or plot(prior3)

## ---- fig.width = 3, fig.height = 3-------------------------------------------
plot_tree_structure(prior3)

## ---- eval = FALSE------------------------------------------------------------
#  posterior3 <- inference_stan(prior3, iter = 15000, warmup = 5000,
#                              seed = 1, init = "0", chains = 1)
#  plot_posterior_stan(posterior3)
#  

## -----------------------------------------------------------------------------
sessionInfo()

