## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
options(rmarkdown.html_vignette.check_title = FALSE)

## ----setup--------------------------------------------------------------------
library(makemyprior)

## -----------------------------------------------------------------------------
# neighborhood graph
graph_path <- paste0(path.package("makemyprior"), "/neonatal.graph")

formula <- y ~ urban + mc(nu) + mc(v) +
  mc(u, model = "besag", graph = graph_path, scale.model = TRUE)


## -----------------------------------------------------------------------------
set.seed(1)
find_pc_prior_param(lower = 0.1, upper = 10, prob = 0.9, N = 2e5)


## -----------------------------------------------------------------------------
prior1 <- make_prior(
  formula, neonatal_data, family = "binomial",
  prior = list(tree = "s1 = (u, v); s2 = (s1, nu)",
               w = list(s1 = list(prior = "pc0", param = 0.25),
                        s2 = list(prior = "pc1", param = 0.75)),
               V = list(s2 = list(prior = "pc",
                                  param = c(3.35, 0.05)))))

prior1


## ---- fig.width = 5, fig.height = 3, eval = FALSE-----------------------------
#  plot_prior(prior1) # or plot(prior1)

## ---- fig.width = 5, fig.height = 3, eval = FALSE-----------------------------
#  plot_tree_structure(prior1)

## ---- eval = FALSE------------------------------------------------------------
#  posterior1 <- inference_stan(prior1, iter = 15000, warmup = 5000,
#                              seed = 1, init = "0", chains = 1)
#  
#  plot_posterior_stan(posterior1, param = "prior", plot_prior = TRUE)
#  

## ---- eval = FALSE------------------------------------------------------------
#  posterior1_inla <- inference_inla(prior1, Ntrials = neonatal_data$Ntrials)
#  plot_posterior_stdev(posterior1_inla)
#  

## -----------------------------------------------------------------------------
prior2 <- make_prior(formula, neonatal_data, family = "binomial")

prior2


## ---- fig.width = 5, fig.height = 3, eval = FALSE-----------------------------
#  plot_prior(prior2)

## ---- fig.width = 5, fig.height = 3, eval = FALSE-----------------------------
#  plot_tree_structure(prior2)

## ---- eval = FALSE------------------------------------------------------------
#  posterior2 <- inference_stan(prior2, iter = 15000, warmup = 5000,
#                              seed = 1, init = "0", chains = 1)
#  
#  plot_posterior_stan(posterior2, param = "prior", plot_prior = TRUE)
#  

## -----------------------------------------------------------------------------
sessionInfo()

