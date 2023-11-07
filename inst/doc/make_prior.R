## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
options(rmarkdown.html_vignette.check_title = FALSE)

## ----setup--------------------------------------------------------------------
library(makemyprior)

## -----------------------------------------------------------------------------
formula <- y ~ x + mc(a) + mc(b)

## -----------------------------------------------------------------------------
p <- 10
m <- 10
n <- m*p

set.seed(1)
data <- list(a = rep(1:p, each = m),
             b = rep(1:m, times = p),
             x = runif(n))

data$y <- data$x + rnorm(p, 0, 0.5)[data$a] +
  rnorm(m, 0, 0.3)[data$b] + rnorm(n, 0, 1)

## -----------------------------------------------------------------------------
prior <- make_prior(formula, data, family = "gaussian",
                    intercept_prior = c(0, 1000),
                    covariate_prior = list(x = c(0, 100)))

## ---- fig.width = 5, fig.height = 2-------------------------------------------
summary(prior)

plot_prior(prior) # or plot(prior)
plot_tree_structure(prior)

## ---- eval = FALSE------------------------------------------------------------
#  new_prior <- makemyprior_gui(prior)

## ---- fig.width = 5, fig.height = 2-------------------------------------------

new_prior <- make_prior(
  formula, data,
  prior = list(
    tree = "s1 = (a, b); s2 = (s1, eps)",
    w = list(s1 = list(prior = "pcM", param = c(0.7, 0.5)),
             s2 = list(prior = "pc1", param = 0.75)),
    V = list(s2 = list(prior = "pc0", param = c(3, 0.05)))
    ),
    covariate_prior = list(x = c(0, 100))
)

summary(new_prior)
plot_prior(new_prior)
plot_tree_structure(new_prior)


## ---- eval = FALSE------------------------------------------------------------
#  compile_stan(save = TRUE)

## ---- eval = FALSE------------------------------------------------------------
#  posterior1 <- inference_stan(new_prior, iter = 1e4, chains = 1, seed = 1)

## ---- fig.width = 5, fig.height = 2, eval = FALSE-----------------------------
#  plot_posterior_stan(posterior1, param = "prior", prior = TRUE) # on the scale of the prior, together with the prior
#  plot_posterior_stan(posterior1, param = "variance") # on variance scale
#  plot_fixed_posterior(posterior1) # fixed effects

## ---- fig.width = 5, fig.height = 2, eval = FALSE-----------------------------
#  prior1 <- inference_stan(new_prior, use_likelihood = FALSE, iter = 1e4, chains = 1, seed = 1)
#  plot_several_posterior_stan(list(Prior = prior1, Posterior = posterior1))

## ---- eval = FALSE------------------------------------------------------------
#  posterior2 <- inference_inla(new_prior)

## ---- fig.width = 5, fig.height = 2, eval = FALSE-----------------------------
#  plot_posterior_variance(posterior2) # on variance scale
#  plot_fixed_posterior(posterior1)

## -----------------------------------------------------------------------------
?makemyprior_models

## ---- fig.width = 5, fig.height = 2-------------------------------------------

prior2 <- make_prior(formula = formula, data = data,
                     prior = list(tree = "(a); (b); (eps)",
                                  V = list(
                                    a = list(prior = "pc", param = c(1, 0.05)),
                                    b = list(prior = "pc", param = c(2, 0.05)),
                                    eps = list(prior = "pc", param = c(3, 0.05))
                                  )))

plot_prior(prior2)
plot_tree_structure(prior2)


## ---- fig.width = 5, fig.height = 2-------------------------------------------

prior3 <- make_prior(formula = formula, data = data,
                     prior = list(tree = "s1 = (a, b); (eps)",
                                  V = list(
                                    s1 = list(prior = "pc", param = c(3, 0.05)),
                                    eps = list(prior = "pc", param = c(3, 0.05))),
                                  w = list(
                                    s1 = list(prior = "pcM", param = c(0.5, 0.8))
                                  )
                                  ))

plot_prior(prior3)
plot_tree_structure(prior3)


## -----------------------------------------------------------------------------
sessionInfo()

