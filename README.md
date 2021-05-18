
# makemyprior

makemyprior is a tool for easy prior construction and visualization. It
helps to formulates joint prior distributions for variance parameters in
latent Gaussian models. The resulting prior is robust and can be created
in an intuitive way. A graphical user interface (GUI) can be used to
choose the joint prior, where the user can click through the model and
select priors. An extensive guide is available in the GUI. The package
allows for direct inference with the specified model and prior. Using a
hierarchical variance decomposition, we formulate a joint variance prior
that takes the whole model structure into account. In this way, existing
knowledge can intuitively be incorporated at the level it applies to.
Alternatively, one can use independent variance priors for each model
components in the latent Gaussian model.

## Installation

You can install the released version of makemyprior from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("makemyprior")
```

## Example

This is an example showing how to implement a prior.

``` r

library(makemyprior)

set.seed(1)
data <- list(
  a = rep(1:10, each = 10),
  b = rep(1:10, times = 10)
)
data$y <- rnorm(10, 0, 0.4)[data$a] + rnorm(10, 0, 0.6)[data$b] + rnorm(100, 0, 1)

formula <- y ~ mc(a) + mc(b)

prior <- make_prior(formula, data, family = "gaussian",
                    prior = list(tree = "s1 = (a, b); s2 = (s1, eps)",
                                 w = list(s2 = list(prior = "pc0", param = 0.25)),
                                 V = list(s2 = list(prior = "pc", param = c(3, 0.05)))),
                    intercept_prior = c(0, 1000))

summary(prior) # printing details
plot(prior) # plotting prior

```
