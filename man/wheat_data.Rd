% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{wheat_data}
\alias{wheat_data}
\title{Genomic wheat breeding model data}
\format{
A list with the following variables
\describe{
  \item{y}{Response}
  \item{a, b, x}{Indexes for the additive, dominance and epistasis genetic effects, respectively}
  \item{Q_a, Q_d, Q_x}{Precision matrices for the genetic effects}
}
}
\usage{
wheat_data
}
\description{
Simulated wheat yield data with 100 observations.
}
\examples{
\dontrun{

vignette("wheat_breeding", package = "makemyprior")
}

if (interactive() && requireNamespace("rstan")){

  wheat_data_scaled <- wheat_data
  wheat_data_scaled$Q_a <- scale_precmat(wheat_data$Q_a)
  wheat_data_scaled$Q_d <- scale_precmat(wheat_data$Q_d)
  wheat_data_scaled$Q_x <- scale_precmat(wheat_data$Q_x)

  formula <- y ~
    mc(a, model = "generic0", Cmatrix = Q_a, constr = TRUE) +
    mc(d, model = "generic0", Cmatrix = Q_d, constr = TRUE) +
    mc(x, model = "generic0", Cmatrix = Q_x, constr = TRUE)

  prior <- make_prior(formula, wheat_data_scaled, prior = list(
    tree = "s1 = (d, x); s2 = (a, s1); s3 = (s2, eps)",
    w = list(s1 = list(prior = "pcM", param = c(0.67, 0.8)),
             s2 = list(prior = "pcM", param = c(0.85, 0.8)),
             s3 = list(prior = "pc0", param = 0.25))))

  posterior <- inference_stan(prior, iter = 150, warmup = 50,
                              chains = 1, seed = 1)
  # Note: For reliable results, increase the number of iterations

  plot(prior)
  plot_tree_structure(prior)
  plot_posterior_fixed(posterior)
  plot_posterior_stan(posterior, param = "prior", prior = TRUE)

}

\dontrun{

posterior <- inference_stan(prior, iter = 150, warmup = 50,
                            chains = 1, seed = 1)

plot(prior)
plot_tree_structure(prior)
plot_posterior_fixed(posterior)
plot_posterior_stan(posterior, param = "prior", prior = TRUE)
}

}
\keyword{datasets}
