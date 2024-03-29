% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_posterior_diagnostics.R
\name{plot_marginal_prior}
\alias{plot_marginal_prior}
\title{Plotting prior for a single parameter (weight or variance (not standard deviation))}
\usage{
plot_marginal_prior(x, obj, param, sd = FALSE)
}
\arguments{
\item{x}{Values to evaluate prior in.}

\item{obj}{An object from \code{make_prior}, \code{inference_stan}, \code{inference_inla}, or \code{makemyprior_gui}}

\item{param}{Name of parameter to plot (see \code{print(obj)} for syntax). Note that only variances
will be plotted, so \code{V[..]} and \code{sigma^2[..]} must be used to indicate those parameters.}

\item{sd}{Whether to plot variance parameters on the standard deviation (\code{TRUE}) or variance (\code{FALSE}, default) scale}
}
\value{
A \link[ggplot2]{ggplot} with the posterior distribution.
See also \link[makemyprior]{makemyprior_plotting}.
}
\description{
Following the parameterization of the prior.
}
\examples{
ex_prior <- makemyprior_example_model()
plot_marginal_prior(seq(0, 1, 0.001), ex_prior, "w[a/a_b]")
plot_marginal_prior(seq(0, 1, 0.001), ex_prior, "w[eps/eps_a_b]")
plot_marginal_prior(seq(0, 5, 0.01), ex_prior, "V[eps_a_b]")

}
\keyword{plot}
