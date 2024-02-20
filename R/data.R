

#' Genomic wheat breeding model data
#'
#' Simulated wheat yield data with 100 observations.
#'
#' @format A list with the following variables
#' \describe{
#'   \item{y}{Response}
#'   \item{a, b, x}{Indexes for the additive, dominance and epistasis genetic effects, respectively}
#'   \item{Q_a, Q_d, Q_x}{Precision matrices for the genetic effects}
#' }
#' @examples
#' \dontrun{
#'
#' vignette("wheat_breeding", package = "makemyprior")
#' }
#'
#' if (interactive() && requireNamespace("rstan")){
#'
#'   wheat_data_scaled <- wheat_data
#'   wheat_data_scaled$Q_a <- scale_precmat(wheat_data$Q_a)
#'   wheat_data_scaled$Q_d <- scale_precmat(wheat_data$Q_d)
#'   wheat_data_scaled$Q_x <- scale_precmat(wheat_data$Q_x)
#'
#'   formula <- y ~
#'     mc(a, model = "generic0", Cmatrix = Q_a, constr = TRUE) +
#'     mc(d, model = "generic0", Cmatrix = Q_d, constr = TRUE) +
#'     mc(x, model = "generic0", Cmatrix = Q_x, constr = TRUE)
#'
#'   prior <- make_prior(formula, wheat_data_scaled, prior = list(
#'     tree = "s1 = (d, x); s2 = (a, s1); s3 = (s2, eps)",
#'     w = list(s1 = list(prior = "pcM", param = c(0.67, 0.8)),
#'              s2 = list(prior = "pcM", param = c(0.85, 0.8)),
#'              s3 = list(prior = "pc0", param = 0.25))))
#'
#'   posterior <- inference_stan(prior, iter = 150, warmup = 50,
#'                               chains = 1, seed = 1)
#'   # Note: For reliable results, increase the number of iterations
#'
#'   plot(prior)
#'   plot_tree_structure(prior)
#'   plot_posterior_fixed(posterior)
#'   plot_posterior_stan(posterior, param = "prior", prior = TRUE)
#'
#' }
#'
#' \dontrun{
#'
#' posterior <- inference_stan(prior, iter = 150, warmup = 50,
#'                             chains = 1, seed = 1)
#'
#' plot(prior)
#' plot_tree_structure(prior)
#' plot_posterior_fixed(posterior)
#' plot_posterior_stan(posterior, param = "prior", prior = TRUE)
#' }
#'
"wheat_data"



#' Latin square experiment data
#'
#' Simulated dataset for latin square experiment with 81 observations.
#'
#' @format A list with the following variables:
#' \describe{
#'   \item{y}{Response}
#'   \item{lin}{Covariate for linear effect of treatment}
#'   \item{row}{Row indexes}
#'   \item{col}{Column indexes}
#'   \item{treat_iid, treat_rw2}{Treatment indexes}
#' }
#' @examples
#' \dontrun{
#'
#' vignette("latin_square", package = "makemyprior")
#' }
#'
#' if (interactive() && requireNamespace("rstan")){
#'
#'   formula <- y ~ lin + mc(row) + mc(col) + mc(iid) +
#'     mc(rw2, model = "rw2", constr = TRUE, lin_constr = TRUE)
#'
#'   prior <- make_prior(
#'     formula, latin_data,
#'     prior = list(tree = "s1 = (rw2, iid);
#'                                  s2 = (row, col, s1); s3 = (s2, eps)",
#'                  w = list(s1 = list(prior = "pc0", param = 0.25),
#'                           s2 = list(prior = "dirichlet"),
#'                           s3 = list(prior = "pc0", param = 0.25))))
#'
#'   posterior <- inference_stan(prior, iter = 150, warmup = 50,
#'                               seed = 1, init = "0", chains = 1)
#'   # Note: For reliable results, increase the number of iterations
#'
#'   plot(prior)
#'   plot_tree_structure(prior)
#'   plot_posterior_fixed(posterior)
#'   plot_posterior_stan(posterior, param = "prior", prior = TRUE)
#' }
#'
#' \dontrun{
#'
#' posterior <- inference_stan(prior, iter = 15000, warmup = 5000,
#'                             seed = 1, init = "0", chains = 1)
#'
#' plot(prior)
#' plot_tree_structure(prior)
#' plot_posterior_fixed(posterior)
#' plot_posterior_stan(posterior, param = "prior", prior = TRUE)
#' }
#'
"latin_data"


#' Neonatal mortality data
#'
#' Simulated neonatal mortality data with 323 observations.
#'
#' @format A list with the following variables:
#' \describe{
#'   \item{y}{Response}
#'   \item{Ntrials}{Number of trials for each cluster}
#'   \item{urban}{Covariate indicating if cluster is urban (1) or rural (0)}
#'   \item{nu}{Cluster effect indexes}
#'   \item{v}{County effect indexes for iid effect}
#'   \item{u}{County effect indexes for Besag effect}
#' }
#' @examples
#' \dontrun{
#'
#' vignette("neonatal_mortality", package = "makemyprior")
#' }
#'
#' if (interactive() && requireNamespace("rstan")){
#'
#'   graph_path <- paste0(path.package("makemyprior"), "/neonatal.graph")
#'
#'   formula <- y ~ urban + mc(nu) + mc(v) +
#'     mc(u, model = "besag", graph = graph_path, scale.model = TRUE)
#'
#'   set.seed(1)
#'   find_pc_prior_param(lower = 0.1, upper = 10, prob = 0.9, N = 2e5)
#'
#'   prior <- make_prior(
#'     formula, neonatal_data, family = "binomial",
#'     prior = list(tree = "s1 = (u, v); s2 = (s1, nu)",
#'                  w = list(s1 = list(prior = "pc0", param = 0.25),
#'                           s2 = list(prior = "pc1", param = 0.75)),
#'                  V = list(s2 = list(prior = "pc",
#'                                     param = c(3.35, 0.05)))))
#'
#'   posterior <- inference_stan(prior, iter = 150, warmup = 50,
#'                               seed = 1, init = "0", chains = 1)
#'   # Note: For reliable results, increase the number of iterations
#'
#'   plot(prior)
#'   plot_tree_structure(prior)
#'   plot_posterior_fixed(posterior)
#'   plot_posterior_stan(posterior, param = "prior", prior = TRUE)
#' }
#'
#' \dontrun{
#'
#' posterior <- inference_stan(prior, iter = 15000, warmup = 5000,
#'                             seed = 1, init = "0", chains = 1)
#'
#' plot(prior)
#' plot_tree_structure(prior)
#' plot_posterior_fixed(posterior)
#' plot_posterior_stan(posterior, param = "prior", prior = TRUE)
#' }
"neonatal_data"



