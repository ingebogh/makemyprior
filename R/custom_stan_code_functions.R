



#' Create a "skeleton" for custom Stan code
#'
#' Makes and saves files with generic code
#' for writing custom Stan code and still use the HD prior.
#' @param save To confirm that files can be saved (default FALSE).
#' @param location Path to location.
#' @return Nothing.
#' @details
#' Must be in an interactive session to store the code.
#' A folder called \code{"my_stan_code"} will be created in \code{location}.
#' If the folder already exists in the
#' specified location, you get an error. The folder contains:
#' \describe{
#'   \item{\code{main_file.stan}}{Main file. Can put all necessary functions here, but for a cleaner code
#'   that is easier to read, we put the functions in separate files.}
#'   \item{\code{jacobian.stan}}{Function that automatically computes the Jacobian, needed
#'   to transform from weights and total variance parameterization to log-variance parameterization.}
#'   \item{\code{prior_distributions.stan}}{Functions for computing the prior distributions.}
#' }
#' The provided code is written so a random intercept model with an intercept, a group effect and
#' a residual effect can be fitted:
#' \describe{
#'   \item{\code{example_custom_stan.R}}{R script showing how one can fit a random intercept model using the provided code.}
#' }
#' The code can be expanded to fit the desired model.
#' This requires some knowledge about Stan.
#' No more documentation is given, as this is merely an offer to users who want to use other models than
#' what are provided in the package already, and will be highly model specific.
#' @examples
#' \dontrun{
#' create_stan_file(TRUE, "")
#' }
#'
#' @export
create_stan_file <- function(save = FALSE, location = ""){

  if (!save) stop("You must confirm that files can be saved in 'location' by adding argument 'save = TRUE'.", call. = FALSE)

  if (interactive()){

  if (substr(location, nchar(location), nchar(location)) != "/") location <- paste0(location, "/")

  tmp <- dir.exists(paste0(location, "my_stan_code/"))
  if (tmp){
    stop(paste0("There is already a folder in ",
                location,
                " called 'my_stan_code'. Rename it, delete it (if you don't need it), ",
                "or choose another location."), call. = FALSE)
  }
  if (!tmp){

    # creating folder
    dir.create(paste0(location, "my_stan_code/"))

    ### adding files ###

    # stan files
    write(main_file_text(), file = paste0(location, "my_stan_code/main_file.stan"))
    write(jacobian_text(), file = paste0(location, "my_stan_code/jacobian.stan"))
    write(prior_dists_text(), file = paste0(location, "my_stan_code/prior_distributions.stan"))

    # R script
    write(make_r_script(), file = paste0(location, "my_stan_code/example_custom_stan.R"))

  }

  cat("\"", paste0(location, "my_stan_code/"), "\" was created with files for customizing a Stan code ",
      "compatible with the prior created with 'make_prior' in the 'makemyprior' package. ",
      "The folder \"my_stan_code\" consists of:\n",
      " * main_file.stan,\n",
      " * jacobian.stan,\n",
      " * prior_distributions.stan,\n",
      " * example_custom_stan.R,\n",
      "Details are included in the files. Some knowledge about rstan is necessary.\n",
      sep = "")

  invisible()

  } else {
    stop("Cannot save the code when the R session is not interactive.", call. = FALSE)
  }

}


main_file_text <- function(){


  my_text <- "


////// INFO //////

// Main file for Stan code that can be customized.
// This require some knowledge on how to write Stan code, see
// https://mc-stan.org/ for details.
// The Jacobian function included is generic and uses some data that can be
// extracted from the object created by 'make_prior' in the R package 'makemyprior'.
// See the file 'example_custom_stan.R' for how this can be done.
// A custom function can also be written, this may increase the speed of the sampling.
// You can parameterize the model in any way you want, using the parameters of
// The desired prior distribution. Then the Jacobian is not necessary.
// We provide code for a parameterization on log-variances.

// Note that we provide very generic functions for computing the HD prior.

// A random intercept model is implemented here, and can be fitted using the code
// in 'example_custom_stan.R'.

functions {

// including a general jacobian-function
#include jacobian.stan
// including the prior distribution functions
#include prior_distributions.stan

}

data {

  int n; // number of observations
  int n_g; // number of observations in each group
  int n_par; // number of variance parameters, in the random intercept model we have 2

  int indexes_group[n]; // indexes for the group effect

  real y[n]; // observations

  int likelihood; // whether or not to sample with likelihood (1 = yes, gives posterior. 0 = no, gives prior)

  // these four are needed for computing PC prior on variance proportions, DO NOT REMOVE
  int n_pc; // number of variance proportions with a PC prior
  int n_knots; // number of knots for splines
  vector[n_knots] knots; // the knots for the splines
  matrix[n_knots,4] prior_coeffs[n_pc]; // coefficients for the splines

  //// input data needed if the generic PC prior functions are used

  int n_hd; // number of variance parameters involved in a HD prior
  int n_cw; // number of variances with individual priors (component-wise priors)
  int n_splits; // number of splits in the model
  int n_rows_in_w_o_u; // number of rows in w_o and w_u

  int which_pc[n_splits]; // which of the HD priors are PC (1) and Dirichlet (0)
  int which_hd[n_hd]; // which variance parameters are involved in HD priors
  int which_cw[n_cw]; // which variances have CW priors

  // matrices keeping track of which variance parameters are involved where in the tree
  int w_o[n_rows_in_w_o_u, n_par];
  int w_u[n_rows_in_w_o_u, n_par];

  int n_totvar; // number of total variances
  matrix[n_totvar,3] totvar_prior_info; // info on the total variance
  int which_theta_in_hd[n_totvar, n_par]; // which of the variance parameters are involved in a hd prior
  int n_splits_each_tree[n_totvar]; // number of splits in each tree

  int row_index_hd_pr[n_totvar, n_rows_in_w_o_u+1]; // which rows in w_o and w_u belongs to each split

  matrix[n_par,3] cw_prior_info; // info on individual priors on variance components

}

parameters{

  real intercept; // intercept
  // add fixed effects (if any) here

  real theta[n_par]; // log variances for group effect (1) and residual effect (2)

  // here you can add other types of parameters, such as a correlation
  // real cor;

  vector[n_g] group_effect; // group effect

}


transformed parameters{

  vector[n] eta;
  vector[n_g] group_new;

  eta = rep_vector(0, n);

  // add intercept and fixed effects (if any) to the linear predictor
  eta += rep_vector(intercept, n);

  // we use a sum-to-zero constraint for the group effect
  group_new = group_effect - mean(group_effect);
  // add random effects to the linear predictor
  eta += group_new[indexes_group] * sqrt(exp(theta[1]));

}




model {

  // sampling intercept and fixed effects (if any) prior
  target += normal_lpdf(intercept | 0, 1000);

  // sampling the joint variance prior
  target += joint_prior_lpdf(theta | which_pc, which_hd, w_o, w_u, likelihood, n_splits_each_tree, row_index_hd_pr, knots, prior_coeffs, n_knots, totvar_prior_info, cw_prior_info, which_cw, n_totvar, which_theta_in_hd);

  // to do this part yourself, you can:
  //    1) reparameterize the model and use the parameterization of the prior directly.
  //       then you can sample from each variance proportion using 'eval_spline_prior()' (for PC prior) here,
  //       or any other distribution you would like to use.
  //    2) use log-variances as internal parameterization, and write your own function to sample
  //       from the joint prior for the log variance parameters and call it here (in the function you write you can use any distribution on
  //       the parameters). here the included jacobian function can be useful, or you can write one yourself.

  // if any non-variance parameters in the model (e.g. a correlation), they can be sampled here, e.g.:
  // target += some_dist(cor);

  // sampling the random effects
  if (likelihood == 1){
    // note that we use the variance when we make the linear predictor eta in the 'transformed parameters' block
    target += normal_lpdf(group_effect | 0, 1);
  } else {
    // if we want to sample from the prior, we do not need to sample the effect with the parameters,
    // and can save time by just sampling the random effects in the model from a Gaussian(0,1) distribution
    target += normal_lpdf(group_effect | 0, 1);
  }

  // sampling from likelihood
  if (likelihood == 1){
    // mean equal to linear predictor, standard deviation equal to residual stdev
    target += normal_lpdf(y | eta, sqrt(exp(theta[2])));
  } else {
    // if we want to sample from the prior, we do not do this
    target += 0;
  }

}


generated quantities {

  // here you can generate other posterior quantities (such as the marginal log-likelihood, or reparamerizations),
  // that you do not need for sampling from the model
  // they can be sampled above as well, but it saves time to do it here instead

  real totvar;
  real varprop;
  vector[n] log_lik;

  // we want to get the variance proportion and total variance returned together with the effect
  // note that this can be computed after the inference as well
  totvar = exp(theta[1]) + exp(theta[2]);
  varprop = exp(theta[1])/totvar;

  // also including example on how to compute the marginal log-likelihood
  for (ind in 1:n){
    log_lik[ind] = normal_lpdf(y[ind] | eta[ind], sqrt(exp(theta[2])));
  }

}

  "

return(my_text)

}


jacobian_text <- function(){


  my_text <- "




// returns the indexes where w_row is 1
// used to pick correct thetas for the logit weights etc
int[] get_indexes(int[] w_row){

  int res[sum(w_row)];
  int ind;

  ind = 1;
  for (i in 1:num_elements(w_row)){
    if (w_row[i] == 1){
      res[ind] = i;
      ind += 1;
    }
  }

  return(res);

}

// returns the indexes where w_row_o is 0, and w_row_u is 1 at the same time
// used to pick correct thetas for the weights in the jacobian
int[] get_indexes2(int[] w_row_o, int[] w_row_u){

  int res[sum(w_row_u)-sum(w_row_o)];
  int ind;

  ind = 1;
  for (i in 1:num_elements(w_row_o)){
    if (w_row_o[i] == 0 && w_row_u[i] == 1){
      res[ind] = i;
      ind += 1;
    }
  }

  return(res);

}


// a function that calculates the jacobian for any structure, from variance proportion and total variance to log variances
// note that the jacobian for total variance is computed elsewhere

// it is not necessary to know if the prior on each weight is PC or Dirichlet, as we only want the expression for the weight

// this works for one tree, as we compute individual jacobians for each tree

// theta = all log variances in this tree (all thetas from the hd_prior_joint_lpdf function)
// w_o = matrix with variances over fraction in each weight
// w_u = matrix with variances under fraction in each weight
// n = number of variance components involved in the tree
// row_index_hd_pr = which rows in w_o and w_u that goes to each weight
real calc_jac_logdet(real[] theta, int[,] w_o, int[,] w_u, int n, int[] row_index_hd_pr, int n_splits){

  matrix[n, n] jac;
  int i_row; // index for row
  real under; // value under fraction bar for each weight
  real over; // value over fraction bar for each weight
  int bool_arg1;

  jac = rep_matrix(0, n, n); // initializing

  // first row in jacobian is for the log total variance
  for (i in 1:n){
    jac[1,i] = exp(theta[i]);
  }

  i_row = 2;

  // for each split in the tree
  for (i_split in 1:n_splits){

    // value under fraction bar is the same for all weights in a split
    under = ( sum(exp(theta[get_indexes(w_u[row_index_hd_pr[i_split],])])) )^2;

    over = 0; // if the parameter is not in the weight, the derivative is 0

    // looping through the weights in this split (not the last one, since that is 1 minus the others)
    // (must use -2 since the next index in row_index_hd_pr gives the index of the next weight (-1 here), and we do not need the last one (-1 more))
    for (ind in row_index_hd_pr[i_split]:(row_index_hd_pr[i_split+1]-2)){

      // looping through each parameter (i.e. each column of jacobian) relevant for this weight to calculate the part above the fraction bar
      for (i_col in get_indexes(w_u[ind,])){

        // testing if the variance corresponding to the index i_col is involved above the fraction bar of this weight
        bool_arg1 = 0;
        for (i in 1:sum(w_o[ind,])){
          if (i_col == get_indexes(w_o[ind,])[i]){
            bool_arg1 = 1;
          }
        }

        // if this variance parameter is above the fraction bar and we are differentiating wrt to it
        if (bool_arg1 == 1){

          over = ( sum(exp(theta[get_indexes2(w_o[ind,], w_u[ind,])])) ) * exp(theta[i_col]);

        } else { // if this variance parameter is not above the fraction bar and we are differentiating wrt to it

          over = -( sum(exp(theta[get_indexes(w_o[ind,])])) ) * exp(theta[i_col]);

        }

        jac[i_row, i_col] = over/under;

      }

      i_row += 1;

    }

  }

  return(log(fabs(determinant(jac))));

}

  "

return(my_text)

}


prior_dists_text <- function(){

  my_text <- "


////// Functions to compute priors. Note that the normalizing constant is not always included (as it is not needed)


// jeffreys prior on variance
// input is log variance
real jeffreys_logvar_lpdf(real x){
  return(0);
}

// inverse gamma on variance (gamma on precision)
// input is log variance
real invgamma_logvar_lpdf(real x, real shape, real rate){
  return(
    -shape*x - 1/(rate*exp(x))
  );
}

// half cauchy on standard deviation
// input is log variance
real halfcauchy_logvar_lpdf(real x, real shape){
  return(
    x/2 - log(pi()) - log(shape + exp(x)/shape)
  );
}

// pc prior (exponential on standard deviation)
// input is log variance
real pc_logvar_lpdf(real x, real U, real alpha){
  real shape = -log(alpha)/U;
  return(
    log(shape/2) + x/2 - shape*exp(x/2)
  );
}


// to get correct prior
// input is log variance, which prior is used, and a vector of parameters
real choose_prior_lpdf(real x, real prior_number, real[] param){

  if (prior_number == 1){
    return(jeffreys_logvar_lpdf(x));
  } else if (prior_number == 2){
    return(pc_logvar_lpdf(x | param[1], param[2]));
  } else if (prior_number == 3){
    return(invgamma_logvar_lpdf(x | param[1], param[2]));
  } else if (prior_number == 4){
    return(halfcauchy_logvar_lpdf(x | param[1]));
  } else if (prior_number == 5) {
    return(normal_lpdf(x | 0, 1));
  } else {
    return(0);
  }

}



// calculates the joint PC prior
// logitw is a value to evaluate the prior in
// knots is knots for the spline
// prior_coeffs is the whole array for the spline coefficients for ALL priors
// idx is which of the priors in prior_coeffs to use
// n_knots is the number of knows (length of knots)
real eval_spline_lpdf(real logitw, vector knots, matrix[] prior_coeffs, int idx, int n_knots){
    real dx;
    real val;
    int idd;
    int low;
    int high;
    int med;

    // Search for correct cubic polynomial
    low = 1;
    high = n_knots;
    while(high-low > 1){
        med = (high+low)/2;
        if(logitw < knots[med]){
            high = med;
        } else{
            low = med;
        }
    }
    idd = low;
    dx = logitw-knots[idd];
    val =  prior_coeffs[idx, idd, 1];
    val += prior_coeffs[idx, idd, 2]*dx;
    val += prior_coeffs[idx, idd, 3]*dx*dx;
    val += prior_coeffs[idx, idd, 4]*dx*dx*dx;

    return(val);
}

real expit(real x){
  return(exp(x)/(1+exp(x)));
}

// evaluates hd prior for a given split (which is identified by the parameter weights)
// includes jacobian from logit weight to weight!
// not including the jacobian from V/w to theta!
// theta = vector with all log variances in the model
// weights_start = index of where in w_o and w_u this weight is located
// w_o = matrix with variances over fraction in each weight
// w_u = matrix with variances under fraction in each weight
// the rest is for the spline
real hd_pc_prior_lpdf(real[] theta, int weights_start, int idx, int[,] w_o, vector knots, matrix[] prior_coeffs, int n_knots){

 real logdens;
 real logitw;
 logdens = 0;

 // calculate logit weight from log variances
 logitw = log(sum(exp(theta[get_indexes(w_o[weights_start,])]))) - log(sum(exp(theta[get_indexes(w_o[weights_start+1,])])));

 logdens += eval_spline_lpdf(logitw | knots, prior_coeffs, idx, n_knots);

 // jacobian from logit weight to weight
 logdens += -log( expit(logitw)*(1-expit(logitw)) );

 return(logdens);

}

// parameter value for the symmetric Dirichlet distribution for a given number of components involved
real get_dirichlet_parameter(int num_comp){

  real values[19] = {1, 0.7535, 0.6834, 0.6485, 0.6274, 0.6132, 0.6029, 0.5951, 0.589, 0.5842, 0.5801, 0.5768, 0.5739, 0.5715, 0.5693, 0.5675, 0.5658, 0.5643, 0.563};

  return(values[num_comp-1]);
}


// dirichlet for a split with variances theta
// do not include the jacobian from V/w to theta!
// theta = vector with all (yes, ALL) log variances in the model
// m = how many variance parameters are involved in this split
// weights = vector (of at least length 2) with which weights are involved in this particular split
// w_o = matrix with variances over fraction in each weight
// w_u = matrix with variances under fraction in each weight
real hd_dirichlet_prior_lpdf(real[] theta, int m, int weights_start, int[,] w_o, int[,] w_u){

  int n_tot; // number of variances in the split
  int weights[m];
  // int m; // number of weights in the split (2 = dual split, 3 = triple split etc.)
  real alpha; // the dirichlet parameter
  real over_part;
  real under_part;

  // if only two weights (that will sum to 1) we have a uniform prior and the log density is 0
  if (m == 2){
    return(0);
  }

  for (i in 1:m){
    weights[i] = weights_start+i-1;
  }

  n_tot = sum(w_u[weights[1],]); // how many variance parameters are involved in this split

  over_part = 0;

  alpha = get_dirichlet_parameter(m);

  for (ind in 1:m){
    over_part += log(sum(exp(theta[get_indexes(w_o[weights[ind],])])));
  }

  // the part under the fraction is the same for all weights in the multisplit,
  // so we can use the same row in w_u for all
  under_part = log(sum(exp(theta[get_indexes(w_u[weights[1],])])));

  return(
    (alpha-1) * ( over_part - m*under_part )
  );


}


// evaluates all priors that has HD priors in a tree (run several times for multiple trees)
// includes total variance
// includes jacobians!
// theta are all log variances involved in this prior tree (may have more than one tree)
real hd_prior_joint_lpdf(real[] theta, int use_likelihood, int[] which_pc, int[,] w_o, int[,] w_u, int n_splits_tree, int[] row_index_hd_pr, vector knots, matrix[] prior_coeffs, int n_knots, matrix totvar_prior_info, int indV){

  real logdens;
  real logV;
  //int indV;
  int pc_prior_ind; // index to keep track of the pc priors

  // if we do not have a HD prior, theta has length 0, and we return 0
  if (num_elements(theta) == 0){
    return(0);
  }

  logV = log(sum(exp(theta)));
  logdens = 0;

  // prior on the total variance
  if (use_likelihood == 0 && totvar_prior_info[indV,1] == 1){ // must have proper variance when we sample from prior
    logdens += normal_lpdf(logV | 0, 1);
  } else {
    logdens += choose_prior_lpdf(logV | totvar_prior_info[indV,1], to_array_1d(totvar_prior_info[indV,2:3]));
  }

  pc_prior_ind = 1;
  // prior on the weights
  for (indw in 1:n_splits_tree){
    if (which_pc[indw] == 1){ // pc prior on weight
      // jacobian from logit(w) to w is included in the function
      logdens += hd_pc_prior_lpdf(theta | row_index_hd_pr[indw], pc_prior_ind, w_o, knots, prior_coeffs, n_knots);
      pc_prior_ind += 1;
      // logdens += hd_pc_prior_lpdf(theta | row_index_hd_pr[indw+1]-(row_index_hd_pr[indw]), indw, w_o, knots, prior_coeffs, n_knots);
    } else { // dirichlet prior on weight(s)
      // second argument used to be: row_index_hd_pr[indw+1]-row_index_hd_pr[indw]-1
      logdens += hd_dirichlet_prior_lpdf(theta | row_index_hd_pr[indw+1]-(row_index_hd_pr[indw]), row_index_hd_pr[indw], w_o, w_u);
    }
  }

  logdens += calc_jac_logdet(theta, w_o, w_u, num_elements(theta), row_index_hd_pr, n_splits_tree); // jacobian from total variance/weights to log variances
  logdens += -logV; // jacobian from log total variance to total variance (not included in jacobian function)

  return(logdens);

}

// computes the prior for individual variance components
real cw_priors_lpdf(real[] theta, matrix cw_prior_info){

  real logdens;

  logdens = 0;

  // prior on the variances with CW priors
  for (indcw in 1:num_elements(theta)){
    if ((cw_prior_info[indcw,1]) > 0){
      logdens += choose_prior_lpdf(theta[indcw] | cw_prior_info[indcw, 1], to_array_1d(cw_prior_info[indcw, 2:3]));
    }
  }

  // no jacobian needed, as the priors are for log variance already

  return(logdens);

}

// computes the whole joint prior for the variance components in the model
real joint_prior_lpdf(real[] theta, int[] which_pc, int[] which_hd, int[,] w_o, int[,] w_u, int use_likelihood, int[] n_splits_each_tree, int[,] row_index_hd_pr, vector knots, matrix[] prior_coeffs, int n_knots, matrix totvar_prior_info, matrix cw_prior_info, int[] which_cw, int n_totvar, int[,] which_theta_in_hd){

  real logdens;

  logdens = 0;

  // we may have more than one tree in the prior
  if (num_elements(which_hd) > 0){
    for (indV in 1:n_totvar){

      logdens += hd_prior_joint_lpdf(theta[get_indexes(which_theta_in_hd[indV,])] | use_likelihood, which_pc, w_o[,get_indexes(which_theta_in_hd[indV,])], w_u[,get_indexes(which_theta_in_hd[indV,])], n_splits_each_tree[indV], get_indexes(row_index_hd_pr[indV,]), knots, prior_coeffs, n_knots, totvar_prior_info, indV);

    }
  }

  logdens += cw_priors_lpdf(theta[which_cw] | cw_prior_info[which_cw,]);

  return(logdens);

}

  "

return(my_text)

}


make_r_script <- function(){

  my_text <- "




# Code showing how to fit the random intercept model with the provided stan-code,
## using the object from 'make_prior'


# The prior for variance parameters can be included automatically using the existing
## functions included in the stan-files, or be coded by the user (may increase speed).
## The priors for any non-variance parameters (e.g. correlations) must be specified
## by the user in the stan-code.
## The model itself must be implemented, only support for the prior is included.
## The random intercept model is included as an example.

library(makemyprior)

ex_form <- y ~ mc(a)

ex_data <- list(a = rep(1:10, each = 10))
set.seed(1); ex_data$y <- 3 + rnorm(10, 0, 1)[ex_data$a] + rnorm(100, 0, 1)

prior <- make_prior(ex_form, ex_data,
                    prior = list(
                      tree = \"s1 = (a, eps)\",
                      w = list(s1 = list(prior = \"pc0\", param = 0.25)),
                      V = list(s1 = list(prior = \"pc\", param = c(3, 0.05)))
                    ))

# compile custom stan-code
library(rstan)
ex_mod <- stan_model(\"my_stan_code/main_file.stan\", auto_write = TRUE)

# Make the whole data object required, see stan-code for explanation on what the do.
# Note that all data specified in the data block in the stan code must be provided, and
# be on the correct format.

## this function extracts and computes the necessary details that the stan code needs
## for the prior just specified (depends on the tree structure, so must be re-run every time the
## prior is updated)
## the function is not documented, and uses some (undocumented) internal functions from makemyprior
get_data_for_hd_prior <- function(prior_obj){

  var_pr_no <- function(prior_names){
    pr_names <- c(\"jeffreys\", \"pc0\", \"invgam\", \"hc\", \"\")
    pr_num <- c()
    for (i in 1:length(prior_names)){
      pr_num[i] <- which(pr_names == prior_names[i])
    }
    pr_num[pr_num == 5] <- 0
    return(pr_num)
  }

  res <- list()

  # here we decide which variance parameter goes where in the theta-vector to stan
  app_to_stan_id <- data.frame(
    stan_id = 1:length(prior_obj$data$random),
    app_id = prior_obj$node_data$orig_nodedata$id,
    app_name = prior_obj$node_data$orig_nodedata$label
  )

  # store information from prior_obj
  res$n_hd <- sum(prior_obj$node_data$orig_nodedata$status == \"attached\") # how many variance parameters involved in hd prior
  res$n_pc <- if (length(prior_obj$prior_data$weights) > 0) sum(sapply(prior_obj$prior_data$weights, function(x) x$prior == \"pc\")) else 0

  # which components are included in a tree and have HD priors
  res$which_hd <- as.array(c(1:(nrow(prior_obj$node_data$orig_nodedata)))[prior_obj$node_data$orig_nodedata$status == \"attached\"])
  res$which_cw <- as.array(c(1:(nrow(prior_obj$node_data$orig_nodedata)))[prior_obj$node_data$orig_nodedata$status == \"detached\"])
  res$n_cw <- length(res$which_cw)

  res$which_pc <- if (res$n_hd > 0) as.array(sapply(prior_obj$prior_data$weights, function(x) if (x$prior == \"pc\") 1 else 0)) else integer()

  res$cw_prior_numbers <- var_pr_no(sapply(prior_obj$prior_data$cw_priors, function(x) x$prior))
  res$cw_prior_parameters <- t(sapply(prior_obj$prior_data$cw_priors,
                                      function(x) {
                                        param <- x$param
                                        return(if (length(param) == 2) param else rep(param, 2))
                                      }))

  res$cw_prior_info <- cbind(res$cw_prior_numbers, res$cw_prior_parameters)

  res$totvar_prior_numbers <- if (res$n_hd > 0) var_pr_no(sapply(prior_obj$prior_data$total_variance, function(x) x$prior)) else array(0, dim = c(1, 1))
  res$totvar_prior_parameters <- if (res$n_hd > 0) {
    t(sapply(prior_obj$prior_data$total_variance,
             function(x) {
               param <- x$param
               return(if (length(param) == 2) param else rep(param, 2))
             }))
  } else array(0, dim = c(1,2))
  res$totvar_prior_info <- cbind(res$totvar_prior_numbers, res$totvar_prior_parameters)

  # n_totvar is the same as the number of separate trees in the prior
  # if no total variances, i.e. no HD priors, we have n_totvar = 1 with only 0's in totvar_prior_info
  res$n_totvar <- nrow(res$totvar_prior_info)

  # which variances are involved in which hd prior
  # length is the same as number of variance parameters, index is 0 for not involved (i.e. CW prior),
  # and then 1, 2, 3, ... for each total variance index (same order as in totvar_prior_*)
  top_nodes <- prior_obj$node_data$nodes$id[prior_obj$node_data$nodes$top_node == TRUE]
  split_nodes <- unique(prior_obj$node_data$edges$from)
  which_theta_in_hd_mat <- matrix(0, nrow = res$n_totvar, ncol = nrow(prior_obj$node_data$orig_nodedata))
  res$n_splits_each_tree <- array(0, dim = 1)
  for (ind in seq_len(res$n_totvar)){
    # does not include top nodes, because they are not in the original nodes (as they are merged nodes)
    in_tree <- app_to_stan_id$stan_id[app_to_stan_id$app_id %in% makemyprior:::nodes_in_tree(prior_obj$node_data, top_nodes[ind])]
    which_theta_in_hd_mat[ind, in_tree] <- 1
    res$n_splits_each_tree[ind] <- sum(split_nodes %in% makemyprior:::nodes_in_tree(prior_obj$node_data, top_nodes[ind]))
  }

  res$n_splits <- length(prior_obj$weight_priors)

  stopifnot(sum(res$n_splits_each_tree) == res$n_splits)

  res$which_theta_in_hd <- which_theta_in_hd_mat

  # make matrices with information on how each weight looks,
  # each row corresponds to a weight, each column to a variance
  weight_info_over <- matrix(0, nrow = length(prior_obj$node_data$edges$from), ncol = nrow(prior_obj$node_data$orig_nodedata))
  weight_info_under <- matrix(0, nrow = length(prior_obj$node_data$edges$from), ncol = nrow(prior_obj$node_data$orig_nodedata))
  res$n_rows_in_w_o_u <- nrow(weight_info_over)

  # each column corresponds to a row in weight_info_over/under, rows to trees
  # (the first is for index 0, so there is one more column than rows in w_i_o/u)
  row_index_hd_pr_mat <- array(0, dim = c(res$n_totvar, length(prior_obj$node_data$edges$from)+1))
  split_nodes <- unique(prior_obj$node_data$edges$from)
  ind3 <- 1 # number of weights (total, w and 1-w counts as two)
  hd_pr_num <- numeric()
  for (ind in seq_len(sum(res$n_splits))){
    tmp <- makemyprior:::get_leaf_node_ids_for_split(prior_obj$node_data, split_nodes[ind])
    under <- app_to_stan_id$stan_id[app_to_stan_id$app_id %in% tmp$under]
    tree_no <- which(which_theta_in_hd_mat[,under[1]] == 1)
    row_index_hd_pr_mat[tree_no, ind3] <- 1
    for (ind2 in 1:(length(tmp$over))){
      over <- app_to_stan_id$stan_id[app_to_stan_id$app_id %in% tmp$over[[ind2]]]
      weight_info_over[ind3, over] <- 1
      weight_info_under[ind3, under] <- 1
      ind3 <- ind3+1
      hd_pr_num <- c(hd_pr_num, ind)
    }
    # find which tree this weight belongs to
    # can find which row in weight_theta_in_hd_mat that is 1 for the first value in under
    row_index_hd_pr_mat[tree_no, ind3] <- 1
  }

  row_index_hd_pr <- array(0, dim = 1)

  # this loop goes over each weight (each HD prior in the prior)
  for (ind2 in seq_len(length(unique(hd_pr_num)))){
    row_index_hd_pr <- c(row_index_hd_pr, which(hd_pr_num == ind2)[length(which(hd_pr_num == ind2))])
  }

  res$row_index_hd_pr <- row_index_hd_pr_mat # which weight has which rows in the weight_info-matrices
  res$row_index_hd_pr_plot <- row_index_hd_pr + 1 # which weight has which rows in the weight_info-matrices

  res$w_o <- weight_info_over
  res$w_u <- weight_info_under

  res$n_knots <- 122 # in the code we have now this is 122, and cannot be changed by the user
  res$knots <- c(-500,
                 seq(-200, -40, length.out = 10)[1:9],
                 seq(-40, -5, length.out = 22)[1:21],
                 seq(-5, -1e-6, length.out = 30),
                 seq(1e-6, 5, length.out = 30),
                 seq(5, 40, length.out = 22)[-1],
                 seq(40, 200, length.out = 10)[-1],
                 500) # these are always the same

  tmp_array <- array(0, dim = c(res$n_pc, 122, 4))
  pc_ind <- if (res$n_hd > 0) which(sapply(prior_obj$weight_priors, function(x) length(x$knots) > 2)) else c()
  for (i in seq_len(length(pc_ind))){
    tmp_array[i,,] <- prior_obj$weight_priors[[pc_ind[i]]]$coeffs
  }
  res$prior_coeffs <- tmp_array

  return(res)

}

### always required
ex_all_data <- list()
ex_all_data$n <- 100
ex_all_data$n_g <- 10
ex_all_data$n_par <- 2
ex_all_data$indexes_group <- ex_data$a
ex_all_data$y <- ex_data$y

### must be included for the HD prior if the generic functions for HD prior are used
# adding the rest of the data needed for the stan-code
ex_all_data <- c(ex_all_data, get_data_for_hd_prior(prior))

# whether to include likelihood (1, sampling from posterior) or not (0, sampling from prior)
ex_all_data$likelihood <- 1


ex_posterior <- sampling(
  ex_mod,
  ex_all_data,
  chains = 1,
  iter = 1e4,
  seed = 1
)



# some results

samps <- rstan::extract(ex_posterior)

hist(samps$intercept, freq = FALSE)
hist(samps$theta[,1], freq = FALSE)
hist(samps$theta[,2], freq = FALSE)
hist(samps$totvar, freq = FALSE)
hist(samps$varprop, freq = FALSE)

  "

return(my_text)

}



