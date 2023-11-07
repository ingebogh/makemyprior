



functions {

// including general jacobian-function
#include include/jacobian_function_general.stan
// including the variance parameter functions
#include include/variance_priors.stan
// functions with the latent component calculation
#include include/latent_models.stan

}

data {

  int n; // number of observations
  int n_r; // number of random effects i.e. variance parameters (does NOT include the residuals)
  int n_f; // number of fixed effects (not including intercept)
  int n_hd; // number of variance parameters involved in a HD prior
  int n_cw; // number of variances with CW prior
  int n_splits; // number of splits in the model
  int<lower=0,upper=n_splits> n_pc; // number of splits with a PC prior
  int n_knots; // number of knots for splines
  int n_rows_in_w_o_u; // number of rows in w_o and w_u
  int add_res; // 1 if gaussian (and we add residuals), 0 else

  real intercept_prior[2]; // prior for intercept
  matrix[n_f, 2] fixed_priors; // priors for fixed effects

  int which_pc[n_splits]; // which of the HD priors are PC (1) and Dirichlet (0)

  int which_hd[n_hd]; // which variance parameters are involved in HD priors
  int which_cw[n_cw]; // which variances have CW priors

  // matrices keeping track of which variance parameters are involved in the weights
  int w_o[n_rows_in_w_o_u, n_r+add_res];
  int w_u[n_rows_in_w_o_u, n_r+add_res];

  int n_totvar; // number of total variances (for now only one is possible, the code is not made for more than 1)
  matrix[n_totvar,3] totvar_prior_info; // info on the total variance
  int which_theta_in_hd[n_totvar, n_r+add_res];
  int n_splits_each_tree[n_totvar]; // number of splits in each tree

  int row_index_hd_pr[n_totvar, n_rows_in_w_o_u+1]; // which rows in w_o and w_u belongs to each split

  matrix[n_r+add_res,3] cw_prior_info; // info on individual priors on variance components

  int use_intercept; // if intercept is included in the model or not

  int effect_sizes[n_r]; // the size of each effect, indicating how much of each vector and matrix will be used (the rest is just ignored)
  int effect_index_start[n_r]; // index of where each effect start in the random effect vector
  int indexes[n,n_r]; // indexes for each of the random effects

  vector[n_knots] knots; // the knots for the splines
  matrix[n_knots,4] prior_coeffs[n_pc]; // coefficients for the splines

  matrix[n,n_f] fixed_effects; // the fixed effects

  real y_cont[n]; // continuous observations
  int y_disc[n]; // discrete observations
  int Ntrials[n]; // number of trials for binomial observations

  int use_index_matrix[n_r]; // whether we use the index matrix of each effect or not
  matrix[n, n] mega_matrix[n_r]; // matrices for structured effects

  int<lower=0> effect_type[n_r]; // which latent model for each of the random effects

  real scaling_factor[n_r]; // scaling factors (only used for besag)

  // indexes for besag effect
  int n_besag;
  int besag_ind1[n_besag];
  int besag_ind2[n_besag];

  int<lower=0,upper=3> likelihood; // which likelihood we use (0 = only prior, 1 = Gaussian, 2 = binomial, 3 = Poisson)

}

parameters{

  real intercept[use_intercept];
  real theta[n_r+add_res]; // log variances, residual variance comes in last
  vector[n_f] coeff; // fixed effect coefficients

  vector[sum(effect_sizes)] random_effects; // the random effects, as a N(0,1)-vector

}


transformed parameters{

  vector[n] eta;

  eta = rep_vector(0, n);

  if (use_intercept) eta =+ rep_vector(intercept[1], n);

  if (n_f > 0) eta += fixed_effects * coeff;

  if (n_r > 0 && likelihood != 0){
    for (ind in 1:n_r){ // adding one effect to the linear predictor at a time
      eta += latent_model_scaling(effect_type[ind], random_effects[(1+effect_index_start[ind]):(effect_sizes[ind]+effect_index_start[ind])], theta[ind], mega_matrix[ind][1:effect_sizes[ind], 1:effect_sizes[ind]])[indexes[,ind]];
    }
  }

}




model {

  if (use_intercept == 1) target += normal_lpdf(intercept | intercept_prior[1], intercept_prior[2]);

  if (n_f > 0) target += normal_lpdf(coeff | to_array_1d(fixed_priors[,1]), to_array_1d(fixed_priors[,2]));

  target += joint_prior_lpdf(theta | which_pc, which_hd, w_o, w_u, likelihood, n_splits_each_tree, row_index_hd_pr, knots, prior_coeffs, n_knots, totvar_prior_info, cw_prior_info, which_cw, n_totvar, which_theta_in_hd);

  // each random effect is handled in the transformed parameters block
  // if not using likelihood, we do not sample the random effects
  if (n_r > 0) {
    if (likelihood != 0){
      for (ind in 1:n_r){ // latent components
        target += latent_model_lpdf(random_effects[(1+effect_index_start[ind]):(effect_sizes[ind]+effect_index_start[ind])] | effect_type[ind], scaling_factor[ind], besag_ind1, besag_ind2);
      }
    } else { // do not sample from random effects if we don't want the posterior
      target += normal_lpdf(random_effects | 0, 1);
    }
  }

  if (likelihood == 1){
    target += normal_lpdf(y_cont | eta, sqrt(exp(theta[n_r+1])));
  } else if (likelihood == 2){
    target += binomial_logit_lpmf(y_disc | Ntrials, eta);
  } else if (likelihood == 3){
    target += poisson_log_lpmf(y_disc | eta);
  } else {
    target += 0; // no likelihood, i.e., sample from prior
  }

}


generated quantities {

  vector[sum(effect_sizes)] effects; // vector with the effects on the correct scale

  if (n_r > 0){
    for (ind in 1:n_r){ // adding one and one effect to the linear predictor

      effects[(1+effect_index_start[ind]):(effect_sizes[ind]+effect_index_start[ind])] = latent_model_scaling(effect_type[ind], random_effects[(1+effect_index_start[ind]):(effect_sizes[ind]+effect_index_start[ind])], theta[ind], mega_matrix[ind][1:effect_sizes[ind], 1:effect_sizes[ind]]);

    }
  }

}













