
// jeffreys prior on variance
// input is log variance
real jeffreys_logvar_lpdf(real x){
  return(0);
}

// inverse gamma on variance (gamma on precision)
// input is log variance
// shape and rate, corresponds to inverse gamma with shape and scale
real invgamma_logvar_lpdf(real x, real shape, real scale){
  return(
    -lgamma(shape) + shape*log(scale) -shape*x - scale/(exp(x))
  );
}

// half cauchy on standard deviation
// input is log variance
real halfcauchy_logvar_lpdf(real x, real scale){
  return(
    x/2 - log(pi()) - log(scale + exp(x)/scale)
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

// THIS FUNCTION IS INCLUDED IN ANOTHER FILE
// returns the indexes where w_row is 1
// used to pick correct thetas for the logit weights
// int[] get_indexes(int[] w_row){
//
//   int res[sum(w_row)];
//   int ind;
//
//   ind = 1;
//   for (i in 1:num_elements(w_row)){
//     if (w_row[i] == 1){
//       res[ind] = i;
//       ind += 1;
//     }
//   }
//
//   return(res);
//
// }

// calculates the joint PC prior
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

  //print("vAluE of m is:", m);
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
    // logdens += normal_lpdf(logV | 0, 1);
    logdens += choose_prior_lpdf(logV | 5, to_array_1d(totvar_prior_info[indV,2:3]));
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
  // logdens += calc_jac_logdet(theta);
  logdens += -logV; // jacobian from log total variance to total variance

  return(logdens);

}

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


real joint_prior_lpdf(real[] theta, int[] which_pc, int[] which_hd, int[,] w_o, int[,] w_u, int use_likelihood, int[] n_splits_each_tree, int[,] row_index_hd_pr, vector knots, matrix[] prior_coeffs, int n_knots, matrix totvar_prior_info, matrix cw_prior_info, int[] which_cw, int n_totvar, int[,] which_theta_in_hd){

  real logdens;

  logdens = 0;

  // we may have more than one tree in the prior
  if (num_elements(which_hd) > 0){
    for (indV in 1:n_totvar){

      logdens += hd_prior_joint_lpdf(theta[get_indexes(which_theta_in_hd[indV,])] | use_likelihood, which_pc, w_o[,get_indexes(which_theta_in_hd[indV,])], w_u[,get_indexes(which_theta_in_hd[indV,])], n_splits_each_tree[indV], get_indexes(row_index_hd_pr[indV,]), knots, prior_coeffs, n_knots, totvar_prior_info, indV);

    }
  }

  //logdens += hd_prior_joint_lpdf(theta[which_hd] | use_likelihood, which_pc, w_o[,which_hd], w_u[,which_hd], n_splits, row_index_hd_pr, knots, prior_coeffs, n_knots, totvar_prior_info, indV);

  logdens += cw_priors_lpdf(theta[which_cw] | cw_prior_info[which_cw,]);

  return(logdens);

}























