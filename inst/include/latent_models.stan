

// scaling/constraints etc
vector latent_model_scaling(int effect_type, vector eff, real theta, matrix mega_matrix){

  vector[num_elements(eff)] eta_part; // this part of the linear predictor
  vector[num_elements(eff)] constr_vec;

  if (effect_type == 1){ // iid no constraints
    eta_part = eff;
  } else if (effect_type == 2){ // iid sum-to-zero constraint
    eta_part = eff - mean(eff);
  } else if (effect_type == 3){ // rw2 no constraints
    eta_part = eff;
  } else if (effect_type == 4){ // rw2 sum-to-zero constraint
    // eta_part = dot_self(eff[1:(num_elements(eff)-2)] - 2*eff[2:(num_elements(eff)-1)] + eff[3:num_elements(eff)]);
    eta_part = eff - mean(eff);
  } else if (effect_type == 5){ // rw2 sum to zero and linear constraint
    // eta_part = dot_self(eff[1:(num_elements(eff)-2)] - 2*eff[2:(num_elements(eff)-1)] + eff[3:num_elements(eff)]);
    for (i in 1:num_elements(eff)) constr_vec[i] = i;
    constr_vec = constr_vec - mean(constr_vec);
    eta_part = eff - mean(eff) - sum(eff .* constr_vec) / (sum(constr_vec .* constr_vec)) * constr_vec;
  } else if (effect_type == 6){ // effect is multiplied by (cholesky) lower triangle matrix, no constraints
    eta_part = mega_matrix * eff;
  } else if (effect_type == 7){ // effect is multiplied by (cholesky) lower triangle matrix, sum-to-zero constraint
    eta_part = mega_matrix * eff;
    eta_part = eta_part - mean(eta_part);
  } else if (effect_type == 8){ // besag, no constraints
    eta_part = eff;
  } else if (effect_type == 9){ // besag, sum-to-zero
    eta_part = eff - mean(eff);
  } else if (effect_type == 10){ // rw1, no constraints
    eta_part = eff;
  } else if (effect_type == 11){ // rw1, sum-to-zero
    eta_part = eff - mean(eff);
  }

  eta_part = eta_part * sqrt(exp(theta));
  return(eta_part);

}


real latent_model_lpdf(vector eff, int effect_type, real scaling_factor, int[] besag_ind1, int[] besag_ind2){

  int n;
  vector[num_elements(eff)] constr_vec;
  real res = 0;
  n = num_elements(eff);

  if (effect_type == 3 || effect_type == 4 || effect_type == 5){ // rw2
    res += -0.5 * scaling_factor * dot_self(eff[1:(n-2)] - 2*eff[2:(n-1)] + eff[3:n]);
    if (effect_type >= 4){
      res += normal_lpdf(sum(eff) | 0, 1*sqrt(n)); // (soft) sum to zero
    }
    if (effect_type == 5){
      for (i in 1:num_elements(eff)) constr_vec[i] = i;
      res += normal_lpdf(dot_product(eff, constr_vec) | 0, 1*sqrt(1*dot_self(constr_vec-mean(constr_vec)))); // (soft) linear sum to zero
    }
  } else if (effect_type == 8 || effect_type == 9){ // besag
    res += -0.5 * scaling_factor * dot_self(eff[besag_ind1] - eff[besag_ind2]);
    if (effect_type == 9) res += normal_lpdf(sum(eff) | 0, 1*sqrt(n)); // (soft) sum to zero
  } else if (effect_type == 10 || effect_type == 11) { // rw1
    res += -0.5 * scaling_factor * dot_self(eff[1:(n-1)]-eff[2:n]);
    if (effect_type == 11) res += normal_lpdf(sum(eff) | 0, 1*sqrt(n)); // (soft) sum to zero
  } else { // 1, 2, 6, 7 (iid and structured effect)
    res += normal_lpdf(eff | 0, 1);
  }

  return(res);

}






