


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

// a function that calculates the jacobian for any structure

// it is not necessary to know if the prior on each weight is PC or Dirichlet, as we only want the expression for the weight

// this works for one tree

// theta = all log variances in this tree (all thetas from the hd_prior_joint_lpdf function)
// w_o = matrix with variances over fraction in each weight
// w_u = matrix with variances under fraction in each weight
// n = number of variance components involved in the tree
// row_index_hd_pr = which rows in w_o and w_u that goes to each weight
real calc_jac_logdet(real[] theta, int[,] w_o, int[,] w_u, int n, int[] row_index_hd_pr, int n_splits){
//matrix calc_jac_logdet2(real[] theta, int n, int[] row_index_hd_pr){
  
  matrix[n, n] jac;
  // int n_splits; // the number of splits (the number of weights is n-1)
  int i_row; // index for row
  real under; // value under fraction bar for each weight
  real over; // value over fraction bar for each weight
  int bool_arg1;
  
  // n_splits = num_elements(row_index_hd_pr)-1; // row_index_hd_pr has one more elements than the number of weights
  
  jac = rep_matrix(0, n, n); // initializing
  
  // first row in jacobian is for the log total variance
  for (i in 1:n){
    jac[1,i] = exp(theta[i]);
  }
  
  i_row = 2;
  
  // for each split in the tree
  for (i_split in 1:n_splits){
    
    //print("i_split: ", i_split);
    
    //print("theas under: ", get_indexes(w_u[row_index_hd_pr[i_split],]));
    
    // value under fraction bar is the same for all weights in a split
    under = ( sum(exp(theta[get_indexes(w_u[row_index_hd_pr[i_split],])])) )^2;
    
    over = 0; // if the parameter is not in the weight, the derivative is 0
    
    // looping through the weights in this split (not the last one, since that is 1 minus the others)
    // (must use -2 since the next index in row_index_hd_pr gives the index of the next weight (-1 here), and we do not need the last one (-1 more))
    for (ind in row_index_hd_pr[i_split]:(row_index_hd_pr[i_split+1]-2)){
      
      //print("ind: ", ind);
      
      //print("hei: ", get_indexes(w_u[ind,]));
      // looping through each parameter (i.e. each column of jacobian) relevant for this weight to calculate the part above the fraction bar
      for (i_col in get_indexes(w_u[ind,])){
        
        //print("i_col: ", i_col);
        
        // testing if the variance corresponding to the index i_col is involved above the fraction bar of this weight
        bool_arg1 = 0;
        for (i in 1:sum(w_o[ind,])){
          if (i_col == get_indexes(w_o[ind,])[i]){
            bool_arg1 = 1;
          }
        }
        
        // if this variance parameter is above the fraction bar and we are differentiating wrt to it
        if (bool_arg1 == 1){
          
          //print("thetas over: ", get_indexes2(w_o[ind,], w_u[ind,]), ", ", i_col);
          
          over = ( sum(exp(theta[get_indexes2(w_o[ind,], w_u[ind,])])) ) * exp(theta[i_col]);
          
        } else { // if this variance parameter is not above the fraction bar and we are differentiating wrt to it
          
          //print("thetas over: ", get_indexes(w_o[ind,]), ", ", i_col);
          
          over = -( sum(exp(theta[get_indexes(w_o[ind,])])) ) * exp(theta[i_col]);
          
        }
        
        jac[i_row, i_col] = over/under;
        
      }
      
      //print("----------------------------------i_row: ", i_row);
      
      i_row += 1;
      
    }
    
  }
  
  return(log(fabs(determinant(jac))));
  
}









