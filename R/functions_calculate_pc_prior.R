


# must first calculate the (perhaps singular) covariance matrices for all effects, of the same side as the data and residuals
# then for each split, we reduce the the matrices to something that is non-singular, so we can calculate the PC prior
# for each split higher up, we just use the basemodels of lower splits to make new covariance matrices



# assume that we need matrices of lower dimension that can lead to non-intrinsic priors
# reduce matrices by removing linearly dependent things, and mapping something
# mats are the two covariance matrices of size nxn
reduce_matrices <- function(mats){

  # this function is only used if we have two matrices
  mat1 <- mats[[1]]
  mat2 <- mats[[2]]

  if (nrow(mat1) != nrow(mat2) || ncol(mat1) != ncol(mat2)) stop("Matrices of different sizes!", call. = FALSE)
  if (!isSymmetric(mat1, tol = 1e-8) || !isSymmetric(mat2, tol = 1e-8)) stop("Matrices not symmetric!", call. = FALSE)

  n <- nrow(mat1)

  qr1_rank <- compute_rank(mat1)
  qr2_rank <- compute_rank(mat2)
  # qr1 <- base::qr(mat1)
  # qr2 <- base::qr(mat2)

  # then they have full rank and do not need transformation
  # TODO: can test this outside the function
  if (qr1_rank == n && qr2_rank == n) return(mats)
  # if (qr1$rank == n && qr2$rank == n) return(mats)

  ## reduce the matrices so they are non-singular using eigenvalue decomposition

  # if the matrices are equal, we must do something clever
  if (sum(mat1 != mat2) == 0){
    stop("The covariance matrices for this split are equal, this is wrong. Reformulate you model.", call. = FALSE)
  }

  mat <- mat1 + mat2

  # make a reduction transformation matrix A
  eig_decomp <- base::eigen(mat)
  not_zero <- abs(eig_decomp$values) > 1e-10
  A <- t(eig_decomp$vectors[,not_zero])

  mat1_red <- A %*% mat1 %*% t(A)
  mat2_red <- A %*% mat2 %*% t(A)

  rank1_red <- compute_rank(mat1_red)
  rank2_red <- compute_rank(mat2_red)
  # rank1_red <- base::qr(mat1_red)$rank
  # rank2_red <- base::qr(mat2_red)$rank

  # if (rank1_red != qr1_rank){
  #   mat1_red <- transf_matrix(mat1_red)
  #   browser()
  # }
  # if (rank2_red != qr2_rank){
  #   mat2_red <- transf_matrix(mat2_red)
  #   browser()
  # }

  return(list(as.matrix(mat1_red), as.matrix(mat2_red)))

}

# not sure what this is doing
transf_matrix <- function(X){

  n <- nrow(X)
  I <- diag(n)
  J <- matrix(1, n, n)

  return(
    (I-J/n) %*% X %*% t(I-J/n)
  )

}


# the index in the matrix_data-list of a given node with id node_id
get_matrix_data_index <- function(matrix_data, node_id) return(which(sapply(matrix_data, function(x) x$id) == node_id))

# multiplies the matrices in the list mats with the weights in weights, and sums up
make_basemodel_covmat <- function(matrix_data, weights){
  mat <- 0
  for (i in 1:length(matrix_data)){
    mat <- mat + weights[i]*matrix_data[[i]]
  }
  return(mat)
}

# calculates the values for the basemodel for each node in a split
# for pc prior, the order is: value for the basemodel-node, 1-value for the basemodel_node
get_basemodel_value <- function(prior_info){

  if (prior_info$prior == "dirichlet"){
    basemodel <- rep(1/prior_info$no_children, prior_info$no_children)
  } else { # pc prior
    # # the order of the basemodel vector must be correct
    # basemodel1 <- prior_info$param$basemodel
    # basemodel_comp <- prior_info$param$basemodel_node
    # # gives an error if basemodel is not the first node, since the variable "basemodel" does not esist
    # basemodel <- if (prior_info$children[1] == basemodel_comp) c(basemodel1, 1-basemodel1) else c(1-basemodel1, basemodel)
    basemodel <- c(prior_info$param$basemodel, 1-prior_info$param$basemodel)
  }

  return(basemodel)

}

# calculate the pc prior for the whole tree
# matrix_data is data-input and is not changed anywhere
calculate_pc_prior <- function(node_data, prior_data, matrix_data, gui = FALSE){

  # if we have no edges, we do not have any splits either
  if (nrow(node_data$edges) == 0) return(list())

  # # get the matrices for all leaf nodes and splits in one list
  # # TODO: make an object with the matrices so we do not need to calculate them every time
  split_matrix_data <- calc_all_covariance_matrices(node_data, prior_data, matrix_data)

  logit_ws <- c(-500,
                seq(-200, -40, length.out = 10)[1:9],
                seq(-40, -5, length.out = 22)[1:21],
                seq(-5, -1e-6, length.out = 30),
                seq(1e-6, 5, length.out = 30),
                seq(5, 40, length.out = 22)[-1],
                seq(40, 200, length.out = 10)[-1],
                500)

  # calculate pc prior for each split node
  split_nodes <- get_split_ids(node_data)

  prior_list <- list() # list for storing prior information

  # is_intrinsic_vec <- c()

  for (ind in 1:length(split_nodes)){

    node_id <- split_nodes[ind]
    split_name <- get_node_name(node_data, node_id)
    base_name <- get_node_name(node_data, prior_data[[get_prior_number(prior_data, node_id)]]$children[1])

    prior_info_tmp <- prior_data[[get_prior_number(prior_data, node_id)]]

    if (prior_info_tmp$prior == "dirichlet"){

      # no pc prior here, just store info on correct format

      prior <- list(
        node_id = node_id,
        n_knots = 2,
        knots = c(0, 0),
        coeffs = array(0, dim = c(2, 4))
      )

      prior_list[[get_prior_number(prior_data, node_id)]] <- prior

    } else { # pc prior

      # only used when basemodel is not 0 or 1
      conc_param <- prior_info_tmp$param$concentration

      children <- prior_info_tmp$children
      above_id <- prior_info_tmp$param$above
      not_above_id <- children[children != above_id]

      # base_mat is the matrix the matrix corresponding to the node that is registered as the basemodel_node in prior_data
      # alt_mat is the other matrix
      base_mat <- as.matrix(split_matrix_data[[get_matrix_data_index(split_matrix_data, above_id)]]$mat)
      alt_mat <- as.matrix(split_matrix_data[[get_matrix_data_index(split_matrix_data, not_above_id)]]$mat)

      stopifnot(all(dim(base_mat) == dim(alt_mat))) # if they have different dimensions, something is wrong

      n <- nrow(base_mat)
      median_to_base <- prior_info_tmp$param$median
      amount_to_base <- prior_info_tmp$param$basemodel # amount of variance that goes to the basemodel node (0, 1 or median)

      qr_base <- list(rank = compute_rank(base_mat))
      qr_alt <- list(rank = compute_rank(alt_mat))
      # qr_base <- base::qr(base_mat)
      # qr_alt <- base::qr(alt_mat)

      ### three cases, definition of base_mat and alt_mat:
      ### base_mat = matrix of node above in the weight, which is not present in basemodel with PC0
      ### alt_mat = matrix of node that is not above in the weight, which is not present in basemodel with PC1
        # base_mat = matrix of above_node, alt_mat = the other matrix. Independent of what basemodel is!!
      ### both_mat = base_mal*w0 + alt_mat*(1-w0)
      # 1: both base_mat and alt_mat non-singular -> ok
      # 2: one of base_mat and alt_mat are singular
        # a: non-singular matrix is basemodel -> ok
        # b: singular matrix is basemodel -> intrinsic
        # c: basemodel is a combination -> ok
      # 3: both base_mat and alt_mat are singular
        # a: basemodel has lower rank than alternative model -> intrinsic
        # b: basemodel has the same rank as both_mat (i.e. same or higher rank than alternative model) OR
            # basemodel is a combination, -> reduce -> ok (we get situation 1)
      if (qr_base$rank == n && qr_alt$rank == n){ # 1 (both are non-singular)

        spline_coeffs <- calc_prior_spline(
          logit_ws = logit_ws,
          intrinsic = FALSE, # not intrinsic
          basemodel = amount_to_base,
          median_to_base = median_to_base,
          base_mat = base_mat,
          alt_mat = alt_mat,
          conc_param = conc_param,
          splitname = split_name,
          basename = base_name,
          gui = gui
        )

      } else if (min(qr_base$rank, qr_alt$rank) < n && max(qr_base$rank, qr_alt$rank) == n){ # 2 (one is singular, one is non-singular)

        if (amount_to_base == 1 && qr_base$rank < n){ # base_mat is basemodel and singular
          is_intrinsic <- TRUE
        } else if (amount_to_base == 0 && qr_alt$rank < n){ # alt_mat is basemodel and singular
          is_intrinsic <- TRUE
        } else if (amount_to_base == 1 && qr_base$rank == n){ # base_mat is basemodel and non-singular
          is_intrinsic <- FALSE
        } else if (amount_to_base == 0 && qr_alt$rank == n){ # alt_mat is basemodel and non-singular
          is_intrinsic <- FALSE
        } else if (!(amount_to_base %in% c(0,1))){ # basemodel is a combination of one singular and one non-singular
          is_intrinsic <- FALSE
        } else {
          stop("Something went wrong in the computation of the PC prior. Change prior or reformulate your model.", call. = FALSE)
        }

        spline_coeffs <- calc_prior_spline(
          logit_ws = logit_ws,
          intrinsic = is_intrinsic,
          basemodel = amount_to_base,
          median_to_base = median_to_base,
          base_mat = base_mat,
          alt_mat = alt_mat,
          conc_param = conc_param,
          splitname = split_name,
          basename = base_name,
          gui = gui
        )

      } else if (max(qr_base$rank, qr_alt$rank) < n) { # 3 (both are singular)

        both_mat <- base_mat + alt_mat
        qr_both <- list(rank = compute_rank(both_mat)) #base::qr(both_mat)

        # the rank of base_mat, alt_mat and both_mat will not change with the reduction,
        # so checking their ranks first to avoids some matrix-operations

        if (amount_to_base == 1 && qr_base$rank < qr_both$rank){ # base_mat is basemodel and has less rank than both_mat
          spline_coeffs <- calc_prior_spline(
            logit_ws = logit_ws,
            intrinsic = TRUE, # intrinsic
            basemodel = amount_to_base,
            median_to_base = median_to_base,
            base_mat = base_mat,
            alt_mat = alt_mat,
            conc_param = conc_param,
            splitname = split_name,
            basename = base_name,
            gui = gui
          )
        } else if (amount_to_base == 0 && qr_alt$rank < qr_both$rank){ # alt_mat is basemodel and has less rank than both_mat
          spline_coeffs <- calc_prior_spline(
            logit_ws = logit_ws,
            intrinsic = TRUE, # intrinsic
            basemodel = amount_to_base,
            median_to_base = median_to_base,
            base_mat = base_mat,
            alt_mat = alt_mat,
            conc_param = conc_param,
            splitname = split_name,
            basename = base_name,
            gui = gui
          )
        } else { # must reduce the dimension of the matrices

          mat_red <- reduce_matrices(list(base_mat, alt_mat))
          base_mat_red <- mat_red[[1]]
          alt_mat_red <- mat_red[[2]]

          base_mat_red <- base_mat_red #/typical_variance(base_mat_red)
          alt_mat_red <- alt_mat_red #/typical_variance(alt_mat_red)

          spline_coeffs <- calc_prior_spline(
            logit_ws = logit_ws,
            intrinsic = FALSE, # not intrinsic
            basemodel = amount_to_base,
            median_to_base = median_to_base,
            base_mat = base_mat_red,
            alt_mat = alt_mat_red,
            conc_param = conc_param,
            splitname = split_name,
            basename = base_name,
            gui = gui
          )
        }
      } else {
        stop("Something went wrong in the computation of the PC prior. Change prior or reformulate your model.", call. = FALSE)
      }

      # if we are in the gui, we need to know if the base-model is intrinsic or not to print a message to the user
      # if (!gui) {
      #   is_intrinsic_vec <- c(is_intrinsic_vec, spline_coeffs$intrinsic)
      # }
      # spline_coeffs <- spline_coeffs$spline_coeffs

      prior <- list(
        node_id = node_id,
        n_knots = length(spline_coeffs$knots),
        knots = spline_coeffs$knots,
        coeffs = spline_coeffs$C
      )

      prior_list[[get_prior_number(prior_data, node_id)]] <- prior

    }

  }

  # if (gui) return(list(prior_list = prior_list, is_intrinsic = is_intrinsic_vec)) else return(prior_list)

  return(prior_list)

}

# calculates the spline coefficients object for a given prior
calc_prior_spline <- function(logit_ws, intrinsic, basemodel, median_to_base, base_mat, alt_mat, conc_param, splitname, basename, gui = FALSE){

  if (!intrinsic){
    y_vals <- pc_dual_logit_stable(
      logitW = logit_ws,
      w0 = basemodel,
      shape = 1,
      median = median_to_base,
      wLeft = 0.5,
      SigA = base_mat,
      SigB = alt_mat,
      conc_param = conc_param,
      splitname = splitname
      # SigA = alt_mat,
      # SigB = base_mat
    )
  } else {

    if (!gui){
      check <- NULL
      if (basemodel == 0 && median_to_base > 0.25){
        check <- paste0("larger than 0.25 for the '", splitname, "' split\n")
      } else if (basemodel == 1 && median_to_base < 0.75){
        check <- paste0("smaller than 0.75 for the '", splitname, "' split\n")
      }
      if (!is.null(check)){
        warning(paste0(
          "The covariance matrix of the base model for '",
          splitname,
          "' is singular. \nThis means that the prior does not change after the median is ",
          "further than 0.25 away from the base model.\n",
          "Here, that means that a median ", check,
          "will not change the prior.",
          "You have used a prior that corresponds to median = '", median_to_base,
          "' where '", basename, "' is the base model \n(median is '",
          median_to_base, "' away from '", basename,
          "', which may be the opposite of that you specified)."
        ), call. = FALSE)
      }
    }

    y_vals <- prior_intrinsic_base(
      lW = logit_ws,
      shape = 1,
      median = median_to_base,
      w0 = basemodel
    )

  }

  spline_coeffs <- getSplinePrior(logit_ws, y_vals)
  spline_coeffs$C[c(1, 121, 122), 3:4] <- 0 # make the spline linear in the tails

  # return(list(spline_coeffs = spline_coeffs, intrinsic = intrinsic))
  return(spline_coeffs)

}


# calculate the matrices for the split nodes and add the leaf node matrices
calc_all_covariance_matrices <- function(node_data, prior_data, matrix_data){

  # calculate pc prior for each split node
  split_nodes <- get_split_ids(node_data)
  # sort by level to simplify
  split_nodes <- split_nodes[order(node_data$nodes$level[node_data$nodes$id %in% split_nodes], decreasing = TRUE)]

  new_matrix_data <- list()

  for (ind in 1:length(matrix_data)){
    new_matrix_data[[ind]] <- list(id = get_node_id(node_data, names(matrix_data)[ind]), mat = matrix_data[[ind]])
  }

  # add the matrices for the split nodes
  for (ind in 1:length(split_nodes)){

    node_id <- split_nodes[ind]

    prior_data_tmp <- prior_data[[get_prior_number(prior_data, node_id)]]

    child_nodes <- prior_data_tmp$children

    mats <- lapply(child_nodes, function(x) new_matrix_data[[get_matrix_data_index(new_matrix_data, x)]])

    basemodel_weights <- get_basemodel_value(prior_data_tmp)

    new_matrix_data <- c(new_matrix_data, list(list(id = node_id, mat = make_basemodel_covmat(lapply(mats, function(x) x$mat), basemodel_weights))))

  }

  return(new_matrix_data)

}

## x = x-values, either on regular scale (logitscale = FALSE), or on logitscale (logitscale = TRUE)
## prior = list with prior information
## logitscale = logical, choose between regular and logit scale
eval_spline_prior <- function(x, prior, logitscale = TRUE){

  eval_spline1 <- function(logitw, knots, coeffs, n_knots){

    # Search for correct cubic polynomial
    low = 1
    high = n_knots
    while(high-low > 1){
      med = (high+low)/2
      if(logitw < knots[med]){
        high = med
      } else{
        low = med
      }
    }
    idd = low
    dx = logitw-knots[idd]
    val = coeffs[idd, 1]
    val = val + coeffs[idd, 2]*dx
    val = val + coeffs[idd, 3]*dx*dx
    val = val + coeffs[idd, 4]*dx*dx*dx

    return(val)
  }

  eval_spline2 <- function(lw) {
    exp(eval_spline1(lw, prior$knots, prior$coeffs, prior$n_knots))
  }

  eval_spline3 <- function(lw) {
    if (length(lw) > 1) {
      sapply(lw, function(t) eval_spline2(t))
    } else {
      return(eval_spline2(lw))
    }
  }

  logit <- function(x) log(x/(1-x))

  if (logitscale){
    return(
      eval_spline3(x)
    )
  } else {
    return(
      eval_spline3(logit(x)) * (x*(1-x))^(-1)
    )
  }

}







