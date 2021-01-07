
# this script contains the functions that are needed to make the prior-object that can be sent to either
# shiny or directly to inference



# functions for specifying the model and running inference


# functions that may be useful: attr, deparse, terms, match.call, match.args, match.fun, update.formula, all.vars
# label <- lapply(as.list(substitute(list(...)))[-1], eval) # finds all arguments in "..."-part of function call

possible_models <- function() return(c("iid", "linear", "rw2", "generic0"))

# make a function similar to f in inla, called f2 (at least for now)
#' Define latent component
#'
#' Function for defining a latent component for the HD prior package. All model components must be specified, and if an HD
#' prior is used, this is specified later
#' @param ... Name of the component (short names is an advantage as they are used in the app), no default (MUST be provided)
#' @param model Type of model, default is "iid" (see list of models SOMEWHERE)
#' @param prior The prior choice for this component if it should have a CW prior (else the prior will be specified when running the model later).
#' Prior should be specified as a list with prior name and parameters (see list of priors SOMEWHERE), e.g. (list(prior = "pc", param = c(3, 0.05))),
#' which is the default.
#' @param constr Sum-to-zero constraints on component (default TRUE)
#' @param constr_lin Linear constraints on component (only for rw2, default = TRUE)
#' @param A Covariance for this component (taken to be identity if not specified)
#' @keywords formula
#' @export
#' @return Should not return anything to the user
#' @examples
#' vignette("create_prior", package = "priorconstruction")
#'
f2 <- function(..., model = "iid", prior = NULL, constr = NULL, constr_lin = FALSE, A = NULL, Cmatrix = NULL){

  label <- as.character(as.list(substitute(list(...)))[[2]])

  if (label == "eps") stop("You cannot use the name 'eps', this is reserved for residuals. Choose another name.")

  # make an object (a list) with the random effect information
  res <- list(label = label, model = model)

  if (!(model %in% possible_models())) stop("Not a valid model type!")

  # if (model == "iid"){
  #   res$A <- A
  #   res$constr <- constr
  #   if (!missing(prior)){
  #     res$prior <- prior
  #   } else {
  #     res$prior <- list(prior = "mw") # if no prior provided, we use an HD prior
  #   }
  # } else if (model == "linear"){
  #   if (is.null(prior)) res$prior <- list(mean.linear = 0, sd.linear = sqrt(1000))
  # } else { 
  #   stop("feil")
  # }
  
  if (model == "linear"){
    if (is.null(prior)) res$prior <- list(mean.linear = 0, sd.linear = sqrt(1000))
  } else { # random models
    # res$constr <- constr
    if (!missing(prior)){
      res$prior <- prior
    } else {
      res$prior <- list(prior = "mw") # if no prior provided, we use an HD prior
    }
    # in case there is some extra data needed for this model type we need later, we have to store this in the object
    res$inference_data <- list()
  }
  
  # linear is FALSE always
  if (!is.null(constr)){
    if (model != "rw2"){
      res$constr <- constr
      res$constr_lin <- FALSE
    } else {
      res$constr <- constr
      res$constr_lin <- if (is.na(constr_lin)) TRUE else constr_lin
    }
  } else {
    if (model != "iid"){
      res$constr <- TRUE
      res$constr_lin <- if (model == "rw2") TRUE else FALSE
    } else {
      res$constr <- TRUE
      res$constr_lin <- FALSE
    }
  }
  
  # adding additional arguments
  if (!missing(A)) res$A <- A
  
  # checking if the user is missing inputs that are required for each model type
  req <- f2_required_input(model)
    
  # for inla formula (only storing names as strings)
  # TODO: make more efficient somehow
  if ("Cmatrix" %in% req){
    if (is.null(substitute(Cmatrix))) stop(paste0("Missing argument(s) ", paste0(req, collapse = ", "), " for model ", model))
    res$inference_data$Cmatrix <- deparse(substitute(Cmatrix))
  }
  
  return(res)

}

# returns which arguments are required for each model type
f2_required_input <- function(model_type){
  
  if (model_type == "generic0") return(c("Cmatrix"))
  
  return(NULL)
  
}


#' Creating a prior object
#'
#' Create a prior object that can either be sent to the app
#' or directly to inference in INLA or Stan.
#' @param formula A formula object, using the function ``f2``.
#' @param data The data used in the model, as a ``data.frame`` or ``list``. All elements must have the same length.
#' @param family A string indicating the likelihood family. ``gaussian`` is the default.
#' @param link Link function. Default is linear link function.
#' @param hd_prior A named list with four elements: ``tree`` is a string indicating the prior tree for all or some of
#' the random model components, ``w`` is a named list with priors for each split, ``V`` is
#' named a list with priors for
#' the total variance(s) (this list must have names for each tree, where the name is the name of the
#' top split),
#' and ``cw`` are component-wise (CW) priors for each.
#' If the lists with priors are not named properly, or is not valid,
#' the function provides default dirichlet
#' priors. If the prior tree is not providided or is not valid, a default prior is given to all components
#' in the tree, in the case where all priors for the components are (correcly) specified the default is a CW prior
#' on all components. If some but not all model components are specified in the tree, the non-specified ones will
#' get CW priors (if which prior is not specified, a PC(3, 0.05)-prior is given to the standard deviation
#' of the component).
#' Do not need to have w and V, and you can provide only one of the if you prefer, but needs tree if the other objects are going to be used.
#' Note that the tree structure overrides any CW priors specified, if they are included in the string for the tree.
#' If no tree is specified, the components with specified priors will be independent, and the others will get a default HD prior
#' @param residual_prior A list indicating the prior to be used on the residuals (if Gaussian likelihood). If
#' this info is provided in the tree structure, it will be ignored. If not provided, and all model components
#' have CW priors, it defaults to PC-SD_0(3, 0.05). If not provided and a prior tree is provided, this will be
#' attached at the top level with rhinkage towards and 75\% variance provided to residuals.
#' @param mean.intercept Mean of the intercept prior (default is 0)
#' @param sd.intercept Standard deviation of the intercept prior (1000)
#' @keywords prior
#' @return
#' @examples
#' vignette("create_prior", package = "priorconstruction")
#'
create_prior <- function(formula, data, family = "gaussian", # link,
                         hd_prior = list(), residual_prior = list(),
                         mean.intercept = 0, sd.intercept = 1000){

  # get information on each component in the model individually
  # (mw-priors are marked with "mw" and dealt with later)
  prior <- decompose_formula(formula, family, data)
  # do not need to find the fixed priors here, can do that when making data-object for inla/stan
  
  # add info on residual prior
  if (family == "gaussian"){
    prior$random_effects$eps$label <- "eps" # always give residuals name "eps", user cannot use this in the data
    tmp_pr <- 0
    # if the prior is specified and the user has eps in the tree, the prior is MW
    if (length(residual_prior) == 0 || is_node_in_tree_string(hd_prior$tree, "eps") == 2){
      tmp_pr <- list(prior = "mw")
    } else {
      tmp_pr <- residual_prior
    }
    prior$random_effects$eps$prior <- tmp_pr
  } else {
    stop("Not implemented yet!")
  }

  # NOTE: tree structure should override prior choices!
  tree <- make_valid_tree(names(prior$random_effects), hd_prior$tree, names(hd_prior$cw))

  # check if the user has specified any CW priors
  for (ind in seq_len(length(hd_prior$cw))){
    # if the name of this prior is in the data, but not in the tree, we give it a CW prior
    nodename <- names(hd_prior$cw)[ind]
    if (nodename %in% names(prior$random_effects)){
        if (is_node_in_tree_string(tree, nodename) == 0) {
          # not in the tree, so we keep the prior and add the node to the tree as a CW node
          tree <- paste0(tree, ";(", nodename, ")")
        } else if (is_node_in_tree_string(tree, nodename) == 1) { # in tree, do not need to add
          # this node in the tree, do nothing
        } else { # is in split, removing the prior
          hd_prior$cw[[ind]] <- list()
        }
    } else {
      warning("You have provided a prior to a non-existing component, ignoring it.", call. = FALSE)
    }
  }

  # remove empty parts of cw prior list
  hd_prior$cw <- hd_prior$cw[!sapply(hd_prior$cw, function(x) length(x) == 0)]

  mw_prior_comps <- names(prior$random_effects)[sapply(prior$random_effects, function(x) x$prior$prior == "mw")]

  node_data <- make_hd_prior_tree(names(prior$random_effects), tree)
  # checking if this is the prior the user asked for, if not we give a warning

  tree_used <- make_hd_prior_string(node_data)
  if (!compare_tree_strings(tree, tree_used, names(prior$random_effects))){
    warning(paste0("Was not able to create the tree requested, using this instead: ",
                   tree_used
                   ), call. = FALSE
    )
    # make node_data with the object that is used
    node_data <- make_hd_prior_tree(names(prior$random_effects), tree_used)
  }

  old_new_names <- node_data$old_new_names
  node_data <- node_data[-4]

  # if the user did not specify eps in the tree, had all other nodes in one tree, and did not specify residual prior,
  # we put the residuals on the top with a PC_0(0.25)-prior with shrinkage to the residuals
  # in this case, the residual node is detached
  # TODO: if not correctly specified tree, use default (always), do not have all these exceptions
  # only have reisduals for gaussian likelihood
  # check if the tree the user sent in had the residuals as a CW prior
  if (family == "gaussian" && !is.null(hd_prior$tree) && !grepl("(eps)", gsub(" ", "", hd_prior$tree))){
  #if (family == "gaussian" && !grepl("(eps)", gsub(" ", "", hd_prior$tree))){
    eps_id <- get_node_id(node_data, "eps")
    top_id <- get_top_nodes(node_data)
    if (length(top_id) == 1 && are_detached(node_data, list(current_node_id = eps_id)) && length(residual_prior) == 0){
      node_data_tmp <- node_data
      # update so the residual node has the right properties
      node_data_tmp$nodes$level <- node_data_tmp$nodes$level+1 # move all nodes down
      node_data_tmp$nodes$level[node_data$nodes$label == "eps"] <- 1 # give the residual node level 1
      old_top_node <- node_data_tmp$nodes[node_data_tmp$nodes$top_node == 1,]
      node_data_tmp$nodes$top_node <- 0
      node_data_tmp$nodes$status <- "attached"

      # make new split node
      new_node_id <- max(node_data$nodes$id)+1
      new_node <- old_top_node
      new_node$id <- new_node_id
      new_node$label <- paste0(new_node$label, "_eps")
      new_node$top_node <- 1
      new_node$level <- 0
      node_data_tmp$nodes <- rbind(node_data_tmp$nodes, new_node)

      # add edges
      new_edges <- data.frame(from = c(new_node_id, new_node_id), to = c(old_top_node$id, eps_id))
      node_data_tmp$edges <- rbind(node_data_tmp$edges, new_edges)

      # update orig_nodedata (top_node and status)
      node_data_tmp$orig_nodedata <- node_data_tmp$nodes[node_data_tmp$nodes$id %in% node_data_tmp$orig_nodedata$id,]

      node_data <- node_data_tmp

      # add prior-info
      new_prior_w <- list(list(prior = "pc0", param = 0.25))
      names(new_prior_w) <- new_node$label
      hd_prior$w <- c(hd_prior$w, new_prior_w)

      # add this node to old_new_name
      old_new_names <- rbind(old_new_names, data.frame(old = new_node$label, new = new_node$label))
      # get a new tree that is used, since we added residuals
      tree_used <- make_hd_prior_string(node_data)

    }
  }

  prior_data <- make_hd_prior_default(node_data, family) # this is the default choices!

  # fill in the priors specified by the user, if any
  # requires that the prior-tree is provided, if not we do not know the structure of the prior
  if (!is.null(hd_prior$tree)){
    if (!is.null(hd_prior$w)){
      # go through prior_data and see if the user has specified a prior
      for (ind in 1:length(prior_data$weights)){
        # check if user has specified a prior for this split node
        which_match <- which(names(hd_prior$w) %in% old_new_names$old[old_new_names$new == prior_data$weight[[ind]]$name])
        if (length(which_match) == 1){
          basemodel_node <- get_basemodel_node_name(prior_data$weights[[ind]], tree, old_new_names)
          prior_data$weights[[ind]] <- make_valid_w_prior(hd_prior$w[[which_match]], prior_data$weights[[ind]], get_node_id(node_data, basemodel_node))
        }
      }
    }
    if (!is.null(hd_prior$V)){
      # TODO: must check that the list is actually named as the tree
      # go through prior_data and see if the user has specified a prior
      for (ind in 1:length(prior_data$total_variance)){
        which_match <- which(names(hd_prior$V) %in% old_new_names$old[old_new_names$new == prior_data$total_variance[[ind]]$name])
        if (length(which_match) == 1){
          prior_data$total_variance[[ind]] <- make_valid_V_prior(hd_prior$V[[which_match]], prior_data$total_variance[[ind]])
        }
      }
    }
  }
  if (!is.null(hd_prior$cw)){
    # TODO: must check that the list is actually named as the tree
    # go through prior_data and see if the user has specified a prior
    for (ind in 1:length(prior_data$cw_priors)){
      which_match <- which(names(hd_prior$cw) %in% old_new_names$old[old_new_names$new == prior_data$cw_priors[[ind]]$name])
      if (length(which_match) == 1){
        prior_data$cw_priors[[ind]] <- make_valid_cw_prior(hd_prior$cw[[which_match]], prior_data$cw_priors[[ind]])
      }
    }
  }

  covmats <- lapply(prior$random_effects, function(x) x$A)
  covmats <- covmats[!sapply(covmats, is.null)]
  
  # arguments needed for the inference in inla
  inference_data <- list()
  model_data <- list()
  for (ind in seq_len(length(prior$random_effects))){
    if (prior$random_effects[[ind]]$label != "eps"){
      model_data[[ind]] <- c(
        list(id = get_node_id(node_data, prior$random_effects[[ind]]$label)), 
        prior$random_effects[[ind]]
      )
      for (ind2 in seq_along(prior$random_effects[[ind]]$inference_data)){
        tmp <- list(data[[prior$random_effects[[ind]]$inference_data[[ind2]]]])
        names(tmp) <- prior$random_effects[[ind]]$inference_data[[ind2]]
        inference_data <- c(inference_data, tmp)
      }
    }
  }
  
  pr <- make_valid_prior_object(
    #data = data,
    data = list(fixed = lapply(prior$fixed_effects, function(x) x$data), random = lapply(prior$random_effects, function(x) x$data)),
    args = list(
      node_data = node_data,
      prior_data = prior_data,
      covmats = covmats,
      response_name = prior$response_name,
      use_intercept = prior$use_intercept,
      family = family,
      model_data = model_data
    )
  )
  

  res <- c(pr, list(response = prior$response,
                    tree = tree_used,
                    use_intercept = prior$use_intercept,
                    prior_intercept = c(mean.intercept, sd.intercept),
                    prior_fixed = t(sapply(prior$fixed_effects, function(x) x$prior)),
                    model_data = model_data,
                    formula = formula,
                    family = family,
                    inference_data = inference_data))
  #class(res) <- "hdprior"

  return(res)

}


# return:
# 0: not in tree at all
# 1: in tree as CW
# 2: in tree in split
is_node_in_tree_string <- function(tree, nodename){

  if (!is.null(tree) && grepl(nodename, tree)){ #
    if (grepl(paste0("\\(", nodename, "\\)"), tree)){ # if CW node
      return(1)
    } else if (grepl(paste0(",", nodename), tree) || grepl(paste0("\\(", nodename), tree) || grepl(paste0(nodename, ","), tree) || grepl(paste0(nodename, "\\)"), tree)) {
      return(2)
    }
  }

  return(0)

}


get_basemodel_node_name <- function(split_info, tree, old_new_names){
  
  oldname <- old_new_names$old[split_info$name == old_new_names$new]
  splits <- strsplit(strsplit(tree, ";")[[1]], "=")
  base_old <- sub("\\(", "", strsplit(splits[[which(sapply(splits, function(x) x[1]) == oldname)]][2], ",")[[1]][1])
  base_new <- old_new_names$new[base_old == old_new_names$old]
  
  return(base_new)
  
}

make_valid_w_prior <- function(new_prior, old_prior, basemodel_id){

  tmp_prior <- old_prior
  # check that we do not put pc prior on multi-split and that the prior-name is valid (else, default choice)
  if (tmp_prior$no_children > 2) return(old_prior)
  if (new_prior$prior == "dirichlet") return(old_prior) # this will be the same, dirichlet is default
  if (new_prior$prior %in% c("pc0", "pc1", "pcM")){
    
    median_value <- new_prior$param
    
    # if the node that is not defined as basemodel-node by the "new" tree is the basemodel,
    # we turn the prior around (PC0(m) = PC1(1-m), PC1(m) = PC0(1-m) and PCM(m) = PCM(1-m))
    if (basemodel_id != tmp_prior$children[1]){
      
      if (new_prior$prior != "pcM"){
        tmp <- "pc0"
        if (new_prior$prior == "pc0") tmp <- "pc1"
        new_prior$prior <- tmp
      }
      basemodel_id <- tmp_prior$children[1]
      median_value <- 1-median_value

    }
    
    tmp_prior$prior <- "pc" #new_prior$prior
    #tmp_prior$param <- new_prior$param
    if (new_prior$prior == "pc0") {
      tmp_prior$param <- data.frame(basemodel_node = basemodel_id,
                                    basemodel = 0,
                                    median = median_value
      )
    } else if (new_prior$prior == "pc1"){
      tmp_prior$param <- data.frame(basemodel_node = basemodel_id,
                                    basemodel = 1,
                                    median = median_value
      )
    } else if (new_prior$prior == "pcM"){
      tmp_prior$param <- data.frame(basemodel_node = basemodel_id,
                                    basemodel = median_value,
                                    median = median_value
      )
    } else {
      # tmp_prior$param <- data.frame(basemodel_node = tmp_prior$children[1],
      #                               basemodel = new_prior$param[1],
      #                               median = new_prior$param[2]
      # )
      stop("Not valid prior name.")
    }
    # warning2("Ingeborg: Have not renamed the old priors yet, so we need to change the parameters a bit here.")
    return(tmp_prior)
  }

  # returning the old prior if the new is not valid
  return(old_prior)

}

make_valid_V_prior <- function(new_prior, old_prior){

  # have not renamed the priors yet, so we change the name of this prior
  # if the user just inputs "pc", it should change the prior to the internal name
  # the old name will be correct, since it is generated by an internal function
  if (new_prior$prior %in% c("pc", "pc0")) {
    new_prior$prior <- "pc0"
  }

  tmp_prior <- old_prior

  # must check that we do not put a jeffreys prior on total variance unless the default has allowed it
  if (new_prior$prior == "jeffreys" && old_prior$prior != "jeffreys") return(old_prior)

  # other than that, we can do "whatever"
  tmp_prior$prior <- new_prior$prior
  if (tmp_prior$prior == "jeffreys"){
    tmp_prior$param <- c(0, 0)
    return(tmp_prior)
  } else {
    if (tmp_prior$prior == "hc"){
      tmp_prior$param <- c(new_prior$param[1], 0)
    } else {
      tmp_prior$param <- new_prior$param
    }
    return(tmp_prior)
  }

  # returning the old prior if the new is not valid
  return(old_prior)

}

make_valid_cw_prior <- function(new_prior, old_prior){

  # have not renamed the priors yet, so we change the name of this prior
  # if the user just inputs "pc", it should change the prior to the internal name
  # the old name will be correct, since it is generated by an internal function
  if (new_prior$prior %in% c("pc", "pc0")) {
    new_prior$prior <- "pc0"
  }

  tmp_prior <- old_prior

  # other than that, we can do "whatever"
  tmp_prior$prior <- new_prior$prior
  if (tmp_prior$prior == "hc"){
    tmp_prior$param <- c(new_prior$param[1], 0)
  } else if (tmp_prior$prior == "invgam") {
    tmp_prior$param <- new_prior$param
  } else if (tmp_prior$prior == "pc0"){
    tmp_prior$param <- new_prior$param
  } else { # if the prior is not valid, we give a PC(3, 0.05)-prior
    tmp_prior$prior <- "pc0"
    tmp_prior$param <- c(3, 0.05)
  }

  return(tmp_prior)

}

# returns info on what user has chosen for each component in the formula
decompose_formula <- function(formula, family, data){

  # decompose formula, extract relevant information from each component
  terms <- terms(formula)

  fixed_effects <- list()
  random_effects <- list()

  term_labs <- attr(terms, "term.labels")
  intercept <- as.logical(attr(terms, "intercept"))
  # not sure if this is usable for more than one response variable (and not sure if we ever have more than one)
  response_name <- if (attr(terms(as.formula(formula)), which = "response")) all.vars(formula)[1] else NULL
  response <- eval(parse(text = response_name), data)
  n <- length(response)

  # separate fixed and random effects
  for (ind in 1:length(term_labs)){

    # f2-function used
    if (substr(term_labs[ind], 1, 3) == "f2("){

      # run f2-function to extract info from formula
      tmp <- eval(parse(text = term_labs[ind]))
      if (tmp$model == "linear"){ # can be used for fixed effect
        fixed_effects[[tmp$label]] <- tmp
        fixed_effects[[tmp$label]]$data <- eval(parse(text = tmp$label), data)
      # } else if (tmp$model == "iid"){
      #   random_effects[[tmp$label]] <- tmp
      #   random_effects[[tmp$label]]$data <- eval(parse(text = tmp$label), data)
      #   random_effects[[tmp$label]]$A <- tmp$A #if (is.null(tmp$A)) diag(n) else tmp$A
      } else {
        random_effects[[tmp$label]] <- tmp
        random_effects[[tmp$label]]$data <- eval(parse(text = tmp$label), data)
        # random_effects[[tmp$label]]$A <- tmp$A #if (is.null(tmp$A)) diag(n) else tmp$A
      }

    } else { # fixed effect
      # Note that fixed effects are added to stan in a separate data frame (can maybe use same as random effects in inla)
      fixed_effects[[term_labs[ind]]] <- list(label = term_labs[ind],
                                              model = "linear",
                                              data = eval(parse(text = term_labs[ind]), data))
      fixed_effects[[term_labs[ind]]]$prior <- list(mean.linear = 0, sd.linear = sqrt(1000))
    }

  }

  # make data-frame for the fixed effects
  if (length(fixed_effects) > 0){
    fixed_effects_data <- as.data.frame(matrix(unlist(sapply(fixed_effects, function(x) x$data), use.names = FALSE), nrow = n))
    names(fixed_effects_data) <- sapply(fixed_effects, function(x) x$label)
    fixed_effects_data <- as.list(fixed_effects_data)
  }

  # make data-frame for the random effects
  if (length(random_effects) > 0){
    random_effects_data <- as.data.frame(matrix(unlist(sapply(random_effects, function(x) x$data), use.names = FALSE), nrow = n))
    names(random_effects_data) <- sapply(random_effects, function(x) x$label)
    random_effects_data <- as.list(random_effects_data)
  }

  if (family == "gaussian") random_effects_data$eps <- 1:length(response)

  return(list(
    fixed_effects = fixed_effects,
    random_effects = random_effects,
    use_intercept = intercept,
    response = response,
    response_name = response_name,
    likelihood = family
  ))

}


# default tree structure for the labels provided (split with all model components)
default_tree <- function(effnames){
  return(
    paste0(paste0(effnames, collapse = "_"), "=(", paste0(effnames, sep = "", collapse = ","), ")", collapse = "")
  )
}

# internal function, creates a prior object that can be sent to shiny/inference
# input will include full node_data and prior_data, so we do not need to make that
make_valid_prior_object <- function(data, args){

  # if (missing(args)) args <- list() # not sure if this is necessary, but keeping it for now

  # load variables from user to send to shiny, putting them in one list
  initial_args <- list()

  initial_args$.nodenames <- as.character(args$node_data$orig_nodedata$label)
  covmats_input <- args$covmats

  random_data <- data$random
  # if we have a gaussian likelihood, we must have a covariance matrix for the residuals as well
  if (args$family == "gaussian"){
    random_data$eps <- 1:length(random_data[[1]])
  }
  covmats <- make_covariance_matrices(random_data) # making iid covariance matrices
  # covmats <- if (length(args$covmats) == 0) make_covariance_matrices(random_data) else args$covmats
  which_covmats_structured <- c()
  # if some of the covariance matrices are structured
  for (ind in seq_along(covmats_input)){
    index <- which(names(covmats) == names(covmats_input)[ind])
    new_covmat <- fix_size_covmat(covmats_input[[ind]], covmats[[index]], names(covmats)[index])
    covmats[[index]] <- new_covmat
    # tmpmat <- prior_obj$covmats[[ind]][1:res$effect_sizes[ind],1:res$effect_sizes[ind]] + 10^(-9)*diag(res$effect_sizes[ind])
    # lower_chol <- t(chol(tmpmat/typical_variance(tmpmat)))
    which_covmats_structured <- c(which_covmats_structured, index)
  }
  
  # if we have any non-iid components, we need to make a covariance matrix for them
  mods <- c(sapply(args$model_data, function(x) x$model))
  non_iids <- which(mods != "iid")
  
  for (ind in seq_along(non_iids)){
    
    if (mods[non_iids[ind]] == "rw2"){
      covmat <- make_rw2_mat(n = length(unique(data$random[[non_iids[ind]]])), type = "cov")
    } else if (mods[non_iids[ind]] %in% possible_models()) {
      covmat <- covmats[[non_iids[ind]]]
    } else {
      stop(paste0(mods[non_iids[ind]], " is not available."), call. = FALSE)
    }
    new_covmat <- fix_size_covmat(covmat, covmats[[non_iids[ind]]],  names(covmats)[non_iids[ind]])
    covmats[[non_iids[ind]]] <- new_covmat/typical_variance(new_covmat)
    which_covmats_structured <- c(which_covmats_structured, non_iids[ind])
    
  }

  initial_args$.initial_nodes <- args$node_data$nodes
  initial_args$.initial_edges <- args$node_data$edges
  initial_args$.no_comp <- nrow(args$node_data$nodes)

  stopifnot(nrow(initial_args$.no_comp) <= 16) # cannot have too many components

  initial_args$.initial_prior_w <- args$prior_data$weights
  initial_args$.initial_prior_V <- args$prior_data$total_variance
  initial_args$.initial_prior_var <- args$prior_data$cw_priors

  initial_args$.response_name <- args$response_name
  initial_args$.use_intercept <- args$use_intercept
  initial_args$.family <- args$family
  initial_args$.model_data <- args$model_data

  # adding the fixed effect info, for visualization only, is not involved in the prior

  initial_args$.fixef_names <- names(data$fixed)

  initial_args$.model_eq <-
    paste("Model equation: ", initial_args$.response_name, " ~ ",
          if (initial_args$.use_intercept) expression("mu + ") else "",
          if (!is.null(initial_args$.fixef_names)) paste0(paste(initial_args$.fixef_names, sep = "", collapse = " + "), " + ") else "",
          paste(args$node_data$orig_nodedata$label,
                sep = "", collapse = " + "),
          sep = "")

  val <- list(
    prior_data = args$prior_data,
    node_data = args$node_data,
    weight_priors = calculate_pc_prior(args$node_data, args$prior_data$weights, covmats),
    .initial_args = initial_args,
    run_inference_directly = FALSE # not sure if we need this
  )

  prior_obj <- c(val, list(data = data,
                           covmats = covmats,
                           which_covmats_structured = which_covmats_structured))
  return(prior_obj)

}


# makes a list with the splits
make_split_list <- function(splits){

  tmp <- strsplit(splits, "=")

  # remove parantheses
  tmp <- lapply(tmp, function(x) {x[2] <- gsub("\\(", "", x[2]); x})
  tmp <- lapply(tmp, function(x) {x[2] <- gsub("\\)", "", x[2]); x})

  return(
    lapply(tmp, function(x) list(parent = x[1], children = strsplit(x[2], ",")[[1]]))
  )

}

# finds the level of the node with name node_name
find_node_level <- function(split_list, node_name){

  tmp_ind <- sapply(split_list, function(x) sum(x$children == node_name))
  lvl <- 0
  while (sum(tmp_ind) == 1){
    lvl <- lvl + 1
    tmp_ind <- sapply(split_list, function(x) sum(x$children == split_list[[which(tmp_ind == 1)]]$parent))
  }

  return(lvl)

}

make_edges <- function(split_list, node_data){

  edges <- data.frame()
  # must match names and id's for the edges
  for (ind in 1:length(split_list)){
    tmp_from <- node_data$id[node_data$label %in% split_list[[ind]]$parent]
    tmp_to <- node_data$id[node_data$label %in% split_list[[ind]]$children]
    edges <- rbind(edges, data.frame(from = tmp_from, to = tmp_to))
  }

  return(edges)

}

# function that makes a valid tree structure
# check in create_prior if this function crashes, in that case set the default prior
make_valid_tree <- function(effnames, tree, cw_comps) {

  # removing all spaces
  tree <- gsub(" ", "", tree)

  # ok to have specified a tree without anything, then we check later if this is default MW or CW
  if (is.null(tree) || length(tree) == 0 || tree == "") {
    if (length(cw_comps) > 0){
      tree <- default_tree(effnames[!(effnames %in% cw_comps)])
      tree <- paste0(tree, ";", paste0("(", cw_comps, ")", sep = "", collapse = ";"))
      warning("Did not find a full tree, using default tree structure for the components without specified priors.", call. = FALSE)
    } else {
      tree <- default_tree(effnames)
      warning("Did not find a tree, using default tree structure instead.", call. = FALSE)
    }
  } else {

    # checking how much ; there are
    no_dividers <- nchar(as.character(tree)) - nchar(gsub(";", "", tree))
    no_splits <- nchar(as.character(tree)) - nchar(gsub("=", "", tree))

    split_names <- sapply(strsplit(strsplit(tree, ";")[[1]][grepl("=", strsplit(tree, ";")[[1]])], "="), function(x) x[1])
    if (sum(split_names %in% effnames) > 0) {
      warning("Cannot use the same name on splits as the random effects! Using default tree structure instead.",
              call. = FALSE)
      tree <- default_tree(effnames)
    }

    effnames_tree <- strsplit(gsub("\\)", "", gsub("\\(", "", sapply(strsplit(strsplit(tree, ";")[[1]], "="), function(x) x[2]))), ",")
    if (sum(table(unlist(effnames_tree)) > 2) > 0){
      warning("One or more effects or splits are found in several splits. Using default tree instead.", call. = FALSE)
      tree <- default_tree(effnames)
    }

    mw_effnames <- unlist(effnames_tree)[unlist(effnames_tree) %in% effnames]
    cw_effnames <- effnames[!(effnames %in% mw_effnames)]
    # some CW components
    if (length(mw_effnames) < length(effnames)){
      for (ind in 1:length(cw_effnames)){
        if (!grepl(paste0("(", cw_effnames[ind], ")"), tree)) {
          tree <- paste0(tree, ";", paste0("(", cw_effnames, ")", collapse = ";", sep = ""))
        }
      }
    }

    if (no_dividers == no_splits && length(cw_effnames) == 0) {
      alternative_tree <- paste0(strsplit(tree, ";")[[1]], sep = "", collapse = ";")
      warning(paste("Assuming you mean", alternative_tree, "."), call. = FALSE)
      tree <- alternative_tree
    }
  }

  return(tree)

}


# makes a string of the prior tree from node_data
make_hd_prior_string <- function(node_data){

  effnames <- as.character(node_data$orig_nodedata$label)
  splits <- get_split_ids(node_data)

  tmp <- c()
  comps_added <- c()
  for (ind in splits){
    tmp <- c(
      tmp,
      paste0(get_node_name(node_data, ind), " = (", paste0(get_node_name(node_data, get_children(node_data, list(current_node_id = ind))),
                                                         sep = "", collapse = ","), ")")
    )
    comps_added <- c(comps_added, get_node_name(node_data, get_children(node_data, list(current_node_id = ind))))
  }
  tmp <- paste0(tmp, sep = "", collapse = "; ")
  if (sum(!(effnames %in% comps_added)) > 0){ # if we have CW components
    # tmp <- paste0(tmp, ";",  paste0("(", effnames[!(effnames %in% comps_added)], ")", sep = "", collapse = ";"))
    if (tmp == "") tmp <- c()
    tmp <- paste0(c(tmp, paste0("(", effnames[!(effnames %in% comps_added)], ")", sep = "", collapse = "; ")), sep = "", collapse = "; ")
  }

  return(tmp)
  # return(paste0(tmp, sep = "", collapse = "; "))

}

# compares two hd prior tree strings
compare_tree_strings <- function(tree_in, tree_out, effnames){

  # assign each effect a number
  eff_name_to_number_orig <- data.frame(name = effnames, number = 1:length(effnames))
  eff_name_to_number_in <- eff_name_to_number_out <- eff_name_to_number_orig

  # remove spaces
  tree_in <- gsub(" ", "", tree_in)
  tree_out <- gsub(" ", "", tree_out)
  # divide the tree into the different splits
  splits_in <- strsplit(tree_in, ";")[[1]]
  splits_out <- strsplit(tree_out, ";")[[1]]

  # WOPS: need to make sure the tree in contains the residual split!
  if (length(splits_in) != length(splits_out)) return(FALSE) # does not have the same number of splits, so they cannot be equal

  names_in <- names_out <- list()
  for (ind in 1:length(splits_in)){
    tmp_in <- collect_split_info(splits_in[[ind]], eff_name_to_number_in)
    eff_name_to_number_in <- tmp_in$eff_num
    names_in[[ind]] <- tmp_in[-1]
    tmp_out <- collect_split_info(splits_out[[ind]], eff_name_to_number_out)
    eff_name_to_number_out <- tmp_out$eff_num
    names_out[[ind]] <- tmp_out[-1]
  }

  # make new name for each split node, for child each node we must check if it is among the original nodes or not,
  # if it is not we must find the children of that node (and so on)
  split_nums_in <- sapply(names_in, function(x) x$num_name)
  split_nums_out <- sapply(names_out, function(x) x$num_name)

  # if these lists have the same elements (order does not matter), the trees are equal
  list_in <- lapply(split_nums_in, function(x) find_leaf_nodes(names_in, x, eff_name_to_number_orig))
  list_out <- lapply(split_nums_out, function(x) find_leaf_nodes(names_out, x, eff_name_to_number_orig))

  match_func <- function(x, y){
    if (length(y) == length(x)){
      if (!sum(sort(y) != sort(x))){
        return(TRUE)
      }
    }
    return(FALSE)
  }

  # loop through one list, find equal element in the other, and remove it
  for (i in 1:length(list_in)){
    match_io <- sapply(list_out, match_func, list_in[[i]])
    if (!sum(match_io)){ # no match, then this is wrong
      return(FALSE)
    } else {
      list_out <- list_out[!match_io]
    }
  }

  if (length(list_out) == 0) return(TRUE)

  return(FALSE)

}

# finding leaf nodes for each split
find_leaf_nodes <- function(split_list, id, eff_num){

  split <- split_list[[which(sapply(split_list, function(x) x$num_name == id))]]
  new_name <- c()
  for (ind in 1:length(split$num_children)){
    if (split$num_children[ind] %in% eff_num$number){
      new_name <- c(new_name, split$num_children[ind])
    } else {
      new_name <- c(new_name, find_leaf_nodes(split_list, split$num_children[ind], eff_num))
    }
  }

  return(new_name)

}

# collect info about this split (name and children)
collect_split_info <- function(split, eff_num){

  split_name <- strsplit(split, "=")[[1]][1]
  children <- strsplit(gsub(")", "", strsplit(split, "\\(")[[1]][2]), ",")[[1]]
  for (i in 1:length(c(split_name, children))){
    if (!(c(split_name, children)[i] %in% eff_num$name)){
      eff_num <- rbind(eff_num, data.frame(name = c(split_name, children)[i], number = nrow(eff_num)+1))
    }
  }

  num_children <- sapply(children, function(x) eff_num$number[eff_num$name == x])

  return(list(eff_num = eff_num, old_name = split_name, children = children,
              num_name = eff_num$num[eff_num$name == split_name], num_children = num_children))

}


# make valid node_data-object
# NOTE: This function always gets something that is valid!
make_hd_prior_tree <- function(effnames, hd_prior_tree, family = "gaussian", mw = TRUE){
  
  # remove spaces
  hd_prior_tree <- gsub(x = hd_prior_tree, pattern = " ", replacement = "")

  n_par <- length(effnames)

  # we may have only CW priors, must check this (in this case, the expression has no equal-signs)
  if (!grepl("=", hd_prior_tree)){
    n_splits <- 0
    splits <- ""
    split_names <- ""
    split_list <- list()
    status_ad <- rep("detached", n_par)

    level <- rep(0, n_par)

  } else {
    # identify splits in the tree

    splits <- strsplit(hd_prior_tree, split = ";", fixed = TRUE)[[1]]
    # remove CW nodes that show up here
    splits <- splits[grepl("=", splits)]
    n_splits <- length(splits)
    stopifnot(n_splits < n_par)
    split_names <- as.character(sapply(splits, function(x) strsplit(x, "=")[[1]][1]))
    split_list <- make_split_list(splits) # list with splits and their children

    # find which nodes are attached and detached
    status_ad <- rep(NA, n_par)
    status_ad[!(effnames %in% unlist(sapply(split_list, function(x) x$children)))] <- "detached"
    status_ad[(effnames %in% unlist(sapply(split_list, function(x) x$children)))] <- "attached"

    level <- sapply(effnames, function(x) find_node_level(split_list, x))
  }

  orig_nodedata <- data.frame(
    id = 1:n_par,
    label = effnames,
    level = level,
    status = status_ad,
    top_node = 0 # none of the original nodes can be top nodes (not defined as top node in tree with one node)
  )

  if (grepl("=", hd_prior_tree)){
    unused_ids <- c(1:(n_par+n_splits))[!(c(1:(n_par+n_splits)) %in% orig_nodedata$id)]
    nodes_split <- data.frame(
      id = unused_ids,
      label = split_names,
      level = sapply(split_names, function(x) find_node_level(split_list, x)),
      status = "attached", # all split nodes are attached, else it is not a split node...
      top_node = as.numeric(sapply(split_names, function(x) find_node_level(split_list, x)) == 0)
    )

    nodes <- rbind(orig_nodedata, nodes_split)
    nodes <- nodes[order(nodes$id),]

    # do not need width and dashes info outside app, as it is only for the graphical visualization, the info is in the prior
    edges <- make_edges(split_list, nodes)

  } else {
    nodes <- orig_nodedata
    edges <- data.frame()
  }

  # the node_data-object for the tree
  node_data <- list(nodes = nodes, edges = edges, orig_nodedata = orig_nodedata)

  if (grepl("=", hd_prior_tree)){
    # adding id to split_list
    split_list <- lapply(split_list, function(x) {x$id <- node_data$nodes$id[node_data$nodes$label == x$parent]; x})
    split_list <- lapply(split_list, function(x) {x$id_children <- node_data$nodes$id[node_data$nodes$label %in% x$children]; x})

    # fix names of nodes
    node_data <- update_node_labels(node_data)
  }

  node_data$nodes$status <- as.character(node_data$nodes$status)
  node_data$orig_nodedata$status <- as.character(node_data$orig_nodedata$status)

  # save old and new names of split nodes
  node_data$old_new_names <- data.frame(
    old = nodes$label,
    new = node_data$nodes$label
  )

  return(node_data)

}


# make valid default prior_data-object
make_hd_prior_default <- function(node_data, family){

  prior_data <- list(weights = list(), total_variance = list(), cw_priors = list())

  # CW priors
  for (ind in 1:nrow(node_data$orig_nodedata)){
    if (node_data$orig_nodedata$status[ind] == "attached"){
      prior_data$cw_priors[[ind]] <- list(id = node_data$orig_nodedata$id[ind],
                                          name = node_data$orig_nodedata$label[ind],
                                          prior = "", param = c(0, 0))
    } else {
      prior_data$cw_priors[[ind]] <- list(id = node_data$orig_nodedata$id[ind],
                                          name = node_data$orig_nodedata$label[ind],
                                          prior = "pc0", param = c(3, 0.05))
    }
  }

  # weight priors
  for (ind in get_split_ids(node_data)){
    prior_data$weights <- c(prior_data$weights,
                            list(list(id = ind,
                                      name = get_node_name(node_data, ind),
                                      prior = "dirichlet",
                                      param = get_diri_param(no_of_children(node_data, list(current_node_id = ind))),
                                      children = get_children(node_data, list(current_node_id = ind)),
                                      no_children = no_of_children(node_data, list(current_node_id = ind))
                            ))
    )
  }

  # total variance(s)
  top_nodes <- get_top_nodes(node_data)
  for (ind in top_nodes){
    prior_data$total_variance <- c(prior_data$total_variance,
                                   list(list(id = ind,
                                             name = get_node_name(node_data, ind),
                                             prior = default_totvar(family, length(top_nodes) + sum(node_data$nodes$status == "detached"))$prior,
                                             param = default_totvar(family, length(top_nodes) + sum(node_data$nodes$status == "detached"))$param
                                   ))
    )
  }

  return(prior_data)

}

# sets default total variance prior
default_totvar <- function(family, n_comp){
  if (family == "gaussian" && n_comp == 1){
    return(list(prior = "jeffreys", param = c(0, 0)))
  } else {
    return(list(prior = "pc0", param = c(3, 0.05)))
  }
}

#' Graphical prior construction
#'
#' This functions opens a shiny app where the specified prior can be seen, and changed.
#' Returns an invisible object,
#'
#' @param prior An object from ``create_prior``, which contains all information about the data, model, and prior (if provided, can
#' also be unspecified from ``create_prior``). Can also be an object from
#' ``inference_stan``. In this case, the samples from the inference will be removed if the prior is changed, but not if
#' the prior is unchanged when the app is closed.
#' @param guide Logical, whether to open the guide directly when the app is started. Default is ``FALSE``. The guide
#' can be opened in the app at any time.
#' #' @keywords prior, shiny
#' @export
#' @return Returns an object that can be sent to ``inference_stan``. Can also be sent to ``run_shiny`` again.
#' @examples
#' vignette("create_prior", package = "priorconstruction")

run_shiny <- function(prior, guide = FALSE){

  inference_results <- list()
  if (length(prior) == 3) {
    inference_results <- list(stan = prior$stan, stan_data = prior$stan_data)
    prior <- prior$prior
  }

  initial_args <- prior$.initial_args

  initial_args$.indata$covmats <- prior$covmats

  initial_args$.guide <- guide

  .GlobalEnv$.initial_args <- initial_args

  # remove input parameters when closing
  on.exit(rm(list = c(".initial_args"), envir = .GlobalEnv))

  # run application
  shinyjs::useShinyjs()
  val <- shiny::runApp(shinyApp(ui, server), quiet = TRUE)

  prior$.initial_args <- update_initial_args(initial_args, val)

  res <- c(val, list(data = prior$data,
                     covmats = prior$covmats,
                     which_covmats_structured = prior$which_covmats_structured,
                     .initial_args = prior$.initial_args,
                     #.initial_args = update_initial_args(initial_args, val),
                     response = prior$response,
                     tree = make_hd_prior_string(val$node_data),
                     use_intercept = prior$use_intercept,
                     prior_intercept = prior$prior_intercept,
                     prior_fixed = prior$prior_fixed,
                     model_data = prior$model_data,
                     formula = prior$formula,
                     family = prior$family
                     ))

  # class(res) <- "hdprior"

  if (res$run_inference_directly){
    res <- inference_inla(prior, use_likelihood = TRUE)
    return(res)
  }

  if (length(inference_results) > 0 && identical(prior, res)){
    res <- c(list(prior = res), inference_results)
  }

  # return(res)
  return(invisible(res))

}

# makes new initial arguments, as they are now supposed to be the same as the output from shiny
update_initial_args <- function(old_initial_args, prior_obj){

  new_initial_args <- old_initial_args

  # these are the ones who change
  new_initial_args$.initial_nodes <- prior_obj$node_data$nodes
  new_initial_args$.initial_edges <- prior_obj$node_data$edges[,names(prior_obj$node_data$edges) %in% c("from", "to")]
  new_initial_args$.initial_prior_w <- prior_obj$prior_data$weights
  new_initial_args$.initial_prior_V <- prior_obj$prior_data$total_variance
  new_initial_args$.initial_prior_var <- prior_obj$prior_data$cw_priors

  # remove guide and indata
  new_initial_args$.indata <- NULL
  new_initial_args$.guide <- NULL

  return(new_initial_args)

}


#' Run inference
#'
#' This function helps you run inference with RStan using a prior object from construct_prior
#' @param prior_obj An object from ``create_prior``, from ``run_shiny``
#' (both with and without the "Close and run" option), from ``inference_stan``, or from ``inference_inla`` (for refitting model)
#' @param use_likelihood Whether to sample from the prior only (0, can be used for e.g. debugging or to look at the priors on variance parameters when using an HD prior), or to use the likelihood and data to get the posterior (1, default).
#' @param print_prior Whether to print a text with the chosen prior or not (default TRUE)
#' @param seed Seed for random number generation, default = 1
#' @param init Initial value of the parameters. All parameters are initialized to 0 (on log-scale for variances) by default
#' @param chains The number of chains used, default = 1
#' @param ... Other arguments to be sent to ``rstan::sampling``
#' @keywords inference
#' @return A named list with a prior object (prior), a stan-object (stan) and some data stan requires (stan_data).
#' @export
#' @examples
#' vignette("create_prior", package = "priorconstruction")

inference_stan <- function(prior_obj, use_likelihood = 1, print_prior = TRUE, seed = 1, init = "0", chains = 1, ...){
  
  if (print_prior) print_prior_choice(prior_obj)

  if (length(prior_obj) == 3) prior_obj <- prior_obj$prior
  stan_data <- make_stan_data_object(prior_obj, use_likelihood)
  stan_mod <- stanmodels$full_file
  
  stan_data$stz_constraint2 <- stan_data$stz_constraint1 # TODO: whyyy?
  # stan_data$constr_vec <- rep(0, 9)

  res_stan <- sampling(
    object = stan_mod,
    data = stan_data,
    seed = seed,
    init = init,
    chains = chains,
    ...
    # seed = 1, # making the simulations reproducable
    # warmup = n_samps,
    # iter = n_samps*2,
    # init = "0", # initialize all parameters to zero
    # control = list(adapt_delta = adapt_delta),
    # chains = 1
  )

  return(list(prior = prior_obj, stan = res_stan, stan_data =
                list(w_o = stan_data$w_o,
                     w_u = stan_data$w_u,
                     row_index_hd_pr_plot = stan_data$row_index_hd_pr_plot,
                     which_hd = stan_data$which_hd,
                     constr = list(c1 = stan_data$stz_constraint1, c2 = stan_data$stz_constraint2, c_vec = stan_data$constr_vec))
              ))

}


# fixed effects can be a data.frame or a list
make_stan_data_object <- function(prior_obj, use_likelihood = 1){

  # in case y has some weird format
  y <- as.numeric(unlist(prior_obj$response))

  effects <- prior_obj$data$random
  res <- list()

  tmp <- get_node_name_id(prior_obj$node_data)

  # here we decide which variance parameter goes where in the theta-vector to stan
  # may not be necessary when the original nodes are initialized from 1, but does not hurt I guess
  app_to_stan_id <- data.frame(
    stan_id = 1:length(effects),
    app_id = tmp$id,
    app_name = tmp$name
  )

  # if gaussian data and we have "eps" in the data, it should be removed after storing info about the residual node
  if (length(effects$eps) == 0){ # for now, the likelihood is always gaussian
    # if (prior_obj$likelihood == "gaussian" && !is.null(effects$eps)){
    effects <- effects[which(names(effects) != "eps")]
  }

  res$n <- length(y)
  res$n_r <- length(effects) # this is the number of non-residual random effects
  res$n_f <- length(prior_obj$data$fixed) #ncol(res$fixed_effects)
  res$fixed_effects <- as.data.frame(prior_obj$data$fixed)
  if (dim(res$fixed_effects)[1] == 0 && dim(res$fixed_effects)[2] == 0){
    res$fixed_effects <- array(0, dim = c(res$n, 0))
  }
  res$n_weights <- sum(table(prior_obj$node_data$edges$from)-1) # the number of weights (triple split = 2 weights)
  # store information from prior_obj
  res$n_hd <- sum(prior_obj$node_data$orig_nodedata$status == "attached") # how many variance parameters involved in hd prior
  res$n_pc <- if (length(prior_obj$prior_data$weights) > 0) sum(sapply(prior_obj$prior_data$weights, function(x) x$prior == "pc")) else 0

  # which components are included in a tree and have HD priors
  res$which_hd <- as.array(c(1:(res$n_r+1))[prior_obj$node_data$orig_nodedata$status == "attached"])
  res$which_cw <- as.array(c(1:(res$n_r+1))[prior_obj$node_data$orig_nodedata$status == "detached"])
  res$n_cw <- length(res$which_cw)

  # # how many nodes are in each tree
  # # for each node with top_status = TRUE, we have a tree
  # # TODO must have some rule for how the different trees are identified in stan
  # res$n_in_tree <- 0
  # res$n_trees <- length(res$n_in_tree)

  res$which_pc <- if (res$n_hd > 0) as.array(sapply(prior_obj$prior_data$weights, function(x) if (x$prior == "pc") 1 else 0)) else integer()

  res$cw_prior_numbers <- get_variance_prior_number(sapply(prior_obj$prior_data$cw_priors, function(x) x$prior))
  res$cw_prior_parameters <- t(sapply(prior_obj$prior_data$cw_priors,
                                      function(x) {
                                        param <- x$param
                                        return(if (length(param) == 2) param else rep(param, 2))
                                      }))

  res$cw_prior_info <- cbind(res$cw_prior_numbers, res$cw_prior_parameters)

  res$totvar_prior_numbers <- if (res$n_hd > 0) get_variance_prior_number(sapply(prior_obj$prior_data$total_variance, function(x) x$prior)) else array(0, dim = c(1, 1))
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
  top_nodes <- get_top_nodes(prior_obj$node_data)
  split_nodes <- get_split_ids(prior_obj$node_data)
  which_theta_in_hd_mat <- matrix(0, nrow = res$n_totvar, ncol = res$n_r+1)
  res$n_splits_each_tree <- array(0, dim = 1)
  for (ind in seq_len(res$n_totvar)){
    # does not include top nodes, because they are not in the original nodes (as they are merged nodes)
    in_tree <- app_to_stan_id$stan_id[app_to_stan_id$app_id %in% nodes_in_tree(prior_obj$node_data, top_nodes[ind])]
    which_theta_in_hd_mat[ind, in_tree] <- 1
    res$n_splits_each_tree[ind] <- sum(split_nodes %in% nodes_in_tree(prior_obj$node_data, top_nodes[ind]))
  }

  res$n_splits <- length(prior_obj$weight_priors)

  stopifnot(sum(res$n_splits_each_tree) == res$n_splits)

  res$which_theta_in_hd <- which_theta_in_hd_mat

  # make matrices with information on how each weight looks,
  # each row corresponds to a weight, each column to a variance
  weight_info_over <- matrix(0, nrow = length(prior_obj$node_data$edges$from), ncol = res$n_r+1)
  weight_info_under <- matrix(0, nrow = length(prior_obj$node_data$edges$from), ncol = res$n_r+1)
  res$n_rows_in_w_o_u <- nrow(weight_info_over)

  # each column corresponds to a row in weight_info_over/under, rows to trees
  # (the first is for index 0, so there is one more column than rows in w_i_o/u)
  row_index_hd_pr_mat <- array(0, dim = c(res$n_totvar, length(prior_obj$node_data$edges$from)+1))
  split_nodes <- get_split_ids(prior_obj$node_data)
  ind3 <- 1 # number of weights (total, w and 1-w counts as two)
  hd_pr_num <- numeric()
  for (ind in seq_len(sum(res$n_splits))){
    tmp <- get_leaf_node_ids_for_split(prior_obj$node_data, split_nodes[ind])
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
    ##tree_no <- which(which_theta_in_hd_mat[,under[1]] == 1)
    row_index_hd_pr_mat[tree_no, ind3] <- 1
  }

  row_index_hd_pr <- array(0, dim = 1)

  # this loop goes over each weight (each HD prior in the prior)
  for (ind2 in seq_len(length(unique(hd_pr_num)))){
    row_index_hd_pr <- c(row_index_hd_pr, which(hd_pr_num == ind2)[length(which(hd_pr_num == ind2))])
  }

  row_index_hd_pr <- row_index_hd_pr

  res$row_index_hd_pr <- row_index_hd_pr_mat # which weight has which rows in the weight_info-matrices
  res$row_index_hd_pr_plot <- row_index_hd_pr + 1 # which weight has which rows in the weight_info-matrices

  res$w_o <- weight_info_over
  res$w_u <- weight_info_under

  res$use_intercept <- as.integer(prior_obj$use_intercept)
  res$use_likelihood <- use_likelihood

  res$effect_sizes <- array(sapply(effects, function(x) length(unique(x))),
                            dim = c(res$n_r))
  res$effect_index_start <- cumsum(res$effect_sizes)-res$effect_sizes
  res$indexes <- as.data.frame(effects)

  res$n_knots <- 122 # in the code we have now this is 122, and cannot be changed by the user
  res$knots <- c(-500,
                 seq(-200, -40, length.out = 10)[1:9],
                 seq(-40, -5, length.out = 22)[1:21],
                 seq(-5, -1e-6, length.out = 30),
                 seq(1e-6, 5, length.out = 30),
                 seq(5, 40, length.out = 22)[-1],
                 seq(40, 200, length.out = 10)[-1],
                 500) # for now these are always the same

  tmp_array <- array(0, dim = c(res$n_pc, 122, 4))
  pc_ind <- if (res$n_hd > 0) which(sapply(prior_obj$weight_priors, function(x) length(x$knots) > 2)) else c()
  # if (length(pc_ind) == 0) pc_ind <- 0
  for (i in seq_len(length(pc_ind))){
    tmp_array[i,,] <- prior_obj$weight_priors[[pc_ind[i]]]$coeffs
  }
  res$prior_coeffs <- tmp_array

  res$y <- y

  # intercept/fixed effects priors
  res$intercept_prior <- array(prior_obj$prior_intercept, dim = c(2))
  res$fixed_priors <- array(as.numeric(prior_obj$prior_fixed), dim = c(res$n_f, 2))
  
  #browser()
  # default: nothing on iid, sum-to-zero on all other effects, linear sum-to-zero on rw2
  # user can specify something else (but if user says constr on rw2, both will have constraints)
  # for rw2: INLA:::inla.rw2 gir precision matrix, lag en sånn funksjon, finn pseudoinverse,
    # legg på noe på diagonalen (kanskje) [brukes til å lage prior], finn cholesky
  
  # which latent model is each random effect
  res$effect_type <- sapply(prior_obj$model_data, model_to_number)
  # res$effect_type <- sapply(pro$model_data, model_to_number)
  
  
  # # info on linear constraints
  # res$stz_constraint1 <- c()
  # res$stz_constraint2 <- c()
  # res$constr_vec <- rep(0, sum(res$effect_sizes)) 
  # for (ind in 1:res$n_r){
  #   res$stz_constraint1[ind] <- as.numeric(prior_obj$model_data[[ind]]$constr1)
  #   res$stz_constraint2[ind] <- as.numeric(prior_obj$model_data[[ind]]$constr2)
  #   if (res$stz_constraint2[ind] == 1){
  #     res$constr_vec[(1+res$effect_index_start[ind]):(res$effect_sizes[ind]+res$effect_index_start[ind])] <- 1:res$effect_sizes[ind] - mean(1:res$effect_sizes[ind])
  #   }
  # }
  # res$stz_constraint1 <- array(res$stz_constraint1)
  # res$stz_constraint2 <- array(res$stz_constraint2)

  which_covmats_structured <- c()
  mega_matrix <- array(0, dim = c(res$n_r, res$n, res$n))
  for (i in seq_along(prior_obj$which_covmats_structured)){
    ind <- prior_obj$which_covmats_structured[i]
    # warning("Dette skal vel ikke stå her?")
    # tmpmat <- prior_obj$covmats[[ind]][1:res$effect_sizes[ind],1:res$effect_sizes[ind]] + 10^(-9)*diag(res$effect_sizes[ind])
    # lower_chol <- t(chol(tmpmat/typical_variance(tmpmat)))
    lower_chol <- t(chol(prior_obj$covmats[[ind]][1:res$effect_sizes[ind],1:res$effect_sizes[ind]] + 10^(-9)*diag(res$effect_sizes[ind])))
    # we need the lower cholesky triangle
    mega_matrix[ind, 1:res$effect_sizes[ind], 1:res$effect_sizes[ind]] <- lower_chol
    which_covmats_structured <- c(which_covmats_structured, ind)
  }
  res$mega_matrix <- mega_matrix
  res$use_index_matrix <- array(0, dim = res$n_r)
  res$use_index_matrix[which_covmats_structured] <- 1

  return(res)

}


# input is model data for one latent component
model_to_number <- function(model_data){
  
  if (model_data$model == "iid"){
    if (!model_data$constr) return(1) else return(2)
  } else if (model_data$model == "rw2"){
    if (!model_data$constr) return(3) else if (!model_data$constr_lin) return(4) else return(5)
  } else if (model_data$model == "z" || model_data$model == "generic0"){ # treated equally in stan
    if (!model_data$constr) return(6) else return(7)
  }
  return(0)
  
}


# making matrix for rw2-effect
# dim = size of effect
# type = precision matrix or covariance matrix ("prec", "cov" or "lower")
make_rw2_mat <- function(n, type = c("prec", "cov", "lower")){
  
  # precision matrix
  prec_mat <- matrix(0, nrow = n, ncol = n)
  
  # diagonal
  diag(prec_mat) <- 1
  diag(prec_mat)[-c(1,2,n-1,n)] <- 6
  diag(prec_mat)[c(2, n-1)] <- 5
  
  # off-diagonal 1
  prec_mat[matrix(c(1:(n-1), 2:n, 2:n, 1:(n-1)), ncol = 2)] <- -4
  prec_mat[matrix(c(1,2,2,1,n-1,n,n,n-1), ncol = 2, byrow = T)] <- -2
  
  # off-diagonal 2
  prec_mat[matrix(c(1:(n-2), 3:n, 3:n, 1:(n-2)), ncol = 2)] <- 1
  
  if (type == "prec") return(prec_mat)
  
  covmat <- MASS::ginv(prec_mat)
  if (type == "cov") return(covmat)
  
  lower <- t(chol(covmat))
  if (type == "lower") return(lower)
  
  return(diag(n))
  
}



# assumes that the first vector in the is the data (NOT ANYMORE???)
# calculate the covariance matrices for each effect, to use to make the PC priors
# all matrices are of the Matrix-class, so they are sparce
make_covariance_matrices <- function(data){

  matrix_data <- lapply(data, make_covariance_matrix)
  matrix_data <- lapply(matrix_data, scale_cov_matrix)
  names(matrix_data) <- names(data)
  # matrix_data <- lapply(data[-1], make_covariance_matrix) # -1 because we do not want a covariance matrix for the data
  # matrix_data <- lapply(matrix_data, scale_cov_matrix)
  # names(matrix_data) <- names(data[-1])

  return(matrix_data)

}


# index is the index vector for each random effect
make_covariance_matrix <- function(index){

  # if this is a vector of data and not indexes, we make an index vector instead
  if (!is.integer(index)) index <- 1:length(index)

  # the size of the dataset
  n <- length(index)
  mat <- matrix(0, nrow = n, ncol = n)

  for (i_r in 1:n){
    mat[i_r, index == index[i_r]] <- 1
  }

  return(Matrix(mat))

}

# generates scaled (so variances are comparable) covariance matrices
# covmats = list, or covmat1 and covmat2 = covariance matrices involved in this split
scale_cov_matrix <- function(covmat) return(covmat/typical_variance(covmat))

# calculates the typical variance
typical_variance <- function(mat) return(exp(mean(log(Matrix::diag(mat)))))

check_covariance_matrices <- function(mats){

  n <- nrow(mats[[1]])

  for (i in 1:length(mats)){
    tmp <- if (base::qr(mats[[i]])$rank < n) "singular" else "not singular"
    cat(names(mats)[i], "is", tmp, "\n")
  }

}

# if the effect has shorter length than the number of observations, it will be
# repeated in a block matrix
fix_size_covmat <- function(new_mat, old_mat, eff_name){

  # old_dim >= new_dim, since old_mat has dimension nxn (n = number of observations)
  new_dim <- dim(new_mat)[1]
  old_dim <- dim(old_mat)[1]

  if (new_dim != old_dim){
    if ((old_dim/new_dim)%%1 != 0){
      warning(paste0(
        "Wrong dimension of input covariance matrix for effect", eff_name,
        ", using one based on index input!"
      ))
      return(old_mat)
    } else {
      return(
        kronecker(diag(old_dim/new_dim), new_mat)
      )
    }
  } else { # if the dimensions are equal, no resizing necessary
    return(new_mat)
  }

}

get_variance_prior_number <- function(prior_names){

  pr_names <- c("jeffreys", "pc0", "invgam", "hc", "")

  pr_num <- c()
  for (i in 1:length(prior_names)){
    pr_num[i] <- which(pr_names == prior_names[i])
  }

  pr_num[pr_num == 5] <- 0

  return(pr_num)

}


# for each split (i.e. for a given split), get id of the leaf-nodes involved in the split,
# and whether it is above or below the fraction in the corresponding proportion
get_leaf_node_ids_for_split <- function(node_data, split_node){

  # if the node is not a split node, return node NULL
  if (is_leaf_node(node_data, list(current_node_id = split_node))) return(NULL)

  # which original nodes (aka components) are involved in the split
  nodes_in_split <- sort(get_leaf_nodes(node_data, split_node))

  # these components are under the fraction in this split for all weights:
  under <- nodes_in_split

  subnodes <- node_data$edges$to[node_data$edges$from == split_node]

  over <- list()
  # this is the same as the number of edges from the node in question
  for (edg in 1:(length(subnodes))){

    over[[edg]] <- get_leaf_nodes(node_data, subnodes[edg])

  }

  return(list(under = under, over = over))

}

# get data-frame with node ids and node names
get_node_name_id <- function(node_data){

  # only the original node labels will be used, since they are the only ones that are variances (new nodes are for "help" in the tree)
  # this does NOT include the top node (but it is not in orig_nodedata anyways)
  node_name_id <- data.frame(
    name = node_data$orig_nodedata$label, # names of original nodes (to use for parameter names)
    id = node_data$orig_nodedata$id # ids of original nodes (for connecting edge and node data)
  )

  return(node_name_id)

}




















## INLA stuff
# TODO rewrite to make more efficient later!!




#' Run inference
#'
#' This function helps you run inference with INLA using a prior object from construct_prior
#' @param prior_obj An object from ``create_prior``, from ``run_shiny``
#' (both with and without the "Close and run" option), from ``inference_stan``, or from ``inference_inla`` (for refitting model)
#' @param use_likelihood Whether to return the prior (0, can be used for e.g. debugging or to look at the priors on variance parameters when using an HD prior), or to use the likelihood and data to get the posterior (1, default).
#' @param print_prior Whether to print a text with the chosen prior or not (default TRUE)
#' @param args A named list with certain arguments for INLA (expert option)
#' @keywords inference, INLA
#' @return A named list with a prior object (prior), an inla-object (inla) and some data inla requires (inla_data).
#' @export
#' @examples
#' vignette("create_prior", package = "priorconstruction")

inference_inla <- function(prior_obj, use_likelihood = 1, 
                           print_prior = TRUE,
                           args = list(initial_res = 0), ...){
  
  if (print_prior) print_prior_choice(prior_obj)
  
  if (length(prior_obj) == 3) prior_obj <- prior_obj$prior
  
  # return only priors (on precision scale)
  if (use_likelihood == 0){
    prior_obj$response <- rep(NA, length(prior_obj$response))
  }
  
  # make a formula for inla
  
  inla_formula <- make_inla_formula(prior_obj)
  
  inla_data <- c(list(y = prior_obj$response), 
                 prior_obj$data$fixed, 
                 prior_obj$data$random,
                 prior_obj$inference_data)
  
  inla_jpr_data <- make_inla_data_object(prior_obj)
  
  inla_res <- inla(
    formula = inla_formula,
    data = inla_data,
    control.family = list(initial = args$initial_res),
    control.fixed = list(mean.intercept = prior_obj$prior_intercept[1], prec.intercept = 1/prior_obj$prior_intercept[2]^2),
    control.expert = list(jp = make_jpr(prior_obj, inla_jpr_data)),
    ...
  )
  
  return(list(
    prior = prior_obj,
    inla = inla_res,
    inla_data = list(
      w_o = inla_jpr_data$w_o,
      w_u = inla_jpr_data$w_u,
      row_index_hd_pr_plot = inla_jpr_data$row_index_hd_pr_plot,
      which_hd = inla_jpr_data$which_hd)
  )
  )
  
}



# create a formula that inla can use
# all variances goes into the jp-part
make_inla_formula <- function(prior_obj){
  
  # the prior object only contains the data from the effects included, so we can make a formula using the
  # existing effects
  
  effnames_fixed <- names(prior_obj$data$fixed)
  effnames_random <- names(prior_obj$data$random)
  
  # do not want the residuals in the name
  if ("eps" %in% effnames_random) effnames_random <- effnames_random[!(effnames_random == "eps")]
  
  inla_formula <- "y ~ "
  # fixed effects
  comps_f <- c()
  for (ind in seq_len(length(effnames_fixed))) {
    ind2 <- which(row.names(prior_obj$prior_fixed) == effnames_fixed[ind])
    comps_f[ind] <- paste0("f(", effnames_fixed[ind], ", model = \"linear\", ",
                           "mean.linear = ", as.numeric(prior_obj$prior_fixed[ind2,1]), ", ",
                           "prec.linear = ", 1/(as.numeric(prior_obj$prior_fixed[ind2,2])^2),
                           ")")
  }
  if (length(comps_f) > 0) inla_formula <- paste0(inla_formula, paste0(comps_f, collapse = " + "))
  
  if (!prior_obj$use_intercept) inla_formula <- paste0(inla_formula, " - 1") # remove intercept if not used
  
  # random effects
  comps_r <- c()
  for (ind in seq_len(length(effnames_random))){
    
    # if there is something else to add to the formula
    extra_arguments <- ""
    
    if (prior_obj$model_data[[ind]]$label == effnames_random[ind]) {
      mod_type <- prior_obj$model_data[[ind]]$model
      # TODO: if this is not correct, results may be weird
      stz_constr <- prior_obj$model_data[[ind]]$constr
      if (mod_type == "rw2" && prior_obj$model_data[[ind]]$constr_lin == T){
        stz_constr <- paste0(
          # stz_constr,
          "FALSE",
          ", rankdef = 2",
          ", extraconstr = list(A = matrix(rep(1:", length(unique(prior_obj$data$random[[ind]])), ", each = 2), nrow = 2), e = matrix(0, 2, 1))"
          # ", extraconstr = list(A = matrix(1:", length(unique(prior_obj$data$random[[ind]])), ", nrow = 1), e = matrix(0, 1, 1))"
        )
      }
      if (mod_type == "generic0"){
        extra_arguments <- paste0(", Cmatrix = ", prior_obj$model_data[[ind]]$inference_data$Cmatrix)
      }
    } else {
      warning(paste0(effnames_random[ind], " did not get correct model type!"))
      mod_type <- "iid"
    }
    scalemodel <- if (mod_type == "rw2") ", scale.model = TRUE" else ""
    
    comps_r[ind] <- paste0("f(",
                           effnames_random[ind],
                           ", ",
                           "model = \"", 
                           mod_type,
                           "\", ", 
                           "constr = ",
                           stz_constr, 
                           ", initial = 0",
                           extra_arguments, # additional arguments for inla formula
                           scalemodel,
                           ")")
    
  }
  
  if (length(comps_r) > 0) {
    if (prior_obj$use_intercept && length(comps_f) > 0){
      inla_formula <- paste0(inla_formula, " + ", paste0(comps_r, collapse = " + "))
    } else { # nothing added to the formula yet
      inla_formula <- paste0(inla_formula, paste0(comps_r, collapse = " + "))
    }
  }
  
  inla_formula <- formula(inla_formula)
  
  return(inla_formula)
  
}



# fixed effects can be a data.frame or a list
make_inla_data_object <- function(prior_obj, use_likelihood = 1){
  
  # in case y has some weird format
  y <- as.numeric(unlist(prior_obj$response))
  
  effects <- prior_obj$data$random
  res <- list()
  
  tmp <- get_node_name_id(prior_obj$node_data)
  
  # here we decide which variance parameter goes where in the theta-vector to stan
  # may not be necessary when the original nodes are initialized from 1, but does not hurt I guess
  app_to_stan_id <- data.frame(
    stan_id = 1:length(effects),
    app_id = tmp$id,
    app_name = tmp$name
  )
  
  # if gaussian data and we have "eps" in the data, it should be removed after storing info about the residual node
  if (length(effects$eps) == 0){ # for now, the likelihood is always gaussian
    # if (prior_obj$likelihood == "gaussian" && !is.null(effects$eps)){
    effects <- effects[which(names(effects) != "eps")]
  }
  
  res$n <- length(y)
  res$n_r <- length(effects) # this is the number of non-residual random effects
  res$n_f <- length(prior_obj$data$fixed) #ncol(res$fixed_effects)
  res$fixed_effects <- as.data.frame(prior_obj$data$fixed)
  if (dim(res$fixed_effects)[1] == 0 && dim(res$fixed_effects)[2] == 0){
    res$fixed_effects <- array(0, dim = c(res$n, 0))
  }
  res$n_weights <- sum(table(prior_obj$node_data$edges$from)-1) # the number of weights (triple split = 2 weights)
  # store information from prior_obj
  res$n_hd <- sum(prior_obj$node_data$orig_nodedata$status == "attached") # how many variance parameters involved in hd prior
  res$n_pc <- if (length(prior_obj$prior_data$weights) > 0) sum(sapply(prior_obj$prior_data$weights, function(x) x$prior == "pc")) else 0
  
  # which components are included in a tree and have HD priors
  res$which_hd <- as.array(c(1:(res$n_r+1))[prior_obj$node_data$orig_nodedata$status == "attached"])
  res$which_cw <- as.array(c(1:(res$n_r+1))[prior_obj$node_data$orig_nodedata$status == "detached"])
  res$n_cw <- length(res$which_cw)
  
  # # how many nodes are in each tree
  # # for each node with top_status = TRUE, we have a tree
  # # TODO must have some rule for how the different trees are identified in stan
  # res$n_in_tree <- 0
  # res$n_trees <- length(res$n_in_tree)
  
  res$which_pc <- if (res$n_hd > 0) as.array(sapply(prior_obj$prior_data$weights, function(x) if (x$prior == "pc") 1 else 0)) else integer()
  
  res$cw_prior_numbers <- get_variance_prior_number(sapply(prior_obj$prior_data$cw_priors, function(x) x$prior))
  res$cw_prior_parameters <- t(sapply(prior_obj$prior_data$cw_priors,
                                      function(x) {
                                        param <- x$param
                                        return(if (length(param) == 2) param else rep(param, 2))
                                      }))
  
  res$cw_prior_info <- cbind(res$cw_prior_numbers, res$cw_prior_parameters)
  
  res$totvar_prior_numbers <- if (res$n_hd > 0) get_variance_prior_number(sapply(prior_obj$prior_data$total_variance, function(x) x$prior)) else array(0, dim = c(1, 1))
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
  top_nodes <- get_top_nodes(prior_obj$node_data)
  split_nodes <- get_split_ids(prior_obj$node_data)
  which_theta_in_hd_mat <- matrix(0, nrow = res$n_totvar, ncol = res$n_r+1)
  res$n_splits_each_tree <- array(0, dim = 1)
  for (ind in seq_len(res$n_totvar)){
    # does not include top nodes, because they are not in the original nodes (as they are merged nodes)
    in_tree <- app_to_stan_id$stan_id[app_to_stan_id$app_id %in% nodes_in_tree(prior_obj$node_data, top_nodes[ind])]
    which_theta_in_hd_mat[ind, in_tree] <- 1
    res$n_splits_each_tree[ind] <- sum(split_nodes %in% nodes_in_tree(prior_obj$node_data, top_nodes[ind]))
  }
  
  res$n_splits <- length(prior_obj$weight_priors)
  
  stopifnot(sum(res$n_splits_each_tree) == res$n_splits)
  
  res$which_theta_in_hd <- which_theta_in_hd_mat
  
  # make matrices with information on how each weight looks,
  # each row corresponds to a weight, each column to a variance
  weight_info_over <- matrix(0, nrow = length(prior_obj$node_data$edges$from), ncol = res$n_r+1)
  weight_info_under <- matrix(0, nrow = length(prior_obj$node_data$edges$from), ncol = res$n_r+1)
  res$n_rows_in_w_o_u <- nrow(weight_info_over)
  
  # each column corresponds to a row in weight_info_over/under, rows to trees
  # (the first is for index 0, so there is one more column than rows in w_i_o/u)
  row_index_hd_pr_mat <- array(0, dim = c(res$n_totvar, length(prior_obj$node_data$edges$from)+1))
  split_nodes <- get_split_ids(prior_obj$node_data)
  ind3 <- 1 # number of weights (total, w and 1-w counts as two)
  hd_pr_num <- numeric()
  for (ind in seq_len(sum(res$n_splits))){
    tmp <- get_leaf_node_ids_for_split(prior_obj$node_data, split_nodes[ind])
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
    ##tree_no <- which(which_theta_in_hd_mat[,under[1]] == 1)
    row_index_hd_pr_mat[tree_no, ind3] <- 1
  }
  
  row_index_hd_pr <- array(0, dim = 1)
  
  # this loop goes over each weight (each HD prior in the prior)
  for (ind2 in seq_len(length(unique(hd_pr_num)))){
    row_index_hd_pr <- c(row_index_hd_pr, which(hd_pr_num == ind2)[length(which(hd_pr_num == ind2))])
  }
  
  row_index_hd_pr <- row_index_hd_pr
  
  res$row_index_hd_pr <- row_index_hd_pr_mat # which weight has which rows in the weight_info-matrices
  res$row_index_hd_pr_plot <- row_index_hd_pr + 1 # which weight has which rows in the weight_info-matrices
  
  res$w_o <- weight_info_over
  res$w_u <- weight_info_under
  
  res$use_intercept <- as.integer(prior_obj$use_intercept)
  res$use_likelihood <- use_likelihood
  
  res$effect_sizes <- array(sapply(effects, function(x) length(unique(x))),
                            dim = c(res$n_r))
  res$effect_index_start <- cumsum(res$effect_sizes)-res$effect_sizes
  res$indexes <- as.data.frame(effects)
  
  res$n_knots <- 122 # in the code we have now this is 122, and cannot be changed by the user
  res$knots <- c(-500,
                 seq(-200, -40, length.out = 10)[1:9],
                 seq(-40, -5, length.out = 22)[1:21],
                 seq(-5, -1e-6, length.out = 30),
                 seq(1e-6, 5, length.out = 30),
                 seq(5, 40, length.out = 22)[-1],
                 seq(40, 200, length.out = 10)[-1],
                 500) # for now these are always the same
  
  tmp_array <- array(0, dim = c(res$n_pc, 122, 4))
  pc_ind <- if (res$n_hd > 0) which(sapply(prior_obj$weight_priors, function(x) length(x$knots) > 2)) else c()
  # if (length(pc_ind) == 0) pc_ind <- 0
  for (i in seq_len(length(pc_ind))){
    tmp_array[i,,] <- prior_obj$weight_priors[[pc_ind[i]]]$coeffs
  }
  res$prior_coeffs <- tmp_array
  
  res$y <- y
  
  # intercept/fixed effects priors
  res$intercept_prior <- array(prior_obj$prior_intercept, dim = c(2))
  res$fixed_priors <- array(as.numeric(prior_obj$prior_fixed), dim = c(res$n_f, 2))
  
  #browser()
  # default: nothing on iid, sum-to-zero on all other effects, linear sum-to-zero on rw2
  # user can specify something else (but if user says constr on rw2, both will have constraints)
  # for rw2: INLA:::inla.rw2 gir precision matrix, lag en sånn funksjon, finn pseudoinverse,
  # legg på noe på diagonalen (kanskje) [brukes til å lage prior], finn cholesky
  
  # which latent model is each random effect
  res$effect_type <- sapply(prior_obj$model_data, model_to_number)
  # res$effect_type <- sapply(pro$model_data, model_to_number)
  
  
  # # info on linear constraints
  # res$stz_constraint1 <- c()
  # res$stz_constraint2 <- c()
  # res$constr_vec <- rep(0, sum(res$effect_sizes)) 
  # for (ind in 1:res$n_r){
  #   res$stz_constraint1[ind] <- as.numeric(prior_obj$model_data[[ind]]$constr1)
  #   res$stz_constraint2[ind] <- as.numeric(prior_obj$model_data[[ind]]$constr2)
  #   if (res$stz_constraint2[ind] == 1){
  #     res$constr_vec[(1+res$effect_index_start[ind]):(res$effect_sizes[ind]+res$effect_index_start[ind])] <- 1:res$effect_sizes[ind] - mean(1:res$effect_sizes[ind])
  #   }
  # }
  # res$stz_constraint1 <- array(res$stz_constraint1)
  # res$stz_constraint2 <- array(res$stz_constraint2)
  
  which_covmats_structured <- c()
  mega_matrix <- array(0, dim = c(res$n_r, res$n, res$n))
  for (i in seq_along(prior_obj$which_covmats_structured)){
    ind <- prior_obj$which_covmats_structured[i]
    # warning("Dette skal vel ikke stå her?")
    # tmpmat <- prior_obj$covmats[[ind]][1:res$effect_sizes[ind],1:res$effect_sizes[ind]] + 10^(-9)*diag(res$effect_sizes[ind])
    # lower_chol <- t(chol(tmpmat/typical_variance(tmpmat)))
    lower_chol <- t(chol(prior_obj$covmats[[ind]][1:res$effect_sizes[ind],1:res$effect_sizes[ind]] + 10^(-9)*diag(res$effect_sizes[ind])))
    # we need the lower cholesky triangle
    mega_matrix[ind, 1:res$effect_sizes[ind], 1:res$effect_sizes[ind]] <- lower_chol
    which_covmats_structured <- c(which_covmats_structured, ind)
  }
  res$mega_matrix <- mega_matrix
  res$use_index_matrix <- array(0, dim = res$n_r)
  res$use_index_matrix[which_covmats_structured] <- 1
  
  return(res)
  
}




# dirichlet for a split with variances theta
# do not include the jacobian from V/w to theta!
# theta = vector with all log variances in this tree (may have several trees in the model, only the ones in this tree is in "theta")
# m = how many variance parameters are involved in this split
# weights = vector (of at least length 2) with which weights are involved in this particular split
# w_o = matrix with variances over fraction in each weight
# w_u = matrix with variances under fraction in each weight
hd_dirichlet_prior_lpdf <- function(theta, m, weights_start, w_o, w_u, get_indexes, get_dirichlet_parameter){
  
  # if only two weights (that will sum to 1) we have a uniform prior and the log density is 0
  if (m == 2){
    return(0);
  }
  
  weights <- c()
  for (i in 1:m){
    weights[i] = weights_start+i-1;
  }
  
  n_tot = sum(w_u[weights[1],]); # how many variance parameters are involved in this split
  
  over_part = 0;
  
  alpha = get_dirichlet_parameter(m);
  
  for (ind in 1:m){
    over_part <- over_part + log(sum(exp(theta[get_indexes(w_o[weights[ind],])])));
  }
  
  norm_const <- 0
  
  # the part under the fraction is the same for all weights in the multisplit,
  # so we can use the same row in w_u for all
  under_part <- log(sum(exp(theta[get_indexes(w_u[weights[1],])])));
  
  return(
    norm_const + (alpha-1) * ( over_part - m*under_part )
  );
  
  
}

expit <- function(x) 1/(1+exp(-x))


# evaluates hd prior for a given split (which is identified by the parameter weights)
# includes jacobian from logit weight to weight!
# not including the jacobian from V/w to theta!
# theta = vector with all log variances in the tree (may have several trees in the model, only the ones in this tree is in "theta")
# weights_start = index of where in w_o and w_u this weight is located
# w_o = matrix with variances over fraction in each weight
# w_u = matrix with variances under fraction in each weight
# the rest is for the spline
hd_pc_prior_lpdf <- function(theta, weights_start, idx, w_o, knots, prior_coeffs, n_knots, get_indexes, eval_spline_lpdf, expit){
  
  logdens <- 0
  
  # expit <- function(x) 1/(1+exp(-x))
  
  # calculate logit weight from log variances
  logitw = log(sum(exp(theta[get_indexes(w_o[weights_start,])]))) - log(sum(exp(theta[get_indexes(w_o[weights_start+1,])])))
  
  logdens <- logdens + eval_spline_lpdf(logitw, knots, prior_coeffs, idx, n_knots)
  
  # jacobian from logit weight to weight
  logdens <- logdens -log( expit(logitw)*(1-expit(logitw)) )
  
  return(logdens)
  
}


# evaluates all priors that has HD priors in a tree (run several times for multiple trees)
# includes total variance
# includes jacobians!
# theta are all log variances involved in this prior tree (may have more than one tree)
hd_prior_joint_lpdf <- function(theta, use_likelihood, which_pc, w_o, w_u, n_splits_tree, row_index_hd_pr, knots, prior_coeffs, n_knots, totvar_prior_info, indV,
                                choose_prior_lpdf, hd_dirichlet_prior_lpdf, hd_pc_prior_lpdf,
                                calc_jac_logdet, get_indexes, get_indexes2, get_dirichlet_parameter,
                                eval_spline_lpdf, expit){
  
  # if we do not have a HD prior, theta has length 0, and we return 0
  if (length(theta) == 0){
    return(0)
  }
  
  logV <- log(sum(exp(theta)))
  logdens <- 0
  
  # prior on the total variance
  if (use_likelihood == 0 && totvar_prior_info[indV,1] == 1){ # must have proper variance when we sample from prior
    logdens <- logdens + dnorm(logV, 0, 1, log = TRUE)
  } else {
    logdens <- logdens + choose_prior_lpdf(logV, totvar_prior_info[indV,1], totvar_prior_info[indV,2:3])
  }
  
  pc_prior_ind <- 1
  # prior on the weights
  for (indw in 1:n_splits_tree){
    if (which_pc[indw] == 1){ # pc prior on weight
      # jacobian from logit(w) to w is included in the function
      logdens <- logdens + hd_pc_prior_lpdf(theta, row_index_hd_pr[indw], pc_prior_ind, w_o, knots, prior_coeffs, n_knots, get_indexes, eval_spline_lpdf, expit)
      pc_prior_ind <- pc_prior_ind + 1
      # logdens += hd_pc_prior_lpdf(theta | row_index_hd_pr[indw+1]-(row_index_hd_pr[indw]), indw, w_o, knots, prior_coeffs, n_knots);
    } else { # dirichlet prior on weight(s)
      # second argument used to be: row_index_hd_pr[indw+1]-row_index_hd_pr[indw]-1
      logdens <- logdens + hd_dirichlet_prior_lpdf(theta, row_index_hd_pr[indw+1]-(row_index_hd_pr[indw]), row_index_hd_pr[indw], w_o, w_u, get_indexes, get_dirichlet_parameter)
    }
  }
  
  logdens <- logdens + calc_jac_logdet(theta, w_o, w_u, row_index_hd_pr, n_splits_tree, get_indexes2) # jacobian from total variance weights to log variance
  logdens <- logdens -logV# jacobian from log total variance to total variance
  
  return(logdens);
  
}



get_dirichlet_parameter <- function(num_comp){
  pars <- c(1, 0.7535, 0.6834, 0.6485, 0.6274,
            0.6132, 0.6029, 0.5951, 0.589, 0.5842,
            0.5801, 0.5768, 0.5739, 0.5715, 0.5693,
            0.5675, 0.5658, 0.5643, 0.563)
  return(pars[num_comp-1])
}

cw_priors_lpdf <- function(theta, cw_prior_info, choose_prior_lpdf){
  
  logdens <- 0
  
  for (indcw in 1:length(theta)){
    if ((cw_prior_info[indcw,1]) > 0){
      logdens <- logdens + choose_prior_lpdf(theta[indcw], cw_prior_info[indcw, 1], cw_prior_info[indcw, 2:3]);
    }
  }
  
  return(logdens)
  
}


# to get correct prior
# input is log variance, which prior is used, and a vector of parameters
choose_prior_lpdf <- function(x, prior_number, param){
  
  if (prior_number == 1){
    return(0)
  } else if (prior_number == 2){
    shape <- -log(param[2])/param[1]; # param = c(U, alpha)
    return(
      log(shape/2) + x/2 - shape*exp(x/2)
    )
  } else if (prior_number == 3){
    return(
      -param[1]*x - 1/(param[2]*exp(x))
    )
  } else if (prior_number == 4){
    return(
      x/2 - log(pi) - log(param[1] + exp(x)/param[1])
    )
  } else if (prior_number == 5) {
    return(dnorm(x, 0, 1, log = TRUE));
  } else {
    return(0)
  }
  
}



joint_prior <- function(theta_prec){
  
  # put residual variance at the end, and transform from log precision to log variance
  theta <- -theta_prec[c(2:length(theta_prec), 1)]
  
  logdens <- 0
  
  # for each tree in the prior:
  if (length(args_123$which_hd) > 0){
    for (indV in seq_len(args_123$n_totvar)){
      logdens <- logdens + hd_prior_joint_lpdf(
        theta[get_indexes(args_123$which_theta_in_hd[indV,])],
        args_123$use_likelihood,
        args_123$which_pc,
        args_123$w_o[,get_indexes(args_123$which_theta_in_hd[indV,])],
        args_123$w_u[,get_indexes(args_123$which_theta_in_hd[indV,])],
        args_123$n_splits_each_tree[indV],
        get_indexes(args_123$row_index_hd_pr[indV,]),
        args_123$knots,
        args_123$prior_coeffs,
        args_123$n_knots,
        args_123$totvar_prior_info,
        indV,
        choose_prior_lpdf,
        hd_dirichlet_prior_lpdf,
        hd_pc_prior_lpdf,
        calc_jac_logdet,
        get_indexes,
        get_indexes2,
        get_dirichlet_parameter,
        eval_spline_lpdf,
        expit
      )
    }
  }
  
  # adding individual priors (CW priors)
  if (length(args_123$which_cw) > 0){
    logdens <- logdens + cw_priors_lpdf(
      theta[args_123$which_cw],
      matrix(args_123$cw_prior_info[args_123$which_cw,], ncol = 3),
      choose_prior_lpdf
    )
  }
  
  return(logdens)
  
}


# evaluating prior for the inla.jp-things, can use splines, but keep this for now
eval_spline_lpdf <- function(logitw, knots, prior_coeffs, idx, n_knots){
  
  # Search for correct cubic polynomial
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
  val <- val + prior_coeffs[idx, idd, 2]*dx;
  val <- val + prior_coeffs[idx, idd, 3]*dx*dx;
  val <- val + prior_coeffs[idx, idd, 4]*dx*dx*dx;
  
  return(val)
  
}

# evaluating prior for the inla.jp-things, can use splines, but keep this for now
eval_spline_lpdf2 <- function(logitw, knots, prior_coeffs, idx, n_knots){
  
  # Search for correct cubic polynomial
  ind <- which(knots >= logitw)[1]
  
  dx <- logitw - knots[ind]
  
  val <- prior_coeffs[idx, ind, 1] + 
    prior_coeffs[idx, ind, 2]*dx +
    prior_coeffs[idx, ind, 3]*dx*dx +
    prior_coeffs[idx, ind, 4]*dx*dx*dx
  
  return(val)
  
}


# making the inla.jp-function to be sent to inla
make_jpr <- function(pr, inla_data){
  
  warning("What is the argument 'pr' doing???????????? Nothing????????")
  
  return(
    inla.jp.define(joint_prior,
                   args_123 = inla_data,
                   hd_prior_joint_lpdf = hd_prior_joint_lpdf,
                   calc_jac_logdet = calc_jac_logdet,
                   choose_prior_lpdf = choose_prior_lpdf,
                   cw_priors_lpdf = cw_priors_lpdf,
                   expit = expit,
                   get_dirichlet_parameter = get_dirichlet_parameter,
                   get_indexes = get_indexes,
                   get_indexes2 = get_indexes2,
                   hd_dirichlet_prior_lpdf = hd_dirichlet_prior_lpdf,
                   hd_pc_prior_lpdf = hd_pc_prior_lpdf,
                   eval_spline_lpdf = eval_spline_lpdf
    )
  )
  
  
}


# jacobian for one tree
### TODO: This function can probably have less inputs and things, but is for now just a rewritten version of the function in stan
# theta = all log variances in this tree (all thetas from the hd_prior_joint_lpdf function)
# w_o = matrix with variances over fraction in each weight
# w_u = matrix with variances under fraction in each weight
# n = number of variance components involved in the tree
# row_index_hd_pr = which rows in w_o and w_u that goes to each weight
# n_splits = number of splits in this tree
calc_jac_logdet <- function(theta, w_o, w_u, row_index_hd_pr, n_splits, get_indexes2){
  
  n <- length(theta)
  jac <- matrix(0, n, n)
  
  # log total variance
  jac[1,] <- exp(theta)
  
  i_row <- 2
  
  # for each split in this tree
  for (i_split in 1:n_splits){
    
    under <- ( sum(exp(theta[which(w_u[row_index_hd_pr[i_split],] == 1)])) )^2
    over <- 0
    
    for (ind in row_index_hd_pr[i_split]:(row_index_hd_pr[i_split+1]-2)){
      
      for (i_col in which(w_u[ind,] == 1)){
        
        # testing if the variance corresponding to the index i_col is involved above the fraction bar of this weight
        
        bool_arg1 = 0
        for (i in 1:sum(w_o[ind,])){
          if (i_col == which(w_o[ind,] == 1)[i]){
            bool_arg1 = 1
          }
        }
        
        # if this variance parameter is above the fraction bar and we are differentiating wrt to it
        if (bool_arg1 == 1){
          
          over <- ( sum(exp(theta[get_indexes2(w_o[ind,], w_u[ind,])])) ) * exp(theta[i_col])
          
        } else { # if this variance parameter is not above the fraction bar and we are differentiating wrt to it
          
          over <- -( sum(exp(theta[which(w_o[ind,] == 1)])) ) * exp(theta[i_col])
          
        }
        
        jac[i_row, i_col] <- over/under
        
      }
      
      i_row <- i_row + 1
      
    }
    
  }
  
  return(log(abs(det(jac))))
  
}





# returns the indexes where w_row is 1
# used to pick correct thetas for the logit weights etc
get_indexes <- function(w_row){
  
  res <- c()
  ind = 1
  for (i in 1:length(w_row)){
    if (w_row[i] == 1){
      res[ind] = i
      ind <- ind + 1
    }
  }
  
  return(res)
  
}

get_indexes2 <- function(w_row_o, w_row_u){
  
  res <- c()
  
  ind = 1
  for (i in 1:length(w_row_o)){
    if (w_row_o[i] == 0 && w_row_u[i] == 1){
      res[ind] = i
      ind <- ind + 1
    }
  }
  
  return(res)
  
}







# create a formula that inla can use
# all variances goes into the jp-part
# for debugging only, makes a formula that can be sent to "original" inla
make_inla_formula2 <- function(prior_obj, use_likelihood = 1, scalemod = T){
  
  # the prior object only contains the data from the effects included, so we can make a formula using the
  # existing effects
  
  # TODO: add possibility of having other models (for now everything is "iid")
  
  effnames_fixed <- names(prior_obj$data$fixed)
  effnames_random <- names(prior_obj$data$random)
  
  # do not want the residuals in the name
  if ("eps" %in% effnames_random) effnames_random <- effnames_random[!(effnames_random == "eps")]
  
  inla_formula <- if (use_likelihood == 1) "y ~ " else "y2 ~ "
  # fixed effects
  comps_f <- c()
  for (ind in seq_len(length(effnames_fixed))) {
    ind2 <- which(row.names(prior_obj$prior_fixed) == effnames_fixed[ind])
    comps_f[ind] <- paste0("f(", effnames_fixed[ind], ", model = \"linear\", ",
                           "mean.linear = ", as.numeric(prior_obj$prior_fixed[ind2,1]), ", ",
                           "prec.linear = ", 1/(as.numeric(prior_obj$prior_fixed[ind2,2])^2),
                           ")")
  }
  if (length(comps_f) > 0) inla_formula <- paste0(inla_formula, paste0(comps_f, collapse = " + "))
  
  if (!prior_obj$use_intercept) inla_formula <- paste0(inla_formula, " - 1") # remove intercept if not used
  
  # random effects
  comps_r <- c()
  for (ind in seq_len(length(effnames_random))){
    #ind2 <- which(sapply(prior_obj$prior_data$cw_priors, function(x) x$name) == effnames_random[ind])
    #pr_name <- convert_to_inla_prior_name(prior_obj$prior_data$cw_priors[[ind2]]$prior)
    #pr_param <- prior_obj$prior_data$cw_priors[[ind2]]$param
    if (prior_obj$model_data[[ind]]$label == effnames_random[ind]) {
      mod_type <- prior_obj$model_data[[ind]]$model
      # for rw2, constr1 and constr2 are equal, for the other models constr2 = 0 always
      # TODO: if this is not correct, results may be weird
      stz_constr <- prior_obj$model_data[[ind]]$constr1
      if (prior_obj$model_data[[ind]]$model == "rw2" && prior_obj$model_data[[ind]]$constr2 == T){
        stz_constr <- paste0(
          stz_constr,
          # ", extraconstr = list(A = matrix(c(rep(1, 9), 1:9), nrow = 2, byrow = T), e = matrix(0, nrow = 2, 1))"
          ", extraconstr = list(A = matrix(1:", length(unique(prior_obj$data$random[[ind]])), ", nrow = 1), e = matrix(0, 1, 1))"
        )
      }
    } else {
      warning(paste0(effnames_random[ind], " did not get correct model type!"))
      mod_type <- "iid"
    }
    scalemodel <- if (scalemod) ", scale.model = TRUE" else ""
    
    comps_r[ind] <- paste0("f(",
                           effnames_random[ind],
                           ", ",
                           "model = \"", 
                           mod_type,
                           "\", ", 
                           "constr = ",
                           stz_constr, 
                           ", initial = 0",
                           scalemodel,
                           ", hyper = pcprior",
                           ")")
    
  }
  
  if (length(comps_r) > 0) {
    if (prior_obj$use_intercept && length(comps_f) > 0){
      inla_formula <- paste0(inla_formula, " + ", paste0(comps_r, collapse = " + "))
    } else { # nothing added to the formula yet
      inla_formula <- paste0(inla_formula, paste0(comps_r, collapse = " + "))
    }
  }
  
  inla_formula <- formula(inla_formula)
  
  return(inla_formula)
  
}




#' Print prior choices in text
#'
#' Prints the chosen prior distribution, with tree structure as a string and prior distributions for the parameters
#' @param prior_obj An object from ``create_prior``, from ``run_shiny``
#' (both with and without the "Close and run" option), from ``inference_stan``, or from ``inference_inla``
#' @keywords print, prior
#' @return None
#' @export
#' @examples
#' vignette("create_prior", package = "priorconstruction")
#' 
print_prior_choice <- function(prior_obj){
  
  if (length(prior_obj) == 3) prior_obj <- prior_obj$prior
  
  cat(prior_obj$tree)
  cat("\n\n")
  
  if (length(prior_obj$prior_data$weights) > 0) cat("Weight priors:\n")
  for (ind in seq_along(prior_obj$prior_data$weights)){
    cat("\t", get_prior_expr_text(prior_obj$prior_data$weights[[ind]], prior_obj$node_data, "weight"), sep = "")
  }
  
  warning("Ingeborg: sjekk om totalvariansen faktisk finnes eller bare er tom")
  if (length(prior_obj$prior_data$total_variance) > 0) cat("Total variance priors:\n")
  for (ind in seq_along(prior_obj$prior_data$total_variance)){
    cat("\t", get_prior_expr_text(prior_obj$prior_data$total_variance[[ind]], prior_obj$node_data, "totvar"), sep = "")
  }
  
  if (any(sapply(prior_obj$prior_data$cw_priors, function(x) nchar(x$prior)))) cat("Independent variance priors:\n")
  
  for (ind in seq_along(prior_obj$prior_data$cw_priors)){
    if (prior_obj$prior_data$cw_priors[[ind]]$prior != "") 
      cat("\t", get_prior_expr_text(prior_obj$prior_data$cw_priors[[ind]], prior_obj$node_data, "cw"), sep = "")
  }
  
  cat("\n")
  
}

# return text for the expression for the prior distribution on the parameter provided (similar to the latex-version)
get_prior_expr_text <- function(prior_data, node_data, param = c("cw", "totvar", "weight")){
  
  param <- match.arg(param)
  
  if (param == "cw" || param == "totvar"){
    
    if (param == "cw") param_type <- if (prior_data$prior %in% c("pc0", "hc")) "sigma" else "sigma^2"
    if (param == "totvar") param_type <- if (prior_data$prior %in% c("pc0", "hc")) "sqrt(V)" else "V"
    param_name <- paste0(param_type, "[", prior_data$name, "]")
    
    
    if (prior_data$prior == "pc0") {
      prior_expr <- paste0("PC0(", prior_data$param[1], ", ", prior_data$param[2], ")")
    } else if (prior_data$prior == "jeffreys"){
      prior_expr <- paste0("Jeffreys'\n")
    } else if (prior_data$prior == "invgam"){
      prior_expr <- paste0("InvGam(", prior_data$param[1], ", ", prior_data$param[2], ")\n")
    } else if (prior_data$prior == "hc"){
      prior_expr <- paste0("HC(", prior_data$param[1], ")\n")
    }
    
  } else { # weight
    
    if (prior_data$prior == "dirichlet"){
      
      param_name <- c()
      for (i in 1:(prior_data$no_children-1)){
        t_under <- as.character(prior_data$name)
        t_over <- as.character(get_node_name(node_data, prior_data$children[i]))
        param_name[i] <- sprintf("w[%s/%s]", t_over, t_under)
      }
      param_name <- paste0(param_name, sep = "", collapse = ", ")
      if (prior_data$no_children > 2) param_name <- paste0("(", param_name, ")", sep = "", collapse = "")
      
      prior_expr <- sprintf("Dirichlet(%s)\n", prior_data$no_children)
    } else if (prior_data$prior == "pc"){
      
      t_under <- as.character(prior_data$name)
      t_over <- as.character(get_node_name(node_data, prior_data$param$basemodel_node))
      param_name <- sprintf("w[%s/%s]", t_over, t_under)
      
      if (prior_data$param$basemodel == 0){
        prior_expr <- sprintf("PC0(%s)\n", prior_data$param$median)
      } else if (prior_data$param$basemodel == 1) {
        prior_expr <- sprintf("PC1(%s)\n", prior_data$param$median)
      } else {
        prior_expr <- sprintf("PCM(%s)\n", prior_data$param$median)
      }
      
    }
    
  }
  
  return(
    paste0(param_name, " ~ ", prior_expr)
  )
  
}




















