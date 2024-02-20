



#' @import ggplot2
#' @import Matrix
#' @importFrom grDevices colorRampPalette gray
#' @importFrom stats as.formula dbeta density dexp dgamma dnorm formula median nlm optimize pgamma qgamma quantile rexp rnorm sd
#' @importFrom utils help
#' @importFrom rlang .data
#' @importFrom shinyjs hidden
#' @importFrom methods is
NULL





possible_models <- function() return(c("iid", "linear", "rw1", "rw2", "generic0", "besag"))

#' Define latent component
#'
#' Function for defining a latent component for the HD prior package. All model components must be specified, and if an HD
#' prior is used, this is specified later. See \link[makemyprior]{make_prior} for more details and examples.
#' @param label Name of the component (short names is an advantage as they are used in the app), no default (MUST be provided)
#' @param model Type of model, default is "iid" (see list of models: \link[makemyprior]{makemyprior_models}, \code{makemyprior_models("latent")}).
#' @param constr Sum-to-zero constraints on component (default TRUE)
#' @param lin_constr Linear sum-to-zero constraint, TRUE/FALSE (only for rw2 and only for Stan)
#' @param Cmatrix Precision for this component when \code{model = "generic0"}.
#' We recommend that the matrixes are scaled to the typical variance of the corresponding covariance matrix is 1.
#' This can be done with \link[makemyprior]{scale_precmat} before sending matrix to \link[makemyprior]{mc}.
#' @param graph Path to graph file for besag effect (see details).
#' @param ... Additional arguments used for inference. For inference with \link[rstan]{rstan} or \link[INLA]{inla},
#' the following is useful (especially for the Besag model):
#' \describe{
#'   \item{\code{scale.model}}{if TRUE, the models are scaled so the geometric mean of the variance (typical variance) is 1}
#' }
#' And some additional arguments that can be used by \link[INLA]{inla}. Useful arguments include:
#' \describe{
#'   \item{\code{rankdef}}{number defining rank deficiency of the model}
#'   \item{\code{extraconstr}}{extra linear constraints (in addition to \code{constr})}
#' }
#' See \link[INLA]{f} for details.
#' @keywords formula
#' @export
#' @details The graph argument is a path to a file describing neighbouring relationship on the following form:
#' First row: number of elements
#' The rest: First number is the index of this element, second is the number of neighbours, the rest is the index numbers
#' of all neighbours for this element.
#' If element 1 and 4 are neighbours, 1 should have 4 in its neighbour list, and 4 should have 1.
#' @return For specifying details on this latent component.
#' @examples
#' \dontrun{
#'
#' vignette("make_prior", package = "makemyprior")
#' }
#'
#' @export
mc <- function(label, model = "iid", constr = NULL, lin_constr = FALSE,
               Cmatrix = NULL, graph = NULL, ...){

  # as.character(as.list(substitute(list(...)))[[2]])
  label <- as.character(substitute(label))

  args <- list(...) # extra arguments that can be sent to inla

  if (!is.null(args$hyper) || !is.null(args$prior)) stop("You cannot specify priors here, use make_prior instead.", call. = FALSE)

  if (label == "eps") stop("You cannot use the name 'eps', this is reserved for residuals. Choose another name.", call. = FALSE)

  # make an object (a list) with the random effect information
  res <- list(label = label, model = model)

  if (!(model %in% possible_models())) stop("Not a valid model type!", call. = FALSE)

  if (model == "linear"){
    res$prior <- list(mean.linear = 0, sd.linear = sqrt(1000))
  } else { # random models
    res$prior <- list(prior = "mw") # if no prior provided, we use an HD prior
  }

  if (is.null(constr)){
    if (model %in% c("iid", "rw1", "rw2", "besag")){
      res$constr <- TRUE
    } else {
      res$constr <- FALSE
    }
  } else res$constr <- constr

  res$lin_constr_stan <- if (is.null(lin_constr)) FALSE else lin_constr

  # checking if the user is missing inputs that are required for each model type
  req <- mc_required_input(model)

  # for inla formula (only storing names as strings)
  # TODO: make more efficient somehow
  if ("Cmatrix" %in% req){
    if (is.null(substitute(Cmatrix))) stop(paste0("Missing argument(s) ", paste0(req, collapse = ", "), " for model ", model), call. = FALSE)
    res$model_data$Cmatrix <- deparse(substitute(Cmatrix))
  }
  if ("graph" %in% req) {
    if (is.null(substitute(graph))) stop(paste0("Missing argument(s) ", paste0(req, collapse = ", "), " for model ", model), call. = FALSE)
    res$model_data$graph <- graph
  }

  # make extraconstr into string:
  if (!is.null(args$extraconstr)) args$extraconstr <- paste0(deparse(substitute(args$extraconstr)), collapse = "")

  # add inla-arguments
  res <- c(res, args)

  return(res)

}



# returns which arguments are required for each model type
mc_required_input <- function(model_type){

  if (model_type == "generic0") return(c("Cmatrix"))
  if (model_type == "besag") return("graph")

  return(NULL)

}




#' Making a prior object
#'
#' Make a prior object with all necessary information about the prior and model.
#' The object can either be sent to \link[makemyprior]{makemyprior_gui}
#' or used directly for inference with Stan (\link[makemyprior]{inference_stan}) or INLA (\link[makemyprior]{inference_inla}).
#' \link[makemyprior]{eval_joint_prior} can be used to evaluate the prior.
#' @param formula A formula object, using the function \link[makemyprior]{mc}.
#' @param data The data used in the model, as a \code{data.frame} or \code{list}. All elements must have the same length.
#' @param family A string indicating the likelihood family. \code{gaussian} with identity link is the default.
#' \code{binomial} with logit link and \code{poisson} with log link are also possible.
#' @param prior
#' \describe{
#'   \item{\code{tree}}{
#'     The tree structure as a string. A split is specified as \code{s1 = (a,b)}, where \code{(s1)}
#'     represents a split node and can be any name except names of the input data in \code{data} and
#'     the reserved \code{(eps)}, which is used for residuals for a Gaussian likelihood. Short names
#'     are recommended. Note that these split names are just used in the initial specification, and
#'     will not be used later as they are changed by the function automatically.
#'     The child nodes for each split are included in parentheses separated by commas, and each
#'     split is separated by semicolons.
#'     Singletons are included as \code{(a)}.
#'   }
#'   \item{\code{V}}{A named list with information on the priors on each top node and singleton,
#'   i.e., all variances.
#'   The names in the list are the top node and singleton names from the \code{tree} argument.
#'   Syntax is \code{V = list(s1 = list(prior = prior_name, param = parameter_vector))}.}
#'   \item{\code{w}}{A named list with information on the priors on each split, i.e., all variance proportions.
#'   The names in the list are the split node names from the \code{tree} argument, and is specified
#'   in the same way as the variance priors in \code{V}.}
#' }
#' Prior on residuals can be defined using \code{eps} in this list in the case of a Gaussian likelihood.
#' @param intercept_prior Parameters for Gaussian prior on intercept, specified as a vector with mean and standard deviation.
#' Default is (0, 1000).
#' @param covariate_prior Parameters for Gaussian prior on coefficients of covariates,
#' specified as named list, each element is a vector with mean and standard deviation. Default is (0, 1000).
#' @keywords prior
#' @return Prior object.
#' @details See \link[makemyprior]{makemyprior_models} for details on available priors and likelihoods.
#' @examples
#' \dontrun{
#'
#' vignette("make_prior", package = "makemyprior")
#' }
#'
#' p <- 10
#' m <- 10
#' n <- m*p
#'
#' set.seed(1)
#' data <- list(a = rep(1:p, each = m),
#'              b = rep(1:m, times = p),
#'              x = runif(n))
#' data$y <- data$x + rnorm(p, 0, 0.5)[data$a] +
#'   rnorm(m, 0, 0.3)[data$b] + rnorm(n, 0, 1)
#'
#' formula <- y ~ x + mc(a) + mc(b)
#'
#' prior <- make_prior(formula, data, family = "gaussian",
#'                     intercept_prior = c(0, 1000),
#'                     covariate_prior = list(x = c(0, 100)))
#' prior
#' plot(prior)
#'
#' @export
make_prior <- function(formula, data, family = "gaussian",
                       prior = list(),
                       intercept_prior = c(),
                       covariate_prior = list()){

  hd_prior <- prior
  # get information on each component in the model individually
  # (mw-priors are marked with "mw" and dealt with later)
  prior <- try(decompose_formula(formula, family, data))

  if (is(prior, "try-error"))
    stop("Could not decompose the formula, check that your input is correct and that the names of the data matches the entries in the formula.", call. = FALSE)

  if (length(intercept_prior) == 0) intercept_prior <- c(0, 1000)

  check_if_valid_input_make_prior(intercept_prior, covariate_prior, hd_prior, prior)

  # add priors on fixed effects
  for (ind in seq_along(prior$fixed_effects)){
    if (names(prior$fixed_effects)[ind] %in% names(covariate_prior)){
      prior$fixed_effects[[ind]]$prior$mean.linear <- covariate_prior[[which(names(prior$fixed_effects)[ind] == names(covariate_prior))]][1]
      prior$fixed_effects[[ind]]$prior$sd.linear <- covariate_prior[[which(names(prior$fixed_effects)[ind] == names(covariate_prior))]][2]
    }
  }

  nodes_specified <- sapply(names(prior$random_effects), function(x) is_node_in_tree_string(hd_prior$tree, x) == 0)

  if (!is.null(hd_prior$tree) && nchar(hd_prior$tree) > 0 && any(nodes_specified)) stop("Specify all components in the tree.", call. = FALSE)

  # add info on residual prior
  if (family == "gaussian"){
    if (nchar(hd_prior$tree) > 0 && is_node_in_tree_string(hd_prior$tree, "eps") == 0) stop("Specify all components in the tree, eps is missing.", call. = FALSE)
    prior$random_effects$eps$label <- "eps" # always give residuals name "eps", user cannot use this in the data
    tmp_pr <- 0
    # if the user has eps in the tree without being alone ("(eps)"), the residuals are in a mw-prior
    if (is_node_in_tree_string(hd_prior$tree, "eps") == 2){
      tmp_pr <- list(prior = "mw")
    } else { # if the user has a cw prior on eps
      # using the specified prior, or default if nothing is specified
      if (length(hd_prior$prior$V$eps) == 0){ # no prior specified at all
        hd_prior$prior$V$eps <- list(prior = "pc0", param = default_pc_prior_param("gaussian"))
      }
      tmp_pr <- hd_prior$prior$V$eps
    }
    prior$random_effects$eps$prior <- tmp_pr
  } else if (!(family %in% c("binomial", "poisson"))) {
    stop("This likelihood family is not available.", call. = FALSE)
  }

  # since the user specifies the cw-prior in the V-part, we must move that over to something called cw
  # (consider rewriting, but now now)
  hd_prior$cw <- hd_prior$V[names(hd_prior$V) %in% names(prior$random_effects)]
  hd_prior$V <- hd_prior$V[!(names(hd_prior$V) %in% names(prior$random_effects))]

  # NOTE: tree structure should override prior choices for individual components!
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
    # warning(paste0("Was not able to create the tree requested, using this instead: ",
    #                tree_used
    #                ), call. = FALSE
    # )
    stop("Not a valid tree structure, provide a correct one, or nothing at all.", call. = FALSE)
    # make node_data with the object that is used
    node_data <- make_hd_prior_tree(names(prior$random_effects), tree_used)
  }

  old_new_names <- node_data$old_new_names
  node_data <- node_data[-4]

  prior_data <- make_hd_prior_default(node_data, family) # this is the default choices! (thus only dirichlet on splits)

  # fill in the priors specified by the user, if any
  # requires that the prior-tree is provided, if not we do not know the structure of the prior and we get default everything
  if (!is.null(hd_prior$tree)){
    if (length(hd_prior$w) > 0){
    # if (!is.null(hd_prior$w)){
      # go through prior_data and see if the user has specified a prior
      for (ind in 1:length(prior_data$weights)){
        # check if user has specified a prior for this split node
        which_match <- which(names(hd_prior$w) %in% old_new_names$old[old_new_names$new == prior_data$weight[[ind]]$name])
        if (length(which_match) == 1){
          above_node <- get_basemodel_node_name(prior_data$weights[[ind]], tree, tree_used, old_new_names)
          prior_user <- hd_prior$w[[which_match]] # from user
          if (prior_user$prior %in% c("pc0", "pc1", "pcM") && prior_data$weights[[ind]]$no_children != 2) stop("Cannot use PC prior on multi-splits.", call. = FALSE)
          if (above_node$old != above_node$new && prior_user$prior != "dirichlet"){
            prior_user$prior <- if (prior_user$prior == "pc0") "pc1" else if (prior_user$prior == "pc1") "pc0" else if (prior_user$prior == "pcM") "pcM" else "dirichlet"
            prior_user$param[1] <- 1-prior_user$param[1] # if we switch the above_nodes, we must also invert median
          }
          prior_data$weights[[ind]] <- make_valid_w_prior(prior_user, prior_data$weights[[ind]]) #, get_node_id(node_data, basemodel_node))
        }
      }
    }
    if (length(hd_prior$V) > 0){
    # if (!is.null(hd_prior$V)){
      # go through prior_data and see if the user has specified a prior
      for (ind in 1:length(prior_data$total_variance)){
        which_match <- which(names(hd_prior$V) %in% old_new_names$old[old_new_names$new == prior_data$total_variance[[ind]]$name])
        if (length(which_match) == 1){
          prior_data$total_variance[[ind]] <- make_valid_V_prior(hd_prior$V[[which_match]], prior_data$total_variance[[ind]], hd_prior$V$label)
        }
      }
    }
  }
  if (length(hd_prior$cw) > 0){
  # if (!is.null(hd_prior$cw)){
    # go through prior_data and see if the user has specified a prior
    for (ind in 1:length(prior_data$cw_priors)){
      which_match <- which(names(hd_prior$cw) %in% old_new_names$old[old_new_names$new == prior_data$cw_priors[[ind]]$name])
      if (length(which_match) == 1){
        prior_data$cw_priors[[ind]] <- make_valid_cw_prior(hd_prior$cw[[which_match]], prior_data$cw_priors[[ind]])
      }
    }
  }

  covmats_us <- lapply(prior$random_effects, function(x) {
    if (!is.null(x$model_data$Cmatrix)){
      prec_mat <- try(eval(parse(text = x$model_data$Cmatrix)), silent = TRUE) # if in environment
      if (is(prec_mat, "try-error")) prec_mat <- eval(parse(text = x$model_data$Cmatrix), data) # if not in environment, check if it is in the data
      covmat <- MASS::ginv(as.matrix(prec_mat))
      #typvar <- typical_variance(covmat)
      covmat <- covmat #/typvar
      return(covmat)
      # return(list(covmat = covmat, scaling = if (is.null(x$model_data$scaleCmatrix) || x$model_data$scaleCmatrix) 1 else typvar))
    }
    return(NULL)
  })
  #scaling <- sapply(covmats, function(x) x$scaling)
  #scaling <- scaling[!sapply(scaling, is.null)]
  #covmats <- sapply(covmats, function(x) x$covmat)
  covmats_us <- covmats_us[!sapply(covmats_us, is.null)]

  lower_chol <- lapply(covmats_us, function(x){
    return(t(chol(x + 1e-9*diag(nrow(x)))))
  })

  covmats <- lapply(covmats_us, function(x) return(x/typical_variance(x)))

  # dim_covmats <- sapply(covmats, nrow)
  # if (any(dim_covmats) >= 500) {
  if (length(prior$response) > 500) {
    warning("Computation of the priors may be slow due to the size of the effects.", call. = FALSE)
  }

  # arguments needed for the latent models
  model_data <- list()
  for (ind in seq_len(length(prior$random_effects))){
    if (prior$random_effects[[ind]]$label != "eps"){
      model_data[[ind]] <- c(
        list(id = get_node_id(node_data, prior$random_effects[[ind]]$label)),
        prior$random_effects[[ind]][names(prior$random_effects[[ind]]) != "model_data"],
        prior$random_effects[[ind]]$model_data
      )
      #model_data[[ind]]$scaling <- scaling[[which(prior$random_effects[[ind]]$label == names(scaling))]]
    }
  }

  inference_data <- NULL
  if (family == "binomial") inference_data <- list(Ntrials = data$Ntrials)
  if (family == "poisson"){
    inference_data <- if (is.null(data$E)) list(E = rep(1, length(prior$response))) else list(E = data$E)
  }

  pr <- make_valid_prior_object(
    data = list(fixed = lapply(prior$fixed_effects, function(x) x$data), random = lapply(prior$random_effects, function(x) x$data)),
    args = list(
      node_data = node_data,
      prior_data = prior_data,
      covmats = covmats, # skalerte matriser!
      response_name = prior$response_name,
      use_intercept = prior$use_intercept,
      family = family,
      latent = model_data
    )
  )

  res <- c(pr, list(response = prior$response,
                    tree = tree_used,
                    use_intercept = prior$use_intercept,
                    prior_intercept = c(intercept_prior[1], intercept_prior[2]),
                    prior_fixed = t(sapply(prior$fixed_effects, function(x) x$prior)),
                    formula = formula,
                    family = family,
                    inference_data = inference_data,
                    lower_cholesky = lower_chol))

  class(res) <- "mmp_prior"

  return(res)

}


check_if_valid_input_make_prior <- function(intp, covp, randp, form){

  if (form$use_intercept) {
  if (!is(intp, "numeric") || length(intp) != 2) {
    stop("Wrong input 'intercept_prior'. Provide a numeric vector of length two or nothing.", call. = FALSE)
  } else if (length(intp) == 2 && !form$use_intercept){
    stop("You have provided a prior for the intercept, but there is no intercept in the model.", call. = FALSE)
  }
  }

  for (ind in seq_along(covp)){
    if (!is(covp[[ind]], "numeric") || length(covp[[ind]]) != 2){
      stop("Wrong input 'covariate_prior'. Provide a numeric vector of length two or nothing for each covariate.
           Put it in a named list.", call. = FALSE)
    }
  }

  if (length(randp) > 3){
    stop("Too many elements provided in input argument 'prior'. Only a named list with 'tree', 'V' and 'w'
         is valid.", call. = FALSE)
  } else if (!all(names(randp) %in% c("tree", "V", "w"))) {
    stop("Argument 'prior' can only have names 'tree', 'V' and 'w'.", call. = FALSE)
  } else if (length(randp) == 1 && names(randp) != "tree"){
    stop("If 'V' or 'w' provided, 'tree' must also be provided.", call. = FALSE)
  }

}

# return:
# 0: not in tree at all
# 1: in tree as CW
# 2: in tree in split
is_node_in_tree_string <- function(tree, nodename){

  if (is.null(tree)) {
    return(-1)
  } else if (!is.null(tree) && grepl(nodename, tree)){ #
    if (grepl(paste0("\\(", nodename, "\\)"), tree)){ # if CW node
      return(1)
    } else if (grepl(paste0(",", nodename), tree) || grepl(paste0("\\(", nodename), tree) || grepl(paste0(nodename, ","), tree) || grepl(paste0(nodename, "\\)"), tree)) {
      return(2)
    }
  }

  return(0)

}


get_basemodel_node_name <- function(split_info, tree_old, tree_new, old_new_names){

  # new basemodel is the node that is first in the split with new names
  tree_old <- gsub(" ", "", tree_old)
  tree_new <- gsub(" ", "", tree_new)
  splits_old <- strsplit(strsplit(tree_old, ";")[[1]], "=")
  splits_new <- strsplit(strsplit(tree_new, ";")[[1]], "=")
  oldname <- old_new_names$old[split_info$name == old_new_names$new]
  newname <- split_info$name
  above_old <- sub("\\(", "", strsplit(splits_old[[which(sapply(splits_old, function(x) x[1]) == oldname)]][2], ",")[[1]][1])
  above_new <- sub("\\(", "", strsplit(splits_new[[which(sapply(splits_new, function(x) x[1]) == newname)]][2], ",")[[1]][1])

  return(list(old = above_old, new = above_new))

}

# make_valid_w_prior <- function(new_prior, old_prior, basemodel_id){
make_valid_w_prior <- function(new_prior, old_prior){

  if (!(names(new_prior)[1] == "prior")) {
    if (new_prior[[1]] != "dirichlet" && names(new_prior)[2] != "param"){
      stop("Priors must be named list with elements 'prior' and, if not a dirichlet prior is used, 'param'.", call. = FALSE)
    }
  }
  tmp_prior <- old_prior
  # check that we do not put pc prior on multi-split and that the prior-name is valid (else, default choice)
  if (tmp_prior$no_children > 2) return(old_prior)
  if (new_prior$prior == "dirichlet") return(old_prior) # this will be the same, dirichlet is default
  if (new_prior$prior %in% c("pc0", "pc1", "pcM")){

    median_value <- new_prior$param[1]
    if (median_value < 0 || median_value > 1) stop("Median must be between 0 and 1.", call. = FALSE)

    # # if the node that is not defined as basemodel-node by the "new" tree is the basemodel,
    # # we turn the prior around (PC0(m) = PC1(1-m), PC1(m) = PC0(1-m) and PCM(m) = PCM(1-m))
    ########### this is fixed before this function is used!
    # if (basemodel_id != tmp_prior$children[1]){
    #
    #   if (new_prior$prior != "pcM"){
    #     tmp <- "pc0"
    #     if (new_prior$prior == "pc0") tmp <- "pc1"
    #     new_prior$prior <- tmp
    #   }
    #   basemodel_id <- tmp_prior$children[1]
    #   median_value <- 1-median_value
    #
    # }

    if (!is.na(new_prior$param[2])) {
      if (new_prior$param[2] < 0.5 || new_prior$param[2] >= 1){
        stop("Concentration parameter must be at least 0.5 and less than 1.", call. = FALSE)
      } else conc <- new_prior$param[2]
    } else {
      conc <- 0.5
    }

    tmp_prior$prior <- "pc" #new_prior$prior
    #tmp_prior$param <- new_prior$param
    if (new_prior$prior == "pc0") {
      tmp_prior$param <- data.frame(#basemodel_node = basemodel_id,
                                    above_node = old_prior$children[1], # node above fraction
                                    basemodel = 0, # variance of above_node is 0 in basemodel
                                    median = median_value,
                                    concentration = 0.5 # not used for pc0 or pc1
      )
    } else if (new_prior$prior == "pc1"){
      tmp_prior$param <- data.frame(#basemodel_node = basemodel_id,
                                    above_node = old_prior$children[1],
                                    basemodel = 1,
                                    median = median_value,
                                    concentration = 0.5 # not used for pc0 or pc1
      )
    } else if (new_prior$prior == "pcM"){
      tmp_prior$param <- data.frame(#basemodel_node = basemodel_id,
                                    above_node = old_prior$children[1],
                                    basemodel = median_value,
                                    median = median_value,
                                    concentration = conc
      )
    }
    return(tmp_prior)
  } else {
    stop("Not valid prior name.", call. = FALSE)
  }

  # returning the old prior if the new is not valid
  return(old_prior)

}

make_valid_V_prior <- function(new_prior, old_prior, lab){

  # have not renamed the priors yet, so we change the name of this prior
  # if the user just inputs "pc", it should change the prior to the internal name
  # the old name will be correct, since it is generated by an internal function
  if (new_prior$prior %in% c("pc", "pc0")) {
    new_prior$prior <- "pc0"
  }

  tmp_prior <- old_prior

  # must check that we do not put a jeffreys prior on total variance unless the default has allowed it
  if (new_prior$prior == "jeffreys" && old_prior$prior != "jeffreys"){
    stop(paste0("Jeffreys' prior can only be used for a prior with a single tree and Gaussian data."), call. = FALSE)
  }

  # other than that, we can do "whatever"
  tmp_prior$prior <- new_prior$prior
  if (tmp_prior$prior == "jeffreys"){
    tmp_prior$param <- c(0, 0)
    return(tmp_prior)
  } else {
    if (tmp_prior$prior %in% c("hc", "hn")){
      if (length(new_prior$param) > 1) warning("You have specified two hyperparameters for half-Cauchy/half-normal, but need only one. Using the first one.", call. = FALSE)
      if (new_prior$param[1] <= 0) stop("The hyperparameter for Hald-Cauchy must be positive.", call. = FALSE)
      tmp_prior$param <- c(new_prior$param[1], 0)
    } else {
      if (!(new_prior$prior %in% c("pc0", "invgam"))) stop("Not valid prior family, use 'jeffreys', 'pc', 'hc' or 'invgam'.", call. = FALSE)
      if (length(new_prior$param) != 2) stop(paste0(
        "You must provide two parameters for ",
        if (new_prior$prior == "pc0") "PC prior." else if (new_prior$prior == "invgam") "inverse gamma.", sep = "", collapse = ""
        ), call. = FALSE)
      if (new_prior$prior == "pc0"){
        if (new_prior$param[1] <= 0 || new_prior$param[2] >= 1 || new_prior$param[2] <= 0) stop("Parameters for PC is wrong. They are (U, alpha), U must be positive and alpha between 0 and 1.", call. = FALSE)
      } else {
        if (new_prior$param[1] <= 0 || new_prior$param[2] <= 0) stop("Parameters for inverge gamma is wrong, both must be positive.", call. = FALSE)
      }
      tmp_prior$param <- new_prior$param
    }
    return(tmp_prior)
  }

  stop(paste0("A valid prior for ", old_prior$name, " was not found. See documentation for how to specify prior for variance correct."), call. = FALSE)
  # returning the old prior if the new is not valid
  # return(old_prior)

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
    if (length(new_prior$param) > 1) warning("You have specified two hyperparameters for Half-Cauchy, but need only one. Using the first one.", call. = FALSE)
    if (new_prior$param[1] <= 0) stop("The hyperparameter for Hald-Cauchy must be positive.", call. = FALSE)
    tmp_prior$param <- c(new_prior$param[1], 0)
  } else if (tmp_prior$prior == "invgam") {
    if (length(new_prior$param) != 2) stop("You must provide two parameters for inverse gamma.", call. = FALSE)
    if (new_prior$param[1] <= 0 || new_prior$param[2] <= 0) stop("Parameters for inverge gamma is wrong, both must be positive.", call. = FALSE)
    tmp_prior$param <- new_prior$param
  } else if (tmp_prior$prior == "pc0"){
    if (length(new_prior$param) != 2) stop("You must provide two parameters for PC prior.", call. = FALSE)
    if (new_prior$param[1] <= 0 || new_prior$param[2] >= 1 || new_prior$param[2] <= 0) stop("Parameters for PC is wrong. They are (U, alpha), U must be positive and alpha between 0 and 1.", call. = FALSE)
    tmp_prior$param <- new_prior$param
  } else if (tmp_prior$prior == "jeffreys"){
    stop("Jeffreys' prior cannot be used for an individual variance parameter.", call. = FALSE)
  } else {
    stop(paste0("A valid prior for ", old_prior$name, " was not found. See documentation for how to specify prior for variance correctly."), call. = FALSE)
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

    if (substr(term_labs[ind], 1, 3) == "mc("){

      # run mc-function to extract info from formula
      tmp <- eval(parse(text = term_labs[ind]))
      if (tmp$label %in% c(names(fixed_effects), names(random_effects))) stop("Use unique names for each effect.", call. = FALSE)
      if (tmp$model == "linear"){ # mc can be used for fixed effect
        fixed_effects[[tmp$label]] <- tmp
        fixed_effects[[tmp$label]]$data <- eval(parse(text = tmp$label), data)
      } else {
        random_effects[[tmp$label]] <- tmp
        tmp2 <- eval(parse(text = tmp$label), data)
        if (any(round(tmp2) != tmp2)){
          # warning(paste0("The indexes for random effects must be integers, but the effect ", tmp$label,
          #                " is transformed to integer indexes. This may be done incorrectly, you should check."
          #                ), call. = FALSE)
          # tmp2 <- as.integer(factor(tmp2, levels = unique(tmp2)))
          stop("The indexes for random effects must be integers.", call. = FALSE)
        }
        random_effects[[tmp$label]]$data <- tmp2
      }

    } else { # fixed effect
      # Note that fixed effects are added to stan in a separate data.frame
      fixed_effects[[term_labs[ind]]] <- list(label = term_labs[ind],
                                              model = "linear",
                                              data = eval(parse(text = term_labs[ind]), data))
      fixed_effects[[term_labs[ind]]]$prior <- list(mean.linear = 0, sd.linear = 1000)
    }

  }

  if (length(random_effects) == 0){
    stop("You cannot use 'makemyprior' for (generalized) linear regression, you need to specify at least one random effects.", call. = FALSE)
  }

  # # make data-frame for the fixed effects
  # if (length(fixed_effects) > 0){
  #   fixed_effects_data <- as.data.frame(matrix(unlist(sapply(fixed_effects, function(x) x$data), use.names = FALSE), nrow = n))
  #   names(fixed_effects_data) <- sapply(fixed_effects, function(x) x$label)
  #   fixed_effects_data <- as.list(fixed_effects_data)
  # }
  #
  # # make data-frame for the random effects
  # if (length(random_effects) > 0){
  #   random_effects_data <- as.data.frame(matrix(unlist(sapply(random_effects, function(x) x$data), use.names = FALSE), nrow = n))
  #   names(random_effects_data) <- sapply(random_effects, function(x) x$label)
  #   random_effects_data <- as.list(random_effects_data)
  # } else { # if we just do (generalized) linear regression
  #   random_effects_data <- data.frame()
  # }
  #
  # if (family == "gaussian") random_effects_data$eps <- 1:length(response)

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
  # covmats <- lapply(random_data, function(x) Diagonal(length(unique(x))))
  # covmats <- if (length(args$covmats) == 0) make_covariance_matrices(random_data) else args$covmats
  which_covmats_structured <- c()
  # if some of the covariance matrices are structured
  # for (ind in seq_along(covmats_input)){
  #   index <- which(names(covmats) == names(covmats_input)[ind])
  #   new_covmat <- fix_size_covmat(covmats_input[[ind]], covmats[[index]], random_data[[index]], names(covmats)[index])
  #   covmats[[index]] <- new_covmat
  #   # tmpmat <- prior_obj$covmats[[ind]][1:res$effect_sizes[ind],1:res$effect_sizes[ind]] + 10^(-9)*diag(res$effect_sizes[ind])
  #   # lower_chol <- t(chol(tmpmat/typical_variance(tmpmat)))
  #   # which_covmats_structured <- c(which_covmats_structured, index)
  # }

  # if we have any non-iid components, we need to make a covariance matrix for them
  mods <- c(sapply(args$latent, function(x) x$model))
  non_iids <- which(mods != "iid")

  #orig_covmats <- lapply(data$random, function(x) Diagonal(length(unique(x))))

  for (ind in seq_along(non_iids)){

    if (mods[non_iids[ind]] == "rw1"){
      covmat <- make_rw1_mat(n = length(unique(data$random[[non_iids[ind]]])), type = "cov")
      args$latent[[non_iids[ind]]]$scaling_factor <- typical_variance(covmat)
    } else if (mods[non_iids[ind]] == "rw2"){
      covmat <- make_rw2_mat(n = length(unique(data$random[[non_iids[ind]]])), type = "cov")
      args$latent[[non_iids[ind]]]$scaling_factor <- typical_variance(covmat)
    } else if (mods[non_iids[ind]] == "generic0"){
      # index <- which(names(covmats) == names(covmats_input)[ind])
      # index <- which(names(covmats_input) == names(covmats)[ind]) 
      index <- which(names(covmats_input) == names(covmats)[non_iids[ind]])
      covmat <- covmats_input[[index]]
    } else if (mods[non_iids[[ind]]] == "besag"){
      precmat <- make_besag_prec_mat(args$latent[[non_iids[[ind]]]]$graph)
      covmat <- MASS::ginv(as.matrix(precmat))
      # put the scaling factor into the model data
      args$latent[[non_iids[ind]]]$scaling_factor <- typical_variance(covmat)
      covmat <- covmat/typical_variance(covmat)
      # indexes for Stan (not necessarily required, but we still include it here)
      precmat[!lower.tri(precmat)] <- 0
      ind_both <- which(as.matrix(precmat) == -1, arr.ind = TRUE)
      args$latent[[non_iids[ind]]]$besag_details$ind1 <- ind_both[,2]
      args$latent[[non_iids[ind]]]$besag_details$ind2 <- ind_both[,1]
    } else {
      stop(paste0(mods[non_iids[ind]], " is not available."), call. = FALSE)
    }
    new_covmat <- fix_size_covmat(covmat, covmats[[non_iids[ind]]], data$random[[non_iids[ind]]], names(covmats)[non_iids[ind]])
    # orig_covmats[[non_iids[ind]]] <- covmat/typical_variance(covmat)
    covmats[[non_iids[ind]]] <- new_covmat/typical_variance(new_covmat) # maybe ind and non non_iids[ind]?
    # covmats[[non_iids[ind]]] <- covmat/typical_variance(covmat)
    if (mods[[non_iids[ind]]] == "generic0") which_covmats_structured <- c(which_covmats_structured, non_iids[ind])
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
  initial_args$.latent <- args$latent

  # adding the fixed effect info, for visualization only, is not involved in the prior

  initial_args$.fixef_names <- names(data$fixed)

  # initial_args$.model_eq <- make_model_eq(initial_args)
    # paste("Model equation: ", initial_args$.response_name, " ~ ",
    #       if (initial_args$.use_intercept) expression("mu + ") else "",
    #       if (!is.null(initial_args$.fixef_names)) paste0(paste(initial_args$.fixef_names, sep = "", collapse = " + "), " + ") else "",
    #       paste(args$node_data$orig_nodedata$label,
    #             sep = "", collapse = " + "),
    #       sep = "")

  val <- list(
    prior_data = args$prior_data,
    node_data = args$node_data,
    weight_priors = calculate_pc_prior(args$node_data, args$prior_data$weights, covmats, gui = FALSE),
    .initial_args = initial_args,
    latent = args$latent,
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
# check in make_prior if this function crashes, in that case set the default prior
make_valid_tree <- function(effnames, tree, cw_comps) {

  # removing all spaces
  tree <- gsub(" ", "", tree)
  tree <- gsub("\n", "", tree)

  # ok to have specified a tree without anything, then we check later if this is default MW or CW
  if (is.null(tree) || length(tree) == 0 || tree == "") {
    if (length(cw_comps) > 0){
      # tree <- default_tree(effnames[!(effnames %in% cw_comps)])
      # tree <- paste0(tree, ";", paste0("(", cw_comps, ")", sep = "", collapse = ";"))
      # warning("Did not find a full tree, using default tree structure for the components without specified priors.", call. = FALSE)
      stop("Could not understand the tree, please provide a correct one, or nothing at all.", call. = FALSE)
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
      # warning("Cannot use the same name on splits as the random effects! Using default tree structure instead.",
      #         call. = FALSE)
      # tree <- default_tree(effnames)
      stop("Use names on the split that differ from the names of the random effects!", call. = FALSE)
    }

    effnames_tree <- unlist(strsplit(gsub("\\)", "", gsub("\\(", "", sapply(strsplit(strsplit(tree, ";")[[1]], "="), function(x) x[length(x)]))), ","))
    if (!all(effnames_tree[!(effnames_tree %in% split_names)] %in% effnames)) {
      eff_miss <- effnames_tree[!(effnames_tree %in% split_names)][which(!(effnames_tree[!(effnames_tree %in% split_names)] %in% effnames))]
      stop(paste0("Found the effect ", eff_miss, " in the tree that is not in the data."), call. = FALSE)
    }
    if (sum(table(effnames_tree) > 1) > 0){
      # warning("One or more effects or splits are found in several splits. Using default tree instead.", call. = FALSE)
      # tree <- default_tree(effnames)
      stop("One or more effects or splits are found several places, this is not possible.", call. = FALSE)
    }

    mw_effnames <- effnames_tree[effnames_tree %in% effnames]
    cw_effnames <- effnames[!(effnames %in% mw_effnames)]
    # some CW components
    if (length(mw_effnames) < length(effnames)){
      for (ind in 1:length(cw_effnames)){
        if (!grepl(paste0("(", cw_effnames[ind], ")"), tree)) {
          tree <- paste0(tree, ";", paste0("(", cw_effnames, ")", collapse = ";", sep = ""))
        }
      }
    }

    # if (no_dividers == no_splits && length(cw_effnames) == 0) {
    #   alternative_tree <- paste0(strsplit(tree, ";")[[1]], sep = "", collapse = ";")
    #   warning(paste("Assuming you mean", alternative_tree, "."), call. = FALSE)
    #   tree <- alternative_tree
    # }
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
                                          prior = "pc0", param = default_pc_prior_param(family))
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
                                      #above_node = get_children(node_data, list(current_node_id = ind))[1], # only relevant for pc prior
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
    return(list(prior = "pc0", param = default_pc_prior_param(family)))
  }
}

#' Graphical prior construction
#'
#' This functions opens a shiny app where the specified prior can be seen, and changed.
#'
#' @param prior An object from \link[makemyprior]{make_prior}.
#' @param guide Logical, whether to open the guide directly when the app is started. Default is \code{FALSE}. The guide
#' can be opened in the app at any time.
#' @param no_pc Turn off computation of the PC prior on splits when using the shiny-app, to avoid slow computations.
#' Upon closing, the PC priors will be computed. Default is \code{FALSE}.
#' @keywords gui
#' @export
#' @return Returns an object that can be sent to \link[makemyprior]{inference_stan} or \link[makemyprior]{inference_inla}.
#' Can also be sent to \link[makemyprior]{makemyprior_gui} again.
#' @examples
#' \dontrun{
#'
#' vignette("make_prior", package = "makemyprior")
#'
#' }
#'
#' if (interactive()){
#'    ex_prior <- makemyprior_example_model()
#'    new_prior <- makemyprior_gui(ex_prior)
#' }
#'
#' @importFrom shinyjs useShinyjs disable disabled enable hideElement reset showElement
#' @import shiny
#' @export
makemyprior_gui <- function(prior, guide = FALSE, no_pc = FALSE){

  if (is(prior, "mmp_inla") || is(prior, "mmp_stan")) {
    stop("This is a posterior object. Provide a prior object with class 'mmp_prior'.", call. = FALSE)
  }

  if (!is(prior, "mmp_prior")) stop("Provide an object with class 'mmp_prior'.", call. = FALSE)

  initial_args <- prior$.initial_args
  initial_args$.indata$covmats <- prior$covmats
  initial_args$.guide <- guide
  initial_args$.no_pc <- no_pc

  if (!interactive()) {
    warning("Can only run GUI in interactive session, and this is not an interactive session. Returning the input.")
    return(prior)
  }

  # run application
  shinyjs::useShinyjs()
  # get initial arguments using shiny::getShinyOption(".initial_args") in the app
  shinyOptions(.initial_args = initial_args)
  val <- runApp(shinyApp(ui, server), quiet = TRUE)

  prior$.initial_args <- update_initial_args(initial_args, val)

  res <- prior
  res$prior_data <- val$prior_data
  res$node_data <- val$node_data
  res$weight_priors <- calculate_pc_prior(val$node_data, val$prior_data$weights, prior$covmats, gui = FALSE)
  res$tree <- make_hd_prior_string(val$node_data)
  # res <- c(val, list(weight_priors = calculate_pc_prior(val$node_data, val$prior_data$weights, prior$covmats, gui = FALSE),
  #                    data = prior$data,
  #                    covmats = prior$covmats,
  #                    which_covmats_structured = prior$which_covmats_structured,
  #                    .initial_args = prior$.initial_args,
  #                    response = prior$response,
  #                    latent = prior$latent,
  #                    tree = make_hd_prior_string(val$node_data),
  #                    use_intercept = prior$use_intercept,
  #                    prior_intercept = prior$prior_intercept,
  #                    prior_fixed = prior$prior_fixed,
  #                    # model_data = prior$model_data,
  #                    inference_data = prior$inference_data,
  #                    formula = prior$formula,
  #                    family = prior$family,
  #                    lower_cholesky = prior$lower_cholesky
  #                    ))

  class(res) <- "mmp_prior"

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
  new_initial_args$.no_pc <- NULL

  return(new_initial_args)

}



#' Run inference
#'
#' This function helps you run inference with \link[rstan]{rstan} using a prior object from \link[makemyprior]{make_prior}.
#' Note that you must install Stan: \code{install.packages("rstan")}, see \href{https://mc-stan.org/}{mc-stan.org}.
#' @param prior_obj An object from \link[makemyprior]{make_prior}, from \link[makemyprior]{makemyprior_gui},
#' from \link[makemyprior]{inference_stan}, or from \link[makemyprior]{inference_inla} (for refitting model)
#' @param use_likelihood Whether to sample from the prior only (\code{FALSE}, can be used for e.g.
#' debugging or to look at the priors on variance parameters when using an HD prior, see also Details), or
#' to use the likelihood and data to get the posterior (\code{TRUE}, default).
#' @param print_prior Whether to print a text with the chosen prior or not (default \code{TRUE})
#' @param path Path to folder. See \link[makemyprior]{compile_stan}. Only necessary if compiled code is
#' stored somewhere else than in \code{tempdir()} or the package directory (checking \code{tempdir()} first).
#' @param ... Other arguments to be sent to \link[rstan]{sampling}. Useful arguments include:
#' \describe{
#'   \item{\code{iter}}{number of iterations for each chain (including burn-in, 2000 is the default)}
#'   \item{\code{warmup}}{number of iterations for the burn-in (default is \code{iter}/2)}
#'   \item{\code{chains}}{number of chains}
#'   \item{\code{init}}{initial value of the model parameters on internal parameterization (log-variances and covariate coefficients)}
#'   \item{\code{seed}}{seed value for random number generation}
#' }
#' See \link[rstan]{sampling} for more details. Note that for inference with \code{stan}, the \code{Ntrials} (the size parameter for binomial likelihood) and 
#' \code{E} (the mean E_i in  E_i*exp(eta_i) for Poisson likelihood) argument must be included in the data object, and
#' cannot be provided to this function.
#' @keywords rstan
#' @return A named list with a prior object (\code{prior}), a stan-object (\code{stan}) and some data stan requires (\code{stan_data}).
#' @details We cannot sample from a Jeffreys' prior since it is improper.
#' If \code{use_likelihood = FALSE} and Jeffreys' prior is used for the total variance, the prior will be changed to a Gaussian(0,1) prior on
#' the log total variance. This means that it does not make sense to look at the variances/standard deviations/precisions,
#' but the variance proportions will be correct. Note that this is only an issue when sampling from the prior
#' (i.e., not using the likelihood).
#' @examples
#' \dontrun{
#'
#' vignette("make_prior", package = "makemyprior")
#' }
#'
#' ex_prior <- makemyprior_example_model()
#' if (interactive() && requireNamespace("rstan")){
#'   posterior <- inference_stan(ex_prior, iter = 100)
#'   # Note: For reliable results, increase the number of iterations (e.g., 'iter = 2000')
#'   plot(posterior)
#' }
#'
#' \dontrun{
#'
#' posterior <- inference_stan(ex_prior, use_likelihood = TRUE, iter = 1e4, chains = 1, seed = 1)
#' plot(posterior)
#' }
#'
#'
#' @export
inference_stan <- function(prior_obj, use_likelihood = TRUE, print_prior = TRUE, path = NULL, ...){

  if (!requireNamespace("rstan", quietly = TRUE)) {
    stop("Package \"rstan\" needed for this function to work. Please install it.", call. = FALSE)
  }

  if (print_prior) {print_prior_choice_ranef(prior_obj); cat("\n")}

  if (!is(prior_obj, "mmp_prior")) stop("Provide a prior object from 'make_prior'.", call. = FALSE)

  stan_data <- make_stan_data_object(prior_obj)

  if (is.null(path) && file.exists(paste0(tempdir(), "/full_file.rds"))){ # default in temporary location
    stan_mod <- try(readRDS(paste0(tempdir(), "/full_file.rds")))
    if (is(stan_mod, "try-error")){
      stop("Compile the Stan-code by running 'compile_stan(save = TRUE).", call. = FALSE)
    }
  } else if (is.null(path) && file.exists(paste0(path.package("makemyprior"), "/full_file.rds"))){
    stan_mod <- try(readRDS(paste0(path.package("makemyprior"), "/full_file.rds")))
    if (is(stan_mod,"try-error")){
      stop("Compile the Stan-code by running 'compile_stan(save = TRUE).", call. = FALSE)
    }
  } else {
    if (!is.null(path)){
      stan_mod <- try(readRDS(paste0(path, "/full_file.rds")))
    } else {
      warning(paste0("Could not find the compiled Stan code, so compiling it now. You can compile the code with 'compile_stan', ",
                     "else the code is re-compiled everytime this function is used (which is a bit slow)."), call. = FALSE, immediate. = TRUE)
      stan_mod <- try(compile_stan(save = FALSE))
    }
    if (is(stan_mod, "try-error")) stop("Something went wrong when compiling the Stan-code. See 'compile_stan'.", call. = FALSE)
  }

  if (!use_likelihood) stan_data$likelihood <- 0

  res_stan <- rstan::sampling(
    object = stan_mod,
    data = stan_data,
    ...
  )

  res <- list(prior = prior_obj, stan = res_stan, stan_data =
                list(
                  # needed
                  w_o = stan_data$w_o,
                  w_u = stan_data$w_u,
                  row_index_hd_pr_plot = stan_data$row_index_hd_pr_plot,
                  which_hd = stan_data$which_hd,
                  effect_index_start = stan_data$effect_index_start,
                  effect_sizes = stan_data$effect_sizes,
                  # debug
                  effect_type = stan_data$effect_type,
                  scaling_factor = stan_data$scaling_factor
                )
  )

  class(res) <- "mmp_stan"

  return(res)

}


#' Compile stan-model
#'
#' Function that compiles the stan-code for inference that is included in the model.
#' The compiled version is stored in a .rds-file, which is by default stored in tempdir().
#' Can also be stored in the package (permanent = TRUE), or in a custom directory.
#' In the latter case, this custom directory must be specified every time
#' \link[makemyprior]{inference_stan} is called.
#' This will also be done by \link[makemyprior]{inference_stan}, but this way can save some time if it is not
#' already pre-compiled.
#' @param save Whether to save stan-file to location in package or not, defaults to \code{FALSE}. Must be in interactive session to save the object to a file.
#' @param permanent If \code{TRUE}, the file is stored in the package directory. If \code{FALSE} (default),
#' the file is saved in \code{tempdir()}.
#' @param path Only used if file cannot be saved in the package folder. This is a path to a folder where the file
#' is stored, do not specify a name for the file! (It will be called \code{"full_file.rds"}, and should not be changed.)
#' This argument makes the \code{permanent} argument being ignored.
#' @details Note that you will get a message saying something about integer division. The PC priors on variance proportions
#' are represented by splines, and to evaluate them in Stan we look up values, and use integer division for this. This does
#' not cause problems.
#' @return Returns the stan-model invisibly.
#' @examples
#'
#' if (interactive() && requireNamespace("rstan")){
#'   compile_stan(save = TRUE) # saving in tempdir()
#' }
#'
#' @export
compile_stan <- function(save = FALSE, permanent = FALSE, path = NULL){

  if (!requireNamespace("rstan", quietly = TRUE)) {
    stop("Package \"rstan\" needed for this function to work. Please install it.", call. = FALSE)
  }

  stan_file <- paste0(path.package("makemyprior"), "/full_file.stan")

  if (!is.null(path)){
    if (save) warning("You are now saving the compiled Stan-object somewhere outside the package. You must provide this path to 'inference_stan' every time you use that function.")
  } else {
    if (permanent){
      path <- path.package("makemyprior")
    } else {
      path <- tempdir()
    }
  }

  stan_mod <- rstan::stan_model(file = stan_file, auto_write = FALSE)

  if (save && interactive()) saveRDS(stan_mod, paste0(path, "/full_file.rds"))
  if (!interactive() && save) stop("Cannot save the compiled stan-file if not in interactive session.", call. = FALSE)

  invisible(stan_mod)

}

make_stan_data_object <- function(prior_obj){

  # in case y has some weird format
  y <- as.numeric(unlist(prior_obj$response))

  add_res <- if (prior_obj$family == "gaussian") 1 else 0

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
  if (length(effects$eps) == 0 && prior_obj$family == "gaussian"){
    # if (prior_obj$likelihood == "gaussian" && !is.null(effects$eps)){
    effects <- effects[which(names(effects) != "eps")]
  }
  # in case the data is not provided as indexes
  for (ind in which(!sapply(effects, is.integer))){
    effects[[ind]] <- as.integer(as.factor(effects[[ind]]))
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
  res$which_hd <- as.array(c(1:(res$n_r+add_res))[prior_obj$node_data$orig_nodedata$status == "attached"])
  res$which_cw <- as.array(c(1:(res$n_r+add_res))[prior_obj$node_data$orig_nodedata$status == "detached"])
  res$n_cw <- length(res$which_cw)

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
  which_theta_in_hd_mat <- matrix(0, nrow = res$n_totvar, ncol = res$n_r+add_res)
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
  weight_info_over <- matrix(0, nrow = length(prior_obj$node_data$edges$from), ncol = res$n_r+add_res)
  weight_info_under <- matrix(0, nrow = length(prior_obj$node_data$edges$from), ncol = res$n_r+add_res)
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
    row_index_hd_pr_mat[tree_no, ind3] <- 1
  }

  row_index_hd_pr <- array(0, dim = 1)

  # this loop goes over each weight (each HD prior in the prior)
  for (ind2 in seq_len(length(unique(hd_pr_num)))){
    row_index_hd_pr <- c(row_index_hd_pr, which(hd_pr_num == ind2)[length(which(hd_pr_num == ind2))])
  }

  row_index_hd_pr <- row_index_hd_pr # TODO: why do we do this? Something missing?

  res$row_index_hd_pr <- row_index_hd_pr_mat # which weight has which rows in the weight_info-matrices
  res$row_index_hd_pr_plot <- row_index_hd_pr + 1 # which weight has which rows in the weight_info-matrices

  res$w_o <- weight_info_over
  res$w_u <- weight_info_under

  res$use_intercept <- as.integer(prior_obj$use_intercept)

  res$effect_sizes <- array(sapply(effects, function(x) length(unique(x))),
                            dim = c(res$n_r))
  res$effect_index_start <- cumsum(res$effect_sizes)-res$effect_sizes
  res$indexes <- as.data.frame(effects)

  res$n_knots <- sapply(prior_obj$weight_priors, function(x) x$n_knots)
  if (!all(res$n_knots %in% c(2, 122))) stop("All splines for the priors must have the same amount of knots.", call. = FALSE)
  res$n_knots <- 122

  res$knots <- c(-500,
                 seq(-200, -40, length.out = 10)[1:9],
                 seq(-40, -5, length.out = 22)[1:21],
                 seq(-5, -1e-6, length.out = 30),
                 seq(1e-6, 5, length.out = 30),
                 seq(5, 40, length.out = 22)[-1],
                 seq(40, 200, length.out = 10)[-1],
                 500) # these are always the same

  for (ind in seq_along(prior_obj$weight_priors)){
    if (prior_obj$weight_priors[[ind]]$n_knots == 122){
      if (any(res$knots != prior_obj$weight_priors[[ind]]$knots))
        stop("The knots differ from what they should be, run again and do not alter the splines for the PC prior!", call. = FALSE)
    }
  }

  tmp_array <- array(0, dim = c(res$n_pc, 122, 4))
  pc_ind <- if (res$n_hd > 0) which(sapply(prior_obj$weight_priors, function(x) length(x$knots) > 2)) else c()
  for (i in seq_len(length(pc_ind))){
    tmp_array[i,,] <- prior_obj$weight_priors[[pc_ind[i]]]$coeffs
  }
  res$prior_coeffs <- tmp_array

  # intercept/fixed effects priors
  res$intercept_prior <- array(prior_obj$prior_intercept, dim = c(2))
  res$fixed_priors <- array(as.numeric(prior_obj$prior_fixed), dim = c(res$n_f, 2))

  ####### data/likelihood info

  if (prior_obj$family == "binomial" && length(prior_obj$inference_data$Ntrials) != res$n) stop("Ntrials must be provided, and have the same length as observations.", call. = FALSE)
  if (prior_obj$family == "poisson" && length(prior_obj$inference_data$E) != res$n) stop("E must have the same length as observations if provided.", call. = FALSE)

  # the format of the data depends on the likelihood, sending it both as integers and as reals
  res$y_cont <- y
  res$y_disc <- as.integer(y)
  res$Ntrials <- if (prior_obj$family == "binomial") as.integer(prior_obj$inference_data$Ntrials) else rep(0, res$n)
  res$e_disc <- if (prior_obj$family == "poisson") as.integer(prior_obj$inference_data$E) else rep(0, res$n)
  res$likelihood <- if (prior_obj$family == "gaussian") 1 else if (prior_obj$family == "binomial") 2 else if (prior_obj$family == "poisson") 3 else stop("Not valid family.", call. = FALSE)
  res$add_res <- add_res


  ############### latent things (not relevant for the prior)

  # which latent model is each random effect
  res$effect_type <- array(sapply(prior_obj$latent, model_to_number))

  # used for generic0
  # which_covmats_structured <- c()
  mega_matrix <- array(0, dim = c(res$n_r, res$n, res$n))
  for (i in seq_along(prior_obj$which_covmats_structured)){
    ind <- prior_obj$which_covmats_structured[i]
    # lower_chol <- t(chol(prior_obj$unscaled_covmats[[ind]][1:res$effect_sizes[ind],1:res$effect_sizes[ind]] + 10^(-9)*diag(res$effect_sizes[ind])))
    mega_matrix[ind, 1:res$effect_sizes[ind], 1:res$effect_sizes[ind]] <- prior_obj$lower_cholesky[[i]] # lower_chol
    # which_covmats_structured <- c(which_covmats_structured, ind)
  }
  res$mega_matrix <- mega_matrix
  res$use_index_matrix <- array(0, dim = res$n_r)
  res$use_index_matrix[prior_obj$which_covmats_structured] <- 1

  res$scaling_factor <- array(sapply(prior_obj$latent, function(x) if (is.null(x$scale.model) || !x$scale.model) 1 else x$scaling_factor))

  res$besag_ind1 <- array(0, 0)
  res$besag_ind2 <- array(0, 0)
  # used for besag
  for (ind in which(res$effect_type %in% c(8, 9))){
    res$besag_ind1 <- prior_obj$latent[[ind]]$besag_details$ind1
    res$besag_ind2 <- prior_obj$latent[[ind]]$besag_details$ind2
  }
  res$n_besag <- length(res$besag_ind1)

  return(res)

}


# input is model data for one latent component (for Stan)
model_to_number <- function(model_data){

  if (model_data$model == "iid"){
    if (!model_data$constr) return(1) else return(2)
  } else if (model_data$model == "rw2"){
    if (!model_data$constr) return(3) else if (!model_data$lin_constr) return(4) else return(5)
  } else if (model_data$model == "generic0"){
    if (!model_data$constr) return(6) else return(7)
  } else if (model_data$model == "besag"){
    if (!model_data$constr) return(8) else return(9)
  } else if (model_data$model == "rw1"){
    if (!model_data$constr) return(10) else return(11)
  }
  return(-1)

}


# making matrix for rw1-effect
# dim = size of effect
# type = precision matrix or covariance matrix ("prec", "cov" or "lower")
make_rw1_mat <- function(n, type = c("prec", "cov", "lower")){

  # precision matrix
  prec_mat <- matrix(0, nrow = n, ncol = n)

  # diagonal
  diag(prec_mat) <- c(1, rep(2, n-2), 1)

  # off-diagonal 1
  prec_mat[matrix(c(1:(n-1), 2:n, 2:n, 1:(n-1)), ncol = 2)] <- -1

  if (type == "prec") return(prec_mat)

  covmat <- MASS::ginv(prec_mat)
  if (type == "cov") return(covmat)

  lower <- t(chol(covmat))
  if (type == "lower") return(lower)

  return(diag(n))

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
  prec_mat[matrix(c(1,2,2,1,n-1,n,n,n-1), ncol = 2, byrow = TRUE)] <- -2

  # off-diagonal 2
  prec_mat[matrix(c(1:(n-2), 3:n, 3:n, 1:(n-2)), ncol = 2)] <- 1

  if (type == "prec") return(prec_mat)

  covmat <- MASS::ginv(prec_mat)
  if (type == "cov") return(covmat)

  lower <- t(chol(covmat))
  if (type == "lower") return(lower)

  return(diag(n))

}


# making precision matrix for besag effect from graph file (must be on same format as required by INLA)
# input is the file path to the graph file
make_besag_prec_mat <- function(graph_path){

  stopifnot(file.exists(graph_path))

  all_lines_char <- readLines(graph_path)
  all_lines <- sapply(all_lines_char, function(x) as.numeric(strsplit(x, split = " ")[[1]]))
  n <- all_lines[[1]]

  mat <- Matrix(0, n, n)

  for (ind in 1:n){

    # how many -1's on this row:
    no_neighb <- all_lines[[ind+1]][2]
    stopifnot(length(all_lines[[ind+1]])-2 == no_neighb)
    if (no_neighb > 0) mat[ind, all_lines[[ind + 1]][c(-1, -2)]] <- -1
    mat[ind, ind] <- no_neighb

  }

  return(mat)

}


# calculate the covariance matrices for each effect, to use to make the PC priors
# all matrices are of the Matrix-class, so they are sparse
make_covariance_matrices <- function(data){

  matrix_data <- lapply(data, make_covariance_matrix)
  matrix_data <- lapply(matrix_data, scale_cov_matrix)
  names(matrix_data) <- names(data)

  return(matrix_data)

}


# index is the index vector for each random effect
make_covariance_matrix <- function(index){

  
  # if this is a vector of data and not indexes, we make an index vector instead. but this should give an error somewhere else
  if (any(round(index) != index)){
    stop("The random effects must have integer indexes.", call. = FALSE)
  }

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


#' Compute the typical variance
#'
#' Computing the typical variance (geometric mean) of a matrix.
#' @param mat Matrix.
#' @return Typical variance.
#' @examples
#' typical_variance(diag(10))
#'
#' @export
typical_variance <- function(mat){

  this_diag <- Matrix::diag(mat)

  if (any(this_diag < 0)){
    stop("Some diagonal elements of this covariance matrix are negative. Change prior or reformulate your model.", call. = FALSE)
  }

  return(exp(mean(log(this_diag))))

}

check_covariance_matrices <- function(mats){

  n <- nrow(mats[[1]])

  for (i in 1:length(mats)){
    tmp <- if (compute_rank(mats[[i]]) < n) "singular" else "not singular"
    cat(names(mats)[i], "is", tmp, "\n")
  }

}

# if the effect has shorter length than the number of observations, it will be
# resized according to the indexes (so dependencies between each element in the
# covariance matrix is still the same)
# should only be for un-structured effects??
fix_size_covmat <- function(new_mat, old_mat, index, eff_name){

  # old_dim >= new_dim, since old_mat has dimension nxn (n = number of observations)
  new_dim <- dim(new_mat)[1]
  old_dim <- dim(old_mat)[1]

  # if (new_dim != old_dim){
  #   if ((old_dim/new_dim)%%1 != 0){
  #     stop(paste0(
  #       "Wrong dimension of input covariance matrix for effect ", eff_name
  #       #", using one based on index input!"
  #     ))
  #     return(old_mat)
  #   } else {
  #     return(
  #       kronecker(diag(old_dim/new_dim), new_mat)
  #     )
  #   }
  # } else { # if the dimensions are equal, no resizing necessary
  #   return(new_mat)
  # }

  new_mat <- as.matrix(new_mat)
  old_mat <- as.matrix(old_mat)

  if (new_dim != old_dim){

    n <- length(index)
    ret_mat <- matrix(0, n, n)

    for (i_r in 1:n){
      for (i_c in 1:n){
        ret_mat[i_r, i_c] <- new_mat[index[i_r], index[i_c]]
      }
    }
    return(ret_mat)
  } else return(new_mat)

}

max_pr_num <- function() return(6)

get_variance_prior_number <- function(prior_names){

  pr_names <- c("jeffreys", "pc0", "invgam", "hc", "hn", "")

  pr_num <- c()
  for (i in 1:length(prior_names)){
    pr_num[i] <- which(pr_names == prior_names[i])
  }

  pr_num[pr_num == max_pr_num()] <- 0

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

# making a tree consisting of only the nodes below the split_id node
make_subtree <- function(node_data, split_id){

  # remove all nodes above the split node, make totally new tree with fewer nodes, but with
  # original indexes and names

  below <- get_nodes_below(node_data, split_id)
  nodes_in_tree <- c(below, split_id)

  node_data_new <- list(
    nodes = node_data$nodes[node_data$nodes$id %in% nodes_in_tree,],
    edges = node_data$edges[node_data$edges$from %in% nodes_in_tree,],
    orig_nodedata = node_data$orig_nodedata[node_data$orig_nodedata$id %in% nodes_in_tree,]
  )

  return(node_data_new)

}













#' Run inference
#'
#' This function helps you run inference with INLA using a prior object from \link[makemyprior]{make_prior}.
#' You must have INLA installed to run this. INLA can be installed with:
#' \code{install.packages("INLA", repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"), dep = TRUE)}.
#' Also see \href{https://www.r-inla.org/}{r-inla.org}.
#' @param prior_obj An object from \link[makemyprior]{make_prior}, from \link[makemyprior]{makemyprior_gui}, from \link[makemyprior]{inference_stan}, or from \link[makemyprior]{inference_inla} (for refitting model)
#' @param use_likelihood Whether to sample from the prior only (FALSE, can be used for e.g. debugging or to look at the priors on variance parameters when using an HD prior,
#' see also Details), or to use the likelihood and data to get the posterior (TRUE, default).
#' @param print_prior Whether to print a text with the chosen prior or not (default TRUE)
#' @param ... Other values to be sent to INLA.
#' Useful arguments include \code{Ntrials} for the binomial likelihood and \code{E} (mean E_i in E_i*exp(eta_i)) for the Poisson likelihood.
#' We set the default value of \code{num.threads} to 1 (this can however be changed).
#' See \link[INLA]{inla} for details.
#' Can be anything sent to \link[INLA]{inla} except for \code{control.expert} and arguments that specify priors.
#' @keywords INLA
#' @return A named list with a prior object (\code{prior}), an inla-object (\code{inla}) and some data inla requires (\code{inla_data}).
#' @details Jeffreys' prior is improper. If \code{use_likelihood = FALSE} and Jeffreys' prior is used for the total variance, the prior will be changed to a Gaussian(0,1) prior on
#' the log total variance. This means that it does not make sense to look at the variances/standard deviations/precisions,
#' but the variance proportions will be correct. Note that this is only an issue when sampling from the prior
#' (i.e., not using the likelihood).
#' @export
#' @examples
#' \dontrun{
#'
#' vignette("make_prior", package = "makemyprior")
#' }
#'
#' ex_prior <- makemyprior_example_model()
#' if (interactive() && requireNamespace("INLA")){
#'   posterior <- inference_inla(ex_prior)
#'   plot(posterior)
#' }
#'
#' @export
inference_inla <- function(prior_obj, use_likelihood = TRUE,
                           print_prior = TRUE, ...){

  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop("Package \"INLA\" needed for this function to work. Please install it, see www.r-inla.org.", call. = FALSE)
  }

  if (print_prior) {print_prior_choice_ranef(prior_obj); cat("\n")}

  if (!is(prior_obj, "mmp_prior")) stop("Provide a prior object from 'make_prior'.", call. = FALSE)

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

  covmats <- NULL
  # if Cmatrix is missing in the environment for a generic0-model, but is in the data-object, add it to the data here
  generic0_mod_mats <- sapply(prior_obj$latent[sapply(prior_obj$latent, function(x) x$model == "generic0")], function(x) c(x$Cmatrix, x$label))
  if (length(generic0_mod_mats) > 0){
    for (ind in seq_len(ncol(generic0_mod_mats))){
      if (is(try(eval(parse(text = generic0_mod_mats[1,ind]))[1], silent = TRUE), "try-error")){
        tmp <- prior_obj$covmats[names(prior_obj$covmats) == generic0_mod_mats[2, ind]] # this is now a list
        names(tmp) <- generic0_mod_mats[1, ind]
        tmp[[1]] <- MASS::ginv(tmp[[1]]) # to make precision matrix
        inla_data <- c(inla_data, tmp)
      }
    }
  }

  input_args <- list(...)

  inla_jpr_data <- make_inla_data_object(prior_obj, input_args)

  # replace jeffreys prior with gaussian prior on total variance if no likelihood
  if (!use_likelihood) {
    inla_jpr_data$totvar_prior_info[inla_jpr_data$totvar_prior_info[,1] == 1,1] <- max_pr_num()
    inla_jpr_data$totvar_prior_numbers[inla_jpr_data$totvar_prior_numbers == 1] <- max_pr_num()
  }

  if (!is.null(input_args$control.expert)){
    stop("You cannot provide the INLA::control.expert argument, it is used for making the joint prior.", call. = FALSE)
  }
  if (!is.null(input_args$control.family$hyper)){
    stop("You cannot specify priors for parameters in this function, use make_prior instead.", call. = FALSE)
  }

  if (!is.null(input_args$control.fixed$mean) || !is.null(input_args$control.fixed$prec) ||
      !is.null(input_args$control.fixed$mean.intercept) || !is.null(input_args$control.fixed$prec.intercept)){
    stop("You cannot specify a prior here, use make_prior instead.", call. = FALSE)
  }

  if (is.null(input_args$num.threads)) input_args$num.threads <- 1

  inla_res <- do.call(
    INLA::inla, c(list(
      formula = inla_formula,
      data = inla_data,
      family = prior_obj$family,
      control.fixed = list(mean.intercept = prior_obj$prior_intercept[1], prec.intercept = 1/prior_obj$prior_intercept[2]^2),
      control.expert = list(jp = make_jpr(inla_jpr_data))
    ), input_args)
  )

  res <- list(
    prior = prior_obj,
    inla = inla_res,
    inla_data = list(
      formula = inla_formula,
      w_o = inla_jpr_data$w_o,
      w_u = inla_jpr_data$w_u,
      row_index_hd_pr_plot = inla_jpr_data$row_index_hd_pr_plot,
      which_hd = inla_jpr_data$which_hd)
  )

  class(res) <- "mmp_inla"

  return(res)

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

    form_constr <- form_extraconstr <- form_rankdef <- extra_arguments <- scalemod <- NULL

    #initial_val <- paste0("initial = ", prior_obj$latent[[ind]]$initial)

    if (prior_obj$latent[[ind]]$label == effnames_random[ind]) {

      mod_type <- prior_obj$latent[[ind]]$model
      if (prior_obj$latent[[ind]]$constr) form_constr <- "constr = TRUE"
      if (length(prior_obj$latent[[ind]]$extraconstr) > 0) form_extraconstr <- paste0("extraconstr = ", prior_obj$latent[[ind]]$extraconstr)
      if (length(prior_obj$latent[[ind]]$rankdef) > 0) form_rankdef <- paste0("rankdef = ", prior_obj$latent[[ind]]$rankdef)
      if (length(prior_obj$latent[[ind]]$scale.model) > 0) scalemod <- paste0("scale.model = ", prior_obj$latent[[ind]]$scale.model)

      if (mod_type == "generic0"){
        extra_arguments <- paste0("Cmatrix = ", prior_obj$latent[[ind]]$Cmatrix)
      }

      if (mod_type == "besag"){
        extra_arguments <- paste0("graph = \"", prior_obj$latent[[ind]]$graph, "\"")
      }

    } else {
      stop(paste0(effnames_random[ind], " did not get correct model type."), call. = FALSE)
      # warning(paste0(effnames_random[ind], " did not get correct model type! Setting 'iid'."))
      # mod_type <- "iid"
    }

    mod_type <- paste0("model = \"", mod_type, "\"")

    comps_r[ind] <- paste0("f(",
                           paste0(c(effnames_random[ind], mod_type, form_constr,
                                  form_extraconstr, form_rankdef, scalemod, extra_arguments), #, initial_val),
                                  collapse = ", ", sep = ""),
                           ")")

  }

  if (length(comps_r) > 0) {
    if (!prior_obj$use_intercept || length(comps_f) > 0){ # either removed intercept or added fixed effects, then the formula is not empty
      inla_formula <- paste0(inla_formula, " + ", paste0(comps_r, collapse = " + "))
    } else { # nothing added to the formula yet
      inla_formula <- paste0(inla_formula, paste0(comps_r, collapse = " + "))
    }
  }

  inla_formula <- formula(inla_formula)

  return(inla_formula)

}



# fixed effects can be a data.frame or a list
make_inla_data_object <- function(prior_obj, input_args){

  if (length(input_args) > 0 && length(input_args$Ntrials) > 0 && length(prior_obj$inference$Ntrials) == 0) prior_obj$inference_data$Ntrials <- input_args$Ntrials
  if (length(input_args) > 0 && length(input_args$E) > 0 && length(prior_obj$inference$E) == 0) prior_obj$inference_data$E <- input_args$E

  return(make_stan_data_object(prior_obj))

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


#' expit
#'
#' Calculates inverse logit, exp(x)/(1+exp(x)) = 1/(1+exp(-x))
#' @param x Real number, or vector of reals
#' @return expit value
#' @examples
#' expit(2)
#' expit(seq(-5, 5, 1))
#'
#' @export
expit <- function(x) 1/(1+exp(-x))

#' logit
#'
#' Calculates logit, log(x/(1-x)) = log(x) - log(1-x)
#' @param x Value between 0 and 1, or vector of such values
#' @return logit value
#' @examples
#' logit(0.5)
#' logit(seq(0, 1, 0.1))
#'
#' @export
logit <- function(x) log(x/(1-x))


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
  # logdens <- logdens -log( expit(logitw)*(1-expit(logitw)) )
  logdens <- logdens -(logitw - 2*log1p(exp(logitw))) # more numerically stable for very small and very large numbers


  return(logdens)

}


# evaluates all priors that has HD priors in a tree (run several times for multiple trees)
# includes total variance
# includes jacobians!
# theta are all log variances involved in this prior tree (may have more than one tree)
hd_prior_joint_lpdf <- function(theta, likelihood, which_pc, w_o, w_u, n_splits_tree, row_index_hd_pr, knots, prior_coeffs, n_knots, totvar_prior_info, indV,
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
  if (likelihood == 0 && totvar_prior_info[indV,1] == 1){ # must have proper variance when we sample from prior
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

  if (prior_number == 1){ # jeffreys
    return(0)
  } else if (prior_number == 2){ # pc
    shape <- -log(param[2])/param[1]; # param = c(U, alpha)
    return(
      log(shape/2) + x/2 - shape*exp(x/2)
    )
  } else if (prior_number == 3){ # inverse gamma(shape,scale)
    return(
      -lgamma(param[1]) + param[1]*log(param[2]) -param[1]*x - param[2]/(exp(x))
    )
  } else if (prior_number == 4){ # hc(scale)
    return(
      x/2 - log(pi) - log(param[1] + exp(x)/param[1])
    )
  } else if (prior_number == 5) { # half-normal(scale) (stdev in normal)
    return(
      -0.5*log(2*pi) - log(param[1]) + x/2 - exp(x)/(2*param[1]^2)
    )
  } else if (prior_number == max_pr_num()) {
    return(dnorm(x, 0, 1, log = TRUE));
  } else {
    return(0)
  }

}


# args_123 (args_123ysg35ovat) contains the data, which lives in the same environment as this function and the
# functions necessary for the prior to be computed
joint_prior <- function(theta_prec){

  args_123ysg35ovat <- args_123ysg35ovat # hack to avoid note about variable not existing
  args_123 <- args_123ysg35ovat # weird name to hopefully avoid problems with global environment

  # transform from log precision to log variance
  # if Gaussian likelihood, put residual variance at the end
  theta <- if (args_123$add_res == 1) -c(theta_prec[-1], theta_prec[1]) else -theta_prec

  logdens <- 0

  # for each tree in the prior:
  if (length(args_123$which_hd) > 0){
    for (indV in seq_len(args_123$n_totvar)){
      logdens <- logdens + hd_prior_joint_lpdf(
        theta[get_indexes(args_123$which_theta_in_hd[indV,])],
        args_123$likelihood,
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
make_jpr <- function(inla_data){

  return(
    INLA::inla.jp.define(joint_prior,
                         # data
                         # using "unique" name to avoid errors if user has a variable with the same name in global environment
                         args_123ysg35ovat = inla_data,
                         # functions
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
# used to pick correct thetas for the logit weights etc (basically just matching functions from stan-code)
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




# return text for the expression for the prior distribution on the parameter provided (similar to the latex-version)
get_prior_expr_text <- function(prior_data, node_data, param = c("cw", "totvar", "weight")){

  param <- match.arg(param)

  if (param == "cw" || param == "totvar"){

    if (param == "cw") param_type <- if (prior_data$prior %in% c("pc0", "hc")) "sigma" else "sigma^2"
    if (param == "totvar") param_type <- if (prior_data$prior %in% c("pc0", "hc", "hn")) "sqrt(V)" else "V"
    param_name <- paste0(param_type, "[", prior_data$name, "]")


    if (prior_data$prior == "pc0") {
      prior_expr <- paste0("PC0(", prior_data$param[1], ", ", prior_data$param[2], ")\n")
    } else if (prior_data$prior == "jeffreys"){
      prior_expr <- paste0("Jeffreys'\n")
    } else if (prior_data$prior == "invgam"){
      prior_expr <- paste0("InvGam(", prior_data$param[1], ", ", prior_data$param[2], ")\n")
    } else if (prior_data$prior == "hc"){
      prior_expr <- paste0("Half-Cauchy(", prior_data$param[1], ")\n")
    } else if (prior_data$prior == "hn"){
      prior_expr <- paste0("Half-normal(", prior_data$param[1], ")\n")
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
      t_over <- as.character(get_node_name(node_data, prior_data$param$above_node))
      param_name <- sprintf("w[%s/%s]", t_over, t_under)

      if (prior_data$param$basemodel == 0){
        prior_expr <- sprintf("PC0(%s)\n", prior_data$param$median)
      } else if (prior_data$param$basemodel == 1) {
        prior_expr <- sprintf("PC1(%s)\n", prior_data$param$median)
      } else {
        prior_expr <- sprintf("PCM(%s, %s)\n", prior_data$param$median, prior_data$param$concentration)
      }

    }

  }

  return(
    paste0(param_name, " ~ ", prior_expr)
  )

}



#' Find suitable PC prior parameters
#'
#' Returns the \code{U} value in P(\code{U} > sigma) = \code{alpha} for a PC prior on standard deviation given an
#' equal-tailed credible interval
#' P(\code{lower} < \code{func}(x) < \code{upper}) = \code{prob} where x is a Gaussian variable with
#' zero mean standard deviation sigma.
#' Note that this function uses sampling.
#' @param lower lower end of credible interval
#' @param upper upper end of credible interval
#' @param alpha tail probability of the PC prior (default = \code{0.05})
#' @param func function to scale Gaussian variables to match desired credible interval scale,
#' default is the exponential function
#' @param N number of samples to use when sampling sigma and x, default is \code{1e4}
#' @param prob amount of mass to put in the credible interval, default is \code{0.95}
#' @return The \code{U}-value to pass to the PC prior. NB! Store result to avoid rerunning this function, as it uses sampling.
#' The function also prints (sampled) quantiles for the U-value that is returned.
#' @examples find_pc_prior_param(0.1, 10)
#'
#' @export
find_pc_prior_param <- function(lower, upper, alpha = 0.05, func = exp, N = 1e4, prob = 0.95){

  optfun <- function(logU){
    stdev <- rexp(N, rate = -log(alpha)/exp(logU))
    samps <- rnorm(N, 0, sd = stdev)
    return(sum(abs(quantile(func(samps), probs = probs) - c(lower, upper))))
  }

  probs <- c((1-prob)/2)
  probs[2] <- 1-probs[1]

  U <- exp(optimize(optfun, lower = -1e6, upper = 10)$minimum)

  resfun <- function(U){
    stdev <- rexp(N, rate = -log(alpha)/U)
    samps <- rnorm(N, 0, sd = stdev)
    return(quantile(func(samps), probs = probs))
  }

  true_quants <- resfun(U)

  cat("U = ", U, "\n", sep = "")
  cat("Prob(", true_quants[1], " < ", deparse(substitute(func)), "(eta) < ", true_quants[2], ") = ", prob, "\n", sep = "")

  return(invisible(U))

}



default_pc_prior_param <- function(family){

  if (family == "gaussian") return(c(3, 0.05))
  if (family == "binomial") return(c(1.6, 0.05))
  if (family == "poisson") return(c(1.6, 0.05))

}


## prints random effect prior choices with text
print_prior_choice_ranef <- function(prior_obj){

  if (is(prior_obj, "mmp_inla") || is(prior_obj, "mmp_stan")) prior_obj <- prior_obj$prior

  cat("Tree structure:", prior_obj$tree)
  cat("\n\n")

  if (length(prior_obj$prior_data$weights) > 0) cat("Weight priors:\n")
  for (ind in seq_along(prior_obj$prior_data$weights)){
    cat("\t", get_prior_expr_text(prior_obj$prior_data$weights[[ind]], prior_obj$node_data, "weight"), sep = "")
  }

  if (length(prior_obj$prior_data$total_variance) > 0) cat("Total variance priors:\n")
  for (ind in seq_along(prior_obj$prior_data$total_variance)){
    cat("\t", get_prior_expr_text(prior_obj$prior_data$total_variance[[ind]], prior_obj$node_data, "totvar"), sep = "")
  }

  if (any(sapply(prior_obj$prior_data$cw_priors, function(x) nchar(x$prior)))) cat("Independent variance priors:\n")

  for (ind in seq_along(prior_obj$prior_data$cw_priors)){
    if (prior_obj$prior_data$cw_priors[[ind]]$prior != "")
      cat("\t", get_prior_expr_text(prior_obj$prior_data$cw_priors[[ind]], prior_obj$node_data, "cw"), sep = "")
  }

}

## prints fixed effect prior choices with text
print_prior_choice_fixef <- function(prior_obj){

  if (is(prior_obj, "mmp_inla") || is(prior_obj, "mmp_stan")) prior_obj <- prior_obj$prior

  cat("Covariate priors: ")
  if (!prior_obj$use_intercept && nrow(prior_obj$prior_fixed) == 0) cat("No covariates.")
  covariates <- row.names(prior_obj$prior_fixed)
  if (prior_obj$use_intercept) {
    cat("intercept ~ N(", prior_obj$prior_intercept[1], ", ", prior_obj$prior_intercept[2], "^2)", sep = "")
    if (length(covariates) > 0) cat(", ")
  }
  for (ind in seq_along(covariates)){
    cat(covariates[ind], " ~ N(", prior_obj$prior_fixed[ind,]$mean.linear, ", ", prior_obj$prior_fixed[ind,]$sd.linear, "^2)", sep = "")
    if (ind < length(covariates)) cat(", ")
  }

}




#' Scaling precision matrix
#'
#' Scaling a precision matrix so the corresponding covariance matrix has typical variance (geometric mean) equal to 1.
#' @param Q Precision matrix
#' @return Precision matrix which is now scaled to have typical variance 1.
#' @examples
#' scale_precmat(diag(10))
#'
#' @export
scale_precmat <- function(Q){

  A <- MASS::ginv(Q)
  A2 <- A/typical_variance(A)
  Q2 <- MASS::ginv(A2)
  return(Q2)

}





######### printing, summary, plot function

#' Print
#'
#' @method print mmp_prior
#' @param x Object of class \code{mmp_prior}, \code{mmp_inla} or \code{mmp_stan}.
#' @param ... For \code{mmp_stan} and \code{mmp_inla}, see \link[base]{print.data.frame}.
#' @return Returns input object invisible.
#' @examples
#' pri <- makemyprior_example_model()
#' pri # or print(pri)
#'
#' if (interactive() && requireNamespace("rstan")){
#'   res_stan <- inference_stan(ex_prior, iter = 100)
#'   # Note: For reliable results, increase the number of iterations (e.g., 'iter = 2000')
#'   res_stan # or print(res_stan)
#' }
#'
#' if (interactive() && requireNamespace("INLA")){
#'   res_inla <- inference_inla(pri)
#'   res_inla # or print(res_inla)
#' }
#'
#' @export
print.mmp_prior <- function(x, ...){

  cat("Model: ")
  print(x$formula)
  print_prior_choice_ranef(x)

  invisible(x)

}

#' @method print mmp_inla
#' @rdname print.mmp_prior
#' @export
print.mmp_inla <- function(x, ...){

  args <- list(...)
  if (is.null(args$row.names)) args$row.names <- FALSE
  if (is.null(args$right)) args$right <- FALSE
  if (is.null(args$digits)) args$digits <- 3

  cat("Model: ")
  print(x$prior$formula)
  cat("Tree structure:", x$prior$tree, "\n")
  cat("\nInference done with INLA\n\n")

  tmp <- make_posterior_summary_inla(x)
  do.call(print.data.frame, c(list(x = tmp), args))

  invisible(tmp)

}

#' @method print mmp_stan
#' @rdname print.mmp_prior
#' @export
print.mmp_stan <- function(x, ...){

  args <- list(...)
  if (is.null(args$row.names)) args$row.names <- FALSE
  if (is.null(args$right)) args$right <- FALSE
  if (is.null(args$digits)) args$digits <- 3

  cat("Model: ")
  print(x$prior$formula)
  cat("Tree structure:", x$prior$tree, "\n")
  cat("\nInference done with Stan.\n\n")

  tmp <- make_posterior_summary_stan(x)
  do.call(print.data.frame, c(list(x = tmp), args))

  invisible(tmp)

}


#' Short summary
#'
#' @method summary mmp_prior
#' @param object Object of class \code{mmp_prior}, \code{mmp_inla} or \code{mmp_stan}.
#' @param ... For \code{mmp_stan} and \code{mmp_inla}, see \link[base]{print.data.frame}.
#' @return Returns summary invisible.
#' @examples
#' pri <- makemyprior_example_model()
#' summary(pri)
#'
#' if (interactive() && requireNamespace("rstan")){
#'   res_stan <- inference_stan(ex_prior, iter = 100)
#'   # Note: For reliable results, increase the number of iterations (e.g., 'iter = 2000')
#'   summary(res_stan)
#' }
#'
#' if (interactive() && requireNamespace("INLA")){
#'   res_inla <- inference_inla(pri)
#'   summary(res_inla)
#' }
#'
#' @export
summary.mmp_prior <- function(object, ...){

  cat("Model: ")
  print(object$formula)
  print_prior_choice_ranef(object)
  cat("\n")
  print_prior_choice_fixef(object)
  cat("\n")

  invisible(object)

}

#' @method summary mmp_inla
#' @rdname summary.mmp_prior
#' @export
summary.mmp_inla <- function(object, ...){

  args <- list(...)
  if (is.null(args$row.names)) args$row.names <- FALSE
  if (is.null(args$right)) args$right <- FALSE
  if (is.null(args$digits)) args$digits <- 3

  cat("Model: ")
  print(object$prior$formula)
  print_prior_choice_ranef(object)
  cat("\n")
  print_prior_choice_fixef(object)
  cat("\n")
  cat("\nInference done with INLA.\n\n")

  tmp <- make_posterior_summary_inla(object)
  do.call(print.data.frame, c(list(x = tmp), args))

  invisible(tmp)

}

#' @method summary mmp_stan
#' @rdname summary.mmp_prior
#' @export
summary.mmp_stan <- function(object, ...){

  args <- list(...)
  if (is.null(args$row.names)) args$row.names <- FALSE
  if (is.null(args$right)) args$right <- FALSE
  if (is.null(args$digits)) args$digits <- 3

  cat("Model: ")
  print(object$prior$formula)
  print_prior_choice_ranef(object)
  cat("\n")
  print_prior_choice_fixef(object)
  cat("\n")
  cat("\nInference done with Stan.\n\n")

  tmp <- make_posterior_summary_stan(object)
  do.call(print.data.frame, c(list(x = tmp), args))

  invisible(tmp)

}


#' Plotting
#'
#' @method plot mmp_prior
#' @param x Object of class \code{mmp_prior}, \code{mmp_inla} or \code{mmp_stan}.
#' @param ... Additional arguments to plotting functions. Varies with what object is sent to function.
#' @return None.
#' @details See \link[makemyprior]{plot_prior} (objects of class \code{mmp_prior}),
#' \link[makemyprior]{plot_posterior_stan} (objects of class \code{mmp_stan}), and
#' \link[makemyprior]{plot_posterior_variance} (objects of class \code{mmp_inla}),
#' @examples
#' pri <- makemyprior_example_model()
#' plot(pri)
#'
#' if (interactive() && requireNamespace("rstan")){
#'   res_stan <- inference_stan(ex_prior, iter = 100)
#'   # Note: For reliable results, increase the number of iterations (e.g., 'iter = 2000')
#'   plot(res_stan)
#' }
#'
#' if (interactive() && requireNamespace("INLA")){
#'   res_inla <- inference_inla(pri)
#'  plot(res_inla)
#' }
#'
#' @export
plot.mmp_prior <- function(x, ...){

  return(plot_prior(x))

}

#' @method plot mmp_inla
#' @rdname plot.mmp_prior
#' @export
plot.mmp_inla <- function(x, ...){

  return(plot_posterior_variance(x))

}

#' @method plot mmp_stan
#' @rdname plot.mmp_prior
#' @export
plot.mmp_stan <- function(x, ...){

  return(plot_posterior_stan(x, ...))

}


#' List available priors, latent models and likelihoods
#'
#' Prints available priors, latent models and likelihoods to use with \link[makemyprior]{make_prior}.
#' @param type Which of priors, latent models and likelihoods to list. Options are
#' \code{prior}, \code{latent} and \code{likelihood}
#' @param select Which in each group to show details about.
#' \code{NULL} gives only list what exists (default),
#' \code{all} gives detailed information about everything in that category.
#' Can also ask for one or more specific priors/latent models/likelihoods
#' @return None.
#' @examples
#' makemyprior_models("prior", c("pc0", "pc1"))
#' makemyprior_models("latent")
#' makemyprior_models("likelihood", "all")
#'
#' @export
makemyprior_models <- function(type = c("prior", "latent", "likelihood"), select = NULL){

  type <- match.arg(type)

  if (type == "prior") {
    makemyprior_options_prior(select)
  } else if (type == "latent"){
    makemyprior_options_latent(select)
  } else if (type == "likelihood"){
    makemyprior_options_likelihood(select)
  } else stop("Only 'prior', 'latent' and 'likelihood' are valid inputs.", call. = FALSE)

}

# internal, listing priors
makemyprior_options_prior <- function(prs = NULL){

  priors_var <- data.frame(
    internal_name = c("pc", "hc", "hn", "invgam", "jeffreys"),
    full_name = c("PC_0", "Half-Cauchy", "Half-normal", "Inverse Gamma", "Jeffreys'"),
    which_param = c("st.dev.", "st.dev.", "st.dev", "variance", "variance"),
    params = c("U and alpha in Prob(stdev > U) = alpha", "scale", "scale", "shape and scale", "none"),
    example = c("'list(prior = \"pc\", param = c(U, alpha))'",
                "'list(prior = \"hc\", param = scale)'",
                "'list(prior = \"hn\", param = scale)'",
                "'list(prior = \"invgam\", param = c(shape, scale))'",
                "'list(prior = \"jeffreys\")'")
  )

  priors_w <- data.frame(
    internal_name = c("pc0", "pc1", "pcM", "dirichlet"),
    full_name = c("PC_0", "PC_1", "PC_M", "Dirichlet"),
    constr = c(rep("Can only be used on a dual split.\n", 3), ""),
    params = c("median", "median", "median and concentration parameter", "none"),
    example = c(
      "'list(prior = \"pc0\", param = median)'",
      "'list(prior = \"pc1\", param = median)'",
      "'list(prior = \"pcM\", param = c(median, concentration))'",
      "'list(prior = \"dirichlet\")'"
    )
  )

  if (is.null(prs)){
    cat(
      "Variance priors: ", paste0("\"", priors_var$internal_name, "\"", collapse = ", "),
      ".\nVariance proportion priors: ", paste0("\"", priors_w$internal_name, "\"", collapse = ", "),
      ".",
      sep = ""
    )
    cat("\n")
  } else {
    if (prs[1] == "all") prs <- c(priors_var$internal_name, priors_w$internal_name)

    if (any(prs %in% priors_var$internal_name)) cat("\n############ Variance priors: ############\n\n")

    for (ind in which(priors_var$internal_name %in% prs)){
      cat(
        "--------  ", priors_var[ind,"full_name"], "  --------\n",
        "Internal name: \"", priors_var[ind, "internal_name"],
        if (priors_var[ind,"internal_name"] == "pc") "\", also possible to use \"pc0\".\n" else "\".\n",
        "This is a prior on a ", priors_var[ind, "which_param"], " parameter.\n",
        "Required hyperparameters: ", priors_var[ind, "params"], ".\n",
        "Example of use: ", priors_var[ind, "example"]
        ,
        sep = "")
      if (priors_var[ind, "internal_name"] == "invgam"){
        cat(
          "\nNote that InvGam(shape1, scale1) on variance corresponds tp Gam(shape2, rate2) on precision, where",
          "shape1 = shape2 and scale1 = rate2."
          )
      }
      cat("\n\n")
    }

    if (any(prs %in% priors_w$internal_name)) cat("\n############ Variance proportion priors: ############\n\n")

    for (ind in which(priors_w$internal_name %in% prs)){
      cat(
        "--------  ", priors_w[ind,"full_name"], "  --------\n",
        "Internal name: \"", priors_w[ind, "internal_name"], "\".\n",
        priors_w[ind, "constr"],
        "Required hyperparameters: ", priors_w[ind, "params"], ".\n",
        "Example of use: ", priors_w[ind, "example"]
        ,
        sep = "")
      cat("\n\n")
    }

    cat("For more details, see 'vignette(\"make_prior\", package = \"makemyprior\")'.\n")

  }

}

# internal, listing latent models
makemyprior_options_latent <- function(lts = NULL){

  latent <- data.frame(
    internal_name = c("iid", "linear", "besag", "rw1", "rw2", "generic0"),
    full_name = c("Independent and identically distributed", "Fixed effect.", "Besag", "Random walk of 1st order",
                  "Random walk of 2nd order", "Structured covariance matrix"),
    options = c("'constr': sum-to-zero, default is TRUE.",
                "",
                "'constr': sum-to-zero, default is TRUE. 'graph': path to graph file, see '?mc' for details.",
                "'constr': sum-to-zero, default is TRUE.",
                "'constr': sum-to-zero, default is TRUE. 'lin_constr': linear constraint, default is FALSE.",
                "'constr': sum-to-zero, default is FALSE. 'Cmatrix': Precision matrix."
                ),
    details = c("Default option. Sum-to-zero constraint is used by default ('constr = TRUE').",
                "For specifying fixed effects (can also just add effect to formula).",
                "",
                "",
                "",
                "We recommend to scale the precision matrix so the typical variance of the effect is 1, see '?mc' and '?scale_precmat'."
                )
  )

  if (is.null(lts)){
    cat("Latent models: ", paste0("\"", latent$internal_name, "\"", collapse = ", "), sep = "")
    cat("\n")
  } else {
    if (lts[1] == "all") lts <- latent$internal_name

    cat("\n############ Latent models: ############\n\n")

    for (ind in which(latent$internal_name %in% lts)){
      cat(
        "--------  ", latent[ind,"full_name"], "  --------\n",
        "Internal name: \"", latent[ind, "internal_name"], "\".\n",
        "Options: ", latent[ind, "options"],
        if (latent[ind, "options"] != "") "\n" else "",
        latent[ind, "details"],
        if (latent[ind, "details"] != "") "\n" else ""
        ,
        sep = "")
      cat("\n")
    }

    # cat("For more details, see 'vignette(\"make_prior\", package = \"makemyprior\")'.\n")

  }

}

# internal, listing likelihoods
makemyprior_options_likelihood <- function(lks = NULL){

  likelihoods <- data.frame(
    internal_name = c("gaussian", "binomial", "poisson"),
    full_name = c("Gaussian", "Binomial", "Poisson"),
    details = c(paste0("Default. A residual effect will be added automatically, do not add this yourself.\n",
                       "Default variance prior is PC(", default_pc_prior_param("gaussian")[1], ", ",
                       default_pc_prior_param("gaussian")[2], "), except for total variance in the case of only one tree, ",
                       "then a Jeffreys' prior is used.\n"),
                paste0("Need to add number of trials with name 'Ntrials' to the data object for",
                "inference with Stan, and to 'inference_inla' for inference with inla.\n",
                "Default variance prior is PC(", default_pc_prior_param("binomial")[1], ", ",
                       default_pc_prior_param("binomial")[2], ").\n"),
                paste0("If the mean E_i (in E_i * exp(eta_i) where eta_i is the linear predictor) is not 1 for all observations, ",
                       "it must be provided to the data object for inference with Stan, and to 'inference_inla' for inference with inla.\n",
                  "Default variance prior is PC(", default_pc_prior_param("poisson")[1], ", ",
                default_pc_prior_param("poisson")[2], ").\n"))
  )

  if (is.null(lks)){
    cat("Likelihoods: ", paste0("\"", likelihoods$internal_name, "\"", collapse = ", "), sep = "")
    cat("\n")
  } else {
    if (lks[1] == "all") lks <- likelihoods$internal_name

    cat("\n############ Likelihoods: ############\n\n")

    for (ind in which(likelihoods$internal_name %in% lks)){
      cat(
        "--------  ", likelihoods[ind,"full_name"], "  --------\n",
        "Internal name: \"", likelihoods[ind, "internal_name"], "\".\n",
        likelihoods[ind, "details"],
        sep = "")
      cat("\n")
    }

    cat("For more details, see 'vignette(\"make_prior\", package = \"makemyprior\")'.\n")

  }

}


#' List of available plotting functions
#'
#' Directs to this help page.
#'
#' @return None.
#' @details The available plotting functions are:
#' \itemize{
#'   \item{\link[makemyprior]{plot_marginal_prior}}
#'   \item{\link[makemyprior]{plot_posterior_fixed}}
#'   \item{\link[makemyprior]{plot_posterior_precision}}
#'   \item{\link[makemyprior]{plot_posterior_stan}}
#'   \item{\link[makemyprior]{plot_posterior_stdev}}
#'   \item{\link[makemyprior]{plot_posterior_variance}}
#'   \item{\link[makemyprior]{plot_prior}}
#'   \item{\link[makemyprior]{plot_several_posterior_stan}}
#'   \item{\link[makemyprior]{plot_tree_structure}}
#' }
#' In addition the following functions can be used to extract posterior samples of
#' effects and parameters:
#' \itemize{
#'   \item{\link[makemyprior]{extract_posterior_effect}}
#'   \item{\link[makemyprior]{extract_posterior_parameter}}
#' }
#'
#' \link[makemyprior]{eval_pc_prior} can be used to evaluate a PC prior for a weight parameter,
#' and \link[makemyprior]{eval_joint_prior} to evaluate the whole joint prior.
#'
#' @export
makemyprior_plotting <- function() help("makemyprior_plotting")


#' Returning a simple example prior object
#'
#' Creating a simple prior object using \link[makemyprior]{make_prior}.
#' Used in examples of other functions in the package.
#' @param seed A seed value for reproducing the data (default \code{seed = 1}).
#' @return An object of class \code{mmp_prior}.
#' @details See the example for what model is made.
#' @examples
#'
#' ex_model <- makemyprior_example_model()
#'
#' \dontrun{
#' # The function corresponds to the following model call:
#'
#' set.seed(1)
#'
#' data <- list(
#'   a = rep(1:10, each = 10),
#'   b = rep(1:10, times = 10)
#' )
#' data$y <- rnorm(10, 0, 0.4)[data$a] + rnorm(10, 0, 0.6)[data$b] + rnorm(100, 0, 1)
#'
#' formula <- y ~ mc(a) + mc(b)
#'
#' prior <- make_prior(formula, data, family = "gaussian",
#'                     prior = list(tree = "s1 = (a, b); s2 = (s1, eps)",
#'                                  w = list(s2 = list(prior = "pc0", param = 0.25)),
#'                                  V = list(s2 = list(prior = "pc", param = c(3, 0.05)))),
#'                     intercept_prior = c(0, 1000))
#' }
#'
#' @export
makemyprior_example_model <- function(seed = 1){

  set.seed(seed)

  data <- list(
    a = rep(1:10, each = 10),
    b = rep(1:10, times = 10)
  )
  data$y <- rnorm(10, 0, 0.4)[data$a] + rnorm(10, 0, 0.6)[data$b] + rnorm(100, 0, 1)

  formula <- y ~ mc(a) + mc(b)

  prior <- make_prior(formula, data, family = "gaussian",
                      prior = list(tree = "s1 = (a, b); s2 = (s1, eps)",
                                   w = list(s2 = list(prior = "pc0", param = 0.25)),
                                   V = list(s2 = list(prior = "pc", param = c(3, 0.05)))),
                      intercept_prior = c(0, 1000))

  return(prior)

}




#' Evaluate the joint variance prior
#'
#' Function for evaluating the joint variance prior stored in \code{prior_obj}. To compute the joint prior, the functions needs
#' to know the transformation from the total variance/variance proportion scale to log-variance scale.
#' This is computed before inference, but is not stored in the \code{mmp_prior}-object.
#' To avoid having to recompute this for every evaluation and thus improve the speed, we make a condensed data object with
#' the function \link[makemyprior]{make_eval_prior_data}.
#' @param prior_obj Object of class \code{mmp_prior}, see \link[makemyprior]{make_prior}.
#' @param theta Vector of log variances. The order of the log variances is
#' the same as specified in the formula, with the residual variance at the end for a Gaussian likelihood. To be sure,
#' you can use \link[makemyprior]{get_parameter_order} to check the order.
#' @param prior_data An object from \link[makemyprior]{make_eval_prior_data}.
#' @return Logarithm of the prior density.
#' @details Note that a Jeffreys' prior is improper and sampling with the prior only will not work when it
#' is used. For sampling from the prior (for example for debugging), use a proper prior for all parameters instead.
#'
#' The
#' \code{make_eval_prior_data} function is used to create a condensed version of the prior object from
#' \code{make_prior}, that only contains what is needed to compute the joint prior. Since the HD prior is chosen on
#' total variances and variance proportions, some additional information is needed
#' to compute the Jacobian for the joint prior. To improve the speed, we do this once before evaluating the prior.
#' 
#' Expert option: \code{make_eval_prior_data} can also be used to extract the prior to be used with 'regular' INLA. See
#' examples for how this can be done.
#' @examples
#'
#' ex_model <- makemyprior_example_model()
#' get_parameter_order(ex_model) # a, b, eps
#' prior_data <- make_eval_prior_data(ex_model)
#' eval_joint_prior(c(0, 0, 0), prior_data)
#' eval_joint_prior(c(-1, 0, 1), prior_data)
#'
#' # a model with only 2 variance parameters
#' if (interactive()){
#'
#'   data <- list(
#'     a = rep(1:10, each = 10)
#'   )
#'   set.seed(1); data$y <- rnorm(10, 0, 0.4)[data$a] + rnorm(100, 0, 1)
#'
#'   # random intercept model
#'   ex_model2 <- make_prior(y ~ mc(a), data, family = "gaussian",
#'                           prior = list(tree = "s2 = (a, eps)",
#'                                        w = list(s2 = list(prior = "pc0", param = 0.25)),
#'                                        V = list(s2 = list(prior = "pc", param = c(3, 0.05)))),
#'                           intercept_prior = c(0, 1000))
#'
#'   prior_data2 <- make_eval_prior_data(ex_model2)
#'   # evaluating the prior in a grid
#'   theta_a <- seq(-8, 4, 0.1)
#'   theta_eps <- seq(-8, 4, 0.1)
#'   res <- matrix(nrow = 0, ncol = 3)
#'   for (ind in 1:length(theta_a)){
#'     for (jnd in 1:length(theta_eps)){
#'       res <- rbind(res, c(theta_a[ind], theta_eps[jnd],
#'                           eval_joint_prior(c(theta_a[ind], theta_eps[jnd]), prior_data2)))
#'     }
#'   }
#'
#'   # graph showing the prior
#'   if (requireNamespace("ggplot2")){
#'     res2 <- as.data.frame(res)
#'     names(res2) <- c("x", "y", "z")
#'     # Note from the "exp(z)" that we use the posterior, and not log posterior, in this plot
#'     ggplot(res2, aes(x = x, y = y, z = exp(z), fill = exp(z))) +
#'       geom_raster() +
#'       geom_contour(color = "black") +
#'       scale_fill_viridis_c(option = "E") +
#'       xlab("Log variance of 'a'") +
#'       ylab("Log residual variance") +
#'       labs(fill = "Density") +
#'       theme_bw()
#'   }
#'
#' }
#' 
#' \dontrun{
#' 
#' # How an HD prior can be computed with \code{makemyprior}, and then sent to regular INLA 
#' # (expert option).
#' # Note the use of the hidden \code{make_jpr}-function.
#' # Also note that the order of the parameters must be the same as in the call to \code{make_prior}.
#' # The residual variance is put in the correct place by \code{make_jpr}.
#' data <- list(
#'   a = rep(1:10, each = 100),
#'   b = rep(1:100, times = 10)
#' )
#' set.seed(1); data$y <- rnorm(100, 0, 0.4)[data$a] + rnorm(100, 0, 0.6)[data$b] + rnorm(1000, 0, 1)
#' prior <- make_prior(y ~ mc(a) + mc(b), data, family = "gaussian",
#'                     prior = list(tree = "s1 = (a, b); s2 = (s1, eps)",
#'                                  w = list(s2 = list(prior = "pc0", param = 0.25)),
#'                                  V = list(s2 = list(prior = "pc", param = c(3, 0.05)))),
#'                     intercept_prior = c(0, 1000))
#' jpr_dat <- make_eval_prior_data(prior)
#' res <- inla(y ~ f(a) + f(b),
#'             data = data,
#'             control.fixed = list(prec.intercept = 1/1000^2),
#'             control.expert = list(jp = makemyprior:::make_jpr(jpr_dat)))
#' }
#' @export
eval_joint_prior <- function(theta, prior_data){

  if (is(prior_data, "mmp_prior")) stop("Use 'make_eval_prior_data()' and send the result to this function. See ?eval_joint_prior.", call. = FALSE)

  logdens <- 0

  for (indV in seq_len(prior_data$n_totvar)){
    logdens <- logdens + hd_prior_joint_lpdf(
      theta[get_indexes(prior_data$which_theta_in_hd[indV,])],
      prior_data$likelihood,
      prior_data$which_pc,
      prior_data$w_o[,get_indexes(prior_data$which_theta_in_hd[indV,])],
      prior_data$w_u[,get_indexes(prior_data$which_theta_in_hd[indV,])],
      prior_data$n_splits_each_tree[indV],
      get_indexes(prior_data$row_index_hd_pr[indV,]),
      prior_data$knots,
      prior_data$prior_coeffs,
      prior_data$n_knots,
      prior_data$totvar_prior_info,
      indV,
      # sending the functions here, since we do this for INLA and re-using the same functions
      # reduces the possibilities for errors at a later stage
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

  # adding individual priors (CW priors)
  if (length(prior_data$which_cw) > 0){
    logdens <- logdens + cw_priors_lpdf(
      theta[prior_data$which_cw],
      matrix(prior_data$cw_prior_info[prior_data$which_cw,], ncol = 3),
      choose_prior_lpdf
    )
  }

  return(logdens)

}


#' @rdname eval_joint_prior
#' @export
make_eval_prior_data <- function(prior_obj){

  dat <- make_stan_data_object(prior_obj)

  res <- list(
    which_theta_in_hd = dat$which_theta_in_hd,
    likelihood = dat$likelihood,
    which_pc = dat$which_pc,
    w_o = dat$w_o,
    w_u = dat$w_u,
    n_splits_each_tree = dat$n_splits_each_tree,
    row_index_hd_pr = dat$row_index_hd_pr,
    knots = dat$knots,
    prior_coeffs = dat$prior_coeffs,
    n_knots = dat$n_knots,
    totvar_prior_info = dat$totvar_prior_info,
    which_cw = dat$which_cw,
    cw_prior_info = dat$cw_prior_info
  )

  return(dat)

}


#' Internal variance parameter order
#'
#' Returns the internal order of the variance parameters related to each random effect in the model.
#' @param prior_obj Object of class \code{mmp_prior}, see \link[makemyprior]{make_prior}.
#' @return Names of the random effects in the model in the order the prior object reads them.
#'
#' @examples
#' ex_model <- makemyprior_example_model()
#' get_parameter_order(ex_model)
#'
#' @export
get_parameter_order <- function(prior_obj) return(prior_obj$node_data$orig_nodedata$label)

















