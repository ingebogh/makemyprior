

# plotting the prior and posterior of each effect

# some of these are also used in the shiny-app

# for one component, either CW or total variance
# all priors are plotted on variance scale
df_prior_variance <- function(prior_data, node_data, var_type){

  if (is.null(prior_data$prior)) return(data.frame())

  if (prior_data$prior == "pc0"){
    xx <- c(seq(0, 20, 0.1))
    yy <- dexp(sqrt(xx), rate = -log(prior_data$param[2])/prior_data$param[1]) * 1/(2*sqrt(xx))
  } else if (prior_data$prior == "invgam"){
    a <- prior_data$param[1]
    b <- prior_data$param[2] # input is rate, but we plot using scale = 1/rate
    xx <- seq(0, 1e4, 100)
    yy <- b^a/(gamma(a)) * (1/xx)^(a+1) * exp(-b/xx)
  } else if (prior_data$prior == "hc"){
    xx <- sqrt(seq(0, 1e4, 100))
    yy <- 2*prior_data$param[1]/(pi*(xx + prior_data$param[1]^2)) * 1/(2*sqrt(xx))
  } else if (prior_data$prior == "jeffreys"){
    xx <- seq(0, 3, 0.01)
    yy <- 1/xx
  } else {
    return(data.frame()) # if this component does not have a CW prior
  }

  title <- if (var_type == "total") bquote(V[.(paste0(
    get_node_name(node_data, get_leaf_nodes(node_data, prior_data$id)), sep = "", collapse = ","))]) else bquote(sigma[.(prior_data$name)]^2)

  df <- data.frame(x = xx, y = yy, param = as.character(as.expression(title)))

  return(df)

}

df_prior_variance_custom <- function(xx, prior_data, use_sd = FALSE){

  if (is.null(prior_data$prior)) return(data.frame())

  if (prior_data$prior == "pc0"){
    if (use_sd){
      yy <- dexp(xx, rate = -log(prior_data$param[2])/prior_data$param[1])
    } else {
      yy <- dexp(sqrt(xx), rate = -log(prior_data$param[2])/prior_data$param[1]) * 1/(2*sqrt(xx))
    }
  } else if (prior_data$prior == "invgam"){
    a <- prior_data$param[1]
    b <- prior_data$param[2] # input is rate, but we plot using scale = 1/rate
    if (use_sd){
      yy <- 2*b^a/(gamma(a)) * (1/xx)^(2*a+1) * exp(-b/xx^2)
    } else {
      yy <- b^a/(gamma(a)) * (1/xx)^(a+1) * exp(-b/xx)
    }
  } else if (prior_data$prior == "hc"){
    if (use_sd){
      yy <- 2*prior_data$param[1]/(pi*(xx^2 + prior_data$param[1]^2))
    } else {
      yy <- 2*prior_data$param[1]/(pi*(xx + prior_data$param[1]^2)) * 1/(2*sqrt(xx))
    }
  } else if (prior_data$prior == "jeffreys"){
    yy <- 1/xx # does not matter if it is stdev or not, prior is still the same
  } else {
    return(data.frame()) # if this component does not have a CW prior
  }

  df <- data.frame(x = xx, y = yy)

  return(df)

}


# for one split
plot_diri_w <- function(prior_data, node_data){

  no_children <- prior_data$no_children
  param <- get_diri_param(no_children)

  titles <- character()
  under <- paste0(get_node_name(node_data, get_leaf_nodes(node_data, prior_data$id)), sep = "", collapse = ",")
  for (ind in 1:(no_children-1)){
    titles[ind] <- as.expression(bquote(omega[frac(.(as.character(paste0(
      get_node_name(node_data, get_leaf_nodes(node_data, prior_data$children[ind])), sep = "", collapse = ","))), .(under))]))
  }

  xx <- seq(0, 1, 0.01)
  df <- data.frame(x = xx, y = dbeta(xx, param, param*(no_children-1)),
                   param = rep(as.character(titles), each = length(xx)))

  return(df)

}

# for one split
plot_pc_w <- function(prior_data, weight_data, node_data){

  # names of the child nodes, so it is clear where the base model for the PC prior is
  under <- paste0(as.character(get_node_name(node_data, get_leaf_nodes(node_data, prior_data$id))), sep = "", collapse = ",")
  over <- paste0(as.character(get_node_name(node_data, get_leaf_nodes(node_data, prior_data$param$above_node))), sep = "", collapse = ",")
  title <- as.expression(bquote(omega[frac(.(over), .(under))]))

  xx <- seq(0, 1, 0.01)
  yy <- eval_spline_prior(xx, weight_data, FALSE)

  df <- data.frame(x = xx, y = yy, param = as.character(title))

  return(df)

}

# return just the name of the weight
plot_pc_w_nocalc <- function(prior_data, node_data){

  # names of the child nodes, so it is clear where the base model for the PC prior is
  under <- paste0(as.character(get_node_name(node_data, get_leaf_nodes(node_data, prior_data$id))), sep = "", collapse = ",")
  over <- paste0(as.character(get_node_name(node_data, get_leaf_nodes(node_data, prior_data$param$above_node))), sep = "", collapse = ",")
  title <- as.expression(bquote(omega[frac(.(over), .(under))]))

  return(data.frame(x = rep(prior_data$param$median, 10), y = seq(0, 1, length.out = 10), param = as.character(title)))

}


# making data-frame for plotting prior or posteriors, or both for Stan
make_dataframe_for_plotting <- function(type = c("prior", "posterior", "both"), res){

  if (is(res, "mmp_stan")){
    prior_obj <- res$prior
    res_stan <- res$stan
    stan_data <- res$stan_data
  } else if (is(res, "mmp_prior")){
    prior_obj <- res
    type <- "prior" # can only plot prior if that is the only object that is provided
  } else {
    stop("Must provide object of class 'mmp_stan' or 'mmp_prior'.", call. = FALSE)
  }

  type <- match.arg(type)

  if (type == "prior" || type == "both"){

    df_pri <- data.frame()

    # for CW priors:
    for (i in 1:length(prior_obj$prior_data$cw_priors)){
      df_pri <- rbind(df_pri, df_prior_variance(prior_obj$prior_data$cw_priors[[i]], prior_obj$node_data, "cw"))
    }

    # for total variances
    for (i in seq_len(length(prior_obj$prior_data$total_variance))){
      tmp <- df_prior_variance(prior_obj$prior_data$total_variance[[i]], prior_obj$node_data, "total")
      if (!(type == "both" && prior_obj$prior_data$total_variance[[i]]$prior == "jeffreys")){
        # don't plot jeffreys prior for posterior total variance
        df_pri <- rbind(df_pri, tmp)
      }
    }

    # for HD priors
    for (i in seq_len(length(prior_obj$prior_data$weights))){
      if (prior_obj$prior_data$weights[[i]]$prior == "pc"){
        df_pri <- rbind(df_pri, plot_pc_w(prior_obj$prior_data$weights[[i]], prior_obj$weight_priors[[i]], prior_obj$node_data))
      } else { # dirichlet
        df_pri <- rbind(df_pri, plot_diri_w(prior_obj$prior_data$weights[[i]], prior_obj$node_data))
      }
    }

  }

  if (type == "posterior" || type == "both"){

    df_post <- data.frame()

    samps <- rstan::extract(res_stan, "theta")$theta

    # for CW priors:
    for (i in 1:length(prior_obj$prior_data$cw_priors)){

      if (prior_obj$prior_data$cw_priors[[i]]$prior != ""){

        xx <- exp(samps[,i])
        title <- as.expression(bquote(sigma[.(prior_obj$prior_data$cw_priors[[i]]$name)]^2))

        df_post <- rbind(df_post, data.frame(x = xx, param = as.character(title)))

      }

    }

    # for total variances
    for (i in seq_len(length(prior_obj$prior_data$total_variance))){

      if (is.null(prior_obj$prior_data$total_variance[[i]]$prior)) next

      xx <- rowSums(exp(samps[,stan_data$which_hd]))
      title <- as.expression(bquote(V[.(paste0(
        get_node_name(prior_obj$node_data,
                      get_leaf_nodes(prior_obj$node_data,prior_obj$prior_data$total_variance[[i]]$id)), sep = "", collapse = ","))]))

      df_post <- rbind(df_post, data.frame(x = xx, param = as.character(title)))

    }

    # for HD priors
    for (i in seq_len(length(prior_obj$prior_data$weights))){

      if (prior_obj$prior_data$weights[[i]]$prior == "pc"){

        # calculate the correct variances over and under fraction bar
        ind2 <- stan_data$row_index_hd_pr_plot[i]
        over <- rowSums(as.matrix(exp(samps[,stan_data$w_o[ind2,] == 1])))
        under <- rowSums(as.matrix(exp(samps[,stan_data$w_u[ind2,] == 1])))

        t_under <- paste0(get_node_name(prior_obj$node_data, get_leaf_nodes(prior_obj$node_data, prior_obj$prior_data$weights[[i]]$id)), sep = "", collapse = ",")
        t_over <- as.character(
          paste0(get_node_name(prior_obj$node_data,
                               get_leaf_nodes(prior_obj$node_data, prior_obj$prior_data$weights[[i]]$param$above_node)),
          sep = "", collapse = ","))
        title <- as.expression(bquote(omega[frac(.(t_over), .(t_under))]))

        df_post <- rbind(df_post, data.frame(x = over/under, param = as.character(title)))

      } else { # dirichlet

        for (ind in 1:(prior_obj$prior_data$weights[[i]]$no_children-1)){

          ind2 <- stan_data$row_index_hd_pr_plot[i]+ind-1

          over <- rowSums(as.matrix(exp(samps[,stan_data$w_o[ind2,] == 1])))
          under <- rowSums(as.matrix(exp(samps[,stan_data$w_u[ind2,] == 1])))

          t_under <- paste0(get_node_name(prior_obj$node_data, get_leaf_nodes(prior_obj$node_data, prior_obj$prior_data$weights[[i]]$id)), sep = "", collapse = ",")
          t_over <- as.character(
            paste0(get_node_name(prior_obj$node_data,
                                 get_leaf_nodes(prior_obj$node_data, prior_obj$prior_data$weights[[i]]$children[ind]))),
            sep = "", collapse = ",")
          title <- as.expression(bquote(omega[frac(.(t_over), .(t_under))]))

          df_post <- rbind(df_post, data.frame(x = over/under, param = as.character(title)))

        }

      }

    }

  }

  df <- if (type == "prior") df_pri else if (type == "posterior") df_post else list(prior = df_pri, posterior = df_post)

  return(df)

}

df_posteriors_inla <- function(res, type = c("variance", "stdev", "precision"), axes = list()){

  type <- match.arg(type)

  if (type == "variance"){
    myfun <- function(x) 1/x
  } else if (type == "stdev") {
    myfun <- function(x) 1/sqrt(x)
  } else {
    myfun <- function(x) x
  }

  randeff_names <- names(res$prior$data$random)

  titles <- c()

  for (ind in seq_len(length(randeff_names))){
    titles[ind] <-
      if (type == "variance") {
        as.expression(bquote(sigma[.(randeff_names[ind])]^2))
      } else if (type == "stdev"){
        as.expression(bquote(sigma[.(randeff_names[ind])]))
      } else {
        as.expression(bquote(1/sigma[.(randeff_names[ind])]^2))
      }
  }

  df_post <- data.frame()

  for (ind in 1:length(res$inla$marginals.hyperpar)){

    thispar <- names(res$inla$marginals.hyperpar)[ind]
    thispar <- strsplit(thispar, "Precision for ")[[1]][2]
    if (res$prior$family == "gaussian" && !is.na(thispar) && thispar == "the Gaussian observations"){
      thispar <- "eps"
    }
    this_title <- titles[randeff_names == thispar]

    title <- if (type == "variance") {
      as.expression(bquote(sigma[.(thispar)]^2))
    } else if (type == "stdev"){
      as.expression(bquote(sigma[.(thispar)]))
    } else {
      as.expression(bquote(1/sigma[.(thispar)]^2))
    }

    tmp <- as.data.frame(INLA::inla.tmarginal(myfun, res$inla$marginals.hyperpar[[ind]]))

    if (length(which(names(axes) == thispar)) > 0){

      tmp_axes <- axes[[which(names(axes) == thispar)]]
      tmp2 <- tmp[tmp$x >= min(tmp_axes) & tmp$x <= max(tmp_axes),]
      if (nrow(tmp2) >= 2) tmp <- tmp2 # only change the axes if we have at least two points

    }

    tmp$param <- as.character(this_title)

    df_post <- rbind(
      df_post,
      tmp
    )

  }

  df_post$param <- factor(df_post$param, levels = as.character(titles))

  return(df_post)

  # gg1 <- ggplot(df_post, aes(x = x, y = y)) +
  #   geom_line() +
  #   facet_wrap(~param, labeller = label_parsed, scales = "free") +
  #   # theme(strip.text = element_text(size = 14, color = "white")) +
  #   ylab("Density") +
  #   theme(axis.title.x = element_blank()) +
  #   theme(panel.background = element_rect(fill ="white", color = gray(0.5)),
  #         panel.grid.major = element_line(color = gray(0.85)),
  #         panel.grid.minor = element_line(color = gray(0.92))) +
  #   theme(strip.background = element_rect(fill = "white", color = gray(0.5)),
  #         strip.text = element_text(color = "black", size = 15))
  #
  # return(gg1)

}

df_posterior_variances <- function(res, param = c("variance", "stdev", "precision")){

  samps <- rstan::extract(res$stan, "theta")$theta

  if (param == "variance"){
    titles <- sapply(as.character(res$prior$node_data$orig_nodedata$label), function(x) as.expression(bquote(sigma[.(x)]^2)))

    df <- data.frame(
      x = exp(unlist(c(samps))),
      param = factor(rep(as.character(titles), each = nrow(samps)), levels = titles)
    )
  } else if (param == "stdev"){
    titles <- sapply(as.character(res$prior$node_data$orig_nodedata$label), function(x) as.expression(bquote(sigma[.(x)])))

    df <- data.frame(
      x = exp(unlist(c(samps))/2),
      param = factor(rep(as.character(titles), each = nrow(samps)), levels = titles)
    )
  } else if (param == "precision"){
    titles <- sapply(as.character(res$prior$node_data$orig_nodedata$label), function(x) as.expression(bquote(1/sigma[.(x)]^2)))

    df <- data.frame(
      x = exp(-unlist(c(samps))),
      param = factor(rep(as.character(titles), each = nrow(samps)), levels = titles)
    )
  }

  return(df)

  # gg1 <- ggplot(df, aes(x = x, y = ..density..)) +
  #   geom_histogram(col = gray(0.5), fill = "#8E8D8A", bins = 40) +
  #   facet_wrap(~name, labeller = label_parsed, scales = "free") +
  #   # theme(strip.text = element_text(size = 14, color = "white")) +
  #   ylab("Density") +
  #   theme(axis.title.x = element_blank()) +
  #   theme(panel.background = element_rect(fill ="white", color = gray(0.5)),
  #         panel.grid.major = element_line(color = gray(0.85)),
  #         panel.grid.minor = element_line(color = gray(0.92))) +
  #   theme(strip.background = element_rect(fill = "white", color = gray(0.5)),
  #         strip.text = element_text(color = "black", size = 15))
  # # theme(strip.background = element_rect(fill = "#8E8D8A", color = gray(0.5)),
  # #       strip.text = element_text(color = plot_text_color, size = 15))
  #
  # return(gg1)

}



#' Plotting the prior tree structure graph
#'
#' Can only be used for visualization in R.
#' @param obj An object from \code{make_prior}, \code{inference_stan}, \code{inference_inla}, or \code{makemyprior_gui}
#' @param nodenames Custom names for each node (optional). Given as a named list with the default names as list names, and the
#' new names as list elements.
#' Do not need to provide all.
#' @keywords plot
#' @return A \link[visNetwork]{visNetwork} with the tree graph
#' See also \link[makemyprior]{makemyprior_plotting}.
#' @examples
#' ex_prior <- makemyprior_example_model()
#' plot_tree_structure(ex_prior)
#'
#' @export
plot_tree_structure <- function(obj, nodenames){

  if (is(obj, "mmp_inla") || is(obj, "mmp_stan")) obj <- obj$prior
  if (!is(obj, "mmp_prior")) stop("Provide an object with class 'mmp_prior', 'mmp_stan' or 'mmp_inla'.", call. = FALSE)

  nd <- list(x = obj$node_data)

  nd$x <- update_basemodel_edges(nd$x, obj$prior_data$weights)

  tmp_nodes <- nd$x$nodes
  tmp_nodes$color.background <- node_palette(calc_variance_proportions(nd$x)$nodes$varprop)
  tmp_nodes$color.background[tmp_nodes$status == "detached"] <- detached_node_color
  tmp_nodes$color.highlight.background <- tmp_nodes$color.background

  if (!missing(nodenames)){
    for (i in seq_len(length(nodenames))){
      tmp_nodes$label[which(tmp_nodes$label == names(nodenames)[i])] <- nodenames[[i]]
    }
  }

  tmp_edges <- nd$x$edges
  tmp_edges$width <- tmp_edges$width*5 # to get a bit thicker edges in the graph
  if (nrow(tmp_edges) > 0) tmp_edges <- cbind(tmp_edges, arrows = "to", hidden = FALSE)

  tmp_edges <- set_detached_edges(tmp_nodes, tmp_edges)
  tmp_edges$label <- ""

  nn <- visNetwork::visNetwork(tmp_nodes, tmp_edges, background = "white")
  nn <- visNetwork::visHierarchicalLayout(nn)
  nn <- visNetwork::visEdges(nn, color = list(highlight = node_highlight_border_color),
                             arrowStrikethrough = FALSE, selectionWidth = 1, smooth = FALSE,
                             font = list(color = node_border_color, size = 14, vadjust = 35,
                                         # background = "#EAE7DC", # background does not move with vadjust
                                         strokeWidth = 0, bold = TRUE)
  )
  nn <- visNetwork::visNodes(nn,
                             color = list(border = node_border_color,
                                          highlight = list(border = node_highlight_border_color)),
                             borderWidth = 1,
                             #borderWidthSelected = 1,
                             font = list(background = "white", color = node_border_color))
  nn <- visNetwork::visOptions(nn, highlightNearest = FALSE)
  nn <- visNetwork::visInteraction(nn, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, zoomView = FALSE)

  return(nn)

}



#' Plotting prior distributions
#'
#' Function for plotting the prior distributions of the random effects on the scale of the parameters chosen
#' @param obj An object from \code{make_prior}, \code{inference_stan}, \code{inference_inla}, or \code{makemyprior_gui}
#' @keywords plot
#' @return A \link[ggplot2]{ggplot} with the prior distributions.
#' See also \link[makemyprior]{makemyprior_plotting}.
#' @examples
#' ex_prior <- makemyprior_example_model()
#' plot_prior(ex_prior)
#'
#' @export
plot_prior <- function(obj){

  if (is(obj, "mmp_inla") || is(obj, "mmp_stan")){
    obj <- obj$prior
  } else if (!is(obj, "mmp_prior")) stop("Invalid input.", call. = FALSE)

  df <- make_dataframe_for_plotting("prior", obj)

  gg <- ggplot(df, aes(x = .data$x, y = .data$y)) + geom_line(na.rm = TRUE) +
    facet_wrap(~param, labeller = label_parsed, scales = "free") +
    # theme(strip.text = element_text(size = 14, color = "white")) +
    ylab("Density") +
    xlim(0, NA) +
    ylim(0, NA) +
    theme(axis.title.x = element_blank()) +
    theme(panel.background = element_rect(fill ="white", color = gray(0.5)),
          panel.grid.major = element_line(color = gray(0.85)),
          panel.grid.minor = element_line(color = gray(0.92))) +
    theme(strip.background = element_rect(fill = "white", color = gray(0.5)),
          strip.text = element_text(color = "black", size = 15))
    # theme(strip.background = element_rect(fill = "#8E8D8A", color = gray(0.5)),
    #       strip.text = element_text(color = plot_text_color, size = 15))

  return(gg)

}


#' Plotting posterior distributions
#'
#' Function for plotting the posterior distributions of the random effect variances
#' on the scale of the tree parameterization.
#' @param obj An object from \code{inference_stan}
#' @param param A string indicating parameterization of plot.
#' \code{"prior"} for scale of parameters,
#' \code{"variance"}, \code{"stdev"} and \code{"precision"} also possible.
#' @param prior Include prior in the plot? Only possible for \code{param = "prior"}.
#' Note that if Jeffreys' prior is used for the total variance, it will not be included in the plot.
#' @keywords plot
#' @return A \link[ggplot2]{ggplot} with the posterior distributions.
#' See also \link[makemyprior]{makemyprior_plotting}.
#' @examples
#'
#'
#' if (interactive() && requireNamespace("rstan")){
#'   ex_prior <- makemyprior_example_model()
#'   res_stan <- inference_stan(ex_prior, iter = 100)
#'   # Note: For reliable results, increase the number of iterations (e.g., 'iter = 2000')
#'   plot_posterior_stan(res_stan)
#' }
#'
#' @export
plot_posterior_stan <- function(obj, param = c("prior", "variance", "stdev", "precision"), prior = FALSE){

  if (is(obj, "mmp_inla")) stop("You cannot use a posterior fitted with inla for this function. Use 'plot_posterior_variance/stdev/precision' instead.", call. = FALSE)
  if (!is(obj, "mmp_stan")) stop("Invalid input.", call. = FALSE)
  param <- match.arg(param)
  if (param != "prior") prior <- FALSE

  if (param == "prior"){
    if (prior){
      df <- make_dataframe_for_plotting("both", obj)
      gg <- ggplot() +
        geom_histogram(data = df$posterior, mapping = aes(x = .data$x, y = .data$..density..), col = gray(0.5), fill = "#8E8D8A", bins = 40) +
        geom_line(data = df$pri, mapping = aes(x = .data$x, y = .data$y), na.rm = TRUE) +
        facet_wrap(~param, labeller = label_parsed, scales = "free") +
        # theme(strip.text = element_text(size = 14, color = "white")) +
        ylab("Density") +
        theme(axis.title.x = element_blank()) +
        theme(panel.background = element_rect(fill ="white", color = gray(0.5)),
              panel.grid.major = element_line(color = gray(0.85)),
              panel.grid.minor = element_line(color = gray(0.92))) +
        theme(strip.background = element_rect(fill = "white", color = gray(0.5)),
              strip.text = element_text(color = "black", size = 15))
    } else {
      df <- make_dataframe_for_plotting("posterior", obj)
      gg <- ggplot(df, aes(x = .data$x, y = .data$..density..)) + geom_histogram(col = gray(0.5), fill = "#8E8D8A", bins = 40) +
        facet_wrap(~param, labeller = label_parsed, scales = "free") +
        # theme(strip.text = element_text(size = 14, color = "white")) +
        ylab("Density") +
        theme(axis.title.x = element_blank()) +
        theme(panel.background = element_rect(fill ="white", color = gray(0.5)),
              panel.grid.major = element_line(color = gray(0.85)),
              panel.grid.minor = element_line(color = gray(0.92))) +
        theme(strip.background = element_rect(fill = "white", color = gray(0.5)),
              strip.text = element_text(color = "black", size = 15))
    }
  } else {
    if (param == "variance"){
      gg <- plot_posterior_variance(obj)
    } else if (param == "stdev"){
      gg <- plot_posterior_stdev(obj)
    } else {
      gg <- plot_posterior_precision(obj)
    }
  }

  return(gg)

}


#' Plotting posterior distributions
#'
#' Function for plotting the posterior distributions of the coefficients of the fixed effects
#' @param obj An object from \code{inference_stan} or \code{inference_inla}
#' @keywords plot
#' @return A \link[ggplot2]{ggplot} with the posterior distributions.
#' See also \link[makemyprior]{makemyprior_plotting}.
#' @examples
#'
#' if (interactive() && requireNamespace("rstan")){
#'   ex_prior <- makemyprior_example_model()
#'   res_stan <- inference_stan(ex_prior, iter = 100)
#'   # Note: For reliable results, increase the number of iterations (e.g., 'iter = 2000')
#'   plot_posterior_fixed(res_stan)
#' }
#'
#' @export
plot_posterior_fixed <- function(obj){

  df_post <- data.frame()

  if (is(obj, "mmp_stan")){

    samps <- rstan::extract(obj$stan, c("intercept", "coeff"))

    if (obj$prior$use_intercept){
      title <- expression(mu)
      # xx <- seq(min(samps$intercept[,1]), max(samps$intercept[,1]), length.out = 200)
      df_post <- rbind(df_post, data.frame(x = samps$intercept[,1], param = as.character(title)))
    }

    for (i in seq_len(length(obj$prior$data$fixed))){
      title <- as.expression(bquote(beta[.(names(obj$prior$data$fixed)[i])]))
      # xx <- seq(min(samps$coeff[,i]), max(samps$coeff[,i]), length.out = 200)
      df_post <- rbind(df_post, data.frame(x = samps$coeff[,i], param = as.character(title)))
    }

    gg <- ggplot() +
      geom_histogram(data = df_post, mapping = aes(x = .data$x, y = .data$..density..), col = gray(0.5), fill = "#8E8D8A", bins = 40) +
      facet_wrap(~param, labeller = label_parsed, scales = "free") +
      # theme(strip.text = element_text(size = 14, color = "white")) +
      ylab("Density") +
      theme(axis.title.x = element_blank()) +
      theme(panel.background = element_rect(fill ="white", color = gray(0.5)),
            panel.grid.major = element_line(color = gray(0.85)),
            panel.grid.minor = element_line(color = gray(0.92))) +
      theme(strip.background = element_rect(fill = "white", color = gray(0.5)),
            strip.text = element_text(color = "black", size = 15))
    # theme(strip.background = element_rect(fill = "#8E8D8A", color = gray(0.5)),
    #       strip.text = element_text(color = plot_text_color, size = 15))

  } else if (is(obj, "mmp_inla")){

    if (obj$prior$use_intercept){
      df_post <- rbind(df_post,
                       data.frame(x = obj$inla$marginals.fixed$`(Intercept)`[,1],
                                  y = obj$inla$marginals.fixed$`(Intercept)`[,2],
                                  param = as.character(expression(mu))))
    }

    for (i in seq_len(length(obj$prior$data$fixed))){
      title <- as.expression(bquote(beta[.(names(obj$prior$data$fixed)[i])]))
      i2 <- which(names(obj$inla$marginals.fixed) == names(obj$prior$data$fixed)[i]) # which effect
      df_post <- rbind(df_post,
                       data.frame(x = obj$inla$marginals.fixed[[i2]][,1],
                                  y = obj$inla$marginals.fixed[[i2]][,2],
                                  param = as.character(title)))
    }

    gg <- ggplot() +
      geom_line(df_post, mapping = aes(x = .data$x, y = .data$y)) +
      facet_wrap(~param, labeller = label_parsed, scales = "free") +
      # theme(strip.text = element_text(size = 14, color = "white")) +
      ylab("Density") +
      theme(axis.title.x = element_blank()) +
      theme(panel.background = element_rect(fill ="white", color = gray(0.5)),
            panel.grid.major = element_line(color = gray(0.85)),
            panel.grid.minor = element_line(color = gray(0.92))) +
      theme(strip.background = element_rect(fill = "white", color = gray(0.5)),
            strip.text = element_text(color = "black", size = 15))


  } else stop("Invalid input.", call. = FALSE)

  return(gg)

}


#' Plotting posterior variances, standard deviations or precisions
#'
#' @param obj An object from \code{inference_stan} or \code{inference_inla}.
#' @keywords plot
#' @return A \link[ggplot2]{ggplot} object with the plot
#' See also \link[makemyprior]{makemyprior_plotting}.
#' @examples
#'
#' if (interactive() && requireNamespace("rstan")){
#'   ex_prior <- makemyprior_example_model()
#'   res_stan <- inference_stan(ex_prior, iter = 100)
#'   # Note: For reliable results, increase the number of iterations (e.g., 'iter = 2000')
#'   plot_posterior_variance(res_stan)
#' }
#'
#' if (interactive() && requireNamespace("INLA")){
#'   ex_prior <- makemyprior_example_model()
#'   res_inla <- inference_inla(ex_prior)
#'   plot_posterior_variance(res_inla)
#' }
#'
#' @export
plot_posterior_variance <- function(obj){

  if (is(obj, "mmp_stan")){

    df <- df_posterior_variances(obj, "variance")
    gg <- ggplot(df, aes(x = .data$x, y = .data$..density..)) + geom_histogram(col = gray(0.5), fill = "#8E8D8A", bins = 40) +
      facet_wrap(~param, labeller = label_parsed, scales = "free") +
      # theme(strip.text = element_text(size = 14, color = "white")) +
      ylab("Density") +
      theme(axis.title.x = element_blank()) +
      theme(panel.background = element_rect(fill ="white", color = gray(0.5)),
            panel.grid.major = element_line(color = gray(0.85)),
            panel.grid.minor = element_line(color = gray(0.92))) +
      theme(strip.background = element_rect(fill = "white", color = gray(0.5)),
            strip.text = element_text(color = "black", size = 15))

  } else if (is(obj, "mmp_inla")){

    df <- df_posteriors_inla(obj, "variance")
    gg <- ggplot() +
      geom_line(df, mapping = aes(x = .data$x, y = .data$y)) +
      facet_wrap(~param, labeller = label_parsed, scales = "free") +
      # theme(strip.text = element_text(size = 14, color = "white")) +
      ylab("Density") +
      theme(axis.title.x = element_blank()) +
      theme(panel.background = element_rect(fill ="white", color = gray(0.5)),
            panel.grid.major = element_line(color = gray(0.85)),
            panel.grid.minor = element_line(color = gray(0.92))) +
      theme(strip.background = element_rect(fill = "white", color = gray(0.5)),
            strip.text = element_text(color = "black", size = 15))

  } else stop("Invalid input.", call. = FALSE)

  return(gg)

}

#' @rdname plot_posterior_variance
#' @export
plot_posterior_stdev <- function(obj){

  if (is(obj, "mmp_stan")){

    df <- df_posterior_variances(obj, "stdev")
    gg <- ggplot(df, aes(x = .data$x, y = .data$..density..)) + geom_histogram(col = gray(0.5), fill = "#8E8D8A", bins = 40) +
      facet_wrap(~param, labeller = label_parsed, scales = "free") +
      # theme(strip.text = element_text(size = 14, color = "white")) +
      ylab("Density") +
      theme(axis.title.x = element_blank()) +
      theme(panel.background = element_rect(fill ="white", color = gray(0.5)),
            panel.grid.major = element_line(color = gray(0.85)),
            panel.grid.minor = element_line(color = gray(0.92))) +
      theme(strip.background = element_rect(fill = "white", color = gray(0.5)),
            strip.text = element_text(color = "black", size = 15))

  } else if (is(obj, "mmp_inla")){

    df <- df_posteriors_inla(obj, "stdev")
    gg <- ggplot() +
      geom_line(df, mapping = aes(x = .data$x, y = .data$y)) +
      facet_wrap(~param, labeller = label_parsed, scales = "free") +
      # theme(strip.text = element_text(size = 14, color = "white")) +
      ylab("Density") +
      theme(axis.title.x = element_blank()) +
      theme(panel.background = element_rect(fill ="white", color = gray(0.5)),
            panel.grid.major = element_line(color = gray(0.85)),
            panel.grid.minor = element_line(color = gray(0.92))) +
      theme(strip.background = element_rect(fill = "white", color = gray(0.5)),
            strip.text = element_text(color = "black", size = 15))

  } else stop("Invalid input.", call. = FALSE)

  return(gg)

}

#' @rdname plot_posterior_variance
#' @export
plot_posterior_precision <- function(obj){

  if (is(obj, "mmp_stan")){

    df <- df_posterior_variances(obj, "precision")
    gg <- ggplot(df, aes(x = .data$x, y = .data$..density..)) + geom_histogram(col = gray(0.5), fill = "#8E8D8A", bins = 40) +
      facet_wrap(~param, labeller = label_parsed, scales = "free") +
      ylab("Density") +
      theme(axis.title.x = element_blank()) +
      theme(panel.background = element_rect(fill ="white", color = gray(0.5)),
            panel.grid.major = element_line(color = gray(0.85)),
            panel.grid.minor = element_line(color = gray(0.92))) +
      theme(strip.background = element_rect(fill = "white", color = gray(0.5)),
            strip.text = element_text(color = "black", size = 15))

  } else if (is(obj, "mmp_inla")){

    df <- df_posteriors_inla(obj, "precision")
    gg <- ggplot() +
      geom_line(df, mapping = aes(x = .data$x, y = .data$y)) +
      facet_wrap(~param, labeller = label_parsed, scales = "free") +
      # theme(strip.text = element_text(size = 14, color = "white")) +
      ylab("Density") +
      theme(axis.title.x = element_blank()) +
      theme(panel.background = element_rect(fill ="white", color = gray(0.5)),
            panel.grid.major = element_line(color = gray(0.85)),
            panel.grid.minor = element_line(color = gray(0.92))) +
      theme(strip.background = element_rect(fill = "white", color = gray(0.5)),
            strip.text = element_text(color = "black", size = 15))

  } else stop("Invalid input.", call. = FALSE)

  return(gg)

}


#' Extract the posterior of a random effect
#'
#' Extract the posterior of a random effect in the model for inference done with Stan
#' @param obj An object from \code{inference_stan}.
#' @param effname Name of the random effect, same name as in the data.
#' @keywords posterior
#' @return Returns a matrix with the posterior samples of the chosen effect
#' @examples
#'
#' if (interactive() && requireNamespace("rstan")){
#'   ex_prior <- makemyprior_example_model()
#'   res_stan <- inference_stan(ex_prior, iter = 100)
#'   # Note: For reliable results, increase the number of iterations (e.g., 'iter = 2000')
#'   extract_posterior_effect(res_stan, "a")
#' }
#'
#' @export
extract_posterior_effect <- function(obj, effname){

  if (!is(obj, "mmp_stan")) stop("Input 'obj' must be of class 'mmp_stan'.", call. = FALSE)

  if (!effname %in% names(obj$prior$data$random)) stop(paste0("The random effect ", effname, " is not in the model."), call. = FALSE)

  samps <- rstan::extract(obj$stan, "effects")$effects

  # find which indexes this effect belongs to
  which_eff <- which(names(obj$prior$data$random) == effname)

  return(
    samps[,(obj$stan_data$effect_index_start[which_eff]) + 1:obj$stan_data$effect_sizes[which_eff]]
  )

}



#' Extract the posterior parameter estimate
#'
#' Extract the posterior parameter estimate of a random effect variance or fixed effect coefficient
#' when inference is done with Stan
#'
#' @param obj An object from \code{inference_stan}.
#' @param param Name of the variance parameter, which is the same as the name of the corresponding
#' fixed or random effect in the data. Intercept is denoted 'intercept', and residual variance is denoted 'eps'.
#' @keywords posterior
#' @return Returns a vector with the posterior samples of the chosen parameter, on variance scale for
#' variances parameters and original (the common) scale for fixed effect coefficients
#' @examples
#'
#' if (interactive() && requireNamespace("rstan")){
#'   ex_prior <- makemyprior_example_model()
#'   res_stan <- inference_stan(ex_prior, iter = 100)
#'   # Note: For reliable results, increase the number of iterations (e.g., 'iter = 2000')
#'   extract_posterior_parameter(res_stan, "intercept")
#'   extract_posterior_parameter(res_stan, "a")
#' }
#'
#' @export
extract_posterior_parameter <- function(obj, param){

  if (!is(obj, "mmp_stan")) stop("Input 'obj' must be of class 'mmp_stan'.", call. = FALSE)

  if (param == "intercept"){
    if (obj$prior$use_intercept){
      return(c(rstan::extract(obj$stan, "intercept")$intercept))
    } else {
      stop("This model has no intercept.", call. = FALSE)
    }
  } else if (param %in% names(obj$prior$data$fixed)){
    samps <- rstan::extract(obj$stan, "coeff")$coeff
    which_eff <- which(names(obj$prior$data$fixed) == param)
    return(samps[,which_eff])
  } else if (param %in% names(obj$prior$data$random)){
    samps <- rstan::extract(obj$stan, "theta")$theta
    which_eff <- which(sapply(obj$prior$prior_data$cw_priors, function(x) x$name) == param)
    return(exp(samps[,which_eff]))
  } else {
    stop(paste0("The effect ", param, " is not in the model."), call. = FALSE)
  }

  stop("Could not return the desired parameter, see if the name is correct.", call. = FALSE)

}


#' Plotting several posterior distributions
#'
#' Function for plotting the posterior distributions of the random effect variances
#' on the scale of the tree parameterization.
#' @param objs A names list with objects of class \code{mmp_stan} from \code{inference_stan}, can be any length
#' (but typically length two for one prior (\code{use_likelihood = FALSE}) and posterior, or two posteriors).
#' @param param A string indicating parameterization of plot.
#' \code{"prior"} for scale of parameters,
#' \code{"variance"}, \code{"stdev"}, \code{"precision"} and \code{"logvariance"} are also possible.
#' @keywords plot
#' @return A \link[ggplot2]{ggplot} with the posterior distributions.
#' @details We cannot sample from a Jeffreys' prior since it is improper.
#' If Jeffreys' prior is used for the total variance, the prior will be changed to a Gaussian(0,1) prior on
#' the log total variance. This means that it does not make sense to look at the variances/standard deviations/precisions,
#' but the variance proportions will be correct.
#' See also \link[makemyprior]{makemyprior_plotting}.
#' @examples
#'
#' if (interactive() && requireNamespace("rstan")){
#'   ex_prior1 <- makemyprior_example_model(seed = 1)
#'   ex_prior2 <- makemyprior_example_model(seed = 2)
#'   # Note: For reliable results, increase the number of iterations (e.g., 'iter = 2000')
#'   res_stan1 <- inference_stan(ex_prior1, iter = 100)
#'   res_stan2 <- inference_stan(ex_prior2, iter = 100)
#'   plot_several_posterior_stan(list(One = res_stan1, Two = res_stan2))
#' }
#'
#' @export
plot_several_posterior_stan <- function(objs, param = c("prior", "variance", "stdev", "precision", "logvariance")){

  if (!is(objs, "list") || any(is.null(names(objs)))) stop("Provide a named list with objects from 'inference_stan'.", call. = FALSE)
  if (!all(sapply(objs, is, class2 = "mmp_stan"))) stop("The list 'objs' must consist of objects from 'inference_stan'.", call. = FALSE)
  param <- match.arg(param)

  df <- data.frame()
  if (param == "prior"){
    for (ind in 1:length(objs)){
      tmp <- make_dataframe_for_plotting("posterior", objs[[ind]])
      tmp$param2 <- names(objs)[ind]
      df <- rbind(df, tmp)
    }
  } else {
    for (ind in 1:length(objs)){
      if (param == "logvariance") {
        tmp <- df_posterior_variances(objs[[ind]], "variance")
        tmp$x <- log(tmp$x)
      } else {
        tmp <- df_posterior_variances(objs[[ind]], param)
      }
      tmp$param2 <- names(objs)[ind]
      df <- rbind(df, tmp)
    }
  }

  df$param2 <- factor(df$param2, levels = names(objs))

  gg <- ggplot() +
    geom_histogram(data = df, mapping = aes(x = .data$x, y = .data$..density.., fill = .data$param2),
                   col = gray(0.5), bins = 40, alpha = 1/(length(objs)+1), position = position_identity()) +
    facet_wrap(~ param, labeller = label_parsed, scales = "free") +
    scale_fill_viridis_d(option = "C", end = 0.7) +
    ylab("Density") +
    theme(axis.title.x = element_blank()) +
    theme(legend.title = element_blank(), legend.position = "bottom") +
    theme(panel.background = element_rect(fill ="white", color = gray(0.5)),
          panel.grid.major = element_line(color = gray(0.85)),
          panel.grid.minor = element_line(color = gray(0.92))) +
    theme(strip.background = element_rect(fill = "white", color = gray(0.5)),
          strip.text = element_text(color = "black", size = 15))

  return(gg)

}





# extracting some posterior values from stan object
make_posterior_summary_stan <- function(res){

  if (!is(res, "mmp_stan")) stop("An object of class 'mmp_stan' must be provided.", call. = FALSE)

  prior_obj <- res$prior

  df <- data.frame(par = c(), mean = c(), median = c(), sd = c())

  samps <- rstan::extract(res$stan, "theta")$theta

  # for CW priors:
  for (i in 1:length(prior_obj$prior_data$cw_priors)){

    if (prior_obj$prior_data$cw_priors[[i]]$prior != ""){

      xx <- exp(samps[,i])
      title <- paste0("sigma^2", "[", prior_obj$prior_data$cw_priors[[i]]$name, "]")

      df <- rbind(df, data.frame(par = title, mean = mean(xx), median = median(xx), sd = sd(xx)))

    }

  }

  # for total variances
  for (i in seq_len(length(prior_obj$prior_data$total_variance))){

    if (is.null(prior_obj$prior_data$total_variance[[i]]$prior)) next

    xx <- rowSums(exp(samps[,res$stan_data$which_hd]))

    title <- paste0(c("V", "[",
                    prior_obj$prior_data$total_variance[[i]]$name,
                    "]"), collapse = "")

    df <- rbind(df, data.frame(par = title, mean = mean(xx), median = median(xx), sd = sd(xx)))

  }

  # for HD priors
  for (i in seq_len(length(prior_obj$prior_data$weights))){

    if (prior_obj$prior_data$weights[[i]]$prior == "pc"){

      # calculate the correct variances over and under fraction bar
      ind2 <- res$stan_data$row_index_hd_pr_plot[i]
      over <- rowSums(as.matrix(exp(samps[,res$stan_data$w_o[ind2,] == 1])))
      under <- rowSums(as.matrix(exp(samps[,res$stan_data$w_u[ind2,] == 1])))

      t_under <- as.character(prior_obj$prior_data$weights[[i]]$name)
      t_over <- as.character(get_node_name(prior_obj$node_data, prior_obj$prior_data$weights[[i]]$param$above_node))
      title <- sprintf("w[%s/%s]", t_over, t_under)

      df <- rbind(df, data.frame(par = title, mean = mean(over/under), median = median(over/under), sd = sd(over/under)))

    } else { # dirichlet

      for (ind in 1:(prior_obj$prior_data$weights[[i]]$no_children-1)){

        ind2 <- res$stan_data$row_index_hd_pr_plot[i]+ind-1

        over <- rowSums(as.matrix(exp(samps[,res$stan_data$w_o[ind2,] == 1])))
        under <- rowSums(as.matrix(exp(samps[,res$stan_data$w_u[ind2,] == 1])))

        t_under <- paste0(get_node_name(prior_obj$node_data, get_leaf_nodes(prior_obj$node_data, prior_obj$prior_data$weights[[i]]$id)), sep = "", collapse = ",")
        t_over <- as.character(get_node_name(prior_obj$node_data, prior_obj$prior_data$weights[[i]]$children[ind]))
        title <- sprintf("w[%s/%s]", t_over, t_under)

        df <- rbind(df, data.frame(par = title, mean = mean(over/under), median = median(over/under), sd = sd(over/under)))

      }

    }

  }


  # intercept
  if (prior_obj$use_intercept){
    samps <- rstan::extract(res$stan, "intercept")$intercept

    df <- rbind(df, data.frame(par = "intercept", mean = mean(samps), median = median(samps), sd = sd(samps)))

  }

  # fixed effects
  if (length(prior_obj$data$fixed) > 0){
    samps <- rstan::extract(res$stan, "coeff")$coeff
  }
  for (ind in seq_along(prior_obj$data$fixed)){

    title <- row.names(prior_obj$prior_fixed)[ind]

    df <- rbind(df, data.frame(par = title, mean = mean(samps[,ind]), median = median(samps[,ind]), sd = sd(samps[,ind])))

  }

  names(df)[1] <- "Param."

  return(df)

}


# extracting some posterior values from inla object
make_posterior_summary_inla <- function(res){

  if (!is(res, "mmp_inla")) stop("An object of class 'mmp_inla' must be provided.", call. = FALSE)

  prior_obj <- res$prior

  tmp <- res$inla$summary.hyperpar[, c("mean", "0.5quant", "sd")]
  names(tmp) <- c("mean", "median", "sd")
  for (ind in seq_len(nrow(tmp))){
    if (row.names(tmp)[ind] == "Precision for the Gaussian observations"){
      row.names(tmp)[ind] <- "1/sigma[eps]^2"
    } else {
      row.names(tmp)[ind] <- paste0("1/sigma[", strsplit(row.names(tmp)[ind], "Precision for ")[[1]][2], "]^2", collapse = "")
    }
  }

  if (nrow(res$inla$summary.fixed) > 0){
    tmp2 <- res$inla$summary.fixed
    row.names(tmp2)[row.names(tmp2) == "(Intercept)"] <- "intercept"
    tmp2 <- tmp2[, c("mean", "0.5quant", "sd")]
    names(tmp2) <- c("mean", "median", "sd")
  } else tmp2 <- data.frame()

  tmp3 <- rbind(tmp, tmp2)

  tmp3 <- cbind(data.frame(par = row.names(tmp3)), tmp3)
  names(tmp3)[1] <- "Param."
  row.names(tmp3) <- 1:nrow(tmp3)

  return(tmp3)

}


#' Plotting prior for a single parameter (weight or variance (not standard deviation))
#'
#' Following the parameterization of the prior.
#' @param x Values to evaluate prior in.
#' @param obj An object from \code{make_prior}, \code{inference_stan}, \code{inference_inla}, or \code{makemyprior_gui}
#' @param param Name of parameter to plot (see \code{print(obj)} for syntax). Note that only variances
#' will be plotted, so \code{V[..]} and \code{sigma^2[..]} must be used to indicate those parameters.
#' @param sd Whether to plot variance parameters on the standard deviation (\code{TRUE}) or variance (\code{FALSE}, default) scale
#' @keywords plot
#' @return A \link[ggplot2]{ggplot} with the posterior distribution.
#' See also \link[makemyprior]{makemyprior_plotting}.
#' @examples
#' ex_prior <- makemyprior_example_model()
#' plot_marginal_prior(seq(0, 1, 0.001), ex_prior, "w[a/a_b]")
#' plot_marginal_prior(seq(0, 1, 0.001), ex_prior, "w[eps/eps_a_b]")
#' plot_marginal_prior(seq(0, 5, 0.01), ex_prior, "V[eps_a_b]")
#'
#' @export
plot_marginal_prior <- function(x, obj, param, sd = FALSE){

  if (is(obj, "mmp_inla") || is(obj, "mmp_stan")){
    obj <- obj$prior
  } else if (!is(obj, "mmp_prior")) stop("Invalid input 'obj'.", call. = FALSE)

  df <- make_dataframe_for_plotting("prior", obj)
  df <- data.frame()

  # make latex-code from input name for matching with data.frame

  if (substr(param, 1, 4) == "sqrt") stop("You must use 'V' to plot total variance, not 'sqrt(V)'.", call. = FALSE)

  if (substr(param, 1, 2) == "w["){

    under <- substr(strsplit(param, "/")[[1]][2], 1, nchar(strsplit(param, "/")[[1]][2])-1)
    ind <- which(sapply(obj$prior_data$weights, function(x) x$name) == under)
    if (length(ind) != 1) stop("Not valid input 'param'. Print prior object for parameter names and format.", call. = FALSE)
    prior_info <- obj$prior_data$weights[[ind]]

    over <- substr(strsplit(param, "/")[[1]][1], 3, nchar(strsplit(param, "/")[[1]][1]))
    under_new <- get_node_strings(get_leaf_nodes(obj$node_data, prior_info$id), obj$node_data)
    new_param <- paste0("omega[frac(\"", over, "\", \"", under_new, "\")]")
    myxlab <- expression(omega)

    df <- data.frame(x = x, y = eval_spline_prior(x, obj$weight_priors[[ind]], FALSE))

  } else if (substr(param, 1, 2) == "V["){

    ind <- which(sapply(obj$prior_data$total_variance, function(x) x$name) == substr(param, 3, nchar(param)-1))
    if (length(ind) != 1) stop("Not valid input 'param'. Print prior object for parameter names and format.", call. = FALSE)
    prior_info <- obj$prior_data$total_variance[[ind]]
    inside <- get_node_strings(get_leaf_nodes(obj$node_data, prior_info$id), obj$node_data)
    new_param <- paste0("V[\"", inside, "\"]")
    myxlab <- if (sd) expression(sqrt(V)) else expression(V)

    df <- df_prior_variance_custom(x, prior_info, sd)

  } else if (substr(param, 1, 5) == "sigma"){

    if (substr(param, 6, 6) == "[") stop("You must use 'sigma^2' to plot variance, not 'sigma'.", call. = FALSE)

    ind <- which(sapply(obj$prior_data$cw_priors, function(x) x$name) == substr(param, 9, nchar(param)-1))
    if (length(ind) != 1) stop("Not valid input 'param'. Print prior object for parameter names and format.", call. = FALSE)
    prior_info <- obj$prior_data$cw_priors[[ind]]
    inside <- get_node_strings(get_leaf_nodes(obj$node_data, prior_info$id), obj$node_data)
    new_param <- paste0("sigma[\"", inside, "\"]^2")
    myxlab <- if (sd) expression(sigma) else expression(sigma^2)

    df <- df_prior_variance_custom(x, prior_info, sd)

  } else stop("Not valid input 'param'. Print prior object for parameter names and format.", call. = FALSE)

  if (nrow(df) == 0) stop("Not valid input 'param'. Print prior object for parameter names and format.", call. = FALSE)

  gg <- ggplot(df, aes(x = .data$x, y = .data$y)) + geom_line(na.rm = TRUE) +
    ylab("Density") +
    xlim(0, NA) +
    ylim(0, NA) +
    xlab(myxlab) +
    theme(panel.background = element_rect(fill ="white", color = gray(0.5)),
          panel.grid.major = element_line(color = gray(0.85)),
          panel.grid.minor = element_line(color = gray(0.92)))

  return(gg)

}


#' Evaluate PC prior for variance proportion
#'
#'Evaluate PC prior for a variance proportion.
#' @param x Values to evaluate prior in.
#' @param obj Prior object from \link[makemyprior]{make_prior}
#' @param param Which weight to plot, indicated using syntax shown when printing (do not need to include the
#' \code{w[..]} part to indicate that it is a variance proportion, but can be included). Print the prior object
#' to see syntax for each weight.
#' @param logitscale Is the input \code{x} on logit-scale? (default \code{FALSE}).
#' @return Returns density for the given variance proportion.
#' @examples
#' ex_prior <- makemyprior_example_model()
#' eval_pc_prior(seq(0, 1, 0.01), ex_prior, "eps/eps_a_b")
#' # or:
#' eval_pc_prior(seq(0, 1, 0.01), ex_prior, "w[eps/eps_a_b]")
#'
#' @export
eval_pc_prior <- function(x, obj, param, logitscale = FALSE){

  if (!is(obj, "mmp_prior")) stop("Invalid input 'obj'.", call. = FALSE)

  if (substr(param, 1, 2) == "w[") param <- substr(param, 3, nchar(param)-1)

  comps <- strsplit(param, "/")[[1]]

  ind <- which(sapply(obj$prior_data$weights, function(x) x$name) == comps[2])
  if (length(ind) != 1) stop("Not valid input 'param'. Print prior object for parameter names and format.", call. = FALSE)
  if (obj$prior_data$weights[[ind]]$prior != "pc") stop("Not a PC prior.", call. = FALSE)

  return(eval_spline_prior(x, obj$weight_priors[[ind]], logitscale = logitscale))

}







