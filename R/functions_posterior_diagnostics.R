


# plotting the prior and posterior of each effect

# some of these are also used in the shiny-app

# for one component, either CW or total variance
# all priors are plotted on variance scale
plot_prior_variance <- function(prior_data, node_data, var_type){

  if (is.null(prior_data$prior)) return(data.frame())

  if (prior_data$prior == "pc0"){
    xx <- c(seq(1, 20, 0.1))
    yy <- dexp(sqrt(xx), rate = -log(prior_data$param[2])/prior_data$param[1]) * 1/(2*sqrt(xx))
  } else if (prior_data$prior == "invgam"){
    xx <- seq(0, 2e5, 400)
    yy <- LaplacesDemon::dinvgamma(xx, shape = prior_data$param[1], scale = 1/prior_data$param[2])
  } else if (prior_data$prior == "hc"){
    xx <- sqrt(seq(0, 10000, 10))
    yy <- LaplacesDemon::dhalfcauchy(sqrt(xx), scale = prior_data$param[1]) * 1/sqrt(xx)
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
  over <- paste0(as.character(get_node_name(node_data, get_leaf_nodes(node_data, prior_data$param$basemodel_node))), sep = "", collapse = ",")
  title <- as.expression(bquote(omega[frac(.(over), .(under))]))

  xx <- seq(0, 1, 0.01)
  yy <- eval_spline_prior(xx, weight_data, F)

  df <- data.frame(x = xx, y = yy, param = as.character(title))

  return(df)

}



# making data-frame for plotting prior or posteriors, or both for Stan
make_dataframe_for_plotting <- function(type = c("prior", "posterior", "both"), res){

  if (length(res) == 3){
    prior_obj <- res$prior
    res_stan <- res$stan
    stan_data <- res$stan_data
  } else {
    prior_obj <- res
    type <- "prior" # can only plot prior if that is the only object that is provided
  }

  type <- match.arg(type)

  if (type == "prior" || type == "both"){

    df_pri <- data.frame()

    # for CW priors:
    for (i in 1:length(prior_obj$prior_data$cw_priors)){
      df_pri <- rbind(df_pri, plot_prior_variance(prior_obj$prior_data$cw_priors[[i]], prior_obj$node_data, "cw"))
    }

    # for total variances
    for (i in seq_len(length(prior_obj$prior_data$total_variance))){
      df_pri <- rbind(df_pri, plot_prior_variance(prior_obj$prior_data$total_variance[[i]], prior_obj$node_data, "total"))
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
                               get_leaf_nodes(prior_obj$node_data, prior_obj$prior_data$weights[[i]]$param$basemodel_node)),
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

#' Plotting prior distributions
#'
#' Function for plotting the prior distributions of the random effects on the scale of the parameters chosen
#' @param res An object from ``create_prior``, ``inference_stan``, or from ``run_shiny``
#' @keywords plot
#' @return A gg-plot with the prior distributions.
#' @examples
#' vignette("plotting", package = "priorconstruction")

# plotting all priors of random effects in the model
plot_priors_ranef <- function(res){

  if (length(res) == 3) res <- res$prior
  df <- make_dataframe_for_plotting("prior", res)

  gg <- ggplot(df, aes(x = x, y = y)) + geom_line(na.rm = TRUE) +
    facet_wrap(~param, labeller = label_parsed, scales = "free") +
    # theme(strip.text = element_text(size = 14, color = "white")) +
    ylab("Density") +
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
#' Function for plotting the posterior distributions of the random effects on the scale of the parameters chosen
#' @param res An object from ``inference_stan``, or from ``run_shiny`` where the "Close and run" option has been used
#' @keywords plot
#' @return A gg-plot with the posterior distributions.
#' @examples
#' vignette("plotting", package = "priorconstruction")

# plotting all posteriors of random effect in the model from stan
plot_posteriors_ranef <- function(res){

  df <- make_dataframe_for_plotting("posterior", res)

  gg <- ggplot(df, aes(x = x, y = ..density..)) + geom_histogram(col = gray(0.5), fill = "#8E8D8A", bins = 40) +
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

  return(gg)

}


#' Plotting prior and posterior distributions
#'
#' Function for plotting the prior and posterior distributions of the random effects together on the scale of the parameters chosen
#' @param res An object from ``inference_stan``, or from ``run_shiny`` where the "Close and run" option has been used
#' @keywords plot
#' @return A gg-plot with the prior and posterior distributions.
#' @examples
#' vignette("plotting", package = "priorconstruction")

plot_priors_posteriors_ranef <- function(res){

  df <- make_dataframe_for_plotting("both", res)

  gg <- ggplot() +
    geom_histogram(data = df$posterior, mapping = aes(x = x, y = ..density..), col = gray(0.5), fill = "#8E8D8A", bins = 40) +
    geom_line(data = df$pri, mapping = aes(x = x, y = y), na.rm = TRUE) +
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

  return(gg)

}


#' Plotting posterior distributions
#'
#' Function for plotting the posterior distributions of the fixed effects
#' @param res An object from ``inference_stan``, or from ``run_shiny`` where the "Close and run" option has been used
#' @keywords plot
#' @return A gg-plot with the posterior distributions.
#' @examples
#' vignette("plotting", package = "priorconstruction")

plot_posteriors_fixef <- function(res){

  df_post <- data.frame()

  samps <- rstan::extract(res$stan, c("intercept", "coeff"))

  if (res$prior$use_intercept){

    title <- expression(mu)

    xx <- seq(min(samps$intercept[,1]), max(samps$intercept[,1]), length.out = 200)

    df_post <- rbind(df_post, data.frame(x = samps$intercept[,1], param = as.character(title)))

  }

  for (i in seq_len(length(res$prior$data$fixed))){

    title <- as.expression(bquote(beta[.(names(res$prior$data$fixed)[i])]))

    xx <- seq(min(samps$coeff[,i]), max(samps$coeff[,i]), length.out = 200)

    df_post <- rbind(df_post, data.frame(x = samps$coeff[,i], param = as.character(title)))

  }

  gg <- ggplot() +
    geom_histogram(data = df_post, mapping = aes(x = x, y = ..density..), col = gray(0.5), fill = "#8E8D8A", bins = 40) +
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

  return(gg)

}


#' Plotting prior and posterior distributions
#'
#' Function for plotting the prior and posterior distributions of the fixed effects
#' @param res An object from ``inference_stan``, or from ``run_shiny`` where the "Close and run" option has been used
#' @keywords plot
#' @return A gg-plot with the prior and posterior distributions.
#' @examples
#' vignette("plotting", package = "priorconstruction")

plot_priors_posteriors_fixef <- function(res){

  df_pri <- data.frame()
  df_post <- data.frame()

  samps <- rstan::extract(res$stan, c("intercept", "coeff"))

  if (res$prior$use_intercept){

    title <- expression(mu)

    xx <- seq(min(samps$intercept[,1]), max(samps$intercept[,1]), length.out = 200)
    df_pri <- rbind(df_pri, data.frame(x = xx, y = dnorm(xx, mean = res$prior$prior_intercept[1],
                                                         sd = res$prior$prior_intercept[2]),
                                       param = as.character(title)))

    df_post <- rbind(df_post, data.frame(x = samps$intercept[,1], param = as.character(title)))

  }

  for (i in seq_len(length(res$prior$data$fixed))){

    title <- as.expression(bquote(beta[.(names(res$prior$data$fixed)[i])]))

    xx <- xx <- seq(min(samps$coeff[,i]), max(samps$coeff[,i]), length.out = 200)
    df_pri <- rbind(df_pri, data.frame(x = xx, y = dnorm(xx, mean = as.numeric(res$prior$prior_fixed[i, 1]),
                                                         sd = as.numeric(res$prior$prior_fixed[i, 2])),
                                       param = as.character(title)))

    df_post <- rbind(df_post, data.frame(x = samps$coeff[,i], param = as.character(title)))

  }

  gg <- ggplot() +
    geom_histogram(data = df_post, mapping = aes(x = x, y = ..density..), col = gray(0.5), fill = "#8E8D8A", bins = 40) +
    geom_line(data = df_pri, mapping = aes(x = x, y = y), na.rm = TRUE) +
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

  return(gg)

}


#' Plotting posterior distributions of the variances
#'
#' Function for plotting the posterior distributions of the variances of each random effect
#' (regardless of the parameterization chosen for the priors)
#' @param res An object from ``inference_stan``, or from ``run_shiny`` where the "Close and run" option has been used
#' @keywords plot
#' @return A gg-plot with the posterior distributions of the variances.
#' @examples
#' vignette("plotting", package = "priorconstruction")

# plot posterior of each variance parameter (internal parameters in the stan code)
plot_posterior_variances <- function(res){

  samps <- rstan::extract(res$stan, "theta")$theta

  titles <- sapply(as.character(res$prior$node_data$orig_nodedata$label), function(x) as.expression(bquote(sigma[.(x)]^2)))

  df <- data.frame(
    x = exp(unlist(c(samps))),
    name = factor(rep(as.character(titles), each = nrow(samps)), levels = titles)
  )

  gg1 <- ggplot(df, aes(x = x, y = ..density..)) +
    geom_histogram(col = gray(0.5), fill = "#8E8D8A", bins = 40) +
    facet_wrap(~name, labeller = label_parsed, scales = "free") +
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

  return(gg1)

}


# making data-frame with priors for weight and 1-weight
make_weight_prior_data_frame <- function(res){

  df_pri <- data.frame()

  for (i in seq_len(length(res$prior_data$weights))){
    if (res$prior_data$weights[[i]]$prior == "pc"){

      # weight in parameterization
      tmp <- plot_pc_w(res$prior_data$weights[[i]], res$weight_priors[[i]], res$node_data)
      df_pri <- rbind(df_pri, cbind(tmp, parent_node = res$prior_data$weights[[i]]$id))
      # 1- weight
      tmp$y <- rev(tmp$y)
      under <- paste0(get_node_name(res$node_data, get_leaf_nodes(res$node_data, res$prior_data$weights[[i]]$id)), sep = "", collapse = ",")
      not_base <- res$prior_data$weights[[i]]$children[!(res$prior_data$weights[[i]]$children %in% res$prior_data$weights[[i]]$param$basemodel_node)]
      over <- paste0(get_node_name(res$node_data, get_leaf_nodes(res$node_data, not_base)), sep = "", collapse = ",")
      tmp$param <- as.character(tmp$param)
      tmp$param <- as.character(as.expression(bquote(omega[frac(.(over), .(under))])))
      tmp$parent_node <- res$prior_data$weights[[i]]$id
      df_pri <- rbind(df_pri, tmp)

    } else { # dirichlet

      # weight in parameterization
      tmp <- plot_diri_w(res$prior_data$weights[[i]], res$node_data)
      tmp <- tmp[tmp$param == tmp$param[1],]

      titles <- character()
      under <- paste0(get_node_name(res$node_data, get_leaf_nodes(res$node_data, res$prior_data$weights[[i]]$id)), sep = "", collapse = ",")
      for (ind in 1:res$prior_data$weight[[i]]$no_children){
        titles[ind] <- as.expression(bquote(omega[frac(.(as.character(
          get_node_name(res$node_data, get_leaf_nodes(res$node_data, res$prior_data$weights[[i]]$children[ind])))), .(under))]))
      }

      df_pri <- rbind(df_pri, data.frame(x = tmp$x, y = tmp$y, param = rep(as.character(titles), each = nrow(tmp)),
                                         parent_node = res$prior_data$weights[[i]]$id))

    }
  }

  return(df_pri)

}

# making data-frame with posteriors for weight and 1-weight
make_weight_posterior_data_frame <- function(res){

  prior_obj <- res$prior
  res_stan <- res$stan
  stan_data <- res$stan_data

  df_post <- data.frame()

  samps <- rstan::extract(res_stan, "theta")$theta

  # for HD priors
  for (i in seq_len(length(prior_obj$prior_data$weights))){

    if (prior_obj$prior_data$weights[[i]]$prior == "pc"){

      # calculate the correct variances over and under fraction bar
      for (i2 in 1:2){

        ind2 <- stan_data$row_index_hd_pr_plot[i]+i2-1

        over <- rowSums(as.matrix(exp(samps[,stan_data$w_o[ind2,] == 1])))
        under <- rowSums(as.matrix(exp(samps[,stan_data$w_u[ind2,] == 1])))

        t_under <- paste0(as.character(get_node_name(prior_obj$node_data, get_leaf_nodes(prior_obj$node_data, prior_obj$prior_data$weights[[i]]$id))), sep = "", collapse = ",")
        if (i2 == 1){
          tmp_id <- prior_obj$prior_data$weights[[i]]$param$basemodel_node
        } else {
          tmp_id <- prior_obj$prior_data$weights[[i]]$children[!(prior_obj$prior_data$weights[[i]]$children %in% prior_obj$prior_data$weights[[i]]$param$basemodel_node)]
        }
        t_over <- paste0(as.character(get_node_name(prior_obj$node_data, get_leaf_nodes(prior_obj$node_data, tmp_id))), sep = "", collapse = ",")
        title <- as.expression(bquote(omega[frac(.(t_over), .(t_under))]))

        df_post <- rbind(df_post, data.frame(x = over/under, param = as.character(title), parent_node = prior_obj$prior_data$weights[[i]]$id))

      }

    } else { # dirichlet

      for (ind in 1:(prior_obj$prior_data$weights[[i]]$no_children)){

        ind2 <- stan_data$row_index_hd_pr_plot[i]+ind-1

        over <- rowSums(as.matrix(exp(samps[,stan_data$w_o[ind2,] == 1])))
        under <- rowSums(as.matrix(exp(samps[,stan_data$w_u[ind2,] == 1])))

        t_under <- paste0(as.character(get_node_name(prior_obj$node_data, get_leaf_nodes(prior_obj$node_data, prior_obj$prior_data$weights[[i]]$id))), sep = "", collapse = ",")
        t_over <- paste0(as.character(get_node_name(prior_obj$node_data, get_leaf_nodes(prior_obj$node_data, prior_obj$prior_data$weights[[i]]$children[ind]))), sep = "", collapse = ",")
        title <- as.expression(bquote(omega[frac(.(t_over), .(t_under))]))

        df_post <- rbind(df_post, data.frame(x = over/under, param = as.character(title), parent_node = prior_obj$prior_data$weights[[i]]$id))

      }

    }

  }

  return(df_post)

}


#' Plotting prior distributions of weights
#'
#' Function for plotting the prior distributions of all weights, also 1-weights for each split (which is not a parameter in the
#' model, but can be useful to see). Total variances and variances of the components with independent priors are not plotted.
#'
#' @param res An object from ``create_prior``, ``inference_stan``, or from ``run_shiny``
#' @keywords plot
#' @return A gg-plot with the posterior distributions of the weights.
#' @examples
#' vignette("plotting", package = "priorconstruction")

plot_all_prior_weights <- function(res){

  if (length(res) == 3) res <- res$prior

  df_pri <- make_weight_prior_data_frame(res)

  # from top level and down means that level number is increasing
  split_level <- get_split_ids_by_level(res$node_data, decreasing = FALSE)
  df_pri$parent_node <- factor(df_pri$parent_node, levels = split_level)

  gg <- list()
  # each split gets its own line in the plot
  for (i in 1:length(split_level)){
    gg[[i]] <- ggplot(subset(df_pri, parent_node == split_level[i]), aes(x = x, y = y)) + geom_line(na.rm = TRUE) +
      facet_wrap(~param, labeller = label_parsed, scales = "free") +
      ylab("Density") +
      ylim(0, NA) +
      theme(axis.title.x = element_blank()) +
      theme(panel.background = element_rect(fill ="white", color = gray(0.5)),
            panel.grid.major = element_line(color = gray(0.85)),
            panel.grid.minor = element_line(color = gray(0.92))) +
      theme(strip.background = element_rect(fill = "white", color = gray(0.5)),
            strip.text = element_text(color = "black", size = 15))

  }

  gg1 <- ggarrange(plotlist = gg, ncol = 1)
  return(gg1)

}


#' Plotting posterior distributions of weights
#'
#' Function for plotting the posterior distributions of all weights, also 1-weights which is not a parameter in the
#' model, but can be useful to see. Total variances and variances of the components with independent priors are not plotted.
#'
#' @param res An object from ``inference_stan``, or from ``run_shiny`` where the "Close and run" option has been used
#' @keywords plot
#' @return A gg-plot with the posterior distributions of the weights.
#' @examples
#' vignette("plotting", package = "priorconstruction")

plot_all_posterior_weights <- function(res){

  df_post <- make_weight_posterior_data_frame(res)

  # from top level and down means that level number is increasing
  split_level <- get_split_ids_by_level(res$prior$node_data, decreasing = FALSE)
  df_post$parent_node <- factor(df_post$parent_node, levels = split_level)

  gg <- list()
  # each split gets its own line in the plot
  for (i in 1:length(split_level)){
    gg[[i]] <- ggplot(subset(df_post, parent_node == split_level[i]), aes(x = x, y = ..density..)) +
      geom_histogram(col = gray(0.5), fill = "#8E8D8A", bins = 40) +
      facet_wrap(~param, labeller = label_parsed, scales = "free") +
      ylab("Density") +
      theme(axis.title.x = element_blank()) +
      theme(panel.background = element_rect(fill ="white", color = gray(0.5)),
            panel.grid.major = element_line(color = gray(0.85)),
            panel.grid.minor = element_line(color = gray(0.92))) +
      theme(strip.background = element_rect(fill = "white", color = gray(0.5)),
            strip.text = element_text(color = "black", size = 15))

  }

  gg1 <- ggarrange(plotlist = gg, ncol = 1)

  return(gg1)

}


#' Plotting prior and posterior distributions of weights
#'
#' Function for plotting the prior and posterior distributions of all weights, also 1-weights which is not a parameter in the
#' model, but can be useful to see. Total variances and variances of the components with independent priors are not plotted.
#'
#' @param res An object from ``inference_stan``, or from ``run_shiny`` where the "Close and run" option has been used
#' @keywords plot
#' @return A gg-plot with the prior and posterior distributions of the weights.
#' @examples
#' vignette("plotting", package = "priorconstruction")

plot_all_prior_posterior_weights <- function(res){

  df_pri <- make_weight_prior_data_frame(res$prior)
  df_post <- make_weight_posterior_data_frame(res)

  # from top level and down means that level number is increasing
  split_level <- get_split_ids_by_level(res$prior$node_data, decreasing = FALSE)
  df_pri$parent_node <- factor(df_pri$parent_node, levels = split_level)
  df_post$parent_node <- factor(df_post$parent_node, levels = split_level)

  gg <- list()
  # each split gets its own line in the plot
  for (i in 1:length(split_level)){
    gg[[i]] <- ggplot() +
      geom_histogram(data = subset(df_post, parent_node == split_level[i]), mapping = aes(x = x, y = ..density..), col = gray(0.5), fill = "#8E8D8A", bins = 40) +
      geom_line(data = subset(df_pri, parent_node == split_level[i]), mapping = aes(x = x, y = y), na.rm = TRUE) +
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

  gg1 <- ggarrange(plotlist = gg, ncol = 1)

  return(gg1)

}



#' Plotting the prior tree graph
#'
#' Function for plotting the prior tree graph (only prior). Requires visNetwork
#'
#' @param res An object from ``create_prior``, ``inference_stan``, or from ``run_shiny``
#' @param nodenames Custom names for each node (optional). Given as a named list with the default names as list names, and the
#' new names as list elements.
#' Do not need to provide all.
#' @keywords plot
#' @return A gg-plot with the tree graph
#' @examples
#' vignette("plotting", package = "priorconstruction")

plot_prior_tree <- function(res, nodenames){

  if (length(res) == 3) res <- res$prior

  nd <- list(x = res$node_data)

  nd$x <- update_basemodel_edges(nd$x, res$prior_data$weights)

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

  nn <- visNetwork(tmp_nodes, tmp_edges, background = "white")
  nn <- visHierarchicalLayout(nn)
  nn <- visEdges(nn, color = list(highlight = node_highlight_border_color),
                 arrowStrikethrough = F, selectionWidth = 1, smooth = FALSE,
                 font = list(color = node_border_color, size = 14, vadjust = 35,
                             # background = "#EAE7DC", # background does not move with vadjust
                             strokeWidth = 0, bold = TRUE)
  )
  nn <- visNodes(nn,
                 color = list(border = node_border_color,
                              highlight = list(border = node_highlight_border_color)),
                 borderWidth = 1,
                 #borderWidthSelected = 1,
                 font = list(background = "white", color = node_border_color))
  nn <- visOptions(nn, highlightNearest = FALSE)
  nn <- visInteraction(nn, F, F, F, F, F, F, F, F, F, F, F, zoomView = F)

  return(nn)

}



#' Plotting the posterior from INLA
#'
#' Function for plotting the prior and posterior distributions of all weights, also 1-weights which is not a parameter in the
#' model, but can be useful to see. Total variances and variances of the components with independent priors are not plotted.
#'
#' @param res An object from ``inference_inla``, or from ``run_shiny`` where the "Close and run" option has been used
#' @param type String indicating whether to plot the variances ("variance", default), standard deviations ("standard deviation") 
#' or precisions ("precision")
#' @param axes Names list with x-axis limits for each variance. Names are given as in the formula, the residuals are called "eps"
#' @keywords plot
#' @return A gg-plot with the posterior distributions of the variances (or precisions).
#' @examples
#' vignette("plotting", package = "priorconstruction")

plot_posteriors_inla <- function(res, type = c("variance", "standard deviation", "precision"), axes = list()){
  
  type <- match.arg(type)
  
  if (type == "variance"){
    myfun <- function(x) 1/x
  } else if (type == "standard deviation") {
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
      } else if (type == "standard deviation"){
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
    } else if (type == "standard deviation"){
      as.expression(bquote(sigma[.(thispar)]))
    } else {
      as.expression(bquote(1/sigma[.(thispar)]^2))
    }
    
    tmp <- as.data.frame(inla.tmarginal(myfun, res$inla$marginals.hyperpar[[ind]]))
    
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
  
  gg1 <- ggplot(df_post, aes(x = x, y = y)) +
    geom_line() +
    facet_wrap(~param, labeller = label_parsed, scales = "free") +
    # theme(strip.text = element_text(size = 14, color = "white")) +
    ylab("Density") +
    theme(axis.title.x = element_blank()) +
    theme(panel.background = element_rect(fill ="white", color = gray(0.5)),
          panel.grid.major = element_line(color = gray(0.85)),
          panel.grid.minor = element_line(color = gray(0.92))) +
    theme(strip.background = element_rect(fill = "white", color = gray(0.5)),
          strip.text = element_text(color = "black", size = 15))
  
  return(gg1)
  
}



### debugging only!

compare_stan_inla_r <- function(res_inla, res_stan){

  df_inla <- data.frame()
  df_stan <- data.frame()
  
  time_used <- paste0("INLA: ", round(res_inla$inla$cpu.used[4], 1), ", Stan: ", round(sum(get_elapsed_time(res_stan$stan)), 1))
  
  samps <- rstan::extract(res_stan$stan, "theta")$theta

  for (ind in seq_len(length(res_inla$prior$data$random))){
    
    tmp_inla <- cbind(
      as.data.frame(res_inla$inla$internal.marginals.hyperpar[[ind]]),
      param = strsplit(names(res_inla$inla$internal.marginals.hyperpar)[ind], "Log precision for ")[[1]][2]
    )
    
    tmp_stan <- data.frame(
      x = -samps[,ind],
      param = names(res_stan$prior$data$random)[ind]
    )
    
    df_inla <- rbind(df_inla, tmp_inla)
    df_stan <- rbind(df_stan, tmp_stan)

  }
  
  df_inla$param[df_inla$param == "the Gaussian observations"] <- "eps"
  
  df_inla$param <- factor(df_inla$param, levels = names(res_inla$prior$data$random))
  df_stan$param <- factor(df_stan$param, levels = names(res_inla$prior$data$random))
  
  ggplot() +
    geom_histogram(data = df_stan, mapping = aes(x = x, y = ..density..), fill = "lightblue", col = "black", bins = 60) +
    geom_line(data = df_inla, mapping = aes(x = x, y = y)) +
    facet_wrap(~ param, scales = "free") +
    ggtitle(time_used)

}

compare_stan_inla_f <- function(res_inla, res_stan){

  df_inla <- data.frame()
  df_stan <- data.frame()
  
  time_used <- paste0("INLA: ", round(res_inla$inla$cpu.used[4], 1), ", Stan: ", round(sum(get_elapsed_time(res_stan$stan)), 1))

  samps <- cbind(rstan::extract(res_stan$stan, "intercept")$intercept, rstan::extract(res_stan$stan, "coeff")$coeff)
  stanpar <- c("(Intercept)", names(res_stan$prior$data$fixed))

  for (ind in seq_len(length(res_inla$prior$data$fixed)+1)){

    tmp_inla <- cbind(
      as.data.frame(res_inla$inla$marginals.fixed[[ind]]),
      param = names(res_inla$inla$marginals.fixed)[ind]
    )

    tmp_stan <- data.frame(
      x = samps[,ind],
      param = stanpar[ind]
    )

    df_inla <- rbind(df_inla, tmp_inla)
    df_stan <- rbind(df_stan, tmp_stan)

  }

  ggplot() +
    geom_histogram(data = df_stan, mapping = aes(x = x, y = ..density..), fill = "lightblue", col = "black", bins = 60) +
    geom_line(data = df_inla, mapping = aes(x = x, y = y)) +
    facet_wrap(~ param, scales = "free") +
    ggtitle(time_used)

}

marginals_random_difference <- function(res_inla, res_stan){
  
  df <- data.frame()
  
  samps_eff <- rstan::extract(res_stan$stan, "random_effects")$random_effects
  samps_theta <- rstan::extract(res_stan$stan, "theta")$theta
  
  n_eff <- length(res_inla$inla$summary.random)
  n_in_eff <- sapply(1:n_eff, function(x) nrow(res_inla$inla$summary.random[[x]]))
  
  for (ind1 in seq_len(n_eff)) {
    for (ind2 in seq_len(n_in_eff[ind1])){
      
      from_inla <- res_inla$inla$summary.random[[ind1]][ind2, 2:6]
      # need to scale the random effects from stan, as we just return N(0,1), and must multiply with the variances and stuff
      
      post_stan <- c(samps_eff[,n_in_eff[ind1]*(ind1-1)+ind2])
      
      from_stan <- c(mean(post_stan), sd(post_stan), quantile(post_stan, c(0.025, 0.5, 0.975)))
      
      tmp <- as.data.frame(rbind(as.numeric(from_inla), as.numeric(from_stan), as.numeric(from_inla - from_stan)))
      names(tmp) <- c("mean", "sd", "low", "median", "high")
      tmp$type <- c("inla", "stan", "diff")
      tmp$eff <- names(res_inla$inla$summary.random)[ind1]
      tmp$effno <- as.character(ind2)
      tmp$eff2 <- paste0(names(res_inla$inla$summary.random)[ind1], ind2)
      
      df <- rbind(df, tmp)
      
    }
  }
  
  df2 <- subset(df, type == "diff")
  df2 <- df2[,-6]

  ggplot(reshape2::melt(df2, id = c("eff", "effno", "eff2")), aes(x = eff2, y = value, col = eff)) +
    geom_hline(yintercept = 0) +
    geom_point() +
    ggtitle("inla-stan") +
    facet_wrap(~ variable, scales = "free")
  
}


marginals_random_difference_plot <- function(res_inla, res_stan, maxeff = 10){
  
  samps_eff <- rstan::extract(res_stan$stan, "effects")$effects
  samps_theta <- rstan::extract(res_stan$stan, "theta")$theta
  
  n_eff <- length(res_inla$inla$summary.random)
  n_in_eff <- sapply(1:n_eff, function(x) nrow(res_inla$inla$summary.random[[x]]))
  
  df_s <- df_i <- data.frame()
  for (ind1 in 1:n_eff){
    for (ind2 in 1:n_in_eff[ind1]){
      
      set.seed(n_in_eff[ind1]*(ind1-1)+ind2); samps_used <- sample(1:nrow(samps_eff), 1000) # just using 1000 samples for efficiency
      tmp_s <- data.frame(
        x = c(samps_eff[,n_in_eff[ind1]*(ind1-1)+ind2])[samps_used],
        eff = names(res_inla$inla$summary.random)[ind1],
        effno = as.character(ind2),
        eff2 = paste0(names(res_inla$inla$summary.random)[ind1], ind2)
      )
      
      tmp_i <- data.frame(
        as.data.frame(res_inla$inla$marginals.random[[ind1]][[ind2]]),
        eff = names(res_inla$inla$summary.random)[ind1],
        effno = as.character(ind2),
        eff2 = paste0(names(res_inla$inla$summary.random)[ind1], ind2)
      )
      
      df_s <- rbind(df_s, tmp_s)
      df_i <- rbind(df_i, tmp_i)
      
    }
  }
  
  if (length(unique(df_s$effno)) > maxeff) {
    df_s <- df_s[df_s$effno %in% 1:maxeff,]
    df_i <- df_i[df_i$effno %in% 1:maxeff,]
  }
  
  ggplot() +
    geom_histogram(data = df_s, mapping = aes(x = x, y = ..density.., fill = eff), bins = 200) +
    geom_line(data = df_i, mapping = aes(x = x, y = y), col = "black") +
    facet_grid(effno ~ eff, scales = "free")
  
}










