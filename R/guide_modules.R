

opening_guide <- function(guide_data){

  gd <- list(x = guide_data)

  if (gd$x$ask_start_guide){

    tmp <- gd$x$counter_next
    gd$x <- initialize_guide(FALSE)
    gd$x$counter_next <- tmp

    showModal(modalDialog(
      paste("You will now be guided though your prior specification. First you specify the model structure using a node tree,",
            "then you choose prior distributions for each model components based on your knowledge about the model,",
            "and at last you get to see the resulting prior distributions based on your choices.",
            "To quit the guide, click 'Quit guide' at any time (your changes may be lost if you do this before the end)."
      ),
      h5(""),
      visNetwork::visNetworkOutput("graph_guide", height = 250),
      actionButton("begin_guide", "Begin guide", class = "button_g"),
      title = "Welcome!",
      footer = actionButton("quit_guide", "Quit guide", class = "button_quit")
    ))

    gd$x$ask_start_guide <- FALSE

  }

  return(list(guide_data = gd$x))

}


starting_guide <- function(guide_data, session){

  gd <- list(x = guide_data)

  tmp <- gd$x$counter_next
  gd$x <- initialize_guide(FALSE)
  gd$x$counter_next <- tmp

  shinyBS::updateCollapse(session, "sidebarpanel",
                 open = c("choose_structure"),#, "closing"), # latest actions and closing should always be visible
                 close = c("choose_priors", "param_expressions", "prior_expressions"))

  gd$x$guide_is_running <- TRUE
  gd$x$ask_start_guide <- FALSE

  # step 1-message
  removeModal(session)
  showModal(modalDialog(
    "We start by choosing the desired model structure. You can do this by clicking the buttons under the 'Modify model structure'-tab to the left. When you are happy with the structure, you can click the 'Next step'-button below the model structure tree. Click a node to choose it, and to choose more than one node, you can use a long-click (click and hold for a bit).",
    h5(""),
    actionButton("step1", "Got it", class = "button_g", value = FALSE),
    title = "Step 1: Modify prior structure",
    footer = actionButton("quit_guide",  "Quit guide", class = "button_quit")
  ))

  return(list(guide_data = gd$x, session = session))

}


go_to_step_2_guide <- function(prior_weight, guide_data, session){

  pd <- list(w = prior_weight)
  gd <- list(x = guide_data)

  #removeModal(session)
  showModal(modalDialog(
    paste("Now that you are happy with the model structure, you can choose prior distributions for each tree and each individual model component.",
          "Robust default priors are already chosen for you",
          "but if you have some knowledge on some of the components you should use this to improve your model.",
          "Do you want to explore the prior choices yourself, or do you want to be guided through all the choices in a step-wise manner?"
    ),
    h5(""),
    # if we have no trees, more straight to the CW part (but we do not need to tell the user about that)
    if (length(pd$w) > 0) actionButton("stepwise_wV", "Step-wise guide", class = "button1") else actionButton("stepwise_cw", "Step-wise guide", class = "button1"),
    actionButton("explore", "Explore myself", class = "button_g"),
    title = "Step 2: Choose priors",
    footer = actionButton("quit_guide", "Quit guide", class = "button_quit")
  ))

  gd$x$step1 <- FALSE
  gd$x$step2 <- TRUE

  shinyBS::updateCollapse(session, "sidebarpanel",
                 open = c("choose_priors"),#, "closing"), # latest actions and closing should always be visible
                 close = c("choose_structure", "param_expressions", "prior_expressions"))

  return(list(prior_weight = pd$w, guide_data = gd$x, session = session))

}


stepwise_wV_guide <- function(){

  #removeModal(session)
  showModal(modalDialog(
    visNetwork::visNetworkOutput("graph_guide", height = 250),
    "Now we will choose the priors for the components in the tree. You will be guided through each level in the node tree(s) from the bottom up. Then you will be guided through how to choose the priors for the total variance(s).",
    HTML("<br>"),
    actionButton("start_guide_splits", "Ok, let's go", class = "button_g"),
    title = "Step 2: Choose prior distributions",
    footer = actionButton("quit_guide", "Quit guide", class = "button_quit")
  ))

}


split_part_guide <- function(node_data, prior_weight, guide_data){

  nd <- list(x = node_data)
  pd <- list(w = prior_weight)
  gd <- list(x = guide_data)

  if (gd$x$stepwise && length(pd$w) > 0 && gd$x$which_split > 0){

    split_nodes <- get_split_ids_by_level(nd$x, TRUE)

    gd$x$pr_ind <- get_prior_number(pd$w, split_nodes[gd$x$which_split])

    # if this is a multisplit, we just give the user a message about which prior is used
    if (!is_dual_split_node(nd$x, list(current_node_id = split_nodes[gd$x$which_split]))){

      if (gd$x$next_modal == 1){

        # techically this should not be necessary
        pd$w[[gd$x$pr_ind]]$prior <- "dirichlet"
        pd$w[[gd$x$pr_ind]]$param <- get_diri_param(pd$w[[gd$x$pr_ind]]$no_children)

        #removeModal(session)
        showModal(modalDialog(
          visNetwork::visNetworkOutput("graph_guide", height = 250),
          paste("This is a multi-split, and this split is given a symmetric Dirichlet prior with equal amount of variance to each child node.", sep = "", collapse = ""),
          HTML("<br>"),
          if (gd$x$which_split < length(split_nodes)){
            actionButton(paste0("next_split", gd$x$counter_next), "Ok, next split", class = "button_g")
          } else {
            actionButton("go_to_step_2_V", "Ok, I am ready for the total variances", class = "button_g")
          },
          title = "Step 2: Choose split prior",
          footer = actionButton("quit_guide", "Quit guide", class = "button_quit")
        ))

      }

    } else {

      if (gd$x$next_modal == 1){

        #removeModal(session)
        showModal(modalDialog(
          visNetwork::visNetworkOutput("graph_guide", height = 250),
          paste("Do you have an intuition on how the variance of ", get_node_name(nd$x, split_nodes[gd$x$which_split]),
                " is distributed, and do you want to use this knowledge?", sep = "", collapse = ""),
          HTML("<br>"),
          actionButton(paste0("use_pc_w", gd$x$counter_next), "Yes", class = "button_g"),
          actionButton(paste0("use_diri", gd$x$counter_next), "No", class = "button_g"),
          title = "Step 2: Choose split prior",
          footer = actionButton("quit_guide", "Quit guide", class = "button_quit")
        ))
        gd$x$next_modal <- 2
      }
    }
  }

  return(list(prior_weight = pd$w, guide_data = gd$x))

}


chosen_dirichlet_guide <- function(node_data, prior_weight, guide_data){

  nd <- list(x = node_data)
  pd <- list(w = prior_weight)
  gd <- list(x = guide_data)

  split_nodes <- get_split_ids_by_level(nd$x, TRUE)

  if (gd$x$next_modal == 2){

    # techically this should not be necessary
    pd$w[[gd$x$pr_ind]]$prior <- "dirichlet"
    pd$w[[gd$x$pr_ind]]$param <- 1

    #removeModal(session)
    showModal(modalDialog(
      visNetwork::visNetworkOutput("graph_guide", height = 250),
      "Then this split gets a symmetric Dirichlet prior with equal amount of variance to each child node.",
      HTML("<br>"),
      if (gd$x$which_split < length(split_nodes)){
        actionButton(paste0("next_split", gd$x$counter_next), "Ok, next split", class = "button_g")
      } else {
        actionButton("go_to_step_2_V", "Ok, I am ready for the total variances", class = "button_g")
      },
      title = "Step 2: Choose split prior",
      footer = actionButton("quit_guide", "Quit guide", class = "button_quit")
    ))

  }

  return(list(prior_weight = pd$w,guide_data = gd$x))

}

# Asks about the basemodel
chosen_pc_guide <- function(node_data, prior_weight, guide_data){

  nd <- list(x = node_data)
  pd <- list(w = prior_weight)
  gd <- list(x = guide_data)

  if (gd$x$next_modal == 2){
    pd$w[[gd$x$pr_ind]]$prior <- "pc"
    pd$w[[gd$x$pr_ind]]$param <- data.frame(above_node = pd$w[[gd$x$pr_ind]]$children[1], median = 0.25, basemodel = 0.25, concentration = 0.5)
    #removeModal(session)
    showModal(modalDialog(
      visNetwork::visNetworkOutput("graph_guide", height = 250),
      HTML(paste("Then we use a PC prior on this split. Which of the effects ",
            paste(get_node_name(nd$x, pd$w[[gd$x$pr_ind]]$children), sep = "", collapse = " and "),
            " should be preferred? (This means that we have shrinkage to this effect,",
            " and is the same as choosing the <b>basemodel</b> of the prior.)",
            sep = "", collapse = ""
      )),
      HTML("<br>"),
      actionButton(paste0("yes_base1", gd$x$counter_next), paste("EFFECT", get_node_name(nd$x, pd$w[[gd$x$pr_ind]]$children[1])), class = "button_g2"),
      actionButton(paste0("yes_base2", gd$x$counter_next), paste("EFFECT", get_node_name(nd$x, pd$w[[gd$x$pr_ind]]$children[2])), class = "button_g2"),
      actionButton(paste0("no_base", gd$x$counter_next), "A combination", class = "button_g"),
      title = "Step 2: Choose split prior",
      footer = actionButton("quit_guide", "Quit guide", class = "button_quit")
    ))
    gd$x$next_modal <- 3
  }

  return(list(prior_weight = pd$w, guide_data = gd$x))

}

median_pc_guide <- function(node_data, prior_weight, guide_data, input, basemod_value){

  nd <- list(x = node_data)
  pd <- list(w = prior_weight)
  gd <- list(x = guide_data)

  if (gd$x$next_modal == 3){

    pd$w[[gd$x$pr_ind]]$param$basemodel <- basemod_value

    #removeModal(session)
    showModal(modalDialog(
      visNetwork::visNetworkOutput("graph_guide", height = 250),
      HTML(paste("How much of the variance of ",
                 pd$w[[gd$x$pr_ind]]$name,
                 " do you think should go to ",
                 get_node_name(nd$x, pd$w[[gd$x$pr_ind]]$children[1]),
                 "? (This is the same as choosing the <b>median</b> of the prior.)",
                 sep = "", collapse = ""
      )),
      HTML("<br>"),
      numericInput("median_guide",
                   #"Median",
                   "",
                   value = 0.5, min = 0, max = 1, step = 0.05, width = "200px"),
      HTML("<br>"),
      if (basemod_value %in% c(0, 1)){
        actionButton(paste0("set_median_01", gd$x$counter_next), "Ok", class = "button_g")
      } else {
        actionButton(paste0("set_median", gd$x$counter_next), "Ok", class = "button_g")
      },
      title = "Step 2: Choose split prior",
      footer = actionButton("quit_guide", "Quit guide", class = "button_quit")
    ))
    gd$x$next_modal <- if (basemod_value %in% c(0, 1)) 4 else 3.5
  }

  return(list(prior_weight = pd$w, guide_data = gd$x))

}

concentration_pc_guide <- function(node_data, prior_weight, guide_data, input){

  nd <- list(x = node_data)
  pd <- list(w = prior_weight)
  gd <- list(x = guide_data)

  if (gd$x$next_modal == 3.5){

    if (!(pd$w[[gd$x$pr_ind]]$param$basemodel %in% c(0, 1))) {

      pd$w[[gd$x$pr_ind]]$param$median <- input$median_guide

      #removeModal(session)
      showModal(modalDialog(
        visNetwork::visNetworkOutput("graph_guide", height = 250),
        HTML(paste("How concentrated do you want the prior for the amount of variance of ",
                   pd$w[[gd$x$pr_ind]]$name, " to ",
                   get_node_name(nd$x, pd$w[[gd$x$pr_ind]]$children[1]),
                   " to be?",
                   sep = "", collapse = ""
        )),
        HTML("<br>"),
        numericInput("conc_param_guide",
                     "",
                     value = 0.5, min = 0.5, max = 1, step = 0.05, width = "200px"),
        HTML("<br>"),
        actionButton(paste0("set_conc_param", gd$x$counter_next), "Ok", class = "button_g"),
        title = "Step 2: Choose split prior",
        footer = actionButton("quit_guide", "Quit guide", class = "button_quit")
      ))
      gd$x$next_modal <- 4
    }
  }

  return(list(prior_weight = pd$w, guide_data = gd$x))

}

chosen_pc_finished <- function(node_data, prior_weight, guide_data, input){

  nd <- list(x = node_data)
  pd <- list(w = prior_weight)
  gd <- list(x = guide_data)

  split_nodes <- get_split_ids_by_level(nd$x, TRUE)

  if (gd$x$next_modal == 4){

    pd$w[[gd$x$pr_ind]]$param$median <- input$median_guide
    pd$w[[gd$x$pr_ind]]$param$concentration <- input$conc_param_guide

    if (pd$w[[gd$x$pr_ind]]$param$basemodel == 1){

      showModal(modalDialog(
        visNetwork::visNetworkOutput("graph_guide", height = 250),
        paste("Then this split gets a PC prior with shrinkage to ", get_node_name(nd$x, pd$w[[gd$x$pr_ind]]$children[1]),
              " and median ", pd$w[[gd$x$pr_ind]]$param$median, ".",
              sep = "", collapse = ""),
        HTML("<br>"),
        if (gd$x$which_split < length(split_nodes)){
          actionButton(paste0("next_split", gd$x$counter_next), "Ok, next split", class = "button_g")
        } else {
          actionButton("go_to_step_2_V", "Ok, I am ready for the total variances", class = "button_g")
        },
        title = "Step 2: Choose split prior",
        footer = actionButton("quit_guide", "Quit guide", class = "button_quit")
      ))


    } else if (pd$w[[gd$x$pr_ind]]$param$basemodel == 0){

      showModal(modalDialog(
        visNetwork::visNetworkOutput("graph_guide", height = 250),
        paste("Then this split gets a PC prior with shrinkage to ", get_node_name(nd$x, pd$w[[gd$x$pr_ind]]$children[2]),
              " and median ", pd$w[[gd$x$pr_ind]]$param$median, ".",
              sep = "", collapse = ""),
        HTML("<br>"),
        if (gd$x$which_split < length(split_nodes)){
          actionButton(paste0("next_split", gd$x$counter_next), "Ok, next split", class = "button_g")
        } else {
          actionButton("go_to_step_2_V", "Ok, I am ready for the total variances", class = "button_g")
        },
        title = "Step 2: Choose split prior",
        footer = actionButton("quit_guide", "Quit guide", class = "button_quit")
      ))

    } else { # basemodel = median

      pd$w[[gd$x$pr_ind]]$param$basemodel <- input$median_guide

      #removeModal(session)
      showModal(modalDialog(
        visNetwork::visNetworkOutput("graph_guide", height = 250),
        paste("Then this split gets a PC prior with shrinkage to the median (which is set to ",
              pd$w[[gd$x$pr_ind]]$param$median, "), and concentration parameter ",
              pd$w[[gd$x$pr_ind]]$param$concentration, ".",
              sep = "", collapse = ""),
        h5(""),
        if (gd$x$which_split < length(split_nodes)){
          actionButton(paste0("next_split", gd$x$counter_next), "Ok, next split", class = "button_g")
        } else {
          actionButton("go_to_step_2_V", "Ok, I am ready for the total variances", class = "button_g")
        },
        title = "Step 2: Choose split prior",
        footer = actionButton("quit_guide", "Quit guide", class = "button_quit")
      ))

    }

  }

  return(list(prior_weight = pd$w))

}



totvar_part_guide <- function(node_data, prior_totvar, guide_data, .initial_args){

  nd <- list(x = node_data)
  pd <- list(V = prior_totvar)
  gd <- list(x = guide_data)

  likl_name <- if (.initial_args$.family == "poisson") "Poisson" else .initial_args$.family

  if (gd$x$stepwise && length(pd$V) > 0 && gd$x$which_totvar > 0){

    top_nodes <- get_top_nodes(nd$x)
    gd$x$pr_ind <- get_prior_number(pd$V, top_nodes[gd$x$which_totvar])

    #removeModal(session)
    showModal(modalDialog(
      visNetwork::visNetworkOutput("graph_guide", height = 250),
      paste("Which prior do you want on ", pd$V[[gd$x$pr_ind]]$name, "? ",
            if (length(pd$V) > 1){
              "Note that you cannot choose Jeffreys' prior for any total variance, because you have more than one tree."
            },
            if (sum(nd$x$nodes$status == "detached") > 0){
              "Note that you cannot choose Jeffreys' prior for any total variance, because you have some model components that are not attached to the tree (they will get independent priors)."
            },
            if (.initial_args$.family != "gaussian"){
              paste0("Note that you cannot choose Jeffreys' prior because you have a ", likl_name, " likelihood.")
            },
            sep = "", collapse = ""),
      HTML("<br>"),
      fluidRow(
        column(6,
               if (length(pd$V) > 1 || sum(nd$x$nodes$status == "detached") > 0 || .initial_args$.family != "gaussian"){
                 radioButtons("var_prior_guide",
                              "",
                              choices = list("PC prior" = "pc0",
                                             "Inverse gamma" = "invgam",
                                             "Half-Cauchy" = "hc"),
                              selected = "pc0")
               } else {
                 radioButtons("var_prior_guide",
                              "",
                              choices = list("Jeffreys'" = "jeffreys",
                                             "PC prior" = "pc0",
                                             "Inverse gamma" = "invgam",
                                             "Half-Cauchy" = "hc"),
                              selected = "jeffreys")
               }),
        div(style = "min-height: 122px; max-height = 122px", column(6,
                                                                    uiOutput("par1_box_guide"),
                                                                    uiOutput("par2_box_guide")
        ))
      ),
      HTML("<br>"),
      if (gd$x$which_totvar < length(top_nodes)){
        actionButton(paste0("next_totvar", gd$x$counter_next), "Ok, next total variance", class = "button_g")
      } else {
        if (sum(nd$x$nodes$status == "detached") > 0){
          actionButton("stepwise_cw", "Ok, I am ready for the independent nodes", class = "button_g")
        } else {
          actionButton("go_to_step_3_1", "Ok, let me see the prior distributions", class = "button_g")
        }
      },
      title = paste("Step 2: Choose total variance prior", if (length(pd$V) > 1) "s", sep = "", collapse = ""),
      footer = actionButton("quit_guide", "Quit guide", class = "button_quit")
    ))

  }

  return(list(guide_data = gd$x))

}


cw_part_guide <- function(node_data, prior_cw, guide_data){

  nd <- list(x = node_data)
  pd <- list(var = prior_cw)
  gd <- list(x = guide_data)

  if (gd$x$stepwise && length(get_detached_nodes(nd$x)) > 0 && gd$x$which_cw > 0){

    cw_nodes <- get_detached_nodes(nd$x)

    gd$x$pr_ind <- get_prior_number(pd$var, cw_nodes[gd$x$which_cw])

    #removeModal(session)
    showModal(modalDialog(
      visNetwork::visNetworkOutput("graph_guide", height = 250),
      paste("Which prior do you want on ", pd$var[[gd$x$pr_ind]]$name, "?", sep = "", collapse = ""),
      HTML("<br>"),
      fluidRow(
        column(6,
               radioButtons("var_prior_guide",
                            "",
                            choices = list("PC prior" = "pc0",
                                           "Inverse gamma" = "invgam",
                                           "Half-Cauchy" = "hc"),
                            selected = "pc0")
        ),
        div(style = "min-height: 122px; max-height = 122px", column(6,
                                                                    uiOutput("par1_box_guide"),
                                                                    uiOutput("par2_box_guide")
        ))
      ),
      if (gd$x$which_cw < length(cw_nodes)){
        actionButton(paste0("next_cw", gd$x$counter_next), "Ok, next variance", class = "button_g")
      } else {
        actionButton("go_to_step_3_2", "Ok, let me see the prior distributions", class = "button_g")
      },
      title = paste("Step 2: Choose variance prior", if (length(pd$var) > 1) "s", sep = "", collapse = ""),
      footer = actionButton("quit_guide", "Quit guide", class = "button_quit")
    ))

  }

  return(list(guide_data = gd$x))

}


finishing_guide <- function(guide_data, session){

  gd <- list(x = guide_data)

  if (gd$x$finished){

    gd$x$stepwise <- FALSE
    gd$x$finished <- FALSE

    removeModal(session)
    showModal(modalDialog(
      paste("These are the priors you have chosen for each component. The guide is now finished. If you want to restart the guide and do everything again,",
            "click 'Restart guide'. You can also explore more yourself, by using the panels to the left. Note that every time you change the tree structure,",
            "the priors for the splits are reset to symmetric Dirichlet priors. You can start the guide again at any time by clicking 'Begin guide' in",
            "the panel to the left. When you are happy with the priors, you close the window and can start the inference."
      ),
      plotOutput("all_plots_guide"),
      #p(""),
      #uiOutput("prior_dists_guide"),
      # p(""),
      actionButton("begin_guide", "Restart guide", class = "button_g"),
      title = "Finished!",
      footer = actionButton("quit_guide", "Quit guide and save prior", class = "button_quit")
    ))

  }

  return(list(guide_data = gd$x, session = session))

}



make_par1_box_guide <- function(id, input){

  renderUI({
    if (input$var_prior_guide == "pc0"){
      numericInput(id, h6("U"), 3, 0, Inf, step = 1)
    } else if (input$var_prior_guide == "invgam"){
      numericInput(id, h6("Shape"), 1, 0, Inf, step = 1)
    } else if (input$var_prior_guide == "hc"){
      numericInput(id, h6("Scale"), 25, 0, Inf, step = 1)
    } else {
      # nothing
    }
  })

}

make_par2_box_guide <- function(id, input){

  renderUI({

    if (input$var_prior_guide == "pc0"){
      numericInput(id, h6("alpha"), 0.05, 0, 1, step = 0.01)
    } else if (input$var_prior_guide == "invgam"){
      numericInput(id, h6("Scale"), 5e-5, 0, Inf, step = 0.01)
    } else if (input$var_prior_guide == "hc"){
      # nothing
    } else {
      # nothing
    }

  })

}



make_graph_guide <- function(node_data, prior_weight, prior_totvar, prior_cw, guide_data){

  nd <- list(x = node_data)
  pd <- list(w = prior_weight, V = prior_totvar, var = prior_cw)
  gd <- list(x = guide_data)

  tmp_nodes <- nd$x$nodes
  tmp_nodes$color.background <- node_palette(0.35)
  #tmp_nodes$color.border <- node_border_color

  tmp_edges <- nd$x$edges
  tmp_edges$width <- if (nrow(tmp_edges) > 0) 2 else NULL
  tmp_edges$dashes <- if (nrow(tmp_edges) > 0) FALSE else NULL
  tmp_edges$label <- NULL
  if (nrow(tmp_edges) > 0) tmp_edges <- cbind(tmp_edges, arrows = "to", hidden = FALSE)

  #tmp_nodes <- set_node_title(tmp_nodes, nd$x)

  # make sure the detached nodes are not all over (does not work great...)
  tmp_edges <- set_detached_edges(tmp_nodes, tmp_edges)

  # for each step in the guide, we highlight relevant areas
  if (gd$x$stepwise){

    tmp_nodes$color.border <- gray(0.5)

    if (gd$x$which_split > 0){ # if we are going through the splits

      split_node <- pd$w[[gd$x$pr_ind]]$id
      children <- pd$w[[gd$x$pr_ind]]$children
      # split_node <- pd$w[[gd$x$which_split]]$id
      # children <- pd$w[[gd$x$which_split]]$children

      tmp_nodes$color.background[tmp_nodes$id == split_node] <- guide_highlight_parent_color
      tmp_nodes$color.background[tmp_nodes$id %in% children] <- guide_highlight_children_color

      tmp_nodes$color.border[tmp_nodes$id %in% c(split_node, children)] <- "black" #node_border_color
      tmp_edges$color.color[tmp_edges$from == split_node] <- "black" #node_border_color
      tmp_edges$color.inherit[!(tmp_edges$from == split_node)] <- FALSE

    } else if (gd$x$which_totvar > 0){ # if we are going through the total variance(s)

      # totvar_node <- pd$V[[gd$x$pr_ind]]$id
      totvar_node <- pd$V[[gd$x$which_totvar]]$id
      tmp_nodes$color.background[tmp_nodes$id == totvar_node] <- guide_highlight_parent_color

      tmp_nodes$color.border[tmp_nodes$id == totvar_node] <- "black" #node_border_color
      tmp_edges$color.color <- gray(0.5)

    } else if (gd$x$which_cw > 0){ # if we are going through the CW priors

      print(gd$x$which_cw)
      # var_node <- pd$var[[gd$x$pr_ind]]$id
      var_node <- pd$var[[gd$x$which_cw]]$id
      tmp_nodes$color.background[tmp_nodes$id == var_node] <- guide_highlight_parent_color

      tmp_nodes$color.border[tmp_nodes$id == var_node] <- "black" #node_border_color

    }

  }

  #nn <- plot_network()(tmp_nodes, tmp_edges, "small")
  nn <- visNetwork::visNetwork(tmp_nodes, tmp_edges)
  nn <- visNetwork::visHierarchicalLayout(nn)
  nn <- visNetwork::visEdges(nn, color = list(color = node_border_color), arrowStrikethrough = FALSE)
  nn <- visNetwork::visNodes(nn,
                 color = list(border = node_border_color),
                 borderWidth = 1)
  nn <- visNetwork::visNodes(nn, size = 30, font = list(size = 20, background = "#FFFFFF", color = node_border_color))
  #nn <- visHierarchicalLayout(nn, nodeSpacing = 50, levelSeparation = 50, treeSpacing = 50)
  nn <- visNetwork::visInteraction(nn,
                       dragNodes = FALSE,
                       dragView = FALSE,
                       zoomView = FALSE,
                       selectable = FALSE)
  return(nn)

}










