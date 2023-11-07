






make_output <- function(input_obj, prior_data, node_data, data, run_inference_directly = FALSE){

  output_list <- list()
  output_list$prior_data <- prior_data # includes all prior data (totvar, weights and independent variance components)
  node_data <- add_remove_AB(node_data, FALSE) # in case some of the node labels have A: and B:
  output_list$node_data <- node_data
  output_list$node_data$orig_nodedata <- get_original_nodes(node_data)
  # output_list$weight_priors <- calculate_pc_prior(node_data, prior_data$weights, data$covmats)
  output_list$run_inference_directly <- run_inference_directly
  return(output_list)

}


initialize_nodes <- function(initial_args){

  tmp <- list(nodes = initial_args$.initial_nodes, edges = initial_args$.initial_edges)
  tmp <- update_basemodel_edges(tmp, initial_args$.initial_prior_w)

  return(tmp)

}

#initialize_prior_w <- function(initial_args) return(initial_args$.initial_prior_w)

initialize_guide <- function(start_guide){

  guide_data <- list(
    ask_start_guide = start_guide,
    guide_is_running = FALSE,
    step1 = FALSE,
    step2 = FALSE,
    stepwise = FALSE,
    which_split = 0,
    which_totvar = 0,
    which_cw = 0,
    counter_next = 0, # to keep track of the next-buttons
    pr_ind = 0,
    next_modal = 0,
    finished = FALSE,
    message_to_user = NULL
  )

  return(guide_data)

}

server <- function(input, output, session) {

  # .initial_args <- reactiveValues(x = NULL)
  # observe({
  #   if (is.null(.initial_args$x)) .initial_args$x <- shiny::getShinyOption(".initial_args")
  # })

  # initializing the prior data based on user input (prior_data)
  pd <- reactiveValues(w = NULL, V = NULL, var = NULL)
  observe({
    if (is.null(pd$w)) pd$w <- shiny::getShinyOption(".initial_args")$.initial_prior_w
    if (is.null(pd$V)) pd$V <- shiny::getShinyOption(".initial_args")$.initial_prior_V
    if (is.null(pd$var)) pd$var <- shiny::getShinyOption(".initial_args")$.initial_prior_var
  })

  # initializing the tree based on user input (node_data)
  nd <- reactiveValues(x = NULL)
  observe({
    if (is.null(nd$x)) nd$x <- initialize_nodes(shiny::getShinyOption(".initial_args")) #(.initial_args)
  })

  td <- reactiveValues(x = NULL) # for message to user
  node_out <- reactiveValues(x = NULL)
  gd <- reactiveValues(x = initialize_guide(shiny::getShinyOption(".initial_args")$.guide)) # keeping track of the guide

  # for debugging
  output$nodes_data_from_shiny <- renderText({
    input$current_node_id
  })
  output$input_out <- renderPrint({
    c(input$weight_prior, input$median, input$basemodel)
  })
  output$node_data_out <- renderPrint({
    nd$x
  })
  output$prior_data_out <- renderPrint({
    list(w = pd$w, V = pd$V, var = pd$var, gd = gd$x)
  })

  observe({
    showNotification(
      paste("This is an early version, which means it probably contains bugs.",
            "Please contact ingeborg.hem@ntnu.no to report bugs and feature suggestions.",
            "If you encounter bugs, try closing the app and starting over."),
      duration = NULL, type = "warning"
    )
  })

  # none of these buttons should be clickable when no nodes are clicked
  observe({
    if (length(input$current_node_id) == 0) disable_everything()
  })

  # things that should happen if/when a node is clicked
  observeEvent(input$current_node_id, {
    clicked_node_action(nd$x, pd$w, pd$V, pd$var, td$x, input, session)
    tmp <- if (length(input$current_node_id) == 0) "" else if (length(input$current_node_id)) "Chosen node: " else "Chosen nodes: "
    model_data_tmp <- shiny::getShinyOption(".initial_args")$.latent
    if (shiny::getShinyOption(".initial_args")$.family == "gaussian"){
      model_data_tmp <- c(model_data_tmp, list(list(id = get_node_id(nd$x, "eps"), label = "eps", model = "residuals")))
    }
    node_out$x <- paste0(tmp,
                         paste0(
                           get_node_name_click(nd$x, input),
                           get_node_model_name_click(nd$x, model_data_tmp, input),
                           sep = "", collapse = ", "))
    if (shiny::getShinyOption(".initial_args")$.family != "gaussian") disable(selector = "[type=radio][value=jeffreys]")
  })

  output$chosen_node <- renderText({
    node_out$x
  })

  # modifying tree structure

  observeEvent(input$detach, {
    tmp <- detach_action(nd$x, pd$w, pd$V, pd$var, td$x, input)
    nd$x <- tmp$node_data
    pd$w <- tmp$prior_weight
    pd$V <- tmp$prior_totvar
    pd$var <- tmp$prior_cw
    node_out$x <- ""

  })

  observeEvent(input$merge, {

    tmp <- merge_action(nd$x, pd$w, pd$V, pd$var, td$x, input)
    nd$x <- tmp$node_data
    pd$w <- tmp$prior_weight
    pd$V <- tmp$prior_totvar
    pd$var <- tmp$prior_cw
    node_out$x <- ""

  })

  observeEvent(input$remove, {

    tmp <- remove_action(nd$x, pd$w, pd$V, pd$var, td$x, input)
    nd$x <- tmp$node_data
    pd$w <- tmp$prior_weight
    pd$V <- tmp$prior_totvar
    pd$var <- tmp$prior_cw
    node_out$x <- ""

  })

  observeEvent(input$attach, {

    tmp <- attach_action(nd$x, pd$w, pd$V, pd$var, td$x, input)
    nd$x <- tmp$node_data
    pd$w <- tmp$prior_weight
    pd$V <- tmp$prior_totvar
    pd$var <- tmp$prior_cw
    node_out$x <- ""

  })

  # change the prior_data-object when the user chooses another prior distribution family for a split
  observeEvent(length(input$weight_prior) & input$median & input$conc_param & is.null(input$basemodel), {
    req(input$median, input$basemodel)
    req(input$median, input$basemodel, input$conc_param)
    pd$w <- update_prior_w(pd$w, nd$x, input)
    if (input$weight_prior == "pc"){
      children <- get_children(nd$x, input)
      median_value <- if (get_prior_info(pd$w, nd$x, input)$prior == "pc") get_prior_info(pd$w, nd$x, input)$param$median else 0.25
      updateRadioButtons(session, "basemodel",
                         choiceNames = c("100% A + 0% B", "0% A + 100% B",
                                         paste0(round(median_value, 2)*100, "% A + ", round(1-median_value, 2)*100, "% B")),
                         choiceValues = c("mod1", "mod2", "comb"),
                         selected = input$basemodel
      )
      if (!gd$x$guide_is_running) showElement("intrinsic")
    } else {
      hideElement("intrinsic")
    }
  })

  observeEvent(length(input$weight_prior) & is.null(input$basemodel), {
    #req(input$basemodel)
    if (input$weight_prior == "pc"){
      showElement("basemodel")
      showElement("median")
      if (input$basemodel == "comb") showElement("conc_param") else hideElement("conc_param")
    } else { # dirichlet
      hideElement("basemodel")
      hideElement("median")
      hideElement("conc_param")
    }
  })

  # if the prior object is changed, we (may) need to update the edges as well
  observeEvent(pd$w, {
    nd$x <- update_basemodel_edges(nd$x, pd$w)
  })

  # update if the parameter values change
  observeEvent(input$par1 & input$par2, {
    if (length(input$current_node_id) == 1 && is_top_node(nd$x, input)){
      pd$V <- update_var_prior(pd$V, input)
    }
    if (length(input$current_node_id) == 1 && are_detached(nd$x, input)){
      pd$var <- update_var_prior(pd$var, input)
    }
  })

  # update the prior family on total variance or CW priors, and if the family is changed, reset the numeric values
  observeEvent(is.null(input$var_prior), {
  # observeEvent(is.null(input$var_prior) & input$par1 & input$par2, {
    if (length(input$current_node_id) == 1 && is_top_node(nd$x, input)){
      # we only reset the parameters if the prior family is changed
      old_family <- get_prior_info(pd$V, nd$x, input)$prior
      pd$V <- update_var_prior(pd$V, input)
      new_family <- get_prior_info(pd$V, nd$x, input)$prior
      if (old_family == new_family){
        new_par1 <- get_prior_info(pd$V, nd$x, input)$param[1]
        new_par2 <- get_prior_info(pd$V, nd$x, input)$param[2]
      } else {
        new_par1 <- get_par1_default_val(new_family)
        new_par2 <- get_par2_default_val(new_family)
      }
      updateNumericInput(session, "par1", value = new_par1,
                         label = get_par1_label(new_family))
      updateNumericInput(session, "par2", value = new_par2,
                         max = get_par2_max(new_family), step = get_par2_step(new_family),
                         label = get_par2_label(new_family))
    }
    if (length(input$current_node_id) == 1 && are_detached(nd$x, input)){
      # we only reset the parameters if the prior family is changed
      old_family <- get_prior_info(pd$var, nd$x, input)$prior
      pd$var <- update_var_prior(pd$var, input)
      new_family <- get_prior_info(pd$var, nd$x, input)$prior
      if (old_family == new_family){
        new_par1 <- get_prior_info(pd$var, nd$x, input)$param[1]
        new_par2 <- get_prior_info(pd$var, nd$x, input)$param[2]
      } else {
        new_par1 <- get_par1_default_val(new_family)
        new_par2 <- get_par2_default_val(new_family)
      }
      updateNumericInput(session, "par1", value = new_par1,
                         label = get_par1_label(new_family))
      updateNumericInput(session, "par2", value = new_par2,
                         max = get_par2_max(new_family), step = get_par2_step(new_family),
                         label = get_par2_label(new_family))
    }
    if (input$var_prior %in% c("pc0", "invgam")){
      showElement("par1")
      showElement("par2")
    } else if (input$var_prior == "hc"){
      showElement("par1")
      hideElement("par2")
    } else {
      hideElement("par1")
      hideElement("par2")
    }
  })

  ##################################################### guide options

  ## is not possible to change tree structure during part 2 (prior choices) of the guide,
  ## so can directly update the prior, do not need to check for tree structure changes

  # make button for beginning guide if it is not running
  output$begin_quit_guide <- renderUI({

    if (!(!is.null(gd$x$guide_is_running) && gd$x$guide_is_running)){
      actionButton("begin_guide_again", "Begin guide", style = "width:140px", class = "button_g", title = "Begin the prior specification guide.", value = NULL)
    }

  })

  # start guide when user opens window
  observe({
    tmp <- opening_guide(gd$x)
    gd$x <- tmp$guide_data
  })

  # when guide is restarted from window
  observeEvent(input$begin_guide_again, {
    tmp <- gd$x$counter_next
    gd$x <- initialize_guide(FALSE)
    gd$x$counter_next <- tmp
    gd$x$ask_start_guide <- TRUE
  })

  # when guide starts
  observeEvent(input$begin_guide, {

    tmp <- starting_guide(gd$x, session)
    gd$x <- tmp$guide_data
    session <- tmp$session

  })


  # when user continues from step 1-message
  observeEvent(input$step1, {
    removeModal(session)
    gd$x$step1 <- TRUE
  })

  # when user moves to step 2
  observeEvent(input$go_to_step_2, {

    tmp <- go_to_step_2_guide(pd$w, gd$x, session)
    pd$w <- tmp$prior_weight
    gd$x <- tmp$guide_data
    session <- tmp$session

  })


  # if the user wants to explore the priors (step 2) with no guide
  observeEvent(input$explore, {
    removeModal(session)
    gd$x$stepwise <- FALSE
    disable("merge")
    disable("remove")
    disable("attach")
    disable("detach")
    shinyBS::updateCollapse(session, "sidebarpanel",
                            open = c("choose_priors"),
                            close = c("choose_structure", "param_expressions", "prior_expressions"))
  })


  #### for the splits

  observeEvent(input$stepwise_wV, {
    stepwise_wV_guide()
  })

  observeEvent(input$start_guide_splits, {
    gd$x$which_split <- 1
    gd$x$next_modal <- 1
    gd$x$stepwise <- TRUE
  })

  observe({

    if (gd$x$guide_is_running && gd$x$which_split > 0){
      tmp <- split_part_guide(nd$x, pd$w, gd$x)
      pd$w <- tmp$prior_weight
      gd$x <- tmp$guide_data
    }

  })

  observeEvent(eval(parse(text = paste0("input$use_diri", gd$x$counter_next))), {
    tmp <- chosen_dirichlet_guide(nd$x, pd$w, gd$x)
    pd$w <- tmp$prior_weight
    gd$x <- tmp$guide_data
  })

  observeEvent(eval(parse(text = paste0("input$use_pc_w", gd$x$counter_next))), {
    tmp <- chosen_pc_guide(nd$x, pd$w, gd$x)
    gd$x <- tmp$guide_data
    pd$w <- tmp$prior_weight
  })

  observeEvent(eval(parse(text = paste0("input$no_base", gd$x$counter_next))), {
    tmp <- median_pc_guide(nd$x, pd$w, gd$x, input, 0.25) # ask about median value
    pd$w <- tmp$prior_weight
    gd$x <- tmp$guide_data
  })

  observeEvent(eval(parse(text = paste0("input$yes_base1", gd$x$counter_next))), {
    tmp <- median_pc_guide(nd$x, pd$w, gd$x, input, 1) # ask about median value
    pd$w <- tmp$prior_weight
    gd$x <- tmp$guide_data
  })

  observeEvent(eval(parse(text = paste0("input$yes_base2", gd$x$counter_next))), {
    tmp <- median_pc_guide(nd$x, pd$w, gd$x, input, 0) # ask about median value
    pd$w <- tmp$prior_weight
    gd$x <- tmp$guide_data
  })

  observeEvent(eval(parse(text = paste0("input$set_median", gd$x$counter_next))), {
    tmp <- concentration_pc_guide(nd$x, pd$w, gd$x, input) # ask about concentration parameter
    pd$w <- tmp$prior_weight
    gd$x <- tmp$guide_data
  })

  observeEvent(eval(parse(text = paste0("input$set_conc_param", gd$x$counter_next))), {
    # inform about the chosen prior and go to next split/totvar
    tmp <- chosen_pc_finished(nd$x, pd$w, gd$x, input)
    pd$w <- tmp$prior_weight
  })

  observeEvent(eval(parse(text = paste0("input$set_median_01", gd$x$counter_next))), {
    # inform about the chosen prior and go to next split/totvar
    tmp <- chosen_pc_finished(nd$x, pd$w, gd$x, input)
    pd$w <- tmp$prior_weight
  })

  observeEvent(eval(parse(text = paste0("input$next_split", gd$x$counter_next))), {
    gd$x$which_split <- gd$x$which_split + 1
    gd$x$counter_next <- gd$x$counter_next + 1
    gd$x$next_modal <- 1
    # print(gd$x$counter_next)
  })



  #### for the total variance(s)

  observeEvent(input$go_to_step_2_V, {
    gd$x$which_totvar <- 1
    gd$x$which_split <- 0
    gd$x$which_cw <- 0
    gd$x$next_modal <- 0
    # gd$x$pr_ind <- 1
  })

  observe({
    if (gd$x$guide_is_running && gd$x$which_totvar > 0){
      tmp <- totvar_part_guide(nd$x, pd$V, gd$x, shiny::getShinyOption(".initial_args"))
      gd$x <- tmp$guide_data
    }
  })

  observeEvent(eval(parse(text = paste0("input$next_totvar", gd$x$counter_next))), {

    pd$V[[gd$x$pr_ind]]$prior <- input$var_prior_guide
    pd$V[[gd$x$pr_ind]]$param <- c(input$par1_guide, input$par2_guide)
    gd$x$which_totvar <- gd$x$which_totvar + 1
    gd$x$counter_next <- gd$x$counter_next + 1

  })



  #### for the independent variance(s)

  observeEvent(input$stepwise_cw, {
    gd$x$stepwise <- TRUE
    gd$x$which_totvar <- 0
    gd$x$which_cw <- 1
  })

  observe({

    if (gd$x$guide_is_running && gd$x$which_cw > 0){
      tmp <- cw_part_guide(nd$x, pd$var, gd$x)
      gd$x <- tmp$guide_data
    }

  })

  observeEvent(eval(parse(text = paste0("input$next_cw", gd$x$counter_next))), {

    pd$var[[gd$x$pr_ind]]$prior <- input$var_prior_guide
    pd$var[[gd$x$pr_ind]]$param <- c(input$par1_guide, input$par2_guide)

    gd$x$counter_next <- gd$x$counter_next + 1
    gd$x$which_cw <- gd$x$which_cw + 1

  })

  output$par1_box_guide <- make_par1_box_guide("par1_guide", input)
  output$par2_box_guide <- make_par2_box_guide("par2_guide", input)

  output$all_plots_guide <- renderPlot({
    if (gd$x$guide_is_running) {
      if (!shiny::getShinyOption(".initial_args")$.no_pc) {
        tmp_w <- calculate_pc_prior(nd$x, pd$w, shiny::getShinyOption(".initial_args")$.indata$covmats, gui = TRUE)
      } else {
        tmp_w <- list()
      }
      plot_priors(pd$w, pd$V, pd$var, tmp_w, nd$x, shiny::getShinyOption(".initial_args")$.no_pc)
    }
  },
  height = function() 150*ceiling((length(shiny::getShinyOption(".initial_args")$.nodenames))/3),
  bg = "transparent")

  observeEvent(input$go_to_step_3, {
    gd$x$finished <- TRUE
  })
  observeEvent(input$go_to_step_3_1, {
    pd$V[[gd$x$pr_ind]]$prior <- input$var_prior_guide
    pd$V[[gd$x$pr_ind]]$param <- c(input$par1_guide, input$par2_guide)
    gd$x$finished <- TRUE
  })
  observeEvent(input$go_to_step_3_2, {
    pd$var[[gd$x$pr_ind]]$prior <- input$var_prior_guide
    pd$var[[gd$x$pr_ind]]$param <- c(input$par1_guide, input$par2_guide)
    gd$x$finished <- TRUE
  })

  observe({
    if (gd$x$guide_is_running && gd$x$finished) {
      tmp <- finishing_guide(gd$x, session)
      gd$x <- tmp$guide_data
      session <- tmp$session
    }
  })


  # update the message to the user
  observe({
    if (!is.null(gd$x$step1) && !is.null(gd$x$step2)){
      if (gd$x$step1){
        gd$x <- guide_make_message_to_user_step1(gd$x, nd$x, input)
      } else if (gd$x$step2 && !gd$x$stepwise){
        gd$x <- guide_make_message_to_user_step2(gd$x, nd$x, list(w = pd$w, V = pd$V, var = pd$var), input)
      } else {
        gd$x$message_to_user <- NULL
      }
    } else {
      gd$x$message_to_user <- NULL
    }
  })

  # when guide is ended from modal window
  observeEvent(input$quit_guide, {
    removeModal(session)
    shinyBS::updateCollapse(session, "sidebarpanel",
                   open = c("choose_structure", "choose_priors"),#, "latest"),#, "closing"), # latest actions and closing should always be visible
                   close = c("param_expressions", "prior_expressions"))
    tmp <- gd$x$counter_next
    gd$x <- initialize_guide(FALSE)
    gd$x$counter_next <- tmp
  })
  # when guide is ended from main window
  observeEvent(input$quit_guide2, {
    removeModal(session)
    shinyBS::updateCollapse(session, "sidebarpanel",
                   open = c("choose_structure", "choose_priors"),#, "latest"),#, "closing"), # latest actions and closing should always be visible
                   close = c("param_expressions", "prior_expressions"))
    tmp <- gd$x$counter_next
    gd$x <- initialize_guide(FALSE)
    gd$x$counter_next <- tmp
  })

  # every time the node object is updated, we reset the message to the user
  observeEvent(nd$x, {
    if (gd$x$guide_is_running){
      gd$x$message_to_user <- "Click a node! To choose more than one node, you can use a long-click (click and hold for a bit). "
    }
  })

  observe({
    if (gd$x$guide_is_running) showElement("guide_message") else hideElement("guide_message")
  })

  output$guide_message <- renderText({
    gd$x$message_to_user
  })

  # the prior distributions, on the parameter scale that is used
  output$prior_dists_guide <- renderUI({
    withMathJax(helpText(
      make_pr_dists_latex(nd$x, pd$w, pd$V, pd$var, NULL)
    ))
  })


  ##################################################### outputs

  output$basemodel_name <- renderText({

    if (is_split_node(nd$x, input)){
      prior_info <- get_prior_info(pd$w, nd$x, input)
      if (prior_info$prior != "pc") {
        res <- paste("1/", prior_info$no_children, "variance to all")
      } else {
        basemod_id <- prior_info$param$above_node
        res <- paste(prior_info$param$basemodel, "variance to", get_node_name(nd$x, basemod_id))
      }
    } else {
      res <- "\n" #No node chosen"
    }

    res

  })

  isValid_median <- reactive({
    req(input$median)
    !is.null(input$median) && input$median >= 0 && input$median <= 1
  })

  isValid_conc_param <- reactive({
    req(input$conc_param)
    !is.null(input$conc_param) && input$conc_param > 0.5 && input$conc_param < 1
  })

  observe({
    if (!isValid_median()){
      updateNumericInput(session, "median", value = 0.5)
    }
  })

  observe({
    if (!isValid_conc_param()){
      updateNumericInput(session, "conc_param", value = 0.5)
    }
  })

  # output$basemodel_choice <- renderUI({
  #   if (isValid_median() && is_split_node(nd$x, input)){
  #     browser()
  #     basemodel_choices(nd$x, input, pd$w)
  #   }
  # })

  output$variance_prior_name <- renderText({

    if (length(input$current_node_id) == 1){
      if (are_detached(nd$x, input)) {
        res <- "CW variance prior"
      } else if (is_top_node(nd$x, input)){
        res <- "Total variance prior"
      } else{
        res <- "No node chosen"
      }
    } else {
      res <- "No node chosen"
    }

    res

  })

  output$graph <- visNetwork::renderVisNetwork({
    make_graph(nd$x)
  })

  # add the network to step 2 (prior choices) of the guide to clearify which nodes are getting a prior in the current step
  output$graph_guide <- visNetwork::renderVisNetwork({
    make_graph_guide(nd$x, pd$w, pd$V, pd$var, gd$x)
  })

  bottom_part <- function(){
    v <- list()
    if (gd$x$guide_is_running){
      if (gd$x$step1){
        v[[1]] <- actionButton("go_to_step_2", "Next step", class = "button1", title = "Next step in choosing prior distributions for your model.")
      } else if (gd$x$step2){
        v[[1]] <- actionButton("go_to_step_3", "Next step", class = "button1", title = "Next step in choosing prior distributions for your model.")
      }
      v[[2]] <- actionButton("quit_guide2", "Quit guide", class = "button_quit", title = "Quit the guide.")
    } else {
      v[[1]] <- actionButton("make_plot", "Plot priors", class = "button1", value = TRUE, title = "Plot the prior distributions you have chosen.")
    }
    return(v)
  }

  plot_priors_pushed <- eventReactive(input$make_plot, {
    tmp_w <- if (!shiny::getShinyOption(".initial_args")$.no_pc) calculate_pc_prior(nd$x, pd$w, shiny::getShinyOption(".initial_args")$.indata$covmats, gui = TRUE) else list()
    plot_priors(pd$w, pd$V, pd$var, tmp_w, nd$x, shiny::getShinyOption(".initial_args")$.no_pc)
    # plot_priors(pd$w, pd$V, pd$var, calculate_pc_prior(nd$x, pd$w, .initial_args$.indata$covmats), nd$x)
  })

  output$plot_or_next_step_button <- renderUI({
    bottom_part()
  })

  output$all_plots <- renderPlot({
    if (!gd$x$guide_is_running) plot_priors_pushed()
  },
  height = function() 150*ceiling((length(shiny::getShinyOption(".initial_args")$.nodenames))/3),
  bg = "transparent")

  # the variance equation, showing how the variance is distributed
  output$param_eq <- renderUI({
    withMathJax(helpText(
      text_latex_param(nd$x, pd$w, pd$V, pd$var, shiny::getShinyOption(".initial_args")$.response_name, input$param_type_eq)
    ))
  })


  # the prior distributions, on the parameter scale that is used
  output$prior_distributions <- renderUI({
    withMathJax(helpText(
      make_pr_dists_latex(nd$x, pd$w, pd$V, pd$var, input$current_node_id)
    ))
  })

  # model equation
  output$model_eq <- renderUI({
    withMathJax(helpText(
      make_model_eq_latex(shiny::getShinyOption(".initial_args"))
    ))
  })

  output$no_pc <- renderText({
    t1 <- ""
    if (shiny::getShinyOption(".initial_args")$.no_pc)
      t1 <- paste0(t1, "\nWill not compute the PC priors on splits before the app closes.")
    t1
  })

  output$intrinsic <- renderText({
    paste0("Note that the base model may have a singular covariance matrix. ",
           "Then nothing happens if the median is changed to be further than 0.25 away from the base model.")
  })

  # resets to input-values
  observeEvent(input$reset, {
    # reset all inputs
    reset("basemodel")
    reset("median")
    reset("conc_param")
    reset("par1")
    reset("par2")
    reset("weight_prior")
    reset("var_prior")
    # reset objects
    pd$w <- shiny::getShinyOption(".initial_args")$.initial_prior_w
    pd$V <- shiny::getShinyOption(".initial_args")$.initial_prior_V
    pd$var <- shiny::getShinyOption(".initial_args")$.initial_prior_var
    nd$x <- initialize_nodes(shiny::getShinyOption(".initial_args"))
    gd$x <- initialize_guide(FALSE)
  })

  # # closes application, sends output back to R and starting inference immediately
  # observeEvent(input$close_and_run, {stopApp(make_output(
  #   input,
  #   isolate(list(weights = pd$w, total_variance = pd$V, cw_priors = pd$var)),
  #   isolate(nd$x),
  #   .initial_args$.indata,
  #   run_inference_directly = TRUE
  # ))}) # with button

  # closes application and sends output back to R
  observeEvent(input$close, {stopApp(make_output(
    input,
    isolate(list(weights = pd$w, total_variance = pd$V, cw_priors = pd$var)),
    isolate(nd$x),
    shiny::getShinyOption(".initial_args")$.indata
    ))}) # with button
  session$onSessionEnded(function() {stopApp(make_output(
    input,
    isolate(list(weights = pd$w, total_variance = pd$V, cw_priors = pd$var)),
    isolate(nd$x),
    shiny::getShinyOption(".initial_args")$.indata
  ))}) # closing window

}


