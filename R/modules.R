



# UI-things

merge_button <- function(){
  actionButton("merge", "Merge", style = "width:80px", class = "button1",
               title = "Make a new split node with the selected nodes as children. You must select at least two nodes on the same level. ")
}

remove_button <- function(){
  actionButton("remove", "Remove", style = "width:80px", class = "button1",
               title = "Remove a split node. The parent of the removed node will become the parent of the children of the removed node. You must select one node (that is not a top node) with at least two children.")
}

detach_button <- function(){
  actionButton("detach", "Detach", style = "width:80px", class = "button1",
               title = "Detach a node from the tree, to give it an independent prior. You must select one node, that is attached to a tree.")
}

attach_button <- function(){
  actionButton("attach", "Attach", style = "width:80px", class = "button1",
               title = "Attach the detached node(s) to the node tree. You must select a node in the tree, and at least one detached node.")
}



# server-things

# modify tree structure buttons

merge_action <- function(node_data, prior_weight, prior_totvar, prior_cw, text_data, input){

    disable("merge") # disable when clicked
    disable("weight_prior")
    text_data <- paste0("Merged ", paste0(get_node_name_click(node_data, input), sep = "", collapse = " and "))
    if (are_detached(node_data, input)){
      node_data <- merge_detached_nodes(node_data, input)
      prior_totvar <- update_V_prior_merged_detached_nodes(prior_totvar, node_data, input)
      prior_cw <- update_variance_prior_attach_detach(prior_cw, node_data, "merge", input)
    } else {
      if (sum(node_data$nodes$top_node[node_data$nodes$id %in% input$current_node_id]) == length(input$current_node_id)){
        node_data <- merge_trees(node_data, input)
        prior_totvar <- update_V_prior_merged_trees(prior_totvar, node_data, input)
      } else {
        node_data <- update_nodes_edges_merge(node_data, input)
      }
    }
    prior_weight <- update_prior_w_merge_remove(prior_weight, node_data, "merge", input)
    node_data <- update_basemodel_edges(node_data, prior_weight)

    return(list(node_data = node_data, prior_weight = prior_weight, prior_totvar = prior_totvar, prior_cw = prior_cw, text_data = text_data))

}

remove_action <- function(node_data, prior_weight, prior_totvar, prior_cw, text_data, input){

  disable("remove") # disable when clicked
  disable("weight_prior")
  text_data <- paste0("Removed ", get_node_name_click(node_data, input))
  node_data <- update_nodes_edges_remove(node_data, input)
  prior_weight<- update_prior_w_merge_remove(prior_weight, node_data, "remove", input)
  node_data <- update_basemodel_edges(node_data, prior_weight)

  return(list(node_data = node_data, prior_weight = prior_weight, prior_totvar = prior_totvar, prior_cw = prior_cw, text_data = text_data))

}

detach_action <- function(node_data, prior_weight, prior_totvar, prior_cw, text_data, input){

  disable("detach") # disable when clicked
  tmp_sibling <- get_siblings(node_data, input$current_node_id)
  tmp_sibling <- tmp_sibling[!(tmp_sibling %in% input$current_node_id)]
  node_data <- detach_node(node_data, input)
  node_data <- update_node_labels(node_data)
  prior_weight <- update_prior_w_attach_detach(prior_weight, node_data, "detach", input)
  prior_cw <- update_variance_prior_attach_detach(prior_cw, node_data, "detach", input)
  prior_totvar <- update_V_prior_attach_detach(prior_totvar, node_data, "detach", tmp_sibling)
  text_data <- paste0("Detached ", get_node_name_click(node_data, input))
  rm(tmp_sibling)

  return(list(node_data = node_data, prior_weight = prior_weight, prior_totvar = prior_totvar, prior_cw = prior_cw, text_data = text_data))

}

attach_action <- function(node_data, prior_weight, prior_totvar, prior_cw, text_data, input){

  disable("attach") # disable when clicked
  node_data <- attach_nodes(node_data, input)
  node_data <- update_node_labels(node_data)
  prior_weight <- update_prior_w_attach_detach(prior_weight, node_data, "attach", input)
  prior_cw <- update_variance_prior_attach_detach(prior_cw, node_data, "attach", input)
  prior_totvar <- update_V_prior_attach_detach(prior_totvar, node_data, "attach")
  text_data <- paste0("Attached ", get_node_name_click(node_data, input))

  return(list(node_data = node_data, prior_weight = prior_weight, prior_totvar = prior_totvar, prior_cw = prior_cw, text_data = text_data))

}


# clicked node action

clicked_node_action <- function(node_data, prior_weight, prior_totvar, prior_cw, text_data, input, session){

  if (length(input$current_node_id) == 1 && is_split_node(node_data, input)) {
    enable("weight_prior")
  } else {
    disable("weight_prior")
    hideElement("basemodel")
    hideElement("median")
    hideElement("conc_param")
  }
  if (!is_dual_split_node(node_data, input)) disable(selector = "[type=radio][value=pc]") # only have the option to use pc prior on dual split

  # remove or include components in the tree, so the user can give CW priors to some components
  if (is_leaf_node(node_data, input) && are_attached(node_data, input)){ # can only detach one node at a time
    enable("detach")
    disable("attach")
  } else if (can_be_attached(node_data, input)){ # can attach several nodes at the same time
    enable("attach")
    disable("detach")
  } else {
    disable("attach")
    disable("detach")
  }

  # everything that is required for merging/removing nodes
  if (enable_merge(node_data, input)) enable("merge") else disable("merge")
  if (enable_remove(node_data, input)) enable("remove") else disable("remove")

  # only have the possibility to choose total variance prior if we click a top-node
  if (length(input$current_node_id) == 1 && (is_top_node(node_data, input) || are_detached(node_data, input))) {
    enable("var_prior")
    showElement("par1")
    showElement("par2")
  } else {
    disable("var_prior")
    hideElement("par1")
    hideElement("par2")
  }

  # only allow jeffreys when we have only one tree, and no cw priors, else the user must use a proper prior
  if (sum(node_data$nodes$status == "detached") > 0 || sum(node_data$nodes$top_node) != 1) {
    disable(selector = "[type=radio][value=jeffreys]")
  }

  # remember prior choices for splits
  if (is_split_node(node_data, input)){
    prior_info <- get_prior_info(prior_weight, node_data, input)
    updateRadioButtons(session, "weight_prior", selected = prior_info$prior)
    children <- get_children(node_data, input)
    if (prior_info$prior == "pc"){
      showElement("median")
      if (input$basemodel == "comb") showElement("conc_param")
      showElement("basemodel")
      updateRadioButtons(session, "basemodel",
                         # choiceNames = c(get_node_name(node_data, children[1]), get_node_name(node_data, children[2]),
                         #                 paste0(round(input$median, 3), " ", get_node_name(node_data, children[1]),
                         #                        " and ",
                         #                        round(1 - input$median, 3), " ", get_node_name(node_data, children[2]))
                         # ),
                         choiceNames = c("100% A + 0% B", "0% A + 100% B",
                                         paste0(round(prior_info$param$median, 2)*100, "% A + ", round(1-prior_info$param$median, 2)*100, "% B")),
                         choiceValues = c("mod1", "mod2", "comb"),
                         selected = if (prior_info$param$basemodel == 1) "mod1" else if (prior_info$param$basemodel == 0) "mod2" else "comb"
      )
      updateNumericInput(session, "median", value = prior_info$param$median)
      updateNumericInput(session, "conc_param", value = prior_info$param$concentration)
    }
  }

  # remember prior choices for total variance
  if (is_top_node(node_data, input)){
    prior_info <- get_prior_info(prior_totvar, node_data, input)
    updateRadioButtons(session, "var_prior", selected = prior_info$prior)
    updateNumericInput(session, "par1", value = prior_info$param[1], label = get_par1_label(prior_info$prior))
    updateNumericInput(session, "par2", value = prior_info$param[2], label = get_par2_label(prior_info$prior),
                       max = get_par2_max(prior_info$prior), step = get_par2_step(prior_info$prior))
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
  }

  # remember prior choices for CW priors
  if (length(input$current_node_id) == 1 && are_detached(node_data, input)){
    prior_info <- get_prior_info(prior_cw, node_data, input)
    updateRadioButtons(session, "var_prior", selected = prior_info$prior)
    showElement("par1")
    showElement("par2")
    updateNumericInput(session, "par1", value = prior_info$param[1], label = get_par1_label(prior_info$prior))
    updateNumericInput(session, "par2", value = prior_info$param[2], label = get_par2_label(prior_info$prior),
                       max = get_par2_max(prior_info$prior), step = get_par2_step(prior_info$prior))
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
  }

  #return(list(session = session))

  # return(list(node_data = node_data, prior_weight = prior_weight, prior_totvar = prior_totvar, prior_cw = prior_cw, text_data = text_data))

}

get_par2_max <- function(prior_name) if (prior_name == "pc0") 1 else Inf
get_par2_step <- function(prior_name) if (prior_name == "pc0") 0.05 else 0.01

get_par1_label <- function(prior_name) if (prior_name == "pc0") "U" else if (prior_name == "invgam") "Shape" else if (prior_name %in% c("hc", "hn")) "Scale"
get_par2_label <- function(prior_name) if (prior_name == "pc0") "alpha" else if (prior_name == "invgam") "Scale"

get_par1_default_val <- function(prior_name) if (prior_name == "pc0") 3 else if (prior_name == "invgam") 1 else if (prior_name == "hc") 25 else if (prior_name == "hn") 1
get_par2_default_val <- function(prior_name) if (prior_name == "pc0") 0.05 else if (prior_name == "invgam") 5e-5


disable_everything <- function(){

  disable("merge")
  disable("remove")
  disable("attach")
  disable("detach")
  disable("weight_prior")
  disable("var_prior")
  hideElement("par1")
  hideElement("par2")
  hideElement("basemodel")
  hideElement("median")
  hideElement("conc_param")
  hideElement("intrinsic")

}


make_graph <- function(node_data){

  nd <- list(x = node_data)

  tmp_nodes <- nd$x$nodes
  tmp_nodes$color.background <- node_palette(calc_variance_proportions(nd$x)$nodes$varprop)
  tmp_nodes$color.background[tmp_nodes$status == "detached"] <- detached_node_color
  #tmp_nodes$color.border <- node_border_color
  tmp_nodes$color.highlight.background <- tmp_nodes$color.background
  #tmp_nodes <- set_node_title(tmp_nodes, nd$x)

  tmp_edges <- nd$x$edges
  tmp_edges$width <- tmp_edges$width*5 # to get a bit thicker edges in the graph
  if (nrow(tmp_edges) > 0) tmp_edges <- cbind(tmp_edges, arrows = "to", hidden = FALSE)

  tmp_edges <- set_detached_edges(tmp_nodes, tmp_edges)

  nn <- visNetwork::visNetwork(tmp_nodes, tmp_edges)
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
                 font = list(background = "#EAE7DC", color = node_border_color))
  nn <- visNetwork::visOptions(nn, highlightNearest = FALSE)
  nn <- visNetwork::visInteraction(nn,
                       dragNodes = FALSE,
                       dragView = FALSE,
                       zoomView = FALSE,
                       multiselect = TRUE)
  nn <- visNetwork::visEvents(nn, select = "function(nodes) {
                    Shiny.onInputChange('current_node_id', nodes.nodes);
                    ;}")

  return(nn)

}


