


### plotting functions

plot_priors <- function(prior_w, prior_tot, prior_cw, prior_weights, node_data, no_pc){

  df_pri <- data.frame()

  # for CW priors:
  for (i in 1:length(prior_cw)){
    if (prior_cw[[i]]$prior != "") {
      tmp_pr <- get_prior_expr(prior_cw[[i]], node_data, "cw")
    }
  }

  # for CW priors:
  for (i in 1:length(prior_cw)){
    tmp <- df_prior_variance(prior_cw[[i]], node_data, "cw")
    df_pri <- rbind(df_pri, tmp)
  }

  # for total variances
  for (i in seq_len(length(prior_tot))){
    tmp <- df_prior_variance(prior_tot[[i]], node_data, "total")
    df_pri <- rbind(df_pri, tmp)
  }

  # for HD priors
  for (i in seq_len(length(prior_w))){
    if (prior_w[[i]]$prior == "pc"){
      if (no_pc){ # if the PC priors should not be computed inside the app
        tmp <- plot_pc_w_nocalc(prior_w[[i]], node_data)
      } else {
        tmp <- plot_pc_w(prior_w[[i]], prior_weights[[i]], node_data)
      }
    } else { # dirichlet
      tmp <- plot_diri_w(prior_w[[i]], node_data)
    }
    df_pri <- rbind(df_pri, tmp)
  }

  gg <- ggplot(df_pri, aes(x = .data$x, y = .data$y)) + geom_line(na.rm = TRUE) +
    facet_wrap(~param, labeller = label_parsed, scales = "free", ncol = 3) +
    theme(strip.text = element_text(size = 14)) +
    #ylab("Density") +
    ylim(0, NA) +
    theme(axis.title = element_blank()) +
    theme(plot.background = element_rect(fill = "transparent", color = NA)) +
    theme(panel.background = element_rect(fill = "white", color = gray(0.5)),
          panel.grid.major = element_line(color = gray(0.85)),
          panel.grid.minor = element_line(color = gray(0.92))) +
    theme(strip.background = element_rect(fill = plot_color, color = gray(0.5)),
          strip.text = element_text(color = plot_text_color, size = 15))

  return(gg)

}


### functions to classify nodes

# checks if clicked node is leaf_node
# if the node has no children, i.e. there are no edges from it, it is a leaf node
# it is NOT a leaf node if it has no edges either to or form it (it has an individual prior and is not in a tree at all)
is_leaf_node <- function(node_data, input) {
  if (length(input$current_node_id) != 1) return(FALSE)
  if (input$current_node_id %in% node_data$edges$from || !(input$current_node_id %in% c(node_data$edges$from, node_data$edges$to))) return(FALSE) else return(TRUE)
}

# checks if several nodes are leaf nodes or not (if one is not leaf node, return FALSE)
are_leaf_nodes <- function(node_data, input){
  if (length(input$current_node_id) == 0) return(FALSE)
  tmp <- sapply(input$current_node_id, function(x) !(x %in% node_data$edges$from))
  if (sum(tmp) != length(tmp)) return(FALSE) else return(TRUE)
}

# checks if clicked node is a node at the top level (can be several top nodes, if we have more than one tree)
# if the node is not a child, i.e. there are no edges to it, it is a top node
# it is not a top node if there are no edges either from or to it (then it has a CW prior)
is_top_node <- function(node_data, input) {
  if (length(input$current_node_id) != 1) return(FALSE)
  #if ((input$current_node_id %in% node_data$edges$to) && node_data$nodes$top_node[node_data$nodes$id == input$current_node_id] == 0) return(FALSE) else return(TRUE)
  if ((node_data$nodes$top_node[node_data$nodes$id == input$current_node_id] == 0) ||
      !(input$current_node_id %in% c(node_data$edges$from, node_data$edges_to))) return(FALSE) else return(TRUE)
}

# counts number of child nodes
no_of_children <- function(node_data, input) {
  if (length(input$current_node_id) != 1) return(FALSE)
  return(sum(node_data$edges$from == input$current_node_id))
}

# returns the id of the child nodes of a given node
get_children <- function(node_data, input){
  if (is_leaf_node(node_data, input)) return(NULL)
  return(node_data$edges$to[node_data$edges$from == input$current_node_id])
}

# tests if the node is a split node
is_split_node <- function(node_data, input) {
  if (length(input$current_node_id) != 1) return(FALSE)
  if (no_of_children(node_data, input) > 1) return(TRUE) else return(FALSE)
}

# tests if the node is a dual split node
is_dual_split_node <- function(node_data, input) {
  if (length(input$current_node_id) != 1) return(FALSE)
  if (no_of_children(node_data, input) == 2) return(TRUE) else return(FALSE)
}

# get name of a node for the node clicked, function for given node id is in another file (called get_node_name)
get_node_name_click <- function(node_data, input) node_data$nodes$label[node_data$nodes$id %in% input$current_node_id]

# get the model type for the node clicked
get_node_model_name_click <- function(node_data, model_data, input){
  sapply(sapply(input$current_node_id, function(y) get_prior_number(model_data, y)),
         function(x) {
           if (is_leaf_node(node_data, list(current_node_id = x)) || are_detached(node_data, list(current_node_id = x))) {
             paste0(" (", model_data[[x]]$model, ")")
           } else {
             ""
           }
         })
}

# check if a node/several nodes is attached to a tree or not (is in a MW prior)
are_attached <- function(node_data, input){
  if (length(input$current_node_id) == 0) return(FALSE)
  if (sum(node_data$nodes$status[node_data$nodes$id %in% input$current_node_id] == "attached") == length(input$current_node_id)) return(TRUE) else return(FALSE)
}

# check if a node/several nodes is standing alone (has a CW prior)
are_detached <- function(node_data, input){
  if (length(input$current_node_id) == 0) return(FALSE)
  if (sum(node_data$nodes$status[node_data$nodes$id %in% input$current_node_id] == "detached") == length(input$current_node_id)) return(TRUE) else return(FALSE)
}

# check if we can attach the chosen nodes in the tree (user must have chosen a parent node for them to be attached)
can_be_attached <- function(node_data, input){
  if (length(input$current_node_id) <= 1) return(FALSE)
  correct_detached <- sum(node_data$nodes$status[node_data$nodes$id %in% input$current_node_id] == "detached") == length(input$current_node_id)-1
  correct_attached <- sum(node_data$nodes$status[node_data$nodes$id %in% input$current_node_id] == "attached") == 1
  if (correct_detached && correct_attached) return(TRUE)
  # if we do not have a tree anymore, we must merge to get a new tree, and not use the attach button
  if (nrow(node_data$edges) == 0) return(FALSE)
  return(FALSE)
}

# get the node_data-object for this tree (returns everything if only one tree)
get_node_data_for_tree <- function(node_data, node_id){

  node_data_this_tree <- node_data
  nodes_in_this_tree <- sort(nodes_in_tree(node_data, node_id))
  node_data_this_tree$nodes <- node_data$nodes[node_data$nodes$id %in% nodes_in_this_tree,]
  node_data_this_tree$edges <- node_data$edges[node_data$edges$from %in% nodes_in_this_tree,]

  return(node_data_this_tree)
}

# find the siblings of a node
get_siblings <- function(node_data, node_id) get_children(node_data, list(current_node_id = get_parent_node_id(node_data, node_id)))

# get the ids of all split nodes by level
get_split_ids_by_level <- function(node_data, decreasing){
  tmp <- node_data$nodes[node_data$nodes$id %in% get_split_ids(node_data),]
  tmp <- tmp[order(tmp$level, decreasing = decreasing),]
  return(tmp$id)
}

get_top_nodes <- function(node_data) node_data$nodes$id[node_data$nodes$top_node == TRUE]

get_detached_nodes <- function(node_data) node_data$nodes$id[node_data$nodes$status == "detached"]


### functions to update the node tree structure

update_nodes_edges_merge <- function(node_data, input){

  nodes_old <- node_data$nodes
  edges_old <- node_data$edges

  new_node_id <- max(nodes_old$id)+1

  # move all merged nodes one level down
  nodes_old$level[nodes_old$id %in% input$current_node_id] <- nodes_old$level[nodes_old$id %in% input$current_node_id] + 1

  # if any of the merged nodes are split nodes, we must move all their children as well
  split_nodes <- input$current_node_id[input$current_node_id %in% node_data$edges$from]
  if (length(split_nodes) > 0){
    for (id in split_nodes){
      child_nodes <- edges_old$to[edges_old$from == id]
      nodes_old$level[nodes_old$id %in% child_nodes] <- nodes_old$level[nodes_old$id %in% child_nodes] + 1
    }
  }

  node_data$nodes <- rbind(nodes_old,
                           data.frame(
                             id = new_node_id,
                             label = paste(nodes_old$label[nodes_old$id %in% input$current_node_id], collapse = "_"),
                             level = nodes_old$level[nodes_old$id %in% input$current_node_id][1]-1,
                             status = "attached",
                             # y = NA,
                             top_node = 0
                           ))

  # remove old edges to the merged nodes

  edges_new <- edges_old[!(edges_old$to %in% input$current_node_id),]
  edges_new$width <- NA

  node_data$edges <- rbind(edges_new, # the old edges to keep
                           data.frame( # the edges to the new node
                             from = edges_old$from[(edges_old$to %in% input$current_node_id)][1],
                             to = new_node_id,
                             width = NA,
                             dashes = FALSE,
                             label = ""
                           ),
                           data.frame( # the edge from the new node to the ones that were merged
                             from = new_node_id,
                             to = input$current_node_id,
                             width = 1/length(input$current_node_id),
                             dashes = FALSE,
                             label = ""
                           )
  )

  node_data$edges$width[is.na(node_data$edges$width)] <- 1/(sum(is.na(node_data$edges$width)))

  return(node_data)

}

update_nodes_edges_remove <- function(node_data, input){

  nodes_old <- node_data$nodes
  edges_old <- node_data$edges

  node_to_remove <- nodes_old[nodes_old$id == input$current_node_id,]

  lvls <- range(nodes_old$level)

  if (is_leaf_node(node_data, input)) { # if the node to remove is on the bottom level, it should not be removed (but this should not be a possible action)
    print("Not possible")
  } else if (node_to_remove$level == min(lvls)) { # if the node is at the top level, nothing should happen/the button should not be usable
    print("Not possible")
  } else { # move all children of the node to be removed up to the level above and move edges from removed node

    nodes_new <- nodes_old[!(nodes_old$id == input$current_node_id),]
    edges_new <- edges_old[!(edges_old$to == input$current_node_id),]
    # the edges from the removed node must be moved to go from the parent of the removed node
    edges_new$from[edges_new$from == input$current_node_id] <- edges_old$from[edges_old$to == input$current_node_id]
    # move all nodes below this (that has this node as some parent) one level up
    nodes_below <- get_nodes_below(node_data, input$current_node_id)
    #nodes_new$level[edges_old$to[(edges_old$from == input$current_node_id)]] <- nodes_old$level[input$current_node_id]
    nodes_new$level[nodes_new$id %in% nodes_below] <- nodes_new$level[nodes_new$id %in% nodes_below] - 1

    # # if the node that is removed is not the one with the highest id, change the id of all nodes with higher id
    # if (nrow(nodes_new) < max(nodes_new$id)){
    #   # change id of all nodes with higher id than the one that was removed
    #   nodes_new$id[nodes_new$id > node_to_remove$id] <- nodes_new$id[nodes_new$id > node_to_remove$id] -1
    #   edges_new$from[edges_new$from > node_to_remove$id] <- edges_new$from[edges_new$from > node_to_remove$id] -1
    #   edges_new$to[edges_new$to > node_to_remove$id] <- edges_new$to[edges_new$to > node_to_remove$id] -1
    # }

    node_data$nodes <- nodes_new
    node_data$edges <- edges_new

  }

  # # if the node that is removed is not the one with the highest id, change the id of all nodes with higher id
  # if (nrow(node_data$nodes) < max(node_data$nodes$id)){
  #   old_new_ids <- data.frame(old = sort(node_data$nodes$id), new = 1:nrow(node_data$nodes))
  #   node_data <- update_node_ids(node_data, old_new_ids)
  # }

  return(node_data)

}

merge_detached_nodes <- function(node_data, input){

  node_ids <- input$current_node_id

  # make a new top node for a new tree
  new_top_node <- data.frame(
    id = max(node_data$nodes$id)+1,
    label = paste(get_node_name(node_data, node_ids), collapse = "_"),
    level = 0,
    status = "attached",
    top_node = 1
  )

  # make new edges from new top node to the merged nodes
  new_edges <- data.frame(
    from = new_top_node$id,
    to = node_ids,
    width = 1/length(node_ids),
    dashes = FALSE,
    label = ""
  )

  node_data$nodes$status[node_data$nodes$id %in% node_ids] <- "attached"
  node_data$nodes$level[node_data$nodes$id %in% node_ids] <- 1

  # update the edges of the detached nodes
  node_data$nodes$level[node_data$nodes$status == "detached"] <- (1:sum(node_data$nodes$status == "detached")-1)%%4

  node_data$nodes <- rbind(node_data$nodes, new_top_node)
  node_data$edges <- rbind(node_data$edges, new_edges)

  return(node_data)

}

merge_trees <- function(node_data, input){

  node_ids <- input$current_node_id

  # make a new top node for a new tree
  new_top_node <- data.frame(
    id = max(node_data$nodes$id)+1,
    label = paste(get_node_name(node_data, node_ids), collapse = "_"),
    level = 0,
    status = "attached",
    top_node = 1
  )

  # make new edges from new top node to the merged nodes
  new_edges <- data.frame(
    from = new_top_node$id,
    to = node_ids,
    width = 1/length(node_ids),
    dashes = FALSE,
    label = ""
  )

  # move all nodes of these trees one level down

  involved_nodes <- unlist(lapply(node_ids, function(x) nodes_in_tree(node_data, x)))

  node_data$nodes$level[node_data$nodes$id %in% involved_nodes] <- node_data$nodes$level[node_data$nodes$id %in% involved_nodes] +1
  node_data$nodes$top_node[node_data$nodes$id %in% involved_nodes] <- 0 # none of the old top nodes are top nodes anymore

  node_data$nodes <- rbind(node_data$nodes, new_top_node)
  node_data$edges <- rbind(node_data$edges, new_edges)

  return(node_data)

}

# for detaching a node from the tree to give them CW priors, can only do this to leaf nodes (button should not work if not only leaf nodes)
# TODO: can only do this for one node at a time now, consider opening for removing more nodes
detach_node <- function(node_data, input){

  old_edges <- node_data$edges
  old_nodes <- node_data$nodes

  # we must check if the removal of the chosen nodes from the tree makes a split be superflous (includes the top split)
  # for each node that will be removed, find the parent node, and see if this parent end up with 1 og 0 nodes
  parent_node <- unique(sapply(input$current_node_id, function(x) get_parent_node_id(node_data, x)))
  # if this parent only has two children, and we are removing one of them, remove the parent and replace with the sibling node of the detached node
  # only do this if there will still be a tree after we remove the node
  nodes_in_this_tree <- nodes_in_tree(node_data, input$current_node_id)
  if (no_of_children(node_data, list(current_node_id = parent_node)) == 2 && length(nodes_in_this_tree) > 2){

    sibling <- get_children(node_data, list(current_node_id = parent_node))
    sibling <- sibling[sibling != input$current_node_id]

    # move edge from parent node that is removed to the sibling node
    old_edges$to[old_edges$to == parent_node] <- sibling

    # move the sibling one level up
    old_nodes$level[old_nodes$id == sibling] <- old_nodes$level[old_nodes$id == sibling] - 1

    # move all nodes below the sibling node in the tree one level up, if any
    nodes_below <- old_nodes$id %in% get_nodes_below(node_data, sibling)
    # if there are no nodes below, the node id is returned, and then we do not move the node up
    # if we have nodes below, we always have two or more
    if (sum(nodes_below) != 1){
      old_nodes$level[nodes_below] <- old_nodes$level[nodes_below] -1
    }

    # remove edges from parent node that is removed
    old_edges <- old_edges[!(old_edges$from == parent_node),]

    # if the parent node is the top node, we need a new top node
    # the new top node is the sibling
    if (is_top_node(node_data, list(current_node_id = parent_node))){
      old_nodes$top_node[old_nodes$id == sibling] <- 1
    }

    old_nodes <- old_nodes[!(old_nodes$id == parent_node),] # remove parent node (that was a split node)

  }

  old_nodes$status[old_nodes$id == input$current_node_id] <- "detached" # update status of detached node

  # if we are left with only one node (in addition to the top node) in this tree, we do not have a tree at all, and remove all edges
  # and set the status of the other node (which means all nodes) that is left to "detached"
  if (sum(old_nodes[old_nodes$id %in% nodes_in_this_tree,]$status == "attached") <= 2){
    old_edges <- old_edges[!(old_edges$from %in% nodes_in_this_tree),] # remove edges we no longer need
    #browser()
    old_nodes[old_nodes$id %in% nodes_in_this_tree,]$status <- "detached" # all nodes are now detached
    old_nodes[old_nodes$id %in% nodes_in_this_tree,]$top_node <- 0 # none of these nodes should be top nodes
    # set correct level of detached nodes, order them by the order of the original nodes
    # (all nodes in this tree should be detached now, but keep "== "detached"" for now)
    old_nodes[old_nodes$id %in% nodes_in_this_tree,]$level <- c(1:sum(old_nodes[old_nodes$id %in% nodes_in_this_tree,]$status == "detached")-1)%%4


  } else {

    # for each node to be attached, we must remove the edges to it, and change the level
    # (for now we put it under the tree, can maybe move it to the right or something later)
    # we must also set the status to "detached"
    for (node_id in input$current_node_id){ # for now we can only detach one node at a time, so this for loop is not really necessary

      old_edges <- old_edges[!(old_edges$to == node_id),]

      # new_level <- max(old_nodes$level[old_nodes$status == "attached"])+1
      # new_level <- max(c(old_nodes$level[old_nodes$status == "detached"], -1), na.rm = TRUE) +1
      # only have four detached nodes in each column, make more columns if necessary (this node is not marked as detached yet)

      old_nodes$status[old_nodes$id == node_id] <- "detached" # this is not necessary when only one node can be detached at a time
      new_level <- (sum(old_nodes$status == "detached")-1)%%4
      old_nodes$level[old_nodes$id == node_id] <- new_level

    }
  }

  node_data$nodes <- old_nodes
  node_data$edges <- old_edges

  # if all nodes are detached, fix levels
  if (sum(node_data$nodes$status == "attached") == 0){
    node_data$nodes$level <- node_data$nodes$id - 2
  }

  return(node_data)

}

# setting the edges of the detached nodes to arrange them in an orderly fashion, and does not change the "global" node-object
# returns the edges only
set_detached_edges <- function(nodes, edges){

  n_det <- sum(nodes$status == "detached")

  if (n_det < 2) return(edges)

  from_numbers <- c(1,2,3,1)+rep(seq(0,15,4), each = 4)
  to_numbers <- 2:16

  tmp_edges <- data.frame(
    from = nodes$id[nodes$status == "detached"][from_numbers[1:(n_det-1)]],
    to = nodes$id[nodes$status == "detached"][to_numbers[1:(n_det-1)]],
    width = 0,
    dashes = FALSE,
    label = "",
    arrows = "to",
    hidden = TRUE
  )

  tmp_edges <- rbind(edges, tmp_edges)

  return(tmp_edges)

}

# for attaching a detached node to another node
# if a leaf node is chosen, a new split with the detach nodes and this leaf node is made
# if a split node is chosen, the detached nodes are added to this split
attach_nodes <- function(node_data, input){

  #detached_nodes_id <- node_data$nodes$id[node_data$nodes$status == "detached"]
  detached_nodes_id <- node_data$nodes$id[node_data$nodes$status == "detached" & node_data$nodes$id %in% input$current_node_id]
  node_in_tree <- input$current_node_id[input$current_node_id %in% node_data$nodes$id[node_data$nodes$status == "attached"]]

  # if the node in the tree is a leaf node, make a new split with this node and the detached nodes
  if (is_leaf_node(node_data, list(current_node_id = node_in_tree)) && nrow(node_data$edges) > 0){

    nodes_old <- node_data$nodes
    edges_old <- node_data$edges

    new_node_id <- max(nodes_old$id)+1 # id of split node

    new_level <- nodes_old$level[nodes_old$id == node_in_tree] + 1

    # move all merged nodes one level down
    nodes_old$level[nodes_old$id %in% input$current_node_id] <- new_level

    # add new split node to tree
    node_data$nodes <- rbind(nodes_old,
                             data.frame(
                               id = new_node_id,
                               label = paste(nodes_old$label[nodes_old$id %in% input$current_node_id], collapse = "_"),
                               level = nodes_old$level[nodes_old$id %in% input$current_node_id][1]-1,
                               status = "attached",
                               #y = NA,
                               top_node = 0
                             ))

    # remove edge to the leaf node
    edges_new <- edges_old[!(edges_old$to %in% node_in_tree),]

    node_data$edges <- rbind(edges_new, # the old edges to keep
                             data.frame( # the edges to the new node
                               from = edges_old$from[(edges_old$to %in% input$current_node_id)][1],
                               to = new_node_id,
                               width = NA,
                               dashes = FALSE,
                               label = ""
                             ),
                             data.frame( # the edge from the new node to the ones that were merged
                               from = new_node_id,
                               to = input$current_node_id,
                               width = 1/length(input$current_node_id),
                               dashes = FALSE,
                               label = ""
                             )
    )

    node_data$edges$width[is.na(node_data$edges$width)] <- 1/(sum(is.na(node_data$edges$width)))

    node_data$nodes$status[node_data$nodes$id == detached_nodes_id] <- "attached"

  } else if (is_split_node(node_data, list(current_node_id = node_in_tree))){ # add to an existing split

    new_edges <- data.frame(from = node_in_tree,
                            to = detached_nodes_id,
                            #level = node_data$edges$level[node_data$edges$from == node_in_tree],
                            width = 1/(no_of_children(node_data, list(current_node_id = node_in_tree))+length(detached_nodes_id)),
                            dashes = FALSE,
                            label = "")
    node_data$edges <- rbind(node_data$edges, new_edges)

    node_data$nodes$level[node_data$nodes$id %in% detached_nodes_id] <- node_data$nodes$level[node_data$nodes$id == node_in_tree]+1

    node_data$nodes$status[node_data$nodes$id %in% detached_nodes_id] <- "attached"

  } else if (nrow(node_data$edges) == 0){ # if we have no tree at all, we can make one
    #browser()
    # make a new tree, the user can later make this be the data-node
    new_node_id <- if (1 %in% node_data$nodes$id) max(node_data$nodes$id)+1 else 1

    new_edges <- data.frame(from = new_node_id, to = detached_nodes_id, width = 1/length(detached_nodes_id), dashes = FALSE, label = "")

    new_top_node <- data.frame(id = new_node_id,
                               label = paste(node_data$nodes$label[node_data$nodes$id %in% detached_nodes_id], collapse = "_"),
                               level = 0,
                               status = "attached",
                               #y = NA,
                               top_node = 1)

    node_data$edges <- new_edges
    node_data$nodes <- rbind(new_top_node, node_data$nodes)

    node_data$nodes$status[node_data$nodes$id %in% detached_nodes_id] <- "attached"
    node_data$nodes$level[node_data$nodes$id %in% detached_nodes_id] <- 1

    print(node_data)

  } else {
    # browser()
    # stop("Should not be possible to end here!!!")
  }

  # if any detached nodes left, update their level
  if (sum(node_data$nodes$status == "detached") > 0){

    node_data$nodes$level[node_data$nodes$status == "detached"] <- c(1:sum(node_data$nodes$status == "detached")-1)%%4

  }

  return(node_data)

}



### functions to update the prior information when the tree structure is changed

# update prior choices for each split based on user input
update_prior_w <- function(prior_data, node_data, input){

  if (length(input$current_node_id) != 1) return(prior_data) # in case this function is called for more than one node (should not be possible)

  # if we have no splits left, we return an empty prior list
  if (length(get_split_ids(node_data)) == 0) return(list())

  node_id <- input$current_node_id

  split_index <- get_prior_number(prior_data, input$current_node_id) # index in prior_data for the id of the split node

  if (input$weight_prior == "dirichlet"){ # dirichlet prior
    prior_data[[split_index]] <- list(id = node_id,
                                      name = get_node_name(node_data,node_id),
                                      prior = "dirichlet",
                                      param = get_diri_param(sum(node_data$edges$from == node_id)),
                                      children = node_data$edges$to[node_data$edges$from == node_id],
                                      no_children = sum(node_data$edges$from == node_id))
  } else { # pc prior
    # choose one of the nodes to be the basemodel (just take the one of them, user can set basemodel to 0 if the other should be basemodel)
    basemod_id <- node_data$edges$to[node_data$edges$from == node_id][1]
    # in the basemodel, we have "basemodel" amount of the variance going to the node "above_node"
    # the median is the median of w_{above_node/both_nodes}
    # (thus basemodel = 0 to above_node means no variance to that node, and basemodel == 1 means all variance to that node)

    #basemodel_value <- input$basemodel

    if (input$basemodel == "mod1"){
      basemodel_value <- 1
    } else if (input$basemodel == "mod2"){
      basemodel_value <- 0
    } else {
      basemodel_value <- input$median
    }

    prior_data[[split_index]] <- list(id = node_id,
                                      name = get_node_name(node_data, node_id),
                                      prior = "pc",
                                      param = data.frame(above_node = basemod_id, basemodel = basemodel_value,
                                                         median = input$median, concentration = input$conc_param),
                                      children = node_data$edges$to[node_data$edges$from == node_id],
                                      no_children = 2) # can only be 2 with pc prior
  }

  return(prior_data)

}


# update prior when the tree changes, default dirichlet on all splits
# for each merge/removal the nodes at that level and one level up are affected
update_prior_w_merge_remove <- function(prior_data, node_data, change, input){

  if (change == "merge"){ # if a new node was made

    # add dirichlet on the new split and on the "old" split

    # get the parent of the new nodes
    new_parent_node <- sapply(input$current_node_id, function(x) get_parent_node_id(node_data, x))
    stopifnot(length(unique(new_parent_node)) == 1) # if the selected nodes do not have the same parent, something is wrong
    new_parent_node <- new_parent_node[1]
    # prior of the new split
    new_prior <- list(id = new_parent_node,
                      #name = as.character(node_data$nodes$label[new_parent_node]),
                      name = as.character(node_data$nodes$label[node_data$nodes$id == new_parent_node]),
                      prior = "dirichlet",
                      param = get_diri_param(sum(node_data$edges$from == new_parent_node)),
                      children = node_data$edges$to[node_data$edges$from == new_parent_node],
                      no_children = sum(node_data$edges$from == new_parent_node))

    prior_data[[length(prior_data)+1]] <- new_prior # adding prior on new split

    # change dirichlet on root note
    # update prior on the old parent of the selected (i.e. merged) nodes, which is now the parent of the new_parent_node
    old_parent_node <- get_parent_node_id(node_data, new_parent_node)
    # which of the elements in prior_data correspronds to this node
    new_old_prior <- list(id = old_parent_node,
                          name = get_node_name(node_data, old_parent_node),
                          prior = "dirichlet",
                          param = get_diri_param(sum(node_data$edges$from == old_parent_node)),
                          children = node_data$edges$to[node_data$edges$from == old_parent_node],
                          no_children = sum(node_data$edges$from == old_parent_node))

    prior_data[[get_prior_number(prior_data, old_parent_node)]] <- new_old_prior # adding new prior for old_parent_node

  } else if (change == "remove") { # if a node was removed

    # get node id of the node that is removed, this may be another node id now!!
    node_id <- input$current_node_id
    # get split id of the

    # find children of the removed node
    children <- prior_data[[get_prior_number(prior_data, node_id)]]$children
    # find the new parent of these nodes (enough to test one of them, since they have the same parent)
    new_parent_node <- node_data$edges$from[node_data$edges$to == children[1]]
    new_parent_split_index <- get_prior_number(prior_data, new_parent_node) # the index in prior_data for this node
    # get all children of this parent node (may be more)
    all_children <- node_data$edges$to[node_data$edges$from == new_parent_node]


    # update the prior of the new parent node
    prior_data[[new_parent_split_index]]$prior <- "dirichlet"
    prior_data[[new_parent_split_index]]$children <- node_data$edges$to[node_data$edges$from == new_parent_node]
    prior_data[[new_parent_split_index]]$no_children <- length(prior_data[[new_parent_split_index]]$children)
    prior_data[[new_parent_split_index]]$param <- get_diri_param(prior_data[[new_parent_split_index]]$no_children)

    # remove object from prior_data that does no longer exist
    prior_data[[get_prior_number(prior_data, node_id)]] <- NULL
    # if (nrow(node_data$nodes) < max(node_data$nodes$id)){ # if the node removed did not have the highest node id
    #   prior_data[[node_id]] <- NULL
    # } else { # if the node removed was the one with the highest id
    #   prior_data[[node_id]] <- list(id = node_id)
    # }

  }

  return(prior_data)

}

update_prior_w_attach_detach <- function(prior_data, node_data, change, input){

  if (change == "detach"){ # if a node was detached

    node_id <- input$current_node_id
    # put a dirichlet prior on the split that was affected by the removal
    # find the split-node that had the detached node as a child, must look in prior_data since it is removed from node_data
    old_parent_node <- unlist(sapply(prior_data, function(x) if (node_id %in% x$children) return(x$id)))

    # if this node still has at least two children (if not, it is already removed from the node_data)
    if (no_of_children(node_data, list(current_node_id = old_parent_node)) >= 2){
      prior_data[[get_prior_number(prior_data, old_parent_node)]]$name <- as.character(node_data$nodes$label[old_parent_node])
      prior_data[[get_prior_number(prior_data, old_parent_node)]]$param <- get_diri_param(sum(node_data$edges$from == old_parent_node))
      prior_data[[get_prior_number(prior_data, old_parent_node)]]$children <- node_data$edges$to[node_data$edges$from == old_parent_node]
      prior_data[[get_prior_number(prior_data, old_parent_node)]]$no_children <- sum(node_data$edges$from == old_parent_node)
    } else {

      # if the parent node does no longer exist, since it was merged with another node, we remove the prior of this node
      # and update the child nodes of the node that is now the parent of the sibling of the detached node (which has the same number of children)

      # store the sibling of the detached node
      sibling_nodes <- prior_data[[get_prior_number(prior_data, old_parent_node)]]$children
      sibling_node <- sibling_nodes[!(sibling_nodes %in% node_id)]

      # if it is only one sibling node and it is a split node, we just remove the prior
      if (length(sibling_node) == 1 && is_split_node(node_data, list(current_node_id = sibling_node))){
        prior_data <- prior_data[-get_prior_number(prior_data, old_parent_node)]
      } else {

        # remove prior of split that does no longer exist
        prior_data <- prior_data[-get_prior_number(prior_data, old_parent_node)]

        # only do this if the tree still has attached nodes, i.e. if there is still a tree here, and if the detached node had children
        # the old_parent_node does not exist anymore, so we must find a node that is in the tree where the detached node used to belong
        if (length(nodes_in_tree(node_data, sibling_node[1])) >= 2){
          #if (sum(node_data$nodes$status == "attached") > 0){
          # must find the parent of the parent node that was removed
          new_parent_node <- get_split_ids(node_data)[sapply(get_split_ids(node_data),
                                                             function(x) old_parent_node %in% prior_data[[get_prior_number(prior_data, x)]]$children)]

          # update the children and the name of this split node
          prior_data[[get_prior_number(prior_data, new_parent_node)]]$children <- node_data$edges$to[node_data$edges$from == new_parent_node]
          prior_data[[get_prior_number(prior_data, new_parent_node)]]$name <- paste(sapply(
            prior_data[[get_prior_number(prior_data, new_parent_node)]]$children,
            function(x) get_node_name(node_data, x)),
            collapse = "_")
        }

      }

    }

  } else { # if nodes are being attached to the tree

    # NOTE: when no attached nodes, we are merging them, not attaching them! So we are in this function only
    # when we have an existing tree

    node_ids <- input$current_node_id

    parent_node <- get_parent_node_id(node_data, node_ids[1]) # they must all have the same parent, not allowed to attach to more than one node

    # if the parent was not a parent before
    if (length(get_prior_number(prior_data, parent_node)) == 0){

      new_prior <- list(id = parent_node,
                        name = as.character(node_data$nodes$label[parent_node]),
                        prior = "dirichlet",
                        param = get_diri_param(sum(node_data$edges$from == parent_node)),
                        children = node_data$edges$to[node_data$edges$from == parent_node],
                        no_children = sum(node_data$edges$from == parent_node))

      # must update the prior for the other as well
      for (ind in 1:length(prior_data)){

        this_parent <- prior_data[[ind]]$id

        prior_data[[ind]]$name <- as.character(node_data$nodes$label[node_data$nodes$id == this_parent])
        prior_data[[ind]]$prior <- "dirichlet"
        prior_data[[ind]]$param <- get_diri_param(sum(node_data$edges$from == this_parent))
        prior_data[[ind]]$children <- node_data$edges$to[node_data$edges$from == this_parent]
        prior_data[[ind]]$no_children <- sum(node_data$edges$from == this_parent)

      }

      prior_data[[length(prior_data)+1]] <- new_prior # adding prior on new split

    } else {
      # if the node was already a parent, we add the new nodes and change the prior to dirichlet
      prior_id <- get_prior_number(prior_data, parent_node)
      prior_data[[prior_id]]$name <- as.character(node_data$nodes$label[parent_node])
      prior_data[[prior_id]]$prior <- "dirichlet"
      prior_data[[prior_id]]$children <- node_data$edges$to[node_data$edges$from == parent_node]
      prior_data[[prior_id]]$no_children <- length(prior_data[[prior_id]]$children)
      prior_data[[prior_id]]$param <- get_diri_param(prior_data[[prior_id]]$no_children)
    }

  }

  return(prior_data)

}

# get the prior-info for node (returns NULL if the node is a leaf node, which has no priors (all other nodes have priors))
get_prior_info <- function(prior_data, node_data, input) {
  if (is_leaf_node(node_data, input)) return(NULL)
  return(prior_data[[get_prior_number(prior_data, input$current_node_id)]])
}

# update the CW priors in the tree
update_variance_prior_attach_detach <- function(prior_data, node_data, change, input){

  # if a node was detached
  if (change == "detach"){

    #node_id <- input$current_node_id

    # for all nodes that has status = "detached", but no variance prior
    all_detached_nodes <- node_data$nodes$id[node_data$nodes$status == "detached"]
    has_no_prior <- unlist(sapply(prior_data, function(x) if (x$prior == "") x$id))

    need_prior_id <- has_no_prior[has_no_prior %in% all_detached_nodes]
    need_prior_number <- which(sapply(prior_data, function(x) x$id %in% need_prior_id))

    #node_number <- which(sapply(prior_data, function(x) x$id) == node_id)
    for (ind in 1:length(need_prior_number)){
      #browser()
      # by default, give a PC(3, 0.05) prior to the component that is detached
      prior_data[[need_prior_number[ind]]] <- list(id = need_prior_id[ind],
                                        name = get_node_name(node_data, need_prior_id[ind]),
                                        prior = "pc0",
                                        param = c(3, 0.05))
    }


  } else { # attach
    # remove the variance prior on the nodes that are now attached to the tree
    # the nodes that has a variance prior, but are attached in a tree
    has_prior <- sapply(prior_data, function(x) if (x$prior != "") x$id else 0)
    for (ind in c(1:length(has_prior))[has_prior %in% node_data$nodes$id[node_data$nodes$status == "attached"]]){
      prior_data[[ind]]$prior <- ""
      prior_data[[ind]]$param <- 0
    }
  }

  return(prior_data)

}

# update total variance for new trees
update_V_prior_merged_detached_nodes <- function(prior_data, node_data, input){

  # check if the first element of the prior_data is null, if so we replace it instead of appending
  if (length(prior_data) == 0 || is.null(prior_data[[1]])) prior_data <- list()
  # add new total variance prior
  prior_data <- c(prior_data, list(list(
    id = get_parent_node_id(node_data, input$current_node_id[1]), # all merged nodes have the same parent
    name = get_node_name(node_data, get_parent_node_id(node_data, input$current_node_id[1])), # name of parent node
    prior = "pc0", # cannot have jeffreys prior if we do not have a single tree with no CW nodes
    param = c(3, 0.05)
  )))

  return(prior_data)

}

# update total variance for merged trees
update_V_prior_merged_trees <- function(prior_data, node_data, input){

  # add new total variance prior (is set to be pc-prior even though we may have only one tree now)
  prior_data <- c(prior_data, list(list(
    id = get_parent_node_id(node_data, input$current_node_id[1]), # all merged nodes have the same parent
    name = get_node_name(node_data, get_parent_node_id(node_data, input$current_node_id[1])), # name of parent node
    prior = "pc0",
    param = c(3, 0.05)
  )))

  # remove old prior data for the nodes that are no longer top nodes
  old_prior_numbers <- sapply(input$current_node_id, function(x) get_prior_number(prior_data, x))
  prior_data <- prior_data[-old_prior_numbers]

  return(prior_data)

}

# updates the names and priors of the total variance priors when we attach or detach nodes
# sibling_removed = (one of the) sibling(s) of the detached node
update_V_prior_attach_detach <- function(prior_data, node_data, change, sibling_removed = 0){

  # if we do not have one single tree (i.e., user detached a node), we do not allow jeffreys prior
  # only done first time the user detached, so only when we have only one detached node, and the prior is jeffreys
  if (length(prior_data) == 1 && sum(node_data$nodes$status == "detached") == 1 && change == "detach"){
    prior_data[[1]]$prior <- "pc0"
    prior_data[[1]]$param <- c(3, 0.05)
  }

  # if there are no attached nodes, remove all total variance priors
  if (sum(node_data$nodes$status == "attached") == 0){
    return(list())
  }

  for (ind in 1:length(prior_data)){
    # if the top node is still a top node, update the name of the total variance
    if (is_top_node(node_data, list(current_node_id = prior_data[[ind]]$id))){
      prior_data[[ind]]$name <- get_node_name(node_data, prior_data[[ind]]$id)
    } else if (sum(nodes_in_tree(node_data, prior_data[[ind]]$id)) == 0){ # if this tree disappeared, remove the total variance
      #prior_data[[ind]] <- list(NULL)
      prior_data[[ind]] <- NULL
    } else { # if we have a top node, but it is the child of the old top node, make new total variance prior for this new top node
      # if the top node is removed, the removed node only had one sibling
      prior_data[[ind]]$id <- sibling_removed
      prior_data[[ind]]$name <- get_node_name(node_data, prior_data[[ind]]$id)
      prior_data[[ind]]$prior <- "pc0"
      prior_data[[ind]]$param <- c(3, 0.05)
    }
  }

  return(prior_data)

}

# update prior on variance (CW or total) when user changes parameters/prior family
update_var_prior <- function(prior_data, input){
  prior_data[[get_prior_number(prior_data, input$current_node_id)]]$prior <- input$var_prior
  prior_data[[get_prior_number(prior_data, input$current_node_id)]]$param <- update_var_prior_params(input$var_prior, c(input$par1, input$par2))
  return(prior_data)
}

update_var_prior_params <- function(prior_name, params){
  if (prior_name %in% c("pc0", "invgam")){
    return(params)
  } else if (prior_name %in% c("hc", "hn")){
    return(c(params[1], 0))
  } else if (prior_name == "jeffreys"){
    return(c(0, 0))
  } else {
    stop("Not valid prior for variance parameter.", call. = FALSE)
  }
}


### other things

# updating the edges of a split
update_basemodel_edges <- function(node_data, prior_data){

  if (nrow(node_data$edges) == 0) return(node_data)

  # if there is no info here (i.e. a merge of detached nodes), add info on dashes and width
  if (ncol(node_data$edges) == 2){
    # set edges to dirichlet-edges (equal width and no dashes)
    new_edges <- cbind(node_data$edges, data.frame(width = rep(NA, nrow(node_data$edges)), dashes = rep(FALSE, nrow(node_data$edges))))
  } else { # remove
    new_edges <- node_data$edges
  }

  new_edges$label <- ""

  # for each of the splits in the tree
  for (split_number in 1:length(prior_data)){

    # prior_data includes id of each split-node and what prior that split has
    split_id <- prior_data[[split_number]]$id
    sub_edges <- new_edges[new_edges$from == split_id,] # the edges from current split

    # if the split has a pc prior, we must change width and basemodel info
    # only dual splits can have a pc prior
    if (prior_data[[split_number]]$prior == "pc"){

      base_place <- sub_edges$to == prior_data[[split_number]]$param$above_node # index of basemodel_node
      not_base_place <- sub_edges$to != prior_data[[split_number]]$param$above_node # the other node in this split

      # the basemodel_node gets the median the user chooses, the other node gets 1-median
      sub_edges$width[base_place] <- prior_data[[split_number]]$param$median
      sub_edges$width[not_base_place] <- 1-prior_data[[split_number]]$param$median

      # indicate where the basemodel is using dashed arrows
      # edge to basemodel_node has dashes if basemodel = 0
      # edge to the other node has dashes if basemodel = 1
      sub_edges$dashes[base_place] <- if (prior_data[[split_number]]$param$basemodel == 0) TRUE else FALSE
      sub_edges$dashes[not_base_place] <- if (prior_data[[split_number]]$param$basemodel == 1) TRUE else FALSE

      sub_edges$label[base_place] <- "A"
      sub_edges$label[not_base_place] <- "B"

    } else {

      sub_edges$width <- 1/prior_data[[split_number]]$no_children
      sub_edges$dashes <- FALSE # no dashes if not pc prior

    }

    new_edges[new_edges$from == split_id,] <- sub_edges

  }

  node_data$edges <- new_edges

  return(node_data)

}

# update the node labels so each split node has name according to the nodes below it
update_node_labels <- function(node_data){

  split_nodes <- get_split_ids(node_data)
  if (is.null(split_nodes)) return(node_data) # if no splits, we do not change anything

  node_data$nodes$label <- as.character(node_data$nodes$label)

  for (node_id in split_nodes){
    new_label <- paste(sapply(get_leaf_nodes(node_data, node_id), function(x) get_node_name(node_data, x)), collapse = "_")
    # new_label <- paste(sapply(get_leaf_nodes(node_data, node_id), function(x) get_node_name(node_data, x)), collapse = "_")
    node_data$nodes$label[node_data$nodes$id == node_id] <- new_label
  }

  return(node_data)

}

# adds or removes the A and B for PC prior basemodel selection
add_remove_AB <- function(node_data, add, node_ids = 0){

  if (nrow(node_data$edges) == 0) return(node_data)
  if (add){
    # remove old names in case this is another node with a PC prior
    node_data <- add_remove_AB(node_data, FALSE)
    node_data$nodes$label[node_data$nodes$id == node_ids[1]] <- paste0("A: ", node_data$nodes$label[node_data$nodes$id == node_ids[1]])
    node_data$nodes$label[node_data$nodes$id == node_ids[2]] <- paste0("B: ", node_data$nodes$label[node_data$nodes$id == node_ids[2]])
  } else {
    node_data$nodes$label[grepl(": ", node_data$nodes$label)] <-
      substr(node_data$nodes$label[grepl(": ", node_data$nodes$label)], 4, nchar(node_data$nodes$label[grepl(": ", node_data$nodes$label)]))
  }

  return(node_data)

}

# calculate the variance proportion on each node according to how much variance goes to each of them
# detached nodes are not included
calc_variance_proportions <- function(node_data){

  tmp_nodes <- node_data$nodes
  tmp_nodes$varprop <- NA
  tmp_edges <- node_data$edges

  node_id_by_level <- tmp_nodes$id[order(tmp_nodes$level)]

  # calculate the amount of variance in each node in the tree, by level
  for (node_id in node_id_by_level){

    if (tmp_nodes$status[tmp_nodes$id == node_id] == "detached") {
      tmp_nodes$varprop[tmp_nodes$id == node_id] <- 0
      next
    }

    # if it is the top node
    #if (node_id==10) browser()
    if (is_top_node(node_data, list(current_node_id = node_id))) {
      tmp_nodes$varprop[tmp_nodes$id == node_id] <- 1
    } else {

      parent_id <- get_parent_node_id(node_data, node_id)

      # must find the width of the edge to this node
      val <- tmp_edges$width[tmp_edges$to == node_id]
      # then we must multiply this with variance proportion that goes to the parent node
      val <- val * tmp_nodes$varprop[tmp_nodes$id == parent_id]

      tmp_nodes$varprop[tmp_nodes$id == node_id] <- val

    }

  }

  # in case something goes wrong here
  if (any(tmp_nodes$varprop < 0) || any(tmp_nodes$varprop > 1)) tmp_nodes$varprop <- rep(1/length(tmp_nodes$varprop), length(tmp_nodes$varprop))

  node_data$nodes$varprop <- tmp_nodes$varprop

  return(node_data)

}

# all nodes get shortnames in the tree, and the full names are available on original nodes (leaf nodes) when hovering
# returns the node-data-frame with new labels and titles
set_node_title <- function(nodes, node_data){

  titles <- c()
  for (i in 1:nrow(nodes)){
    titles[i] <- if (is_leaf_node(node_data, list(current_node_id = nodes$id[i]))) as.character(nodes$label[i]) else NA
  }

  nodes$title <- titles
  nodes$label <- shorten_names(nodes)

  return(nodes)

}

shorten_names <- function(nodes){

  new_nodenames <- c()
  nodenames <- nodes$label
  # max three letters for each node, additional number if necessary
  for (i in 1:length(nodenames)){
    new_nodename <- substr(nodenames[i], 1, 3)
    if (new_nodename %in% new_nodenames){
      which_equal_names <- which(sapply(new_nodenames, function(x) substr(x, 1, 3)) %in% new_nodename)
      new_nodename <- paste0(new_nodename, length(which_equal_names))
    }
    new_nodenames[i] <- new_nodename
  }

  return(new_nodenames)

}

# returns a hex color with intensity given by num (num = 1 is 100%, num = 0 gives another color)
node_palette <- function(num) {
  cols <- colorRampPalette(c("#FFFFFF", node_color))(101)[round(num*100, 0)+1]
  return(cols)
}

node_color <- "#E85A4F" # col1
detached_node_color <- "white" # col5
node_border_color <- gray(0.35)
node_highlight_border_color <- "black"
guide_highlight_parent_color <- node_color # col4
guide_highlight_children_color <- node_color #colorRampPalette(c("#FFFFFF", guide_highlight_parent_color))(3)[2]
plot_color <- "#8E8D8A"
plot_text_color <- "white"



### functions for enabling/disabling merge/remove buttons and numeric inputs

# attach/detach is just one test, do not need separate function
# make for the rest if we need an object keeping track of enabled/disabled buttons etc

enable_merge <- function(node_data, input){

  node_ids <- input$current_node_id

  # must choose at least two nodes
  if (length(node_ids) < 2) return(FALSE)

  # can only merge two or more nodes at the same level/that has the same parent, and when that parent has at least three children,
  # where at least one is not chosen
  same_level <- length(unique(sapply(node_ids, function(id) get_parent_node_id(node_data, id)))) == 1
  same_parent <- length(node_data$edges$to[node_data$edges$from == get_parent_node_id(node_data, node_ids[1])]) > 2
  more_siblings <- sum(node_data$nodes$level == node_data$nodes$level[node_data$nodes$id == node_ids[1]]) > length(node_ids)

  if (same_level && same_parent && more_siblings) return(TRUE)

  # if all chosen nodes are detached
  if (are_detached(node_data, input)) return(TRUE)

  # merge trees
  all_top_nodes <- sum(node_data$nodes$top_node[node_data$nodes$id %in% node_ids]) == length(node_ids)
  if (all_top_nodes) return(TRUE)

  # all other cases
  return(FALSE)

}

enable_remove <- function(node_data, input){

  node_id <- input$current_node_id

  # must be one node
  if (length(node_id) != 1) return(FALSE)

  # must be split node that is not a top node
  if (is_split_node(node_data, input) && !is_top_node(node_data, input)) return(TRUE)

  # all other cases
  return(FALSE)

}


### making the variance equations

text_latex_param <- function(node_data, prior_weight, prior_totvar, prior_cw, response_name, func_of_var){

  if (!func_of_var) {
    make_var_eq_latex(node_data, prior_weight, prior_totvar, prior_cw, response_name)
  } else {
    make_par_eq_latex(node_data, prior_weight, prior_totvar, prior_cw, response_name)
  }

}

# variances as function of parameters
make_var_eq_latex <- function(node_data, prior_weight, prior_totvar, prior_cw, response_name){

  # for each split in each tree, we store the name of the weight(s) of that split

  comps <- c()

  for (tree in seq_len(length(prior_totvar))){

    # in this case, we have no total variance and thus no HD prior
    if (is.null(prior_totvar[[tree]])) break

    top_node <- prior_totvar[[tree]]$id
    weight_frame <- data.frame()

    totvar <- sprintf("V_{%s}", get_node_strings(get_leaf_nodes(node_data, top_node), node_data))

    node_data_this_tree <- get_node_data_for_tree(node_data, top_node)

    split_nodes <- get_split_ids(node_data_this_tree)

    for (split_id in split_nodes){

      children <- get_children(node_data_this_tree, list(current_node_id = split_id))

      under <- get_node_strings(get_leaf_nodes(node_data_this_tree, split_id), node_data_this_tree)
      over <- sapply(children, function(x) get_node_strings(get_leaf_nodes(node_data_this_tree, x), node_data_this_tree))

      weight_names <- sapply(1:(length(children)-1), function(x) sprintf("\\omega_{\\frac{%s}{%s}}", over[x], under))
      weight_names <- c(weight_names, paste(c("(1-", paste(weight_names, sep = "", collapse = "-"), ")"), sep = "", collapse = ""))

      weight_frame <- rbind(weight_frame, data.frame(id = children, wn = weight_names))

    }

    # for each leaf node in this tree, we must multiply the correct weights and total variances

    leaf_nodes <- get_leaf_nodes(node_data_this_tree, top_node)

    full_weights <- c()

    for (id in leaf_nodes){

        full_weight <- paste(weight_frame[weight_frame$id %in% c(id, get_nodes_above(node_data_this_tree, id)),]$wn, sep = "", collapse = "")

        # multiply by total variance
        full_weight <- paste(totvar, full_weight, sep = "", collapse = "")

        # this is equal to the variance of some component
        full_weight <- paste(sprintf("\\sigma_{%s}^2 =", get_node_name(node_data_this_tree, id)), full_weight)

        comps <- c(comps, full_weight)

    }

  }

  # if any CW priors, just write sigma = sigma

  for (ind in 1:length(prior_cw)){
    if (prior_cw[[ind]]$prior != ""){
      comps <- c(comps, sprintf("\\sigma_{%s}^2 = \\sigma_{%s}^2", prior_cw[[ind]]$name, prior_cw[[ind]]$name))
    }
  }

  return(
    withMathJax(
      helpText(
        HTML(paste(
          paste("\\(", comps, "\\)"),
          collapse = " <br/> "
        ))
      )
    )
  )


}

# parameters as function of variances
make_par_eq_latex <- function(node_data, prior_weight, prior_totvar, prior_cw, response_name){

  comps <- c()

  # for each total variance, find the variances it is a sum of.
  # then find the weights belonging to this total variance and which variances that are involved in that weight.
  # at last include CW priors (if any).

  for (tree in seq_len(length(prior_totvar))){

    top_node <- prior_totvar[[tree]]$id
    leaf_nodes <- get_leaf_nodes(node_data, top_node)

    totvar <- paste0(
      "V_{", paste0(gsub("_", "\\_", get_node_name(node_data, leaf_nodes), fixed = TRUE), sep = "", collapse = ","),
      "}",
      "=",
      paste0("\\sigma_{", gsub("_", "\\_", get_node_name(node_data, leaf_nodes), fixed = TRUE), "}^2", sep = "", collapse = "+")
    )

    weights <- c()

    for (we in seq_len(length(prior_weight))){

      split_node <- prior_weight[[we]]$id
      leaf_nodes_in_split <- get_leaf_nodes(node_data, split_node)
      child_nodes <- get_children(node_data, list(current_node_id = split_node))

      under <- gsub("_", "\\_", get_node_name(node_data, leaf_nodes_in_split), fixed = TRUE)

      if (prior_weight[[we]]$no_children == 2){ # dual split

        weights <- c(weights, paste0(
          "\\omega_{\\frac{",
          paste0(gsub("_", "\\_", get_node_name(node_data, child_nodes[1]), fixed = TRUE), sep = "", collapse = ","),
          "}{",
          paste0(under, sep = "", collapse = ","),
          "}} = \\frac{",
          paste0("\\sigma_{", gsub("_", "\\_", get_node_name(node_data, child_nodes[1]), fixed = TRUE), "}^2", sep = "", collapse = "+"),
          "}{",
          paste0("\\sigma_{", under, "}^2", sep = "", collapse = "+"),
          "}"
        )
        )

      } else { # multi-split

        for (ind in 1:(prior_weight[[we]]$no_children-1)){
          weights <- c(weights, paste0(
            "\\omega_{\\frac{",
            paste0(gsub("_", "\\_", get_node_name(node_data, child_nodes[ind]), fixed = TRUE), sep = "", collapse = ","),
            "}{",
            paste0(under, sep = "", collapse = ","),
            "}} = \\frac{",
            paste0("\\sigma_{", gsub("_", "\\_", get_node_name(node_data, child_nodes[ind]), fixed = TRUE), "}^2", sep = "", collapse = "+"),
            "}{",
            paste0("\\sigma_{", under, "}^2", sep = "", collapse = "+"),
            "}"
          ))
        }

      }

    }

    comps <- c(comps, c(totvar, weights))

  }

  cw <- c()
  for (cw_ind in seq_len(length(prior_cw))){
    if (prior_cw[[cw_ind]]$prior != ""){
      cw <- c(cw, paste0(
        "\\sigma_{", gsub("_", "\\_", prior_cw[[cw_ind]]$name, fixed = TRUE), "}^2", "=",
        "\\sigma_{", gsub("_", "\\_", prior_cw[[cw_ind]]$name, fixed = TRUE), "}^2"
      ))
    }
  }


  return(
    withMathJax(
      helpText(
        HTML(paste(
          paste("\\(", c(comps, cw), "\\)"),
          collapse = " <br/> "
        ))
      )
    )
  )

}

# return latex-code for the expression for the prior distribution on the parameter provided
get_prior_expr <- function(prior_data, node_data, param = c("cw", "totvar", "weight")){

  param <- match.arg(param)

  if (param == "cw"){

    if (prior_data$prior == "pc0") {
      param_name <- sprintf("\\sigma_{%s}", gsub("_", "\\_", prior_data$name, fixed = TRUE)) #, "}")
      prior_expr <- sprintf("\\mathrm{PC}_{\\mathrm{0}}(%s, %s)", prior_data$param[1], prior_data$param[2])
    } else if (prior_data$prior == "jeffreys"){
      param_name <- sprintf("\\sigma_{%s}^2", gsub("_", "\\_", prior_data$name, fixed = TRUE))
      prior_expr <- paste0("\\mathrm{Jeffreys'}")
    } else if (prior_data$prior == "invgam"){
      param_name <- sprintf("\\sigma_{%s}^2", gsub("_", "\\_", prior_data$name, fixed = TRUE))
      prior_expr <- paste0("\\mathrm{InvGam}(", prior_data$param[1], ",", prior_data$param[2], ")")
    } else if (prior_data$prior == "hc"){
      param_name <- sprintf("\\sigma_{%s}", gsub("_", "\\_", prior_data$name, fixed = TRUE)) #, "}")
      prior_expr <- paste0("\\mathrm{HC}(", prior_data$param[1], ")")
    }

  } else if (param == "totvar") {

    if (prior_data$prior == "pc0") {
      param_name <- sprintf("\\sqrt{V_{%s}}", get_node_strings(get_leaf_nodes(node_data, prior_data$id), node_data))
      prior_expr <- sprintf("\\mathrm{PC}_{\\mathrm{0}}(%s, %s)", prior_data$param[1], prior_data$param[2])
    } else if (prior_data$prior == "jeffreys"){
      param_name <- sprintf("V_{%s}", get_node_strings(get_leaf_nodes(node_data, prior_data$id), node_data))
      prior_expr <- paste0("\\mathrm{Jeffreys'}")
    } else if (prior_data$prior == "invgam"){
      param_name <- sprintf("V_{%s}", get_node_strings(get_leaf_nodes(node_data, prior_data$id), node_data))
      prior_expr <- sprintf("\\mathrm{InvGam}(%s, %s)", prior_data$param[1], prior_data$param[2])
    } else if (prior_data$prior == "hc"){
      param_name <- sprintf("\\sqrt{V_{%s}}", get_node_strings(get_leaf_nodes(node_data, prior_data$id), node_data))
      prior_expr <- sprintf("\\mathrm{HC}(%s)", prior_data$param[1])
    }

  } else { # weight

    if (prior_data$prior == "dirichlet"){

      param_name <- c()
      for (i in 1:(prior_data$no_children-1)){
        t_under <- get_node_strings(get_leaf_nodes(node_data, prior_data$id), node_data)
        t_over <- get_node_strings(get_leaf_nodes(node_data, prior_data$children[i]), node_data)
        param_name[i] <- sprintf("\\omega_{\\frac{%s}{%s}}", t_over, t_under)
      }
      param_name <- paste0(param_name, sep = "", collapse = ", ")
      if (prior_data$no_children > 2) param_name <- paste0("(", param_name, ")", sep = "", collapse = "")

      prior_expr <- sprintf("\\mathrm{Dirichlet}(%s)", prior_data$no_children)
    } else if (prior_data$prior == "pc"){

      t_under <- get_node_strings(get_leaf_nodes(node_data, prior_data$id), node_data)
      t_over <- get_node_strings(get_leaf_nodes(node_data, prior_data$param$above_node), node_data)
      param_name <- sprintf("\\omega_{\\frac{%s}{%s}}", t_over, t_under)

      if (prior_data$param$basemodel == 0){
        prior_expr <- sprintf("\\mathrm{PC}_{\\mathrm{0}}(%s)", prior_data$param$median)
      } else if (prior_data$param$basemodel == 1) {
        prior_expr <- sprintf("\\mathrm{PC}_{\\mathrm{1}}(%s)", prior_data$param$median)
      } else {
        prior_expr <- sprintf("\\mathrm{PC}_{\\mathrm{M}}(%s, %s)", prior_data$param$median, prior_data$param$concentration)
      }

    }

  }

  return(
    paste0(param_name, " \\sim ", prior_expr)
  )

}

# make proper strings for the node-names, with "," instead of "_" (not for the node_data-objects, only latex-things)
get_node_strings <- function(nodes, node_data){

  this_name <- paste0(get_node_name(node_data, nodes), sep = "", collapse = ",")
  this_name <- gsub("_", "\\_", this_name, fixed = TRUE)

  return(this_name)

}

# making expressions for the prior of each model parameter
make_pr_dists_latex <- function(node_data, prior_weight, prior_totvar, prior_cw, input_id){

  pr_exprs <- c()
  bold <- if (length(input_id) == 1) TRUE else FALSE

  # for CW priors:
  for (i in 1:length(prior_cw)){
    if (prior_cw[[i]]$prior != "") {
      tmp_pr <- get_prior_expr(prior_cw[[i]], node_data, "cw")
      if (bold && prior_cw[[i]]$id == input_id){
        tmp_pr <- paste0(
          "{\\color{black}{", tmp_pr, "}}"
        )
      }
      pr_exprs <- c(pr_exprs, tmp_pr)
    }
  }

  # for total variances
  for (i in seq_len(length(prior_totvar))){
    tmp_pr <- get_prior_expr(prior_totvar[[i]], node_data, "totvar")
    if (bold && prior_totvar[[i]]$id == input_id){
      tmp_pr <- paste0(
        "{\\color{black}{", tmp_pr, "}}"
      )
    }
    pr_exprs <- c(pr_exprs, tmp_pr)
  }

  # for weight priors
  for (i in seq_len(length(prior_weight))){
    tmp_pr <- get_prior_expr(prior_weight[[i]], node_data, "weight")
    if (bold && prior_weight[[i]]$id == input_id){
      tmp_pr <- paste0(
        "{\\color{black}{", tmp_pr, "}}"
      )
    }
    pr_exprs <- c(pr_exprs, tmp_pr)
  }

  return(
    withMathJax(
      helpText(
        HTML(paste(
          paste("\\(", pr_exprs, "\\)"),
          collapse = " <br/> "
        ))
      )
    )
  )


}

# make model equation in latex
make_model_eq_latex <- function(initial_args){

  if (initial_args$.family == "gaussian"){
    meq <- "Gaussian likelihood. \\(y = "
  } else if (initial_args$.family == "binomial") {
    meq <- "Binomial likelihood. \\(\\eta = "
  } else {
    meq <- "Poisson likelihood. \\(\\eta = "
  }

  if (initial_args$.use_intercept) meq <- paste0(meq, "\\mu +", collapse = "")
  for (ind in seq_along(initial_args$.fixef_names))
    meq <- paste0(meq, paste0(c(initial_args$.fixef_names[ind], "\\beta_{", gsub("_", "\\_", initial_args$.fixef_names[ind], fixed = TRUE),
                                "} ",
                                " + "), collapse = ""), collapse = "")
  for (ind in seq_along(initial_args$.nodenames)){
    # if (initial_args$.nodenames[ind] == "eps"){ # if this is the residual term, use epsilon instead of "eps"
    #   meq <- paste0(c("(+ ", meq, gsub("_", "\\_", initial_args$.nodenames[ind], fixed = TRUE, ")")), collapse = "")
    # } else {
      meq <- paste0(c(meq, gsub("_", "\\_", initial_args$.nodenames[ind], fixed = TRUE)), collapse = "")
    # }

    if (ind < length(initial_args$.nodenames)) meq <- paste0(c(meq, " + "), collapse = "")
  }
  meq <- paste0(meq, ", \\ \\eta = g(\\mu) = ", collapse = "")
  meq <- if (initial_args$.family == "gaussian"){
    paste0(meq, "\\mu", collapse = "")
  } else if (initial_args$.family == "poisson") {
    paste0(meq, "\\mathrm{log}(\\mu)", collapse = "")
  } else if (initial_args$.family == "binomial"){
    paste0(meq, "\\mathrm{logit}(\\mu)", collapse = "")
  }

  return(paste0(c(meq, "\\)"), collapse = ""))
  # return(paste0(c("\\(", meq, "\\)"), collapse = ""))

}

# make model equation, in r mathematical expression (not used)
make_model_eq_rexpr <- function(initial_args){
  meq <- c(expression(eta ~ "="))
  if (initial_args$.use_intercept) meq <- c(meq, expression(mu ~ "+"))
  for (ind in seq_along(initial_args$.fixef_names)) meq <- c(meq, bquote(beta[.(initial_args$.fixef_names[ind])] ~ "+"))
  for (ind in seq_along(initial_args$.nodenames)){
    if (ind == length(initial_args$.nodenames)){ # no plus on the last element, as it is the last
      meq <- c(meq, bquote(.(initial_args$.nodenames[ind]) * ","))
    } else {
      meq <- c(meq, bquote(.(initial_args$.nodenames[ind]) ~ "+"))
    }
  }
  meq <- c(meq, expression(eta ~ "=" ~ g * "(" * mu * ")" ~ "=" ))
  meq <- if (initial_args$.family == "gaussian"){
    c(meq, expression(mu))
  } else if (initial_args$.family == "poisson"){
    c(meq, expression(log(mu)))
  } else if (initial_args$.family == "binomial"){
    c(meq, expression(logit(mu)))
  }

  if (initial_args$.no_pc){
    c(meq, "\nPC priors on splits are not computed before the app closes. This means that the PC priors on splits are not plotted in the GUI.")
  }

  return(meq)
}


### for the guide

# give a message to the user to explain what he/she should to at step 1 of the guide
guide_make_message_to_user_step1 <- function(guide_data, node_data, input){

  msg <- "All model structure settings are available in the panel 'Modify model structure' to the left.\n"

  # if the user is at step 1 (modify model structure), give message based on the number of nodes clicked

  # if something happened, reset message
  # if the user has not chosen any nodes, or the chosen node has been removed
  if (length(input$current_node_id) == 0 || !(input$current_node_id %in% node_data$nodes$id) || (input$merge || input$remove || input$attach || input$detach)){
    msg <- paste(msg, "Click a node! To choose more than one node, you can use a long-click (click and hold for a bit). ")
    # end no nodes clicked
  } else if (length(input$current_node_id) == 1){

    msg <- paste(
      msg,
      "You have chosen one node (node ", get_node_name_click(node_data, input), "). ",
      sep = ""
    )

    if (are_attached(node_data, input)){

      msg <- paste(
        msg,
        "This node is attached to a tree. "
      )

      if (is_top_node(node_data, input)){
        if (sum(node_data$nodes$top_node) == 1){
          msg <- paste(
            msg,
            "This is a top node. You cannot do anything with this, but you can change the structure of its children. Try clicking one of them! ",
            sep = ""
          )
        } else if (sum(node_data$nodes$top_node) > 1){
          msg <- paste(
            msg,
            "This is a top node. You can merge it with another top node. ",
            sep = ""
          )
        }
      } else if (is_split_node(node_data, input)){
        msg <- paste(
          msg,
          "This is a split node. You can remove it to move its children up to the level of this node. ",
          sep = ""
        )
      } else if (is_split_node(node_data, input) && length(get_siblings(node_data, input$current_node_id)) > 2) {
        msg <- paste(
          msg,
          "This is a split node. You can remove it to move its children up to the level of this node, or click one or more additional nodes to merge it with another node at this level to make a new split. ",
          sep = ""
        )
      } else if (length(get_siblings(node_data, input$current_node_id)) > 2){
        msg <- paste(
          msg,
          "You can choose to detach this node from the tree, or to click one or more additional nodes to combine them in a new split. ",
          sep = ""
        )
      } else if (length(get_siblings(node_data, input$current_node_id)) <= 2){
        msg <- paste(
          msg,
          "You can choose to detach this node from the tree. ",
          sep = ""
        )
      }

    } else if (are_detached(node_data, input)){

      msg <- paste(
        msg,
        "This node is not attached to a tree. Click on a node in a tree to attach it, or choose another detached node to make a new tree by merging the nodes. ",
        sep = ""
      )

    } else {browser(); stop("Should not be possible")}
    # end 1 node clicked
  } else if (length(input$current_node_id) > 1){

    msg <- paste(
      msg,
      "You have chosen ", length(input$current_node_id), " nodes (nodes ", paste(get_node_name_click(node_data, input), collapse = ", "), "). ",
      sep = ""
    )

    if (are_attached(node_data, input)){

      if (length(unique(node_data$nodes$level[node_data$nodes$id %in% input$current_node_id])) > 1){
        msg <- paste(
          msg,
          "The nodes you have chosen are at different levels in the tree, and then you cannot do anything to change the model structure. "
        )
      } else if (sum(sapply(input$current_node_id, function(x) is_top_node(node_data, list(current_node_id = x)))) == length(input$current_node_id)) {
        msg <- paste(
          msg,
          "You can click 'Merge' to merge these trees into one single tree. "
        )
      } else if (length(unique(sapply(input$current_node_id, function(x) get_parent_node_id(node_data, x)))) > 1){
        msg <- paste(
          msg,
          "You have chosen non-top nodes in different trees, and then you cannot do anything to change the model structure. ",
          sep = ""
        )
      } else if (length(input$current_node_id) == length(get_siblings(node_data, input$current_node_id[1]))){
        msg <- paste(
          msg,
          "You have chosen all nodes at this level, and then you cannot do anything to change the model structure. "
        )
      } else {
        msg <- paste(
          msg,
          "You can merge these nodes to make a new split node in the tree. "
        )
      }

    } else if (are_detached(node_data, input)){

      msg <- paste(
        msg,
        "If you want to merge these nodes and make a new tree, you can click 'Merge'. If you want to attach them to an existing tree, choose a node in that tree as well and click 'Attach'. ",
        sep = ""
      )

    } else if (can_be_attached(node_data, input)){

      det_nodes <- node_data$nodes$id[node_data$nodes$status == "detached"]
      det_nodes <- det_nodes[det_nodes %in% input$current_node_id]
      msg <- paste(
        msg,
        "You can attach ", paste(get_node_name(node_data, det_nodes), sep = "", collapse = ", "), " to ",
        get_node_name(node_data, input$current_node_id[!(input$current_node_id %in% det_nodes)]), ". ",
        sep = "", collapse = ""
      )

    } else {
      msg <- paste(
        msg,
        "You cannot do anything with this combination of nodes. Choose some other combination of nodes to change the model structure. "
      )
    }
    # end more than 1 node clicked
  }



  guide_data$message_to_user <- msg

  return(guide_data)

}

# give a message to the user to explain what he/she should to at step 2 of the guide
guide_make_message_to_user_step2 <- function(guide_data, node_data, prior_data, input){

  msg <- "All prior settings are available in the panel 'Change priors' to the left.\n"

  if (length(input$current_node_id) != 1){
    msg <- paste(
      msg,
      "Choose one (and only one) node to choose prior for it."
    )
  } else if (is_leaf_node(node_data, input)){
    msg <- paste(
      msg,
      "You must choose a split node, a top node, or a detached node to choose a prior."
    )
  } else if (is_top_node(node_data, input)){
    msg <- paste(
      msg,
      "This is a top node. You can choose a prior for the total variance of this tree. Note that the Half-Cauchy prior is a prior on standard deviations, while the others are on variances. "
    )
    if (sum(node_data$nodes$top_node) == 1 && sum(node_data$nodes$status == "detached") == 0){
      msg <- paste(
        msg,
        "Since you have only one tree and no individual components, you can choose a Jeffreys' prior on the total variance. "
      )
    }
  } else if (are_detached(node_data, input)){
    msg <- paste(
      msg,
      "This is a detached node. The variance of the effect this node represents can be given a prior that is independent of the other components in the model. Note that the Half-Cauchy prior is a prior on standard deviations, while the others are on variances."
    )
  }
  if (is_split_node(node_data, input)){
    if (is_dual_split_node(node_data, input)){
      msg <- paste(
        msg,
        "This is a dual-split node. If you have no preference on how the variance in the part of the model this node (",
        get_node_name_click(node_data, input),
        ") is distributed to each of the children (",
        paste(get_node_name(node_data, get_children(node_data, input)), sep = "", collapse = ", "),
        "), you can choose a Dirichlet distribution on the weight this node represents. ",
        "If you on the other hand has a preference on how the variance is distributed, or you want shrinkage to one",
        "of the effects these nodes represent, you can choose a PC prior and choose suitable parameters.",
        sep = "", collapse = ""
      )
    } else {
      msg <- paste(
        msg,
        paste(
          "This is a multi-split node. The variance in the part of the model this node (", get_node_name_click(node_data, input),
          ") represents is distributed equally to each of the children (",
          paste(get_node_name(node_data, get_children(node_data, input)), sep = "", collapse = ", "), ") through a symmetric Dirichlet distribution. ",
          sep = "", collapse = ""
        )
      )
    }

  }

  guide_data$message_to_user <- msg

  return(guide_data)

}

