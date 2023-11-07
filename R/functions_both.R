

## functions that are used in both shiny and elsewhere

# returns the dirichlet parameter for a given number of weights
# for more than 20 nodes this does not work yet, but that is too many parameters anyways
get_diri_param <- function(no_param){

  if (no_param == 1 | no_param > 20) stop("Too many or few children for the Dirichlet prior.", call. = FALSE)

  values <- c(1, 0.7535, 0.6834, 0.6485, 0.6274,
              0.6132, 0.6029, 0.5951, 0.589, 0.5842,
              0.5801, 0.5768, 0.5739, 0.5715, 0.5693,
              0.5675, 0.5658, 0.5643, 0.563)

  return(values[no_param-1])

}

# for a given node, find which leaf nodes that contributes to the variance in that node
get_leaf_nodes <- function(node_data, node_id){

  orig_nodes <- get_original_nodes(node_data)$id

  # then the node is already a leaf-node and has no children, unless it is a top node (level == 0 is top node)
  # TODO: top node is no longer a part of the original_node-data
  #if ((node_id %in% orig_nodes) && node_data$nodes$level[node_data$nodes$id == node_id] != 0) {
  if ((node_id %in% orig_nodes)) {
    return(node_id)
  } else {

    # the children of the node
    kids <- node_data$edges$to[node_data$edges$from == node_id]

    orig_kids <- kids[(kids %in% orig_nodes)]
    not_orig_kids <- kids[!(kids %in% orig_nodes)]

    # for each of the children that are not original nodes, call the function again
    for (id in not_orig_kids){
      orig_kids <- c(orig_kids, get_leaf_nodes(node_data, id))
    }

    return(orig_kids)

  }

  return(NULL)

}

# returns all nodes below this node in the tree that has this node as parent/grand-parent/...
get_nodes_below <- function(node_data, node_id){

  # if the node is a leaf node, return the node-id
  if (!(node_id %in% node_data$edges$from)) return(node_id)

  # in all other cases, find all decendents
  kids <- node_data$edges$to[node_data$edges$from == node_id]

  for (id in kids){
    kids <- c(kids, get_nodes_below(node_data, id))
  }

  return(sort(unique(kids)))

}

# returns all nodes above this node in the tree (direct parents, not uncle nodes)
get_nodes_above <- function(node_data, node_id){

  # if the node is a top node, return the node_id
  if (!(node_id %in% node_data$edges$to)) return(node_id)

  parent_id <- c(get_parent_node_id(node_data, node_id))

  parent_id <- c(parent_id, get_nodes_above(node_data, parent_id))

  return(sort(unique(parent_id)))

}

# get the name for a node with given node id
get_node_name <- function(node_data, node_id) node_data$nodes$label[node_data$nodes$id %in% node_id]

# get the id of a node with a given name
get_node_id <- function(node_data, node_name) node_data$nodes$id[node_data$nodes$label %in% node_name]

get_original_nodes <- function(node_data){

  # which nodes are the top node (may be only 1)
  top_node_id <- node_data$nodes$id[!(node_data$nodes$id %in% node_data$edges$to)]

  # identify the nodes that are not in the from-list
  not_in_from <- node_data$nodes$id[!(node_data$nodes$id %in% node_data$edges$from)]

  #return(node_data$nodes[node_data$nodes$id %in% c(top_node_id, not_in_from),])
  return(node_data$nodes[node_data$nodes$id %in% not_in_from,])

  # return(node_data$orig_nodedata)

}

# get id of parent node
# if the node has no edges to it, return the id of the input node (node_id)
get_parent_node_id <- function(node_data, node_id){
  if (sum(node_data$edges$to == node_id) == 0) return(node_id)
  return(node_data$edges$from[node_data$edges$to == node_id])
}

# get the ids of all split nodes
get_split_ids <- function(node_data) return(unique(node_data$edges$from))

# which index in the prior_data-list corresponds to a given node id
# get_split_number <- function(prior_data, node_id) {print("use get_prior_number instead!"); which(sapply(prior_data, function(x) x$id) == node_id)}
get_split_number <- function(prior_data, node_id) which(sapply(prior_data, function(x) x$id) == node_id)

# which index in the prior_data-list corresponds to a given node id (split, top node or CW)
get_prior_number <- function(prior_data, node_id) which(sapply(prior_data, function(x) x$id) == node_id)

# which nodes are in the tree where the node with node_id belongs/is attached
nodes_in_tree <- function(node_data, node_id){

  if (sum(node_data$nodes$top_node) == 0) return(0)

  # find all nodes above in the tree
  parents <- get_nodes_above(node_data, node_id)

  # for each parent, find the child nodes
  kids <- c()
  for (id in parents){
    kids <- c(kids, get_nodes_below(node_data, id))
  }

  return(sort(unique(c(kids, parents, node_id))))

}

# which nodes are the first generation children of this node
first_gen_children <- function(node_data, node_id) node_data$edges$to[node_data$edges$from == node_id]


