# Deterministic connected-component operations for zero-based CSR graphs.
#
# Graph nodes follow id2patch order and each adjacency row is sorted. Component
# discovery uses ascending-node breadth-first traversal; induced graphs retain
# that ordering so patch/PU IDs and serialized graph lists remain reproducible.


find_csr_components <- function(row_ptr, col_idx, active_nodes = NULL) {
  row_ptr <- as.integer(row_ptr)
  col_idx <- as.integer(col_idx)
  n_nodes <- length(row_ptr) - 1L

  if (is.null(active_nodes)) {
    active_nodes <- rep(TRUE, n_nodes)
  } else {
    active_nodes <- as.logical(active_nodes)
  }

  component_id_by_node <- integer(n_nodes)
  queue <- integer(n_nodes)
  next_component_id <- 0L

  for (start_node in which(active_nodes)) {
    if (component_id_by_node[start_node] != 0L) next

    next_component_id <- next_component_id + 1L
    component_id_by_node[start_node] <- next_component_id
    queue_head <- 1L
    queue_tail <- 1L
    queue[1L] <- start_node

    while (queue_head <= queue_tail) {
      current_node <- queue[queue_head]
      queue_head <- queue_head + 1L

      first_edge_index <- row_ptr[current_node] + 1L
      last_edge_index <- row_ptr[current_node + 1L]
      if (last_edge_index < first_edge_index) next

      for (edge_index in first_edge_index:last_edge_index) {
        neighbor_node <- col_idx[edge_index]
        if (active_nodes[neighbor_node] &&
            component_id_by_node[neighbor_node] == 0L) {
          component_id_by_node[neighbor_node] <- next_component_id
          queue_tail <- queue_tail + 1L
          queue[queue_tail] <- neighbor_node
        }
      }
    }
  }

  component_id_by_node
}

component_areas <- function(component_id_by_node, patch_area) {
  component_nodes <- component_id_by_node > 0L
  if (!any(component_nodes)) return(numeric(0L))

  as.numeric(rowsum(
    as.numeric(patch_area)[component_nodes],
    component_id_by_node[component_nodes],
    reorder = FALSE
  ))
}

build_csr_subgraph <- function(row_ptr, col_idx, kept_node_indices) {
  row_ptr <- as.integer(row_ptr)
  col_idx <- as.integer(col_idx)
  kept_node_indices <- as.integer(kept_node_indices)

  n_original_nodes <- length(row_ptr) - 1L
  n_kept_nodes <- length(kept_node_indices)
  old_to_new_node_id <- integer(n_original_nodes)
  old_to_new_node_id[kept_node_indices] <- seq_len(n_kept_nodes)

  new_row_ptr <- integer(n_kept_nodes + 1L)

  for (new_node_index in seq_len(n_kept_nodes)) {
    original_node_index <- kept_node_indices[new_node_index]
    first_edge_index <- row_ptr[original_node_index] + 1L
    last_edge_index <- row_ptr[original_node_index + 1L]

    if (last_edge_index >= first_edge_index) {
      original_neighbors <- col_idx[first_edge_index:last_edge_index]
      remapped_neighbors <- old_to_new_node_id[original_neighbors]
      kept_neighbors <- remapped_neighbors[remapped_neighbors > 0L]
      new_row_ptr[new_node_index + 1L] <-
        new_row_ptr[new_node_index] + length(kept_neighbors)
    } else {
      new_row_ptr[new_node_index + 1L] <- new_row_ptr[new_node_index]
    }
  }

  new_col_idx <- integer(new_row_ptr[n_kept_nodes + 1L])
  write_position <- 1L

  for (new_node_index in seq_len(n_kept_nodes)) {
    original_node_index <- kept_node_indices[new_node_index]
    first_edge_index <- row_ptr[original_node_index] + 1L
    last_edge_index <- row_ptr[original_node_index + 1L]

    if (last_edge_index >= first_edge_index) {
      original_neighbors <- col_idx[first_edge_index:last_edge_index]
      remapped_neighbors <- old_to_new_node_id[original_neighbors]
      kept_neighbors <- remapped_neighbors[remapped_neighbors > 0L]
      n_kept_neighbors <- length(kept_neighbors)

      if (n_kept_neighbors > 0L) {
        new_col_idx[write_position:(write_position + n_kept_neighbors - 1L)] <-
          kept_neighbors
        write_position <- write_position + n_kept_neighbors
      }
    }
  }

  list(
    row_ptr = as.integer(new_row_ptr),
    col_idx = as.integer(new_col_idx)
  )
}

find_alive_components <- function(row_ptr, col_idx, alive_nodes, patch_area) {
  component_id_by_node <- find_csr_components(row_ptr, col_idx, active_nodes = alive_nodes)
  list(
    component_id_by_node = component_id_by_node,
    component_area = component_areas(component_id_by_node, patch_area)
  )
}
