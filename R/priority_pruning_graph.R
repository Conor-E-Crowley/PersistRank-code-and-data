# Deterministic PU repair after patch loss during pruning.
#
# Alive graph nodes are traversed in ascending order, component areas are
# summed in node order, and strict PU-area `>` decisions preserve the original
# PU ID for the first survivor before allocating new IDs sequentially.




# Rebuild one PU after some of its patches have died.

rebuild_pu_after_patch_loss <- function(
  pu_graph,
  alive_nodes,
  patch_area,
  pu_area_threshold,
  next_available_pu_id
) {
  pu_area_threshold <- validate_area_threshold(
    pu_area_threshold,
    "pu_area_threshold"
  )

  # Count currently alive graph nodes.
  n_alive_nodes <- sum(alive_nodes)

  # If no nodes survive, return no surviving PU graphs.
  if (n_alive_nodes == 0L) {
    return(list(
      surviving_pu_graphs = list(),
      dropped_patch_ids = integer(0L),
      next_available_pu_id = next_available_pu_id
    ))
  }

  # If exactly one patch survives, handle it without BFS.
  if (n_alive_nodes == 1L) {
    # Read the surviving node index.
    surviving_node_index <- which(alive_nodes)[1L]

    # Read the surviving patch area.
    surviving_patch_area <- patch_area[surviving_node_index]

    # If the single surviving patch is below the PU threshold, drop it.
    if (!area_exceeds_threshold(surviving_patch_area, pu_area_threshold)) {
      return(list(
        surviving_pu_graphs = list(),
        dropped_patch_ids = pu_graph$id2patch[surviving_node_index],
        next_available_pu_id = next_available_pu_id
      ))
    }

    # Build a one-node PU graph with no edges.
    one_node_graph <- list(
      pu_id = pu_graph$pu_id,
      id2patch = pu_graph$id2patch[surviving_node_index],
      row_ptr = c(0L, 0L),
      col_idx = integer(0L)
    )

    # Return the one-node surviving graph.
    return(list(
      surviving_pu_graphs = list(one_node_graph),
      dropped_patch_ids = integer(0L),
      next_available_pu_id = next_available_pu_id
    ))
  }

  # Find connected components among surviving graph nodes.
  component_info <- find_alive_components(
    row_ptr = pu_graph$row_ptr,
    col_idx = pu_graph$col_idx,
    alive_nodes = alive_nodes,
    patch_area = patch_area
  )

  # Read node -> component labels.
  component_id_by_node <- component_info$component_id_by_node

  # Read component areas.
  component_area <- component_info$component_area

  # Keep component labels whose area is above the PU threshold.
  surviving_component_labels <- which(
    area_exceeds_threshold(component_area, pu_area_threshold)
  )

  # Identify component labels whose area is below the PU threshold.
  dropped_component_labels <- which(
    !area_exceeds_threshold(component_area, pu_area_threshold)
  )

  # Allocate list for surviving component graphs.
  surviving_pu_graphs <- vector("list", length(surviving_component_labels))

  # Assign PU IDs to surviving components if any survive.
  if (length(surviving_component_labels) > 0L) {
    # If only one component survives, preserve the original PU ID.
    if (length(surviving_component_labels) == 1L) {
      # Reuse the original PU ID.
      surviving_pu_ids <- pu_graph$pu_id
    } else {
      # Give the first surviving component the original PU ID and others new IDs.
      surviving_pu_ids <- c(
        pu_graph$pu_id,
        next_available_pu_id + seq_len(length(surviving_component_labels) - 1L)
      )

      # Advance the next available PU ID to the last assigned ID.
      next_available_pu_id <- surviving_pu_ids[length(surviving_pu_ids)]
    }

    # Build one standalone graph per surviving component.
    for (j in seq_along(surviving_component_labels)) {
      # Read this surviving component label.
      component_label <- surviving_component_labels[j]

      # Identify original node indices in this component.
      component_node_indices <- which(component_id_by_node == component_label)

      # Build the induced CSR subgraph for this component.
      component_subgraph <- build_csr_subgraph(
        row_ptr = pu_graph$row_ptr,
        col_idx = pu_graph$col_idx,
        kept_node_indices = component_node_indices
      )

      # Store this component graph.
      surviving_pu_graphs[[j]] <- list(
        pu_id = surviving_pu_ids[j],
        id2patch = pu_graph$id2patch[component_node_indices],
        row_ptr = component_subgraph$row_ptr,
        col_idx = component_subgraph$col_idx
      )
    }
  }

  # Identify nodes belonging to dropped subthreshold components.
  dropped_node_indices <- which(component_id_by_node %in% dropped_component_labels)

  # Convert dropped component nodes to dropped patch IDs.
  dropped_patch_ids <- pu_graph$id2patch[dropped_node_indices]

  # Return surviving graphs, dropped patch IDs, and updated next PU ID.
  list(
    surviving_pu_graphs = surviving_pu_graphs,
    dropped_patch_ids = dropped_patch_ids,
    next_available_pu_id = next_available_pu_id
  )
}
