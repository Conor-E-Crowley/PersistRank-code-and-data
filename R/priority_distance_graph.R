# Population-unit graph repair after targeted distance predicates.
#
# Only edges incident to the current recheck set may change. Surviving
# components preserve canonical node/adjacency order, strict PU-area `>`
# semantics, and the established PU-ID allocation sequence.




# Rebuild one PU after distance-invalid edges are removed.

rebuild_pu_after_edge_filter <- function(
  pu_graph,              # filtered CSR graph for one PU
  patch_area_by_patch,   # named numeric vector: patch_id -> area
  pu_area_threshold,     # species PU-area threshold
  next_available_pu_id   # next available PU ID if splitting occurs
) {
  pu_area_threshold <- validate_area_threshold(
    pu_area_threshold,
    "pu_area_threshold"
  )

  # Read graph node -> patch IDs.
  id2patch <- as.integer(pu_graph$id2patch)

  # Read CSR row pointer.
  row_ptr <- as.integer(pu_graph$row_ptr)

  # Read CSR column index.
  col_idx <- as.integer(pu_graph$col_idx)

  # Count graph nodes.
  n_nodes <- length(id2patch)

  # Return no surviving graphs if the PU has no nodes.
  if (n_nodes == 0L) {
    return(list(
      surviving_pu_graphs = list(),
      dropped_patch_ids = integer(0L),
      next_available_pu_id = next_available_pu_id
    ))
  }

  # Handle one-node graph without component search.
  if (n_nodes == 1L) {
    # Read the single patch area.
    patch_area <- as.numeric(patch_area_by_patch[as.character(id2patch[1L])])

    # Drop the single patch if its area is missing or below threshold.
    if (!area_exceeds_threshold(patch_area, pu_area_threshold)) {
      return(list(
        surviving_pu_graphs = list(),
        dropped_patch_ids = id2patch,
        next_available_pu_id = next_available_pu_id
      ))
    }

    # Otherwise return a one-node surviving PU graph.
    return(list(
      surviving_pu_graphs = list(list(
        pu_id = pu_graph$pu_id,
        id2patch = id2patch,
        row_ptr = c(0L, 0L),
        col_idx = integer(0L)
      )),
      dropped_patch_ids = integer(0L),
      next_available_pu_id = next_available_pu_id
    ))
  }

  # Label connected components in the filtered graph.
  component_id_by_node <- find_csr_components(row_ptr, col_idx)

  # Read patch area for each graph node.
  patch_area_by_node <- as.numeric(patch_area_by_patch[as.character(id2patch)])

  # Sum patch area by connected component.
  component_area <- as.numeric(
    rowsum(
      patch_area_by_node,
      component_id_by_node,
      reorder = FALSE
    )
  )

  # Components strictly above threshold survive.
  surviving_component_labels <- which(
    area_exceeds_threshold(component_area, pu_area_threshold)
  )

  # Components at or below threshold are dropped.
  dropped_component_labels <- which(
    !area_exceeds_threshold(component_area, pu_area_threshold)
  )

  # Allocate surviving graph list.
  surviving_pu_graphs <- vector("list", length(surviving_component_labels))

  # Assign PU IDs to surviving components.
  if (length(surviving_component_labels) > 0L) {
    # If only one component survives, keep the original PU ID.
    if (length(surviving_component_labels) == 1L) {
      surviving_pu_ids <- pu_graph$pu_id
    } else {
      # First surviving component keeps the old ID; others receive new IDs.
      surviving_pu_ids <- c(
        pu_graph$pu_id,
        next_available_pu_id + seq_len(length(surviving_component_labels) - 1L)
      )

      # Advance the next available PU ID.
      next_available_pu_id <- surviving_pu_ids[length(surviving_pu_ids)]
    }

    # Build one CSR subgraph per surviving component.
    for (j in seq_along(surviving_component_labels)) {
      # Read this surviving component label.
      component_label <- surviving_component_labels[j]

      # Identify original node indices in this component.
      component_node_indices <- which(component_id_by_node == component_label)

      # Build CSR subgraph for the component.
      component_subgraph <- build_csr_subgraph(
        row_ptr = row_ptr,
        col_idx = col_idx,
        kept_node_indices = component_node_indices
      )

      # Store this surviving PU graph.
      surviving_pu_graphs[[j]] <- list(
        pu_id = as.integer(surviving_pu_ids[j]),
        id2patch = as.integer(id2patch[component_node_indices]),
        row_ptr = as.integer(component_subgraph$row_ptr),
        col_idx = as.integer(component_subgraph$col_idx)
      )
    }
  }

  # Convert dropped component nodes to dropped patch IDs.
  dropped_patch_ids <- id2patch[
    component_id_by_node %in% dropped_component_labels
  ]

  # Return surviving graphs, dropped patches, and next available PU ID.
  list(
    surviving_pu_graphs = surviving_pu_graphs,
    dropped_patch_ids = as.integer(dropped_patch_ids),
    next_available_pu_id = as.integer(next_available_pu_id)
  )
}
# The supplied sparse predicate contains only edges incident to patches selected
# for rechecking. Unrelated existing edges are retained without another spatial
# predicate.
distance_predicate_rows <- function(distance_predicate_lookup) {
  if (all(c("patch_u", "patch_v", "edge_valid") %in%
          names(distance_predicate_lookup))) {
    return(data.frame(
      patch_u = as.integer(distance_predicate_lookup$patch_u),
      patch_v = as.integer(distance_predicate_lookup$patch_v),
      edge_valid = as.logical(distance_predicate_lookup$edge_valid)
    ))
  }
  values <- distance_predicate_lookup$edge_valid_by_key
  keys <- names(values)
  if (!length(values)) {
    return(data.frame(
      patch_u = integer(), patch_v = integer(), edge_valid = logical()
    ))
  }
  endpoints <- strsplit(keys, "|", fixed = TRUE)
  data.frame(
    patch_u = as.integer(vapply(endpoints, `[[`, character(1L), 1L)),
    patch_v = as.integer(vapply(endpoints, `[[`, character(1L), 2L)),
    edge_valid = as.logical(values)
  )
}

filter_distance_invalid_edges_for_pu <- function(
  pu_graph,
  distance_predicate_lookup,
  return_diagnostics = FALSE
) {
  if (!exists("stage6_filter_distance_graph_cpp", mode = "function")) {
    stop("The compiled distance-filter kernel is not loaded.", call. = FALSE)
  }
  id2patch <- get_pu_graph_id2patch(pu_graph)
  predicate <- distance_predicate_rows(distance_predicate_lookup)
  in_graph <- predicate$patch_u %in% id2patch & predicate$patch_v %in% id2patch
  predicate <- predicate[in_graph, , drop = FALSE]
  invalid <- predicate[!predicate$edge_valid, , drop = FALSE]
  lazy <- is_lazy_added_node_overlay_graph(pu_graph)
  input_adjacency <- if (lazy) {
    length(pu_graph$base_graph$col_idx) + 2L * nrow(pu_graph$overlay_edges)
  } else {
    length(pu_graph$col_idx)
  }

  # A standard CSR graph is already the exact required output when every
  # rechecked edge remains valid. Avoid scanning or copying its adjacency.
  if (!lazy && !nrow(invalid)) {
    graph <- list(
      pu_id = as.integer(pu_graph$pu_id),
      id2patch = as.integer(pu_graph$id2patch),
      row_ptr = as.integer(pu_graph$row_ptr),
      col_idx = as.integer(pu_graph$col_idx)
    )
    diagnostics <- c(
      invalid_undirected_edges = 0,
      input_adjacency_entries = input_adjacency,
      output_adjacency_entries = length(graph$col_idx),
      adjacency_entries_scanned = 0,
      adjacency_entries_avoided = input_adjacency,
      unchanged = 1,
      lazy_materialized = 0
    )
    return(if (isTRUE(return_diagnostics)) {
      list(graph = graph, diagnostics = diagnostics)
    } else graph)
  }

  if (lazy) {
    base_graph <- pu_graph$base_graph
    base_node_count <- as.integer(pu_graph$base_node_count)
    overlay <- pu_graph$overlay_edges
    overlay_u <- if (is.null(overlay) || !nrow(overlay)) integer() else as.integer(overlay$u)
    overlay_v <- if (is.null(overlay) || !nrow(overlay)) integer() else as.integer(overlay$v)
  } else {
    base_graph <- pu_graph
    base_node_count <- length(id2patch)
    overlay_u <- overlay_v <- integer()
  }
  result <- stage6_filter_distance_graph_cpp(
    id2patch = as.integer(id2patch),
    base_node_count = as.integer(base_node_count),
    base_row_ptr = as.integer(base_graph$row_ptr),
    base_col_idx = as.integer(base_graph$col_idx),
    overlay_u = overlay_u,
    overlay_v = overlay_v,
    invalid_patch_u = as.integer(invalid$patch_u),
    invalid_patch_v = as.integer(invalid$patch_v),
    pu_id = as.integer(pu_graph$pu_id)
  )
  graph <- result[c("pu_id", "id2patch", "row_ptr", "col_idx")]
  graph$pu_id <- as.integer(graph$pu_id)
  graph$id2patch <- as.integer(graph$id2patch)
  graph$row_ptr <- as.integer(graph$row_ptr)
  graph$col_idx <- as.integer(graph$col_idx)
  diagnostics <- c(
    invalid_undirected_edges = as.numeric(result$invalid_undirected_edges),
    input_adjacency_entries = as.numeric(result$input_adjacency_entries),
    output_adjacency_entries = as.numeric(result$output_adjacency_entries),
    adjacency_entries_scanned = as.numeric(result$input_adjacency_entries),
    adjacency_entries_avoided = 0,
    unchanged = 0,
    lazy_materialized = as.numeric(lazy)
  )
  if (isTRUE(return_diagnostics)) {
    list(graph = graph, diagnostics = diagnostics)
  } else graph
}
