# Provisional graph updates after patch fragmentation.
#
# Additions-only changes use a lazy overlay over the immutable base CSR graph.
# Hard cases use the compiled provisional-graph kernel. Fragment siblings
# inherit origin neighbours until the distance stage removes invalid edges.


# Make a lazy overlay graph for additions-only fragmentation.
#
# This is the fast path for cases where fragmentation added new patch nodes
# but did not remove any old patch nodes.
#
# Adding nodes cannot disconnect an already connected graph, so this helper
# preserves the old CSR graph and stores only the provisional edges involving
# newly added fragment nodes. The distance stage can then consume this lazy
# graph and materialize or update it as needed.
#

make_lazy_added_node_overlay_pu_graph <- function(
  current_pu_patch_rows,
  old_pu_graph,
  new_patch_nodes
) {
  # A lazy overlay graph must wrap an existing old graph.
  if (is.null(old_pu_graph)) {
    stop("make_lazy_added_node_overlay_pu_graph() received NULL old_pu_graph.")
  }

  # Read current patch IDs.
  current_patch_ids <- as.integer(current_pu_patch_rows$patch_id)

  # Read current origin patch IDs.
  current_origin_patch_ids <- as.integer(current_pu_patch_rows$origin_patch_id)

  # Read old graph patch IDs.
  old_patch_ids <- as.integer(old_pu_graph$id2patch)

  # Read old CSR row pointer.
  old_row_ptr <- as.integer(old_pu_graph$row_ptr)

  # Read old CSR column index.
  old_col_idx <- as.integer(old_pu_graph$col_idx)

  # Stop if current patch IDs are duplicated.
  if (anyDuplicated(current_patch_ids)) {
    duplicated_patch_ids <- unique(current_patch_ids[duplicated(current_patch_ids)])

    stop(
      "Duplicate patch_id values found inside current_pu_patch_rows: ",
      paste(utils::head(duplicated_patch_ids, 20L), collapse = ", "),
      if (length(duplicated_patch_ids) > 20L) " ..." else ""
    )
  }

  # Stop if required current patch identity fields contain missing values.
  if (any(is.na(current_patch_ids)) || any(is.na(current_origin_patch_ids))) {
    stop("current_pu_patch_rows contains NA patch_id or origin_patch_id values.")
  }

  # Confirm that no old patch nodes are missing.
  missing_old_patch_ids <- setdiff(old_patch_ids, current_patch_ids)

  # Stop if this additions-only helper was called after old nodes disappeared.
  if (length(missing_old_patch_ids)) {
    stop(
      "Lazy overlay graph was requested even though old patch nodes are missing: ",
      paste(utils::head(missing_old_patch_ids, 20L), collapse = ", "),
      if (length(missing_old_patch_ids) > 20L) " ..." else ""
    )
  }

  # Preserve current patch-row order for added patch nodes.
  added_patch_ids <- unique(as.integer(
    current_patch_ids[current_patch_ids %in% as.integer(new_patch_nodes)]
  ))

  # Stop if no added nodes were supplied.
  if (!length(added_patch_ids)) {
    stop("make_lazy_added_node_overlay_pu_graph() was called with no added patch nodes.")
  }

  # Count original graph nodes.
  base_node_count <- length(old_patch_ids)

  # Build old patch_id -> old graph row lookup.
  old_row_by_patch_id <- stats::setNames(
    seq_along(old_patch_ids),
    as.character(old_patch_ids)
  )

  # Build added patch_id -> new graph row lookup.
  added_row_by_patch_id <- stats::setNames(
    base_node_count + seq_along(added_patch_ids),
    as.character(added_patch_ids)
  )

  # Combine old and added patch row lookups.
  full_row_by_patch_id <- c(
    old_row_by_patch_id,
    added_row_by_patch_id
  )

  # Build patch_id -> origin_patch_id lookup for current rows.
  current_origin_by_patch_id <- stats::setNames(
    current_origin_patch_ids,
    as.character(current_patch_ids)
  )

  # Split current patch IDs by origin patch.
  current_patches_by_origin <- split(
    current_patch_ids,
    current_origin_patch_ids
  )

  # Allocate overlay-edge table list.
  overlay_edge_tables <- vector("list", length(added_patch_ids))

  # Initialize overlay-edge table write index.
  overlay_edge_table_index <- 0L

  # Loop over newly added patch nodes.
  for (added_patch_id in added_patch_ids) {
    # Convert added patch ID to character key.
    added_patch_key <- as.character(added_patch_id)

    # Read graph row for this added patch.
    added_row <- as.integer(added_row_by_patch_id[added_patch_key])

    # Read the origin patch for this added patch.
    origin_patch_id <- as.integer(current_origin_by_patch_id[added_patch_key])

    # Stop if origin lookup failed.
    if (is.na(origin_patch_id)) {
      stop("Could not find origin_patch_id for added patch ", added_patch_id, ".")
    }

    # Read the old graph row corresponding to the origin patch.
    origin_row <- as.integer(old_row_by_patch_id[as.character(origin_patch_id)])

    # Stop if the origin patch is not present in the old graph.
    if (is.na(origin_row)) {
      stop(
        "Added patch ",
        added_patch_id,
        " has origin_patch_id ",
        origin_patch_id,
        ", but that origin is not present in the old PU graph."
      )
    }

    # Compute first old CSR edge position for the origin patch.
    first_edge_index <- old_row_ptr[origin_row] + 1L

    # Compute last old CSR edge position for the origin patch.
    last_edge_index <- old_row_ptr[origin_row + 1L]

    # Initialize neighbor origin patch IDs.
    neighbor_origin_patch_ids <- integer(0L)

    # If the origin node has old graph neighbors, read their origin patch IDs.
    if (last_edge_index >= first_edge_index) {
      # Read old neighboring graph rows.
      neighbor_rows <- as.integer(old_col_idx[first_edge_index:last_edge_index])

      # Convert old neighboring graph rows to old patch IDs.
      neighbor_origin_patch_ids <- as.integer(old_patch_ids[neighbor_rows])
    }

    # Added nodes provisionally connect to same-origin fragments and old-neighbor descendants.
    target_origin_ids <- unique(as.integer(c(
      origin_patch_id,
      neighbor_origin_patch_ids
    )))

    # Keep target origins that currently have descendants.
    target_origin_names <- intersect(
      as.character(target_origin_ids),
      names(current_patches_by_origin)
    )

    # Read all current patch IDs descended from those target origins.
    target_patch_ids <- unique(as.integer(unlist(
      current_patches_by_origin[target_origin_names],
      use.names = FALSE
    )))

    # Do not connect the added patch to itself.
    target_patch_ids <- target_patch_ids[target_patch_ids != added_patch_id]

    # Skip this added patch if it has no provisional targets.
    if (!length(target_patch_ids)) {
      next
    }

    # Map target patch IDs to graph rows.
    target_rows <- as.integer(
      full_row_by_patch_id[as.character(target_patch_ids)]
    )

    # Stop if any target row failed to map.
    if (any(is.na(target_rows))) {
      stop(
        "Failed to map one or more target patches to graph rows while ",
        "building lazy overlay edges for added patch ",
        added_patch_id,
        "."
      )
    }

    # Canonicalize overlay edges by graph-row ID.
    edge_low <- pmin(added_row, target_rows)

    # Canonicalize overlay edges by graph-row ID.
    edge_high <- pmax(added_row, target_rows)

    # Advance overlay-edge table write index.
    overlay_edge_table_index <- overlay_edge_table_index + 1L

    # Store overlay edges for this added patch.
    overlay_edge_tables[[overlay_edge_table_index]] <- data.table::data.table(
      u = as.integer(edge_low),
      v = as.integer(edge_high)
    )
  }

  # Combine overlay edge tables if any were created.
  overlay_edges <- if (overlay_edge_table_index > 0L) {
    # Bind only initialized edge-table slots.
    overlay_edges_raw <- data.table::rbindlist(
      overlay_edge_tables[seq_len(overlay_edge_table_index)],
      use.names = TRUE
    )

    # Drop self-edges.
    overlay_edges_raw <- overlay_edges_raw[u != v]

    # Deduplicate overlay edges.
    overlay_edges_raw <- unique(overlay_edges_raw, by = c("u", "v"))

    # Sort overlay edges deterministically.
    data.table::setorder(overlay_edges_raw, u, v)

    # Return the overlay edge table.
    overlay_edges_raw
  } else {
    # Return an empty overlay edge table.
    data.table::data.table(
      u = integer(),
      v = integer()
    )
  }

  # Return the lazy overlay graph object.
  list(
    graph_type = "lazy_added_node_overlay",
    species = old_pu_graph$species,
    pu_id = as.integer(old_pu_graph$pu_id),
    id2patch = as.integer(c(old_patch_ids, added_patch_ids)),
    base_graph = old_pu_graph,
    base_node_count = as.integer(base_node_count),
    added_patch_ids = as.integer(added_patch_ids),
    added_origin_patch_ids = as.integer(
      current_origin_by_patch_id[as.character(added_patch_ids)]
    ),
    overlay_edges = overlay_edges
  )
}


# Build a provisional fragmentation graph with the validated compiled kernel.
# The allocation-heavy R oracle is confined to the test helper.

build_provisional_pu_graph <- function(
  current_pu_patch_rows,
  old_pu_graph,
  return_diagnostics = FALSE
) {
  n_current_patches <- nrow(current_pu_patch_rows)
  if (n_current_patches == 0L) {
    graph <- list(
      pu_id = integer(0L),
      id2patch = integer(0L),
      row_ptr = 0L,
      col_idx = integer(0L)
    )
    if (!isTRUE(return_diagnostics)) return(graph)
    return(list(
      graph = graph,
      diagnostics = c(
        candidate_undirected_edges = 0,
        unique_undirected_edges = 0,
        output_adjacency_entries = 0
      )
    ))
  }
  if (is.null(old_pu_graph)) {
    stop("build_provisional_pu_graph() received NULL old_pu_graph.")
  }
  if (!exists(
    "stage6_build_provisional_fragment_graph_cpp",
    mode = "function",
    inherits = TRUE
  )) {
    stop(
      "The compiled provisional-fragment graph kernel is not loaded.",
      call. = FALSE
    )
  }

  current_patch_ids <- as.integer(current_pu_patch_rows$patch_id)
  current_origin_patch_ids <- as.integer(current_pu_patch_rows$origin_patch_id)
  pu_id_value <- as.integer(current_pu_patch_rows$pu_id[1L])
  kernel_result <- stage6_build_provisional_fragment_graph_cpp(
    current_patch_ids = current_patch_ids,
    current_origin_patch_ids = current_origin_patch_ids,
    old_patch_ids = as.integer(old_pu_graph$id2patch),
    old_row_ptr = as.integer(old_pu_graph$row_ptr),
    old_col_idx = as.integer(old_pu_graph$col_idx),
    pu_id = pu_id_value
  )
  graph <- kernel_result[c("pu_id", "id2patch", "row_ptr", "col_idx")]
  graph$pu_id <- as.integer(graph$pu_id)
  graph$id2patch <- as.integer(graph$id2patch)
  graph$row_ptr <- as.integer(graph$row_ptr)
  graph$col_idx <- as.integer(graph$col_idx)
  if (!isTRUE(return_diagnostics)) return(graph)

  list(
    graph = graph,
    diagnostics = c(
      candidate_undirected_edges = as.numeric(
        kernel_result$candidate_undirected_edges
      ),
      unique_undirected_edges = as.numeric(
        kernel_result$unique_undirected_edges
      ),
      output_adjacency_entries = as.numeric(
        kernel_result$output_adjacency_entries
      )
    )
  )
}


# Return compact graph-storage counts without materializing edges. Lazy
# overlays retain their base graph and store one row per undirected overlay
# edge, so each overlay row contributes two adjacency entries.
fragmentation_graph_storage_counts <- function(pu_graph) {
  if (is.null(pu_graph) || is.null(pu_graph$id2patch)) {
    stop("Cannot profile an invalid PU graph.", call. = FALSE)
  }
  node_count <- length(pu_graph$id2patch)
  if (identical(pu_graph$graph_type, "lazy_added_node_overlay")) {
    base <- fragmentation_graph_storage_counts(pu_graph$base_graph)
    overlay_edges <- pu_graph$overlay_edges
    overlay_count <- if (is.null(overlay_edges)) 0 else nrow(overlay_edges)
    return(c(
      nodes = as.numeric(node_count),
      adjacency_entries = as.numeric(base[["adjacency_entries"]]) +
        2 * as.numeric(overlay_count)
    ))
  }
  adjacency <- if (is.null(pu_graph$col_idx)) integer() else pu_graph$col_idx
  c(
    nodes = as.numeric(node_count),
    adjacency_entries = as.numeric(length(adjacency))
  )
}
