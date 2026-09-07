# Distance-connectivity stage driver.
#
# Affected species are processed sequentially. Candidate cells and graph edges
# come from compact indexes and targeted CSR traversal; commits update mutable
# graph, patch, compact-index, and dense-cell state only after local validation.
# The returned profiles contain counters already collected by scientific work.

# The fragmentation stage identifies patches whose dispersal links may have
# changed. This driver rechecks only those local candidate patches, removes
# links beyond the species dispersal threshold, and rebuilds affected
# persistence units.

log_distance_connectivity_stage <- function(
  stage_index,
  touched_species_count,
  input_recheck_patch_count,
  removed_cell_count,
  remaining_patch_count,
  remaining_pu_count,
  remaining_alive_cell_count
) {
  runtime_log_event(
    "distance_done",
    stage = as.integer(stage_index),
    touched_species = as.integer(touched_species_count),
    input_recheck_patches = as.integer(input_recheck_patch_count),
    removed_cells = as.integer(removed_cell_count),
    remaining_patches = as.integer(remaining_patch_count),
    remaining_pus = as.integer(remaining_pu_count),
    remaining_alive = as.integer(remaining_alive_cell_count)
  )

  invisible(NULL)
}

empty_distance_timing <- function() {
  stats::setNames(
    numeric(17L),
    c(
      "candidate_edge_setup_seconds",
      "candidate_index_extract_seconds",
      "candidate_extent_seconds",
      "local_raster_seconds",
      "polygonize_seconds",
      "sf_conversion_seconds",
      "geometry_validation_seconds",
      "geometry_repair_seconds",
      "predicate_spatial_seconds",
      "predicate_mapping_seconds",
      "pu_index_seconds",
      "pu_repair_seconds",
      "graph_commit_seconds",
      "compact_index_seconds",
      "cell_state_seconds",
      "state_commit_seconds",
      "other_seconds"
    )
  )
}

empty_distance_geometry_counts <- function() {
  stats::setNames(
    numeric(7L),
    c(
      "candidate_patches",
      "candidate_cells",
      "local_raster_cells",
      "max_expansion_ratio",
      "geometry_features_checked",
      "invalid_geometry_features",
      "geometry_repair_invoked"
    )
  )
}

empty_distance_counts <- function() {
  stats::setNames(numeric(14L), c(
    "affected_pus", "edge_filtered_pus", "unchanged_distance_pus",
    "compiled_filtered_pus", "lazy_materialized_pus",
    "invalid_edges_removed", "adjacency_entries_scanned",
    "adjacency_entries_avoided", "graphs_inspected",
    "graph_nodes_scanned", "recheck_rows",
    "csr_adjacency_entries_visited", "overlay_edges_scanned",
    "unique_recheck_edges"
  ))
}

empty_distance_species_workload <- function() {
  stats::setNames(numeric(29L), c(
    "recheck_patches", "candidate_patches", "candidate_cells",
    "local_raster_cells", "raster_expansion_ratio", "rechecked_edges",
    "affected_pus", "edge_filtered_pus", "input_graph_nodes",
    "input_graph_adjacency_entries", "output_graph_nodes",
    "output_graph_adjacency_entries", "surviving_components",
    "dropped_patches", "unchanged_distance_pus", "compiled_filtered_pus",
    "lazy_materialized_pus", "invalid_edges_removed",
    "adjacency_entries_scanned", "adjacency_entries_avoided",
    "geometry_features_checked", "invalid_geometry_features",
    "geometry_repair_invoked", "graphs_inspected", "graph_nodes_scanned",
    "recheck_rows", "csr_adjacency_entries_visited",
    "overlay_edges_scanned", "unique_recheck_edges"
  ))
}

empty_distance_species_timing <- function() {
  stats::setNames(numeric(19L), c(
    "candidate_edge_setup_seconds", "candidate_index_extract_seconds",
    "candidate_extent_seconds",
    "local_raster_seconds", "polygonize_seconds", "sf_conversion_seconds",
    "geometry_validation_seconds", "geometry_repair_seconds",
    "predicate_spatial_seconds", "predicate_mapping_seconds",
    "pu_index_seconds",
    "edge_filter_seconds", "component_rebuild_seconds", "pu_assembly_seconds",
    "graph_commit_seconds", "compact_index_seconds", "cell_state_seconds",
    "other_seconds", "elapsed_seconds"
  ))
}

new_distance_species_profile <- function() {
  list(
    workload = empty_distance_species_workload(),
    timing = empty_distance_species_timing()
  )
}

validate_distance_species_profile <- function(profile) {
  valid <- is.list(profile) &&
    is.numeric(profile$workload) &&
    identical(names(profile$workload), names(empty_distance_species_workload())) &&
    all(is.finite(profile$workload) & profile$workload >= 0) &&
    is.numeric(profile$timing) &&
    identical(names(profile$timing), names(empty_distance_species_timing())) &&
    all(is.finite(profile$timing) & profile$timing >= 0)
  if (!isTRUE(valid)) {
    stop("Distance species profile fields must be finite named nonnegative values.")
  }

  classified <- setdiff(
    names(profile$timing),
    c("other_seconds", "elapsed_seconds")
  )
  if (!isTRUE(all.equal(
    sum(profile$timing[c(classified, "other_seconds")]),
    as.numeric(profile$timing[["elapsed_seconds"]]),
    tolerance = 1e-8,
    check.attributes = FALSE
  ))) {
    stop("Distance species profile timings are not mutually exclusive.")
  }
  profile
}

finalize_distance_species_profile <- function(profile, elapsed_seconds) {
  elapsed_seconds <- as.numeric(elapsed_seconds)
  if (length(elapsed_seconds) != 1L || !is.finite(elapsed_seconds) ||
      elapsed_seconds < 0) {
    stop("Distance species elapsed time must be finite and nonnegative.")
  }
  classified <- setdiff(
    names(profile$timing),
    c("other_seconds", "elapsed_seconds")
  )
  classified_seconds <- sum(profile$timing[classified])
  profile$timing[["other_seconds"]] <- max(
    0,
    elapsed_seconds - classified_seconds
  )
  profile$timing[["elapsed_seconds"]] <-
    classified_seconds + profile$timing[["other_seconds"]]
  validate_distance_species_profile(profile)
}

new_distance_geometry_profile <- function(n_species) {
  n_species <- as.integer(n_species)
  if (length(n_species) != 1L || is.na(n_species) || n_species < 0L) {
    stop("Distance geometry profile size must be one nonnegative integer.")
  }
  matrix(
    0,
    nrow = n_species,
    ncol = 8L,
    dimnames = list(NULL, c(
      "expansion_ratio",
      "candidate_cells",
      "local_raster_cells",
      "local_raster_seconds",
      "polygonize_seconds",
      "sf_conversion_seconds",
      "geometry_validation_seconds",
      "geometry_repair_seconds"
    ))
  )
}

summarize_distance_geometry_profile <- function(profile) {
  expected_names <- colnames(new_distance_geometry_profile(0L))
  if (!is.matrix(profile) || !is.numeric(profile) ||
      !identical(colnames(profile), expected_names) ||
      any(!is.finite(profile)) || any(profile < 0)) {
    stop("Distance geometry profile values must be finite and nonnegative.")
  }
  if (!nrow(profile)) {
    return(NULL)
  }

  expansion <- profile[, "expansion_ratio"]
  local_cells <- profile[, "local_raster_cells"]
  expansion_quantiles <- as.numeric(stats::quantile(
    expansion,
    probs = c(0.5, 0.9, 0.99),
    names = FALSE,
    type = 7
  ))
  total_local_cells <- sum(local_cells)
  top_five_share <- if (total_local_cells > 0) {
    100 * sum(utils::head(sort(local_cells, decreasing = TRUE), 5L)) /
      total_local_cells
  } else {
    0
  }

  actionable_seconds <- rowSums(profile[, c(
    "local_raster_seconds",
    "polygonize_seconds",
    "sf_conversion_seconds",
    "geometry_validation_seconds",
    "geometry_repair_seconds"
  ), drop = FALSE])
  compact <- expansion < 10
  sparse <- expansion >= 10 & expansion < 100
  extreme <- expansion >= 100
  displayed_cost <- round(c(
    compact_lt10_seconds = sum(actionable_seconds[compact]),
    sparse_10_100_seconds = sum(actionable_seconds[sparse]),
    extreme_ge100_seconds = sum(actionable_seconds[extreme])
  ), 2L)

  list(
    sparsity = list(
      profiled_species = as.integer(nrow(profile)),
      expansion_p50 = round(expansion_quantiles[1L], 2L),
      expansion_p90 = round(expansion_quantiles[2L], 2L),
      expansion_p99 = round(expansion_quantiles[3L], 2L),
      expansion_max = round(max(expansion), 2L),
      species_ge10 = as.integer(sum(expansion >= 10)),
      species_ge100 = as.integer(sum(expansion >= 100)),
      species_ge1000 = as.integer(sum(expansion >= 1000)),
      max_local_raster_cells = max(local_cells),
      top5_local_raster_cells_pct = round(top_five_share, 2L)
    ),
    cost = list(
      compact_lt10_species = as.integer(sum(compact)),
      sparse_10_100_species = as.integer(sum(sparse)),
      extreme_ge100_species = as.integer(sum(extreme)),
      compact_lt10_seconds = displayed_cost[["compact_lt10_seconds"]],
      sparse_10_100_seconds = displayed_cost[["sparse_10_100_seconds"]],
      extreme_ge100_seconds = displayed_cost[["extreme_ge100_seconds"]],
      actionable_geometry_seconds = sum(displayed_cost)
    )
  )
}

log_distance_geometry_profile <- function(stage_index, profile) {
  summary <- summarize_distance_geometry_profile(profile)
  if (is.null(summary)) {
    return(invisible(NULL))
  }
  do.call(
    runtime_log_event,
    c(list(event = "distance_sparsity", stage = as.integer(stage_index)), summary$sparsity)
  )
  do.call(
    runtime_log_event,
    c(
      list(event = "distance_geometry_cost", stage = as.integer(stage_index)),
      summary$cost
    )
  )
  invisible(NULL)
}

log_distance_timing <- function(
  stage_index,
  timing_seconds,
  path_counts,
  geometry_counts = empty_distance_geometry_counts()
) {
  if (!is.numeric(timing_seconds) ||
      !identical(names(timing_seconds), names(empty_distance_timing())) ||
      any(!is.finite(timing_seconds)) || any(timing_seconds < 0)) {
    stop("Distance timings must be finite named nonnegative values.")
  }
  if (!is.numeric(path_counts) ||
      !identical(names(path_counts), names(empty_distance_counts())) ||
      any(!is.finite(path_counts)) || any(path_counts < 0) ||
      any(path_counts != floor(path_counts))) {
    stop("Distance path counts must be named nonnegative integers.")
  }
  if (!is.numeric(geometry_counts) ||
      !identical(names(geometry_counts), names(empty_distance_geometry_counts())) ||
      any(!is.finite(geometry_counts)) || any(geometry_counts < 0) ||
      any(geometry_counts[c(
        "candidate_patches", "candidate_cells", "local_raster_cells",
        "geometry_features_checked", "invalid_geometry_features",
        "geometry_repair_invoked"
      )] != floor(geometry_counts[c(
        "candidate_patches", "candidate_cells", "local_raster_cells",
        "geometry_features_checked", "invalid_geometry_features",
        "geometry_repair_invoked"
      )]))) {
    stop("Distance geometry counts must be valid named nonnegative values.")
  }
  displayed_timing <- round(timing_seconds, 2L)
  candidate_cell_seconds <- sum(displayed_timing[c(
    "candidate_edge_setup_seconds", "candidate_index_extract_seconds",
    "candidate_extent_seconds"
  )])
  predicate_seconds <- sum(displayed_timing[c(
    "predicate_spatial_seconds", "predicate_mapping_seconds"
  )])
  geometry_fields <- c(
    "candidate_edge_setup_seconds", "candidate_index_extract_seconds",
    "candidate_extent_seconds", "local_raster_seconds", "polygonize_seconds",
    "sf_conversion_seconds", "geometry_validation_seconds",
    "geometry_repair_seconds"
  )
  geometry_seconds <- sum(displayed_timing[geometry_fields])
  expansion_ratio <- if (geometry_counts[["candidate_cells"]] > 0) {
    geometry_counts[["local_raster_cells"]] / geometry_counts[["candidate_cells"]]
  } else {
    0
  }
  runtime_log_event(
    "distance_timing",
    stage = as.integer(stage_index),
    geometry_seconds = geometry_seconds,
    candidate_cell_seconds = candidate_cell_seconds,
    candidate_edge_setup_seconds = displayed_timing[["candidate_edge_setup_seconds"]],
    candidate_index_extract_seconds = displayed_timing[["candidate_index_extract_seconds"]],
    candidate_extent_seconds = displayed_timing[["candidate_extent_seconds"]],
    local_raster_seconds = displayed_timing[["local_raster_seconds"]],
    polygonize_seconds = displayed_timing[["polygonize_seconds"]],
    sf_conversion_seconds = displayed_timing[["sf_conversion_seconds"]],
    geometry_validation_seconds = displayed_timing[["geometry_validation_seconds"]],
    geometry_repair_seconds = displayed_timing[["geometry_repair_seconds"]],
    predicate_seconds = predicate_seconds,
    predicate_spatial_seconds = displayed_timing[["predicate_spatial_seconds"]],
    predicate_mapping_seconds = displayed_timing[["predicate_mapping_seconds"]],
    pu_index_seconds = displayed_timing[["pu_index_seconds"]],
    pu_repair_seconds = displayed_timing[["pu_repair_seconds"]],
    graph_commit_seconds = displayed_timing[["graph_commit_seconds"]],
    compact_index_seconds = displayed_timing[["compact_index_seconds"]],
    cell_state_seconds = displayed_timing[["cell_state_seconds"]],
    state_commit_seconds = displayed_timing[["state_commit_seconds"]],
    other_seconds = displayed_timing[["other_seconds"]],
    affected_pus = as.integer(path_counts[["affected_pus"]]),
    edge_filtered_pus = as.integer(path_counts[["edge_filtered_pus"]]),
    unchanged_distance_pus = as.integer(path_counts[["unchanged_distance_pus"]]),
    compiled_filtered_pus = as.integer(path_counts[["compiled_filtered_pus"]]),
    lazy_materialized_pus = as.integer(path_counts[["lazy_materialized_pus"]]),
    invalid_edges_removed = as.integer(path_counts[["invalid_edges_removed"]]),
    adjacency_entries_scanned = as.numeric(path_counts[["adjacency_entries_scanned"]]),
    adjacency_entries_avoided = as.numeric(path_counts[["adjacency_entries_avoided"]]),
    graphs_inspected = as.integer(path_counts[["graphs_inspected"]]),
    graph_nodes_scanned = as.numeric(path_counts[["graph_nodes_scanned"]]),
    recheck_rows = as.integer(path_counts[["recheck_rows"]]),
    csr_adjacency_entries_visited = as.numeric(
      path_counts[["csr_adjacency_entries_visited"]]
    ),
    overlay_edges_scanned = as.numeric(path_counts[["overlay_edges_scanned"]]),
    unique_recheck_edges = as.numeric(path_counts[["unique_recheck_edges"]]),
    candidate_patches = as.integer(geometry_counts[["candidate_patches"]]),
    candidate_cells = as.integer(geometry_counts[["candidate_cells"]]),
    local_raster_cells = as.integer(geometry_counts[["local_raster_cells"]]),
    geometry_features_checked = as.integer(
      geometry_counts[["geometry_features_checked"]]
    ),
    invalid_geometry_features = as.integer(
      geometry_counts[["invalid_geometry_features"]]
    ),
    geometry_repair_invoked = as.integer(
      geometry_counts[["geometry_repair_invoked"]]
    ),
    raster_expansion_ratio = round(expansion_ratio, 2L),
    max_species_expansion_ratio = round(geometry_counts[["max_expansion_ratio"]], 2L),
    total_seconds = sum(displayed_timing)
  )
  invisible(NULL)
}

log_distance_predicate_profile <- function(stage_index, profile) {
  expected_names <- names(empty_distance_predicate_profile())
  if (!is.numeric(profile) ||
      !identical(names(profile), expected_names) ||
      any(!is.finite(profile)) || any(profile < 0) ||
      any(profile[seq_len(11L)] != floor(profile[seq_len(11L)]))) {
    stop("Distance predicate profile must contain finite nonnegative values.")
  }
  displayed <- round(profile, 2L)
  graph_valid_pct <- if (profile[["graph_edges"]] > 0) {
    100 * profile[["valid_graph_edges"]] / profile[["graph_edges"]]
  } else {
    0
  }
  non_graph_pct <- if (profile[["unique_spatial_keys"]] > 0) {
    100 * profile[["non_graph_spatial_keys"]] /
      profile[["unique_spatial_keys"]]
  } else {
    0
  }
  runtime_log_event(
    "distance_predicate_profile",
    stage = as.integer(stage_index),
    predicate_calls = as.integer(profile[["predicate_calls"]]),
    recheck_patches = as.numeric(profile[["recheck_patches"]]),
    candidate_patches = as.numeric(profile[["candidate_patches"]]),
    candidate_pair_upper_bound = as.numeric(
      profile[["candidate_pair_upper_bound"]]
    ),
    graph_edges = as.numeric(profile[["graph_edges"]]),
    directed_predicate_hits = as.numeric(
      profile[["directed_predicate_hits"]]
    ),
    unique_spatial_keys = as.numeric(profile[["unique_spatial_keys"]]),
    self_spatial_keys = as.numeric(profile[["self_spatial_keys"]]),
    valid_graph_edges = as.numeric(profile[["valid_graph_edges"]]),
    invalid_graph_edges = as.numeric(profile[["invalid_graph_edges"]]),
    non_graph_spatial_keys = as.numeric(
      profile[["non_graph_spatial_keys"]]
    ),
    graph_edge_valid_pct = round(graph_valid_pct, 2L),
    non_graph_spatial_key_pct = round(non_graph_pct, 2L),
    predicate_spatial_seconds = displayed[["predicate_spatial_seconds"]],
    predicate_mapping_seconds = displayed[["predicate_mapping_seconds"]],
    residual_seconds = displayed[["residual_seconds"]],
    total_seconds = displayed[["total_seconds"]]
  )
  invisible(NULL)
}


# Run the distance-based connectivity stage.

run_distance_connectivity_stage <- function(
  stage_index,
  patches_requiring_distance_recheck,
  patch_table,
  pu_graphs_by_key,
  alive_species_count_by_cell,
  patch_id_by_species_env,
  patch_cell_index_by_species_env,
  species_params,
  template_raster,
  state_representation = c("dense", "compact"),
  predicate_strategy = c("matrix", "edge_pairs"),
  collect_runtime_diagnostics = FALSE,
  emit_logs = TRUE,
  progress_callback = NULL,
  completion_callback = NULL
) {
  distance_stage_started <- proc.time()[["elapsed"]]
  state_representation <- match.arg(state_representation)
  predicate_strategy <- match.arg(predicate_strategy)
  compact_state <- identical(state_representation, "compact")
  if (!is.logical(collect_runtime_diagnostics) ||
      length(collect_runtime_diagnostics) != 1L ||
      is.na(collect_runtime_diagnostics)) {
    stop("collect_runtime_diagnostics must be TRUE or FALSE.", call. = FALSE)
  }
  # -------------------------------------------------------------------
  # 1. Standardize table inputs
  # -------------------------------------------------------------------

  patch_table <- data.table::as.data.table(patch_table)

  species_params <- data.table::as.data.table(species_params)

  patches_requiring_distance_recheck <-
    data.table::as.data.table(patches_requiring_distance_recheck)

  input_recheck_patch_count <- nrow(patches_requiring_distance_recheck)


  # -------------------------------------------------------------------
  # 2. Validate package dependency
  # -------------------------------------------------------------------

  # Stop early if sf is unavailable, because local distance predicates require sf.
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("The distance stage requires the 'sf' package.")
  }


  # -------------------------------------------------------------------
  # 3. Initialize distance-stage accumulator
  # -------------------------------------------------------------------

  # Store cells that become globally empty during distance-stage updates.
  removed_cells_in_distance <- integer(0L)
  distance_timing_seconds <- empty_distance_timing()
  distance_path_counts <- empty_distance_counts()
  distance_geometry_counts <- empty_distance_geometry_counts()
  distance_predicate_profile <- empty_distance_predicate_profile()
  candidate_ids_by_species <- if (collect_runtime_diagnostics) list() else NULL
  distance_workload <- if (collect_runtime_diagnostics) stats::setNames(numeric(11L), c(
    "recheck_patches", "candidate_patches_before_filter",
    "edge_endpoint_patches", "isolated_rechecks_avoided",
    "zero_edge_species", "all_pair_candidate_upper_bound",
    "graph_edge_pairs_evaluated", "valid_edges", "invalid_edges",
    "predicate_spatial_seconds", "predicate_mapping_seconds"
  )) else NULL
  rechecked_edges <- 0L
  profiling_enabled <- is.function(completion_callback)
  if (!is.null(completion_callback) && !profiling_enabled) {
    stop("completion_callback must be NULL or a function.")
  }


  # -------------------------------------------------------------------
  # 4. Fast return if no patches were flagged for distance recheck
  # -------------------------------------------------------------------

  # If there are no recheck patches, there is no distance-stage work.
  if (nrow(patches_requiring_distance_recheck) == 0L) {
    # Count current patch-table rows.
    remaining_patch_count <- nrow(patch_table)

    # Count current unique species-PU combinations.
    remaining_pu_count <- data.table::uniqueN(patch_table, by = c("species", "pu_id"))

    # Count globally alive cells.
    remaining_alive_cell_count <- sum(alive_species_count_by_cell > 0L)

    if (isTRUE(emit_logs)) log_distance_connectivity_stage(
      stage_index = stage_index,
      touched_species_count = 0L,
      input_recheck_patch_count = 0L,
      removed_cell_count = 0L,
      remaining_patch_count = remaining_patch_count,
      remaining_pu_count = remaining_pu_count,
      remaining_alive_cell_count = remaining_alive_cell_count
    )
    if (isTRUE(emit_logs)) {
      log_distance_timing(stage_index, distance_timing_seconds, distance_path_counts)
      log_distance_predicate_profile(stage_index, distance_predicate_profile)
    }

    # Return unchanged state.
    return(list(
      patch_table = patch_table,
      pu_graphs = pu_graphs_by_key,
      alive_species_count_by_cell = alive_species_count_by_cell,
      removed_cells_in_distance = integer(0L),
      touched_species = character(),
      timing_seconds = distance_timing_seconds,
      path_counts = distance_path_counts,
      geometry_counts = distance_geometry_counts,
      rechecked_edges = rechecked_edges,
      workload = distance_workload,
      candidate_ids_by_species = candidate_ids_by_species
    ))
  }


  # -------------------------------------------------------------------
  # 5. Restrict recheck patches to species available in the live env
  # -------------------------------------------------------------------

  # Read species currently represented by live cell -> patch vectors.
  species_available_in_env <- sort(ls(
    envir = if (compact_state) {
      patch_cell_index_by_species_env
    } else {
      patch_id_by_species_env
    },
    all.names = TRUE
  ))

  # Keep only recheck rows whose species has a live species vector.
  patches_requiring_distance_recheck <- patches_requiring_distance_recheck[
    species %in% species_available_in_env
  ]

  # Identify touched species after filtering.
  touched_species <- sort(unique(as.character(
    patches_requiring_distance_recheck$species
  )))

  # If no touched species remain, exit cleanly.
  if (!length(touched_species)) {
    # Count current patch-table rows.
    remaining_patch_count <- nrow(patch_table)

    # Count current unique species-PU combinations.
    remaining_pu_count <- data.table::uniqueN(patch_table, by = c("species", "pu_id"))

    # Count globally alive cells.
    remaining_alive_cell_count <- sum(alive_species_count_by_cell > 0L)

    if (isTRUE(emit_logs)) log_distance_connectivity_stage(
      stage_index = stage_index,
      touched_species_count = 0L,
      input_recheck_patch_count = input_recheck_patch_count,
      removed_cell_count = 0L,
      remaining_patch_count = remaining_patch_count,
      remaining_pu_count = remaining_pu_count,
      remaining_alive_cell_count = remaining_alive_cell_count
    )
    if (isTRUE(emit_logs)) {
      log_distance_timing(stage_index, distance_timing_seconds, distance_path_counts)
      log_distance_predicate_profile(stage_index, distance_predicate_profile)
    }

    # Return unchanged state.
    return(list(
      patch_table = patch_table,
      pu_graphs = pu_graphs_by_key,
      alive_species_count_by_cell = alive_species_count_by_cell,
      removed_cells_in_distance = integer(0L),
      touched_species = character(),
      timing_seconds = distance_timing_seconds,
      path_counts = distance_path_counts,
      geometry_counts = distance_geometry_counts,
      rechecked_edges = rechecked_edges,
      workload = distance_workload,
      candidate_ids_by_species = candidate_ids_by_species
    ))
  }


  # ===================================================================
  # 6. Index immutable phase inputs and process touched species
  # ===================================================================

  distance_geometry_profile <- new_distance_geometry_profile(length(touched_species))
  distance_geometry_profile_count <- 0L

  # Index the immutable patch table once instead of scanning it per species.
  patch_rows_by_species <- split(
    seq_len(nrow(patch_table)),
    as.character(patch_table$species)
  )

  # Index the compact distance-recheck input once for direct species lookup.
  recheck_patch_ids_by_species <- split(
    patches_requiring_distance_recheck$patch_id,
    as.character(patches_requiring_distance_recheck$species)
  )

  # Collect final rows for species that complete distance reconstruction.
  # The global patch table is replaced once after the species loop, avoiding
  # a full-table filter and bind for every touched species.
  replacement_species <- character(0L)
  replacement_rows_by_species <- vector("list", 0L)
  template_geometry <- prepare_distance_template_geometry(template_raster)

  # Loop over species that have at least one surviving recheck patch.
  for (species_index in seq_along(touched_species)) {
    species_name <- touched_species[[species_index]]
    if (is.function(progress_callback)) {
      progress_callback(species_index, length(touched_species), species_name)
    }
    # Read recheck patch IDs for this species.
    recheck_patch_ids_for_species <- sort(unique(as.integer(
      recheck_patch_ids_by_species[[species_name]]
    )))

    # Read current patch rows for this species.
    species_patch_row_ids <- patch_rows_by_species[[species_name]]
    if (is.null(species_patch_row_ids)) {
      species_patch_row_ids <- integer(0L)
    }
    species_patch_rows_before <- patch_table[species_patch_row_ids]

    # Skip species with no surviving patch rows.
    if (nrow(species_patch_rows_before) == 0L) {
      next
    }

    # Keep only recheck patches still present in the current patch table.
    recheck_patch_ids_for_species <- intersect(
      recheck_patch_ids_for_species,
      species_patch_rows_before$patch_id
    )

    # Skip this species if no recheck patches survive.
    if (!length(recheck_patch_ids_for_species)) {
      next
    }

    # Identify current PUs containing at least one recheck patch.
    affected_pu_ids <- sort(unique(
      species_patch_rows_before$pu_id[
        species_patch_rows_before$patch_id %in% recheck_patch_ids_for_species
      ]
    ))

    # Skip this species if no affected PUs exist.
    if (!length(affected_pu_ids)) {
      next
    }

    species_profile_started <- if (profiling_enabled) {
      proc.time()[["elapsed"]]
    } else {
      NA_real_
    }
    species_profile <- if (profiling_enabled) {
      new_distance_species_profile()
    } else {
      NULL
    }
    if (profiling_enabled) {
      species_profile$workload[["recheck_patches"]] <-
        length(recheck_patch_ids_for_species)
      species_profile$workload[["affected_pus"]] <- length(affected_pu_ids)
    }
    finish_species_profile <- if (profiling_enabled) {
      function() {
        completed <- finalize_distance_species_profile(
          species_profile,
          proc.time()[["elapsed"]] - species_profile_started
        )
        completion_callback(
          species_index,
          length(touched_species),
          species_name,
          completed
        )
        invisible(NULL)
      }
    } else {
      NULL
    }

    # Read this species' dispersal distance in kilometers.
    dispersal_threshold_km <- get_species_dispersal_threshold_km(
      species_name = species_name,
      species_params = species_params
    )

    subphase_started <- proc.time()[["elapsed"]]

    # Extract only existing graph edges incident to a recheck patch. Stage 6
    # retains isolated recheck nodes in its candidate matrix; Stage 7.2 drops
    # them because its pairwise predicate needs only actual edge endpoints.
    distance_edge_spec <- extract_distance_recheck_edges(
      species_name = species_name,
      recheck_patch_ids_for_species = recheck_patch_ids_for_species,
      affected_pu_ids = affected_pu_ids,
      pu_graphs_by_key = pu_graphs_by_key,
      return_diagnostics = TRUE
    )
    edge_diagnostic_fields <- names(distance_edge_spec$diagnostics)
    distance_path_counts[edge_diagnostic_fields] <-
      distance_path_counts[edge_diagnostic_fields] +
      distance_edge_spec$diagnostics
    rechecked_edges <- rechecked_edges + nrow(distance_edge_spec$edges)
    candidate_patch_ids <- distance_edge_spec$candidate_patch_ids
    candidate_patch_count_before_filter <- if (collect_runtime_diagnostics) {
      length(candidate_patch_ids)
    } else {
      0L
    }
    edge_endpoint_patch_ids <- if (
      identical(predicate_strategy, "edge_pairs") ||
        collect_runtime_diagnostics
    ) {
      sort(unique(c(
        distance_edge_spec$edges$patch_u,
        distance_edge_spec$edges$patch_v
      )))
    } else {
      integer()
    }
    if (collect_runtime_diagnostics) {
      distance_workload[["recheck_patches"]] <-
        distance_workload[["recheck_patches"]] +
        length(recheck_patch_ids_for_species)
      distance_workload[["candidate_patches_before_filter"]] <-
        distance_workload[["candidate_patches_before_filter"]] +
        candidate_patch_count_before_filter
      distance_workload[["edge_endpoint_patches"]] <-
        distance_workload[["edge_endpoint_patches"]] +
        length(edge_endpoint_patch_ids)
      distance_workload[["all_pair_candidate_upper_bound"]] <-
        distance_workload[["all_pair_candidate_upper_bound"]] +
        as.numeric(length(recheck_patch_ids_for_species)) *
        as.numeric(candidate_patch_count_before_filter)
      distance_workload[["graph_edge_pairs_evaluated"]] <-
        distance_workload[["graph_edge_pairs_evaluated"]] +
        nrow(distance_edge_spec$edges)
      if (!nrow(distance_edge_spec$edges)) {
        distance_workload[["zero_edge_species"]] <-
          distance_workload[["zero_edge_species"]] + 1
      }
    }
    if (identical(predicate_strategy, "edge_pairs")) {
      if (collect_runtime_diagnostics) {
        distance_workload[["isolated_rechecks_avoided"]] <-
          distance_workload[["isolated_rechecks_avoided"]] +
          length(setdiff(recheck_patch_ids_for_species, edge_endpoint_patch_ids))
      }
      candidate_patch_ids <- edge_endpoint_patch_ids
    }

    # Keep only candidate patches still present in the patch table.
    candidate_patch_ids <- intersect(
      candidate_patch_ids,
      species_patch_rows_before$patch_id
    )
    if (collect_runtime_diagnostics) {
      candidate_ids_by_species[[species_name]] <-
        as.integer(candidate_patch_ids)
    }
    if (profiling_enabled) {
      species_profile$workload[["rechecked_edges"]] <-
        nrow(distance_edge_spec$edges)
      species_profile$workload[["candidate_patches"]] <-
        length(candidate_patch_ids)
      species_profile$workload[edge_diagnostic_fields] <-
        distance_edge_spec$diagnostics
    }

    # Skip if no candidate patches remain.
    if (!length(candidate_patch_ids)) {
      candidate_preparation_seconds <- proc.time()[["elapsed"]] - subphase_started
      distance_timing_seconds[["candidate_edge_setup_seconds"]] <-
        distance_timing_seconds[["candidate_edge_setup_seconds"]] +
        candidate_preparation_seconds
      if (profiling_enabled) {
        species_profile$timing[["candidate_edge_setup_seconds"]] <-
          candidate_preparation_seconds
        finish_species_profile()
      }
      next
    }

    # Read this species' compact live patch-to-cell index.
    species_patch_index <- get(
      species_name,
      envir = patch_cell_index_by_species_env,
      inherits = FALSE
    )
    candidate_preparation_seconds <- proc.time()[["elapsed"]] - subphase_started
    if (profiling_enabled) {
      species_profile$timing[["candidate_edge_setup_seconds"]] <-
        candidate_preparation_seconds
    }

    # Build local polygons only for candidate patches.
    geometry_result <- build_local_patch_polygons(
      species_name = species_name,
      template_raster = template_raster,
      species_patch_index = species_patch_index,
      candidate_patch_ids = candidate_patch_ids,
      template_geometry = template_geometry,
      return_diagnostics = TRUE
    )
    candidate_patch_polygons <- geometry_result$polygons

    if (profiling_enabled) {
      species_profile$timing[names(geometry_result$timing_seconds)] <-
        geometry_result$timing_seconds
      species_profile$workload[["candidate_cells"]] <-
        geometry_result$counts[["candidate_cells"]]
      species_profile$workload[["local_raster_cells"]] <-
        geometry_result$counts[["local_raster_cells"]]
      species_profile$workload[["raster_expansion_ratio"]] <-
        geometry_result$counts[["expansion_ratio"]]
      species_profile$workload[["geometry_features_checked"]] <-
        geometry_result$counts[["geometry_features_checked"]]
      species_profile$workload[["invalid_geometry_features"]] <-
        geometry_result$counts[["invalid_geometry_features"]]
      species_profile$workload[["geometry_repair_invoked"]] <-
        geometry_result$counts[["geometry_repair_invoked"]]
    }

    geometry_fields <- names(geometry_result$timing_seconds)
    distance_timing_seconds[geometry_fields] <-
      distance_timing_seconds[geometry_fields] + geometry_result$timing_seconds
    distance_timing_seconds[["candidate_edge_setup_seconds"]] <-
      distance_timing_seconds[["candidate_edge_setup_seconds"]] +
      candidate_preparation_seconds
    distance_geometry_counts[["candidate_patches"]] <-
      distance_geometry_counts[["candidate_patches"]] +
      geometry_result$counts[["candidate_patches"]]
    distance_geometry_counts[["candidate_cells"]] <-
      distance_geometry_counts[["candidate_cells"]] +
      geometry_result$counts[["candidate_cells"]]
    distance_geometry_counts[["local_raster_cells"]] <-
      distance_geometry_counts[["local_raster_cells"]] +
      geometry_result$counts[["local_raster_cells"]]
    distance_geometry_counts[["geometry_features_checked"]] <-
      distance_geometry_counts[["geometry_features_checked"]] +
      geometry_result$counts[["geometry_features_checked"]]
    distance_geometry_counts[["invalid_geometry_features"]] <-
      distance_geometry_counts[["invalid_geometry_features"]] +
      geometry_result$counts[["invalid_geometry_features"]]
    distance_geometry_counts[["geometry_repair_invoked"]] <-
      distance_geometry_counts[["geometry_repair_invoked"]] +
      geometry_result$counts[["geometry_repair_invoked"]]
    distance_geometry_counts[["max_expansion_ratio"]] <- max(
      distance_geometry_counts[["max_expansion_ratio"]],
      geometry_result$counts[["expansion_ratio"]]
    )
    distance_geometry_profile_count <- distance_geometry_profile_count + 1L
    distance_geometry_profile[distance_geometry_profile_count, ] <- c(
      expansion_ratio = geometry_result$counts[["expansion_ratio"]],
      candidate_cells = geometry_result$counts[["candidate_cells"]],
      local_raster_cells = geometry_result$counts[["local_raster_cells"]],
      local_raster_seconds = geometry_result$timing_seconds[["local_raster_seconds"]],
      polygonize_seconds = geometry_result$timing_seconds[["polygonize_seconds"]],
      sf_conversion_seconds = geometry_result$timing_seconds[["sf_conversion_seconds"]],
      geometry_validation_seconds = geometry_result$timing_seconds[["geometry_validation_seconds"]],
      geometry_repair_seconds = geometry_result$timing_seconds[["geometry_repair_seconds"]]
    )

    # Skip if no candidate geometry could be reconstructed.
    if (nrow(candidate_patch_polygons) == 0L) {
      if (profiling_enabled) finish_species_profile()
      next
    }
    # Build sparse Boolean within-distance lookup for candidate polygons.
    distance_predicate_result <- build_distance_predicate_lookup(
      candidate_patch_polygons = candidate_patch_polygons,
      dispersal_threshold_km = dispersal_threshold_km,
      recheck_patch_ids = recheck_patch_ids_for_species,
      recheck_edges = distance_edge_spec$edges,
      predicate_strategy = predicate_strategy,
      return_diagnostics = TRUE
    )
    distance_predicate_lookup <- distance_predicate_result$lookup
    distance_predicate_profile <- distance_predicate_profile +
      distance_predicate_result$profile
    predicate_fields <- names(distance_predicate_result$timing_seconds)
    distance_timing_seconds[predicate_fields] <-
      distance_timing_seconds[predicate_fields] +
      distance_predicate_result$timing_seconds
    if (profiling_enabled) {
      species_profile$timing[predicate_fields] <-
        distance_predicate_result$timing_seconds
    }
    subphase_started <- proc.time()[["elapsed"]]


    # -----------------------------------------------------------------
    # 6A. Prepare affected/unaffected patch rows and graph collectors
    # -----------------------------------------------------------------

    # Keep patch rows from unaffected PUs unchanged.
    unaffected_species_patch_rows <- species_patch_rows_before[
      !(pu_id %in% affected_pu_ids),
      list(
        species = species,
        patch_id = patch_id,
        pu_id = pu_id,
        patch_area_km2 = patch_area_km2
      )
    ]

    # Allocate rebuilt patch-row collector for affected PUs.
    rebuilt_patch_rows_for_affected_pus <- vector("list", 0L)

    # Allocate rebuilt graph collector for affected PUs.
    rebuilt_graphs_for_affected_pus <- list()

    # Standard graphs with no invalid distance edge remain bound in place.
    passthrough_graph_keys <- character()

    # Track patches dropped after distance-edge filtering.
    dropped_patch_ids_from_distance <- integer(0L)

    # Initialize next available PU ID for possible PU splits.
    next_available_pu_id <- if (nrow(species_patch_rows_before)) {
      max(species_patch_rows_before$pu_id)
    } else {
      0L
    }

    # Read this species' PU-area threshold.
    pu_area_threshold <- species_params[
      species == species_name
    ]$min_population_area_km2[1L]

    # Index affected rows once so the PU loop does not repeatedly scan all
    # current rows for this species.
    affected_patch_rows_by_pu <- index_affected_patch_rows_by_pu(
      species_patch_rows_before,
      affected_pu_ids
    )
    distance_path_counts[["affected_pus"]] <-
      distance_path_counts[["affected_pus"]] + length(affected_pu_ids)
    distance_timing_seconds[["pu_index_seconds"]] <-
      distance_timing_seconds[["pu_index_seconds"]] +
      proc.time()[["elapsed"]] - subphase_started
    if (profiling_enabled) {
      species_profile$timing[["pu_index_seconds"]] <-
        proc.time()[["elapsed"]] - subphase_started
    }
    subphase_started <- proc.time()[["elapsed"]]


    # -----------------------------------------------------------------
    # 6B. Filter and rebuild each affected PU
    # -----------------------------------------------------------------

    # Loop over affected PUs.
    for (affected_pu_id in affected_pu_ids) {
      # Read current patch rows for this PU.
      current_rows_in_this_pu <- species_patch_rows_before[
        affected_patch_rows_by_pu[[as.character(affected_pu_id)]]
      ]

      # Build graph key for this species-PU.
      graph_key <- paste0(species_name, "|", affected_pu_id)

      # Read current graph object for this PU.
      old_pu_graph <- priority_graph_get(pu_graphs_by_key, graph_key)

      if (is.null(old_pu_graph)) {
        stop(
          "Missing PU graph during distance-stage update for species '",
          species_name,
          "' and pu_id ",
          affected_pu_id,
          ". This indicates inconsistent evolving state between patch_table ",
          "and pu_graphs_by_key."
        )
      }
      if (profiling_enabled) {
        input_storage <- fragmentation_graph_storage_counts(old_pu_graph)
        species_profile$workload[["input_graph_nodes"]] <-
          species_profile$workload[["input_graph_nodes"]] +
          input_storage[["nodes"]]
        species_profile$workload[["input_graph_adjacency_entries"]] <-
          species_profile$workload[["input_graph_adjacency_entries"]] +
          input_storage[["adjacency_entries"]]
      }

      # Restrict recheck patches to this PU.
      recheck_patch_ids_in_pu <- intersect(
        recheck_patch_ids_for_species,
        current_rows_in_this_pu$patch_id
      )

      # If this affected PU has no recheck patches after restriction, keep it.
      if (!length(recheck_patch_ids_in_pu)) {
        profile_started <- if (profiling_enabled) {
          proc.time()[["elapsed"]]
        } else {
          NA_real_
        }
        # Materialize lazy graphs with the compiled CSR kernel; standard graphs
        # pass through without scanning their adjacency.
        kept_result <- filter_distance_invalid_edges_for_pu(
          old_pu_graph,
          distance_predicate_lookup = list(
            patch_u = integer(), patch_v = integer(), edge_valid = logical(),
            edge_valid_by_key = stats::setNames(logical(), character())
          ),
          return_diagnostics = TRUE
        )
        graph_to_keep <- kept_result$graph
        graph_to_keep$species <- species_name
        keep_diagnostics <- kept_result$diagnostics
        distance_path_counts[["unchanged_distance_pus"]] <-
          distance_path_counts[["unchanged_distance_pus"]] + 1
        distance_path_counts[["compiled_filtered_pus"]] <-
          distance_path_counts[["compiled_filtered_pus"]] +
          as.numeric(keep_diagnostics[["adjacency_entries_scanned"]] > 0)
        distance_path_counts[["lazy_materialized_pus"]] <-
          distance_path_counts[["lazy_materialized_pus"]] +
          keep_diagnostics[["lazy_materialized"]]
        distance_path_counts[["adjacency_entries_scanned"]] <-
          distance_path_counts[["adjacency_entries_scanned"]] +
          keep_diagnostics[["adjacency_entries_scanned"]]
        distance_path_counts[["adjacency_entries_avoided"]] <-
          distance_path_counts[["adjacency_entries_avoided"]] +
          keep_diagnostics[["adjacency_entries_avoided"]]
        if (profiling_enabled) {
          species_profile$timing[["edge_filter_seconds"]] <-
            species_profile$timing[["edge_filter_seconds"]] +
            proc.time()[["elapsed"]] - profile_started
          species_profile$workload[["unchanged_distance_pus"]] <-
            species_profile$workload[["unchanged_distance_pus"]] + 1
          species_profile$workload[["compiled_filtered_pus"]] <-
            species_profile$workload[["compiled_filtered_pus"]] +
            as.numeric(keep_diagnostics[["adjacency_entries_scanned"]] > 0)
          species_profile$workload[["lazy_materialized_pus"]] <-
            species_profile$workload[["lazy_materialized_pus"]] +
            keep_diagnostics[["lazy_materialized"]]
          species_profile$workload[["adjacency_entries_scanned"]] <-
            species_profile$workload[["adjacency_entries_scanned"]] +
            keep_diagnostics[["adjacency_entries_scanned"]]
          species_profile$workload[["adjacency_entries_avoided"]] <-
            species_profile$workload[["adjacency_entries_avoided"]] +
            keep_diagnostics[["adjacency_entries_avoided"]]
          profile_started <- proc.time()[["elapsed"]]
        }

        if (isTRUE(keep_diagnostics[["unchanged"]] == 1) &&
            !is_lazy_added_node_overlay_graph(old_pu_graph)) {
          passthrough_graph_keys <- c(passthrough_graph_keys, graph_key)
        } else {
          rebuilt_graphs_for_affected_pus[[graph_key]] <- graph_to_keep
        }

        rebuilt_patch_rows_for_affected_pus[[
          length(rebuilt_patch_rows_for_affected_pus) + 1L
        ]] <- current_rows_in_this_pu

        if (profiling_enabled) {
          output_storage <- fragmentation_graph_storage_counts(graph_to_keep)
          species_profile$workload[["output_graph_nodes"]] <-
            species_profile$workload[["output_graph_nodes"]] +
            output_storage[["nodes"]]
          species_profile$workload[["output_graph_adjacency_entries"]] <-
            species_profile$workload[["output_graph_adjacency_entries"]] +
            output_storage[["adjacency_entries"]]
          species_profile$workload[["surviving_components"]] <-
            species_profile$workload[["surviving_components"]] + 1
          species_profile$timing[["pu_assembly_seconds"]] <-
            species_profile$timing[["pu_assembly_seconds"]] +
            proc.time()[["elapsed"]] - profile_started
        }

        # Continue to the next affected PU.
        next
      }

      # Filter distance-invalid edges from this PU graph.
      distance_path_counts[["edge_filtered_pus"]] <-
        distance_path_counts[["edge_filtered_pus"]] + 1L
      profile_started <- if (profiling_enabled) {
        proc.time()[["elapsed"]]
      } else {
        NA_real_
      }
      filtered_result <- filter_distance_invalid_edges_for_pu(
        pu_graph = old_pu_graph,
        distance_predicate_lookup = distance_predicate_lookup,
        return_diagnostics = TRUE
      )
      filtered_pu_graph <- filtered_result$graph
      filter_diagnostics <- filtered_result$diagnostics
      distance_path_counts[["unchanged_distance_pus"]] <-
        distance_path_counts[["unchanged_distance_pus"]] +
        as.numeric(filter_diagnostics[["invalid_undirected_edges"]] == 0)
      distance_path_counts[["compiled_filtered_pus"]] <-
        distance_path_counts[["compiled_filtered_pus"]] +
        as.numeric(filter_diagnostics[["adjacency_entries_scanned"]] > 0)
      distance_path_counts[["lazy_materialized_pus"]] <-
        distance_path_counts[["lazy_materialized_pus"]] +
        filter_diagnostics[["lazy_materialized"]]
      distance_path_counts[["invalid_edges_removed"]] <-
        distance_path_counts[["invalid_edges_removed"]] +
        filter_diagnostics[["invalid_undirected_edges"]]
      distance_path_counts[["adjacency_entries_scanned"]] <-
        distance_path_counts[["adjacency_entries_scanned"]] +
        filter_diagnostics[["adjacency_entries_scanned"]]
      distance_path_counts[["adjacency_entries_avoided"]] <-
        distance_path_counts[["adjacency_entries_avoided"]] +
        filter_diagnostics[["adjacency_entries_avoided"]]
      if (profiling_enabled) {
        species_profile$workload[["edge_filtered_pus"]] <-
          species_profile$workload[["edge_filtered_pus"]] + 1
        species_profile$timing[["edge_filter_seconds"]] <-
          species_profile$timing[["edge_filter_seconds"]] +
          proc.time()[["elapsed"]] - profile_started
        for (field in c(
          "unchanged_distance_pus", "compiled_filtered_pus",
          "lazy_materialized_pus", "invalid_edges_removed",
          "adjacency_entries_scanned", "adjacency_entries_avoided"
        )) {
          value <- switch(
            field,
            unchanged_distance_pus = as.numeric(
              filter_diagnostics[["invalid_undirected_edges"]] == 0
            ),
            compiled_filtered_pus = as.numeric(
              filter_diagnostics[["adjacency_entries_scanned"]] > 0
            ),
            lazy_materialized_pus = filter_diagnostics[["lazy_materialized"]],
            invalid_edges_removed = filter_diagnostics[["invalid_undirected_edges"]],
            adjacency_entries_scanned = filter_diagnostics[["adjacency_entries_scanned"]],
            adjacency_entries_avoided = filter_diagnostics[["adjacency_entries_avoided"]]
          )
          species_profile$workload[[field]] <-
            species_profile$workload[[field]] + value
        }
        profile_started <- proc.time()[["elapsed"]]
      }

      # With no removed edge, topology and PU thresholds are unchanged. Retain
      # the standard CSR result directly and avoid component reconstruction.
      unchanged_pu_above_threshold <- area_exceeds_threshold(
        sum(current_rows_in_this_pu$patch_area_km2), pu_area_threshold
      )
      if (filter_diagnostics[["invalid_undirected_edges"]] == 0 &&
          unchanged_pu_above_threshold) {
        if (!is_lazy_added_node_overlay_graph(old_pu_graph)) {
          graph_to_keep <- old_pu_graph
          passthrough_graph_keys <- c(passthrough_graph_keys, graph_key)
        } else {
          graph_to_keep <- list(
            species = species_name,
            pu_id = as.integer(affected_pu_id),
            id2patch = as.integer(filtered_pu_graph$id2patch),
            row_ptr = as.integer(filtered_pu_graph$row_ptr),
            col_idx = as.integer(filtered_pu_graph$col_idx)
          )
          rebuilt_graphs_for_affected_pus[[graph_key]] <- graph_to_keep
        }
        rebuilt_patch_rows_for_affected_pus[[
          length(rebuilt_patch_rows_for_affected_pus) + 1L
        ]] <- current_rows_in_this_pu
        if (profiling_enabled) {
          species_profile$workload[["output_graph_nodes"]] <-
            species_profile$workload[["output_graph_nodes"]] +
            length(graph_to_keep$id2patch)
          species_profile$workload[["output_graph_adjacency_entries"]] <-
            species_profile$workload[["output_graph_adjacency_entries"]] +
            length(graph_to_keep$col_idx)
          species_profile$workload[["surviving_components"]] <-
            species_profile$workload[["surviving_components"]] + 1
          species_profile$timing[["pu_assembly_seconds"]] <-
            species_profile$timing[["pu_assembly_seconds"]] +
            proc.time()[["elapsed"]] - profile_started
        }
        next
      }

      # Build named patch_id -> patch_area vector for this PU.
      patch_area_by_patch <- setNames(
        object = current_rows_in_this_pu$patch_area_km2,
        nm = as.character(current_rows_in_this_pu$patch_id)
      )

      if (profiling_enabled) {
        species_profile$timing[["pu_assembly_seconds"]] <-
          species_profile$timing[["pu_assembly_seconds"]] +
          proc.time()[["elapsed"]] - profile_started
        profile_started <- proc.time()[["elapsed"]]
      }

      rebuilt_pu <- rebuild_pu_after_edge_filter(
        pu_graph = filtered_pu_graph,
        patch_area_by_patch = patch_area_by_patch,
        pu_area_threshold = pu_area_threshold,
        next_available_pu_id = next_available_pu_id
      )
      if (profiling_enabled) {
        species_profile$timing[["component_rebuild_seconds"]] <-
          species_profile$timing[["component_rebuild_seconds"]] +
          proc.time()[["elapsed"]] - profile_started
        profile_started <- proc.time()[["elapsed"]]
      }

      # Carry forward the next available PU ID.
      next_available_pu_id <- rebuilt_pu$next_available_pu_id

      # Accumulate dropped patch IDs from subthreshold split components.
      if (length(rebuilt_pu$dropped_patch_ids)) {
        dropped_patch_ids_from_distance <- c(
          dropped_patch_ids_from_distance,
          rebuilt_pu$dropped_patch_ids
        )
      }

      # Store surviving PU graph objects.
      for (surviving_graph in rebuilt_pu$surviving_pu_graphs) {
        # Build graph key for this surviving component.
        surviving_graph_key <- paste0(species_name, "|", surviving_graph$pu_id)

        # Store standard CSR graph for this surviving component.
        rebuilt_graphs_for_affected_pus[[surviving_graph_key]] <- list(
          species = species_name,
          pu_id = as.integer(surviving_graph$pu_id),
          id2patch = as.integer(surviving_graph$id2patch),
          row_ptr = as.integer(surviving_graph$row_ptr),
          col_idx = as.integer(surviving_graph$col_idx)
        )
        if (profiling_enabled) {
          species_profile$workload[["output_graph_nodes"]] <-
            species_profile$workload[["output_graph_nodes"]] +
            length(surviving_graph$id2patch)
          species_profile$workload[["output_graph_adjacency_entries"]] <-
            species_profile$workload[["output_graph_adjacency_entries"]] +
            length(surviving_graph$col_idx)
        }
      }

      for (surviving_graph in rebuilt_pu$surviving_pu_graphs) {
        # Read patch IDs in this surviving component.
        component_patch_ids <- as.integer(surviving_graph$id2patch)

        # Extract component patch rows in graph patch order.
        component_rows <- current_rows_in_this_pu[
          match(component_patch_ids, patch_id)
        ]

        # Assign this component's PU ID.
        component_rows[, pu_id := as.integer(surviving_graph$pu_id)]

        # Store component rows.
        rebuilt_patch_rows_for_affected_pus[[
          length(rebuilt_patch_rows_for_affected_pus) + 1L
        ]] <- component_rows
      }
      if (profiling_enabled) {
        species_profile$workload[["surviving_components"]] <-
          species_profile$workload[["surviving_components"]] +
          length(rebuilt_pu$surviving_pu_graphs)
        species_profile$timing[["pu_assembly_seconds"]] <-
          species_profile$timing[["pu_assembly_seconds"]] +
          proc.time()[["elapsed"]] - profile_started
      }
    }


    # -----------------------------------------------------------------
    # 6C. Apply patch drops from distance-stage PU rebuilding
    # -----------------------------------------------------------------

    # Deduplicate dropped patch IDs for this species.
    dropped_patch_ids_from_distance <- sort(unique(as.integer(
      dropped_patch_ids_from_distance
    )))
    profile_started <- if (profiling_enabled) {
      proc.time()[["elapsed"]]
    } else {
      NA_real_
    }
    if (profiling_enabled) {
      species_profile$workload[["dropped_patches"]] <-
        length(dropped_patch_ids_from_distance)
    }

    # -----------------------------------------------------------------
    # 6D. Rebuild this species' patch rows
    # -----------------------------------------------------------------

    # Combine rebuilt affected-PU rows, or create an empty table.
    rebuilt_species_patch_rows <- if (length(rebuilt_patch_rows_for_affected_pus)) {
      data.table::rbindlist(
        rebuilt_patch_rows_for_affected_pus,
        use.names = TRUE,
        fill = TRUE
      )
    } else {
      data.table::data.table(
        species = character(),
        patch_id = integer(),
        pu_id = integer(),
        patch_area_km2 = numeric()
      )
    }

    # Assemble this species' complete replacement rows. An empty result is a
    # valid replacement when distance repair removes every surviving patch.
    final_species_patch_rows <- data.table::rbindlist(
      list(
        unaffected_species_patch_rows,
        rebuilt_species_patch_rows
      ),
      use.names = TRUE,
      fill = TRUE
    )

    replacement_species <- c(replacement_species, species_name)
    replacement_rows_by_species[[length(replacement_rows_by_species) + 1L]] <-
      final_species_patch_rows
    if (profiling_enabled) {
      species_profile$timing[["pu_assembly_seconds"]] <-
        species_profile$timing[["pu_assembly_seconds"]] +
        proc.time()[["elapsed"]] - profile_started
    }
    distance_timing_seconds[["pu_repair_seconds"]] <-
      distance_timing_seconds[["pu_repair_seconds"]] +
      proc.time()[["elapsed"]] - subphase_started


    # -----------------------------------------------------------------
    # 6E. Replace affected PU graphs
    # -----------------------------------------------------------------

    # Commit affected graph removals and replacements as one hashed batch.
    graph_keys_to_replace <- setdiff(
      paste0(species_name, "|", affected_pu_ids),
      passthrough_graph_keys
    )
    graph_actions <- c(
      lapply(
        graph_keys_to_replace,
        function(graph_key) list(kind = "remove", key = graph_key, graph = NULL)
      ),
      lapply(
        names(rebuilt_graphs_for_affected_pus),
        function(graph_key) list(
          kind = "set",
          key = graph_key,
          graph = rebuilt_graphs_for_affected_pus[[graph_key]]
        )
      )
    )
    if (length(rebuilt_graphs_for_affected_pus) && any(vapply(
      rebuilt_graphs_for_affected_pus,
      is_lazy_added_node_overlay_graph,
      logical(1L)
    ))) {
      stop(
        "Distance repair produced a lazy overlay graph for species: ",
        species_name,
        call. = FALSE
      )
    }
    subphase_started <- proc.time()[["elapsed"]]
    pu_graphs_by_key <- priority_graph_apply_actions(
      pu_graphs_by_key,
      graph_actions
    )
    graph_commit_elapsed <- proc.time()[["elapsed"]] - subphase_started
    distance_timing_seconds[["graph_commit_seconds"]] <-
      distance_timing_seconds[["graph_commit_seconds"]] +
      graph_commit_elapsed
    if (profiling_enabled) {
      species_profile$timing[["graph_commit_seconds"]] <-
        graph_commit_elapsed
    }
    subphase_started <- proc.time()[["elapsed"]]


    # -----------------------------------------------------------------
    # 6F. Update alive-species counts for cells where this species vanished
    # -----------------------------------------------------------------

    # Leave live cell/index bindings untouched unless distance repair actually
    # drops patches. When it does, remove only their compact indexed cells.
    species_cells_removed_in_distance <- integer()
    if (length(dropped_patch_ids_from_distance)) {
      index_update <- update_species_patch_index_incremental(
        patch_index = species_patch_index,
        replacement_rows_by_origin = list(),
        expected_patch_ids = sort(unique(as.integer(final_species_patch_rows$patch_id))),
        pre_update_patch_ids = species_patch_rows_before$patch_id
      )
      compact_index_elapsed <- proc.time()[["elapsed"]] - subphase_started
      distance_timing_seconds[["compact_index_seconds"]] <-
        distance_timing_seconds[["compact_index_seconds"]] +
        compact_index_elapsed
      if (profiling_enabled) {
        species_profile$timing[["compact_index_seconds"]] <-
          compact_index_elapsed
      }
      subphase_started <- proc.time()[["elapsed"]]
      species_patch_ids <- if (compact_state) {
        NULL
      } else {
        get(
          species_name,
          envir = patch_id_by_species_env,
          inherits = FALSE
        )
      }
      species_cells_removed_in_distance <- if (compact_state) {
        as.integer(index_update$removed_cells)
      } else {
        index_update$removed_cells[
          !is.na(species_patch_ids[index_update$removed_cells])
        ]
      }
      if (!length(species_cells_removed_in_distance)) {
        stop(
          "Distance-dropped patches have no live indexed cells for species: ",
          species_name,
          call. = FALSE
        )
      }
      if (!compact_state) {
        species_patch_ids[species_cells_removed_in_distance] <- NA_integer_
        assign(species_name, species_patch_ids, envir = patch_id_by_species_env)
      }
      assign(
        species_name,
        index_update$index,
        envir = patch_cell_index_by_species_env
      )
    }

    # Read global alive-species counts before removing this species.
    count_before_species_loss <- alive_species_count_by_cell[
      species_cells_removed_in_distance
    ]

    # Identify cells that remain globally alive after losing this species.
    cells_still_alive_after_species_loss <- species_cells_removed_in_distance[
      count_before_species_loss > 0L
    ]

    # Identify cells that become globally empty after losing this species.
    cells_becoming_empty <- species_cells_removed_in_distance[
      count_before_species_loss == 1L
    ]

    # Decrement alive-species counts for cells that were still globally alive.
    if (length(cells_still_alive_after_species_loss)) {
      alive_species_count_by_cell[cells_still_alive_after_species_loss] <-
        alive_species_count_by_cell[cells_still_alive_after_species_loss] - 1L
    }

    # Add newly empty cells to the distance-stage accumulator.
    if (length(cells_becoming_empty)) {
      removed_cells_in_distance <- c(
        removed_cells_in_distance,
        cells_becoming_empty
      )
    }

    distance_timing_seconds[["cell_state_seconds"]] <-
      distance_timing_seconds[["cell_state_seconds"]] +
      proc.time()[["elapsed"]] - subphase_started

    if (profiling_enabled) {
      species_profile$timing[["cell_state_seconds"]] <-
        proc.time()[["elapsed"]] - subphase_started
      finish_species_profile()
    }


  }

  rm(patch_rows_by_species, recheck_patch_ids_by_species)


  # -------------------------------------------------------------------
  # 7. Commit all species patch-table replacements once
  # -------------------------------------------------------------------

  if (length(replacement_species)) {
    subphase_started <- proc.time()[["elapsed"]]
    unchanged_patch_rows <- patch_table[
      !(species %chin% replacement_species)
    ]

    patch_table <- data.table::rbindlist(
      c(list(unchanged_patch_rows), replacement_rows_by_species),
      use.names = TRUE,
      fill = TRUE
    )
    distance_timing_seconds[["state_commit_seconds"]] <-
      distance_timing_seconds[["state_commit_seconds"]] +
      proc.time()[["elapsed"]] - subphase_started
  }


  # -------------------------------------------------------------------
  # 8. Finalize distance-stage removed cells
  # -------------------------------------------------------------------

  # Deduplicate and sort globally empty cells created by the distance stage.
  removed_cells_in_distance <- sort(
    unique(as.integer(removed_cells_in_distance))
  )


  # -------------------------------------------------------------------
  # 9. Compute final stage summaries and log one line
  # -------------------------------------------------------------------

  # Count remaining patch-table rows.
  remaining_patch_count <- nrow(patch_table)

  # Count remaining unique species-PU combinations.
  remaining_pu_count <- data.table::uniqueN(
    patch_table,
    by = c("species", "pu_id")
  )

  # Count remaining globally alive cells.
  remaining_alive_cell_count <- sum(alive_species_count_by_cell > 0L)

  classified_fields <- setdiff(
    names(distance_timing_seconds),
    "other_seconds"
  )
  distance_timing_seconds[["other_seconds"]] <- max(
    0,
    proc.time()[["elapsed"]] - distance_stage_started -
      sum(distance_timing_seconds[classified_fields])
  )
  if (isTRUE(emit_logs)) log_distance_connectivity_stage(
    stage_index = stage_index,
    touched_species_count = length(touched_species),
    input_recheck_patch_count = input_recheck_patch_count,
    removed_cell_count = length(removed_cells_in_distance),
    remaining_patch_count = remaining_patch_count,
    remaining_pu_count = remaining_pu_count,
    remaining_alive_cell_count = remaining_alive_cell_count
  )
  if (isTRUE(emit_logs)) log_distance_timing(
    stage_index,
    distance_timing_seconds,
    distance_path_counts,
    distance_geometry_counts
  )
  if (isTRUE(emit_logs)) {
    log_distance_predicate_profile(stage_index, distance_predicate_profile)
  }
  if (isTRUE(emit_logs) && distance_geometry_profile_count > 0L) {
    log_distance_geometry_profile(
      stage_index,
      distance_geometry_profile[seq_len(distance_geometry_profile_count), , drop = FALSE]
    )
  }
  if (collect_runtime_diagnostics) {
    distance_workload[["valid_edges"]] <-
      distance_predicate_profile[["valid_graph_edges"]]
    distance_workload[["invalid_edges"]] <-
      distance_predicate_profile[["invalid_graph_edges"]]
    distance_workload[["predicate_spatial_seconds"]] <-
      distance_predicate_profile[["predicate_spatial_seconds"]]
    distance_workload[["predicate_mapping_seconds"]] <-
      distance_predicate_profile[["predicate_mapping_seconds"]]
  }


  # -------------------------------------------------------------------
  # 10. Return updated evolving state
  # -------------------------------------------------------------------

  # Return the updated state needed by the outer pipeline.
  list(
    patch_table = patch_table,
    pu_graphs = pu_graphs_by_key,
    alive_species_count_by_cell = alive_species_count_by_cell,
    removed_cells_in_distance = removed_cells_in_distance,
    touched_species = touched_species,
    timing_seconds = distance_timing_seconds,
    path_counts = distance_path_counts,
    geometry_counts = distance_geometry_counts,
    rechecked_edges = rechecked_edges,
    workload = distance_workload,
    candidate_ids_by_species = candidate_ids_by_species
  )
}
