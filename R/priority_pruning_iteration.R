# One sparse-frontier pruning iteration.
#
# The function selects a deterministic cell batch, updates affected patch/PU
# state, invokes compiled batched component repair, commits validated graph
# actions, and records compact timing/work counters. Large dense vectors and
# graph storage are mutated only at established synchronization points.

# Helper functions for one frontier-cell pruning iteration: scoring,
# partial selection, patch loss, CSR graph repair, and threshold bookkeeping.

log_pruning_iteration <- function(
  prune_iteration,
  frontier_cell_count,
  selected_cell_count = NULL,
  newly_empty_cell_count = NULL,
  touched_patch_count = NULL,
  dropped_patch_threshold_count = NULL,
  dropped_component_threshold_count = NULL,
  dropped_pu_threshold_count = NULL,
  dropped_pu_count = NULL,
  new_pu_count = NULL,
  frontier_exhausted = FALSE,
  stage_index = NA_integer_
) {
  # If the frontier is exhausted, print the short stopping message.
  if (isTRUE(frontier_exhausted)) {
    # Print exactly one frontier-exhausted log line.
    runtime_log_event(
      "prune_step",
      console = FALSE,
      stage = as.integer(stage_index),
      iter = as.integer(prune_iteration),
      status = "frontier_exhausted",
      frontier = as.integer(frontier_cell_count)
    )

    return(invisible(NULL))
  }

  # Print exactly one normal pruning-iteration summary line.
  runtime_log_event(
    "prune_step",
    console = FALSE,
    stage = as.integer(stage_index),
    iter = as.integer(prune_iteration),
    frontier = as.integer(frontier_cell_count),
    selected = as.integer(selected_cell_count),
    new_empty = as.integer(newly_empty_cell_count),
    touched_patches = as.integer(touched_patch_count),
    repair_patch = as.integer(dropped_patch_threshold_count),
    repair_component = as.integer(dropped_component_threshold_count),
    repair_pu_patch = as.integer(dropped_pu_threshold_count),
    dropped_pus = as.integer(dropped_pu_count),
    new_pus = as.integer(new_pu_count)
  )

  invisible(NULL)
}


# Compact pruning timing and batched patch-table repair.

empty_pruning_timing <- function() {
  stats::setNames(
    numeric(12L),
    c(
      "frontier_fetch_seconds",
      "score_cache_seconds",
      "frontier_scoring_seconds",
      "batch_selection_seconds",
      "ecology_diagnostics_seconds",
      "patch_update_seconds",
      "pu_repair_seconds",
      "drop_cell_lookup_seconds",
      "alive_count_update_seconds",
      "frontier_update_seconds",
      "iteration_bookkeeping_seconds",
      "deferred_cell_clear_seconds"
    )
  )
}

empty_pruning_work_counts <- function() {
  stats::setNames(
    integer(3L),
    c(
      "threshold_species_updates_deferred",
      "threshold_cell_entries_deferred",
      "cascade_empty_cells"
    )
  )
}

empty_pruning_pu_repair_timing <- function() {
  stats::setNames(numeric(4L), c(
    "setup_seconds", "component_rebuild_seconds",
    "transaction_seconds", "pu_threshold_seconds"
  ))
}

empty_pruning_pu_repair_counts <- function() {
  stats::setNames(numeric(11L), c(
    "species_pu_repairs", "graph_nodes_examined",
    "adjacency_entries_examined", "dead_patch_nodes",
    "surviving_components", "split_pus",
    "component_dropped_patches", "threshold_candidate_pus",
    "threshold_dropped_pus", "compiled_species_batches",
    "compiled_pu_workloads"
  ))
}

format_pruning_assignment_duplicates <- function(component_assignments) {
  duplicated_keys <- component_assignments[
    ,
    .N,
    by = .(species, patch_id)
  ][N > 1L, .(species, patch_id)]

  details <- component_assignments[
    duplicated_keys,
    on = .(species, patch_id),
    nomatch = 0L
  ][
    ,
    .(pu_ids = paste(unique(as.integer(pu_id)), collapse = ",")),
    by = .(species, patch_id)
  ]

  paste(
    sprintf(
      "species=%s, patch_id=%d, assigned_pu_ids=%s",
      details$species,
      details$patch_id,
      details$pu_ids
    )[seq_len(min(5L, nrow(details)))],
    collapse = "; "
  )
}

pruning_pre_step_pu_high_water <- function(cache_entry, species_name) {
  pre_step_pu_lookup <- cache_entry$pu_id_by_patch_id
  represented_pu_ids <- as.numeric(
    pre_step_pu_lookup[!is.na(pre_step_pu_lookup)]
  )

  if (!length(represented_pu_ids) ||
      any(!is.finite(represented_pu_ids)) ||
      any(represented_pu_ids < 1) ||
      any(represented_pu_ids != floor(represented_pu_ids))) {
    stop("Invalid pre-removal PU-ID lookup for species: ", species_name)
  }

  as.integer(max(represented_pu_ids))
}

find_pruning_patch_threshold_drops <- function(
  patch_table,
  touched_patches,
  species_params
) {
  touched_patches <- unique(
    touched_patches[, .(
      species = as.character(species),
      patch_id = as.integer(patch_id)
    )],
    by = c("species", "patch_id")
  )
  if (!nrow(touched_patches)) {
    return(data.table::data.table(species = character(), patch_id = integer()))
  }

  touched_patch_rows <- patch_table[
    touched_patches,
    on = .(species, patch_id),
    nomatch = 0L
  ]
  if (nrow(touched_patch_rows) != nrow(touched_patches)) {
    stop("A directly touched patch is missing from the current patch table.")
  }

  # Put the compact touched-patch table on the i side so species that were not
  # touched cannot create unmatched NA patch rows.
  threshold_rows <- species_params[
    touched_patch_rows,
    on = .(species),
    nomatch = 0L
  ]
  if (nrow(threshold_rows) != nrow(touched_patch_rows)) {
    stop("A directly touched patch has no species threshold parameters.")
  }

  threshold_rows[
    !area_exceeds_threshold(patch_area_km2, min_patch_area_km2),
    .(species = as.character(species), patch_id = as.integer(patch_id))
  ]
}

snapshot_pruning_affected_graphs <- function(
  affected_patch_to_pu_before_step,
  pu_graphs_by_key
) {
  if (!nrow(affected_patch_to_pu_before_step)) {
    return(list(
      affected_pus = data.table::data.table(
        species = character(),
        pu_id = integer()
      ),
      affected_keys = character(),
      graph_snapshot = list()
    ))
  }

  affected_pus <- unique(
    affected_patch_to_pu_before_step[, .(species, pu_id)],
    by = c("species", "pu_id")
  )
  affected_keys <- paste0(affected_pus$species, "|", affected_pus$pu_id)
  graph_present <- priority_graph_has(pu_graphs_by_key, affected_keys)
  if (any(!graph_present)) {
    stop(
      "Missing pre-iteration graph(s) for affected PU key(s): ",
      paste(utils::head(affected_keys[!graph_present], 5L), collapse = ", ")
    )
  }

  graph_snapshot <- priority_graph_snapshot(
    pu_graphs_by_key,
    affected_keys,
    require_all = TRUE
  )
  mapping_keys <- paste0(
    affected_patch_to_pu_before_step$species,
    "|",
    affected_patch_to_pu_before_step$pu_id
  )
  expected_patch_ids_by_key <- split(
    affected_patch_to_pu_before_step$patch_id,
    mapping_keys
  )

  for (i in seq_along(affected_keys)) {
    graph_key <- affected_keys[[i]]
    species_name <- affected_pus$species[[i]]
    pu_id_value <- as.integer(affected_pus$pu_id[[i]])
    graph <- graph_snapshot[[i]]
    graph_patch_ids <- as.integer(graph$id2patch)
    expected_patch_ids <- as.integer(expected_patch_ids_by_key[[graph_key]])

    if (!is.list(graph) ||
        !identical(as.character(graph$species), species_name) ||
        !identical(as.integer(graph$pu_id), pu_id_value)) {
      stop("Affected pre-iteration graph metadata differs from its key: ", graph_key)
    }
    if (!length(graph_patch_ids) || anyNA(graph_patch_ids) ||
        any(graph_patch_ids < 1L) || anyDuplicated(graph_patch_ids)) {
      stop("Affected pre-iteration graph contains invalid patch IDs: ", graph_key)
    }
    if (length(graph_patch_ids) != length(expected_patch_ids) ||
        !setequal(graph_patch_ids, expected_patch_ids)) {
      stop("Affected pre-iteration graph and patch mapping differ for: ", graph_key)
    }
  }

  list(
    affected_pus = affected_pus,
    affected_keys = affected_keys,
    graph_snapshot = graph_snapshot
  )
}

validate_pruning_graph_transaction <- function(
  affected_patch_to_pu_before_step,
  direct_patch_drops,
  component_drops,
  component_assignments,
  graph_actions,
  pu_graphs_before,
  pre_step_max_pu_by_species
) {
  direct_patch_drops <- unique(
    data.table::as.data.table(direct_patch_drops)[, .(species, patch_id)],
    by = c("species", "patch_id")
  )
  component_drops <- data.table::as.data.table(component_drops)
  component_assignments <- data.table::as.data.table(component_assignments)

  if (anyDuplicated(component_assignments, by = c("species", "patch_id"))) {
    stop(
      "Component assignments contain duplicate species/patch IDs: ",
      format_pruning_assignment_duplicates(component_assignments)
    )
  }

  original_patches <- unique(
    affected_patch_to_pu_before_step[, .(species, patch_id)],
    by = c("species", "patch_id")
  )
  if (nrow(direct_patch_drops[
    !original_patches,
    on = .(species, patch_id)
  ])) {
    stop("Direct patch-threshold drops fall outside the affected original PUs.")
  }

  surviving_after_direct <- original_patches[
    !direct_patch_drops,
    on = .(species, patch_id)
  ]
  terminal_rows <- data.table::rbindlist(
    list(
      component_assignments[, .(species, patch_id)],
      component_drops[, .(species, patch_id)]
    ),
    use.names = TRUE
  )
  if (anyDuplicated(terminal_rows, by = c("species", "patch_id"))) {
    stop("An affected patch is both retained and dropped during component repair.")
  }
  if (nrow(surviving_after_direct[
    !terminal_rows,
    on = .(species, patch_id)
  ]) || nrow(terminal_rows[
    !surviving_after_direct,
    on = .(species, patch_id)
  ])) {
    stop("Component repair does not exactly partition the surviving affected patches.")
  }

  action_keys <- vapply(graph_actions, `[[`, character(1L), "key")
  action_kinds <- vapply(graph_actions, `[[`, character(1L), "kind")
  allowed_action_kinds <- c("remove_original", "replace_original", "add_new")
  if (anyNA(action_keys) || any(!nzchar(action_keys)) ||
      anyNA(action_kinds) || any(!action_kinds %in% allowed_action_kinds)) {
    stop("A staged component graph action has an invalid key or action type.")
  }

  original_action_keys <- action_keys[action_kinds != "add_new"]
  expected_original_pus <- unique(
    affected_patch_to_pu_before_step[, .(species, pu_id)],
    by = c("species", "pu_id")
  )
  expected_original_keys <- if (nrow(expected_original_pus)) {
    paste0(
      expected_original_pus$species,
      "|",
      expected_original_pus$pu_id
    )
  } else {
    character()
  }
  if (anyDuplicated(original_action_keys) ||
      !setequal(original_action_keys, expected_original_keys)) {
    stop("Each affected original PU must receive exactly one graph action.")
  }

  new_action_keys <- action_keys[action_kinds == "add_new"]
  if (anyDuplicated(new_action_keys)) {
    stop("New component graph keys are duplicated within one pruning iteration.")
  }
  colliding_new_keys <- new_action_keys[
    priority_graph_has(pu_graphs_before, new_action_keys)
  ]
  if (length(colliding_new_keys)) {
    stop(
      "New component graph key collides with a pre-iteration graph: ",
      paste(utils::head(colliding_new_keys, 5L), collapse = ", ")
    )
  }

  graph_membership_rows <- vector("list", 0L)
  new_pu_ids_by_species <- list()
  for (action in graph_actions) {
    if (identical(action$kind, "remove_original")) {
      next
    }

    graph <- action$graph
    graph_patch_ids <- as.integer(graph$id2patch)
    expected_key <- paste0(graph$species, "|", as.integer(graph$pu_id))
    if (!identical(action$key, expected_key) ||
        !length(graph_patch_ids) || anyNA(graph_patch_ids) ||
        any(graph_patch_ids < 1L) || anyDuplicated(graph_patch_ids)) {
      stop("A staged component graph has invalid metadata or patch IDs: ", action$key)
    }

    if (identical(action$kind, "add_new")) {
      species_name <- as.character(graph$species)
      pu_id_value <- as.integer(graph$pu_id)
      species_high_water <- as.integer(pre_step_max_pu_by_species[[species_name]])
      if (length(species_high_water) != 1L || is.na(species_high_water) ||
          pu_id_value <= species_high_water) {
        stop("A new PU ID does not exceed its species pre-removal maximum: ", action$key)
      }
      new_pu_ids_by_species[[species_name]] <- c(
        new_pu_ids_by_species[[species_name]],
        pu_id_value
      )
    }

    graph_membership_rows[[length(graph_membership_rows) + 1L]] <-
      data.table::data.table(
        species = as.character(graph$species),
        patch_id = graph_patch_ids,
        pu_id = as.integer(graph$pu_id)
      )
  }

  for (species_name in names(new_pu_ids_by_species)) {
    new_ids <- as.integer(new_pu_ids_by_species[[species_name]])
    if (anyDuplicated(new_ids) || (length(new_ids) > 1L && any(diff(new_ids) <= 0L))) {
      stop("New PU IDs are not strictly increasing for species: ", species_name)
    }
  }

  graph_membership <- if (length(graph_membership_rows)) {
    data.table::rbindlist(graph_membership_rows, use.names = TRUE)
  } else {
    data.table::data.table(species = character(), patch_id = integer(), pu_id = integer())
  }
  if (nrow(graph_membership) != nrow(component_assignments) ||
      nrow(graph_membership[
        !component_assignments,
        on = .(species, patch_id, pu_id)
      ]) || nrow(component_assignments[
        !graph_membership,
        on = .(species, patch_id, pu_id)
      ])) {
    stop("Staged graph membership differs from component patch assignments.")
  }

  invisible(TRUE)
}

validate_pruning_component_commit <- function(
  patch_table,
  pu_graphs_by_key,
  component_assignments,
  component_drops,
  graph_actions
) {
  if (nrow(component_drops) && nrow(patch_table[
    component_drops,
    on = .(species, patch_id),
    nomatch = 0L
  ])) {
    stop("A component-threshold patch remains in the committed patch table.")
  }

  if (nrow(component_assignments)) {
    committed_assignments <- patch_table[
      component_assignments,
      on = .(species, patch_id),
      nomatch = 0L,
      .(
        species,
        patch_id,
        expected_pu_id = as.integer(i.pu_id),
        committed_pu_id = as.integer(pu_id)
      )
    ]
    if (nrow(committed_assignments) != nrow(component_assignments) ||
        any(committed_assignments$expected_pu_id != committed_assignments$committed_pu_id)) {
      stop("Committed patch-table PU assignments differ from the staged transaction.")
    }
  }

  for (action in graph_actions) {
    committed_graph <- priority_graph_get(pu_graphs_by_key, action$key)
    if (identical(action$kind, "remove_original")) {
      if (!is.null(committed_graph)) {
        stop("A removed original PU graph remains after commit: ", action$key)
      }
    } else if (is.null(committed_graph) ||
               !identical(committed_graph, action$graph)) {
      stop("A committed PU graph differs from its staged replacement: ", action$key)
    }
  }

  invisible(TRUE)
}

apply_pruning_component_changes <- function(
  patch_table,
  component_assignments,
  component_drops
) {
  component_assignments <- data.table::as.data.table(component_assignments)
  component_drops <- data.table::as.data.table(component_drops)

  if (nrow(component_assignments)) {
    required_assignment_columns <- c("species", "patch_id", "pu_id")
    if (!identical(names(component_assignments), required_assignment_columns)) {
      stop("Component assignments have an invalid schema.")
    }
    if (anyDuplicated(component_assignments, by = c("species", "patch_id"))) {
      stop(
        "Component assignments contain duplicate species/patch IDs: ",
        format_pruning_assignment_duplicates(component_assignments)
      )
    }
  }

  if (nrow(component_drops)) {
    required_drop_columns <- c("species", "patch_id")
    if (!identical(names(component_drops), required_drop_columns)) {
      stop("Component drops have an invalid schema.")
    }
    if (anyDuplicated(component_drops, by = c("species", "patch_id"))) {
      stop("Component drops contain duplicate species/patch IDs.")
    }
  }

  if (nrow(component_assignments) && nrow(component_drops) &&
      nrow(component_assignments[
        component_drops,
        on = .(species, patch_id),
        nomatch = 0L
      ])) {
    stop("A component-repair patch cannot be both retained and dropped.")
  }

  if (nrow(component_drops)) {
    patch_table <- patch_table[
      !component_drops,
      on = .(species, patch_id)
    ]
  }

  if (nrow(component_assignments)) {
    matched_assignments <- patch_table[
      component_assignments,
      on = .(species, patch_id),
      nomatch = 0L,
      .N
    ]
    if (!identical(as.integer(matched_assignments), nrow(component_assignments))) {
      stop("Component assignments reference missing patch-table rows.")
    }

    patch_table[
      component_assignments,
      on = .(species, patch_id),
      pu_id := i.pu_id
    ]
  }

  patch_table
}

collect_pruning_pu_drop_rows <- function(patch_table, dropped_pus) {
  dropped_pus <- data.table::as.data.table(dropped_pus)
  if (!nrow(dropped_pus)) {
    return(data.table::data.table(
      species = character(),
      patch_id = integer(),
      pu_id = integer()
    ))
  }
  if (!identical(names(dropped_pus), c("species", "pu_id")) ||
      anyDuplicated(dropped_pus, by = c("species", "pu_id"))) {
    stop("Dropped PUs have an invalid or duplicated species/PU contract.")
  }

  dropped_rows <- patch_table[
    dropped_pus,
    on = .(species, pu_id),
    nomatch = 0L,
    .(species, patch_id, pu_id)
  ]
  represented_pus <- unique(
    dropped_rows[, .(species, pu_id)],
    by = c("species", "pu_id")
  )
  if (nrow(represented_pus) != nrow(dropped_pus) ||
      nrow(dropped_pus[
        !represented_pus,
        on = .(species, pu_id)
      ])) {
    stop("A PU selected for threshold removal has no patch-table rows.")
  }

  dropped_rows
}


# Reconstruct pre-step mappings only for directly affected PUs.

build_affected_patch_to_pu_mapping <- function(
  touched_patches,
  patch_score_cache,
  pu_graphs_by_key
) {
  empty_mapping <- data.table::data.table(
    species = character(),
    patch_id = integer(),
    pu_id = integer()
  )
  if (!nrow(touched_patches)) {
    return(empty_mapping)
  }

  touched_patches <- unique(
    touched_patches[, .(
      species = as.character(species),
      patch_id = as.integer(patch_id)
    )],
    by = c("species", "patch_id")
  )

  touched_by_species <- split(touched_patches$patch_id, touched_patches$species)
  touched_mapping_rows <- vector("list", length(touched_by_species))
  names(touched_mapping_rows) <- names(touched_by_species)

  for (species_name in names(touched_by_species)) {
    patch_ids <- touched_by_species[[species_name]]
    cache_entry <- patch_score_cache[[species_name]]
    pu_lookup <- cache_entry$pu_id_by_patch_id

    if (is.null(cache_entry) || is.null(pu_lookup)) {
      stop("Missing pre-step patch-to-PU cache for species: ", species_name)
    }
    if (anyNA(patch_ids) || any(patch_ids < 1L) ||
        any(patch_ids > length(pu_lookup))) {
      stop("Touched patch IDs cannot index the pre-step PU cache for species: ", species_name)
    }

    pu_ids <- pu_lookup[patch_ids]
    if (anyNA(pu_ids) || any(pu_ids < 1L)) {
      stop("Touched patches have no valid pre-step PU for species: ", species_name)
    }

    touched_mapping_rows[[species_name]] <- data.table::data.table(
      species = species_name,
      patch_id = patch_ids,
      pu_id = as.integer(pu_ids)
    )
  }

  touched_mapping <- data.table::rbindlist(
    touched_mapping_rows,
    use.names = TRUE
  )
  affected_pus <- unique(
    touched_mapping[, .(species, pu_id)],
    by = c("species", "pu_id")
  )

  # Inspect each affected species lookup once, then validate and materialize
  # only the original PU graphs that can participate in cascading repair.
  affected_pus_by_species <- split(affected_pus$pu_id, affected_pus$species)
  mapping_rows <- vector("list", nrow(affected_pus))
  output_index <- 0L

  for (species_name in names(affected_pus_by_species)) {
    pu_ids <- affected_pus_by_species[[species_name]]
    pu_lookup <- patch_score_cache[[species_name]]$pu_id_by_patch_id
    represented_patch_ids <- which(!is.na(pu_lookup) & pu_lookup %in% pu_ids)
    represented_by_pu <- split(
      represented_patch_ids,
      pu_lookup[represented_patch_ids]
    )

    for (pu_id in pu_ids) {
      graph_key <- paste0(species_name, "|", pu_id)
      pu_graph <- priority_graph_get(pu_graphs_by_key, graph_key)
      if (is.null(pu_graph)) {
        stop("Missing pre-step PU graph for affected key: ", graph_key)
      }

      graph_patch_ids <- as.integer(pu_graph$id2patch)
      expected_patch_ids <- as.integer(represented_by_pu[[as.character(pu_id)]])
      if (!length(graph_patch_ids) || anyNA(graph_patch_ids) ||
          any(graph_patch_ids < 1L) || anyDuplicated(graph_patch_ids)) {
        stop("Invalid pre-step patch IDs in affected PU graph: ", graph_key)
      }
      if (length(expected_patch_ids) != length(graph_patch_ids) ||
          !setequal(expected_patch_ids, graph_patch_ids)) {
        stop("Pre-step cache and PU graph patch sets differ for: ", graph_key)
      }

      output_index <- output_index + 1L
      mapping_rows[[output_index]] <- data.table::data.table(
        species = species_name,
        patch_id = graph_patch_ids,
        pu_id = as.integer(pu_id)
      )
    }
  }

  mapping <- data.table::rbindlist(
    mapping_rows[seq_len(output_index)],
    use.names = TRUE
  )
  if (anyDuplicated(mapping, by = c("species", "patch_id"))) {
    stop("Affected pre-step patch mapping contains duplicate species/patch IDs.")
  }

  mapping
}

run_pruning_iteration <- function(
  prune_iteration,
  cells_to_remove_per_step,
  patch_table,
  pu_graphs_by_key,
  alive_species_count_by_cell,
  cell_area_by_cell,
  species_params,
  patch_id_by_species_env,
  patch_cell_index_by_species_env,
  cell_species_index = NULL,
  frontier_state,
  rook_neighbor_index,
  species_names,
  patch_score_cache,
  dirty_species_for_scoring,
  stage_index = NA_integer_,
  ecology_log_every_iterations = NULL,
  ecology_previous_snapshot = NULL,
  initial_alive_cell_count = NULL,
  current_alive_cell_count = NULL
) {
  pruning_timing_seconds <- empty_pruning_timing()
  pruning_work_counts <- empty_pruning_work_counts()
  frontier_scoring_counts <- empty_frontier_scoring_counts()
  frontier_scoring_timing <- empty_frontier_scoring_timing()
  pu_repair_profile_timing <- empty_pruning_pu_repair_timing()
  pu_repair_profile_counts <- empty_pruning_pu_repair_counts()
  pu_repair_graph_actions <- 0
  phase_started <- proc.time()[["elapsed"]]

  # -------------------------------------------------------------------
  # STEP 1. Identify the current frontier
  # -------------------------------------------------------------------

  # Read all currently removable frontier cells.
  frontier_cells <- get_frontier_cells(frontier_state)
  pruning_timing_seconds[["frontier_fetch_seconds"]] <-
    proc.time()[["elapsed"]] - phase_started

  # Stop this pruning stage if the frontier is too small for a full batch.
  if (length(frontier_cells) < cells_to_remove_per_step) {
    # Print the single frontier-exhausted log line.
    log_pruning_iteration(
      stage_index = stage_index,
      prune_iteration = prune_iteration,
      frontier_cell_count = length(frontier_cells),
      frontier_exhausted = TRUE
    )

    # Return unchanged state and empty change bookkeeping.
    return(list(
      patch_table = patch_table,
      pu_graphs = pu_graphs_by_key,
      alive_species_count_by_cell = alive_species_count_by_cell,
      frontier_state = frontier_state,
      newly_empty_cells = integer(0L),
      changed_patches = data.table::data.table(
        species = character(),
        patch_id = integer()
      ),
      threshold_dropped_patches = data.table::data.table(
        species = character(),
        patch_id = integer()
      ),
      patch_score_cache = patch_score_cache,
      dirty_species_for_scoring = dirty_species_for_scoring,
      ecology_previous_snapshot = ecology_previous_snapshot,
      current_alive_cell_count = current_alive_cell_count,
      pruning_timing_seconds = pruning_timing_seconds,
      pruning_work_counts = pruning_work_counts,
      frontier_scoring_counts = frontier_scoring_counts,
      frontier_scoring_timing = frontier_scoring_timing,
      pu_repair_profile_timing = pu_repair_profile_timing,
      pu_repair_profile_counts = pu_repair_profile_counts,
      pu_repair_graph_actions = pu_repair_graph_actions
    ))
  }

  # -------------------------------------------------------------------
  # STEP 2. Score the frontier and choose cells to remove
  # -------------------------------------------------------------------

  # Update cached patch scores only for species whose state changed.
  phase_started <- proc.time()[["elapsed"]]
  patch_score_cache <- update_patch_score_cache(
    patch_table = patch_table,
    species_params = species_params,
    patch_score_cache = patch_score_cache,
    dirty_species = dirty_species_for_scoring,
    include_ecology = !is.null(ecology_log_every_iterations)
  )
  pruning_timing_seconds[["score_cache_seconds"]] <-
    proc.time()[["elapsed"]] - phase_started

  # Score every current frontier cell.
  phase_started <- proc.time()[["elapsed"]]
  frontier_scoring_result <- score_frontier_cells(
    frontier_cells = frontier_cells,
    patch_scores_by_species = patch_score_cache,
    patch_id_by_species_env = patch_id_by_species_env,
    cell_species_index = cell_species_index,
    return_diagnostics = TRUE
  )
  frontier_scores <- frontier_scoring_result$scores
  frontier_scoring_counts <- frontier_scoring_result$counts
  frontier_scoring_timing <- frontier_scoring_result$timing_seconds
  pruning_timing_seconds[["frontier_scoring_seconds"]] <-
    proc.time()[["elapsed"]] - phase_started

  phase_started <- proc.time()[["elapsed"]]
  selected_cells <- choose_frontier_cells_to_remove(
    frontier_cells = frontier_cells,
    frontier_scores = frontier_scores,
    cells_to_remove_per_step = cells_to_remove_per_step
  )
  if (length(selected_cells) != cells_to_remove_per_step || anyDuplicated(selected_cells)) {
    stop("A normal pruning iteration must select exactly one unique full batch of frontier cells.")
  }
  if (!all(selected_cells %in% frontier_cells) ||
      any(alive_species_count_by_cell[selected_cells] <= 0L)) {
    stop("Selected pruning cells must be currently alive members of the frontier.")
  }
  pruning_timing_seconds[["batch_selection_seconds"]] <-
    proc.time()[["elapsed"]] - phase_started

  # Print optional ecological diagnostics from the exact pre-removal state
  # used to score this batch. No additional frontier or patch-table pass is
  # performed.
  phase_started <- proc.time()[["elapsed"]]
  if (ecology_log_due(prune_iteration, ecology_log_every_iterations)) {
    if (is.null(initial_alive_cell_count) || is.null(current_alive_cell_count)) {
      stop("Ecological diagnostics require initial and current alive-cell counts.")
    }

    current_snapshot <- ecology_snapshot_from_cache(
      patch_score_cache = patch_score_cache,
      stage_index = stage_index,
      prune_iteration = prune_iteration,
      current_alive_cell_count = current_alive_cell_count,
      initial_alive_cell_count = initial_alive_cell_count
    )
    cutoff <- ecology_cutoff_cell(
      selected_cells = selected_cells,
      frontier_cells = frontier_cells,
      frontier_scores = frontier_scores,
      patch_score_cache = patch_score_cache,
      patch_id_by_species_env = patch_id_by_species_env
    )
    extremes <- ecology_extreme_pus(patch_score_cache)

    log_ecology_diagnostics(
      current_snapshot = current_snapshot,
      previous_snapshot = ecology_previous_snapshot,
      cutoff = cutoff,
      extremes = extremes
    )
    ecology_previous_snapshot <- current_snapshot
  }
  pruning_timing_seconds[["ecology_diagnostics_seconds"]] <-
    proc.time()[["elapsed"]] - phase_started
  phase_started <- proc.time()[["elapsed"]]


  # -------------------------------------------------------------------
  # STEP 3. Remove selected cells at the global cell level
  # -------------------------------------------------------------------

  # Store selected cells as integer cell IDs.
  newly_empty_cells <- as.integer(selected_cells)

  # Mark selected cells as globally removed.
  alive_species_count_by_cell[selected_cells] <- 0L

  # Read exact cell areas for patch-area decrementing.
  removed_cell_area <- cell_area_by_cell[selected_cells]


  # -------------------------------------------------------------------
  # STEP 4. Translate removed cells into removed area by patch
  # -------------------------------------------------------------------

  # Summarize selected-cell area by species-specific patch.
  removed_area_by_patch <- summarize_removed_patch_area(
    selected_cells = selected_cells,
    removed_cell_area = removed_cell_area,
    species_names = species_names,
    patch_id_by_species_env = patch_id_by_species_env,
    patch_score_cache = patch_score_cache
  )

  # Reconstruct the pre-step mapping only for original PUs containing a
  # directly touched patch. All later threshold cascades are confined to this
  # compact set of PUs.
  affected_patch_to_pu_before_step <- build_affected_patch_to_pu_mapping(
    touched_patches = removed_area_by_patch[, .(species, patch_id)],
    patch_score_cache = patch_score_cache,
    pu_graphs_by_key = pu_graphs_by_key
  )

  # Do not rewrite every full species vector here. Selected cells are already
  # globally dead, cannot re-enter the frontier, and cannot contribute to a
  # later score. Their stale patch IDs are cleared once at the pruning-stage
  # boundary. Threshold losses still update tables, graphs, counts, and the
  # frontier immediately; only their dense-vector clearing is deferred.


  # -------------------------------------------------------------------
  # STEP 5. Decrement patch areas
  # -------------------------------------------------------------------

  # Join removed-area rows to the pre-step patch -> PU mapping.
  patch_area_decrements <- removed_area_by_patch[
    affected_patch_to_pu_before_step,
    on = .(species, patch_id),
    nomatch = 0L
  ][
    ,
    .(species, patch_id, pu_id, area_removed)
  ]

  # Subtract removed area from affected patch-table rows.
  patch_table[
    patch_area_decrements,
    patch_area_km2 := patch_area_km2 - i.area_removed,
    on = .(species, patch_id, pu_id)
  ]

  # Record directly touched patches.
  touched_patches <- unique(
    patch_area_decrements[, .(species, patch_id)]
  )


  # -------------------------------------------------------------------
  # STEP 6. Drop patches below the patch-level threshold
  # -------------------------------------------------------------------

  # Identify touched patches that now fall below their species patch threshold.
  patches_dropped_by_patch_threshold <- find_pruning_patch_threshold_drops(
    patch_table = patch_table,
    touched_patches = touched_patches,
    species_params = species_params
  )

  # Remove patch-threshold failures from the patch table.
  patch_table <- patch_table[
    !patches_dropped_by_patch_threshold,
    on = .(species, patch_id)
  ]

  pruning_timing_seconds[["patch_update_seconds"]] <-
    proc.time()[["elapsed"]] - phase_started
  phase_started <- proc.time()[["elapsed"]]


  # -------------------------------------------------------------------
  # STEP 7. Repair affected PUs after patch losses
  # -------------------------------------------------------------------

  pu_repair_profile_started <- phase_started
  pu_repair_subphase_started <- phase_started

  # Identify original species-PU pairs containing newly dropped patches.
  affected_species_pu_pairs <- merge(
    patches_dropped_by_patch_threshold,
    affected_patch_to_pu_before_step,
    by = c("species", "patch_id")
  )

  component_repair_pus <- unique(
    affected_species_pu_pairs[, .(species, pu_id)],
    by = c("species", "pu_id")
  )
  component_repair_mapping <- affected_patch_to_pu_before_step[
    component_repair_pus,
    on = .(species, pu_id),
    nomatch = 0L
  ]

  # Freeze affected original graphs before any graph action is staged.
  affected_graph_context <- snapshot_pruning_affected_graphs(
    affected_patch_to_pu_before_step = component_repair_mapping,
    pu_graphs_by_key = pu_graphs_by_key
  )

  # Allocate collector for component-threshold patch drops.
  component_threshold_drop_records <- vector("list", 0L)

  # Allocate collector for final patch-to-PU assignments from rebuilt components.
  component_assignment_records <- vector("list", 0L)

  # Stage ordered graph removals and replacements without mutating live graphs.
  graph_actions <- vector("list", 0L)

  # Store each affected species' immutable pre-removal PU-ID high-water mark.
  pre_step_max_pu_by_species <- list()

  # Initialize count of new PUs created by splitting.
  new_pus_created <- 0L

  # Index current patch rows once for all species requiring component repair.
  repair_species <- unique(as.character(affected_species_pu_pairs$species))
  repair_row_mask <- patch_table$species %chin% repair_species
  repair_rows_by_species <- split(
    which(repair_row_mask),
    as.character(patch_table$species[repair_row_mask])
  )
  pu_repair_profile_counts[["species_pu_repairs"]] <-
    nrow(component_repair_pus)
  pu_repair_profile_timing[["setup_seconds"]] <-
    proc.time()[["elapsed"]] - pu_repair_subphase_started
  pu_repair_subphase_started <- proc.time()[["elapsed"]]

  # Process each species with patch-threshold drops.
  for (species_name in repair_species) {
    # Mark rows for this species.
    species_rows <- affected_species_pu_pairs$species == species_name

    # Split dropped patch IDs by original PU.
    dropped_patch_ids_by_pu <- split(
      affected_species_pu_pairs$patch_id[species_rows],
      affected_species_pu_pairs$pu_id[species_rows]
    )

    # Read this species' current rows once for all repaired PUs.
    species_patch_row_ids <- repair_rows_by_species[[species_name]]
    if (is.null(species_patch_row_ids)) {
      species_patch_row_ids <- integer(0L)
    }
    species_patch_rows <- patch_table[species_patch_row_ids]

    # Initialize new IDs above the immutable pre-removal species maximum.
    next_available_pu_id <- pruning_pre_step_pu_high_water(
      cache_entry = patch_score_cache[[species_name]],
      species_name = species_name
    )
    pre_step_max_pu_by_species[[species_name]] <- next_available_pu_id

    # Read this species' PU-area threshold.
    pu_area_threshold <- species_params[.(species_name), min_population_area_km2]

    # Build current patch_id -> patch_area lookup for this species.
    patch_area_lookup <- species_patch_rows[
      ,
      .(patch_id, patch_area_km2)
    ]

    # Assemble all affected PUs for this species once. The compiled kernel
    # preserves this exact PU order, node order, BFS order, and PU-ID sequence.
    pu_id_strings <- names(dropped_patch_ids_by_pu)
    species_repair_graphs <- lapply(
      pu_id_strings,
      function(pu_id_string) {
        affected_graph_context$graph_snapshot[[
          paste0(species_name, "|", pu_id_string)
        ]]
      }
    )
    species_repair_areas <- lapply(
      species_repair_graphs,
      function(pu_graph) {
        area_match <- fastmatch::fmatch(
          pu_graph$id2patch,
          patch_area_lookup$patch_id
        )
        node_patch_area <- numeric(length(pu_graph$id2patch))
        valid_area <- !is.na(area_match)
        node_patch_area[valid_area] <-
          patch_area_lookup$patch_area_km2[area_match[valid_area]]
        node_patch_area
      }
    )
    species_repair <- stage6_rebuild_species_pus_cpp(
      pu_graphs = species_repair_graphs,
      dead_patch_ids_by_pu = unname(dropped_patch_ids_by_pu[pu_id_strings]),
      patch_area_by_pu = species_repair_areas,
      pu_area_threshold = pu_area_threshold,
      next_available_pu_id = next_available_pu_id
    )
    next_available_pu_id <- species_repair$next_available_pu_id
    compiled_counts <- species_repair$profiling_counts
    pu_repair_profile_counts[names(compiled_counts)] <-
      pu_repair_profile_counts[names(compiled_counts)] + compiled_counts
    pu_repair_profile_counts[["compiled_species_batches"]] <-
      pu_repair_profile_counts[["compiled_species_batches"]] + 1

    # Process each affected original PU.
    for (pu_position in seq_along(pu_id_strings)) {
      pu_id_string <- pu_id_strings[[pu_position]]
      # Build graph key for this species-PU pair.
      graph_key <- paste0(species_name, "|", pu_id_string)

      # Read the compiled result corresponding to this original PU.
      pu_update <- species_repair$updates[[pu_position]]

      # Read patches dropped because surviving components were too small.
      component_dropped_patch_ids <- pu_update$dropped_patch_ids

      # Record component-threshold failures for one batched table removal.
      if (length(component_dropped_patch_ids)) {
        # Record component-threshold patch drops.
        component_threshold_drop_records[[length(component_threshold_drop_records) + 1L]] <-
          data.table(species = species_name, patch_id = component_dropped_patch_ids)

      }

      # Read surviving component graphs.
      surviving_pu_graphs <- pu_update$surviving_pu_graphs

      # Stage removal of the original graph if no component survives.
      if (!length(surviving_pu_graphs)) {
        graph_actions[[length(graph_actions) + 1L]] <- list(
          kind = "remove_original",
          key = graph_key,
          graph = NULL
        )
      } else {
        # Count new PUs if this PU split into multiple surviving components.
        if (length(surviving_pu_graphs) > 1L) {
          new_pus_created <- new_pus_created + (length(surviving_pu_graphs) - 1L)
        }

        # Store all surviving component graphs.
        for (j in seq_along(surviving_pu_graphs)) {
          # Read one surviving component graph.
          component_graph <- surviving_pu_graphs[[j]]

          staged_graph <- list(
            species = species_name,
            pu_id = as.integer(component_graph$pu_id),
            id2patch = as.integer(component_graph$id2patch),
            row_ptr = as.integer(component_graph$row_ptr),
            col_idx = as.integer(component_graph$col_idx)
          )

          # Stage the first component under the original graph key.
          if (j == 1L) {
            graph_actions[[length(graph_actions) + 1L]] <- list(
              kind = "replace_original",
              key = graph_key,
              graph = staged_graph
            )
          } else {
            # Build graph key for split-off component.
            new_key <- paste0(species_name, "|", component_graph$pu_id)

            # Stage the split-off component under its collision-free new key.
            graph_actions[[length(graph_actions) + 1L]] <- list(
              kind = "add_new",
              key = new_key,
              graph = staged_graph
            )
          }

          # Record final component membership for one batched update join.
          component_assignment_records[[length(component_assignment_records) + 1L]] <-
            data.table(
              species = species_name,
              patch_id = as.integer(component_graph$id2patch),
              pu_id = as.integer(component_graph$pu_id)
            )
        }
      }
    }
  }

  # Collapse component-threshold patch drops into one unique table.
  patches_dropped_by_component_threshold <- if (length(component_threshold_drop_records)) {
    unique(
      rbindlist(component_threshold_drop_records, use.names = TRUE, fill = TRUE),
      by = c("species", "patch_id")
    )
  } else {
    data.table(species = character(), patch_id = integer())
  }

  # Collapse final component assignments without changing their processing order.
  component_assignments <- if (length(component_assignment_records)) {
    rbindlist(component_assignment_records, use.names = TRUE, fill = TRUE)
  } else {
    data.table(species = character(), patch_id = integer(), pu_id = integer())
  }
  pu_repair_profile_timing[["component_rebuild_seconds"]] <-
    proc.time()[["elapsed"]] - pu_repair_subphase_started
  pu_repair_subphase_started <- proc.time()[["elapsed"]]
  pu_repair_graph_actions <- length(graph_actions)

  # Validate the entire affected graph/table transaction before live mutation.
  validate_pruning_graph_transaction(
    affected_patch_to_pu_before_step = component_repair_mapping,
    direct_patch_drops = patches_dropped_by_patch_threshold,
    component_drops = patches_dropped_by_component_threshold,
    component_assignments = component_assignments,
    graph_actions = graph_actions,
    pu_graphs_before = pu_graphs_by_key,
    pre_step_max_pu_by_species = pre_step_max_pu_by_species
  )

  # Apply all component removals and PU-ID updates to the global table once.
  patch_table <- apply_pruning_component_changes(
    patch_table = patch_table,
    component_assignments = component_assignments,
    component_drops = patches_dropped_by_component_threshold
  )

  # Commit the validated graph actions in their original processing order.
  pu_graphs_by_key <- priority_graph_apply_actions(
    pu_graphs_by_key, graph_actions
  )

  # Check only the affected rows and graphs after the in-memory commit.
  validate_pruning_component_commit(
    patch_table = patch_table,
    pu_graphs_by_key = pu_graphs_by_key,
    component_assignments = component_assignments,
    component_drops = patches_dropped_by_component_threshold,
    graph_actions = graph_actions
  )

  rm(
    repair_rows_by_species,
    repair_row_mask,
    component_assignments,
    component_repair_mapping,
    component_repair_pus,
    affected_graph_context,
    graph_actions,
    pre_step_max_pu_by_species
  )
  pu_repair_profile_timing[["transaction_seconds"]] <-
    proc.time()[["elapsed"]] - pu_repair_subphase_started
  pu_repair_subphase_started <- proc.time()[["elapsed"]]


  # -------------------------------------------------------------------
  # STEP 8. Drop whole PUs that now fail the PU-level threshold
  # -------------------------------------------------------------------

  # Build candidate patches whose area loss could make a PU subthreshold.
  candidate_patches <- unique(
    rbindlist(
      list(
        removed_area_by_patch[, .(species, patch_id)],
        patches_dropped_by_patch_threshold
      ),
      use.names = TRUE,
      fill = TRUE
    ),
    by = c("species", "patch_id")
  )

  # Identify original PUs containing those candidate patches.
  affected_pus_before <- unique(
    affected_patch_to_pu_before_step[
      candidate_patches,
      on = .(species, patch_id),
      nomatch = 0L
    ][
      ,
      .(species, pu_id)
    ],
    by = c("species", "pu_id")
  )

  # Recover original patches from affected original PUs, keeping only current survivors.
  affected_patches_that_still_exist <- affected_patch_to_pu_before_step[
    affected_pus_before,
    on = .(species, pu_id),
    nomatch = 0L
  ][
    ,
    .(species, patch_id)
  ][
    patch_table,
    on = .(species, patch_id),
    nomatch = 0L
  ]

  candidate_pus_after <- unique(
    affected_patches_that_still_exist[, .(species, pu_id)],
    by = c("species", "pu_id")
  )
  pu_repair_profile_counts[["threshold_candidate_pus"]] <-
    nrow(candidate_pus_after)

  # Recompute total area for only candidate PUs.
  pu_area_after_update <- patch_table[
    candidate_pus_after,
    .(pu_area_km2 = sum(patch_area_km2)),
    on = .(species, pu_id),
    by = .EACHI
  ]

  # Identify candidate PUs now below the PU-area threshold.
  pus_dropped_by_pu_threshold <- species_params[
    pu_area_after_update,
    on = .(species)
  ][
    !area_exceeds_threshold(pu_area_km2, min_population_area_km2),
    .(species, pu_id)
  ]
  pu_repair_profile_counts[["threshold_dropped_pus"]] <-
    nrow(pus_dropped_by_pu_threshold)
  pu_repair_graph_actions <-
    pu_repair_graph_actions + nrow(pus_dropped_by_pu_threshold)

  # Collect every patch belonging to a failed PU in one indexed join.
  pu_threshold_drop_rows <- collect_pruning_pu_drop_rows(
    patch_table = patch_table,
    dropped_pus = pus_dropped_by_pu_threshold
  )

  # Remove all failed PU graphs through one hashed-store transaction.
  if (nrow(pus_dropped_by_pu_threshold)) {
    pu_drop_actions <- lapply(
      paste0(
        pus_dropped_by_pu_threshold$species,
        "|",
        pus_dropped_by_pu_threshold$pu_id
      ),
      function(graph_key) list(kind = "remove", key = graph_key, graph = NULL)
    )
    pu_graphs_by_key <- priority_graph_apply_actions(
      pu_graphs_by_key,
      pu_drop_actions
    )
  }

  # Remove all failed PUs from the global patch table once.
  if (nrow(pus_dropped_by_pu_threshold)) {
    patch_table <- patch_table[
      !pus_dropped_by_pu_threshold,
      on = .(species, pu_id)
    ]
  }

  # Retain the unchanged two-column dropped-patch contract.
  patches_dropped_by_pu_threshold <- unique(
    pu_threshold_drop_rows[, .(species, patch_id)],
    by = c("species", "patch_id")
  )

  # Build one compact pending-drop list after every table-level repair is known.
  all_pending_patch_drops <- unique(
    rbindlist(
      list(
        patches_dropped_by_patch_threshold,
        patches_dropped_by_component_threshold,
        patches_dropped_by_pu_threshold
      ),
      use.names = TRUE,
      fill = TRUE
    ),
    by = c("species", "patch_id")
  )
  pending_patch_drops_by_species <- split(
    all_pending_patch_drops$patch_id,
    all_pending_patch_drops$species
  )
  pu_repair_profile_timing[["pu_threshold_seconds"]] <-
    proc.time()[["elapsed"]] - pu_repair_subphase_started

  pruning_timing_seconds[["pu_repair_seconds"]] <-
    proc.time()[["elapsed"]] - pu_repair_profile_started


  # -------------------------------------------------------------------
  # STEP 9. Apply accumulated patch drops to cell-level state
  # -------------------------------------------------------------------

  # Loop over species with possible pending dropped patches.
  for (species_name in names(pending_patch_drops_by_species)) {
    subphase_started <- proc.time()[["elapsed"]]

    # Deduplicate dropped patch IDs for this species.
    patch_ids_to_drop <- unique(pending_patch_drops_by_species[[species_name]])

    # Skip species with no pending dropped patches.
    if (!length(patch_ids_to_drop)) {
      pruning_timing_seconds[["drop_cell_lookup_seconds"]] <-
        pruning_timing_seconds[["drop_cell_lookup_seconds"]] +
        proc.time()[["elapsed"]] - subphase_started
      next
    }

    patch_ids_to_drop <- sort(unique(as.integer(patch_ids_to_drop)))

    # Keep only valid positive patch IDs.
    patch_ids_to_drop <- patch_ids_to_drop[
      !is.na(patch_ids_to_drop) &
        patch_ids_to_drop >= 1L
    ]

    # Skip if no valid patch IDs remain.
    if (!length(patch_ids_to_drop)) {
      pruning_timing_seconds[["drop_cell_lookup_seconds"]] <-
        pruning_timing_seconds[["drop_cell_lookup_seconds"]] +
        proc.time()[["elapsed"]] - subphase_started
      next
    }

    # Read species-specific patch_id -> cells index.
    patch_index <- get(
      species_name,
      envir = patch_cell_index_by_species_env,
      inherits = FALSE
    )

    # A dropped patch must remain represented in the compact index until the
    # pruning-stage boundary, when its dense cells are cleared once.
    if (is.null(patch_index)) {
      stop("Missing compact patch-cell index for threshold-dropped species: ", species_name)
    }

    # Look up all cells belonging to dropped patches for this species.
    dropped_patch_cells <- get_patch_cells_batch(
      patch_index = patch_index,
      patch_ids = patch_ids_to_drop
    )

    represented_patch_ids <- sort(unique(as.integer(dropped_patch_cells$patch_id)))
    missing_patch_ids <- setdiff(patch_ids_to_drop, represented_patch_ids)
    if (length(missing_patch_ids)) {
      stop(
        "Threshold-dropped patch IDs are missing from the compact index for species ",
        species_name, ": ", paste(utils::head(missing_patch_ids, 5L), collapse = ", ")
      )
    }

    # Deduplicate valid cell IDs and fail if internal state is structurally impossible.
    cells_to_clear <- validate_priority_cell_ids(
      dropped_patch_cells$cell,
      n_cells = length(alive_species_count_by_cell),
      label = paste0("Dropped-patch cells for species ", species_name)
    )

    # Skip if no valid cells remain.
    if (!length(cells_to_clear)) {
      pruning_timing_seconds[["drop_cell_lookup_seconds"]] <-
        pruning_timing_seconds[["drop_cell_lookup_seconds"]] +
        proc.time()[["elapsed"]] - subphase_started
      next
    }

    pruning_timing_seconds[["drop_cell_lookup_seconds"]] <-
      pruning_timing_seconds[["drop_cell_lookup_seconds"]] +
      proc.time()[["elapsed"]] - subphase_started
    subphase_started <- proc.time()[["elapsed"]]

    # Read global alive-species counts before removing this species.
    count_before_species_loss <- alive_species_count_by_cell[cells_to_clear]

    # Identify cells still globally alive before this species is removed.
    still_alive_cells <- cells_to_clear[count_before_species_loss > 0L]

    # Identify cells that become globally empty after this species is removed.
    cells_becoming_empty <- cells_to_clear[count_before_species_loss == 1L]

    # Decrement alive-species counts for cells that were still globally alive.
    if (length(still_alive_cells)) {
      alive_species_count_by_cell[still_alive_cells] <-
        alive_species_count_by_cell[still_alive_cells] - 1L
    }

    # Add cascading newly empty cells to the iteration-level list.
    if (length(cells_becoming_empty)) {
      newly_empty_cells <- c(newly_empty_cells, cells_becoming_empty)
    }

    pruning_work_counts[["cascade_empty_cells"]] <-
      pruning_work_counts[["cascade_empty_cells"]] +
      length(cells_becoming_empty)
    pruning_timing_seconds[["alive_count_update_seconds"]] <-
      pruning_timing_seconds[["alive_count_update_seconds"]] +
      proc.time()[["elapsed"]] - subphase_started

    # Physical dense-vector clearing is deferred to the stage boundary. The
    # refreshed score cache masks these dropped patch IDs in later iterations.
    pruning_work_counts[["threshold_species_updates_deferred"]] <-
      pruning_work_counts[["threshold_species_updates_deferred"]] + 1L
    pruning_work_counts[["threshold_cell_entries_deferred"]] <-
      pruning_work_counts[["threshold_cell_entries_deferred"]] +
      length(cells_to_clear)
  }


  # -------------------------------------------------------------------
  # STEP 10. Final bookkeeping for the stage-level orchestrator
  # -------------------------------------------------------------------

  subphase_started <- proc.time()[["elapsed"]]

  # Deduplicate and sort all cells that became globally empty in this iteration.
  newly_empty_cells <- sort(unique(as.integer(newly_empty_cells)))

  if (!is.null(current_alive_cell_count)) {
    current_alive_cell_count <- as.integer(
      current_alive_cell_count - length(newly_empty_cells)
    )
    if (current_alive_cell_count < 0L) {
      stop("Incremental alive-cell count became negative during pruning.")
    }
  }

  # Update incremental frontier state using only newly empty cells.
  frontier_state <- update_frontier_state_after_empty_cells(
    frontier_state = frontier_state,
    newly_empty_cells = newly_empty_cells,
    rook_neighbor_index = rook_neighbor_index
  )

  pruning_timing_seconds[["frontier_update_seconds"]] <-
    proc.time()[["elapsed"]] - subphase_started
  subphase_started <- proc.time()[["elapsed"]]

  # Collect every patch changed or removed during this pruning iteration.
  changed_patches <- unique(
    rbindlist(
      list(
        removed_area_by_patch[, .(species, patch_id)],
        patches_dropped_by_patch_threshold,
        patches_dropped_by_component_threshold,
        patches_dropped_by_pu_threshold
      ),
      use.names = TRUE,
      fill = TRUE
    ),
    by = c("species", "patch_id")
  )

  dirty_species_for_next_iteration <-
    sort(unique(as.character(changed_patches$species)))

  deferred_drop_species <- sort(unique(as.character(all_pending_patch_drops$species)))
  missing_dirty_species <- setdiff(
    deferred_drop_species,
    dirty_species_for_next_iteration
  )
  if (length(missing_dirty_species)) {
    stop(
      "Threshold-dropped species were not marked dirty for score-cache refresh: ",
      paste(utils::head(missing_dirty_species, 5L), collapse = ", ")
    )
  }


  # -------------------------------------------------------------------
  # STEP 11. Print exactly one pruning-iteration log line
  # -------------------------------------------------------------------

  # Print the single concise production log message for this iteration.
  log_pruning_iteration(
    stage_index = stage_index,
    prune_iteration = prune_iteration,
    frontier_cell_count = length(frontier_cells),
    selected_cell_count = length(selected_cells),
    newly_empty_cell_count = length(newly_empty_cells),
    touched_patch_count = nrow(touched_patches),
    dropped_patch_threshold_count = nrow(patches_dropped_by_patch_threshold),
    dropped_component_threshold_count = nrow(patches_dropped_by_component_threshold),
    dropped_pu_threshold_count = nrow(patches_dropped_by_pu_threshold),
    dropped_pu_count = nrow(pus_dropped_by_pu_threshold),
    new_pu_count = as.integer(new_pus_created),
    frontier_exhausted = FALSE
  )

  pruning_timing_seconds[["iteration_bookkeeping_seconds"]] <-
    proc.time()[["elapsed"]] - subphase_started


  # -------------------------------------------------------------------
  # STEP 12. Return updated state
  # -------------------------------------------------------------------

  # Return exactly the objects needed by the pruning-stage orchestrator.
  list(
    patch_table = patch_table,
    pu_graphs = pu_graphs_by_key,
    alive_species_count_by_cell = alive_species_count_by_cell,
    frontier_state = frontier_state,
    newly_empty_cells = newly_empty_cells,
    changed_patches = changed_patches,
    threshold_dropped_patches = all_pending_patch_drops[, .(species, patch_id)],
    patch_score_cache = patch_score_cache,
    dirty_species_for_scoring = dirty_species_for_next_iteration,
    ecology_previous_snapshot = ecology_previous_snapshot,
    current_alive_cell_count = current_alive_cell_count,
    pruning_timing_seconds = pruning_timing_seconds,
    pruning_work_counts = pruning_work_counts,
    frontier_scoring_counts = frontier_scoring_counts,
    frontier_scoring_timing = frontier_scoring_timing,
    pu_repair_profile_timing = pu_repair_profile_timing,
    pu_repair_profile_counts = pu_repair_profile_counts,
    pu_repair_graph_actions = pu_repair_graph_actions
  )
}
