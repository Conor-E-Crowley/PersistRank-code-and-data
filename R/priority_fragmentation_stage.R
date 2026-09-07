# Fragmentation-stage driver.
#
# Only pruning-changed patches are relabelled. Fragment IDs, retained original
# IDs, strict patch/PU threshold decisions, lazy graph overlays, and distance
# recheck IDs are committed in canonical species and patch order.

# Keep diagnostic products above R's integer limit finite. Both operands must
# be promoted before multiplication; converting an already-overflowed integer
# result to numeric would preserve NA.
fragmentation_dense_cell_workload <- function(n_species, n_cells) {
  as.numeric(n_species) * as.numeric(n_cells)
}

# After each pruning stage, only changed patches are relabeled. Newly split
# fragments are filtered by the minimum patch-size rule, affected persistence
# units are rebuilt, and patches that remain eligible for distance filtering are
# returned to the distance-connectivity stage.

log_fragmentation_stage <- function(
  stage_index,
  touched_species_count,
  input_changed_patch_count,
  removed_cell_count,
  recheck_patch_count,
  remaining_patch_count,
  remaining_pu_count,
  remaining_alive_cell_count
) {
  runtime_log_event(
    "fragmentation_done",
    stage = as.integer(stage_index),
    touched_species = as.integer(touched_species_count),
    input_changed_patches = as.integer(input_changed_patch_count),
    removed_cells = as.integer(removed_cell_count),
    recheck_patches = as.integer(recheck_patch_count),
    remaining_patches = as.integer(remaining_patch_count),
    remaining_pus = as.integer(remaining_pu_count),
    remaining_alive = as.integer(remaining_alive_cell_count)
  )

  invisible(NULL)
}

empty_fragmentation_timing <- function() {
  stats::setNames(
    numeric(9L),
    c(
      "component_relabel_seconds",
      "patch_threshold_seconds",
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

empty_fragmentation_counts <- function() {
  stats::setNames(
    integer(4L),
    c(
      "affected_pus",
      "unchanged_graph_pus",
      "added_node_only_pus",
      "full_rebuild_pus"
    )
  )
}

log_fragmentation_timing <- function(stage_index, timing_seconds, path_counts) {
  if (!is.numeric(timing_seconds) ||
      !identical(names(timing_seconds), names(empty_fragmentation_timing())) ||
      any(!is.finite(timing_seconds)) || any(timing_seconds < 0)) {
    stop("Fragmentation timings must be finite named nonnegative values.")
  }
  if (!is.numeric(path_counts) ||
      !identical(names(path_counts), names(empty_fragmentation_counts())) ||
      any(!is.finite(path_counts)) || any(path_counts < 0) ||
      any(path_counts != floor(path_counts))) {
    stop("Fragmentation path counts must be named nonnegative integers.")
  }
  displayed_timing <- round(timing_seconds, 2L)
  runtime_log_event(
    "fragmentation_timing",
    stage = as.integer(stage_index),
    component_relabel_seconds = displayed_timing[["component_relabel_seconds"]],
    patch_threshold_seconds = displayed_timing[["patch_threshold_seconds"]],
    pu_index_seconds = displayed_timing[["pu_index_seconds"]],
    pu_repair_seconds = displayed_timing[["pu_repair_seconds"]],
    graph_commit_seconds = displayed_timing[["graph_commit_seconds"]],
    compact_index_seconds = displayed_timing[["compact_index_seconds"]],
    cell_state_seconds = displayed_timing[["cell_state_seconds"]],
    state_commit_seconds = displayed_timing[["state_commit_seconds"]],
    other_seconds = displayed_timing[["other_seconds"]],
    affected_pus = as.integer(path_counts[["affected_pus"]]),
    unchanged_graph_pus = as.integer(path_counts[["unchanged_graph_pus"]]),
    added_node_only_pus = as.integer(path_counts[["added_node_only_pus"]]),
    full_rebuild_pus = as.integer(path_counts[["full_rebuild_pus"]]),
    total_seconds = sum(displayed_timing)
  )
  invisible(NULL)
}

log_fragmentation_index_profile <- function(stage_index, counts, timing_seconds) {
  if (!is.numeric(counts) ||
      !identical(names(counts), names(empty_fragmentation_index_counts())) ||
      any(!is.finite(counts)) || any(counts < 0) ||
      any(counts != floor(counts))) {
    stop("Fragmentation index profile counts must be named nonnegative integers.")
  }
  if (!is.numeric(timing_seconds) ||
      !identical(names(timing_seconds), names(empty_fragmentation_index_timing())) ||
      any(!is.finite(timing_seconds)) || any(timing_seconds < 0)) {
    stop("Fragmentation index profile timings must be finite named nonnegative values.")
  }
  displayed <- round(timing_seconds, 2L)
  classified <- sum(timing_seconds[c(
    "flatten_seconds", "kernel_seconds", "result_packaging_seconds"
  )])
  ratio <- if (counts[["replacement_entries"]] > 0) {
    counts[["input_index_entries"]] / counts[["replacement_entries"]]
  } else {
    0
  }
  runtime_log_event(
    "fragmentation_index_profile",
    stage = as.integer(stage_index),
    species_indexes_updated = as.integer(counts[["species_indexes_updated"]]),
    input_index_entries = as.numeric(counts[["input_index_entries"]]),
    changed_origin_segments = as.numeric(counts[["changed_origin_segments"]]),
    replacement_entries = as.numeric(counts[["replacement_entries"]]),
    output_index_entries = as.numeric(counts[["output_index_entries"]]),
    removed_entries = as.numeric(counts[["removed_entries"]]),
    relabel_entries = as.numeric(counts[["relabel_entries"]]),
    existing_to_replacement_ratio = round(ratio, 2L),
    flatten_seconds = displayed[["flatten_seconds"]],
    kernel_seconds = displayed[["kernel_seconds"]],
    result_packaging_seconds = displayed[["result_packaging_seconds"]],
    residual_seconds = round(max(
      0, timing_seconds[["total_seconds"]] - classified
    ), 2L),
    total_seconds = displayed[["total_seconds"]]
  )
  invisible(NULL)
}


# Replace every successfully processed species in one table operation. A NULL
# slot means that the species was skipped; a zero-row table is an intentional
# replacement for a species that lost all surviving patches.
replace_fragmentation_patch_rows <- function(
  patch_table,
  touched_species,
  replacement_rows_by_species
) {
  replacement_indices <- which(!vapply(
    replacement_rows_by_species,
    is.null,
    logical(1L)
  ))

  if (!length(replacement_indices)) {
    return(patch_table)
  }

  replaced_species <- touched_species[replacement_indices]
  untouched_rows <- patch_table[!(species %chin% replaced_species)]
  replacement_rows <- data.table::rbindlist(
    replacement_rows_by_species[replacement_indices],
    use.names = TRUE,
    fill = TRUE
  )

  data.table::rbindlist(
    list(untouched_rows, replacement_rows),
    use.names = TRUE,
    fill = TRUE
  )
}


# Validate the fused compiled component analysis before mutating live state.
validate_fragmentation_component_analysis <- function(
  result,
  changed_patch_ids,
  species_name
) {
  required <- c(
    "origin_patch_ids", "origin_offsets", "cells", "component_id_by_cell",
    "component_origin_patch_id", "component_id", "component_area_km2",
    "component_n_cells", "component_first_cell", "component_count_by_origin",
    "timing_seconds"
  )
  if (!is.list(result) || !identical(names(result), required)) {
    stop("Invalid compiled fragmentation-analysis result for species: ",
         species_name, call. = FALSE)
  }
  changed_patch_ids <- as.integer(changed_patch_ids)
  valid <- identical(as.integer(result$origin_patch_ids), changed_patch_ids) &&
    length(result$origin_offsets) == length(changed_patch_ids) + 1L &&
    identical(as.integer(result$origin_offsets[[1L]]), 0L) &&
    all(diff(result$origin_offsets) > 0L) &&
    tail(result$origin_offsets, 1L) == length(result$cells) &&
    length(result$cells) == length(result$component_id_by_cell) &&
    length(result$component_count_by_origin) == length(changed_patch_ids) &&
    all(result$component_count_by_origin > 0L)
  component_n <- sum(result$component_count_by_origin)
  valid <- valid && all(vapply(
    result[c(
      "component_origin_patch_id", "component_id", "component_area_km2",
      "component_n_cells", "component_first_cell"
    )],
    length,
    integer(1L)
  ) == component_n) &&
    all(is.finite(result$component_area_km2)) &&
    all(result$component_area_km2 > 0) &&
    all(result$component_n_cells > 0L) &&
    identical(
      as.integer(result$component_origin_patch_id),
      rep.int(changed_patch_ids, result$component_count_by_origin)
    ) &&
    identical(
      as.integer(result$component_id),
      sequence(result$component_count_by_origin)
    ) &&
    is.numeric(result$timing_seconds) &&
    identical(
      names(result$timing_seconds),
      c(
        "patch_cell_lookup_seconds", "component_label_seconds",
        "fragment_area_seconds"
      )
    ) &&
    all(is.finite(result$timing_seconds) & result$timing_seconds >= 0)
  if (!isTRUE(valid)) {
    stop("Invalid compiled fragmentation-analysis result for species: ",
         species_name, call. = FALSE)
  }
  result
}


analyze_fragmentation_components <- function(
  patch_index,
  changed_patch_ids,
  current_patch_id_by_cell,
  cell_area_by_cell,
  rook_neighbor_index,
  species_name,
  profile = FALSE
) {
  if (!exists(
    "stage6_analyze_fragmentation_components_cpp",
    mode = "function",
    inherits = TRUE
  )) {
    stop(
      "The compiled Stage 6 fragmentation-analysis kernel is not loaded.",
      call. = FALSE
    )
  }
  result <- stage6_analyze_fragmentation_components_cpp(
    index_pid = as.integer(patch_index$pid),
    index_cell = as.integer(patch_index$cell),
    changed_patch_ids = as.integer(changed_patch_ids),
    current_patch_id_by_cell = as.integer(current_patch_id_by_cell),
    cell_area_by_cell = as.numeric(cell_area_by_cell),
    rook_start = as.integer(rook_neighbor_index$start),
    rook_end = as.integer(rook_neighbor_index$end),
    rook_to = as.integer(rook_neighbor_index$to),
    n_cells = as.integer(rook_neighbor_index$n_cells),
    profile = isTRUE(profile)
  )
  validate_fragmentation_component_analysis(
    result,
    changed_patch_ids,
    species_name
  )
}

analyze_fragmentation_components_compact <- function(
  patch_index,
  changed_patch_ids,
  cell_area_by_cell,
  rook_neighbor_index,
  species_name,
  profile = FALSE
) {
  if (!exists(
    "stage72_analyze_fragmentation_components_cpp",
    mode = "function",
    inherits = TRUE
  )) {
    stop(
      "The compiled compact Stage 7.2 fragmentation kernel is not loaded.",
      call. = FALSE
    )
  }
  result <- stage72_analyze_fragmentation_components_cpp(
    index_pid = as.integer(patch_index$pid),
    index_cell = as.integer(patch_index$cell),
    changed_patch_ids = as.integer(changed_patch_ids),
    cell_area_by_cell = as.numeric(cell_area_by_cell),
    rook_start = as.integer(rook_neighbor_index$start),
    rook_end = as.integer(rook_neighbor_index$end),
    rook_to = as.integer(rook_neighbor_index$to),
    n_cells = as.integer(rook_neighbor_index$n_cells),
    profile = isTRUE(profile)
  )
  validate_fragmentation_component_analysis(
    result,
    changed_patch_ids,
    species_name
  )
}


# Run fragmentation for only the species and patches changed during the
# preceding pruning stage.
# It:
#
#   1. checks changed patches for rook-contiguity fragmentation
#   2. relabels newly disconnected fragments
#   3. drops fragments below the patch-size threshold
#   4. rebuilds only affected PU graph structures
#   5. drops PU components below the PU-size threshold
#   6. updates live environment-based species state
#   7. returns cells removed during fragmentation
#   8. returns surviving patches needing distance-connectivity recheck
#
# It does not write rasters or masks.
#

run_fragmentation_stage <- function(
  stage_index,
  changed_patches_in_stage,
  patch_table,
  pu_graphs_by_key,
  alive_species_count_by_cell,
  patch_id_by_species_env,
  patch_cell_index_by_species_env,
  species_params,
  cell_area_by_cell,
  rook_neighbor_index,
  state_representation = c("dense", "compact"),
  collect_runtime_diagnostics = FALSE,
  emit_logs = TRUE,
  progress_callback = NULL,
  completion_callback = NULL
) {
  fragmentation_stage_started <- proc.time()[["elapsed"]]
  state_representation <- match.arg(state_representation)
  compact_state <- identical(state_representation, "compact")
  if (!is.logical(collect_runtime_diagnostics) ||
      length(collect_runtime_diagnostics) != 1L ||
      is.na(collect_runtime_diagnostics)) {
    stop("collect_runtime_diagnostics must be TRUE or FALSE.", call. = FALSE)
  }
  # -------------------------------------------------------------------
  # 1. Standardize inputs
  # -------------------------------------------------------------------

  patch_table <- data.table::as.data.table(patch_table)

  species_params <- data.table::as.data.table(species_params)

  changed_patches_in_stage <- data.table::as.data.table(changed_patches_in_stage)


  # -------------------------------------------------------------------
  # 2. Validate the required rook-neighbor index
  # -------------------------------------------------------------------

  if (
    is.null(rook_neighbor_index) ||
      is.null(rook_neighbor_index$start) ||
      is.null(rook_neighbor_index$end) ||
      is.null(rook_neighbor_index$to) ||
      is.null(rook_neighbor_index$n_cells)
  ) {
    stop(
      "run_fragmentation_stage() requires a valid rook_neighbor_index. ",
      "Build it once with build_rook_neighbor_index() and pass it into ",
      "the fragmentation stage."
    )
  }


  # -------------------------------------------------------------------
  # 3. Initialize fragmentation-stage accumulators
  # -------------------------------------------------------------------

  # Store cells that become globally empty during fragmentation.
  removed_cells_in_fragmentation <- integer(0L)
  fragmentation_timing_seconds <- empty_fragmentation_timing()
  fragmentation_path_counts <- empty_fragmentation_counts()
  fragmentation_index_counts <- empty_fragmentation_index_counts()
  fragmentation_index_timing <- empty_fragmentation_index_timing()

  input_changed_patch_count <- nrow(changed_patches_in_stage)


  # -------------------------------------------------------------------
  # 4. Fast return if no patches changed during pruning
  # -------------------------------------------------------------------

  # If no changed patches were supplied, fragmentation has nothing to update.
  if (nrow(changed_patches_in_stage) == 0L) {
    # Count current patch rows.
    remaining_patch_count <- nrow(patch_table)

    # Count current unique species-PU combinations.
    remaining_pu_count <- data.table::uniqueN(patch_table, by = c("species", "pu_id"))

    # Count current globally alive cells.
    remaining_alive_cell_count <- sum(alive_species_count_by_cell > 0L)

    if (isTRUE(emit_logs)) log_fragmentation_stage(
      stage_index = stage_index,
      touched_species_count = 0L,
      input_changed_patch_count = 0L,
      removed_cell_count = 0L,
      recheck_patch_count = 0L,
      remaining_patch_count = remaining_patch_count,
      remaining_pu_count = remaining_pu_count,
      remaining_alive_cell_count = remaining_alive_cell_count
    )
    if (isTRUE(emit_logs)) log_fragmentation_timing(
      stage_index,
      fragmentation_timing_seconds,
      fragmentation_path_counts
    )
    if (isTRUE(emit_logs)) log_fragmentation_index_profile(
      stage_index,
      fragmentation_index_counts,
      fragmentation_index_timing
    )

    # Return unchanged state and empty fragmentation bookkeeping.
    return(list(
      patch_table = patch_table,
      pu_graphs = pu_graphs_by_key,
      alive_species_count_by_cell = alive_species_count_by_cell,
      removed_cells_in_fragmentation = integer(0L),
      patches_requiring_distance_recheck = data.table::data.table(
        species = character(),
        patch_id = integer()
      ),
      touched_species = character(),
      timing_seconds = fragmentation_timing_seconds,
      path_counts = fragmentation_path_counts,
      index_timing_seconds = fragmentation_index_timing
    ))
  }


  # -------------------------------------------------------------------
  # 5. Restrict changed patches to species available in the live env
  # -------------------------------------------------------------------

  # Read species currently represented in the live cell -> patch environment.
  species_available_in_env <- sort(ls(
    envir = if (compact_state) {
      patch_cell_index_by_species_env
    } else {
      patch_id_by_species_env
    },
    all.names = TRUE
  ))

  # Keep only changed patches whose species has a live environment vector.
  changed_patches_in_stage <- changed_patches_in_stage[
    species %in% species_available_in_env
  ]

  # Identify species touched by the changed patches that remain after filtering.
  touched_species <- sort(unique(as.character(changed_patches_in_stage$species)))

  # If filtering removed all touched species, no fragmentation update is possible.
  if (!length(touched_species)) {
    # Count current patch rows.
    remaining_patch_count <- nrow(patch_table)

    # Count current unique species-PU combinations.
    remaining_pu_count <- data.table::uniqueN(patch_table, by = c("species", "pu_id"))

    # Count current globally alive cells.
    remaining_alive_cell_count <- sum(alive_species_count_by_cell > 0L)

    if (isTRUE(emit_logs)) log_fragmentation_stage(
      stage_index = stage_index,
      touched_species_count = 0L,
      input_changed_patch_count = input_changed_patch_count,
      removed_cell_count = 0L,
      recheck_patch_count = 0L,
      remaining_patch_count = remaining_patch_count,
      remaining_pu_count = remaining_pu_count,
      remaining_alive_cell_count = remaining_alive_cell_count
    )
    if (isTRUE(emit_logs)) log_fragmentation_timing(
      stage_index,
      fragmentation_timing_seconds,
      fragmentation_path_counts
    )
    if (isTRUE(emit_logs)) log_fragmentation_index_profile(
      stage_index,
      fragmentation_index_counts,
      fragmentation_index_timing
    )

    # Return unchanged state and empty fragmentation bookkeeping.
    return(list(
      patch_table = patch_table,
      pu_graphs = pu_graphs_by_key,
      alive_species_count_by_cell = alive_species_count_by_cell,
      removed_cells_in_fragmentation = integer(0L),
      patches_requiring_distance_recheck = data.table::data.table(
        species = character(),
        patch_id = integer()
      ),
      touched_species = character(),
      timing_seconds = fragmentation_timing_seconds,
      path_counts = fragmentation_path_counts,
      index_timing_seconds = fragmentation_index_timing
    ))
  }


  # -------------------------------------------------------------------
  # 6. Index immutable phase inputs and allocate outputs
  # -------------------------------------------------------------------

  # Index the immutable patch table once instead of scanning it per species.
  patch_rows_by_species <- split(
    seq_len(nrow(patch_table)),
    as.character(patch_table$species)
  )

  # Index the compact changed-patch input once for direct species lookup.
  changed_patch_ids_by_species <- split(
    changed_patches_in_stage$patch_id,
    as.character(changed_patches_in_stage$species)
  )

  # Create one output slot per touched species.
  recheck_patch_rows_by_species <- vector("list", length(touched_species))

  # Store updated patch rows for one batched global-table replacement. NULL
  # means that a species was skipped; a zero-row table is a valid replacement.
  replacement_patch_rows_by_species <- vector("list", length(touched_species))
  profiling_enabled <- is.function(completion_callback)


  # ===================================================================
  # 7. Process each touched species
  # ===================================================================

  # Loop over species touched by pruning-stage patch changes.
  for (species_index in seq_along(touched_species)) {
    species_name <- touched_species[species_index]
    if (is.function(progress_callback)) {
      progress_callback(species_index, length(touched_species), species_name)
    }
    species_started <- if (profiling_enabled) proc.time()[["elapsed"]] else NA_real_
    species_timing_before <- if (profiling_enabled) {
      fragmentation_timing_seconds
    } else {
      NULL
    }
    changed_patch_cells_inspected <- 0L
    profile_unchanged_pus <- 0L
    profile_added_only_pus <- 0L
    profile_full_rebuild_pus <- 0L
    profile_dropped_pus <- 0L
    profile_graph_nodes <- 0
    profile_graph_adjacency_entries <- 0
    profile_overlay_edges <- 0
    profile_overlay_seconds <- 0
    profile_provisional_graph_seconds <- 0
    profile_component_rebuild_seconds <- 0
    profile_full_rebuild_input_nodes <- 0
    profile_full_rebuild_input_adjacency_entries <- 0
    profile_provisional_candidate_edges <- 0
    profile_provisional_unique_edges <- 0
    profile_provisional_adjacency_entries <- 0
    profile_patch_cell_lookup_seconds <- 0
    profile_component_label_seconds <- 0
    profile_fragment_area_seconds <- 0
    profile_fragment_assignment_seconds <- 0
    profile_component_row_assembly_seconds <- 0
    profile_components_found <- 0L
    profile_split_origin_patches <- 0L
    profile_largest_changed_patch_cells <- 0L
    profile_max_components_in_one_patch <- 0L

    # Read changed patch IDs for this species.
    changed_patch_ids_for_species <- sort(unique(as.integer(
      changed_patch_ids_by_species[[species_name]]
    )))

    species_patch_ids_before <- if (compact_state) {
      NULL
    } else {
      get(
        species_name,
        envir = patch_id_by_species_env,
        inherits = FALSE
      )
    }

    # Read current patch-table rows for this species.
    species_patch_row_ids <- patch_rows_by_species[[species_name]]
    if (is.null(species_patch_row_ids)) {
      species_patch_row_ids <- integer(0L)
    }
    species_patch_rows_before <- patch_table[species_patch_row_ids]

    # Skip this species if it has no surviving patch rows.
    if (nrow(species_patch_rows_before) == 0L) {
      next
    }

    # Keep only changed patches that still survive in this species' patch rows.
    changed_patch_ids_for_species <- intersect(
      changed_patch_ids_for_species,
      species_patch_rows_before$patch_id
    )

    # Skip this species if none of its changed patches survive.
    if (!length(changed_patch_ids_for_species)) {
      next
    }


    # -----------------------------------------------------------------
    # 7A. Initialize patch origin map for this species
    # -----------------------------------------------------------------

    # Build a table mapping current patches to their origin patch and PU.
    patch_origin_map <- species_patch_rows_before[
      ,
      list(
        patch_id = as.integer(patch_id),
        origin_patch_id = as.integer(patch_id),
        pu_id = as.integer(pu_id)
      )
    ]

    # Initialize new patch IDs above the current species-specific maximum.
    next_available_patch_id <- max(species_patch_rows_before$patch_id)

    # Read this species' patch -> cells index from the live environment.
    species_patch_cell_index_before <- get(
      species_name,
      envir = patch_cell_index_by_species_env,
      inherits = FALSE
    )

    if (is.null(species_patch_cell_index_before)) {
      stop(
        "Missing patch-cell index for species '",
        species_name,
        "' during fragmentation repair."
      )
    }

    # Accumulate exact fragment rows and compact-index replacements from the
    # component cells already visited by the scientific repair.
    fragment_patch_row_records <- vector("list", 0L)
    replacement_index_rows_by_origin <- list()


    # -----------------------------------------------------------------
    # 7B. Pass 1: relabel changed patches by rook contiguity
    # -----------------------------------------------------------------

    subphase_started <- proc.time()[["elapsed"]]
    component_analysis <- if (compact_state) {
      analyze_fragmentation_components_compact(
        patch_index = species_patch_cell_index_before,
        changed_patch_ids = changed_patch_ids_for_species,
        cell_area_by_cell = cell_area_by_cell,
        rook_neighbor_index = rook_neighbor_index,
        species_name = species_name,
        profile = profiling_enabled
      )
    } else {
      analyze_fragmentation_components(
        patch_index = species_patch_cell_index_before,
        changed_patch_ids = changed_patch_ids_for_species,
        current_patch_id_by_cell = species_patch_ids_before,
        cell_area_by_cell = cell_area_by_cell,
        rook_neighbor_index = rook_neighbor_index,
        species_name = species_name,
        profile = profiling_enabled
      )
    }
    if (profiling_enabled) {
      profile_patch_cell_lookup_seconds <-
        component_analysis$timing_seconds[["patch_cell_lookup_seconds"]]
      profile_component_label_seconds <-
        component_analysis$timing_seconds[["component_label_seconds"]]
      profile_fragment_area_seconds <-
        component_analysis$timing_seconds[["fragment_area_seconds"]]
      profile_components_found <- sum(
        component_analysis$component_count_by_origin
      )
      profile_split_origin_patches <- sum(
        component_analysis$component_count_by_origin > 1L
      )
      profile_largest_changed_patch_cells <- max(diff(
        component_analysis$origin_offsets
      ))
      profile_max_components_in_one_patch <- max(
        component_analysis$component_count_by_origin
      )
    }
    changed_patch_cells_inspected <- length(component_analysis$cells)

    patch_row_positions <- match(
      changed_patch_ids_for_species,
      species_patch_rows_before$patch_id
    )
    if (anyNA(patch_row_positions)) {
      stop("Changed fragmentation patches are missing species rows for: ",
           species_name, call. = FALSE)
    }
    expected_origin_areas <- species_patch_rows_before$patch_area_km2[
      patch_row_positions
    ]
    inherited_pu_ids <- species_patch_rows_before$pu_id[patch_row_positions]
    component_cursor <- 0L

    # Allocate final fragment IDs in the established origin/component order.
    for (origin_index in seq_along(changed_patch_ids_for_species)) {
      fine_started <- if (profiling_enabled) proc.time()[["elapsed"]] else NA_real_
      original_patch_id <- changed_patch_ids_for_species[[origin_index]]
      cell_positions <- seq.int(
        component_analysis$origin_offsets[[origin_index]] + 1L,
        component_analysis$origin_offsets[[origin_index + 1L]]
      )
      patch_cells <- component_analysis$cells[cell_positions]
      fragment_labels_in_patch_cells <-
        component_analysis$component_id_by_cell[cell_positions]
      component_count <-
        component_analysis$component_count_by_origin[[origin_index]]
      component_positions <- component_cursor + seq_len(component_count)
      component_cursor <- component_cursor + component_count
      unique_fragment_labels <- seq_len(component_count)
      fragment_area_table <- data.table::data.table(
        fragment_label = as.integer(component_analysis$component_id[
          component_positions
        ]),
        fragment_area_km2 = as.numeric(
          component_analysis$component_area_km2[component_positions]
        )
      )

      expected_origin_area <- expected_origin_areas[[origin_index]]
      observed_origin_area <- sum(fragment_area_table$fragment_area_km2)
      area_tolerance <- 1e-8 * max(1, abs(expected_origin_area))
      if (!is.finite(expected_origin_area) ||
          abs(observed_origin_area - expected_origin_area) > area_tolerance) {
        stop(
          "Fragment cell areas do not conserve current patch area for ",
          species_name, "|", original_patch_id,
          ": cells=", format(observed_origin_area, digits = 16),
          ", patch_table=", format(expected_origin_area, digits = 16), ".",
          call. = FALSE
        )
      }

      keep_fragment_label <- fragment_area_table[
        order(-fragment_area_km2, fragment_label)
      ]$fragment_label[1L]
      inherited_pu_id <- inherited_pu_ids[[origin_index]]
      if (is.na(inherited_pu_id)) {
        stop(
          "Could not find inherited PU id for species '", species_name,
          "', patch_id ", original_patch_id,
          " during fragmentation repair."
        )
      }
      if (profiling_enabled) {
        profile_fragment_area_seconds <-
          profile_fragment_area_seconds +
          proc.time()[["elapsed"]] - fine_started
        fine_started <- proc.time()[["elapsed"]]
      }

      final_patch_id_by_component <- rep.int(
        as.integer(original_patch_id),
        component_count
      )
      for (fragment_label in unique_fragment_labels) {
        fragment_patch_id <- as.integer(original_patch_id)
        if (fragment_label != keep_fragment_label) {
          next_available_patch_id <- next_available_patch_id + 1L
          fragment_patch_id <- as.integer(next_available_patch_id)
          final_patch_id_by_component[[fragment_label]] <- fragment_patch_id
        }
        fragment_patch_row_records[[length(fragment_patch_row_records) + 1L]] <-
          data.table::data.table(
            species = species_name,
            patch_id = fragment_patch_id,
            pu_id = as.integer(inherited_pu_id),
            patch_area_km2 = fragment_area_table$fragment_area_km2[[fragment_label]],
            origin_patch_id = as.integer(original_patch_id),
            first_cell = as.integer(
              component_analysis$component_first_cell[
                component_positions[[fragment_label]]
              ]
            )
          )
      }

      replacement_index_rows_by_origin[[as.character(original_patch_id)]] <-
        data.table::data.table(
          origin_patch_id = as.integer(original_patch_id),
          pid = as.integer(final_patch_id_by_component[
            fragment_labels_in_patch_cells
          ]),
          cell = as.integer(patch_cells)
        )
      if (profiling_enabled) {
        profile_fragment_assignment_seconds <-
          profile_fragment_assignment_seconds +
          proc.time()[["elapsed"]] - fine_started
      }
    }

    # -----------------------------------------------------------------
    # 7C. Reconstruct patch rows from the component results
    # -----------------------------------------------------------------

    fine_started <- if (profiling_enabled) proc.time()[["elapsed"]] else NA_real_
    fragment_patch_rows <- data.table::rbindlist(
      fragment_patch_row_records,
      use.names = TRUE,
      fill = TRUE
    )
    species_patch_rows_current <- data.table::copy(species_patch_rows_before)
    species_patch_rows_current[, origin_patch_id := as.integer(patch_id)]
    original_fragment_rows <- fragment_patch_rows[patch_id == origin_patch_id]
    species_patch_rows_current[
      original_fragment_rows,
      patch_area_km2 := i.patch_area_km2,
      on = .(patch_id)
    ]
    new_fragment_rows <- fragment_patch_rows[patch_id != origin_patch_id]
    if (nrow(new_fragment_rows)) {
      species_patch_rows_current <- data.table::rbindlist(
        list(
          species_patch_rows_current,
          new_fragment_rows[, .(
            species,
            patch_id,
            pu_id,
            patch_area_km2,
            origin_patch_id
          )]
        ),
        use.names = TRUE,
        fill = TRUE
      )
    }
    patch_origin_map <- species_patch_rows_current[
      ,
      .(patch_id, origin_patch_id, pu_id)
    ]
    if (profiling_enabled) {
      profile_component_row_assembly_seconds <-
        profile_component_row_assembly_seconds +
        proc.time()[["elapsed"]] - fine_started
    }
    fragmentation_timing_seconds[["component_relabel_seconds"]] <-
      fragmentation_timing_seconds[["component_relabel_seconds"]] +
      proc.time()[["elapsed"]] - subphase_started


    # -----------------------------------------------------------------
    # 7D. Pass 2: enforce minimum patch-size threshold
    # -----------------------------------------------------------------

    subphase_started <- proc.time()[["elapsed"]]

    # Read this species' minimum patch-size threshold.
    patch_area_threshold <- species_params[
      species == species_name
    ]$min_patch_area_km2[1L]

    # Identify current patches below the patch-size threshold.
    dropped_patch_ids_by_fragmentation <- species_patch_rows_current[
      !area_exceeds_threshold(patch_area_km2, patch_area_threshold)
    ]$patch_id

    # Track patch losses in preallocated per-PU slots and concatenate once.
    dropped_patch_id_records <- vector("list", 1L)
    dropped_patch_id_records[[1L]] <- as.integer(dropped_patch_ids_by_fragmentation)
    dropped_patch_record_cursor <- 1L

    # Remove below-threshold rows directly. Their exact cells remain available
    # in the compact pre-fragmentation index.
    if (length(dropped_patch_ids_by_fragmentation)) {
      species_patch_rows_current <- species_patch_rows_current[
        !(patch_id %in% dropped_patch_ids_by_fragmentation)
      ]
    }
    fragmentation_timing_seconds[["patch_threshold_seconds"]] <-
      fragmentation_timing_seconds[["patch_threshold_seconds"]] +
      proc.time()[["elapsed"]] - subphase_started
    subphase_started <- proc.time()[["elapsed"]]


    # -----------------------------------------------------------------
    # 7E. Identify affected PUs
    # -----------------------------------------------------------------

    # Affected PUs are those containing changed origin patches.
    affected_pu_ids <- sort(unique(
      patch_origin_map$pu_id[
        patch_origin_map$origin_patch_id %in% changed_patch_ids_for_species
      ]
    ))

    # Resize the drop collector now that the exact affected-PU count is known.
    length(dropped_patch_id_records) <- length(affected_pu_ids) + 1L

    # Keep existing species patch rows from PUs unaffected by fragmentation.
    unaffected_species_patch_rows <- species_patch_rows_before[
      !(pu_id %in% affected_pu_ids),
      list(
        species = species,
        patch_id = patch_id,
        pu_id = pu_id,
        patch_area_km2 = patch_area_km2
      )
    ]

    # Keep current patch rows from affected PUs only.
    affected_species_patch_rows_current <- species_patch_rows_current[
      pu_id %in% affected_pu_ids
    ]
    affected_species_patch_rows_current[, final_order__ := NA_integer_]

    # Index affected rows once so the PU loop does not repeatedly scan this
    # species' complete affected patch table.
    affected_patch_rows_by_pu <- index_affected_patch_rows_by_pu(
      affected_species_patch_rows_current,
      affected_pu_ids
    )

    # Collect only real graph changes. Unchanged graphs remain bound in the
    # hashed store; their historical remove/re-add ordering is reproduced by a
    # lightweight order-only commit after repair.
    graph_keys_to_remove <- character(length(affected_pu_ids))
    graph_remove_cursor <- 0L
    rebuilt_graphs_for_affected_pus <- list()
    final_affected_graph_keys <- character(nrow(affected_species_patch_rows_current))
    final_graph_key_cursor <- 0L
    final_row_order_cursor <- 0L

    next_available_pu_id <- if (nrow(species_patch_rows_before)) {
      max(species_patch_rows_before$pu_id)
    } else {
      0L
    }

    pu_area_threshold <- species_params[
      species == species_name
    ]$min_population_area_km2[1L]
    fragmentation_path_counts[["affected_pus"]] <-
      fragmentation_path_counts[["affected_pus"]] + length(affected_pu_ids)
    fragmentation_timing_seconds[["pu_index_seconds"]] <-
      fragmentation_timing_seconds[["pu_index_seconds"]] +
      proc.time()[["elapsed"]] - subphase_started
    subphase_started <- proc.time()[["elapsed"]]


    # -----------------------------------------------------------------
    # 7F. Pass 3: rebuild only affected PUs
    # -----------------------------------------------------------------

    # Loop over affected original PUs.
    for (affected_pu_id in affected_pu_ids) {
      # Read integer row positions without materializing a PU-sized table.
      current_row_positions_in_this_pu <-
        affected_patch_rows_by_pu[[as.character(affected_pu_id)]]
      current_patch_ids <- as.integer(
        affected_species_patch_rows_current$patch_id[current_row_positions_in_this_pu]
      )

      # Build the old graph key.
      old_graph_key <- paste0(species_name, "|", affected_pu_id)

      # If all patches in this old PU disappeared, there is no graph to rebuild.
      if (!length(current_row_positions_in_this_pu)) {
        if (profiling_enabled) profile_dropped_pus <- profile_dropped_pus + 1L
        graph_remove_cursor <- graph_remove_cursor + 1L
        graph_keys_to_remove[[graph_remove_cursor]] <- old_graph_key
        next
      }

      old_pu_graph <- priority_graph_get(pu_graphs_by_key, old_graph_key)

      # A surviving PU must have a corresponding old graph.
      if (is.null(old_pu_graph)) {
        stop(
          "Missing PU graph during fragmentation-stage update for species '",
          species_name,
          "' and pu_id ",
          affected_pu_id,
          ". This indicates inconsistent evolving state between patch_table ",
          "and pu_graphs_by_key."
        )
      }

      if (profiling_enabled) {
        graph_storage <- fragmentation_graph_storage_counts(old_pu_graph)
        profile_graph_nodes <- profile_graph_nodes + graph_storage[["nodes"]]
        profile_graph_adjacency_entries <-
          profile_graph_adjacency_entries +
          graph_storage[["adjacency_entries"]]
      }

      # Read old graph patch IDs.
      old_patch_ids <- as.integer(old_pu_graph$id2patch)

      if (anyNA(old_patch_ids) || any(old_patch_ids < 1L) ||
          anyDuplicated(old_patch_ids)) {
        stop("Invalid patch IDs in fragmentation graph: ", old_graph_key)
      }
      if (anyNA(current_patch_ids) || any(current_patch_ids < 1L) ||
          anyDuplicated(current_patch_ids)) {
        stop("Invalid current patch IDs during fragmentation: ", old_graph_key)
      }

      # Classify topology directly from indexed vectors. Sorting is deferred
      # until a genuinely changed graph needs the historical ordered sets.
      current_in_old <- match(current_patch_ids, old_patch_ids, nomatch = 0L)
      old_in_current <- match(old_patch_ids, current_patch_ids, nomatch = 0L)
      new_patch_nodes_unsorted <- current_patch_ids[current_in_old == 0L]
      removed_old_patch_nodes_unsorted <- old_patch_ids[old_in_current == 0L]

      # Determine whether the PU graph node set changed.
      pu_patch_node_set_changed <-
        length(new_patch_nodes_unsorted) > 0L ||
        length(removed_old_patch_nodes_unsorted) > 0L

      # Compute current total area for this PU.
      current_pu_area_km2 <- sum(
        affected_species_patch_rows_current$patch_area_km2[
          current_row_positions_in_this_pu
        ]
      )


      # ---------------------------------------------------------------
      # Fast path 1: no patch nodes were added or removed
      # ---------------------------------------------------------------

      # If the node set did not change, graph topology can be reused.
      if (!pu_patch_node_set_changed) {
        if (profiling_enabled) {
          profile_unchanged_pus <- profile_unchanged_pus + 1L
        }
        fragmentation_path_counts[["unchanged_graph_pus"]] <-
          fragmentation_path_counts[["unchanged_graph_pus"]] + 1L
        # If the PU is now below threshold, remove all of its current patches.
        if (!area_exceeds_threshold(current_pu_area_km2, pu_area_threshold)) {
          dropped_patch_record_cursor <- dropped_patch_record_cursor + 1L
          dropped_patch_id_records[[dropped_patch_record_cursor]] <- current_patch_ids
          graph_remove_cursor <- graph_remove_cursor + 1L
          graph_keys_to_remove[[graph_remove_cursor]] <- old_graph_key
        } else {
          # Keep the graph binding and current rows untouched. Record only the
          # final row/key order that the former remove-and-readd path produced.
          row_count <- length(current_row_positions_in_this_pu)
          affected_species_patch_rows_current$final_order__[
            current_row_positions_in_this_pu
          ] <- final_row_order_cursor + seq_len(row_count)
          final_row_order_cursor <- final_row_order_cursor + row_count
          final_graph_key_cursor <- final_graph_key_cursor + 1L
          final_affected_graph_keys[[final_graph_key_cursor]] <- old_graph_key
        }

        # Continue to the next affected PU.
        next
      }


      # ---------------------------------------------------------------
      # Fast path 2: added nodes only, no old nodes removed
      # ---------------------------------------------------------------

      # Check whether fragmentation only added patch nodes.
      only_added_nodes_no_old_nodes_removed <-
        length(new_patch_nodes_unsorted) > 0L &&
        length(removed_old_patch_nodes_unsorted) == 0L

      # If only new nodes were added, the old graph cannot have disconnected.
      if (only_added_nodes_no_old_nodes_removed) {
        if (profiling_enabled) {
          profile_added_only_pus <- profile_added_only_pus + 1L
        }
        fragmentation_path_counts[["added_node_only_pus"]] <-
          fragmentation_path_counts[["added_node_only_pus"]] + 1L
        # If the PU is below threshold, remove all current patches in this PU.
        if (!area_exceeds_threshold(current_pu_area_km2, pu_area_threshold)) {
          dropped_patch_record_cursor <- dropped_patch_record_cursor + 1L
          dropped_patch_id_records[[dropped_patch_record_cursor]] <- current_patch_ids
          graph_remove_cursor <- graph_remove_cursor + 1L
          graph_keys_to_remove[[graph_remove_cursor]] <- old_graph_key
        } else {
          current_rows_in_this_pu <- affected_species_patch_rows_current[
            current_row_positions_in_this_pu
          ]
          new_patch_nodes <- sort(as.integer(new_patch_nodes_unsorted))
          # Build a lazy overlay graph containing only edges for added nodes.
          overlay_started <- if (profiling_enabled) {
            proc.time()[["elapsed"]]
          } else {
            NA_real_
          }
          lazy_overlay_graph <- make_lazy_added_node_overlay_pu_graph(
            current_pu_patch_rows = current_rows_in_this_pu,
            old_pu_graph = old_pu_graph,
            new_patch_nodes = new_patch_nodes
          )
          if (profiling_enabled) {
            profile_overlay_seconds <- profile_overlay_seconds +
              proc.time()[["elapsed"]] - overlay_started
            profile_overlay_edges <- profile_overlay_edges +
              nrow(lazy_overlay_graph$overlay_edges)
          }

          rebuilt_graphs_for_affected_pus[[old_graph_key]] <- lazy_overlay_graph
          row_count <- length(current_row_positions_in_this_pu)
          affected_species_patch_rows_current$final_order__[
            current_row_positions_in_this_pu
          ] <- final_row_order_cursor + seq_len(row_count)
          final_row_order_cursor <- final_row_order_cursor + row_count
          final_graph_key_cursor <- final_graph_key_cursor + 1L
          final_affected_graph_keys[[final_graph_key_cursor]] <- old_graph_key
        }

        # Continue to the next affected PU.
        next
      }


      # ---------------------------------------------------------------
      # Full path: old patch nodes were removed
      # ---------------------------------------------------------------

      # Build the full provisional graph for current patches in this PU.
      if (profiling_enabled) {
        profile_full_rebuild_pus <- profile_full_rebuild_pus + 1L
      }
      fragmentation_path_counts[["full_rebuild_pus"]] <-
        fragmentation_path_counts[["full_rebuild_pus"]] + 1L
      if (profiling_enabled) {
        profile_full_rebuild_input_nodes <-
          profile_full_rebuild_input_nodes + length(old_pu_graph$id2patch)
        profile_full_rebuild_input_adjacency_entries <-
          profile_full_rebuild_input_adjacency_entries +
          length(old_pu_graph$col_idx)
      }
      current_rows_in_this_pu <- affected_species_patch_rows_current[
        current_row_positions_in_this_pu
      ]
      graph_remove_cursor <- graph_remove_cursor + 1L
      graph_keys_to_remove[[graph_remove_cursor]] <- old_graph_key
      provisional_graph_started <- if (profiling_enabled) {
        proc.time()[["elapsed"]]
      } else {
        NA_real_
      }
      provisional_result <- build_provisional_pu_graph(
        current_pu_patch_rows = current_rows_in_this_pu,
        old_pu_graph = old_pu_graph,
        return_diagnostics = profiling_enabled
      )
      if (profiling_enabled) {
        profile_provisional_graph_seconds <-
          profile_provisional_graph_seconds +
          proc.time()[["elapsed"]] - provisional_graph_started
        provisional_pu_graph <- provisional_result$graph
        provisional_diagnostics <- provisional_result$diagnostics
        profile_provisional_candidate_edges <-
          profile_provisional_candidate_edges +
          provisional_diagnostics[["candidate_undirected_edges"]]
        profile_provisional_unique_edges <-
          profile_provisional_unique_edges +
          provisional_diagnostics[["unique_undirected_edges"]]
        profile_provisional_adjacency_entries <-
          profile_provisional_adjacency_entries +
          provisional_diagnostics[["output_adjacency_entries"]]
      } else {
        provisional_pu_graph <- provisional_result
      }

      # Rebuild this PU after patch loss and possible graph disconnection.
      component_rebuild_started <- if (profiling_enabled) {
        proc.time()[["elapsed"]]
      } else {
        NA_real_
      }
      rebuilt_pu <- rebuild_pu_after_patch_loss(
        pu_graph = provisional_pu_graph,
        alive_nodes = rep(TRUE, nrow(current_rows_in_this_pu)),
        patch_area = current_rows_in_this_pu$patch_area_km2,
        pu_area_threshold = pu_area_threshold,
        next_available_pu_id = next_available_pu_id
      )
      if (profiling_enabled) {
        profile_component_rebuild_seconds <-
          profile_component_rebuild_seconds +
          proc.time()[["elapsed"]] - component_rebuild_started
      }

      # Carry forward the next available PU ID after possible splitting.
      next_available_pu_id <- rebuilt_pu$next_available_pu_id

      # If any components are too small, clear their patches from the species vector.
      if (length(rebuilt_pu$dropped_patch_ids)) {
        dropped_patch_record_cursor <- dropped_patch_record_cursor + 1L
        dropped_patch_id_records[[dropped_patch_record_cursor]] <-
          as.integer(rebuilt_pu$dropped_patch_ids)
      }

      # Store surviving PU graph objects.
      for (surviving_graph in rebuilt_pu$surviving_pu_graphs) {
        # Build graph key for this surviving PU.
        surviving_graph_key <- paste0(species_name, "|", surviving_graph$pu_id)

        rebuilt_graphs_for_affected_pus[[surviving_graph_key]] <- list(
          species = species_name,
          pu_id = as.integer(surviving_graph$pu_id),
          id2patch = as.integer(surviving_graph$id2patch),
          row_ptr = as.integer(surviving_graph$row_ptr),
          col_idx = as.integer(surviving_graph$col_idx)
        )
        final_graph_key_cursor <- final_graph_key_cursor + 1L
        final_affected_graph_keys[[final_graph_key_cursor]] <- surviving_graph_key
      }

      # Store patch rows for surviving PU components with updated PU IDs.
      for (surviving_graph in rebuilt_pu$surviving_pu_graphs) {
        # Read patch IDs in this surviving component.
        component_patch_ids <- as.integer(surviving_graph$id2patch)

        # Extract current rows for these component patches.
        component_row_positions <- match(
          component_patch_ids,
          current_rows_in_this_pu$patch_id
        )
        if (anyNA(component_row_positions)) {
          stop(
            "A rebuilt fragmentation component references a missing patch row for ",
            species_name, "|", affected_pu_id, ".",
            call. = FALSE
          )
        }
        component_global_positions <- current_row_positions_in_this_pu[
          component_row_positions
        ]
        if (length(component_global_positions)) {
          affected_species_patch_rows_current$pu_id[component_global_positions] <-
            as.integer(surviving_graph$pu_id)
          row_count <- length(component_global_positions)
          affected_species_patch_rows_current$final_order__[
            component_global_positions
          ] <- final_row_order_cursor + seq_len(row_count)
          final_row_order_cursor <- final_row_order_cursor + row_count
        }
      }
    }

    graph_keys_to_remove <- graph_keys_to_remove[seq_len(graph_remove_cursor)]
    final_affected_graph_keys <-
      final_affected_graph_keys[seq_len(final_graph_key_cursor)]
    all_dropped_patch_ids <- unique(as.integer(unlist(
      dropped_patch_id_records[seq_len(dropped_patch_record_cursor)],
      use.names = FALSE
    )))


    # -----------------------------------------------------------------
    # 7G. Combine unaffected and rebuilt patch rows for this species
    # -----------------------------------------------------------------

    # Materialize affected rows once in the exact historical PU/component
    # order. Unchanged PUs reached this point without per-PU row copies.
    rebuilt_species_patch_rows <- affected_species_patch_rows_current[
      !is.na(final_order__)
    ][order(final_order__)]
    rebuilt_species_patch_rows[, final_order__ := NULL]

    # Combine unaffected rows and rebuilt affected rows, retaining origin IDs.
    final_species_patch_rows_with_origin <- data.table::rbindlist(
      list(
        unaffected_species_patch_rows[
          ,
          list(
            species = species,
            patch_id = patch_id,
            pu_id = pu_id,
            patch_area_km2 = patch_area_km2,
            origin_patch_id = patch_id
          )
        ],
        rebuilt_species_patch_rows[
          ,
          list(
            species = species,
            patch_id = patch_id,
            pu_id = pu_id,
            patch_area_km2 = patch_area_km2,
            origin_patch_id = origin_patch_id
          )
        ]
      ),
      use.names = TRUE,
      fill = TRUE
    )
    fragmentation_timing_seconds[["pu_repair_seconds"]] <-
      fragmentation_timing_seconds[["pu_repair_seconds"]] +
      proc.time()[["elapsed"]] - subphase_started


    # -----------------------------------------------------------------
    # 7H. Stage this species' final patch rows
    # -----------------------------------------------------------------

    # Store the final four-column table for one global replacement after all
    # touched species have been processed. Storing an empty table intentionally
    # removes a species that lost every surviving patch.
    replacement_patch_rows_by_species[[species_index]] <-
      final_species_patch_rows_with_origin[
        ,
        list(
          species = species,
          patch_id = patch_id,
          pu_id = pu_id,
          patch_area_km2 = patch_area_km2
        )
      ]


    # -----------------------------------------------------------------
    # 7I. Replace affected PU graphs in the global graph list
    # -----------------------------------------------------------------

    if (anyDuplicated(final_affected_graph_keys)) {
      stop("Fragmentation produced duplicate final affected graph keys.", call. = FALSE)
    }

    # Commit only genuinely changed graph bindings. A separate order-only
    # update preserves the established named-list serialization contract.
    graph_actions <- c(
      lapply(
        unique(graph_keys_to_remove),
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
    subphase_started <- proc.time()[["elapsed"]]


    # -----------------------------------------------------------------
    # 7J. Update alive-species counts for cells where this species vanished
    # -----------------------------------------------------------------

    # Update the compact index using only changed/stale patch segments. This
    # returns the exact cells removed without scanning the full species vector.
    final_species_patch_ids <- sort(unique(as.integer(
      final_species_patch_rows_with_origin$patch_id
    )))
    retained_dropped_patch_ids <- intersect(
      all_dropped_patch_ids,
      final_species_patch_ids
    )
    if (length(retained_dropped_patch_ids)) {
      stop(
        "A fragmentation-dropped patch remains in the final species rows: ",
        paste(utils::head(retained_dropped_patch_ids, 5L), collapse = ", "),
        call. = FALSE
      )
    }
    index_update <- update_species_patch_index_incremental(
      patch_index = species_patch_cell_index_before,
      replacement_rows_by_origin = replacement_index_rows_by_origin,
      expected_patch_ids = final_species_patch_ids,
      pre_update_patch_ids = species_patch_rows_before$patch_id,
      return_diagnostics = TRUE
    )
    fragmentation_index_counts <- fragmentation_index_counts +
      index_update$diagnostic_counts
    fragmentation_index_timing <- fragmentation_index_timing +
      index_update$diagnostic_timing_seconds
    fragmentation_timing_seconds[["compact_index_seconds"]] <-
      fragmentation_timing_seconds[["compact_index_seconds"]] +
      proc.time()[["elapsed"]] - subphase_started

    # All table, graph, and compact-index changes are now validated. Commit
    # only the affected graph bindings and reproduce established key order.
    subphase_started <- proc.time()[["elapsed"]]
    pu_graphs_by_key <- priority_graph_apply_actions(
      pu_graphs_by_key,
      graph_actions
    )
    pu_graphs_by_key <- priority_graph_reposition_keys_to_end(
      pu_graphs_by_key,
      final_affected_graph_keys
    )
    fragmentation_timing_seconds[["graph_commit_seconds"]] <-
      fragmentation_timing_seconds[["graph_commit_seconds"]] +
      proc.time()[["elapsed"]] - subphase_started
    subphase_started <- proc.time()[["elapsed"]]
    species_cells_removed_in_fragmentation <- if (compact_state) {
      as.integer(index_update$removed_cells)
    } else {
      index_update$removed_cells[
        !is.na(species_patch_ids_before[index_update$removed_cells])
      ]
    }

    # Apply all local relabels and removals in one copy-on-write update.
    if (!compact_state && (
      nrow(index_update$relabel_rows) ||
        length(species_cells_removed_in_fragmentation)
    )) {
      species_patch_ids_after <- species_patch_ids_before
      if (nrow(index_update$relabel_rows)) {
        species_patch_ids_after[index_update$relabel_rows$cell] <-
          index_update$relabel_rows$patch_id
      }
      if (length(index_update$removed_cells)) {
        species_patch_ids_after[index_update$removed_cells] <- NA_integer_
      }
      assign(
        species_name,
        species_patch_ids_after,
        envir = patch_id_by_species_env
      )
    }

    # Read global alive-species counts before removing this species from those cells.
    count_before_species_loss <- alive_species_count_by_cell[
      species_cells_removed_in_fragmentation
    ]

    # Identify cells that remain globally alive after losing this species.
    cells_still_alive_after_species_loss <- species_cells_removed_in_fragmentation[
      count_before_species_loss > 0L
    ]

    # Identify cells that become globally empty after losing this species.
    cells_becoming_empty <- species_cells_removed_in_fragmentation[
      count_before_species_loss == 1L
    ]

    # Decrement alive-species counts for cells that were still globally alive.
    if (length(cells_still_alive_after_species_loss)) {
      alive_species_count_by_cell[cells_still_alive_after_species_loss] <-
        alive_species_count_by_cell[cells_still_alive_after_species_loss] - 1L
    }

    # Add newly empty cells to the fragmentation-stage accumulator.
    if (length(cells_becoming_empty)) {
      removed_cells_in_fragmentation <- c(
        removed_cells_in_fragmentation,
        cells_becoming_empty
      )
    }


    # -----------------------------------------------------------------
    # 7K. Write updated species state back to live environments
    # -----------------------------------------------------------------

    assign(
      species_name,
      index_update$index,
      envir = patch_cell_index_by_species_env
    )


    # -----------------------------------------------------------------
    # 7L. Record surviving patches requiring distance recheck
    # -----------------------------------------------------------------

    # Store patches whose geometry or identity changed and survived this stage.
    recheck_patch_rows_by_species[[species_index]] <- unique(
      final_species_patch_rows_with_origin[
        origin_patch_id %in% changed_patch_ids_for_species |
          patch_id != origin_patch_id,
        list(
          species = species,
          patch_id = patch_id
        )
      ],
      by = c("species", "patch_id")
    )
    fragmentation_timing_seconds[["cell_state_seconds"]] <-
      fragmentation_timing_seconds[["cell_state_seconds"]] +
      proc.time()[["elapsed"]] - subphase_started

    if (profiling_enabled) {
      species_elapsed <- proc.time()[["elapsed"]] - species_started
      species_timing <- fragmentation_timing_seconds - species_timing_before
      broad_component_seconds <- species_timing[["component_relabel_seconds"]]
      component_classified_seconds <- sum(c(
        profile_patch_cell_lookup_seconds,
        profile_component_label_seconds,
        profile_fragment_area_seconds,
        profile_fragment_assignment_seconds,
        profile_component_row_assembly_seconds
      ))
      profile_component_other_seconds <- max(
        0,
        broad_component_seconds - component_classified_seconds
      )
      component_seconds <- component_classified_seconds +
        profile_component_other_seconds
      patch_threshold_seconds <- species_timing[["patch_threshold_seconds"]]
      pu_index_seconds <- species_timing[["pu_index_seconds"]]
      pu_repair_seconds <- species_timing[["pu_repair_seconds"]]
      compact_index_seconds <- species_timing[["compact_index_seconds"]]
      graph_commit_seconds <- species_timing[["graph_commit_seconds"]]
      cell_state_seconds <- species_timing[["cell_state_seconds"]]
      profile_full_rebuild_seconds <-
        profile_provisional_graph_seconds +
        profile_component_rebuild_seconds
      pu_other_seconds <- max(
        0,
        pu_repair_seconds - profile_overlay_seconds -
          profile_full_rebuild_seconds
      )
      classified_seconds <- sum(c(
        component_seconds,
        patch_threshold_seconds,
        pu_index_seconds,
        profile_overlay_seconds,
        profile_provisional_graph_seconds,
        profile_component_rebuild_seconds,
        pu_other_seconds,
        compact_index_seconds,
        graph_commit_seconds,
        cell_state_seconds
      ))
      other_seconds <- max(0, species_elapsed - classified_seconds)
      reported_elapsed <- classified_seconds + other_seconds
      completion_callback(
        species_index,
        length(touched_species),
        species_name,
        changed_patches = length(changed_patch_ids_for_species),
        changed_patch_cells = changed_patch_cells_inspected,
        affected_pus = length(affected_pu_ids),
        unchanged_pus = profile_unchanged_pus,
        added_only_pus = profile_added_only_pus,
        full_rebuild_pus = profile_full_rebuild_pus,
        dropped_pus = profile_dropped_pus,
        affected_graph_nodes = profile_graph_nodes,
        affected_graph_adjacency_entries = profile_graph_adjacency_entries,
        overlay_edges_created = profile_overlay_edges,
        full_rebuild_input_nodes = profile_full_rebuild_input_nodes,
        full_rebuild_input_adjacency_entries =
          profile_full_rebuild_input_adjacency_entries,
        provisional_candidate_edges = profile_provisional_candidate_edges,
        provisional_unique_edges = profile_provisional_unique_edges,
        provisional_adjacency_entries =
          profile_provisional_adjacency_entries,
        new_fragments = nrow(new_fragment_rows),
        components_found = profile_components_found,
        split_origin_patches = profile_split_origin_patches,
        largest_changed_patch_cells = profile_largest_changed_patch_cells,
        max_components_in_one_patch = profile_max_components_in_one_patch,
        patch_cell_lookup_seconds = profile_patch_cell_lookup_seconds,
        component_label_seconds = profile_component_label_seconds,
        fragment_area_seconds = profile_fragment_area_seconds,
        fragment_assignment_seconds = profile_fragment_assignment_seconds,
        component_row_assembly_seconds =
          profile_component_row_assembly_seconds,
        component_other_seconds = profile_component_other_seconds,
        component_seconds = component_seconds,
        patch_threshold_seconds = patch_threshold_seconds,
        pu_index_seconds = pu_index_seconds,
        overlay_seconds = profile_overlay_seconds,
        provisional_graph_seconds = profile_provisional_graph_seconds,
        component_rebuild_seconds = profile_component_rebuild_seconds,
        full_rebuild_seconds = profile_full_rebuild_seconds,
        pu_other_seconds = pu_other_seconds,
        compact_index_seconds = compact_index_seconds,
        graph_commit_seconds = graph_commit_seconds,
        cell_state_seconds = cell_state_seconds,
        other_seconds = other_seconds,
        elapsed_seconds = reported_elapsed
      )
    }
  }

  rm(patch_rows_by_species, changed_patch_ids_by_species)


  # -------------------------------------------------------------------
  # 8. Replace all successfully processed species in one table operation
  # -------------------------------------------------------------------

  # Preserve the prior sequential result order: untouched original rows first,
  # followed by replacement rows in sorted touched-species processing order.
  subphase_started <- proc.time()[["elapsed"]]
  patch_table <- replace_fragmentation_patch_rows(
    patch_table = patch_table,
    touched_species = touched_species,
    replacement_rows_by_species = replacement_patch_rows_by_species
  )
  fragmentation_timing_seconds[["state_commit_seconds"]] <-
    fragmentation_timing_seconds[["state_commit_seconds"]] +
    proc.time()[["elapsed"]] - subphase_started


  # -------------------------------------------------------------------
  # 9. Finalize removed cells and distance-recheck patch table
  # -------------------------------------------------------------------

  # Deduplicate and sort cells that became globally empty in fragmentation.
  removed_cells_in_fragmentation <- sort(
    unique(as.integer(removed_cells_in_fragmentation))
  )

  # Keep non-empty distance-recheck tables.
  non_empty_recheck_tables <- Filter(
    f = function(x) is.data.frame(x) && nrow(x) > 0L,
    x = recheck_patch_rows_by_species
  )

  # Combine distance-recheck tables, or return an empty table.
  patches_requiring_distance_recheck <- if (length(non_empty_recheck_tables)) {
    # Combine and deduplicate recheck rows.
    unique(
      data.table::rbindlist(non_empty_recheck_tables, use.names = TRUE, fill = TRUE),
      by = c("species", "patch_id")
    )
  } else {
    # Return an empty recheck table.
    data.table::data.table(
      species = character(),
      patch_id = integer()
    )
  }


  # -------------------------------------------------------------------
  # 10. Log one summary line for this fragmentation stage
  # -------------------------------------------------------------------

  # Count remaining patch-table rows.
  remaining_patch_count <- nrow(patch_table)

  # Count remaining unique species-PU combinations.
  remaining_pu_count <- data.table::uniqueN(patch_table, by = c("species", "pu_id"))

  # Count cells still containing at least one surviving species.
  remaining_alive_cell_count <- sum(alive_species_count_by_cell > 0L)

  classified_fields <- setdiff(
    names(fragmentation_timing_seconds),
    "other_seconds"
  )
  fragmentation_timing_seconds[["other_seconds"]] <- max(
    0,
    proc.time()[["elapsed"]] - fragmentation_stage_started -
      sum(fragmentation_timing_seconds[classified_fields])
  )
  if (isTRUE(emit_logs)) log_fragmentation_stage(
    stage_index = stage_index,
    touched_species_count = length(touched_species),
    input_changed_patch_count = input_changed_patch_count,
    removed_cell_count = length(removed_cells_in_fragmentation),
    recheck_patch_count = nrow(patches_requiring_distance_recheck),
    remaining_patch_count = remaining_patch_count,
    remaining_pu_count = remaining_pu_count,
    remaining_alive_cell_count = remaining_alive_cell_count
  )
  if (isTRUE(emit_logs)) log_fragmentation_timing(
    stage_index,
    fragmentation_timing_seconds,
    fragmentation_path_counts
  )
  if (isTRUE(emit_logs)) log_fragmentation_index_profile(
    stage_index,
    fragmentation_index_counts,
    fragmentation_index_timing
  )


  # -------------------------------------------------------------------
  # 11. Return updated live state
  # -------------------------------------------------------------------

  # Return the updated state needed by the outer pipeline and distance stage.
  list(
    patch_table = patch_table,
    pu_graphs = pu_graphs_by_key,
    alive_species_count_by_cell = alive_species_count_by_cell,
    removed_cells_in_fragmentation = removed_cells_in_fragmentation,
    patches_requiring_distance_recheck = patches_requiring_distance_recheck,
    touched_species = touched_species,
    timing_seconds = fragmentation_timing_seconds,
    path_counts = fragmentation_path_counts,
    index_timing_seconds = fragmentation_index_timing,
    workload = if (collect_runtime_diagnostics) c(
      state_representation_compact = as.numeric(compact_state),
      affected_species = as.numeric(length(touched_species)),
      changed_patches = as.numeric(input_changed_patch_count),
      changed_index_entries = as.numeric(
        fragmentation_index_counts[["replacement_entries"]]
      ),
      full_index_entries_scanned = if (compact_state) 0 else as.numeric(
        sum(vapply(
          touched_species,
          function(species_name) {
            length(get(
              species_name,
              envir = patch_cell_index_by_species_env,
              inherits = FALSE
            )$cell)
          },
          numeric(1L)
        ))
      ),
      dense_vectors_materialized = if (compact_state) 0 else as.numeric(
        length(touched_species)
      ),
      dense_cells_materialized = if (compact_state) 0 else
        fragmentation_dense_cell_workload(
          length(touched_species), rook_neighbor_index$n_cells
        ),
      dense_cells_avoided = if (compact_state)
        fragmentation_dense_cell_workload(
          length(touched_species), rook_neighbor_index$n_cells
        ) else 0,
      wrapper_seconds = as.numeric(
        sum(fragmentation_timing_seconds)
      )
    ) else NULL
  )
}
