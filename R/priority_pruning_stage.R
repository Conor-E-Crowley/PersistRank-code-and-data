# Repeated sparse-frontier pruning within one global stage.
#
# Dense cell state is updated incrementally and empty-cell clearing is deferred
# until the safe synchronization points already established by the optimized
# implementation. Returned changed patches are the only fragmentation input.


# Print one pruning-stage log message.
#
# This helper prints exactly one summary message per completed pruning stage.
#
# It does not perform any expensive spatial reconstruction. All values passed
# to it are already available from the final stage state or from the stage
# accumulators.
#

log_pruning_stage <- function(
  stage_index,                      # global stage index
  completed_pruning_iterations,      # number of pruning iterations actually completed
  removed_cell_count,                # number of unique cells removed during this pruning stage
  changed_patch_count,               # number of unique patches changed during this pruning stage
  frontier_exhausted,                # whether pruning stopped because the frontier was exhausted
  remaining_patch_count,             # number of patch-table rows remaining after pruning
  remaining_pu_count,                # number of species-PU combinations remaining after pruning
  remaining_alive_cell_count         # number of globally alive cells remaining after pruning
) {
  runtime_log_event(
    "pruning_done",
    stage = as.integer(stage_index),
    iters = as.integer(completed_pruning_iterations),
    removed_cells = as.integer(removed_cell_count),
    changed_patches = as.integer(changed_patch_count),
    frontier_exhausted = isTRUE(frontier_exhausted),
    remaining_patches = as.integer(remaining_patch_count),
    remaining_pus = as.integer(remaining_pu_count),
    remaining_alive = as.integer(remaining_alive_cell_count)
  )

  invisible(NULL)
}


log_pruning_timing <- function(
  stage_index,
  completed_pruning_iterations,
  timing_seconds,
  work_counts = empty_pruning_work_counts(),
  species_vectors_flushed = 0L
) {
  expected_names <- names(empty_pruning_timing())
  if (!is.numeric(timing_seconds) ||
      !identical(names(timing_seconds), expected_names) ||
      any(!is.finite(timing_seconds)) ||
      any(timing_seconds < 0)) {
    stop("Pruning subphase timings must be finite named nonnegative values.")
  }
  if (!is.numeric(work_counts) ||
      !identical(names(work_counts), names(empty_pruning_work_counts())) ||
      any(!is.finite(work_counts)) || any(work_counts < 0L) ||
      any(work_counts != floor(work_counts))) {
    stop("Pruning workload counts must be named nonnegative integers.")
  }

  displayed_timing <- round(timing_seconds, 2)
  score_select_seconds <- sum(displayed_timing[c(
    "frontier_fetch_seconds",
    "score_cache_seconds",
    "frontier_scoring_seconds",
    "batch_selection_seconds",
    "ecology_diagnostics_seconds"
  )])
  cell_update_fields <- c(
    "drop_cell_lookup_seconds",
    "alive_count_update_seconds",
    "frontier_update_seconds",
    "iteration_bookkeeping_seconds"
  )
  cell_update_seconds <- sum(displayed_timing[cell_update_fields])

  runtime_log_event(
    "pruning_timing",
    stage = as.integer(stage_index),
    iters = as.integer(completed_pruning_iterations),
    score_select_seconds = score_select_seconds,
    frontier_fetch_seconds = displayed_timing[["frontier_fetch_seconds"]],
    score_cache_seconds = displayed_timing[["score_cache_seconds"]],
    frontier_scoring_seconds = displayed_timing[["frontier_scoring_seconds"]],
    batch_selection_seconds = displayed_timing[["batch_selection_seconds"]],
    ecology_diagnostics_seconds = displayed_timing[["ecology_diagnostics_seconds"]],
    patch_update_seconds = displayed_timing[["patch_update_seconds"]],
    pu_repair_seconds = displayed_timing[["pu_repair_seconds"]],
    cell_update_seconds = cell_update_seconds,
    drop_cell_lookup_seconds = displayed_timing[["drop_cell_lookup_seconds"]],
    alive_count_update_seconds = displayed_timing[["alive_count_update_seconds"]],
    frontier_update_seconds = displayed_timing[["frontier_update_seconds"]],
    iteration_bookkeeping_seconds = displayed_timing[["iteration_bookkeeping_seconds"]],
    deferred_cell_clear_seconds = displayed_timing[["deferred_cell_clear_seconds"]],
    threshold_species_updates_deferred = as.integer(
      work_counts[["threshold_species_updates_deferred"]]
    ),
    threshold_cell_entries_deferred = as.integer(
      work_counts[["threshold_cell_entries_deferred"]]
    ),
    cascade_empty_cells = as.integer(work_counts[["cascade_empty_cells"]]),
    species_vectors_flushed = as.integer(species_vectors_flushed),
    total_seconds = sum(displayed_timing)
  )

  invisible(NULL)
}

log_pruning_profiles <- function(
  stage_index,
  completed_pruning_iterations,
  frontier_counts,
  frontier_timing,
  frontier_total_seconds,
  pu_repair_counts,
  pu_repair_timing,
  pu_repair_total_seconds,
  pu_repair_graph_actions
) {
  if (!is.numeric(frontier_counts) ||
      !identical(names(frontier_counts), names(empty_frontier_scoring_counts())) ||
      any(!is.finite(frontier_counts)) || any(frontier_counts < 0) ||
      any(frontier_counts != floor(frontier_counts))) {
    stop("Frontier scoring profile counts must be named nonnegative integers.")
  }
  if (!is.numeric(frontier_timing) ||
      !identical(names(frontier_timing), names(empty_frontier_scoring_timing())) ||
      any(!is.finite(frontier_timing)) || any(frontier_timing < 0) ||
      length(frontier_total_seconds) != 1L ||
      !is.finite(frontier_total_seconds) || frontier_total_seconds < 0) {
    stop("Frontier scoring profile timings must be finite named nonnegative values.")
  }
  if (!is.numeric(pu_repair_counts) ||
      !identical(names(pu_repair_counts), names(empty_pruning_pu_repair_counts())) ||
      any(!is.finite(pu_repair_counts)) || any(pu_repair_counts < 0) ||
      any(pu_repair_counts != floor(pu_repair_counts))) {
    stop("PU repair profile counts must be named nonnegative integers.")
  }
  if (!is.numeric(pu_repair_timing) ||
      !identical(names(pu_repair_timing), names(empty_pruning_pu_repair_timing())) ||
      any(!is.finite(pu_repair_timing)) || any(pu_repair_timing < 0) ||
      length(pu_repair_total_seconds) != 1L ||
      !is.finite(pu_repair_total_seconds) || pu_repair_total_seconds < 0 ||
      length(pu_repair_graph_actions) != 1L ||
      !is.finite(pu_repair_graph_actions) || pu_repair_graph_actions < 0 ||
      pu_repair_graph_actions != floor(pu_repair_graph_actions)) {
    stop("PU repair profile timings and graph-action count are invalid.")
  }

  frontier_display <- round(frontier_timing, 2L)
  frontier_classified <- sum(frontier_timing[c(
    "input_setup_seconds", "kernel_seconds", "result_finalize_seconds"
  )])
  if (frontier_classified > frontier_total_seconds + 1e-8) {
    stop("Frontier scoring profile exceeds the measured scoring total.")
  }
  kernel_measured <- sum(frontier_timing[c(
    "kernel_scan_seconds", "kernel_accumulation_seconds",
    "kernel_output_seconds"
  )])
  mean_contributions <- if (frontier_counts[["scored_frontier_cells"]] > 0) {
    frontier_counts[["finite_contributions"]] /
      frontier_counts[["scored_frontier_cells"]]
  } else {
    0
  }
  runtime_log_event(
    "frontier_scoring_profile",
    stage = as.integer(stage_index),
    iters = as.integer(completed_pruning_iterations),
    score_calls = as.integer(frontier_counts[["score_calls"]]),
    frontier_cells = as.numeric(frontier_counts[["frontier_cells"]]),
    active_species = as.numeric(frontier_counts[["active_species"]]),
    species_frontier_probes = as.numeric(
      frontier_counts[["species_frontier_probes"]]
    ),
    sparse_membership_probes = as.numeric(
      frontier_counts[["sparse_membership_probes"]]
    ),
    dense_probes_avoided = as.numeric(
      frontier_counts[["dense_probes_avoided"]]
    ),
    finite_contributions = as.numeric(frontier_counts[["finite_contributions"]]),
    scored_frontier_cells = as.numeric(
      frontier_counts[["scored_frontier_cells"]]
    ),
    mean_contributions_per_scored_cell = round(mean_contributions, 2L),
    input_setup_seconds = frontier_display[["input_setup_seconds"]],
    kernel_seconds = frontier_display[["kernel_seconds"]],
    result_finalize_seconds = frontier_display[["result_finalize_seconds"]],
    kernel_scan_seconds = frontier_display[["kernel_scan_seconds"]],
    kernel_accumulation_seconds =
      frontier_display[["kernel_accumulation_seconds"]],
    kernel_output_seconds = frontier_display[["kernel_output_seconds"]],
    kernel_residual_seconds = round(max(
      0, frontier_timing[["kernel_seconds"]] - kernel_measured
    ), 2L),
    wrapper_residual_seconds = round(max(
      0, frontier_total_seconds - frontier_classified
    ), 2L),
    total_seconds = round(frontier_total_seconds, 2L)
  )

  pu_display <- round(pu_repair_timing, 2L)
  if (sum(pu_repair_timing) > pu_repair_total_seconds + 1e-8) {
    stop("PU repair profile exceeds the measured PU-repair total.")
  }
  pu_residual <- max(0, pu_repair_total_seconds - sum(pu_repair_timing))
  runtime_log_event(
    "pruning_pu_repair_profile",
    stage = as.integer(stage_index),
    iters = as.integer(completed_pruning_iterations),
    species_pu_repairs = as.numeric(pu_repair_counts[["species_pu_repairs"]]),
    graph_nodes_examined = as.numeric(
      pu_repair_counts[["graph_nodes_examined"]]
    ),
    adjacency_entries_examined = as.numeric(
      pu_repair_counts[["adjacency_entries_examined"]]
    ),
    dead_patch_nodes = as.numeric(pu_repair_counts[["dead_patch_nodes"]]),
    surviving_components = as.numeric(
      pu_repair_counts[["surviving_components"]]
    ),
    split_pus = as.numeric(pu_repair_counts[["split_pus"]]),
    component_dropped_patches = as.numeric(
      pu_repair_counts[["component_dropped_patches"]]
    ),
    threshold_candidate_pus = as.numeric(
      pu_repair_counts[["threshold_candidate_pus"]]
    ),
    threshold_dropped_pus = as.numeric(
      pu_repair_counts[["threshold_dropped_pus"]]
    ),
    compiled_species_batches = as.numeric(
      pu_repair_counts[["compiled_species_batches"]]
    ),
    compiled_pu_workloads = as.numeric(
      pu_repair_counts[["compiled_pu_workloads"]]
    ),
    graph_actions_staged = as.numeric(pu_repair_graph_actions),
    setup_seconds = pu_display[["setup_seconds"]],
    component_rebuild_seconds = pu_display[["component_rebuild_seconds"]],
    transaction_seconds = pu_display[["transaction_seconds"]],
    pu_threshold_seconds = pu_display[["pu_threshold_seconds"]],
    residual_seconds = round(pu_residual, 2L),
    total_seconds = round(pu_repair_total_seconds, 2L)
  )
  invisible(NULL)
}


# Clear ordinary removals and threshold-dropped patches from each affected
# species vector in one stage-boundary assignment. Compact patch indexes remain
# unchanged here and are synchronized by the subsequent fragmentation splice.
flush_pruning_pending_cells <- function(
  patch_id_by_species_env,
  patch_cell_index_by_species_env,
  removed_cells,
  changed_species,
  threshold_dropped_patches,
  patch_table,
  alive_species_count_by_cell
) {
  removed_cells <- sort(unique(as.integer(removed_cells)))
  changed_species <- unique(as.character(changed_species))
  if (anyNA(changed_species) || any(!nzchar(changed_species))) {
    stop("Deferred pruning clears contain invalid changed-species names.")
  }

  if (!data.table::is.data.table(threshold_dropped_patches)) {
    threshold_dropped_patches <- data.table::as.data.table(threshold_dropped_patches)
  }
  required_drop_columns <- c("species", "patch_id")
  if (!all(required_drop_columns %in% names(threshold_dropped_patches))) {
    stop("Deferred threshold drops must contain species and patch_id columns.")
  }
  threshold_dropped_patches <- unique(
    threshold_dropped_patches[, .(
      species = as.character(species),
      patch_id = as.integer(patch_id)
    )],
    by = c("species", "patch_id")
  )
  if (nrow(threshold_dropped_patches) &&
      (anyNA(threshold_dropped_patches$species) ||
       any(!nzchar(threshold_dropped_patches$species)) ||
       anyNA(threshold_dropped_patches$patch_id) ||
       any(threshold_dropped_patches$patch_id < 1L))) {
    stop("Deferred threshold drops contain invalid species or patch IDs.")
  }

  if (nrow(threshold_dropped_patches)) {
    still_live <- patch_table[
      threshold_dropped_patches,
      on = .(species, patch_id),
      nomatch = 0L
    ]
    if (nrow(still_live)) {
      stop("Deferred threshold-dropped patches remain in the final pruning patch table.")
    }
  }

  affected_species <- sort(unique(c(
    changed_species,
    threshold_dropped_patches$species
  )))

  if (!length(removed_cells) && !nrow(threshold_dropped_patches)) {
    return(0L)
  }
  if (!length(affected_species)) {
    stop("Deferred pruning clears contain pending cells but no affected species.")
  }
  if (anyNA(removed_cells) || any(removed_cells < 1L) ||
      any(removed_cells > length(alive_species_count_by_cell))) {
    stop("Deferred pruning clears contain invalid raster-cell IDs.")
  }
  if (length(removed_cells) && any(alive_species_count_by_cell[removed_cells] != 0L)) {
    stop("Deferred pruning clears include cells that remain globally alive.")
  }

  flushed <- 0L
  for (species_name in affected_species) {
    if (!exists(
      species_name,
      envir = patch_id_by_species_env,
      inherits = FALSE
    )) {
      stop("Missing species vector during deferred pruning clear: ", species_name)
    }
    patch_ids <- get(
      species_name,
      envir = patch_id_by_species_env,
      inherits = FALSE
    )
    if (length(patch_ids) != length(alive_species_count_by_cell)) {
      stop("Species vector has the wrong cell count during deferred clear: ", species_name)
    }

    ordinary_cells <- integer()
    if (length(removed_cells)) {
      ordinary_cells <- removed_cells[!is.na(patch_ids[removed_cells])]
    }

    species_drops <- threshold_dropped_patches[species == species_name]
    threshold_cells <- integer()
    if (nrow(species_drops)) {
      if (!exists(
        species_name,
        envir = patch_cell_index_by_species_env,
        inherits = FALSE
      )) {
        stop("Missing compact patch-cell index during deferred clear: ", species_name)
      }
      patch_index <- get(
        species_name,
        envir = patch_cell_index_by_species_env,
        inherits = FALSE
      )
      indexed <- get_patch_cells_batch(
        patch_index = patch_index,
        patch_ids = species_drops$patch_id
      )
      if (anyDuplicated(indexed$cell)) {
        stop("Deferred threshold-drop index contains duplicate cells for species: ", species_name)
      }
      represented <- sort(unique(as.integer(indexed$patch_id)))
      missing_patch_ids <- setdiff(sort(species_drops$patch_id), represented)
      if (length(missing_patch_ids)) {
        stop(
          "Deferred threshold-dropped patches are missing from the compact index for species ",
          species_name, ": ", paste(utils::head(missing_patch_ids, 5L), collapse = ", ")
        )
      }
      threshold_cells <- validate_priority_cell_ids(
        indexed$cell,
        n_cells = length(patch_ids),
        label = paste0("Deferred threshold-drop cells for species ", species_name)
      )
      expected_by_cell <- indexed$patch_id[match(threshold_cells, indexed$cell)]
      observed_by_cell <- patch_ids[threshold_cells]
      mismatch <- is.na(observed_by_cell) | observed_by_cell != expected_by_cell
      if (any(mismatch)) {
        first <- which(mismatch)[1L]
        stop(
          "Deferred threshold-clear dense/index mismatch for species ", species_name,
          ", patch ", expected_by_cell[first], ", cell ", threshold_cells[first], "."
        )
      }
    }

    cells_to_clear <- sort(unique(c(ordinary_cells, threshold_cells)))
    if (length(cells_to_clear)) {
      patch_ids[cells_to_clear] <- NA_integer_
      assign(species_name, patch_ids, envir = patch_id_by_species_env)
      flushed <- flushed + 1L
    }
    if (length(cells_to_clear) && any(!is.na(patch_ids[cells_to_clear]))) {
      stop("Deferred pruning clear failed for species: ", species_name)
    }
  }

  as.integer(flushed)
}


# Run one full pruning stage.
#
# This function is the stage-level driver for the pruning portion of the
# spatial prioritization pipeline.
#
# It repeatedly calls run_pruning_iteration(), updating the live state after
# each call. It does not write rasters or masks. Instead, it returns the
# per-iteration removed-cell vectors so the outer pipeline can record removal
# order in memory and write one final removal-order raster later.
#
run_pruning_stage <- function(
  stage_index,                      # global stage index for this pruning stage
  pruning_iterations_per_stage,     # maximum pruning iterations to run in this stage
  cells_to_remove_per_iteration,    # number of frontier cells to remove per iteration
  patch_table,                      # current patch table
  pu_graphs_by_key,                 # current PU graph objects keyed by "species|pu_id"
  alive_species_count_by_cell,      # current number of surviving species in each cell
  patch_id_by_species_env,          # live environment: species -> cell-to-patch vector
  patch_cell_index_by_species_env,  # live environment: species -> patch-to-cells index
  cell_species_index = NULL,        # immutable cell -> possible species CSR index
  cell_area_by_cell,                # area of each raster cell in km^2
  frontier_state,                   # current incremental frontier state
  rook_neighbor_index,              # compact rook-neighbor lookup
  species_params,                   # species-level thresholds and persistence parameters
  ecology_log_every_iterations = NULL,
  ecology_previous_snapshot = NULL,
  initial_alive_cell_count = NULL,
  current_alive_cell_count = NULL
) {
  # -------------------------------------------------------------------
  # 1) Validate scalar stage-level arguments
  # -------------------------------------------------------------------

  if (!is.numeric(stage_index) || length(stage_index) != 1L || is.na(stage_index)) {
    stop("stage_index must be one non-missing numeric value.")
  }

  if (!is.numeric(pruning_iterations_per_stage) ||
      length(pruning_iterations_per_stage) != 1L ||
      is.na(pruning_iterations_per_stage) ||
      pruning_iterations_per_stage < 1) {
    stop("pruning_iterations_per_stage must be a single integer >= 1.")
  }

  if (!is.numeric(cells_to_remove_per_iteration) ||
      length(cells_to_remove_per_iteration) != 1L ||
      is.na(cells_to_remove_per_iteration) ||
      cells_to_remove_per_iteration < 1) {
    stop("cells_to_remove_per_iteration must be a single integer >= 1.")
  }

  pruning_iterations_per_stage <- as.integer(pruning_iterations_per_stage)

  cells_to_remove_per_iteration <- as.integer(cells_to_remove_per_iteration)


  # -------------------------------------------------------------------
  # 2) Prepare pruning-stage state shared across iterations
  # -------------------------------------------------------------------

  # Build the species-name vector once for this pruning stage.
  species_names <- sort(ls(envir = patch_id_by_species_env, all.names = TRUE))

  # Start this pruning stage with no patch-score cache.
  patch_score_cache <- NULL

  # Force the first pruning iteration to score all species currently in patch_table.
  dirty_species_for_scoring <- sort(unique(as.character(patch_table$species)))


  # -------------------------------------------------------------------
  # 3) Allocate per-iteration accumulators
  # -------------------------------------------------------------------

  # Allocate one removed-cell vector slot per possible pruning iteration.
  removed_cells_each_iteration <- vector("list", pruning_iterations_per_stage)

  # Allocate one changed-patch table slot per possible pruning iteration.
  changed_patches_each_iteration <- vector("list", pruning_iterations_per_stage)

  # Retain only compact species/patch IDs whose dense-cell clearing is deferred.
  threshold_drops_each_iteration <- vector("list", pruning_iterations_per_stage)

  # Count how many pruning iterations actually completed.
  completed_pruning_iterations <- 0L

  # Track whether the stage ended because the frontier was exhausted.
  frontier_exhausted <- FALSE

  # Accumulate compact timing totals without retaining per-iteration records.
  pruning_timing_seconds <- empty_pruning_timing()
  pruning_work_counts <- empty_pruning_work_counts()
  frontier_scoring_counts <- empty_frontier_scoring_counts()
  frontier_scoring_timing <- empty_frontier_scoring_timing()
  pu_repair_profile_timing <- empty_pruning_pu_repair_timing()
  pu_repair_profile_counts <- empty_pruning_pu_repair_counts()
  pu_repair_graph_actions <- 0


  # -------------------------------------------------------------------
  # 4) Run repeated pruning iterations
  # -------------------------------------------------------------------

  # Loop over the maximum allowed pruning iterations for this stage.
  for (iteration_in_stage in seq_len(pruning_iterations_per_stage)) {
    # Run one pruning iteration using the current live state.
    pruning_result <- run_pruning_iteration(
      prune_iteration                 = iteration_in_stage,
      cells_to_remove_per_step        = cells_to_remove_per_iteration,
      patch_table                     = patch_table,
      pu_graphs_by_key                = pu_graphs_by_key,
      alive_species_count_by_cell     = alive_species_count_by_cell,
      cell_area_by_cell               = cell_area_by_cell,
      species_params                  = species_params,
      patch_id_by_species_env         = patch_id_by_species_env,
      patch_cell_index_by_species_env = patch_cell_index_by_species_env,
      cell_species_index              = cell_species_index,
      frontier_state                  = frontier_state,
      rook_neighbor_index             = rook_neighbor_index,
      species_names                   = species_names,
      patch_score_cache               = patch_score_cache,
      dirty_species_for_scoring       = dirty_species_for_scoring,
      stage_index                     = stage_index,
      ecology_log_every_iterations    = ecology_log_every_iterations,
      ecology_previous_snapshot       = ecology_previous_snapshot,
      initial_alive_cell_count        = initial_alive_cell_count,
      current_alive_cell_count        = current_alive_cell_count
    )

    # Replace the stage patch table with the updated patch table.
    patch_table <- pruning_result$patch_table

    # Replace the stage PU graph list with the updated graph list.
    pu_graphs_by_key <- pruning_result$pu_graphs

    # Replace the stage alive-species count vector with the updated vector.
    alive_species_count_by_cell <- pruning_result$alive_species_count_by_cell

    # Replace the stage frontier state with the updated frontier state.
    frontier_state <- pruning_result$frontier_state

    # Carry forward the updated patch-score cache.
    patch_score_cache <- pruning_result$patch_score_cache

    # Carry forward the dirty-species set for the next pruning iteration.
    dirty_species_for_scoring <- pruning_result$dirty_species_for_scoring

    # Carry the tiny logging-only state across scheduled diagnostics.
    ecology_previous_snapshot <- pruning_result$ecology_previous_snapshot
    current_alive_cell_count <- pruning_result$current_alive_cell_count

    iteration_timing <- pruning_result$pruning_timing_seconds
    if (!is.numeric(iteration_timing) ||
        !identical(names(iteration_timing), names(pruning_timing_seconds)) ||
        any(!is.finite(iteration_timing)) ||
        any(iteration_timing < 0)) {
      stop("A pruning iteration returned invalid subphase timings.")
    }
    pruning_timing_seconds <- pruning_timing_seconds + iteration_timing

    iteration_counts <- pruning_result$pruning_work_counts
    if (!is.numeric(iteration_counts) ||
        !identical(names(iteration_counts), names(pruning_work_counts)) ||
        any(!is.finite(iteration_counts)) || any(iteration_counts < 0L) ||
        any(iteration_counts != floor(iteration_counts))) {
      stop("A pruning iteration returned invalid workload counts.")
    }
    pruning_work_counts <- pruning_work_counts + iteration_counts

    iteration_frontier_counts <- pruning_result$frontier_scoring_counts
    iteration_frontier_timing <- pruning_result$frontier_scoring_timing
    iteration_pu_timing <- pruning_result$pu_repair_profile_timing
    iteration_pu_counts <- pruning_result$pu_repair_profile_counts
    if (!is.numeric(iteration_frontier_counts) ||
        !identical(names(iteration_frontier_counts), names(frontier_scoring_counts)) ||
        any(!is.finite(iteration_frontier_counts)) ||
        any(iteration_frontier_counts < 0) ||
        !is.numeric(iteration_frontier_timing) ||
        !identical(names(iteration_frontier_timing), names(frontier_scoring_timing)) ||
        any(!is.finite(iteration_frontier_timing)) ||
        any(iteration_frontier_timing < 0) ||
        !is.numeric(iteration_pu_timing) ||
        !identical(names(iteration_pu_timing), names(pu_repair_profile_timing)) ||
        any(!is.finite(iteration_pu_timing)) || any(iteration_pu_timing < 0) ||
        !is.numeric(iteration_pu_counts) ||
        !identical(names(iteration_pu_counts), names(pu_repair_profile_counts)) ||
        any(!is.finite(iteration_pu_counts)) || any(iteration_pu_counts < 0)) {
      stop("A pruning iteration returned invalid profiling diagnostics.")
    }
    frontier_scoring_counts <-
      frontier_scoring_counts + iteration_frontier_counts
    frontier_scoring_timing <-
      frontier_scoring_timing + iteration_frontier_timing
    pu_repair_profile_timing <-
      pu_repair_profile_timing + iteration_pu_timing
    pu_repair_profile_counts <-
      pu_repair_profile_counts + iteration_pu_counts
    pu_repair_graph_actions <-
      pu_repair_graph_actions + pruning_result$pu_repair_graph_actions

    # Extract cells that became newly empty during this pruning iteration.
    removed_cells_this_iteration <- pruning_result$newly_empty_cells

    # Extract patches changed or removed during this pruning iteration.
    changed_patches_this_iteration <- pruning_result$changed_patches

    threshold_drops_this_iteration <- pruning_result$threshold_dropped_patches
    if (!is.data.frame(threshold_drops_this_iteration) ||
        !all(c("species", "patch_id") %in% names(threshold_drops_this_iteration))) {
      stop("A pruning iteration returned invalid deferred threshold-drop bookkeeping.")
    }

    # Store this iteration's newly empty cells.
    removed_cells_each_iteration[[iteration_in_stage]] <- removed_cells_this_iteration

    # Store this iteration's changed patches.
    changed_patches_each_iteration[[iteration_in_stage]] <- changed_patches_this_iteration

    threshold_drops_each_iteration[[iteration_in_stage]] <-
      data.table::as.data.table(threshold_drops_this_iteration)[, .(
        species = as.character(species),
        patch_id = as.integer(patch_id)
      )]

    # Mark this pruning iteration as completed.
    completed_pruning_iterations <- as.integer(iteration_in_stage)

    # If this iteration removed no cells, the frontier is exhausted.
    if (length(removed_cells_this_iteration) == 0L) {
      # Record that the pruning stage ended by frontier exhaustion.
      frontier_exhausted <- TRUE

      # Stop the pruning-stage loop early.
      break
    }
  }


  # -------------------------------------------------------------------
  # 5) Trim per-iteration accumulators to completed iterations
  # -------------------------------------------------------------------

  # If at least one pruning iteration completed, keep only completed slots.
  if (completed_pruning_iterations > 0L) {
    # Keep removed-cell vectors for completed iterations only.
    removed_cells_each_iteration <-
      removed_cells_each_iteration[seq_len(completed_pruning_iterations)]

    # Keep changed-patch tables for completed iterations only.
    changed_patches_each_iteration <-
      changed_patches_each_iteration[seq_len(completed_pruning_iterations)]

    threshold_drops_each_iteration <-
      threshold_drops_each_iteration[seq_len(completed_pruning_iterations)]
  } else {
    # Use an empty list if no pruning iterations completed.
    removed_cells_each_iteration <- list()

    # Use an empty list if no pruning iterations completed.
    changed_patches_each_iteration <- list()

    threshold_drops_each_iteration <- list()
  }


  # -------------------------------------------------------------------
  # 6) Consolidate removed cells across the full pruning stage
  # -------------------------------------------------------------------

  # If no iterations completed, no cells were removed.
  if (!length(removed_cells_each_iteration)) {
    # Store an empty integer vector for stage-level removed cells.
    removed_cells_in_stage <- integer(0L)
  } else {
    # Combine, deduplicate, sort, and coerce removed cells to integer IDs.
    removed_cells_in_stage <- sort(
      unique(
        as.integer(
          unlist(
            removed_cells_each_iteration,
            use.names = FALSE
          )
        )
      )
    )
  }


  # -------------------------------------------------------------------
  # 7) Consolidate changed patches across the full pruning stage
  # -------------------------------------------------------------------

  # Keep only changed-patch objects that are data frames with at least one row.
  non_empty_changed_patch_tables <- Filter(
    f = function(x) {
      is.data.frame(x) && nrow(x) > 0L
    },
    x = changed_patches_each_iteration
  )

  # If no patches changed, return an empty changed-patch table.
  if (!length(non_empty_changed_patch_tables)) {
    # Create an empty two-column changed-patch table.
    changed_patches_in_stage <- data.table::data.table(
      species = character(),
      patch_id = integer()
    )
  } else {
    # Combine all changed-patch tables and deduplicate by species and patch ID.
    changed_patches_in_stage <- unique(
      data.table::rbindlist(
        non_empty_changed_patch_tables,
        use.names = TRUE,
        fill = TRUE
      ),
      by = c("species", "patch_id")
    )
  }

  non_empty_threshold_drop_tables <- Filter(
    f = function(x) is.data.frame(x) && nrow(x) > 0L,
    x = threshold_drops_each_iteration
  )
  if (!length(non_empty_threshold_drop_tables)) {
    threshold_dropped_patches_in_stage <- data.table::data.table(
      species = character(),
      patch_id = integer()
    )
  } else {
    threshold_dropped_patches_in_stage <- unique(
      data.table::rbindlist(
        non_empty_threshold_drop_tables,
        use.names = TRUE,
        fill = FALSE
      ),
      by = c("species", "patch_id")
    )
  }


  # -------------------------------------------------------------------
  # 8) Synchronize dense species vectors before spatial repair
  # -------------------------------------------------------------------

  deferred_clear_started <- proc.time()[["elapsed"]]
  species_vectors_flushed <- flush_pruning_pending_cells(
    patch_id_by_species_env = patch_id_by_species_env,
    patch_cell_index_by_species_env = patch_cell_index_by_species_env,
    removed_cells = removed_cells_in_stage,
    changed_species = changed_patches_in_stage$species,
    threshold_dropped_patches = threshold_dropped_patches_in_stage,
    patch_table = patch_table,
    alive_species_count_by_cell = alive_species_count_by_cell
  )
  pruning_timing_seconds[["deferred_cell_clear_seconds"]] <-
    pruning_timing_seconds[["deferred_cell_clear_seconds"]] +
    proc.time()[["elapsed"]] - deferred_clear_started


  # -------------------------------------------------------------------
  # 9) Compute final stage summaries for logging
  # -------------------------------------------------------------------

  # Count remaining patch-table rows.
  remaining_patch_count <- nrow(patch_table)

  # Count remaining unique species-PU combinations.
  remaining_pu_count <- data.table::uniqueN(
    patch_table,
    by = c("species", "pu_id")
  )

  # Count cells that still contain at least one surviving species.
  remaining_alive_cell_count <- if (is.null(current_alive_cell_count)) {
    sum(alive_species_count_by_cell > 0L)
  } else {
    as.integer(current_alive_cell_count)
  }

  # Count unique cells removed during this pruning stage.
  removed_cell_count <- length(removed_cells_in_stage)

  # Count unique patches changed during this pruning stage.
  changed_patch_count <- nrow(changed_patches_in_stage)


  # -------------------------------------------------------------------
  # 10) Print exactly one pruning-stage summary log message
  # -------------------------------------------------------------------

  # Print the concise pruning-stage summary.
  log_pruning_stage(
    stage_index                  = stage_index,
    completed_pruning_iterations = completed_pruning_iterations,
    removed_cell_count           = removed_cell_count,
    changed_patch_count          = changed_patch_count,
    frontier_exhausted           = frontier_exhausted,
    remaining_patch_count        = remaining_patch_count,
    remaining_pu_count           = remaining_pu_count,
    remaining_alive_cell_count   = remaining_alive_cell_count
  )

  log_pruning_timing(
    stage_index = stage_index,
    completed_pruning_iterations = completed_pruning_iterations,
    timing_seconds = pruning_timing_seconds,
    work_counts = pruning_work_counts,
    species_vectors_flushed = species_vectors_flushed
  )
  log_pruning_profiles(
    stage_index = stage_index,
    completed_pruning_iterations = completed_pruning_iterations,
    frontier_counts = frontier_scoring_counts,
    frontier_timing = frontier_scoring_timing,
    frontier_total_seconds =
      pruning_timing_seconds[["frontier_scoring_seconds"]],
    pu_repair_counts = pu_repair_profile_counts,
    pu_repair_timing = pu_repair_profile_timing,
    pu_repair_total_seconds = pruning_timing_seconds[["pu_repair_seconds"]],
    pu_repair_graph_actions = pu_repair_graph_actions
  )


  # -------------------------------------------------------------------
  # 11) Return updated state and pruning-stage bookkeeping
  # -------------------------------------------------------------------

  # Return the updated live state and compact stage bookkeeping.
  list(
    patch_table                   = patch_table,
    pu_graphs                     = pu_graphs_by_key,
    alive_species_count_by_cell   = alive_species_count_by_cell,
    frontier_state                = frontier_state,
    removed_cells_each_iteration  = removed_cells_each_iteration,
    removed_cells_in_stage        = removed_cells_in_stage,
    changed_patches_in_stage      = changed_patches_in_stage,
    completed_pruning_iterations  = completed_pruning_iterations,
    frontier_exhausted            = frontier_exhausted,
    ecology_previous_snapshot     = ecology_previous_snapshot,
    current_alive_cell_count      = remaining_alive_cell_count,
    timing_seconds                = pruning_timing_seconds
  )
}
