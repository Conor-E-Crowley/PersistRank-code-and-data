# Stage-wise reverse-removal pipeline orchestration.
#
# Mutable patch, graph, cell, frontier, and checkpoint state advances in
# canonical stage order. Scientific work remains in the pruning,
# fragmentation, and distance modules; this function coordinates those phases
# and writes validated stage/final outputs.


run_priority_pipeline <- function(
  cells_to_remove_per_iteration,    # frontier cells removed per pruning iteration
  pruning_iterations_per_stage,     # maximum pruning iterations per global stage
  patch_table,                      # current patch table
  pu_graphs_by_key,                 # current PU graphs keyed by "species|pu_id"
  alive_species_count_by_cell,      # current surviving-species count by raster cell
  patch_id_by_species_env,          # live env: species -> cell-to-patch vector
  patch_cell_index_by_species_env,  # live env: species -> patch-to-cells index
  cell_area_by_cell,                # area of each raster cell in km^2
  rook_neighbor_pairs,              # rook-neighbor cell pairs for frontier tracking
  species_params,                   # species thresholds, persistence params, dispersal distance
  template_raster,                  # shared raster geometry for the distance stage
  mask_template_raster,             # single-layer grid template for removal_order.tif
  output_dir,                       # root output directory
  max_stages = Inf,                 # maximum number of full global stages
  resume_checkpoint = NULL,         # optional completed-stage checkpoint
  checkpoint_every_stages = NULL,   # NULL disables stage-boundary checkpoints
  checkpoint_keep = 2L,             # number of recent checkpoints to retain
  checkpoint_metadata = NULL,       # metadata written to future checkpoints
  ecology_log_every_iterations = NULL # NULL disables ecological diagnostics
) {
  # -------------------------------------------------------------------
  # 1. Validate scalar control inputs
  # -------------------------------------------------------------------

  cells_to_remove_per_iteration <- validate_priority_count(
    cells_to_remove_per_iteration,
    "cells_to_remove_per_iteration"
  )
  pruning_iterations_per_stage <- validate_priority_count(
    pruning_iterations_per_stage,
    "pruning_iterations_per_stage"
  )
  max_stages <- validate_max_stages(max_stages)

  if (!is.null(ecology_log_every_iterations)) {
    ecology_log_every_iterations <- validate_priority_count(
      ecology_log_every_iterations,
      "ecology_log_every_iterations"
    )
  }

  checkpoint_enabled <- !is.null(checkpoint_every_stages)

  if (isTRUE(checkpoint_enabled)) {
    if (is.null(checkpoint_keep)) {
      checkpoint_keep <- 2L
    }
    checkpoint_every_stages <- validate_priority_count(
      checkpoint_every_stages,
      "checkpoint_every_stages"
    )
    checkpoint_keep <- validate_priority_count(checkpoint_keep, "checkpoint_keep")

    if (is.null(checkpoint_metadata) || !is.list(checkpoint_metadata)) {
      stop("checkpoint_metadata must be supplied when checkpointing is enabled.")
    }
  } else {
    checkpoint_every_stages <- NULL
  }

  resume_requested <- !is.null(resume_checkpoint)

  if (isTRUE(resume_requested) && !is.list(resume_checkpoint)) {
    stop("resume_checkpoint must be a priority checkpoint list.")
  }

  if (isTRUE(resume_requested) && is.null(checkpoint_metadata)) {
    checkpoint_metadata <- resume_checkpoint$metadata
  }

  # -------------------------------------------------------------------
  # 2. Standardize table inputs
  # -------------------------------------------------------------------

  patch_table <- data.table::as.data.table(patch_table)

  species_params <- data.table::as.data.table(species_params)


  # -------------------------------------------------------------------
  # 3. Validate required columns and template objects
  # -------------------------------------------------------------------

  validate_priority_species_parameters(
    species_params,
    contract = checkpoint_metadata$analysis_contract %||% canonical_analysis_contract()
  )

  if (!is.character(output_dir) || length(output_dir) != 1L || !nzchar(output_dir)) {
    stop("output_dir must be one non-empty character path.")
  }

  if (is.null(mask_template_raster)) {
    stop("mask_template_raster must be supplied for writing removal_order.tif.")
  }

  if (is.null(template_raster)) {
    stop("template_raster must be supplied for the distance stage.")
  }


  # -------------------------------------------------------------------
  # 4. Create the minimal production output directories
  # -------------------------------------------------------------------

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  patch_lookup_output_dir <- file.path(output_dir, "patch_lookup_tables")

  dir.create(patch_lookup_output_dir, recursive = TRUE, showWarnings = FALSE)

  removal_order_output_path <- file.path(output_dir, "removal_order.tif")

  rankmap_output_path <- file.path(output_dir, "rankmap.tif")
  
  removal_events_output_path <- file.path(output_dir, "removal_events.csv")


  # -------------------------------------------------------------------
  # 5. Initialize or restore global stage-loop state
  # -------------------------------------------------------------------

  if (isTRUE(resume_requested)) {
    validate_priority_checkpoint_schema(
      checkpoint = resume_checkpoint,
      expected_metadata = checkpoint_metadata,
      n_cells = length(cell_area_by_cell)
    )

    patch_table <- data.table::as.data.table(resume_checkpoint$patch_table)
    pu_graphs_by_key <- resume_checkpoint$pu_graphs_by_key
    alive_species_count_by_cell <- resume_checkpoint$alive_species_count_by_cell
    patch_id_by_species_env <- resume_checkpoint$patch_id_by_species_env
    removal_order_by_cell <- resume_checkpoint$removal_order_by_cell
    removal_event_rows <- resume_checkpoint$removal_event_rows
    removal_step_counter <- as.integer(resume_checkpoint$removal_step_counter)
    completed_stages <- as.integer(resume_checkpoint$completed_stages)
    frontier_exhausted <- isTRUE(resume_checkpoint$frontier_exhausted)
    stage_index <- as.integer(resume_checkpoint$stage_index)

    validate_no_stage_outputs_after_checkpoint(
      patch_lookup_output_dir = patch_lookup_output_dir,
      checkpoint_stage = completed_stages
    )

    patch_cell_index_by_species_env <-
      rebuild_patch_cell_index_env_from_patch_ids(patch_id_by_species_env)

    gc(FALSE)
  } else {
    stage_index <- 0L

    # Allocate the removal-order vector; NA means not removed or not initially alive.
    removal_order_by_cell <- rep(NA_integer_, length(alive_species_count_by_cell))

    removal_step_counter <- 1L

    # Store one small metadata row per non-empty removal event.
    removal_event_rows <- vector("list", 0L)

    # Count completed full global stages.
    completed_stages <- 0L

    # Track whether pruning has exhausted the frontier.
    frontier_exhausted <- FALSE
  }

  # Use hashed graph lookup and mutation throughout the iterative pipeline.
  # Bundles and checkpoints are validated in their established list form first.
  pu_graphs_by_key <- new_priority_graph_store(pu_graphs_by_key)

  cell_state <- validate_priority_cell_state(
    alive_species_count_by_cell = alive_species_count_by_cell,
    cell_area_by_cell = cell_area_by_cell,
    removal_order_by_cell = removal_order_by_cell,
    retained_species_count = nrow(species_params)
  )
  alive_species_count_by_cell <- cell_state$alive_species_count_by_cell
  removal_order_by_cell <- cell_state$removal_order_by_cell
  initial_alive_by_cell <- cell_state$initial_alive_by_cell
  initial_alive_cell_count <- cell_state$initial_alive_cell_count
  current_alive_cell_count <- cell_state$current_alive_cell_count
  initial_alive_area_km2 <- cell_state$initial_alive_area_km2
  rm(cell_state)
  n_cells <- length(alive_species_count_by_cell)

  # Build one immutable sparse membership superset for frontier scoring. Species
  # presence can only disappear during prioritization, so current dense patch
  # vectors remain the authority for every candidate contribution.
  cell_species_index <- build_cell_species_index(
    patch_cell_index_by_species_env = patch_cell_index_by_species_env,
    species_names = sort(ls(envir = patch_id_by_species_env, all.names = TRUE)),
    n_cells = n_cells,
    alive_species_count_by_cell = alive_species_count_by_cell
  )
  # Keep logging-only state outside checkpoints and scientific outputs. A
  # resumed process intentionally begins with a fresh diagnostic baseline.
  ecology_previous_snapshot <- NULL

  log_priority_start(
    cells_to_remove_per_iteration = cells_to_remove_per_iteration,
    pruning_iterations_per_stage = pruning_iterations_per_stage,
    initial_alive_cell_count = initial_alive_cell_count,
    species_count = length(unique(as.character(species_params$species)))
  )
  log_priority_runtime_strategy()

  if (isTRUE(resume_requested)) {
    log_priority_resume(
      checkpoint_path = attr(resume_checkpoint, "checkpoint_path", exact = TRUE),
      stage_index = stage_index,
      removal_step_counter = removal_step_counter,
      alive_cell_count = current_alive_cell_count
    )
  }


  # -------------------------------------------------------------------
  # 6. Define local helper to record removal order in memory
  # -------------------------------------------------------------------
  
  record_removed_cells <- function(
    removed_cells,
    stage_index,
    event_type,
    pruning_iteration = NA_integer_
  ) {
    # Return immediately if this event removed no cells.
    if (!length(removed_cells)) {
      return(invisible(0L))
    }
  
    removed_cells <- validate_priority_cell_ids(
      removed_cells,
      n_cells = n_cells,
      label = paste0("Stage 6 ", event_type, " removal cells")
    )
  
    # Return if no valid cell IDs remain.
    if (!length(removed_cells)) {
      return(invisible(0L))
    }
  
    # Keep only cells that:
    #   1. were alive at the start of the run, and
    #   2. have not already been assigned a removal order.
    eligible_cells <- removed_cells[
      initial_alive_by_cell[removed_cells] &
        is.na(removal_order_by_cell[removed_cells])
    ]
  
    # Return if this event contains no newly recordable cells.
    if (!length(eligible_cells)) {
      return(invisible(0L))
    }
  
    this_step <- as.integer(removal_step_counter)
  
    # Assign this removal step to all eligible cells in the removal-order raster vector.
    removal_order_by_cell[eligible_cells] <<- this_step
  
    # Compute how many cells were actually recorded in this event.
    cells_removed_this_event <- length(eligible_cells)
  
    # Compute the area removed in this event.
    area_removed_this_event_km2 <- sum(
      cell_area_by_cell[eligible_cells],
      na.rm = TRUE
    )
  
    # Append one metadata row for this removal event.
    removal_event_rows[[length(removal_event_rows) + 1L]] <<-
      data.table::data.table(
        removal_step = this_step,
        stage = as.integer(stage_index),
        event_type = as.character(event_type),
        pruning_iteration = as.integer(pruning_iteration),
        cells_removed = as.integer(cells_removed_this_event),
        area_removed_km2 = as.numeric(area_removed_this_event_km2)
      )
  
    # Advance the removal-event counter once per non-empty recorded event.
    removal_step_counter <<- removal_step_counter + 1L
  
    # Return the number of cells recorded in this event.
    invisible(cells_removed_this_event)
  }
  
  
  
  record_final_retained_cells <- function(stage_index) {
    # Identify initially alive cells that never received a removal-order value.
    final_retained_cells <- which(
      initial_alive_by_cell &
        is.na(removal_order_by_cell)
    )
  
    # Assign one final removal-order step to all retained-to-end cells. If all
    # cells have already been removed, keep a zero-cell terminal event so the
    # event table still has an explicit endpoint.
    this_step <- as.integer(removal_step_counter)
  
    if (length(final_retained_cells)) {
      removal_order_by_cell[final_retained_cells] <<- this_step
    }
  
    # Compute retained-to-end area.
    final_retained_area_km2 <- sum(
      cell_area_by_cell[final_retained_cells],
      na.rm = TRUE
    )
  
    # Add a final event row so removal_order.tif and removal_events.csv stay consistent.
    removal_event_rows[[length(removal_event_rows) + 1L]] <<-
      data.table::data.table(
        removal_step = this_step,
        stage = as.integer(stage_index),
        event_type = "final_retained",
        pruning_iteration = NA_integer_,
        cells_removed = as.integer(length(final_retained_cells)),
        area_removed_km2 = as.numeric(final_retained_area_km2)
      )
  
    # Advance the event counter.
    removal_step_counter <<- removal_step_counter + 1L
  
    invisible(length(final_retained_cells))
  }


  # -------------------------------------------------------------------
  # 7. Build frontier-tracking objects once
  # -------------------------------------------------------------------

  # Build the compact rook-neighbor index once from the input neighbor pairs.
  rook_neighbor_index <- build_rook_neighbor_index(
    rook_neighbor_pairs = rook_neighbor_pairs,
    n_cells = n_cells,
    force_symmetric = TRUE
  )

  frontier_state <- initialize_frontier_state(
    alive_species_count_by_cell = alive_species_count_by_cell,
    rook_neighbor_index = rook_neighbor_index
  )


  # -------------------------------------------------------------------
  # 8. Define optional stage-checkpoint writer
  # -------------------------------------------------------------------

  checkpoint_output_dir <- priority_checkpoint_dir(output_dir)
  last_checkpoint_stage <- if (isTRUE(resume_requested)) {
    as.integer(completed_stages)
  } else {
    NA_integer_
  }

  save_stage_checkpoint <- function(force = FALSE) {
    if (!isTRUE(checkpoint_enabled) || completed_stages < 1L) {
      return(invisible(NULL))
    }

    if (!isTRUE(force) && completed_stages %% checkpoint_every_stages != 0L) {
      return(invisible(NULL))
    }

    if (!is.na(last_checkpoint_stage) &&
        as.integer(last_checkpoint_stage) == as.integer(completed_stages)) {
      return(invisible(NULL))
    }

    checkpoint_started <- proc.time()[["elapsed"]]
    checkpoint <- make_priority_checkpoint(
      metadata = checkpoint_metadata,
      stage_index = completed_stages,
      completed_stages = completed_stages,
      frontier_exhausted = frontier_exhausted,
      patch_table = patch_table,
      pu_graphs_by_key = pu_graphs_by_key,
      alive_species_count_by_cell = alive_species_count_by_cell,
      patch_id_by_species_env = patch_id_by_species_env,
      removal_order_by_cell = removal_order_by_cell,
      removal_event_rows = removal_event_rows,
      removal_step_counter = removal_step_counter
    )

    checkpoint_path <- write_priority_checkpoint(
      checkpoint = checkpoint,
      checkpoint_dir = checkpoint_output_dir,
      checkpoint_keep = checkpoint_keep
    )

    last_checkpoint_stage <<- as.integer(completed_stages)
    runtime_log_event(
      "priority_phase_elapsed",
      stage = completed_stages,
      phase = "checkpoint",
      elapsed_seconds = round(proc.time()[["elapsed"]] - checkpoint_started, 2)
    )

    invisible(checkpoint_path)
  }


  # -------------------------------------------------------------------
  # 9. Run the global stage loop
  # -------------------------------------------------------------------

  # Run until the frontier is exhausted or max_stages is reached.
  while (!isTRUE(frontier_exhausted) && stage_index < max_stages) {
    # Advance to the next global stage.
    stage_index <- stage_index + 1L
    stage_started <- proc.time()[["elapsed"]]

    runtime_log_event(
      "priority_stage_start",
      stage = as.integer(stage_index),
      max_stages = as.character(max_stages),
      alive = as.integer(current_alive_cell_count),
      patches = as.integer(nrow(patch_table)),
      pus = as.integer(data.table::uniqueN(patch_table, by = c("species", "pu_id")))
    )


    # ---------------------------------------------------------------
    # 8A. Run the pruning stage
    # ---------------------------------------------------------------

    # Run repeated pruning iterations for this global stage.
    pruning_started <- proc.time()[["elapsed"]]
    pruning_stage_result <- run_pruning_stage(
      stage_index = stage_index,
      pruning_iterations_per_stage = pruning_iterations_per_stage,
      cells_to_remove_per_iteration = cells_to_remove_per_iteration,
      patch_table = patch_table,
      pu_graphs_by_key = pu_graphs_by_key,
      alive_species_count_by_cell = alive_species_count_by_cell,
      patch_id_by_species_env = patch_id_by_species_env,
      patch_cell_index_by_species_env = patch_cell_index_by_species_env,
      cell_species_index              = cell_species_index,
      cell_area_by_cell = cell_area_by_cell,
      frontier_state = frontier_state,
      rook_neighbor_index = rook_neighbor_index,
      species_params = species_params,
      ecology_log_every_iterations = ecology_log_every_iterations,
      ecology_previous_snapshot = ecology_previous_snapshot,
      initial_alive_cell_count = initial_alive_cell_count,
      current_alive_cell_count = current_alive_cell_count
    )
    pruning_elapsed <- proc.time()[["elapsed"]] - pruning_started
    log_priority_phase_elapsed(
      stage_index = stage_index,
      phase = "pruning",
      elapsed_seconds = pruning_elapsed,
      timing_seconds = pruning_stage_result$timing_seconds
    )

    patch_table <- pruning_stage_result$patch_table

    pu_graphs_by_key <- pruning_stage_result$pu_graphs

    alive_species_count_by_cell <- pruning_stage_result$alive_species_count_by_cell

    frontier_state <- pruning_stage_result$frontier_state

    ecology_previous_snapshot <- pruning_stage_result$ecology_previous_snapshot
    current_alive_cell_count <- pruning_stage_result$current_alive_cell_count

    # Record whether pruning exhausted the frontier in this stage.
    frontier_exhausted <- isTRUE(pruning_stage_result$frontier_exhausted)
    
    # Record removal order for each pruning iteration in chronological order.
    if (length(pruning_stage_result$removed_cells_each_iteration)) {
      for (iteration_index in seq_along(pruning_stage_result$removed_cells_each_iteration)) {
        record_removed_cells(
          removed_cells = pruning_stage_result$removed_cells_each_iteration[[iteration_index]],
          stage_index = stage_index,
          event_type = "prune",
          pruning_iteration = iteration_index
        )
      }
    }


    # ---------------------------------------------------------------
    # 8B. Run the fragmentation stage
    # ---------------------------------------------------------------

    log_priority_pipeline_boundary_start(
      stage_index = stage_index,
      boundary_step = "fragmentation_start",
      workload_label = "input_changed_patches",
      workload_count = nrow(pruning_stage_result$changed_patches_in_stage)
    )

    # Run fragmentation repair on patches changed by pruning.
    fragmentation_hotspots <- list()
    fragmentation_hotspot_callback <- function(
      species_index,
      species_count,
      species_name,
      ...
    ) {
      fields <- list(...)
      fragmentation_hotspots[[length(fragmentation_hotspots) + 1L]] <<-
        c(list(species = species_name), fields)
      invisible(NULL)
    }
    fragmentation_started <- proc.time()[["elapsed"]]
    fragmentation_stage_result <- run_fragmentation_stage(
      stage_index = stage_index,
      changed_patches_in_stage = pruning_stage_result$changed_patches_in_stage,
      patch_table = patch_table,
      pu_graphs_by_key = pu_graphs_by_key,
      alive_species_count_by_cell = alive_species_count_by_cell,
      patch_id_by_species_env = patch_id_by_species_env,
      patch_cell_index_by_species_env = patch_cell_index_by_species_env,
      species_params = species_params,
      cell_area_by_cell = cell_area_by_cell,
      rook_neighbor_index = rook_neighbor_index,
      completion_callback = fragmentation_hotspot_callback
    )
    fragmentation_elapsed <- proc.time()[["elapsed"]] - fragmentation_started
    log_priority_phase_elapsed(
      stage_index = stage_index,
      phase = "fragmentation",
      elapsed_seconds = fragmentation_elapsed,
      timing_seconds = fragmentation_stage_result$timing_seconds
    )
    log_fragmentation_species_hotspots(
      stage_index,
      fragmentation_hotspots
    )

    patch_table <- fragmentation_stage_result$patch_table

    pu_graphs_by_key <- fragmentation_stage_result$pu_graphs

    alive_species_count_by_cell <- fragmentation_stage_result$alive_species_count_by_cell

    current_alive_cell_count <- as.integer(
      current_alive_cell_count -
        length(fragmentation_stage_result$removed_cells_in_fragmentation)
    )

    frontier_state <- update_frontier_state_after_empty_cells(
      frontier_state = frontier_state,
      newly_empty_cells = fragmentation_stage_result$removed_cells_in_fragmentation,
      rook_neighbor_index = rook_neighbor_index
    )

    # Record fragmentation removals as one removal-order event.
    record_removed_cells(
      removed_cells = fragmentation_stage_result$removed_cells_in_fragmentation,
      stage_index = stage_index,
      event_type = "fragmentation",
      pruning_iteration = NA_integer_
    )


    # ---------------------------------------------------------------
    # 8C. Run the distance-connectivity stage
    # ---------------------------------------------------------------

    log_priority_pipeline_boundary_start(
      stage_index = stage_index,
      boundary_step = "distance_start",
      workload_label = "input_recheck_patches",
      workload_count = nrow(fragmentation_stage_result$patches_requiring_distance_recheck)
    )

    # Run distance-based connectivity repair on patches flagged by fragmentation.
    distance_hotspots <- list()
    distance_hotspot_callback <- function(
      species_index,
      species_count,
      species_name,
      profile
    ) {
      distance_hotspots[[length(distance_hotspots) + 1L]] <<- list(
        species = species_name,
        elapsed_seconds = profile$timing[["elapsed_seconds"]],
        profile = profile
      )
      invisible(NULL)
    }
    distance_started <- proc.time()[["elapsed"]]
    distance_stage_result <- run_distance_connectivity_stage(
      stage_index = stage_index,
      patches_requiring_distance_recheck =
        fragmentation_stage_result$patches_requiring_distance_recheck,
      patch_table = patch_table,
      pu_graphs_by_key = pu_graphs_by_key,
      alive_species_count_by_cell = alive_species_count_by_cell,
      patch_id_by_species_env = patch_id_by_species_env,
      patch_cell_index_by_species_env = patch_cell_index_by_species_env,
      species_params = species_params,
      template_raster = template_raster,
      completion_callback = distance_hotspot_callback
    )
    distance_elapsed <- proc.time()[["elapsed"]] - distance_started
    log_priority_phase_elapsed(
      stage_index = stage_index,
      phase = "distance",
      elapsed_seconds = distance_elapsed,
      timing_seconds = distance_stage_result$timing_seconds
    )
    log_distance_species_hotspots(stage_index, distance_hotspots)

    patch_table <- distance_stage_result$patch_table

    pu_graphs_by_key <- distance_stage_result$pu_graphs

    alive_species_count_by_cell <- distance_stage_result$alive_species_count_by_cell

    current_alive_cell_count <- as.integer(
      current_alive_cell_count - length(distance_stage_result$removed_cells_in_distance)
    )
    if (current_alive_cell_count < 0L) {
      stop("Incremental alive-cell count became negative during stage repair.")
    }

    frontier_state <- update_frontier_state_after_empty_cells(
      frontier_state = frontier_state,
      newly_empty_cells = distance_stage_result$removed_cells_in_distance,
      rook_neighbor_index = rook_neighbor_index
    )

    # Record distance-stage removals as one removal-order event.
    record_removed_cells(
      removed_cells = distance_stage_result$removed_cells_in_distance,
      stage_index = stage_index,
      event_type = "distance",
      pruning_iteration = NA_integer_
    )


    # ---------------------------------------------------------------
    # 8D. Write the patch lookup table for this completed full stage
    # ---------------------------------------------------------------

    write_stage_patch_lookup_table(
      patch_table = patch_table,
      stage_index = stage_index,
      patch_lookup_output_dir = patch_lookup_output_dir
    )


    # ---------------------------------------------------------------
    # 8E. Print one end-of-stage pipeline summary
    # ---------------------------------------------------------------

    # Count remaining patch-table rows.
    remaining_patch_count <- nrow(patch_table)

    # Count remaining unique species-PU combinations.
    remaining_pu_count <- data.table::uniqueN(
      patch_table,
      by = c("species", "pu_id")
    )

    # Count remaining globally alive cells.
    remaining_alive_cell_count <- current_alive_cell_count

    log_priority_pipeline_stage(
      stage_index = stage_index,
      frontier_exhausted = frontier_exhausted,
      removal_step_count = removal_step_counter - 1L,
      remaining_patch_count = remaining_patch_count,
      remaining_pu_count = remaining_pu_count,
      remaining_alive_cell_count = remaining_alive_cell_count,
      initial_alive_cell_count = initial_alive_cell_count
    )
    runtime_log_event(
      "priority_stage_elapsed",
      stage = stage_index,
      elapsed_seconds = round(proc.time()[["elapsed"]] - stage_started, 2)
    )


    # ---------------------------------------------------------------
    # 8F. Record completed stage count
    # ---------------------------------------------------------------

    completed_stages <- as.integer(stage_index)

    save_stage_checkpoint(force = FALSE)
  }

  save_stage_checkpoint(force = TRUE)

  runtime_log_event(
    "priority_final_outputs_start",
    completed_stages = as.integer(completed_stages),
    frontier_exhausted = isTRUE(frontier_exhausted),
    removal_steps = as.integer(removal_step_counter) - 1L
  )


  # -------------------------------------------------------------------
  # 9. Assign the terminal retained layer after the stopping condition
  # -------------------------------------------------------------------
  
  # Zonation rank maps assign every analysis-domain cell a priority rank.
  # Therefore, regardless of whether this run stopped because the frontier was
  # exhausted or because max_stages was reached, assign every still-unordered
  # initially alive cell to one terminal removal-order layer. Use a synthetic
  # stage after the last completed stage so stage_meta rows remain aligned with
  # actual patch_lookup_tables/stage_patch_lookup_stage_####.csv files.
  final_retained_cell_count <- record_final_retained_cells(
    stage_index = as.integer(completed_stages + 1L)
  )


  # -------------------------------------------------------------------
  # 10. Stage and atomically commit the final output set
  # -------------------------------------------------------------------
  staged_removal_order_path <- tempfile(
    "removal_order_", tmpdir = output_dir, fileext = ".tif"
  )
  staged_rankmap_path <- tempfile(
    "rankmap_", tmpdir = output_dir, fileext = ".tif"
  )
  staged_removal_events_path <- tempfile(
    "removal_events_", tmpdir = output_dir, fileext = ".csv"
  )
  staged_final_paths <- c(
    staged_removal_order_path,
    staged_rankmap_path,
    staged_removal_events_path
  )
  on.exit(unlink(staged_final_paths, force = TRUE), add = TRUE)

  # Combine one-row event tables, or create an empty event table if no cells were removed.
  removal_events <- if (length(removal_event_rows)) {
    data.table::rbindlist(
      removal_event_rows,
      use.names = TRUE,
      fill = TRUE
    )
  } else {
    data.table::data.table(
      removal_step = integer(),
      stage = integer(),
      event_type = character(),
      pruning_iteration = integer(),
      cells_removed = integer(),
      area_removed_km2 = numeric()
    )
  }
  
  # Add cumulative summaries if at least one removal event occurred.
  if (nrow(removal_events)) {
    data.table::setorder(removal_events, removal_step)
  
    # Cumulative number of initially alive cells removed.
    removal_events[, cum_cells_removed := cumsum(cells_removed)]
  
    # Cumulative area removed.
    removal_events[, cum_area_removed_km2 := cumsum(area_removed_km2)]
  
    # Cumulative proportion of initially alive cells removed.
    removal_events[, cum_prop_cells_removed :=
      cum_cells_removed / as.numeric(initial_alive_cell_count)
    ]
  
    # Cumulative proportion of initially alive area removed.
    removal_events[, cum_prop_area_removed :=
      cum_area_removed_km2 / as.numeric(initial_alive_area_km2)
    ]
  
    # Remaining initially alive cells after this event.
    removal_events[, cells_retained :=
      as.integer(initial_alive_cell_count - cum_cells_removed)
    ]
  
    # Remaining initially alive area after this event.
    removal_events[, area_retained_km2 :=
      as.numeric(initial_alive_area_km2 - cum_area_removed_km2)
    ]
  } else {
    # Add the same columns to the empty table so downstream code has a stable schema.
    removal_events[, `:=`(
      cum_cells_removed = integer(),
      cum_area_removed_km2 = numeric(),
      cum_prop_cells_removed = numeric(),
      cum_prop_area_removed = numeric(),
      cells_retained = integer(),
      area_retained_km2 = numeric()
    )]
  }

  validate_removal_events_table(removal_events)

  # Validate the complete in-memory raster/event contract once, then reuse the
  # returned integer vector for both output rasters.
  removal_order_by_cell <- validate_priority_final_cell_state(
    removal_order_by_cell = removal_order_by_cell,
    initial_alive_by_cell = initial_alive_by_cell,
    removal_events = removal_events
  )

  write_removal_order_raster(
    removal_order_by_cell = removal_order_by_cell,
    template_raster = mask_template_raster,
    output_path = staged_removal_order_path,
    validate_values = FALSE
  )

  data.table::fwrite(
    x = removal_events,
    file = staged_removal_events_path
  )
  
  # proportions. The final retained layer maps to 1. Cells outside the initial
  # alive domain remain NA.
  write_rankmap_raster(
    removal_order_by_cell = removal_order_by_cell,
    removal_events = removal_events,
    template_raster = mask_template_raster,
    output_path = staged_rankmap_path,
    value_col = "cum_prop_cells_removed",
    validate_values = FALSE
  )

  validate_staged_priority_raster(
    staged_removal_order_path,
    mask_template_raster,
    "removal-order raster",
    expected_datatype = "INT4S"
  )
  validate_staged_priority_raster(
    staged_rankmap_path,
    mask_template_raster,
    "rank-map raster",
    expected_datatype = "FLT4S"
  )
  runtime_log_event(
    "priority_final_raster_datatypes_validated",
    removal_order = "INT4S",
    rankmap = "FLT4S"
  )
  assert(
    file.exists(staged_removal_events_path) &&
      file.info(staged_removal_events_path)$size[[1L]] > 0,
    "Staged removal-events CSV is missing or empty."
  )
  priority_commit_file_set(
    staged_paths = staged_final_paths,
    target_paths = c(
      removal_order_output_path,
      rankmap_output_path,
      removal_events_output_path
    ),
    overwrite = isTRUE(resume_requested),
    label = "Stage 6 final outputs"
  )

  runtime_log_event(
    "priority_final_outputs_done",
    removal_order = normalizePath(removal_order_output_path, mustWork = FALSE),
    rankmap = normalizePath(rankmap_output_path, mustWork = FALSE),
    removal_events = normalizePath(removal_events_output_path, mustWork = FALSE),
    final_retained_cells = as.integer(final_retained_cell_count)
  )


  # -------------------------------------------------------------------
  # 11. Compute compact final-run summaries
  # -------------------------------------------------------------------

  # Count initially alive cells that received a removal order.
  removed_initial_alive_cell_count <- sum(
    initial_alive_by_cell &
      !is.na(removal_order_by_cell)
  )

  # Count initially alive cells that still have no removal order.
  unremoved_initial_alive_cell_count <- sum(
    initial_alive_by_cell &
      is.na(removal_order_by_cell)
  )


  # -------------------------------------------------------------------
  # 12. Return compact output metadata
  # -------------------------------------------------------------------

  # The durable artifacts are the scientific contract. Do not return the large
  # mutable graph, patch, cell, or frontier objects after they are published.
  list(
    completed_stages = completed_stages,
    frontier_exhausted = frontier_exhausted,
    removal_order_path = removal_order_output_path,
    rankmap_path = rankmap_output_path,
    removal_events_path = removal_events_output_path,
    patch_lookup_output_dir = patch_lookup_output_dir,
    removal_step_count = removal_step_counter - 1L,
    initial_alive_cell_count = initial_alive_cell_count,
    removed_initial_alive_cell_count = removed_initial_alive_cell_count,
    unremoved_initial_alive_cell_count = unremoved_initial_alive_cell_count,
    final_remaining_patches = nrow(patch_table),
    final_remaining_pus = data.table::uniqueN(
      patch_table, by = c("species", "pu_id")
    ),
    final_alive_cells = sum(alive_species_count_by_cell > 0L)
  )
}
