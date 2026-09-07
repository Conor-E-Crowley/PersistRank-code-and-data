# Stage 6 pipeline-level console diagnostics.
#
# These helpers format existing append-only events. They receive compact
# counters and top-five species profiles only; no scientific state or raster
# data is copied for logging.


format_pipeline_pct <- function(numerator, denominator, digits = 1L) {
  if (!is.finite(numerator) || !is.finite(denominator) || denominator <= 0) {
    return("n/a")
  }

  sprintf(paste0("%.", as.integer(digits), "f%%"), 100 * numerator / denominator)
}

log_priority_start <- function(
  cells_to_remove_per_iteration,
  pruning_iterations_per_stage,
  initial_alive_cell_count,
  species_count
) {
  runtime_log_event(
    "priority_start",
    objective = "per_area_marginal_persistence_loss",
    loss10 = "log10_loss",
    removal_order = "lower_loss_removed_first",
    batch = as.integer(cells_to_remove_per_iteration),
    stage_iters = as.integer(pruning_iterations_per_stage),
    initial_alive = as.integer(initial_alive_cell_count),
    species = as.integer(species_count)
  )

  invisible(NULL)
}

log_priority_runtime_strategy <- function() {
  runtime_log_event(
    "priority_runtime_strategy",
    patch_index = "flat_cpp",
    distance_geometry = "single_crop",
    frontier_scoring = "sparse_csr",
    pu_repair = "compiled_batch",
    species_hotspot_limit = 5L
  )
  invisible(NULL)
}

log_priority_phase_elapsed <- function(
  stage_index,
  phase,
  elapsed_seconds,
  timing_seconds = NULL
) {
  fields <- list(
    event = "priority_phase_elapsed",
    stage = as.integer(stage_index),
    phase = as.character(phase),
    elapsed_seconds = round(as.numeric(elapsed_seconds), 2L)
  )
  if (!is.null(timing_seconds)) {
    classified_seconds <- sum(as.numeric(timing_seconds))
    accounting_gap_seconds <- max(
      0,
      as.numeric(elapsed_seconds) - classified_seconds
    )
    coverage_pct <- if (elapsed_seconds > 0) {
      min(100, 100 * classified_seconds / elapsed_seconds)
    } else {
      100
    }
    fields <- c(fields, list(
      classified_seconds = round(classified_seconds, 2L),
      accounting_gap_seconds = round(accounting_gap_seconds, 2L),
      timer_coverage_pct = round(coverage_pct, 1L)
    ))
  }
  do.call(runtime_log_event, fields)
  invisible(NULL)
}

priority_top_hotspots <- function(records, limit = 5L) {
  if (!length(records)) return(list())
  elapsed <- vapply(
    records,
    function(record) as.numeric(record$elapsed_seconds),
    numeric(1L)
  )
  species <- vapply(records, `[[`, character(1L), "species")
  records[utils::head(order(-elapsed, species), as.integer(limit))]
}

log_fragmentation_species_hotspots <- function(
  stage_index,
  records,
  limit = 5L
) {
  hotspots <- priority_top_hotspots(records, limit)
  for (rank_index in seq_along(hotspots)) {
    record <- hotspots[[rank_index]]
    runtime_log_event(
      "fragmentation_species_hotspot",
      stage = as.integer(stage_index),
      rank = as.integer(rank_index),
      species = record$species,
      elapsed_seconds = round(record$elapsed_seconds, 3L),
      changed_patches = as.integer(record$changed_patches),
      changed_patch_cells = as.numeric(record$changed_patch_cells),
      affected_pus = as.integer(record$affected_pus),
      full_rebuild_pus = as.integer(record$full_rebuild_pus),
      component_seconds = round(record$component_seconds, 3L),
      patch_threshold_seconds = round(
        record$patch_threshold_seconds,
        3L
      ),
      pu_index_seconds = round(record$pu_index_seconds, 3L),
      full_rebuild_seconds = round(record$full_rebuild_seconds, 3L),
      compact_index_seconds = round(record$compact_index_seconds, 3L),
      graph_commit_seconds = round(record$graph_commit_seconds, 3L),
      cell_state_seconds = round(record$cell_state_seconds, 3L),
      other_seconds = round(record$other_seconds, 3L)
    )
  }
  invisible(NULL)
}

log_distance_species_hotspots <- function(
  stage_index,
  records,
  limit = 5L
) {
  hotspots <- priority_top_hotspots(records, limit)
  for (rank_index in seq_along(hotspots)) {
    record <- hotspots[[rank_index]]
    workload <- record$profile$workload
    timing <- record$profile$timing
    runtime_log_event(
      "distance_species_hotspot",
      stage = as.integer(stage_index),
      rank = as.integer(rank_index),
      species = record$species,
      elapsed_seconds = round(timing[["elapsed_seconds"]], 3L),
      candidate_patches = as.integer(workload[["candidate_patches"]]),
      candidate_cells = as.numeric(workload[["candidate_cells"]]),
      local_raster_cells = as.numeric(workload[["local_raster_cells"]]),
      raster_expansion_ratio = round(
        workload[["raster_expansion_ratio"]],
        2L
      ),
      rechecked_edges = as.numeric(workload[["rechecked_edges"]]),
      candidate_edge_setup_seconds = round(
        timing[["candidate_edge_setup_seconds"]],
        3L
      ),
      candidate_index_extract_seconds = round(
        timing[["candidate_index_extract_seconds"]],
        3L
      ),
      candidate_extent_seconds = round(
        timing[["candidate_extent_seconds"]],
        3L
      ),
      local_raster_seconds = round(timing[["local_raster_seconds"]], 3L),
      polygonize_seconds = round(timing[["polygonize_seconds"]], 3L),
      sf_conversion_seconds = round(timing[["sf_conversion_seconds"]], 3L),
      geometry_validation_seconds = round(
        timing[["geometry_validation_seconds"]],
        3L
      ),
      predicate_seconds = round(
        timing[["predicate_spatial_seconds"]] +
          timing[["predicate_mapping_seconds"]],
        3L
      ),
      edge_filter_seconds = round(timing[["edge_filter_seconds"]], 3L),
      pu_index_seconds = round(timing[["pu_index_seconds"]], 3L),
      component_rebuild_seconds = round(
        timing[["component_rebuild_seconds"]],
        3L
      ),
      pu_assembly_seconds = round(timing[["pu_assembly_seconds"]], 3L),
      graph_commit_seconds = round(timing[["graph_commit_seconds"]], 3L),
      compact_index_seconds = round(timing[["compact_index_seconds"]], 3L),
      cell_state_seconds = round(timing[["cell_state_seconds"]], 3L),
      other_seconds = round(timing[["other_seconds"]], 3L)
    )
  }
  invisible(NULL)
}

log_priority_pipeline_stage <- function(
  stage_index,                 # completed global stage index
  frontier_exhausted,          # whether pruning exhausted the frontier
  removal_step_count,          # number of removal-order events recorded so far
  remaining_patch_count,       # number of patch-table rows remaining
  remaining_pu_count,          # number of species-PU combinations remaining
  remaining_alive_cell_count,  # number of cells with at least one species remaining
  initial_alive_cell_count     # number of initially alive cells
) {
  status <- if (isTRUE(frontier_exhausted)) "frontier_exhausted" else "running"

  runtime_log_event(
    "priority_stage_done",
    stage = as.integer(stage_index),
    removal_steps = as.integer(removal_step_count),
    remaining_alive = as.integer(remaining_alive_cell_count),
    remaining_alive_pct = format_pipeline_pct(
      remaining_alive_cell_count,
      initial_alive_cell_count,
      1L
    ),
    remaining_patches = as.integer(remaining_patch_count),
    remaining_pus = as.integer(remaining_pu_count),
    status = status
  )

  invisible(NULL)
}


# Print one boundary-start pipeline message.

log_priority_pipeline_boundary_start <- function(
  stage_index,                 # global stage index about to enter a boundary step
  boundary_step,               # boundary step label, e.g. fragmentation_start
  workload_label,              # name of the input workload count
  workload_count               # size of the input workload
) {
  fields <- list(stage = as.integer(stage_index))
  fields[[as.character(workload_label)]] <- as.integer(workload_count)

  do.call(
    runtime_log_event,
    c(list(event = as.character(boundary_step)), fields)
  )

  invisible(NULL)
}


# Print one completed-stage pipeline message.
