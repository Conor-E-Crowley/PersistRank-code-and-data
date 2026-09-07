# Incremental Stage 7.2 benchmark reconstruction.
#
# Input state contains the current patch table, hashed CSR graph store, flat
# patch-to-cell indexes, alive-species counts, cell areas, species thresholds,
# and rook-neighbour index. One call removes the requested cells, applies the
# shared strict patch/PU thresholds, repairs fragmentation and distance
# connectivity, and returns updated state plus canonical LUT rows only for
# touched species. Species, patch, PU, graph-node, adjacency, and raster-cell
# ordering follow the Stage 6 contracts.
#
# The retained benchmark landscape is nested, so unchanged species are never
# reconstructed. Compact indexes and graph stores are the authoritative mutable
# runtime representation; Stage 7.2 never materializes full species rasters.
# Stages and species remain sequential to preserve deterministic IDs and bound
# peak memory.


# Stable aggregate diagnostics avoid growing event tables inside hot loops and
# are formatted only once per completed stage.
stage72_incremental_timing <- function() {
  stats::setNames(
    numeric(26L),
    c(
      "rank_cell_lookup_seconds", "area_threshold_seconds",
      "fragmentation_seconds", "fragmentation_graph_repair_seconds",
      "distance_candidate_edge_setup_seconds",
      "distance_candidate_index_extract_seconds",
      "distance_candidate_extent_seconds", "distance_local_raster_seconds",
      "distance_polygonize_seconds", "distance_sf_conversion_seconds",
      "distance_geometry_validation_seconds",
      "distance_geometry_repair_seconds", "distance_predicate_spatial_seconds",
      "distance_predicate_mapping_seconds",
      "distance_pu_index_seconds", "distance_pu_repair_seconds",
      "distance_graph_commit_seconds", "distance_compact_index_seconds",
      "distance_cell_state_seconds", "distance_state_commit_seconds",
      "distance_other_seconds",
      "state_commit_seconds", "lut_bookkeeping_seconds",
      "lut_canonicalization_seconds", "lut_bind_seconds",
      "validation_write_seconds"
    )
  )
}

stage72_incremental_timing_totals <- function(timing) {
  assert(
    identical(names(timing), names(stage72_incremental_timing())) &&
      all(is.finite(timing) & timing >= 0),
    "Stage 7.2 incremental timings must be finite named non-negative values."
  )
  distance_non_graph_fields <- setdiff(
    grep("^distance_.*_seconds$", names(timing), value = TRUE),
    c("distance_pu_repair_seconds", "distance_graph_commit_seconds")
  )
  c(
    distance_geometry_seconds = sum(timing[distance_non_graph_fields]),
    graph_repair_seconds = sum(timing[c(
      "fragmentation_graph_repair_seconds",
      "distance_pu_repair_seconds", "distance_graph_commit_seconds"
    )]),
    lut_assembly_seconds = sum(timing[c(
      "lut_bookkeeping_seconds", "lut_canonicalization_seconds",
      "lut_bind_seconds"
    )]),
    total_seconds = sum(timing)
  )
}

stage72_add_distance_timing <- function(timing, distance_timing) {
  expected <- empty_distance_timing()
  assert(
    is.numeric(distance_timing) && identical(names(distance_timing), names(expected)) &&
      all(is.finite(distance_timing) & distance_timing >= 0),
    "Stage 7.2 distance timings must match the shared distance timing contract."
  )
  destination <- paste0("distance_", names(distance_timing))
  timing[destination] <- timing[destination] + distance_timing
  timing
}

stage72_touched_patch_rows <- function(patch_table, touched_species) {
  touched_species <- as.character(touched_species)
  if (!length(touched_species)) return(list())
  species <- as.character(patch_table$species)
  positions <- which(data.table::`%chin%`(species, touched_species))
  grouped <- split(
    positions,
    factor(species[positions], levels = touched_species),
    drop = FALSE
  )
  stats::setNames(lapply(grouped, as.integer), touched_species)
}

stage72_incremental_counts <- function() {
  stats::setNames(
    numeric(26L),
    c(
      "removed_rank_cells", "touched_species", "touched_patches",
      "fragmented_patches", "new_fragments", "candidate_geometry_patches",
      "rechecked_edges", "dropped_patches", "dropped_pus",
      "reused_species", "lut_rows", "unchanged_distance_pus",
      "compiled_filtered_pus", "lazy_materialized_pus",
      "invalid_edges_removed", "adjacency_entries_scanned",
      "adjacency_entries_avoided", "geometry_features_checked",
      "invalid_geometry_features", "geometry_repair_invocations",
      "graphs_inspected", "graph_nodes_scanned", "recheck_rows",
      "csr_adjacency_entries_visited", "overlay_edges_scanned",
      "unique_recheck_edges"
    )
  )
}

stage72_log_incremental_timing <- function(stage, timing, counts) {
  assert(
    identical(names(timing), names(stage72_incremental_timing())) &&
      all(is.finite(timing) & timing >= 0),
    "Stage 7.2 incremental timings must be finite named non-negative values."
  )
  assert(
    identical(names(counts), names(stage72_incremental_counts())) &&
      all(is.finite(counts) & counts >= 0),
    "Stage 7.2 incremental counters must be finite named non-negative values."
  )
  totals <- stage72_incremental_timing_totals(timing)
  candidate_cell_seconds <- sum(timing[c(
    "distance_candidate_edge_setup_seconds",
    "distance_candidate_index_extract_seconds",
    "distance_candidate_extent_seconds"
  )])
  predicate_seconds <- sum(timing[c(
    "distance_predicate_spatial_seconds",
    "distance_predicate_mapping_seconds"
  )])
  runtime_log_event(
    "stage72_incremental_timing", console = FALSE, stage = as.integer(stage),
    rank_cell_lookup_seconds = sprintf("%.3f", timing[["rank_cell_lookup_seconds"]]),
    area_threshold_seconds = sprintf("%.3f", timing[["area_threshold_seconds"]]),
    fragmentation_seconds = sprintf("%.3f", timing[["fragmentation_seconds"]]),
    fragmentation_graph_repair_seconds = sprintf(
      "%.3f", timing[["fragmentation_graph_repair_seconds"]]
    ),
    distance_candidate_cell_seconds = sprintf("%.3f", candidate_cell_seconds),
    distance_candidate_edge_setup_seconds = sprintf(
      "%.3f", timing[["distance_candidate_edge_setup_seconds"]]
    ),
    distance_candidate_index_extract_seconds = sprintf(
      "%.3f", timing[["distance_candidate_index_extract_seconds"]]
    ),
    distance_candidate_extent_seconds = sprintf(
      "%.3f", timing[["distance_candidate_extent_seconds"]]
    ),
    distance_local_raster_seconds = sprintf("%.3f", timing[["distance_local_raster_seconds"]]),
    distance_polygonize_seconds = sprintf("%.3f", timing[["distance_polygonize_seconds"]]),
    distance_sf_conversion_seconds = sprintf("%.3f", timing[["distance_sf_conversion_seconds"]]),
    distance_geometry_validation_seconds = sprintf(
      "%.3f", timing[["distance_geometry_validation_seconds"]]
    ),
    distance_geometry_repair_seconds = sprintf("%.3f", timing[["distance_geometry_repair_seconds"]]),
    distance_predicate_seconds = sprintf("%.3f", predicate_seconds),
    distance_predicate_spatial_seconds = sprintf(
      "%.3f", timing[["distance_predicate_spatial_seconds"]]
    ),
    distance_predicate_mapping_seconds = sprintf(
      "%.3f", timing[["distance_predicate_mapping_seconds"]]
    ),
    distance_pu_index_seconds = sprintf("%.3f", timing[["distance_pu_index_seconds"]]),
    distance_pu_repair_seconds = sprintf("%.3f", timing[["distance_pu_repair_seconds"]]),
    distance_graph_commit_seconds = sprintf("%.3f", timing[["distance_graph_commit_seconds"]]),
    distance_compact_index_seconds = sprintf("%.3f", timing[["distance_compact_index_seconds"]]),
    distance_cell_state_seconds = sprintf("%.3f", timing[["distance_cell_state_seconds"]]),
    distance_state_commit_seconds = sprintf("%.3f", timing[["distance_state_commit_seconds"]]),
    distance_other_seconds = sprintf("%.3f", timing[["distance_other_seconds"]]),
    distance_geometry_seconds = sprintf("%.3f", totals[["distance_geometry_seconds"]]),
    graph_repair_seconds = sprintf("%.3f", totals[["graph_repair_seconds"]]),
    state_commit_seconds = sprintf("%.3f", timing[["state_commit_seconds"]]),
    lut_bookkeeping_seconds = sprintf("%.3f", timing[["lut_bookkeeping_seconds"]]),
    lut_canonicalization_seconds = sprintf("%.3f", timing[["lut_canonicalization_seconds"]]),
    lut_bind_seconds = sprintf("%.3f", timing[["lut_bind_seconds"]]),
    lut_assembly_seconds = sprintf("%.3f", totals[["lut_assembly_seconds"]]),
    validation_write_seconds = sprintf("%.3f", timing[["validation_write_seconds"]]),
    total_seconds = sprintf("%.3f", totals[["total_seconds"]]),
    removed_rank_cells = as.integer(counts[["removed_rank_cells"]]),
    touched_species = as.integer(counts[["touched_species"]]),
    touched_patches = as.integer(counts[["touched_patches"]]),
    fragmented_patches = as.integer(counts[["fragmented_patches"]]),
    new_fragments = as.integer(counts[["new_fragments"]]),
    candidate_geometry_patches = as.integer(counts[["candidate_geometry_patches"]]),
    rechecked_edges = as.integer(counts[["rechecked_edges"]]),
    dropped_patches = as.integer(counts[["dropped_patches"]]),
    dropped_pus = as.integer(counts[["dropped_pus"]]),
    reused_species = as.integer(counts[["reused_species"]]),
    lut_rows = as.integer(counts[["lut_rows"]]),
    unchanged_distance_pus = as.integer(counts[["unchanged_distance_pus"]]),
    compiled_filtered_pus = as.integer(counts[["compiled_filtered_pus"]]),
    lazy_materialized_pus = as.integer(counts[["lazy_materialized_pus"]]),
    invalid_edges_removed = as.integer(counts[["invalid_edges_removed"]]),
    adjacency_entries_scanned = sprintf(
      "%.0f", counts[["adjacency_entries_scanned"]]
    ),
    adjacency_entries_avoided = sprintf(
      "%.0f", counts[["adjacency_entries_avoided"]]
    ),
    graphs_inspected = as.integer(counts[["graphs_inspected"]]),
    graph_nodes_scanned = sprintf("%.0f", counts[["graph_nodes_scanned"]]),
    recheck_rows = as.integer(counts[["recheck_rows"]]),
    csr_adjacency_entries_visited = sprintf(
      "%.0f", counts[["csr_adjacency_entries_visited"]]
    ),
    overlay_edges_scanned = sprintf(
      "%.0f", counts[["overlay_edges_scanned"]]
    ),
    unique_recheck_edges = sprintf(
      "%.0f", counts[["unique_recheck_edges"]]
    ),
    geometry_features_checked = sprintf(
      "%.0f", counts[["geometry_features_checked"]]
    ),
    invalid_geometry_features = as.integer(
      counts[["invalid_geometry_features"]]
    ),
    geometry_repair_invocations = as.integer(
      counts[["geometry_repair_invocations"]]
    )
  )
  invisible(totals[["total_seconds"]])
}

stage72_progress_callback <- function(stage, phase, cadence = 10L) {
  cadence <- as.integer(cadence)
  assert(
    length(cadence) == 1L && !is.na(cadence) && cadence > 0L,
    "Stage 7.2 progress cadence must be a positive integer."
  )
  phase_started <- proc.time()[["elapsed"]]
  function(index, total, scientific_name) {
    index <- as.integer(index)
    total <- as.integer(total)
    if (index == 1L || index == total || index %% cadence == 0L) {
      runtime_log_event(
        "stage72_progress", stage = as.integer(stage),
        phase = phase,
        species = paste0(index, "/", total),
        scientific_name = scientific_name,
        phase_elapsed_s = sprintf(
          "%.1f", proc.time()[["elapsed"]] - phase_started
        )
      )
    }
    invisible(NULL)
  }
}

stage72_log_phase <- function(stage, phase, status, ..., started = NULL) {
  fields <- c(
    list(
      event = "stage72_phase", stage = as.integer(stage),
      phase = phase, status = status
    ),
    list(...)
  )
  if (!is.null(started)) {
    fields$elapsed_seconds <- sprintf(
      "%.1f", proc.time()[["elapsed"]] - started
    )
  }
  do.call(runtime_log_event, fields)
  invisible(NULL)
}

stage72_log_fragmentation_workload <- function(
  stage, workload, total_live_index_entries, timing_seconds,
  index_timing_seconds = NULL
) {
  workload_names <- names(workload)
  workload <- as.numeric(workload)
  names(workload) <- workload_names
  assert(all(is.finite(workload) & workload >= 0),
         "Stage 7.2 fragmentation workload counters must be finite and non-negative.")
  total_live_index_entries <- as.numeric(total_live_index_entries)
  assert(length(total_live_index_entries) == 1L &&
           is.finite(total_live_index_entries) && total_live_index_entries >= 0,
         "Stage 7.2 live-index workload must be finite and non-negative.")
  value <- function(field) {
    if (field %in% names(workload)) workload[[field]] else 0
  }
  # `timing_seconds` is the same named vector run_fragmentation_stage() already
  # computes every stage (see empty_fragmentation_timing() in
  # priority_fragmentation_stage.R); this just surfaces the full per-phase
  # breakdown that previously collapsed into "component_analysis_seconds" and
  # "wrapper_seconds" alone. wrapper_seconds remains the sum of every phase
  # below it, kept as an explicit cross-check total.
  timing_names <- names(timing_seconds)
  timing_seconds <- as.numeric(timing_seconds)
  names(timing_seconds) <- timing_names
  assert(all(is.finite(timing_seconds) & timing_seconds >= 0),
         "Stage 7.2 fragmentation timings must be finite and non-negative.")
  timing_value <- function(field) {
    if (field %in% names(timing_seconds)) timing_seconds[[field]] else 0
  }
  # `index_timing_seconds` is the compact-index splice sub-breakdown that
  # run_fragmentation_stage() already accumulates every stage (see
  # empty_fragmentation_index_timing() / update_species_patch_index_incremental()
  # in priority_fragmentation_patches.R). It splits compact_index_seconds above
  # into: time building the R-side replacement rows (index_flatten_seconds),
  # time inside the compiled splice kernel itself (index_kernel_seconds), and
  # time packaging the kernel's results back into R (index_packaging_seconds).
  # index_total_seconds is kept as an explicit cross-check total, analogous to
  # wrapper_seconds above.
  if (is.null(index_timing_seconds)) {
    index_timing_seconds <- empty_fragmentation_index_timing()
  }
  index_timing_names <- names(index_timing_seconds)
  index_timing_seconds <- as.numeric(index_timing_seconds)
  names(index_timing_seconds) <- index_timing_names
  assert(all(is.finite(index_timing_seconds) & index_timing_seconds >= 0),
         "Stage 7.2 compact-index timings must be finite and non-negative.")
  index_timing_value <- function(field) {
    if (field %in% names(index_timing_seconds)) {
      index_timing_seconds[[field]]
    } else {
      0
    }
  }
  runtime_log_event(
    "stage72_fragmentation_workload", stage = as.integer(stage),
    affected_species = sprintf("%.0f", value("affected_species")),
    changed_patches = sprintf("%.0f", value("changed_patches")),
    total_live_index_entries = sprintf("%.0f", total_live_index_entries),
    changed_index_entries = sprintf(
      "%.0f", value("changed_index_entries")
    ),
    full_index_entries_scanned = sprintf(
      "%.0f", value("full_index_entries_scanned")
    ),
    dense_vectors_materialized = sprintf(
      "%.0f", value("dense_vectors_materialized")
    ),
    dense_cells_materialized = sprintf(
      "%.0f", value("dense_cells_materialized")
    ),
    dense_cells_avoided = sprintf(
      "%.0f", value("dense_cells_avoided")
    ),
    component_analysis_seconds = sprintf(
      "%.3f", timing_value("component_relabel_seconds")
    ),
    patch_threshold_seconds = sprintf(
      "%.3f", timing_value("patch_threshold_seconds")
    ),
    pu_index_seconds = sprintf("%.3f", timing_value("pu_index_seconds")),
    pu_repair_seconds = sprintf("%.3f", timing_value("pu_repair_seconds")),
    graph_commit_seconds = sprintf(
      "%.3f", timing_value("graph_commit_seconds")
    ),
    compact_index_seconds = sprintf(
      "%.3f", timing_value("compact_index_seconds")
    ),
    index_flatten_seconds = sprintf(
      "%.3f", index_timing_value("flatten_seconds")
    ),
    index_kernel_seconds = sprintf(
      "%.3f", index_timing_value("kernel_seconds")
    ),
    index_packaging_seconds = sprintf(
      "%.3f", index_timing_value("result_packaging_seconds")
    ),
    index_total_seconds = sprintf(
      "%.3f", index_timing_value("total_seconds")
    ),
    cell_state_seconds = sprintf("%.3f", timing_value("cell_state_seconds")),
    state_commit_seconds = sprintf(
      "%.3f", timing_value("state_commit_seconds")
    ),
    other_seconds = sprintf("%.3f", timing_value("other_seconds")),
    wrapper_seconds = sprintf(
      "%.3f", value("wrapper_seconds")
    )
  )
  invisible(NULL)
}

stage72_log_predicate_workload <- function(stage, workload, candidate_patches) {
  workload_names <- names(workload)
  workload <- as.numeric(workload)
  names(workload) <- workload_names
  candidate_patches <- as.numeric(candidate_patches)
  assert(all(is.finite(workload) & workload >= 0) &&
           length(candidate_patches) == 1L && is.finite(candidate_patches) &&
           candidate_patches >= 0,
         "Stage 7.2 predicate workload counters must be finite and non-negative.")
  value <- function(field) {
    if (field %in% names(workload)) workload[[field]] else 0
  }
  all_pairs_upper_bound <- value("all_pair_candidate_upper_bound")
  graph_pairs <- value("graph_edge_pairs_evaluated")
  pair_reduction <- if (all_pairs_upper_bound > 0) {
    100 * (1 - graph_pairs / all_pairs_upper_bound)
  } else {
    0
  }
  runtime_log_event(
    "stage72_predicate_workload", stage = as.integer(stage),
    recheck_patches = sprintf("%.0f", value("recheck_patches")),
    candidate_patches_before_filter = sprintf(
      "%.0f", value("candidate_patches_before_filter")
    ),
    candidate_patches_after_filter = sprintf("%.0f", candidate_patches),
    edge_endpoint_patches = sprintf(
      "%.0f", value("edge_endpoint_patches")
    ),
    isolated_rechecks_avoided = sprintf(
      "%.0f", value("isolated_rechecks_avoided")
    ),
    zero_edge_species = sprintf(
      "%.0f", value("zero_edge_species")
    ),
    all_pair_candidate_upper_bound = sprintf("%.0f", all_pairs_upper_bound),
    graph_edge_pairs_evaluated = sprintf("%.0f", graph_pairs),
    pair_reduction_pct = sprintf("%.2f", max(0, pair_reduction)),
    valid_edges = sprintf("%.0f", value("valid_edges")),
    invalid_edges = sprintf("%.0f", value("invalid_edges")),
    predicate_spatial_seconds = sprintf(
      "%.3f", value("predicate_spatial_seconds")
    ),
    predicate_mapping_seconds = sprintf(
      "%.3f", value("predicate_mapping_seconds")
    )
  )
  invisible(NULL)
}

stage72_geometry_reuse_diagnostic <- function(
  stage, current_candidates, previous_candidates, recheck_patches
) {
  started <- proc.time()[["elapsed"]]
  current_candidates <- current_candidates %||% list()
  previous_available <- !is.null(previous_candidates)
  previous_candidates <- previous_candidates %||% list()
  recheck_patches <- data.table::as.data.table(recheck_patches)
  species_names <- sort(names(current_candidates))
  totals <- stats::setNames(numeric(7L), c(
    "current_candidates", "previous_candidates", "overlap",
    "invalidated_overlap", "potential_hits", "potential_misses",
    "candidate_species"
  ))
  totals[["candidate_species"]] <- length(species_names)

  for (scientific_name in species_names) {
    current <- sort(unique(as.integer(current_candidates[[scientific_name]])))
    previous <- sort(unique(as.integer(
      previous_candidates[[scientific_name]] %||% integer()
    )))
    overlap <- intersect(current, previous)
    invalidated <- intersect(
      overlap,
      recheck_patches[species == scientific_name, patch_id]
    )
    hits <- setdiff(overlap, invalidated)
    totals[["current_candidates"]] <-
      totals[["current_candidates"]] + length(current)
    totals[["previous_candidates"]] <-
      totals[["previous_candidates"]] + length(previous)
    totals[["overlap"]] <- totals[["overlap"]] + length(overlap)
    totals[["invalidated_overlap"]] <-
      totals[["invalidated_overlap"]] + length(invalidated)
    totals[["potential_hits"]] <-
      totals[["potential_hits"]] + length(hits)
    totals[["potential_misses"]] <-
      totals[["potential_misses"]] + length(current) - length(hits)
  }
  hit_pct <- if (totals[["current_candidates"]] > 0) {
    100 * totals[["potential_hits"]] / totals[["current_candidates"]]
  } else {
    0
  }
  runtime_log_event(
    "stage72_geometry_reuse", stage = as.integer(stage),
    previous_stage = if (previous_available) as.integer(stage) - 1L else "NA",
    previous_stage_available = tolower(as.character(previous_available)),
    candidate_species = sprintf("%.0f", totals[["candidate_species"]]),
    recheck_species = data.table::uniqueN(recheck_patches$species),
    recheck_patches = nrow(recheck_patches),
    current_candidates = sprintf("%.0f", totals[["current_candidates"]]),
    previous_candidates = sprintf("%.0f", totals[["previous_candidates"]]),
    overlap = sprintf("%.0f", totals[["overlap"]]),
    invalidated_overlap = sprintf("%.0f", totals[["invalidated_overlap"]]),
    potential_hits = sprintf("%.0f", totals[["potential_hits"]]),
    potential_misses = sprintf("%.0f", totals[["potential_misses"]]),
    potential_hit_pct = sprintf("%.2f", hit_pct),
    diagnostic_seconds = sprintf(
      "%.3f", proc.time()[["elapsed"]] - started
    )
  )
  invisible(totals)
}

# Compact-index and LUT contracts. Compact entries are patch-ID sorted, cells
# are unique, and compiled canonicalization restores stable output row order.
stage72_validate_compact_index <- function(index, species_name, n_cells) {
  assert(
    is.list(index) && all(c("pid", "cell") %in% names(index)),
    paste0("Missing compact patch index for species: ", species_name)
  )
  pid <- as.integer(index$pid)
  cell <- as.integer(index$cell)
  assert(
    length(pid) == length(cell) && length(pid) > 0L &&
      !anyNA(pid) && !anyNA(cell) && all(pid > 0L) &&
      all(cell > 0L & cell <= n_cells) && !anyDuplicated(cell) &&
      !is.unsorted(pid),
    paste0("Invalid compact patch index for species: ", species_name)
  )
  list(pid = pid, cell = cell)
}

stage72_canonical_species_lut_from_rows <- function(
  scientific_name, safe_species, species_rows, patch_index
) {
  rows <- data.table::as.data.table(species_rows)
  if (!nrow(rows)) {
    assert(
      is.list(patch_index) && !length(patch_index$pid) && !length(patch_index$cell),
      paste0("Patch table and compact index disagree for species: ", scientific_name)
    )
    return(empty_rank_lut())
  }
  need_cols(
    rows, c("patch_id", "pu_id", "patch_area_km2"),
    paste0("Stage 7.2 patch rows for ", scientific_name)
  )
  index <- stage72_validate_compact_index(
    patch_index, scientific_name, max(patch_index$cell)
  )
  mapping <- tryCatch(
    stage6_canonicalize_rank_lut_cpp(
      patch_ids = as.integer(rows$patch_id),
      pu_ids = as.integer(rows$pu_id),
      index_pid = index$pid,
      index_cell = index$cell
    ),
    error = function(error) {
      stop(
        "Canonical LUT construction failed for species ", scientific_name,
        ": ", conditionMessage(error), call. = FALSE
      )
    }
  )
  row_order <- as.integer(mapping$row_order)
  assert(
    identical(sort(row_order), seq_len(nrow(rows))),
    paste0("Canonical LUT row mapping is invalid for species: ", scientific_name)
  )
  data.table::data.table(
    stage = rep(NA_integer_, nrow(rows)),
    method = rep("rank", nrow(rows)),
    scientificName = rep(scientific_name, nrow(rows)),
    species = rep(safe_species, nrow(rows)),
    patch_id = as.integer(mapping$patch_id),
    pu_id = as.integer(mapping$pu_id),
    patch_area_km2 = as.numeric(rows$patch_area_km2[row_order]),
    patch_n_cells = as.integer(mapping$patch_n_cells)
  )
}

stage72_bundle_compact_indexes <- function(bundle, sp, n_cells) {
  scientific_names <- as.character(sp$scientificName)
  indexes <- bundle$patch_cell_index_by_species_list
  stats::setNames(lapply(scientific_names, function(scientific_name) {
    index <- indexes[[scientific_name]] %||% indexes[[species_id(scientific_name)]]
    stage72_validate_compact_index(index, scientific_name, n_cells)
  }), scientific_names)
}

# Fresh state retains only fields required by consecutive-stage repair.
stage72_initialize_incremental_state <- function(bundle, sp, n_cells,
                                                 indexes = NULL) {
  scientific_names <- as.character(sp$scientificName)
  if (is.null(indexes)) {
    indexes <- stage72_bundle_compact_indexes(bundle, sp, n_cells)
  }
  patch_table <- data.table::copy(data.table::as.data.table(bundle$patch_table))
  if (!"species" %in% names(patch_table) && "scientificName" %in% names(patch_table)) {
    patch_table[, species := as.character(scientificName)]
  }
  patch_table <- patch_table[species %in% scientific_names, .(
    species = as.character(species),
    patch_id = as.integer(patch_id),
    pu_id = as.integer(pu_id),
    patch_area_km2 = as.numeric(patch_area_km2)
  )]
  params <- data.table::copy(data.table::as.data.table(bundle$species_params))[
    species %in% scientific_names
  ]
  data.table::setkey(params, species)
  validate_priority_species_parameters(params, "Stage 7.2 incremental species parameters")
  rook_index <- build_rook_neighbor_index(bundle$rook_neighbor_pairs, n_cells)
  list(
    patch_table = patch_table,
    graphs = new_priority_graph_store(bundle$pu_graphs_by_key),
    indexes = indexes,
    alive_counts = as.integer(bundle$alive_species_count_by_cell),
    cell_area = as.numeric(bundle$cell_area_by_cell),
    species_params = params,
    rook_index = rook_index,
    n_cells = as.integer(n_cells)
  )
}

# Apply direct patch losses and any resulting PU component/area losses before
# fragmentation. Graph actions are validated as one transaction before commit.
stage72_repair_direct_thresholds <- function(state, touched_patches) {
  touched_patches <- unique(data.table::as.data.table(touched_patches)[, .(
    species = as.character(species), patch_id = as.integer(patch_id)
  )], by = c("species", "patch_id"))
  empty_drops <- data.table::data.table(species = character(), patch_id = integer())
  empty_pus <- data.table::data.table(species = character(), pu_id = integer())
  if (!nrow(touched_patches)) {
    return(list(state = state, dropped_patches = empty_drops, dropped_pus = empty_pus))
  }

  pre_rows <- state$patch_table[touched_patches, on = .(species, patch_id), nomatch = 0L]
  assert(
    nrow(pre_rows) == nrow(touched_patches),
    "A benchmark-touched patch is missing from the incremental patch table."
  )
  direct_drops <- find_pruning_patch_threshold_drops(
    state$patch_table, touched_patches, state$species_params
  )
  affected_original_pus <- unique(pre_rows[, .(species, pu_id)])
  affected_mapping <- state$patch_table[
    affected_original_pus, on = .(species, pu_id), nomatch = 0L,
    .(species, patch_id, pu_id)
  ]

  component_drops <- empty_drops
  component_assignments <- data.table::data.table(
    species = character(), patch_id = integer(), pu_id = integer()
  )
  graph_actions <- list()
  pre_step_max <- list()

  if (nrow(direct_drops)) {
    repair_pus <- unique(merge(
      direct_drops, affected_mapping, by = c("species", "patch_id")
    )[, .(species, pu_id)])
    data.table::setorder(repair_pus, species, pu_id)
    repair_mapping <- affected_mapping[
      repair_pus, on = .(species, pu_id), nomatch = 0L
    ]
    graph_context <- snapshot_pruning_affected_graphs(repair_mapping, state$graphs)
    component_drop_records <- list()
    assignment_records <- list()

    for (species_name in unique(repair_pus$species)) {
      species_mapping <- repair_mapping[species == species_name]
      pu_ids <- unique(species_mapping$pu_id)
      all_species_mapping <- state$patch_table[
        species == species_name, .(patch_id, pu_id)
      ]
      # Sparse patch IDs are expected after fragmentation. Lookup gaps must be
      # NA because zero is a genuinely invalid represented PU ID.
      pre_lookup <- rep.int(
        NA_integer_,
        max(all_species_mapping$patch_id)
      )
      pre_lookup[all_species_mapping$patch_id] <- all_species_mapping$pu_id
      cache_entry <- list(pu_id_by_patch_id = pre_lookup)
      next_pu_id <- pruning_pre_step_pu_high_water(cache_entry, species_name)
      pre_step_max[[species_name]] <- next_pu_id
      patch_area <- state$patch_table[species == species_name, .(patch_id, patch_area_km2)]
      pu_threshold <- state$species_params[
        species == species_name, min_population_area_km2
      ]
      assert(
        length(pu_threshold) == 1L,
        paste0("Missing population-area threshold for species: ", species_name)
      )

      for (pu_id in pu_ids) {
        key <- paste0(species_name, "|", pu_id)
        graph <- graph_context$graph_snapshot[[key]]
        dead <- direct_drops[species == species_name, patch_id]
        alive_nodes <- !graph$id2patch %in% dead
        area_match <- match(graph$id2patch, patch_area$patch_id)
        node_area <- numeric(length(graph$id2patch))
        node_area[!is.na(area_match)] <- patch_area$patch_area_km2[area_match[!is.na(area_match)]]
        rebuilt <- rebuild_pu_after_patch_loss(
          pu_graph = graph,
          alive_nodes = alive_nodes,
          patch_area = node_area,
          pu_area_threshold = pu_threshold,
          next_available_pu_id = next_pu_id
        )
        next_pu_id <- rebuilt$next_available_pu_id
        if (length(rebuilt$dropped_patch_ids)) {
          component_drop_records[[length(component_drop_records) + 1L]] <-
            data.table::data.table(
              species = species_name,
              patch_id = as.integer(rebuilt$dropped_patch_ids)
            )
        }
        survivors <- rebuilt$surviving_pu_graphs
        if (!length(survivors)) {
          graph_actions[[length(graph_actions) + 1L]] <- list(
            kind = "remove_original", key = key, graph = NULL
          )
        } else {
          for (j in seq_along(survivors)) {
            survivor <- survivors[[j]]
            staged <- list(
              species = species_name,
              pu_id = as.integer(survivor$pu_id),
              id2patch = as.integer(survivor$id2patch),
              row_ptr = as.integer(survivor$row_ptr),
              col_idx = as.integer(survivor$col_idx)
            )
            graph_actions[[length(graph_actions) + 1L]] <- list(
              kind = if (j == 1L) "replace_original" else "add_new",
              key = if (j == 1L) key else paste0(species_name, "|", staged$pu_id),
              graph = staged
            )
            assignment_records[[length(assignment_records) + 1L]] <-
              data.table::data.table(
                species = species_name,
                patch_id = staged$id2patch,
                pu_id = staged$pu_id
              )
          }
        }
      }
    }
    if (length(component_drop_records)) {
      component_drops <- unique(data.table::rbindlist(component_drop_records),
                                by = c("species", "patch_id"))
    }
    if (length(assignment_records)) {
      component_assignments <- data.table::rbindlist(assignment_records)
    }
    validate_pruning_graph_transaction(
      affected_patch_to_pu_before_step = repair_mapping,
      direct_patch_drops = direct_drops,
      component_drops = component_drops,
      component_assignments = component_assignments,
      graph_actions = graph_actions,
      pu_graphs_before = state$graphs,
      pre_step_max_pu_by_species = pre_step_max
    )
  }

  all_component_drops <- unique(data.table::rbindlist(
    list(direct_drops, component_drops), use.names = TRUE
  ), by = c("species", "patch_id"))
  candidate_table <- apply_pruning_component_changes(
    state$patch_table,
    component_assignments = component_assignments,
    component_drops = all_component_drops
  )

  candidate_pus <- unique(candidate_table[
    affected_mapping[, .(species, patch_id)],
    on = .(species, patch_id), nomatch = 0L,
    .(species, pu_id)
  ])
  pu_areas <- candidate_table[
    candidate_pus, on = .(species, pu_id), nomatch = 0L,
    .(pu_area_km2 = sum(patch_area_km2)), by = .EACHI
  ]
  dropped_pus <- if (nrow(pu_areas)) {
    state$species_params[pu_areas, on = .(species)][
      !area_exceeds_threshold(pu_area_km2, min_population_area_km2),
      .(species, pu_id)
    ]
  } else empty_pus
  pu_drop_rows <- if (nrow(dropped_pus)) {
    candidate_table[dropped_pus, on = .(species, pu_id), nomatch = 0L,
                    .(species, patch_id)]
  } else empty_drops
  if (nrow(dropped_pus)) {
    for (key in paste0(dropped_pus$species, "|", dropped_pus$pu_id)) {
      graph_actions[[length(graph_actions) + 1L]] <- list(
        kind = "remove", key = key, graph = NULL
      )
    }
    candidate_table <- candidate_table[!dropped_pus, on = .(species, pu_id)]
  }

  all_drops <- unique(data.table::rbindlist(
    list(all_component_drops, pu_drop_rows), use.names = TRUE
  ), by = c("species", "patch_id"))
  state$patch_table <- candidate_table
  state$graphs <- priority_graph_apply_actions(state$graphs, graph_actions)
  list(state = state, dropped_patches = all_drops, dropped_pus = dropped_pus)
}

# Advance one completed benchmark stage. The four visible scientific phases
# are rank/area update, threshold repair, fragmentation, and distance repair.
stage72_apply_incremental_stage <- function(
  state, removed_cells_by_species, safe_species_by_name, stage_index,
  template_raster, build_luts = TRUE
) {
  timing <- stage72_incremental_timing()
  counts <- stage72_incremental_counts()
  stage_started <- proc.time()[["elapsed"]]
  index_environment <- new.env(parent = emptyenv())
  changed_records <- list()
  before_patch_ids <- list()
  before_pu_ids <- list()
  touched_species <- character()
  working_patch_table <- NULL

  rank_area_started <- proc.time()[["elapsed"]]
  stage72_log_phase(
    stage_index, "rank_area", "start",
    candidate_species = length(removed_cells_by_species)
  )
  rank_area_progress <- stage72_progress_callback(stage_index, "rank_area")
  requested_species <- names(removed_cells_by_species)
  for (species_index in seq_along(requested_species)) {
    scientific_name <- requested_species[[species_index]]
    rank_area_progress(
      species_index, length(requested_species), scientific_name
    )
    lookup_started <- proc.time()[["elapsed"]]
    index <- state$indexes[[scientific_name]]
    positions <- match(as.integer(removed_cells_by_species[[scientific_name]]), index$cell)
    live_positions <- positions[!is.na(positions)]
    timing[["rank_cell_lookup_seconds"]] <- timing[["rank_cell_lookup_seconds"]] +
      proc.time()[["elapsed"]] - lookup_started
    if (!length(live_positions)) next

    live_cells <- as.integer(index$cell[live_positions])
    live_patch_ids <- as.integer(index$pid[live_positions])
    if (is.null(working_patch_table)) {
      working_patch_table <- data.table::copy(state$patch_table)
    }
    touched_species <- c(touched_species, scientific_name)
    before_rows <- working_patch_table[species == scientific_name]
    before_patch_ids[[scientific_name]] <- before_rows$patch_id
    before_pu_ids[[scientific_name]] <- before_rows$pu_id

    threshold_started <- proc.time()[["elapsed"]]
    decrements <- data.table::data.table(
      patch_id = live_patch_ids,
      area_removed = state$cell_area[live_cells]
    )[, .(area_removed = sum(area_removed)), by = patch_id]
    matched <- working_patch_table[
      data.table::data.table(species = scientific_name, patch_id = decrements$patch_id),
      on = .(species, patch_id), nomatch = 0L, .N
    ]
    assert(
      identical(as.integer(matched), nrow(decrements)),
      paste0("Removed benchmark cells reference missing live patches for species: ", scientific_name)
    )
    working_patch_table[
      data.table::data.table(
        species = scientific_name,
        patch_id = decrements$patch_id,
        area_removed = decrements$area_removed
      ),
      on = .(species, patch_id),
      patch_area_km2 := patch_area_km2 - i.area_removed
    ]
    keep <- rep(TRUE, length(index$cell))
    keep[live_positions] <- FALSE
    updated_index <- list(pid = index$pid[keep], cell = index$cell[keep])
    exact_areas <- data.table::data.table(
      patch_id = updated_index$pid,
      cell = updated_index$cell
    )[patch_id %in% decrements$patch_id,
      .(exact_area = sum(state$cell_area[cell])), by = patch_id]
    exact_areas <- merge(
      data.table::data.table(patch_id = decrements$patch_id),
      exact_areas, by = "patch_id", all.x = TRUE, sort = FALSE
    )
    exact_areas[is.na(exact_area), exact_area := 0]
    observed_after <- working_patch_table[
      data.table::data.table(species = scientific_name, patch_id = exact_areas$patch_id),
      on = .(species, patch_id), patch_area_km2
    ]
    tolerance <- 1e-8 * pmax(1, abs(exact_areas$exact_area))
    assert(
      all(abs(observed_after - exact_areas$exact_area) <= tolerance),
      paste0("Incremental patch-area conservation failed for species: ", scientific_name)
    )
    working_patch_table[
      data.table::data.table(
        species = scientific_name,
        patch_id = exact_areas$patch_id,
        exact_area = exact_areas$exact_area
      ),
      on = .(species, patch_id),
      patch_area_km2 := i.exact_area
    ]
    assert(
      all(working_patch_table[species == scientific_name, patch_area_km2] >= -1e-9),
      paste0("Incremental patch-area underflow for species: ", scientific_name)
    )
    timing[["area_threshold_seconds"]] <- timing[["area_threshold_seconds"]] +
      proc.time()[["elapsed"]] - threshold_started

    assign(scientific_name, updated_index, envir = index_environment)
    changed_records[[length(changed_records) + 1L]] <- data.table::data.table(
      species = scientific_name,
      patch_id = sort(unique(live_patch_ids))
    )
  }

  if (!length(touched_species)) {
    stage72_log_phase(
      stage_index, "rank_area", "complete",
      " touched_species=0", started = rank_area_started
    )
    stage72_log_fragmentation_workload(
      stage = stage_index,
      workload = c(
        affected_species = 0, changed_patches = 0,
        changed_index_entries = 0, full_index_entries_scanned = 0,
        dense_vectors_materialized = 0, dense_cells_materialized = 0,
        dense_cells_avoided = 0, wrapper_seconds = 0
      ),
      total_live_index_entries = 0,
      timing_seconds = empty_fragmentation_timing(),
      index_timing_seconds = empty_fragmentation_index_timing()
    )
    stage72_log_predicate_workload(
      stage = stage_index,
      workload = stats::setNames(
        numeric(11L),
        c(
          "recheck_patches", "candidate_patches_before_filter",
          "edge_endpoint_patches", "isolated_rechecks_avoided",
          "zero_edge_species", "all_pair_candidate_upper_bound",
          "graph_edge_pairs_evaluated", "valid_edges", "invalid_edges",
          "predicate_spatial_seconds", "predicate_mapping_seconds"
        )
      ),
      candidate_patches = 0
    )
    stage72_geometry_reuse_diagnostic(
      stage = stage_index,
      current_candidates = list(),
      previous_candidates =
        state$runtime_previous_distance_candidates %||% NULL,
      recheck_patches = data.table::data.table(
        species = character(), patch_id = integer()
      )
    )
    state$runtime_previous_distance_candidates <- list()
    return(list(
      state = state,
      luts = list(),
      touched_species = character(),
      timing = timing,
      counts = counts,
      runtime_workload = list(fragmentation = numeric(), distance = numeric())
    ))
  }
  stage72_log_phase(
    stage_index, "rank_area", "complete",
    touched_species = length(touched_species),
    started = rank_area_started
  )
  state$patch_table <- working_patch_table
  changed <- unique(data.table::rbindlist(changed_records), by = c("species", "patch_id"))
  counts[["touched_species"]] <- length(touched_species)
  counts[["touched_patches"]] <- nrow(changed)

  threshold_started <- proc.time()[["elapsed"]]
  stage72_log_phase(
    stage_index, "threshold_repair", "start",
    touched_patches = nrow(changed)
  )
  repaired <- stage72_repair_direct_thresholds(state, changed)
  state <- repaired$state
  timing[["area_threshold_seconds"]] <- timing[["area_threshold_seconds"]] +
    proc.time()[["elapsed"]] - threshold_started
  counts[["dropped_patches"]] <- nrow(repaired$dropped_patches)
  counts[["dropped_pus"]] <- nrow(repaired$dropped_pus)
  stage72_log_phase(
    stage_index, "threshold_repair", "complete",
    dropped_patches = nrow(repaired$dropped_patches),
    dropped_pus = nrow(repaired$dropped_pus),
    started = threshold_started
  )

  for (scientific_name in touched_species) {
    index <- get(scientific_name, envir = index_environment, inherits = FALSE)
    dropped_ids <- repaired$dropped_patches[species == scientific_name, patch_id]
    if (length(dropped_ids)) {
      dropped_cells <- index$cell[index$pid %in% dropped_ids]
      if (length(dropped_cells)) {
        positive <- state$alive_counts[dropped_cells] > 0L
        state$alive_counts[dropped_cells[positive]] <-
          state$alive_counts[dropped_cells[positive]] - 1L
      }
      keep <- !index$pid %in% dropped_ids
      index <- list(pid = index$pid[keep], cell = index$cell[keep])
    }
    assign(scientific_name, index, envir = index_environment)
  }
  changed_survivors <- changed[state$patch_table, on = .(species, patch_id), nomatch = 0L,
                               .(species, patch_id)]
  fragmentation_started <- proc.time()[["elapsed"]]
  total_live_index_entries <- 0
  for (scientific_name in touched_species) {
    total_live_index_entries <- total_live_index_entries + length(get(
      scientific_name,
      envir = index_environment,
      inherits = FALSE
    )$cell)
  }
  stage72_log_phase(
    stage_index, "fragmentation", "start",
    changed_patches = nrow(changed_survivors)
  )
  fragmented <- run_fragmentation_stage(
    stage_index = stage_index,
    changed_patches_in_stage = changed_survivors,
    patch_table = state$patch_table,
    pu_graphs_by_key = state$graphs,
    alive_species_count_by_cell = state$alive_counts,
    patch_id_by_species_env = NULL,
    patch_cell_index_by_species_env = index_environment,
    species_params = state$species_params,
    cell_area_by_cell = state$cell_area,
    rook_neighbor_index = state$rook_index,
    state_representation = "compact",
    collect_runtime_diagnostics = TRUE,
    emit_logs = FALSE,
    progress_callback = stage72_progress_callback(stage_index, "fragmentation")
  )
  fragmentation_elapsed <- proc.time()[["elapsed"]] - fragmentation_started
  state$patch_table <- fragmented$patch_table
  state$graphs <- fragmented$pu_graphs
  state$alive_counts <- fragmented$alive_species_count_by_cell
  stage72_log_phase(
    stage_index, "fragmentation", "complete",
    recheck_patches = nrow(fragmented$patches_requiring_distance_recheck),
    started = fragmentation_started
  )
  stage72_log_fragmentation_workload(
    stage = stage_index,
    workload = fragmented$workload %||% numeric(),
    total_live_index_entries = total_live_index_entries,
    timing_seconds = fragmented$timing_seconds %||% empty_fragmentation_timing(),
    index_timing_seconds = fragmented$index_timing_seconds %||%
      empty_fragmentation_index_timing()
  )

  distance_started <- proc.time()[["elapsed"]]
  stage72_log_phase(
    stage_index, "distance", "start",
    recheck_patches = nrow(fragmented$patches_requiring_distance_recheck)
  )
  distance <- run_distance_connectivity_stage(
    stage_index = stage_index,
    patches_requiring_distance_recheck = fragmented$patches_requiring_distance_recheck,
    patch_table = state$patch_table,
    pu_graphs_by_key = state$graphs,
    alive_species_count_by_cell = state$alive_counts,
    patch_id_by_species_env = NULL,
    patch_cell_index_by_species_env = index_environment,
    species_params = state$species_params,
    template_raster = template_raster,
    state_representation = "compact",
    predicate_strategy = "edge_pairs",
    collect_runtime_diagnostics = TRUE,
    emit_logs = FALSE,
    progress_callback = stage72_progress_callback(stage_index, "distance")
  )
  state$patch_table <- distance$patch_table
  state$graphs <- distance$pu_graphs
  state$alive_counts <- distance$alive_species_count_by_cell
  stage72_log_phase(
    stage_index, "distance", "complete",
    rechecked_edges = distance$rechecked_edges %||% 0L,
    started = distance_started
  )
  stage72_log_predicate_workload(
    stage = stage_index,
    workload = distance$workload %||% numeric(),
    candidate_patches =
      distance$geometry_counts[["candidate_patches"]] %||% 0
  )
  stage72_geometry_reuse_diagnostic(
    stage = stage_index,
    current_candidates = distance$candidate_ids_by_species,
    previous_candidates = state$runtime_previous_distance_candidates %||% NULL,
    recheck_patches = fragmented$patches_requiring_distance_recheck
  )
  state$runtime_previous_distance_candidates <-
    distance$candidate_ids_by_species

  fragmentation_graph <- sum(fragmented$timing_seconds[
    intersect(c("pu_repair_seconds", "graph_commit_seconds"),
              names(fragmented$timing_seconds))
  ])
  timing[["fragmentation_graph_repair_seconds"]] <- fragmentation_graph
  timing <- stage72_add_distance_timing(timing, distance$timing_seconds)
  timing[["fragmentation_seconds"]] <- max(0, fragmentation_elapsed - fragmentation_graph)
  counts[["fragmented_patches"]] <- nrow(changed_survivors)
  counts[["candidate_geometry_patches"]] <-
    distance$geometry_counts[["candidate_patches"]] %||%
    nrow(fragmented$patches_requiring_distance_recheck)
  counts[["geometry_features_checked"]] <-
    distance$geometry_counts[["geometry_features_checked"]] %||% 0
  counts[["invalid_geometry_features"]] <-
    distance$geometry_counts[["invalid_geometry_features"]] %||% 0
  counts[["geometry_repair_invocations"]] <-
    distance$geometry_counts[["geometry_repair_invoked"]] %||% 0
  counts[["rechecked_edges"]] <- distance$rechecked_edges %||% 0L
  if (!is.null(distance$path_counts)) {
    for (field in c(
      "unchanged_distance_pus", "compiled_filtered_pus",
      "lazy_materialized_pus", "invalid_edges_removed",
      "adjacency_entries_scanned", "adjacency_entries_avoided",
      "graphs_inspected", "graph_nodes_scanned", "recheck_rows",
      "csr_adjacency_entries_visited", "overlay_edges_scanned",
      "unique_recheck_edges"
    )) {
      counts[[field]] <- distance$path_counts[[field]] %||% 0
    }
  }

  bookkeeping_started <- proc.time()[["elapsed"]]
  patch_rows <- stage72_touched_patch_rows(state$patch_table, touched_species)
  total_new_fragments <- 0L
  total_dropped_patches <- 0L
  total_dropped_pus <- 0L
  for (scientific_name in touched_species) {
    index <- get0(scientific_name, envir = index_environment, inherits = FALSE)
    if (is.null(index)) index <- list(pid = integer(), cell = integer())
    state$indexes[[scientific_name]] <- index
    row_positions <- patch_rows[[scientific_name]]
    after_patch_ids <- state$patch_table$patch_id[row_positions]
    after_pu_ids <- state$patch_table$pu_id[row_positions]
    total_new_fragments <- total_new_fragments +
      length(setdiff(after_patch_ids, before_patch_ids[[scientific_name]]))
    total_dropped_patches <- total_dropped_patches +
      length(setdiff(before_patch_ids[[scientific_name]], after_patch_ids))
    total_dropped_pus <- total_dropped_pus +
      length(setdiff(before_pu_ids[[scientific_name]], after_pu_ids))
  }
  counts[["new_fragments"]] <- counts[["new_fragments"]] + total_new_fragments
  counts[["dropped_patches"]] <- max(counts[["dropped_patches"]], total_dropped_patches)
  counts[["dropped_pus"]] <- max(counts[["dropped_pus"]], total_dropped_pus)
  timing[["lut_bookkeeping_seconds"]] <-
    proc.time()[["elapsed"]] - bookkeeping_started

  canonicalization_started <- proc.time()[["elapsed"]]
  luts <- list()
  # Stage 7.2 persists patch LUTs and therefore requests canonical rows. Stage
  # 7.4 consumes the repaired state directly and disables this output-only
  # materialization; all scientific repair steps above are identical.
  if (isTRUE(build_luts)) {
    luts <- stats::setNames(vector("list", length(touched_species)), touched_species)
    for (scientific_name in touched_species) {
      rows <- state$patch_table[
        patch_rows[[scientific_name]],
        .(patch_id, pu_id, patch_area_km2)
      ]
      luts[[scientific_name]] <- stage72_canonical_species_lut_from_rows(
        scientific_name, safe_species_by_name[[scientific_name]],
        rows, state$indexes[[scientific_name]]
      )
    }
  }
  timing[["lut_canonicalization_seconds"]] <-
    proc.time()[["elapsed"]] - canonicalization_started
  accounted <- sum(timing[
    setdiff(names(timing), c("state_commit_seconds", "validation_write_seconds"))
  ])
  timing[["state_commit_seconds"]] <- max(
    0, proc.time()[["elapsed"]] - stage_started - accounted
  )
  list(
    state = state,
    luts = luts,
    touched_species = touched_species,
    timing = timing,
    counts = counts,
    # Optional consumers such as Stage 7.4 aggregate counters already produced
    # by the repair stages. Returning them performs no additional traversal.
    runtime_workload = list(
      fragmentation = fragmented$workload %||% numeric(),
      distance = distance$workload %||% numeric()
    )
  )
}
