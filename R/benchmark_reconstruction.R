# Shared exact-target benchmark reconstruction and recovery.
#
# The mutable state is the same compact patch/index/CSR representation used by
# Stage 7.2. Targets are visited once in descending retained-cell order. Each
# exact retained-cell state is committed once. Persistence is
# intentionally absent from this module so reporting can never initialize the
# rank raster, compiled kernels, or recovery checkpoint.


# Load the one scenario-owned initialization and hand it back unmodified.
#
# species_params remains curve-neutral: it contains one alpha/beta pair per
# named curve and no bare a_pred/b_pred columns. Downstream consumers select
# their required curves only after the shared schema has been validated.
load_benchmark_initialization <- function(config) {
  reference <- readRDS(config$paths$initialization_bundle)
  validate_priority_initialization_schema(reference, "Stage 6 shared initialization")
  validate_priority_initialization_metadata(
    reference, config$taxa_tag, config$sdm,
    "Stage 6 shared initialization",
    contract = config$contract
  )
  reference
}

benchmark_runtime_species <- function(bundle) {
  names <- as.character(bundle$metadata$retained_species)
  data.table::data.table(
    scientificName = names,
    species = species_id(names),
    species_order = seq_along(names)
  )
}

benchmark_checkpoint_metadata <- function(config, plan, retained_species, n_cells) {
  list(
    schema = benchmark_reconstruction_schema(),
    library_manifest = expected_benchmark_library_manifest(config),
    union_target_id = as.integer(plan$union$union_target_id),
    target_keep_n = as.integer(plan$union$keep_n),
    retained_species = as.character(retained_species),
    n_cells = as.integer(n_cells)
  )
}

validate_benchmark_mutable_state <- function(state, metadata) {
  required <- c("patch_table", "indexes", "graphs", "alive_counts")
  assert(is.list(state) && all(required %in% names(state)),
         "Benchmark checkpoint mutable state is incomplete.")
  patches <- data.table::as.data.table(state$patch_table)
  need_cols(patches, c("species", "patch_id", "pu_id", "patch_area_km2"),
            "benchmark checkpoint patch table")
  assert(!anyDuplicated(patches[, .(species, patch_id)]) &&
           all(patches$patch_id > 0L) && all(patches$pu_id > 0L) &&
           all(is.finite(patches$patch_area_km2) & patches$patch_area_km2 > 0),
         "Benchmark checkpoint patch rows are invalid.")
  indexes <- state$indexes
  assert(is.list(indexes) && identical(names(indexes), metadata$retained_species),
         "Benchmark checkpoint compact-index species/order is invalid.")
  expected_alive <- integer(metadata$n_cells)
  for (scientific_name in metadata$retained_species) {
    index <- indexes[[scientific_name]]
    assert(is.list(index) && all(c("pid", "cell") %in% names(index)),
           paste0("Missing benchmark compact index for ", scientific_name, "."))
    index <- list(pid = as.integer(index$pid), cell = as.integer(index$cell))
    assert(length(index$pid) == length(index$cell) && !anyNA(index$pid) &&
             !anyNA(index$cell) && all(index$pid > 0L) &&
             all(index$cell > 0L & index$cell <= metadata$n_cells) &&
             !anyDuplicated(index$cell) && !is.unsorted(index$pid),
           paste0("Invalid benchmark compact index for ", scientific_name, "."))
    represented <- sort(unique(index$pid))
    expected <- sort(patches[species == scientific_name, unique(as.integer(patch_id))])
    assert(identical(represented, expected), paste0(
      "Benchmark checkpoint patch table/index mismatch for ", scientific_name, "."
    ))
    expected_alive[index$cell] <- expected_alive[index$cell] + 1L
  }
  alive <- as.integer(state$alive_counts)
  assert(identical(alive, expected_alive),
         "Benchmark checkpoint alive counts disagree with compact indexes.")
  graphs <- if (is_priority_graph_store(state$graphs)) {
    priority_graph_as_list(state$graphs)
  } else state$graphs
  expected_keys <- unique(paste0(patches$species, "|", patches$pu_id))
  assert(is.list(graphs) && !is.null(names(graphs)) &&
           setequal(names(graphs), expected_keys) && !anyDuplicated(names(graphs)),
         "Benchmark checkpoint graph keys disagree with the patch table.")
  for (key in names(graphs)) {
    graph <- graphs[[key]]
    assert(!is_lazy_added_node_overlay_graph(graph) && all(
      c("pu_id", "id2patch", "row_ptr", "col_idx") %in% names(graph)
    ), paste0("Invalid benchmark checkpoint graph: ", key))
    ids <- as.integer(graph$id2patch)
    row_ptr <- as.integer(graph$row_ptr)
    col_idx <- as.integer(graph$col_idx)
    assert(length(ids) > 0L && !anyNA(ids) && !anyDuplicated(ids) &&
             length(row_ptr) == length(ids) + 1L && row_ptr[[1L]] == 0L &&
             !is.unsorted(row_ptr) && tail(row_ptr, 1L) == length(col_idx) &&
             !anyNA(col_idx) && all(col_idx > 0L & col_idx <= length(ids)),
           paste0("Invalid benchmark checkpoint graph CSR: ", key))
    split <- strsplit(key, "|", fixed = TRUE)[[1L]]
    graph_species <- paste(split[-length(split)], collapse = "|")
    graph_pu <- as.integer(tail(split, 1L))
    expected_ids <- patches[species == graph_species & pu_id == graph_pu,
                            sort(as.integer(patch_id))]
    assert(identical(sort(ids), expected_ids), paste0(
      "Benchmark checkpoint graph membership mismatch: ", key
    ))
  }
  invisible(TRUE)
}

make_benchmark_checkpoint <- function(state, metadata, completed_target, counters) {
  serial_state <- list(
    patch_table = state$patch_table,
    indexes = state$indexes,
    graphs = priority_graph_as_list(state$graphs),
    alive_counts = as.integer(state$alive_counts)
  )
  validate_benchmark_mutable_state(serial_state, metadata)
  list(
    metadata = metadata,
    completed_union_target = as.integer(completed_target),
    state = serial_state,
    counters = counters
  )
}

validate_benchmark_checkpoint <- function(checkpoint, metadata, plan) {
  assert(is.list(checkpoint) && identical(checkpoint$metadata, metadata),
         "Benchmark reconstruction checkpoint is incompatible.")
  target <- as.integer(checkpoint$completed_union_target)
  assert(length(target) == 1L && target %in% plan$union$union_target_id,
         "Benchmark checkpoint completed target is invalid.")
  validate_benchmark_mutable_state(checkpoint$state, metadata)
  checkpoint
}

restore_benchmark_state <- function(checkpoint_state, bundle, sp, n_cells) {
  state <- stage72_initialize_incremental_state(bundle, sp, n_cells)
  state$patch_table <- data.table::copy(data.table::as.data.table(checkpoint_state$patch_table))
  state$indexes <- checkpoint_state$indexes
  state$graphs <- new_priority_graph_store(checkpoint_state$graphs)
  state$alive_counts <- as.integer(checkpoint_state$alive_counts)
  state
}

write_benchmark_checkpoint <- function(checkpoint, path) {
  ensure_writable_dir(dirname(path), "benchmark checkpoint directory")
  staged <- tempfile("benchmark_checkpoint_", tmpdir = dirname(path), fileext = ".rds")
  on.exit(unlink(staged, force = TRUE), add = TRUE)
  saveRDS(checkpoint, staged, compress = FALSE, version = 3L)
  # Avoid deserializing a potentially multi-GB checkpoint after every write.
  # saveRDS() is atomic with respect to this staged path; the transaction below
  # exposes the file only after serialization has completed successfully.
  assert(file.exists(staged) && file.info(staged)$size > 0,
         "Failed to stage the benchmark checkpoint.")
  project_file_set_transaction(
    staged, path, overwrite = file.exists(path), label = "benchmark checkpoint"
  )
  invisible(path)
}

empty_benchmark_checkpoint_metrics <- function() {
  stats::setNames(
    numeric(3L),
    c("prepare_seconds", "write_seconds", "total_seconds")
  )
}

# Build, validate, serialize, and report one deterministic recovery checkpoint.
# Runtime metrics remain outside the checkpoint so its serialized schema and
# scientific compatibility contract are unchanged.
save_benchmark_checkpoint <- function(
  state, metadata, completed_target, counters, path,
  checkpoint_interval_targets, reason = c("initial", "interval", "final")
) {
  reason <- match.arg(reason)
  checkpoint_started <- proc.time()[["elapsed"]]
  prepare_started <- proc.time()[["elapsed"]]
  checkpoint <- make_benchmark_checkpoint(
    state, metadata, completed_target, counters
  )
  prepare_seconds <- proc.time()[["elapsed"]] - prepare_started
  patch_rows <- nrow(checkpoint$state$patch_table)
  graph_count <- length(checkpoint$state$graphs)
  index_entries <- sum(vapply(
    checkpoint$state$indexes,
    function(index) length(index$cell),
    numeric(1L)
  ))
  write_started <- proc.time()[["elapsed"]]
  write_benchmark_checkpoint(checkpoint, path)
  write_seconds <- proc.time()[["elapsed"]] - write_started
  checkpoint_bytes <- as.numeric(file.info(path)$size)
  total_seconds <- proc.time()[["elapsed"]] - checkpoint_started
  metrics <- c(
    prepare_seconds = prepare_seconds,
    write_seconds = write_seconds,
    total_seconds = total_seconds
  )
  assert(all(is.finite(metrics) & metrics >= 0) &&
           is.finite(checkpoint_bytes) && checkpoint_bytes > 0,
         "Benchmark checkpoint diagnostics must be finite and non-negative.")
  runtime_log_event(
    "benchmark_checkpoint",
    union_target_id = as.integer(completed_target),
    reason = reason,
    checkpoint_interval_targets = as.integer(checkpoint_interval_targets),
    patch_rows = patch_rows,
    graph_count = graph_count,
    index_entries = sprintf("%.0f", index_entries),
    checkpoint_bytes = sprintf("%.0f", checkpoint_bytes),
    prepare_seconds = sprintf("%.3f", prepare_seconds),
    write_seconds = sprintf("%.3f", write_seconds),
    total_seconds = sprintf("%.3f", total_seconds)
  )
  rm(checkpoint)
  metrics
}

prepare_benchmark_rank_state <- function(config, bundle, plan, template) {
  rank <- terra::rast(config$paths$zonation_rankmap)
  assert(terra::compareGeom(rank, template, stopOnError = FALSE),
         "Zonation rank map geometry differs from the Stage 6 raster template.")
  alive0 <- as.integer(bundle$alive_species_count_by_cell) > 0L
  rank_values <- terra::values(rank, mat = FALSE)
  assert(length(rank_values) == length(alive0),
         "Zonation rank map cell count differs from the shared Stage 6 state.")
  rank_table <- build_stage7_rank_table(
    rank_values, alive0, bundle$cell_area_by_cell
  )
  assert(nrow(rank_table) == plan$union$keep_n[[1L]],
         "The target plan initial cell count differs from the ranked Stage 6 domain.")
  list(rank_table = rank_table, rank_raster = rank)
}

benchmark_rank_position <- function(ranked_cells, total_cells) {
  ranked_cells <- as.integer(ranked_cells)
  assert(!anyNA(ranked_cells) && !anyDuplicated(ranked_cells) &&
           all(ranked_cells >= 1L & ranked_cells <= total_cells),
         "Ranked cells must be unique valid template-cell IDs.")
  position <- integer(total_cells)
  position[ranked_cells] <- seq_along(ranked_cells)
  position
}

benchmark_species_prefix <- function(cells, rank_position, target_counts,
                                     scientific_name) {
  cells <- as.integer(cells)
  assert(length(cells) > 0L && !anyNA(cells) && !anyDuplicated(cells) &&
           all(cells >= 1L & cells <= length(rank_position)),
         paste0("Invalid initial cells for ", scientific_name, "."))
  positions <- rank_position[cells]
  assert(all(positions > 0L), paste0(
    scientific_name, " has cells outside the ranked initial domain."
  ))
  order_index <- order(positions)
  positions <- positions[order_index]
  list(
    cells = cells[order_index],
    retained_counts = vapply(
      target_counts, function(count) findInterval(count, positions), integer(1L)
    )
  )
}

prepare_benchmark_species_prefixes <- function(bundle, sp, rank_cells, keep_n,
                                             n_cells) {
  rank_position <- benchmark_rank_position(rank_cells, n_cells)
  prefixes <- stats::setNames(vector("list", nrow(sp)), sp$scientificName)
  counts <- matrix(0L, nrow = nrow(sp), ncol = length(keep_n),
                   dimnames = list(sp$scientificName, as.character(seq_along(keep_n) - 1L)))
  for (i in seq_len(nrow(sp))) {
    scientific_name <- sp$scientificName[[i]]
    index <- bundle$patch_cell_index_by_species_list[[scientific_name]] %||%
      bundle$patch_cell_index_by_species_list[[species_id(scientific_name)]]
    ordered <- benchmark_species_prefix(
      index$cell, rank_position, keep_n, scientific_name
    )
    prefixes[[scientific_name]] <- ordered$cells
    counts[i, ] <- ordered$retained_counts
  }
  list(cells = prefixes, counts = counts)
}

# Advance the shared benchmark state and checkpoint at configured target
# intervals, plus the final unique target.
reconstruct_benchmark_targets <- function(config, plan, bundle, template) {
  setup_started <- proc.time()[["elapsed"]]
  sp <- benchmark_runtime_species(bundle)
  n_cells <- terra::ncell(template)
  metadata <- benchmark_checkpoint_metadata(
    config, plan, sp$scientificName, n_cells
  )
  library <- benchmark_library_from_config(config)
  if (file.exists(library$checkpoint)) {
    candidate <- tryCatch(readRDS(library$checkpoint), error = function(e) NULL)
    compatible <- is.list(candidate) && identical(candidate$metadata, metadata)
    if (compatible) {
      completed_prefix <- plan$union[
        union_target_id <= as.integer(candidate$completed_union_target), keep_n
      ]
      compatible <- all(file.exists(benchmark_lookup_file(library, completed_prefix)))
    }
    if (!compatible) {
      if (file.exists(library$checkpoint_backup)) unlink(library$checkpoint_backup)
      assert(file.rename(library$checkpoint, library$checkpoint_backup),
             "Could not archive the incompatible benchmark checkpoint.")
      runtime_log_event(
        "benchmark_checkpoint",
        action = "restart_from_initial",
        completed_lookups = "preserved"
      )
    }
  }
  rank_state_started <- proc.time()[["elapsed"]]
  rank_state <- prepare_benchmark_rank_state(config, bundle, plan, template)
  rank_state_seconds <- proc.time()[["elapsed"]] - rank_state_started
  species_prefix_started <- proc.time()[["elapsed"]]
  prefixes <- prepare_benchmark_species_prefixes(
    bundle, sp, rank_state$rank_table$cell, plan$union$keep_n, n_cells
  )
  species_prefix_seconds <- proc.time()[["elapsed"]] - species_prefix_started
  counters <- c(
    rank_raster_reads = 1, kernel_loads = 1, compact_entries_touched = 0,
    dense_cells_avoided = 0, graph_edge_pairs = 0, geometry_endpoints = 0
  )
  state_initialize_seconds <- 0
  checkpoint_restore_seconds <- 0
  stage_zero_lookup_seconds <- 0
  stage_zero_checkpoint_metrics <- empty_benchmark_checkpoint_metrics()
  resumed <- file.exists(config$paths$checkpoint)

  if (resumed) {
    restore_started <- proc.time()[["elapsed"]]
    checkpoint <- validate_benchmark_checkpoint(
      readRDS(config$paths$checkpoint), metadata, plan
    )
    state <- restore_benchmark_state(checkpoint$state, bundle, sp, n_cells)
    completed <- checkpoint$completed_union_target
    counters <- checkpoint$counters
    assert(is.numeric(counters) && all(is.finite(counters) & counters >= 0),
           "Benchmark checkpoint runtime counters must be finite and non-negative.")
    checkpoint_restore_seconds <- proc.time()[["elapsed"]] - restore_started
    runtime_log_event(
      "benchmark_resume",
      completed_union_target = completed,
      next_union_target = completed + 1L
    )
  } else {
    initialize_started <- proc.time()[["elapsed"]]
    state <- stage72_initialize_incremental_state(bundle, sp, n_cells)
    state_initialize_seconds <- proc.time()[["elapsed"]] - initialize_started
    completed <- 0L
    stage_zero_lookup_started <- proc.time()[["elapsed"]]
    write_benchmark_lookup(
      config, plan$union[union_target_id == 0L, keep_n][[1L]], state$patch_table
    )
    stage_zero_lookup_seconds <-
      proc.time()[["elapsed"]] - stage_zero_lookup_started
    stage_zero_checkpoint_metrics <- save_benchmark_checkpoint(
      state, metadata, 0L, counters, config$paths$checkpoint,
      config$checkpoint_interval_targets, reason = "initial"
    )
  }

  setup_total_seconds <- proc.time()[["elapsed"]] - setup_started
  setup_accounted_seconds <- sum(c(
    rank_state_seconds, species_prefix_seconds, state_initialize_seconds,
    checkpoint_restore_seconds, stage_zero_lookup_seconds,
    stage_zero_checkpoint_metrics[["total_seconds"]]
  ))
  setup_residual_seconds <- max(0, setup_total_seconds - setup_accounted_seconds)
  setup_diagnostics <- c(
    rank_state_seconds = rank_state_seconds,
    species_prefix_seconds = species_prefix_seconds,
    state_initialize_seconds = state_initialize_seconds,
    checkpoint_restore_seconds = checkpoint_restore_seconds,
    stage_zero_lookup_seconds = stage_zero_lookup_seconds,
    stage_zero_checkpoint_seconds = stage_zero_checkpoint_metrics[["total_seconds"]],
    accounted_seconds = setup_accounted_seconds,
    residual_seconds = setup_residual_seconds,
    total_seconds = setup_total_seconds
  )
  assert(all(is.finite(setup_diagnostics) & setup_diagnostics >= 0),
         "Benchmark reconstruction-setup diagnostics must be finite and non-negative.")
  runtime_log_event(
    "benchmark_reconstruction_setup",
    rank_state_seconds = sprintf("%.3f", rank_state_seconds),
    species_prefix_seconds = sprintf("%.3f", species_prefix_seconds),
    state_initialize_seconds = sprintf("%.3f", state_initialize_seconds),
    checkpoint_restore_seconds = sprintf("%.3f", checkpoint_restore_seconds),
    stage_zero_lookup_seconds = sprintf("%.3f", stage_zero_lookup_seconds),
    stage_zero_checkpoint_seconds = sprintf(
      "%.3f", stage_zero_checkpoint_metrics[["total_seconds"]]
    ),
    accounted_seconds = sprintf("%.3f", setup_accounted_seconds),
    residual_seconds = sprintf("%.3f", setup_residual_seconds),
    total_seconds = sprintf("%.3f", setup_total_seconds),
    resumed = resumed,
    completed_union_target = as.integer(completed)
  )

  previous_keep <- plan$union[union_target_id == completed, keep_n][[1L]]
  previous_counts <- prefixes$counts[, as.character(completed)]
  remaining <- plan$union[union_target_id > completed]
  for (i in seq_len(nrow(remaining))) {
    target <- remaining[i]
    started <- proc.time()[["elapsed"]]
    target_prepare_started <- proc.time()[["elapsed"]]
    current_id <- as.integer(target$union_target_id)
    keep_n <- as.integer(target$keep_n)
    removed_global <- if (previous_keep > keep_n) {
      rank_state$rank_table$cell[seq.int(keep_n + 1L, previous_keep)]
    } else integer()
    if (length(removed_global)) state$alive_counts[removed_global] <- 0L
    current_counts <- prefixes$counts[, as.character(current_id)]
    removed_by_species <- list()
    for (j in seq_len(nrow(sp))) {
      if (previous_counts[[j]] > current_counts[[j]]) {
        scientific_name <- sp$scientificName[[j]]
        removed_by_species[[scientific_name]] <- prefixes$cells[[scientific_name]][
          seq.int(current_counts[[j]] + 1L, previous_counts[[j]])
        ]
      }
    }
    target_prepare_seconds <- proc.time()[["elapsed"]] - target_prepare_started
    updated <- stage72_apply_incremental_stage(
      state, removed_by_species,
      stats::setNames(sp$species, sp$scientificName), current_id, template,
      build_luts = FALSE
    )
    state <- updated$state
    write_benchmark_lookup(config, keep_n, state$patch_table)
    workload_value <- function(workload, field) {
      if (field %in% names(workload)) as.numeric(workload[[field]]) else 0
    }
    counters[["compact_entries_touched"]] <- counters[["compact_entries_touched"]] +
      workload_value(updated$runtime_workload$fragmentation, "changed_index_entries")
    counters[["dense_cells_avoided"]] <- counters[["dense_cells_avoided"]] +
      workload_value(updated$runtime_workload$fragmentation, "dense_cells_avoided")
    counters[["graph_edge_pairs"]] <- counters[["graph_edge_pairs"]] +
      workload_value(updated$runtime_workload$distance, "graph_edge_pairs_evaluated")
    counters[["geometry_endpoints"]] <- counters[["geometry_endpoints"]] +
      as.numeric(updated$counts[["candidate_geometry_patches"]] %||% 0)
    incremental_seconds <- sum(updated$timing)
    detailed_log_started <- proc.time()[["elapsed"]]
    stage72_log_incremental_timing(current_id, updated$timing, updated$counts)
    detailed_log_seconds <- proc.time()[["elapsed"]] - detailed_log_started

    consumers <- plan$consumers[union_target_id == current_id]
    # Full recovery checkpoints (validation + disk write) are the dominant
    # per-target cost, and re-running them on every target buys no additional
    # safety: a crash between checkpoints re-derives the skipped targets
    # deterministically from the cached prefixes on resume (see
    # validate_benchmark_checkpoint()/restore_benchmark_state()). Skip both the
    # validation and the write except every checkpoint_interval_targets
    # targets, and always checkpoint the final target so a completed run
    # always leaves a checkpoint reflecting its true end state.
    checkpoint_due <- (current_id %% config$checkpoint_interval_targets == 0L) ||
      i == nrow(remaining)
    checkpoint_metrics <- empty_benchmark_checkpoint_metrics()
    if (checkpoint_due) {
      checkpoint_reason <- if (i == nrow(remaining)) "final" else "interval"
      checkpoint_metrics <- save_benchmark_checkpoint(
        state, metadata, current_id, counters,
        config$paths$checkpoint, config$checkpoint_interval_targets,
        reason = checkpoint_reason
      )
    }
    touched_species_count <- length(updated$touched_species)
    fragmentation_seconds <- updated$timing[["fragmentation_seconds"]] %||% 0
    distance_seconds <- sum(updated$timing[
      grep("^distance_.*_seconds$", names(updated$timing))
    ])
    consumer_count <- nrow(consumers)
    removed_rank_cells <- length(removed_global)
    previous_keep <- keep_n
    previous_counts <- current_counts
    completed <- current_id
    cleanup_started <- proc.time()[["elapsed"]]
    rm(updated)
    gc(FALSE)
    cleanup_gc_seconds <- proc.time()[["elapsed"]] - cleanup_started
    accounted_seconds <- sum(c(
      target_prepare_seconds, incremental_seconds, detailed_log_seconds,
      checkpoint_metrics[["prepare_seconds"]],
      checkpoint_metrics[["write_seconds"]], cleanup_gc_seconds
    ))
    total_seconds <- proc.time()[["elapsed"]] - started
    residual_seconds <- max(0, total_seconds - accounted_seconds)
    target_diagnostics <- c(
      target_prepare_seconds = target_prepare_seconds,
      incremental_seconds = incremental_seconds,
      detailed_log_seconds = detailed_log_seconds,
      checkpoint_prepare_seconds = checkpoint_metrics[["prepare_seconds"]],
      checkpoint_write_seconds = checkpoint_metrics[["write_seconds"]],
      cleanup_gc_seconds = cleanup_gc_seconds,
      accounted_seconds = accounted_seconds,
      residual_seconds = residual_seconds,
      total_seconds = total_seconds
    )
    assert(all(is.finite(target_diagnostics) & target_diagnostics >= 0),
           "Benchmark union-target diagnostics must be finite and non-negative.")
    runtime_log_event(
      "benchmark_union_target",
      console = FALSE,
      union_target_id = current_id,
      keep_n = keep_n,
      removed_rank_cells = removed_rank_cells,
      consumer_count = consumer_count,
      touched_species = touched_species_count,
      fragmentation_seconds = sprintf("%.3f", fragmentation_seconds),
      distance_seconds = sprintf("%.3f", distance_seconds),
      target_prepare_seconds = sprintf("%.3f", target_prepare_seconds),
      incremental_seconds = sprintf("%.3f", incremental_seconds),
      detailed_log_seconds = sprintf("%.3f", detailed_log_seconds),
      checkpoint_prepare_seconds = sprintf(
        "%.3f", checkpoint_metrics[["prepare_seconds"]]
      ),
      checkpoint_write_seconds = sprintf(
        "%.3f", checkpoint_metrics[["write_seconds"]]
      ),
      checkpoint_seconds = sprintf(
        "%.3f", checkpoint_metrics[["total_seconds"]]
      ),
      checkpoint_written = checkpoint_due,
      cleanup_gc_seconds = sprintf("%.3f", cleanup_gc_seconds),
      accounted_seconds = sprintf("%.3f", accounted_seconds),
      residual_seconds = sprintf("%.3f", residual_seconds),
      total_seconds = sprintf("%.3f", total_seconds)
    )
  }
  runtime_log_event(
    "benchmark_reconstruction_workload",
    rank_raster_reads = counters[["rank_raster_reads"]],
    kernel_loads = counters[["kernel_loads"]],
    compact_entries_touched = counters[["compact_entries_touched"]],
    dense_cells_avoided = counters[["dense_cells_avoided"]],
    graph_edge_pairs = counters[["graph_edge_pairs"]],
    geometry_endpoints = counters[["geometry_endpoints"]]
  )
  finalize_benchmark_library(config, plan)
  list(counters = counters, completed_union_target = completed)
}
