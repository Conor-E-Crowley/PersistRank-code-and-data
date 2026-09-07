# Deterministic compiled kernels for measured Stage 6/7.2 hot paths.
#
# The source-aware Rcpp cache is supplied by the caller. Contract v11 and its
# startup self-tests must pass before run output changes. Missing kernels are
# fatal; exact R oracles live only in the Stage 6 test helper.

stage6_runtime_kernel_contract_expected <- function() {
  "stage6_runtime_kernels_v11"
}

stage6_runtime_kernels_available <- function() {
  all(vapply(
    c(
      "stage6_runtime_kernel_contract",
      "stage6_build_cell_species_index_cpp",
      "stage6_rebuild_species_pus_cpp",
      "stage6_score_frontier_cpp",
      "stage6_score_frontier_sparse_cpp",
      "stage6_splice_patch_index_cpp",
      "stage6_label_patch_components_cpp",
      "stage6_analyze_fragmentation_components_cpp",
      "stage72_analyze_fragmentation_components_cpp",
      "stage6_extract_candidate_patch_cells_cpp",
      "stage6_extract_distance_recheck_edges_cpp",
      "stage6_build_provisional_fragment_graph_cpp",
      "stage6_filter_distance_graph_cpp",
      "stage6_canonicalize_rank_lut_cpp"
    ),
    exists,
    logical(1L),
    mode = "function",
    inherits = TRUE
  ))
}

validate_stage6_runtime_kernels <- function() {
  if (!stage6_runtime_kernels_available()) {
    stop("The compiled Stage 6 runtime kernels are not loaded.", call. = FALSE)
  }
  contract <- stage6_runtime_kernel_contract()
  if (!identical(as.character(contract), stage6_runtime_kernel_contract_expected())) {
    stop(
      "The compiled Stage 6 runtime-kernel contract is incompatible. Expected '",
      stage6_runtime_kernel_contract_expected(), "' but received '",
      as.character(contract), "'.",
      call. = FALSE
    )
  }

  frontier_state <- stage6_score_frontier_cpp(
    frontier_cells = 1:3,
    patch_ids_by_species = list(
      c(1L, 1L, NA_integer_),
      c(1L, NA_integer_, 1L)
    ),
    scores_by_species = list(log(0.2), log(0.3)),
    return_diagnostics = TRUE
  )
  frontier <- -(frontier_state$max_log_term + log(frontier_state$sum_exp_terms))
  expected_frontier <- c(-log(0.2 + 0.3), -log(0.2), -log(0.3))
  if (!isTRUE(all.equal(frontier, expected_frontier, tolerance = 1e-15))) {
    stop("The compiled Stage 6 frontier kernel failed its deterministic self-test.", call. = FALSE)
  }
  if (!identical(
    names(frontier_state$diagnostic_counts),
    c(
      "score_calls", "frontier_cells", "active_species",
      "species_frontier_probes", "sparse_membership_probes",
      "dense_probes_avoided", "finite_contributions",
      "scored_frontier_cells"
    )
  ) || !identical(
    names(frontier_state$diagnostic_timing_seconds),
    c(
      "kernel_scan_seconds", "kernel_accumulation_seconds",
      "kernel_output_seconds"
    )
  ) || any(!is.finite(c(
    frontier_state$diagnostic_counts,
    frontier_state$diagnostic_timing_seconds
  ))) || any(c(
    frontier_state$diagnostic_counts,
    frontier_state$diagnostic_timing_seconds
  ) < 0)) {
    stop("The compiled Stage 6 frontier diagnostics failed their contract self-test.", call. = FALSE)
  }

  cell_species <- stage6_build_cell_species_index_cpp(
    cells_by_species = list(c(1L, 2L), c(1L, 3L)),
    n_cells = 3L,
    use_raw_species_ids = TRUE
  )
  sparse_frontier_state <- stage6_score_frontier_sparse_cpp(
    frontier_cells = 1:3,
    cell_offsets = cell_species$offsets,
    species_ids_input = cell_species$species_ids,
    canonical_species_count = 2L,
    active_species_ids = 1:2,
    patch_ids_by_species = list(
      c(1L, 1L, NA_integer_),
      c(1L, NA_integer_, 1L)
    ),
    scores_by_species = list(log(0.2), log(0.3)),
    return_diagnostics = TRUE
  )
  if (!identical(
    sparse_frontier_state[c("max_log_term", "sum_exp_terms", "has_score")],
    frontier_state[c("max_log_term", "sum_exp_terms", "has_score")]
  ) || !identical(
    names(sparse_frontier_state$diagnostic_counts),
    names(frontier_state$diagnostic_counts)
  )) {
    stop("The compiled sparse frontier kernel failed its deterministic self-test.", call. = FALSE)
  }

  repair <- stage6_rebuild_species_pus_cpp(
    pu_graphs = list(list(
      pu_id = 2L,
      id2patch = 1:3,
      row_ptr = c(0L, 1L, 3L, 4L),
      col_idx = c(2L, 1L, 3L, 2L)
    )),
    dead_patch_ids_by_pu = list(2L),
    patch_area_by_pu = list(c(2, 0, 3)),
    pu_area_threshold = 1,
    next_available_pu_id = 7L
  )
  if (!identical(repair$next_available_pu_id, 8L) ||
      !identical(
        lapply(repair$updates[[1L]]$surviving_pu_graphs, `[[`, "id2patch"),
        list(1L, 3L)
      ) ||
      !identical(
        vapply(
          repair$updates[[1L]]$surviving_pu_graphs,
          `[[`,
          integer(1L),
          "pu_id"
        ),
        c(2L, 8L)
      )) {
    stop("The compiled Stage 6 PU-repair kernel failed its deterministic self-test.", call. = FALSE)
  }

  splice <- stage6_splice_patch_index_cpp(
    old_pid = c(1L, 1L, 2L),
    old_cell = c(1L, 3L, 4L),
    expected_patch_ids_input = c(1L, 2L),
    pre_update_patch_ids_input = 1L,
    replacement_origin_ids_input = 1L,
    replacement_origin = c(1L, 1L),
    replacement_pid = c(1L, 2L),
    replacement_cell = c(1L, 3L)
  )
  if (!identical(splice$pid, c(1L, 2L)) ||
      !identical(splice$cell, c(1L, 3L)) ||
      !identical(splice$removed_cells, 4L)) {
    stop("The compiled Stage 6 compact-index kernel failed its deterministic self-test.", call. = FALSE)
  }

  components <- stage6_label_patch_components_cpp(
    patch_cells_input = c(4L, 1L, 2L),
    rook_start = c(1L, 2L, 4L, 5L),
    rook_end = c(1L, 3L, 4L, 5L),
    rook_to = c(2L, 1L, 3L, 2L, 4L),
    n_cells = 4L
  )
  if (!identical(components$cells, c(1L, 2L, 4L)) ||
      !identical(components$component_id, c(1L, 1L, 2L)) ||
      !identical(components$component_count, 2L)) {
    stop("The compiled Stage 6 component kernel failed its deterministic self-test.", call. = FALSE)
  }

  fragmentation <- stage6_analyze_fragmentation_components_cpp(
    index_pid = c(1L, 1L, 1L, 2L),
    index_cell = c(1L, 2L, 4L, 3L),
    changed_patch_ids = c(1L, 2L),
    current_patch_id_by_cell = c(1L, 1L, 2L, 1L),
    cell_area_by_cell = c(1, 2, 4, 8),
    rook_start = c(1L, 2L, 4L, 5L),
    rook_end = c(1L, 3L, 4L, 5L),
    rook_to = c(2L, 1L, 3L, 2L, 4L),
    n_cells = 4L,
    profile = FALSE
  )
  if (!identical(fragmentation$origin_patch_ids, c(1L, 2L)) ||
      !identical(fragmentation$origin_offsets, c(0L, 3L, 4L)) ||
      !identical(fragmentation$cells, c(1L, 2L, 4L, 3L)) ||
      !identical(fragmentation$component_id_by_cell, c(1L, 1L, 2L, 1L)) ||
      !identical(fragmentation$component_origin_patch_id, c(1L, 1L, 2L)) ||
      !identical(fragmentation$component_id, c(1L, 2L, 1L)) ||
      !identical(fragmentation$component_area_km2, c(3, 8, 4)) ||
      !identical(fragmentation$component_n_cells, c(2L, 1L, 1L)) ||
      !identical(fragmentation$component_first_cell, c(1L, 4L, 3L)) ||
      !identical(fragmentation$component_count_by_origin, c(2L, 1L))) {
    stop(
      "The compiled Stage 6 fragmentation-analysis kernel failed its deterministic self-test.",
      call. = FALSE
    )
  }

  compact_fragmentation <- stage72_analyze_fragmentation_components_cpp(
    index_pid = c(1L, 1L, 1L, 2L),
    index_cell = c(1L, 2L, 4L, 3L),
    changed_patch_ids = c(1L, 2L),
    cell_area_by_cell = c(1, 2, 4, 8),
    rook_start = c(1L, 2L, 4L, 5L),
    rook_end = c(1L, 3L, 4L, 5L),
    rook_to = c(2L, 1L, 3L, 2L, 4L),
    n_cells = 4L,
    profile = FALSE
  )
  if (!identical(compact_fragmentation, fragmentation)) {
    stop(
      "The compact Stage 7.2 fragmentation kernel failed its deterministic self-test.",
      call. = FALSE
    )
  }

  candidate_cells <- stage6_extract_candidate_patch_cells_cpp(
    index_pid = c(1L, 1L, 3L, 7L, 7L),
    index_cell = c(2L, 5L, 8L, 11L, 12L),
    candidate_patch_ids_input = c(7L, 1L, 7L)
  )
  if (!identical(candidate_cells$patch_id, c(1L, 1L, 7L, 7L)) ||
      !identical(candidate_cells$cell, c(2L, 5L, 11L, 12L))) {
    stop(
      "The compiled Stage 6 candidate-cell kernel failed its deterministic self-test.",
      call. = FALSE
    )
  }

  standard_recheck_graph <- list(
    id2patch = c(10L, 20L, 30L),
    row_ptr = c(0L, 1L, 3L, 4L),
    col_idx = c(2L, 1L, 3L, 2L)
  )
  lazy_recheck_graph <- list(
    graph_type = "lazy_added_node_overlay",
    id2patch = c(40L, 50L, 60L),
    base_node_count = 2L,
    base_graph = list(
      id2patch = c(40L, 50L),
      row_ptr = c(0L, 1L, 2L),
      col_idx = c(2L, 1L)
    ),
    overlay_edges = data.frame(
      u = c(1L, 1L, 2L),
      v = c(2L, 3L, 3L)
    )
  )
  recheck_edges <- stage6_extract_distance_recheck_edges_cpp(
    affected_graphs = list(standard_recheck_graph, lazy_recheck_graph),
    recheck_patch_ids_input = c(99L, 20L, 40L)
  )
  if (!identical(recheck_edges$patch_u, c(10L, 20L, 40L, 40L)) ||
      !identical(recheck_edges$patch_v, c(20L, 30L, 50L, 60L)) ||
      !identical(
        recheck_edges$candidate_patch_ids,
        c(10L, 20L, 30L, 40L, 50L, 60L, 99L)
      ) ||
      !identical(as.numeric(recheck_edges$unique_recheck_edges), 4)) {
    stop(
      "The compiled Stage 6 distance recheck-edge kernel failed its deterministic self-test.",
      call. = FALSE
    )
  }

  provisional <- stage6_build_provisional_fragment_graph_cpp(
    current_patch_ids = c(21L, 10L, 22L),
    current_origin_patch_ids = c(20L, 10L, 20L),
    old_patch_ids = c(10L, 20L, 30L),
    old_row_ptr = c(0L, 1L, 3L, 4L),
    old_col_idx = c(2L, 1L, 3L, 2L),
    pu_id = 4L
  )
  if (!identical(as.integer(provisional$pu_id), 4L) ||
      !identical(provisional$id2patch, c(21L, 10L, 22L)) ||
      !identical(provisional$row_ptr, c(0L, 2L, 4L, 6L)) ||
      !identical(provisional$col_idx, c(2L, 3L, 1L, 3L, 1L, 2L)) ||
      !identical(as.numeric(provisional$candidate_undirected_edges), 3) ||
      !identical(as.numeric(provisional$unique_undirected_edges), 3) ||
      !identical(as.numeric(provisional$output_adjacency_entries), 6)) {
    stop(
      "The compiled Stage 6 provisional-graph kernel failed its deterministic self-test.",
      call. = FALSE
    )
  }

  filtered <- stage6_filter_distance_graph_cpp(
    id2patch = c(10L, 20L, 30L, 40L),
    base_node_count = 3L,
    base_row_ptr = c(0L, 1L, 3L, 4L),
    base_col_idx = c(2L, 1L, 3L, 2L),
    overlay_u = c(3L, 4L),
    overlay_v = c(4L, 1L),
    invalid_patch_u = 20L,
    invalid_patch_v = 30L,
    pu_id = 7L
  )
  if (!identical(as.integer(filtered$pu_id), 7L) ||
      !identical(filtered$id2patch, c(10L, 20L, 30L, 40L)) ||
      !identical(filtered$row_ptr, c(0L, 2L, 3L, 4L, 6L)) ||
      !identical(filtered$col_idx, c(2L, 4L, 1L, 4L, 1L, 3L)) ||
      !identical(as.numeric(filtered$invalid_undirected_edges), 1) ||
      !identical(as.numeric(filtered$input_adjacency_entries), 8) ||
      !identical(as.numeric(filtered$output_adjacency_entries), 6)) {
    stop(
      "The compiled Stage 6 distance-filter kernel failed its deterministic self-test.",
      call. = FALSE
    )
  }

  canonical <- stage6_canonicalize_rank_lut_cpp(
    patch_ids = c(9L, 3L, 7L),
    pu_ids = c(8L, 4L, 8L),
    index_pid = c(3L, 7L, 7L, 9L),
    index_cell = c(2L, 5L, 6L, 9L)
  )
  if (!identical(canonical$row_order, c(2L, 3L, 1L)) ||
      !identical(canonical$patch_id, 1:3) ||
      !identical(canonical$pu_id, c(1L, 2L, 2L)) ||
      !identical(canonical$patch_n_cells, c(1L, 2L, 1L))) {
    stop(
      "The compiled Stage 7.2 canonical-LUT kernel failed its deterministic self-test.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

load_stage6_runtime_kernels <- function(cache_dir) {
  cache_dir <- validate_path_param(cache_dir, "cache_dir")

  # Rcpp wrappers can remain in a long-lived interactive session after their
  # shared library has been unloaded. Existence alone is therefore not proof
  # that the native symbols are callable. Reuse only a validated kernel set;
  # otherwise reload it from the source-aware cache below.
  if (stage6_runtime_kernels_available()) {
    kernels_valid <- tryCatch(
      {
        validate_stage6_runtime_kernels()
        TRUE
      },
      error = function(error) FALSE
    )
    if (kernels_valid) return(invisible(TRUE))
  }
  if (!requireNamespace("Rcpp", quietly = TRUE)) {
    stop(
      "Active Stage 6 prioritization requires Rcpp and a working C++ toolchain. ",
      "Install Rcpp before starting or resuming the pipeline.",
      call. = FALSE
    )
  }

  ensure_writable_dir(cache_dir, "Stage 6 Rcpp cache directory")
  source_path <- file.path("src", "priority_stage6_kernels.cpp")
  need_file(source_path, "Stage 6 runtime-kernel source")

  load_from_source <- function(rebuild) {
    Rcpp::sourceCpp(
      file = source_path,
      cacheDir = cache_dir,
      rebuild = rebuild,
      showOutput = FALSE,
      verbose = FALSE,
      env = globalenv()
    )
  }

  load_error <- tryCatch(
    {
      load_from_source(FALSE)
      validate_stage6_runtime_kernels()
      NULL
    },
    error = identity
  )

  # A cached shared library can itself be stale after an interrupted or
  # upgraded session. Rebuild once before reporting an actionable failure.
  if (!is.null(load_error)) {
    load_error <- tryCatch(
      {
        load_from_source(TRUE)
        validate_stage6_runtime_kernels()
        NULL
      },
      error = identity
    )
  }
  if (!is.null(load_error)) {
    stop(
      "Active Stage 6 prioritization requires the deterministic runtime ",
      "kernels, but they could not be compiled or loaded. Verify the Rcpp ",
      "compiler toolchain. Original error: ", conditionMessage(load_error),
      call. = FALSE
    )
  }
  invisible(TRUE)
}
