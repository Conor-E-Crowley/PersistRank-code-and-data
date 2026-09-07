# Candidate geometry and targeted WGS84/S2 distance predicates.
#
# Candidate cells come from flat patch indexes, are polygonized through the
# existing single-crop Terra path, and retain canonical patch-ID order.
# Stage 6 retains its established candidate-matrix predicate. Stage 7.2 uses
# the same validated polygons and threshold semantics but evaluates only the
# requested graph-edge endpoint pairs; neither path uses a planar approximation
# for WGS84 data.


# Read the species-specific dispersal threshold.

get_species_dispersal_threshold_km <- function(
  species_name,        # character scalar naming the focal species
  species_params       # data.table containing species-level parameters
) {
  # Read this species' dispersal threshold in kilometers.
  dispersal_threshold_km <- species_params[
    species == species_name
  ]$dispersal_distance_km[1L]

  # Stop if the dispersal threshold is missing, non-finite, or non-positive.
  if (!is.finite(dispersal_threshold_km) || dispersal_threshold_km <= 0) {
    stop("Invalid dispersal_distance_km for species: ", species_name)
  }

  # Return the valid dispersal threshold.
  dispersal_threshold_km
}



# Identify lazy added-node overlay PU graphs.

is_lazy_added_node_overlay_graph <- function(
  pu_graph              # one PU graph object
) {
  # Return TRUE only for list objects explicitly marked as lazy overlay graphs.
  is.list(pu_graph) &&
    !is.null(pu_graph$graph_type) &&
    identical(pu_graph$graph_type, "lazy_added_node_overlay")
}



# Get effective patch IDs from a standard or lazy graph.

get_pu_graph_id2patch <- function(
  pu_graph              # standard CSR graph or lazy overlay graph
) {
  # Lazy overlay graphs already expose the full effective id2patch vector.
  if (is_lazy_added_node_overlay_graph(pu_graph)) {
    return(as.integer(pu_graph$id2patch))
  }

  # Standard CSR graphs also store id2patch directly.
  as.integer(pu_graph$id2patch)
}
# Copy only requested compact-index segments with the validated runtime kernel.
# The exact R oracle is confined to the test helper.
extract_candidate_patch_cells <- function(
  species_name,
  species_patch_index,
  candidate_patch_ids
) {
  candidate_patch_ids_numeric <- suppressWarnings(as.numeric(candidate_patch_ids))
  if (
    anyNA(candidate_patch_ids_numeric) ||
      any(!is.finite(candidate_patch_ids_numeric)) ||
      any(candidate_patch_ids_numeric < 1) ||
      any(candidate_patch_ids_numeric != floor(candidate_patch_ids_numeric))
  ) {
    stop(
      "Candidate distance patch IDs must be positive finite integers for species: ",
      species_name
    )
  }
  candidate_patch_ids <- sort(unique(as.integer(candidate_patch_ids_numeric)))
  if (!length(candidate_patch_ids)) {
    return(data.table::data.table(
      patch_id = integer(),
      cell = integer()
    ))
  }
  if (
    !is.list(species_patch_index) ||
      !all(c("pid", "cell") %in% names(species_patch_index)) ||
      length(species_patch_index$pid) != length(species_patch_index$cell) ||
      !is.integer(species_patch_index$pid) ||
      !is.integer(species_patch_index$cell) ||
      !length(species_patch_index$pid)
  ) {
    stop("Invalid compact patch-to-cell index for species: ", species_name)
  }
  if (!exists(
    "stage6_extract_candidate_patch_cells_cpp",
    mode = "function",
    inherits = TRUE
  )) {
    stop(
      "The compiled Stage 6 candidate-cell kernel is not loaded for species: ",
      species_name,
      call. = FALSE
    )
  }

  extracted <- tryCatch(
    stage6_extract_candidate_patch_cells_cpp(
      index_pid = species_patch_index$pid,
      index_cell = species_patch_index$cell,
      candidate_patch_ids_input = candidate_patch_ids
    ),
    error = function(error) {
      stop(
        "Could not extract candidate distance cells for species '",
        species_name,
        "': ",
        conditionMessage(error),
        call. = FALSE
      )
    }
  )
  data.table::data.table(
    patch_id = extracted$patch_id,
    cell = extracted$cell
  )
}



# Extract only graph edges incident to recheck patches in compiled code.
extract_distance_recheck_edges <- function(
  species_name,
  recheck_patch_ids_for_species,
  affected_pu_ids,
  pu_graphs_by_key,
  return_diagnostics = FALSE
) {
  recheck_numeric <- suppressWarnings(as.numeric(
    recheck_patch_ids_for_species
  ))
  if (
    anyNA(recheck_numeric) || any(!is.finite(recheck_numeric)) ||
      any(recheck_numeric < 1) || any(recheck_numeric != floor(recheck_numeric))
  ) {
    stop(
      "Distance-recheck patch IDs must be positive finite integers for species: ",
      species_name,
      call. = FALSE
    )
  }
  recheck_patch_ids_for_species <- sort(unique(as.integer(recheck_numeric)))

  if (!exists(
    "stage6_extract_distance_recheck_edges_cpp",
    mode = "function",
    inherits = TRUE
  )) {
    stop(
      "The compiled distance recheck-edge kernel is not loaded.",
      call. = FALSE
    )
  }

  graph_keys <- paste0(species_name, "|", affected_pu_ids)
  affected_graphs <- priority_graph_snapshot(
    pu_graphs_by_key,
    graph_keys,
    require_all = FALSE
  )
  kernel_result <- stage6_extract_distance_recheck_edges_cpp(
    affected_graphs = unname(affected_graphs),
    recheck_patch_ids_input = recheck_patch_ids_for_species
  )

  edges <- data.table::data.table(
    patch_u = as.integer(kernel_result$patch_u),
    patch_v = as.integer(kernel_result$patch_v)
  )
  edges[, edge_key := paste0(patch_u, "|", patch_v)]
  result <- list(
    edges = edges,
    candidate_patch_ids = as.integer(kernel_result$candidate_patch_ids)
  )
  if (!isTRUE(return_diagnostics)) return(result)

  result$diagnostics <- c(
    graphs_inspected = as.numeric(kernel_result$graphs_inspected),
    graph_nodes_scanned = as.numeric(kernel_result$graph_nodes_scanned),
    recheck_rows = as.numeric(kernel_result$recheck_rows),
    csr_adjacency_entries_visited = as.numeric(
      kernel_result$csr_adjacency_entries_visited
    ),
    overlay_edges_scanned = as.numeric(kernel_result$overlay_edges_scanned),
    unique_recheck_edges = as.numeric(kernel_result$unique_recheck_edges)
  )
  result
}



# Store sparse Boolean dispersal decisions; never build a dense distance matrix.

empty_distance_predicate_profile <- function() {
  stats::setNames(numeric(15L), c(
    "predicate_calls", "recheck_patches", "candidate_patches",
    "candidate_pair_upper_bound", "graph_edges",
    "directed_predicate_hits", "unique_spatial_keys",
    "self_spatial_keys", "valid_graph_edges", "invalid_graph_edges",
    "non_graph_spatial_keys", "predicate_spatial_seconds",
    "predicate_mapping_seconds", "residual_seconds", "total_seconds"
  ))
}

build_distance_predicate_lookup <- function(
  candidate_patch_polygons,
  dispersal_threshold_km,
  recheck_patch_ids,
  recheck_edges,
  predicate_strategy = c("matrix", "edge_pairs"),
  return_diagnostics = FALSE
) {
  predicate_strategy <- match.arg(predicate_strategy)
  predicate_started <- proc.time()[["elapsed"]]
  predicate_timing <- c(
    predicate_spatial_seconds = 0,
    predicate_mapping_seconds = 0
  )
  predicate_profile <- empty_distance_predicate_profile()
  predicate_profile[["predicate_calls"]] <- 1
  finish_predicate <- function(lookup) {
    if (!isTRUE(return_diagnostics)) return(lookup)
    predicate_profile[["predicate_spatial_seconds"]] <-
      predicate_timing[["predicate_spatial_seconds"]]
    predicate_profile[["predicate_mapping_seconds"]] <-
      predicate_timing[["predicate_mapping_seconds"]]
    total_seconds <- proc.time()[["elapsed"]] - predicate_started
    predicate_profile[["residual_seconds"]] <- max(
      0, total_seconds - sum(predicate_timing)
    )
    predicate_profile[["total_seconds"]] <-
      sum(predicate_timing) + predicate_profile[["residual_seconds"]]
    list(
      lookup = lookup,
      timing_seconds = predicate_timing,
      profile = predicate_profile
    )
  }
  recheck_edges <- data.table::as.data.table(recheck_edges)
  if (!all(c("patch_u", "patch_v", "edge_key") %in% names(recheck_edges))) {
    stop("recheck_edges must contain patch_u, patch_v, and edge_key.", call. = FALSE)
  }
  patch_ids <- as.integer(candidate_patch_polygons$patch_id)
  if (anyDuplicated(patch_ids)) {
    stop("candidate_patch_polygons contains duplicated patch_id values.", call. = FALSE)
  }
  recheck_patch_ids <- intersect(
    sort(unique(as.integer(recheck_patch_ids))),
    patch_ids
  )
  predicate_profile[["recheck_patches"]] <- length(recheck_patch_ids)
  predicate_profile[["candidate_patches"]] <- length(patch_ids)
  predicate_profile[["candidate_pair_upper_bound"]] <-
    as.numeric(length(recheck_patch_ids)) * as.numeric(length(patch_ids))
  predicate_profile[["graph_edges"]] <- nrow(recheck_edges)
  if (!nrow(recheck_edges)) {
    return(finish_predicate(list(
      edge_valid_by_key = stats::setNames(logical(), character()),
      patch_u = integer(), patch_v = integer(), edge_valid = logical()
    )))
  }
  missing_endpoints <- setdiff(
    unique(c(recheck_edges$patch_u, recheck_edges$patch_v)),
    patch_ids
  )
  if (length(missing_endpoints)) {
    stop(
      "Distance predicate geometry is missing graph-edge endpoint patch IDs: ",
      paste(utils::head(missing_endpoints, 10L), collapse = ", "),
      call. = FALSE
    )
  }

  subphase_started <- proc.time()[["elapsed"]]
  if (identical(predicate_strategy, "edge_pairs")) {
    endpoint_u <- match(recheck_edges$patch_u, patch_ids)
    endpoint_v <- match(recheck_edges$patch_v, patch_ids)
    if (sf::st_is_longlat(candidate_patch_polygons)) {
      if (!requireNamespace("s2", quietly = TRUE)) {
        stop("The edge-pair WGS84 predicate requires the 's2' package.")
      }
      geographies <- sf::st_as_s2(candidate_patch_polygons)
      # The vectorized pairwise predicate avoids the candidate matrix. The
      # prepared variant is slower for this one-query-per-edge workload.
      edge_valid <- s2::s2_dwithin(
        geographies[endpoint_u],
        geographies[endpoint_v],
        distance = dispersal_threshold_km * 1000
      )
    } else {
      edge_distance <- sf::st_distance(
        candidate_patch_polygons[endpoint_u, ],
        candidate_patch_polygons[endpoint_v, ],
        by_element = TRUE
      )
      threshold <- units::set_units(
        dispersal_threshold_km * 1000,
        "m"
      )
      edge_valid <- edge_distance <= threshold
    }
  } else {
    x_rows <- match(recheck_patch_ids, patch_ids)
    within_rows <- sf::st_is_within_distance(
      x = candidate_patch_polygons[x_rows, ],
      y = candidate_patch_polygons,
      dist = dispersal_threshold_km * 1000,
      sparse = TRUE
    )
  }
  predicate_timing[["predicate_spatial_seconds"]] <-
    proc.time()[["elapsed"]] - subphase_started
  subphase_started <- proc.time()[["elapsed"]]
  if (identical(predicate_strategy, "edge_pairs")) {
    edge_valid <- as.logical(edge_valid)
    true_edge_keys <- recheck_edges$edge_key[edge_valid]
    predicate_profile[["directed_predicate_hits"]] <- sum(edge_valid)
    predicate_profile[["unique_spatial_keys"]] <- sum(edge_valid)
    predicate_profile[["self_spatial_keys"]] <- 0
    predicate_profile[["non_graph_spatial_keys"]] <- 0
  } else {
    true_edge_keys <- unique(unlist(lapply(seq_along(within_rows), function(i) {
      if (!length(within_rows[[i]])) return(character())
      patch_from <- rep.int(recheck_patch_ids[[i]], length(within_rows[[i]]))
      patch_to <- patch_ids[within_rows[[i]]]
      paste0(pmin.int(patch_from, patch_to), "|", pmax.int(patch_from, patch_to))
    }), use.names = FALSE))
    edge_valid <- recheck_edges$edge_key %in% true_edge_keys
    self_edge_keys <- intersect(
      true_edge_keys,
      paste0(recheck_patch_ids, "|", recheck_patch_ids)
    )
    predicate_profile[["directed_predicate_hits"]] <- sum(lengths(within_rows))
    predicate_profile[["unique_spatial_keys"]] <- length(true_edge_keys)
    predicate_profile[["self_spatial_keys"]] <- length(self_edge_keys)
    predicate_profile[["non_graph_spatial_keys"]] <-
      length(setdiff(true_edge_keys, recheck_edges$edge_key))
  }
  predicate_profile[["valid_graph_edges"]] <- sum(edge_valid)
  predicate_profile[["invalid_graph_edges"]] <- sum(!edge_valid)
  predicate_timing[["predicate_mapping_seconds"]] <-
    proc.time()[["elapsed"]] - subphase_started
  finish_predicate(list(
    edge_valid_by_key = stats::setNames(edge_valid, recheck_edges$edge_key),
    patch_u = as.integer(recheck_edges$patch_u),
    patch_v = as.integer(recheck_edges$patch_v),
    edge_valid = as.logical(edge_valid)
  ))
}


prepare_distance_template_geometry <- function(template_raster) {
  list(
    raster = terra::rast(
      terra::ext(template_raster),
      nrows = terra::nrow(template_raster),
      ncols = terra::ncol(template_raster),
      crs = terra::crs(template_raster, proj = TRUE)
    ),
    resolution = terra::res(template_raster),
    crs = terra::crs(template_raster, proj = TRUE),
    ncell = terra::ncell(template_raster)
  )
}



# Build local patch polygons for one species.
#
# This helper polygonizes only a local crop around candidate patches instead
# of polygonizing the full species raster.
#

validate_and_repair_patch_polygons <- function(
  patch_polygons_sf,
  species_name,
  validity_fn = sf::st_is_valid,
  repair_fn = sf::st_make_valid
) {
  validation_started <- proc.time()[["elapsed"]]
  validity <- validity_fn(patch_polygons_sf)
  validation_seconds <- proc.time()[["elapsed"]] - validation_started
  if (!is.logical(validity) || length(validity) != nrow(patch_polygons_sf)) {
    stop(
      "Geometry validation returned an invalid result for species: ",
      species_name,
      call. = FALSE
    )
  }
  invalid <- is.na(validity) | !validity
  invalid_count <- sum(invalid)
  repair_seconds <- 0
  repair_invoked <- invalid_count > 0L

  if (repair_invoked) {
    repair_started <- proc.time()[["elapsed"]]
    patch_polygons_sf <- repair_fn(patch_polygons_sf)
    repaired_validity <- validity_fn(patch_polygons_sf)
    repair_seconds <- proc.time()[["elapsed"]] - repair_started
    if (!is.logical(repaired_validity) ||
        length(repaired_validity) != nrow(patch_polygons_sf)) {
      stop(
        "Geometry repair validation returned an invalid result for species: ",
        species_name,
        call. = FALSE
      )
    }
    failed <- which(is.na(repaired_validity) | !repaired_validity)
    if (length(failed)) {
      stop(
        "Geometry repair left invalid feature(s) for species '",
        species_name,
        "' at row(s): ",
        paste(utils::head(failed, 20L), collapse = ", "),
        if (length(failed) > 20L) " ..." else "",
        ".",
        call. = FALSE
      )
    }
  }

  list(
    polygons = patch_polygons_sf,
    timing_seconds = c(
      geometry_validation_seconds = validation_seconds,
      geometry_repair_seconds = repair_seconds
    ),
    counts = c(
      geometry_features_checked = as.integer(length(validity)),
      invalid_geometry_features = as.integer(invalid_count),
      geometry_repair_invoked = as.integer(repair_invoked)
    )
  )
}

build_local_patch_polygons <- function(
  species_name,               # focal species, used in validation errors
  template_raster,            # shared single-layer raster template
  species_patch_index,        # compact sorted patch_id -> cell index
  candidate_patch_ids,        # patch IDs needing local geometry
  template_geometry = NULL,
  return_diagnostics = FALSE
) {
  geometry_timing <- stats::setNames(
    numeric(7L),
    c(
      "candidate_index_extract_seconds", "candidate_extent_seconds",
      "local_raster_seconds", "polygonize_seconds",
      "sf_conversion_seconds", "geometry_validation_seconds",
      "geometry_repair_seconds"
    )
  )
  geometry_validation_counts <- c(
    geometry_features_checked = 0L,
    invalid_geometry_features = 0L,
    geometry_repair_invoked = 0L
  )
  finish_geometry <- function(polygons, candidate_cells = 0L, local_cells = 0L) {
    if (!isTRUE(return_diagnostics)) {
      return(polygons)
    }
    list(
      polygons = polygons,
      timing_seconds = geometry_timing,
      counts = c(
        candidate_patches = as.integer(length(candidate_patch_ids)),
        candidate_cells = as.integer(candidate_cells),
        local_raster_cells = as.integer(local_cells),
        expansion_ratio = if (candidate_cells > 0L) local_cells / candidate_cells else 0,
        geometry_validation_counts
      )
    )
  }

  if (is.null(template_raster) || terra::nlyr(template_raster) != 1L) {
    stop("template_raster must be one raster layer.")
  }

  subphase_started <- proc.time()[["elapsed"]]
  # Retrieve candidate cells directly from the compact live index.
  candidate_patch_cells <- extract_candidate_patch_cells(
    species_name = species_name,
    species_patch_index = species_patch_index,
    candidate_patch_ids = candidate_patch_ids
  )
  geometry_timing[["candidate_index_extract_seconds"]] <-
    proc.time()[["elapsed"]] - subphase_started
  subphase_started <- proc.time()[["elapsed"]]

  if (is.null(template_geometry)) {
    template_geometry <- prepare_distance_template_geometry(template_raster)
  }
  if (any(candidate_patch_cells$cell > template_geometry$ncell)) {
    stop("Candidate distance cells exceed the template bounds for species: ", species_name)
  }

  # Return an empty sf object if no candidate cells are present.
  if (!nrow(candidate_patch_cells)) {
    empty_polygons <- sf::st_sf(
      patch_id = integer(),
      geometry = sf::st_sfc(
        crs = sf::st_crs(terra::crs(template_raster, proj = TRUE))
      )
    )
    geometry_timing[["candidate_extent_seconds"]] <-
      proc.time()[["elapsed"]] - subphase_started
    return(finish_geometry(empty_polygons))
  }

  # Read x-y coordinates for candidate cells.
  candidate_xy <- terra::xyFromCell(
    template_raster,
    candidate_patch_cells$cell
  )

  # Read raster resolution.
  raster_resolution <- template_geometry$resolution

  # Build a local crop extent with a one-cell margin.
  local_extent <- terra::ext(
    min(candidate_xy[, 1]) - raster_resolution[1],
    max(candidate_xy[, 1]) + raster_resolution[1],
    min(candidate_xy[, 2]) - raster_resolution[2],
    max(candidate_xy[, 2]) + raster_resolution[2]
  )
  geometry_timing[["candidate_extent_seconds"]] <-
    proc.time()[["elapsed"]] - subphase_started

  # Recreate only the shared template geometry, without its source values, and
  # crop that empty geometry to the candidate extent.
  subphase_started <- proc.time()[["elapsed"]]
  local_raster <- terra::crop(
    template_geometry$raster,
    local_extent,
    snap = "out"
  )

  # Map global candidate-cell centers into the cropped raster.
  local_candidate_cells <- terra::cellFromXY(local_raster, candidate_xy)

  if (anyNA(local_candidate_cells)) {
    stop(
      "Could not map all candidate distance cells into the local raster for species: ",
      species_name
    )
  }

  if (anyDuplicated(local_candidate_cells)) {
    stop("Candidate distance cells map to duplicate local cells for species: ", species_name)
  }

  # Populate only candidate cells. Numeric values avoid Terra's integer
  # value-type dispatch while preserving exact integer-valued patch IDs.
  local_patch_ids <- rep(NA_real_, terra::ncell(local_raster))
  local_patch_ids[local_candidate_cells] <- as.numeric(
    candidate_patch_cells$patch_id
  )
  local_raster <- terra::setValues(local_raster, local_patch_ids)
  local_raster_cell_count <- terra::ncell(local_raster)
  geometry_timing[["local_raster_seconds"]] <-
    proc.time()[["elapsed"]] - subphase_started

  # Polygonize the local candidate-patch raster.
  subphase_started <- proc.time()[["elapsed"]]
  patch_polygons_vect <- terra::as.polygons(
    local_raster,
    dissolve = TRUE,
    na.rm = TRUE
  )
  geometry_timing[["polygonize_seconds"]] <-
    proc.time()[["elapsed"]] - subphase_started

  # Convert terra vector polygons to sf.
  subphase_started <- proc.time()[["elapsed"]]
  patch_polygons_sf <- sf::st_as_sf(patch_polygons_vect)

  # Read the active sf geometry column name.
  geometry_column_name <- attr(patch_polygons_sf, "sf_column")

  # Identify non-geometry attribute columns.
  non_geometry_columns <- setdiff(names(patch_polygons_sf), geometry_column_name)

  # Stop if no patch-ID attribute column was created.
  if (!length(non_geometry_columns)) {
    stop("Could not identify the patch-ID attribute column after polygonization.")
  }

  # Rename the first non-geometry column to patch_id.
  names(patch_polygons_sf)[
    names(patch_polygons_sf) == non_geometry_columns[1L]
  ] <- "patch_id"

  # Coerce polygon patch IDs to integer.
  patch_polygons_sf$patch_id <- as.integer(patch_polygons_sf$patch_id)
  geometry_timing[["sf_conversion_seconds"]] <-
    proc.time()[["elapsed"]] - subphase_started

  # Preserve valid Terra-derived geometry unchanged. If any feature is
  # invalid, retain the established full-object st_make_valid() behavior.
  geometry_result <- validate_and_repair_patch_polygons(
    patch_polygons_sf = patch_polygons_sf,
    species_name = species_name
  )
  patch_polygons_sf <- geometry_result$polygons
  geometry_timing[names(geometry_result$timing_seconds)] <-
    geometry_result$timing_seconds
  geometry_validation_counts <- geometry_result$counts

  # Return the local candidate patch polygons.
  finish_geometry(
    patch_polygons_sf,
    candidate_cells = nrow(candidate_patch_cells),
    local_cells = local_raster_cell_count
  )
}
