# Deterministic single-species Stage 5 scientific kernel.
#
# Inputs are one validated species row plus shared raster context. The kernel
# returns either a completed filtered outcome or validated lookup, CSR graph,
# and integer patch raster information. Recovery and final transactions belong
# to stage5_checkpoints.R and stage5_workflow.R.


# ---- Diagnostics and atomic raster publication ------------------------------

format_patch_log_metric <- function(name, value) {
  if (is.null(value) || !length(value) || all(is.na(value))) return("NA")
  if (is.character(value)) return(paste(value, collapse = ","))
  if (is.logical(value)) return(paste(ifelse(value, "TRUE", "FALSE"), collapse = ","))

  value <- as.numeric(value[[1L]])
  if (!is.finite(value)) return("NA")
  if (grepl("_pct$", name)) {
    return(formatC(value, format = "f", digits = 1, big.mark = ","))
  }
  if (grepl("_km2$", name)) {
    return(formatC(value, format = "f", digits = 3, big.mark = ","))
  }
  if (grepl("(^n_|_cells$|_layers$|remaining_cells$)", name)) {
    return(formatC(value, format = "f", digits = 0, big.mark = ","))
  }
  format(value, trim = TRUE, scientific = FALSE, big.mark = ",")
}

log_patch_step <- function(species,
                           status,
                           ...,
                           species_index = NULL,
                           species_total = NULL) {
  metrics <- list(...)
  details <- if (length(metrics)) {
    values <- mapply(
      format_patch_log_metric,
      names(metrics),
      metrics,
      USE.NAMES = FALSE
    )
    paste0(" | ", paste0(names(metrics), "=", values, collapse = " | "))
  } else {
    ""
  }
  progress <- if (
    length(species_index) == 1L && is.finite(species_index) &&
      length(species_total) == 1L && is.finite(species_total)
  ) {
    sprintf("[%d/%d] ", as.integer(species_index), as.integer(species_total))
  } else {
    ""
  }

  message(sprintf(
    "[%s] %s[%s] %s%s",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    progress,
    species,
    status,
    details
  ))
}

clump_patch_raster <- function(x, backend = c("terra", "fasterRaster"), diagonal = FALSE) {
  backend <- match.arg(backend)
  if (identical(backend, "terra")) {
    return(terra::patches(x, directions = if (isTRUE(diagonal)) 8L else 4L))
  }
  tryCatch(
    terra::rast(fasterRaster::clump(fasterRaster::fast(x), diagonal = isTRUE(diagonal))),
    error = function(e) patch_abort("fasterRaster clumping failed: ", conditionMessage(e))
  )
}


write_patch_raster_atomic <- function(raster, path, patch_ids, template,
                                      overwrite = FALSE,
                                      rename_file = file.rename) {
  assert(!file.exists(path) || isTRUE(overwrite), paste0("Refusing to overwrite patch raster: ", path))
  ensure_writable_dir(dirname(path), "Stage 5 patch raster directory")
  temporary <- tempfile("stage5_patch_", tmpdir = dirname(path), fileext = ".tif")
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  terra::writeRaster(
    raster,
    filename = temporary,
    overwrite = TRUE,
    wopt = list(
      datatype = "INT4U",
      gdal = c("TILED=YES", "BLOCKXSIZE=512", "BLOCKYSIZE=512", "COMPRESS=LZW", "BIGTIFF=IF_SAFER")
    )
  )
  check <- terra::rast(temporary)
  assert(terra::nlyr(check) == 1L, "Temporary patch raster must contain exactly one layer.")
  assert(terra::compareGeom(check, template, stopOnError = FALSE),
         "Temporary patch raster geometry differs from the Stage 5 template.")
  values <- as.numeric(terra::unique(check, na.rm = TRUE)[[1L]])
  values <- sort(unique(values[is.finite(values)]))
  assert(all(values > 0 & values == floor(values)), "Temporary patch raster contains invalid patch IDs.")
  assert(identical(as.integer(values), sort(unique(as.integer(patch_ids)))),
         "Temporary patch raster IDs do not match the species patch lookup.")
  project_file_set_transaction(
    staged_paths = temporary,
    target_paths = path,
    overwrite = overwrite,
    rename_file = rename_file,
    label = "Stage 5 patch raster"
  )
  invisible(path)
}

# ---- Canonical within-PU connectivity ----------------------------------------

# Each PU graph is emitted in deterministic patch order as zero-based CSR
# adjacency; retained patch and PU identifiers are already canonical here.
build_pu_connectivity <- function(nb, patches_sf, patch_to_pu_kept, final_patch_ids, species_name) {
  upper_neighbors <- lapply(seq_along(nb), function(i) {
    neighbors <- as.integer(nb[[i]])
    neighbors[neighbors > i]
  })
  edge_i <- unlist(Map(rep.int, seq_along(upper_neighbors), lengths(upper_neighbors)), use.names = FALSE)
  edge_j <- unlist(upper_neighbors, use.names = FALSE)
  edges <- data.frame(
    patch_i = patches_sf$patch_id[edge_i],
    patch_j = patches_sf$patch_id[edge_j]
  )

  keep_prefilter_ids <- as.integer(names(final_patch_ids))
  edges <- edges[
    edges$patch_i %in% keep_prefilter_ids & edges$patch_j %in% keep_prefilter_ids,
    ,
    drop = FALSE
  ]

  if (nrow(edges)) {
    edges$patch_i <- unname(as.integer(final_patch_ids[as.character(edges$patch_i)]))
    edges$patch_j <- unname(as.integer(final_patch_ids[as.character(edges$patch_j)]))
    edges <- edges[stats::complete.cases(edges), , drop = FALSE]
  }

  patch_to_pu <- stats::setNames(
    as.integer(patch_to_pu_kept$pu_id),
    as.character(patch_to_pu_kept$patch_id)
  )
  if (nrow(edges)) {
    edge_pu_i <- unname(patch_to_pu[as.character(edges$patch_i)])
    edge_pu_j <- unname(patch_to_pu[as.character(edges$patch_j)])
    edges <- edges[!is.na(edge_pu_i) & edge_pu_i == edge_pu_j, , drop = FALSE]
    edge_pu <- edge_pu_i[!is.na(edge_pu_i) & edge_pu_i == edge_pu_j]
    edges_by_pu <- split(edges, edge_pu)
  } else {
    edges_by_pu <- list()
  }

  out <- list()

  for (pu in sort(unique(patch_to_pu_kept$pu_id))) {
    patch_ids <- sort(unique(patch_to_pu_kept$patch_id[patch_to_pu_kept$pu_id == pu]))
    patch_ids <- patch_ids[!is.na(patch_ids)]
    if (!length(patch_ids)) next

    idx_map <- seq_along(patch_ids)
    names(idx_map) <- as.character(patch_ids)

    sub_edges <- edges_by_pu[[as.character(pu)]] %||% edges[0, , drop = FALSE]

    if (!nrow(sub_edges)) {
      row_ptr <- as.integer(c(0L, rep.int(0L, length(patch_ids))))
      col_idx <- integer(0)
    } else {
      symmetric_edges <- rbind(
        data.frame(i = sub_edges$patch_i, j = sub_edges$patch_j),
        data.frame(i = sub_edges$patch_j, j = sub_edges$patch_i)
      )
      symmetric_edges <- unique(symmetric_edges[symmetric_edges$i != symmetric_edges$j, , drop = FALSE])
      symmetric_edges <- symmetric_edges[order(symmetric_edges$i, symmetric_edges$j), , drop = FALSE]

      row_local <- as.integer(idx_map[as.character(symmetric_edges$i)])
      col_local <- as.integer(idx_map[as.character(symmetric_edges$j)])
      ok <- !is.na(row_local) & !is.na(col_local)
      row_local <- row_local[ok]
      col_local <- col_local[ok]

      split_neighbors <- split(col_local - 1L, row_local)
      adj_list <- vector("list", length(patch_ids))
      for (i in seq_along(adj_list)) {
        adj_list[[i]] <- as.integer(sort(unique(split_neighbors[[as.character(i)]] %||% integer(0))))
      }

      lens <- vapply(adj_list, length, integer(1))
      row_ptr <- as.integer(c(0L, cumsum(lens)))
      col_idx <- as.integer(unlist(adj_list, use.names = FALSE))
    }

    key <- paste0(species_name, "|", pu)
    out[[key]] <- list(
      species   = species_name,
      pu_id     = as.integer(pu),
      patch_ids = as.integer(patch_ids),
      row_ptr   = as.integer(row_ptr),
      col_idx   = as.integer(col_idx)
    )
  }

  out
}

# ---- One-species patch and population-unit construction ---------------------

# This kernel owns the large raster, polygon, and graph intermediates for one
# species at a time; callers can release them before beginning the next species.
patch_components_for_species <- function(row,
                                         mapped_habitat,
                                         species_name = row$scientificName,
                                         mapped_habitat_cells = NULL,
                                         mapped_habitat_area_km2 = NULL,
                                         cell_area_km2 = NULL,
                                         clump_backend = "terra",
                                         retain_intermediates = FALSE) {
  if (is.null(mapped_habitat_cells)) {
    mapped_habitat_cells <- non_na_cell_count(mapped_habitat)
  }
  if (is.null(mapped_habitat_area_km2)) {
    mapped_habitat_area_km2 <- raster_area_km2(mapped_habitat)
  }

  component_result <- function(status, stop_stage = NA_character_, stop_reason = NA_character_, ...) {
    c(
      list(
        status = status,
        species = species_name,
        stop_stage = stop_stage,
        stop_reason = stop_reason,
        mapped_habitat_cells = mapped_habitat_cells,
        mapped_habitat_area_km2 = mapped_habitat_area_km2
      ),
      list(...)
    )
  }

  # Habitat intersection is already complete. Rook-connected clumping
  # (diagonal = FALSE) defines the raw patches on the trimmed occupied extent.
  mapped_habitat_trim <- terra::trim(mapped_habitat)
  clump_started <- Sys.time()
  patch_raw <- clump_patch_raster(mapped_habitat_trim, backend = clump_backend, diagonal = FALSE)
  clump_elapsed <- as.numeric(difftime(Sys.time(), clump_started, units = "secs"))
  names(patch_raw) <- "patch_id_raw"

  cell_area <- if (is.null(cell_area_km2)) {
    terra::cellSize(mapped_habitat_trim, unit = "km")
  } else {
    terra::crop(cell_area_km2, mapped_habitat_trim, snap = "near")
  }
  names(cell_area) <- "patch_area_km2"
  cell_count <- terra::ifel(is.na(mapped_habitat_trim), NA_integer_, 1L)
  names(cell_count) <- "patch_cells"
  patch_area <- as.data.frame(terra::zonal(
    c(cell_area, cell_count),
    patch_raw,
    "sum",
    na.rm = TRUE
  ))
  if (!nrow(patch_area)) {
    return(component_result(
      status = "skipped_no_patch_ids",
      stop_stage = "patch_clumping",
      stop_reason = "Clumping produced no patch identifiers.",
      n_raw_patches = 0L,
      raw_patch_cells = 0,
      raw_patch_area_km2 = 0
    ))
  }

  names(patch_area)[1:3] <- c("patch_id_raw", "patch_area_km2", "patch_cells")
  patch_area <- patch_area[!is.na(patch_area$patch_id_raw), , drop = FALSE]
  if (!nrow(patch_area)) {
    return(component_result(
      status = "skipped_no_valid_patch_ids",
      stop_stage = "patch_clumping",
      stop_reason = "Clumping produced no valid patch identifiers.",
      n_raw_patches = 0L,
      raw_patch_cells = 0,
      raw_patch_area_km2 = 0
    ))
  }

  # Patch areas use per-cell geodesic area, and the minimum-area predicate is
  # intentionally strict; threshold equality does not retain a patch.
  n_raw_patches <- nrow(patch_area)
  raw_patch_cells <- sum(patch_area$patch_cells, na.rm = TRUE)
  raw_patch_area_km2 <- sum(patch_area$patch_area_km2, na.rm = TRUE)
  keep_raw_patches <- patch_area$patch_id_raw[
    area_exceeds_threshold(patch_area$patch_area_km2, row$min_patch_km2)
  ]
  retained_patch_rows <- patch_area$patch_id_raw %in% keep_raw_patches
  retained_patch_cells <- sum(patch_area$patch_cells[retained_patch_rows], na.rm = TRUE)
  retained_patch_area_km2 <- sum(patch_area$patch_area_km2[retained_patch_rows], na.rm = TRUE)

  if (!length(keep_raw_patches)) {
    return(component_result(
      status = "skipped_no_patches_above_min_patch",
      stop_stage = "patch_area_filter",
      stop_reason = "No raw patches exceeded the minimum patch-area threshold.",
      n_raw_patches = n_raw_patches,
      raw_patch_cells = raw_patch_cells,
      raw_patch_area_km2 = raw_patch_area_km2,
      n_patches_after_patch_filter = 0L,
      retained_patch_cells = 0,
      retained_patch_area_km2 = 0
    ))
  }

  patch_prefilter_trim <- terra::classify(
    patch_raw,
    rcl    = cbind(as.integer(keep_raw_patches), seq_along(keep_raw_patches)),
    others = NA_integer_
  )
  names(patch_prefilter_trim) <- "patch_id"

  patch_prefilter <- terra::extend(patch_prefilter_trim, mapped_habitat)
  names(patch_prefilter) <- "patch_id"

  patch_table_prefilter <- data.frame(
    patch_id_prefilter = seq_along(keep_raw_patches),
    patch_area_km2 = patch_area$patch_area_km2[match(keep_raw_patches, patch_area$patch_id_raw)],
    patch_cells = patch_area$patch_cells[match(keep_raw_patches, patch_area$patch_id_raw)]
  )

  polygon_started <- Sys.time()
  patch_polygons <- tryCatch(
    sf::st_as_sf(terra::as.polygons(patch_prefilter, values = TRUE, dissolve = TRUE, na.rm = TRUE)),
    error = function(e) NULL
  )
  polygon_elapsed <- as.numeric(difftime(Sys.time(), polygon_started, units = "secs"))
  if (is.null(patch_polygons) || !nrow(patch_polygons)) {
    return(component_result(
      status = "skipped_polygonization_failed",
      stop_stage = "patch_polygonization",
      stop_reason = "Retained patch cells could not be polygonized.",
      n_raw_patches = n_raw_patches,
      raw_patch_cells = raw_patch_cells,
      raw_patch_area_km2 = raw_patch_area_km2,
      n_patches_after_patch_filter = length(keep_raw_patches),
      retained_patch_cells = retained_patch_cells,
      retained_patch_area_km2 = retained_patch_area_km2
    ))
  }

  value_col <- setdiff(names(patch_polygons), attr(patch_polygons, "sf_column"))[1]
  patch_polygons$patch_id <- as.integer(patch_polygons[[value_col]])
  patch_polygons <- patch_polygons[!is.na(patch_polygons$patch_id), ]
  if (!nrow(patch_polygons)) {
    return(component_result(
      status = "skipped_no_patch_polygons",
      stop_stage = "patch_polygonization",
      stop_reason = "Polygonization produced no valid patch polygons.",
      n_raw_patches = n_raw_patches,
      raw_patch_cells = raw_patch_cells,
      raw_patch_area_km2 = raw_patch_area_km2,
      n_patches_after_patch_filter = length(keep_raw_patches),
      retained_patch_cells = retained_patch_cells,
      retained_patch_area_km2 = retained_patch_area_km2
    ))
  }

  # Patches within the species dispersal distance form candidate population
  # units through connected components of the undirected distance graph.
  distance_started <- Sys.time()
  neighbors <- sf::st_is_within_distance(patch_polygons, patch_polygons, dist = row$disp_km * 1000)
  patch_graph <- igraph::graph_from_adj_list(neighbors, mode = "out") |>
    igraph::as_undirected("collapse") |>
    igraph::simplify(remove.multiple = TRUE, remove.loops = TRUE)

  component_id <- igraph::components(patch_graph)$membership
  distance_elapsed <- as.numeric(difftime(Sys.time(), distance_started, units = "secs"))
  patch_polygons$pu_id_prefilter <- as.integer(match(component_id, sort(unique(component_id))))

  patch_to_pu_prefilter <- data.frame(
    patch_id_prefilter = patch_polygons$patch_id,
    pu_id_prefilter = patch_polygons$pu_id_prefilter
  )

  patch_area_with_pu <- merge(patch_table_prefilter, patch_to_pu_prefilter, by = "patch_id_prefilter")
  pu_area <- aggregate(
    cbind(patch_area_km2, patch_cells) ~ pu_id_prefilter,
    patch_area_with_pu,
    sum,
    na.rm = TRUE
  )
  names(pu_area) <- c("pu_id_prefilter", "pu_area_km2", "pu_cells")
  candidate_pu_area_km2 <- sum(pu_area$pu_area_km2, na.rm = TRUE)
  candidate_pu_cells <- sum(pu_area$pu_cells, na.rm = TRUE)

  # Apply the strict population-area filter before canonical renumbering, so
  # final PU and patch IDs are dense, stable, and independent of discarded IDs.
  keep_prefilter_pu <- sort(pu_area$pu_id_prefilter[
    area_exceeds_threshold(pu_area$pu_area_km2, row$min_pop_km2)
  ])
  if (!length(keep_prefilter_pu)) {
    return(component_result(
      status = "skipped_no_pus_above_min_pop",
      stop_stage = "population_unit_area_filter",
      stop_reason = "No candidate population units exceeded the minimum population-area threshold.",
      n_raw_patches = n_raw_patches,
      raw_patch_cells = raw_patch_cells,
      raw_patch_area_km2 = raw_patch_area_km2,
      n_patches_after_patch_filter = length(keep_raw_patches),
      retained_patch_cells = retained_patch_cells,
      retained_patch_area_km2 = retained_patch_area_km2,
      n_candidate_pus = nrow(pu_area),
      candidate_pu_cells = candidate_pu_cells,
      candidate_pu_area_km2 = candidate_pu_area_km2,
      n_final_pus = 0L,
      n_final_patches = 0L,
      final_cells = 0,
      final_area_km2 = 0
    ))
  }

  patch_to_pu_prefilter$pu_id <- match(patch_to_pu_prefilter$pu_id_prefilter, keep_prefilter_pu)
  patch_to_pu_kept <- patch_to_pu_prefilter[!is.na(patch_to_pu_prefilter$pu_id), , drop = FALSE]

  keep_prefilter_patch_ids <- sort(unique(patch_to_pu_kept$patch_id_prefilter))
  final_patch_ids <- seq_along(keep_prefilter_patch_ids)
  names(final_patch_ids) <- as.character(keep_prefilter_patch_ids)

  patch_final <- terra::classify(
    patch_prefilter,
    rcl    = cbind(as.integer(keep_prefilter_patch_ids), final_patch_ids),
    others = NA_integer_
  )
  names(patch_final) <- "patch_id"

  patch_to_pu_kept$patch_id <- as.integer(final_patch_ids[as.character(patch_to_pu_kept$patch_id_prefilter)])

  patch_lookup <- data.frame(
    scientificName = species_name,
    patch_id = patch_to_pu_kept$patch_id,
    pu_id = patch_to_pu_kept$pu_id,
    patch_area_km2 = patch_table_prefilter$patch_area_km2[
      match(patch_to_pu_kept$patch_id_prefilter, patch_table_prefilter$patch_id_prefilter)
    ],
    stringsAsFactors = FALSE
  )

  connectivity <- build_pu_connectivity(
    nb = neighbors,
    patches_sf = patch_polygons,
    patch_to_pu_kept = patch_to_pu_kept[, c("patch_id", "pu_id")],
    final_patch_ids = final_patch_ids,
    species_name = species_name
  )

  final_pu_rows <- pu_area$pu_id_prefilter %in% keep_prefilter_pu
  final_cells <- sum(pu_area$pu_cells[final_pu_rows], na.rm = TRUE)
  final_area_km2 <- sum(pu_area$pu_area_km2[final_pu_rows], na.rm = TRUE)

  # Large diagnostic rasters and geometries remain optional; the default result
  # transfers only the objects needed for publication and downstream analysis.
  result <- component_result(
    status = "retained",
    stop_stage = "completed",
    stop_reason = NA_character_,
    n_raw_patches = n_raw_patches,
    raw_patch_cells = raw_patch_cells,
    raw_patch_area_km2 = raw_patch_area_km2,
    n_patches_after_patch_filter = length(keep_raw_patches),
    retained_patch_cells = retained_patch_cells,
    retained_patch_area_km2 = retained_patch_area_km2,
    n_candidate_pus = nrow(pu_area),
    candidate_pu_cells = candidate_pu_cells,
    candidate_pu_area_km2 = candidate_pu_area_km2,
    n_final_pus = length(keep_prefilter_pu),
    n_final_patches = length(unique(patch_lookup$patch_id)),
    final_cells = final_cells,
    final_area_km2 = final_area_km2,
    clump_elapsed_seconds = clump_elapsed,
    polygon_elapsed_seconds = polygon_elapsed,
    distance_elapsed_seconds = distance_elapsed,
    patch_lookup = patch_lookup,
    connectivity = connectivity,
    patch_final = patch_final
  )
  if (isTRUE(retain_intermediates)) {
    patch_filter_status_trim <- terra::ifel(
      !is.na(patch_prefilter_trim), 2L,
      terra::ifel(!is.na(patch_raw), 1L, NA_integer_)
    )
    patch_filter_status <- terra::extend(patch_filter_status_trim, mapped_habitat)
    names(patch_filter_status) <- "patch_status"
    pu_candidate <- terra::classify(
      patch_prefilter,
      rcl = as.matrix(patch_to_pu_prefilter[, c("patch_id_prefilter", "pu_id_prefilter")]),
      others = NA_integer_
    )
    names(pu_candidate) <- "pu_prefilter"
    pu_final <- terra::classify(
      patch_prefilter,
      rcl = as.matrix(patch_to_pu_kept[, c("patch_id_prefilter", "pu_id")]),
      others = NA_integer_
    )
    names(pu_final) <- "pu_id"
    final_units_sf <- patch_polygons |>
      dplyr::mutate(
        pu_id = patch_to_pu_kept$pu_id[match(patch_id, patch_to_pu_kept$patch_id_prefilter)]
      ) |>
      dplyr::filter(!is.na(pu_id)) |>
      dplyr::group_by(pu_id) |>
      dplyr::summarise(geometry = sf::st_union(geometry), .groups = "drop")
    result <- c(result, list(
      patch_filter_status = patch_filter_status,
      pu_candidate = pu_candidate,
      pu_final = pu_final,
      final_units_sf = final_units_sf,
      kept_candidate_pus = keep_prefilter_pu
    ))
  }
  result
}

# ---- Species orchestration and publication ----------------------------------

# Construct and validate a complete single-species result before publishing its
# patch raster atomically; no later species shares ownership of intermediates.
process_patch_species <- function(row,
                                  paths,
                                  run_options,
                                  habitat_masks,
                                  template,
                                  cell_area_km2 = NULL,
                                  clump_backend = "terra",
                                  verbose = FALSE,
                                  label_map = habitat_label_map(),
                                  species_index = NULL,
                                  species_total = NULL) {
  species_name <- row$scientificName
  species_started <- Sys.time()
  log_step <- function(status, ...) {
    log_patch_step(
      species_name,
      status,
      ...,
      species_index = species_index,
      species_total = species_total
    )
  }
  base_result <- list(
    species = species_name,
    class_lc = row$class_lc,
    sdm_method = row$sdm_method
  )
  finish_filtered <- function(result) {
    result$elapsed_seconds <- as.numeric(difftime(Sys.time(), species_started, units = "secs"))
    log_step(
      "FILTERED_OUT",
      stop_stage = result$stop_stage,
      remaining_cells = result$remaining_cells %||% NA_real_,
      reason = result$stop_reason,
      elapsed_seconds = result$elapsed_seconds
    )
    result
  }

  log_step(
    "START",
    taxon = row$class_lc,
    sdm_method = row$sdm_method
  )

  # Intersect the aligned SDM presence with the species habitat masks before
  # any patch topology or area threshold is evaluated.
  mapped <- mapped_habitat_for_species(
    row, habitat_masks, template,
    cell_area_km2 = cell_area_km2,
    label_map = label_map
  )
  if (!identical(mapped$status, "ok")) {
    return(finish_filtered(utils::modifyList(base_result, list(
      status = mapped$status,
      stop_stage = "habitat_label_mapping",
      stop_reason = "No supplied habitat labels mapped to available land-cover masks.",
      habitat_labels = mapped$labels,
      habitat_mask_names = mapped$mask_names,
      remaining_cells = NA_real_
    ))))
  }

  base_result <- utils::modifyList(base_result, list(
    source_sdm_layers = mapped$source_sdm_layers,
    source_sdm_total_cells = mapped$source_sdm_total_cells,
    source_sdm_non_na_cells = mapped$source_sdm_non_na_cells,
    source_sdm_presence_cells = mapped$source_sdm_presence_cells,
    alignment_action = mapped$alignment_action,
    aligned_presence_cells = mapped$aligned_presence_cells,
    habitat_labels = mapped$labels,
    habitat_mask_names = mapped$mask_names,
    mapped_habitat_cells = mapped$mapped_habitat_cells
  ))
  mapped_habitat_area_km2 <- mapped$mapped_habitat_area_km2
  base_result$mapped_habitat_area_km2 <- mapped_habitat_area_km2
  log_step(
    "HABITAT",
    alignment = mapped$alignment_action,
    mapped_cells = mapped$mapped_habitat_cells,
    mapped_area_km2 = mapped_habitat_area_km2
  )
  if (isTRUE(verbose)) {
    log_step(
      "INPUT_DETAIL",
      source_cells = mapped$source_sdm_total_cells,
      source_non_na_cells = mapped$source_sdm_non_na_cells,
      source_presence_cells = mapped$source_sdm_presence_cells,
      habitat_masks = paste(mapped$mask_names, collapse = ",")
    )
  }

  if (mapped$aligned_presence_cells <= 0) {
    return(finish_filtered(utils::modifyList(base_result, list(
      status = "skipped_no_aligned_presence",
      stop_stage = "aligned_sdm_presence",
      stop_reason = "The aligned SDM contained no presence cells.",
      remaining_cells = 0
    ))))
  }
  if (mapped$mapped_habitat_cells <= 0) {
    return(finish_filtered(utils::modifyList(base_result, list(
      status = "skipped_no_habitat_intersection",
      stop_stage = "habitat_sdm_intersection",
      stop_reason = "The habitat and aligned SDM rasters had no intersecting cells.",
      remaining_cells = 0
    ))))
  }

  components <- patch_components_for_species(
    row,
    mapped$mapped_habitat,
    species_name,
    mapped_habitat_cells = mapped$mapped_habitat_cells,
    mapped_habitat_area_km2 = mapped_habitat_area_km2,
    cell_area_km2 = cell_area_km2,
    clump_backend = clump_backend
  )
  components <- utils::modifyList(base_result, components)

  has_patch_summary <- is.finite(components$n_raw_patches %||% NA_real_)
  has_pu_summary <- is.finite(components$n_candidate_pus %||% NA_real_)
  if (has_patch_summary || has_pu_summary) {
    log_step(
      "STRUCTURE",
      patches = if (has_patch_summary) {
        paste0(
          components$n_raw_patches,
          "->",
          components$n_patches_after_patch_filter %||% NA_integer_
        )
      } else {
        NA_character_
      },
      pus = if (has_pu_summary) {
        paste0(
          components$n_candidate_pus,
          "->",
          components$n_final_pus %||% NA_integer_
        )
      } else {
        NA_character_
      },
      final_area_km2 = components$final_area_km2 %||% NA_real_
    )
  }
  if (isTRUE(verbose) && identical(components$status, "retained")) {
    log_step(
      "TIMING",
      clump_seconds = components$clump_elapsed_seconds,
      polygon_seconds = components$polygon_elapsed_seconds,
      distance_seconds = components$distance_elapsed_seconds
    )
  }

  if (!identical(components$status, "retained")) {
    components$remaining_cells <- switch(
      components$stop_stage,
      patch_polygonization = components$retained_patch_cells %||% NA_real_,
      0
    )
    return(finish_filtered(components))
  }

  out_file <- file.path(paths$patch_dir, patch_filename_from_scientific(species_name))
  if (file.exists(out_file) && !run_options$overwrite_patch_rasters) {
    patch_abort(paste0("Output patch raster exists and overwrite_patch_rasters is FALSE: ", out_file))
  }

  write_patch_raster_atomic(
    components$patch_final,
    path = out_file,
    patch_ids = components$patch_lookup$patch_id,
    template = template,
    overwrite = run_options$overwrite_patch_rasters
  )

  components$patch_final <- NULL
  components$class_lc <- row$class_lc
  components$sdm_method <- row$sdm_method
  components$patch_raster <- out_file
  components$elapsed_seconds <- as.numeric(difftime(Sys.time(), species_started, units = "secs"))
  log_step(
    "RETAINED",
    patch_raster = basename(out_file),
    elapsed_seconds = components$elapsed_seconds
  )
  components
}

patch_status_summary <- function(results) {
  dplyr::bind_rows(lapply(results, function(x) {
    data.frame(
      species = x$species %||% NA_character_,
      class_lc = x$class_lc %||% NA_character_,
      sdm_method = x$sdm_method %||% NA_character_,
      status = x$status %||% NA_character_,
      stop_stage = x$stop_stage %||% NA_character_,
      stop_reason = x$stop_reason %||% NA_character_,
      source_sdm_layers = x$source_sdm_layers %||% NA_integer_,
      source_sdm_total_cells = x$source_sdm_total_cells %||% NA_real_,
      source_sdm_non_na_cells = x$source_sdm_non_na_cells %||% NA_real_,
      source_sdm_presence_cells = x$source_sdm_presence_cells %||% NA_real_,
      alignment_action = x$alignment_action %||% NA_character_,
      aligned_presence_cells = x$aligned_presence_cells %||% NA_real_,
      mapped_habitat_cells = x$mapped_habitat_cells %||% NA_real_,
      mapped_habitat_area_km2 = x$mapped_habitat_area_km2 %||% NA_real_,
      n_raw_patches = x$n_raw_patches %||% NA_integer_,
      raw_patch_cells = x$raw_patch_cells %||% NA_real_,
      raw_patch_area_km2 = x$raw_patch_area_km2 %||% NA_real_,
      n_patches_after_patch_filter = x$n_patches_after_patch_filter %||% NA_integer_,
      retained_patch_cells = x$retained_patch_cells %||% NA_real_,
      retained_patch_area_km2 = x$retained_patch_area_km2 %||% NA_real_,
      n_candidate_pus = x$n_candidate_pus %||% NA_integer_,
      candidate_pu_cells = x$candidate_pu_cells %||% NA_real_,
      candidate_pu_area_km2 = x$candidate_pu_area_km2 %||% NA_real_,
      n_final_pus = x$n_final_pus %||% NA_integer_,
      n_final_patches = x$n_final_patches %||% NA_integer_,
      final_cells = x$final_cells %||% NA_real_,
      final_area_km2 = x$final_area_km2 %||% NA_real_,
      remaining_cells = x$remaining_cells %||% NA_real_,
      patch_raster = x$patch_raster %||% NA_character_,
      clump_elapsed_seconds = x$clump_elapsed_seconds %||% NA_real_,
      polygon_elapsed_seconds = x$polygon_elapsed_seconds %||% NA_real_,
      distance_elapsed_seconds = x$distance_elapsed_seconds %||% NA_real_,
      elapsed_seconds = x$elapsed_seconds %||% NA_real_,
      error = x$error %||% NA_character_,
      stringsAsFactors = FALSE
    )
  }))
}
