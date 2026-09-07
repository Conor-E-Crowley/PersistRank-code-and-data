# Output writers for the priority-removal sequence.
#
# Stage lookup and event rows preserve their serialized columns and canonical
# ordering. Final rasters retain template geometry, established datatypes, NA
# treatment, and cell values. Multi-file commitment uses same-directory
# staging and rollback without hashing or rereading large rasters.


# Stage lookup and event tables are scientific outputs. Validation is
# centralized here so the pipeline loop only orchestrates state transitions.
validate_priority_output_integer <- function(x, label, positive = FALSE) {
  values <- suppressWarnings(as.numeric(x))
  lower_ok <- if (isTRUE(positive)) values > 0 else values >= 0
  assert(
    length(values) == length(x) &&
      all(is.finite(values)) &&
      all(lower_ok) &&
      all(values == floor(values)),
    paste0(
      label,
      if (isTRUE(positive)) {
        " must contain positive integers."
      } else {
        " must contain non-negative integers."
      }
    )
  )
  as.integer(values)
}

validate_priority_output_number <- function(x, label, positive = FALSE) {
  values <- suppressWarnings(as.numeric(x))
  lower_ok <- if (isTRUE(positive)) values > 0 else values >= 0
  assert(
    length(values) == length(x) && all(is.finite(values)) && all(lower_ok),
    paste0(
      label,
      if (isTRUE(positive)) {
        " must contain positive finite values."
      } else {
        " must contain non-negative finite values."
      }
    )
  )
  values
}

validate_stage_patch_lookup_table <- function(
  patch_table,
  label = "stage patch lookup table"
) {
  assert(
    is.data.frame(patch_table),
    paste0(label, " must be a data.frame or data.table.")
  )
  need_cols(
    patch_table,
    c("species", "patch_id", "pu_id", "patch_area_km2"),
    label
  )
  if (!nrow(patch_table)) return(invisible(TRUE))

  species <- trimws(as.character(patch_table$species))
  assert(
    all(!is.na(species)) && all(nzchar(species)),
    paste0(label, " contains blank species values.")
  )
  validate_priority_output_integer(
    patch_table$patch_id,
    paste0(label, "$patch_id"),
    positive = TRUE
  )
  validate_priority_output_integer(
    patch_table$pu_id,
    paste0(label, "$pu_id"),
    positive = TRUE
  )
  validate_priority_output_number(
    patch_table$patch_area_km2,
    paste0(label, "$patch_area_km2"),
    positive = TRUE
  )
  invisible(TRUE)
}

write_stage_patch_lookup_table <- function(
  patch_table,
  stage_index,
  patch_lookup_output_dir
) {
  stage_index <- validate_priority_count(stage_index, "stage_index")
  validate_stage_patch_lookup_table(patch_table)
  output_path <- file.path(
    patch_lookup_output_dir,
    sprintf("stage_patch_lookup_stage_%04d.csv", stage_index)
  )
  tmp_path <- tempfile(
    pattern = paste0(basename(output_path), "."),
    tmpdir = patch_lookup_output_dir,
    fileext = ".tmp"
  )
  on.exit(unlink(tmp_path, force = TRUE), add = TRUE)
  data.table::fwrite(patch_table, tmp_path)
  priority_commit_file_set(
    staged_paths = tmp_path,
    target_paths = output_path,
    overwrite = FALSE,
    label = "Stage 6 patch lookup"
  )
  runtime_log_event(
    "stage_patch_lookup_written",
    stage = stage_index,
    rows = as.integer(nrow(patch_table)),
    path = normalizePath(output_path, mustWork = FALSE)
  )
  invisible(output_path)
}

validate_removal_events_table <- function(
  removal_events,
  label = "removal_events.csv"
) {
  assert(
    is.data.frame(removal_events),
    paste0(label, " must be a data.frame or data.table.")
  )
  need_cols(removal_events, c(
    "removal_step", "stage", "event_type", "pruning_iteration",
    "cells_removed", "area_removed_km2", "cum_cells_removed",
    "cum_area_removed_km2", "cum_prop_cells_removed",
    "cum_prop_area_removed", "cells_retained", "area_retained_km2"
  ), label)
  if (!nrow(removal_events)) return(invisible(TRUE))

  removal_step <- validate_priority_output_integer(
    removal_events$removal_step,
    paste0(label, "$removal_step"),
    positive = TRUE
  )
  stage <- validate_priority_output_integer(
    removal_events$stage,
    paste0(label, "$stage")
  )
  cells_removed <- validate_priority_output_integer(
    removal_events$cells_removed,
    paste0(label, "$cells_removed")
  )
  cum_cells_removed <- validate_priority_output_integer(
    removal_events$cum_cells_removed,
    paste0(label, "$cum_cells_removed")
  )
  cells_retained <- validate_priority_output_integer(
    removal_events$cells_retained,
    paste0(label, "$cells_retained")
  )

  event_type <- as.character(removal_events$event_type)
  valid_event_types <- c("prune", "fragmentation", "distance", "final_retained")
  assert(
    all(!is.na(event_type)) && all(event_type %in% valid_event_types),
    paste0(label, " contains unknown event_type values.")
  )
  assert(
    sum(event_type == "final_retained") == 1L,
    paste0(label, " must contain exactly one final_retained event.")
  )

  pruning_iteration <- suppressWarnings(as.numeric(removal_events$pruning_iteration))
  prune_rows <- event_type == "prune"
  assert(
    all(!prune_rows | (
      is.finite(pruning_iteration) &
        pruning_iteration > 0 &
        pruning_iteration == floor(pruning_iteration)
    )) &&
      all(prune_rows | is.na(pruning_iteration)),
    paste0(
      label,
      "$pruning_iteration must be positive integers for prune events and NA otherwise."
    )
  )

  area_removed <- validate_priority_output_number(
    removal_events$area_removed_km2,
    paste0(label, "$area_removed_km2")
  )
  cum_area_removed <- validate_priority_output_number(
    removal_events$cum_area_removed_km2,
    paste0(label, "$cum_area_removed_km2")
  )
  area_retained <- validate_priority_output_number(
    removal_events$area_retained_km2,
    paste0(label, "$area_retained_km2")
  )

  assert(
    identical(removal_step, seq_len(nrow(removal_events))),
    paste0(label, "$removal_step must be contiguous from 1.")
  )
  assert(
    identical(cum_cells_removed, as.integer(cumsum(cells_removed))),
    paste0(label, "$cum_cells_removed does not match cumulative cells_removed.")
  )
  assert(
    all(abs(cum_area_removed - cumsum(area_removed)) <= 1e-6),
    paste0(
      label,
      "$cum_area_removed_km2 does not match cumulative area_removed_km2."
    )
  )
  assert(
    all(cells_retained >= 0L) && all(area_retained >= 0) && all(stage >= 0L),
    paste0(label, " contains invalid retained values or stage indices.")
  )

  for (column in c("cum_prop_cells_removed", "cum_prop_area_removed")) {
    values <- suppressWarnings(as.numeric(removal_events[[column]]))
    assert(
      all(is.finite(values)) && all(values >= 0) && all(values <= 1 + 1e-9),
      paste0(label, "$", column, " must contain values in [0, 1].")
    )
  }
  invisible(TRUE)
}

# Rename a staged file set into place without hashing or rereading large raster
# payloads. Same-directory backups make rollback a metadata operation.
priority_commit_file_set <- function(staged_paths,
                                     target_paths,
                                     overwrite = FALSE,
                                     rename_file = file.rename,
                                     label = "Stage 6") {
  staged_paths <- as.character(staged_paths)
  target_paths <- as.character(target_paths)
  assert(
    length(staged_paths) > 0L && length(staged_paths) == length(target_paths),
    paste0(label, " commit requires equal, nonempty path sets.")
  )
  assert(
    !anyDuplicated(staged_paths) && !anyDuplicated(target_paths),
    paste0(label, " commit paths must be unique.")
  )
  assert(
    all(file.exists(staged_paths)),
    paste0(label, " staged output is missing: ", paste(staged_paths[!file.exists(staged_paths)], collapse = ", "))
  )

  existed <- file.exists(target_paths)
  assert(
    isTRUE(overwrite) || !any(existed),
    paste0(label, " output already exists: ", paste(target_paths[existed], collapse = ", "))
  )
  backups <- paste0(target_paths, ".stage6_backup")
  assert(!any(file.exists(backups)), paste0(label, " has unresolved backup files."))

  backed_up <- committed <- logical(length(target_paths))
  rollback <- function() {
    unlink(target_paths[committed], force = TRUE)
    for (i in which(backed_up)) {
      if (!isTRUE(rename_file(backups[[i]], target_paths[[i]]))) {
        stop("Could not restore Stage 6 output: ", target_paths[[i]], call. = FALSE)
      }
    }
    unlink(c(staged_paths, backups), force = TRUE)
  }

  failure <- NULL
  tryCatch({
    for (i in which(existed)) {
      if (!isTRUE(rename_file(target_paths[[i]], backups[[i]]))) {
        stop("Could not back up Stage 6 output: ", target_paths[[i]], call. = FALSE)
      }
      backed_up[[i]] <- TRUE
    }
    for (i in seq_along(target_paths)) {
      if (!isTRUE(rename_file(staged_paths[[i]], target_paths[[i]]))) {
        stop("Could not commit Stage 6 output: ", target_paths[[i]], call. = FALSE)
      }
      committed[[i]] <- TRUE
    }
  }, error = function(e) failure <<- e)

  if (!is.null(failure)) {
    rollback()
    stop(failure)
  }

  unlink(backups, force = TRUE)
  invisible(target_paths)
}

validate_staged_priority_raster <- function(path,
                                            template_raster,
                                            label,
                                            expected_datatype) {
  expected_datatype <- validate_scalar_string(expected_datatype, "expected_datatype")
  assert(file.exists(path), paste0("Missing staged ", label, ": ", path))
  assert(file.info(path)$size[[1L]] > 0, paste0("Empty staged ", label, ": ", path))
  staged <- terra::rast(path)
  assert(terra::nlyr(staged) == 1L, paste0(label, " must have one layer."))
  assert(
    isTRUE(terra::compareGeom(staged, template_raster[[1L]], stopOnError = FALSE)),
    paste0(label, " geometry differs from the Stage 6 template.")
  )
  found_datatype <- terra::datatype(staged)
  assert(
    length(found_datatype) == 1L && identical(found_datatype, expected_datatype),
    paste0(
      label, " must use datatype ", expected_datatype,
      "; found ", paste(found_datatype, collapse = ", "), "."
    )
  )
  invisible(path)
}

validate_removal_order_values <- function(removal_order_by_cell, label = "removal_order_by_cell") {
  values <- suppressWarnings(as.numeric(removal_order_by_cell))
  ok <- is.na(values) | (is.finite(values) & values > 0 & values == floor(values))
  assert(
    length(values) == length(removal_order_by_cell) && all(ok),
    paste0(label, " must contain positive integer removal steps or NA.")
  )
  as.integer(values)
}

validate_priority_final_cell_state <- function(
  removal_order_by_cell,
  initial_alive_by_cell,
  removal_events,
  label = "Stage 6 final cell state"
) {
  values <- validate_removal_order_values(
    removal_order_by_cell,
    paste0(label, "$removal_order_by_cell")
  )
  assert(
    is.logical(initial_alive_by_cell) &&
      length(initial_alive_by_cell) == length(values) &&
      all(!is.na(initial_alive_by_cell)),
    paste0(label, "$initial_alive_by_cell must be a complete logical vector of matching length.")
  )
  assert(is.data.frame(removal_events), paste0(label, "$removal_events must be a table."))
  need_cols(removal_events, c("removal_step", "cells_removed"), paste0(label, "$removal_events"))
  assert(nrow(removal_events) > 0L, paste0(label, "$removal_events is empty."))

  missing_inside <- initial_alive_by_cell & is.na(values)
  ranked_outside <- !initial_alive_by_cell & !is.na(values)
  assert(
    !any(missing_inside),
    paste0(
      label, " leaves ", sum(missing_inside),
      " initially occupied cell(s) without a removal step."
    )
  )
  assert(
    !any(ranked_outside),
    paste0(
      label, " assigns removal steps to ", sum(ranked_outside),
      " cell(s) outside the initial domain."
    )
  )

  event_steps <- suppressWarnings(as.numeric(removal_events$removal_step))
  expected_counts <- suppressWarnings(as.numeric(removal_events$cells_removed))
  assert(
    all(is.finite(event_steps)) &&
      all(event_steps == floor(event_steps)) &&
      identical(as.integer(event_steps), seq_len(nrow(removal_events))),
    paste0(label, "$removal_events$removal_step must be contiguous positive integers.")
  )
  assert(
    all(is.finite(expected_counts)) &&
      all(expected_counts >= 0) &&
      all(expected_counts == floor(expected_counts)),
    paste0(label, "$removal_events$cells_removed must contain nonnegative integers.")
  )

  domain_steps <- values[initial_alive_by_cell]
  assert(
    all(domain_steps <= nrow(removal_events)),
    paste0(label, " contains a removal step absent from removal_events.")
  )
  observed_counts <- tabulate(domain_steps, nbins = nrow(removal_events))
  assert(
    identical(observed_counts, as.integer(expected_counts)),
    paste0(label, " raster cell counts do not match removal_events$cells_removed.")
  )

  values
}

validate_rankmap_event_values <- function(removal_events, value_col) {
  assert(is.data.frame(removal_events), "removal_events must be a data.frame or data.table.")
  need_cols(removal_events, c("removal_step", value_col), "removal_events")

  steps <- suppressWarnings(as.numeric(removal_events$removal_step))
  assert(
    length(steps) == nrow(removal_events) &&
      all(is.finite(steps)) &&
      all(steps > 0) &&
      all(steps == floor(steps)),
    "removal_events$removal_step must contain positive integers."
  )
  assert(
    anyDuplicated(as.integer(steps)) == 0L,
    "removal_events$removal_step contains duplicate values."
  )

  rank_values <- suppressWarnings(as.numeric(removal_events[[value_col]]))
  assert(
    length(rank_values) == nrow(removal_events) &&
      all(is.finite(rank_values)) &&
      all(rank_values >= 0) &&
      all(rank_values <= 1),
    paste0("removal_events$", value_col, " must contain finite values in [0, 1].")
  )

  list(steps = as.integer(steps), rank_values = rank_values)
}

write_removal_order_raster <- function(removal_order_by_cell,
                                       template_raster,
                                       output_path,
                                       validate_values = TRUE) {
  assert(!is.null(template_raster), "template_raster must be supplied.")
  assert(
    is.logical(validate_values) && length(validate_values) == 1L && !is.na(validate_values),
    "validate_values must be TRUE or FALSE."
  )
  if (isTRUE(validate_values)) {
    removal_order_by_cell <- validate_removal_order_values(removal_order_by_cell)
  }

  output_raster <- terra::rast(template_raster[[1L]])
  assert(
    terra::ncell(output_raster) == length(removal_order_by_cell),
    "Length of removal_order_by_cell does not match template_raster."
  )

  terra::values(output_raster) <- as.integer(removal_order_by_cell)
  terra::writeRaster(
    x = output_raster,
    filename = output_path,
    datatype = "INT4S",
    overwrite = TRUE
  )

  invisible(output_path)
}

write_rankmap_raster <- function(
  removal_order_by_cell,
  removal_events,
  template_raster,
  output_path,
  value_col = "cum_prop_cells_removed",
  validate_values = TRUE
) {
  assert(!is.null(template_raster), "template_raster must be supplied.")
  assert(
    is.logical(validate_values) && length(validate_values) == 1L && !is.na(validate_values),
    "validate_values must be TRUE or FALSE."
  )
  if (isTRUE(validate_values)) {
    removal_order_by_cell <- validate_removal_order_values(removal_order_by_cell)
  }
  event_values <- validate_rankmap_event_values(removal_events, value_col)

  output_raster <- terra::rast(template_raster[[1L]])
  assert(
    terra::ncell(output_raster) == length(removal_order_by_cell),
    "Length of removal_order_by_cell does not match template_raster."
  )

  rank_values <- rep(NA_real_, length(removal_order_by_cell))
  ranked_cells <- which(!is.na(removal_order_by_cell))

  if (length(ranked_cells)) {
    step_to_rank <- stats::setNames(
      event_values$rank_values,
      as.character(event_values$steps)
    )
    rank_values[ranked_cells] <- step_to_rank[as.character(removal_order_by_cell[ranked_cells])]
    assert(
      all(is.finite(rank_values[ranked_cells])),
      "rankmap values could not be assigned for one or more ranked cells."
    )
  }

  finite_rank_values <- rank_values[is.finite(rank_values)]
  assert(
    !length(finite_rank_values) ||
      (min(finite_rank_values) >= 0 && max(finite_rank_values) <= 1),
    "rankmap values must be in [0, 1]."
  )

  terra::values(output_raster) <- rank_values
  terra::writeRaster(
    x = output_raster,
    filename = output_path,
    overwrite = TRUE,
    wopt = list(
      datatype = "FLT4S",
      gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=6", "TILED=YES", "BIGTIFF=IF_SAFER")
    )
  )

  invisible(output_path)
}
