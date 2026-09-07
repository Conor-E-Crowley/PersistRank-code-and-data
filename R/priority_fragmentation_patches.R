# Flat patch-to-cell indexes used by fragmentation and distance repair.
#
# Each species index is sorted by patch ID and raster-cell ID. Production
# updates use the compiled splice kernel to replace only changed/stale
# segments. Resume reconstruction is owned separately by the checkpoint
# module.

# Index one species' affected patch rows by PU.

index_affected_patch_rows_by_pu <- function(patch_rows, affected_pu_ids) {
  patch_rows <- data.table::as.data.table(patch_rows)
  if (!"pu_id" %in% names(patch_rows)) {
    stop("patch_rows must contain pu_id.", call. = FALSE)
  }

  affected_pu_ids <- as.integer(affected_pu_ids)
  if (anyNA(affected_pu_ids) || any(affected_pu_ids < 1L) ||
      anyDuplicated(affected_pu_ids)) {
    stop("affected_pu_ids must be unique positive integers.", call. = FALSE)
  }
  if (!length(affected_pu_ids)) {
    return(stats::setNames(vector("list", 0L), character()))
  }

  row_is_affected <- patch_rows$pu_id %in% affected_pu_ids
  affected_rows <- which(row_is_affected)
  rows_by_pu <- split(
    affected_rows,
    as.character(patch_rows$pu_id[affected_rows])
  )
  output <- stats::setNames(
    rep(list(integer()), length(affected_pu_ids)),
    as.character(affected_pu_ids)
  )
  represented_names <- intersect(names(rows_by_pu), names(output))
  output[represented_names] <- lapply(rows_by_pu[represented_names], as.integer)
  output
}


# Update the flat sorted compact index from local patch changes in compiled code.
empty_fragmentation_index_counts <- function() {
  stats::setNames(numeric(7L), c(
    "species_indexes_updated", "input_index_entries",
    "changed_origin_segments", "replacement_entries",
    "output_index_entries", "removed_entries", "relabel_entries"
  ))
}

empty_fragmentation_index_timing <- function() {
  stats::setNames(numeric(4L), c(
    "flatten_seconds", "kernel_seconds",
    "result_packaging_seconds", "total_seconds"
  ))
}

update_species_patch_index_incremental <- function(
  patch_index,
  replacement_rows_by_origin,
  expected_patch_ids,
  pre_update_patch_ids,
  return_diagnostics = FALSE
) {
  update_started <- proc.time()[["elapsed"]]
  if (!exists("stage6_splice_patch_index_cpp", mode = "function", inherits = TRUE)) {
    stop("The compiled Stage 6 compact-index kernel is not loaded.", call. = FALSE)
  }
  if (is.null(patch_index) || is.null(patch_index$pid) || is.null(patch_index$cell)) {
    stop("patch_index must contain pid and cell vectors.", call. = FALSE)
  }

  flatten_started <- proc.time()[["elapsed"]]
  replacement_rows_by_origin <- replacement_rows_by_origin[!vapply(
    replacement_rows_by_origin,
    is.null,
    logical(1L)
  )]
  replacement_origin_ids <- if (length(replacement_rows_by_origin)) {
    suppressWarnings(as.integer(names(replacement_rows_by_origin)))
  } else {
    integer()
  }

  flattened <- lapply(seq_along(replacement_rows_by_origin), function(i) {
    rows <- data.table::as.data.table(replacement_rows_by_origin[[i]])
    if (!all(c("pid", "cell") %in% names(rows))) {
      stop("Compact-index replacement rows must contain pid and cell.", call. = FALSE)
    }
    data.table::data.table(
      origin = rep.int(replacement_origin_ids[[i]], nrow(rows)),
      pid = as.integer(rows$pid),
      cell = as.integer(rows$cell)
    )
  })
  flat <- if (length(flattened)) {
    data.table::rbindlist(flattened, use.names = TRUE)
  } else {
    data.table::data.table(origin = integer(), pid = integer(), cell = integer())
  }
  flatten_seconds <- proc.time()[["elapsed"]] - flatten_started

  kernel_started <- proc.time()[["elapsed"]]
  result <- stage6_splice_patch_index_cpp(
    old_pid = as.integer(patch_index$pid),
    old_cell = as.integer(patch_index$cell),
    expected_patch_ids_input = as.integer(expected_patch_ids),
    pre_update_patch_ids_input = as.integer(pre_update_patch_ids),
    replacement_origin_ids_input = as.integer(replacement_origin_ids),
    replacement_origin = flat$origin,
    replacement_pid = flat$pid,
    replacement_cell = flat$cell
  )
  kernel_seconds <- proc.time()[["elapsed"]] - kernel_started

  packaging_started <- proc.time()[["elapsed"]]
  output <- list(
    index = if (length(result$pid)) {
      list(pid = as.integer(result$pid), cell = as.integer(result$cell))
    } else {
      NULL
    },
    removed_cells = as.integer(result$removed_cells),
    relabel_rows = data.table::data.table(
      cell = as.integer(result$relabel_cells),
      patch_id = as.integer(result$relabel_patch_ids)
    )
  )
  packaging_seconds <- proc.time()[["elapsed"]] - packaging_started
  if (!isTRUE(return_diagnostics)) return(output)

  output$diagnostic_counts <- stats::setNames(c(
    1,
    length(patch_index$pid),
    length(replacement_origin_ids),
    nrow(flat),
    length(result$pid),
    length(result$removed_cells),
    length(result$relabel_cells)
  ), names(empty_fragmentation_index_counts()))
  output$diagnostic_timing_seconds <- stats::setNames(c(
    flatten_seconds,
    kernel_seconds,
    packaging_seconds,
    proc.time()[["elapsed"]] - update_started
  ), names(empty_fragmentation_index_timing()))
  output
}
