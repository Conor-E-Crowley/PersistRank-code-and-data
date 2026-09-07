# Compact patch-cell access, rook-neighbour indexing, and frontier state.
#
# Patch IDs and cell IDs remain in ascending canonical order. The frontier is
# initialized once, then updated only around newly empty cells so pruning never
# rescans the complete raster solely to rediscover boundary cells.

get_patch_cells_batch <- function(patch_index, patch_ids) {
  patch_ids <- sort(unique(as.integer(patch_ids)))
  patch_ids <- patch_ids[!is.na(patch_ids) & patch_ids >= 1L]

  if (!length(patch_ids) ||
      is.null(patch_index) ||
      is.null(patch_index$pid) ||
      is.null(patch_index$cell) ||
      !length(patch_index$pid)) {
    return(data.table::data.table(
      patch_id = integer(),
      cell = integer()
    ))
  }

  # The flat index is sorted by patch ID. Two binary interval lookups locate
  # every requested segment without scanning unrelated patch memberships.
  lo <- findInterval(patch_ids - 0.5, patch_index$pid) + 1L
  hi <- findInterval(patch_ids + 0.5, patch_index$pid)
  patch_present <- lo <= hi

  if (!any(patch_present)) {
    return(data.table::data.table(
      patch_id = integer(),
      cell = integer()
    ))
  }

  patch_ids <- patch_ids[patch_present]
  lo <- lo[patch_present]
  hi <- hi[patch_present]
  n_cells_by_patch <- as.integer(hi - lo + 1L)
  patch_index_rows <- sequence(n_cells_by_patch) +
    rep.int(lo - 1L, n_cells_by_patch)

  data.table::data.table(
    patch_id = rep.int(patch_ids, n_cells_by_patch),
    cell = as.integer(patch_index$cell[patch_index_rows])
  )
}

build_rook_neighbor_index <- function(rook_neighbor_pairs, n_cells, force_symmetric = TRUE) {
  if (is.null(rook_neighbor_pairs) || ncol(rook_neighbor_pairs) < 2L) {
    stop("rook_neighbor_pairs must have at least two columns.")
  }

  n_cells <- as.integer(n_cells)
  from <- as.integer(rook_neighbor_pairs[, 1L])
  to <- as.integer(rook_neighbor_pairs[, 2L])

  valid <- !is.na(from) & !is.na(to) &
    from >= 1L & from <= n_cells &
    to >= 1L & to <= n_cells &
    from != to
  from <- from[valid]
  to <- to[valid]

  # Canonical CSR-like storage is directed, deduplicated, and ordered by source
  # then target. Production rook pairs are made symmetric here once.
  if (isTRUE(force_symmetric)) {
    pairs <- data.table::data.table(
      from = c(from, to),
      to = c(to, from)
    )
  } else {
    pairs <- data.table::data.table(
      from = from,
      to = to
    )
  }
  pairs <- unique(pairs, by = c("from", "to"))
  data.table::setorder(pairs, from, to)

  neighbor_count <- tabulate(pairs$from, nbins = n_cells)
  end_index <- cumsum(neighbor_count)
  start_index <- end_index - neighbor_count + 1L
  start_index[neighbor_count == 0L] <- NA_integer_
  end_index[neighbor_count == 0L] <- NA_integer_

  list(
    from = as.integer(pairs$from),
    to = as.integer(pairs$to),
    start = as.integer(start_index),
    end = as.integer(end_index),
    n_cells = n_cells
  )
}

initialize_frontier_state <- function(alive_species_count_by_cell, rook_neighbor_index) {
  alive_by_cell <- alive_species_count_by_cell > 0L
  alive_pair <- alive_by_cell[rook_neighbor_index$from] &
    alive_by_cell[rook_neighbor_index$to]
  alive_neighbor_count_by_cell <- tabulate(
    rook_neighbor_index$from[alive_pair],
    nbins = rook_neighbor_index$n_cells
  )

  # A live cell is on the removal frontier when at least one of its four rook
  # neighbours is absent or empty.
  frontier_flag_by_cell <- alive_by_cell & alive_neighbor_count_by_cell < 4L

  list(
    alive_by_cell = alive_by_cell,
    alive_neighbor_count_by_cell = as.integer(alive_neighbor_count_by_cell),
    frontier_flag_by_cell = frontier_flag_by_cell
  )
}

get_frontier_cells <- function(frontier_state) {
  which(frontier_state$frontier_flag_by_cell)
}

update_frontier_state_after_empty_cells <- function(
  frontier_state,
  newly_empty_cells,
  rook_neighbor_index
) {
  newly_empty_cells <- unique(as.integer(newly_empty_cells))
  newly_empty_cells <- newly_empty_cells[
    !is.na(newly_empty_cells) &
      newly_empty_cells >= 1L &
      newly_empty_cells <= rook_neighbor_index$n_cells
  ]
  newly_empty_cells <- newly_empty_cells[
    frontier_state$alive_by_cell[newly_empty_cells]
  ]

  if (!length(newly_empty_cells)) {
    return(frontier_state)
  }

  frontier_state$alive_by_cell[newly_empty_cells] <- FALSE
  frontier_state$frontier_flag_by_cell[newly_empty_cells] <- FALSE

  # Each raster cell has at most four rook neighbours. Preallocation avoids
  # growing this vector inside the hot update loop.
  affected_neighbors <- integer(length(newly_empty_cells) * 4L)
  write_pos <- 0L

  for (cell in newly_empty_cells) {
    start_pos <- rook_neighbor_index$start[cell]
    if (is.na(start_pos)) {
      next
    }

    end_pos <- rook_neighbor_index$end[cell]
    neighbors <- rook_neighbor_index$to[start_pos:end_pos]
    live_neighbors <- neighbors[frontier_state$alive_by_cell[neighbors]]
    if (!length(live_neighbors)) {
      next
    }

    frontier_state$alive_neighbor_count_by_cell[live_neighbors] <-
      frontier_state$alive_neighbor_count_by_cell[live_neighbors] - 1L

    n_live <- length(live_neighbors)
    affected_neighbors[(write_pos + 1L):(write_pos + n_live)] <- live_neighbors
    write_pos <- write_pos + n_live
  }

  # Reclassify only live neighbours adjacent to the newly empty cells.
  if (write_pos > 0L) {
    affected_neighbors <- unique(affected_neighbors[seq_len(write_pos)])
    frontier_state$frontier_flag_by_cell[affected_neighbors] <-
      frontier_state$alive_by_cell[affected_neighbors] &
      frontier_state$alive_neighbor_count_by_cell[affected_neighbors] < 4L
  }

  frontier_state
}
