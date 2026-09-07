# Persistence-loss scoring for pruning.
#
# PU and patch scores remain on the log scale. Frontier accumulation visits the
# immutable sparse cell-to-species membership index in canonical species order,
# preserving floating-point reduction order and lower-cell-ID tie breaking.

compute_pu_log_scores <- function(
  pu_area_table,
  species_params,
  include_diagnostics = FALSE
) {
  # Join species-level parameters onto current PU-area rows.
  pu_data <- species_params[pu_area_table, on = .(species)]

  # Read current PU areas in km^2.
  pu_area_km2 <- pu_data$pu_area_km2

  # Read species density values.
  species_density <- pu_data$density

  # Read fitted Gompertz alpha parameters.
  persistence_a <- pu_data$a_pred

  # Read fitted Gompertz beta parameters.
  persistence_b <- pu_data$b_pred

  # Read the canonical species-specific minimum population-area thresholds.
  pu_area_threshold_km2 <- pu_data$min_population_area_km2

  # Read species names aligned with PU rows.
  species_name <- pu_data$species

  # Read PU IDs aligned with PU rows.
  pu_id <- pu_data$pu_id

  # Compute each PU's area above the quasi-extinction area threshold.
  area_above_threshold <- pu_area_km2 - pu_area_threshold_km2
  if (any(!is.finite(area_above_threshold) | area_above_threshold <= 0)) {
    bad <- which(!is.finite(area_above_threshold) | area_above_threshold <= 0)
    stop(
      "PU area must remain strictly above its population-area threshold before scoring. First affected: ",
      species_name[bad[[1L]]], "|", pu_id[bad[[1L]]],
      call. = FALSE
    )
  }

  # Take log area-above-threshold once for reuse.
  log_area_above_threshold <- log(area_above_threshold)

  # Compute log(a * density^(-b)).
  log_a_density_term <-
    log(persistence_a) - persistence_b * log(species_density)

  # Compute log(a * density^(-b) * area_above_threshold^(-b)).
  log_inner_term <-
    log_a_density_term -
    persistence_b * log_area_above_threshold

  # Convert the inner term back to ordinary scale.
  inner_term <- exp(log_inner_term)

  # Compute current PU persistence.
  pu_persistence <- exp(-inner_term)

  # Compute log(1 - exp(-inner_term)) without first rounding persistence to
  # exactly one when inner_term is positive but smaller than machine precision.
  log_one_minus_persistence <- log(-expm1(-inner_term))

  # Sum log failure probabilities across PUs within each species.
  species_sum_log_failure <- tapply(log_one_minus_persistence, species_name, sum)

  # Remove each PU's own failure term to get the product over all other PUs.
  log_other_pu_term <-
    species_sum_log_failure[species_name] - log_one_minus_persistence

  # Compute log derivative of PU persistence with respect to PU area.
  log_persistence_derivative <-
    log_a_density_term +
    log(persistence_b) -
    (persistence_b + 1) * log_area_above_threshold -
    inner_term

  # Combine derivative and redundancy terms into the final PU log marginal score.
  log_score <- log_persistence_derivative + log_other_pu_term

  # Stop if any PU score is non-finite.
  if (any(!is.finite(log_score))) {
    stop("compute_pu_log_scores() produced non-finite PU scores.")
  }

  # Build the compact score table used by the pruning cache.
  score_table <- data.table::data.table(
    species = species_name,
    pu_id = pu_id,
    log_score = log_score
  )

  # Ordinary scoring retains the original minimal three-column contract.
  if (!isTRUE(include_diagnostics)) {
    return(score_table)
  }

  # Species persistence is one minus the product of all current PU failures.
  species_persistence <- -expm1(species_sum_log_failure[species_name])
  if (any(!is.finite(species_persistence) |
          species_persistence < 0 | species_persistence > 1)) {
    stop("compute_pu_log_scores() produced invalid species persistence values.")
  }

  # Retain only intermediates already calculated for scientific scoring. These
  # compact PU rows support optional ecological console diagnostics without
  # rescoring PUs or expanding the large patch table.
  score_table[, `:=`(
    pu_area_km2 = as.numeric(pu_area_km2),
    area_threshold_km2 = as.numeric(pu_area_threshold_km2),
    pu_persistence = as.numeric(pu_persistence),
    species_persistence = as.numeric(species_persistence),
    log_local_sensitivity = as.numeric(log_persistence_derivative),
    log_redundancy_multiplier = as.numeric(log_other_pu_term)
  )]

  score_table
}



# Update the patch-score cache for dirty species only.

update_patch_score_cache <- function(
  patch_table,
  species_params,
  patch_score_cache = NULL,
  dirty_species = NULL,
  include_ecology = FALSE
) {
  # Stop if the patch table is unexpectedly empty.
  if (!nrow(patch_table)) {
    stop("update_patch_score_cache() received an empty patch_table.")
  }

  # Normalize an explicitly supplied dirty-species set once.
  requested_dirty_species <- if (is.null(dirty_species)) {
    NULL
  } else {
    sort(unique(as.character(dirty_species)))
  }

  # Initialization and conservative refreshes must inspect every represented
  # species. Ordinary incremental refreshes instead scan the species column
  # only once to extract the requested dirty rows.
  if (is.null(patch_score_cache)) {
    patch_score_cache <- list()
    current_species <- sort(unique(as.character(patch_table$species)))
    dirty_species <- current_species
    dirty_patch_table <- patch_table
  } else if (is.null(requested_dirty_species)) {
    current_species <- sort(unique(as.character(patch_table$species)))
    stale_species <- setdiff(names(patch_score_cache), current_species)
    if (length(stale_species)) {
      patch_score_cache[stale_species] <- NULL
    }
    dirty_species <- current_species
    dirty_patch_table <- patch_table
  } else {
    dirty_row <- data.table::`%chin%`(
      as.character(patch_table$species),
      requested_dirty_species
    )
    dirty_patch_table <- patch_table[dirty_row]
    dirty_species <- sort(unique(as.character(dirty_patch_table$species)))

    # Requested dirty species absent from the current table have gone extinct;
    # remove only those known stale entries without rescanning all species.
    extinct_dirty_species <- setdiff(requested_dirty_species, dirty_species)
    if (length(extinct_dirty_species)) {
      patch_score_cache[extinct_dirty_species] <- NULL
    }
  }

  # Recompute score vectors only for dirty species.
  if (length(dirty_species)) {
    # Stop if dirty species were requested but no rows exist for them.
    if (!nrow(dirty_patch_table)) {
      stop("Dirty species were requested, but dirty_patch_table is empty.")
    }

    # Aggregate dirty patch rows to current PU areas.
    dirty_pu_area_table <- dirty_patch_table[
      ,
      list(pu_area_km2 = sum(patch_area_km2)),
      by = .(species, pu_id)
    ]

    # Compute current PU-level log scores for dirty species.
    dirty_pu_score_table <- compute_pu_log_scores(
      pu_area_table = dirty_pu_area_table,
      species_params = species_params,
      include_diagnostics = include_ecology
    )

    # Attach each dirty patch to its current PU log score.
    dirty_patch_score_table <- dirty_patch_table[
      ,
      list(species, patch_id, pu_id)
    ][
      dirty_pu_score_table[, .(species, pu_id, log_score)],
      on = .(species, pu_id),
      nomatch = 0L
    ][
      ,
      list(
        species = as.character(species),
        patch_id = as.integer(patch_id),
        pu_id = as.integer(pu_id),
        score = log_score
      )
    ]

    # Stop if no patch scores were produced for dirty species.
    if (!nrow(dirty_patch_score_table)) {
      stop("update_patch_score_cache() produced no scores for dirty species.")
    }

    # Identify dirty species that actually received score rows.
    scored_dirty_species <- sort(unique(dirty_patch_score_table$species))

    # Identify dirty species missing from the score table.
    missing_dirty_species <- setdiff(dirty_species, scored_dirty_species)

    # Stop if any dirty surviving species did not receive scores.
    if (length(missing_dirty_species)) {
      stop(
        "No patch scores were produced for dirty species: ",
        paste(missing_dirty_species, collapse = ", ")
      )
    }

    # Stop if patch IDs cannot be used as positive integer vector indices.
    if (any(is.na(dirty_patch_score_table$patch_id) |
            dirty_patch_score_table$patch_id < 1L)) {
      stop("Patch IDs must be positive integers for direct patch-score indexing.")
    }

    # Stop if any patch score is non-finite.
    if (any(!is.finite(dirty_patch_score_table$score))) {
      stop("update_patch_score_cache() produced non-finite patch scores.")
    }

    # Split row numbers by species for per-species score-vector construction.
    row_index_by_species <- split(
      seq_len(nrow(dirty_patch_score_table)),
      dirty_patch_score_table$species
    )

    # Rebuild one direct patch_id -> score vector per dirty species.
    for (species_name in names(row_index_by_species)) {
      # Read row positions for this species.
      row_index <- row_index_by_species[[species_name]]

      # Read patch IDs for this species.
      patch_ids <- dirty_patch_score_table$patch_id[row_index]

      # Read log scores for this species.
      scores <- dirty_patch_score_table$score[row_index]

      # Read PU IDs aligned with patch IDs for direct pre-step lookup.
      pu_ids <- dirty_patch_score_table$pu_id[row_index]

      # Stop if the score table contains duplicate patch IDs for one species.
      if (anyDuplicated(patch_ids)) {
        stop("Duplicate patch IDs found while building score cache for: ", species_name)
      }

      # Allocate direct lookup vector, filling absent patch IDs with -Inf.
      score_by_patch_id <- rep(-Inf, max(patch_ids))

      # Store finite patch scores at their patch-ID positions.
      score_by_patch_id[patch_ids] <- scores

      # Build the direct patch_id -> pu_id lookup used by pruning repair. This
      # compact integer vector replaces a full patch-table snapshot each step.
      pu_id_by_patch_id <- rep(NA_integer_, max(patch_ids))
      pu_id_by_patch_id[patch_ids] <- pu_ids

      # Store this species' direct lookup vectors in the cache.
      cache_entry <- list(
        score_by_patch_id = score_by_patch_id,
        pu_id_by_patch_id = pu_id_by_patch_id
      )

      if (isTRUE(include_ecology)) {
        pu_diagnostics <- dirty_pu_score_table[
          species == species_name,
          .(
            species = as.character(species),
            pu_id = as.integer(pu_id),
            log_score = as.numeric(log_score),
            pu_area_km2 = as.numeric(pu_area_km2),
            area_threshold_km2 = as.numeric(area_threshold_km2),
            pu_persistence = as.numeric(pu_persistence),
            species_persistence = as.numeric(species_persistence),
            log_local_sensitivity = as.numeric(log_local_sensitivity),
            log_redundancy_multiplier = as.numeric(log_redundancy_multiplier)
          )
        ]

        if (!nrow(pu_diagnostics) || any(!is.finite(pu_diagnostics$log_score))) {
          stop("Missing or invalid PU diagnostics for species: ", species_name)
        }

        # Determine per-species extrema with deterministic PU-ID tie-breaking.
        maximum_score <- max(pu_diagnostics$log_score)
        minimum_score <- min(pu_diagnostics$log_score)
        maximum_candidates <- which(pu_diagnostics$log_score == maximum_score)
        minimum_candidates <- which(pu_diagnostics$log_score == minimum_score)
        maximum_index <- maximum_candidates[[which.min(
          pu_diagnostics$pu_id[maximum_candidates]
        )]]
        minimum_index <- minimum_candidates[[which.min(
          pu_diagnostics$pu_id[minimum_candidates]
        )]]

        cache_entry$pu_diagnostics <- pu_diagnostics
        cache_entry$species_summary <- data.table::data.table(
          species = species_name,
          pu_count = as.integer(nrow(pu_diagnostics)),
          total_pu_area_km2 = as.numeric(sum(pu_diagnostics$pu_area_km2))
        )
        cache_entry$maximum_value_pu <- pu_diagnostics[maximum_index]
        cache_entry$minimum_value_pu <- pu_diagnostics[minimum_index]
      }

      patch_score_cache[[species_name]] <- cache_entry
    }
  }

  # Stop if the score cache is unexpectedly empty.
  if (!length(patch_score_cache)) {
    stop("update_patch_score_cache() produced an empty patch-score cache.")
  }

  # Return the updated score cache.
  patch_score_cache
}



# Dispatch frontier scoring to the validated compiled Stage 6 kernel.

empty_frontier_scoring_counts <- function() {
  stats::setNames(numeric(8L), c(
    "score_calls", "frontier_cells", "active_species",
    "species_frontier_probes", "sparse_membership_probes",
    "dense_probes_avoided", "finite_contributions",
    "scored_frontier_cells"
  ))
}

empty_frontier_scoring_timing <- function() {
  stats::setNames(numeric(6L), c(
    "input_setup_seconds", "kernel_seconds", "result_finalize_seconds",
    "kernel_scan_seconds", "kernel_accumulation_seconds",
    "kernel_output_seconds"
  ))
}

score_frontier_cells <- function(
  frontier_cells,
  patch_scores_by_species,
  patch_id_by_species_env,
  cell_species_index = NULL,
  return_diagnostics = FALSE
) {
  if (!exists("stage6_score_frontier_cpp", mode = "function", inherits = TRUE)) {
    stop("The compiled Stage 6 frontier kernel is not loaded.", call. = FALSE)
  }
  if (!length(frontier_cells)) {
    scores <- numeric()
    if (!isTRUE(return_diagnostics)) return(scores)
    return(list(
      scores = scores,
      counts = stats::setNames(
        c(1, 0, length(patch_scores_by_species), 0, 0, 0, 0, 0),
        names(empty_frontier_scoring_counts())
      ),
      timing_seconds = empty_frontier_scoring_timing()
    ))
  }

  input_started <- proc.time()[["elapsed"]]
  species_names <- names(patch_scores_by_species)
  patch_ids_by_species <- lapply(species_names, function(species_name) {
    get(species_name, envir = patch_id_by_species_env, inherits = FALSE)
  })
  scores_by_species <- lapply(species_names, function(species_name) {
    patch_scores_by_species[[species_name]]$score_by_patch_id
  })
  input_setup_seconds <- proc.time()[["elapsed"]] - input_started
  kernel_started <- proc.time()[["elapsed"]]
  accumulator <- if (is.null(cell_species_index)) {
    stage6_score_frontier_cpp(
      frontier_cells = as.integer(frontier_cells),
      patch_ids_by_species = patch_ids_by_species,
      scores_by_species = scores_by_species,
      return_diagnostics = return_diagnostics
    )
  } else {
    if (!is.list(cell_species_index) ||
        !all(c("offsets", "species_ids", "species_names") %in%
          names(cell_species_index))) {
      stop("cell_species_index has an invalid runtime contract.")
    }
    active_species_ids <- match(species_names, cell_species_index$species_names)
    if (anyNA(active_species_ids)) {
      stop("Frontier score species are missing from cell_species_index.")
    }
    stage6_score_frontier_sparse_cpp(
      frontier_cells = as.integer(frontier_cells),
      cell_offsets = cell_species_index$offsets,
      species_ids_input = cell_species_index$species_ids,
      canonical_species_count = length(cell_species_index$species_names),
      active_species_ids = as.integer(active_species_ids),
      patch_ids_by_species = patch_ids_by_species,
      scores_by_species = scores_by_species,
      return_diagnostics = return_diagnostics
    )
  }
  kernel_seconds <- proc.time()[["elapsed"]] - kernel_started
  finalize_started <- proc.time()[["elapsed"]]
  frontier_scores <- rep(-Inf, length(frontier_cells))
  scored_position <- accumulator$has_score &
    accumulator$sum_exp_terms > 0 &
    is.finite(accumulator$max_log_term)
  frontier_scores[scored_position] <- -(
    accumulator$max_log_term[scored_position] +
      log(accumulator$sum_exp_terms[scored_position])
  )
  if (any(!is.finite(frontier_scores))) {
    bad <- which(!is.finite(frontier_scores))
    stop(
      "score_frontier_cells() produced non-finite frontier scores for ",
      length(bad),
      " frontier cells. This indicates missing/stale patch-score lookups. ",
      "Example raster cell ids: ",
      paste(head(frontier_cells[bad], 10L), collapse = ", ")
    )
  }
  result_finalize_seconds <- proc.time()[["elapsed"]] - finalize_started
  if (!isTRUE(return_diagnostics)) return(frontier_scores)

  counts <- stats::setNames(
    as.numeric(accumulator$diagnostic_counts),
    names(empty_frontier_scoring_counts())
  )
  kernel_timing <- as.numeric(accumulator$diagnostic_timing_seconds)
  timing_seconds <- stats::setNames(
    c(
      input_setup_seconds,
      kernel_seconds,
      result_finalize_seconds,
      kernel_timing
    ),
    names(empty_frontier_scoring_timing())
  )
  if (any(!is.finite(c(counts, timing_seconds))) ||
      any(c(counts, timing_seconds) < 0)) {
    stop("Frontier scoring diagnostics must be finite and nonnegative.")
  }
  list(
    scores = frontier_scores,
    counts = counts,
    timing_seconds = timing_seconds
  )
}


# Build the immutable cell-to-possible-species runtime index.

build_cell_species_index <- function(
  patch_cell_index_by_species_env,
  species_names,
  n_cells,
  alive_species_count_by_cell = NULL
) {
  if (!exists(
    "stage6_build_cell_species_index_cpp",
    mode = "function",
    inherits = TRUE
  )) {
    stop("The compiled Stage 6 cell-species index kernel is not loaded.")
  }
  species_names <- as.character(species_names)
  if (!length(species_names) || anyNA(species_names) ||
      any(!nzchar(species_names)) || anyDuplicated(species_names)) {
    stop("Cell-species index requires unique nonblank species names.")
  }
  if (!identical(species_names, sort(species_names))) {
    stop("Cell-species index species names must use canonical sorted order.")
  }
  n_cells <- as.integer(n_cells)
  if (length(n_cells) != 1L || is.na(n_cells) || n_cells < 1L) {
    stop("Cell-species index requires a positive cell count.")
  }
  cells_by_species <- lapply(species_names, function(species_name) {
    patch_index <- get(
      species_name,
      envir = patch_cell_index_by_species_env,
      inherits = FALSE
    )
    if (is.null(patch_index)) return(integer())
    as.integer(patch_index$cell)
  })
  built <- stage6_build_cell_species_index_cpp(
    cells_by_species = cells_by_species,
    n_cells = n_cells,
    use_raw_species_ids = length(species_names) <= 255L
  )
  if (!is.null(alive_species_count_by_cell) &&
      !identical(
        as.integer(built$membership_count_by_cell),
        as.integer(alive_species_count_by_cell)
      )) {
    stop("Cell-species memberships differ from alive_species_count_by_cell.")
  }
  list(
    offsets = as.integer(built$offsets),
    species_ids = built$species_ids,
    species_names = species_names,
    entry_count = as.numeric(built$entry_count)
  )
}



# Choose which frontier cells to remove.

choose_frontier_cells_to_remove <- function(
  frontier_cells,
  frontier_scores,
  cells_to_remove_per_step
) {
  # Count available frontier cells.
  n_frontier <- length(frontier_cells)

  # Stop if the requested removal batch size is invalid.
  if (cells_to_remove_per_step < 1L) {
    stop("cells_to_remove_per_step must be >= 1.")
  }

  # If the request is at least the full frontier, return the full frontier deterministically.
  if (cells_to_remove_per_step >= n_frontier) {
    # Order higher score first, then lower cell index.
    ordered_positions <- order(-frontier_scores, frontier_cells)

    # Return ordered frontier cell IDs.
    return(frontier_cells[ordered_positions])
  }

  # Find the cutoff score for the top requested number of frontier cells.
  cutoff_score <- Rfast::nth(
    frontier_scores,
    k = cells_to_remove_per_step,
    descending = TRUE,
    index.return = FALSE
  )

  # Keep frontier positions whose score is at least the cutoff score.
  selected_positions <- which(frontier_scores >= cutoff_score)

  # If ties over-selected cells, break ties deterministically.
  if (length(selected_positions) > cells_to_remove_per_step) {
    # Order tied candidates by higher score, then lower cell index.
    tie_break_order <- order(
      -frontier_scores[selected_positions],
      frontier_cells[selected_positions]
    )

    # Keep exactly the requested number of selected positions.
    selected_positions <-
      selected_positions[tie_break_order][seq_len(cells_to_remove_per_step)]
  }

  # Convert selected positions back to global raster-cell IDs.
  frontier_cells[selected_positions]
}



# Summarize removed area by species-specific patch.

summarize_removed_patch_area <- function(
  selected_cells,
  removed_cell_area,
  species_names,
  patch_id_by_species_env,
  patch_score_cache
) {
  # Combine one per-species summary table into a single data.table.
  data.table::rbindlist(
    # Loop over all species present in the current pruning stage.
    lapply(species_names, function(species_name) {
      # A missing cache entry means that the species has no surviving patches.
      # Its dense vector may still contain deferred stale IDs, but those IDs are
      # not part of the current ecological state.
      cache_entry <- patch_score_cache[[species_name]]
      if (is.null(cache_entry)) {
        return(NULL)
      }
      pu_id_by_patch_id <- cache_entry$pu_id_by_patch_id
      if (is.null(pu_id_by_patch_id) || !is.atomic(pu_id_by_patch_id)) {
        stop("Patch-score cache is missing the live-PU lookup for species: ", species_name)
      }

      # Read this species' patch IDs in the selected cells.
      patch_ids <- get(
        species_name,
        envir = patch_id_by_species_env,
        inherits = FALSE
      )[selected_cells]

      # Dense patch IDs can remain temporarily at cells belonging to patches
      # removed by an earlier threshold cascade. Only IDs with a current,
      # positive PU lookup are live in this pre-removal state.
      valid <- !is.na(patch_ids) &
        patch_ids >= 1L &
        patch_ids <= length(pu_id_by_patch_id)
      if (any(valid)) {
        selected_pu_ids <- pu_id_by_patch_id[patch_ids[valid]]
        live <- !is.na(selected_pu_ids) &
          is.finite(selected_pu_ids) &
          selected_pu_ids >= 1 &
          selected_pu_ids == floor(selected_pu_ids)
        valid[which(valid)] <- live
      }

      # Return no rows if this species is absent from all selected cells.
      if (!any(valid)) {
        return(NULL)
      }

      # Sum removed cell area by patch ID for this species.
      area_sum <- rowsum(
        removed_cell_area[valid],
        patch_ids[valid],
        reorder = FALSE
      )

      # Return one row per affected species-specific patch.
      data.table::data.table(
        species = species_name,
        patch_id = as.integer(rownames(area_sum)),
        area_removed = as.numeric(area_sum[, 1])
      )
    }),
    use.names = TRUE,
    fill = TRUE
  )
}
