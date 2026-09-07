# Lightweight ecological diagnostics for Stage 6 pruning.
#
# These helpers operate only on compact score-cache entries and one selected
# cell. They do not read rasters, traverse graphs, aggregate the patch table, or
# rescore the frontier.

ecology_log_due <- function(iteration, every_iterations) {
  if (is.null(every_iterations)) {
    return(FALSE)
  }

  iteration <- as.integer(iteration)
  iteration == 1L || iteration %% as.integer(every_iterations) == 0L
}

ecology_round <- function(x, digits = 4L) {
  round(as.numeric(x), digits = as.integer(digits))
}

ecology_snapshot_from_cache <- function(
  patch_score_cache,
  stage_index,
  prune_iteration,
  current_alive_cell_count,
  initial_alive_cell_count
) {
  species_names <- sort(names(patch_score_cache))
  if (!length(species_names)) {
    stop("Cannot build ecological diagnostics from an empty score cache.")
  }

  species_rows <- lapply(species_names, function(species_name) {
    summary_row <- patch_score_cache[[species_name]]$species_summary
    if (is.null(summary_row) || nrow(summary_row) != 1L) {
      stop("Missing ecological species summary for: ", species_name)
    }
    summary_row
  })

  species_state <- data.table::rbindlist(species_rows, use.names = TRUE)
  data.table::setorder(species_state, species)

  list(
    stage = as.integer(stage_index),
    iteration = as.integer(prune_iteration),
    retained_cells_pct = 100 * as.numeric(current_alive_cell_count) /
      as.numeric(initial_alive_cell_count),
    represented_species = as.integer(nrow(species_state)),
    pu_count = as.integer(sum(species_state$pu_count)),
    single_pu_species = as.integer(sum(species_state$pu_count == 1L)),
    species_state = species_state
  )
}

format_ecology_snapshot_id <- function(snapshot) {
  paste0("stage", snapshot$stage, ":iter", snapshot$iteration)
}

format_lost_species <- function(species_names, max_names = 5L) {
  species_names <- sort(unique(as.character(species_names)))
  if (!length(species_names)) {
    return("none")
  }

  shown <- utils::head(species_names, max_names)
  suffix <- if (length(species_names) > max_names) {
    paste0(" (+", length(species_names) - max_names, " more)")
  } else {
    ""
  }
  paste0(paste(shown, collapse = "; "), suffix)
}

ecology_area_change <- function(previous_snapshot, current_snapshot, top_n = 3L) {
  if (is.null(previous_snapshot)) {
    return(list(
      status = "baseline",
      previous = NA_character_,
      current = format_ecology_snapshot_id(current_snapshot),
      largest_losses = "none",
      species_losing_all_pus_count = 0L,
      species_losing_all_pus = "none"
    ))
  }

  previous <- previous_snapshot$species_state[, .(
    species,
    previous_area_km2 = total_pu_area_km2,
    previous_pu_count = pu_count
  )]
  current <- current_snapshot$species_state[, .(
    species,
    current_area_km2 = total_pu_area_km2,
    current_pu_count = pu_count
  )]

  comparison <- merge(previous, current, by = "species", all = TRUE)
  comparison[is.na(previous_area_km2), `:=`(
    previous_area_km2 = 0,
    previous_pu_count = 0L
  )]
  comparison[is.na(current_area_km2), `:=`(
    current_area_km2 = 0,
    current_pu_count = 0L
  )]
  comparison[, area_loss_km2 := pmax(previous_area_km2 - current_area_km2, 0)]
  comparison[, area_loss_pct := data.table::fifelse(
    previous_area_km2 > 0,
    100 * area_loss_km2 / previous_area_km2,
    0
  )]

  ranked <- comparison[area_loss_km2 > 0]
  data.table::setorder(ranked, -area_loss_km2, species)
  ranked <- utils::head(ranked, as.integer(top_n))

  largest_losses <- if (nrow(ranked)) {
    paste(
      sprintf(
        "%s:%.3f km2 (%.2f%%)",
        ranked$species,
        ranked$area_loss_km2,
        ranked$area_loss_pct
      ),
      collapse = "; "
    )
  } else {
    "none"
  }

  lost_species <- comparison[
    previous_pu_count > 0L & current_pu_count == 0L,
    species
  ]

  list(
    status = "change",
    previous = format_ecology_snapshot_id(previous_snapshot),
    current = format_ecology_snapshot_id(current_snapshot),
    largest_losses = largest_losses,
    species_losing_all_pus_count = as.integer(length(lost_species)),
    species_losing_all_pus = format_lost_species(lost_species)
  )
}

ecology_extreme_pus <- function(patch_score_cache) {
  species_names <- sort(names(patch_score_cache))
  maximum_rows <- lapply(species_names, function(species_name) {
    patch_score_cache[[species_name]]$maximum_value_pu
  })
  minimum_rows <- lapply(species_names, function(species_name) {
    patch_score_cache[[species_name]]$minimum_value_pu
  })

  if (any(vapply(maximum_rows, is.null, logical(1L))) ||
      any(vapply(minimum_rows, is.null, logical(1L)))) {
    stop("Ecological PU extrema are missing from the score cache.")
  }

  maximum_rows <- data.table::rbindlist(maximum_rows, use.names = TRUE)
  minimum_rows <- data.table::rbindlist(minimum_rows, use.names = TRUE)

  list(
    # Species names were assembled in lexical order and each cached row already
    # uses the lowest PU ID to break within-species ties. which.max()/which.min()
    # therefore provide the complete deterministic tie rule without sorting.
    most = maximum_rows[which.max(maximum_rows$log_score)],
    least = minimum_rows[which.min(minimum_rows$log_score)]
  )
}

ecology_cutoff_cell <- function(
  selected_cells,
  frontier_cells,
  frontier_scores,
  patch_score_cache,
  patch_id_by_species_env
) {
  selected_position <- fastmatch::fmatch(selected_cells, frontier_cells)
  if (anyNA(selected_position)) {
    stop("Selected ecology-diagnostic cells are absent from the frontier.")
  }

  selected_frontier_scores <- frontier_scores[selected_position]
  selected_log_loss <- -selected_frontier_scores
  selected_log10_loss <- selected_log_loss / log(10)

  # The cutoff cell has the greatest admitted loss, equivalently the smallest
  # stored negative-log frontier score. Break exact ties by lower cell ID.
  cutoff_order <- order(-selected_log_loss, selected_cells)
  cutoff_index <- cutoff_order[[1L]]
  cutoff_cell <- as.integer(selected_cells[[cutoff_index]])
  cutoff_log_loss <- as.numeric(selected_log_loss[[cutoff_index]])

  contribution_rows <- lapply(sort(names(patch_score_cache)), function(species_name) {
    cache_entry <- patch_score_cache[[species_name]]
    patch_ids <- get(species_name, envir = patch_id_by_species_env, inherits = FALSE)
    patch_id <- patch_ids[[cutoff_cell]]
    if (is.na(patch_id) || patch_id < 1L ||
        patch_id > length(cache_entry$score_by_patch_id)) {
      return(NULL)
    }

    log_contribution <- cache_entry$score_by_patch_id[[patch_id]]
    pu_id <- cache_entry$pu_id_by_patch_id[[patch_id]]
    if (!is.finite(log_contribution) || is.na(pu_id)) {
      return(NULL)
    }

    data.table::data.table(
      species = species_name,
      pu_id = as.integer(pu_id),
      log_contribution = as.numeric(log_contribution)
    )
  })
  contribution_rows <- Filter(Negate(is.null), contribution_rows)
  if (!length(contribution_rows)) {
    stop("No species contribution was found for the selected cutoff cell.")
  }

  contributions <- data.table::rbindlist(contribution_rows, use.names = TRUE)
  driver_order <- order(
    -contributions$log_contribution,
    contributions$species,
    contributions$pu_id
  )
  driver <- contributions[driver_order[[1L]]]

  driver_cache <- patch_score_cache[[driver$species[[1L]]]]
  driver_pu <- driver_cache$pu_diagnostics[
    pu_id == driver$pu_id[[1L]]
  ]
  if (nrow(driver_pu) != 1L) {
    stop("The cutoff-cell driver does not resolve to one PU diagnostic row.")
  }

  list(
    selected_log10_loss_median = stats::median(selected_log10_loss),
    selected_log10_loss_cutoff = max(selected_log10_loss),
    cutoff_cell = cutoff_cell,
    driver = driver_pu,
    contribution_share = exp(driver$log_contribution[[1L]] - cutoff_log_loss)
  )
}

log_ecology_diagnostics <- function(
  current_snapshot,
  previous_snapshot,
  cutoff,
  extremes
) {
  runtime_log_event(
    "ecology_state",
    stage = current_snapshot$stage,
    iter = current_snapshot$iteration,
    retained_cells_pct = ecology_round(current_snapshot$retained_cells_pct, 1L),
    represented_species = current_snapshot$represented_species,
    pus = current_snapshot$pu_count,
    single_pu_species = current_snapshot$single_pu_species
  )

  driver <- cutoff$driver
  runtime_log_event(
    "ecology_score",
    stage = current_snapshot$stage,
    iter = current_snapshot$iteration,
    selected_log10_loss_median = ecology_round(cutoff$selected_log10_loss_median),
    selected_log10_loss_cutoff = ecology_round(cutoff$selected_log10_loss_cutoff),
    cutoff_cell = cutoff$cutoff_cell,
    driver = paste0(driver$species[[1L]], "|PU", driver$pu_id[[1L]]),
    contribution_share = ecology_round(cutoff$contribution_share),
    pu_persistence = ecology_round(driver$pu_persistence[[1L]]),
    species_persistence = ecology_round(driver$species_persistence[[1L]]),
    a_over_a0 = ecology_round(
      driver$pu_area_km2[[1L]] / driver$area_threshold_km2[[1L]]
    ),
    log10_sensitivity = ecology_round(
      driver$log_local_sensitivity[[1L]] / log(10)
    ),
    redundancy_multiplier = ecology_round(
      exp(driver$log_redundancy_multiplier[[1L]])
    )
  )

  most <- extremes$most
  least <- extremes$least
  runtime_log_event(
    "ecology_pu_extremes",
    stage = current_snapshot$stage,
    iter = current_snapshot$iteration,
    most = paste0(most$species[[1L]], "|PU", most$pu_id[[1L]]),
    most_log10_m = ecology_round(most$log_score[[1L]] / log(10)),
    most_pu_persistence = ecology_round(most$pu_persistence[[1L]]),
    most_species_persistence = ecology_round(most$species_persistence[[1L]]),
    most_a_over_a0 = ecology_round(
      most$pu_area_km2[[1L]] / most$area_threshold_km2[[1L]]
    ),
    most_log10_sensitivity = ecology_round(
      most$log_local_sensitivity[[1L]] / log(10)
    ),
    most_redundancy = ecology_round(exp(most$log_redundancy_multiplier[[1L]])),
    least = paste0(least$species[[1L]], "|PU", least$pu_id[[1L]]),
    least_log10_m = ecology_round(least$log_score[[1L]] / log(10)),
    least_pu_persistence = ecology_round(least$pu_persistence[[1L]]),
    least_species_persistence = ecology_round(least$species_persistence[[1L]]),
    least_a_over_a0 = ecology_round(
      least$pu_area_km2[[1L]] / least$area_threshold_km2[[1L]]
    ),
    least_log10_sensitivity = ecology_round(
      least$log_local_sensitivity[[1L]] / log(10)
    ),
    least_redundancy = ecology_round(exp(least$log_redundancy_multiplier[[1L]]))
  )

  area_change <- ecology_area_change(previous_snapshot, current_snapshot)
  runtime_log_event(
    "ecology_area_change",
    stage = current_snapshot$stage,
    iter = current_snapshot$iteration,
    status = area_change$status,
    previous = area_change$previous,
    current = area_change$current,
    largest_losses = area_change$largest_losses,
    species_losing_all_pus_count = area_change$species_losing_all_pus_count,
    species_losing_all_pus = area_change$species_losing_all_pus
  )

  invisible(current_snapshot)
}
