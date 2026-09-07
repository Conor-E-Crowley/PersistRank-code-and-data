# Stage 7.1 report-table preparation.
#
# The Stage 7.1 workflow owns loading this definition-only module after the
# scientific persistence implementation. It converts an in-memory result into
# the five small tables displayed by the Rmd; it performs no reads or writes.

# Return the manuscript-facing label for a Stage 7.1 assemblage scope.
stage71_scope_label <- function(x) {
  unname(c(
    all = "All selected species",
    mammals = "Mammals",
    birds = "Birds"
  )[as.character(x)])
}

# Prepare every table displayed by the Stage 7.1 Rmd.
#
# `result` is the validated value returned by `run_stage71()`. The returned
# list preserves the former Rmd transformations, row order, column names, and
# column types. This function neither mutates `result` nor accesses files.
stage71_report_tables <- function(result) {
  assert(is.list(result) && is.list(result$persistence),
         "Stage 7.1 report tables require a persistence result.")

  summary_table <- data.table::copy(result$persistence$summary)
  first_order <- min(summary_table$stage_order)
  last_order <- max(summary_table$stage_order)
  state_summary <- summary_table[
    stage_order %in% c(first_order, last_order)
  ]
  state_summary[, state := ifelse(
    stage_order == first_order,
    "Initial",
    "Latest"
  )]
  state_summary[, scope := stage71_scope_label(scope)]
  state_summary <- state_summary[, .(
    State = state,
    Scope = scope,
    `Cells removed (%)` = pct_cells_removed_end,
    `Area removed (%)` = pct_area_removed_end,
    Species = n_species,
    `Species with PUs` = n_species_with_pu,
    Mean = mean_persist,
    Median = median_persist,
    `10th percentile` = q10_persist,
    `90th percentile` = q90_persist,
    `Fraction > 0.5` = frac_gt_05
  )]

  largest_losses <- utils::head(
    result$persistence$changes[order(-persistence_loss, scientificName)],
    10L
  )
  largest_losses <- largest_losses[, .(
    Species = scientificName,
    Taxon = stage71_scope_label(taxon),
    SDM = sdm_method,
    Initial = initial_persistence,
    Latest = latest_persistence,
    Loss = persistence_loss,
    `Latest PUs` = latest_n_pu,
    `Latest PU area (km2)` = latest_pu_area_km2
  )]

  smallest_auc <- utils::head(
    result$persistence$changes[
      order(normalized_persistence_auc, persistence_auc, scientificName)
    ],
    12L
  )
  smallest_auc <- smallest_auc[, .(
    Species = scientificName,
    Taxon = stage71_scope_label(taxon),
    SDM = sdm_method,
    `Persistence AUC` = persistence_auc,
    `Normalized AUC` = normalized_persistence_auc,
    Initial = initial_persistence,
    Latest = latest_persistence,
    `Latest PUs` = latest_n_pu,
    `Latest PU area (km2)` = latest_pu_area_km2
  )]

  lowest_latest <- utils::head(
    result$persistence$changes[order(latest_persistence, scientificName)],
    10L
  )
  lowest_latest <- lowest_latest[, .(
    Species = scientificName,
    Taxon = stage71_scope_label(taxon),
    SDM = sdm_method,
    Initial = initial_persistence,
    Latest = latest_persistence,
    Loss = persistence_loss,
    `Latest PUs` = latest_n_pu,
    `Latest PU area (km2)` = latest_pu_area_km2
  )]

  without_pus <- result$persistence$changes[
    latest_n_pu == 0L
  ][order(taxon, scientificName), .(
    Species = scientificName,
    Taxon = stage71_scope_label(taxon),
    SDM = sdm_method,
    Initial = initial_persistence,
    Latest = latest_persistence,
    Loss = persistence_loss
  )]

  list(
    state_summary = state_summary,
    largest_losses = largest_losses,
    smallest_auc = smallest_auc,
    lowest_latest = lowest_latest,
    without_pus = without_pus
  )
}
