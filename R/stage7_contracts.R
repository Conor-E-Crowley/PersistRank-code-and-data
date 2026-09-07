# Authoritative Stage 6 state indexing for Stages 7.1-8.1.
#
# Reads the shared Stage 6 initialization once plus completed removal ledgers
# and patch-lookup filenames; returns one row per completed state. Writes
# nothing. The initialization is released before returning, so multiple curves
# are indexed without retaining a large spatial object. `stage_meta.csv` is
# validated only as a Stage 7.1 export.


stage7_lookup_files <- function(directory) {
  need_dir(directory, "Stage 6 patch lookup directory")
  files <- list.files(
    directory, pattern = "^stage_patch_lookup_stage_[0-9]{4}[.]csv$",
    full.names = TRUE
  )
  assert(length(files) > 0L, paste0("No Stage 6 patch lookups in: ", directory))
  table <- data.table::data.table(path = files)
  table[, stage := as.integer(sub(
    "^stage_patch_lookup_stage_([0-9]{4})[.]csv$", "\\1", basename(path)
  ))]
  data.table::setorder(table, stage)
  assert(!anyDuplicated(table$stage) &&
           identical(table$stage, seq_len(max(table$stage))),
         "Stage 6 patch lookups must be unique and contiguous from stage one.")
  table[]
}

read_stage6_state_index <- function(run, optimization_curve, config,
                                    bundle = NULL, removal_events = NULL) {
  if (is.null(bundle)) {
    need_file(
      config$paths$initialization_bundle,
      "Stage 6 shared initialization"
    )
  }
  if (is.null(removal_events)) {
    need_file(run$removal_events, paste0(optimization_curve, " removal events"))
  }
  lookup_files <- stage7_lookup_files(run$patch_lookup_dir)
  if (is.null(bundle)) bundle <- readRDS(config$paths$initialization_bundle)
  validate_priority_initialization_schema(bundle, "Stage 6 shared initialization")
  validate_priority_initialization_metadata(
    bundle, config$taxa_tag, config$sdm,
    "Stage 6 shared initialization", contract = config$contract
  )
  alive <- as.integer(bundle$alive_species_count_by_cell) > 0L
  initial_cells <- as.integer(sum(alive))
  initial_area <- as.numeric(sum(bundle$cell_area_by_cell[alive]))
  events <- if (is.null(removal_events)) {
    data.table::fread(run$removal_events)
  } else {
    data.table::copy(data.table::as.data.table(removal_events))
  }
  required <- c(
    "removal_step", "stage", "event_type", "pruning_iteration", "cells_removed",
    "area_removed_km2", "cum_cells_removed", "cum_area_removed_km2",
    "cum_prop_cells_removed", "cum_prop_area_removed", "cells_retained",
    "area_retained_km2"
  )
  need_cols(events, required, paste0(optimization_curve, " removal_events.csv"))
  data.table::setorder(events, removal_step)
  assert(sum(events$event_type == "final_retained") == 1L,
         paste0(optimization_curve, " requires one final_retained event."))
  integer_fields <- c(
    "removal_step", "stage", "cells_removed", "cum_cells_removed", "cells_retained"
  )
  valid_integers <- vapply(integer_fields, function(field) {
    values <- suppressWarnings(as.numeric(events[[field]]))
    all(is.finite(values) & values >= 0 & values == floor(values))
  }, logical(1L))
  assert(all(valid_integers),
         paste0(optimization_curve, " removal-event integer fields are invalid."))
  assert(identical(as.integer(events$removal_step), seq_len(nrow(events))) &&
           !anyDuplicated(events$removal_step) && all(events$cells_removed >= 0L) &&
           identical(as.numeric(events$cum_cells_removed),
                     cumsum(as.numeric(events$cells_removed))),
         paste0(optimization_curve, " removal-event sequence is invalid."))
  assert(all(is.finite(events$area_removed_km2) & events$area_removed_km2 >= 0) &&
           all(diff(events$cum_area_removed_km2) >= -1e-6) &&
           all(abs(events$cum_area_removed_km2 -
                     cumsum(events$area_removed_km2)) <= 1e-6 * max(1, initial_area)),
         paste0(optimization_curve, " removal-event areas are invalid."))
  assert(all(events$cum_prop_cells_removed >= 0 &
               events$cum_prop_cells_removed <= 1) &&
           all(events$cum_prop_area_removed >= 0 &
                 events$cum_prop_area_removed <= 1),
         paste0(optimization_curve, " removal-event proportions are invalid."))
  assert(all(events$cells_retained + events$cum_cells_removed == initial_cells),
         paste0(optimization_curve, " cell counts disagree with its bundle."))
  area_tolerance <- 1e-5 * max(1, initial_area)
  assert(all(abs(events$area_retained_km2 + events$cum_area_removed_km2 -
                   initial_area) <= area_tolerance),
         paste0(optimization_curve, " areas disagree with its bundle."))

  ecological <- events[event_type != "final_retained"]
  last_by_stage <- ecological[, .SD[.N], by = stage]
  data.table::setorder(last_by_stage, stage)
  within_stage <- ecological[, .(
    cells_removed_stage = sum(cells_removed),
    area_removed_stage_km2 = sum(area_removed_km2),
    events_in_stage = .N
  ), by = stage]
  positive <- data.table::data.table(
    stage = lookup_files$stage,
    patch_lookup_path = normalizePath(
      lookup_files$path, winslash = "/", mustWork = TRUE
    )
  )
  previous <- findInterval(positive$stage, last_by_stage$stage)
  assert(all(previous > 0L),
         paste0(optimization_curve, " lookup has no preceding removal event."))
  final <- last_by_stage[previous]
  positive[, `:=`(
    last_removal_step = as.integer(final$removal_step),
    alive_end = as.integer(final$cells_retained),
    removed_end = as.integer(final$cum_cells_removed),
    area_retained_km2 = as.numeric(final$area_retained_km2),
    area_removed_km2 = as.numeric(final$cum_area_removed_km2),
    prop_cells_removed_end = as.numeric(final$cum_prop_cells_removed),
    prop_area_removed_end = as.numeric(final$cum_prop_area_removed)
  )]
  positive <- merge(positive, within_stage, by = "stage", all.x = TRUE, sort = FALSE)
  positive[is.na(cells_removed_stage), `:=`(
    cells_removed_stage = 0, area_removed_stage_km2 = 0, events_in_stage = 0L
  )]
  state <- data.table::rbindlist(list(
    data.table::data.table(
      stage = 0L, patch_lookup_path = NA_character_, last_removal_step = 0L,
      alive_end = initial_cells, removed_end = 0L,
      area_retained_km2 = initial_area, area_removed_km2 = 0,
      cells_removed_stage = 0L, area_removed_stage_km2 = 0,
      events_in_stage = 0L, prop_cells_removed_end = 0,
      prop_area_removed_end = 0
    ), positive
  ), use.names = TRUE)
  state[, `:=`(
    pct_cells_removed_end = 100 * prop_cells_removed_end,
    pct_area_removed_end = 100 * prop_area_removed_end,
    stage_order = seq_len(.N) - 1L,
    optimization_curve = optimization_curve, run_id = run$run_id,
    keep_n = as.integer(alive_end),
    target_cells_retained = as.integer(alive_end),
    target_area_retained_km2 = as.numeric(area_retained_km2),
    pct_cells_removed = as.numeric(100 * prop_cells_removed_end)
  )]
  assert(all(diff(state$keep_n) <= 0L) &&
           all(diff(state$area_retained_km2) <= area_tolerance),
         paste0(optimization_curve, " retained states must be non-increasing."))
  state[]
}

read_completed_stage6_states <- function(config) {
  need_file(config$paths$initialization_bundle, "Stage 6 shared initialization")
  bundle <- readRDS(config$paths$initialization_bundle)
  states <- stats::setNames(lapply(config$curves, function(curve) {
    read_stage6_state_index(
      config$paths$runs[[curve]], curve, config, bundle = bundle
    )
  }), config$curves)
  rm(bundle)
  initial <- vapply(states, function(state) state$keep_n[[1L]], integer(1L))
  assert(length(unique(initial)) == 1L && initial[[1L]] > 0L,
         "Selected Stage 6 curves must share one non-empty initial domain.")
  states
}
