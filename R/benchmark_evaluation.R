# Shared benchmark persistence evaluation and compact trajectory contracts.
#
# Patch rows are aggregated to population-unit area once per state, then the
# same PU state is evaluated under all five fitted coefficient sets. Only one
# fixed-denominator species row per curve/function/stage is retained; absent
# species are represented explicitly with zero area, PUs, and persistence.


prepare_benchmark_parameters <- function(reference_bundle, config) {
  reference_curve <- if ("q50" %in% config$curves) "q50" else config$curves[[1L]]
  primary_config <- list(curve = reference_curve, taxa_tag = config$taxa_tag,
                         sdm = config$sdm, contract = config$contract)
  base <- prepare_stage71_bundle_parameters(reference_bundle, primary_config)
  curves <- prepare_stage71_curve_parameters(
    base, reference_bundle$species_params, reference_curve
  )
  list(base = base, curves = curves)
}

aggregate_benchmark_population_units <- function(lookup, base_params, stage, label) {
  patches <- standardize_stage71_patch_lookup(lookup, base_params, stage, label)
  if (!nrow(patches)) {
    return(data.table::data.table(
      scientificName = character(), species = character(), pu_id = integer(),
      pu_area_km2 = numeric()
    ))
  }
  pu <- patches[, .(pu_area_km2 = sum(patch_area_km2)),
                by = .(scientificName, species, pu_id)]
  thresholds <- base_params[, .(scientificName, species, threshold_area_km2)]
  pu <- merge(pu, thresholds, by = c("scientificName", "species"),
              all.x = TRUE, sort = FALSE)
  assert(all(is.finite(pu$pu_area_km2) & pu$pu_area_km2 > pu$threshold_area_km2),
         paste0(label, " contains a represented PU at or below its strict area threshold."))
  pu[, threshold_area_km2 := NULL]
  data.table::setorder(pu, scientificName, pu_id)
  pu[]
}

evaluate_benchmark_population_units <- function(pu, parameter_sets, optimization_curve, method,
                                 stage, stage_order, keep_n,
                                 pct_cells_removed) {
  base <- parameter_sets$base
  pieces <- lapply(persistence_curves(), function(evaluation_curve) {
    parameters <- parameter_sets$curves[curve == evaluation_curve, .(
      scientificName, species, density, threshold_area_km2,
      alpha, beta, parameter_order
    )]
    represented <- if (nrow(pu)) {
      joined <- merge(pu, parameters, by = c("scientificName", "species"),
                      all.x = TRUE, sort = FALSE)
      assert(all(!is.na(joined$alpha)),
             "A benchmark PU could not be joined to its persistence parameters.")
      joined[, pu_persistence := stage7_pu_persistence(
        pu_area_km2, density, alpha, beta, threshold_area_km2,
        paste0("Benchmark ", method, " PU persistence")
      )]
      joined[, .(
        species_persistence = stage7_species_persistence(pu_persistence),
        n_pu = .N,
        total_pu_area_km2 = sum(pu_area_km2)
      ), by = .(scientificName, species)]
    } else {
      data.table::data.table(
        scientificName = character(), species = character(),
        species_persistence = numeric(), n_pu = integer(),
        total_pu_area_km2 = numeric()
      )
    }
    full <- merge(
      base[, .(scientificName, species, taxon, parameter_order)], represented,
      by = c("scientificName", "species"), all.x = TRUE, sort = FALSE
    )
    full[is.na(species_persistence), `:=`(
      species_persistence = 0, n_pu = 0L, total_pu_area_km2 = 0
    )]
    data.table::setorder(full, parameter_order)
    full[, .(
      optimization_curve = as.character(optimization_curve),
      evaluation_curve = evaluation_curve,
      method = as.character(method),
      stage = as.integer(stage),
      stage_order = as.integer(stage_order),
      keep_n = as.integer(keep_n),
      pct_cells_removed = as.numeric(pct_cells_removed),
      scientificName, species, taxon,
      species_persistence = as.numeric(species_persistence),
      n_pu = as.integer(n_pu),
      total_pu_area_km2 = as.numeric(total_pu_area_km2)
    )]
  })
  out <- data.table::rbindlist(pieces, use.names = TRUE)
  assert(
    nrow(out) == nrow(base) * length(persistence_curves()) &&
      !anyDuplicated(out[, .(evaluation_curve, species)]) &&
      all(is.finite(out$species_persistence) & out$species_persistence >= 0 &
            out$species_persistence <= 1),
    paste0("Invalid compact ", method, " persistence state at stage ", stage, ".")
  )
  out[]
}

benchmark_trajectory_columns <- function() c(
  "optimization_curve", "evaluation_curve", "method", "stage", "stage_order",
  "keep_n", "pct_cells_removed", "scientificName", "species", "taxon",
  "species_persistence", "n_pu", "total_pu_area_km2"
)

validate_benchmark_trajectories <- function(x, plan, n_species, method,
                                          label = "benchmark trajectories") {
  dt <- data.table::as.data.table(x)
  required <- benchmark_trajectory_columns()
  assert(identical(names(dt), required), paste0(
    label, " must contain exactly: ", paste(required, collapse = ", "), "."
  ))
  assert(nrow(dt) == nrow(plan$consumers) * length(persistence_curves()) * n_species,
         paste0(label, " has incomplete curve/stage/species coverage."))
  assert(all(dt$method == method) &&
           !anyDuplicated(dt[, .(optimization_curve, evaluation_curve, stage, species)]),
         paste0(label, " contains invalid method labels or duplicate rows."))
  assert(all(dt$optimization_curve %in% persistence_curves()) &&
           all(dt$evaluation_curve %in% persistence_curves()) &&
           all(is.finite(dt$species_persistence) & dt$species_persistence >= 0 &
                 dt$species_persistence <= 1) &&
           all(!is.na(dt$n_pu) & dt$n_pu >= 0L) &&
           all(is.finite(dt$total_pu_area_km2) & dt$total_pu_area_km2 >= 0),
         paste0(label, " contains invalid values."))

  # Validate target identity separately from the repeated species/function
  # rows. This prevents a structurally plausible cache from silently pairing
  # a stage with a different retained-cell target.
  target_fields <- c(
    "optimization_curve", "stage", "stage_order", "keep_n",
    "pct_cells_removed"
  )
  observed_targets <- unique(dt[, target_fields, with = FALSE])
  expected_targets <- plan$consumers[, target_fields, with = FALSE]
  data.table::setorder(observed_targets, optimization_curve, stage)
  data.table::setorder(expected_targets, optimization_curve, stage)
  assert(benchmark_target_plans_equal(observed_targets, expected_targets),
         paste0(label, " does not match the current retained-cell target plan."))
  coverage <- dt[, .N, by = .(optimization_curve, evaluation_curve, stage)]
  species_coverage <- dt[, .N, by = species]
  assert(
    data.table::uniqueN(dt[, .(optimization_curve, evaluation_curve)]) ==
      data.table::uniqueN(dt$optimization_curve) * length(persistence_curves()) &&
      all(coverage$N == n_species) &&
      data.table::uniqueN(dt$species) == n_species &&
      all(species_coverage$N == nrow(plan$consumers) * length(persistence_curves())),
    paste0(label, " does not use a fixed complete species denominator.")
  )
  invisible(dt)
}

expand_benchmark_state <- function(state_rows, consumers) {
  data.table::rbindlist(lapply(seq_len(nrow(consumers)), function(index) {
    consumer <- consumers[index]
    expanded <- data.table::copy(state_rows)
    expanded[, `:=`(
      optimization_curve = consumer$optimization_curve,
      stage = as.integer(consumer$stage),
      stage_order = as.integer(consumer$stage_order),
      keep_n = as.integer(consumer$keep_n),
      pct_cells_removed = as.numeric(consumer$pct_cells_removed)
    )]
    expanded
  }), use.names = TRUE)
}

# Stream one curve's Stage 6 patch sequence. No stage-level patch or PU table
# survives the iteration, which bounds memory independently of stage count.
evaluate_pipeline_curve <- function(curve, bundle, run_meta,
                                          parameter_sets) {
  curve_started <- proc.time()[["elapsed"]]
  lookup_read_seconds <- 0
  pu_aggregation_seconds <- 0
  persistence_seconds <- 0
  stages <- vector("list", nrow(run_meta))
  for (i in seq_len(nrow(run_meta))) {
    row <- run_meta[i]
    stage <- as.integer(row$stage)
    lookup_started <- proc.time()[["elapsed"]]
    lookup <- if (stage == 0L) {
      bundle$patch_table
    } else {
      data.table::fread(row$patch_lookup_path)
    }
    lookup_read_seconds <- lookup_read_seconds +
      proc.time()[["elapsed"]] - lookup_started
    aggregation_started <- proc.time()[["elapsed"]]
    pu <- aggregate_benchmark_population_units(
      lookup, parameter_sets$base, stage,
      paste0(curve, " pipeline stage ", stage)
    )
    pu_aggregation_seconds <- pu_aggregation_seconds +
      proc.time()[["elapsed"]] - aggregation_started
    persistence_started <- proc.time()[["elapsed"]]
    stages[[i]] <- evaluate_benchmark_population_units(
      pu, parameter_sets, curve, "pipeline", stage,
      row$stage_order, row$keep_n, row$pct_cells_removed
    )
    persistence_seconds <- persistence_seconds +
      proc.time()[["elapsed"]] - persistence_started
    rm(lookup, pu)
  }
  bind_started <- proc.time()[["elapsed"]]
  trajectories <- data.table::rbindlist(stages, use.names = TRUE)
  curve_bind_seconds <- proc.time()[["elapsed"]] - bind_started
  total_seconds <- proc.time()[["elapsed"]] - curve_started
  timing <- c(
    lookup_read_seconds = lookup_read_seconds,
    pu_aggregation_seconds = pu_aggregation_seconds,
    persistence_seconds = persistence_seconds,
    curve_bind_seconds = curve_bind_seconds,
    total_seconds = total_seconds
  )
  assert(all(is.finite(timing) & timing >= 0),
         "Benchmark pipeline-curve timings must be finite and non-negative.")
  list(trajectories = trajectories, timing = timing)
}

evaluate_pipeline_trajectories <- function(
  config,
  run_meta,
  parameter_sets,
  initialization = NULL
) {
  initialization_read_started <- proc.time()[["elapsed"]]
  owns_initialization <- is.null(initialization)
  if (owns_initialization) initialization <- load_benchmark_initialization(config)
  initialization_read_seconds <- if (owns_initialization) {
    proc.time()[["elapsed"]] - initialization_read_started
  } else 0
  pieces <- vector("list", length(config$curves))
  for (i in seq_along(config$curves)) {
    curve_started <- proc.time()[["elapsed"]]
    curve <- config$curves[[i]]
    curve_result <- evaluate_pipeline_curve(
      curve, initialization, run_meta[[curve]], parameter_sets
    )
    pieces[[i]] <- curve_result$trajectories
    curve_timing <- curve_result$timing
    trajectory_rows <- nrow(pieces[[i]])
    cleanup_started <- proc.time()[["elapsed"]]
    rm(curve_result)
    gc(FALSE)
    cleanup_gc_seconds <- proc.time()[["elapsed"]] - cleanup_started
    accounted_seconds <- sum(c(
      curve_timing[c(
        "lookup_read_seconds", "pu_aggregation_seconds",
        "persistence_seconds", "curve_bind_seconds"
      )],
      cleanup_gc_seconds
    ))
    total_seconds <- proc.time()[["elapsed"]] - curve_started
    residual_seconds <- max(0, total_seconds - accounted_seconds)
    diagnostics <- c(
      initialization_read_seconds = if (i == 1L) {
        initialization_read_seconds
      } else 0,
      curve_timing[c(
        "lookup_read_seconds", "pu_aggregation_seconds",
        "persistence_seconds", "curve_bind_seconds"
      )],
      cleanup_gc_seconds = cleanup_gc_seconds,
      accounted_seconds = accounted_seconds,
      residual_seconds = residual_seconds,
      total_seconds = total_seconds
    )
    assert(all(is.finite(diagnostics) & diagnostics >= 0),
           "Benchmark pipeline diagnostics must be finite and non-negative.")
    runtime_log_event(
      "benchmark_pipeline_curve",
      curve = curve,
      stages = nrow(run_meta[[curve]]),
      initialization_read_seconds = sprintf(
        "%.3f", if (i == 1L) initialization_read_seconds else 0
      ),
      lookup_read_seconds = sprintf(
        "%.3f", curve_timing[["lookup_read_seconds"]]
      ),
      pu_aggregation_seconds = sprintf(
        "%.3f", curve_timing[["pu_aggregation_seconds"]]
      ),
      persistence_seconds = sprintf(
        "%.3f", curve_timing[["persistence_seconds"]]
      ),
      curve_bind_seconds = sprintf(
        "%.3f", curve_timing[["curve_bind_seconds"]]
      ),
      cleanup_gc_seconds = sprintf("%.3f", cleanup_gc_seconds),
      accounted_seconds = sprintf("%.3f", accounted_seconds),
      residual_seconds = sprintf("%.3f", residual_seconds),
      total_seconds = sprintf("%.3f", total_seconds),
      trajectory_rows = trajectory_rows
    )
  }
  if (owns_initialization) rm(initialization)
  data.table::rbindlist(pieces, use.names = TRUE)
}

# Re-evaluate durable exact-target lookups without reading the rank raster or
# running patch/PU repair. This is the normal application-storage report path.
evaluate_exact_lookup_library <- function(config, plan, parameter_sets) {
  evaluation_started <- proc.time()[["elapsed"]]
  lookup_read_seconds <- 0
  pu_aggregation_seconds <- 0
  persistence_seconds <- 0
  consumer_expansion_seconds <- 0
  cleanup_seconds <- 0
  pieces <- vector("list", nrow(plan$union))
  for (i in seq_len(nrow(plan$union))) {
    target <- plan$union[i]
    retained <- as.integer(target$keep_n)
    phase_started <- proc.time()[["elapsed"]]
    object <- read_benchmark_lookup(config, retained, validate_manifest = FALSE)
    lookup_read_seconds <- lookup_read_seconds +
      proc.time()[["elapsed"]] - phase_started
    consumers <- plan$consumers[union_target_id == target$union_target_id]
    phase_started <- proc.time()[["elapsed"]]
    pu <- aggregate_benchmark_population_units(
      object$patch_table, parameter_sets$base, target$union_target_id,
      paste0("benchmark lookup retained_cells=", retained)
    )
    pu_aggregation_seconds <- pu_aggregation_seconds +
      proc.time()[["elapsed"]] - phase_started
    phase_started <- proc.time()[["elapsed"]]
    one <- evaluate_benchmark_population_units(
      pu, parameter_sets, consumers$optimization_curve[[1L]], "benchmark",
      consumers$stage[[1L]], consumers$stage_order[[1L]], retained,
      consumers$pct_cells_removed[[1L]]
    )
    persistence_seconds <- persistence_seconds +
      proc.time()[["elapsed"]] - phase_started
    phase_started <- proc.time()[["elapsed"]]
    pieces[[i]] <- expand_benchmark_state(one, consumers)
    consumer_expansion_seconds <- consumer_expansion_seconds +
      proc.time()[["elapsed"]] - phase_started
    phase_started <- proc.time()[["elapsed"]]
    rm(object, pu, one)
    cleanup_seconds <- cleanup_seconds +
      proc.time()[["elapsed"]] - phase_started
  }
  bind_started <- proc.time()[["elapsed"]]
  trajectories <- data.table::rbindlist(pieces, use.names = TRUE)
  bind_seconds <- proc.time()[["elapsed"]] - bind_started
  total_seconds <- proc.time()[["elapsed"]] - evaluation_started
  accounted_seconds <- sum(c(
    lookup_read_seconds, pu_aggregation_seconds, persistence_seconds,
    consumer_expansion_seconds, bind_seconds, cleanup_seconds
  ))
  diagnostics <- c(
    lookup_read_seconds = lookup_read_seconds,
    pu_aggregation_seconds = pu_aggregation_seconds,
    persistence_seconds = persistence_seconds,
    consumer_expansion_seconds = consumer_expansion_seconds,
    bind_seconds = bind_seconds,
    cleanup_seconds = cleanup_seconds,
    accounted_seconds = accounted_seconds,
    residual_seconds = max(0, total_seconds - accounted_seconds),
    total_seconds = total_seconds
  )
  assert(all(is.finite(diagnostics) & diagnostics >= 0),
         "Benchmark lookup-evaluation timings must be finite and non-negative.")
  runtime_log_event(
    "benchmark_lookup_evaluation",
    benchmark = config$rank_method,
    exact_targets = as.integer(nrow(plan$union)),
    consumer_states = as.integer(nrow(plan$consumers)),
    lookup_read_seconds = sprintf("%.3f", lookup_read_seconds),
    pu_aggregation_seconds = sprintf("%.3f", pu_aggregation_seconds),
    persistence_seconds = sprintf("%.3f", persistence_seconds),
    consumer_expansion_seconds = sprintf("%.3f", consumer_expansion_seconds),
    bind_seconds = sprintf("%.3f", bind_seconds),
    cleanup_seconds = sprintf("%.3f", cleanup_seconds),
    accounted_seconds = sprintf("%.3f", accounted_seconds),
    residual_seconds = sprintf("%.3f", max(0, total_seconds - accounted_seconds)),
    total_seconds = sprintf("%.3f", total_seconds),
    trajectory_rows = as.integer(nrow(trajectories))
  )
  trajectories
}
