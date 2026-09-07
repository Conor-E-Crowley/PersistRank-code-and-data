# Shared benchmark manuscript statistics.
#
# Pipeline and benchmark states are paired at exact retained-cell targets.
# Complete-sequence differences use the pipeline run's removal-percentage axis;
# retained-area differences are integrated exactly over the combined knots of
# both piecewise-linear PU-area curves. Direction classifications use the fixed
# 1e-9 near-zero tolerance.

benchmark_tolerance <- function() 1e-9

benchmark_trapezoid <- function(x, y) {
  dt <- data.table::data.table(x = as.numeric(x), y = as.numeric(y))
  dt <- dt[is.finite(x) & is.finite(y), .(y = mean(y)), by = x]
  data.table::setorder(dt, x)
  if (nrow(dt) < 2L) return(NA_real_)
  sum(diff(dt$x) * (head(dt$y, -1L) + tail(dt$y, -1L)) / 2)
}

benchmark_absolute_trapezoid <- function(x, y) {
  dt <- data.table::data.table(x = as.numeric(x), y = as.numeric(y))
  dt <- dt[is.finite(x) & is.finite(y), .(y = mean(y)), by = x]
  data.table::setorder(dt, x)
  if (nrow(dt) < 2L) return(NA_real_)

  x1 <- head(dt$x, -1L)
  x2 <- tail(dt$x, -1L)
  y1 <- head(dt$y, -1L)
  y2 <- tail(dt$y, -1L)
  widths <- x2 - x1
  areas <- widths * (abs(y1) + abs(y2)) / 2
  crossing <- (y1 < 0 & y2 > 0) | (y1 > 0 & y2 < 0)
  if (any(crossing)) {
    zero_x <- x1[crossing] -
      y1[crossing] * widths[crossing] / (y2[crossing] - y1[crossing])
    areas[crossing] <-
      (zero_x - x1[crossing]) * abs(y1[crossing]) / 2 +
      (x2[crossing] - zero_x) * abs(y2[crossing]) / 2
  }
  sum(areas)
}

benchmark_exact_area_comparison <- function(
    pipeline_area,
    pipeline_persistence,
    benchmark_area,
    benchmark_persistence
) {
  invalid_result <- function() {
    list(
      valid = FALSE,
      area_knots = numeric(),
      pipeline_persistence = numeric(),
      benchmark_persistence = numeric(),
      difference = numeric(),
      lower = NA_real_,
      upper = NA_real_,
      span = NA_real_,
      signed_auc = NA_real_,
      mean_difference = NA_real_
    )
  }

  prepare <- function(area, persistence) {
    x <- data.table::data.table(
      area = c(as.numeric(area), 0),
      persistence = c(as.numeric(persistence), 0)
    )
    x <- x[
      is.finite(area) & is.finite(persistence) & area >= 0,
      .(persistence = mean(persistence)),
      by = area
    ]
    data.table::setorder(x, area)
    x
  }

  pipeline <- prepare(pipeline_area, pipeline_persistence)
  benchmark <- prepare(benchmark_area, benchmark_persistence)
  if (nrow(pipeline) < 2L || nrow(benchmark) < 2L) {
    return(invalid_result())
  }

  lower <- max(min(pipeline$area), min(benchmark$area))
  upper <- min(max(pipeline$area), max(benchmark$area))
  if (!is.finite(lower) || !is.finite(upper) || upper <= lower) {
    return(invalid_result())
  }

  area_knots <- sort(unique(c(
    lower,
    pipeline$area[pipeline$area >= lower & pipeline$area <= upper],
    benchmark$area[benchmark$area >= lower & benchmark$area <= upper],
    upper
  )))
  if (length(area_knots) < 2L) {
    return(invalid_result())
  }

  pipeline_at_knots <- stats::approx(
    pipeline$area, pipeline$persistence, area_knots,
    ties = "ordered", rule = 1
  )$y
  benchmark_at_knots <- stats::approx(
    benchmark$area, benchmark$persistence, area_knots,
    ties = "ordered", rule = 1
  )$y
  difference <- pipeline_at_knots - benchmark_at_knots
  if (any(!is.finite(difference))) {
    return(invalid_result())
  }

  span <- upper - lower
  signed_auc <- benchmark_trapezoid(area_knots, difference)
  list(
    valid = is.finite(signed_auc),
    area_knots = area_knots,
    pipeline_persistence = pipeline_at_knots,
    benchmark_persistence = benchmark_at_knots,
    difference = difference,
    lower = lower,
    upper = upper,
    span = span,
    signed_auc = signed_auc,
    mean_difference = signed_auc / span
  )
}

benchmark_area_difference <- function(pipeline_area, pipeline_persistence,
                                    benchmark_area, benchmark_persistence) {
  benchmark_exact_area_comparison(
    pipeline_area, pipeline_persistence,
    benchmark_area, benchmark_persistence
  )$mean_difference
}

pair_benchmark_trajectories <- function(pipeline, benchmark) {
  keys <- c(
    "optimization_curve", "evaluation_curve", "stage", "stage_order", "keep_n",
    "pct_cells_removed", "scientificName", "species", "taxon"
  )
  p <- data.table::copy(data.table::as.data.table(pipeline))
  b <- data.table::copy(data.table::as.data.table(benchmark))
  p[, method := NULL]
  b[, method := NULL]
  data.table::setnames(p, c("species_persistence", "n_pu", "total_pu_area_km2"),
                       paste0("pipeline_", c("persistence", "n_pu", "area_km2")))
  data.table::setnames(b, c("species_persistence", "n_pu", "total_pu_area_km2"),
                       paste0("benchmark_", c("persistence", "n_pu", "area_km2")))
  paired <- merge(p, b, by = keys, all = FALSE, sort = FALSE)
  assert(nrow(paired) == nrow(pipeline) && nrow(paired) == nrow(benchmark) &&
           !anyDuplicated(paired[, c(keys), with = FALSE]),
         "Pipeline and benchmark trajectories do not pair exactly.")
  paired[, difference := pipeline_persistence - benchmark_persistence]
  paired[]
}

benchmark_late_stage_keys <- function(paired, targets, curves) {
  stages <- unique(paired[, .(
    optimization_curve, stage, stage_order, pct_cells_removed
  )])
  data.table::rbindlist(lapply(curves, function(curve) {
    available <- stages[optimization_curve == curve]
    data.table::rbindlist(lapply(targets, function(target) {
      out <- available[which.min(abs(pct_cells_removed - target))]
      out[, target := as.numeric(target)]
      out
    }))
  }))
}

calculate_benchmark_statistics <- function(pipeline, benchmark, config) {
  paired <- pair_benchmark_trajectories(pipeline, benchmark)
  n_species <- data.table::uniqueN(paired$species)
  eps <- benchmark_tolerance()
  late_keys <- benchmark_late_stage_keys(
    paired, config$late_stage_removal_percentages, config$curves
  )
  late_species <- merge(
    paired, late_keys[, .(optimization_curve, stage, target)],
    by = c("optimization_curve", "stage"), all = FALSE, sort = FALSE
  )
  late <- late_species[, .(
    actual_removal_percentage = unique(pct_cells_removed),
    n_species = .N,
    mean_pipeline_persistence = mean(pipeline_persistence),
    mean_benchmark_persistence = mean(benchmark_persistence),
    mean_difference = mean(difference),
    n_below_pipeline = sum(pipeline_persistence < config$persistence_threshold),
    n_below_benchmark = sum(benchmark_persistence < config$persistence_threshold),
    n_below_only_benchmark = sum(
      benchmark_persistence < config$persistence_threshold &
        pipeline_persistence >= config$persistence_threshold
    ),
    n_below_only_pipeline = sum(
      pipeline_persistence < config$persistence_threshold &
        benchmark_persistence >= config$persistence_threshold
    )
  ), by = .(optimization_curve, evaluation_curve, target, stage)]
  assert(all(late$n_species == n_species),
         "Benchmark late-stage statistics do not use the fixed species denominator.")

  sequence_species <- paired[order(stage_order), {
    span <- max(pct_cells_removed) - min(pct_cells_removed)
    list(
      n_stages = .N,
      mean_difference = benchmark_trapezoid(pct_cells_removed, difference) / span
    )
  }, by = .(optimization_curve, evaluation_curve, scientificName, species)]
  assert(all(is.finite(sequence_species$mean_difference)),
         "Benchmark sequence differences are invalid.")
  sequence <- sequence_species[, .(
    n_species = .N,
    mean_difference = mean(mean_difference),
    n_pipeline = sum(mean_difference > eps),
    n_benchmark = sum(mean_difference < -eps),
    n_near_zero = sum(abs(mean_difference) <= eps)
  ), by = .(optimization_curve, evaluation_curve)]

  area_species <- paired[, .(
    mean_difference = benchmark_area_difference(
      pipeline_area_km2, pipeline_persistence,
      benchmark_area_km2, benchmark_persistence
    )
  ), by = .(optimization_curve, evaluation_curve, scientificName, species)]
  assert(all(is.finite(area_species$mean_difference)),
         "Benchmark retained-area differences are invalid.")
  area <- area_species[, .(
    n_species = .N,
    mean_difference = mean(mean_difference),
    n_pipeline = sum(mean_difference > eps),
    n_benchmark = sum(mean_difference < -eps),
    n_near_zero = sum(abs(mean_difference) <= eps)
  ), by = .(optimization_curve, evaluation_curve)]
  assert(all(sequence$n_pipeline + sequence$n_benchmark + sequence$n_near_zero == n_species) &&
           all(area$n_pipeline + area$n_benchmark + area$n_near_zero == n_species),
         "Benchmark direction categories do not partition the species set.")

  species_stability <- sequence_species[, .(
    positive = sum(mean_difference > eps),
    negative = sum(mean_difference < -eps),
    near_zero = sum(abs(mean_difference) <= eps)
  ), by = .(optimization_curve, scientificName, species)]
  stability <- species_stability[, .(
    n_species = .N,
    pipeline_all_five = sum(positive == 5L),
    benchmark_all_five = sum(negative == 5L),
    near_zero_all_five = sum(near_zero == 5L),
    direction_changes = sum(positive > 0L & negative > 0L)
  ), by = optimization_curve]
  late_stability <- late[, .(
    positive_evaluation_curves = sum(mean_difference > eps),
    benchmark_evaluation_curves = sum(mean_difference < -eps),
    near_zero_evaluation_curves = sum(abs(mean_difference) <= eps)
  ), by = .(optimization_curve, target)]

  list(
    late_keys = late_keys,
    late = late,
    sequence_species = sequence_species,
    sequence = sequence,
    area_species = area_species,
    area = area,
    species_stability = species_stability,
    stability = stability,
    late_stability = late_stability,
    primary_late = late[optimization_curve == config$primary_curve &
                          evaluation_curve == config$primary_curve],
    primary_sequence = sequence[optimization_curve == config$primary_curve &
                                  evaluation_curve == config$primary_curve],
    primary_area = area[optimization_curve == config$primary_curve &
                          evaluation_curve == config$primary_curve],
    primary_robustness = late[optimization_curve == config$primary_curve]
  )
}
