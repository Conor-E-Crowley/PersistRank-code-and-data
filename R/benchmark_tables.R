# Shared format-neutral benchmark table preparation.
#
# Reporting is deliberately separate from scientific calculations. Tables use
# canonical ascending curve order and contain the comparisons emphasized in
# the manuscript. No HTML, CSS, render state, or derived cache is produced.

benchmark_curve_label <- function(x) {
  labels <- benchmark_curve_labels()
  out <- unname(labels[as.character(x)])
  assert(!anyNA(out), "Unsupported curve in benchmark report.")
  out
}

benchmark_result_matrix <- function(x, value, digits = 3L, formatter = NULL) {
  dt <- data.table::copy(data.table::as.data.table(x))
  if (is.null(formatter)) {
    dt[, display_value := round(get(value), digits)]
  } else {
    dt[, display_value := formatter(.SD), .SDcols = names(dt)]
  }
  wide <- data.table::dcast(
    dt, optimization_curve ~ evaluation_curve, value.var = "display_value"
  )
  wide[, optimization_curve := benchmark_curve_label(optimization_curve)]
  selected <- persistence_curves()[persistence_curves() %in% unique(dt$optimization_curve)]
  expected_rows <- benchmark_curve_labels()[selected]
  expected <- benchmark_curve_labels()[persistence_curves()]
  wide <- wide[match(expected_rows, optimization_curve)]
  data.table::setcolorder(wide, c("optimization_curve", persistence_curves()))
  data.table::setnames(wide, c("Optimization curve", unname(expected)))
  wide[]
}

benchmark_direction_cell <- function(dt) {
  sprintf("%d / %d / %d", dt$n_pipeline, dt$n_benchmark, dt$n_near_zero)
}

build_benchmark_report_tables <- function(statistics, config) {
  primary_late <- data.table::copy(statistics$primary_late)[, .(
    `Target removal` = paste0(format(target, trim = TRUE), "%"),
    `Completed stage` = stage,
    `Actual removal` = sprintf("%.1f%%", actual_removal_percentage),
    `Mean pipeline` = round(mean_pipeline_persistence, 3),
    `Mean benchmark` = round(mean_benchmark_persistence, 3),
    `Mean difference` = round(mean_difference, 3),
    `Below pipeline` = n_below_pipeline,
    `Below benchmark` = n_below_benchmark,
    `Below only benchmark` = n_below_only_benchmark,
    `Below only pipeline` = n_below_only_pipeline
  )]
  primary_sequence <- statistics$primary_sequence[, .(
    `Mean sequence difference` = round(mean_difference, 4),
    `Pipeline / benchmark / near zero` = sprintf(
      "%d / %d / %d", n_pipeline, n_benchmark, n_near_zero
    )
  )]
  primary_area <- statistics$primary_area[, .(
    `Mean retained-area difference` = round(mean_difference, 4),
    `Pipeline / benchmark / near zero` = sprintf(
      "%d / %d / %d", n_pipeline, n_benchmark, n_near_zero
    )
  )]
  robustness <- data.table::copy(statistics$primary_robustness)
  robustness[, `:=`(
    `Evaluation curve` = benchmark_curve_label(evaluation_curve),
    `Target removal` = paste0(format(target, trim = TRUE), "%"),
    `Mean difference` = round(mean_difference, 3)
  )]
  robustness <- robustness[, .(
    `Evaluation curve`, `Target removal`, `Mean difference`,
    `Below only benchmark` = n_below_only_benchmark,
    `Below only pipeline` = n_below_only_pipeline
  )]

  late_matrices <- stats::setNames(lapply(config$late_stage_removal_percentages, function(target_i) {
    benchmark_result_matrix(statistics$late[target == target_i], "mean_difference", 3L)
  }), paste0("late_", config$late_stage_removal_percentages))

  stability <- data.table::copy(statistics$stability)
  stability[, `Optimization curve` := benchmark_curve_label(optimization_curve)]
  stability <- stability[, .(
    `Optimization curve`,
    `Species favoring pipeline under all five` = pipeline_all_five,
    `Species favoring benchmark under all five` = benchmark_all_five,
    `Near zero under all five` = near_zero_all_five,
    `Direction changes among functions` = direction_changes
  )]
  list(
    primary_late = primary_late,
    primary_sequence = primary_sequence,
    primary_area = primary_area,
    robustness = robustness,
    late_matrices = late_matrices,
    sequence_mean = benchmark_result_matrix(statistics$sequence, "mean_difference", 4L),
    sequence_directions = benchmark_result_matrix(
      statistics$sequence, "mean_difference", formatter = benchmark_direction_cell
    ),
    area_mean = benchmark_result_matrix(statistics$area, "mean_difference", 4L),
    area_directions = benchmark_result_matrix(
      statistics$area, "mean_difference", formatter = benchmark_direction_cell
    ),
    stability = stability
  )
}
