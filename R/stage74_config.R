# Stage 7.4 report-only configuration.
#
# Reads and validates the canonical application, scenario, run, and persistence
# manifests. Returns one benchmark-method configuration augmented only with
# report settings. Writes nothing and performs no spatial reconstruction.


stage74_config <- function(params, paths = project_paths()) {
  allowed <- c(
    "mode", "application", "quasi_extinction_abundance",
    "persistence_horizon_years", "cells_to_remove_per_iteration",
    "pruning_iterations_per_stage", "optimization_curves", "main_curve",
    "rank_method", "late_stage_removal_percentages", "persistence_threshold"
  )
  validate_parameter_names(params, allowed, "Stage 7.4")
  context <- benchmark_run_context(params, paths)
  config <- benchmark_method_config(list(
    taxa = context$taxa, sdm = context$sdm,
    cells_to_remove_per_iteration = context$cells_to_remove_per_iteration,
    pruning_iterations_per_stage = context$pruning_iterations_per_stage,
    rank_method = params$rank_method,
    optimization_curves = context$curves
  ), context$project_paths, context$contract)
  report <- benchmark_report_settings(params, context$curves)
  config[names(report)] <- report
  config$mode <- validate_scalar_choice(
    params$mode, c("inspect", "report"), "params$mode"
  )
  config$scenario <- context$scenario
  config
}
