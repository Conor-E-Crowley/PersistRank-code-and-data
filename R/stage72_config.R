# Stage 7.2 single-method exact-lookup configuration.
#
# This wrapper adds only the inspect/resume/restart operation to the shared
# reconstruction configuration. Stage 7.2 is the only workflow allowed to
# request restart, which deletes only the selected method library.


stage72_lookup_config <- function(params, paths = project_paths()) {
  allowed <- c(
    "mode", "application", "quasi_extinction_abundance",
    "persistence_horizon_years", "cells_to_remove_per_iteration",
    "pruning_iterations_per_stage", "optimization_curves", "rank_method"
  )
  validate_parameter_names(params, allowed, "Stage 7.2")
  mode <- validate_scalar_choice(
    params$mode, c("inspect", "resume", "restart"), "params$mode"
  )
  context <- benchmark_run_context(params, paths)
  config <- benchmark_method_config(list(
    taxa = context$taxa, sdm = context$sdm,
    cells_to_remove_per_iteration = context$cells_to_remove_per_iteration,
    pruning_iterations_per_stage = context$pruning_iterations_per_stage,
    rank_method = params$rank_method,
    optimization_curves = context$curves
  ), context$project_paths, context$contract)
  config$mode <- mode
  config$scenario <- context$scenario
  config
}
