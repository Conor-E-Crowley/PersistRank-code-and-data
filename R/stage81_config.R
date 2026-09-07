# Stage 8.1 canonical multi-benchmark configuration.
#
# Reads only the canonical application/run/model manifests through the shared
# benchmark context. Returns one compact configuration for four fixed methods;
# writes nothing and performs no lookup or spatial work.


stage81_config <- function(params, paths = project_paths()) {
  allowed <- c(
    "mode", "application", "quasi_extinction_abundance",
    "persistence_horizon_years", "cells_to_remove_per_iteration",
    "pruning_iterations_per_stage", "optimization_curves", "main_curve",
    "late_stage_removal_percentages", "persistence_threshold", "focal_species",
    "assemblage_point_style", "focal_point_style", "assemblage_central_stat"
  )
  validate_parameter_names(params, allowed, "Stage 8.1")
  context <- benchmark_run_context(params, paths)
  report <- benchmark_report_settings(params, context$curves)
  presentation <- report_presentation_settings(params)
  list(
    mode = validate_scalar_choice(params$mode, c("inspect", "resume"), "params$mode"),
    application = context$application,
    cells_to_remove_per_iteration = context$cells_to_remove_per_iteration,
    pruning_iterations_per_stage = context$pruning_iterations_per_stage,
    taxa = context$taxa, sdm = context$sdm, curves = context$curves,
    main_curve = report$primary_curve,
    late_stage_removal_percentages = report$late_stage_removal_percentages,
    persistence_threshold = report$persistence_threshold,
    methods = benchmark_rank_methods(),
    method_labels = benchmark_rank_labels(),
    focal_species = presentation$focal_species,
    assemblage_point_style = presentation$assemblage_point_style,
    focal_point_style = presentation$focal_point_style,
    assemblage_central_stat = presentation$assemblage_central_stat,
    runtime_log = file.path(
      context$project_paths$priority_run, "Logs", "stage81_report.log"
    ),
    contract = context$contract, scenario = context$scenario,
    project_paths = context$project_paths
  )
}

stage81_benchmark_config <- function(config, method) {
  benchmark_method_config(list(
    taxa = config$taxa, sdm = config$sdm,
    cells_to_remove_per_iteration = config$cells_to_remove_per_iteration,
    pruning_iterations_per_stage = config$pruning_iterations_per_stage,
    rank_method = method, optimization_curves = config$curves
  ), config$project_paths, config$contract)
}
