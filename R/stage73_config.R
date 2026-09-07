# Canonical Stage 7.3 detailed ABF/main-curve figure configuration.
#
# Returns paths and presentation settings for one curve and one completed exact
# lookup library. Reads no scientific data, writes nothing, and supports only
# inspection or in-memory report regeneration. Stage 7.3 owns the seven-panel figure and
# focal-species results; cross-benchmark statistics belong to Stage 7.4.


stage73_figure_id <- function(rank_method) {
  method <- validate_benchmark_rank_method(rank_method, "rank_method")
  if (identical(method, "abf")) {
    return("persistence_comparison_assemblage_focal")
  }
  paste0("persistence_comparison_assemblage_focal_", method)
}

stage73_figure_path <- function(project, rank_method) {
  file.path(
    project$comparison_figures,
    paste0("persistence_comparison_", toupper(rank_method), ".png")
  )
}

stage73_artifact_paths <- function(project, curve, rank_method) {
  run <- application_curve_paths(project, curve)
  library <- benchmark_lookup_library_paths(project, rank_method)
  list(
    run_id = run$run_id,
    initialization_bundle = project$priority_initialization,
    removal_events = run$removal_events,
    pipe_lookup_dir = run$patch_lookup_dir,
    lookup_library = library,
    figure = stage73_figure_path(project, rank_method)
  )
}

stage73_config <- function(params, paths, contract) {
  allowed <- c(
    "mode", "curve", "taxa", "sdm", "cells_to_remove_per_iteration",
    "pruning_iterations_per_stage", "rank_method", "persistence_threshold",
    "focal_species", "assemblage_point_style", "focal_point_style",
    "assemblage_central_stat"
  )
  validate_parameter_names(params, allowed, "Stage 7.3")
  taxon <- validate_priority_public_taxa(params$taxa, "params$taxa")
  presentation <- report_presentation_settings(params)
  curve <- unname(as.character(validate_persistence_curve(params$curve, "params$curve")))
  method <- validate_benchmark_rank_method(params$rank_method, "params$rank_method")
  artifacts <- stage73_artifact_paths(paths, curve, method)
  list(
    mode = validate_scalar_choice(params$mode, c("inspect", "report"), "params$mode"),
    curve = curve, taxa = taxon$taxa, taxa_tag = taxon$taxa_tag,
    sdm = validate_priority_sdm_tag(params$sdm, "params$sdm"),
    cells_to_remove_per_iteration = validate_priority_count(
      params$cells_to_remove_per_iteration, "params$cells_to_remove_per_iteration"
    ),
    pruning_iterations_per_stage = validate_priority_count(
      params$pruning_iterations_per_stage, "params$pruning_iterations_per_stage"
    ),
    rank_method = method, rank_label = benchmark_rank_label(method),
    persistence_threshold = presentation$persistence_threshold,
    focal_species = presentation$focal_species,
    assemblage_point_style = presentation$assemblage_point_style,
    focal_point_style = presentation$focal_point_style,
    assemblage_central_stat = presentation$assemblage_central_stat,
    contract = validate_analysis_contract(contract, "Stage 7.3 contract"),
    run_id = artifacts$run_id, paths = artifacts, project_paths = paths
  )
}

validate_stage73_inputs <- function(config) {
  need_file(config$paths$initialization_bundle, "Stage 6 shared initialization")
  need_file(config$paths$removal_events, "Stage 6 removal events")
  need_dir(config$paths$pipe_lookup_dir, "Stage 6 patch lookups")
  need_file(config$paths$lookup_library$manifest, "benchmark lookup manifest")
  need_dir(config$paths$lookup_library$lookups, "benchmark exact lookups")
  invisible(TRUE)
}

inspect_stage73 <- function(config) {
  paths <- c(
    initialization = config$paths$initialization_bundle,
    removal_events = config$paths$removal_events,
    pipeline_lookups = config$paths$pipe_lookup_dir,
    benchmark_manifest = config$paths$lookup_library$manifest,
    benchmark_lookups = config$paths$lookup_library$lookups,
    figure = config$paths$figure
  )
  status <- project_artifact_status(paths)
  status$artifact <- names(paths)
  status[, c("artifact", "path", "exists", "type", "size_bytes", "modified")]
}

stage73_application_config <- function(params, paths = project_paths()) {
  context_params <- params
  context_params$optimization_curves <- params$main_curve %||% "q50"
  context <- benchmark_run_context(context_params, paths)
  stage73_config(list(
    mode = params$mode, curve = params$main_curve %||% "q50",
    taxa = context$taxa, sdm = context$sdm,
    cells_to_remove_per_iteration = params$cells_to_remove_per_iteration,
    pruning_iterations_per_stage = params$pruning_iterations_per_stage,
    rank_method = params$rank_method %||% "abf",
    persistence_threshold = params$persistence_threshold,
    focal_species = params$focal_species,
    assemblage_point_style = params$assemblage_point_style,
    focal_point_style = params$focal_point_style,
    assemblage_central_stat = params$assemblage_central_stat
  ), context$project_paths, context$contract)
}
