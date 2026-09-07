# Canonical paths for reusable models and spatial applications.
#
# Loaded by model, application, and reporting workflows after shared path and
# analysis contracts. Sourcing defines deterministic schema, identity, and path
# helpers only: it performs no filesystem access, package loading, or writes.


application_storage_schema <- function() 1L

validate_application_name <- function(x, label = "application") {
  value <- tolower(validate_scalar_string(x, label))
  assert(
    grepl("^[a-z][a-z0-9]*(?:_[a-z0-9]+)*$", value),
    paste0(
      label, " must be a short lowercase name using letters, numbers, and ",
      "single underscores (for example: madagascar, africa, or global)."
    )
  )
  value
}

storage_threshold_tag <- function(quasi_extinction_abundance,
                                  persistence_horizon_years = 1000L) {
  threshold <- validate_scalar_integer(
    quasi_extinction_abundance, "quasi_extinction_abundance"
  )
  horizon <- validate_scalar_integer(
    persistence_horizon_years, "persistence_horizon_years"
  )
  tag <- paste0("qe_", threshold)
  if (horizon != 1000L) tag <- paste0(tag, "_horizon_", horizon, "yr")
  tag
}

storage_run_tag <- function(cells_to_remove_per_iteration,
                            pruning_iterations_per_stage) {
  cells <- validate_scalar_integer(
    cells_to_remove_per_iteration, "cells_to_remove_per_iteration"
  )
  iterations <- validate_scalar_integer(
    pruning_iterations_per_stage, "pruning_iterations_per_stage"
  )
  paste0("remove_", cells, "_stage_", iterations)
}

demography_model_paths <- function(base = project_paths()) {
  root <- file.path(base$models, "Demography")
  figures <- file.path(base$model_figures, "Demography")
  list(
    root = root,
    outputs = file.path(root, "Outputs"),
    manifest = file.path(root, "model_manifest.rds"),
    figures = figures,
    project = list(
      root = base$root, data = root, raw = base$raw,
      clean = file.path(root, "Outputs"), results = file.path(root, "Outputs"),
      checkpoints = file.path(root, "Checkpoints"), figures = figures,
      si_figures = file.path(figures, "SI"), rcpp_cache = base$rcpp_cache
    )
  )
}

persistence_model_paths <- function(base = project_paths(),
                                    quasi_extinction_abundance,
                                    persistence_horizon_years = 1000L) {
  tag <- storage_threshold_tag(
    quasi_extinction_abundance, persistence_horizon_years
  )
  root <- file.path(base$persistence_models, tag)
  stage2 <- file.path(root, "Stage2")
  stage3 <- file.path(root, "Stage3")
  demography <- demography_model_paths(base)
  list(
    root = root,
    tag = tag,
    manifest = file.path(root, "model_manifest.rds"),
    run_state = file.path(root, "run_state.rds"),
    stage2 = stage2,
    stage3 = stage3,
    project = list(
      root = base$root,
      data = root,
      raw = base$raw,
      clean = stage3,
      results = stage2,
      checkpoints = file.path(stage2, "Checkpoints"),
      figures = file.path(base$model_figures, "Persistence", tag),
      si_figures = file.path(base$model_figures, "Persistence", tag, "SI"),
      rcpp_cache = base$rcpp_cache
    ),
    source_project = list(
      root = base$root,
      data = demography$root,
      raw = base$raw,
      clean = demography$outputs,
      results = demography$outputs,
      checkpoints = file.path(demography$root, "Checkpoints"),
      figures = demography$figures,
      si_figures = file.path(demography$figures, "SI"),
      rcpp_cache = base$rcpp_cache
    )
  )
}

application_root_paths <- function(base = project_paths(), application) {
  application <- validate_application_name(application)
  root <- file.path(base$applications, application)
  list(
    application = application,
    root = root,
    manifest = file.path(root, "application_manifest.rds"),
    inputs = file.path(root, "application_inputs.rds")
  )
}

application_scenario_paths <- function(base = project_paths(), application,
                                       quasi_extinction_abundance,
                                       persistence_horizon_years = 1000L,
                                       cells_to_remove_per_iteration = 1000L,
                                       pruning_iterations_per_stage = 50L) {
  app <- application_root_paths(base, application)
  threshold_tag <- storage_threshold_tag(
    quasi_extinction_abundance, persistence_horizon_years
  )
  run_tag <- storage_run_tag(
    cells_to_remove_per_iteration, pruning_iterations_per_stage
  )
  root <- file.path(app$root, threshold_tag)
  species <- file.path(root, "Species")
  spatial <- file.path(root, "Spatial")
  priority_initialization <- file.path(spatial, "priority_initialization.rds")
  zonation <- file.path(root, "Zonation")
  runs <- file.path(root, "Runs")
  run_root <- file.path(runs, run_tag)
  figure_root <- file.path(
    base$application_figures, app$application, threshold_tag, run_tag
  )
  list(
    application = app$application,
    application_root = app$root,
    application_manifest = app$manifest,
    application_inputs = app$inputs,
    threshold_tag = threshold_tag,
    root = root,
    manifest = file.path(root, "scenario_manifest.rds"),
    species = species,
    spatial = spatial,
    priority_initialization = priority_initialization,
    zonation = zonation,
    runs = runs,
    run_tag = run_tag,
    run_root = run_root,
    run_manifest = file.path(run_root, "run_manifest.rds"),
    run_state = file.path(run_root, "run_state.rds"),
    benchmark_lookups = file.path(run_root, "BenchmarkLookups"),
    figure_root = figure_root,
    project = list(
      root = base$root,
      data = root,
      raw = base$raw,
      clean = species,
      results = root,
      checkpoints = file.path(root, "Checkpoints"),
      figures = file.path(figure_root, "Priority"),
      si_figures = file.path(figure_root, "Comparisons"),
      rcpp_cache = base$rcpp_cache,
      patches = file.path(spatial, "Patches"),
      patch_lookup = file.path(spatial, "all_patch_lookup.rds"),
      connectivity = file.path(spatial, "all_connectivity.rds"),
      stage5_metadata = file.path(spatial, "stage5_build_metadata.rds"),
      priority_initialization = priority_initialization,
      binary_patches = file.path(zonation, "Patches_binary"),
      zonation = zonation,
      zonation_feature_list = file.path(zonation, "feature_list.txt"),
      stage5_checkpoints = file.path(spatial, "Checkpoints"),
      priority_runs = file.path(run_root, "Curves"),
      priority_run = run_root,
      benchmark_lookups = file.path(run_root, "BenchmarkLookups"),
      application = app$application,
      threshold_tag = threshold_tag,
      run_tag = run_tag,
      scenario_root = root,
      application_manifest = app$manifest,
      application_inputs = app$inputs,
      run_manifest = file.path(run_root, "run_manifest.rds"),
      run_state = file.path(run_root, "run_state.rds"),
      spatial_figures = file.path(figure_root, "Spatial"),
      priority_figures = file.path(figure_root, "Priority"),
      comparison_figures = file.path(figure_root, "Comparisons")
    )
  )
}

application_curve_paths <- function(project, curve) {
  curve <- unname(as.character(validate_persistence_curve(curve, "curve")))
  assert(is.list(project) && !is.null(project$priority_run),
         "application_curve_paths() requires canonical application paths.")
  root <- file.path(project$priority_run, "Curves", curve)
  run <- file.path(root, "Run")
  list(
    run_id = curve,
    out_curve = run,
    ana = file.path(run, "Analysis"),
    removal_events = file.path(run, "removal_events.csv"),
    removal_order = file.path(run, "removal_order.tif"),
    patch_lookup_dir = file.path(run, "patch_lookup_tables")
  )
}

benchmark_lookup_library_paths <- function(project, rank_method) {
  method <- validate_benchmark_rank_method(rank_method, "rank_method")
  assert(!is.null(project$benchmark_lookups),
         "Benchmark lookup paths require an application priority run.")
  root <- file.path(project$benchmark_lookups, toupper(method))
  list(
    rank_method = method,
    rank_label = benchmark_rank_label(method),
    root = root,
    manifest = file.path(root, "lookup_manifest.rds"),
    checkpoint = file.path(root, "checkpoint.rds"),
    checkpoint_backup = file.path(root, "checkpoint.rds.bak"),
    lookups = file.path(root, "Lookups")
  )
}

benchmark_lookup_file <- function(library, retained_cells) {
  assert(
    is.numeric(retained_cells) && length(retained_cells) > 0L &&
      all(is.finite(retained_cells)) &&
      all(retained_cells == floor(retained_cells)) &&
      all(retained_cells >= 0) && all(retained_cells <= .Machine$integer.max),
    "retained_cells must contain nonnegative integers within the R integer range."
  )
  retained_cells <- as.integer(retained_cells)
  file.path(library$lookups, paste0("retained_cells_", retained_cells, ".rds"))
}
