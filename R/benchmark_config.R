# Canonical configuration for shared benchmark methods and reports.
#
# Inputs: validated application/run settings plus the manifest-based project
# paths. Returns: one small configuration containing scientific settings and
# canonical paths. Reads: application and run manifests only in the public
# application wrapper. Writes nothing. Cost is negligible. Stages 7.2, 7.4,
# and 8.1 share this contract.


benchmark_curve_labels <- function() {
  c(q025 = "2.5%", q16 = "16%", q50 = "50%", q84 = "84%", q975 = "97.5%")
}
BENCHMARK_CHECKPOINT_INTERVAL_TARGETS <- 5L

benchmark_selected_curves <- function(value, completed = persistence_curves()) {
  if (is.null(value) || !length(value)) value <- completed
  if (is.list(value) && identical(names(value), "value")) value <- value$value
  value <- as.character(unlist(value, use.names = FALSE))
  assert(length(value) > 0L && !anyDuplicated(value) && all(value %in% persistence_curves()),
         "optimization_curves must be a unique non-empty canonical subset.")
  persistence_curves()[persistence_curves() %in% value]
}

benchmark_artifact_paths <- function(project, rank_method, curves) {
  assert(is.list(project) && !is.null(project$priority_run),
         "Benchmark workflows require canonical application paths.")
  runs <- stats::setNames(lapply(curves, function(curve) {
    application_curve_paths(project, curve)
  }), curves)
  reference <- if ("q50" %in% curves) "q50" else curves[[1L]]
  list(
    runs = runs,
    initialization_bundle = project$priority_initialization,
    template_raster = runs[[reference]]$removal_order,
    zonation_rankmap = zonation_rankmap_path(project$zonation, rank_method),
    lookup_library = benchmark_lookup_library_paths(project, rank_method),
    checkpoint = benchmark_lookup_library_paths(project, rank_method)$checkpoint,
    checkpoint_dir = benchmark_lookup_library_paths(project, rank_method)$root,
    reconstruction_log = file.path(
      project$priority_run, "Logs",
      paste0("benchmark_", rank_method, "_reconstruction.log")
    ),
    report_log = file.path(
      project$priority_run, "Logs",
      paste0("benchmark_", rank_method, "_report.log")
    ),
    rcpp_cache = project$rcpp_cache
  )
}

# Validate one method configuration. Failure is immediate for unsupported
# settings, missing canonical identity, or a non-canonical curve selection.
benchmark_method_config <- function(params, paths, contract) {
  allowed <- c(
    "taxa", "sdm", "cells_to_remove_per_iteration",
    "pruning_iterations_per_stage", "rank_method", "optimization_curves"
  )
  validate_parameter_names(params, allowed, "Benchmark configuration")
  taxa <- validate_priority_public_taxa(params$taxa, "params$taxa")
  sdm <- validate_priority_sdm_tag(params$sdm, "params$sdm")
  cells <- validate_priority_count(params$cells_to_remove_per_iteration,
                                   "params$cells_to_remove_per_iteration")
  iterations <- validate_priority_count(params$pruning_iterations_per_stage,
                                        "params$pruning_iterations_per_stage")
  method <- validate_benchmark_rank_method(params$rank_method, "params$rank_method")
  curves <- benchmark_selected_curves(params$optimization_curves)
  contract <- validate_analysis_contract(contract, "benchmark persistence contract")
  list(
    taxa = taxa$taxa, taxa_tag = taxa$taxa_tag, sdm = sdm,
    cells_to_remove_per_iteration = cells,
    pruning_iterations_per_stage = iterations,
    rank_method = method, rank_label = benchmark_rank_label(method), curves = curves,
    checkpoint_interval_targets = BENCHMARK_CHECKPOINT_INTERVAL_TARGETS,
    contract = contract,
    paths = benchmark_artifact_paths(paths, method, curves),
    project_paths = paths
  )
}

# Add benchmark curve selection to the application-owned completed-run context.
# Manifest reading and validation remain entirely in application_context.R.
benchmark_run_context <- function(params, paths = project_paths()) {
  context <- read_application_run_context(application_run_identity(params, paths))
  curves <- benchmark_selected_curves(
    params$optimization_curves, context$completed_curves
  )
  assert(all(curves %in% context$completed_curves),
         "Every selected optimization curve must be complete in the run manifest.")
  context$curves <- curves
  context$horizon <- context$persistence_horizon_years
  context$threshold <- context$quasi_extinction_abundance
  context
}

benchmark_report_settings <- function(params, curves) {
  primary <- unname(as.character(validate_persistence_curve(
    params$main_curve %||% "q50", "params$main_curve"
  )))
  assert(primary %in% curves, "main_curve must be selected in optimization_curves.")
  late <- as.numeric(unlist(params$late_stage_removal_percentages, use.names = FALSE))
  assert(length(late) > 0L && all(is.finite(late) & late > 0 & late < 100) &&
           !anyDuplicated(late) && !is.unsorted(late),
         "late_stage_removal_percentages must be unique increasing values in (0,100).")
  persistence_threshold <- validate_report_persistence_threshold(
    params$persistence_threshold
  )
  list(
    curves = curves, primary_curve = primary,
    late_stage_removal_percentages = late,
    persistence_threshold = persistence_threshold
  )
}

validate_benchmark_reconstruction_inputs <- function(config) {
  need_file(config$paths$zonation_rankmap, paste0(config$rank_label, " rank map"))
  need_file(config$project_paths$stage5_metadata, "Stage 5 spatial source metadata")
  need_file(config$paths$initialization_bundle, "Stage 6 shared initialization")
  need_file(config$paths$template_raster, "Stage 6 removal-order template")
  for (curve in config$curves) {
    run <- config$paths$runs[[curve]]
    need_file(run$removal_events, paste0(curve, " Stage 6 removal events"))
    need_dir(run$patch_lookup_dir, paste0(curve, " Stage 6 patch lookups"))
  }
  invisible(TRUE)
}
