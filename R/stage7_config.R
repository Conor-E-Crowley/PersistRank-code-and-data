# Stage 7.1 configuration.
#
# Stage 7.1 deterministically reports one completed Stage 6 run. Its selectors
# identify the run; the caller supplies either the application-centered project
# context or the supported low-level canonical context.


stage7_config <- function(params, paths, contract) {
  assert(is.list(params), "Stage 7.1 parameters must be supplied as a list.")
  assert(is.list(paths), "paths must be a project_paths() result.")
  assert(!is.null(paths$priority_run),
         "Stage 7.1 requires canonical application paths.")
  contract <- validate_analysis_contract(contract, "Stage 7.1 abundance contract")

  public_names <- c(
    "curve", "taxa", "sdm", "cells_to_remove_per_iteration",
    "pruning_iterations_per_stage"
  )
  validate_parameter_names(params, public_names, "Stage 7.1")

  curve <- validate_persistence_curve(params$curve, "params$curve")
  taxon_selection <- validate_priority_public_taxa(params$taxa, "params$taxa")
  sdm <- validate_priority_sdm_tag(params$sdm, "params$sdm")
  cells <- validate_priority_count(
    params$cells_to_remove_per_iteration,
    "params$cells_to_remove_per_iteration"
  )
  iterations <- validate_priority_count(
    params$pruning_iterations_per_stage,
    "params$pruning_iterations_per_stage"
  )
  run <- application_curve_paths(paths, curve)
  artifact_paths <- c(
    run,
    list(
      initialization_bundle = paths$priority_initialization,
      patch_dir = paths$patches,
      removal_order = file.path(run$out_curve, "removal_order.tif"),
      rankmap = file.path(run$out_curve, "rankmap.tif"),
      removal_events = file.path(run$out_curve, "removal_events.csv"),
      pipe_lookup_dir = file.path(run$out_curve, "patch_lookup_tables"),
      stage_meta = file.path(run$ana, "stage_meta.csv"),
      si_figures = c(
        pipeline_persistence = file.path(
          paths$si_figures, "S7_pipeline_persistence.png"
        ),
        curve_uncertainty = file.path(
          paths$si_figures, "S7_curve_uncertainty.png"
        ),
        species_persistence = file.path(
          paths$si_figures, "S7_species_persistence.png"
        )
      ),
      cell_removal_order_figure = file.path(
        paths$figures, "cell_removal_order.png"
      )
    )
  )
  list(
    curve = unname(as.character(curve)),
    taxa = taxon_selection$taxa,
    taxa_tag = taxon_selection$taxa_tag,
    sdm = sdm,
    cells_to_remove_per_iteration = cells,
    pruning_iterations_per_stage = iterations,
    run_id = run$run_id,
    contract = contract,
    paths = artifact_paths
  )
}

# Build Stage 7.1 reporting configuration from one completed application run.
stage71_application_config <- function(params, context) {
  stage_params <- list(
    curve = params$curve,
    taxa = context$taxa,
    sdm = context$sdm,
    cells_to_remove_per_iteration = context$cells_to_remove_per_iteration,
    pruning_iterations_per_stage = context$pruning_iterations_per_stage
  )
  report_paths <- context$paths
  report_paths$figures <- file.path(context$paths$priority_figures, params$curve)
  report_paths$si_figures <- file.path(
    context$paths$priority_figures, params$curve, "SI"
  )
  stage7_config(stage_params, report_paths, context$contract)
}

validate_stage7_active_inputs <- function(config) {
  p <- config$paths
  need_dir(p$out_curve, "Stage 6 run directory")
  need_file(p$initialization_bundle, "Stage 6 shared initialization")
  need_file(p$removal_order, "Stage 6 removal_order.tif")
  need_file(p$rankmap, "Stage 6 rankmap.tif")
  need_file(p$removal_events, "Stage 6 removal_events.csv")
  need_dir(p$pipe_lookup_dir, "Stage 6 patch lookup directory")
  invisible(TRUE)
}

log_stage7_config <- function(config) {
  log_msg(
    "CONFIG | stage=7.1",
    " run_id=", config$run_id,
    " curve=", config$curve,
    " taxa=", config$taxa,
    " sdm=", config$sdm
  )
  invisible(config)
}
