# Shared benchmark reconstruction and report operations.
#
# Reconstruction and reporting are deliberately disjoint. The resume operation
# may read rasters, load compiled kernels, and write checkpoints/lookups. The
# evaluation operations are read-only and can regenerate every trajectory,
# statistic, and table in a fresh R session from completed exact lookups.
# Missing-target reconstruction appends method-specific operational diagnostics;
# a complete resume returns before opening the log or spatial runtime.

if (!exists("project_source", mode = "function")) {
  source(file.path("R", "project_loader.R"))
}
project_source(c(
  # Shared contracts and benchmark configuration.
  "R/project_utils.R", "R/project_artifacts.R", "R/project_rmd.R",
  "R/project_runtime_log.R", "R/project_transactions.R",
  "R/project_paths.R", "R/demographic_contract.R",
  "R/analysis_contract.R", "R/application_storage.R", "R/storage_artifacts.R",
  "R/application_context.R",
  "R/report_config.R",
  "R/priority_run_config.R", "R/benchmark_config.R", "R/stage7_contracts.R",
  # Read-only planning, evaluation, statistics, and reconstruction interfaces.
  "R/benchmark_targets.R", "R/benchmark_library.R",
  "R/stage7_persistence.R", "R/benchmark_evaluation.R",
  "R/benchmark_statistics.R", "R/benchmark_tables.R",
  "R/benchmark_reconstruction.R"
))

load_benchmark_spatial_runtime <- function(cache_directory) {
  # Spatial reconstruction is intentionally absent from inspect/report modes.
  project_source(c(
    "R/stage7_rank_inputs.R", "R/priority_runtime_kernels.R",
    "R/priority_graph_store.R", "R/priority_pruning_frontier.R",
    "R/priority_csr_graph.R", "R/priority_pruning_graph.R",
    "R/priority_pruning_scoring.R", "R/priority_ecology_logging.R",
    "R/priority_pruning_iteration.R", "R/priority_fragmentation_patches.R",
    "R/priority_fragmentation_graph.R", "R/priority_fragmentation_stage.R",
    "R/priority_distance_geometry.R", "R/priority_distance_graph.R",
    "R/priority_distance_stage.R", "R/stage7_rank_incremental.R"
  ))
  load_stage6_runtime_kernels(cache_directory)
  invisible(TRUE)
}

# Build only missing exact targets. A complete library returns before terra,
# the rank raster, compiled kernels, or checkpoint code is initialized.
resume_benchmark_library_active <- function(config, plan) {
  validate_benchmark_reconstruction_inputs(config)
  initialize_benchmark_library(config)

  # Validate every already committed required target before expensive spatial
  # initialization. A corrupt exact state must never be skipped by a checkpoint.
  benchmark_library_index(config, plan, validate = TRUE)

  started <- proc.time()[["elapsed"]]
  reference <- load_benchmark_initialization(config)
  assert(requireNamespace("terra", quietly = TRUE),
         "Benchmark reconstruction requires package 'terra'.")
  template <- terra::rast(config$paths$template_raster)
  assert(terra::ncell(template) == length(reference$alive_species_count_by_cell),
         "Stage 6 template and bundle cell counts differ.")
  load_benchmark_spatial_runtime(config$paths$rcpp_cache)
  before <- benchmark_library_index(config, plan, validate = FALSE)$exists
  reconstruction <- reconstruct_benchmark_targets(config, plan, reference, template)
  after <- benchmark_library_index(config, plan, validate = TRUE)$exists
  rm(reference, template)
  invisible(gc(FALSE))
  assert(all(after), "Benchmark reconstruction ended with missing exact targets.")
  list(
    benchmark = config$rank_method, status = "complete",
    created = sum(!before & after), reused = sum(before),
    completed_union_target = reconstruction$completed_union_target,
    counters = reconstruction$counters,
    elapsed_seconds = proc.time()[["elapsed"]] - started
  )
}

resume_benchmark_library <- function(config, plan) {
  if (benchmark_library_complete(config, plan)) {
    finalize_benchmark_library(config, plan)
    return(list(
      benchmark = config$rank_method, status = "complete",
      created = 0L, reused = nrow(plan$union),
      completed_union_target = max(plan$union$union_target_id),
      counters = c(rank_raster_reads = 0, kernel_loads = 0),
      runtime_log_path = config$paths$reconstruction_log
    ))
  }
  result <- with_runtime_log(
    path = config$paths$reconstruction_log,
    stage = "7.2",
    operation = "benchmark_reconstruction",
    mode = config$mode %||% "resume",
    context = list(
      application = config$project_paths$application,
      run = config$project_paths$run_tag,
      benchmark = config$rank_method,
      curves = paste(config$curves, collapse = ","),
      required_targets = as.integer(nrow(plan$union))
    ),
    code = function() resume_benchmark_library_active(config, plan)
  )
  result$runtime_log_path <- config$paths$reconstruction_log
  result
}

# Read each required exact lookup once, evaluate it, and return benchmark
# trajectories. This is read-only and fails on an incomplete/incompatible
# library or an invalid exact patch/PU table.
evaluate_benchmark_library <- function(config, plan, parameters) {
  assert(benchmark_library_complete(config, plan), paste0(
    config$rank_label, " lookup library is incomplete. Run Stage 7.2 or the ",
    "Stage 8.1 lookup-libraries chunk first."
  ))
  trajectories <- evaluate_exact_lookup_library(config, plan, parameters)
  validate_benchmark_trajectories(
    trajectories, plan, nrow(parameters$base), "benchmark",
    paste0(config$rank_label, " benchmark trajectories")
  )
  trajectories
}

# Evaluate authoritative Stage 6 states against the same species parameters.
# Supplying run_states avoids rereading state indexes already used for planning.
evaluate_pipeline_states <- function(config, plan, parameters,
                                     run_states = read_completed_stage6_states(config),
                                     initialization = NULL) {
  trajectories <- evaluate_pipeline_trajectories(
    config, run_states, parameters, initialization
  )
  validate_benchmark_trajectories(
    trajectories, plan, nrow(parameters$base), "pipeline",
    "pipeline trajectories"
  )
  trajectories
}
