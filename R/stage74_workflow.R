# Stage 7.4 report-only entry point.
#
# Sources the shared benchmark implementation deterministically. Stage 7.4
# reads completed Stage 6 states and exact benchmark lookups, returns in-memory
# statistics and plain tables, and writes no scientific or recovery artifact.
# Its sole disk side effect is the method-specific append-only runtime log.

if (!exists("project_source", mode = "function")) {
  source(file.path("R", "project_loader.R"))
}
project_source(c(
  # Shared benchmark report definitions; no spatial runtime is loaded here.
  "R/project_utils.R", "R/project_artifacts.R", "R/project_rmd.R",
  "R/project_runtime_log.R", "R/project_transactions.R",
  "R/project_paths.R", "R/project_config.R", "R/demographic_contract.R",
  "R/analysis_contract.R", "R/application_storage.R", "R/storage_artifacts.R",
  "R/application_context.R",
  "R/report_config.R",
  "R/priority_run_config.R", "R/benchmark_config.R", "R/stage7_contracts.R",
  "R/benchmark_targets.R", "R/benchmark_library.R",
  "R/stage7_persistence.R", "R/benchmark_evaluation.R",
  "R/benchmark_statistics.R", "R/benchmark_tables.R",
  "R/benchmark_reconstruction.R", "R/benchmark_workflow.R",
  "R/stage74_config.R"
))

# Regenerate one complete report in memory. Only the operational runtime log is
# written; large pipeline and benchmark trajectories are released before return.
run_stage74_report <- function(config, run_states, plan) {
  started <- proc.time()[["elapsed"]]
  phase_started <- proc.time()[["elapsed"]]
  reference <- load_benchmark_initialization(config)
  parameters <- prepare_benchmark_parameters(reference, config)
  parameter_preparation_seconds <- proc.time()[["elapsed"]] - phase_started

  phase_started <- proc.time()[["elapsed"]]
  pipeline <- evaluate_pipeline_states(
    config, plan, parameters, run_states, reference
  )
  rm(reference)
  invisible(gc(FALSE))
  pipeline_evaluation_seconds <- proc.time()[["elapsed"]] - phase_started

  phase_started <- proc.time()[["elapsed"]]
  benchmark <- evaluate_benchmark_library(config, plan, parameters)
  benchmark_evaluation_seconds <- proc.time()[["elapsed"]] - phase_started

  phase_started <- proc.time()[["elapsed"]]
  statistics <- calculate_benchmark_statistics(pipeline, benchmark, config)
  statistics_seconds <- proc.time()[["elapsed"]] - phase_started

  phase_started <- proc.time()[["elapsed"]]
  tables <- build_benchmark_report_tables(statistics, config)
  table_preparation_seconds <- proc.time()[["elapsed"]] - phase_started

  phase_started <- proc.time()[["elapsed"]]
  rm(pipeline, benchmark, parameters, run_states)
  invisible(gc(FALSE))
  cleanup_seconds <- proc.time()[["elapsed"]] - phase_started
  total_seconds <- proc.time()[["elapsed"]] - started
  accounted_seconds <- sum(c(
    parameter_preparation_seconds, pipeline_evaluation_seconds,
    benchmark_evaluation_seconds, statistics_seconds,
    table_preparation_seconds, cleanup_seconds
  ))
  runtime_log_event(
    "benchmark_report_timing",
    benchmark = config$rank_method,
    parameter_preparation_seconds = sprintf("%.3f", parameter_preparation_seconds),
    pipeline_evaluation_seconds = sprintf("%.3f", pipeline_evaluation_seconds),
    benchmark_evaluation_seconds = sprintf("%.3f", benchmark_evaluation_seconds),
    statistics_seconds = sprintf("%.3f", statistics_seconds),
    table_preparation_seconds = sprintf("%.3f", table_preparation_seconds),
    cleanup_seconds = sprintf("%.3f", cleanup_seconds),
    accounted_seconds = sprintf("%.3f", accounted_seconds),
    residual_seconds = sprintf("%.3f", max(0, total_seconds - accounted_seconds)),
    total_seconds = sprintf("%.3f", total_seconds)
  )
  list(
    configuration = list(
      application = config$project_paths$application,
      threshold = config$contract$quasi_extinction_abundance,
      horizon_years = config$contract$persistence_horizon_years,
      run_tag = config$project_paths$run_tag,
      benchmark_method = config$rank_method,
      optimization_curves = config$curves,
      main_curve = config$primary_curve,
      late_stage_removal_percentages = config$late_stage_removal_percentages,
      persistence_threshold = config$persistence_threshold
    ),
    plan = plan, statistics = statistics,
    report_tables = tables, elapsed = total_seconds,
    runtime_log_path = config$paths$report_log
  )
}

run_stage74 <- function(config) {
  run_states <- read_completed_stage6_states(config)
  plan <- build_benchmark_target_plan(run_states)
  if (identical(config$mode, "inspect")) {
    return(inspect_benchmark_library(config, plan))
  }
  assert(identical(config$mode, "report"),
         "run_stage74() requires mode = 'report'.")
  need_file(config$paths$species_csv, "persistence species table")
  assert(benchmark_library_complete(config, plan), paste0(
    config$rank_label, " lookup library is incomplete. Complete Stage 7.2 first."
  ))
  with_runtime_log(
    path = config$paths$report_log,
    stage = "7.4",
    operation = "benchmark_report",
    mode = config$mode,
    context = list(
      application = config$project_paths$application,
      run = config$project_paths$run_tag,
      benchmark = config$rank_method,
      curves = paste(config$curves, collapse = ",")
    ),
    code = function() run_stage74_report(config, run_states, plan)
  )
}
