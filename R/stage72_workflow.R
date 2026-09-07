# Stage 7.2 chronological single-method lookup workflow.

if (!exists("project_source", mode = "function")) {
  source(file.path("R", "project_loader.R"))
}
project_source(c(
  # Shared benchmark reporting definitions; spatial reconstruction stays lazy.
  "R/project_utils.R", "R/project_artifacts.R", "R/project_rmd.R",
  "R/project_runtime_log.R", "R/project_transactions.R",
  "R/project_paths.R", "R/project_config.R", "R/demographic_contract.R",
  "R/analysis_contract.R", "R/application_storage.R", "R/storage_artifacts.R",
  "R/application_context.R",
  "R/report_config.R",
  "R/priority_run_config.R", "R/benchmark_config.R", "R/stage72_config.R",
  "R/stage7_contracts.R", "R/benchmark_targets.R", "R/benchmark_library.R",
  "R/stage7_persistence.R", "R/benchmark_evaluation.R",
  "R/benchmark_statistics.R", "R/benchmark_tables.R",
  "R/benchmark_reconstruction.R", "R/benchmark_workflow.R"
))

stage72_target_plan <- function(config) {
  build_benchmark_target_plan(read_completed_stage6_states(config))
}

inspect_stage72 <- function(config, plan) {
  inspect_benchmark_library(config, plan)
}

run_stage72 <- function(config, plan) {
  if (identical(config$mode, "restart")) restart_benchmark_library(config)
  resume_benchmark_library(config, plan)
}
