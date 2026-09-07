# Stage 7.3 orchestration.
#
# The workflow loads each dependency once, constructs or validates the four
# reusable comparison tables, derives all report statistics, and writes one
# seven-panel figure. Large method-specific objects are released before the
# function returns to keep peak memory bounded.

if (!exists("project_source", mode = "function")) {
  source(file.path("R", "project_loader.R"))
}
project_source(c(
  # Shared contracts, benchmark readers, and presentation.
  "R/project_utils.R", "R/project_artifacts.R", "R/project_rmd.R",
  "R/project_runtime_log.R", "R/project_transactions.R",
  "R/figure_utils.R", "R/project_paths.R", "R/project_config.R",
  "R/demographic_contract.R", "R/analysis_contract.R",
  "R/application_storage.R", "R/storage_artifacts.R", "R/application_context.R",
  "R/report_config.R",
  "R/priority_run_config.R", "R/stage73_config.R",
  "R/benchmark_config.R", "R/stage7_contracts.R", "R/benchmark_targets.R",
  "R/benchmark_library.R", "R/stage7_persistence.R",
  "R/benchmark_evaluation.R", "R/benchmark_statistics.R",
  "R/benchmark_tables.R", "R/benchmark_reconstruction.R",
  "R/benchmark_workflow.R",
  # Stage 7.3 comparison, statistics, figures, and report.
  "R/stage73_comparison.R", "R/stage73_statistics.R",
  "R/stage73_figure_data.R", "R/stage73_focal_figures.R",
  "R/stage73_assemblage_figures.R", "R/stage73_figures.R",
  "R/stage73_report.R", "R/stage73_report_pipeline.R"
))

write_stage73_figure_atomic <- function(figure, config) {
  target <- validate_path_param(config$paths$figure, "Stage 7.3 figure path")
  figure_id <- stage73_figure_id(config$rank_method)
  target_existed <- file.exists(target)
  ensure_writable_dir(dirname(target), "Stage 7.3 figure directory")

  staging_dir <- tempfile("stage73_figure_", tmpdir = dirname(target))
  dir.create(staging_dir)
  on.exit(unlink(staging_dir, recursive = TRUE, force = TRUE), add = TRUE)
  staged <- suppressMessages(save_manuscript_figure(
    figure,
    figure_id,
    figure_dir = staging_dir
  ))
  project_file_set_transaction(
    staged,
    target,
    overwrite = TRUE,
    label = "Stage 7.3 combined figure"
  )
  log_msg(
    "FIGURE | written=", target,
    " method=", config$rank_method,
    " existing_replaced=", target_existed
  )
  invisible(target)
}

run_stage73 <- function(config) {
  assert(
    identical(config$mode, "report"),
    "run_stage73() requires mode = report."
  )
  started <- proc.time()[["elapsed"]]
  report_result <- build_stage73_report_result(config)
  figure_path <- write_stage73_figure_atomic(report_result$figure, config)
  statistics <- report_result$statistics
  report_tables <- report_result$report_tables
  rm(report_result)
  gc(FALSE)

  elapsed <- proc.time()[["elapsed"]] - started
  log_msg(
    "STAGE 7.3 COMPLETE | figure=1",
    " method=", config$rank_method,
    " focal_species=", paste(config$focal_species, collapse = ";"),
    " elapsed_s=", sprintf("%.1f", elapsed)
  )
  list(
    statistics = statistics,
    report_tables = report_tables,
    focal_species = config$focal_species,
    figure_written = TRUE,
    figure_path = figure_path,
    elapsed = elapsed
  )
}
