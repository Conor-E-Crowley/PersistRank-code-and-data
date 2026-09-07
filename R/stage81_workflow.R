# Stage 8.1 multi-benchmark orchestration.
#
# The public operations mirror the Rmd chunks. Target planning and inspection
# are read-only. Lookup resumption is the only spatially expensive operation.
# Results and the optional detailed figure are independently regenerated from
# disk, so neither depends on objects created by the lookup chunk. Missing-target
# reconstruction uses method logs; multi-method reporting uses one run-level
# Stage 8.1 operational log.

if (!exists("project_source", mode = "function")) {
  source(file.path("R", "project_loader.R"))
}
project_source(c(
  # Shared benchmark planning and read-only reporting definitions.
  "R/project_utils.R", "R/project_artifacts.R", "R/project_rmd.R",
  "R/project_runtime_log.R", "R/project_transactions.R",
  "R/project_paths.R", "R/project_config.R", "R/demographic_contract.R",
  "R/analysis_contract.R", "R/application_storage.R", "R/storage_artifacts.R",
  "R/application_context.R",
  "R/report_config.R",
  "R/priority_run_config.R", "R/benchmark_config.R", "R/stage81_config.R",
  "R/stage7_contracts.R", "R/benchmark_targets.R", "R/benchmark_library.R",
  "R/stage7_persistence.R", "R/benchmark_evaluation.R",
  "R/benchmark_statistics.R", "R/benchmark_tables.R",
  "R/benchmark_reconstruction.R", "R/benchmark_workflow.R"
))

# Read completed Stage 6 state indexes and return the shared consumer/target plan.
stage81_target_plan <- function(config) {
  shared <- stage81_benchmark_config(config, "abf")
  states <- read_completed_stage6_states(shared)
  build_benchmark_target_plan(states)
}

# Return one lightweight readiness row per fixed benchmark method. No exact
# lookup table or reconstruction source is opened.
inspect_stage81 <- function(config, plan) {
  rows <- lapply(config$methods, function(method) {
    method_config <- stage81_benchmark_config(config, method)
    inspect_benchmark_library(method_config, plan)
  })
  do.call(rbind, rows)
}

# Non-destructively complete missing exact targets method by method. The
# function fails unless mode is resume and never exposes a restart operation.
resume_stage81_lookups <- function(config, plan) {
  assert(identical(config$mode, "resume"),
         "Stage 8.1 lookup reconstruction requires mode = resume.")
  rows <- lapply(config$methods, function(method) {
    method_config <- stage81_benchmark_config(config, method)
    result <- resume_benchmark_library(method_config, plan)
    data.frame(
      benchmark = method,
      benchmark_label = unname(config$method_labels[[method]]),
      status = result$status,
      created = result$created,
      reused = result$reused,
      elapsed_seconds = result$elapsed_seconds %||% 0,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

stage81_primary_rows <- function(method, label, statistics) {
  add <- function(table) data.frame(
    benchmark = method, benchmark_label = label,
    as.data.frame(table), check.names = FALSE, stringsAsFactors = FALSE
  )
  list(
    late = add(statistics$primary_late),
    sequence = add(statistics$primary_sequence),
    retained_area = add(statistics$primary_area)
  )
}

stage81_bind_primary <- function(method_results, field) {
  do.call(rbind, lapply(method_results, function(result) result$primary[[field]]))
}

# Regenerate trajectories, statistics, and plain tables without touching the
# reconstruction runtime. Only the shared pipeline trajectory remains resident
# while benchmark methods are read and released one at a time.
build_stage81_results_active <- function(config, plan) {
  # calculate_benchmark_statistics() (shared with Stage 7.4) reads the primary
  # comparison curve as config$primary_curve. Stage 8.1's own config exposes
  # that same value as config$main_curve instead, to match the Rmd-facing
  # parameter name shared with Stages 7.3 and 9. Alias it here, on this local
  # copy only, so the calculate_benchmark_statistics() call below finds it.
  config$primary_curve <- config$main_curve

  report_started <- proc.time()[["elapsed"]]
  setup_started <- proc.time()[["elapsed"]]
  shared <- stage81_benchmark_config(config, "abf")
  run_states <- read_completed_stage6_states(shared)
  expected <- build_benchmark_target_plan(run_states)
  assert(benchmark_target_plans_equal(plan$consumers, expected$consumers),
         "The supplied Stage 8.1 target plan is stale.")
  for (method in config$methods) {
    method_config <- stage81_benchmark_config(config, method)
    assert(benchmark_library_complete(method_config, plan), paste0(
      unname(config$method_labels[[method]]),
      " lookup library is incomplete. Run the lookup-libraries chunk first."
    ))
  }
  reference <- load_benchmark_initialization(shared)
  parameters <- prepare_benchmark_parameters(reference, shared)
  shared_setup_seconds <- proc.time()[["elapsed"]] - setup_started

  pipeline_started <- proc.time()[["elapsed"]]
  pipeline <- evaluate_pipeline_states(
    shared, plan, parameters, run_states, reference
  )
  rm(reference)
  invisible(gc(FALSE))
  pipeline_evaluation_seconds <- proc.time()[["elapsed"]] - pipeline_started

  method_results <- stats::setNames(vector("list", length(config$methods)), config$methods)
  method_total_seconds <- 0
  for (method in config$methods) {
    method_started <- proc.time()[["elapsed"]]
    method_config <- stage81_benchmark_config(config, method)
    phase_started <- proc.time()[["elapsed"]]
    benchmark <- evaluate_benchmark_library(method_config, plan, parameters)
    evaluation_seconds <- proc.time()[["elapsed"]] - phase_started
    phase_started <- proc.time()[["elapsed"]]
    statistics <- calculate_benchmark_statistics(pipeline, benchmark, config)
    statistics_seconds <- proc.time()[["elapsed"]] - phase_started
    phase_started <- proc.time()[["elapsed"]]
    tables <- build_benchmark_report_tables(statistics, config)
    table_preparation_seconds <- proc.time()[["elapsed"]] - phase_started
    phase_started <- proc.time()[["elapsed"]]
    method_results[[method]] <- list(
      benchmark = method,
      benchmark_label = unname(config$method_labels[[method]]),
      statistics = statistics,
      report_tables = tables,
      primary = stage81_primary_rows(
        method, unname(config$method_labels[[method]]), statistics
      )
    )
    result_preparation_seconds <- proc.time()[["elapsed"]] - phase_started
    cleanup_started <- proc.time()[["elapsed"]]
    rm(benchmark, statistics, tables)
    invisible(gc(FALSE))
    cleanup_seconds <- proc.time()[["elapsed"]] - cleanup_started
    total_seconds <- proc.time()[["elapsed"]] - method_started
    accounted_seconds <- sum(c(
      evaluation_seconds, statistics_seconds, table_preparation_seconds,
      result_preparation_seconds, cleanup_seconds
    ))
    method_total_seconds <- method_total_seconds + total_seconds
    runtime_log_event(
      "stage81_method_timing",
      benchmark = method,
      evaluation_seconds = sprintf("%.3f", evaluation_seconds),
      statistics_seconds = sprintf("%.3f", statistics_seconds),
      table_preparation_seconds = sprintf("%.3f", table_preparation_seconds),
      result_preparation_seconds = sprintf("%.3f", result_preparation_seconds),
      cleanup_seconds = sprintf("%.3f", cleanup_seconds),
      accounted_seconds = sprintf("%.3f", accounted_seconds),
      residual_seconds = sprintf("%.3f", max(0, total_seconds - accounted_seconds)),
      total_seconds = sprintf("%.3f", total_seconds)
    )
  }
  final_started <- proc.time()[["elapsed"]]
  rm(pipeline, parameters, run_states)
  invisible(gc(FALSE))
  result <- list(
    methods = method_results,
    cross_benchmark = list(
      late = stage81_bind_primary(method_results, "late"),
      sequence = stage81_bind_primary(method_results, "sequence"),
      retained_area = stage81_bind_primary(method_results, "retained_area")
    ),
    target_diagnostics = plan$diagnostics,
    runtime_log_path = config$runtime_log
  )
  final_seconds <- proc.time()[["elapsed"]] - final_started
  total_seconds <- proc.time()[["elapsed"]] - report_started
  accounted_seconds <- sum(c(
    shared_setup_seconds, pipeline_evaluation_seconds,
    method_total_seconds, final_seconds
  ))
  runtime_log_event(
    "stage81_report_timing",
    methods = paste(config$methods, collapse = ","),
    shared_setup_seconds = sprintf("%.3f", shared_setup_seconds),
    pipeline_evaluation_seconds = sprintf("%.3f", pipeline_evaluation_seconds),
    method_seconds = sprintf("%.3f", method_total_seconds),
    final_preparation_seconds = sprintf("%.3f", final_seconds),
    accounted_seconds = sprintf("%.3f", accounted_seconds),
    residual_seconds = sprintf("%.3f", max(0, total_seconds - accounted_seconds)),
    total_seconds = sprintf("%.3f", total_seconds)
  )
  result
}

build_stage81_results <- function(config, plan) {
  assert(identical(config$mode, "resume"),
         "Stage 8.1 result generation requires mode = resume.")
  with_runtime_log(
    path = config$runtime_log,
    stage = "8.1",
    operation = "multi_benchmark_report",
    mode = config$mode,
    context = list(
      application = config$application,
      run = config$scenario$run_tag,
      benchmarks = paste(config$methods, collapse = ","),
      curves = paste(config$curves, collapse = ",")
    ),
    code = function() build_stage81_results_active(config, plan)
  )
}

# Regenerate the detailed ABF/main-curve figure separately from report tables.
build_stage81_figure <- function(config, plan) {
  # The optional detailed figure is loaded only when explicitly requested.
  project_source(c(
    "R/figure_utils.R", "R/report_config.R", "R/stage73_config.R",
    "R/stage73_comparison.R",
    "R/stage73_statistics.R", "R/stage73_figure_data.R",
    "R/stage73_focal_figures.R", "R/stage73_assemblage_figures.R",
    "R/stage73_figures.R", "R/stage73_report.R",
    "R/stage73_report_pipeline.R",
    "R/stage73_workflow.R"
  ))
  figure_config <- stage73_config(list(
    mode = "report", curve = config$main_curve, taxa = config$taxa,
    sdm = config$sdm,
    cells_to_remove_per_iteration = config$cells_to_remove_per_iteration,
    pruning_iterations_per_stage = config$pruning_iterations_per_stage,
    rank_method = "abf", persistence_threshold = config$persistence_threshold,
    focal_species = config$focal_species,
    assemblage_point_style = config$assemblage_point_style,
    focal_point_style = config$focal_point_style,
    assemblage_central_stat = config$assemblage_central_stat
  ), config$project_paths, config$contract)
  states <- read_completed_stage6_states(stage81_benchmark_config(config, "abf"))
  expected_plan <- build_benchmark_target_plan(states)
  assert(benchmark_target_plans_equal(plan$consumers, expected_plan$consumers),
         "The supplied Stage 8.1 target plan is stale.")
  assert(benchmark_library_complete(
    stage81_benchmark_config(config, "abf"), plan
  ), "The detailed figure requires a complete ABF lookup library.")
  rm(states, expected_plan)
  result <- run_stage73(figure_config)
  focal_results <- data.table::as.data.table(result$statistics$auc$species)[
    scientificName %in% config$focal_species
  ]
  list(
    focal_results = focal_results,
    focal_species = result$focal_species,
    figure_path = result$figure_path,
    elapsed_seconds = result$elapsed
  )
}

stage81_artifact_summary <- function(config, plan, figure = NULL) {
  rows <- do.call(rbind, lapply(config$methods, function(method) {
    method_config <- stage81_benchmark_config(config, method)
    index <- benchmark_library_index(method_config, plan, validate = FALSE)
    data.frame(
      benchmark = method,
      exact_lookup_files = sum(index$exists),
      lookup_bytes = sum(index$bytes, na.rm = TRUE),
      lookup_directory = method_config$paths$lookup_library$lookups,
      stringsAsFactors = FALSE
    )
  }))
  reconstruction_logs <- do.call(rbind, lapply(config$methods, function(method) {
    path <- stage81_benchmark_config(config, method)$paths$reconstruction_log
    data.frame(
      operation = "lookup reconstruction", benchmark = method,
      path = path, exists = file.exists(path), stringsAsFactors = FALSE
    )
  }))
  report_logs <- data.frame(
    operation = "report", benchmark = "all",
    path = config$runtime_log, exists = file.exists(config$runtime_log),
    stringsAsFactors = FALSE
  )
  list(
    lookups = rows,
    runtime_logs = rbind(reconstruction_logs, report_logs),
    figure = if (is.null(figure)) data.frame() else data.frame(
      artifact = "detailed ABF figure", path = figure$figure_path,
      exists = file.exists(figure$figure_path), stringsAsFactors = FALSE
    )
  )
}
