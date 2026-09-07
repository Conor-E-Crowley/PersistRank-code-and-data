# Reusable Stage 2 orchestration shared by the standalone report and Stage 8.
# Inspection is read-only. Active simulation appends setup, curve, and trait
# performance diagnostics to the derived model-level runtime log while retaining
# the existing scientific point/checkpoint transactions.

if (!exists("project_source", mode = "function")) {
  source(file.path("R", "project_loader.R"))
}
project_source(c(
  # Shared contracts and storage.
  "R/project_utils.R", "R/project_artifacts.R", "R/project_rmd.R",
  "R/project_runtime_log.R", "R/project_transactions.R",
  "R/project_paths.R", "R/project_config.R", "R/demographic_contract.R",
  "R/analysis_contract.R", "R/application_storage.R", "R/storage_artifacts.R",
  "R/model_lifecycle.R",
  # Stage 2 configuration, inputs, recovery, and simulation.
  "R/stage2_config.R", "R/stage2_runtime.R", "R/stage2_inputs.R",
  "R/stage2_checkpoints.R",
  "R/stage2_simulation.R"
))

run_stage2_impl <- function(config) {
  assert(is.list(config) && !is.null(config$sim), "run_stage2() requires a stage2_config() result.")
  started <- Sys.time()
  setup_started <- unname(proc.time()[["elapsed"]])
  check_stage2_preflight(config)
  load_packages(required_stage2_packages())
  old_omp <- Sys.getenv("OMP_NUM_THREADS", unset = NA_character_)
  on.exit({
    if (is.na(old_omp)) Sys.unsetenv("OMP_NUM_THREADS") else Sys.setenv(OMP_NUM_THREADS = old_omp)
  }, add = TRUE)

  input_started <- unname(proc.time()[["elapsed"]])
  inputs <- load_persist_inputs(config)
  if (!is.null(inputs$birds)) config$bird_model <- inputs$birds$bird_model
  input_seconds <- unname(proc.time()[["elapsed"]] - input_started)
  cpp_info <- NULL
  crn_ctx <- NULL
  kernel_seconds <- 0
  crn_seconds <- 0
  if (config$active) {
    if (!is.null(config$threads)) Sys.setenv(OMP_NUM_THREADS = as.character(config$threads))
    kernel_started <- unname(proc.time()[["elapsed"]])
    cpp_info <- compile_persist_cpp(config$paths, rebuild_cpp = FALSE)
    kernel_seconds <- unname(proc.time()[["elapsed"]] - kernel_started)
    crn_started <- unname(proc.time()[["elapsed"]])
    run_crn_self_test(config$sim, q_levels = config$grid$q_levels)
    crn_ctx <- create_crn_context(config$sim)
    crn_seconds <- unname(proc.time()[["elapsed"]] - crn_started)
  }

  traits <- list(mammals = NULL, birds = NULL)
  metadata <- list(mammals = NULL, birds = NULL)
  if (config$selected_mammals) {
    traits$mammals <- mammal_persistence_trait_table(inputs$mammals$mass)
    metadata$mammals <- stage2_taxon_metadata(inputs$mammals, config)
  }
  if (config$selected_birds) {
    traits$birds <- bird_persistence_trait_table(inputs$birds$genlength, inputs$birds$bird_model)
    metadata$birds <- stage2_taxon_metadata(inputs$birds, config)
  }
  branch_summary <- data.frame(
    taxon = c(if (config$selected_mammals) "mammals", if (config$selected_birds) "birds"),
    trait_branch_rows = c(
      if (config$selected_mammals) nrow(traits$mammals),
      if (config$selected_birds) nrow(traits$birds)
    ),
    stringsAsFactors = FALSE
  )
  inspect_artifacts <- function() dplyr::bind_rows(
    if (config$selected_mammals) dplyr::bind_rows(lapply(
      c(config$paths$mammals_points, config$paths$mammals_partial),
      inspect_stage2_artifact, trait_table = traits$mammals, config = config,
      metadata = metadata$mammals, idx_col = "mass_idx", value_col = "Mass_g"
    )),
    if (config$selected_birds) dplyr::bind_rows(lapply(
      c(config$paths$birds_points, config$paths$birds_partial),
      inspect_stage2_artifact, trait_table = traits$birds, config = config,
      metadata = metadata$birds, idx_col = "genlength_idx", value_col = "GenLength",
      extra_cols = "bird_sigma_model_group"
    ))
  )
  artifact_status <- if (!config$active) inspect_artifacts() else NULL

  if (config$active) {
    setup_total_seconds <- unname(proc.time()[["elapsed"]] - setup_started)
    accounted_seconds <- sum(c(input_seconds, kernel_seconds, crn_seconds))
    runtime_log_event(
      "stage2_setup_timing",
      taxa = config$taxa,
      threads = config$threads %||% "default",
      draws = as.integer(config$sim$n_draws),
      replicates = as.integer(config$sim$reps),
      chunk_size = as.integer(config$sim$chunk_size),
      input_seconds = sprintf("%.3f", input_seconds),
      kernel_seconds = sprintf("%.3f", kernel_seconds),
      crn_context_seconds = sprintf("%.3f", crn_seconds),
      accounted_seconds = sprintf("%.3f", accounted_seconds),
      residual_seconds = sprintf("%.3f", max(0, setup_total_seconds - accounted_seconds)),
      total_seconds = sprintf("%.3f", setup_total_seconds)
    )
  }

  run_taxon <- function(taxon) {
    mammal <- identical(taxon, "mammals")
    run_one_scenario(
      trait_table = traits[[taxon]],
      sampler = if (mammal) make_mammal_sampler(inputs$mammals, config$sim) else make_bird_sampler(inputs$birds, config$sim),
      out_file = config$paths[[paste0(taxon, "_points")]],
      partial_file = config$paths[[paste0(taxon, "_partial")]],
      checkpoint_file = config$paths[[paste0(taxon, "_checkpoint")]],
      taxon = taxon, crn_ctx = crn_ctx, sim = config$sim, grid = config$grid,
      metadata = metadata[[taxon]],
      idx_col = if (mammal) "mass_idx" else "genlength_idx",
      value_col = if (mammal) "Mass_g" else "GenLength",
      extra_cols = if (mammal) character() else "bird_sigma_model_group",
      existing_output = config$mode, label = paste(taxon, "persistence points"),
      verbose = config$verbose
    )
  }
  if (config$active && config$selected_mammals) run_taxon("mammals")
  if (config$active && config$selected_birds) run_taxon("birds")
  if (config$active) artifact_status <- inspect_artifacts()
  outputs <- c(
    if (config$selected_mammals) config$paths$mammals_points,
    if (config$selected_birds) config$paths$birds_points
  )
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  rm(inputs, traits, metadata, crn_ctx)
  list(
    action = if (config$active) "complete" else "inspect",
    paths = outputs, branch_summary = branch_summary,
    artifact_status = artifact_status, cpp_info = cpp_info,
    grid_signature = config$grid$signature, elapsed_seconds = elapsed
  )
}

run_stage2 <- function(config) {
  assert(is.list(config) && !is.null(config$sim),
         "run_stage2() requires a stage2_config() result.")
  if (!isTRUE(config$active)) return(run_stage2_impl(config))
  result <- with_runtime_log(
    path = config$paths$runtime_log,
    stage = "2",
    operation = "persistence_simulation",
    mode = config$mode,
    context = list(
      taxa = config$taxa,
      threshold = as.integer(config$sim$quasi_extinction_abundance),
      horizon_years = as.integer(config$sim$years),
      curves = paste(config$curves, collapse = ",")
    ),
    code = function() run_stage2_impl(config)
  )
  result$runtime_log_path <- config$paths$runtime_log
  result
}

# Run Stage 2 and, for active modes, publish its recoverable model boundary.
# Standalone Rmds and Stage 8 call this same wrapper so the simulation itself
# and its lifecycle commit cannot drift into two orchestration paths.
run_stage2_model <- function(config, model) {
  assert(
    is.list(model) && !is.null(model$root) && !is.null(model$manifest),
    "run_stage2_model() requires canonical persistence_model_paths()."
  )
  result <- run_stage2(config)
  if (isTRUE(config$active)) {
    update_persistence_model_state(
      model, "stage_2", result$paths,
      list(grid_signature = result$grid_signature)
    )
  }
  result
}
