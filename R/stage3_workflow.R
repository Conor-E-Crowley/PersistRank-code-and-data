# Reusable Stage 3 orchestration shared by the standalone report and Stage 8.

if (!exists("project_source", mode = "function")) {
  source(file.path("R", "project_loader.R"))
}
project_source(c(
  # Shared contracts, storage, and presentation.
  "R/project_utils.R", "R/project_artifacts.R", "R/project_rmd.R",
  "R/project_runtime_log.R", "R/project_transactions.R",
  "R/project_paths.R", "R/project_config.R", "R/demographic_contract.R",
  "R/analysis_contract.R", "R/application_storage.R", "R/storage_artifacts.R",
  "R/model_lifecycle.R", "R/figure_utils.R",
  # Stage 3 fitting, smoothing, contracts, presentation, and output.
  "R/stage3_config.R", "R/stage3_inputs.R", "R/stage3_gompertz.R",
  "R/stage3_smoothing.R",
  "R/stage3_models.R", "R/stage3_figures.R",
  "R/stage3_parameter_figures.R", "R/stage3_persistence_figures.R",
  "R/stage3_outputs.R"
))

run_stage3 <- function(config) {
  assert(is.list(config) && !is.null(config$loess), "run_stage3() requires a stage3_config() result.")
  started <- Sys.time()
  check_stage3_preflight(config)
  load_gompertz_packages(config$write_figures)
  inputs <- load_stage3_inputs(config)
  config <- resolve_stage3_input_settings(config, inputs)
  persistence_points <- NULL
  loess_predictions <- NULL
  input_summary <- dplyr::bind_rows(lapply(names(Filter(Negate(is.null), inputs)), function(taxon) {
    item <- inputs[[taxon]]
    tibble::tibble(taxon = taxon, rows = nrow(item$points), traits = item$trait_count,
                   blocks = item$block_count, md5 = item$md5)
  }))
  if (config$mode == "fit") {
    persistence_points <- bind_stage3_points(inputs)
    gompertz_parameters <- fit_all_gompertz_parameters(persistence_points, config$k0)
    assert_all_gompertz_fits(gompertz_parameters)
    persistence_models <- build_gompertz_loess_models(gompertz_parameters, inputs, config)
  } else {
    persistence_models <- readRDS(config$paths$loess_models_rds)
    validate_gompertz_loess_models(
      persistence_models, config$selected_mammals, config$selected_birds,
      config$curves, config$k0
    )
    validate_stage3_reuse_provenance(persistence_models, inputs, config)
    gompertz_parameters <- persistence_models$gompertz_fits
  }
  bird_branches <- if (config$selected_birds) {
    persistence_models$meta$bird_model_spec$model_groups
  } else {
    character()
  }
  figures <- list()
  figure_targets <- character()
  if (isTRUE(config$write_figures)) {
    loess_predictions <- predict_saved_stage3_models(persistence_models, config)
    figure_targets <- c(config$paths$selected_figures)
    if (config$selected_mammals) {
      if (is.null(persistence_points)) persistence_points <- bind_stage3_points(inputs)
      figures$gompertz_wolff <- build_gompertz_wolff_figure(
        persistence_points, gompertz_parameters, k0 = config$k0, curve = config$main_curve
      )
    }
    figures$gompertz_parameters <- build_gompertz_loess_combined_figure(
      gompertz_parameters, loess_predictions, config$curves
    )
    figures$abundance_persistence <- build_abundance_persistence_figure(
      loess_predictions,
      k0 = config$k0,
      curve = config$heatmap_curve,
      include_mammals = config$selected_mammals,
      bird_model_groups = bird_branches
    )
  }
  if (identical(config$mode, "fit") || length(figures)) {
    write_stage3_outputs(persistence_models, figures, figure_targets, config)
  }
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  rm(inputs, persistence_points, gompertz_parameters, loess_predictions, persistence_models, figures)
  list(
    action = config$mode, model = config$paths$loess_models_rds,
    figures = unname(figure_targets), input_summary = input_summary,
    quasi_extinction_abundance = config$k0,
    persistence_horizon_years = config$expected_horizon,
    bird_branches = bird_branches, elapsed_seconds = elapsed
  )
}

# Run Stage 3, publish its recovery boundary, and install the reusable model
# manifest as one lifecycle operation shared by the standalone Rmd and Stage 8.
run_stage3_model <- function(config, model) {
  assert(
    is.list(model) && !is.null(model$root) && !is.null(model$manifest),
    "run_stage3_model() requires canonical persistence_model_paths()."
  )
  result <- run_stage3(config)
  update_persistence_model_state(
    model, "stage_3", result$model,
    list(
      persistence_horizon_years = result$persistence_horizon_years,
      quasi_extinction_abundance = result$quasi_extinction_abundance
    )
  )
  write_persistence_model_manifest_from_artifacts(
    model, result$persistence_horizon_years,
    result$quasi_extinction_abundance, config$loess$span
  )
  result
}
