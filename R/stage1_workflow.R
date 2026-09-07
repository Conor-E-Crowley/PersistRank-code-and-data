# Stage 1 demographic-model workflow.
#
# Upstream: raw calibration tables and a validated stage1_config().
# Downstream: Stage 2 consumes the installed posterior tables and trait grids.
# Sourcing this module only loads definitions. run_stage1() owns preflight,
# fitting/reuse, figure construction, the atomic output transaction, and the
# demographic-model manifest; it returns compact summaries for the Rmd.

if (!exists("project_source", mode = "function")) {
  source(file.path("R", "project_loader.R"))
}
project_source(c(
  # Shared contracts and storage.
  "R/project_utils.R", "R/project_artifacts.R", "R/project_rmd.R",
  "R/project_runtime_log.R", "R/project_transactions.R",
  "R/project_paths.R", "R/project_config.R", "R/demographic_contract.R",
  "R/analysis_contract.R", "R/application_storage.R", "R/storage_artifacts.R",
  "R/model_lifecycle.R", "R/figure_utils.R",
  # Stage 1 matching, calibration, fitting, presentation, and output.
  "R/bird_generation_lengths.R", "R/stage1_config.R", "R/stage1_artifacts.R",
  "R/stage1_bird_matching.R", "R/stage1_data.R", "R/stage1_models.R",
  "R/stage1_figures.R", "R/stage1_posterior_figures.R",
  "R/stage1_outputs.R"
))

# Run the complete Stage 1 operation without exposing intermediate scientific
# objects to an interactive document. Reuse mode still reconstructs the trait
# grids and figures from validated posterior artifacts, exactly as fit mode
# does, but never invokes JAGS.
run_stage1 <- function(config, model = demography_model_paths()) {
  assert(
    is.list(config) && is.list(config$paths) && is.list(config$bayes),
    "run_stage1() requires a stage1_config() result."
  )
  assert(
    is.list(model) && !is.null(model$root) && !is.null(model$manifest),
    "run_stage1() requires canonical demography_model_paths()."
  )

  started <- Sys.time()
  paths <- config$paths
  bayes <- config$bayes
  check_stage1_preflight(config)
  load_packages(required_stage1_packages())
  assert_jags_available(config$fit_models)

  log_msg(
    "Stage 1 | start | mode=", config$mode,
    " | separate bird intercepts=",
    if (length(config$bird_model$separate_intercepts)) {
      paste(config$bird_model$separate_intercepts, collapse = ",")
    } else {
      "none"
    },
    " | chains=", bayes$n_chains,
    " | posterior draws=",
    if (config$fit_models) config$expected_fit_draws else "validate saved files"
  )

  inputs <- read_calibration_inputs(paths)
  calibration <- prepare_calibration_data(
    inputs,
    bird_model = config$bird_model,
    verbose = config$verbose
  )
  diet_counts <- calibration$bird_sigma_data |>
    dplyr::count(diet5_group, bird_sigma_model_group, name = "n") |>
    dplyr::arrange(match(diet5_group, bird_diet5_categories()))

  bird_range <- range(inputs$bird_generation_lengths$GenLength, finite = TRUE)
  trait_grids <- list(
    mammal = tibble::tibble(
      Mass_g = log_space(
        config$grid$mammal_mass_min_g,
        config$grid$mammal_mass_max_g,
        config$grid$points
      )
    ),
    bird = tibble::tibble(
      GenLength = log_space(
        bird_range[[1L]], bird_range[[2L]], config$grid$points
      )
    )
  )

  if (config$fit_models) {
    fits <- fit_all_demographic_models(
      calibration,
      bird_model = config$bird_model,
      bayes = bayes
    )
    posterior_draws <- lapply(fits, `[[`, "draws")
    posterior_draw_count <- config$expected_fit_draws
  } else {
    posterior_draws <- read_stage1_posteriors(paths, config$bird_model)
    posterior_draw_count <- attr(posterior_draws, "posterior_draws")
  }

  log_posterior_intervals(posterior_draws)
  stage1_r2 <- descriptive_stage1_r2(calibration, config$bird_model)
  figures <- list(
    demographic_calibration = build_demographic_calibration_figure(
      calibration,
      trait_grids,
      posterior_draws,
      bird_model = config$bird_model,
      stage1_r2 = stage1_r2
    ),
    demographic_posteriors = build_stage1_coefficient_figure(posterior_draws)
  )

  if (config$fit_models) {
    write_stage1_fit_outputs(
      posterior_draws,
      trait_grids,
      figures = figures,
      paths = paths,
      bird_model = config$bird_model,
      expected_draws = config$expected_fit_draws
    )
  } else {
    write_stage1_reuse_outputs(config, prepared = trait_grids, figures = figures)
  }

  manifest <- write_demography_model_manifest(
    model,
    unlist(paths$clean, use.names = FALSE)
  )
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))

  # Large calibration, posterior, fit, and plot objects are intentionally not
  # returned. Their validated durable representations are the Stage 1 contract.
  list(
    action = config$mode,
    calibration_summary = calibration$summary,
    bird_diet_counts = diet_counts,
    posterior_draws_per_model = as.integer(posterior_draw_count),
    bird_model_branches = config$bird_model$model_groups,
    output_files = unname(unlist(paths$clean, use.names = FALSE)),
    figure_files = unname(unlist(paths$figures, use.names = FALSE)),
    manifest = manifest,
    elapsed_seconds = elapsed
  )
}
