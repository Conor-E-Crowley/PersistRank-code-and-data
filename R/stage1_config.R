# Pure Stage 1 configuration and path contracts.
#
# Loaded by stage1_workflow.R before operational input/artifact modules. Public
# Rmd parameters are converted to one validated list; sourcing and configuration
# construction perform no filesystem access, package loading, fitting, or writes.


stage1_paths <- function(project = project_paths(), params) {
  selected_input <- function(name, default) {
    project_input_path(
      params[[name]], default, project$root, paste0("params$", name)
    )
  }
  list(
    input_dir = project$raw,
    data_dir = project$data,
    figure_dir = project$figures,
    si_dir = project$si_figures,
    raw = list(
      mammal_rmax = selected_input(
        "mammal_rmax_file", file.path(project$raw, "mammal_rmax.txt")
      ),
      sigma = selected_input("sigma_file", file.path(project$raw, "sigma.csv")),
      bird_traits = selected_input(
        "bird_traits_file", file.path(project$raw, "bird_data.txt")
      ),
      bird_generation_lengths = selected_input(
        "bird_generation_lengths_file",
        file.path(project$raw, "cobi13486-sup-0004-tables4.xlsx")
      ),
      bird_growth = selected_input(
        "bird_growth_file",
        file.path(project$raw, "bird_growth_niel_lebreton.csv")
      ),
      iucn_bird_synonyms = selected_input(
        "iucn_bird_synonyms_file",
        file.path(project$raw, "iucn_bird_synonyms.csv")
      ),
      input_synonyms = selected_input(
        "curated_synonyms_file", file.path(project$raw, "input_synonyms.csv")
      )
    ),
    clean = list(
      mammal_growth_posterior =
        file.path(project$clean, "mammal_growth_posterior.csv"),
      mammal_environmental_variation_posterior =
        file.path(project$clean, "mammal_environmental_variation_posterior.csv"),
      bird_growth_posterior =
        file.path(project$clean, "bird_growth_posterior.csv"),
      bird_environmental_variation_posterior =
        file.path(project$clean, "bird_environmental_variation_posterior.csv"),
      mammal_mass_grid = file.path(project$clean, "mammal_mass_grid.csv"),
      bird_generation_length_grid =
        file.path(project$clean, "bird_generation_length_grid.csv")
    ),
    figures = list(
      demographic_calibration = file.path(project$figures, "demographic_calibration.png"),
      demographic_posteriors =
        file.path(project$si_figures, "S1_demographic_posteriors.png")
    )
  )
}

stage1_config <- function(params, paths = project_paths()) {
  mode <- validate_scalar_choice(params$mode, c("reuse", "fit"), "params$mode")
  fit_models <- identical(mode, "fit")
  bird_model <- bird_demographic_model_spec(
    normalize_bird_separate_intercepts(
      params$bird_sigma_separate_intercepts,
      "params$bird_sigma_separate_intercepts"
    )
  )

  n_chains <- validate_scalar_integer(
    params$n_chains, "params$n_chains", if (fit_models) 2L else 1L
  )
  n_adapt <- validate_scalar_integer(params$n_adapt, "params$n_adapt")
  n_iter <- validate_scalar_integer(params$n_iter, "params$n_iter")
  thin <- validate_scalar_integer(params$thin, "params$thin")
  assert(
    n_iter %% thin == 0L,
    "params$n_iter must be divisible by params$thin."
  )

  prior_sd <- validate_scalar_number(
    params$coefficient_prior_sd,
    "params$coefficient_prior_sd",
    minimum = 0,
    minimum_open = TRUE
  )
  residual_min <- validate_scalar_number(
    params$residual_sd_min, "params$residual_sd_min", 0
  )
  residual_max <- validate_scalar_number(
    params$residual_sd_max,
    "params$residual_sd_max",
    minimum = 0,
    minimum_open = TRUE
  )
  assert(
    residual_max > residual_min,
    "params$residual_sd_max must exceed params$residual_sd_min."
  )
  mass_min <- validate_scalar_number(
    params$mammal_mass_min_g,
    "params$mammal_mass_min_g",
    minimum = 0,
    minimum_open = TRUE
  )
  mass_max <- validate_scalar_number(
    params$mammal_mass_max_g,
    "params$mammal_mass_max_g",
    minimum = 0,
    minimum_open = TRUE
  )
  assert(
    mass_max > mass_min,
    "params$mammal_mass_max_g must exceed params$mammal_mass_min_g."
  )

  bayes <- list(
    n_chains = n_chains,
    n_adapt = n_adapt,
    n_iter = n_iter,
    thin = thin,
    seed = validate_scalar_integer(params$seed, "params$seed", 0L),
    prior_mean = validate_scalar_number(
      params$coefficient_prior_mean, "params$coefficient_prior_mean"
    ),
    prior_sd = prior_sd,
    tau_beta = 1 / prior_sd^2,
    sigma_min = residual_min,
    sigma_max = residual_max,
    rhat_max = validate_scalar_number(
      params$rhat_max, "params$rhat_max", minimum = 1
    ),
    ess_min = validate_scalar_number(
      params$ess_min,
      "params$ess_min",
      minimum = 0,
      minimum_open = TRUE
    )
  )

  list(
    mode = mode,
    fit_models = fit_models,
    verbose = validate_scalar_logical(params$verbose, "params$verbose"),
    bird_model = bird_model,
    paths = stage1_paths(paths, params),
    bayes = bayes,
    expected_fit_draws = as.integer(n_chains * n_iter / thin),
    grid = list(
      points = validate_scalar_integer(
        params$trait_grid_points, "params$trait_grid_points", 5L
      ),
      mammal_mass_min_g = mass_min,
      mammal_mass_max_g = mass_max
    )
  )
}

required_stage1_packages <- function() {
  c(
    "dplyr", "readr", "readxl", "tibble", "stringr",
    "ggplot2", "cowplot", "scales"
  )
}
