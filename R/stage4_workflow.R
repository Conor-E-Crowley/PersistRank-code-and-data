# Reusable Stage 4 species-table workflow with a shared transactional IUCN cache.

if (!exists("project_source", mode = "function")) {
  source(file.path("R", "project_loader.R"))
}
project_source(c(
  # Shared contracts and application adapters.
  "R/project_utils.R", "R/project_artifacts.R", "R/project_rmd.R",
  "R/project_runtime_log.R", "R/project_transactions.R",
  "R/project_paths.R", "R/project_config.R", "R/demographic_contract.R",
  "R/analysis_contract.R", "R/priority_run_config.R",
  "R/species_table_config.R", "R/patch_contract.R", "R/patch_config.R",
  "R/application_storage.R", "R/storage_artifacts.R", "R/model_lifecycle.R",
  "R/application_state.R",
  "R/application_context.R", "R/application_rmd.R", "R/figure_utils.R",
  # Stage 4 inputs, matching, traits, build, and figures.
  "R/bird_generation_lengths.R", "R/species_table_inputs.R",
  "R/species_name_candidates.R", "R/species_name_diagnostics.R",
  "R/species_table_names.R", "R/species_trait_covariates.R",
  "R/species_gompertz_parameters.R", "R/species_table_traits.R",
  "R/species_table_iucn.R", "R/species_table_build.R",
  "R/species_table_figures.R"
))

run_stage4 <- function(config) {
  assert(is.list(config) && !is.null(config$contract), "run_stage4() requires a stage4_config() result.")
  started <- Sys.time()
  check_stage4_preflight(config)
  if (config$build_figure_only) {
    load_species_table_figure_packages()
    rebuild_area_curve_figure(
      species_csv = config$paths$clean$species_table,
      focal_species = config$figure_focal_species,
      curve = config$main_curve,
      figure_dir = config$paths$figure_dir,
      overwrite = TRUE
    )
    return(list(
      action = "figure", rows = NA_integer_, outputs = config$paths$area_curve,
      figure = config$paths$area_curve,
      elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))
    ))
  }

  load_species_table_packages()
  core <- read_species_core_inputs(config)
  resolved <- resolve_species_inputs(core, verbose = config$verbose)
  species_traits <- attach_bird_generation_lengths(
    resolved$species_inputs,
    bird_generation_lengths = core$bird_generation_lengths,
    synonyms = core$synonyms,
    input_synonyms = core$input_synonyms,
    verbose = config$verbose
  )
  species_parameters <- add_trait_covariates(
    species_traits,
    random_effects = core$random_effects,
    models = core$models,
    curves = config$curves,
    use_gompertz_mammals = config$fit_mammal_curves,
    use_gompertz_birds = config$fit_bird_curves,
    contract = config$contract
  )
  validate_stage4_species_parameters(
    species_parameters, curves = config$curves,
    use_gompertz_mammals = config$fit_mammal_curves,
    use_gompertz_birds = config$fit_bird_curves
  )
  habitats <- if (!is.null(config$habitat_rows)) {
    rows <- tibble::as_tibble(config$habitat_rows)
    need_cols(rows, c(
      "genusName", "speciesName", "habitat_codes_suitable",
      "habitats_level1", "habitats_mixed"
    ), "frozen application IUCN rows")
    list(rows = rows, queried = 0L)
  } else resolve_species_habitats(
    species_parameters,
    cache_path = config$paths$clean$iucn_cache,
    mode = config$iucn_mode,
    pause_seconds = config$iucn_pause_seconds,
    verbose = config$verbose
  )
  species_table <- assemble_species_table(
    species_parameters, habitats$rows, curves = config$curves,
    use_gompertz_mammals = config$fit_mammal_curves,
    use_gompertz_birds = config$fit_bird_curves
  )
  write_species_table(
    species_table,
    table_path = config$paths$clean$species_table,
    curves = config$curves,
    use_gompertz_mammals = config$fit_mammal_curves,
    use_gompertz_birds = config$fit_bird_curves,
    overwrite = TRUE
  )
  class_counts <- species_table |>
    dplyr::count(className, sdm_method, name = "n")
  rows <- nrow(species_table)
  queried_iucn_species <- habitats$queried
  habitat_rows <- if (is.null(config$habitat_rows)) habitats$rows else NULL
  outputs <- config$paths$clean$species_table
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  rm(core, resolved, species_traits, species_parameters, habitats, species_table)
  list(
    action = "build", rows = rows, class_counts = class_counts,
    outputs = outputs,
    queried_iucn_species = queried_iucn_species,
    habitat_rows = habitat_rows,
    elapsed_seconds = elapsed
  )
}
