# Reusable Stage 5.1 spatial-process figure workflow.

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
  # Stage 5.1 spatial and presentation definitions.
  "R/patch_species_inputs.R", "R/patch_habitat.R",
  "R/patch_processing.R", "R/patch_process_figure.R"
))

run_stage51 <- function(config) {
  if (is.null(config$target_species)) {
    log_msg("Stage 5.1 | status=not_applicable | process_figure_species is null")
    return(list(status = "not_applicable", species = NA_character_,
                figure = character(), elapsed_seconds = 0))
  }
  check_stage51_preflight(config)
  load_patch_packages(include_figures = TRUE)
  started <- Sys.time()
  species_table <- read_patch_species_base(
    config$paths$species_csv, contract = config$contract
  )
  target_row <- species_table[
    species_table$scientificName == config$target_species, , drop = FALSE
  ]
  assert(nrow(target_row) == 1L, paste0(
    "Stage 5.1 focal species must match exactly one species-table row: ",
    config$target_species, ". Found ", nrow(target_row), "."
  ))
  # Load optional GRASS/fasterRaster definitions only for applicable active work.
  project_source("R/patch_runtime.R")
  resolved_backend <- resolve_patch_clump_backend(config$clump_backend, config$grass_dir, "Stage 5.1")
  layers <- build_single_species_layers(
    target_row, paths = config$paths, roi = config$study_area,
    clump_backend = resolved_backend
  )
  process_figure <- make_single_species_process_figure(layers)
  write_stage51_figure_atomic(process_figure, config)
  summary <- list(
    status = "complete", species = target_row$scientificName,
    mapped_area_km2 = layers$mapped_habitat_area_km2,
    raw_patches = layers$n_raw_patches, kept_patches = layers$n_kept_patches,
    candidate_population_units = layers$n_candidate_pus,
    final_population_units = layers$n_final_pus,
    source_raster = target_row$raster_path,
    figure = config$paths$figure,
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))
  )
  rm(species_table, target_row, layers, process_figure)
  summary
}
