# Standalone-Rmd adapter for shared application contexts and recovery records.
#
# Stage workflows own loading this file. Sourcing it defines adapters only and
# writes nothing; application state changes occur through explicit record and
# finalize calls. It deliberately does not load the Stage 8 workflow.

if (!exists("project_source", mode = "function")) {
  source(file.path("R", "project_loader.R"))
}
project_source(c(
  # Shared contracts and application architecture.
  "R/project_utils.R",
  "R/project_artifacts.R",
  "R/project_rmd.R",
  "R/project_runtime_log.R",
  "R/project_transactions.R",
  "R/project_paths.R",
  "R/demographic_contract.R",
  "R/analysis_contract.R",
  "R/priority_run_config.R",
  "R/species_table_config.R",
  "R/patch_contract.R",
  "R/patch_config.R",
  "R/application_storage.R",
  "R/storage_artifacts.R",
  "R/model_lifecycle.R",
  "R/application_state.R",
  "R/application_context.R",
  "R/application_lifecycle.R"
))

# Convert standalone Rmd parameters to the shared lightweight application
# context. Optional manifest/model validation is read-only and fails on mismatch.
application_rmd_context <- function(params, require_application = TRUE,
                                    require_scenario = FALSE,
                                    paths = project_paths()) {
  identity <- application_run_identity(params, paths)
  application <- identity$application
  threshold <- identity$quasi_extinction_abundance
  horizon <- identity$persistence_horizon_years
  cells <- identity$cells_to_remove_per_iteration
  iterations <- identity$pruning_iterations_per_stage
  scenario <- identity$scenario
  model <- identity$model
  app_manifest <- if (file.exists(scenario$application_manifest)) {
    readRDS(scenario$application_manifest)
  } else NULL
  if (isTRUE(require_application)) {
    assert(!is.null(app_manifest), paste0(
      "Missing application manifest: ", scenario$application_manifest,
      "\nRun Stage 4 in build mode for this application first."
    ))
    assert(identical(app_manifest$schema, application_storage_schema()) &&
             identical(app_manifest$type, "spatial_application"),
           "Unsupported application manifest.")
  }
  minimum_patch <- app_manifest$contract$minimum_patch_abundance %||%
    params$minimum_patch_abundance %||% 10L
  contract <- new_analysis_contract(horizon, threshold, minimum_patch)
  context <- list(
    schema = application_storage_schema(), mode = "resume",
    application = application,
    contract = contract,
    taxa = app_manifest$contract$taxa %||% params$taxa %||% "both",
    sdm = app_manifest$contract$sdm %||% params$sdm %||% "ppm_rangebag",
    cells_to_remove_per_iteration = cells,
    pruning_iterations_per_stage = iterations,
    model = model, scenario = scenario, base_paths = identity$base_paths,
    paths = identity$paths
  )
  if (isTRUE(require_scenario)) {
    validate_persistence_model(context)
    application_validate_scenario_handoff(context)
  }
  context
}

# Build the Stage 4-owned context from requested selectors rather than silently
# substituting values from an existing application manifest.
application_rmd_stage4_context <- function(params, paths = project_paths()) {
  context <- application_rmd_context(
    params, require_application = FALSE, paths = paths
  )
  checksum_inputs <- !identical(params$mode, "figure")
  # Stage 4 owns these inputs. Use the requested values here so an existing
  # application with different inputs is rejected rather than silently reused.
  context$taxa <- validate_taxa_selector(
    params$taxa %||% "both", "params$taxa"
  )$value
  context$sdm <- validate_priority_sdm_tag(
    params$sdm %||% "ppm_rangebag", "params$sdm"
  )
  context$contract <- new_analysis_contract(
    context$contract$persistence_horizon_years,
    context$contract$quasi_extinction_abundance,
    params$minimum_patch_abundance %||% 10L
  )
  context$inputs <- list(
    species_list = application_selected_file(
      params$species_list_file, file.path(paths$raw, "simple_summary.csv"),
      paths, "species list", checksum = checksum_inputs
    ),
    synonyms = application_selected_file(
      params$synonyms_file, file.path(paths$raw, "synonyms.csv"),
      paths, "synonym list", checksum = checksum_inputs
    ),
    curated_synonyms = application_selected_file(
      params$curated_synonyms_file, file.path(paths$raw, "input_synonyms.csv"),
      paths, "curated synonym list", checksum = checksum_inputs
    ),
    mammal_traits = application_selected_file(
      params$mammal_traits_file, file.path(paths$raw, "mammal_data.txt"),
      paths, "mammal traits", checksum = checksum_inputs
    ),
    bird_traits = application_selected_file(
      params$bird_traits_file, file.path(paths$raw, "bird_data.txt"),
      paths, "bird traits", checksum = checksum_inputs
    ),
    bird_generation_lengths = application_selected_file(
      params$bird_generation_lengths_file,
      file.path(paths$raw, "cobi13486-sup-0004-tables4.xlsx"),
      paths, "bird generation lengths", checksum = checksum_inputs
    ),
    random_effects = application_selected_file(
      params$random_effects_file, file.path(paths$raw, "random_effects.csv"),
      paths, "trait random effects", checksum = checksum_inputs
    ),
    landcover = application_selected_file(
      params$landcover_file, file.path(paths$raw, "esacci_2022_pfts.tif"),
      paths, "land-cover raster", checksum = checksum_inputs
    )
  )
  sdm_input <- application_sdm_input(params, paths, checksum = checksum_inputs)
  context$sdm_input_mode <- sdm_input$mode
  context$sdm_parent_dir <- project_input_path(
    params$sdm_parent_dir, ".", paths$root, "params$sdm_parent_dir"
  )
  context$sdm_index <- sdm_input$index
  context$study_area <- application_study_area(
    params, paths, checksum = checksum_inputs
  )
  context$roi_bounds <- context$study_area$bounds
  context
}

# Prepare manifest/state objects for an active standalone Stage 4 transaction.
# Any disk changes occur only through the called atomic application helpers.
prepare_application_rmd_stage4 <- function(context) {
  validate_persistence_model(context)
  requested <- application_manifest(context, active = TRUE)
  manifest <- application_write_or_validate_manifest(
    context$scenario$application_manifest, requested,
    "application manifest",
    "Use a new application name for different species, SDMs, land cover, or study area."
  )
  list(
    manifest = manifest,
    frozen_inputs = read_application_inputs(context)
  )
}

# Resolve and validate standalone spatial inputs before scientific work begins.
# The returned context contains contracts only, not loaded raster cell values.
validate_application_rmd_spatial_inputs <- function(context, params) {
  app <- readRDS(context$scenario$application_manifest)
  landcover <- application_selected_file(
    params$landcover_file, file.path(context$base_paths$raw, "esacci_2022_pfts.tif"),
    context$base_paths, "land-cover raster"
  )
  assert(isTRUE(landcover$exists), paste0(
    "Missing land-cover raster: ", landcover$path
  ))
  study_area <- application_study_area(params, context$base_paths)
  requested <- list(
    landcover = list(bytes = landcover$bytes, md5 = landcover$md5),
    study_area = application_study_area_contract(
      list(study_area = study_area, inputs = list(landcover = landcover)),
      spatial = identical(study_area$mode, "vector")
    )
  )
  expected <- app$contract[c("landcover", "study_area")]
  differences <- application_manifest_differences(expected, requested)
  assert(!length(differences), paste0(
    "The selected spatial inputs do not match application '",
    context$application, "':\n- ", paste(differences, collapse = "\n- "),
    "\nUse the original inputs or a new application name."
  ))
  invisible(landcover)
}
