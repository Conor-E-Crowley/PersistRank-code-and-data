# Shared context and manifest contracts for named spatial applications.
#
# Standalone Rmd adapters and Stage 8 load this definition-only module after
# path, storage, analysis, species, and patch contracts. Reading and validating
# context is lightweight; spatial libraries are loaded only when a caller asks
# for canonical vector geometry during an active application build.

# Scientific handoffs and optional presentation products deliberately have
# separate contracts. Stage 6 consumes only Stage 4--5 science; benchmark and
# figure products may be produced later without changing scenario identity.
application_scenario_core_stages <- function() c("stage_4", "stage_5")

application_scenario_optional_stages <- function() {
  c("stage_5_1", "stage_5_2", "stage_5_3")
}

application_scenario_shared_stages <- function() "stage_6_initialization"

application_curve_core_stages <- function() "stage_6"

application_curve_optional_stages <- function() "stage_7_1"

application_boundary_stages <- function() c(
  application_scenario_core_stages(),
  application_scenario_optional_stages(),
  application_scenario_shared_stages(),
  application_curve_core_stages(),
  application_curve_optional_stages()
)

# Resolve the small application/run identity shared by standalone stages,
# Stage 8, and completed-run reports. This performs validation and path
# derivation only: it reads no manifests, hashes no files, and writes nothing.
application_run_identity <- function(params, paths = project_paths()) {
  assert(is.list(params), "Application parameters must be supplied as a list.")
  application <- validate_application_name(
    params$application %||% "madagascar", "params$application"
  )
  horizon <- validate_scalar_integer(
    params$persistence_horizon_years %||% 1000L,
    "params$persistence_horizon_years"
  )
  threshold <- validate_scalar_integer(
    params$quasi_extinction_abundance %||% 500L,
    "params$quasi_extinction_abundance"
  )
  cells <- validate_priority_count(
    params$cells_to_remove_per_iteration %||% 1000L,
    "params$cells_to_remove_per_iteration"
  )
  iterations <- validate_priority_count(
    params$pruning_iterations_per_stage %||% 50L,
    "params$pruning_iterations_per_stage"
  )
  model <- persistence_model_paths(paths, threshold, horizon)
  scenario <- application_scenario_paths(
    paths, application, threshold, horizon, cells, iterations
  )
  list(
    application = application,
    persistence_horizon_years = horizon,
    quasi_extinction_abundance = threshold,
    cells_to_remove_per_iteration = cells,
    pruning_iterations_per_stage = iterations,
    model = model,
    scenario = scenario,
    base_paths = paths,
    paths = scenario$project
  )
}

# Read and validate the four manifests defining one completed application run.
# Each manifest is read once. Manifest-owned taxa, SDM, and minimum patch
# abundance are authoritative; the function performs no writes or spatial work.
read_application_run_context <- function(identity) {
  assert(is.list(identity) && !is.null(identity$scenario) && !is.null(identity$model),
         "identity must be an application_run_identity() result.")
  required <- c(
    application = identity$scenario$application_manifest,
    scenario = identity$scenario$manifest,
    run = identity$scenario$run_manifest,
    persistence_model = identity$model$manifest
  )
  for (label in names(required)) need_file(required[[label]], paste(label, "manifest"))
  manifests <- lapply(required, readRDS)
  application <- manifests$application
  scenario <- manifests$scenario
  run <- manifests$run
  model <- manifests$persistence_model

  assert(identical(application$schema, application_storage_schema()) &&
           identical(application$type, "spatial_application"),
         "Unsupported application manifest.")
  assert(identical(scenario$schema, application_storage_schema()) &&
           identical(scenario$type, "threshold_scenario") &&
           identical(scenario$contract$threshold_tag, identity$scenario$threshold_tag),
         "Incompatible threshold-scenario manifest.")
  expected_schedule <- list(
    cells_to_remove_per_iteration = identity$cells_to_remove_per_iteration,
    pruning_iterations_per_stage = identity$pruning_iterations_per_stage
  )
  assert(identical(run$schema, application_storage_schema()) &&
           identical(run$type, "priority_run") &&
           identical(run$contract$schedule, expected_schedule),
         "Incompatible priority-run manifest or removal schedule.")
  assert(identical(model$schema, application_storage_schema()) &&
           identical(model$type, "persistence_model") &&
           identical(model$contract$taxa, "both") &&
           identical(
             model$contract$persistence_horizon_years,
             identity$persistence_horizon_years
           ) &&
           identical(
             model$contract$quasi_extinction_abundance,
             identity$quasi_extinction_abundance
           ) && identical(model$contract$curves, persistence_curves()),
         "Incompatible persistence-model manifest.")
  completed <- if (!is.list(run$curves)) character() else names(run$curves)[vapply(
    run$curves, function(curve) identical(curve$status, "complete"), logical(1L)
  )]
  if (length(completed)) {
    initialization <- scenario$stages[[application_scenario_shared_stages()]]
    assert(
      is.list(initialization) && length(initialization$outputs) == 1L,
      paste0(
        "Completed Stage 6 curves require the shared Stage 6 initialization ",
        "record. Reinitialize Stage 6 and regenerate the selected runs."
      )
    )
    initialization_record <- initialization$outputs[[1L]]
    expected_path <- storage_relative_path(
      identity$scenario$priority_initialization,
      identity$base_paths$root
    )
    for (curve in completed) {
      entry <- run$curves[[curve]]
      coefficient <- entry$coefficient_identity
      shared <- entry$initialization
      assert(
        is.list(shared) &&
          identical(shared$path, expected_path) &&
          identical(shared$md5, initialization_record$md5) &&
          is.list(coefficient) &&
          identical(coefficient$curve, curve) &&
          identical(coefficient$alpha_column, paste0("alpha_", curve)) &&
          identical(coefficient$beta_column, paste0("beta_", curve)) &&
          is.character(coefficient$md5) && length(coefficient$md5) == 1L &&
          !is.na(coefficient$md5) && grepl("^[[:xdigit:]]{32}$", coefficient$md5),
        paste0(
          "Completed curve ", curve,
          " does not identify the current shared Stage 6 initialization and ",
          "coefficient pair. Regenerate that Stage 6 run."
        )
      )
    }
  }
  completed <- persistence_curves()[persistence_curves() %in% completed]
  app_contract <- application$contract
  list(
    application = identity$application,
    persistence_horizon_years = identity$persistence_horizon_years,
    quasi_extinction_abundance = identity$quasi_extinction_abundance,
    cells_to_remove_per_iteration = identity$cells_to_remove_per_iteration,
    pruning_iterations_per_stage = identity$pruning_iterations_per_stage,
    taxa = app_contract$taxa,
    sdm = app_contract$sdm,
    contract = new_analysis_contract(
      identity$persistence_horizon_years,
      identity$quasi_extinction_abundance,
      app_contract$minimum_patch_abundance
    ),
    completed_curves = completed,
    project_paths = identity$paths,
    scenario = identity$scenario,
    model = identity$model,
    manifests = manifests
  )
}

# Resolve a repository-relative input to a normalized lightweight fingerprint.
# Missing files return an explicit unavailable record so inspect mode stays usable.
application_selected_file <- function(value, default, paths, label,
                                      checksum = TRUE) {
  checksum <- validate_scalar_logical(checksum, "checksum")
  value <- project_input_path(value, default, paths$root, label)
  if (!file.exists(value) || dir.exists(value)) {
    return(list(
      path = value, bytes = NA_real_, md5 = NA_character_, exists = FALSE,
      checksum_status = "unavailable"
    ))
  }
  if (!isTRUE(checksum)) {
    return(list(
      path = value, bytes = unname(as.numeric(file.info(value)$size)),
      md5 = NA_character_, exists = TRUE, checksum_status = "not_evaluated"
    ))
  }
  c(
    storage_file_fingerprint(value, label),
    list(exists = TRUE, checksum_status = "verified")
  )
}

# Select indexed or folder-based SDM discovery without opening spatial data.
application_sdm_input <- function(params, paths, checksum = TRUE) {
  value <- params$sdm_index_file
  if (is.null(value) || !length(value)) {
    return(list(mode = "folders", index = NULL))
  }
  list(
    mode = "index",
    index = application_selected_file(
      value, NULL, paths, "SDM index CSV", checksum = checksum
    )
  )
}

# Normalize bounds/vector study-area parameters. Vector mode fingerprints the
# source file but deliberately defers spatial geometry work.
application_study_area <- function(params, paths, checksum = TRUE) {
  area <- patch_study_area_from_params(params, paths$root)
  if (identical(area$mode, "vector")) {
    area$fingerprint <- if (!isTRUE(checksum)) {
      list(path = area$file, bytes = if (file.exists(area$file)) {
        unname(as.numeric(file.info(area$file)$size))
      } else NA_real_, md5 = NA_character_)
    } else if (file.exists(area$file) && !dir.exists(area$file)) {
      storage_file_fingerprint(area$file, "study-area vector")
    } else {
      list(path = area$file, bytes = NA_real_, md5 = NA_character_)
    }
  }
  area
}

# Build the deterministic SDM inventory used by active application manifests.
# It reads metadata and hashes rasters but never loads raster cells into memory.
application_sdm_inventory <- function(config) {
  if (identical(config$sdm_input_mode, "index")) {
    index <- data.table::fread(config$sdm_index$path)
    need_cols(
      index,
      c("scientificName", "taxon_class", "sdm_method", "raster_path"),
      "SDM index CSV"
    )
    index$raster_path <- resolve_sdm_index_raster_paths(
      index$raster_path, config$sdm_index$path
    )
    missing <- index$raster_path[!file.exists(index$raster_path)]
    assert(!length(missing), paste0(
      "SDM index references missing raster(s):\n", paste(missing, collapse = "\n")
    ))
    files <- data.frame(
      scientificName = as.character(index$scientificName),
      taxon_class = as.character(index$taxon_class),
      sdm_method = as.character(index$sdm_method),
      bytes = as.numeric(file.info(index$raster_path)$size),
      md5 = unname(tools::md5sum(index$raster_path)),
      stringsAsFactors = FALSE
    )
    files <- files[do.call(
      order, files[c("taxon_class", "sdm_method", "scientificName")]
    ), ]
    rownames(files) <- NULL
    return(list(mode = "index", files = files))
  }

  folders <- selected_sdm_directories(
    config$sdm_parent_dir, config$taxa, config$sdm
  )
  roots <- folders$dir_path
  missing <- roots[!dir.exists(roots)]
  assert(!length(missing), paste0(
    "Selected SDM folder(s) are missing:\n", paste(missing, collapse = "\n")
  ))
  raster_files <- lapply(
    roots, list.files, pattern = "_bin[.]tif$", full.names = TRUE
  )
  empty <- roots[lengths(raster_files) == 0L]
  assert(!length(empty), paste0(
    "Selected SDM folder(s) contain no *_bin.tif files:\n",
    paste(empty, collapse = "\n")
  ))
  file_rows <- lapply(seq_along(roots), function(i) {
    paths <- sort(raster_files[[i]])
    data.frame(
      taxon_class = folders$taxon_class[[i]],
      sdm_method = folders$sdm_method[[i]],
      file = basename(paths),
      bytes = as.numeric(file.info(paths)$size),
      md5 = unname(tools::md5sum(paths)),
      stringsAsFactors = FALSE
    )
  })
  files <- do.call(rbind, file_rows)
  assert(nrow(files) > 0L, "Selected SDM folders contain no *_bin.tif files.")
  list(mode = "folders", files = files)
}

# Return the stable study-area contract. `spatial = TRUE` additionally computes
# canonical geometry identity and is reserved for active Stage 4 validation.
application_study_area_contract <- function(config, spatial = FALSE) {
  area <- config$study_area
  out <- list(mode = area$mode)
  if (identical(area$mode, "bounds")) {
    out$bounds <- unname(as.numeric(area$bounds))
  }
  if (identical(area$mode, "vector")) {
    out$file_md5 <- area$fingerprint$md5
    out$layer <- area$layer
    if (isTRUE(spatial)) {
      load_packages(c("sf", "terra"))
      vector <- sf::st_read(area$file, layer = area$layer, quiet = TRUE)
      assert(nrow(vector) > 0L, "Study-area vector is empty.")
      landcover <- terra::rast(config$inputs$landcover$path)
      vector <- sf::st_transform(
        sf::st_make_valid(sf::st_union(vector)),
        terra::crs(landcover, proj = TRUE)
      )
      wkb <- sf::st_as_binary(sf::st_geometry(vector), EWKB = TRUE)
      temporary <- tempfile("study_area_", fileext = ".rds")
      on.exit(unlink(temporary, force = TRUE), add = TRUE)
      saveRDS(wkb, temporary, version = 3L)
      out$canonical_geometry_md5 <- unname(tools::md5sum(temporary))
    }
  }
  out
}

# Construct the immutable application input contract. Active mode verifies
# complete fingerprints and inventories; inspect mode remains metadata-only.
application_contract <- function(config, active = FALSE) {
  if (isTRUE(active)) {
    landcover <- config$inputs$landcover
    assert(
      isTRUE(landcover$exists) && is.finite(landcover$bytes) &&
        landcover$bytes > 0 && is.character(landcover$md5) &&
        length(landcover$md5) == 1L && !is.na(landcover$md5) &&
        grepl("^[[:xdigit:]]{32}$", landcover$md5),
      "Active application validation requires a verified land-cover fingerprint."
    )
  }
  list(
    schema = application_storage_schema(),
    input_files = lapply(
      config$inputs[c(
        "species_list", "synonyms", "curated_synonyms", "mammal_traits",
        "bird_traits", "bird_generation_lengths", "random_effects"
      )],
      function(x) list(bytes = x$bytes, md5 = x$md5)
    ),
    taxa = config$taxa,
    sdm = config$sdm,
    sdm_inventory = if (isTRUE(active)) {
      application_sdm_inventory(config)
    } else {
      list(
        mode = config$sdm_input_mode,
        index_md5 = config$sdm_index$md5 %||% NULL
      )
    },
    landcover = list(
      bytes = config$inputs$landcover$bytes,
      md5 = config$inputs$landcover$md5
    ),
    study_area = application_study_area_contract(config, spatial = active),
    minimum_patch_abundance = as.integer(config$contract$minimum_patch_abundance),
    clump_contract = "strict_patch_and_population_unit_area_thresholds_v1"
  )
}

# Wrap the application contract in the unchanged durable manifest schema.
application_manifest <- function(config, active = TRUE) {
  list(
    schema = application_storage_schema(), type = "spatial_application",
    application = config$application,
    contract = application_contract(config, active),
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
}

# Construct the threshold-scenario handoff from completed Stage 4–5 records.
# The result is in memory; the caller owns its atomic write transaction.
application_scenario_manifest <- function(config, state, model_manifest,
                                          application_manifest) {
  list(
    schema = application_storage_schema(), type = "threshold_scenario",
    application = config$application,
    contract = list(
      threshold_tag = config$scenario$threshold_tag,
      persistence_model_manifest_md5 = unname(tools::md5sum(config$model$manifest)),
      application_manifest_md5 = unname(tools::md5sum(
        config$scenario$application_manifest
      )),
      abundance_contract = unclass(config$contract)
    ),
    stages = state$stages[intersect(
      c(application_scenario_core_stages(),
        application_scenario_optional_stages(),
        application_scenario_shared_stages()),
      names(state$stages)
    )],
    persistence_model = storage_relative_path(
      config$model$root, config$base_paths$root
    ),
    application_manifest = application_manifest,
    model_manifest = model_manifest
  )
}

# Validate the persisted Stage 4–5 handoff, its model/application identities,
# abundance contract, and every recorded artifact. This function is read-only.
application_validate_scenario_handoff <- function(
  config, required_stages = application_scenario_core_stages()
) {
  required_stages <- unique(as.character(required_stages))
  allowed <- c(
    application_scenario_core_stages(),
    application_scenario_optional_stages()
  )
  assert(
    length(required_stages) > 0L && !anyNA(required_stages) &&
      all(required_stages %in% allowed),
    "required_stages must contain known scenario-stage names."
  )
  need_file(config$scenario$manifest, "threshold-scenario manifest")
  manifest <- readRDS(config$scenario$manifest)
  expected <- list(
    threshold_tag = config$scenario$threshold_tag,
    persistence_model_manifest_md5 = unname(tools::md5sum(config$model$manifest)),
    application_manifest_md5 = unname(tools::md5sum(
      config$scenario$application_manifest
    )),
    abundance_contract = unclass(config$contract)
  )
  differences <- application_manifest_differences(manifest$contract, expected)
  assert(!length(differences), paste0(
    "The reused Stage 4–5 products are incompatible:\n- ",
    paste(differences, collapse = "\n- "),
    "\nRestart from stage_4 to rebuild this scenario."
  ))
  missing_stages <- required_stages[!required_stages %in% names(manifest$stages)]
  present <- required_stages[required_stages %in% names(manifest$stages)]
  invalid_stages <- present[!vapply(
    manifest$stages[present],
    application_record_valid,
    logical(1L),
    root = config$base_paths$root
  )]
  assert(!length(missing_stages) && !length(invalid_stages), paste0(
    "The reused Stage 4–5 handoff is incomplete or has changed artifacts.",
    if (length(missing_stages)) {
      paste0(" Missing: ", paste(missing_stages, collapse = ", "), ".")
    } else "",
    if (length(invalid_stages)) {
      paste0(" Invalid: ", paste(invalid_stages, collapse = ", "), ".")
    } else "",
    " Restart from stage_4 to rebuild this scenario."
  ))
  manifest
}

# Return the immutable Stage 6 schedule/algorithm identity for a run manifest.
application_run_contract <- function(config) list(
  schedule = list(
    cells_to_remove_per_iteration = config$cells_to_remove_per_iteration,
    pruning_iterations_per_stage = config$pruning_iterations_per_stage
  ),
  algorithm_schema = priority_initialization_schema_version()
)

# Read frozen Stage 4 habitat rows, returning NULL when none exist. Invalid
# schema/content fails without changing the file.
read_application_inputs <- function(config) {
  if (!file.exists(config$scenario$application_inputs)) return(NULL)
  saved <- readRDS(config$scenario$application_inputs)
  assert(
    identical(saved$schema, application_storage_schema()) &&
      is.data.frame(saved$habitat_rows),
    "Unsupported application_inputs.rds file."
  )
  saved
}

# Atomically freeze the Stage 4 habitat rows used by reproducible reruns.
# Returns the installed object invisibly; it does not retain spatial objects.
write_application_inputs <- function(config, habitat_rows) {
  object <- list(
    schema = application_storage_schema(),
    habitat_rows = as.data.frame(habitat_rows)
  )
  storage_atomic_save_rds(
    object, config$scenario$application_inputs, "frozen application IUCN inputs"
  )
  object
}
