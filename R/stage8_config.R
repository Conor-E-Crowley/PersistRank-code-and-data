# Stage 8 application-centered configuration.
#
# Human-readable paths are derived from the few settings people commonly vary:
# application name, quasi-extinction threshold, and Stage 6 removal schedule.
# All remaining scientific settings are recorded in manifests and validated at
# the appropriate resume boundary; none are encoded in opaque directory names.


stage8_curve_order <- function() c("q025", "q16", "q50", "q84", "q975")

# Public execution order. Shared Stage 6 initialization remains an internal
# boundary of stage_6 and is therefore intentionally absent.
stage8_stage_order <- function() c(
  "stage_2", "stage_3", "stage_4", "stage_5",
  "stage_5_1", "stage_5_2", "stage_5_3", "stage_6", "stage_7_1"
)

stage8_core_stage_order <- function() {
  c("stage_2", "stage_3", "stage_4", "stage_5", "stage_6")
}

stage8_selected_stages <- function(params) {
  supplied <- params$stages
  if (is.list(supplied) && identical(names(supplied), "value")) {
    supplied <- supplied$value
  }
  supplied <- as.character(unlist(supplied, use.names = FALSE))
  assert(length(supplied) > 0L, "params$stages must contain at least one stage.")
  assert(
    !anyNA(supplied) && all(nzchar(trimws(supplied))),
    "params$stages must not contain missing or blank values."
  )
  supplied <- trimws(supplied)
  assert(!anyDuplicated(supplied), "params$stages contains duplicate stages.")
  invalid <- setdiff(supplied, stage8_stage_order())
  assert(!length(invalid), paste0(
    "Unsupported Stage 8 stage(s): ", paste(invalid, collapse = ", "), "."
  ))

  core <- stage8_core_stage_order()
  selected_core <- core[core %in% supplied]
  if (length(selected_core) > 1L) {
    positions <- match(selected_core, core)
    assert(
      identical(positions, seq.int(min(positions), max(positions))),
      paste0(
        "Selected core stages must form one contiguous range: ",
        paste(core, collapse = " -> "), "."
      )
    )
  }

  parents <- c(
    stage_5_1 = "stage_5", stage_5_2 = "stage_5",
    stage_5_3 = "stage_5", stage_7_1 = "stage_6"
  )
  for (stage in intersect(names(parents), supplied)) {
    parent <- unname(parents[[stage]])
    earlier <- core[seq_len(match(parent, core) - 1L)]
    assert(
      parent %in% supplied || !any(earlier %in% supplied),
      paste0(
        stage, " cannot reuse omitted ", parent,
        " when an earlier selected core stage may invalidate it. Include ",
        parent, " in params$stages."
      )
    )
  }
  stage8_stage_order()[stage8_stage_order() %in% supplied]
}

stage8_stage_selected <- function(config, stage) {
  stage %in% config$stages
}

stage8_restart_boundary <- function(params, mode, stages) {
  value <- params$restart_from
  if (is.null(value) || !length(value)) value <- NULL
  if (!identical(mode, "restart")) {
    assert(is.null(value), "params$restart_from must be null unless mode is restart.")
    return(NULL)
  }
  value <- validate_scalar_choice(
    value, c("stage_2", "stage_4", "stage_6"), "params$restart_from"
  )
  assert(
    value %in% stages,
    "params$restart_from must also be present in params$stages."
  )
  value
}

stage8_priority_curves <- function(params) {
  supplied <- params$priority_curves
  if (is.list(supplied) && identical(names(supplied), "value")) supplied <- supplied$value
  supplied <- as.character(unlist(supplied, use.names = FALSE))
  assert(length(supplied) > 0L, "params$priority_curves must contain at least one curve.")
  assert(!anyDuplicated(supplied), "params$priority_curves contains duplicate curves.")
  invalid <- setdiff(supplied, stage8_curve_order())
  assert(!length(invalid), paste0(
    "Unsupported priority curve(s): ", paste(invalid, collapse = ", "), "."
  ))
  stage8_curve_order()[stage8_curve_order() %in% supplied]
}

# Validate the public Stage 8 parameters and derive its lightweight execution
# context. Construction may read one existing application manifest to recover
# its authoritative contract. It loads no spatial/native runtime and writes
# nothing; active workflow operations own those effects. Settings and input
# records owned only by omitted stages are NULL rather than synthetic defaults.
stage8_config <- function(params, paths = project_paths()) {
  assert(is.list(params), "Stage 8 parameters must be supplied as a list.")
  allowed <- c(
    "mode", "stages", "restart_from", "stage3_figures", "application",
    "species_list_file", "synonyms_file", "curated_synonyms_file",
    "mammal_traits_file", "bird_traits_file",
    "bird_generation_lengths_file", "random_effects_file",
    "zonation_settings_template",
    "landcover_file", "roi_xmin", "roi_xmax", "roi_ymin", "roi_ymax",
    "study_area_mode", "study_area_file", "study_area_layer",
    "sdm_index_file", "sdm_parent_dir", "process_figure_species",
    "clump_backend", "grass_dir",
    "priority_curves", "taxa", "sdm",
    "persistence_horizon_years", "quasi_extinction_abundance",
    "minimum_patch_abundance", "demographic_uncertainty", "curves",
    "population_cap_factor", "growth_rate_buffer", "n_posterior_draws",
    "replicates_per_draw", "posterior_chunk_size", "base_seed",
    "anchor_probabilities", "persistence_grid_min", "persistence_grid_max",
    "persistence_grid_step", "additional_persistence_probabilities",
    "k_search_start", "k_search_relative_tolerance", "k_search_max",
    "k_rounding", "loess_span", "threads",
    "cells_to_remove_per_iteration", "pruning_iterations_per_stage",
    "max_stages", "checkpoint_every_stages"
  )
  validate_parameter_names(params, allowed, "Stage 8")

  mode <- validate_scalar_choice(params$mode, c("inspect", "resume", "restart"), "params$mode")
  stages <- stage8_selected_stages(params)
  restart_from <- stage8_restart_boundary(params, mode, stages)
  stage3_figures <- if ("stage_3" %in% stages) {
    validate_scalar_logical(params$stage3_figures %||% TRUE, "params$stage3_figures")
  } else NULL
  needs_priority_curves <- any(c("stage_5_3", "stage_6", "stage_7_1") %in% stages)
  curves <- if (needs_priority_curves) stage8_priority_curves(params) else character()
  # q50 is the natural summary curve. For a subset without q50, use its first
  # canonical curve so the driver never needs a second curve selector.
  report_curve <- if (!length(curves)) NULL else if ("q50" %in% curves) {
    "q50"
  } else curves[[1L]]

  identity <- application_run_identity(params, paths)
  application <- identity$application
  cells <- identity$cells_to_remove_per_iteration
  iterations <- identity$pruning_iterations_per_stage
  model <- identity$model
  scenario <- identity$scenario

  application_stages <- c(
    application_scenario_core_stages(),
    application_scenario_optional_stages(), "stage_6", "stage_7_1"
  )
  needs_application <- any(application_stages %in% stages)
  existing_application <- if (needs_application && !"stage_4" %in% stages &&
      file.exists(scenario$application_manifest)) {
    readRDS(scenario$application_manifest)
  } else NULL
  if (!is.null(existing_application)) {
    assert(
      is.list(existing_application) &&
        identical(existing_application$schema, application_storage_schema()) &&
        identical(existing_application$type, "spatial_application") &&
        is.list(existing_application$contract) &&
        all(c("taxa", "sdm", "minimum_patch_abundance") %in%
              names(existing_application$contract)),
      "Unsupported application manifest."
    )
  }
  existing_contract <- existing_application$contract %||% list()
  uses_requested_application <- "stage_4" %in% stages || is.null(existing_application)
  taxa <- if (!uses_requested_application && !is.null(existing_contract$taxa)) {
    validate_taxa_selector(existing_contract$taxa, "application manifest taxa")$value
  } else if (needs_application) {
    validate_taxa_selector(params$taxa, "params$taxa")$value
  } else NULL
  sdm <- if (!uses_requested_application && !is.null(existing_contract$sdm)) {
    validate_priority_sdm_tag(existing_contract$sdm, "application manifest SDM")
  } else if (needs_application) {
    validate_priority_sdm_tag(params$sdm, "params$sdm")
  } else NULL
  minimum_patch <- if (!uses_requested_application &&
      !is.null(existing_contract$minimum_patch_abundance)) {
    existing_contract$minimum_patch_abundance
  } else if (needs_application) params$minimum_patch_abundance else 10L
  contract <- new_analysis_contract(
    identity$persistence_horizon_years,
    identity$quasi_extinction_abundance,
    minimum_patch
  )

  needs_stage4_inputs <- "stage_4" %in% stages
  needs_landcover <- any(c("stage_4", "stage_5", "stage_5_1") %in% stages)
  evaluate_checksums <- !identical(mode, "inspect")

  selected_input <- function(selected, value, default, label) {
    if (!selected) return(NULL)
    application_selected_file(
      value, default, paths, label, checksum = evaluate_checksums
    )
  }
  inputs <- list(
    species_list = selected_input(
      needs_stage4_inputs, params$species_list_file,
      file.path(paths$raw, "simple_summary.csv"), "species list"
    ),
    synonyms = selected_input(
      needs_stage4_inputs, params$synonyms_file,
      file.path(paths$raw, "synonyms.csv"), "synonym list"
    ),
    curated_synonyms = selected_input(
      needs_stage4_inputs, params$curated_synonyms_file,
      file.path(paths$raw, "input_synonyms.csv"), "curated synonym list"
    ),
    mammal_traits = selected_input(
      needs_stage4_inputs, params$mammal_traits_file,
      file.path(paths$raw, "mammal_data.txt"), "mammal traits"
    ),
    bird_traits = selected_input(
      needs_stage4_inputs, params$bird_traits_file,
      file.path(paths$raw, "bird_data.txt"), "bird traits"
    ),
    bird_generation_lengths = selected_input(
      needs_stage4_inputs, params$bird_generation_lengths_file,
      file.path(paths$raw, "cobi13486-sup-0004-tables4.xlsx"),
      "bird generation lengths"
    ),
    random_effects = selected_input(
      needs_stage4_inputs, params$random_effects_file,
      file.path(paths$raw, "random_effects.csv"), "trait random effects"
    ),
    landcover = selected_input(
      needs_landcover, params$landcover_file,
      file.path(paths$raw, "esacci_2022_pfts.tif"), "land-cover raster"
    )
  )
  sdm_input <- if (needs_stage4_inputs) {
    application_sdm_input(params, paths, checksum = evaluate_checksums)
  } else NULL
  study_area <- if (needs_landcover) {
    application_study_area(params, paths, checksum = evaluate_checksums)
  } else NULL
  needs_stage2 <- "stage_2" %in% stages
  k_start <- if (needs_stage2) {
    resolve_stage2_k_search_start(
      params$k_search_start, contract$quasi_extinction_abundance
    )
  } else NULL
  # Persistence models are generic project assets and always contain both
  # taxonomic branches. Stage 8 uses the same configured scientific controls as
  # the standalone Stage 2 adapter; there is no second set of hidden defaults.
  stage2_params <- if (needs_stage2) list(
    mode = if (identical(mode, "inspect")) "inspect" else if (
      identical(restart_from, "stage_2")
    ) "restart" else "resume",
    taxa = "both", verbose = FALSE,
    demographic_uncertainty = params$demographic_uncertainty,
    curves = params$curves,
    persistence_horizon_years = contract$persistence_horizon_years,
    quasi_extinction_abundance = contract$quasi_extinction_abundance,
    population_cap_factor = params$population_cap_factor,
    growth_rate_buffer = params$growth_rate_buffer,
    n_posterior_draws = params$n_posterior_draws,
    replicates_per_draw = params$replicates_per_draw,
    posterior_chunk_size = params$posterior_chunk_size,
    base_seed = params$base_seed,
    anchor_probabilities = params$anchor_probabilities,
    persistence_grid_min = params$persistence_grid_min,
    persistence_grid_max = params$persistence_grid_max,
    persistence_grid_step = params$persistence_grid_step,
    additional_persistence_probabilities = params$additional_persistence_probabilities,
    k_search_start = k_start,
    k_search_relative_tolerance = params$k_search_relative_tolerance,
    k_search_max = params$k_search_max,
    k_rounding = params$k_rounding,
    threads = params$threads
  ) else NULL
  preview <- if (needs_stage2) {
    stage2_config(stage2_params, model$project, model$source_project)
  } else NULL
  needs_stage6 <- "stage_6" %in% stages
  checkpoint_every <- if (needs_stage6) params$checkpoint_every_stages else NULL
  if (needs_stage6 && !is.null(checkpoint_every)) {
    checkpoint_every <- validate_priority_count(
      checkpoint_every, "params$checkpoint_every_stages"
    )
  }
  needs_patch_runtime <- any(c("stage_5", "stage_5_1") %in% stages)
  grass_dir <- if (needs_patch_runtime) {
    value <- validate_patch_grass_param(params$grass_dir)
    if (is.na(value)) NULL else value
  } else NULL

  list(
    schema = application_storage_schema(), mode = mode, stages = stages,
    restart_from = restart_from, stage3_figures = stage3_figures,
    application = application, contract = contract, taxa = taxa, sdm = sdm,
    priority_curves = curves, report_curve = report_curve,
    inputs = inputs,
    sdm_input_mode = sdm_input$mode %||% NULL,
    sdm_index = sdm_input$index %||% NULL,
    sdm_parent_dir = if (needs_stage4_inputs) project_input_path(
      params$sdm_parent_dir, ".", paths$root, "params$sdm_parent_dir"
    ) else NULL,
    zonation_settings_template = if ("stage_5_2" %in% stages) {
      project_input_path(
        params$zonation_settings_template,
        file.path(paths$raw, "settings.z5.txt"), paths$root,
        "params$zonation_settings_template"
      )
    } else NULL,
    study_area = study_area, roi_bounds = study_area$bounds %||% NULL,
    process_figure_species = if (!"stage_5_1" %in% stages ||
      is.null(params$process_figure_species)) NULL else
      validate_scalar_string(params$process_figure_species, "params$process_figure_species"),
    verbose = FALSE,
    loess_span = if ("stage_3" %in% stages) {
      validate_scalar_number(
        params$loess_span, "params$loess_span",
        minimum = 0, maximum = 1, minimum_open = TRUE
      )
    } else NULL,
    clump_backend = if (needs_patch_runtime) {
      validate_patch_backend(params$clump_backend %||% "auto")
    } else NULL,
    grass_dir = grass_dir,
    cells_to_remove_per_iteration = cells,
    pruning_iterations_per_stage = iterations,
    max_stages = if (needs_stage6) {
      validate_max_stages(params$max_stages, "params$max_stages")
    } else NULL,
    ecology_log_every_iterations = if (needs_stage6) 10L else NULL,
    checkpoint_every_stages = checkpoint_every,
    checkpoint_keep = if (needs_stage6) 2L else NULL,
    threads = preview$threads %||% NULL, stage2_params = stage2_params,
    model = model, scenario = scenario, base_paths = identity$base_paths,
    paths = identity$paths,
    existing_application_manifest = existing_application
  )
}

# Validate only the spatial inputs consumed by a selected Stage 5/5.1 when
# Stage 4 is omitted. This avoids rescanning SDM/species inputs while preventing
# a different land cover or study area from entering an existing application.
stage8_validate_spatial_handoff <- function(config, manifest) {
  assert(
    is.list(manifest) && identical(manifest$schema, application_storage_schema()) &&
      identical(manifest$type, "spatial_application"),
    "Unsupported application manifest."
  )
  found_landcover <- manifest$contract$landcover
  requested_landcover <- config$inputs$landcover
  assert(
    isTRUE(requested_landcover$exists) &&
      identical(requested_landcover$bytes, found_landcover$bytes) &&
      identical(requested_landcover$md5, found_landcover$md5),
    "Selected Stage 5 spatial work requires the land cover recorded by the application manifest."
  )
  requested_area <- application_study_area_contract(config, spatial = FALSE)
  found_area <- manifest$contract$study_area
  fields <- intersect(
    c("mode", "bounds", "file_md5", "layer"),
    union(names(requested_area), names(found_area))
  )
  differences <- application_manifest_differences(
    found_area[fields], requested_area[fields]
  )
  assert(!length(differences), paste0(
    "Selected Stage 5 spatial work has a different study area:\n- ",
    paste(differences, collapse = "\n- ")
  ))
  invisible(TRUE)
}

# Stage-specific translations used by Stage 8. These return public-shaped
# parameter lists only; the application-aware stage adapters perform the final
# low-level config construction shared with standalone Rmds.
stage8_stage4_params <- function(config, frozen_inputs = NULL) list(
  mode = "build", verbose = config$verbose, taxa = config$taxa,
  sdm = config$sdm, sdm_index_file = config$sdm_index$path %||% NULL,
  iucn_mode = if (is.null(frozen_inputs)) "cache_or_query" else "cache_only",
  iucn_pause_seconds = 2, sdm_parent_dir = config$sdm_parent_dir,
  species_list_file = config$inputs$species_list$path,
  synonyms_file = config$inputs$synonyms$path,
  curated_synonyms_file = config$inputs$curated_synonyms$path,
  mammal_traits_file = config$inputs$mammal_traits$path,
  bird_traits_file = config$inputs$bird_traits$path,
  bird_generation_lengths_file = config$inputs$bird_generation_lengths$path,
  random_effects_file = config$inputs$random_effects$path
)

stage8_stage5_params <- function(config, mode = "resume") list(
  mode = mode, taxa = config$taxa, verbose = config$verbose,
  clump_backend = config$clump_backend, grass_dir = config$grass_dir,
  landcover_file = config$inputs$landcover$path,
  study_area_mode = config$study_area$mode,
  study_area_file = config$study_area$file,
  study_area_layer = config$study_area$layer,
  roi_xmin = config$roi_bounds[["xmin"]] %||% 0,
  roi_xmax = config$roi_bounds[["xmax"]] %||% 1,
  roi_ymin = config$roi_bounds[["ymin"]] %||% 0,
  roi_ymax = config$roi_bounds[["ymax"]] %||% 1
)

stage8_stage51_params <- function(config) c(
  stage8_stage5_params(config),
  list(process_figure_species = config$process_figure_species)
)

stage8_stage52_params <- function(config, mode = "inputs") list(
  mode = mode,
  zonation_settings_template = config$zonation_settings_template
)

stage8_stage53_params <- function(config) list(
  taxa = config$taxa, curve = config$report_curve
)

stage8_stage6_params <- function(config, mode, curve = NULL) {
  if (identical(mode, "initialize")) return(list(mode = mode))
  list(
    mode = mode, curve = curve, taxa = config$taxa, sdm = config$sdm,
    cells_to_remove_per_iteration = config$cells_to_remove_per_iteration,
    pruning_iterations_per_stage = config$pruning_iterations_per_stage,
    max_stages = config$max_stages,
    ecology_log_every_iterations = config$ecology_log_every_iterations,
    checkpoint_every_stages = config$checkpoint_every_stages,
    checkpoint_keep = config$checkpoint_keep, resume_stage = NULL
  )
}

stage8_stage71_params <- function(config, curve) list(
  curve = curve, taxa = config$taxa, sdm = config$sdm,
  cells_to_remove_per_iteration = config$cells_to_remove_per_iteration,
  pruning_iterations_per_stage = config$pruning_iterations_per_stage
)

describe_stage8_config <- function(config) {
  curves <- if (length(config$priority_curves)) {
    paste(config$priority_curves, collapse = ",")
  } else "not selected"
  log_msg(
    "STAGE 8 | mode=", config$mode,
    " stages=", paste(config$stages, collapse = ","),
    " restart_from=", config$restart_from %||% "none",
    " application=", config$application,
    " threshold=", config$contract$quasi_extinction_abundance,
    " model=", config$model$root,
    " scenario=", config$scenario$root,
    " run=", config$scenario$run_root,
    " curves=", curves
  )
  invisible(config)
}

stage8_format_setting <- function(value) {
  if (is.null(value) || !length(value)) return("automatic")
  if (is.list(value) && identical(names(value), "value")) value <- value$value
  if (is.logical(value)) return(paste(tolower(as.character(value)), collapse = ", "))
  paste(as.character(value), collapse = ", ")
}

stage8_settings_table <- function(config) {
  values <- list(
    mode = config$mode,
    stages = config$stages,
    restart_from = config$restart_from,
    application = config$application,
    quasi_extinction_abundance = config$contract$quasi_extinction_abundance,
    persistence_horizon_years = config$contract$persistence_horizon_years,
    cells_to_remove_per_iteration = config$cells_to_remove_per_iteration,
    pruning_iterations_per_stage = config$pruning_iterations_per_stage
  )
  if (any(c(
    application_scenario_core_stages(), application_scenario_optional_stages(),
    "stage_6", "stage_7_1"
  ) %in% config$stages)) values <- c(values, list(
    minimum_patch_abundance = config$contract$minimum_patch_abundance,
    taxa = config$taxa, sdm = config$sdm
  ))
  if (stage8_stage_selected(config, "stage_4")) values <- c(values, list(
    sdm_input = if (identical(config$sdm_input_mode, "index")) {
      "index file"
    } else "standard folders"
  ))
  if (any(c("stage_4", "stage_5", "stage_5_1") %in% config$stages)) {
    values <- c(values, list(study_area_mode = config$study_area$mode))
  }
  if (any(c("stage_5", "stage_5_1") %in% config$stages)) values <- c(
    values, list(
      clump_backend = config$clump_backend, grass_dir = config$grass_dir
    )
  )
  if (stage8_stage_selected(config, "stage_3")) {
    values <- c(values, list(stage3_figures = config$stage3_figures))
  }
  if (length(config$priority_curves)) values <- c(values, list(
    priority_curves = config$priority_curves,
    report_curve = config$report_curve
  ))
  if (stage8_stage_selected(config, "stage_2")) values <- c(values, list(
    k_search_start = config$stage2_params$k_search_start,
    threads = config$threads
  ))
  if (stage8_stage_selected(config, "stage_6")) values <- c(values, list(
    max_stages = config$max_stages,
    checkpoint_every_stages = config$checkpoint_every_stages
  ))
  data.frame(
    setting = names(values),
    value = vapply(values, stage8_format_setting, character(1L)),
    stringsAsFactors = FALSE, row.names = NULL
  )
}
