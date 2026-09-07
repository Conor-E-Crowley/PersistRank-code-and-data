# Human-readable project configuration.
#
# A configuration YAML owns scientific settings and input paths. Numbered Rmds
# supply only operation selectors (for example mode, curve, or benchmark). This
# module performs no scientific reads and writes nothing: it parses one small
# YAML file, rejects ambiguous keys, and returns stage-shaped parameter lists.

project_config_section_keys <- function() {
  list(
    demography = c(
      "bird_sigma_separate_intercepts", "seed", "n_chains", "n_adapt",
      "n_iter", "thin", "coefficient_prior_mean", "coefficient_prior_sd",
      "residual_sd_min", "residual_sd_max", "rhat_max", "ess_min",
      "trait_grid_points", "mammal_mass_min_g", "mammal_mass_max_g"
    ),
    persistence = c(
      "demographic_uncertainty", "curves", "persistence_horizon_years",
      "quasi_extinction_abundance", "population_cap_factor",
      "growth_rate_buffer", "n_posterior_draws", "replicates_per_draw",
      "posterior_chunk_size", "base_seed", "anchor_probabilities",
      "persistence_grid_min", "persistence_grid_max",
      "persistence_grid_step", "additional_persistence_probabilities",
      "k_search_start", "k_search_relative_tolerance", "k_search_max",
      "k_rounding", "loess_span"
    ),
    application = c(
      "application", "minimum_patch_abundance", "taxa", "sdm"
    ),
    inputs = c(
      "mammal_rmax_file", "sigma_file", "bird_growth_file",
      "mammal_traits_file", "bird_traits_file",
      "bird_generation_lengths_file", "iucn_bird_synonyms_file",
      "species_list_file", "synonyms_file", "curated_synonyms_file",
      "random_effects_file", "landcover_file", "sdm_parent_dir",
      "sdm_index_file", "zonation_settings_template"
    ),
    spatial = c(
      "study_area_mode", "study_area_file", "study_area_layer",
      "roi_xmin", "roi_xmax", "roi_ymin", "roi_ymax"
    ),
    priority = c(
      "cells_to_remove_per_iteration", "pruning_iterations_per_stage"
    ),
    reporting = c(
      "late_stage_removal_percentages", "persistence_threshold",
      "focal_species", "assemblage_point_style", "focal_point_style",
      "assemblage_central_stat"
    )
  )
}

project_config_stage_keys <- function(stage) {
  shared_identity <- c(
    "application", "quasi_extinction_abundance",
    "persistence_horizon_years", "cells_to_remove_per_iteration",
    "pruning_iterations_per_stage"
  )
  spatial <- c(
    "landcover_file", "study_area_mode", "study_area_file",
    "study_area_layer", "roi_xmin", "roi_xmax", "roi_ymin", "roi_ymax"
  )
  report <- c(
    "persistence_threshold", "focal_species", "assemblage_point_style",
    "focal_point_style", "assemblage_central_stat"
  )
  keys <- list(
    stage1 = c(
      "bird_sigma_separate_intercepts", "seed", "n_chains", "n_adapt",
      "n_iter", "thin", "coefficient_prior_mean", "coefficient_prior_sd",
      "residual_sd_min", "residual_sd_max", "rhat_max", "ess_min",
      "trait_grid_points", "mammal_mass_min_g", "mammal_mass_max_g",
      "mammal_rmax_file", "sigma_file", "bird_traits_file",
      "bird_generation_lengths_file", "bird_growth_file",
      "iucn_bird_synonyms_file", "curated_synonyms_file"
    ),
    stage2 = c(
      "demographic_uncertainty", "curves", "persistence_horizon_years",
      "quasi_extinction_abundance", "population_cap_factor",
      "growth_rate_buffer", "n_posterior_draws", "replicates_per_draw",
      "posterior_chunk_size", "base_seed", "anchor_probabilities",
      "persistence_grid_min", "persistence_grid_max",
      "persistence_grid_step", "additional_persistence_probabilities",
      "k_search_start", "k_search_relative_tolerance", "k_search_max",
      "k_rounding"
    ),
    stage3 = c(
      "loess_span", "persistence_horizon_years",
      "quasi_extinction_abundance"
    ),
    stage4 = c(
      shared_identity, "minimum_patch_abundance", "taxa", "sdm",
      "sdm_index_file", "sdm_parent_dir", "species_list_file",
      "synonyms_file", "curated_synonyms_file", "mammal_traits_file",
      "bird_traits_file", "bird_generation_lengths_file",
      "random_effects_file", spatial
    ),
    stage5 = c(shared_identity, spatial),
    stage51 = c(shared_identity, spatial),
    stage52 = c(shared_identity, "zonation_settings_template"),
    stage53 = shared_identity,
    stage6 = shared_identity,
    stage71 = shared_identity,
    stage72 = shared_identity,
    stage73 = c(shared_identity, report),
    stage74 = c(shared_identity, "late_stage_removal_percentages", report[1L]),
    stage8 = c(
      shared_identity, "minimum_patch_abundance", "taxa", "sdm",
      "sdm_index_file", "sdm_parent_dir", "species_list_file",
      "synonyms_file", "curated_synonyms_file", "mammal_traits_file",
      "bird_traits_file", "bird_generation_lengths_file",
      "random_effects_file", "zonation_settings_template", spatial,
      "demographic_uncertainty", "curves", "population_cap_factor",
      "growth_rate_buffer", "n_posterior_draws", "replicates_per_draw",
      "posterior_chunk_size", "base_seed", "anchor_probabilities",
      "persistence_grid_min", "persistence_grid_max",
      "persistence_grid_step", "additional_persistence_probabilities",
      "k_search_start", "k_search_relative_tolerance", "k_search_max",
      "k_rounding", "loess_span"
    ),
    stage81 = c(
      shared_identity, "late_stage_removal_percentages", report
    )
  )
  assert(stage %in% names(keys), paste0("Unknown project configuration stage: ", stage))
  unique(keys[[stage]])
}

project_config_operational_keys <- function(stage) {
  keys <- list(
    stage1 = c("mode", "verbose"),
    stage2 = c("mode", "verbose", "threads"),
    stage3 = c("mode", "write_figures"),
    stage4 = c("mode", "verbose", "iucn_mode", "iucn_pause_seconds"),
    stage5 = c("mode", "verbose", "clump_backend", "grass_dir"),
    stage51 = c("process_figure_species", "clump_backend", "grass_dir"),
    stage52 = "mode",
    stage53 = "curve",
    stage6 = c(
      "mode", "curve", "max_stages", "ecology_log_every_iterations",
      "checkpoint_every_stages", "checkpoint_keep", "resume_stage"
    ),
    stage71 = "curve",
    stage72 = c("mode", "optimization_curves", "rank_method"),
    stage73 = c("mode", "main_curve", "rank_method"),
    stage74 = c("mode", "optimization_curves", "main_curve", "rank_method"),
    stage8 = c(
      "mode", "stages", "restart_from", "stage3_figures",
      "process_figure_species", "clump_backend", "grass_dir",
      "priority_curves", "threads", "max_stages", "checkpoint_every_stages"
    ),
    stage81 = c("mode", "optimization_curves", "main_curve")
  )
  assert(stage %in% names(keys), paste0("Unknown project configuration stage: ", stage))
  keys[[stage]]
}

read_project_config <- function(path = "config/madagascar.yml",
                                root = project_paths()$root) {
  assert(requireNamespace("yaml", quietly = TRUE), paste0(
    "Package 'yaml' is required to read the project configuration."
  ))
  path <- project_input_path(path, root = root, label = "config_file")
  need_file(path, "project configuration")
  document <- yaml::read_yaml(path)
  assert(is.list(document), "Project configuration must be a YAML mapping.")
  assert(
    identical(as.integer(document$schema_version), 1L),
    "Project configuration schema_version must be 1."
  )
  section_keys <- project_config_section_keys()
  sections <- names(section_keys)
  unknown <- setdiff(names(document), c("schema_version", sections))
  missing <- setdiff(sections, names(document))
  assert(!length(unknown), paste0(
    "Unknown project configuration section(s): ", paste(unknown, collapse = ", ")
  ))
  assert(!length(missing), paste0(
    "Missing project configuration section(s): ", paste(missing, collapse = ", ")
  ))

  values <- list()
  for (section in sections) {
    entries <- document[[section]]
    assert(
      is.list(entries) && length(entries) > 0L &&
        !is.null(names(entries)) && all(nzchar(names(entries))),
      paste0("Configuration section '", section, "' must be a named mapping.")
    )
    unexpected <- setdiff(names(entries), section_keys[[section]])
    absent <- setdiff(section_keys[[section]], names(entries))
    assert(!length(unexpected), paste0(
      "Unknown key(s) in configuration section '", section, "': ",
      paste(unexpected, collapse = ", ")
    ))
    assert(!length(absent), paste0(
      "Missing key(s) from configuration section '", section, "': ",
      paste(absent, collapse = ", ")
    ))
    for (name in names(entries)) values[name] <- list(entries[[name]])
  }
  structure(
    list(path = path, root = normalizePath(root, winslash = "/", mustWork = FALSE),
         values = values),
    class = "project_configuration"
  )
}

# Cheap, read-only input diagnostics for a copied application configuration.
# This resolves paths and counts filenames only; it never opens tabular or
# spatial data, calculates checksums, contacts IUCN, or loads native runtimes.
inspect_project_config <- function(configuration) {
  assert(
    inherits(configuration, "project_configuration"),
    "configuration must come from read_project_config()."
  )
  input_keys <- c(
    "mammal_rmax_file", "sigma_file", "bird_growth_file",
    "mammal_traits_file", "bird_traits_file",
    "bird_generation_lengths_file", "iucn_bird_synonyms_file",
    "species_list_file", "synonyms_file", "curated_synonyms_file",
    "random_effects_file", "landcover_file", "sdm_parent_dir",
    "sdm_index_file", "zonation_settings_template", "study_area_file"
  )
  directory_keys <- "sdm_parent_dir"
  optional_keys <- c("sdm_index_file", "study_area_file")
  rows <- lapply(input_keys, function(key) {
    value <- configuration$values[[key]]
    configured <- !is.null(value) && length(value) > 0L
    path <- if (configured) {
      project_input_path(value, root = configuration$root, label = key)
    } else NA_character_
    is_directory <- key %in% directory_keys
    exists <- if (!configured) NA else if (is_directory) dir.exists(path) else
      file.exists(path) && !dir.exists(path)
    data.frame(
      key = key, required = !key %in% optional_keys,
      kind = if (is_directory) "directory" else "file",
      configured = configured, exists = exists, path = path,
      stringsAsFactors = FALSE
    )
  })
  inputs <- do.call(rbind, rows)
  rownames(inputs) <- NULL

  parent <- inputs$path[inputs$key == "sdm_parent_dir"]
  definitions <- data.frame(
    taxon = c("mammals", "mammals", "birds", "birds"),
    method = c("PPM", "RangeBag", "PPM", "RangeBag"),
    directory = c(
      "mammal_ppm_bin", "mammal_rangebag_bin",
      "bird_ppm_bin", "bird_rangebag_bin"
    ),
    stringsAsFactors = FALSE
  )
  definitions$path <- file.path(parent, definitions$directory)
  definitions$exists <- dir.exists(definitions$path)
  definitions$raster_count <- vapply(
    seq_len(nrow(definitions)),
    function(i) if (!definitions$exists[[i]]) 0L else length(list.files(
      definitions$path[[i]], pattern = "_bin[.]tif$",
      recursive = FALSE, ignore.case = TRUE, no.. = TRUE
    )),
    integer(1L)
  )
  list(inputs = inputs, sdm_directories = definitions)
}

# Select only the settings owned by one stage, then overlay the Rmd's explicit
# operational choices. This single precedence rule avoids hidden defaults and
# prevents unrelated settings from reaching stage-specific validators.
project_config_params <- function(configuration, stage, operational = list()) {
  assert(
    inherits(configuration, "project_configuration"),
    "configuration must come from read_project_config()."
  )
  assert(is.list(operational), "operational parameters must be a list.")
  allowed_operations <- project_config_operational_keys(stage)
  unknown_operations <- setdiff(names(operational), allowed_operations)
  assert(!length(unknown_operations), paste0(
    "Unsupported ", stage, " operational parameter(s): ",
    paste(unknown_operations, collapse = ", ")
  ))
  required <- project_config_stage_keys(stage)
  missing <- setdiff(required, names(configuration$values))
  assert(!length(missing), paste0(
    "Project configuration is missing ", stage, " key(s): ",
    paste(missing, collapse = ", ")
  ))
  params <- configuration$values[required]
  for (name in names(operational)) params[name] <- list(operational[[name]])
  params
}

rmd_project_config <- function(rmd_params, stage,
                               root = project_paths()$root) {
  assert(is.list(rmd_params), "Rmd parameters must be a list.")
  config_file <- rmd_params$config_file %||% "config/madagascar.yml"
  operational <- rmd_params[setdiff(names(rmd_params), "config_file")]
  configuration <- read_project_config(config_file, root)
  list(
    configuration = configuration,
    params = project_config_params(configuration, stage, operational)
  )
}
