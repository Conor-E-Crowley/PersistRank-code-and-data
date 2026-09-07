# Stage 4 configuration and mode-specific preflight.
#
# This module is the sole interpreter of Stage 4 YAML. It derives every
# project-owned path from project_paths(), validates the external SDM parent,
# and keeps figure mode independent of all build-only inputs and services.
# Scientific calculations and file writes belong to the downstream modules.


species_abort <- function(...) stop(paste0(...), call. = FALSE)

species_table_curve_columns <- function(curves) {
  curves <- as.character(curves)
  assert(
    length(curves) > 0L &&
      all(!is.na(curves)) &&
      all(curves %in% persistence_curves()),
    paste0(
      "species-table curve labels must be drawn from: ",
      paste(persistence_curves(), collapse = ", "), "."
    )
  )

  as.vector(rbind(paste0("alpha_", curves), paste0("beta_", curves)))
}

validate_species_table_curve_parameters <- function(species_table,
                                                    curves,
                                                    rows = NULL,
                                                    label = "species_table.csv",
                                                    rows_label = "selected species") {
  assert(is.data.frame(species_table), paste0(label, " must be a data.frame or tibble."))
  curve_cols <- species_table_curve_columns(curves)
  need_cols(species_table, curve_cols, label)

  if (is.null(rows)) {
    rows <- rep(TRUE, nrow(species_table))
  }
  assert(
    is.logical(rows) &&
      length(rows) == nrow(species_table) &&
      !any(is.na(rows)),
    paste0(label, " validation row mask must be TRUE/FALSE for each row.")
  )

  if (!any(rows)) return(invisible(TRUE))

  species_names <- if ("scientificName" %in% names(species_table)) {
    as.character(species_table$scientificName)
  } else {
    paste0("row ", seq_len(nrow(species_table)))
  }

  bad_details <- character(0)
  for (col in curve_cols) {
    values <- suppressWarnings(as.numeric(species_table[[col]]))
    bad <- rows & (!is.finite(values) | values <= 0)
    if (any(bad)) {
      examples <- utils::head(species_names[bad], 8L)
      bad_details <- c(
        bad_details,
        paste0(
          "  ", col, ": ", sum(bad), " invalid ",
          rows_label, " (examples: ", paste(examples, collapse = "; "),
          if (sum(bad) > length(examples)) "; ..." else "",
          ")"
        )
      )
    }
  }

  if (length(bad_details)) {
    species_abort(
      label,
      " contains missing, non-finite, or non-positive Gompertz curve parameters for ",
      rows_label,
      ".\n",
      paste(bad_details, collapse = "\n"),
      "\nThese alpha_*/beta_* columns are required for regenerated canonical species_table.csv outputs. ",
      "If this is the checked-in incomplete species table, regenerate Stages 2-4 before Stage 4 figures or Stage 6."
    )
  }

  invisible(TRUE)
}

required_species_table_packages <- function() {
  c("readr", "readxl", "dplyr", "stringr", "tibble", "tools")
}

load_species_table_packages <- function() {
  load_packages(required_species_table_packages())
}

required_species_table_figure_packages <- function() {
  c("readr", "dplyr", "stringr", "tibble", "ggplot2", "scales")
}

load_species_table_figure_packages <- function() {
  load_packages(required_species_table_figure_packages())
}

species_table_sdm_methods <- function(sdm) {
  sdm <- validate_priority_sdm_tag(sdm)
  switch(
    sdm,
    ppm = "PPM",
    rangebag = "RangeBag",
    ppm_rangebag = c("PPM", "RangeBag")
  )
}

# One canonical mapping owns the conventional SDM folder contract used by
# Stage 4 discovery, application manifests, and Stage 8 inspection.
selected_sdm_directories <- function(parent, taxa, sdm) {
  selection <- validate_taxa_selector(taxa, "taxa")
  methods <- species_table_sdm_methods(sdm)
  classes <- c(
    if (selection$selected_mammals) "Mammalia",
    if (selection$selected_birds) "Aves"
  )
  directories <- data.frame(
    taxon_class = c("Mammalia", "Mammalia", "Aves", "Aves"),
    sdm_method = c("PPM", "RangeBag", "PPM", "RangeBag"),
    directory = c(
      "mammal_ppm_bin", "mammal_rangebag_bin",
      "bird_ppm_bin", "bird_rangebag_bin"
    ),
    stringsAsFactors = FALSE
  )
  directories <- directories[
    directories$taxon_class %in% classes &
      directories$sdm_method %in% methods,
    , drop = FALSE
  ]
  directories$dir_path <- file.path(parent, directories$directory)
  directories
}

# Raster paths in an SDM index are portable: absolute paths are kept, while
# relative paths are interpreted beside the index CSV.
resolve_sdm_index_raster_paths <- function(x, index_file) {
  x <- trimws(as.character(x))
  assert(length(x) > 0L && all(!is.na(x) & nzchar(x)),
         "SDM index raster_path values must be non-empty.")
  absolute <- grepl("^/|^[A-Za-z]:[/\\\\]", x)
  x[!absolute] <- file.path(dirname(index_file), x[!absolute])
  normalizePath(x, winslash = "/", mustWork = FALSE)
}

stage4_config <- function(params,
                          paths = project_paths(),
                          contract = canonical_analysis_contract(),
                          iucn_mode = "refresh",
                          shared_paths = project_paths(),
                          gompertz_models_file = NULL,
                          habitat_rows = NULL) {
  contract <- validate_analysis_contract(contract, "Stage 4 abundance contract")
  iucn_mode <- validate_scalar_choice(
    iucn_mode, c("cache_only", "cache_or_query", "refresh"), "iucn_mode"
  )
  mode <- validate_scalar_choice(
    params$mode,
    c("figure", "build"),
    "params$mode"
  )
  selection <- validate_taxa_selector(params$taxa, "params$taxa")
  build_figure_only <- mode == "figure"

  sdm <- validate_priority_sdm_tag(params$sdm, "params$sdm")
  iucn_pause_seconds <- validate_scalar_number(
    params$iucn_pause_seconds %||% 2,
    "params$iucn_pause_seconds",
    minimum = 0
  )
  sdm_parent_dir <- project_input_path(
    params$sdm_parent_dir, ".", paths$root, "params$sdm_parent_dir"
  )
  rasters <- selected_sdm_directories(sdm_parent_dir, selection$value, sdm)
  sdm_input_mode <- if (is.null(params$sdm_index_file) ||
                         !length(params$sdm_index_file)) "folders" else "index"
  sdm_index_file <- NULL
  if (identical(sdm_input_mode, "index")) {
    assert(!is.null(params$sdm_index_file) && length(params$sdm_index_file),
           "params$sdm_index_file is required when an SDM index is selected.")
    sdm_index_file <- project_input_path(
      params$sdm_index_file, root = paths$root, label = "params$sdm_index_file"
    )
  }

  curves <- persistence_curves()

  selected_file <- function(value, default, label) {
    project_input_path(value, default, paths$root, label)
  }
  paths <- list(
    raw = list(
      summary = selected_file(
        params$species_list_file, file.path(paths$raw, "simple_summary.csv"),
        "params$species_list_file"
      ),
      synonyms = selected_file(
        params$synonyms_file, file.path(paths$raw, "synonyms.csv"),
        "params$synonyms_file"
      ),
      mammal_traits = selected_file(
        params$mammal_traits_file, file.path(paths$raw, "mammal_data.txt"),
        "params$mammal_traits_file"
      ),
      bird_traits = selected_file(
        params$bird_traits_file, file.path(paths$raw, "bird_data.txt"),
        "params$bird_traits_file"
      ),
      bird_generation_lengths = selected_file(
        params$bird_generation_lengths_file,
        file.path(paths$raw, "cobi13486-sup-0004-tables4.xlsx"),
        "params$bird_generation_lengths_file"
      ),
      input_synonyms = selected_file(
        params$curated_synonyms_file, file.path(paths$raw, "input_synonyms.csv"),
        "params$curated_synonyms_file"
      ),
      random_effects = selected_file(
        params$random_effects_file, file.path(paths$raw, "random_effects.csv"),
        "params$random_effects_file"
      )
    ),
    rasters = rasters,
    clean = list(
      gompertz_models = selected_file(
        gompertz_models_file, file.path(paths$clean, "persistence_curve_models.rds"),
        "gompertz_models_file"
      ),
      species_table = file.path(paths$clean, "species_table.csv"),
      iucn_cache = shared_paths$iucn_cache %||%
        file.path(shared_paths$cache %||% shared_paths$clean, "iucn_habitat_cache.rds")
    ),
    figure_dir = paths$figures,
    area_curve = file.path(paths$figures, "area_curve.png")
  )
  paths$intended_outputs <- if (build_figure_only) {
    c(area_curve = paths$area_curve)
  } else {
    c(species_table = paths$clean$species_table)
  }

  list(
    mode = mode,
    build_figure_only = build_figure_only,
    taxa = selection$value,
    fit_mammal_curves = selection$selected_mammals,
    fit_bird_curves = selection$selected_birds,
    verbose = validate_scalar_logical(params$verbose, "params$verbose"),
    sdm = sdm,
    sdm_input_mode = sdm_input_mode,
    sdm_index_file = sdm_index_file,
    curves = curves,
    main_curve = main_persistence_curve(),
    iucn_pause_seconds = iucn_pause_seconds,
    iucn_mode = iucn_mode,
    habitat_rows = habitat_rows,
    contract = contract,
    figure_focal_species = "Fossa fossana",
    paths = paths
  )
}

# Translate a shared application context to the Stage 4 configuration. The
# application taxa selector owns both inventory and coefficient selection.
stage4_application_config <- function(params, context, frozen_inputs = NULL,
                                      iucn_mode = NULL) {
  assert(is.list(context) && !is.null(context$contract) && !is.null(context$paths),
         "context must be a shared application context.")
  stage_params <- params[c(
    "mode", "verbose", "taxa", "sdm", "sdm_index_file",
    "iucn_pause_seconds", "sdm_parent_dir", "species_list_file",
    "synonyms_file", "curated_synonyms_file", "mammal_traits_file",
    "bird_traits_file", "bird_generation_lengths_file", "random_effects_file"
  )]
  stage_params$taxa <- context$taxa
  stage_params$sdm <- context$sdm
  stage_params$sdm_index_file <- context$sdm_index$path %||%
    params$sdm_index_file %||% NULL
  stage_params$sdm_parent_dir <- context$sdm_parent_dir %||%
    params$sdm_parent_dir %||% context$base_paths$root
  stage_params$species_list_file <- context$inputs$species_list$path %||%
    params$species_list_file
  stage_params$synonyms_file <- context$inputs$synonyms$path %||% params$synonyms_file
  stage_params$curated_synonyms_file <- context$inputs$curated_synonyms$path %||%
    params$curated_synonyms_file
  stage_params$mammal_traits_file <- params$mammal_traits_file
  stage_params$bird_traits_file <- params$bird_traits_file
  stage_params$bird_generation_lengths_file <- params$bird_generation_lengths_file
  stage_params$random_effects_file <- params$random_effects_file
  requested_iucn_mode <- iucn_mode %||% params$iucn_mode %||% "refresh"
  stage4_config(
    stage_params, paths = context$paths, contract = context$contract,
    iucn_mode = requested_iucn_mode,
    shared_paths = context$base_paths,
    gompertz_models_file = file.path(
      context$model$stage3, "persistence_curve_models.rds"
    ),
    habitat_rows = frozen_inputs$habitat_rows %||% NULL
  )
}

check_stage4_preflight <- function(config) {
  if (config$build_figure_only) {
    need_file(
      config$paths$clean$species_table,
      "existing Stage 4 species table for figure-only mode"
    )
    ensure_writable_dir(config$paths$figure_dir, "Stage 4 figure directory")
    return(invisible(TRUE))
  }

  ensure_writable_dir(dirname(config$paths$clean$species_table), "species table output directory")

  invisible(lapply(config$paths$raw, need_file))
    if (identical(config$sdm_input_mode, "folders")) {
      invisible(lapply(config$paths$rasters$dir_path, need_dir))
    } else {
      need_file(config$sdm_index_file, "SDM index CSV")
    }
  if (config$fit_mammal_curves || config$fit_bird_curves) {
    need_file(config$paths$clean$gompertz_models, "Stage 3 Gompertz LOESS model")
  }
  invisible(TRUE)
}

describe_stage4_config <- function(config) {
  log_msg(
    "Stage 4 | mode=", config$mode,
    "| taxa=", config$taxa,
    "| sdm=", config$sdm
  )
  if (config$build_figure_only) {
    log_msg(
      "Stage 4 | action=figure-only",
      "| source=", config$paths$clean$species_table,
      "| curve=", config$main_curve,
      "| focal_species=", config$figure_focal_species
    )
  } else {
    log_msg(
      "Stage 4 | curves=", paste(config$curves, collapse = ","),
      "| fitted mammal curves=", config$fit_mammal_curves,
      "| fitted bird curves=", config$fit_bird_curves,
      "| IUCN pause=", config$iucn_pause_seconds, "s",
      "| IUCN mode=", config$iucn_mode
    )
  }
  log_msg(
    "Stage 4 | intended outputs=", paste(config$paths$intended_outputs, collapse = ", ")
  )
  if (config$build_figure_only) {
    log_msg(
      "Stage 4 | retained input=", config$paths$clean$species_table,
      "| no species-table or IUCN rebuild"
    )
    return(invisible(config))
  }
  invisible(config)
}
