# Stage 5 configuration, canonical paths, and clumping-backend discovery.
#
# Project-owned locations and scientific thresholds are intentionally derived
# from shared contracts. Only a machine-specific GRASS location is accepted
# from YAML.

required_patch_packages <- function(include_figures = FALSE) {
  pkgs <- c("terra", "sf", "igraph", "dplyr", "readr", "stringr", "tools")
  if (isTRUE(include_figures)) {
    pkgs <- c(pkgs, "ggplot2", "tidyterra", "cowplot", "rnaturalearth", "scales")
  }
  unique(pkgs)
}

load_patch_packages <- function(include_figures = FALSE) {
  pkgs <- required_patch_packages(include_figures)
  load_packages(pkgs)

  sf::sf_use_s2(TRUE)
  invisible(pkgs)
}

patch_paths <- function(project = project_paths(), landcover_file = NULL) {
  landcover_file <- project_input_path(
    landcover_file, file.path(project$raw, "esacci_2022_pfts.tif"),
    project$root, "landcover_file"
  )
  list(
    species_csv = file.path(project$clean, "species_table.csv"),
    landcover_tif = landcover_file,
    patch_dir = project$patches,
    patch_lookup_rds = project$patch_lookup %||% file.path(project$clean, "all_patch_lookup.rds"),
    connectivity_rds = project$connectivity %||% file.path(project$clean, "all_connectivity.rds"),
    metadata_rds = project$stage5_metadata %||% file.path(project$clean, "stage5_build_metadata.rds"),
    checkpoint_dir = project$stage5_checkpoints,
    runtime_log = file.path(project$data, "Logs", "stage5.log"),
    figure_dir = project$spatial_figures %||% project$figures,
    figure = file.path(project$spatial_figures %||% project$figures, "spatial_process.png")
  )
}

patch_roi_bounds <- function(xmin = 43.18, xmax = 50.56,
                             ymin = -25.64, ymax = -11.89) {
  values <- c(
    xmin = validate_scalar_number(xmin, "roi_xmin"),
    xmax = validate_scalar_number(xmax, "roi_xmax"),
    ymin = validate_scalar_number(ymin, "roi_ymin"),
    ymax = validate_scalar_number(ymax, "roi_ymax")
  )
  assert(values[["xmin"]] < values[["xmax"]], "roi_xmin must be less than roi_xmax.")
  assert(values[["ymin"]] < values[["ymax"]], "roi_ymin must be less than roi_ymax.")
  structure(
    values,
    class = c("patch_roi_bounds", "numeric")
  )
}

patch_roi_from_params <- function(params) {
  patch_roi_bounds(
    params$roi_xmin %||% 43.18,
    params$roi_xmax %||% 50.56,
    params$roi_ymin %||% -25.64,
    params$roi_ymax %||% -11.89
  )
}

patch_roi <- function(bounds = patch_roi_bounds()) {
  assert(inherits(bounds, "patch_roi_bounds"), "bounds must be the shared Stage 5 ROI.")

  terra::ext(
    bounds[["xmin"]],
    bounds[["xmax"]],
    bounds[["ymin"]],
    bounds[["ymax"]]
  )
}

patch_study_area_from_params <- function(params, root = ".") {
  mode <- validate_scalar_choice(
    params$study_area_mode %||% "bounds",
    c("bounds", "vector", "full_raster"),
    "params$study_area_mode"
  )
  if (identical(mode, "bounds")) {
    return(list(mode = mode, bounds = patch_roi_from_params(params), file = NULL, layer = NULL))
  }
  if (identical(mode, "full_raster")) {
    return(list(mode = mode, bounds = NULL, file = NULL, layer = NULL))
  }
  file <- project_input_path(
    params$study_area_file, root = root, label = "params$study_area_file"
  )
  layer <- params$study_area_layer
  if (!is.null(layer) && length(layer)) {
    layer <- validate_scalar_string(layer, "params$study_area_layer")
  } else layer <- NULL
  list(mode = mode, bounds = NULL, file = file, layer = layer)
}

patch_study_area_label <- function(study_area) {
  if (identical(study_area$mode, "bounds")) {
    return(paste(unname(study_area$bounds), collapse = ","))
  }
  if (identical(study_area$mode, "vector")) return(study_area$file)
  "complete land-cover raster"
}

validate_patch_backend <- function(x, label = "params$clump_backend") {
  value <- tolower(validate_scalar_string(x, label))
  assert(
    value %in% c("auto", "fasterraster", "terra"),
    paste0(label, " must be one of: auto, fasterRaster, terra.")
  )
  if (identical(value, "fasterraster")) "fasterRaster" else value
}

validate_patch_grass_param <- function(x, label = "params$grass_dir") {
  value <- optional_path(x, label)
  if (!is.na(value)) need_dir(value, label)
  value
}

stage5_config <- function(params, paths = project_paths(),
                          contract = canonical_analysis_contract(),
                          species_table_file = NULL,
                          landcover_fingerprint = NULL) {
  contract <- validate_analysis_contract(contract, "Stage 5 abundance contract")
  mode <- validate_scalar_choice(
    params$mode,
    c("inspect", "resume", "restart"),
    "params$mode"
  )
  selection <- validate_taxa_selector(params$taxa, "params$taxa")
  verbose <- validate_scalar_logical(params$verbose, "params$verbose")
  artifact_paths <- patch_paths(paths, params$landcover_file %||% NULL)
  if (!is.null(species_table_file)) {
    artifact_paths$species_csv <- normalizePath(
      validate_path_param(species_table_file, "species_table_file"),
      winslash = "/", mustWork = FALSE
    )
  }
  assert(
    !identical(
      normalizePath(artifact_paths$patch_lookup_rds, winslash = "/", mustWork = FALSE),
      normalizePath(artifact_paths$connectivity_rds, winslash = "/", mustWork = FALSE)
    ),
    "Stage 5 patch lookup and connectivity paths must be different files."
  )
  list(
    mode = mode,
    taxa = selection$value,
    selected_mammals = selection$selected_mammals,
    selected_birds = selection$selected_birds,
    verbose = verbose,
    contract = contract,
    clump_backend = validate_patch_backend(params$clump_backend),
    grass_dir = validate_patch_grass_param(params$grass_dir),
    landcover_fingerprint = landcover_fingerprint,
    study_area = patch_study_area_from_params(params, paths$root),
    roi_bounds = if (identical(params$study_area_mode %||% "bounds", "bounds")) {
      patch_roi_from_params(params)
    } else NULL,
    paths = artifact_paths
  )
}

stage51_config <- function(params, paths = project_paths(),
                           contract = canonical_analysis_contract(),
                           species_table_file = NULL,
                           target_species = "Cryptoprocta ferox") {
  contract <- validate_analysis_contract(contract, "Stage 5.1 abundance contract")
  artifact_paths <- patch_paths(paths, params$landcover_file %||% NULL)
  if (!is.null(species_table_file)) {
    artifact_paths$species_csv <- normalizePath(
      validate_path_param(species_table_file, "species_table_file"),
      winslash = "/", mustWork = FALSE
    )
  }
  list(
    target_species = if (is.null(target_species)) NULL else
      validate_scalar_string(target_species, "target_species"),
    contract = contract,
    clump_backend = validate_patch_backend(params$clump_backend),
    grass_dir = validate_patch_grass_param(params$grass_dir),
    study_area = patch_study_area_from_params(params, paths$root),
    roi_bounds = if (identical(params$study_area_mode %||% "bounds", "bounds")) {
      patch_roi_from_params(params)
    } else NULL,
    paths = list(
      species_csv = artifact_paths$species_csv,
      landcover_tif = artifact_paths$landcover_tif,
      figure_dir = artifact_paths$figure_dir,
      figure = artifact_paths$figure
    )
  )
}

# Translate application-owned selectors and paths to Stage 5 configuration.
stage5_application_config <- function(params, context,
                                      landcover_fingerprint = NULL) {
  stage_params <- params[c(
    "mode", "verbose", "clump_backend", "grass_dir", "landcover_file",
    "study_area_mode", "study_area_file", "study_area_layer",
    "roi_xmin", "roi_xmax", "roi_ymin", "roi_ymax"
  )]
  stage_params$taxa <- context$taxa
  stage5_config(
    stage_params, context$paths, context$contract,
    landcover_fingerprint = landcover_fingerprint
  )
}

# Translate the same application spatial settings to the explanatory Stage 5.1
# configuration without reading rasters or attaching spatial packages.
stage51_application_config <- function(params, context) {
  stage_params <- params[c(
    "clump_backend", "grass_dir", "landcover_file", "study_area_mode",
    "study_area_file", "study_area_layer", "roi_xmin", "roi_xmax",
    "roi_ymin", "roi_ymax"
  )]
  stage51_config(
    stage_params, context$paths, context$contract,
    target_species = params$process_figure_species %||%
      context$process_figure_species %||% NULL
  )
}

check_stage5_preflight <- function(config) {
  need_file(config$paths$species_csv, "Stage 4 species table")
  need_file(config$paths$landcover_tif, "ESA CCI land-cover raster")
  if (!identical(config$mode, "inspect")) {
    ensure_writable_dir(dirname(config$paths$patch_dir), "Stage 5 clean-data directory")
    ensure_writable_dir(dirname(config$paths$checkpoint_dir), "Stage 5 checkpoint parent")
  }
  invisible(TRUE)
}

check_stage51_preflight <- function(config) {
  need_file(config$paths$species_csv, "Stage 4 species table")
  need_file(config$paths$landcover_tif, "ESA CCI land-cover raster")
  ensure_writable_dir(config$paths$figure_dir, "Stage 5.1 figure directory")
  invisible(TRUE)
}

resolve_patch_clump_backend <- function(requested, grass_dir = NULL, label = "Stage 5") {
  label <- validate_scalar_string(label, "backend log label")
  requested <- validate_patch_backend(requested, "clump_backend")
  if (identical(requested, "terra")) {
    log_msg(label, "| clump backend=terra (explicit)")
    return("terra")
  }

  initialized <- suppressWarnings(init_faster_raster(grass_dir))
  if (!is.na(initialized)) {
    log_msg(label, "| clump backend=fasterRaster | GRASS=", initialized)
    return("fasterRaster")
  }
  if (identical(requested, "fasterRaster")) {
    patch_abort(
      "clump_backend=fasterRaster was requested, but fasterRaster/GRASS could not be initialized."
    )
  }
  warning(
    paste0(label, " is using terra::patches() because fasterRaster/GRASS is unavailable; clumping may be substantially slower."),
    call. = FALSE
  )
  log_msg(label, "| clump backend=terra (automatic fallback)")
  "terra"
}

describe_stage5_config <- function(config, label = "Stage 5") {
  if (!is.null(config$paths$patch_lookup_rds)) {
    log_msg(
      label, "| mode=", config$mode,
      "| backend=", config$clump_backend,
      "| study_area=", patch_study_area_label(config$study_area),
      "| taxa=", config$taxa,
      "| mammals=", config$selected_mammals,
      "| birds=", config$selected_birds,
      "| checkpoint_dir=", config$paths$checkpoint_dir,
      "| patch_dir=", config$paths$patch_dir
    )
  } else {
    log_msg(
      label, "| backend=", config$clump_backend,
      "| study_area=", patch_study_area_label(config$study_area),
      "| target_species=", config$target_species,
      "| figure=", config$paths$figure
    )
  }
  invisible(config)
}
