# Stage 5.2 configuration and read-only preflight.
#
# The Stage 5.2 workflow owns loading this definition-only module. Configuration
# derives canonical Zonation paths and validates operation prerequisites without
# opening rasters, attaching packages, or writing artifacts.

stage52_config <- function(params, paths = project_paths()) {
  validate_parameter_names(
    params, c("mode", "zonation_settings_template"), "Stage 5.2"
  )
  mode <- validate_scalar_choice(
    params$mode, c("inspect", "inputs", "feature_list"), "params$mode"
  )
  methods <- benchmark_rank_methods()
  artifact_paths <- list(
    patch_raster_dir = paths$patches,
    patch_lookup_rds = paths$patch_lookup %||%
      file.path(paths$clean, "all_patch_lookup.rds"),
    binary_patch_raster_dir = paths$binary_patches,
    zonation_parent_dir = paths$zonation,
    settings_template = project_input_path(
      params$zonation_settings_template,
      file.path(paths$raw, "settings.z5.txt"),
      paths$root,
      "params$zonation_settings_template"
    ),
    settings_file = file.path(paths$zonation, "settings.z5.txt")
  )
  artifact_paths$feature_list <- paths$zonation_feature_list
  artifact_paths$zonation_method_dirs <- stats::setNames(
    file.path(artifact_paths$zonation_parent_dir, toupper(methods)), methods
  )
  artifact_paths$rankmap_sources <- stats::setNames(
    file.path(artifact_paths$zonation_method_dirs, "rankmap.tif"), methods
  )
  list(mode = mode, paths = artifact_paths)
}

stage52_application_config <- function(params, context) {
  stage52_config(
    params[c("mode", "zonation_settings_template")], context$paths
  )
}

# Discover the authoritative top-level feature set for feature-list mode.
# This metadata-only scan excludes directories and nested work products.
stage52_existing_feature_paths <- function(config) {
  need_dir(
    config$paths$binary_patch_raster_dir, "Stage 5.2 binary feature directory"
  )
  paths <- list.files(
    config$paths$binary_patch_raster_dir, pattern = "[.]tif$",
    full.names = TRUE, recursive = FALSE, ignore.case = TRUE, no.. = TRUE
  )
  if (length(paths)) {
    info <- file.info(paths)
    paths <- paths[!is.na(info$isdir) & !info$isdir]
  }
  assert(length(paths) > 0L, paste0(
    "No top-level .tif feature rasters were found in: ",
    config$paths$binary_patch_raster_dir
  ))
  sort(normalizePath(paths, winslash = "/", mustWork = TRUE))
}

describe_stage52_config <- function(config) {
  log_msg("Stage 5.2 | mode=", config$mode)
  log_msg(
    "Stage 5.2 | patch_dir=", config$paths$patch_raster_dir,
    "| binary_dir=", config$paths$binary_patch_raster_dir,
    "| feature_list=", config$paths$feature_list,
    "| settings=", config$paths$settings_file
  )
  log_msg("Stage 5.2 | zonation_outputs=", config$paths$zonation_parent_dir)
  invisible(config)
}

check_stage52_preflight <- function(config) {
  if (identical(config$mode, "inspect")) return(invisible(TRUE))
  need_file(config$paths$settings_template, "Zonation settings template")
  if (identical(config$mode, "feature_list")) {
    stage52_existing_feature_paths(config)
    return(invisible(TRUE))
  }
  need_dir(config$paths$patch_raster_dir, "Stage 5 patch raster directory")
  need_file(config$paths$patch_lookup_rds, "Stage 5 patch lookup")
  invisible(TRUE)
}
