# Stage 5.2 Zonation input preparation workflow.
#
# `inputs` rebuilds every binary feature raster from the finalized Stage 5
# artifacts and commits the complete set with its feature list and settings
# file. `feature_list` only refreshes those two text files from the binary
# rasters currently on disk; it never opens or changes a raster. `inspect` is
# read-only.

if (!exists("project_source", mode = "function")) {
  source(file.path("R", "project_loader.R"))
}
project_source(c(
  # Shared contracts and application adapters; raster work remains inside run.
  "R/project_utils.R", "R/project_artifacts.R", "R/project_rmd.R",
  "R/project_runtime_log.R", "R/project_transactions.R",
  "R/project_paths.R", "R/project_config.R", "R/demographic_contract.R",
  "R/analysis_contract.R", "R/priority_run_config.R",
  "R/species_table_config.R", "R/patch_contract.R", "R/patch_config.R",
  "R/application_storage.R", "R/storage_artifacts.R", "R/model_lifecycle.R",
  "R/application_state.R",
  "R/application_context.R", "R/application_rmd.R", "R/stage52_config.R"
))

stage52_normalize_output_paths <- function(paths) {
  file.path(
    normalizePath(dirname(paths), winslash = "/", mustWork = TRUE),
    basename(paths)
  )
}


inspect_stage52 <- function(config) {
  check_stage52_preflight(config)
  artifact_paths <- c(
    "patch lookup" = config$paths$patch_lookup_rds,
    "patch raster directory" = config$paths$patch_raster_dir,
    "binary raster directory" = config$paths$binary_patch_raster_dir,
    "feature list" = config$paths$feature_list,
    "settings template" = config$paths$settings_template,
    "settings file" = config$paths$settings_file,
    stats::setNames(
      config$paths$rankmap_sources,
      paste0("Zonation ", toupper(names(config$paths$rankmap_sources)), " rankmap")
    )
  )
  status <- project_artifact_status(artifact_paths)
  status$artifact <- names(artifact_paths)
  status[, c("artifact", "path", "exists", "type", "size_bytes", "modified")]
}

stage52_raster_values <- function(raster, label) {
  assert(terra::nlyr(raster) == 1L, paste0(label, " must contain exactly one layer."))
  values <- as.numeric(terra::unique(raster, na.rm = TRUE)[[1L]])
  values <- values[is.finite(values)]
  assert(length(values) > 0L, paste0(label, " contains no retained patch cells."))
  assert(all(values > 0 & values == floor(values)), paste0(label, " contains invalid patch IDs."))
  sort(unique(as.integer(values)))
}

write_binary_patch_raster <- function(input_path, output_path, template = NULL,
                                      expected_patch_ids = NULL) {
  input <- terra::rast(input_path)
  source_ids <- stage52_raster_values(input, paste0("Stage 5 patch raster ", input_path))
  if (!is.null(expected_patch_ids)) {
    assert(
      identical(source_ids, sort(unique(as.integer(expected_patch_ids)))),
      paste0("Stage 5 patch raster IDs do not match the lookup: ", input_path)
    )
  }
  if (!is.null(template)) {
    assert(
      terra::compareGeom(input, template, stopOnError = FALSE),
      paste0("Stage 5 patch raster geometry differs from the common template: ", input_path)
    )
  }
  assert(!file.exists(output_path),
         paste0("Staged binary patch raster already exists: ", output_path))
  ensure_writable_dir(dirname(output_path), "Stage 5.2 raster staging directory")
  binary <- terra::ifel(is.na(input), NA, 1L)
  names(binary) <- names(input)
  terra::writeRaster(
    binary, output_path, overwrite = FALSE, datatype = "INT2S", NAflag = -32768
  )
  check <- terra::rast(output_path)
  assert(terra::compareGeom(check, input, stopOnError = FALSE),
         "Temporary binary patch raster changed geometry.")
  values <- as.numeric(terra::unique(check, na.rm = TRUE)[[1L]])
  assert(identical(sort(unique(values[is.finite(values)])), 1),
         "Temporary binary patch raster contains values other than 1 and NA.")
  rm(check)
  invisible(output_path)
}

write_zonation_feature_list <- function(feature_paths, validation_paths, output_path) {
  assert(!file.exists(output_path),
         paste0("Staged Zonation feature list already exists: ", output_path))
  assert(length(feature_paths) == length(validation_paths) &&
           all(file.exists(validation_paths)),
         "Every feature-list entry must have a validated staged raster.")
  ensure_writable_dir(
    dirname(output_path), "Stage 5.2 feature-list staging directory"
  )
  normalized <- stage52_normalize_output_paths(feature_paths)
  writeLines(
    c(
      '"weight" "filename"',
      sprintf("1 %s", shQuote(normalized, type = "cmd"))
    ),
    output_path
  )
  check <- readLines(output_path, warn = FALSE)
  assert(
    length(check) == length(feature_paths) + 1L,
    "Temporary Zonation feature list has an invalid row count."
  )
  invisible(output_path)
}

# Preserve the supplied Zonation settings and replace only its feature-list
# pointer. Requiring exactly one key prevents silently producing an ambiguous
# or structurally different settings file.
stage52_settings_lines <- function(template_path, feature_list_path) {
  need_file(template_path, "Zonation settings template")
  lines <- readLines(template_path, warn = FALSE)
  key <- grep(
    "^[[:space:]]*feature list file[[:space:]]*=",
    lines,
    ignore.case = TRUE
  )
  assert(length(key) == 1L,
         "Zonation settings template must contain exactly one 'feature list file =' line.")
  lines[[key]] <- paste0(
    "feature list file = ",
    stage52_normalize_output_paths(feature_list_path)
  )
  lines
}

write_zonation_settings <- function(template_path, feature_list_path, output_path) {
  assert(!file.exists(output_path),
         paste0("Staged Zonation settings file already exists: ", output_path))
  lines <- stage52_settings_lines(template_path, feature_list_path)
  ensure_writable_dir(dirname(output_path), "Stage 5.2 settings staging directory")
  writeLines(lines, output_path)
  assert(identical(readLines(output_path, warn = FALSE), lines),
         "Temporary Zonation settings file validation failed.")
  invisible(output_path)
}

# Atomically refresh feature_list.txt and the settings file that points to it.
# The existing binary rasters are treated as inputs and are never opened.
run_stage52_feature_list <- function(config, started, .rename_file = file.rename) {
  feature_paths <- stage52_existing_feature_paths(config)
  ensure_writable_dir(config$paths$zonation_parent_dir, "Zonation parent directory")

  staging_dir <- tempfile(
    "stage52_text_staging_",
    tmpdir = config$paths$zonation_parent_dir
  )
  dir.create(staging_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(staging_dir, recursive = TRUE, force = TRUE), add = TRUE)
  staged_feature_path <- file.path(staging_dir, basename(config$paths$feature_list))
  staged_settings_path <- file.path(staging_dir, basename(config$paths$settings_file))
  write_zonation_feature_list(
    feature_paths = feature_paths,
    validation_paths = feature_paths,
    output_path = staged_feature_path
  )
  write_zonation_settings(
    template_path = config$paths$settings_template,
    feature_list_path = config$paths$feature_list,
    output_path = staged_settings_path
  )
  transaction <- stage5_file_set_transaction(
    c(staged_feature_path, staged_settings_path),
    c(config$paths$feature_list, config$paths$settings_file),
    overwrite = TRUE,
    rename_file = .rename_file
  )
  transaction$finalize()

  expected_lines <- c(
    '"weight" "filename"',
    sprintf("1 %s", shQuote(feature_paths, type = "cmd"))
  )
  assert(
    identical(readLines(config$paths$feature_list, warn = FALSE), expected_lines),
    "Committed Stage 5.2 feature_list.txt does not match the current binary features."
  )
  assert(
    identical(
      readLines(config$paths$settings_file, warn = FALSE),
      stage52_settings_lines(
        config$paths$settings_template,
        config$paths$feature_list
      )
    ),
    "Committed Stage 5.2 settings file does not match its template and feature list."
  )
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  log_msg(
    "Stage 5.2 complete | action=feature_list | features=", length(feature_paths),
    "| elapsed=", round(elapsed, 2), "s"
  )
  list(
    action = "feature_list",
    feature_list = config$paths$feature_list,
    settings_file = config$paths$settings_file,
    feature_paths = feature_paths,
    feature_count = length(feature_paths),
    elapsed_seconds = elapsed
  )
}

run_stage52 <- function(config, .rename_file = file.rename) {
  started <- Sys.time()
  assert(!identical(config$mode, "inspect"),
         "run_stage52() does not run in inspect mode; use inspect_stage52().")
  check_stage52_preflight(config)
  if (identical(config$mode, "feature_list")) {
    return(run_stage52_feature_list(config, started, .rename_file))
  }

  assert(identical(config$mode, "inputs"), "Unsupported Stage 5.2 run mode.")
  load_packages("terra")
  patch_lookup <- readRDS(config$paths$patch_lookup_rds)
  validate_patch_lookup_object(patch_lookup, "Stage 5.2 source patch lookup")
  species <- sort(unique(patch_squish(patch_lookup$scientificName)))
  assert(length(species) > 0L, "Stage 5 patch lookup contains no retained species.")
  input_paths <- file.path(
    config$paths$patch_raster_dir,
    vapply(species, patch_filename_from_scientific, character(1L))
  )
  missing <- input_paths[!file.exists(input_paths)]
  assert(!length(missing), paste0(
    "Missing Stage 5 patch raster(s):\n",
    paste(missing, collapse = "\n")
  ))
  ensure_writable_dir(config$paths$binary_patch_raster_dir, "binary patch raster directory")
  ensure_writable_dir(config$paths$zonation_parent_dir, "Zonation parent directory")
  for (path in config$paths$zonation_method_dirs) {
    ensure_writable_dir(path, "Zonation method directory")
  }
  output_paths <- file.path(config$paths$binary_patch_raster_dir, basename(input_paths))

  # Build beside the final directory, then promote the complete species set
  # and feature list together. Existing final artifacts remain usable if any
  # conversion or validation fails.
  staging_dir <- tempfile(
    "stage52_binary_staging_",
    tmpdir = config$paths$binary_patch_raster_dir
  )
  dir.create(staging_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(staging_dir, recursive = TRUE, force = TRUE), add = TRUE)
  staged_binary_paths <- file.path(staging_dir, basename(output_paths))
  staged_feature_path <- file.path(staging_dir, basename(config$paths$feature_list))
  staged_settings_path <- file.path(staging_dir, basename(config$paths$settings_file))

  template <- terra::rast(input_paths[[1L]])
  for (i in seq_along(input_paths)) {
    log_msg(
      "Stage 5.2 | raster", paste0(i, "/", length(input_paths)),
      "|", basename(input_paths[[i]])
    )
    write_binary_patch_raster(
      input_paths[[i]], staged_binary_paths[[i]], template,
      expected_patch_ids = patch_lookup$patch_id[
        patch_squish(patch_lookup$scientificName) == species[[i]]
      ]
    )
  }
  write_zonation_feature_list(
    feature_paths = output_paths,
    validation_paths = staged_binary_paths,
    output_path = staged_feature_path
  )
  write_zonation_settings(
    template_path = config$paths$settings_template,
    feature_list_path = config$paths$feature_list,
    output_path = staged_settings_path
  )
  transaction <- stage5_file_set_transaction(
    c(staged_binary_paths, staged_feature_path, staged_settings_path),
    c(output_paths, config$paths$feature_list, config$paths$settings_file),
    overwrite = TRUE,
    rename_file = .rename_file
  )
  transaction$finalize()
  log_msg(
    "Stage 5.2 complete | species=", length(species),
    "| elapsed=", round(as.numeric(difftime(Sys.time(), started, units = "secs")), 2), "s"
  )
  list(
    action = "inputs",
    species = species,
    input_paths = input_paths,
    binary_paths = output_paths,
    feature_list = config$paths$feature_list,
    settings_file = config$paths$settings_file,
    feature_paths = output_paths,
    feature_count = length(output_paths),
    method_dirs = config$paths$zonation_method_dirs,
    rankmaps_present = file.exists(config$paths$rankmap_sources),
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))
  )
}
