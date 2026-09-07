# Stage 5 species-level recovery and final artifact promotion.
#
# Checkpoints contain completed scientific species results only. They never
# contain the shared land-cover raster, habitat masks, or live terra objects;
# those inexpensive-to-reconstruct runtime objects are rebuilt after restart.


stage5_checkpoint_schema <- function() 2L

stage5_object_md5 <- function(x) {
  path <- tempfile("stage5_signature_", fileext = ".rds")
  on.exit(unlink(path, force = TRUE), add = TRUE)
  saveRDS(x, path, version = 3)
  unname(tools::md5sum(path))
}

stage5_file_manifest <- function(paths, label) {
  paths <- normalizePath(paths, winslash = "/", mustWork = TRUE)
  assert(!anyDuplicated(paths), paste0(label, " contains duplicate file paths."))
  info <- file.info(paths)
  assert(all(!is.na(info$size) & info$size > 0), paste0(label, " contains an empty file."))
  data.frame(
    path = paths,
    size = as.numeric(info$size),
    md5 = unname(tools::md5sum(paths)),
    stringsAsFactors = FALSE
  )
}

stage5_landcover_file_manifest <- function(path, fingerprint = NULL) {
  path <- normalizePath(
    validate_path_param(path, "Stage 5 land-cover input"),
    winslash = "/", mustWork = TRUE
  )
  if (is.null(fingerprint)) {
    return(stage5_file_manifest(path, "Stage 5 land-cover input"))
  }

  assert(is.list(fingerprint), "landcover_fingerprint must be a verified file record.")
  fingerprint_path <- normalizePath(
    validate_path_param(fingerprint$path, "landcover_fingerprint$path"),
    winslash = "/", mustWork = FALSE
  )
  assert(
    identical(fingerprint_path, path),
    "The trusted land-cover fingerprint refers to a different file."
  )
  current_bytes <- unname(as.numeric(file.info(path)$size))
  assert(
    length(fingerprint$bytes) == 1L && is.finite(fingerprint$bytes) &&
      fingerprint$bytes > 0 &&
      identical(as.numeric(fingerprint$bytes), current_bytes),
    "The land-cover file size changed after its fingerprint was verified."
  )
  assert(
    is.character(fingerprint$md5) && length(fingerprint$md5) == 1L &&
      !is.na(fingerprint$md5) && grepl("^[[:xdigit:]]{32}$", fingerprint$md5),
    "landcover_fingerprint$md5 must be a verified MD5 checksum."
  )
  data.frame(
    path = path,
    size = current_bytes,
    md5 = as.character(fingerprint$md5),
    stringsAsFactors = FALSE
  )
}

stage5_landcover_contract <- function(path, fingerprint = NULL) {
  file <- stage5_landcover_file_manifest(path, fingerprint)
  raster <- terra::rast(file$path[[1L]])
  list(
    file = file,
    nrow = terra::nrow(raster),
    ncol = terra::ncol(raster),
    nlyr = terra::nlyr(raster),
    resolution = unname(terra::res(raster)),
    extent = unname(as.vector(terra::ext(raster))),
    crs = terra::crs(raster, proj = TRUE),
    names = names(raster),
    categorical = unname(vapply(
      seq_len(terra::nlyr(raster)),
      function(i) isTRUE(terra::is.factor(raster[[i]])),
      logical(1)
    ))
  )
}

stage5_helper_contract <- function() {
  helper_paths <- file.path(
    "R",
    c(
      "analysis_contract.R", "patch_config.R", "patch_species_inputs.R",
      "patch_habitat.R", "patch_contract.R", "patch_processing.R",
      "stage5_checkpoints.R",
      "stage5_workflow.R"
    )
  )
  stage5_file_manifest(helper_paths, "Stage 5 helper contract")
}

# These files transport configuration, coordinate recovery, and emit logs.
# Their complete hashes remain in every manifest for provenance, but changes
# to them alone do not invalidate scientific outputs. Scientific changes in
# any of these modules must be accompanied by a Stage 5 schema increment.
stage5_runtime_helper_files <- function() {
  c("patch_config.R", "stage5_checkpoints.R", "stage5_workflow.R")
}

stage5_manifest_scientific_view <- function(manifest) {
  if (!is.list(manifest)) return(manifest)
  out <- manifest
  if (is.data.frame(out$helpers) && "path" %in% names(out$helpers)) {
    keep <- !basename(out$helpers$path) %in% stage5_runtime_helper_files()
    out$helpers <- out$helpers[keep, , drop = FALSE]
    rownames(out$helpers) <- NULL
  }
  out
}

stage5_package_contract <- function(backend) {
  packages <- c("terra", "sf", "igraph", "dplyr", "readr", "stringr")
  if (identical(backend, "fasterRaster")) {
    packages <- c(packages, "fasterRaster")
  }
  stats::setNames(
    vapply(packages, function(package) as.character(utils::packageVersion(package)), character(1)),
    packages
  )
}

build_stage5_manifest <- function(config, selected_species, backend) {
  sdm_paths <- unique(as.character(selected_species$raster_path))
  list(
    schema = stage5_checkpoint_schema(),
    taxa = config$taxa,
    species = as.character(selected_species$scientificName),
    species_table = stage5_file_manifest(
      config$paths$species_csv,
      "Stage 5 species table"
    ),
    sdm = stage5_file_manifest(sdm_paths, "Stage 5 SDM inputs"),
    landcover = stage5_landcover_contract(
      config$paths$landcover_tif,
      config$landcover_fingerprint %||% NULL
    ),
    habitat_crosswalk = stage5_object_md5(list(
      classes = landcover_class_map(),
      labels = habitat_label_map()
    )),
    study_area = config$study_area,
    thresholds = list(
      persistence_horizon_years = persistence_horizon_years(config$contract),
      quasi_extinction_abundance = quasi_extinction_abundance(config$contract),
      minimum_patch_abundance = minimum_patch_abundance(config$contract)
    ),
    backend = backend,
    helpers = stage5_helper_contract(),
    packages = stage5_package_contract(backend)
  )
}

stage5_manifest_compatible <- function(found, expected) {
  is.list(found) && identical(
    stage5_manifest_scientific_view(found),
    stage5_manifest_scientific_view(expected)
  )
}

stage5_manifest_path <- function(config) {
  file.path(config$paths$checkpoint_dir, "manifest.rds")
}

stage5_completed_dir <- function(config) {
  file.path(config$paths$checkpoint_dir, "completed")
}

stage5_species_checkpoint_name <- function(index, species) {
  stem <- tools::file_path_sans_ext(patch_filename_from_scientific(species))
  sprintf("%04d_%s", as.integer(index), stem)
}

stage5_species_checkpoint_dir <- function(config, index, species) {
  file.path(
    stage5_completed_dir(config),
    stage5_species_checkpoint_name(index, species)
  )
}

initialize_stage5_checkpoint <- function(config, manifest, restart = FALSE) {
  checkpoint_dir <- config$paths$checkpoint_dir
  if (isTRUE(restart) && dir.exists(checkpoint_dir)) {
    unlink(checkpoint_dir, recursive = TRUE, force = TRUE)
  }
  if (dir.exists(checkpoint_dir)) {
    path <- stage5_manifest_path(config)
    need_file(path, "Stage 5 checkpoint manifest")
    found <- readRDS(path)
    assert(
      stage5_manifest_compatible(found, manifest),
      paste0(
        "The Stage 5 checkpoint is incompatible with the current inputs or settings: ",
        checkpoint_dir,
        "\nUse mode='restart' to discard only this partial state."
      )
    )
    ensure_writable_dir(stage5_completed_dir(config), "Stage 5 completed-species directory")
    return(invisible(FALSE))
  }

  ensure_writable_dir(dirname(checkpoint_dir), "Stage 5 checkpoint parent")
  staging <- tempfile("stage5_checkpoint_", tmpdir = dirname(checkpoint_dir))
  assert(dir.create(staging), "Could not create the Stage 5 checkpoint staging directory.")
  on.exit(unlink(staging, recursive = TRUE, force = TRUE), add = TRUE)
  dir.create(file.path(staging, "completed"))
  saveRDS(manifest, file.path(staging, "manifest.rds"), version = 3)
  assert(
    identical(readRDS(file.path(staging, "manifest.rds")), manifest),
    "Temporary Stage 5 checkpoint manifest failed validation."
  )
  assert(file.rename(staging, checkpoint_dir), "Could not install the Stage 5 checkpoint.")
  invisible(TRUE)
}

validate_stage5_checkpoint_result <- function(result, directory, template = NULL) {
  assert(is.list(result), "Stage 5 checkpoint result must be a list.")
  assert(
    is.character(result$species) && length(result$species) == 1L &&
      !is.na(result$species) && nzchar(result$species),
    "Stage 5 checkpoint result has an invalid species name."
  )
  assert(
    is.character(result$status) && length(result$status) == 1L &&
      !is.na(result$status) && result$status != "error",
    "Stage 5 checkpoint result must be a completed scientific outcome."
  )

  if (!identical(result$status, "retained")) return(invisible(TRUE))

  validate_patch_lookup_object(
    result$patch_lookup,
    paste0("Stage 5 checkpoint lookup for ", result$species)
  )
  connectivity <- set_patch_connectivity_attrs(result$connectivity)
  validate_patch_output_consistency(
    result$patch_lookup,
    connectivity,
    paste0("Stage 5 checkpoint for ", result$species)
  )
  raster_path <- file.path(directory, basename(result$patch_raster))
  need_file(raster_path, paste0("Stage 5 checkpoint raster for ", result$species))
  raster <- terra::rast(raster_path)
  if (!is.null(template)) {
    assert(
      terra::compareGeom(raster, template, stopOnError = FALSE),
      paste0("Stage 5 checkpoint raster geometry changed for ", result$species, ".")
    )
  }
  raster_ids <- as.numeric(terra::unique(raster, na.rm = TRUE)[[1L]])
  raster_ids <- sort(unique(as.integer(raster_ids[is.finite(raster_ids)])))
  expected_ids <- sort(unique(as.integer(result$patch_lookup$patch_id)))
  assert(
    identical(raster_ids, expected_ids),
    paste0("Stage 5 checkpoint raster IDs differ from the lookup for ", result$species, ".")
  )
  invisible(TRUE)
}

write_stage5_species_checkpoint <- function(result, config, index, template,
                                            rename_file = file.rename) {
  assert(is.function(rename_file), "rename_file must be a function.")
  destination <- stage5_species_checkpoint_dir(config, index, result$species)
  assert(!dir.exists(destination), paste0("Stage 5 species checkpoint already exists: ", destination))
  parent <- stage5_completed_dir(config)
  ensure_writable_dir(parent, "Stage 5 completed-species directory")
  staging <- tempfile(
    paste0(".", stage5_species_checkpoint_name(index, result$species), "_"),
    tmpdir = parent
  )
  assert(dir.create(staging), "Could not create a Stage 5 species checkpoint staging directory.")
  on.exit(unlink(staging, recursive = TRUE, force = TRUE), add = TRUE)

  stored <- result
  if (identical(stored$status, "retained")) {
    source_raster <- stored$patch_raster
    target_raster <- file.path(staging, basename(source_raster))
    assert(
      file.copy(source_raster, target_raster, overwrite = FALSE, copy.mode = TRUE),
      paste0("Could not stage the patch raster for ", stored$species, ".")
    )
    stored$patch_raster <- basename(target_raster)
  }
  saveRDS(stored, file.path(staging, "result.rds"), version = 3)
  check <- readRDS(file.path(staging, "result.rds"))
  validate_stage5_checkpoint_result(check, staging, template)
  assert(
    rename_file(staging, destination),
    paste0("Could not atomically install the checkpoint for ", stored$species, ".")
  )
  need_file(
    file.path(destination, "result.rds"),
    paste0("installed Stage 5 checkpoint result for ", stored$species)
  )
  if (identical(check$status, "retained")) {
    installed_raster <- file.path(destination, basename(check$patch_raster))
    need_file(installed_raster, paste0(
      "installed Stage 5 checkpoint raster for ", stored$species
    ))
    check$patch_raster <- installed_raster
  }
  invisible(check)
}

read_stage5_species_checkpoint <- function(config, index, species, template = NULL) {
  directory <- stage5_species_checkpoint_dir(config, index, species)
  if (!dir.exists(directory)) return(NULL)
  path <- file.path(directory, "result.rds")
  need_file(path, paste0("Stage 5 checkpoint result for ", species))
  result <- readRDS(path)
  assert(
    identical(result$species, species),
    paste0("Stage 5 checkpoint species does not match its canonical slot: ", directory)
  )
  validate_stage5_checkpoint_result(result, directory, template)
  if (identical(result$status, "retained")) {
    result$patch_raster <- file.path(directory, basename(result$patch_raster))
  }
  result
}

stage5_final_artifacts_present <- function(config) {
  c(
    patch_dir = dir.exists(config$paths$patch_dir),
    patch_lookup = file.exists(config$paths$patch_lookup_rds),
    connectivity = file.exists(config$paths$connectivity_rds),
    metadata = file.exists(config$paths$metadata_rds)
  )
}

validate_stage5_final_outputs <- function(config) {
  present <- stage5_final_artifacts_present(config)
  required <- present[c("patch_dir", "patch_lookup", "connectivity")]
  if (!all(required)) {
    assert(
      !any(required),
      "Stage 5 final artifacts are incomplete; patch directory, lookup, and connectivity must exist together."
    )
    return(NULL)
  }
  patch_lookup <- readRDS(config$paths$patch_lookup_rds)
  connectivity <- readRDS(config$paths$connectivity_rds)
  validate_patch_output_consistency(patch_lookup, connectivity, "Stage 5 final outputs")
  species <- sort(unique(patch_squish(patch_lookup$scientificName)))
  expected <- file.path(
    config$paths$patch_dir,
    vapply(species, patch_filename_from_scientific, character(1))
  )
  missing <- expected[!file.exists(expected)]
  assert(!length(missing), paste0("Stage 5 final patch raster(s) are missing:\n", paste(missing, collapse = "\n")))
  extras <- setdiff(
    normalizePath(
      list.files(config$paths$patch_dir, pattern = "\\.tif$", full.names = TRUE),
      winslash = "/",
      mustWork = FALSE
    ),
    normalizePath(expected, winslash = "/", mustWork = FALSE)
  )
  assert(!length(extras), paste0("Stage 5 final patch directory contains orphan raster(s):\n", paste(extras, collapse = "\n")))
  list(
    patch_lookup = patch_lookup,
    connectivity = connectivity,
    metadata = if (present[["metadata"]]) readRDS(config$paths$metadata_rds) else NULL
  )
}

promote_stage5_outputs <- function(config, manifest, results) {
  retained <- results[vapply(results, function(x) identical(x$status, "retained"), logical(1))]
  assert(length(retained) > 0L, "No species were retained; Stage 5 final outputs were not replaced.")
  patch_lookup <- dplyr::bind_rows(lapply(retained, `[[`, "patch_lookup"))
  connectivity <- unlist(lapply(retained, `[[`, "connectivity"), recursive = FALSE)
  connectivity <- set_patch_connectivity_attrs(connectivity)
  validate_patch_output_consistency(patch_lookup, connectivity, "assembled Stage 5 outputs")

  clean_dir <- dirname(config$paths$patch_lookup_rds)
  staging_root <- tempfile("stage5_final_", tmpdir = clean_dir)
  assert(dir.create(staging_root), "Could not create Stage 5 final staging directory.")
  on.exit(unlink(staging_root, recursive = TRUE, force = TRUE), add = TRUE)
  staged_patch_dir <- file.path(staging_root, "Patches")
  dir.create(staged_patch_dir)

  for (result in retained) {
    target <- file.path(staged_patch_dir, patch_filename_from_scientific(result$species))
    assert(
      file.copy(result$patch_raster, target, overwrite = FALSE, copy.mode = TRUE),
      paste0("Could not stage final patch raster for ", result$species, ".")
    )
  }
  staged_lookup <- file.path(staging_root, "all_patch_lookup.rds")
  staged_connectivity <- file.path(staging_root, "all_connectivity.rds")
  staged_metadata <- file.path(staging_root, "stage5_build_metadata.rds")
  metadata <- list(
    schema = 1L,
    manifest = manifest,
    completed_species = vapply(results, `[[`, character(1), "species"),
    retained_species = vapply(retained, `[[`, character(1), "species"),
    status = patch_status_summary(results)
  )
  saveRDS(patch_lookup, staged_lookup, version = 3)
  saveRDS(connectivity, staged_connectivity, version = 3)
  saveRDS(metadata, staged_metadata, version = 3)
  validate_patch_output_consistency(
    readRDS(staged_lookup),
    readRDS(staged_connectivity),
    "staged Stage 5 outputs"
  )

  targets <- c(
    config$paths$patch_dir,
    config$paths$patch_lookup_rds,
    config$paths$connectivity_rds,
    config$paths$metadata_rds
  )
  staged <- c(staged_patch_dir, staged_lookup, staged_connectivity, staged_metadata)
  backups <- paste0(targets, ".stage5_backup")
  assert(!any(file.exists(backups) | dir.exists(backups)), "A stale Stage 5 backup blocks final promotion.")
  existed <- file.exists(targets) | dir.exists(targets)

  committed <- logical(length(targets))
  tryCatch({
    for (i in seq_along(targets)) {
      if (existed[[i]]) {
        assert(file.rename(targets[[i]], backups[[i]]), paste0("Could not back up ", targets[[i]], "."))
      }
      assert(file.rename(staged[[i]], targets[[i]]), paste0("Could not install ", targets[[i]], "."))
      committed[[i]] <- TRUE
    }
  }, error = function(e) {
    for (i in rev(seq_along(targets))) {
      if (committed[[i]] && (file.exists(targets[[i]]) || dir.exists(targets[[i]]))) {
        unlink(targets[[i]], recursive = TRUE, force = TRUE)
      }
      if (existed[[i]] && (file.exists(backups[[i]]) || dir.exists(backups[[i]]))) {
        file.rename(backups[[i]], targets[[i]])
      }
    }
    stop(e)
  })
  unlink(backups[existed], recursive = TRUE, force = TRUE)
  validate_stage5_final_outputs(config)
  unlink(config$paths$checkpoint_dir, recursive = TRUE, force = TRUE)
  invisible(list(
    patch_lookup = patch_lookup,
    connectivity = connectivity,
    metadata = metadata
  ))
}
