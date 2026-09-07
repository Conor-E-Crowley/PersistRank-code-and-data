# Exact retained-cell benchmark lookup libraries shared by Stages 7.2, 7.4,
# and 8.1.
#
# Inputs: one canonical benchmark configuration and a target plan derived from
# completed Stage 6 states. Returns: regenerated status/index tables or validated
# lookup objects. Reads: rank map, the shared Stage 6 initialization, manifest,
# and exact lookup files. Writes: one immutable library manifest and missing
# exact lookup files. Cost: cheap except when called by reconstruction. Side
# effects are limited to atomic file commits; compatible exact targets are
# never overwritten.


benchmark_library_schema <- function() 2L
benchmark_lookup_schema <- function() 2L
benchmark_reconstruction_schema <- function() 2L

benchmark_library_from_config <- function(config) {
  library <- config$paths$lookup_library
  assert(is.list(library), "Configuration does not define a benchmark lookup library.")
  library
}

benchmark_manifest_fingerprint <- function(path, label, root) {
  record <- storage_file_fingerprint(path, label, root = root)
  list(
    path = record$relative_path %||% record$path,
    bytes = record$bytes,
    md5 = record$md5
  )
}

expected_benchmark_library_identity <- function(config) {
  list(
    schema = benchmark_library_schema(),
    benchmark_method = config$rank_method,
    application = config$project_paths$application,
    threshold = unclass(config$contract),
    removal_run = list(
      run_tag = config$project_paths$run_tag,
      cells_to_remove_per_iteration = config$cells_to_remove_per_iteration,
      pruning_iterations_per_stage = config$pruning_iterations_per_stage
    ),
    reconstruction_schema = benchmark_reconstruction_schema()
  )
}

expected_benchmark_library_manifest <- function(config) {
  identity <- expected_benchmark_library_identity(config)
  list(
    schema = identity$schema,
    benchmark_method = identity$benchmark_method,
    rank_map = benchmark_manifest_fingerprint(
      config$paths$zonation_rankmap, paste0(config$rank_label, " rank map"),
      config$project_paths$root
    ),
    application = identity$application,
    threshold = identity$threshold,
    removal_run = identity$removal_run,
    spatial_source = benchmark_manifest_fingerprint(
      config$project_paths$stage5_metadata, "Stage 5 spatial source",
      config$project_paths$root
    ),
    reconstruction_schema = identity$reconstruction_schema
  )
}

# Report-time validation deliberately ignores reconstruction source files. The
# exact lookups remain usable after transfer to a machine without rank rasters.
validate_benchmark_library_manifest <- function(config, verify_sources = FALSE) {
  library <- benchmark_library_from_config(config)
  need_file(library$manifest, paste0(config$rank_label, " lookup manifest"))
  observed <- readRDS(library$manifest)
  required_names <- c(
    "schema", "benchmark_method", "rank_map", "application", "threshold",
    "removal_run", "spatial_source", "reconstruction_schema"
  )
  assert(identical(names(observed), required_names), paste0(
    config$rank_label, " lookup manifest has an unsupported schema."
  ))
  identity <- observed[c(
    "schema", "benchmark_method", "application", "threshold",
    "removal_run", "reconstruction_schema"
  )]
  assert(identical(identity, expected_benchmark_library_identity(config)), paste0(
    config$rank_label,
    " lookup library is incompatible with the scientific run identity. ",
    "Rebuild only this method with Stage 7.2."
  ))
  if (isTRUE(verify_sources)) {
    expected <- expected_benchmark_library_manifest(config)
    assert(identical(observed, expected), paste0(
      config$rank_label,
      " lookup library is incompatible with the rank map or spatial source. ",
      "Rebuild only this method with Stage 7.2."
    ))
  }
  invisible(observed)
}

initialize_benchmark_library <- function(config) {
  library <- benchmark_library_from_config(config)
  expected <- expected_benchmark_library_manifest(config)
  if (file.exists(library$manifest)) {
    validate_benchmark_library_manifest(config, verify_sources = TRUE)
    return(invisible(expected))
  }
  existing <- if (dir.exists(library$lookups)) list.files(
    library$lookups, pattern = "^retained_cells_[0-9]+[.]rds$", full.names = TRUE
  ) else character()
  assert(!length(existing), paste0(
    "Exact lookup files exist without a library manifest: ", library$lookups,
    ". Rebuild this method with Stage 7.2."
  ))
  ensure_writable_dir(library$root, "benchmark lookup-library directory")
  storage_atomic_save_rds(expected, library$manifest, "benchmark lookup manifest")
  invisible(expected)
}

benchmark_lookup_object <- function(retained_cells, patch_table) {
  list(
    schema = benchmark_lookup_schema(),
    retained_cells = as.integer(retained_cells),
    patch_table = data.table::copy(data.table::as.data.table(patch_table))
  )
}

validate_benchmark_lookup <- function(object, retained_cells,
                                             label = "benchmark lookup") {
  assert(is.list(object) && identical(object$schema, benchmark_lookup_schema()),
         paste0(label, " has an unsupported schema."))
  assert(identical(as.integer(object$retained_cells), as.integer(retained_cells)),
         paste0(label, " belongs to a different retained-cell target."))
  patches <- data.table::as.data.table(object$patch_table)
  need_cols(patches, c("species", "patch_id", "pu_id", "patch_area_km2"), label)
  assert(!anyDuplicated(patches[, .(species, patch_id)]) &&
           all(patches$patch_id > 0L) && all(patches$pu_id > 0L) &&
           all(is.finite(patches$patch_area_km2) & patches$patch_area_km2 > 0),
         paste0(label, " contains invalid patch/PU state."))
  invisible(object)
}

read_benchmark_lookup <- function(config, retained_cells,
                                  validate_manifest = TRUE) {
  if (isTRUE(validate_manifest)) validate_benchmark_library_manifest(config)
  path <- benchmark_lookup_file(
    benchmark_library_from_config(config), retained_cells
  )
  need_file(path, paste0(config$rank_label, " retained-cell lookup ", retained_cells))
  object <- readRDS(path)
  validate_benchmark_lookup(object, retained_cells, basename(path))
  object
}

write_benchmark_lookup <- function(config, retained_cells, patch_table) {
  library <- benchmark_library_from_config(config)
  path <- benchmark_lookup_file(library, retained_cells)
  if (file.exists(path)) {
    validate_benchmark_lookup(readRDS(path), retained_cells, basename(path))
    return(list(path = path, status = "reused"))
  }
  ensure_writable_dir(library$lookups, "benchmark exact-lookup directory")
  object <- benchmark_lookup_object(retained_cells, patch_table)
  staged <- tempfile("lookup_", tmpdir = library$lookups, fileext = ".rds")
  on.exit(unlink(staged, force = TRUE), add = TRUE)
  saveRDS(object, staged, version = 3L)
  validate_benchmark_lookup(readRDS(staged), retained_cells, basename(path))
  project_file_set_transaction(
    staged, path, overwrite = FALSE,
    label = paste0("benchmark lookup retained_cells=", retained_cells)
  )
  list(path = path, status = "created")
}

benchmark_library_index <- function(config, plan, validate = TRUE) {
  library <- benchmark_library_from_config(config)
  targets <- as.integer(plan$union$keep_n)
  paths <- benchmark_lookup_file(library, targets)
  rows <- data.table::data.table(
    retained_cells = targets,
    lookup_file = paths,
    exists = file.exists(paths)
  )
  rows[, `:=`(rows = NA_integer_, bytes = NA_real_, modified = as.POSIXct(NA))]
  for (i in which(rows$exists)) {
    if (isTRUE(validate)) {
      object <- readRDS(paths[[i]])
      validate_benchmark_lookup(object, targets[[i]], basename(paths[[i]]))
      rows[i, rows := nrow(object$patch_table)]
    }
    info <- file.info(paths[[i]])
    rows[i, `:=`(
      bytes = as.numeric(info$size), modified = info$mtime
    )]
  }
  rows[]
}

benchmark_library_complete <- function(config, plan, validate_lookups = FALSE) {
  library <- benchmark_library_from_config(config)
  if (!file.exists(library$manifest)) return(FALSE)
  validate_benchmark_library_manifest(config)
  index <- benchmark_library_index(config, plan, validate = validate_lookups)
  all(index$exists)
}

inspect_benchmark_library <- function(config, plan) {
  library <- benchmark_library_from_config(config)
  required <- nrow(plan$union)
  paths <- benchmark_lookup_file(library, plan$union$keep_n)
  manifest_status <- if (!file.exists(library$manifest)) {
    "missing"
  } else if (inherits(try(validate_benchmark_library_manifest(config), silent = TRUE),
                         "try-error")) {
    "incompatible"
  } else {
    "compatible"
  }
  data.frame(
    benchmark = config$rank_method,
    benchmark_label = config$rank_label,
    required = required,
    existing = sum(file.exists(paths)),
    missing = sum(!file.exists(paths)),
    checkpoint = file.exists(library$checkpoint),
    compatibility = manifest_status,
    stringsAsFactors = FALSE
  )
}

finalize_benchmark_library <- function(config, plan) {
  validate_benchmark_library_manifest(config)
  index <- benchmark_library_index(config, plan, validate = TRUE)
  assert(all(index$exists), "Cannot finalize an incomplete benchmark lookup library.")
  library <- benchmark_library_from_config(config)
  checkpoints <- c(library$checkpoint, library$checkpoint_backup)
  for (checkpoint in checkpoints[file.exists(checkpoints)]) {
    unlink(checkpoint, force = TRUE)
    assert(!file.exists(checkpoint),
           "Could not remove a completed reconstruction checkpoint.")
  }
  invisible(index)
}

# Stage 7.2 alone may deliberately delete one method library.
restart_benchmark_library <- function(config) {
  library <- benchmark_library_from_config(config)
  expected_parent <- normalizePath(
    config$project_paths$benchmark_lookups, winslash = "/", mustWork = FALSE
  )
  target <- normalizePath(library$root, winslash = "/", mustWork = FALSE)
  assert(dirname(target) == expected_parent &&
           basename(target) == toupper(config$rank_method),
         "Refusing to reset an unexpected benchmark lookup directory.")
  if (dir.exists(target)) unlink(target, recursive = TRUE, force = TRUE)
  assert(!dir.exists(target), paste0("Could not reset benchmark library: ", target))
  invisible(target)
}
