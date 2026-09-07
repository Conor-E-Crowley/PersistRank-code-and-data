# Filesystem artifact helpers loaded after project_utils.R by numbered workflows.
#
# This definition-only module owns file/directory assertions, writable output
# preparation, read-only artifact status, and cache-freshness validation. Calls
# may inspect or prepare explicitly supplied paths; sourcing performs no I/O,
# package loading, compilation, or expensive runtime initialization.

need_file <- function(path, label = path) {
  assert(file.exists(path), paste0("Missing file: ", label, "\nPath: ", path))
}

need_dir <- function(path, label = path) {
  assert(dir.exists(path), paste0("Missing directory: ", label, "\nPath: ", path))
}

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  need_dir(path)
  invisible(path)
}

ensure_writable_dir <- function(path, label = path) {
  ensure_dir(path)
  probe <- tempfile("write_test_", tmpdir = path, fileext = ".tmp")
  writable <- tryCatch(file.create(probe), error = function(e) FALSE)
  if (isTRUE(writable)) unlink(probe)
  assert(writable, paste0("Directory is not writable: ", label, "\nPath: ", path))
  invisible(path)
}

# Return lightweight metadata without opening or deserializing an artifact.
# Directories deliberately have no byte size; missing paths remain in output.
project_artifact_status <- function(paths) {
  paths <- unname(as.character(paths))
  exists <- file.exists(paths) | dir.exists(paths)
  info <- file.info(paths)
  data.frame(
    path = paths,
    exists = exists,
    type = ifelse(
      !exists, "missing", ifelse(dir.exists(paths), "directory", "file")
    ),
    size_bytes = ifelse(exists & !dir.exists(paths), info$size, NA_real_),
    modified = ifelse(
      exists, format(info$mtime, "%Y-%m-%d %H:%M:%S"), NA_character_
    ),
    stringsAsFactors = FALSE
  )
}

# Validate that a cached artifact exists and is no older than every dependency.
# The operation is read-only and returns TRUE invisibly; missing/unreadable paths
# and stale caches fail with path-specific regeneration guidance.
validate_cached_artifact_freshness <- function(
  cached_path,
  dependency_paths,
  label = "cached artifact",
  rebuild_hint = "Regenerate the cached artifact before reusing it."
) {
  cached_path <- validate_path_param(cached_path, "cached_path")
  label <- validate_scalar_string(label, "label")
  rebuild_hint <- validate_scalar_string(rebuild_hint, "rebuild_hint")

  assert(
    is.character(dependency_paths),
    "dependency_paths must be a character vector of source paths."
  )

  dependency_paths <- trimws(dependency_paths)
  dependency_paths <- unique(path.expand(
    dependency_paths[!is.na(dependency_paths) & nzchar(dependency_paths)]
  ))

  need_file(cached_path, label)

  if (!length(dependency_paths)) {
    return(invisible(TRUE))
  }

  preview_paths <- function(paths, max_paths = 8L) {
    paths <- normalizePath(paths, winslash = "/", mustWork = FALSE)
    if (length(paths) <= max_paths) {
      return(paths)
    }
    c(
      head(paths, max_paths),
      paste0("... and ", length(paths) - max_paths, " more")
    )
  }

  missing_dependencies <- dependency_paths[!file.exists(dependency_paths)]
  assert(
    !length(missing_dependencies),
    paste0(
      "Cannot validate freshness for ", label,
      " because source artifact(s) are missing:\n",
      paste(preview_paths(missing_dependencies), collapse = "\n")
    )
  )

  cached_info <- file.info(cached_path)
  assert(
    !is.na(cached_info$mtime[[1L]]),
    paste0("Cannot read modification time for ", label, ": ", cached_path)
  )

  dependency_info <- file.info(dependency_paths)
  missing_mtime <- is.na(dependency_info$mtime)
  assert(
    !any(missing_mtime),
    paste0(
      "Cannot read modification time for source artifact(s) used by ",
      label,
      ":\n",
      paste(preview_paths(dependency_paths[missing_mtime]), collapse = "\n")
    )
  )

  newer_dependencies <- dependency_paths[
    dependency_info$mtime > cached_info$mtime[[1L]]
  ]

  assert(
    !length(newer_dependencies),
    paste0(
      "Cached ", label, " is older than source artifact(s):\n",
      paste(preview_paths(newer_dependencies), collapse = "\n"),
      "\n",
      rebuild_hint
    )
  )

  invisible(TRUE)
}
