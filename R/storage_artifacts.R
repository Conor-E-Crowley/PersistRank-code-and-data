# Relocatable artifact records and atomic metadata persistence.
#
# Loaded after project artifact and transaction infrastructure by workflows that
# validate or write model/application metadata. Sourcing defines functions only;
# checksums, filesystem reads, and writes occur solely on explicit calls.

storage_relative_path <- function(path, root) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  root <- normalizePath(root, winslash = "/", mustWork = FALSE)
  prefix <- paste0(root, "/")
  inside <- startsWith(path, prefix)
  path[inside] <- substring(path[inside], nchar(prefix) + 1L)
  path
}

storage_file_fingerprint <- function(path, label = "input file", root = NULL) {
  path <- normalizePath(validate_path_param(path, label), winslash = "/", mustWork = TRUE)
  info <- file.info(path)
  assert(!isTRUE(info$isdir), paste0(label, " must be a file: ", path))
  md5 <- unname(tools::md5sum(path))
  assert(!is.na(md5), paste0("Cannot checksum ", label, ": ", path))
  record <- list(path = path, bytes = unname(as.numeric(info$size)), md5 = md5)
  if (!is.null(root)) {
    relative <- storage_relative_path(path, root)
    if (!identical(relative, path)) record$relative_path <- relative
  }
  record
}

# Resolve a canonical project-relative record or an explicitly external path.
storage_record_path <- function(record, root, must_work = FALSE) {
  assert(is.list(record) && is.character(record$path) && length(record$path) == 1L,
         "Artifact record must contain one path.")
  root <- normalizePath(root, winslash = "/", mustWork = FALSE)
  resolved <- if (!is.null(record$relative_path)) {
    file.path(root, record$relative_path)
  } else if (grepl("^/|^[A-Za-z]:[/\\\\]", record$path)) {
    record$path
  } else {
    file.path(root, record$path)
  }
  resolved <- normalizePath(resolved, winslash = "/", mustWork = FALSE)
  if (isTRUE(must_work)) {
    assert(file.exists(resolved) || dir.exists(resolved), paste0(
      "Recorded artifact is missing: ", resolved
    ))
  }
  resolved
}

storage_fingerprint_equal <- function(x, y) {
  identical(x$bytes, y$bytes) && identical(x$md5, y$md5)
}

storage_fingerprint_lists_equal <- function(x, y) {
  identical(names(x), names(y)) && length(x) == length(y) &&
    all(Map(storage_fingerprint_equal, x, y) |> unlist(use.names = FALSE))
}

storage_atomic_save_rds <- function(object, path, label = "project metadata") {
  ensure_writable_dir(dirname(path), paste0(label, " directory"))
  staged <- tempfile("metadata_", tmpdir = dirname(path), fileext = ".rds")
  on.exit(unlink(staged, force = TRUE), add = TRUE)
  saveRDS(object, staged, version = 3L)
  assert(identical(readRDS(staged), object), paste0("Failed to validate staged ", label, "."))
  project_file_set_transaction(
    staged, path, overwrite = file.exists(path), label = label
  )
  invisible(path)
}
