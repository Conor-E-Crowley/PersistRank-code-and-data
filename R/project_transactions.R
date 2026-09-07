# Atomic filesystem transactions loaded after project_artifacts.R.
#
# Workflows use these helpers to promote validated staged files or replace a
# directory without exposing partial output. Sourcing is definition-only. Calls
# may move/remove only their explicit staged, target, backup, or replacement
# paths and retain byte-identical rollback validation.

# Commit a validated set of staged files as one unit. Existing targets are
# restored byte-for-byte if any backup or replacement fails.
project_file_set_transaction <- function(staged_paths,
                                         target_paths,
                                         overwrite = FALSE,
                                         rename_file = file.rename,
                                         defer_finalize = FALSE,
                                         label = "Project") {
  staged_paths <- as.character(staged_paths)
  target_paths <- as.character(target_paths)
  assert(
    length(staged_paths) > 0L && length(staged_paths) == length(target_paths),
    paste0(label, " file transaction requires equal, nonempty staged and target path sets.")
  )
  assert(
    !anyDuplicated(staged_paths) && !anyDuplicated(target_paths),
    paste0(label, " file transaction paths must be unique.")
  )

  missing_staged <- staged_paths[!file.exists(staged_paths)]
  assert(
    !length(missing_staged),
    paste0(label, " file transaction has missing staged file(s):\n", paste(missing_staged, collapse = "\n"))
  )
  conflicts <- target_paths[file.exists(target_paths)]
  assert(
    isTRUE(overwrite) || !length(conflicts),
    paste0(
      label, " file transaction would replace existing file(s) without authorization:\n",
      paste(conflicts, collapse = "\n")
    )
  )
  invisible(lapply(unique(dirname(target_paths)), ensure_writable_dir, label = paste0(label, " output directory")))

  existed <- file.exists(target_paths)
  original_md5 <- rep(NA_character_, length(target_paths))
  original_md5[existed] <- unname(tools::md5sum(target_paths[existed]))
  backups <- vapply(seq_along(target_paths), function(i) {
    tempfile(
      "project_set_backup_",
      tmpdir = dirname(target_paths[[i]]),
      fileext = paste0(".", tools::file_ext(target_paths[[i]]))
    )
  }, character(1L))
  backed_up <- logical(length(target_paths))
  committed <- logical(length(target_paths))
  active <- TRUE

  rollback <- function() {
    if (!active) return(invisible(FALSE))
    unlink(target_paths[committed], force = TRUE)
    restore_failed <- character(0L)
    for (i in which(backed_up)) {
      if (!file.exists(backups[[i]]) || !isTRUE(rename_file(backups[[i]], target_paths[[i]]))) {
        restore_failed <- c(restore_failed, target_paths[[i]])
      }
    }
    if (length(restore_failed)) {
      stop(
        label, " rollback could not restore file(s):\n",
        paste(restore_failed, collapse = "\n"),
        call. = FALSE
      )
    }
    restored <- which(existed)
    if (length(restored)) {
      restored_md5 <- unname(tools::md5sum(target_paths[restored]))
      if (!identical(restored_md5, original_md5[restored])) {
        stop(label, " rollback restored files with changed checksums.", call. = FALSE)
      }
    }
    unlink(c(staged_paths, backups), force = TRUE)
    active <<- FALSE
    invisible(FALSE)
  }
  finalize <- function() {
    if (active) unlink(backups[backed_up], force = TRUE)
    unlink(staged_paths[file.exists(staged_paths)], force = TRUE)
    active <<- FALSE
    invisible(TRUE)
  }

  failure <- NULL
  tryCatch({
    for (i in which(existed)) {
      if (!isTRUE(rename_file(target_paths[[i]], backups[[i]]))) {
        stop("Could not back up ", label, " output: ", target_paths[[i]], call. = FALSE)
      }
      backed_up[[i]] <- TRUE
    }
    for (i in seq_along(target_paths)) {
      if (!isTRUE(rename_file(staged_paths[[i]], target_paths[[i]]))) {
        stop("Could not commit ", label, " output: ", target_paths[[i]], call. = FALSE)
      }
      committed[[i]] <- TRUE
    }
  }, error = function(e) failure <<- e)
  if (!is.null(failure)) {
    rollback()
    stop(failure)
  }

  transaction <- list(
    targets = target_paths,
    backups = backups,
    rollback = rollback,
    finalize = finalize
  )
  if (isTRUE(defer_finalize)) return(transaction)
  transaction$finalize()
  invisible(target_paths)
}

# Return the fixed sibling path used while a directory replacement is active.
directory_replacement_backup_path <- function(path) {
  paste0(validate_path_param(path, "replacement directory"), ".restart_backup")
}

# Move an existing directory aside until the caller finalizes or rolls back.
# The returned closures are idempotent; rollback removes a partial replacement
# and restores the original, while finalize removes only the explicit backup.
begin_directory_replacement <- function(
  path,
  label = "Directory",
  rename_path = file.rename,
  remove_path = function(x) {
    unlink(x, recursive = TRUE, force = TRUE)
    !file.exists(x)
  }
) {
  path <- validate_path_param(path, paste0(label, " path"))
  backup <- directory_replacement_backup_path(path)
  assert(
    !file.exists(backup),
    paste0("Unresolved ", label, " replacement backup: ", backup)
  )

  had_previous <- file.exists(path)
  if (had_previous) {
    assert(
      isTRUE(rename_path(path, backup)),
      paste0("Could not back up existing ", label, ": ", path)
    )
  }
  active <- TRUE

  rollback <- function() {
    if (!active) return(invisible(FALSE))
    if (file.exists(path)) {
      assert(
        isTRUE(remove_path(path)) && !file.exists(path),
        paste0("Could not remove partial ", label, " during rollback: ", path)
      )
    }
    if (had_previous) {
      assert(
        file.exists(backup) && isTRUE(rename_path(backup, path)),
        paste0("Could not restore previous ", label, ". Backup retained at: ", backup)
      )
    }
    active <<- FALSE
    invisible(FALSE)
  }

  finalize <- function() {
    if (!active) return(invisible(TRUE))
    if (file.exists(backup)) {
      assert(
        isTRUE(remove_path(backup)) && !file.exists(backup),
        paste0("Could not remove completed ", label, " backup: ", backup)
      )
    }
    active <<- FALSE
    invisible(TRUE)
  }

  list(
    path = path,
    backup = backup,
    had_previous = had_previous,
    rollback = rollback,
    finalize = finalize
  )
}
