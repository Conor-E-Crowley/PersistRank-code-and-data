# Deterministic source loader for repository R modules.
#
# Workflows own dependency order and call project_source() with repository-
# relative files. Definitions are installed in the global environment to retain
# the existing Rmd interface. The private registry stores paths and loading
# states only; sourcing this file again deliberately preserves that registry.

if (!exists(".project_loader_state", envir = .GlobalEnv, inherits = FALSE) ||
    !is.environment(get(".project_loader_state", envir = .GlobalEnv))) {
  .project_loader_state <- new.env(parent = emptyenv())
  .project_loader_state$loaded <- new.env(parent = emptyenv())
  .project_loader_state$loading <- new.env(parent = emptyenv())
  .project_loader_state$stack <- character()
}

# Source repository modules once per R session.
#
# `files` is an ordered character vector. Paths are normalized before registry
# lookup, so equivalent spellings identify the same module. `reload = TRUE`
# intentionally replaces definitions for interactive development. A module is
# marked loaded only after sys.source() succeeds, and recursive dependency
# cycles report the complete active chain.
project_source <- function(files, reload = FALSE) {
  if (!is.character(files)) {
    stop("files must be a character vector of repository R-module paths.", call. = FALSE)
  }
  if (!is.logical(reload) || length(reload) != 1L || is.na(reload)) {
    stop("reload must be TRUE or FALSE.", call. = FALSE)
  }
  if (!length(files)) return(invisible(character()))

  invalid <- which(is.na(files) | !nzchar(trimws(files)))
  if (length(invalid)) {
    stop(
      "Invalid source path at files[[", invalid[[1L]], "]]: paths must be non-empty and non-missing.",
      call. = FALSE
    )
  }
  files <- unique(files)
  for (path in files) {
    if (!file.exists(path)) {
      stop("Missing source file: ", path, call. = FALSE)
    }
    if (dir.exists(path)) {
      stop("Source path is a directory, not a file: ", path, call. = FALSE)
    }
  }
  normalized <- vapply(
    files,
    normalizePath,
    character(1L),
    winslash = "/",
    mustWork = TRUE,
    USE.NAMES = FALSE
  )
  normalized <- unique(normalized)

  state <- .project_loader_state
  for (path in normalized) {
    if (exists(path, envir = state$loading, inherits = FALSE)) {
      start <- match(path, state$stack)
      cycle <- c(state$stack[seq.int(start, length(state$stack))], path)
      stop(
        "Circular source dependency: ", paste(cycle, collapse = " -> "),
        call. = FALSE
      )
    }
    if (!reload && exists(path, envir = state$loaded, inherits = FALSE)) next

    assign(path, TRUE, envir = state$loading)
    state$stack <- c(state$stack, path)
    completed <- FALSE
    tryCatch(
      {
        sys.source(path, envir = .GlobalEnv)
        completed <- TRUE
      },
      finally = {
        rm(list = path, envir = state$loading)
        state$stack <- head(state$stack, -1L)
      }
    )
    if (completed) assign(path, TRUE, envir = state$loaded)
  }
  invisible(normalized)
}
