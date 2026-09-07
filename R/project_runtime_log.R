# Append-only operational runtime logging loaded after project_utils.R.
#
# Workflows call these helpers only around active scientific work. Sourcing
# attaches no package and writes no file; it preserves an existing session state
# environment so an explicit development reload cannot silently orphan a log.

format_runtime_log_value <- function(x) {
  if (is.null(x)) {
    return("NULL")
  }

  if (!length(x)) {
    return("empty")
  }

  value <- paste(as.character(x), collapse = ",")
  if (is.na(value)) {
    value <- "NA"
  }

  quote_value <- grepl("[[:space:]]", value)
  value <- gsub("\\", "\\\\", value, fixed = TRUE)
  value <- gsub("\r", "\\r", value, fixed = TRUE)
  value <- gsub("\n", "\\n", value, fixed = TRUE)
  value <- gsub("\t", "\\t", value, fixed = TRUE)
  value <- gsub('"', '\\"', value, fixed = TRUE)
  if (quote_value) {
    value <- paste0('"', value, '"')
  }

  value
}

# Re-sourcing deterministic workflow dependencies must replace function
# definitions without discarding a session already opened by the caller.
if (!exists(".runtime_log_state", inherits = FALSE)) {
  .runtime_log_state <- new.env(parent = emptyenv())
  .runtime_log_state$path <- NULL
  .runtime_log_state$session_id <- NULL
  .runtime_log_state$started <- NULL
  .runtime_log_state$serial <- 0L
}
if (is.null(.runtime_log_state$serial)) .runtime_log_state$serial <- 0L

runtime_log_active <- function() {
  path <- .runtime_log_state$path
  is.character(path) && length(path) == 1L && !is.na(path) && nzchar(path)
}

runtime_log_path <- function() {
  if (!runtime_log_active()) return(NULL)
  .runtime_log_state$path
}

clear_runtime_log <- function() {
  .runtime_log_state$path <- NULL
  .runtime_log_state$session_id <- NULL
  .runtime_log_state$started <- NULL
  invisible(NULL)
}

runtime_log_timestamp <- function(time = Sys.time()) {
  format(time, "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC")
}

runtime_log_session_id <- function(time = Sys.time()) {
  .runtime_log_state$serial <- .runtime_log_state$serial + 1L
  stamp <- format(time, "%Y%m%dT%H%M%OS6Z", tz = "UTC")
  paste0(stamp, "-pid", Sys.getpid(), "-s", .runtime_log_state$serial)
}

append_runtime_log_line <- function(line) {
  path <- runtime_log_path()
  assert(!is.null(path), "A runtime log session is not active.")
  tryCatch(
    cat(line, "\n", file = path, append = TRUE, sep = ""),
    error = function(error) {
      stop(
        "Failed to append runtime log file ", path, ": ",
        conditionMessage(error),
        call. = FALSE
      )
    }
  )
  invisible(path)
}

runtime_log_event <- function(event, ..., console = TRUE) {
  assert(runtime_log_active(), "runtime_log_event() requires an active runtime log session.")
  event <- validate_scalar_string(event, "runtime log event")
  fields <- list(...)
  field_names <- names(fields)

  formatted_fields <- character(0L)

  if (length(fields)) {
    formatted_fields <- vapply(
      seq_along(fields),
      function(i) {
        key <- if (length(field_names)) field_names[[i]] else ""
        value <- fields[[i]]

        if (!nzchar(key)) {
          return(as.character(value))
        }

        paste0(key, "=", format_runtime_log_value(value))
      },
      character(1L)
    )
  }

  suffix <- if (length(formatted_fields)) {
    paste0(" | ", paste(formatted_fields, collapse = " "))
  } else {
    ""
  }

  line <- sprintf(
    "[%s] %s | session_id=%s%s",
    runtime_log_timestamp(),
    as.character(event),
    .runtime_log_state$session_id,
    suffix
  )

  # Persist first so a selected diagnostic can never appear only on screen.
  append_runtime_log_line(line)
  if (isTRUE(console)) {
    cat(line, "\n", sep = "")
    flush.console()
  }

  invisible(NULL)
}

# Execute one zero-argument operation inside an append-only log session. The
# directory/file is prepared only when called; nested sessions fail. The active
# state is always cleared, and completion/error events retain their old order.
with_runtime_log <- function(path, stage, operation, mode, context = list(), code) {
  assert(!runtime_log_active(), "Nested runtime log sessions are not supported.")
  path <- validate_path_param(path, "runtime log path")
  stage <- validate_scalar_string(stage, "runtime log stage")
  operation <- validate_scalar_string(operation, "runtime log operation")
  mode <- validate_scalar_string(mode, "runtime log mode")
  assert(is.list(context), "runtime log context must be a list.")
  assert(is.function(code) && length(formals(code)) == 0L,
         "runtime log code must be a zero-argument function.")
  if (length(context)) {
    assert(!is.null(names(context)) && all(nzchar(names(context))),
           "runtime log context fields must all be named.")
  }

  directory <- dirname(path)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  assert(dir.exists(directory), paste0("Could not create runtime log directory: ", directory))
  tryCatch({
    connection <- file(path, open = "a")
    close(connection)
  }, error = function(error) {
    stop(
      "Could not open runtime log file ", path, ": ",
      conditionMessage(error),
      call. = FALSE
    )
  })

  started_time <- Sys.time()
  .runtime_log_state$path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  .runtime_log_state$session_id <- runtime_log_session_id(started_time)
  .runtime_log_state$started <- unname(proc.time()[["elapsed"]])
  on.exit(clear_runtime_log(), add = TRUE)

  do.call(runtime_log_event, c(list(
    event = "runtime_session_start",
    stage = stage,
    operation = operation,
    mode = mode,
    log_file = .runtime_log_state$path
  ), context))

  result <- tryCatch(
    code(),
    error = function(error) {
      elapsed <- unname(proc.time()[["elapsed"]] - .runtime_log_state$started)
      runtime_log_event(
        "runtime_session_end",
        stage = stage,
        operation = operation,
        mode = mode,
        status = "error",
        condition_class = class(error)[[1L]],
        condition_message = conditionMessage(error),
        elapsed_seconds = sprintf("%.3f", elapsed)
      )
      stop(error)
    }
  )

  elapsed <- unname(proc.time()[["elapsed"]] - .runtime_log_state$started)
  runtime_log_event(
    "runtime_session_end",
    stage = stage,
    operation = operation,
    mode = mode,
    status = "complete",
    elapsed_seconds = sprintf("%.3f", elapsed)
  )
  result
}
