# Recovery state and manifest transactions for named spatial applications.
#
# Application workflows load this definition-only module after storage and
# analysis contracts. It preserves the existing schema and artifact records;
# writes occur only when an explicit state/manifest function is called.

application_state_schema <- function() 3L

# Describe one durable file or directory relative to `root`. The returned
# record includes sizes and checksums; missing paths fail before any state write.
application_file_record <- function(path, root = project_paths()$root) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  relative <- storage_relative_path(path, root)
  relative <- if (identical(relative, path)) NULL else relative
  if (dir.exists(path)) {
    files <- list.files(
      path, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE
    )
    files <- files[file.exists(files) & !dir.exists(files)]
    entries <- data.frame(
      path = substring(
        normalizePath(files, winslash = "/", mustWork = TRUE),
        nchar(path) + 2L
      ),
      size = as.numeric(file.info(files)$size),
      md5 = unname(tools::md5sum(files)),
      stringsAsFactors = FALSE
    )
    record <- list(
      path = path, type = "directory", files = length(files),
      size = sum(file.info(files)$size, na.rm = TRUE), md5 = NA_character_,
      entries = entries
    )
    if (!is.null(relative)) record$relative_path <- relative
    return(record)
  }
  info <- file.info(path)
  record <- list(
    path = path, type = "file", files = 1L,
    size = unname(as.numeric(info$size)), md5 = unname(tools::md5sum(path))
  )
  if (!is.null(relative)) record$relative_path <- relative
  record
}

# Convert an ordered path set to validated artifact records. This is read-only
# and fails atomically if any requested output is absent.
application_artifact_records <- function(paths, root = project_paths()$root) {
  paths <- unique(as.character(paths))
  paths <- paths[!is.na(paths) & nzchar(paths)]
  missing <- paths[!file.exists(paths) & !dir.exists(paths)]
  assert(!length(missing), paste0(
    "Completed Stage 8 artifact(s) are missing:\n", paste(missing, collapse = "\n")
  ))
  lapply(paths, application_file_record, root = root)
}

# Recheck a saved record against current files without modifying either one.
# Directory validation checks every recorded member checksum.
application_record_valid <- function(record, root = project_paths()$root) {
  if (!is.list(record)) return(FALSE)
  outputs <- record$outputs %||% list()
  upstream <- record$upstream %||% list()
  if (!length(outputs) && !identical(record$status, "not_applicable")) {
    return(FALSE)
  }
  records <- c(outputs, upstream)
  all(vapply(records, function(x) {
    path <- storage_record_path(x, root)
    if (!file.exists(path) && !dir.exists(path)) return(FALSE)
    if (identical(x$type, "directory")) {
      files <- file.path(path, x$entries$path)
      if (!all(file.exists(files) & !dir.exists(files))) return(FALSE)
      return(identical(unname(tools::md5sum(files)), x$entries$md5))
    }
    identical(unname(tools::md5sum(path)), x$md5)
  }, logical(1L)))
}

# Create the in-memory recovery envelope for one application/scenario/run.
new_application_state <- function(config) {
  list(
    schema = application_state_schema(), application = config$application,
    threshold_tag = config$scenario$threshold_tag,
    run_tag = config$scenario$run_tag,
    stages = list()
  )
}

# Read the run recovery state and import scenario-owned records when safe.
# New schedules inherit all published scenario stages. An existing schedule
# imports a subsequently published shared Stage 6 initialization only when its
# own record is missing; an existing record may be newer if publication was
# interrupted between the run-state and scenario-manifest transactions.
# Unsupported schemas and identities fail without changing disk state.
read_application_state <- function(config) {
  path <- config$scenario$run_state
  if (!file.exists(path)) {
    state <- new_application_state(config)
    # Stage 4--5 products and the shared Stage 6 initialization belong to the
    # threshold scenario and can seed a new removal schedule. Curve execution
    # and Stage 7.1 records remain isolated to their run.
    if (file.exists(config$scenario$manifest)) {
      scenario <- readRDS(config$scenario$manifest)
      keys <- intersect(
        c(application_scenario_core_stages(),
          application_scenario_optional_stages(),
          application_scenario_shared_stages()),
        names(scenario$stages)
      )
      state$stages[keys] <- scenario$stages[keys]
    }
    return(state)
  }
  state <- readRDS(path)
  assert(
    identical(state$schema, application_state_schema()),
    paste0(
      "Unsupported Stage 8 recovery-state schema. Inspect the folder and ",
      "restart the selected boundary."
    )
  )
  assert(
    identical(state$application, config$application) &&
      identical(state$threshold_tag, config$scenario$threshold_tag) &&
      identical(state$run_tag, config$scenario$run_tag),
    paste0(
      "Stage 8 recovery state belongs to a different application, threshold, ",
      "or removal schedule."
    )
  )

  shared_key <- application_scenario_shared_stages()
  if (is.null(state$stages[[shared_key]]) &&
      file.exists(config$scenario$manifest)) {
    scenario <- readRDS(config$scenario$manifest)
    if (!is.null(scenario$stages[[shared_key]])) {
      state$stages[[shared_key]] <- scenario$stages[[shared_key]]
    }
  }
  state
}

# Atomically replace the run recovery file and invisibly return `state`.
write_application_state <- function(state, config) {
  storage_atomic_save_rds(state, config$scenario$run_state, "Stage 8 recovery state")
  invisible(state)
}

# Construct the stable recovery key for a stage, optionally scoped to a curve.
application_stage_key <- function(stage, curve = NULL) {
  if (is.null(curve)) stage else paste0(stage, "_", curve)
}

# Mark a valid completed boundary reused in memory; return NULL when its record
# or artifacts cannot be reused. This function does not persist the result.
application_reuse_boundary <- function(state, key, root = project_paths()$root,
                                       expected_summary = list()) {
  assert(is.list(expected_summary), "expected_summary must be a list.")
  if (length(expected_summary)) {
    assert(
      !is.null(names(expected_summary)) &&
        all(nzchar(names(expected_summary))) && !anyDuplicated(names(expected_summary)),
      "expected_summary must have unique non-empty names."
    )
  }
  record <- state$stages[[key]]
  if (is.null(record) || !application_record_valid(record, root)) return(NULL)
  matches <- vapply(names(expected_summary), function(field) {
    field %in% names(record$summary) &&
      identical(record$summary[[field]], expected_summary[[field]])
  }, logical(1L))
  if (length(matches) && !all(matches)) return(NULL)
  state$stages[[key]]$status <- "reused"
  state
}

# Validate output/upstream artifacts, commit one completion record, and then
# atomically persist the full state. No record is written for missing outputs.
application_finish_boundary <- function(state, key, status, outputs, summary,
                                        upstream = character(), config) {
  assert(
    status %in% c("complete", "reused", "not_applicable"),
    "Invalid Stage 8 completion status."
  )
  state$stages[[key]] <- list(
    status = status,
    outputs = if (identical(status, "not_applicable") && !length(outputs)) {
      list()
    } else {
      application_artifact_records(outputs, config$base_paths$root)
    },
    upstream = if (length(upstream)) {
      application_artifact_records(upstream, config$base_paths$root)
    } else list(),
    summary = summary
  )
  write_application_state(state, config)
  state
}

# Return compact stage/status rows without reading or hashing artifacts.
application_status_table <- function(state) {
  if (!length(state$stages)) return(data.frame(
    stage = character(), status = character(), stringsAsFactors = FALSE
  ))
  data.frame(
    stage = names(state$stages),
    status = vapply(state$stages, `[[`, character(1L), "status"),
    stringsAsFactors = FALSE, row.names = NULL
  )
}

# Flatten recorded outputs for display, resolving relocatable paths under root.
application_artifact_table <- function(state, root = project_paths()$root) {
  rows <- lapply(names(state$stages), function(stage) {
    outputs <- state$stages[[stage]]$outputs
    if (!length(outputs)) return(NULL)
    data.frame(
      stage = stage,
      path = vapply(outputs, storage_record_path, character(1L), root = root),
      type = vapply(outputs, `[[`, character(1L), "type"),
      size = vapply(outputs, `[[`, numeric(1L), "size"),
      stringsAsFactors = FALSE, row.names = NULL
    )
  })
  do.call(rbind, Filter(Negate(is.null), rows)) %||% data.frame(
    stage = character(), path = character(), type = character(), size = numeric()
  )
}

# Aggregate recorded byte counts by stage for presentation only.
summarize_application_artifacts <- function(x) {
  if (!nrow(x)) return(x)
  aggregate(size ~ stage, x, function(z) sum(z, na.rm = TRUE)) |>
    stats::setNames(c("stage", "total_bytes"))
}

# Recursively report contract-field differences in stable field order.
application_manifest_differences <- function(found, expected, prefix = "") {
  fields <- union(names(found), names(expected))
  unlist(lapply(fields, function(field) {
    label <- if (nzchar(prefix)) paste0(prefix, ".", field) else field
    if (!field %in% names(found)) {
      return(paste0(label, ": missing from existing manifest"))
    }
    if (!field %in% names(expected)) {
      return(paste0(label, ": not expected by current code"))
    }
    a <- found[[field]]
    b <- expected[[field]]
    if (is.list(a) && is.list(b) && !is.data.frame(a) && !is.data.frame(b)) {
      return(application_manifest_differences(a, b, label))
    }
    if (!identical(a, b)) {
      paste0(label, ": existing and requested values differ")
    } else character()
  }), use.names = FALSE)
}

# Create a manifest atomically or return the compatible existing manifest.
# Contract differences fail with `change_hint`; existing files are never patched.
application_write_or_validate_manifest <- function(path, manifest, label,
                                                    change_hint) {
  if (!file.exists(path)) {
    storage_atomic_save_rds(manifest, path, label)
    return(manifest)
  }
  found <- readRDS(path)
  differences <- application_manifest_differences(found$contract, manifest$contract)
  assert(!length(differences), paste0(
    label, " already exists for different inputs:\n- ",
    paste(differences, collapse = "\n- "), "\n", change_hint
  ))
  found
}

# Scenario inputs are immutable, but stage records may change when a completed
# presentation step is regenerated. Validate the contract before refreshing.
application_write_scenario_manifest <- function(path, manifest) {
  if (file.exists(path)) {
    found <- readRDS(path)
    differences <- application_manifest_differences(found$contract, manifest$contract)
    assert(!length(differences), paste0(
      "threshold-scenario manifest already exists for different inputs:\n- ",
      paste(differences, collapse = "\n- "),
      "\nRestart from stage_4 to rebuild this threshold scenario."
    ))
  }
  storage_atomic_save_rds(manifest, path, "threshold-scenario manifest")
  manifest
}
