# Stage 5 chronological orchestration.
#
# The scientific per-species kernel remains in patch_processing.R. This module
# owns status inspection, checkpoint reuse, progress reporting, and final
# promotion while keeping existing final artifacts safe during a replacement.
# Active work appends setup and per-attempt species timings to the scenario-level
# runtime log; inspection creates no log session.

if (!exists("project_source", mode = "function")) {
  source(file.path("R", "project_loader.R"))
}
project_source(c(
  # Shared contracts and application adapters.
  "R/project_utils.R", "R/project_artifacts.R", "R/project_rmd.R",
  "R/project_runtime_log.R", "R/project_transactions.R",
  "R/project_paths.R", "R/project_config.R", "R/demographic_contract.R",
  "R/analysis_contract.R", "R/priority_run_config.R",
  "R/species_table_config.R", "R/patch_contract.R", "R/patch_config.R",
  "R/application_storage.R", "R/storage_artifacts.R", "R/model_lifecycle.R",
  "R/application_state.R",
  "R/application_context.R", "R/application_rmd.R",
  # Stage 5 spatial inputs, processing, and recovery.
  "R/patch_species_inputs.R", "R/patch_habitat.R",
  "R/patch_processing.R", "R/stage5_checkpoints.R"
))

stage5_elapsed_seconds <- function(started) {
  elapsed <- unname(proc.time()[["elapsed"]] - started)
  assert(length(elapsed) == 1L && is.finite(elapsed) && elapsed >= 0,
         "Stage 5 measured an invalid elapsed time.")
  elapsed
}

stage5_validate_timings <- function(values, label) {
  values <- unlist(values, use.names = TRUE)
  assert(
    length(values) > 0L && all(is.finite(values)) && all(values >= 0),
    paste0(label, " timing values must be finite and nonnegative.")
  )
  values
}

stage5_log_setup_timing <- function(values, status) {
  values <- stage5_validate_timings(values, "Stage 5 setup")
  expected <- c(
    "final_validation_seconds", "backend_setup_seconds", "manifest_seconds",
    "checkpoint_initialize_seconds", "raster_context_seconds", "total_seconds"
  )
  assert(identical(names(values), expected),
         "Stage 5 setup timing fields are incomplete or out of order.")
  status <- validate_scalar_choice(status, c("processing", "final_reused"), "setup status")
  do.call(runtime_log_event, c(
    list(event = "stage5_setup_timing", status = status),
    as.list(stats::setNames(sprintf("%.6f", values), names(values)))
  ))
  invisible(values)
}

stage5_log_species_timing <- function(species, taxon, sdm, status, values) {
  species <- validate_scalar_string(species, "species")
  taxon <- validate_scalar_string(taxon, "taxon")
  sdm <- validate_scalar_string(sdm, "sdm")
  status <- validate_scalar_string(status, "status")
  values <- stage5_validate_timings(values, "Stage 5 species")
  expected <- c(
    "scientific_seconds", "clump_seconds", "polygon_seconds",
    "distance_seconds", "other_scientific_seconds", "checkpoint_seconds",
    "cleanup_gc_seconds", "orchestration_seconds", "total_seconds"
  )
  assert(identical(names(values), expected),
         "Stage 5 species timing fields are incomplete or out of order.")
  do.call(runtime_log_event, c(
    list(
      event = "stage5_species_timing", species = species,
      taxon = taxon, sdm = sdm, status = status
    ),
    as.list(stats::setNames(sprintf("%.6f", values), names(values)))
  ))
  invisible(values)
}

stage5_status_row <- function(config, final, final_compatible,
                              checkpoint_available, checkpoint_compatible) {
  present <- stage5_final_artifacts_present(config)
  data.frame(
    taxa = config$taxa,
    final_complete = !is.null(final),
    final_compatible = final_compatible,
    final_provenance = if (is.null(final)) {
      "none"
    } else if (is.null(final$metadata)) {
      "unavailable"
    } else {
      "recorded"
    },
    patch_directory = present[["patch_dir"]],
    patch_lookup = present[["patch_lookup"]],
    connectivity = present[["connectivity"]],
    checkpoint = checkpoint_available,
    checkpoint_compatible = checkpoint_compatible,
    stringsAsFactors = FALSE
  )
}

# Inspect only inexpensive filesystem state. Exact compatibility requires
# reading large final objects and hashing every selected SDM, so it belongs to
# active resume/restart preflight rather than the read-only inspection path.
stage5_inspection_status <- function(config, selected_species) {
  present <- c(
    patch_dir = dir.exists(config$paths$patch_dir),
    patch_lookup = file.exists(config$paths$patch_lookup_rds),
    connectivity = file.exists(config$paths$connectivity_rds),
    metadata = file.exists(config$paths$metadata_rds)
  )
  required <- present[c("patch_dir", "patch_lookup", "connectivity")]
  final_state <- if (all(required)) {
    "complete"
  } else if (any(required)) {
    "partial"
  } else {
    "missing"
  }

  checkpoint_available <- dir.exists(config$paths$checkpoint_dir)
  completed_dir <- file.path(config$paths$checkpoint_dir, "completed")
  completed_species <- if (dir.exists(completed_dir)) {
    length(list.dirs(completed_dir, recursive = FALSE, full.names = TRUE))
  } else {
    0L
  }

  data.frame(
    taxa = config$taxa,
    selected_species = nrow(selected_species),
    requested_backend = config$clump_backend,
    final_state = final_state,
    final_complete = identical(final_state, "complete"),
    final_provenance = if (present[["metadata"]]) {
      "recorded"
    } else if (identical(final_state, "complete")) {
      "unavailable"
    } else {
      "none"
    },
    patch_directory = present[["patch_dir"]],
    patch_lookup = present[["patch_lookup"]],
    connectivity = present[["connectivity"]],
    checkpoint = checkpoint_available,
    checkpoint_manifest = file.exists(file.path(
      config$paths$checkpoint_dir, "manifest.rds"
    )),
    completed_checkpoint_species = completed_species,
    compatibility = "validated only by resume/restart",
    stringsAsFactors = FALSE
  )
}

stage5_final_compatible <- function(final, manifest) {
  if (is.null(final) || is.null(final$metadata)) return(FALSE)
  stage5_manifest_compatible(final$metadata$manifest, manifest)
}

stage5_attempt_species <- function(row, index, total, config, raster_context, backend,
                                   work_dir) {
  assert(dir.exists(work_dir), "Stage 5 species work directory does not exist.")
  processing_paths <- config$paths
  processing_paths$patch_dir <- work_dir
  started <- Sys.time()
  tryCatch(
    process_patch_species(
      row,
      processing_paths,
      run_options = list(overwrite_patch_rasters = TRUE),
      habitat_masks = raster_context$habitat_masks,
      template = raster_context$template,
      cell_area_km2 = raster_context$cell_area_km2,
      clump_backend = backend,
      verbose = config$verbose,
      label_map = habitat_label_map(),
      species_index = index,
      species_total = total
    ),
    error = function(e) {
      elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
      log_patch_step(
        row$scientificName,
        "ERROR",
        stop_stage = "unhandled_error",
        reason = conditionMessage(e),
        elapsed_seconds = elapsed,
        species_index = index,
        species_total = total
      )
      list(
        status = "error",
        species = row$scientificName,
        class_lc = row$class_lc,
        sdm_method = row$sdm_method,
        stop_stage = "unhandled_error",
        stop_reason = conditionMessage(e),
        elapsed_seconds = elapsed,
        error = conditionMessage(e)
      )
    }
  )
}

run_stage5_workflow_impl <- function(config) {
  started <- Sys.time()
  species_table <- read_patch_species_base(
    config$paths$species_csv, contract = config$contract
  )
  selected_species <- filter_patch_species_taxa(
    species_table,
    do_mammals = config$selected_mammals,
    do_birds = config$selected_birds
  )
  selected_counts <- selected_species |>
    dplyr::count(class_lc, sdm_method, name = "n") |>
    dplyr::arrange(class_lc, sdm_method)
  log_msg("Stage 5 | selected species by taxon and SDM")
  print(selected_counts)

  # Inspection is intentionally metadata-only. Exact scientific compatibility
  # requires hashing every SDM and deserializing the large final graph object;
  # active modes perform those checks immediately before reuse or replacement.
  if (identical(config$mode, "inspect")) {
    return(list(
      action = "inspect",
      selected_species_count = nrow(selected_species),
      status = stage5_inspection_status(config, selected_species),
      status_summary = data.frame(),
      outputs_committed = FALSE,
      commit_blocked = FALSE,
      elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))
    ))
  }

  setup_started <- unname(proc.time()[["elapsed"]])
  phase_started <- unname(proc.time()[["elapsed"]])
  final <- validate_stage5_final_outputs(config)
  final_validation_seconds <- stage5_elapsed_seconds(phase_started)

  phase_started <- unname(proc.time()[["elapsed"]])
  backend <- resolve_patch_clump_backend(config$clump_backend, config$grass_dir)
  backend_setup_seconds <- stage5_elapsed_seconds(phase_started)

  phase_started <- unname(proc.time()[["elapsed"]])
  manifest <- build_stage5_manifest(config, selected_species, backend)
  manifest_path <- stage5_manifest_path(config)
  checkpoint_available <- dir.exists(config$paths$checkpoint_dir)
  checkpoint_compatible <- FALSE
  if (checkpoint_available && file.exists(manifest_path)) {
    checkpoint_compatible <- stage5_manifest_compatible(readRDS(manifest_path), manifest)
  }
  final_compatible <- stage5_final_compatible(final, manifest)
  status <- stage5_status_row(
    config,
    final,
    final_compatible,
    checkpoint_available,
    checkpoint_compatible
  )
  manifest_seconds <- stage5_elapsed_seconds(phase_started)

  if (
    identical(config$mode, "resume") &&
      !checkpoint_available &&
      final_compatible
  ) {
    stage5_log_setup_timing(list(
      final_validation_seconds = final_validation_seconds,
      backend_setup_seconds = backend_setup_seconds,
      manifest_seconds = manifest_seconds,
      checkpoint_initialize_seconds = 0,
      raster_context_seconds = 0,
      total_seconds = stage5_elapsed_seconds(setup_started)
    ), "final_reused")
    return(list(
      action = "final_reused",
      selected_species_count = nrow(selected_species),
      status = status,
      status_summary = patch_status_summary(list()),
      patch_rows = nrow(final$patch_lookup),
      population_units = length(final$connectivity),
      outputs_committed = FALSE,
      commit_blocked = FALSE,
      backend = backend,
      elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))
    ))
  }

  phase_started <- unname(proc.time()[["elapsed"]])
  initialize_stage5_checkpoint(
    config,
    manifest,
    restart = identical(config$mode, "restart")
  )
  checkpoint_initialize_seconds <- stage5_elapsed_seconds(phase_started)

  phase_started <- unname(proc.time()[["elapsed"]])
  raster_context <- load_patch_raster_context(
    config$paths$landcover_tif, config$study_area
  )
  raster_context_seconds <- stage5_elapsed_seconds(phase_started)
  log_msg(
    "Stage 5 | raster context ready | habitat_layers=", terra::nlyr(raster_context$habitat_masks),
    "| cells=", terra::ncell(raster_context$template),
    "| elapsed=", round(raster_context_seconds, 2), "s"
  )
  stage5_log_setup_timing(list(
    final_validation_seconds = final_validation_seconds,
    backend_setup_seconds = backend_setup_seconds,
    manifest_seconds = manifest_seconds,
    checkpoint_initialize_seconds = checkpoint_initialize_seconds,
    raster_context_seconds = raster_context_seconds,
    total_seconds = stage5_elapsed_seconds(setup_started)
  ), "processing")

  results <- vector("list", nrow(selected_species))
  completed_before <- 0L
  durations <- numeric()
  for (i in seq_len(nrow(selected_species))) {
    species <- selected_species$scientificName[[i]]
    cached <- read_stage5_species_checkpoint(
      config,
      i,
      species,
      template = raster_context$template
    )
    if (!is.null(cached)) {
      results[[i]] <- cached
      completed_before <- completed_before + 1L
      log_patch_step(
        species,
        "CHECKPOINT_REUSED",
        species_index = i,
        species_total = nrow(selected_species)
      )
      next
    }

    species_started <- unname(proc.time()[["elapsed"]])
    scientific_seconds <- 0
    checkpoint_seconds <- 0
    work_dir <- tempfile(
      paste0("stage5_species_", sprintf("%04d", i), "_"),
      tmpdir = dirname(config$paths$checkpoint_dir)
    )
    assert(dir.create(work_dir), "Could not create the Stage 5 species work directory.")
    result <- tryCatch({
      scientific_started <- unname(proc.time()[["elapsed"]])
      candidate <- stage5_attempt_species(
        selected_species[i, , drop = FALSE],
        i,
        nrow(selected_species),
        config,
        raster_context,
        backend,
        work_dir
      )
      scientific_seconds <- stage5_elapsed_seconds(scientific_started)
      if (!identical(candidate$status, "error")) {
        checkpoint_started <- unname(proc.time()[["elapsed"]])
        candidate <- write_stage5_species_checkpoint(
          candidate,
          config,
          i,
          template = raster_context$template
        )
        checkpoint_seconds <- stage5_elapsed_seconds(checkpoint_started)
      }
      candidate
    }, finally = {
      unlink(work_dir, recursive = TRUE, force = TRUE)
    })
    results[[i]] <- result
    cleanup_started <- unname(proc.time()[["elapsed"]])
    gc(FALSE)
    cleanup_gc_seconds <- stage5_elapsed_seconds(cleanup_started)
    phase_seconds <- function(field) {
      value <- suppressWarnings(as.numeric(result[[field]] %||% 0))
      if (length(value) != 1L || !is.finite(value) || value < 0) 0 else value
    }
    clump_seconds <- phase_seconds("clump_elapsed_seconds")
    polygon_seconds <- phase_seconds("polygon_elapsed_seconds")
    distance_seconds <- phase_seconds("distance_elapsed_seconds")
    classified_scientific <- sum(c(
      clump_seconds, polygon_seconds, distance_seconds
    ))
    scientific_seconds <- max(scientific_seconds, classified_scientific)
    other_scientific_seconds <- max(
      0, scientific_seconds - classified_scientific
    )
    total_seconds <- stage5_elapsed_seconds(species_started)
    accounted <- scientific_seconds + checkpoint_seconds + cleanup_gc_seconds
    # A timer-resolution edge can make separately sampled phases exceed the
    # enclosing measurement by a few microseconds. Preserve nonnegative,
    # internally accountable diagnostics without altering any workflow step.
    total_seconds <- max(total_seconds, accounted)
    orchestration_seconds <- max(0, total_seconds - accounted)
    stage5_log_species_timing(
      species = species,
      taxon = as.character(selected_species$class_lc[[i]]),
      sdm = as.character(selected_species$sdm_method[[i]]),
      status = as.character(result$status),
      values = list(
      scientific_seconds = scientific_seconds,
      clump_seconds = clump_seconds,
      polygon_seconds = polygon_seconds,
      distance_seconds = distance_seconds,
      other_scientific_seconds = other_scientific_seconds,
      checkpoint_seconds = checkpoint_seconds,
      cleanup_gc_seconds = cleanup_gc_seconds,
      orchestration_seconds = orchestration_seconds,
      total_seconds = total_seconds
      )
    )
    durations <- c(durations, total_seconds)
    if (identical(result$status, "error")) next
    remaining <- nrow(selected_species) - i
    eta <- if (length(durations)) remaining * mean(utils::tail(durations, 20L)) else NA_real_
    log_msg(
      "Stage 5 | checkpoint committed | completed=", i, "/", nrow(selected_species),
      "| species=", species,
      "| ETA=", if (is.finite(eta)) paste0(round(eta / 60, 1), " min") else "NA"
    )
  }

  status_summary <- patch_status_summary(results)
  errors <- which(status_summary$status == "error")
  if (length(errors)) {
    return(list(
      action = "checkpointed_with_errors",
      selected_species_count = nrow(selected_species),
      completed_before = completed_before,
      status_summary = status_summary,
      patch_rows = 0L,
      population_units = 0L,
      outputs_committed = FALSE,
      commit_blocked = TRUE,
      backend = backend,
      elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))
    ))
  }

  promoted <- promote_stage5_outputs(config, manifest, results)
  list(
    action = "committed",
    selected_species_count = nrow(selected_species),
    completed_before = completed_before,
    status_summary = status_summary,
    patch_rows = nrow(promoted$patch_lookup),
    population_units = length(promoted$connectivity),
    outputs_committed = TRUE,
    commit_blocked = FALSE,
    backend = backend,
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))
  )
}

run_stage5 <- function(config) {
  assert(is.list(config) && !is.null(config$paths$runtime_log),
         "run_stage5() requires a stage5_config() result.")
  check_stage5_preflight(config)
  if (identical(config$mode, "inspect")) {
    return(run_stage5_workflow_impl(config))
  }
  # Runtime discovery remains absent from workflow sourcing and inspection.
  project_source("R/patch_runtime.R")
  load_patch_packages(include_figures = FALSE)
  result <- with_runtime_log(
    path = config$paths$runtime_log,
    stage = "5",
    operation = "patch_and_connectivity_construction",
    mode = config$mode,
    context = list(
      taxa = config$taxa,
      clump_backend = config$clump_backend,
      threshold = as.integer(config$contract$quasi_extinction_abundance)
    ),
    code = function() run_stage5_workflow_impl(config)
  )
  result$runtime_log_path <- config$paths$runtime_log
  result
}

report_stage5 <- function(result) {
  if (identical(result$action, "inspect")) {
    print(result$status)
    log_msg(
      "Stage 5 | inspection complete | selected=", result$selected_species_count,
      "| files_written=none",
      "| elapsed=", round(result$elapsed_seconds, 2), "s"
    )
    return(invisible(result))
  }

  if (identical(result$action, "final_reused")) {
    log_msg(
      "Stage 5 | existing compatible final outputs retained",
      "| selected=", result$selected_species_count,
      "| patch_rows=", result$patch_rows,
      "| population_units=", result$population_units,
      "| files_written=none",
      "| elapsed=", round(result$elapsed_seconds, 2), "s"
    )
    return(invisible(result))
  }

  status <- result$status_summary
  log_msg(
    "Stage 5 | summary | action=", result$action,
    "| selected=", result$selected_species_count,
    "| retained=", sum(status$status == "retained", na.rm = TRUE),
    "| filtered=", sum(grepl("^skipped_", status$status), na.rm = TRUE),
    "| errors=", sum(status$status == "error", na.rm = TRUE)
  )
  if (nrow(status)) {
    print(status |>
      dplyr::count(class_lc, sdm_method, status, stop_stage, name = "n") |>
      dplyr::arrange(class_lc, sdm_method, status, stop_stage))
  }
  log_msg(
    "Stage 5 | outputs_committed=", result$outputs_committed,
    "| patch_rows=", result$patch_rows,
    "| population_units=", result$population_units,
    "| elapsed=", round(result$elapsed_seconds, 2), "s"
  )
  if (isTRUE(result$commit_blocked)) {
    errors <- status[status$status == "error", c("species", "error"), drop = FALSE]
    patch_abort(
      "Stage 5 encountered technical species errors. Completed checkpoints were retained and final outputs were not replaced.\n",
      paste(apply(errors, 1, paste, collapse = " | "), collapse = "\n")
    )
  }
  invisible(result)
}
