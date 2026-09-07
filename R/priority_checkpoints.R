# Stage-boundary checkpoints for the spatial-prioritization pipeline.
#
# Schema 4 stores completed scientific state in canonical species, patch, PU,
# graph, event, and cell order. Runtime-only compact and sparse indexes are
# reconstructed after resume. Checkpoints are written uncompressed and
# atomically; retention is applied only after a new checkpoint validates.


# ---- Schema and deterministic checkpoint discovery --------------------------

# Filenames encode the completed stage; discovery sorts by that numeric stage
# so "latest" has one unambiguous meaning on every filesystem.
priority_checkpoint_schema_version <- function() {
  4L
}

priority_checkpoint_dir <- function(output_dir) {
  file.path(output_dir, "checkpoints")
}

priority_checkpoint_metadata_fields <- function() {
  c(
    "checkpoint_schema_version",
    "initialization_schema_version",
    "initialization_created_at",
    "curve_label",
    "do_mammals",
    "do_birds",
    "do_ppm",
    "do_rangebag",
    "retained_species",
    "n_cells",
    "cells_to_remove_per_iteration",
    "pruning_iterations_per_stage",
    "analysis_contract"
  )
}

priority_checkpoint_state_fields <- function() {
  c(
    "stage_index",
    "completed_stages",
    "frontier_exhausted",
    "patch_table",
    "pu_graphs_by_key",
    "alive_species_count_by_cell",
    "patch_id_by_species_env",
    "removal_order_by_cell",
    "removal_event_rows",
    "removal_step_counter"
  )
}

priority_checkpoint_required_fields <- function() {
  c("metadata", priority_checkpoint_state_fields())
}

validate_checkpoint_nonnegative_integer <- function(x, label) {
  assert(
    is.numeric(x) &&
      length(x) == 1L &&
      is.finite(x) &&
      x >= 0 &&
      x == floor(x),
    paste0(label, " must be a single non-negative integer.")
  )
  as.integer(x)
}

validate_checkpoint_retained_species <- function(x, label = "retained_species") {
  assert(
    is.character(x),
    paste0(label, " must be a character vector of species names.")
  )
  x <- trimws(unname(x))
  assert(
    length(x) > 0L &&
      all(!is.na(x)) &&
      all(nzchar(trimws(x))),
    paste0(label, " must contain at least one non-blank species name.")
  )
  x
}

priority_checkpoint_filename <- function(stage_index) {
  stage_index <- validate_priority_count(stage_index, "stage_index")
  sprintf("priority_checkpoint_stage_%04d.rds", stage_index)
}

priority_checkpoint_path <- function(checkpoint_dir, stage_index) {
  file.path(checkpoint_dir, priority_checkpoint_filename(stage_index))
}

priority_checkpoint_stage_from_path <- function(path) {
  file_name <- basename(path)
  match <- regexec("^priority_checkpoint_stage_([0-9]+)\\.rds$", file_name)
  parsed <- regmatches(file_name, match)[[1L]]

  if (length(parsed) != 2L) {
    return(NA_integer_)
  }

  as.integer(parsed[[2L]])
}

list_priority_checkpoint_files <- function(checkpoint_dir) {
  if (!dir.exists(checkpoint_dir)) {
    return(data.frame(
      path = character(),
      stage = integer(),
      stringsAsFactors = FALSE
    ))
  }

  files <- list.files(
    checkpoint_dir,
    pattern = "^priority_checkpoint_stage_[0-9]+\\.rds$",
    full.names = TRUE
  )

  if (!length(files)) {
    return(data.frame(
      path = character(),
      stage = integer(),
      stringsAsFactors = FALSE
    ))
  }

  stages <- vapply(files, priority_checkpoint_stage_from_path, integer(1L))
  keep <- !is.na(stages)

  files <- files[keep]
  stages <- stages[keep]
  ordering <- order(stages, files)

  data.frame(
    path = files[ordering],
    stage = stages[ordering],
    stringsAsFactors = FALSE
  )
}

resolve_priority_checkpoint_path <- function(resume_checkpoint_path, checkpoint_dir) {
  if (is.null(resume_checkpoint_path)) {
    return(NULL)
  }

  resume_checkpoint_path <- validate_scalar_string(
    resume_checkpoint_path,
    "resume_checkpoint_path"
  )

  if (identical(resume_checkpoint_path, "latest")) {
    checkpoints <- list_priority_checkpoint_files(checkpoint_dir)

    if (!nrow(checkpoints)) {
      stop("No priority checkpoints found in: ", checkpoint_dir)
    }

    return(checkpoints$path[[nrow(checkpoints)]])
  }

  resume_checkpoint_path <- path.expand(resume_checkpoint_path)

  if (!file.exists(resume_checkpoint_path)) {
    stop("Requested priority checkpoint does not exist: ", resume_checkpoint_path)
  }

  resume_checkpoint_path
}

stage_patch_lookup_stage_from_path <- function(path) {
  file_name <- basename(path)
  match <- regexec("^stage_patch_lookup_stage_([0-9]+)\\.csv$", file_name)
  parsed <- regmatches(file_name, match)[[1L]]

  if (length(parsed) != 2L) {
    return(NA_integer_)
  }

  as.integer(parsed[[2L]])
}

list_stage_patch_lookup_files <- function(patch_lookup_output_dir) {
  if (!dir.exists(patch_lookup_output_dir)) {
    return(data.frame(
      path = character(),
      stage = integer(),
      stringsAsFactors = FALSE
    ))
  }

  files <- list.files(
    patch_lookup_output_dir,
    pattern = "^stage_patch_lookup_stage_[0-9]+\\.csv$",
    full.names = TRUE
  )

  if (!length(files)) {
    return(data.frame(
      path = character(),
      stage = integer(),
      stringsAsFactors = FALSE
    ))
  }

  stages <- vapply(files, stage_patch_lookup_stage_from_path, integer(1L))
  keep <- !is.na(stages)

  files <- files[keep]
  stages <- stages[keep]
  ordering <- order(stages, files)

  data.frame(
    path = files[ordering],
    stage = stages[ordering],
    stringsAsFactors = FALSE
  )
}

# ---- Completed-stage output boundary ----------------------------------------

# A resumable checkpoint requires every lookup through its stage and forbids
# later lookups that could otherwise be mistaken for products of the resumed run.
validate_no_stage_outputs_after_checkpoint <- function(
  patch_lookup_output_dir,
  checkpoint_stage
) {
  checkpoint_stage <- validate_checkpoint_nonnegative_integer(
    checkpoint_stage,
    "checkpoint_stage"
  )
  stage_files <- list_stage_patch_lookup_files(patch_lookup_output_dir)

  if (checkpoint_stage > 0L) {
    expected_stages <- seq_len(checkpoint_stage)
    missing_stages <- setdiff(expected_stages, stage_files$stage)

    if (length(missing_stages)) {
      stop(
        "Stage patch lookup files are missing for checkpointed stage(s): ",
        paste(sprintf("%04d", missing_stages), collapse = ", "),
        ". Resume from this checkpoint requires the completed-stage lookup ",
        "tables in: ",
        patch_lookup_output_dir
      )
    }
  }

  stale <- stage_files[stage_files$stage > checkpoint_stage, , drop = FALSE]
  if (nrow(stale)) {
    stop(
      "Stage patch lookup files exist after the requested checkpoint stage. ",
      "Remove or archive these stale outputs before resuming from stage ",
      checkpoint_stage,
      ": ",
      paste(basename(stale$path), collapse = ", ")
    )
  }

  invisible(TRUE)
}

# ---- Scientific identity and compatibility metadata -------------------------

# These fields identify the exact initialized analysis. Resume validation is an
# equality check; metadata is never repaired or inferred from partial state.
priority_checkpoint_metadata <- function(
  initialization_schema_version,
  initialization_created_at,
  curve_label,
  do_mammals,
  do_birds,
  do_ppm,
  do_rangebag,
  retained_species,
  n_cells,
  cells_to_remove_per_iteration,
  pruning_iterations_per_stage,
  analysis_contract = canonical_analysis_contract()
) {
  initialization_schema_version <- validate_priority_count(
    initialization_schema_version,
    "initialization_schema_version"
  )
  initialization_created_at <- validate_scalar_string(
    initialization_created_at,
    "initialization_created_at"
  )
  curve_label <- validate_persistence_curve(curve_label, "curve_label")
  do_mammals <- validate_scalar_logical(do_mammals, "do_mammals")
  do_birds <- validate_scalar_logical(do_birds, "do_birds")
  do_ppm <- validate_scalar_logical(do_ppm, "do_ppm")
  do_rangebag <- validate_scalar_logical(do_rangebag, "do_rangebag")
  retained_species <- validate_checkpoint_retained_species(retained_species)
  n_cells <- validate_priority_count(n_cells, "n_cells")
  cells_to_remove_per_iteration <- validate_priority_count(
    cells_to_remove_per_iteration,
    "cells_to_remove_per_iteration"
  )
  pruning_iterations_per_stage <- validate_priority_count(
    pruning_iterations_per_stage,
    "pruning_iterations_per_stage"
  )
  analysis_contract <- validate_analysis_contract(
    analysis_contract, "checkpoint analysis_contract"
  )

  list(
    checkpoint_schema_version = priority_checkpoint_schema_version(),
    initialization_schema_version = initialization_schema_version,
    initialization_created_at = initialization_created_at,
    curve_label = unname(as.character(curve_label)),
    do_mammals = do_mammals,
    do_birds = do_birds,
    do_ppm = do_ppm,
    do_rangebag = do_rangebag,
    retained_species = retained_species,
    n_cells = n_cells,
    cells_to_remove_per_iteration = cells_to_remove_per_iteration,
    pruning_iterations_per_stage = pruning_iterations_per_stage,
    analysis_contract = unclass(analysis_contract)
  )
}

validate_priority_checkpoint_metadata <- function(checkpoint, expected_metadata) {
  assert(is.list(checkpoint), "Priority checkpoint must be a list.")
  assert(
    is.list(expected_metadata),
    "Expected priority checkpoint metadata must be a list."
  )

  if (is.null(checkpoint$metadata) || !is.list(checkpoint$metadata)) {
    stop("Priority checkpoint does not contain a metadata list.")
  }

  actual <- checkpoint$metadata

  required_fields <- priority_checkpoint_metadata_fields()

  missing_fields <- setdiff(required_fields, names(actual))
  if (length(missing_fields)) {
    stop(
      "Priority checkpoint metadata is missing field(s): ",
      paste(missing_fields, collapse = ", ")
      )
  }

  missing_expected_fields <- setdiff(required_fields, names(expected_metadata))
  if (length(missing_expected_fields)) {
    stop(
      "Expected priority checkpoint metadata is missing field(s): ",
      paste(missing_expected_fields, collapse = ", ")
    )
  }

  actual_version <- tryCatch(
    validate_priority_count(
      actual$checkpoint_schema_version,
      "checkpoint_schema_version"
    ),
    error = function(e) NA_integer_
  )
  expected_version <- priority_checkpoint_schema_version()
  if (!identical(actual_version, expected_version)) {
    actual_version_label <- if (length(actual_version) == 1L && !is.na(actual_version)) {
      actual_version
    } else {
      "NA"
    }

    stop(
      "Unsupported priority checkpoint schema. Expected version ",
      expected_version,
      "; found ",
      actual_version_label,
      ". Restart Stage 6 from the validated shared initialization."
    )
  }

  compare_scalar <- function(field, actual_value, expected_value) {
    if (!identical(actual_value, expected_value)) {
      stop("Priority checkpoint metadata mismatch for ", field, ".")
    }
  }

  compare_scalar(
    "initialization_schema_version",
    validate_priority_count(
      actual$initialization_schema_version,
      "checkpoint initialization_schema_version"
    ),
    validate_priority_count(
      expected_metadata$initialization_schema_version,
      "expected initialization_schema_version"
    )
  )
  compare_scalar(
    "initialization_created_at",
    validate_scalar_string(
      actual$initialization_created_at,
      "checkpoint initialization_created_at"
    ),
    validate_scalar_string(
      expected_metadata$initialization_created_at,
      "expected initialization_created_at"
    )
  )
  compare_scalar(
    "curve_label",
    validate_persistence_curve(
      actual$curve_label,
      "checkpoint metadata$curve_label"
    ),
    validate_persistence_curve(
      expected_metadata$curve_label,
      "expected metadata$curve_label"
    )
  )
  compare_scalar(
    "do_mammals",
    validate_scalar_logical(actual$do_mammals, "checkpoint do_mammals"),
    validate_scalar_logical(expected_metadata$do_mammals, "expected do_mammals")
  )
  compare_scalar(
    "do_birds",
    validate_scalar_logical(actual$do_birds, "checkpoint do_birds"),
    validate_scalar_logical(expected_metadata$do_birds, "expected do_birds")
  )
  compare_scalar(
    "do_ppm",
    validate_scalar_logical(actual$do_ppm, "checkpoint do_ppm"),
    validate_scalar_logical(expected_metadata$do_ppm, "expected do_ppm")
  )
  compare_scalar(
    "do_rangebag",
    validate_scalar_logical(actual$do_rangebag, "checkpoint do_rangebag"),
    validate_scalar_logical(expected_metadata$do_rangebag, "expected do_rangebag")
  )
  compare_scalar(
    "retained_species",
    validate_checkpoint_retained_species(actual$retained_species, "checkpoint retained_species"),
    validate_checkpoint_retained_species(expected_metadata$retained_species, "expected retained_species")
  )
  compare_scalar(
    "n_cells",
    validate_priority_count(actual$n_cells, "checkpoint n_cells"),
    validate_priority_count(expected_metadata$n_cells, "expected n_cells")
  )
  compare_scalar(
    "analysis_contract",
    validate_analysis_contract(actual$analysis_contract, "checkpoint analysis_contract"),
    validate_analysis_contract(expected_metadata$analysis_contract, "expected analysis_contract")
  )
  compare_scalar(
    "cells_to_remove_per_iteration",
    validate_priority_count(
      actual$cells_to_remove_per_iteration,
      "checkpoint cells_to_remove_per_iteration"
    ),
    validate_priority_count(
      expected_metadata$cells_to_remove_per_iteration,
      "expected cells_to_remove_per_iteration"
    )
  )
  compare_scalar(
    "pruning_iterations_per_stage",
    validate_priority_count(
      actual$pruning_iterations_per_stage,
      "checkpoint pruning_iterations_per_stage"
    ),
    validate_priority_count(
      expected_metadata$pruning_iterations_per_stage,
      "expected pruning_iterations_per_stage"
    )
  )

  invisible(TRUE)
}

# ---- Durable scientific-state validation ------------------------------------

# Runtime-only indexes are deliberately absent here. The stored objects must be
# complete enough to reconstruct them without changing scientific state.
validate_priority_checkpoint_state <- function(checkpoint, n_cells) {
  assert(is.list(checkpoint), "Priority checkpoint must be a list.")

  required_fields <- priority_checkpoint_state_fields()

  missing_fields <- setdiff(required_fields, names(checkpoint))
  if (length(missing_fields)) {
    stop(
      "Priority checkpoint is missing field(s): ",
      paste(missing_fields, collapse = ", ")
    )
  }

  stage_index <- validate_checkpoint_nonnegative_integer(
    checkpoint$stage_index,
    "checkpoint stage_index"
  )
  completed_stages <- validate_checkpoint_nonnegative_integer(
    checkpoint$completed_stages,
    "checkpoint completed_stages"
  )
  validate_priority_count(
    checkpoint$removal_step_counter,
    "checkpoint removal_step_counter"
  )

  assert(
    completed_stages == stage_index,
    "Checkpoint completed_stages must equal checkpoint stage_index."
  )

  assert(
    is.logical(checkpoint$frontier_exhausted) &&
      length(checkpoint$frontier_exhausted) == 1L,
    "Checkpoint frontier_exhausted must be one logical value."
  )

  assert(
    is.data.frame(checkpoint$patch_table),
    "Checkpoint patch_table must be a data.frame or data.table."
  )

  need_cols(
    checkpoint$patch_table,
    c("species", "patch_id", "pu_id", "patch_area_km2"),
    "checkpoint patch_table"
  )

  assert(
    is.list(checkpoint$pu_graphs_by_key),
    "Checkpoint pu_graphs_by_key must be a list."
  )

  assert(
    length(checkpoint$alive_species_count_by_cell) == n_cells,
    "Checkpoint alive_species_count_by_cell length does not match n_cells."
  )

  assert(
    is.environment(checkpoint$patch_id_by_species_env),
    "Checkpoint patch_id_by_species_env must be an environment."
  )

  assert(
    length(checkpoint$removal_order_by_cell) == n_cells,
    "Checkpoint removal_order_by_cell length does not match n_cells."
  )

  assert(
    is.list(checkpoint$removal_event_rows),
    "Checkpoint removal_event_rows must be a list."
  )

  invisible(TRUE)
}

validate_priority_checkpoint_schema <- function(
  checkpoint,
  expected_metadata = NULL,
  n_cells = NULL
) {
  assert(is.list(checkpoint), "Priority checkpoint must be a list.")

  missing_fields <- setdiff(priority_checkpoint_required_fields(), names(checkpoint))
  if (length(missing_fields)) {
    stop(
      "Priority checkpoint is missing field(s): ",
      paste(missing_fields, collapse = ", ")
    )
  }

  if (is.null(expected_metadata)) {
    expected_metadata <- checkpoint$metadata
  }

  validate_priority_checkpoint_metadata(
    checkpoint = checkpoint,
    expected_metadata = expected_metadata
  )

  if (is.null(n_cells)) {
    n_cells <- validate_priority_count(
      checkpoint$metadata$n_cells,
      "checkpoint metadata$n_cells"
    )
  }

  validate_priority_checkpoint_state(
    checkpoint = checkpoint,
    n_cells = n_cells
  )

  invisible(TRUE)
}

# ---- Runtime-index reconstruction -------------------------------------------

# Rebuild the omitted patch-to-cell lookup in sorted patch-id order. This order
# is part of deterministic resume behavior, not a new scientific computation.
rebuild_patch_cell_index_env_from_patch_ids <- function(patch_id_by_species_env) {
  assert(
    is.environment(patch_id_by_species_env),
    "patch_id_by_species_env must be an environment."
  )

  species_names <- sort(ls(envir = patch_id_by_species_env, all.names = TRUE))
  patch_cell_index_by_species_env <- new.env(parent = emptyenv())

  for (species_name in species_names) {
    species_patch_ids <- get(
      species_name,
      envir = patch_id_by_species_env,
      inherits = FALSE
    )

    occupied_cells <- which(!is.na(species_patch_ids))

    if (!length(occupied_cells)) {
      assign(species_name, NULL, envir = patch_cell_index_by_species_env)
      next
    }

    occupied_patch_ids <- species_patch_ids[occupied_cells]
    ordering <- order(occupied_patch_ids)

    assign(
      species_name,
      list(
        pid = as.integer(occupied_patch_ids[ordering]),
        cell = as.integer(occupied_cells[ordering])
      ),
      envir = patch_cell_index_by_species_env
    )
  }

  patch_cell_index_by_species_env
}

# ---- Checkpoint construction, atomic publication, and retention -------------

make_priority_checkpoint <- function(
  metadata,
  stage_index,
  completed_stages,
  frontier_exhausted,
  patch_table,
  pu_graphs_by_key,
  alive_species_count_by_cell,
  patch_id_by_species_env,
  removal_order_by_cell,
  removal_event_rows,
  removal_step_counter
) {
  pu_graphs_by_key <- priority_graph_as_list(pu_graphs_by_key)
  stage_index <- validate_checkpoint_nonnegative_integer(stage_index, "stage_index")
  completed_stages <- validate_checkpoint_nonnegative_integer(
    completed_stages,
    "completed_stages"
  )
  removal_step_counter <- validate_priority_count(
    removal_step_counter,
    "removal_step_counter"
  )
  frontier_exhausted <- validate_scalar_logical(
    frontier_exhausted,
    "frontier_exhausted"
  )
  assert(
    completed_stages == stage_index,
    "completed_stages must equal stage_index when writing a checkpoint."
  )
  validate_priority_checkpoint_metadata(list(metadata = metadata), metadata)

  list(
    metadata = metadata,
    stage_index = stage_index,
    completed_stages = completed_stages,
    frontier_exhausted = frontier_exhausted,
    patch_table = patch_table,
    pu_graphs_by_key = pu_graphs_by_key,
    alive_species_count_by_cell = alive_species_count_by_cell,
    patch_id_by_species_env = patch_id_by_species_env,
    removal_order_by_cell = removal_order_by_cell,
    removal_event_rows = removal_event_rows,
    removal_step_counter = removal_step_counter
  )
}

prune_priority_checkpoints <- function(checkpoint_dir, checkpoint_keep) {
  checkpoint_keep <- validate_priority_count(checkpoint_keep, "checkpoint_keep")

  checkpoints <- list_priority_checkpoint_files(checkpoint_dir)
  if (nrow(checkpoints) <= checkpoint_keep) {
    return(invisible(character()))
  }

  old_paths <- checkpoints$path[seq_len(nrow(checkpoints) - checkpoint_keep)]
  unlink(old_paths, force = TRUE)

  invisible(old_paths)
}

log_priority_checkpoint_saved <- function(stage_index, checkpoint_path, checkpoint_keep) {
  runtime_log_event(
    "priority_checkpoint_saved",
    stage = as.integer(stage_index),
    path = normalizePath(checkpoint_path, mustWork = FALSE),
    keep = as.integer(checkpoint_keep)
  )

  invisible(NULL)
}

write_priority_checkpoint <- function(checkpoint, checkpoint_dir, checkpoint_keep = 2L) {
  checkpoint_keep <- validate_priority_count(checkpoint_keep, "checkpoint_keep")
  validate_priority_checkpoint_schema(checkpoint)

  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(checkpoint_dir)) {
    stop("Could not create priority checkpoint directory: ", checkpoint_dir)
  }

  checkpoint_path <- priority_checkpoint_path(
    checkpoint_dir = checkpoint_dir,
    stage_index = checkpoint$stage_index
  )
  tmp_path <- tempfile(
    pattern = paste0(basename(checkpoint_path), "."),
    tmpdir = checkpoint_dir,
    fileext = ".tmp"
  )
  backup_path <- paste0(checkpoint_path, ".bak")

  on.exit({
    if (file.exists(tmp_path)) {
      unlink(tmp_path, force = TRUE)
    }
  }, add = TRUE)

  # Publish only a complete temporary serialization. An existing checkpoint is
  # retained as a backup until the replacement has been promoted successfully.
  saveRDS(checkpoint, tmp_path, compress = FALSE)
  if (!file.exists(tmp_path)) {
    stop("Failed to write temporary priority checkpoint: ", tmp_path)
  }

  if (file.exists(checkpoint_path)) {
    if (file.exists(backup_path)) {
      stop("Unresolved priority checkpoint backup exists: ", backup_path)
    }

    backup_created <- file.rename(checkpoint_path, backup_path)
    if (!isTRUE(backup_created) || !file.exists(backup_path)) {
      stop(
        "Failed to preserve existing priority checkpoint before replacement: ",
        checkpoint_path
      )
    }
  }

  renamed <- file.rename(tmp_path, checkpoint_path)
  if (!isTRUE(renamed)) {
    if (file.exists(backup_path)) {
      restored <- file.rename(backup_path, checkpoint_path)
      if (!isTRUE(restored)) {
        stop(
          "Failed to move the new checkpoint and restore the previous checkpoint. Backup retained at: ",
          backup_path
        )
      }
    }
    stop("Failed to move priority checkpoint into place: ", checkpoint_path)
  }

  if (file.exists(backup_path)) {
    unlink(backup_path, force = TRUE)
  }

  # Retention runs last, after the new canonical checkpoint is safely in place.
  prune_priority_checkpoints(
    checkpoint_dir = checkpoint_dir,
    checkpoint_keep = checkpoint_keep
  )

  log_priority_checkpoint_saved(
    stage_index = checkpoint$stage_index,
    checkpoint_path = checkpoint_path,
    checkpoint_keep = checkpoint_keep
  )

  gc(FALSE)

  invisible(checkpoint_path)
}

log_priority_resume <- function(
  checkpoint_path,
  stage_index,
  removal_step_counter,
  alive_cell_count
) {
  checkpoint_label <- if (
    is.null(checkpoint_path) ||
      !length(checkpoint_path) ||
      is.na(checkpoint_path[1L]) ||
      !nzchar(checkpoint_path[1L])
  ) {
    "in_memory"
  } else {
    normalizePath(checkpoint_path, mustWork = FALSE)
  }

  runtime_log_event(
    "priority_resume",
    checkpoint = checkpoint_label,
    stage = as.integer(stage_index),
    next_stage = as.integer(stage_index) + 1L,
    removal_steps = as.integer(removal_step_counter) - 1L,
    alive = as.integer(alive_cell_count)
  )

  invisible(NULL)
}
