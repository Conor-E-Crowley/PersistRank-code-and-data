# Stage 2 recovery helpers.
#
# A checkpoint contains only the active trait/curve identity and completed
# K-to-quantile evaluations. Trajectory random numbers remain owned by the C++
# CRN context and are deterministically reconstructed after a restart.


stage2_checkpoint_schema <- function() 1L

stage2_checkpoint_signature <- function(sim, grid, metadata) {
  paste(
    "schema", stage2_checkpoint_schema(),
    "years", sim$years,
    "threshold", sim$quasi_extinction_abundance,
    "cap", format(sim$cap_factor, digits = 17),
    "buffer", format(sim$r_buffer, digits = 17),
    "draws", sim$n_draws,
    "reps", sim$reps,
    "chunk", sim$chunk_size,
    "seed", sim$base_seed,
    "uncertainty", sim$demographic_uncertainty,
    "grid", grid$signature,
    "rm", metadata$posterior_rm_md5,
    "sigma", metadata$posterior_sigma_md5,
    "kernel", metadata$simulator_contract,
    sep = "|"
  )
}

checkpoint_cache_as_list <- function(cache) {
  keys <- sort(ls(cache, all.names = TRUE))
  stats::setNames(lapply(keys, function(key) cache[[key]]), keys)
}

write_stage2_checkpoint <- function(
  path,
  taxon,
  trait_key,
  curve,
  cache,
  sim,
  grid,
  metadata
) {
  state <- list(
    schema_version = stage2_checkpoint_schema(),
    saved_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    signature = stage2_checkpoint_signature(sim, grid, metadata),
    taxon = taxon,
    trait_key = trait_key,
    curve = curve,
    evaluations = checkpoint_cache_as_list(cache)
  )
  ensure_writable_dir(dirname(path), "Stage 2 checkpoint directory")
  temporary <- tempfile(".stage2_checkpoint_", tmpdir = dirname(path), fileext = ".rds")
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(state, temporary, version = 3L)
  check <- readRDS(temporary)
  assert(
    identical(check$signature, state$signature) &&
      identical(check$evaluations, state$evaluations),
    "Temporary Stage 2 checkpoint validation failed."
  )
  project_file_set_transaction(
    temporary, path, overwrite = TRUE, label = "Stage 2 checkpoint"
  )
  invisible(state)
}

read_stage2_checkpoint <- function(path, sim, grid, metadata) {
  if (!file.exists(path)) return(NULL)
  state <- tryCatch(readRDS(path), error = function(e) {
    stop("Stage 2 checkpoint cannot be read: ", conditionMessage(e), call. = FALSE)
  })
  assert(
    is.list(state) &&
      identical(state$schema_version, stage2_checkpoint_schema()),
    "Stage 2 checkpoint has an unsupported schema."
  )
  expected <- stage2_checkpoint_signature(sim, grid, metadata)
  assert(
    identical(state$signature, expected),
    paste0(
      "Stage 2 checkpoint is incompatible with the current scientific ",
      "settings or source posteriors. Use mode: 'restart'."
    )
  )
  assert(
    is.list(state$evaluations) &&
      all(vapply(
        state$evaluations,
        function(x) is.numeric(x) && all(is.finite(x)),
        logical(1)
      )),
    "Stage 2 checkpoint contains invalid cached evaluations."
  )
  state
}

restore_stage2_evaluation_cache <- function(cache, checkpoint, trait_key, curve) {
  if (is.null(checkpoint) ||
      !identical(checkpoint$trait_key, trait_key) ||
      !identical(checkpoint$curve, curve)) {
    return(0L)
  }
  for (key in names(checkpoint$evaluations)) {
    cache[[key]] <- checkpoint$evaluations[[key]]
  }
  length(checkpoint$evaluations)
}

remove_stage2_checkpoint <- function(path) {
  if (file.exists(path)) {
    assert(file.remove(path), paste0("Could not remove completed checkpoint: ", path))
  }
  invisible(TRUE)
}

inspect_stage2_artifact <- function(
  path,
  trait_table,
  config,
  metadata,
  idx_col,
  value_col,
  extra_cols = character()
) {
  if (!file.exists(path)) {
    return(data.frame(
      path = path, exists = FALSE, compatible = NA,
      completed_blocks = 0L, remaining_blocks = NA_integer_,
      detail = "missing", stringsAsFactors = FALSE
    ))
  }
  result <- tryCatch(
    inspect_existing_persist_points(
      path = path,
      trait_table = trait_table,
      sim = config$sim,
      grid = config$grid,
      metadata = metadata,
      idx_col = idx_col,
      value_col = value_col,
      extra_cols = extra_cols,
      label = basename(path),
      require_complete = FALSE
    ),
    error = identity
  )
  if (inherits(result, "error")) {
    return(data.frame(
      path = path, exists = TRUE, compatible = FALSE,
      completed_blocks = NA_integer_, remaining_blocks = NA_integer_,
      detail = conditionMessage(result), stringsAsFactors = FALSE
    ))
  }
  data.frame(
    path = path,
    exists = TRUE,
    compatible = TRUE,
    completed_blocks = length(result$completed_keys),
    remaining_blocks = length(result$missing_keys),
    detail = if (length(result$incomplete_keys)) {
      paste(length(result$incomplete_keys), "partial block(s) will be dropped")
    } else {
      "valid"
    },
    stringsAsFactors = FALSE
  )
}
