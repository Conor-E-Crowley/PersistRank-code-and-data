# Reusable demographic- and persistence-model lifecycle contracts.
#
# Loaded after canonical storage paths and artifact primitives by Stages 1-8
# and the standalone application adapter. Sourcing has no side effects; explicit
# validation/state/manifest calls perform the documented reads or atomic writes.

# Install the canonical Stage 1 model manifest from the complete output set.
# Publication is the trust boundary: every declared artifact must exist here.
write_demography_model_manifest <- function(model, files) {
  files <- as.character(files)
  manifest <- list(
    schema = application_storage_schema(),
    type = "demographic_model",
    contract = list(calibration = "canonical_stage1"),
    files = lapply(files, storage_file_fingerprint, root = model$root)
  )
  storage_atomic_save_rds(
    manifest, model$manifest, "demographic-model manifest"
  )
  invisible(manifest)
}

update_persistence_model_state <- function(model, stage, files, summary = list()) {
  assert(stage %in% c("stage_2", "stage_3"),
         "Persistence-model state supports only stage_2 and stage_3.")
  state <- if (file.exists(model$run_state)) readRDS(model$run_state) else list(
    schema = application_storage_schema(), type = "persistence_model_state",
    stages = list()
  )
  assert(identical(state$schema, application_storage_schema()) &&
           identical(state$type, "persistence_model_state"),
         "Unsupported persistence-model recovery state.")
  state$stages[[stage]] <- list(
    status = "complete",
    files = lapply(
      as.character(files), storage_file_fingerprint, root = model$root
    ),
    summary = summary
  )
  storage_atomic_save_rds(state, model$run_state, "persistence-model recovery state")
  invisible(state)
}

validate_demography_model <- function(base = project_paths()) {
  paths <- demography_model_paths(base)
  need_file(paths$manifest, "demographic-model manifest")
  manifest <- readRDS(paths$manifest)
  assert(identical(manifest$schema, application_storage_schema()) &&
           identical(manifest$type, "demographic_model") &&
           is.list(manifest$files) && length(manifest$files) > 0L,
         "Unsupported or incomplete demographic-model manifest.")
  for (record in manifest$files) {
    path <- storage_record_path(record, paths$root, must_work = TRUE)
    assert(identical(unname(tools::md5sum(path)), record$md5),
           paste0("Demographic-model artifact checksum mismatch: ", path))
  }
  manifest
}

persistence_model_contract <- function(config, stage2_config = NULL) {
  grid_signature <- if (is.null(stage2_config)) NULL else stage2_config$grid$signature
  model_taxa <- if (is.null(stage2_config)) "both" else stage2_config$taxa
  assert(
    identical(model_taxa, "both"),
    "A reusable persistence model must contain both mammal and bird branches."
  )
  list(
    taxa = model_taxa,
    persistence_horizon_years = as.integer(config$contract$persistence_horizon_years),
    quasi_extinction_abundance = as.integer(config$contract$quasi_extinction_abundance),
    curves = persistence_curves(),
    grid_signature = grid_signature,
    loess_span = as.numeric(config$loess_span)
  )
}

validate_persistence_model <- function(config, require_files = TRUE) {
  need_file(config$model$manifest, "persistence-model manifest")
  manifest <- readRDS(config$model$manifest)
  assert(identical(manifest$schema, application_storage_schema()) &&
           identical(manifest$type, "persistence_model"),
         "Unsupported persistence-model manifest.")
  expected <- persistence_model_contract(config, NULL)
  for (field in c("taxa", "persistence_horizon_years",
                  "quasi_extinction_abundance", "curves")) {
    assert(identical(manifest$contract[[field]], expected[[field]]), paste0(
      "Persistence model has a different ", field, ". Start from stage_2."
    ))
  }
  if (isTRUE(require_files)) {
    for (record in manifest$files) {
      path <- storage_record_path(record, config$model$root, must_work = TRUE)
      assert(identical(unname(tools::md5sum(path)), record$md5),
             paste0("Persistence-model artifact checksum mismatch: ", path))
    }
  }
  manifest
}

write_persistence_model_manifest_from_artifacts <- function(
  model, persistence_horizon_years, quasi_extinction_abundance, loess_span = 0.5
) {
  points <- c(
    mammals = file.path(model$stage2, "persistence_points_mammals.csv"),
    birds = file.path(model$stage2, "persistence_points_birds.csv")
  )
  points <- points[file.exists(points)]
  model_file <- file.path(model$stage3, "persistence_curve_models.rds")
  need_file(model_file, "Stage 3 persistence model")
  assert(length(points) == 2L, paste0(
    "A reusable persistence model requires completed mammal and bird Stage 2 point tables."
  ))
  signatures <- unique(unlist(lapply(points, function(path) {
    unique(data.table::fread(path, select = "grid_signature")$grid_signature)
  }), use.names = FALSE))
  assert(length(signatures) == 1L, "Stage 2 point tables have inconsistent grid signatures.")
  manifest <- list(
    schema = application_storage_schema(), type = "persistence_model",
    contract = list(
      taxa = "both",
      persistence_horizon_years = as.integer(persistence_horizon_years),
      quasi_extinction_abundance = as.integer(quasi_extinction_abundance),
      curves = persistence_curves(), grid_signature = signatures[[1L]],
      loess_span = as.numeric(loess_span)
    ),
    files = lapply(
      c(points, model = model_file), storage_file_fingerprint, root = model$root
    ),
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
  if (file.exists(model$manifest)) {
    existing <- readRDS(model$manifest)
    assert(identical(existing$contract, manifest$contract) &&
             storage_fingerprint_lists_equal(existing$files, manifest$files), paste0(
      "A different persistence-model manifest already exists at ", model$root,
      ". Remove that model directory only when deliberately rebuilding it."
    ))
    return(existing)
  }
  storage_atomic_save_rds(manifest, model$manifest, "persistence-model manifest")
  manifest
}
