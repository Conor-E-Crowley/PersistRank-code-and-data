# Application recovery lifecycle shared by standalone Rmds and Stage 8.
#
# Workflows load this definition-only module after application context/state.
# Boundary construction is purely in memory. Commits retain the existing
# artifact hashing, state schema, and atomic manifest transactions; this module
# never loads scientific, spatial, reporting, or compiled runtime code.

application_boundary <- function(stage, status, outputs, summary = list(),
                                 upstream = character(), curve = NULL) {
  stage <- validate_scalar_choice(stage, application_boundary_stages(), "stage")
  status <- validate_scalar_choice(
    status, c("complete", "reused", "not_applicable"), "status"
  )
  curve_stages <- c("stage_6", "stage_7_1")
  if (stage %in% curve_stages) {
    assert(!is.null(curve), paste0("curve is required for ", stage, "."))
    curve <- unname(as.character(validate_persistence_curve(curve, "curve")))
  } else {
    assert(is.null(curve), paste0("curve is not accepted for ", stage, "."))
  }
  normalize_paths <- function(x, label, allow_empty = TRUE) {
    assert(is.character(x), paste0(label, " must be a character vector."))
    x <- trimws(x)
    assert(!anyNA(x) && all(nzchar(x)), paste0(
      label, " must not contain missing or blank paths."
    ))
    x <- unique(path.expand(x))
    assert(isTRUE(allow_empty) || length(x), paste0(
      label, " must contain at least one path."
    ))
    x
  }
  outputs <- normalize_paths(
    outputs, "outputs", allow_empty = identical(status, "not_applicable")
  )
  upstream <- normalize_paths(upstream, "upstream")
  assert(is.list(summary), "summary must be a list.")
  list(
    key = application_stage_key(stage, curve), status = status,
    outputs = outputs, upstream = upstream, summary = summary
  )
}

# Validate and atomically persist one prepared boundary. Supplying state avoids
# rereading the recovery file, which is important for Stage 8's ordered run.
application_commit_boundary <- function(context, boundary, state = NULL) {
  assert(is.list(boundary) && identical(
    names(boundary), c("key", "status", "outputs", "upstream", "summary")
  ), "boundary must be an application_boundary() result.")
  if (is.null(state)) state <- read_application_state(context)
  application_finish_boundary(
    state, boundary$key, boundary$status, boundary$outputs,
    boundary$summary, boundary$upstream, context
  )
}

application_stage4_boundary <- function(result, context) {
  application_boundary(
    "stage_4", "complete", result$outputs,
    result[setdiff(names(result), c("outputs", "habitat_rows"))],
    c(context$model$manifest, context$scenario$application_inputs)
  )
}

application_stage5_boundary <- function(config, result) {
  application_boundary(
    "stage_5",
    if (identical(result$action, "final_reused")) "reused" else "complete",
    c(config$paths$patch_dir, config$paths$patch_lookup_rds,
      config$paths$connectivity_rds, config$paths$metadata_rds),
    list(
      action = result$action, status = result$status,
      elapsed_seconds = result$elapsed_seconds
    ),
    c(config$paths$species_csv, config$paths$landcover_tif)
  )
}

application_stage51_boundary <- function(config, result) {
  application_boundary(
    "stage_5_1", result$status, result$figure %||% character(),
    result[setdiff(names(result), "figure")],
    c(config$paths$species_csv, config$paths$landcover_tif)
  )
}

application_stage52_boundary <- function(config, result) {
  if (!identical(config$mode, "inputs")) return(NULL)
  application_boundary(
    "stage_5_2", "complete",
    c(result$binary_paths, result$feature_list, result$settings_file,
      result$method_dirs),
    list(
      species = length(result$species),
      rankmaps_present = result$rankmaps_present
    ),
    c(config$paths$patch_lookup_rds, config$paths$settings_template)
  )
}

application_stage53_boundary <- function(config, result) {
  application_boundary(
    "stage_5_3", "complete", config$paths$figure,
    list(
      coverage = result$coverage, pu_summary = result$pu_summary,
      species_summary = result$species_summary,
      correlations = result$correlations, curve = result$curve,
      figure = config$paths$figure
    ),
    c(config$paths$species_csv, config$paths$patch_lookup_rds)
  )
}

application_stage6_initialization_boundary <- function(config, result) {
  application_boundary(
    "stage_6_initialization",
    "complete",
    config$paths$initialization_bundle,
    result[setdiff(names(result), "initialization_bundle_path")],
    c(config$paths$species_table, config$paths$patch_lookup,
      config$paths$connectivity, config$paths$patch_dir)
  )
}

application_stage6_run_boundary <- function(config, result) {
  pipeline <- result$priority_pipeline_result
  if (!isTRUE(pipeline$frontier_exhausted)) return(NULL)
  application_boundary(
    "stage_6", "complete", stage6_scientific_outputs(config$paths),
    list(
      completed_stages = pipeline$completed_stages,
      frontier_exhausted = TRUE,
      initialization_created_at = result$initialization_created_at,
      coefficient_identity = result$coefficient_identity
    ),
    character(), config$curve
  )
}

application_stage71_boundary <- function(config, result) {
  summary <- list(
    surface = result$surface, outputs = unname(result$outputs),
    elapsed = result$elapsed
  )
  application_boundary(
    "stage_7_1", "complete", unname(result$outputs), summary,
    stage6_scientific_outputs(config$paths),
    config$curve
  )
}

# Validate completed Stage 4--5 records and atomically publish their unchanged
# threshold-scenario handoff. Optional supplied objects avoid duplicate reads.
application_finalize_scenario <- function(
  context, state = NULL, model_manifest = NULL, application_manifest = NULL,
  require_complete = TRUE
) {
  if (is.null(state)) state <- read_application_state(context)
  required <- application_scenario_core_stages()
  missing <- required[!required %in% names(state$stages)]
  present <- required[required %in% names(state$stages)]
  invalid <- present[!vapply(
    state$stages[present], application_record_valid, logical(1L),
    root = context$base_paths$root
  )]
  message <- paste0(
    "Cannot finalize the Stage 4–5 scenario.",
    if (length(missing)) paste0(
      " Missing stages: ", paste(missing, collapse = ", "), "."
    ) else "",
    if (length(invalid)) paste0(
      " Changed stages: ", paste(invalid, collapse = ", "), "."
    ) else "",
    " Run the missing numbered Rmds first."
  )
  if (length(missing) || length(invalid)) {
    if (isTRUE(require_complete)) stop(message, call. = FALSE)
    log_msg("APPLICATION SCENARIO | status=incomplete | ", message)
    return(invisible(NULL))
  }
  if (is.null(application_manifest)) {
    application_manifest <- readRDS(context$scenario$application_manifest)
  }
  if (is.null(model_manifest)) model_manifest <- validate_persistence_model(context)
  manifest <- application_scenario_manifest(
    context, state, model_manifest, application_manifest
  )
  application_write_scenario_manifest(context$scenario$manifest, manifest)
}

# Return the compatible priority-run manifest or a new in-memory envelope.
application_priority_run_manifest <- function(context) {
  manifest <- list(
    schema = application_storage_schema(), type = "priority_run",
    contract = application_run_contract(context), curves = list()
  )
  if (!file.exists(context$scenario$run_manifest)) return(manifest)
  existing <- readRDS(context$scenario$run_manifest)
  differences <- application_manifest_differences(
    existing$contract, manifest$contract
  )
  assert(!length(differences), paste0(
    "Existing run manifest does not match this removal schedule."
  ))
  existing
}

# Report whether one canonical curve is scientifically complete in a validated
# priority-run manifest. This read-only query keeps curve-entry ownership out of
# Stage 8 sequencing and report adapters.
application_curve_manifest_complete <- function(manifest, curve) {
  curve <- unname(as.character(validate_persistence_curve(curve, "curve")))
  entry <- manifest$curves[[curve]]
  is.list(entry) && identical(entry$status, "complete")
}

# Validate the shared initialization and curve-scoped execution records, then
# atomically mark the curve complete in the priority-run manifest.
application_mark_curve_complete <- function(
  context, curve, state = NULL, manifest = NULL
) {
  curve <- unname(as.character(validate_persistence_curve(curve, "curve")))
  if (is.null(state)) state <- read_application_state(context)
  keys <- c(
    application_scenario_shared_stages(),
    vapply(
      application_curve_core_stages(), application_stage_key,
      character(1L), curve = curve
    )
  )
  missing <- keys[!keys %in% names(state$stages)]
  curve_keys <- setdiff(keys, application_scenario_shared_stages())
  present_curve_keys <- curve_keys[curve_keys %in% names(state$stages)]
  # The shared record is checksum-validated once before the multi-curve loop.
  # Rehashing the same large immutable file while committing every curve would
  # add repeated I/O without strengthening that already-established boundary.
  invalid <- present_curve_keys[!vapply(
    state$stages[present_curve_keys], application_record_valid, logical(1L),
    root = context$base_paths$root
  )]
  assert(!length(missing) && !length(invalid), paste0(
    "Cannot mark curve ", curve, " complete.",
    if (length(missing)) paste0(
      " Missing records: ", paste(missing, collapse = ", "), "."
    ) else "",
    if (length(invalid)) paste0(
      " Changed records: ", paste(invalid, collapse = ", "), "."
    ) else ""
  ))
  if (is.null(manifest)) manifest <- application_priority_run_manifest(context)
  initialization <- state$stages[[application_scenario_shared_stages()]]
  execution <- state$stages[[application_stage_key("stage_6", curve)]]
  initialization_record <- initialization$outputs[[1L]]
  coefficient <- execution$summary$coefficient_identity
  assert(
    length(initialization$outputs) == 1L && is.list(coefficient) &&
      identical(coefficient$curve, curve),
    paste0("Incomplete Stage 6 identity records for curve ", curve, ".")
  )
  manifest$curves[[curve]] <- list(
    status = "complete",
    path = storage_relative_path(
      application_curve_paths(context$paths, curve)$out_curve,
      context$scenario$run_root
    ),
    initialization = list(
      path = storage_relative_path(
        context$scenario$priority_initialization,
        context$base_paths$root
      ),
      md5 = initialization_record$md5
    ),
    coefficient_identity = coefficient
  )
  storage_atomic_save_rds(
    manifest, context$scenario$run_manifest, "priority-run manifest"
  )
  manifest
}
