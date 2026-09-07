# Stage 6 public configuration and read-only preflight.
#
# Scientific selectors and Stage 4–5 artifacts determine the shared
# initialization. Curve and removal schedule determine only execution outputs.
# Logging, checkpoint retention, and max_stages control execution only.


stage6_artifact_paths <- function(project, curve = NULL) {
  shared <- list(
    scenario_root = project$scenario_root,
    species_table = file.path(project$clean, "species_table.csv"),
    patch_dir = project$patches,
    patch_lookup = project$patch_lookup %||%
      file.path(project$clean, "all_patch_lookup.rds"),
    connectivity = project$connectivity %||%
      file.path(project$clean, "all_connectivity.rds"),
    initialization_bundle = project$priority_initialization,
    initialization_runtime_log = file.path(
      project$scenario_root, "Logs", "stage6_initialization.log"
    )
  )
  if (is.null(curve)) return(shared)

  run <- application_curve_paths(project = project, curve = curve)
  c(shared, list(
    output_root = project$priority_runs,
    rcpp_cache = project$rcpp_cache,
    checkpoint_dir = file.path(run$out_curve, "checkpoints"),
    runtime_log = file.path(
      project$priority_run, "Logs", paste0("stage6_", curve, ".log")
    ),
    replacement_backup = directory_replacement_backup_path(run$out_curve),
    patch_lookup_output_dir = run$patch_lookup_dir,
    rankmap = file.path(run$out_curve, "rankmap.tif")
  ), run)
}

stage6_config <- function(params, paths, preflight = TRUE,
                          contract = canonical_analysis_contract()) {
  assert(is.list(params), "Stage 6 parameters must be supplied as a list.")
  assert(is.list(paths), "paths must be a project_paths() result.")
  contract <- validate_analysis_contract(contract, "Stage 6 abundance contract")

  mode <- validate_scalar_choice(
    params$mode,
    c("inspect", "initialize", "run", "resume"),
    "params$mode"
  )
  taxon_selection <- validate_priority_public_taxa(params$taxa, "params$taxa")
  sdm <- validate_priority_sdm_tag(params$sdm, "params$sdm")

  # Shared initialization has no curve or execution schedule. Returning here
  # prevents irrelevant run/checkpoint options from becoming operational
  # preconditions for a scenario-owned artifact.
  if (identical(mode, "initialize")) {
    artifact_paths <- stage6_artifact_paths(project = paths)
    config <- list(
      mode = mode,
      curve = NULL,
      taxa = taxon_selection$taxa,
      taxa_tag = taxon_selection$taxa_tag,
      sdm = sdm,
      contract = contract,
      paths = artifact_paths
    )
    if (!isTRUE(preflight)) return(config)

    need_file(artifact_paths$species_table, "Stage 4 species table")
    need_dir(artifact_paths$patch_dir, "Stage 5 patch raster directory")
    need_file(artifact_paths$patch_lookup, "Stage 5 patch lookup")
    need_file(artifact_paths$connectivity, "Stage 5 connectivity")
    ensure_writable_dir(
      dirname(artifact_paths$initialization_bundle),
      "Stage 6 initialization directory"
    )
    return(config)
  }

  curve <- validate_persistence_curve(params$curve, "params$curve")
  cells <- validate_priority_count(
    params$cells_to_remove_per_iteration,
    "params$cells_to_remove_per_iteration"
  )
  iterations <- validate_priority_count(
    params$pruning_iterations_per_stage,
    "params$pruning_iterations_per_stage"
  )
  max_stages <- validate_max_stages(params$max_stages, "params$max_stages")

  ecology_every <- params$ecology_log_every_iterations
  if (!is.null(ecology_every)) {
    ecology_every <- validate_priority_count(
      ecology_every,
      "params$ecology_log_every_iterations"
    )
  }
  checkpoint_every <- params$checkpoint_every_stages
  if (!is.null(checkpoint_every)) {
    checkpoint_every <- validate_priority_count(
      checkpoint_every,
      "params$checkpoint_every_stages"
    )
  }
  checkpoint_keep <- validate_priority_count(
    params$checkpoint_keep,
    "params$checkpoint_keep"
  )
  resume_stage <- params$resume_stage
  if (!is.null(resume_stage)) {
    resume_stage <- validate_priority_count(resume_stage, "params$resume_stage")
  }
  assert(
    identical(mode, "resume") || is.null(resume_stage),
    "params$resume_stage is used only with params$mode = 'resume'."
  )

  artifact_paths <- stage6_artifact_paths(project = paths, curve = curve)

  config <- list(
    mode = mode,
    curve = unname(as.character(curve)),
    taxa = taxon_selection$taxa,
    taxa_tag = taxon_selection$taxa_tag,
    sdm = sdm,
    cells_to_remove_per_iteration = cells,
    pruning_iterations_per_stage = iterations,
    max_stages = max_stages,
    ecology_log_every_iterations = ecology_every,
    checkpoint_every_stages = checkpoint_every,
    checkpoint_keep = checkpoint_keep,
    resume_stage = resume_stage,
    contract = contract,
    paths = artifact_paths
  )

  if (!isTRUE(preflight) || identical(mode, "inspect")) return(config)

  need_file(artifact_paths$initialization_bundle, "Stage 6 shared initialization")
  need_dir(artifact_paths$patch_dir, "Stage 5 patch raster directory")
  ensure_writable_dir(artifact_paths$output_root, "Stage 6 run-output root")
  assert(
    !file.exists(artifact_paths$replacement_backup),
    paste0(
      "An unresolved Stage 6 replacement backup exists:\n",
      artifact_paths$replacement_backup,
      "\nResolve it before starting or resuming this run."
    )
  )

  if (identical(mode, "resume")) {
    need_dir(artifact_paths$out_curve, "Stage 6 run directory to resume")
    if (is.null(resume_stage)) {
      checkpoints <- if (dir.exists(artifact_paths$checkpoint_dir)) {
        list.files(
          artifact_paths$checkpoint_dir,
          pattern = "^priority_checkpoint_stage_[0-9]+[.]rds$",
          full.names = TRUE
        )
      } else {
        character()
      }
      assert(
        length(checkpoints) > 0L,
        paste0(
          "No Stage 6 checkpoints are available in:\n",
          artifact_paths$checkpoint_dir
        )
      )
    } else {
      need_file(
        file.path(
          artifact_paths$checkpoint_dir,
          sprintf("priority_checkpoint_stage_%04d.rds", resume_stage)
        ),
        "Stage 6 resume checkpoint"
      )
    }
  }

  config
}

# Translate public Stage 6 settings plus manifest-authoritative application
# selectors into the existing low-level Stage 6 configuration.
stage6_application_config <- function(params, context, preflight = TRUE) {
  owned <- if (identical(params$mode, "initialize")) {
    "mode"
  } else {
    c(
      "mode", "curve", "cells_to_remove_per_iteration",
      "pruning_iterations_per_stage", "max_stages",
      "ecology_log_every_iterations", "checkpoint_every_stages",
      "checkpoint_keep", "resume_stage"
    )
  }
  stage_params <- params[owned]
  stage_params$taxa <- context$taxa
  stage_params$sdm <- context$sdm
  stage6_config(stage_params, context$paths, preflight, context$contract)
}

describe_stage6_config <- function(config) {
  log_msg(
    "CONFIG | mode=", config$mode,
    " curve=", config$curve %||% "all (initialization)",
    " taxa=", config$taxa,
    if (identical(config$taxa, "both")) {
      paste0(" (artifact tag ", config$taxa_tag, ")")
    } else {
      ""
    },
    " sdm=", config$sdm
  )
  log_msg("INITIALIZATION | ", config$paths$initialization_bundle)
  if (!is.null(config$curve)) {
    log_msg(
      "SCHEDULE | cells/iteration=", config$cells_to_remove_per_iteration,
      " iterations/stage=", config$pruning_iterations_per_stage,
      " max_stages=", as.character(config$max_stages),
      " ecology_log_every=", config$ecology_log_every_iterations %||% "disabled",
      " checkpoint_every=", config$checkpoint_every_stages %||% "disabled",
      " checkpoint_keep=", config$checkpoint_keep
    )
    log_msg("RUN | id=", config$paths$run_id, " path=", config$paths$out_curve)
  }
  invisible(config)
}
# Return the complete durable output set that identifies exhausted Stage 6
# science. This is pure path selection: it performs no existence or checksum
# work and preserves the established artifact order.
stage6_scientific_outputs <- function(paths) {
  c(
    paths$patch_lookup_output_dir,
    paths$removal_events,
    paths$removal_order,
    paths$rankmap
  )
}
