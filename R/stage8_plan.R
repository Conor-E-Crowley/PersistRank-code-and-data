# Read-only Stage 8 handoff discovery and execution planning.
#
# Stage 8 loads this definition-only module after lightweight configuration,
# storage, context, and state contracts. Selection is literal: omitted stages
# are never added implicitly, and only their required durable handoffs are
# checked. Planning loads no scientific, spatial, or native runtime.

stage8_model_files <- function(config) {
  model_files <- if (file.exists(config$model$manifest)) {
    manifest <- tryCatch(readRDS(config$model$manifest), error = function(e) NULL)
    if (is.list(manifest$files)) {
      vapply(
        manifest$files, storage_record_path, character(1L),
        root = config$model$root
      )
    } else character()
  } else character()
  unique(c(
    config$model$manifest,
    file.path(config$model$stage3, "persistence_curve_models.rds"),
    model_files
  ))
}

stage8_demography_files <- function(config) {
  paths <- demography_model_paths(config$base_paths)
  c(
    paths$manifest,
    file.path(paths$outputs, c(
      "mammal_mass_grid.csv", "bird_generation_length_grid.csv",
      "mammal_growth_posterior.csv",
      "mammal_environmental_variation_posterior.csv",
      "bird_growth_posterior.csv",
      "bird_environmental_variation_posterior.csv"
    ))
  )
}

stage8_stage2_point_files <- function(config) {
  file.path(config$model$stage2, c(
    "persistence_points_mammals.csv", "persistence_points_birds.csv"
  ))
}

# Return only prerequisites whose producing stages were deliberately omitted.
# Selected upstream stages satisfy downstream requirements during this run.
stage8_required_handoffs <- function(config) {
  selected <- config$stages
  required <- character()
  add <- function(paths) required <<- c(required, paths)
  # Several selected stages can need the same model handoff. Resolve it lazily
  # and at most once so unrelated selections do not read its manifest.
  model_files <- NULL
  get_model_files <- function() {
    if (is.null(model_files)) model_files <<- stage8_model_files(config)
    model_files
  }

  if ("stage_2" %in% selected) add(stage8_demography_files(config))
  if ("stage_3" %in% selected && !"stage_2" %in% selected) {
    add(stage8_stage2_point_files(config))
  }
  if ("stage_4" %in% selected && !"stage_3" %in% selected) {
    add(get_model_files())
  }

  application_stages <- c(
    "stage_4", "stage_5", "stage_5_1", "stage_5_2", "stage_5_3",
    "stage_6", "stage_7_1"
  )
  if (any(application_stages %in% selected) &&
      !"stage_3" %in% selected) {
    add(get_model_files())
  }

  if ("stage_4" %in% selected) {
    add(vapply(config$inputs, `[[`, character(1L), "path"))
    add(config$sdm_index$path %||% character())
    add(config$study_area$file %||% character())
  }

  if ("stage_5" %in% selected && !"stage_4" %in% selected) {
    add(c(
      config$scenario$application_manifest,
      file.path(config$paths$clean, "species_table.csv"),
      config$inputs$landcover$path,
      config$study_area$file %||% character()
    ))
  }

  if ("stage_5_1" %in% selected && !"stage_5" %in% selected) {
    add(c(
      config$scenario$manifest,
      file.path(config$paths$clean, "species_table.csv"),
      config$inputs$landcover$path,
      config$study_area$file %||% character()
    ))
  }
  if ("stage_5_2" %in% selected && !"stage_5" %in% selected) {
    add(c(
      config$scenario$manifest, config$paths$patches,
      config$paths$patch_lookup, config$zonation_settings_template
    ))
  }
  if ("stage_5_3" %in% selected && !"stage_5" %in% selected) {
    add(c(
      config$scenario$manifest,
      file.path(config$paths$clean, "species_table.csv"),
      config$paths$patch_lookup
    ))
  }

  if ("stage_6" %in% selected && !"stage_5" %in% selected) {
    add(c(
      get_model_files(), config$scenario$application_manifest,
      config$scenario$manifest,
      file.path(config$paths$clean, "species_table.csv"),
      config$paths$patches, config$paths$patch_lookup,
      config$paths$connectivity
    ))
  }

  if ("stage_7_1" %in% selected && !"stage_6" %in% selected) {
    add(c(
      get_model_files(), config$scenario$application_manifest,
      config$scenario$manifest, config$scenario$run_manifest,
      config$paths$priority_initialization
    ))
    for (curve in config$priority_curves) {
      add(stage6_scientific_outputs(stage6_artifact_paths(config$paths, curve)))
    }
  }

  unique(required[!is.na(required) & nzchar(required)])
}

stage8_plan <- function(config, state = NULL) {
  required <- stage8_required_handoffs(config)
  exists <- file.exists(required) | dir.exists(required)
  missing <- required[!exists]
  if (is.null(state)) state <- read_application_state(config)
  # Active execution validates each selected recovery record immediately before
  # reuse. Avoid hashing the same large outputs and upstream directories once
  # for planning and again for execution. Inspect mode has no execution pass,
  # so its plan performs the authoritative checksum validation here.
  inspect_mode <- identical(config$mode, "inspect")
  record_reusable <- function(record) {
    !is.null(record) && inspect_mode &&
      application_record_valid(record, config$base_paths$root)
  }
  stage4_record_ready <- !is.null(state$stages$stage_4) &&
    (!inspect_mode || record_reusable(state$stages$stage_4))
  missing_state <- if (
    "stage_5" %in% config$stages && !"stage_4" %in% config$stages &&
      !stage4_record_ready
  ) "valid Stage 4 recovery record" else character()
  rows <- list()
  add <- function(role, block, stage, artifact, action, reason,
                  exists = file.exists(artifact) || dir.exists(artifact)) {
    rows[[length(rows) + 1L]] <<- data.frame(
      role = role, block = block, stage = stage, source = "project",
      artifact = artifact, exists = exists,
      compatible = !identical(action, "blocked"), action = action,
      reason = reason, stringsAsFactors = FALSE
    )
  }

  for (i in seq_along(required)) {
    add(
      "prerequisite", "handoff", "handoff", required[[i]],
      if (exists[[i]]) "reuse" else "blocked",
      if (exists[[i]]) "omitted producer handoff is present" else
        "omitted producer handoff is missing",
      exists[[i]]
    )
  }
  if (length(missing_state)) {
    add(
      "prerequisite", "handoff", "stage_4", config$scenario$run_state,
      "blocked", missing_state, exists = FALSE
    )
  }

  representative <- list(
    stage_2 = config$model$stage2,
    stage_3 = config$model$manifest,
    stage_4 = file.path(config$paths$clean, "species_table.csv"),
    stage_5 = config$paths$patch_lookup,
    stage_5_1 = file.path(config$paths$spatial_figures, "spatial_process.png"),
    stage_5_2 = config$paths$zonation_feature_list,
    stage_5_3 = file.path(config$paths$spatial_figures, "initial_persistence.png")
  )
  for (stage in setdiff(stage8_stage_order(), c("stage_6", "stage_7_1"))) {
    selected <- stage %in% config$stages
    artifact <- representative[[stage]]
    record <- state$stages[[stage]]
    reusable <- record_reusable(record)
    present <- file.exists(artifact) || dir.exists(artifact)
    add(
      if (selected) "selected" else "omitted", "stage", stage, artifact,
      if (!selected) "skip" else if (reusable) "reuse" else if (present) {
        "resume"
      } else "run",
      if (!selected) "not selected" else if (reusable) {
        "recorded artifacts are valid"
      } else if (present) "existing work requires validation" else
        "selected output is absent",
      present
    )
  }

  if (!length(config$priority_curves)) {
    add(
      "omitted", "stage_6", "stage_6", config$paths$priority_runs,
      "skip", "not selected"
    )
    add(
      "omitted", "stage_7_1", "stage_7_1", config$paths$priority_figures,
      "skip", "not selected"
    )
  }
  for (curve in config$priority_curves) {
    curve_paths <- application_curve_paths(config$paths, curve)
    run_key <- application_stage_key("stage_6", curve)
    report_key <- application_stage_key("stage_7_1", curve)
    run_record <- state$stages[[run_key]]
    run_reusable <- record_reusable(run_record)
    run_selected <- "stage_6" %in% config$stages
    add(
      if (run_selected) "selected" else "omitted", "stage_6", run_key,
      curve_paths$out_curve,
      if (!run_selected) "skip" else if (run_reusable) "reuse" else if (
        dir.exists(curve_paths$out_curve)
      ) "resume" else "run",
      if (!run_selected) "not selected" else "curve is independently recoverable"
    )

    report_record <- state$stages[[report_key]]
    report_reusable <- record_reusable(report_record)
    report_selected <- "stage_7_1" %in% config$stages
    report_artifact <- file.path(
      config$paths$priority_figures, curve, "cell_removal_order.png"
    )
    add(
      if (report_selected) "selected" else "omitted", "stage_7_1",
      report_key, report_artifact,
      if (!report_selected) "skip" else if (report_reusable) "reuse" else if (
        run_selected && !run_reusable
      ) "deferred" else "run",
      if (!report_selected) "not selected" else if (
        run_selected && !run_reusable
      ) "awaiting selected Stage 6 curve" else "report is independently recoverable"
    )
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  attr(out, "missing") <- c(missing, missing_state)
  out
}

stage8_assert_plan_ready <- function(plan) {
  missing <- attr(plan, "missing") %||% character()
  assert(!length(missing), paste0(
    "Stage 8 cannot start because omitted-stage handoffs are missing:\n",
    paste(missing, collapse = "\n")
  ))
  invisible(TRUE)
}
