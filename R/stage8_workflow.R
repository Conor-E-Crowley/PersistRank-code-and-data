# Selective Stage 8 sequencing for reusable models and named applications.
#
# Sourcing defines lightweight configuration, planning, and execution helpers
# only. The selected stage vector is literal; omitted prerequisites are reused
# and validated but never run. Scientific workflows and expensive runtimes are
# loaded only when their selected execution boundary is reached.

if (!exists("project_source", mode = "function")) {
  source(file.path("R", "project_loader.R"))
}
project_source(c(
  # Shared contracts and application architecture.
  "R/project_utils.R", "R/project_artifacts.R", "R/project_rmd.R",
  "R/project_runtime_log.R", "R/project_transactions.R", "R/project_paths.R",
  "R/project_config.R",
  "R/demographic_contract.R", "R/analysis_contract.R",
  "R/priority_run_config.R", "R/stage2_config.R",
  "R/species_table_config.R", "R/patch_contract.R", "R/patch_config.R",
  "R/application_storage.R", "R/storage_artifacts.R", "R/model_lifecycle.R",
  "R/application_state.R", "R/application_context.R",
  # Lightweight Stage 6 paths and Stage 8 orchestration.
  "R/stage6_config.R", "R/application_lifecycle.R", "R/stage8_config.R",
  "R/stage8_plan.R"
))

stage8_restart_transaction <- function(config) {
  boundary <- config$restart_from
  assert(!is.null(boundary), "restart mode requires config$restart_from.")
  if (identical(boundary, "stage_2") && file.exists(config$model$manifest) &&
      dir.exists(config$base_paths$applications)) {
    manifests <- list.files(
      config$base_paths$applications,
      pattern = "^scenario_manifest[.]rds$", recursive = TRUE, full.names = TRUE
    )
    manifests <- setdiff(
      normalizePath(manifests, winslash = "/", mustWork = FALSE),
      normalizePath(config$scenario$manifest, winslash = "/", mustWork = FALSE)
    )
    model_md5 <- unname(tools::md5sum(config$model$manifest))
    referenced <- manifests[vapply(manifests, function(path) {
      manifest <- tryCatch(readRDS(path), error = function(e) NULL)
      identical(manifest$contract$persistence_model_manifest_md5, model_md5)
    }, logical(1L))]
    assert(!length(referenced), paste0(
      "This persistence model is referenced by another application and cannot be restarted in place:\n",
      paste(referenced, collapse = "\n"),
      "\nUse a different threshold/horizon or deliberately remove dependent applications first."
    ))
  }
  targets <- switch(
    boundary,
    stage_2 = c(config$model$root, config$scenario$root),
    stage_4 = config$scenario$root,
    stage_6 = vapply(config$priority_curves, function(curve) {
      application_curve_paths(config$paths, curve)$out_curve
    }, character(1L))
  )
  targets <- unique(targets[file.exists(targets) | dir.exists(targets)])
  if (!length(targets)) return(list())
  transactions <- list()
  for (target in targets) {
    transactions[[length(transactions) + 1L]] <- tryCatch(
      begin_directory_replacement(target, paste0("Stage 8 ", boundary, " output")),
      error = function(e) {
        invisible(lapply(rev(transactions), function(x) x$rollback()))
        stop(e)
      }
    )
  }
  transactions
}

# Execute selected model stages. Stage 2 alone updates recoverable model state;
# Stage 3 is the boundary that publishes the complete persistence manifest.
stage8_execute_model_stages <- function(config) {
  selected2 <- stage8_stage_selected(config, "stage_2")
  selected3 <- stage8_stage_selected(config, "stage_3")
  figures <- character()
  completed <- character()

  if (selected2) {
    validate_demography_model(config$base_paths)
    project_source("R/stage2_workflow.R")
    s2 <- stage2_config(
      config$stage2_params, config$model$project, config$model$source_project
    )
    result2 <- run_stage2_model(s2, config$model)
    completed <- c(completed, "stage_2")
    rm(result2, s2); gc(FALSE)
  }

  if (selected3) {
    project_source("R/stage3_workflow.R")
    s3 <- stage3_config(
      list(
        mode = if (file.exists(file.path(
          config$model$stage3, "persistence_curve_models.rds"
        ))) "reuse" else "fit",
        taxa = "both", loess_span = config$loess_span,
        write_figures = config$stage3_figures
      ),
      paths = config$model$project
    )
    result3 <- run_stage3_model(s3, config$model)
    figures <- result3$figures
    completed <- c(completed, "stage_3")
    rm(result3, s3); gc(FALSE)
  }

  needs_model <- any(c(
    "stage_4", "stage_5", "stage_5_1", "stage_5_2", "stage_5_3",
    "stage_6", "stage_7_1"
  ) %in% config$stages)
  manifest <- if (selected3 || (!selected2 && needs_model)) {
    validate_persistence_model(config)
  } else NULL
  list(
    manifest = manifest, figures = figures, completed = completed,
    complete = !is.null(manifest)
  )
}

# Execute Stage 4--5 science and independently selected optional products.
# Core scenario publication requires only Stage 4 and Stage 5 records.
stage8_execute_scenario_stages <- function(config, model_manifest, state) {
  selected <- intersect(
    c(application_scenario_core_stages(),
      application_scenario_optional_stages()),
    config$stages
  )
  if (!length(selected)) {
    return(list(
      state = state, manifest = NULL, application_manifest = NULL,
      stage53 = NULL, figures = character(), completed = character(),
      complete = file.exists(config$scenario$manifest)
    ))
  }

  ensure_writable_dir(config$scenario$root, "application threshold directory")
  app_manifest <- if (stage8_stage_selected(config, "stage_4")) {
    requested <- application_manifest(config, active = TRUE)
    application_write_or_validate_manifest(
      config$scenario$application_manifest, requested, "application manifest",
      "Use a new application name for different species, SDMs, land cover, or study area."
    )
  } else {
    config$existing_application_manifest %||% {
      need_file(config$scenario$application_manifest, "application manifest")
      readRDS(config$scenario$application_manifest)
    }
  }
  completed <- character(); figures <- character(); stage53_summary <- NULL

  if (!stage8_stage_selected(config, "stage_4") &&
      any(c("stage_5", "stage_5_1") %in% selected)) {
    stage8_validate_spatial_handoff(config, app_manifest)
  }

  if (stage8_stage_selected(config, "stage_5") &&
      !stage8_stage_selected(config, "stage_4")) {
    reused_stage4 <- application_reuse_boundary(
      state, "stage_4", config$base_paths$root
    )
    assert(
      !is.null(reused_stage4),
      "Selected Stage 5 requires a valid omitted Stage 4 recovery record."
    )
    state <- reused_stage4
  }

  optional_only <- any(application_scenario_optional_stages() %in% selected) &&
    !any(application_scenario_core_stages() %in% selected)
  if (optional_only) application_validate_scenario_handoff(config)

  if (stage8_stage_selected(config, "stage_4")) {
    frozen <- read_application_inputs(config)
    reused <- application_reuse_boundary(state, "stage_4", config$base_paths$root)
    if (!is.null(reused)) {
      state <- reused
    } else {
      project_source("R/stage4_workflow.R")
      s4 <- stage4_application_config(
        stage8_stage4_params(config, frozen), config,
        frozen_inputs = frozen,
        iucn_mode = if (is.null(frozen)) "cache_or_query" else "cache_only"
      )
      result <- run_stage4(s4)
      if (is.null(frozen)) write_application_inputs(config, result$habitat_rows)
      state <- application_commit_boundary(
        config, application_stage4_boundary(result, config), state
      )
      rm(result, s4); gc(FALSE)
    }
    completed <- c(completed, "stage_4")
  }

  if (stage8_stage_selected(config, "stage_5")) {
    reused <- application_reuse_boundary(state, "stage_5", config$base_paths$root)
    if (!is.null(reused)) {
      state <- reused
    } else {
      project_source("R/stage5_workflow.R")
      s5 <- stage5_application_config(
        stage8_stage5_params(config), config,
        landcover_fingerprint = config$inputs$landcover
      )
      result <- run_stage5(s5)
      state <- application_commit_boundary(
        config, application_stage5_boundary(s5, result), state
      )
      rm(result, s5); gc(FALSE)
    }
    completed <- c(completed, "stage_5")
  }

  if (stage8_stage_selected(config, "stage_5_1")) {
    reused <- application_reuse_boundary(
      state, "stage_5_1", config$base_paths$root,
      expected_summary = list(
        species = config$process_figure_species %||% NA_character_
      )
    )
    if (!is.null(reused)) {
      state <- reused
      figures <- c(figures, vapply(
        state$stages$stage_5_1$outputs, storage_record_path, character(1L),
        root = config$base_paths$root
      ))
    } else {
      project_source("R/stage51_workflow.R")
      s51 <- stage51_application_config(stage8_stage51_params(config), config)
      result <- run_stage51(s51)
      boundary <- application_stage51_boundary(s51, result)
      state <- application_commit_boundary(config, boundary, state)
      figures <- c(figures, boundary$outputs)
      rm(result, s51); gc(FALSE)
    }
    completed <- c(completed, "stage_5_1")
  }

  if (stage8_stage_selected(config, "stage_5_2")) {
    reused <- application_reuse_boundary(state, "stage_5_2", config$base_paths$root)
    if (!is.null(reused)) {
      state <- reused
    } else {
      project_source("R/stage52_workflow.R")
      s52 <- stage52_application_config(stage8_stage52_params(config), config)
      result <- run_stage52(s52)
      state <- application_commit_boundary(
        config, application_stage52_boundary(s52, result), state
      )
      rm(result, s52); gc(FALSE)
    }
    completed <- c(completed, "stage_5_2")
  }

  if (stage8_stage_selected(config, "stage_5_3")) {
    reused <- application_reuse_boundary(
      state, "stage_5_3", config$base_paths$root,
      expected_summary = list(curve = config$report_curve)
    )
    if (!is.null(reused)) {
      state <- reused
      stage53_summary <- state$stages$stage_5_3$summary
      figures <- c(figures, stage53_summary$figure %||% character())
    } else {
      project_source("R/stage53_workflow.R")
      s53 <- stage53_application_config(stage8_stage53_params(config), config)
      result <- run_stage53(s53)
      boundary <- application_stage53_boundary(s53, result)
      stage53_summary <- boundary$summary
      figures <- c(figures, boundary$outputs)
      state <- application_commit_boundary(config, boundary, state)
      rm(result, s53); gc(FALSE)
    }
    completed <- c(completed, "stage_5_3")
  }

  scenario_manifest <- application_finalize_scenario(
    config, state, model_manifest, app_manifest,
    require_complete = FALSE
  )
  list(
    state = state, manifest = scenario_manifest,
    application_manifest = app_manifest, stage53 = stage53_summary,
    figures = figures, completed = completed,
    complete = !is.null(scenario_manifest) || file.exists(config$scenario$manifest)
  )
}

# Execute selected Stage 6 curves and optional Stage 7.1 reports. Scientific
# curve completion is committed before presentation so reports can be omitted
# or regenerated later without changing the run handoff.
stage8_execute_priority_stages <- function(config, state) {
  run_selected <- stage8_stage_selected(config, "stage_6")
  report_selected <- stage8_stage_selected(config, "stage_7_1")
  if (!run_selected && !report_selected) return(list(
    state = state, incomplete = list(), stage71 = list(),
    figures = character(), completed = character(), completed_curves = character()
  ))

  application_validate_scenario_handoff(config)
  # Read or initialize this once. The returned object is threaded through curve
  # commits so neither report-only runs nor multi-curve runs reread it per curve.
  run_manifest <- application_priority_run_manifest(config)
  summaries <- list(); incomplete <- list(); figures <- character()
  completed <- character(); completed_curves <- character()

  # Stage 6 initialization belongs to the threshold scenario, not a curve or
  # removal schedule. Validate or construct it exactly once before entering the
  # curve loop, then publish its recovery record back to the scenario manifest
  # so compatible alternative schedules inherit the same artifact.
  if (run_selected) {
    project_source("R/stage6_workflow.R")
    initialization_key <- application_scenario_shared_stages()
    s6i <- stage6_application_config(
      stage8_stage6_params(config, "initialize"), config, TRUE
    )
    reused <- application_reuse_boundary(
      state, initialization_key, config$base_paths$root
    )
    if (!is.null(reused)) {
      state <- reused
    } else {
      result <- build_stage6_initialization(s6i)
      state <- application_commit_boundary(
        config,
        application_stage6_initialization_boundary(s6i, result),
        state
      )
      application_finalize_scenario(config, state)
      rm(result)
    }
  }

  for (curve in config$priority_curves) {
    if (run_selected) {
      run_key <- application_stage_key("stage_6", curve)
      reused <- application_reuse_boundary(state, run_key, config$base_paths$root)
      if (!is.null(reused)) {
        state <- reused
      } else {
        inspect <- stage6_application_config(
          stage8_stage6_params(config, "inspect", curve), config, FALSE
        )
        status <- inspect_stage6(inspect)
        run_mode <- if (nrow(status$checkpoints)) "resume" else "run"
        s6 <- stage6_application_config(
          stage8_stage6_params(config, run_mode, curve), config, TRUE
        )
        result <- run_stage6(
          s6,
          initialization_record_validated = TRUE
        )
        pipeline <- result$priority_pipeline_result
        if (!isTRUE(pipeline$frontier_exhausted)) {
          incomplete[[curve]] <- list(completed_stages = pipeline$completed_stages)
          rm(result, pipeline, s6, inspect, status); gc(FALSE)
          next
        }
        state <- application_commit_boundary(
          config, application_stage6_run_boundary(s6, result), state
        )
        rm(result, pipeline, s6, inspect, status); gc(FALSE)
      }
      run_manifest <- application_mark_curve_complete(
        config, curve, state, run_manifest
      )
      completed_curves <- c(completed_curves, curve)
    } else {
      assert(
        application_curve_manifest_complete(run_manifest, curve),
        paste0("Stage 7.1 requires completed Stage 6 curve: ", curve, ".")
      )
      completed_curves <- c(completed_curves, curve)
    }

    if (report_selected) {
      report_key <- application_stage_key("stage_7_1", curve)
      reused <- application_reuse_boundary(state, report_key, config$base_paths$root)
      if (!is.null(reused)) {
        state <- reused
        summaries[[curve]] <- state$stages[[report_key]]$summary
        figures <- c(
          figures,
          summaries[[curve]]$outputs %||% character()
        )
      } else {
        project_source("R/stage71_workflow.R")
        s71 <- stage71_application_config(
          stage8_stage71_params(config, curve), config
        )
        result <- run_stage71(s71)
        boundary <- application_stage71_boundary(s71, result)
        summaries[[curve]] <- boundary$summary
        figures <- c(figures, boundary$outputs)
        state <- application_commit_boundary(config, boundary, state)
        rm(result, s71); gc(FALSE)
      }
    }
  }
  if (exists("s6i", inherits = FALSE)) rm(s6i)
  if (run_selected && length(completed_curves)) completed <- c(completed, "stage_6")
  if (report_selected && !length(incomplete)) completed <- c(completed, "stage_7_1")
  list(
    state = state, incomplete = incomplete, stage71 = summaries,
    figures = figures, completed = completed,
    completed_curves = unique(completed_curves)
  )
}

# Inspect or execute the exact selected application stages.
run_stage8 <- function(config) {
  describe_stage8_config(config)
  state <- read_application_state(config)
  plan <- stage8_plan(config, state)
  if (identical(config$mode, "inspect")) {
    return(list(
      mode = "inspect", application = config$application,
      selected_stages = config$stages,
      completed_stages = character(),
      skipped_stages = setdiff(stage8_stage_order(), config$stages),
      model_root = config$model$root, scenario_root = config$scenario$root,
      run_root = config$scenario$run_root,
      model_complete = file.exists(config$model$manifest),
      scenario_complete = file.exists(config$scenario$manifest),
      completed_curves = character(), complete = FALSE, plan = plan,
      inputs = config$inputs, study_area = config$study_area,
      status = application_status_table(state),
      artifacts = application_artifact_table(state, config$base_paths$root),
      manifest = if (file.exists(config$scenario$manifest)) {
        readRDS(config$scenario$manifest)
      } else NULL,
      stage53 = NULL, stage71 = NULL, stage71_by_curve = list(),
      stage6_incomplete = list(), figures = character()
    ))
  }
  stage8_assert_plan_ready(plan)

  transaction <- if (identical(config$mode, "restart")) {
    stage8_restart_transaction(config)
  } else NULL
  success <- FALSE
  on.exit({
    if (length(transaction)) invisible(lapply(rev(transaction), function(x) {
      if (success) x$finalize() else x$rollback()
    }))
  }, add = TRUE)

  if (identical(config$restart_from, "stage_2") ||
      identical(config$restart_from, "stage_4")) {
    state <- new_application_state(config)
  }

  model <- stage8_execute_model_stages(config)
  scenario <- stage8_execute_scenario_stages(config, model$manifest, state)
  state <- scenario$state
  priority <- stage8_execute_priority_stages(config, state)
  state <- priority$state

  complete <- !length(priority$incomplete)
  completed <- unique(c(model$completed, scenario$completed, priority$completed))
  figures <- unique(c(model$figures, scenario$figures, priority$figures))
  figures <- figures[grepl("[.](png|jpe?g|svg|pdf)$", figures, ignore.case = TRUE)]
  success <- TRUE
  if (length(transaction)) {
    invisible(lapply(transaction, function(x) x$finalize()))
    transaction <- list()
  }
  list(
    mode = config$mode, application = config$application,
    selected_stages = config$stages, completed_stages = completed,
    skipped_stages = setdiff(stage8_stage_order(), config$stages),
    model_root = config$model$root, scenario_root = config$scenario$root,
    run_root = config$scenario$run_root,
    model_complete = model$complete,
    scenario_complete = file.exists(config$scenario$manifest),
    completed_curves = priority$completed_curves,
    complete = complete,
    stop_reason = if (complete) NULL else paste0(
      "Stage 6 reached max_stages before frontier exhaustion for: ",
      paste(names(priority$incomplete), collapse = ", "),
      ". Resume after increasing max_stages."
    ),
    plan = plan, inputs = config$inputs, study_area = config$study_area,
    status = application_status_table(state),
    artifacts = application_artifact_table(state, config$base_paths$root),
    manifest = if (file.exists(config$scenario$manifest)) {
      readRDS(config$scenario$manifest)
    } else NULL,
    stage53 = scenario$stage53,
    stage71 = if (is.null(config$report_curve)) NULL else
      priority$stage71[[config$report_curve]] %||% NULL,
    stage71_by_curve = priority$stage71,
    stage6_incomplete = priority$incomplete,
    figures = figures
  )
}
