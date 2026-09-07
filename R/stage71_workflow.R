# Stage 7.1: summarize a completed Stage 6 run and create its result figures.
#
# Stage zero is read from the immutable shared initialization; positive stages
# are streamed from the lookup CSVs. The initial patch table is consumed once
# and released before surface validation and figure construction. All public
# outputs are staged and promoted together so an incomplete run cannot leave a
# mixture of old and new results.

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
  "R/application_context.R", "R/stage6_config.R",
  "R/application_lifecycle.R", "R/application_rmd.R",
  # Stage 7.1 report definitions.
  "R/stage7_config.R", "R/stage7_contracts.R", "R/stage7_persistence.R",
  "R/stage71_report.R", "R/figure_utils.R", "R/stage71_figures.R",
  "R/stage7_removal_order_map.R"
))

write_stage71_outputs <- function(
  stage_meta,
  persistence,
  rankmap,
  config
) {
  assert(requireNamespace("ggplot2", quietly = TRUE), "ggplot2 is required for Stage 7.1 figures.")
  targets <- c(
    stage_meta = config$paths$stage_meta,
    config$paths$si_figures,
    cell_removal_order = config$paths$cell_removal_order_figure
  )
  assert(
    identical(
      names(targets),
      c(
        "stage_meta", "pipeline_persistence", "curve_uncertainty",
        "species_persistence", "cell_removal_order"
      )
    ) &&
      !anyDuplicated(targets),
    "Stage 7.1 output paths are incomplete or duplicated."
  )

  ensure_writable_dir(config$paths$ana, "Stage 7.1 analysis directory")
  staging_root <- tempfile("stage71_outputs_", tmpdir = config$paths$ana)
  dir.create(staging_root)
  on.exit(unlink(staging_root, recursive = TRUE, force = TRUE), add = TRUE)

  staged_meta <- file.path(staging_root, "stage_meta.csv")
  data.table::fwrite(stage_meta, staged_meta)
  reread_meta <- data.table::fread(staged_meta)
  assert(
    identical(names(reread_meta), names(stage_meta)) &&
      nrow(reread_meta) == nrow(stage_meta),
    "Staged Stage 7.1 metadata failed readback validation."
  )

  save_staged_figure <- function(plot, figure_id) {
    path <- suppressMessages(save_manuscript_figure(
      plot,
      figure_id,
      figure_dir = staging_root
    ))
    need_file(path, paste0("staged Stage 7.1 figure ", figure_id))
    assert(
      file.info(path)$size[[1L]] > 0,
      paste0("Staged Stage 7.1 figure is empty: ", path)
    )
    path
  }

  removal_plot <- build_stage71_removal_order_map(rankmap)
  staged_map <- save_staged_figure(removal_plot, "cell_removal_order")
  rm(removal_plot, rankmap)
  invisible(gc(FALSE))

  pipeline_plot <- build_stage71_persistence_summary_plot(persistence$summary)
  staged_pipeline <- save_staged_figure(
    pipeline_plot,
    "stage7_pipeline_persistence"
  )
  rm(pipeline_plot)
  invisible(gc(FALSE))

  uncertainty_plot <- build_stage71_mean_curve_uncertainty_plot(
    persistence$curve_summary
  )
  staged_uncertainty <- save_staged_figure(
    uncertainty_plot,
    "stage7_curve_uncertainty"
  )
  rm(uncertainty_plot)
  invisible(gc(FALSE))

  species_plot <- build_stage71_persistence_heatmap(persistence$species)
  staged_species <- save_staged_figure(
    species_plot,
    "stage7_species_persistence"
  )
  rm(species_plot)
  invisible(gc(FALSE))

  staged <- c(
    stage_meta = staged_meta,
    pipeline_persistence = staged_pipeline,
    curve_uncertainty = staged_uncertainty,
    species_persistence = staged_species,
    cell_removal_order = staged_map
  )
  project_file_set_transaction(
    staged_paths = unname(staged),
    target_paths = unname(targets),
    overwrite = TRUE,
    label = "Stage 7.1 result set"
  )
  log_msg(
    "STAGE 7.1 OUTPUTS | files=", length(targets),
    " bytes=", sum(file.info(targets)$size),
    " stage_meta=", targets[["stage_meta"]],
    " figure_dir=", dirname(config$paths$cell_removal_order_figure)
  )
  targets
}

run_stage71 <- function(config) {
  start <- proc.time()[["elapsed"]]
  validate_stage7_active_inputs(config)
  check_stage71_removal_map_packages()
  log_msg(
    "STAGE 7.1 START | run_id=", config$run_id,
    " initialization=", config$paths$initialization_bundle,
    " run_dir=", config$paths$out_curve
  )
  suppressPackageStartupMessages({
    library(data.table)
    library(terra)
  })
  lookup_files <- stage7_lookup_files(config$paths$pipe_lookup_dir)
  bundle <- readRDS(config$paths$initialization_bundle)
  persistence_params <- prepare_stage71_bundle_parameters(
    bundle,
    config
  )
  curve_params <- prepare_stage71_curve_parameters(
    persistence_params,
    bundle$species_params,
    configured_curve = config$curve
  )
  curves <- unique(curve_params$curve)
  baseline_state <- new.env(parent = emptyenv())
  baseline_state$patch_table <- bundle$patch_table
  bundle_for_index <- bundle
  alive0 <- as.integer(bundle$alive_species_count_by_cell) > 0L
  cell_area <- as.numeric(bundle$cell_area_by_cell)
  log_msg(
    "STAGE 7.1 INPUTS | lookup_stages=", nrow(lookup_files),
    " retained_species=", nrow(persistence_params),
    " primary_curve=", config$curve,
    " coefficients=shared Stage 6 initialization",
    " baseline_patches=", nrow(bundle$patch_table),
    " raster_cells=", length(alive0)
  )
  rm(bundle)
  template <- terra::rast(config$paths$removal_order)
  assert(length(alive0) == terra::ncell(template), "Bundle alive-domain length differs from the Stage 6 raster.")
  assert(length(cell_area) == terra::ncell(template), "Bundle cell-area length differs from the Stage 6 raster.")
  assert(all(is.finite(cell_area[alive0]) & cell_area[alive0] > 0), "Initially alive cells require positive finite areas.")
  alive0_n <- sum(alive0)
  alive0_area <- sum(cell_area[alive0])

  events <- data.table::fread(config$paths$removal_events)
  state_index <- read_stage6_state_index(list(
    run_id = config$run_id,
    removal_events = config$paths$removal_events,
    patch_lookup_dir = config$paths$pipe_lookup_dir
  ), config$curve, config, bundle = bundle_for_index,
  removal_events = events)
  rm(bundle_for_index)
  stage_meta_columns <- c(
    "stage", "patch_lookup_path", "last_removal_step", "alive_end", "removed_end",
    "area_retained_km2", "area_removed_km2", "cells_removed_stage",
    "area_removed_stage_km2", "events_in_stage", "prop_cells_removed_end",
    "prop_area_removed_end", "pct_cells_removed_end", "pct_area_removed_end"
  )
  stage_meta <- state_index[, stage_meta_columns, with = FALSE]
  rm(state_index)
  log_msg(
    "STAGE 7.1 METADATA | ecological_states=", nrow(stage_meta),
    " target=", config$paths$stage_meta
  )
  rm(cell_area, template)
  invisible(gc(FALSE))

  load_pipeline_lookup <- function(stage, path) {
    if (as.integer(stage) == 0L) {
      assert(
        exists("patch_table", envir = baseline_state, inherits = FALSE),
        "Stage 7.1 baseline patch table was requested more than once."
      )
      out <- get("patch_table", envir = baseline_state, inherits = FALSE)
      rm(list = "patch_table", envir = baseline_state)
      return(out)
    }
    need_file(path, paste0("Stage 6 patch lookup for stage ", stage))
    data.table::fread(
      path,
      select = c("species", "patch_id", "pu_id", "patch_area_km2")
    )
  }
  log_persistence_stage <- function(x) {
    log_msg(
      "persistence_stage | stage=", x$stage,
      " lookup_rows=", x$lookup_rows,
      " pus=", x$pu_count,
      " represented_species=", x$represented_species,
      " mean=", sprintf("%.6f", x$mean_persist),
      " median=", sprintf("%.6f", x$median_persist),
      " elapsed_seconds=", sprintf("%.2f", x$elapsed_seconds)
    )
  }
  log_msg(
    "STAGE 7.1 PERSISTENCE START | ecological_states=", nrow(stage_meta),
    " retained_species=", nrow(persistence_params),
    " primary_curve=", config$curve,
    " uncertainty_curves=", paste(curves, collapse = ","),
    " uncertainty_coefficients=species_table.csv",
    " stage0_source=shared_initialization positive_stage_source=stage_lookup_csv"
  )
  persistence <- compute_stage71_pipeline_persistence(
    stage_meta = stage_meta,
    params = persistence_params,
    lookup_loader = load_pipeline_lookup,
    progress = log_persistence_stage,
    curve_params = curve_params
  )
  assert(
    !exists("patch_table", envir = baseline_state, inherits = FALSE),
    "Stage 7.1 did not release the baseline patch table after stage-zero persistence."
  )
  log_msg(
    "STAGE 7.1 PERSISTENCE COMPLETE | species_stage_rows=", nrow(persistence$species),
    " summary_rows=", nrow(persistence$summary),
    " auc_species=", nrow(persistence$changes),
    " elapsed_seconds=", sprintf("%.1f", proc.time()[["elapsed"]] - start)
  )

  surface <- validate_stage71_removal_surface(
    rankmap = terra::rast(config$paths$rankmap),
    removal_order = terra::rast(config$paths$removal_order),
    events = events,
    alive0 = alive0
  )
  log_msg(
    "STAGE 7.1 SURFACE | ranked_cells=", surface$summary$n_ranked,
    " range=[", sprintf("%.6f", surface$summary$min_rank),
    ",", sprintf("%.6f", surface$summary$max_rank), "]",
    " terminal_cells=", surface$summary$n_terminal,
    " path=", config$paths$rankmap
  )
  rm(alive0, events)
  invisible(gc(FALSE))

  outputs <- write_stage71_outputs(
    stage_meta = stage_meta,
    persistence = persistence,
    rankmap = surface$rankmap,
    config = config
  )
  surface$rankmap <- NULL
  invisible(gc(FALSE))

  main_scope <- if (config$taxa == "both") "all" else config$taxa
  main_summary <- persistence$summary[scope == main_scope]
  data.table::setorder(main_summary, stage_order)
  initial_summary <- main_summary[1L]
  latest_summary <- main_summary[.N]
  lowest_auc <- persistence$changes[which.min(normalized_persistence_auc)]
  elapsed <- proc.time()[["elapsed"]] - start
  log_msg("STAGE 7.1 COMPLETE | stages=", nrow(stage_meta),
          " alive_cells=", alive0_n,
          " alive_area_km2=", sprintf("%.3f", alive0_area),
          " ranked_cells=", surface$summary$n_ranked,
          " priority_range=[", sprintf("%.6f", surface$summary$min_rank),
          ",", sprintf("%.6f", surface$summary$max_rank), "]",
          " retained_species=", nrow(persistence_params),
          " primary_curve=", config$curve,
          " uncertainty_curves=", paste(curves, collapse = ","),
          " persistence_rows=", nrow(persistence$species),
          " initial_mean=", sprintf("%.6f", initial_summary$mean_persist),
          " latest_mean=", sprintf("%.6f", latest_summary$mean_persist),
          " initial_median=", sprintf("%.6f", initial_summary$median_persist),
          " latest_median=", sprintf("%.6f", latest_summary$median_persist),
          " latest_species_with_pu=", latest_summary$n_species_with_pu,
          " latest_frac_gt_05=", sprintf("%.6f", latest_summary$frac_gt_05),
          " lowest_auc_species=", shQuote(lowest_auc$scientificName),
          " lowest_normalized_auc=", sprintf("%.6f", lowest_auc$normalized_persistence_auc),
          " persistence_range=[", sprintf("%.6f", min(persistence$species$sp_persist)),
          ",", sprintf("%.6f", max(persistence$species$sp_persist)), "]",
          " stage_meta=", outputs[["stage_meta"]],
          " outputs=", length(outputs),
          " elapsed_s=", sprintf("%.1f", elapsed))
  list(
    config = config,
    stage_meta = stage_meta,
    surface = surface$summary,
    persistence = persistence,
    outputs = outputs,
    elapsed = elapsed
  )
}
