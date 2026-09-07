# Shared Stage 6 initialization, inspection, and curve-run lifecycle.
#
# This module keeps orchestration outside the optimized scientific stages. A
# run validates its immutable initialization before changing output state, reconstructs
# transient indexes after resume, and commits replacement runs only after all
# scientific outputs validate. Shared initialization and active curve runs use
# separate scenario-level and run/curve runtime logs.

if (!exists("project_source", mode = "function")) {
  source(file.path("R", "project_loader.R"))
}
project_source(c(
  # Lightweight contracts and application adapters used by every mode.
  "R/project_utils.R", "R/project_artifacts.R", "R/project_rmd.R",
  "R/project_runtime_log.R", "R/project_transactions.R",
  "R/project_paths.R", "R/project_config.R", "R/demographic_contract.R",
  "R/analysis_contract.R", "R/priority_run_config.R",
  "R/species_table_config.R", "R/patch_contract.R", "R/patch_config.R",
  "R/application_storage.R", "R/storage_artifacts.R", "R/model_lifecycle.R",
  "R/application_state.R",
  "R/application_context.R", "R/stage6_config.R",
  "R/application_lifecycle.R", "R/application_rmd.R"
))

stage6_module_contract <- function(kind = c("initialization", "runtime")) {
  kind <- match.arg(kind)
  initialization <- c(
    patch_paths = "R/patch_config.R",
    validate_patch_output_consistency = "R/patch_contract.R",
    initialize_priority_inputs = "R/priority_inputs.R"
  )
  if (identical(kind, "initialization")) return(initialization)

  c(
    patch_paths = "R/patch_config.R",
    new_priority_graph_store = "R/priority_graph_store.R",
    load_stage6_runtime_kernels = "R/priority_runtime_kernels.R",
    find_csr_components = "R/priority_csr_graph.R",
    build_rook_neighbor_index = "R/priority_pruning_frontier.R",
    score_frontier_cells = "R/priority_pruning_scoring.R",
    log_ecology_diagnostics = "R/priority_ecology_logging.R",
    rebuild_pu_after_patch_loss = "R/priority_pruning_graph.R",
    run_pruning_iteration = "R/priority_pruning_iteration.R",
    run_pruning_stage = "R/priority_pruning_stage.R",
    update_species_patch_index_incremental = "R/priority_fragmentation_patches.R",
    build_provisional_pu_graph = "R/priority_fragmentation_graph.R",
    run_fragmentation_stage = "R/priority_fragmentation_stage.R",
    build_distance_predicate_lookup = "R/priority_distance_geometry.R",
    filter_distance_invalid_edges_for_pu = "R/priority_distance_graph.R",
    run_distance_connectivity_stage = "R/priority_distance_stage.R",
    write_rankmap_raster = "R/priority_outputs.R",
    write_priority_checkpoint = "R/priority_checkpoints.R",
    log_priority_start = "R/priority_logging.R",
    run_priority_pipeline = "R/priority_pipeline.R"
  )
}

load_stage6_modules <- function(kind = c("initialization", "runtime")) {
  kind <- match.arg(kind)
  target <- environment(load_stage6_modules)
  modules <- stage6_module_contract(kind)
  for (required_function in names(modules)) {
    if (!exists(
      required_function,
      envir = target,
      mode = "function",
      inherits = TRUE
    )) {
      sys.source(modules[[required_function]], envir = target)
    }
  }
  invisible(modules)
}

inspect_stage6 <- function(config) {
  assert(is.list(config), "config must be a validated Stage 6 configuration.")
  assert(identical(config$mode, "inspect"), "inspect_stage6() requires mode = 'inspect'.")

  checkpoints <- if (dir.exists(config$paths$checkpoint_dir)) {
    files <- list.files(
      config$paths$checkpoint_dir,
      pattern = "^priority_checkpoint_stage_[0-9]+[.]rds$",
      full.names = TRUE
    )
    stages <- suppressWarnings(as.integer(sub(
      "^priority_checkpoint_stage_([0-9]+)[.]rds$",
      "\\1",
      basename(files)
    )))
    order_index <- order(stages, files)
    data.frame(path = files[order_index], stage = stages[order_index])
  } else {
    data.frame(path = character(), stage = integer())
  }

  artifact_paths <- c(
    initialization = config$paths$initialization_bundle,
    run = config$paths$out_curve,
    replacement_backup = config$paths$replacement_backup,
    runtime_log = config$paths$runtime_log,
    patch_lookups = config$paths$patch_lookup_output_dir,
    removal_events = config$paths$removal_events,
    removal_order = config$paths$removal_order,
    rankmap = config$paths$rankmap
  )
  artifacts <- project_artifact_status(artifact_paths)
  artifacts$artifact <- names(artifact_paths)
  artifacts <- artifacts[, c("artifact", "path", "exists", "type", "size_bytes", "modified")]

  log_msg(
    "STATUS | initialization=", if (artifacts$exists[artifacts$artifact == "initialization"]) "present" else "missing",
    " run=", if (artifacts$exists[artifacts$artifact == "run"]) "present" else "missing",
    " checkpoints=", nrow(checkpoints)
  )
  invisible(list(artifacts = artifacts, checkpoints = checkpoints))
}

build_stage6_initialization_active <- function(config) {
  assert(is.list(config), "config must be a validated Stage 6 configuration.")
  assert(
    identical(config$mode, "initialize"),
    "build_stage6_initialization() requires mode = 'initialize'."
  )
  load_stage6_modules("initialization")
  result <- initialize_priority_inputs(config)

  expected_path <- normalizePath(
    config$paths$initialization_bundle,
    winslash = "/",
    mustWork = FALSE
  )
  actual_path <- normalizePath(
    result$initialization_bundle_path,
    winslash = "/",
    mustWork = FALSE
  )
  assert(
    identical(actual_path, expected_path) && file.exists(expected_path),
    paste0(
      "Stage 6 initialization was not installed at its configured path:\n",
      expected_path
    )
  )

  invisible(result)
}

build_stage6_initialization <- function(config) {
  assert(is.list(config), "config must be a validated Stage 6 configuration.")
  assert(identical(config$mode, "initialize"),
         "build_stage6_initialization() requires mode = 'initialize'.")
  result <- with_runtime_log(
    path = config$paths$initialization_runtime_log,
    stage = "6",
    operation = "priority_initialization",
    mode = config$mode,
    context = list(
      scenario = config$paths$scenario_root,
      taxa = config$taxa_tag,
      sdm = config$sdm
    ),
    code = function() build_stage6_initialization_active(config)
  )
  result$runtime_log_path <- config$paths$initialization_runtime_log
  invisible(result)
}

stage6_initialization_dependency_paths <- function(retained_species,
                                           species_table_csv_path,
                                           patch_lookup_rds_path,
                                           connectivity_rds_path,
                                           patch_raster_dir) {
  retained_species <- trimws(as.character(retained_species))
  assert(
    length(retained_species) > 0L &&
      all(!is.na(retained_species)) &&
      all(nzchar(retained_species)),
    "retained_species must contain at least one nonblank species name."
  )

  c(
    validate_path_param(species_table_csv_path, "species_table_csv_path"),
    validate_path_param(patch_lookup_rds_path, "patch_lookup_rds_path"),
    validate_path_param(connectivity_rds_path, "connectivity_rds_path"),
    file.path(
      validate_path_param(patch_raster_dir, "patch_raster_dir"),
      vapply(retained_species, patch_filename_from_scientific, character(1L))
    )
  )
}

validate_stage6_initialization_freshness <- function(initialization_path,
                                             retained_species,
                                             species_table_csv_path,
                                             patch_lookup_rds_path,
                                             connectivity_rds_path,
                                             patch_raster_dir) {
  dependencies <- stage6_initialization_dependency_paths(
    retained_species = retained_species,
    species_table_csv_path = species_table_csv_path,
    patch_lookup_rds_path = patch_lookup_rds_path,
    connectivity_rds_path = connectivity_rds_path,
    patch_raster_dir = patch_raster_dir
  )
  validate_cached_artifact_freshness(
    cached_path = initialization_path,
    dependency_paths = dependencies,
    label = "Stage 6 shared initialization",
    rebuild_hint = "Run Stage 6 with mode: 'initialize' before starting the pipeline."
  )
  invisible(dependencies)
}

run_stage6_active <- function(config, initialization_record_validated = FALSE) {
  assert(is.list(config), "config must be a validated Stage 6 configuration.")
  assert(config$mode %in% c("run", "resume"), "run_stage6() requires mode = 'run' or 'resume'.")
  load_stage6_modules("runtime")
  load_stage6_runtime_kernels(cache_dir = config$paths$rcpp_cache)
  run_started <- proc.time()[["elapsed"]]

  flags <- priority_flags_from_tags(config$taxa_tag, config$sdm)
  initialization_bundle_path <- config$paths$initialization_bundle
  pipeline_output_dir <- config$paths$out_curve
  checkpoint_output_dir <- config$paths$checkpoint_dir
  resume_requested <- identical(config$mode, "resume")
  resume_checkpoint_file <- if (resume_requested) {
    resolve_priority_checkpoint_path(
      resume_checkpoint_path = if (is.null(config$resume_stage)) {
        "latest"
      } else {
        priority_checkpoint_path(checkpoint_output_dir, config$resume_stage)
      },
      checkpoint_dir = checkpoint_output_dir
    )
  } else {
    NULL
  }
  if (resume_requested) {
    assert(
      identical(
        normalizePath(dirname(resume_checkpoint_file), winslash = "/", mustWork = TRUE),
        normalizePath(checkpoint_output_dir, winslash = "/", mustWork = TRUE)
      ),
      "The resume checkpoint must belong to the selected Stage 6 run directory."
    )
  }
  if (resume_requested) {
    need_dir(pipeline_output_dir, "Stage 6 run directory to resume")
  }

  assert(
    !file.exists(config$paths$replacement_backup),
    paste0("An unresolved Stage 6 replacement backup exists: ", config$paths$replacement_backup)
  )
  replacement_transaction <- NULL
  run_succeeded <- FALSE
  need_file(initialization_bundle_path, "saved Stage 6 shared initialization")
  need_dir(config$paths$patch_dir, "Stage 5 patch raster directory")

  # ---------------------------------------------------------------------
  # Load the shared curve-neutral initialization
  # ---------------------------------------------------------------------
  initialization_bundle <- readRDS(initialization_bundle_path)

  if (is.null(initialization_bundle$metadata)) {
    stop("Saved Stage 6 initialization does not contain metadata.")
  }


  validate_priority_initialization_schema(initialization_bundle)
  retained_species <- validate_priority_initialization_metadata(
    initialization_bundle,
    taxa = config$taxa_tag,
    sdm = config$sdm,
    contract = config$contract,
    label = "saved Stage 6 initialization"
  )
  initialization_created_at <- validate_scalar_string(
    initialization_bundle$metadata$created_at,
    "saved Stage 6 initialization metadata$created_at"
  )

  if (!initialization_record_validated) {
    initialization_dependencies <- validate_stage6_initialization_freshness(
      initialization_path = initialization_bundle_path,
      retained_species = retained_species,
      species_table_csv_path = config$paths$species_table,
      patch_lookup_rds_path = config$paths$patch_lookup,
      connectivity_rds_path = config$paths$connectivity,
      patch_raster_dir = config$paths$patch_dir
    )
    runtime_log_event(
      "priority_initialization_freshness_validated",
      source_files = length(initialization_dependencies),
      initialization = normalizePath(initialization_bundle_path, mustWork = FALSE)
    )
    rm(initialization_dependencies)
  } else {
    # The Rmd and Stage 8 validate the scenario recovery record—its output and
    # every Stage 4/5 source checksum—once before execution. Do not repeat even
    # the cheaper all-raster timestamp scan for every curve in that operation.
    runtime_log_event(
      "priority_initialization_recovery_record_reused",
      initialization = normalizePath(initialization_bundle_path, mustWork = FALSE)
    )
  }

  # Package loading is deferred until initialization identity and source freshness
  # have passed their inexpensive metadata checks.
  load_packages(c("data.table", "terra", "sf", "Rfast", "fastmatch"))

  # ---------------------------------------------------------------------
  # Restore static objects and either fresh or checkpointed mutable state
  # ---------------------------------------------------------------------

  cell_area_by_cell <- initialization_bundle$cell_area_by_cell
  species_params <- priority_species_parameters_for_curve(
    initialization_bundle$species_params,
    config$curve,
    contract = config$contract
  )
  coefficient_identity <- priority_curve_coefficient_identity(
    initialization_bundle$species_params,
    config$curve,
    contract = config$contract
  )
  rook_neighbor_pairs <- initialization_bundle$rook_neighbor_pairs


  # ---------------------------------------------------------------------
  # Validate required static objects
  # ---------------------------------------------------------------------

  if (!length(cell_area_by_cell)) {
    stop("Loaded cell_area_by_cell is empty.")
  }

  if (!nrow(species_params)) {
    stop("Loaded species_params is empty.")
  }

  if (!nrow(rook_neighbor_pairs)) {
    stop("Loaded rook_neighbor_pairs is empty.")
  }

  required_species_param_columns <- c(
    "species",
    "density",
    "dispersal_distance_km",
    "min_patch_area_km2",
    "min_population_area_km2",
    "a_pred",
    "b_pred",
    "taxon_class",
    "sdm_method",
    "redlist_category"
  )

  missing_species_param_columns <- setdiff(required_species_param_columns, names(species_params))
  if (length(missing_species_param_columns)) {
    stop(
      "Loaded species_params is missing required column(s): ",
      paste(missing_species_param_columns, collapse = ", ")
    )
  }

  validate_priority_species_parameters(
    species_params, "Loaded species_params", contract = config$contract
  )

  checkpoint_metadata <- priority_checkpoint_metadata(
    initialization_schema_version = initialization_bundle$metadata$schema_version,
    initialization_created_at = initialization_created_at,
    curve_label = config$curve,
    do_mammals = flags$do_mammals,
    do_birds = flags$do_birds,
    do_ppm = flags$do_ppm,
    do_rangebag = flags$do_rangebag,
    retained_species = retained_species,
    n_cells = length(cell_area_by_cell),
    cells_to_remove_per_iteration = config$cells_to_remove_per_iteration,
    pruning_iterations_per_stage = config$pruning_iterations_per_stage,
    analysis_contract = config$contract
  )

  resume_checkpoint <- NULL

  if (isTRUE(resume_requested)) {
    initialization_bundle$patch_id_by_species_list <- NULL
    initialization_bundle$patch_cell_index_by_species_list <- NULL
    initialization_bundle$patch_table <- NULL
    initialization_bundle$pu_graphs_by_key <- NULL
    initialization_bundle$alive_species_count_by_cell <- NULL
    gc(FALSE)

    resume_checkpoint <- readRDS(resume_checkpoint_file)
    attr(resume_checkpoint, "checkpoint_path") <- resume_checkpoint_file

    validate_priority_checkpoint_schema(
      checkpoint = resume_checkpoint,
      expected_metadata = checkpoint_metadata,
      n_cells = length(cell_area_by_cell)
    )

    patch_table <- as.data.table(resume_checkpoint$patch_table)
    pu_graphs_by_key <- resume_checkpoint$pu_graphs_by_key
    alive_species_count_by_cell <- resume_checkpoint$alive_species_count_by_cell
    patch_id_by_species_env <- resume_checkpoint$patch_id_by_species_env
    patch_cell_index_by_species_env <- NULL

    gc(FALSE)
  } else {
    patch_table <- as.data.table(initialization_bundle$patch_table)
    pu_graphs_by_key <- initialization_bundle$pu_graphs_by_key
    alive_species_count_by_cell <- initialization_bundle$alive_species_count_by_cell

    if (!nrow(patch_table)) {
      stop("Loaded patch_table is empty.")
    }

    if (!is.list(pu_graphs_by_key)) {
      stop("Loaded pu_graphs_by_key is not a list.")
    }

    if (!length(alive_species_count_by_cell)) {
      stop("Loaded alive_species_count_by_cell is empty.")
    }

    if (is.null(initialization_bundle$patch_id_by_species_list) ||
        !is.list(initialization_bundle$patch_id_by_species_list)) {
      stop("Bundle does not contain a valid patch_id_by_species_list.")
    }

    patch_id_by_species_env <- new.env(parent = emptyenv())

    for (species_name in names(initialization_bundle$patch_id_by_species_list)) {
      assign(
        species_name,
        initialization_bundle$patch_id_by_species_list[[species_name]],
        envir = patch_id_by_species_env
      )
    }

    if (is.null(initialization_bundle$patch_cell_index_by_species_list) ||
        !is.list(initialization_bundle$patch_cell_index_by_species_list)) {
      stop("Bundle does not contain a valid patch_cell_index_by_species_list.")
    }

    patch_cell_index_by_species_env <- new.env(parent = emptyenv())

    for (species_name in names(initialization_bundle$patch_cell_index_by_species_list)) {
      assign(
        species_name,
        initialization_bundle$patch_cell_index_by_species_list[[species_name]],
        envir = patch_cell_index_by_species_env
      )
    }
  }

  rm(initialization_bundle)
  invisible(gc(FALSE))

  required_patch_columns <- c("species", "patch_id", "pu_id", "patch_area_km2")
  missing_patch_columns <- setdiff(required_patch_columns, names(patch_table))
  if (length(missing_patch_columns)) {
    stop(
      "Loaded patch_table is missing required column(s): ",
      paste(missing_patch_columns, collapse = ", ")
    )
  }


  # ---------------------------------------------------------------------
  # Reconstruct the shared raster template
  # ---------------------------------------------------------------------
  template_species <- retained_species[[1]]

  mask_template_raster_path <- file.path(
    config$paths$patch_dir,
    patch_filename_from_scientific(template_species)
  )

  if (!file.exists(mask_template_raster_path)) {
    stop(
      "Could not reconstruct mask template raster. Missing file: ",
      mask_template_raster_path
    )
  }

  mask_template_raster <- terra::rast(mask_template_raster_path)
  assert(terra::nlyr(mask_template_raster) == 1L, "Stage 6 template raster must have one layer.")

  missing_patch_rasters <- retained_species[!file.exists(file.path(
    config$paths$patch_dir,
    vapply(retained_species, patch_filename_from_scientific, character(1L))
  ))]
  assert(
    !length(missing_patch_rasters),
    paste0("Missing patch rasters for retained species: ", paste(missing_patch_rasters, collapse = ", "))
  )


  # ---------------------------------------------------------------------
  # Optional post-load consistency checks
  # ---------------------------------------------------------------------
  species_in_env <- ls(envir = patch_id_by_species_env, all.names = TRUE)

  if (!setequal(retained_species, species_in_env)) {
    stop(
      "Retained species in metadata do not match species stored in patch_id_by_species_list."
    )
  }

  if (!all(unique(patch_table$species) %in% species_params$species)) {
    stop("Some species in patch_table are missing from species_params.")
  }

  if (!all(unique(patch_table$species) %in% retained_species)) {
    stop("Some species in patch_table are missing from retained_species metadata.")
  }

  if (length(alive_species_count_by_cell) != terra::ncell(mask_template_raster)) {
    stop(
      "alive_species_count_by_cell length does not match the number of cells in the mask template raster."
    )
  }

  if (length(cell_area_by_cell) != terra::ncell(mask_template_raster)) {
    stop(
      "cell_area_by_cell length does not match the number of cells in the mask template raster."
    )
  }


  # A fresh run replaces any prior run through a rollback-capable transaction.
  # Resume modifies the existing run in place and therefore never enters this
  # transaction.
  if (identical(config$mode, "run")) {
    replacement_transaction <- begin_directory_replacement(
      pipeline_output_dir,
      label = "Stage 6 run directory"
    )
  }

  on.exit({
    if (!is.null(replacement_transaction) && !run_succeeded) {
      replacement_transaction$rollback()
    }
  }, add = TRUE)

  dir.create(pipeline_output_dir, recursive = TRUE, showWarnings = FALSE)

  runtime_log_event(
    "priority_run_context",
    run_id = config$paths$run_id,
    output_dir = pipeline_output_dir,
    resume = isTRUE(resume_requested),
    checkpoint = if (isTRUE(resume_requested)) {
      normalizePath(resume_checkpoint_file, mustWork = FALSE)
    } else {
      "none"
    }
  )

  runtime_log_event(
    "priority_initialization_loaded",
    initialization = normalizePath(initialization_bundle_path, mustWork = FALSE),
    species = as.integer(length(retained_species)),
    cells = as.integer(length(cell_area_by_cell)),
    resume = isTRUE(resume_requested)
  )


  # ---------------------------------------------------------------------
  # Run the full priority pipeline
  # ---------------------------------------------------------------------
  priority_pipeline_result <- run_priority_pipeline(
    cells_to_remove_per_iteration   = config$cells_to_remove_per_iteration,
    pruning_iterations_per_stage    = config$pruning_iterations_per_stage,
    patch_table                     = patch_table,
    pu_graphs_by_key                = pu_graphs_by_key,
    alive_species_count_by_cell     = alive_species_count_by_cell,
    patch_id_by_species_env         = patch_id_by_species_env,
    patch_cell_index_by_species_env = patch_cell_index_by_species_env,
    cell_area_by_cell               = cell_area_by_cell,
    rook_neighbor_pairs             = rook_neighbor_pairs,
    species_params                  = species_params,
    template_raster                 = mask_template_raster,
    mask_template_raster            = mask_template_raster,
    output_dir                      = pipeline_output_dir,
    max_stages                      = config$max_stages,
    resume_checkpoint               = resume_checkpoint,
    checkpoint_every_stages         = config$checkpoint_every_stages,
    checkpoint_keep                 = config$checkpoint_keep,
    checkpoint_metadata             = checkpoint_metadata,
    ecology_log_every_iterations    = config$ecology_log_every_iterations
  )


  # ---------------------------------------------------------------------
  # Small summary prints
  # ---------------------------------------------------------------------
  if (!is.null(replacement_transaction)) {
    replacement_transaction$finalize()
  }
  run_succeeded <- TRUE

  runtime_log_event(
    "priority_run_done",
    initialization = initialization_bundle_path,
    output_dir = pipeline_output_dir,
    log_file = runtime_log_path(),
    removal_order = normalizePath(
      priority_pipeline_result$removal_order_path,
      mustWork = FALSE
    ),
    rankmap = normalizePath(priority_pipeline_result$rankmap_path, mustWork = FALSE),
    removal_events = normalizePath(
      priority_pipeline_result$removal_events_path,
      mustWork = FALSE
    ),
    patch_lookup_dir = normalizePath(
      priority_pipeline_result$patch_lookup_output_dir,
      mustWork = FALSE
    ),
    curve = config$curve,
    taxa = config$taxa_tag,
    sdm = config$sdm,
    completed_stages = priority_pipeline_result$completed_stages,
    frontier_exhausted = priority_pipeline_result$frontier_exhausted,
    removal_steps = priority_pipeline_result$removal_step_count,
    initial_alive_cells = priority_pipeline_result$initial_alive_cell_count,
    removed_initial_alive_cells =
      priority_pipeline_result$removed_initial_alive_cell_count,
    unremoved_initial_alive_cells =
      priority_pipeline_result$unremoved_initial_alive_cell_count,
    final_remaining_patches = priority_pipeline_result$final_remaining_patches,
    final_remaining_pus = priority_pipeline_result$final_remaining_pus,
    final_alive_cells = priority_pipeline_result$final_alive_cells,
    removal_order_bytes = file.info(priority_pipeline_result$removal_order_path)$size[[1L]],
    rankmap_bytes = file.info(priority_pipeline_result$rankmap_path)$size[[1L]],
    removal_events_bytes = file.info(priority_pipeline_result$removal_events_path)$size[[1L]],
    elapsed_seconds = round(proc.time()[["elapsed"]] - run_started, 2)
  )

  # Public callers need only completion state and the three durable removal
  # surfaces. Large mutable pipeline state is released inside the pipeline.
  public_pipeline_result <- priority_pipeline_result[c(
    "completed_stages", "frontier_exhausted", "removal_order_path",
    "rankmap_path", "removal_events_path"
  )]
  invisible(list(
    initialization_bundle_path = initialization_bundle_path,
    initialization_created_at = initialization_created_at,
    coefficient_identity = coefficient_identity,
    pipeline_output_dir = pipeline_output_dir,
    runtime_log_path = runtime_log_path(),
    priority_pipeline_result = public_pipeline_result
  ))
}

run_stage6 <- function(config, initialization_record_validated = FALSE) {
  assert(is.list(config), "config must be a validated Stage 6 configuration.")
  assert(config$mode %in% c("run", "resume"),
         "run_stage6() requires mode = 'run' or 'resume'.")
  initialization_record_validated <- validate_scalar_logical(
    initialization_record_validated,
    "initialization_record_validated"
  )
  with_runtime_log(
    path = config$paths$runtime_log,
    stage = "6",
    operation = "priority_pipeline",
    mode = config$mode,
    context = list(
      application = config$paths$application %||% "unspecified",
      run_id = config$paths$run_id,
      curve = config$curve,
      taxa = config$taxa_tag,
      sdm = config$sdm,
      cells_per_iteration = as.integer(config$cells_to_remove_per_iteration),
      iterations_per_stage = as.integer(config$pruning_iterations_per_stage),
      ecology_log_every = config$ecology_log_every_iterations %||% "disabled",
      max_stages = as.character(config$max_stages)
    ),
    code = function() run_stage6_active(
      config,
      initialization_record_validated = initialization_record_validated
    )
  )
}
