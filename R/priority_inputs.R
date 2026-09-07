# Shared Stage 6 prioritization-initialization construction.
#
# The builder combines the selected Stage 4 species parameters with finalized
# Stage 5 patch rasters, lookup rows, and CSR connectivity. It reads species
# rasters sequentially, retains flat patch-sorted cell indexes, and writes one
# validated schema-5 bundle through an atomic same-directory transaction.
# Species, patch, PU, graph, adjacency, and raster-cell order are canonical;
# this file deliberately avoids raster stacks and species-by-cell matrices.


validate_priority_patch_values <- function(values, n_cells, species_name) {
  species_name <- validate_scalar_string(species_name, "species_name")
  n_cells <- validate_priority_count(n_cells, "n_cells")
  values <- suppressWarnings(as.numeric(values))

  assert(
    length(values) == n_cells,
    paste0(
      "Patch raster cell count differs from the common Stage 5 grid for species: ",
      species_name
    )
  )

  present <- !is.na(values)
  invalid <- is.nan(values) |
    (present & (!is.finite(values) | values <= 0 | values != floor(values)))
  if (any(invalid)) {
    bad <- which(invalid)
    preview <- utils::head(bad, 10L)
    details <- paste0(
      "cell=", preview,
      " value=", format(values[preview], trim = TRUE, scientific = FALSE)
    )
    stop(
      "Patch raster contains ", length(bad),
      " invalid value(s) for species ", species_name,
      ". Expected positive finite integer patch IDs or NA. First invalid cells: ",
      paste(details, collapse = "; "),
      call. = FALSE
    )
  }

  as.integer(values)
}

validate_priority_patch_id_match <- function(
  raster_patch_ids,
  lookup_patch_ids,
  species_name
) {
  species_name <- validate_scalar_string(species_name, "species_name")
  raster_patch_ids <- sort(unique(as.integer(raster_patch_ids)))
  lookup_patch_ids <- sort(unique(as.integer(lookup_patch_ids)))

  if (!identical(raster_patch_ids, lookup_patch_ids)) {
    # These set differences are intentionally calculated only on failure.
    raster_only <- setdiff(raster_patch_ids, lookup_patch_ids)
    lookup_only <- setdiff(lookup_patch_ids, raster_patch_ids)
    preview_ids <- function(x) {
      if (!length(x)) "none" else paste(utils::head(x, 10L), collapse = ", ")
    }
    stop(
      "Patch ID mismatch for species '", species_name,
      "': raster-only IDs=", preview_ids(raster_only),
      "; lookup-only IDs=", preview_ids(lookup_only), ".",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

convert_priority_connectivity_graphs <- function(connectivity, retained_species) {
  assert(is.list(connectivity), "Stage 5 connectivity must be a list.")
  retained_species <- trimws(as.character(retained_species))
  assert(
    length(retained_species) > 0L &&
      all(!is.na(retained_species)) &&
      all(nzchar(retained_species)),
    "retained_species must contain at least one nonblank species name."
  )

  keep <- vapply(
    connectivity,
    function(entry) {
      is.list(entry) &&
        !is.null(entry$species) &&
        as.character(entry$species) %in% retained_species
    },
    logical(1L)
  )
  selected <- which(keep)
  graphs <- vector("list", length(selected))
  graph_keys <- character(length(selected))

  for (i in seq_along(selected)) {
    entry <- connectivity[[selected[[i]]]]
    species <- validate_scalar_string(as.character(entry$species), "connectivity species")
    pu_id <- validate_priority_count(entry$pu_id, "connectivity pu_id")
    patch_ids <- as.integer(entry$patch_ids)
    row_ptr <- as.integer(entry$row_ptr)
    stored_col_idx <- as.integer(entry$col_idx)
    n_nodes <- length(patch_ids)

    assert(n_nodes > 0L, "Stage 5 connectivity entries must contain patch_ids.")
    assert(
      length(row_ptr) == n_nodes + 1L &&
        row_ptr[[1L]] == 0L &&
        all(diff(row_ptr) >= 0L) &&
        row_ptr[[length(row_ptr)]] == length(stored_col_idx),
      paste0("Invalid stored CSR offsets for ", species, " PU ", pu_id, ".")
    )
    assert(
      !length(stored_col_idx) ||
        all(stored_col_idx >= 0L & stored_col_idx < n_nodes),
      paste0("Stored CSR neighbor index is out of bounds for ", species, " PU ", pu_id, ".")
    )

    graph_keys[[i]] <- paste0(species, "|", pu_id)
    graphs[[i]] <- list(
      species = species,
      pu_id = pu_id,
      id2patch = patch_ids,
      row_ptr = row_ptr,
      col_idx = stored_col_idx + 1L
    )
  }

  assert(!anyDuplicated(graph_keys), "Stage 5 connectivity contains duplicate species/PU keys.")
  stats::setNames(graphs, graph_keys)
}


normalize_taxon_class <- function(x) {
  x0 <- tolower(trimws(as.character(x)))
  x0 <- gsub("[^a-z]+", "", x0)

  out <- rep(NA_character_, length(x0))
  out[x0 %in% c("mammalia", "mammal", "mammals")] <- "mammalia"
  out[x0 %in% c("aves", "bird", "birds")]        <- "aves"

  out
}

normalize_sdm_method <- function(x) {
  x0 <- tolower(trimws(as.character(x)))
  x0 <- gsub("[^a-z0-9]+", "", x0)

  out <- rep(NA_character_, length(x0))
  out[x0 %in% c("ppm")] <- "ppm"
  out[x0 %in% c("rangebag", "rangebags", "range")] <- "rangebag"

  out
}

# Standardize Red List category labels from species_table.csv.
normalize_redlist_category <- function(x) {
  x0 <- tolower(trimws(as.character(x)))
  x0 <- gsub("[^a-z0-9]+", "", x0)

  out <- rep(NA_character_, length(x0))

  out[x0 %in% c("lc", "leastconcern")] <- "least_concern"
  out[x0 %in% c("nt", "nearthreatened")] <- "near_threatened"
  out[x0 %in% c("vu", "vulnerable")] <- "vulnerable"
  out[x0 %in% c("en", "endangered")] <- "endangered"
  out[x0 %in% c("cr", "criticallyendangered")] <- "critically_endangered"
  out[x0 %in% c("dd", "datadeficient")] <- "data_deficient"

  out
}

format_species_preview <- function(x, n = 20L) {
  x <- sort(unique(trimws(as.character(x))))
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) return("")

  shown <- utils::head(x, n)
  suffix <- if (length(x) > n) {
    paste0("\n... and ", length(x) - n, " more")
  } else {
    ""
  }

  paste0(paste(shown, collapse = "\n"), suffix)
}

initialize_priority_inputs <- function(config) {
  assert(is.list(config), "config must be a validated Stage 6 configuration.")
  assert(
    identical(config$mode, "initialize"),
    "initialize_priority_inputs() requires mode = 'initialize'."
  )

  taxa_tag <- validate_priority_taxa_tag(config$taxa_tag, "config$taxa_tag")
  sdm <- validate_priority_sdm_tag(config$sdm, "config$sdm")
  flags <- priority_flags_from_tags(taxa_tag, sdm)
  contract <- validate_analysis_contract(config$contract, "Stage 6 abundance contract")

  species_table_csv_path <- validate_path_param(
    config$paths$species_table,
    "config$paths$species_table"
  )
  patch_raster_dir <- validate_path_param(config$paths$patch_dir, "config$paths$patch_dir")
  patch_lookup_rds_path <- validate_path_param(
    config$paths$patch_lookup,
    "config$paths$patch_lookup"
  )
  connectivity_rds_path <- validate_path_param(
    config$paths$connectivity,
    "config$paths$connectivity"
  )
  # The validated scenario configuration owns the one shared output path. It
  # deliberately contains neither a curve nor a removal-schedule component.
  initialization_bundle_path <- validate_path_param(
    config$paths$initialization_bundle,
    "config$paths$initialization_bundle"
  )
  initialization_bundle_dir <- dirname(initialization_bundle_path)
  # ---------------------------------------------------------------------
  # Existence checks for upstream outputs
  # ---------------------------------------------------------------------
  if (!file.exists(species_table_csv_path)) {
    stop("Missing species table output from Stage 4: ", species_table_csv_path)
  }

  if (!dir.exists(patch_raster_dir)) {
    stop("Missing patch raster directory from Stage 5: ", patch_raster_dir)
  }

  if (!file.exists(patch_lookup_rds_path)) {
    stop("Missing patch lookup output from Stage 5: ", patch_lookup_rds_path)
  }

  if (!file.exists(connectivity_rds_path)) {
    stop("Missing connectivity output from Stage 5: ", connectivity_rds_path)
  }


  # ---------------------------------------------------------------------
  # Package checks
  # ---------------------------------------------------------------------
  load_packages(c("data.table", "terra"))


  # ---------------------------------------------------------------------
  # Read species table from Stage 4 and build species_params
  # ---------------------------------------------------------------------
  species_table_raw <- data.table::fread(species_table_csv_path)

  required_species_columns <- c(
    "scientificName",
    "className",
    "redlistCategory",
    "sdm_method",
    "density",
    "dispersal_dist",
    "min_patch_size",
    "min_pop_size",
    persistence_coefficient_columns()
  )

  missing_species_columns <- setdiff(required_species_columns, names(species_table_raw))
  if (length(missing_species_columns)) {
    stop(
      "species_table.csv is missing required column(s): ",
      paste(missing_species_columns, collapse = ", ")
    )
  }

  # Standardize taxonomic class and SDM-method labels before filtering.
  species_table_raw[
    ,
    taxon_class_std := normalize_taxon_class(className)
  ]

  species_table_raw[
    ,
    sdm_method_std := normalize_sdm_method(sdm_method)
  ]

  species_table_raw[
    ,
    redlist_category_std := normalize_redlist_category(redlistCategory)
  ]

  unknown_taxon_values <- sort(unique(
    as.character(species_table_raw$className[
      is.na(species_table_raw$taxon_class_std) &
        !is.na(species_table_raw$className) &
        nzchar(trimws(as.character(species_table_raw$className)))
    ])
  ))

  if (length(unknown_taxon_values)) {
    stop(
      "Unrecognized className value(s) in species_table.csv: ",
      paste(unknown_taxon_values, collapse = ", ")
    )
  }

  unknown_sdm_values <- sort(unique(
    as.character(species_table_raw$sdm_method[
      is.na(species_table_raw$sdm_method_std) &
        !is.na(species_table_raw$sdm_method) &
        nzchar(trimws(as.character(species_table_raw$sdm_method)))
    ])
  ))

  unknown_redlist_values <- sort(unique(
    as.character(species_table_raw$redlistCategory[
      is.na(species_table_raw$redlist_category_std) &
        !is.na(species_table_raw$redlistCategory) &
        nzchar(trimws(as.character(species_table_raw$redlistCategory)))
    ])
  ))

  if (length(unknown_redlist_values)) {
    stop(
      "Unrecognized redlistCategory value(s) in species_table.csv: ",
      paste(unknown_redlist_values, collapse = ", ")
    )
  }

  if (length(unknown_sdm_values)) {
    stop(
      "Unrecognized sdm_method value(s) in species_table.csv: ",
      paste(unknown_sdm_values, collapse = ", "),
      ". Expected values equivalent to PPM or RangeBag."
    )
  }

  selected_taxon_rows <- (
    (flags$do_mammals & species_table_raw$taxon_class_std == "mammalia") |
    (flags$do_birds   & species_table_raw$taxon_class_std == "aves")
  )

  selected_sdm_rows <- (
    (flags$do_ppm      & species_table_raw$sdm_method_std == "ppm") |
    (flags$do_rangebag & species_table_raw$sdm_method_std == "rangebag")
  )

  selected_species_rows <- selected_taxon_rows & selected_sdm_rows

  species_table_selected <- species_table_raw[selected_species_rows]

  if (!nrow(species_table_selected)) {
    stop(
      "After applying taxon and SDM-source gates, no species remain. ",
      "Current selectors: taxa=", taxa_tag, ", sdm=", sdm, "."
    )
  }

  duplicate_selected_species <- sort(unique(as.character(
    species_table_selected$scientificName[duplicated(species_table_selected$scientificName)]
  )))
  if (length(duplicate_selected_species)) {
    stop(
      "Selected species_table.csv rows contain duplicate scientificName values after applying taxon and SDM gates. ",
      "For mixed PPM+RangeBag runs, rebuild Stage 4 with params$sdm = \"ppm_rangebag\" ",
      "so ambiguous PPM/RangeBag raster matches are diagnosed before Stage 6. ",
      "First duplicate species:\n",
      format_species_preview(duplicate_selected_species),
      call. = FALSE
    )
  }

  selection_summary <- species_table_selected[
    ,
    .N,
    by = .(taxon_class_std, sdm_method_std, redlist_category_std)
  ][order(taxon_class_std, sdm_method_std, redlist_category_std)]

  message(
    "Species selected before spatial-output filtering: ",
    nrow(species_table_selected),
    " species across ",
    nrow(selection_summary),
    " taxon/SDM/Red List metadata groups."
  )

  species_table_selected <- add_canonical_area_thresholds(
    species_table_selected,
    validate = TRUE,
    label = "species_params area thresholds",
    contract = contract
  )

  species_params <- species_table_selected[
    ,
    .(
      species               = as.character(scientificName),
      taxon_class           = as.character(taxon_class_std),
      redlist_category      = as.character(redlist_category_std),
      sdm_method            = as.character(sdm_method_std),
      density               = as.numeric(density),
      dispersal_distance_km = as.numeric(dispersal_dist),
      min_patch_area_km2      = as.numeric(min_patch_area_km2),
      min_population_area_km2 = as.numeric(min_population_area_km2)
    )
  ]
  coefficient_columns <- persistence_coefficient_columns()
  for (column in coefficient_columns) {
    species_params[[column]] <- suppressWarnings(as.numeric(
      species_table_selected[[column]]
    ))
  }

  if (any(!is.finite(species_params$dispersal_distance_km) |
          species_params$dispersal_distance_km <= 0)) {
    stop("species_params contains non-positive or non-finite dispersal_distance_km values.")
  }

  if (any(!is.finite(species_params$density) | species_params$density <= 0)) {
    stop("species_params contains non-positive or non-finite density values.")
  }

  if (any(!is.finite(species_params$min_patch_area_km2) |
          species_params$min_patch_area_km2 <= 0)) {
    stop("species_params contains non-positive or non-finite min_patch_area_km2 values.")
  }

  if (any(!is.finite(species_params$min_population_area_km2) |
          species_params$min_population_area_km2 <= 0)) {
    stop("species_params contains non-positive or non-finite min_population_area_km2 values.")
  }

  invalid_coefficients <- vapply(coefficient_columns, function(column) {
    values <- species_params[[column]]
    any(!is.finite(values) | values <= 0)
  }, logical(1L))
  if (any(invalid_coefficients)) {
    stop(
      "Selected species have missing or non-positive Gompertz parameters in: ",
      paste(coefficient_columns[invalid_coefficients], collapse = ", "),
      ". Regenerate Stages 2–4 before initializing Stage 6."
    )
  }

  setkey(species_params, species)

  rm(
    species_table_raw,
    species_table_selected,
    selection_summary,
    selected_species_rows,
    selected_taxon_rows,
    selected_sdm_rows
  )
  gc(FALSE)


  # ---------------------------------------------------------------------
  # Read patch lookup from Stage 5 and standardize it to patch_table
  # ---------------------------------------------------------------------
  patch_lookup_raw <- readRDS(patch_lookup_rds_path)
  all_connectivity_raw <- readRDS(connectivity_rds_path)

  validate_patch_output_consistency(
    patch_lookup_raw,
    all_connectivity_raw,
    label = "Stage 5 outputs used for priority-input initialization"
  )

  patch_lookup_dt <- as.data.table(patch_lookup_raw)
  patch_table <- patch_lookup_dt[
    ,
    .(
      species        = as.character(scientificName),
      patch_id       = as.integer(patch_id),
      pu_id          = as.integer(pu_id),
      patch_area_km2 = as.numeric(patch_area_km2)
    )
  ]

  patch_table <- patch_table[species %in% species_params$species]

  if (!nrow(patch_table)) {
    stop("patch_table is empty after filtering to the selected species.")
  }

  setkey(patch_table, species, patch_id)

  rm(patch_lookup_raw, patch_lookup_dt)
  gc(FALSE)


  # ---------------------------------------------------------------------
  # Restrict to species that actually have final spatial outputs
  # ---------------------------------------------------------------------
  candidate_species <- species_params$species

  patch_raster_file_by_species <- setNames(
    file.path(
      patch_raster_dir,
      vapply(candidate_species, patch_filename_from_scientific, character(1L))
    ),
    candidate_species
  )

  species_with_patch_raster <- names(patch_raster_file_by_species)[file.exists(patch_raster_file_by_species)]
  species_with_patch_lookup <- unique(patch_table$species)

  missing_patch_raster_species <- setdiff(candidate_species, species_with_patch_raster)
  missing_patch_lookup_species <- setdiff(candidate_species, species_with_patch_lookup)

  if (length(missing_patch_raster_species)) {
    message(
      "Selected species without final Stage 5 patch rasters are excluded from the priority bundle. ",
      "If this is unexpected, rerun Stage 5 with the current species table and SDM settings. First affected species:\n",
      format_species_preview(missing_patch_raster_species)
    )
  }

  if (length(missing_patch_lookup_species)) {
    message(
      "Selected species without Stage 5 patch lookup rows are excluded from the priority bundle. ",
      "If this is unexpected, rerun Stage 5 with the current species table and SDM settings. First affected species:\n",
      format_species_preview(missing_patch_lookup_species)
    )
  }

  retained_spatial_species <- intersect(species_with_patch_raster, species_with_patch_lookup)

  if (!length(retained_spatial_species)) {
    stop("No species remain after restricting to those with both final patch rasters and patch lookup rows.")
  }

  species_params <- species_params[species %in% retained_spatial_species]
  patch_table    <- patch_table[species %in% retained_spatial_species]

  patch_raster_file_by_species <- patch_raster_file_by_species[species_params$species]

  final_selection_summary <- species_params[
    ,
    .N,
    by = .(taxon_class, sdm_method, redlist_category)
  ][order(taxon_class, sdm_method, redlist_category)]

  message(
    "Species retained after requiring patch rasters and patch lookup rows: ",
    nrow(species_params),
    " species across ",
    nrow(final_selection_summary),
    " taxon/SDM/Red List metadata groups."
  )

  rm(final_selection_summary)


  # ---------------------------------------------------------------------
  # Read each Stage 5 patch raster once and build the cell indexes
  # ---------------------------------------------------------------------
  initialization_started <- proc.time()[["elapsed"]]
  species_names <- names(patch_raster_file_by_species)
  template_raster <- terra::rast(patch_raster_file_by_species[[1L]])
  assert(terra::nlyr(template_raster) == 1L, "Stage 5 patch rasters must be single-layer.")

  n_cells <- terra::ncell(template_raster)
  alive_species_count_by_cell <- integer(n_cells)
  patch_id_by_species_env <- new.env(parent = emptyenv())
  patch_cell_index_by_species_env <- new.env(parent = emptyenv())

  for (i in seq_along(species_names)) {
    species_started <- proc.time()[["elapsed"]]
    species_name <- species_names[[i]]
    species_raster <- if (i == 1L) {
      template_raster
    } else {
      terra::rast(patch_raster_file_by_species[[species_name]])
    }

    assert(
      terra::nlyr(species_raster) == 1L,
      paste0("Patch raster must be single-layer for species: ", species_name)
    )
    assert(
      isTRUE(terra::compareGeom(template_raster, species_raster, stopOnError = FALSE)),
      paste0("Patch raster geometry differs from the common Stage 5 grid for species: ", species_name)
    )

    patch_ids <- validate_priority_patch_values(
      terra::values(species_raster, mat = FALSE),
      n_cells = n_cells,
      species_name = species_name
    )
    valid_cells <- !is.na(patch_ids)

    raster_patch_ids <- sort(unique(patch_ids[valid_cells]))
    table_patch_ids <- sort(unique(as.integer(patch_table[.(species_name), patch_id])))
    validate_priority_patch_id_match(
      raster_patch_ids = raster_patch_ids,
      lookup_patch_ids = table_patch_ids,
      species_name = species_name
    )

    assign(species_name, patch_ids, envir = patch_id_by_species_env)
    alive_species_count_by_cell[valid_cells] <-
      alive_species_count_by_cell[valid_cells] + 1L

    occupied_cells <- which(valid_cells)
    occupied_patch_ids <- patch_ids[occupied_cells]
    ordering <- order(occupied_patch_ids)
    assign(
      species_name,
      list(
        pid = as.integer(occupied_patch_ids[ordering]),
        cell = as.integer(occupied_cells[ordering])
      ),
      envir = patch_cell_index_by_species_env
    )

    runtime_log_event(
      "priority_initialization_raster_loaded",
      species_index = paste0(i, "/", length(species_names)),
      species = species_name,
      occupied_cells = sum(valid_cells),
      patches = length(raster_patch_ids),
      elapsed_seconds = round(proc.time()[["elapsed"]] - species_started, 2)
    )

    rm(species_raster, patch_ids, valid_cells, occupied_cells, occupied_patch_ids, ordering)
  }

  # Global grid quantities are identical for every retained species and therefore
  # need to be calculated only from the common template.
  cell_area_by_cell <- as.numeric(terra::values(
    terra::cellSize(template_raster, unit = "km"),
    mat = FALSE
  ))
  rook_neighbor_pairs <- terra::adjacent(
    template_raster,
    cells = seq_len(n_cells),
    directions = 4,
    pairs = TRUE
  )


  # ---------------------------------------------------------------------
  # Convert Stage 5 connectivity graphs to priority-bundle CSR graphs
  # ---------------------------------------------------------------------
  csr_version   <- attr(all_connectivity_raw, "csr_version", exact = TRUE)
  col_idx_base  <- attr(all_connectivity_raw, "col_idx_base", exact = TRUE)
  col_idx_space <- attr(all_connectivity_raw, "col_idx_space", exact = TRUE)

  if (!is.null(csr_version) && csr_version != 2L) {
    stop("Unexpected csr_version attribute in all_connectivity.rds: ", csr_version)
  }

  if (!is.null(col_idx_base) && col_idx_base != "0-based") {
    stop("Unexpected col_idx_base attribute in all_connectivity.rds: ", col_idx_base)
  }

  if (!is.null(col_idx_space) && col_idx_space != "pu_local_index") {
    stop("Unexpected col_idx_space attribute in all_connectivity.rds: ", col_idx_space)
  }

  pu_graphs_by_key <- convert_priority_connectivity_graphs(
    all_connectivity_raw,
    retained_species = species_params$species
  )

  rm(all_connectivity_raw)
  gc(FALSE)


  # ---------------------------------------------------------------------
  # Internal validation checks
  # ---------------------------------------------------------------------
  if (anyDuplicated(patch_table[, .(species, patch_id)]) > 0L) {
    stop("patch_table contains duplicate (species, patch_id) rows.")
  }

  if (!all(unique(patch_table$species) %in% species_params$species)) {
    stop("Some species in patch_table are missing from species_params.")
  }

  if (!all(unique(patch_table$species) %in% species_names)) {
    stop("Some species in patch_table are missing from the retained patch rasters.")
  }

  # ---------------------------------------------------------------------
  # Save the reusable curve-neutral initialization before any pruning call
  # ---------------------------------------------------------------------
  #
  # This bundle captures the clean starting state before any cell-removal step.
  # The two environments are saved as named lists, then rebuilt on load. That is
  # simpler and more robust than saving the environments directly.
  #
  dir.create(initialization_bundle_dir, recursive = TRUE, showWarnings = FALSE)

  patch_id_species_names <- ls(envir = patch_id_by_species_env, all.names = TRUE)
  patch_id_by_species_list <- setNames(
    lapply(patch_id_species_names, function(species_name) {
      get(species_name, envir = patch_id_by_species_env, inherits = FALSE)
    }),
    patch_id_species_names
  )

  patch_index_species_names <- ls(envir = patch_cell_index_by_species_env, all.names = TRUE)
  patch_cell_index_by_species_list <- setNames(
    lapply(patch_index_species_names, function(species_name) {
      get(species_name, envir = patch_cell_index_by_species_env, inherits = FALSE)
    }),
    patch_index_species_names
  )

  initialization_bundle <- list(
    metadata = list(
      schema_version            = priority_initialization_schema_version(),
      analysis_contract         = unclass(contract),
      created_at                = format(Sys.time(), "%Y-%m-%d %H:%M:%OS6 %z"),
      taxa_tag                  = taxa_tag,
      sdm_source_tag            = sdm,
      retained_species          = species_params$species,
      retained_species_table    = species_params[
        ,
        .(
          species          = species,
          taxon_class      = taxon_class,
          redlist_category = redlist_category,
          sdm_method       = sdm_method
        )
      ]
    ),
    patch_table                      = patch_table,
    pu_graphs_by_key                 = pu_graphs_by_key,
    alive_species_count_by_cell      = alive_species_count_by_cell,
    cell_area_by_cell                = cell_area_by_cell,
    species_params                   = species_params,
    rook_neighbor_pairs              = rook_neighbor_pairs,
    patch_id_by_species_list         = patch_id_by_species_list,
    patch_cell_index_by_species_list = patch_cell_index_by_species_list
  )

  tmp_bundle_path <- tempfile(
    pattern = paste0(basename(initialization_bundle_path), "."),
    tmpdir = dirname(initialization_bundle_path),
    fileext = ".tmp"
  )
  backup_bundle_path <- paste0(initialization_bundle_path, ".bak")
  on.exit(unlink(tmp_bundle_path, force = TRUE), add = TRUE)
  assert(
    !file.exists(backup_bundle_path),
    paste0("Unresolved Stage 6 initialization backup exists: ", backup_bundle_path)
  )

  saveRDS(initialization_bundle, tmp_bundle_path)
  assert(
    file.exists(tmp_bundle_path) && file.info(tmp_bundle_path)$size[[1L]] > 0,
    paste0("Failed to write temporary Stage 6 initialization: ", tmp_bundle_path)
  )

  had_bundle <- file.exists(initialization_bundle_path)
  if (had_bundle) {
    assert(
      file.rename(initialization_bundle_path, backup_bundle_path),
      paste0("Could not back up existing Stage 6 initialization: ", initialization_bundle_path)
    )
  }

  committed <- file.rename(tmp_bundle_path, initialization_bundle_path)
  if (!isTRUE(committed)) {
    if (had_bundle && file.exists(backup_bundle_path)) {
      restored <- file.rename(backup_bundle_path, initialization_bundle_path)
      if (!isTRUE(restored)) {
        stop(
          "Could not commit or restore the Stage 6 initialization. Backup retained at: ",
          backup_bundle_path,
          call. = FALSE
        )
      }
    }
    stop("Could not commit Stage 6 initialization: ", initialization_bundle_path, call. = FALSE)
  }
  unlink(backup_bundle_path, force = TRUE)

  initialization_elapsed <- proc.time()[["elapsed"]] - initialization_started
  runtime_log_event(
    "priority_initialization_saved",
    path = normalizePath(initialization_bundle_path, mustWork = FALSE),
    species = nrow(species_params),
    cells = length(cell_area_by_cell),
    patches = nrow(patch_table),
    pus = data.table::uniqueN(patch_table, by = c("species", "pu_id")),
    graphs = length(pu_graphs_by_key),
    size_bytes = file.info(initialization_bundle_path)$size[[1L]],
    elapsed_seconds = round(initialization_elapsed, 2)
  )

  invisible(list(
    initialization_bundle_path = initialization_bundle_path,
    retained_species = nrow(species_params),
    initial_alive_cells = sum(alive_species_count_by_cell > 0L),
    patch_count = nrow(patch_table),
    pu_count = data.table::uniqueN(patch_table, by = c("species", "pu_id")),
    elapsed_seconds = initialization_elapsed
  ))
}
