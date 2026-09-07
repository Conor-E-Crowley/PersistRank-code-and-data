# Shared Stage 6–7 selectors, artifact identities, and serialized contracts.
#
# Taxon and SDM tags identify the shared initialization; curve, removal-batch,
# and pruning-stage tags identify only run outputs. Initialization metadata and
# rank-LUT validators remain here because Stage 7 consumes the same immutable
# identities and ordering contracts.


priority_flags_from_tags <- function(taxa, sdm) {
  taxa <- validate_priority_taxa_tag(taxa)
  sdm <- validate_priority_sdm_tag(sdm)
  list(
    do_mammals = taxa %in% c("mammals", "mammals_birds"),
    do_birds = taxa %in% c("birds", "mammals_birds"),
    do_ppm = sdm %in% c("ppm", "ppm_rangebag"),
    do_rangebag = sdm %in% c("rangebag", "ppm_rangebag")
  )
}

validate_priority_count <- function(x, label) {
  assert(
    is.numeric(x) &&
      length(x) == 1L &&
      is.finite(x) &&
      x > 0 &&
      x == floor(x),
    paste0(label, " must be a single positive integer.")
  )
  as.integer(x)
}

validate_priority_cell_ids <- function(x, n_cells, label = "cell IDs") {
  n_cells <- validate_priority_count(n_cells, "n_cells")
  if (!length(x)) return(integer())

  values <- suppressWarnings(as.numeric(x))
  invalid <- !is.finite(values) |
    values < 1 |
    values > n_cells |
    values != floor(values)

  if (any(invalid)) {
    bad <- which(invalid)
    preview <- utils::head(bad, 10L)
    details <- paste0(
      "position=", preview,
      " value=", format(values[preview], trim = TRUE, scientific = FALSE)
    )
    stop(
      label, " contains ", length(bad),
      " invalid raster-cell ID(s); expected finite integers in 1:", n_cells,
      ". First invalid values: ", paste(details, collapse = "; "),
      call. = FALSE
    )
  }

  sort(unique(as.integer(values)))
}

validate_priority_cell_state <- function(
  alive_species_count_by_cell,
  cell_area_by_cell,
  removal_order_by_cell,
  retained_species_count,
  label = "Stage 6 cell state"
) {
  retained_species_count <- validate_priority_count(
    retained_species_count,
    "retained_species_count"
  )
  n_cells <- length(alive_species_count_by_cell)
  assert(n_cells > 0L, paste0(label, " is empty."))
  assert(
    length(cell_area_by_cell) == n_cells && length(removal_order_by_cell) == n_cells,
    paste0(label, " vectors must have equal lengths.")
  )
  assert(
    is.numeric(alive_species_count_by_cell),
    paste0(label, "$alive_species_count_by_cell must be numeric.")
  )
  assert(
    is.numeric(cell_area_by_cell),
    paste0(label, "$cell_area_by_cell must be numeric.")
  )
  assert(
    is.numeric(removal_order_by_cell),
    paste0(label, "$removal_order_by_cell must be numeric or integer.")
  )

  invalid_alive <- !is.finite(alive_species_count_by_cell) |
    alive_species_count_by_cell < 0 |
    alive_species_count_by_cell > retained_species_count |
    alive_species_count_by_cell != floor(alive_species_count_by_cell)
  assert(
    !any(invalid_alive),
    paste0(
      label,
      "$alive_species_count_by_cell must contain finite integer counts in [0, ",
      retained_species_count, "]."
    )
  )

  invalid_order <- is.nan(removal_order_by_cell) |
    (!is.na(removal_order_by_cell) & (
      !is.finite(removal_order_by_cell) |
        removal_order_by_cell <= 0 |
        removal_order_by_cell != floor(removal_order_by_cell)
    ))
  assert(
    !any(invalid_order),
    paste0(label, "$removal_order_by_cell must contain positive integers or NA.")
  )

  current_alive_by_cell <- alive_species_count_by_cell > 0
  assert(
    !any(current_alive_by_cell & !is.na(removal_order_by_cell)),
    paste0(label, " contains cells that are simultaneously alive and already removed.")
  )
  initial_alive_by_cell <- current_alive_by_cell | !is.na(removal_order_by_cell)
  assert(any(initial_alive_by_cell), paste0(label, " contains no initially alive cells."))

  initial_cell_area <- cell_area_by_cell[initial_alive_by_cell]
  assert(
    all(is.finite(initial_cell_area) & initial_cell_area > 0),
    paste0(label, " contains non-positive or non-finite area within the initial domain.")
  )

  list(
    alive_species_count_by_cell = if (is.integer(alive_species_count_by_cell)) {
      alive_species_count_by_cell
    } else {
      as.integer(alive_species_count_by_cell)
    },
    removal_order_by_cell = if (is.integer(removal_order_by_cell)) {
      removal_order_by_cell
    } else {
      as.integer(removal_order_by_cell)
    },
    initial_alive_by_cell = initial_alive_by_cell,
    initial_alive_cell_count = as.integer(sum(initial_alive_by_cell)),
    current_alive_cell_count = as.integer(sum(current_alive_by_cell)),
    initial_alive_area_km2 = as.numeric(sum(initial_cell_area))
  )
}

validate_priority_species_parameters <- function(
  species_params,
  label = "Stage 6 species_params",
  contract = canonical_analysis_contract()
) {
  contract <- validate_analysis_contract(contract, paste0(label, " abundance contract"))
  assert(is.data.frame(species_params), paste0(label, " must be a data.frame or data.table."))
  required <- c(
    "species", "density", "dispersal_distance_km", "min_patch_area_km2",
    "min_population_area_km2", "a_pred", "b_pred"
  )
  need_cols(species_params, required, label)
  assert(nrow(species_params) > 0L, paste0(label, " is empty."))

  species <- trimws(as.character(species_params$species))
  assert(
    all(!is.na(species)) && all(nzchar(species)) && !anyDuplicated(species),
    paste0(label, "$species must contain unique nonblank names.")
  )

  numeric_columns <- required[-1L]
  invalid_by_column <- lapply(numeric_columns, function(column) {
    values <- suppressWarnings(as.numeric(species_params[[column]]))
    !is.finite(values) | values <= 0
  })
  names(invalid_by_column) <- numeric_columns
  invalid_rows <- Reduce(`|`, invalid_by_column)
  if (any(invalid_rows)) {
    bad_species <- utils::head(species[invalid_rows], 20L)
    bad_columns <- names(invalid_by_column)[vapply(invalid_by_column, any, logical(1L))]
    stop(
      label, " contains non-positive or non-finite values in: ",
      paste(bad_columns, collapse = ", "),
      ". First affected species: ", paste(bad_species, collapse = ", "),
      call. = FALSE
    )
  }

  validate_density_area_thresholds(
    density = species_params$density,
    min_patch_area_km2 = species_params$min_patch_area_km2,
    min_population_area_km2 = species_params$min_population_area_km2,
    label = paste0(label, " area thresholds"),
    contract = contract
  )
  invisible(TRUE)
}

# Materialize the compact curve-specific view consumed by the unchanged Stage 6
# optimizer. The shared artifact retains all coefficients; only this small
# species table is copied and receives the established a_pred/b_pred names.
priority_species_parameters_for_curve <- function(
  species_params,
  curve,
  label = "Stage 6 shared species_params",
  contract = canonical_analysis_contract()
) {
  curve <- unname(as.character(validate_persistence_curve(curve, "curve")))
  alpha_column <- paste0("alpha_", curve)
  beta_column <- paste0("beta_", curve)
  need_cols(
    species_params,
    c(alpha_column, beta_column),
    paste0(label, " coefficients")
  )
  selected <- data.table::copy(data.table::as.data.table(species_params))
  selected[, `:=`(
    a_pred = suppressWarnings(as.numeric(get(alpha_column))),
    b_pred = suppressWarnings(as.numeric(get(beta_column)))
  )]
  validate_priority_species_parameters(
    selected,
    paste0(label, " [", curve, "]"),
    contract = contract
  )
  data.table::setkey(selected, species)
  selected[]
}

# Fingerprint only the selected compact coefficient pair. This identity is
# stored in the curve manifest/checkpoints without duplicating initialized
# spatial state or rehashing the large shared artifact.
priority_curve_coefficient_identity <- function(
  species_params,
  curve,
  contract = canonical_analysis_contract()
) {
  curve <- unname(as.character(validate_persistence_curve(curve, "curve")))
  selected <- priority_species_parameters_for_curve(
    species_params, curve, contract = contract
  )
  ordered <- order(as.character(selected$species))
  identity <- list(
    curve = curve,
    species = as.character(selected$species[ordered]),
    alpha = as.numeric(selected$a_pred[ordered]),
    beta = as.numeric(selected$b_pred[ordered])
  )
  path <- tempfile("stage6_coefficients_", fileext = ".rds")
  on.exit(unlink(path, force = TRUE), add = TRUE)
  saveRDS(identity, path, version = 3L)
  list(
    curve = curve,
    alpha_column = paste0("alpha_", curve),
    beta_column = paste0("beta_", curve),
    md5 = unname(tools::md5sum(path))
  )
}

validate_max_stages <- function(x, label = "max_stages") {
  assert(
    is.numeric(x) &&
      length(x) == 1L &&
      !is.na(x) &&
      x > 0 &&
      (is.infinite(x) || x == floor(x)),
    paste0(label, " must be a single positive integer or Inf.")
  )
  if (is.finite(x)) as.integer(x) else Inf
}

validate_priority_taxa_tag <- function(x, label = "taxa") {
  x <- validate_scalar_string(x, label)
  assert(
    x %in% c("mammals", "birds", "mammals_birds"),
    paste0(label, " must be one of: mammals, birds, mammals_birds.")
  )
  x
}

validate_priority_public_taxa <- function(x, label = "taxa") {
  selection <- validate_taxa_selector(x, label)
  list(
    taxa = selection$value,
    taxa_tag = switch(
      selection$value,
      mammals = "mammals",
      birds = "birds",
      both = "mammals_birds"
    ),
    selected_mammals = selection$selected_mammals,
    selected_birds = selection$selected_birds
  )
}

validate_priority_sdm_tag <- function(x, label = "sdm") {
  validate_scalar_choice(x, c("ppm", "rangebag", "ppm_rangebag"), label)
}

validate_priority_initialization_metadata <- function(
  bundle,
  taxa,
  sdm,
  label = "Stage 6 shared initialization",
  contract = canonical_analysis_contract()
) {
  metadata <- bundle$metadata
  assert(!is.null(metadata), paste0(label, " is missing metadata."))

  taxa <- validate_priority_taxa_tag(taxa, "expected taxa")
  sdm <- validate_priority_sdm_tag(sdm, "expected sdm")
  contract <- validate_analysis_contract(contract, paste0(label, " expected abundance contract"))

  required_fields <- c(
    "taxa_tag",
    "sdm_source_tag",
    "retained_species"
  )
  missing_fields <- setdiff(required_fields, names(metadata))
  assert(
    !length(missing_fields),
    paste0(
      label,
      " metadata is missing field(s): ",
      paste(missing_fields, collapse = ", "),
      ". Regenerate the shared initialization with 6_spatial_prioritization_pipeline.Rmd."
    )
  )

  found_taxa <- validate_priority_taxa_tag(
    metadata$taxa_tag,
    paste0(label, " metadata$taxa_tag")
  )
  assert(
    identical(found_taxa, taxa),
    paste0(
      label,
      " metadata$taxa_tag mismatch. Expected ",
      taxa,
      " but initialization contains ",
      metadata$taxa_tag,
      "."
    )
  )

  found_sdm <- validate_priority_sdm_tag(metadata$sdm_source_tag, paste0(label, " metadata$sdm_source_tag"))
  assert(
    identical(found_sdm, sdm),
    paste0(
      label,
      " metadata$sdm_source_tag mismatch. Expected ",
      sdm,
      " but the shared initialization contains ",
      metadata$sdm_source_tag,
      "."
    )
  )
  found_contract <- validate_analysis_contract(
    metadata$analysis_contract,
    paste0(label, " metadata$analysis_contract")
  )
  assert(
    identical(found_contract, contract),
    paste0(label, " abundance contract does not match the requested analysis.")
  )

  retained_species <- trimws(as.character(metadata$retained_species))
  assert(
    length(retained_species) > 0L &&
      all(!is.na(retained_species)) &&
      all(nzchar(retained_species)),
    paste0(label, " metadata$retained_species must contain at least one species name.")
  )

  for (curve in persistence_curves()) {
    priority_species_parameters_for_curve(
      bundle$species_params,
      curve,
      label = paste0(label, " species_params"),
      contract = contract
    )
  }

  invisible(retained_species)
}

zonation_rankmap_path <- function(zonation_root, rank_method) {
  rank_method <- validate_benchmark_rank_method(rank_method)
  zonation_root <- validate_path_param(zonation_root, "zonation_root")
  file.path(zonation_root, toupper(rank_method), "rankmap.tif")
}

empty_rank_lut <- function() {
  data.table::data.table(
    stage = integer(),
    method = character(),
    scientificName = character(),
    species = character(),
    patch_id = integer(),
    pu_id = integer(),
    patch_area_km2 = numeric(),
    patch_n_cells = integer()
  )
}

normalized_rank_lut_columns <- function() {
  c("method", "stage", "scientificName", "species", "patch_id", "pu_id", "patch_area_km2")
}

empty_normalized_rank_lut <- function() {
  empty_rank_lut()[, normalized_rank_lut_columns(), with = FALSE]
}

standardize_rank_lut <- function(x, method, stage, species_params) {
  dt <- data.table::as.data.table(x)

  if (!nrow(dt)) {
    return(empty_normalized_rank_lut())
  }

  if (!"scientificName" %in% names(dt)) {
    assert("species" %in% names(dt),
           "LUT must contain either scientificName or species.")

    species_raw <- as.character(dt$species)

    if (all(species_raw %in% species_params$scientificName)) {
      dt[, scientificName := species_raw]
    } else if (all(species_raw %in% species_params$species)) {
      map <- species_params[, .(species, scientificName)]
      dt[, species := species_raw]
      dt <- merge(dt, map, by = "species", all.x = TRUE, sort = FALSE)
      assert(all(!is.na(dt$scientificName)),
             "Could not map all species IDs in LUT to scientificName.")
    } else {
      stop(
        "LUT species column is neither scientific names nor sid(scientificName) IDs.",
        call. = FALSE
      )
    }
  }

  assert("pu_id" %in% names(dt), "LUT is missing pu_id.")
  assert("patch_area_km2" %in% names(dt), "LUT is missing patch_area_km2.")

  if (!"patch_id" %in% names(dt)) {
    dt[, patch_id := seq_len(.N)]
  }

  dt[, scientificName := trimws(as.character(scientificName))]
  dt[, species := species_id(scientificName)]
  dt <- dt[species %in% species_params$species]

  if (!nrow(dt)) {
    return(empty_normalized_rank_lut())
  }

  patch_id_num <- suppressWarnings(as.numeric(dt$patch_id))
  pu_id_num <- suppressWarnings(as.numeric(dt$pu_id))
  bad_ids <- !is.finite(patch_id_num) |
    patch_id_num <= 0 |
    patch_id_num != floor(patch_id_num) |
    !is.finite(pu_id_num) |
    pu_id_num <= 0 |
    pu_id_num != floor(pu_id_num)
  assert(
    !any(bad_ids),
    paste0(
      "LUT for method=", method, ", stage=", stage,
      " contains invalid patch_id or pu_id; both must be positive integers."
    )
  )

  dt[, patch_id := as.integer(patch_id_num)]
  dt[, pu_id := as.integer(pu_id_num)]
  dt[, patch_area_km2 := as.numeric(patch_area_km2)]

  bad <- dt[
    !is.finite(patch_area_km2) |
      patch_area_km2 <= 0
  ]

  assert(nrow(bad) == 0L,
         paste0("LUT for method=", method, ", stage=", stage,
                " contains invalid patch_area_km2."))

  dt[, .(
    method = as.character(method),
    stage = as.integer(stage),
    scientificName,
    species,
    patch_id,
    pu_id,
    patch_area_km2
  )]
}
