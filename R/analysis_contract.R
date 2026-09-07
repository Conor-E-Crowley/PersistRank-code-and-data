# Cross-stage scientific and artifact contracts.
#
# The numeric constants below define the canonical manuscript analysis used by
# the spatial pipeline. Stage 2 may deliberately use another simulation horizon
# or quasi-extinction threshold through its public YAML. Those values are
# written to every persistence-point row and must be passed explicitly by
# consumers; this module validates mismatches but never silently overrides an
# artifact's recorded settings.


new_analysis_contract <- function(
  persistence_horizon_years = 1000L,
  quasi_extinction_abundance = 500L,
  minimum_patch_abundance = 10L
) {
  positive_integer <- function(value, label) {
    value <- suppressWarnings(as.numeric(value))
    assert(
      length(value) == 1L && is.finite(value) && value > 0 &&
        value == floor(value),
      paste0(label, " must be one positive integer.")
    )
    as.integer(value)
  }
  structure(
    list(
      persistence_horizon_years = positive_integer(
        persistence_horizon_years, "persistence_horizon_years"
      ),
      quasi_extinction_abundance = positive_integer(
        quasi_extinction_abundance, "quasi_extinction_abundance"
      ),
      minimum_patch_abundance = positive_integer(
        minimum_patch_abundance, "minimum_patch_abundance"
      )
    ),
    class = c("analysis_contract", "list")
  )
}

validate_analysis_contract <- function(contract, label = "analysis contract") {
  assert(is.list(contract), paste0(label, " must be a list."))
  required <- c(
    "persistence_horizon_years", "quasi_extinction_abundance",
    "minimum_patch_abundance"
  )
  assert(
    all(required %in% names(contract)),
    paste0(label, " must contain: ", paste(required, collapse = ", "), ".")
  )
  out <- new_analysis_contract(
    contract$persistence_horizon_years,
    contract$quasi_extinction_abundance,
    contract$minimum_patch_abundance
  )
  out
}

canonical_analysis_contract <- function() new_analysis_contract()

quasi_extinction_abundance <- function(contract = NULL) {
  if (is.null(contract)) contract <- canonical_analysis_contract()
  validate_analysis_contract(contract)$quasi_extinction_abundance
}

minimum_patch_abundance <- function(contract = NULL) {
  if (is.null(contract)) contract <- canonical_analysis_contract()
  validate_analysis_contract(contract)$minimum_patch_abundance
}

persistence_horizon_years <- function(contract = NULL) {
  if (is.null(contract)) contract <- canonical_analysis_contract()
  validate_analysis_contract(contract)$persistence_horizon_years
}

priority_initialization_schema_version <- function() {
  5L
}

persistence_quantiles <- function() {
  c(q50 = 0.50, q16 = 0.16, q025 = 0.025, q84 = 0.84, q975 = 0.975)
}

persistence_curves <- function() {
  names(persistence_quantiles())
}

persistence_coefficient_columns <- function(curves = persistence_curves()) {
  curves <- as.character(curves)
  unlist(lapply(
    curves,
    function(curve) c(paste0("alpha_", curve), paste0("beta_", curve))
  ), use.names = FALSE)
}

main_persistence_curve <- function() {
  "q50"
}

persistence_point_metadata_columns <- function() {
  c(
    "curve_probability",
    "persistence_horizon_years",
    "quasi_extinction_abundance",
    "cap_factor",
    "r_buffer",
    "n_draws",
    "reps",
    "chunk_size",
    "base_seed",
    "posterior_rm_seed",
    "posterior_sigma_seed",
    "demographic_uncertainty",
    "residual_growth_seed",
    "residual_environmental_variation_seed",
    "bird_model_signature",
    "grid_signature",
    "posterior_sampling",
    "posterior_rm_md5",
    "posterior_sigma_md5",
    "simulator_contract"
  )
}

validate_persistence_point_metadata <- function(
  x,
  label = "persistence-point table",
  expected_horizon = NULL,
  expected_quasi_extinction = NULL
) {
  need_cols(x, c("curve", persistence_point_metadata_columns()), label)
  assert(NROW(x) > 0L, paste0(label, " contains no rows."))

  validate_expected_integer <- function(value, arg_label) {
    if (is.null(value)) return(NULL)
    value <- suppressWarnings(as.numeric(value))
    assert(
      length(value) == 1L &&
        is.finite(value) &&
        value > 0 &&
        value == floor(value),
      paste0(arg_label, " must be a single positive integer.")
    )
    as.integer(value)
  }

  expected_horizon <- validate_expected_integer(expected_horizon, "expected_horizon")
  expected_quasi_extinction <- validate_expected_integer(
    expected_quasi_extinction,
    "expected_quasi_extinction"
  )

  metadata_number <- function(col) {
    suppressWarnings(as.numeric(x[[col]]))
  }

  validate_integer_metadata <- function(col, expected = NULL, minimum = NULL) {
    value <- metadata_number(col)
    ok <- is.finite(value) & value == floor(value)
    if (!is.null(expected)) ok <- ok & value == expected
    if (!is.null(minimum)) ok <- ok & value >= minimum
    expectation <- if (!is.null(expected)) {
      paste0("integer value ", expected)
    } else if (!is.null(minimum)) {
      paste0("integer value >= ", minimum)
    } else {
      "integer value"
    }
    assert(
      all(ok),
      paste0(label, " contains invalid ", col, " metadata; expected ", expectation, ".")
    )
    invisible(value)
  }

  curve <- as.character(x$curve)
  expected_prob <- unname(persistence_quantiles()[curve])
  curve_probability <- metadata_number("curve_probability")

  bad_prob <- !is.finite(expected_prob) |
    !is.finite(curve_probability) |
    abs(curve_probability - expected_prob) > sqrt(.Machine$double.eps)
  assert(
    !any(bad_prob),
    paste0(label, " has curve_probability values inconsistent with analysis_contract.R.")
  )

  validate_integer_metadata(
    "persistence_horizon_years",
    expected = expected_horizon,
    minimum = 1L
  )
  validate_integer_metadata(
    "quasi_extinction_abundance",
    expected = expected_quasi_extinction,
    minimum = 1L
  )

  positive_cols <- c("cap_factor", "r_buffer")
  for (col in positive_cols) {
    value <- metadata_number(col)
    assert(
      all(is.finite(value) & value > 0),
      paste0(label, " contains invalid ", col, " metadata.")
    )
  }

  for (col in c("n_draws", "reps", "chunk_size")) {
    validate_integer_metadata(col, minimum = 1L)
  }

  for (col in c("base_seed", "posterior_rm_seed", "posterior_sigma_seed")) {
    validate_integer_metadata(col, minimum = 0L)
  }
  uncertainty <- trimws(as.character(x$demographic_uncertainty))
  assert(
    all(uncertainty %in% c("coefficients_only", "posterior_predictive")),
    paste0(label, " contains invalid demographic_uncertainty metadata.")
  )
  for (col in c(
    "residual_growth_seed",
    "residual_environmental_variation_seed"
  )) {
    value <- suppressWarnings(as.numeric(x[[col]]))
    expected_missing <- uncertainty == "coefficients_only"
    assert(
      all((expected_missing & is.na(value)) |
            (!expected_missing & is.finite(value) & value >= 0 & value == floor(value))),
      paste0(label, " contains invalid ", col, " metadata.")
    )
  }

  validate_text_metadata <- function(col, expected = NULL) {
    value <- trimws(as.character(x[[col]]))
    ok <- !is.na(value) & nzchar(value)
    if (!is.null(expected)) ok <- ok & value == expected
    assert(all(ok), paste0(label, " contains invalid ", col, " metadata."))
  }
  for (col in c(
    "bird_model_signature", "grid_signature", "posterior_rm_md5",
    "posterior_sigma_md5", "simulator_contract"
  )) {
    validate_text_metadata(col)
  }
  validate_text_metadata("posterior_sampling", "independent_without_replacement")

  invisible(TRUE)
}

validate_persistence_curve <- function(curve, label = "curve") {
  curve <- validate_scalar_string(curve, label)
  assert(
    curve %in% persistence_curves(),
    paste0(
      label, " must be one of: ",
      paste(persistence_curves(), collapse = ", "), "."
    )
  )
  curve
}

benchmark_rank_labels <- function() {
  c(
    abf = "Zonation ABF",
    caz1 = "Zonation CAZ1",
    caz2 = "Zonation CAZ2",
    cazmax = "Zonation CAZMAX"
  )
}

benchmark_rank_methods <- function() {
  names(benchmark_rank_labels())
}

validate_benchmark_rank_method <- function(method, label = "rank_method") {
  validate_scalar_choice(method, benchmark_rank_methods(), label)
}

benchmark_rank_label <- function(method) {
  method <- validate_benchmark_rank_method(method)
  unname(benchmark_rank_labels()[method])
}

area_exceeds_threshold <- function(area, threshold) {
  is.finite(area) & is.finite(threshold) & area > threshold
}

validate_area_threshold <- function(threshold, label = "area threshold") {
  assert(
    length(threshold) == 1L && is.finite(threshold) && threshold > 0,
    paste0(label, " must be one finite positive value.")
  )
  as.numeric(threshold)
}

validate_density_area_thresholds <- function(
  density,
  min_patch_area_km2,
  min_population_area_km2,
  label = "species thresholds",
  tolerance = sqrt(.Machine$double.eps),
  contract = NULL
) {
  if (is.null(contract)) contract <- canonical_analysis_contract()
  contract <- validate_analysis_contract(contract)
  expected_patch <- minimum_patch_abundance(contract) / density
  expected_population <- quasi_extinction_abundance(contract) / density

  ok <- all(is.finite(density) & density > 0) &&
    all(is.finite(min_patch_area_km2) & min_patch_area_km2 > 0) &&
    all(is.finite(min_population_area_km2) & min_population_area_km2 > 0) &&
    isTRUE(all.equal(
      as.numeric(min_patch_area_km2),
      as.numeric(expected_patch),
      tolerance = tolerance,
      check.attributes = FALSE
    )) &
    isTRUE(all.equal(
      as.numeric(min_population_area_km2),
      as.numeric(expected_population),
      tolerance = tolerance,
      check.attributes = FALSE
    ))

  assert(
    ok,
    paste0(
      label, " are inconsistent with density and the requested abundance ",
      "thresholds (", minimum_patch_abundance(contract), " and ",
      quasi_extinction_abundance(contract), " individuals)."
    )
  )
  invisible(TRUE)
}

add_canonical_area_thresholds <- function(
  x,
  density_col = "density",
  min_patch_col = "min_patch_size",
  min_population_col = "min_pop_size",
  validate = TRUE,
  label = "species-table area thresholds",
  contract = NULL
) {
  need_cols(x, c(density_col, min_patch_col, min_population_col), label)

  x$min_patch_area_km2 <- suppressWarnings(as.numeric(x[[min_patch_col]]))
  x$min_population_area_km2 <- suppressWarnings(as.numeric(x[[min_population_col]]))

  if (isTRUE(validate)) {
    validate_density_area_thresholds(
      density = suppressWarnings(as.numeric(x[[density_col]])),
      min_patch_area_km2 = x$min_patch_area_km2,
      min_population_area_km2 = x$min_population_area_km2,
      label = label,
      contract = contract
    )
  }

  x
}

validate_priority_initialization_schema <- function(
  bundle,
  label = "Stage 6 shared initialization"
) {
  assert(is.list(bundle), paste0(label, " must be an R list."))

  required_fields <- c(
    "metadata",
    "patch_table",
    "pu_graphs_by_key",
    "alive_species_count_by_cell",
    "cell_area_by_cell",
    "species_params",
    "rook_neighbor_pairs",
    "patch_id_by_species_list",
    "patch_cell_index_by_species_list"
  )
  missing_fields <- setdiff(required_fields, names(bundle))
  assert(
    !length(missing_fields),
    paste0(
      label,
      " is missing required top-level field(s): ",
      paste(missing_fields, collapse = ", "),
      ". Regenerate the shared initialization with 6_spatial_prioritization_pipeline.Rmd."
    )
  )

  metadata <- bundle$metadata
  assert(is.list(metadata), paste0(label, " metadata must be an R list."))

  found <- suppressWarnings(as.numeric(metadata$schema_version))
  expected <- priority_initialization_schema_version()
  found_ok <- length(found) == 1L && is.finite(found) &&
    found == floor(found) && found == expected
  assert(
    found_ok,
    paste0(
      "Unsupported ", label, " schema. Expected schema_version ", expected,
      "; found ",
      if (length(found) == 1L && is.finite(found)) found else "missing or non-numeric schema_version",
      ". Regenerate the shared initialization with 6_spatial_prioritization_pipeline.Rmd."
    )
  )

  required_metadata_fields <- c(
    "schema_version",
    "created_at",
    "taxa_tag",
    "sdm_source_tag",
    "retained_species",
    "retained_species_table"
  )
  missing_metadata_fields <- setdiff(required_metadata_fields, names(metadata))
  assert(
    !length(missing_metadata_fields),
    paste0(
      label,
      " metadata is missing required field(s): ",
      paste(missing_metadata_fields, collapse = ", "),
      ". Regenerate the shared initialization with 6_spatial_prioritization_pipeline.Rmd."
    )
  )

  assert(
    "analysis_contract" %in% names(metadata),
    paste0(label, " schema ", expected, " requires metadata$analysis_contract.")
  )
  validate_analysis_contract(
    metadata$analysis_contract,
    paste0(label, " metadata$analysis_contract")
  )
  species_columns <- c(
    "species", "taxon_class", "redlist_category", "sdm_method", "density",
    "dispersal_distance_km", "min_patch_area_km2",
    "min_population_area_km2", persistence_coefficient_columns()
  )
  missing_species_columns <- setdiff(species_columns, names(bundle$species_params))
  assert(
    !length(missing_species_columns),
    paste0(
      label, " species_params is missing shared coefficient/base field(s): ",
      paste(missing_species_columns, collapse = ", "), "."
    )
  )
  assert(
    !any(c("a_pred", "b_pred") %in% names(bundle$species_params)),
    paste0(
      label,
      " contains obsolete curve-selected species parameters. Reinitialize Stage 6."
    )
  )
  invisible(TRUE)
}
