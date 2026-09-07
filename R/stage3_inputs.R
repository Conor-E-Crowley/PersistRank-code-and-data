# Stage 3 persistence-point handoff and operational prerequisites.
#
# Loaded by stage3_workflow.R after pure configuration. Sourcing is definition
# only; preflight, package attachment, Stage 2 artifact reads, provenance checks,
# and input-derived settings occur only when their functions are called.

load_gompertz_packages <- function(include_figures = TRUE) {
  load_packages(required_gompertz_packages(include_figures))
}

check_stage3_preflight <- function(config) {
  if (config$selected_mammals) need_file(config$paths$mammals_points, "mammal persistence points")
  if (config$selected_birds) need_file(config$paths$birds_points, "bird persistence points")
  if (isTRUE(config$write_figures)) {
    ensure_writable_dir(config$paths$figure_dir, "Stage 3 figure directory")
    ensure_writable_dir(config$paths$si_dir, "Stage 3 SI figure directory")
  }

  ensure_writable_dir(dirname(config$paths$loess_models_rds), "Stage 3 model directory")
  if (config$mode == "reuse") {
    need_file(config$paths$loess_models_rds, "saved Stage 3 persistence-curve model")
  }
  invisible(TRUE)
}

stage3_numeric_equal <- function(x, y, tolerance = sqrt(.Machine$double.eps)) {
  x <- suppressWarnings(as.numeric(x)); y <- suppressWarnings(as.numeric(y))
  length(x) == length(y) && all(is.finite(x)) && all(is.finite(y)) &&
    all(abs(x - y) <= tolerance * pmax(1, abs(x), abs(y)))
}

stage3_unique_metadata <- function(x, columns, label) {
  out <- list()
  for (column in columns) {
    values <- unique(x[[column]])
    values <- values[!is.na(values)]
    assert(length(values) <= 1L, paste0(label, " must contain one consistent ", column, " value."))
    out[[column]] <- if (length(values)) values[[1L]] else NA
  }
  out
}

stage3_shared_provenance_columns <- function() {
  c(
    "persistence_horizon_years", "quasi_extinction_abundance", "cap_factor", "r_buffer",
    "n_draws", "reps", "chunk_size", "base_seed", "grid_signature", "posterior_sampling",
    "demographic_uncertainty", "simulator_contract"
  )
}

load_stage3_points <- function(path, taxon = c("mammal", "bird"), config) {
  taxon <- match.arg(taxon)
  label <- paste0(taxon, " persistence points")
  need_file(path, label)
  raw <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)

  idx_col <- if (taxon == "mammal") "mass_idx" else "genlength_idx"
  value_col <- if (taxon == "mammal") "Mass_g" else "GenLength"
  extra_cols <- if (taxon == "bird") "bird_sigma_model_group" else character()
  required <- unique(c(
    idx_col, value_col, extra_cols, "curve", "curve_probability", "p_target", "K", "p_eval",
    persistence_point_metadata_columns()
  ))
  need_cols(raw, required, label)
  assert(nrow(raw) > 0L, paste0(label, " contains no rows."))

  input_horizon <- unique(suppressWarnings(as.numeric(raw$persistence_horizon_years)))
  input_k0 <- unique(suppressWarnings(as.numeric(raw$quasi_extinction_abundance)))
  assert(
    length(input_horizon) == 1L && is.finite(input_horizon) &&
      input_horizon > 0 && input_horizon == floor(input_horizon),
    paste0(label, " must contain one positive integer simulation horizon.")
  )
  assert(
    length(input_k0) == 1L && is.finite(input_k0) &&
      input_k0 > 0 && input_k0 == floor(input_k0),
    paste0(label, " must contain one positive integer quasi-extinction abundance.")
  )
  validate_persistence_point_metadata(
    raw, label, expected_horizon = input_horizon,
    expected_quasi_extinction = input_k0
  )

  points <- raw |>
    dplyr::transmute(
      group = if (taxon == "mammal") "Mammals" else "Birds",
      predictor = if (taxon == "mammal") "Mass_g" else "GenLength",
      bird_sigma_model_group = if (taxon == "bird") trimws(as.character(.data[[extra_cols]])) else NA_character_,
      trait_idx = suppressWarnings(as.integer(.data[[idx_col]])),
      trait_value = suppressWarnings(as.numeric(.data[[value_col]])),
      curve = trimws(as.character(.data$curve)),
      curve_probability = suppressWarnings(as.numeric(.data$curve_probability)),
      p_target = suppressWarnings(as.numeric(.data$p_target)),
      K = suppressWarnings(as.numeric(.data$K)),
      p_eval = suppressWarnings(as.numeric(.data$p_eval))
    )

  assert(!anyDuplicated(points), paste0(label, " contains duplicate persistence-point rows."))
  assert(all(is.finite(points$trait_idx) & points$trait_idx > 0 & points$trait_idx == floor(points$trait_idx)),
         paste0(label, " contains invalid trait indices."))
  assert(all(is.finite(points$trait_value) & points$trait_value > 0), paste0(label, " contains invalid trait values."))
  assert(identical(unique(points$curve), config$curves) || setequal(unique(points$curve), config$curves),
         paste0(label, " must contain exactly the five Stage 3 curves: ", paste(config$curves, collapse = ", "), "."))
  expected_probability <- unname(config$curve_probabilities[points$curve])
  assert(stage3_numeric_equal(points$curve_probability, expected_probability),
         paste0(label, " has curve probabilities inconsistent with the curve contract."))
  assert(all(is.finite(points$K) & points$K > input_k0), paste0(label, " contains K <= Nq or nonfinite K."))
  assert(all(is.finite(points$p_target) & points$p_target > 0 & points$p_target < 1),
         paste0(label, " contains invalid p_target values."))
  assert(all(is.finite(points$p_eval) & points$p_eval >= 0 & points$p_eval <= 1),
         paste0(label, " contains invalid evaluated persistence values."))

  if (taxon == "bird") {
    signature <- unique(trimws(as.character(raw$bird_model_signature)))
    assert(length(signature) == 1L, paste0(label, " has inconsistent bird model signatures."))
    bird_model <- bird_model_spec_from_signature(signature)
    validate_bird_model_group(points$bird_sigma_model_group, bird_model,
                              paste0(label, " bird_sigma_model_group"))
    assert(
      identical(sort(unique(points$bird_sigma_model_group)), sort(bird_model$model_groups)),
      paste0(label, " must contain every modeled bird group: ",
             paste(bird_model$model_groups, collapse = ", "), ".")
    )
  } else {
    bird_model <- NULL
  }

  key_cols <- c("trait_idx", "bird_sigma_model_group", "curve")
  ordered <- points |> dplyr::arrange(dplyr::across(dplyr::all_of(key_cols)), .data$p_target)
  block_check <- ordered |>
    dplyr::group_by(dplyr::across(dplyr::all_of(key_cols))) |>
    dplyr::summarise(
      rows = dplyr::n(), targets = dplyr::n_distinct(.data$p_target),
      target_grid = paste(format(.data$p_target, digits = 17), collapse = "|"),
      target_increases = all(diff(.data$p_target) > 0),
      k_increases = all(diff(.data$K) > 0),
      persistence_nondecreasing = all(diff(.data$p_eval) >= -sqrt(.Machine$double.eps)),
      .groups = "drop"
    )
  assert(all(block_check$rows >= 3L), paste0(label, " must contain at least three p_target rows per trait/curve block."))
  assert(all(block_check$targets == block_check$rows), paste0(label, " repeats p_target within a block."))
  assert(length(unique(block_check$target_grid)) == 1L, paste0(label, " uses inconsistent p_target grids across blocks."))
  assert(all(block_check$target_increases & block_check$k_increases & block_check$persistence_nondecreasing),
         paste0(label, " must be monotone in p_target, K, and evaluated persistence within every block."))

  trait_map <- points |> dplyr::distinct(.data$trait_idx, .data$bird_sigma_model_group, .data$trait_value)
  assert(!anyDuplicated(trait_map[c("trait_idx", "bird_sigma_model_group")]),
         paste0(label, " maps a trait index to multiple values."))
  trait_checks <- trait_map |>
    dplyr::arrange(.data$bird_sigma_model_group, .data$trait_idx) |>
    dplyr::group_by(.data$bird_sigma_model_group) |>
    dplyr::summarise(
      n_traits = dplyr::n(),
      contiguous_indices = identical(as.integer(.data$trait_idx), seq_len(dplyr::n())),
      increasing_values = all(diff(.data$trait_value) > 0), .groups = "drop"
    )
  assert(all(trait_checks$n_traits >= 5L), paste0(label, " must contain at least five trait values per model group."))
  assert(all(trait_checks$contiguous_indices), paste0(label, " trait indices must be contiguous and start at one."))
  assert(all(trait_checks$increasing_values), paste0(label, " trait values must increase with trait index."))

  minimum_loess_traits <- ceiling((config$loess$degree + 3) / config$loess$span)
  assert(all(trait_checks$n_traits >= minimum_loess_traits),
         paste0(label, " has too few trait values for finite LOESS uncertainty bands with span ",
                config$loess$span, ". At least ", minimum_loess_traits, " are required."))

  if (taxon == "bird") {
    reference_grid <- trait_map |>
      dplyr::filter(.data$bird_sigma_model_group == bird_model$model_groups[[1L]]) |>
      dplyr::arrange(.data$trait_idx)
    for (model_group in bird_model$model_groups[-1L]) {
      comparison <- trait_map |>
        dplyr::filter(.data$bird_sigma_model_group == .env$model_group) |>
        dplyr::arrange(.data$trait_idx)
      assert(identical(reference_grid$trait_idx, comparison$trait_idx) &&
               stage3_numeric_equal(reference_grid$trait_value, comparison$trait_value),
             paste0(label, " must use one generation-length grid for every bird branch."))
    }
  }

  metadata <- stage3_unique_metadata(
    raw, setdiff(persistence_point_metadata_columns(), "curve_probability"), label
  )
  provenance <- metadata[stage3_shared_provenance_columns()]
  list(
    points = points, metadata = metadata, provenance = provenance,
    source = normalizePath(path, winslash = "/", mustWork = TRUE),
    md5 = unname(tools::md5sum(path)), p_target = sort(unique(points$p_target)),
    trait_count = length(unique(points$trait_idx)), block_count = nrow(block_check),
    bird_model = bird_model
  )
}

validate_stage3_shared_provenance <- function(inputs) {
  if (is.null(inputs$mammals) || is.null(inputs$birds)) return(invisible(TRUE))
  columns <- stage3_shared_provenance_columns()
  mismatched <- columns[!vapply(columns, function(column) {
    x <- inputs$mammals$provenance[[column]]
    y <- inputs$birds$provenance[[column]]
    if (is.numeric(x) || is.numeric(y)) stage3_numeric_equal(x, y) else identical(as.character(x), as.character(y))
  }, logical(1))]
  assert(!length(mismatched),
         paste0("Mammal and bird Stage 2 inputs have mismatched scientific provenance: ",
                paste(mismatched, collapse = ", "), "."))
  invisible(TRUE)
}

resolve_stage3_input_settings <- function(config, inputs) {
  active <- Filter(Negate(is.null), inputs)
  assert(length(active) > 0L, "Stage 3 requires at least one validated taxon input.")
  provenance <- active[[1L]]$provenance
  config$expected_horizon <- as.integer(provenance$persistence_horizon_years)
  config$k0 <- as.integer(provenance$quasi_extinction_abundance)
  config
}

load_stage3_inputs <- function(config) {
  inputs <- list(
    mammals = if (config$selected_mammals) load_stage3_points(config$paths$mammals_points, "mammal", config) else NULL,
    birds = if (config$selected_birds) load_stage3_points(config$paths$birds_points, "bird", config) else NULL
  )
  validate_stage3_shared_provenance(inputs)
  inputs
}
