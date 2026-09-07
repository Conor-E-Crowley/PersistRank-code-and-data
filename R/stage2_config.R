# Pure Stage 2 configuration and probability-grid construction.
#
# Loaded first by stage2_workflow.R. Parameter validation and path derivation
# perform no filesystem access, package loading, compilation, or logging;
# operational preflight belongs to stage2_inputs.R.


required_stage2_packages <- function() c("data.table", "Rcpp", "dplyr")

stage2_probabilities <- function(x, label, minimum_length = 1L) {
  values <- suppressWarnings(as.numeric(unlist(
    x, recursive = TRUE, use.names = FALSE
  )))
  assert(
    length(values) >= minimum_length &&
      all(is.finite(values) & values > 0 & values < 1),
    paste0(
      label, " must contain at least ", minimum_length,
      " finite probabilities in (0, 1)."
    )
  )
  assert(!anyDuplicated(values), paste0(label, " must not contain duplicates."))
  assert(all(diff(values) > 0), paste0(label, " must be strictly increasing."))
  values
}

normalize_stage2_curves <- function(x, label = "params$curves") {
  curves <- unlist(x, recursive = TRUE, use.names = FALSE)
  assert(is.character(curves) && length(curves), paste0(label, " must be a character vector."))
  curves <- trimws(curves)
  assert(all(!is.na(curves) & nzchar(curves)), paste0(label, " contains blank values."))
  assert(!anyDuplicated(curves), paste0(label, " contains duplicates."))
  invalid <- setdiff(curves, persistence_curves())
  assert(
    !length(invalid),
    paste0(label, " contains unsupported curve(s): ", paste(invalid, collapse = ", "), ".")
  )
  assert(main_persistence_curve() %in% curves, paste0(label, " must include q50."))
  persistence_curves()[persistence_curves() %in% curves]
}

persist_numeric_close <- function(x, expected, tolerance = 1e-10) {
  x <- suppressWarnings(as.numeric(x))
  expected <- as.numeric(expected)
  is.finite(x) & is.finite(expected) &
    abs(x - expected) <= tolerance * pmax(1, abs(expected))
}

persist_values_in_grid <- function(x, grid, tolerance = 1e-10) {
  values <- suppressWarnings(as.numeric(x))
  vapply(
    values,
    function(value) any(persist_numeric_close(value, grid, tolerance)),
    logical(1)
  )
}

format_invalid_persist_grid_values <- function(x, grid, max_values = 5L) {
  values <- suppressWarnings(as.numeric(x))
  bad <- unique(values[!persist_values_in_grid(values, grid)])
  bad <- utils::head(bad, max_values)
  paste0(
    "Invalid p_target values: ",
    paste(format(bad, digits = 17), collapse = ", "),
    "."
  )
}

build_stage2_probability_grid <- function(
  minimum,
  maximum,
  step,
  additional
) {
  minimum <- validate_scalar_number(
    minimum, "params$persistence_grid_min",
    minimum = 0, minimum_open = TRUE
  )
  maximum <- validate_scalar_number(
    maximum, "params$persistence_grid_max",
    minimum = 0, minimum_open = TRUE
  )
  step <- validate_scalar_number(
    step, "params$persistence_grid_step",
    minimum = 0, minimum_open = TRUE
  )
  assert(maximum < 1 && maximum > minimum, "The regular probability range must lie inside (0, 1).")

  raw_steps <- (maximum - minimum) / step
  n_steps <- round(raw_steps)
  tolerance <- 1e-10 * max(1, abs(raw_steps))
  assert(
    abs(raw_steps - n_steps) <= tolerance,
    "params$persistence_grid_step must divide the configured interval exactly."
  )
  regular <- minimum + seq.int(0L, n_steps) * step
  regular[[length(regular)]] <- maximum

  additional <- if (!length(unlist(additional, recursive = TRUE))) {
    numeric()
  } else {
    stage2_probabilities(
      as.numeric(unlist(additional, recursive = TRUE, use.names = FALSE)),
      "params$additional_persistence_probabilities"
    )
  }
  assert(
    !any(vapply(
      additional,
      function(value) any(abs(regular - value) <= 1e-12),
      logical(1)
    )),
    "Regular and additional persistence probabilities must not overlap."
  )
  final <- sort(c(regular, additional))
  stage2_probabilities(final, "final persistence grid", minimum_length = 3L)
}

stage2_signature_number <- function(x) {
  format(as.numeric(x), digits = 17, scientific = FALSE, trim = TRUE)
}

# Resolve NULL before constructing the grid so every artifact records the
# effective integer, keeping automatic and equivalent explicit settings exact.
resolve_stage2_k_search_start <- function(requested, quasi_extinction_abundance) {
  threshold <- validate_scalar_integer(
    quasi_extinction_abundance, "quasi_extinction_abundance"
  )
  if (is.null(requested) || !length(requested)) {
    requested <- 2^ceiling(log2(as.numeric(threshold) + 1))
    assert(
      is.finite(requested) && requested <= .Machine$integer.max,
      paste0(
        "The automatic k_search_start for quasi_extinction_abundance=",
        threshold, " exceeds the supported integer range. Supply an explicit ",
        "k_search_start within that range or use a smaller threshold."
      )
    )
  }
  start_k <- validate_scalar_integer(requested, "params$k_search_start")
  assert(
    start_k > threshold,
    "params$k_search_start must exceed the quasi-extinction abundance."
  )
  start_k
}

stage2_grid_signature <- function(grid) {
  paste0(
    "curves=", paste(
      names(grid$q_levels), stage2_signature_number(grid$q_levels),
      sep = ":", collapse = ","
    ),
    "|anchors=", paste(stage2_signature_number(grid$p_anchor), collapse = ","),
    "|targets=", paste(stage2_signature_number(grid$p_grid), collapse = ","),
    "|start_k=", stage2_signature_number(grid$k_search$start_k),
    "|rel_tol=", stage2_signature_number(grid$k_search$rel_tol),
    "|rounding=", grid$k_search$rounding,
    "|max_k=", stage2_signature_number(grid$k_search$max_k)
  )
}

stage2_paths <- function(project = project_paths(), source_project = project) {
  list(
    data_dir = project$data,
    clean_dir = project$clean,
    results_dir = project$results,
    checkpoint_dir = project$checkpoints,
    runtime_log = file.path(project$data, "Logs", "stage2.log"),
    mammal_mass_grid = file.path(source_project$clean, "mammal_mass_grid.csv"),
    bird_generation_length_grid =
      file.path(source_project$clean, "bird_generation_length_grid.csv"),
    mammal_growth_posterior =
      file.path(source_project$clean, "mammal_growth_posterior.csv"),
    mammal_environmental_variation_posterior =
      file.path(source_project$clean, "mammal_environmental_variation_posterior.csv"),
    bird_growth_posterior =
      file.path(source_project$clean, "bird_growth_posterior.csv"),
    bird_environmental_variation_posterior =
      file.path(source_project$clean, "bird_environmental_variation_posterior.csv"),
    cpp_file = file.path(project$root, "src", "simulate_persist_probs_cpp.cpp"),
    rcpp_cache = project$rcpp_cache,
    mammals_points = file.path(project$results, "persistence_points_mammals.csv"),
    birds_points = file.path(project$results, "persistence_points_birds.csv"),
    mammals_partial =
      file.path(project$results, "persistence_points_mammals.partial.csv"),
    birds_partial =
      file.path(project$results, "persistence_points_birds.partial.csv"),
    mammals_checkpoint = file.path(project$checkpoints, "stage2_mammals.rds"),
    birds_checkpoint = file.path(project$checkpoints, "stage2_birds.rds")
  )
}

stage2_config <- function(params, paths = project_paths(), source_paths = paths) {
  mode <- validate_scalar_choice(
    params$mode, c("inspect", "resume", "restart"), "params$mode"
  )
  selection <- validate_taxa_selector(params$taxa, "params$taxa")
  uncertainty <- validate_scalar_choice(
    params$demographic_uncertainty,
    c("coefficients_only", "posterior_predictive"),
    "params$demographic_uncertainty"
  )
  selected_curves <- normalize_stage2_curves(params$curves)

  anchors <- stage2_probabilities(
    params$anchor_probabilities,
    "params$anchor_probabilities",
    minimum_length = 3L
  )
  assert(
    any(anchors == 0.5),
    "params$anchor_probabilities must include 0.50 exactly."
  )
  p_grid <- build_stage2_probability_grid(
    params$persistence_grid_min,
    params$persistence_grid_max,
    params$persistence_grid_step,
    params$additional_persistence_probabilities
  )

  sim <- list(
    years = validate_scalar_integer(
      params$persistence_horizon_years,
      "params$persistence_horizon_years"
    ),
    quasi_extinction_abundance = validate_scalar_integer(
      params$quasi_extinction_abundance,
      "params$quasi_extinction_abundance"
    ),
    cap_factor = validate_scalar_number(
      params$population_cap_factor, "params$population_cap_factor", 1
    ),
    r_buffer = validate_scalar_number(
      params$growth_rate_buffer,
      "params$growth_rate_buffer",
      minimum = 0,
      maximum = 1,
      minimum_open = TRUE
    ),
    n_draws = validate_scalar_integer(
      params$n_posterior_draws, "params$n_posterior_draws"
    ),
    reps = validate_scalar_integer(
      params$replicates_per_draw, "params$replicates_per_draw"
    ),
    chunk_size = validate_scalar_integer(
      params$posterior_chunk_size, "params$posterior_chunk_size"
    ),
    base_seed = validate_scalar_integer(
      params$base_seed, "params$base_seed", 0L
    ),
    demographic_uncertainty = uncertainty
  )
  assert(
    sim$chunk_size <= sim$n_draws,
    "params$posterior_chunk_size cannot exceed params$n_posterior_draws."
  )

  rounding <- validate_scalar_choice(
    params$k_rounding, c("round", "ceiling"), "params$k_rounding"
  )
  k_search <- list(
    start_k = resolve_stage2_k_search_start(
      params$k_search_start, sim$quasi_extinction_abundance
    ),
    rel_tol = validate_scalar_number(
      params$k_search_relative_tolerance,
      "params$k_search_relative_tolerance",
      minimum = 0,
      maximum = 1,
      maximum_open = TRUE,
      minimum_open = TRUE
    ),
    max_k = validate_scalar_number(
      params$k_search_max,
      "params$k_search_max",
      minimum = 0,
      minimum_open = TRUE
    ),
    rounding = rounding,
    round_k = identical(rounding, "round")
  )
  assert(
    k_search$max_k > k_search$start_k,
    "params$k_search_max must exceed params$k_search_start."
  )

  threads <- params$threads
  if (is.null(threads) || !length(threads)) {
    threads <- NULL
  } else {
    threads <- validate_scalar_integer(threads, "params$threads")
  }

  grid <- list(
    q_levels = persistence_quantiles()[selected_curves],
    p_anchor = anchors,
    p_grid = p_grid,
    k_search = k_search
  )
  grid$signature <- stage2_grid_signature(grid)

  list(
    mode = mode,
    active = mode != "inspect",
    taxa = selection$value,
    selected_mammals = selection$selected_mammals,
    selected_birds = selection$selected_birds,
    verbose = validate_scalar_logical(params$verbose, "params$verbose"),
    curves = selected_curves,
    threads = threads,
    sim = sim,
    grid = grid,
    paths = stage2_paths(paths, source_paths)
  )
}

describe_stage2_config <- function(config) {
  paste0(
    "Stage 2 | mode=", config$mode,
    " | taxa=", config$taxa,
    " | uncertainty=", config$sim$demographic_uncertainty,
    " | curves=", paste(config$curves, collapse = ","),
    " | draws=", config$sim$n_draws,
    " | replicates=", config$sim$reps,
    " | years=", config$sim$years,
    " | quasi_extinction_abundance=", config$sim$quasi_extinction_abundance,
    " | k_search_start=", config$grid$k_search$start_k,
    " | targets=", length(config$grid$p_grid),
    " | threads=", config$threads %||% "runtime-default"
  )
}
