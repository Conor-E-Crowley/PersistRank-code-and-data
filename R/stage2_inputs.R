# Stage 1 handoff and deterministic input preparation for Stage 2.
#
# Loaded after pure configuration and the deferred runtime contract. Sourcing
# is definition-only; explicit preflight/read calls validate artifacts and
# prepare directories. Native compilation remains deferred to stage2_runtime.R.

check_stage2_preflight <- function(config) {
  paths <- config$paths
  selected <- list()
  if (config$selected_mammals) {
    selected$mammals <- c(
      paths$mammal_mass_grid,
      paths$mammal_growth_posterior,
      paths$mammal_environmental_variation_posterior
    )
  }
  if (config$selected_birds) {
    selected$birds <- c(
      paths$bird_generation_length_grid,
      paths$bird_growth_posterior,
      paths$bird_environmental_variation_posterior
    )
  }
  invisible(lapply(unlist(selected, use.names = FALSE), need_file))

  if (config$active) {
    need_file(paths$cpp_file, "Stage 2 C++ simulator")
    ensure_writable_dir(paths$results_dir, "Stage 2 results directory")
    ensure_writable_dir(paths$checkpoint_dir, "Stage 2 checkpoint directory")
    ensure_writable_dir(paths$rcpp_cache, "Stage 2 Rcpp cache directory")
  }
  invisible(TRUE)
}


validate_stage2_trait_grid <- function(values, label, minimum_points = 5L) {
  values <- as.numeric(values)
  assert(
    length(values) >= minimum_points,
    paste0(label, " must contain at least ", minimum_points, " trait values.")
  )
  assert(all(is.finite(values) & values > 0), paste0(label, " must contain positive finite values."))
  assert(!anyDuplicated(values), paste0(label, " must not contain duplicate values."))
  assert(all(diff(values) > 0), paste0(label, " must be strictly increasing."))
  values
}

# ---- Stage 1 handoff -----------------------------------------------------------

read_stage2_numeric_table <- function(path, cols, label) {
  need_file(path, label)
  out <- data.table::fread(path)
  assert(
    identical(names(out), cols),
    paste0(
      label, " must contain exactly these columns in order: ",
      paste(cols, collapse = ", "), "."
    )
  )

  bad_cols <- cols[
    !vapply(out[, cols, with = FALSE], function(x) is.numeric(x) && all(is.finite(x)), logical(1))
  ]
  assert(
    length(bad_cols) == 0L,
    paste0(label, " contains non-numeric or non-finite values in column(s): ",
           paste(bad_cols, collapse = ", "))
  )

  out[, cols, with = FALSE]
}

load_mammal_demographic_inputs <- function(paths) {
  mass <- validate_stage2_trait_grid(
    read_stage2_numeric_table(paths$mammal_mass_grid, "Mass_g", "mammal_mass_grid.csv")[["Mass_g"]],
    "mammal_mass_grid.csv"
  )

  list(
    mass = mass,
    post_rm = read_stage2_numeric_table(
      paths$mammal_growth_posterior,
      c("alpha", "beta_logM", "residual_sd"),
      "mammal_growth_posterior.csv"
    ),
    post_sigma = read_stage2_numeric_table(
      paths$mammal_environmental_variation_posterior,
      c("alpha", "beta_logM", "residual_sd"),
      "mammal_environmental_variation_posterior.csv"
    ),
    provenance = list(
      posterior_rm_md5 = unname(tools::md5sum(paths$mammal_growth_posterior)),
      posterior_sigma_md5 = unname(tools::md5sum(
        paths$mammal_environmental_variation_posterior
      ))
    )
  )
}

load_bird_demographic_inputs <- function(paths) {
  genlength <- validate_stage2_trait_grid(
    read_stage2_numeric_table(
      paths$bird_generation_length_grid,
      "GenLength",
      "bird_generation_length_grid.csv"
    )[["GenLength"]],
    "bird_generation_length_grid.csv"
  )
  posterior_header <- names(data.table::fread(
    paths$bird_environmental_variation_posterior,
    nrows = 0L
  ))
  bird_model <- bird_model_spec_from_posterior(posterior_header)
  sigma_cols <- c(
    "alpha", "beta_logGenLength",
    bird_model$coefficient_names, "residual_sd"
  )

  list(
    genlength = genlength,
    bird_model = bird_model,
    diet_groups = bird_model$model_groups,
    post_rm = read_stage2_numeric_table(
      paths$bird_growth_posterior,
      c("alpha", "beta_logGenLength", "residual_sd"),
      "bird_growth_posterior.csv"
    ),
    post_sigma = read_stage2_numeric_table(
      paths$bird_environmental_variation_posterior,
      sigma_cols,
      "bird_environmental_variation_posterior.csv"
    ),
    provenance = list(
      posterior_rm_md5 = unname(tools::md5sum(paths$bird_growth_posterior)),
      posterior_sigma_md5 = unname(tools::md5sum(
        paths$bird_environmental_variation_posterior
      ))
    )
  )
}

mammal_persistence_trait_table <- function(mass) {
  data.frame(
    mass_idx = seq_along(mass),
    Mass_g = as.numeric(mass),
    stringsAsFactors = FALSE
  )
}

bird_persistence_trait_table <- function(genlength, bird_model) {
  diet_groups <- bird_model$model_groups
  validate_bird_model_group(
    diet_groups, bird_model, "bird persistence model group"
  )
  grid <- expand.grid(
    genlength_idx = seq_along(genlength),
    bird_sigma_model_group = diet_groups,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  grid$GenLength <- as.numeric(genlength[grid$genlength_idx])
  grid[, c("genlength_idx", "GenLength", "bird_sigma_model_group")]
}

sample_posterior_indices <- function(n_rm, n_sigma, n_draws, base_seed, rm_offset, sigma_offset) {
  assert(
    n_draws <= n_rm && n_draws <= n_sigma,
    paste0(
      "n_draws must be no larger than available posterior rows ",
      "(rm=", n_rm, ", sigma=", n_sigma, ")."
    )
  )

  rm_seed <- as.integer(base_seed + rm_offset)
  sigma_seed <- as.integer(base_seed + sigma_offset)

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(rm_seed)
  idx_rm <- sample.int(n_rm, n_draws, replace = FALSE)
  set.seed(sigma_seed)
  idx_sigma <- sample.int(n_sigma, n_draws, replace = FALSE)
  list(rm = idx_rm, sigma = idx_sigma, rm_seed = rm_seed, sigma_seed = sigma_seed)
}

# Posterior row sampling and residual deviates use independent, documented
# seed offsets. Residual vectors are constructed once and reused across traits
# and bird branches; the C++ population trajectories use a separate CRN stream.
deterministic_standard_normals <- function(n, seed) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed))
  stats::rnorm(n)
}

# ---- Runtime inputs and metadata ----------------------------------------------

stage2_residual_streams <- function(sim, taxon) {
  if (identical(sim$demographic_uncertainty, "coefficients_only")) {
    return(list(growth = NULL, environmental_variation = NULL, seeds = integer()))
  }
  offsets <- switch(
    taxon,
    mammals = c(growth = 11001L, environmental_variation = 12001L),
    birds = c(growth = 13001L, environmental_variation = 14001L),
    stop("Unknown Stage 2 taxon.", call. = FALSE)
  )
  seeds <- as.integer(sim$base_seed + offsets)
  names(seeds) <- names(offsets)
  list(
    growth = deterministic_standard_normals(sim$n_draws, seeds[["growth"]]),
    environmental_variation = deterministic_standard_normals(
      sim$n_draws, seeds[["environmental_variation"]]
    ),
    seeds = seeds
  )
}

stage2_taxon_metadata <- function(inputs, config) {
  residual_seed <- function(name) {
    seeds <- inputs$residual$seeds
    if (length(seeds) && name %in% names(seeds)) {
      as.integer(seeds[[name]])
    } else {
      NA_integer_
    }
  }
  list(
    posterior_rm_seed = inputs$idx$rm_seed,
    posterior_sigma_seed = inputs$idx$sigma_seed,
    demographic_uncertainty = config$sim$demographic_uncertainty,
    residual_growth_seed = residual_seed("growth"),
    residual_environmental_variation_seed =
      residual_seed("environmental_variation"),
    bird_model_signature = if (is.null(inputs$bird_model)) {
      "not_applicable"
    } else {
      bird_model_signature(inputs$bird_model)
    },
    grid_signature = config$grid$signature,
    posterior_sampling = inputs$provenance$posterior_sampling,
    posterior_rm_md5 = inputs$provenance$posterior_rm_md5,
    posterior_sigma_md5 = inputs$provenance$posterior_sigma_md5,
    simulator_contract = stage2_simulator_contract()
  )
}

load_persist_inputs <- function(config) {
  paths <- config$paths
  sim <- config$sim
  out <- list(mammals = NULL, birds = NULL)

  if (config$selected_mammals) {
    m <- load_mammal_demographic_inputs(paths)
    m$idx <- sample_posterior_indices(
      n_rm = nrow(m$post_rm),
      n_sigma = nrow(m$post_sigma),
      n_draws = sim$n_draws,
      base_seed = sim$base_seed,
      rm_offset = 7001L,
      sigma_offset = 8001L
    )
    m$provenance$posterior_sampling <- "independent_without_replacement"
    m$residual <- stage2_residual_streams(sim, "mammals")
    m$bird_model <- NULL
    out$mammals <- m
  }

  if (config$selected_birds) {
    b <- load_bird_demographic_inputs(paths)
    b$idx <- sample_posterior_indices(
      n_rm = nrow(b$post_rm),
      n_sigma = nrow(b$post_sigma),
      n_draws = sim$n_draws,
      base_seed = sim$base_seed,
      rm_offset = 9001L,
      sigma_offset = 10001L
    )
    b$provenance$posterior_sampling <- "independent_without_replacement"
    b$residual <- stage2_residual_streams(sim, "birds")
    out$birds <- b
  }

  out
}
