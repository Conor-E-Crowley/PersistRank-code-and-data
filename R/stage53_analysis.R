# Stage 5.3 initial-state persistence calculations.
#
# stage53_workflow.R owns loading this definition-only module after the shared
# analysis and patch contracts. Calls read the configured species/lookup inputs
# and return transient tables; sourcing has no side effect or expensive runtime.

stage53_selected_classes <- function(config) {
  c(
    if (config$selected_mammals) "mammalia",
    if (config$selected_birds) "aves"
  )
}

read_stage53_species <- function(config) {
  alpha_col <- paste0("alpha_", config$curve)
  beta_col <- paste0("beta_", config$curve)
  required <- c(
    "scientificName", "className", "redlistCategory", "sdm_method",
    "density", "min_pop_size", alpha_col, beta_col
  )
  species <- readr::read_csv(
    config$paths$species_csv,
    show_col_types = FALSE,
    progress = FALSE,
    col_select = dplyr::all_of(required)
  ) |>
    dplyr::mutate(
      scientificName = patch_squish(scientificName),
      class_lc = stringr::str_to_lower(stringr::str_squish(className)),
      className = dplyr::case_when(
        class_lc == "mammalia" ~ "MAMMALIA",
        class_lc == "aves" ~ "AVES",
        TRUE ~ NA_character_
      ),
      sdm_method = stringr::str_squish(as.character(sdm_method)),
      density = suppressWarnings(as.numeric(density)),
      min_pop_size = suppressWarnings(as.numeric(min_pop_size)),
      alpha = suppressWarnings(as.numeric(.data[[alpha_col]])),
      beta = suppressWarnings(as.numeric(.data[[beta_col]]))
    ) |>
    dplyr::filter(class_lc %in% stage53_selected_classes(config))

  assert(
    nrow(species) > 0L,
    paste0("Stage 5.3 selected no species for taxa='", config$taxa, "'.")
  )
  assert(
    all(!is.na(species$scientificName) & nzchar(species$scientificName)),
    "Stage 5.3 selected species contain blank scientific names."
  )
  assert(
    !anyDuplicated(species$scientificName),
    "Stage 5.3 species table contains duplicate names after whitespace normalization."
  )
  assert(
    all(!is.na(species$sdm_method) & nzchar(species$sdm_method) &
          species$sdm_method %in% c("PPM", "RangeBag")),
    "Stage 5.3 species table sdm_method must contain only PPM or RangeBag."
  )
  assert(
    all(is.finite(species$density) & species$density > 0 &
          is.finite(species$min_pop_size) & species$min_pop_size > 0),
    "Stage 5.3 selected species require positive finite density and min_pop_size values."
  )
  expected_min_pop <- quasi_extinction_abundance(config$contract) / species$density
  threshold_close <- abs(species$min_pop_size - expected_min_pop) <=
    sqrt(.Machine$double.eps) * pmax(1, abs(expected_min_pop))
  assert(
    all(threshold_close),
    paste0(
      "Stage 5.3 min_pop_size values are inconsistent with density and the requested ",
      quasi_extinction_abundance(config$contract), "-individual threshold."
    )
  )
  assert(
    all(is.finite(species$alpha) & species$alpha > 0 &
          is.finite(species$beta) & species$beta > 0),
    "Stage 5.3 selected species have invalid Gompertz parameters."
  )

  species
}

stage53_pu_persistence <- function(area_km2, density, alpha, beta, threshold_km2) {
  delta <- as.numeric(area_km2) - as.numeric(threshold_km2)
  density <- as.numeric(density); alpha <- as.numeric(alpha); beta <- as.numeric(beta)
  assert(
    all(is.finite(delta) & delta > 0 & is.finite(density) & density > 0 &
          is.finite(alpha) & alpha > 0 & is.finite(beta) & beta > 0),
    "Stage 5.3 PU persistence inputs must be finite, positive, and above the population threshold."
  )
  log_inner <- log(alpha) - beta * log(density) - beta * log(delta)
  out <- exp(-exp(log_inner))
  assert(all(is.finite(out)), "Stage 5.3 produced nonfinite PU persistence values.")
  pmin(pmax(out, 0), 1)
}

stage53_species_persistence <- function(x) {
  x <- as.numeric(x)
  assert(length(x) > 0L && all(is.finite(x) & x >= 0 & x <= 1),
         "Species persistence requires finite PU probabilities in [0, 1].")
  if (any(x >= 1)) return(1)
  as.numeric(-expm1(sum(log1p(-x))))
}

calculate_stage53_persistence <- function(species_selected, lookup) {
  validate_patch_lookup_object(lookup, "Stage 5.3 patch lookup")
  represented_names <- intersect(
    species_selected$scientificName,
    unique(patch_squish(lookup$scientificName))
  )
  absent <- species_selected |>
    dplyr::filter(!scientificName %in% represented_names) |>
    dplyr::select(scientificName, className, redlistCategory)
  assert(
    length(represented_names) > 0L,
    "No selected species are represented in the Stage 5 patch lookup."
  )

  pu <- lookup |>
    dplyr::filter(scientificName %in% represented_names) |>
    dplyr::group_by(scientificName, pu_id) |>
    dplyr::summarise(
      pu_area_km2 = sum(patch_area_km2),
      n_patches = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::left_join(
      species_selected |>
        dplyr::select(scientificName, className, density, min_pop_size, alpha, beta),
      by = "scientificName"
    ) |>
    dplyr::mutate(
      P_pu = stage53_pu_persistence(pu_area_km2, density, alpha, beta, min_pop_size)
    )
  species <- pu |>
    dplyr::group_by(scientificName, className) |>
    dplyr::summarise(
      n_pu = dplyr::n(),
      total_pu_area_km2 = sum(pu_area_km2),
      species_persistence = stage53_species_persistence(P_pu),
      .groups = "drop"
    ) |>
    dplyr::left_join(
      species_selected |>
        dplyr::select(scientificName, sdm_method),
      by = "scientificName"
    )

  list(
    represented_names = represented_names,
    absent_species = absent,
    pu = pu,
    species = species
  )
}
