# Stage 1 calibration-input assembly and validation.
#
# Loaded by stage1_workflow.R after deterministic bird matching. It owns the
# published exclusions, joins, and final calibration tables. Sourcing defines
# constants and functions only; file reads and expensive work occur on call.

MAMMAL_MARINE_RMAX_SPECIES <- c(
  "Phoca fasciata", "Phoca groenlandica", "Balaena mysticetus", "Eubalaena glacialis",
  "Balaenoptera acutorostrata", "Balaenoptera borealis", "Balaenoptera musculus",
  "Balaenoptera physalus", "Megaptera novaeangliae", "Globicephala melas",
  "Lagenorhynchus acutus", "Orcinus orca", "Stenella attenuata", "Stenella coeruleoalba",
  "Stenella longirostris", "Tursiops truncatus", "Eschrichtius robustus",
  "Delphinapterus leucas", "Monodon monoceros", "Phocoena phocoena",
  "Phocoenoides dalli", "Physeter catodon", "Pontoporia blainvillei", "Berardius bairdii",
  "Trichechus manatus"
)

MAMMAL_MARINE_SIGMA_SPECIES <- c(
  "Phoca groenlandica", "Stenella attenuata", "Eschrichtius robustus", "Trichechus manatus"
)

# ---- Source-table readers and basic contracts --------------------------------

read_calibration_inputs <- function(paths) {
  raw_mammal_rmax <- readr::read_tsv(paths$raw$mammal_rmax, show_col_types = FALSE)
  raw_sigma <- readr::read_csv(
    paths$raw$sigma,
    show_col_types = FALSE,
    col_select = dplyr::any_of(c("Class", "Genus", "Species", "Vr", "Mass"))
  )
  bird_elton_traits <- read_bird_elton_traits(paths$raw$bird_traits)
  bird_generation_lengths <- read_bird_generation_lengths(paths$raw$bird_generation_lengths)
  raw_iucn_bird_synonyms <- readr::read_csv(paths$raw$iucn_bird_synonyms, show_col_types = FALSE)
  raw_input_synonyms <- readr::read_csv(paths$raw$input_synonyms, show_col_types = FALSE)
  bird_growth <- readr::read_csv(
    paths$raw$bird_growth,
    show_col_types = FALSE,
    progress = FALSE
  )

  validate_raw_inputs(
    raw_mammal_rmax,
    raw_sigma,
    bird_elton_traits,
    bird_generation_lengths,
    raw_iucn_bird_synonyms,
    raw_input_synonyms,
    bird_growth
  )

  list(
    mammal_rmax = raw_mammal_rmax,
    sigma = raw_sigma,
    bird_elton_traits = bird_elton_traits,
    bird_generation_lengths = bird_generation_lengths,
    bird_growth = validate_niel_rmax_table(bird_growth),
    iucn_bird_synonyms = raw_iucn_bird_synonyms,
    input_synonyms = raw_input_synonyms
  )
}

validate_raw_inputs <- function(raw_mammal_rmax, raw_sigma, bird_elton_traits,
                                bird_generation_lengths,
                                raw_iucn_bird_synonyms, raw_input_synonyms,
                                bird_growth) {
  need_cols(raw_mammal_rmax, c("log10(M)", "log10(rm)", "Order", "Species"), "mammal_rmax.txt")
  need_cols(raw_sigma, c("Genus", "Species", "Class", "Mass", "Vr"), "sigma.csv")
  need_cols(bird_elton_traits, c("elton_scientific", "elton_key", "Diet_5Cat"), "EltonTraits bird data")
  need_cols(bird_generation_lengths, c("bird_scientific", "GenLength", "bird_key"), "Bird et al. generation lengths")
  need_cols(raw_iucn_bird_synonyms, c("scientificName", "genusName", "speciesName"), "iucn_bird_synonyms.csv")
  need_cols(raw_input_synonyms, c("scientificName", "synonym_scientificName"), "input_synonyms.csv")
  need_cols(bird_growth, c("Species", "lambda"), "bird_growth_niel_lebreton.csv")

  invisible(TRUE)
}

build_bird_rmax_join <- function(niel_bird_rmax, bird_generation_lengths,
                                 raw_iucn_bird_synonyms, raw_input_synonyms,
                                 verbose = FALSE) {
  match_result <- match_bird_generation_lengths(
    niel_bird_rmax$Species,
    bird_generation_lengths,
    raw_iucn_bird_synonyms,
    raw_input_synonyms,
    query_label = "Niel & Lebreton bird growth"
  )
  report_bird_generation_matches(
    match_result,
    "Niel & Lebreton bird growth records",
    include_all_matches = TRUE,
    verbose = verbose
  )

  if (nrow(match_result$unmatched) > 0) {
    stop(
      "Niel & Lebreton bird growth records must all match Bird et al. generation lengths before fitting.",
      call. = FALSE
    )
  }

  niel_bird_rmax |>
    dplyr::left_join(
      match_result$matches |>
        dplyr::transmute(
          Species = query_scientific,
          bird_scientific,
          GenLength,
          match_source
        ),
      by = "Species"
    )
}

build_bird_sigma_join <- function(bird_sigma, bird_generation_lengths,
                                  raw_iucn_bird_synonyms, raw_input_synonyms,
                                  verbose = FALSE) {
  match_result <- match_bird_generation_lengths(
    bird_sigma$sigma_scientific,
    bird_generation_lengths,
    raw_iucn_bird_synonyms,
    raw_input_synonyms,
    query_label = "Brook bird sigma"
  )
  report_bird_generation_matches(
    match_result,
    "Brook bird sigma records",
    verbose = verbose
  )

  if (nrow(match_result$unmatched) > 0L || nrow(match_result$ambiguous) > 0L ||
      nrow(match_result$manual_excluded) > 0L) {
    stop(
      "Brook bird sigma records must all resolve to one Bird et al. generation length before fitting.",
      call. = FALSE
    )
  }

  bird_sigma |>
    dplyr::inner_join(
      match_result$matches |>
        dplyr::transmute(
          sigma_scientific = query_scientific,
          bird_scientific,
          GenLength,
          match_source
        ),
      by = "sigma_scientific"
    )
}

prepare_calibration_data <- function(
  inputs,
  bird_model,
  niel_bird_rmax = inputs$bird_growth,
  verbose = FALSE
) {
  niel_bird_rmax <- validate_niel_rmax_table(niel_bird_rmax)

  raw_sigma_parsed <- inputs$sigma |>
    dplyr::mutate(Mass = as_num(Mass), Vr = as_num(Vr))

  invalid_sigma <- !is.finite(raw_sigma_parsed$Mass) | raw_sigma_parsed$Mass <= 0 |
    !is.finite(raw_sigma_parsed$Vr) | raw_sigma_parsed$Vr <= 0
  raw_sigma <- raw_sigma_parsed[!invalid_sigma, , drop = FALSE]

  mammal_rmax_source <- inputs$mammal_rmax |>
    dplyr::rename(log10_M = `log10(M)`, log10_rm = `log10(rm)`) |>
    dplyr::mutate(
      exclusion_reason = dplyr::case_when(
        Order == "Chiroptera" ~ "bat",
        Species %in% MAMMAL_MARINE_RMAX_SPECIES ~ "marine_or_aquatic",
        !is.finite(log10_M) | !is.finite(log10_rm) ~ "invalid_numeric_value",
        TRUE ~ "retained"
      )
    )

  mammal_rmax <- mammal_rmax_source |>
    dplyr::filter(
      exclusion_reason == "retained"
    ) |>
    dplyr::transmute(Mass_g = 10^log10_M, rm = 10^log10_rm)

  sigma_std <- raw_sigma |>
    dplyr::transmute(
      Class = Class,
      Species = paste(Genus, Species),
      Mass_g = Mass,
      sigma = sqrt(Vr)
    )

  mammal_sigma <- sigma_std |>
    dplyr::filter(Class == "Mammalia", !(Species %in% MAMMAL_MARINE_SIGMA_SPECIES)) |>
    dplyr::transmute(Mass_g, sigma)

  bird_sigma <- sigma_std |>
    dplyr::filter(Class == "Aves") |>
    dplyr::transmute(sigma_scientific = Species, sigma)

  bird_rmax <- build_bird_rmax_join(
    niel_bird_rmax,
    inputs$bird_generation_lengths,
    inputs$iucn_bird_synonyms,
    inputs$input_synonyms,
    verbose = verbose
  ) |>
    dplyr::mutate(rm = log(lambda)) |>
    dplyr::transmute(
      Species,
      bird_scientific,
      GenLength,
      rm,
      match_source
    )

  bird_sigma_data <- build_bird_sigma_join(
    bird_sigma,
    inputs$bird_generation_lengths,
    inputs$iucn_bird_synonyms,
    inputs$input_synonyms,
    verbose = verbose
  ) |>
    dplyr::transmute(
      sigma_scientific,
      bird_scientific,
      GenLength,
      sigma,
      match_source
    ) |>
    attach_bird_sigma_model_groups(
      inputs$bird_elton_traits,
      inputs$iucn_bird_synonyms,
      inputs$input_synonyms,
      bird_model = bird_model,
      verbose = verbose
    )

  validate_calibration_data(
    mammal_rmax, mammal_sigma, bird_rmax, bird_sigma_data, bird_model
  )

  summary <- tibble::tribble(
    ~dataset, ~raw_records, ~excluded_records, ~retained_records,
    "mammal_rm", nrow(mammal_rmax_source), sum(mammal_rmax_source$exclusion_reason != "retained"), nrow(mammal_rmax),
    "mammal_sigma", sum(raw_sigma_parsed$Class == "Mammalia", na.rm = TRUE),
      sum(raw_sigma_parsed$Class == "Mammalia" & invalid_sigma, na.rm = TRUE) +
        sum(raw_sigma$Class == "Mammalia" & paste(raw_sigma$Genus, raw_sigma$Species) %in% MAMMAL_MARINE_SIGMA_SPECIES),
      nrow(mammal_sigma),
    "bird_rm", nrow(niel_bird_rmax), nrow(niel_bird_rmax) - nrow(bird_rmax), nrow(bird_rmax),
    "bird_sigma", sum(raw_sigma_parsed$Class == "Aves", na.rm = TRUE),
      sum(raw_sigma_parsed$Class == "Aves" & invalid_sigma, na.rm = TRUE) +
        sum(raw_sigma$Class == "Aves", na.rm = TRUE) - nrow(bird_sigma_data),
      nrow(bird_sigma_data)
  )

  exclusion_counts <- mammal_rmax_source |>
    dplyr::count(exclusion_reason, name = "n") |>
    dplyr::arrange(exclusion_reason)

  list(
    mammal_rmax = mammal_rmax,
    mammal_sigma = mammal_sigma,
    bird_rmax = bird_rmax,
    bird_sigma_data = bird_sigma_data,
    summary = summary,
    mammal_rmax_exclusions = exclusion_counts
  )
}

# ---- Final scientific invariants ----------------------------------------------

assert_positive_finite_columns <- function(data, columns, label) {
  bad <- columns[
    !vapply(
      data[columns],
      function(x) is.numeric(x) && all(is.finite(x) & x > 0),
      logical(1L)
    )
  ]
  assert(
    length(bad) == 0L,
    paste0(label, " must contain positive finite values in: ", paste(bad, collapse = ", "), ".")
  )
  invisible(TRUE)
}

validate_calibration_data <- function(
  mammal_rmax,
  mammal_sigma,
  bird_rmax,
  bird_sigma_data,
  bird_model
) {
  for (item in list(mammal_rmax, mammal_sigma, bird_rmax, bird_sigma_data)) {
    assert(nrow(item) > 1L, "Each calibration dataset requires at least two records.")
  }
  assert_positive_finite_columns(mammal_rmax, c("Mass_g", "rm"), "Mammal rm calibration data")
  assert_positive_finite_columns(mammal_sigma, c("Mass_g", "sigma"), "Mammal sigma calibration data")
  assert_positive_finite_columns(bird_rmax, c("GenLength", "rm"), "Bird rm calibration data")
  assert_positive_finite_columns(
    bird_sigma_data,
    c("GenLength", "sigma"),
    "Bird sigma calibration data"
  )

  validate_bird_sigma_diet_matches(bird_sigma_data, bird_model)

  # Validate selected categories against the actual calibration data here,
  # rather than baking manuscript-specific sample counts into reusable code.
  bird_demographic_model_spec(
    bird_model$separate_intercepts,
    observed_diet5 = bird_sigma_data$diet5_group,
    require_observed = TRUE
  )

  invisible(TRUE)
}

