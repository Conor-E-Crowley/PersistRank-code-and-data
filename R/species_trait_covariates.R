# Stage 4 ecological trait, density, and dispersal covariates.
#
# Loaded by stage4_workflow.R before Gompertz enrichment. It owns published
# coefficients and deterministic covariate calculations. Sourcing defines small
# constants and functions only; no package is attached and no artifact is read.

diet_cols <- list(
  inv = "Diet-Inv",
  plant = c("Diet-Fruit", "Diet-Nect", "Diet-Seed", "Diet-PlantO"),
  vert = c("Diet-Vend", "Diet-Vect", "Diet-Vfish", "Diet-Vunk", "Diet-Scav")
)

# Published allometric coefficients. Body mass is supplied in grams for
# density models and kilograms for home-range models; areas are km^2.
mammal_density_coefficients <- list(
  intercept = 2.4814753048,
  mass = -0.3296510182,
  diet = c(herbivore = -0.1003301887, omnivore = -0.0552981064, carnivore = 0)
)
bird_density_coefficients <- list(
  intercept = 1.30947404816817,
  mass = -0.341277579477041,
  diet = c(
    FruiNect = 0.42760614375972505,
    Invertebrate = 0.181241393478166,
    Omnivore = 0.28471583669890205,
    PlantSeed = 0.34312864063468,
    VertFishScav = 0
  )
)

diet_sum <- function(df, cols) {
  present <- intersect(cols, names(df))
  if (!length(present)) return(rep(NA_real_, nrow(df)))
  mat <- dplyr::mutate(df[present], dplyr::across(dplyr::everything(), as_number))
  rowSums(mat, na.rm = TRUE)
}

diet_aggregates <- function(df) {
  inv <- if (diet_cols$inv %in% names(df)) as_number(df[[diet_cols$inv]]) else rep(NA_real_, nrow(df))
  fruit_nectar <- diet_sum(df, c("Diet-Fruit", "Diet-Nect"))
  plant_seed <- diet_sum(df, c("Diet-Seed", "Diet-PlantO"))
  plant <- fruit_nectar + plant_seed
  vert <- diet_sum(df, diet_cols$vert)

  tibble::tibble(
    Diet_Inv = inv,
    Diet_FruiNect = fruit_nectar,
    Diet_PlantSeed = plant_seed,
    Diet_AllPlants = plant,
    Diet_VertFishScav = vert,
    Diet_AllAnimal = dplyr::coalesce(inv, 0) + dplyr::coalesce(vert, 0)
  )
}

classify_mammal_diet <- function(all_animal, all_plants, mass_g) {
  dplyr::case_when(
    !is.finite(mass_g) ~ NA_character_,
    all_animal >= 80 ~ "carnivore",
    all_plants >= 80 ~ "herbivore",
    TRUE ~ "omnivore"
  )
}

get_bird_diet5 <- function(df) {
  preferred <- c("Diet.5Cat", "Diet-5Cat", "Diet_5Cat", "diet.5cat", "diet-5cat", "diet_5cat")
  hit <- intersect(preferred, names(df))
  if (length(hit)) {
    x <- stringr::str_squish(as.character(df[[hit[1]]]))
    x[x == ""] <- NA_character_
    return(x)
  }

  agg <- diet_aggregates(df)
  dom <- pmax(agg$Diet_Inv, agg$Diet_VertFishScav, agg$Diet_AllPlants, na.rm = TRUE)
  out <- rep(NA_character_, nrow(df))
  out[is.finite(dom) & dom >= 50 & dom == agg$Diet_VertFishScav] <- "VertFishScav"
  out[is.finite(dom) & dom >= 50 & dom == agg$Diet_Inv] <- "Invertebrate"
  out[is.finite(dom) & dom >= 50 & dom == agg$Diet_AllPlants & agg$Diet_FruiNect >= agg$Diet_PlantSeed] <- "FruiNect"
  out[is.finite(dom) & dom >= 50 & dom == agg$Diet_AllPlants & agg$Diet_PlantSeed > agg$Diet_FruiNect] <- "PlantSeed"
  out[is.na(out) & is.finite(dom)] <- "Omnivore"
  out
}

attach_random_effects <- function(x, random_effects) {
  n_before <- nrow(x)
  out <- x |>
    dplyr::mutate(
      class_key = clean_name(taxon_class),
      order_key = clean_name(orderName),
      family_key = clean_name(familyName),
      species_key = clean_name(scientificName)
    ) |>
    dplyr::left_join(
      dplyr::distinct(dplyr::select(random_effects, Class, Order, Effect_Order)),
      by = c("class_key" = "Class", "order_key" = "Order")
    ) |>
    dplyr::left_join(
      dplyr::distinct(dplyr::select(random_effects, Class, Family, Effect_Family)),
      by = c("class_key" = "Class", "family_key" = "Family")
    ) |>
    dplyr::left_join(
      dplyr::distinct(dplyr::select(random_effects, Class, Species, Effect_Species)),
      by = c("class_key" = "Class", "species_key" = "Species")
    ) |>
    dplyr::mutate(
      Effect_Order = dplyr::coalesce(Effect_Order, 0),
      Effect_Family = dplyr::coalesce(Effect_Family, 0),
      Effect_Species = dplyr::coalesce(Effect_Species, 0)
    ) |>
    dplyr::select(-class_key, -order_key, -family_key, -species_key)
  assert(nrow(out) == n_before, "Random-effect joins unexpectedly changed the number of species rows.")
  out
}

compute_density_and_dispersal <- function(x, contract = canonical_analysis_contract()) {
  contract <- validate_analysis_contract(contract, "Stage 4 abundance contract")
  mass_g <- x$BodyMass.Value
  mass_kg <- mass_g / 1000
  mass_ok <- is.finite(mass_g) & mass_g > 0

  mammal_density <- rep(NA_real_, nrow(x))
  bird_density <- rep(NA_real_, nrow(x))
  mammal_hr <- rep(NA_real_, nrow(x))
  bird_hr <- rep(NA_real_, nrow(x))

  is_mammal <- x$taxon_class == "Mammalia" & mass_ok
  is_bird <- x$taxon_class == "Aves" & mass_ok

  mammal_diet_effect <- dplyr::case_when(
    x$Diet == "herbivore" ~ mammal_density_coefficients$diet[["herbivore"]],
    x$Diet == "omnivore" ~ mammal_density_coefficients$diet[["omnivore"]],
    x$Diet == "carnivore" ~ mammal_density_coefficients$diet[["carnivore"]],
    TRUE ~ 0
  )
  mammal_density[is_mammal] <- 10^(
    mammal_density_coefficients$intercept +
      x$Effect_Order[is_mammal] + x$Effect_Family[is_mammal] + x$Effect_Species[is_mammal] +
      mammal_diet_effect[is_mammal] + mammal_density_coefficients$mass * log10(mass_g[is_mammal])
  )
  mammal_hr[is_mammal] <- dplyr::case_when(
    x$Diet[is_mammal] == "carnivore" ~ 0.38 * mass_kg[is_mammal]^1.13,
    TRUE ~ 0.054 * mass_kg[is_mammal]
  )

  bird_diet_effect <- unname(bird_density_coefficients$diet[x$Diet])
  ok_bird <- is_bird & is.finite(bird_diet_effect)
  bird_density[ok_bird] <- 10^(
    bird_density_coefficients$intercept +
      x$Effect_Order[ok_bird] + x$Effect_Family[ok_bird] + x$Effect_Species[ok_bird] +
      bird_diet_effect[ok_bird] + bird_density_coefficients$mass * log10(mass_g[ok_bird])
  )
  bird_hr[is_bird] <- dplyr::case_when(
    x$Diet[is_bird] == "VertFishScav" ~ 210 * mass_kg[is_bird]^1.13,
    TRUE ~ 37 * mass_kg[is_bird]
  )

  x |>
    dplyr::mutate(
      density = dplyr::case_when(
        taxon_class == "Mammalia" ~ mammal_density,
        taxon_class == "Aves" ~ bird_density,
        TRUE ~ NA_real_
      ),
      home_range_size = dplyr::case_when(
        taxon_class == "Mammalia" ~ mammal_hr,
        taxon_class == "Aves" ~ bird_hr,
        TRUE ~ NA_real_
      ),
      dispersal_dist = dplyr::case_when(
        taxon_class == "Mammalia" ~ 5.6 * sqrt(mammal_hr),
        taxon_class == "Aves" ~ 12 * sqrt(bird_hr),
        TRUE ~ NA_real_
      ),
      min_patch_size = minimum_patch_abundance(contract) / density,
      min_pop_size = quasi_extinction_abundance(contract) / density
    )
}

predict_positive_trait_loess <- function(fit, trait_value) {
  tv <- as.numeric(trait_value)
  out <- rep(NA_real_, length(tv))
  ok <- is.finite(tv) & tv > 0
  if (!any(ok)) return(out)

  x_train <- fit$x
  if (is.matrix(x_train) || is.data.frame(x_train)) x_train <- x_train[, 1]
  rng <- range(as.numeric(x_train), finite = TRUE)
  logTrait <- pmin(pmax(log10(tv[ok]), rng[1]), rng[2])
  out[ok] <- exp(as.numeric(stats::predict(fit, newdata = data.frame(logTrait = logTrait))))
  out
}

report_loess_trait_range <- function(fit, trait_value, label) {
  training_x <- fit$x
  if (is.matrix(training_x) || is.data.frame(training_x)) training_x <- training_x[, 1]
  limits <- range(10^as.numeric(training_x), finite = TRUE)
  values <- as.numeric(trait_value)
  log_msg(
    "Stage 4 | LOESS range |", label,
    "| below=", sum(values < limits[1]),
    "| within=", sum(values >= limits[1] & values <= limits[2]),
    "| above=", sum(values > limits[2]),
    "| training=", paste(signif(limits, 5), collapse = "-")
  )
  invisible(limits)
}

