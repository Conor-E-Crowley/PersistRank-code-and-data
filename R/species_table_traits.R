# High-level Stage 4 trait-enrichment coordination.
#
# Loaded by stage4_workflow.R after covariate and Gompertz definitions. It
# assembles the existing enriched species table from validated inputs. Sourcing
# is side-effect free; all scientific work occurs only when called.

add_trait_covariates <- function(species_inputs, random_effects, models, curves,
                                 use_gompertz_mammals, use_gompertz_birds,
                                 contract = canonical_analysis_contract()) {
  contract <- validate_analysis_contract(contract, "Stage 4 abundance contract")
  diet <- diet_aggregates(species_inputs)

  out <- species_inputs |>
    dplyr::mutate(
      Diet = dplyr::case_when(
        taxon_class == "Mammalia" ~ classify_mammal_diet(diet$Diet_AllAnimal, diet$Diet_AllPlants, BodyMass.Value),
        taxon_class == "Aves" ~ get_bird_diet5(species_inputs),
        TRUE ~ NA_character_
      ),
      diet5_group = dplyr::if_else(
        taxon_class == "Aves",
        Diet,
        NA_character_
      ),
      bird_sigma_model_group = NA_character_
    ) |>
    attach_random_effects(random_effects) |>
    compute_density_and_dispersal(contract = contract)

  bad_diet <- is.na(out$Diet) | !nzchar(as.character(out$Diet))
  if (any(bad_diet)) {
    species_abort(
      "Retained species lack a recognized diet category:\n",
      paste(utils::head(out$scientificName[bad_diet], 25L), collapse = "\n"),
      if (sum(bad_diet) > 25L) "\n..." else ""
    )
  }

  if (!is.null(models)) {
    add_gompertz_parameters(
      out,
      models = models,
      curves = curves,
      use_mammals = use_gompertz_mammals,
      use_birds = use_gompertz_birds,
      contract = contract
    )
  } else {
    for (curve in curves) {
      out[[paste0("alpha_", curve)]] <- NA_real_
      out[[paste0("beta_", curve)]] <- NA_real_
    }
    out
  }
}

