# Stage 4 Gompertz-model validation and parameter enrichment.
#
# Loaded by stage4_workflow.R after trait covariates. It requires a validated
# Stage 3 model contract and explicit taxon/curve selectors. Sourcing has no
# side effects; prediction occurs only when enrichment is requested.

validate_gompertz_model_contract <- function(models, curves, use_mammals = TRUE,
                                             use_birds = TRUE,
                                             contract = canonical_analysis_contract()) {
  contract <- validate_analysis_contract(contract, "Stage 4 abundance contract")
  if (is.null(models)) {
    species_abort("Gompertz LOESS models are required but were not loaded.")
  }

  expected_version <- "persistence_curve_models_v6_dynamic_bird"
  found_version <- models$meta$schema_version %||% NA_character_
  if (!identical(found_version, expected_version)) {
    species_abort(
      "Unexpected persistence_curve_models.rds schema. Expected meta$schema_version == ",
      expected_version,
      "; found ",
      if (length(found_version) && !is.na(found_version)) found_version else "missing",
      ". Regenerate Stage 3 before building species_table.csv."
    )
  }

  expected_curves <- persistence_curves()
  expected_probabilities <- persistence_quantiles()[expected_curves]
  if (!identical(as.character(models$meta$curves), expected_curves) ||
      !identical(as.character(curves), expected_curves)) {
    species_abort("Gompertz model metadata must declare the ordered five-curve analysis contract.")
  }
  found_probabilities <- suppressWarnings(as.numeric(models$meta$curve_probabilities[expected_curves]))
  if (length(found_probabilities) != length(expected_probabilities) ||
      any(!is.finite(found_probabilities)) ||
      any(abs(found_probabilities - unname(expected_probabilities)) > sqrt(.Machine$double.eps))) {
    species_abort("Gompertz model metadata has invalid persistence-curve probabilities.")
  }
  if (!identical(models$meta$main_curve, main_persistence_curve())) {
    species_abort("Gompertz model metadata does not declare q50 as the main curve.")
  }
  model_k0 <- suppressWarnings(as.numeric(models$meta$K0))
  if (length(model_k0) != 1L || !is.finite(model_k0) ||
      model_k0 != quasi_extinction_abundance(contract)) {
    species_abort("Gompertz model metadata does not use the requested quasi-extinction abundance.")
  }
  if (!identical(models$meta$x_transform, "log10(trait_value)") ||
      !identical(models$meta$y_transform, "log(parameter)")) {
    species_abort("Gompertz model metadata has unexpected trait or parameter transformations.")
  }
  created <- as.character(models$meta$created_utc %||% NA_character_)
  if (length(created) != 1L || is.na(created) || !nzchar(trimws(created))) {
    species_abort("Gompertz model metadata lacks a valid creation time.")
  }
  loess_span <- suppressWarnings(as.numeric(models$meta$loess_span))
  loess_degree <- suppressWarnings(as.numeric(models$meta$loess_degree))
  loess_z <- suppressWarnings(as.numeric(models$meta$loess_z_display))
  if (length(loess_span) != 1L || !is.finite(loess_span) || loess_span <= 0 || loess_span > 1 ||
      length(loess_degree) != 1L || !is.finite(loess_degree) || loess_degree != 2 ||
      !identical(as.character(models$meta$loess_family), "gaussian") ||
      length(loess_z) != 1L || !is.finite(loess_z) || loess_z <= 0) {
    species_abort("Gompertz model metadata has invalid LOESS settings.")
  }
  required_provenance <- c(
    "persistence_horizon_years", "quasi_extinction_abundance", "cap_factor", "r_buffer",
    "n_draws", "reps", "chunk_size", "base_seed", "grid_signature", "posterior_sampling",
    "demographic_uncertainty", "simulator_contract"
  )
  if (!is.list(models$meta$stage2_provenance) ||
      !all(required_provenance %in% names(models$meta$stage2_provenance))) {
    species_abort("Gompertz model metadata lacks complete Stage 2 scientific provenance.")
  }

  present_taxa <- as.character(models$meta$present_taxa %||% character())
  selected_taxa <- as.character(models$meta$selected_taxa %||% character())
  if (!identical(selected_taxa, present_taxa) || !length(present_taxa) ||
      anyDuplicated(present_taxa) || !all(present_taxa %in% c("Mammals", "Birds"))) {
    species_abort("Gompertz model metadata has inconsistent selected_taxa/present_taxa values.")
  }

  source_files <- models$meta$source_files %||% character()
  source_md5 <- models$meta$source_md5 %||% character()
  expected_source_names <- tolower(present_taxa)
  if (!identical(names(source_files), expected_source_names) ||
      !identical(names(source_md5), expected_source_names) ||
      any(is.na(source_files) | !nzchar(trimws(as.character(source_files)))) ||
      any(is.na(source_md5) | !grepl("^[[:xdigit:]]{32}$", as.character(source_md5)))) {
    species_abort("Gompertz model metadata must provide one named source path and MD5 checksum per present taxon.")
  }

  provenance <- models$meta$stage2_provenance
  positive_integer <- c(
    "persistence_horizon_years", "quasi_extinction_abundance", "n_draws", "reps", "chunk_size"
  )
  for (field in positive_integer) {
    value <- suppressWarnings(as.numeric(provenance[[field]]))
    if (length(value) != 1L || !is.finite(value) || value <= 0 || value != floor(value)) {
      species_abort("Gompertz model metadata has invalid Stage 2 provenance: ", field, ".")
    }
  }
  base_seed <- suppressWarnings(as.numeric(provenance$base_seed))
  cap_factor <- suppressWarnings(as.numeric(provenance$cap_factor))
  r_buffer <- suppressWarnings(as.numeric(provenance$r_buffer))
  nonblank_scalar <- function(x) {
    x <- as.character(x)
    length(x) == 1L && !is.na(x) && nzchar(trimws(x))
  }
  if (length(base_seed) != 1L || !is.finite(base_seed) || base_seed < 0 || base_seed != floor(base_seed) ||
      length(cap_factor) != 1L || !is.finite(cap_factor) || cap_factor < 1 ||
      length(r_buffer) != 1L || !is.finite(r_buffer) || r_buffer <= 0 || r_buffer > 1 ||
      !identical(as.character(provenance$posterior_sampling), "independent_without_replacement") ||
      !nonblank_scalar(provenance$grid_signature) ||
      !nonblank_scalar(provenance$simulator_contract)) {
    species_abort("Gompertz model metadata contains invalid Stage 2 scientific provenance values.")
  }
  if (as.numeric(provenance$quasi_extinction_abundance) != as.numeric(models$meta$K0)) {
    species_abort("Gompertz model Stage 2 provenance and model metadata use different quasi-extinction abundances.")
  }

  if (isTRUE(use_mammals) && !"Mammals" %in% present_taxa) {
    species_abort("The Stage 3 model does not contain the requested mammal branch.")
  }
  if (isTRUE(use_birds) && !"Birds" %in% present_taxa) {
    species_abort("The Stage 3 model does not contain the requested bird branch.")
  }
  if (!"Mammals" %in% present_taxa && !is.null(models$mammals)) {
    species_abort("The Stage 3 model contains a mammal branch absent from present_taxa metadata.")
  }
  if (!"Birds" %in% present_taxa && !is.null(models$birds)) {
    species_abort("The Stage 3 model contains a bird branch absent from present_taxa metadata.")
  }

  predictors <- models$meta$predictors %||% character()
  if (isTRUE(use_mammals) && !identical(unname(predictors[["Mammals"]]), "Mass_g")) {
    species_abort("Gompertz model metadata does not declare Mammals predictor as Mass_g.")
  }
  if (isTRUE(use_birds) && !identical(unname(predictors[["Birds"]]), "GenLength")) {
    species_abort("Gompertz model metadata does not declare Birds predictor as GenLength.")
  }

  check_curve_set <- function(x, label) {
    found <- names(x)
    if (!is.list(x) || is.null(found) || !setequal(found, curves)) {
      species_abort(
        "Gompertz model object has an invalid ", label, " curve structure. Expected direct curve names: ",
        paste(curves, collapse = ", "),
        ". Found: ",
        if (length(found)) paste(found, collapse = ", ") else "none",
        ". Expected the current Stage 3 model structure."
      )
    }
    for (curve in curves) {
      fit <- x[[curve]]
      if (!is.list(fit) || !inherits(fit$alpha, "loess") || !inherits(fit$beta, "loess")) {
        species_abort("Gompertz model object has unusable alpha/beta LOESS fits for ", label, " ", curve, ".")
      }
      for (parameter in c("alpha", "beta")) {
        trait_range <- suppressWarnings(as.numeric(fit$trait_range))
        if (length(trait_range) != 2L || any(!is.finite(trait_range)) ||
            any(trait_range <= 0) || trait_range[1] >= trait_range[2]) {
          species_abort("Gompertz model has an invalid trait range for ", label, " ", curve, ".")
        }
        trait_grid <- 10^seq(log10(trait_range[1]), log10(trait_range[2]), length.out = 101L)
        prediction <- exp(as.numeric(stats::predict(
          fit[[parameter]],
          newdata = data.frame(logTrait = log10(trait_grid))
        )))
        if (any(!is.finite(prediction) | prediction <= 0)) {
          species_abort("Gompertz model has nonpositive or nonfinite ", parameter,
                        " predictions for ", label, " ", curve, ".")
        }
      }
    }
  }

  if (isTRUE(use_mammals)) check_curve_set(models$mammals, "mammal")
  if (isTRUE(use_birds)) {
    bird_model <- models$meta$bird_model_spec
    validate_bird_model_spec(bird_model)
    if (!is.list(models$birds) || !identical(names(models$birds), bird_model$model_groups)) {
      species_abort(
        "Gompertz model object has an invalid bird model-group structure. Expected: ",
        paste(bird_model$model_groups, collapse = ", "), "."
      )
    }
    for (diet_group in bird_model$model_groups) {
      check_curve_set(models$birds[[diet_group]], paste0("bird ", diet_group))
    }
  }

  invisible(TRUE)
}

add_gompertz_parameters <- function(x, models, curves, use_mammals = TRUE,
                                    use_birds = TRUE,
                                    contract = canonical_analysis_contract()) {
  for (curve in curves) {
    x[[paste0("alpha_", curve)]] <- NA_real_
    x[[paste0("beta_", curve)]] <- NA_real_
  }

  if (isTRUE(use_mammals) || isTRUE(use_birds)) {
    validate_gompertz_model_contract(
      models, curves, use_mammals, use_birds, contract = contract
    )
  }

  if (isTRUE(use_mammals)) {
    ii <- x$taxon_class == "Mammalia" & is.finite(x$BodyMass.Value) & x$BodyMass.Value > 0
    if (any(ii)) {
      report_loess_trait_range(
        models$mammals[[curves[[1L]]]]$alpha,
        x$BodyMass.Value[ii],
        "mammals"
      )
    }
    for (curve in curves) {
      x[[paste0("alpha_", curve)]][ii] <- predict_positive_trait_loess(models$mammals[[curve]]$alpha, x$BodyMass.Value[ii])
      x[[paste0("beta_", curve)]][ii] <- predict_positive_trait_loess(models$mammals[[curve]]$beta, x$BodyMass.Value[ii])
    }
  }

  if (isTRUE(use_birds)) {
    ii <- x$taxon_class == "Aves" & is.finite(x$GenLength) & x$GenLength > 0
    bird_model <- models$meta$bird_model_spec
    validate_bird_diet5(
      x$diet5_group[ii],
      "species_table diet5_group"
    )
    bird_model_group <- rep(NA_character_, nrow(x))
    bird_model_group[ii] <- map_bird_diet_to_model_group(
      x$diet5_group[ii],
      bird_model
    )
    x$bird_sigma_model_group <- bird_model_group
    for (diet_group in bird_model$model_groups) {
      jj <- ii & bird_model_group == diet_group
      if (!any(jj)) next
      report_loess_trait_range(
        models$birds[[diet_group]][[curves[[1L]]]]$alpha,
        x$GenLength[jj],
        paste("birds", diet_group)
      )
      for (curve in curves) {
        x[[paste0("alpha_", curve)]][jj] <- predict_positive_trait_loess(
          models$birds[[diet_group]][[curve]]$alpha,
          x$GenLength[jj]
        )
        x[[paste0("beta_", curve)]][jj] <- predict_positive_trait_loess(
          models$birds[[diet_group]][[curve]]$beta,
          x$GenLength[jj]
        )
      }
    }
  }

  x
}

