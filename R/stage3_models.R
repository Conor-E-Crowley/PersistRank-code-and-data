# Saved-model contracts and provenance for Stage 3.
#
# Loaded by stage3_workflow.R after fitting and smoothing definitions. It builds
# and validates the persisted model object and reuse handoff. Sourcing is
# definition-only; model reads remain owned by the workflow.

validate_stage3_model_provenance <- function(provenance, k0, label) {
  required <- stage3_shared_provenance_columns()
  assert(is.list(provenance) && all(required %in% names(provenance)),
         paste0(label, " lacks complete Stage 2 scientific provenance."))

  positive_integer <- c("persistence_horizon_years", "quasi_extinction_abundance", "n_draws", "reps", "chunk_size")
  for (field in positive_integer) {
    value <- suppressWarnings(as.numeric(provenance[[field]]))
    assert(length(value) == 1L && is.finite(value) && value > 0 && value == floor(value),
           paste0(label, " has invalid Stage 2 provenance: ", field, "."))
  }
  base_seed <- suppressWarnings(as.numeric(provenance$base_seed))
  assert(length(base_seed) == 1L && is.finite(base_seed) && base_seed >= 0 && base_seed == floor(base_seed),
         paste0(label, " has invalid Stage 2 provenance: base_seed."))

  cap_factor <- suppressWarnings(as.numeric(provenance$cap_factor))
  r_buffer <- suppressWarnings(as.numeric(provenance$r_buffer))
  assert(length(cap_factor) == 1L && is.finite(cap_factor) && cap_factor >= 1,
         paste0(label, " has invalid Stage 2 provenance: cap_factor."))
  assert(length(r_buffer) == 1L && is.finite(r_buffer) && r_buffer > 0 && r_buffer <= 1,
         paste0(label, " has invalid Stage 2 provenance: r_buffer."))
  assert(stage3_numeric_equal(provenance$quasi_extinction_abundance, k0),
         paste0(label, " has Stage 2 provenance inconsistent with Nq."))

  for (field in c("grid_signature", "simulator_contract")) {
    value <- trimws(as.character(provenance[[field]]))
    assert(length(value) == 1L && !is.na(value) && nzchar(value),
           paste0(label, " has invalid Stage 2 provenance: ", field, "."))
  }
  assert(identical(as.character(provenance$posterior_sampling), "independent_without_replacement"),
         paste0(label, " has invalid Stage 2 provenance: posterior_sampling."))
  invisible(TRUE)
}

validate_stage3_gompertz_fit_table <- function(fits, models, present, curves, label) {
  required <- c(
    "group", "predictor", "bird_sigma_model_group", "trait_idx", "trait_value", "curve",
    "alpha", "beta", "fit_ok", "n_pts", "rmse_lin", "rmse_prob", "r2_prob",
    "max_abs_error", "fit_msg"
  )
  assert(is.data.frame(fits), paste0(label, " must contain a gompertz_fits table."))
  need_cols(fits, required, paste0(label, " gompertz_fits"))
  assert(nrow(fits) > 0L, paste0(label, " gompertz_fits is empty."))

  group <- as.character(fits$group)
  predictor <- as.character(fits$predictor)
  bird_group <- as.character(fits$bird_sigma_model_group)
  curve <- as.character(fits$curve)
  trait_idx <- suppressWarnings(as.numeric(fits$trait_idx))
  trait_value <- suppressWarnings(as.numeric(fits$trait_value))
  assert(setequal(unique(group), present) && all(group %in% present),
         paste0(label, " gompertz_fits does not match the present taxa."))
  assert(all(is.finite(trait_idx) & trait_idx > 0 & trait_idx == floor(trait_idx)) &&
           all(is.finite(trait_value) & trait_value > 0),
         paste0(label, " gompertz_fits contains invalid trait metadata."))
  assert(all(curve %in% curves), paste0(label, " gompertz_fits contains an invalid curve."))
  assert(all(predictor[group == "Mammals"] == "Mass_g") &&
           all(predictor[group == "Birds"] == "GenLength"),
         paste0(label, " gompertz_fits contains invalid predictor metadata."))
  assert(all(is.na(fits$bird_sigma_model_group[group == "Mammals"])),
         paste0(label, " mammal Gompertz fits must not have bird-group values."))
  if (any(group == "Birds")) {
    assert(all(!is.na(bird_group[group == "Birds"]) & nzchar(bird_group[group == "Birds"])),
           paste0(label, " gompertz_fits contains invalid bird groups."))
  }

  assert(is.logical(fits$fit_ok) && all(!is.na(fits$fit_ok) & fits$fit_ok),
         paste0(label, " gompertz_fits contains unsuccessful fits."))
  numeric_fields <- c("alpha", "beta", "n_pts", "rmse_lin", "rmse_prob", "r2_prob", "max_abs_error")
  assert(all(vapply(fits[numeric_fields], function(x) all(is.finite(suppressWarnings(as.numeric(x)))), logical(1))),
         paste0(label, " gompertz_fits contains nonfinite coefficients or diagnostics."))
  assert(all(fits$alpha > 0 & fits$beta > 0) &&
           all(fits$n_pts >= 3 & fits$n_pts == floor(fits$n_pts)) &&
           all(fits$rmse_lin >= 0 & fits$rmse_prob >= 0 & fits$max_abs_error >= 0),
         paste0(label, " gompertz_fits contains invalid coefficients or diagnostics."))
  assert(is.character(fits$fit_msg) && all(!is.na(fits$fit_msg)),
         paste0(label, " gompertz_fits contains invalid fit messages."))

  bird_key <- ifelse(is.na(bird_group), "", bird_group)
  keys <- paste(group, bird_key, trait_idx, curve, sep = "\r")
  assert(!anyDuplicated(keys), paste0(label, " gompertz_fits contains duplicate model blocks."))
  trait_keys <- paste(group, bird_key, trait_idx, sep = "\r")
  trait_map <- unique(data.frame(key = trait_keys, value = trait_value, stringsAsFactors = FALSE))
  assert(!anyDuplicated(trait_map$key), paste0(label, " gompertz_fits maps a trait index to multiple values."))
  coverage <- split(curve, trait_keys)
  assert(all(vapply(coverage, function(x) identical(sort(unique(x)), sort(curves)), logical(1))),
         paste0(label, " gompertz_fits lacks complete five-curve coverage."))

  validate_branch_counts <- function(row_ids, branch, branch_label) {
    for (curve_name in curves) {
      block_ids <- row_ids[curve[row_ids] == curve_name]
      block <- fits[block_ids, , drop = FALSE]
      stored <- branch[[curve_name]]
      assert(nrow(block) == stored$n &&
               stage3_numeric_equal(range(block$trait_value), stored$trait_range),
             paste0(label, " gompertz_fits is inconsistent with the ", branch_label,
                    " LOESS metadata for ", curve_name, "."))
    }
  }
  if ("Mammals" %in% present) validate_branch_counts(which(group == "Mammals"), models$mammals, "mammal")
  if ("Birds" %in% present) {
    for (bird_name in sort(unique(bird_group[group == "Birds"]))) {
      validate_branch_counts(which(group == "Birds" & bird_group == bird_name), models$birds[[bird_name]], paste0("bird ", bird_name))
    }
  }
  invisible(TRUE)
}

validate_gompertz_loess_models <- function(models, run_mammals = NULL, run_birds = NULL,
                                           curves = persistence_curves(), k0 = NULL,
                                           label = "persistence_curve_models.rds") {
  assert(is.list(models) && is.list(models$meta), paste0(label, " must contain a meta list."))
  meta <- models$meta
  assert(identical(meta$schema_version, gompertz_model_schema()),
         paste0(label, " schema must be ", gompertz_model_schema(), "."))
  assert(identical(as.character(meta$curves), as.character(curves)), paste0(label, " has an invalid curve order."))
  assert(identical(names(meta$curve_probabilities), curves) &&
           stage3_numeric_equal(unname(meta$curve_probabilities), unname(persistence_quantiles()[curves])),
         paste0(label, " has invalid curve probabilities."))
  assert(identical(meta$main_curve, main_persistence_curve()), paste0(label, " has an invalid main curve."))
  expected_k0 <- if (is.null(k0)) meta$K0 else k0
  assert(
    length(expected_k0) == 1L &&
      is.finite(expected_k0) &&
      expected_k0 > 0 &&
      expected_k0 == floor(expected_k0) &&
      stage3_numeric_equal(meta$K0, expected_k0),
    paste0(label, " has an invalid Nq.")
  )
  assert(identical(unname(meta$predictors[["Mammals"]]), "Mass_g") &&
           identical(unname(meta$predictors[["Birds"]]), "GenLength"),
         paste0(label, " has invalid predictor metadata."))
  assert(identical(meta$x_transform, "log10(trait_value)") && identical(meta$y_transform, "log(parameter)"),
         paste0(label, " has invalid transformation metadata."))
  assert(length(meta$created_utc) == 1L && !is.na(meta$created_utc) && nzchar(trimws(as.character(meta$created_utc))),
         paste0(label, " has invalid creation-time metadata."))
  assert(length(meta$loess_span) == 1L && is.finite(meta$loess_span) && meta$loess_span > 0 && meta$loess_span <= 1 &&
           identical(as.integer(meta$loess_degree), 2L) && identical(as.character(meta$loess_family), "gaussian"),
         paste0(label, " has invalid LOESS metadata."))
  assert(length(meta$loess_z_display) == 1L && is.finite(meta$loess_z_display) && meta$loess_z_display > 0,
         paste0(label, " has invalid display-band metadata."))
  present <- as.character(meta$present_taxa)
  selected <- as.character(meta$selected_taxa)
  allowed_taxa <- c("Mammals", "Birds")
  assert(length(present) > 0L && identical(selected, present) &&
           identical(present, allowed_taxa[allowed_taxa %in% present]) && !anyDuplicated(present),
         paste0(label, " has invalid present_taxa metadata."))
  expected_sources <- tolower(present)
  assert(is.character(meta$source_files) && identical(names(meta$source_files), expected_sources) &&
           all(!is.na(meta$source_files) & nzchar(trimws(meta$source_files))),
         paste0(label, " has invalid source-file metadata."))
  assert(is.character(meta$source_md5) && identical(names(meta$source_md5), expected_sources) &&
           all(grepl("^[[:xdigit:]]{32}$", meta$source_md5)),
         paste0(label, " has invalid source-checksum metadata."))
  validate_stage3_model_provenance(meta$stage2_provenance, meta$K0, label)

  if (is.null(run_mammals)) run_mammals <- "Mammals" %in% present
  if (is.null(run_birds)) run_birds <- "Birds" %in% present
  assert(is.logical(run_mammals) && length(run_mammals) == 1L && !is.na(run_mammals) &&
           is.logical(run_birds) && length(run_birds) == 1L && !is.na(run_birds),
         "Requested Stage 3 taxa must be TRUE/FALSE.")
  assert(is.null(models$mammals) == !("Mammals" %in% present) && is.null(models$birds) == !("Birds" %in% present),
         paste0(label, " has model branches inconsistent with present_taxa."))
  if ("Birds" %in% present) {
    validate_bird_model_spec(meta$bird_model_spec)
    bird_groups <- meta$bird_model_spec$model_groups
  } else {
    bird_groups <- character()
    assert(
      is.null(meta$bird_model_spec),
      paste0(label, " has an unexpected bird model specification.")
    )
  }
  if (isTRUE(run_mammals)) {
    assert("Mammals" %in% present, paste0(label, " does not contain the requested mammal branch."))
  }
  if ("Mammals" %in% present) {
    validate_loess_curve_set(models$mammals, curves, meta$loess_z_display, paste0(label, " mammal branch"))
  }
  if (isTRUE(run_birds)) {
    assert("Birds" %in% present, paste0(label, " does not contain the requested bird branch."))
  }
  if ("Birds" %in% present) {
    assert(is.list(models$birds) && identical(names(models$birds), bird_groups),
           paste0(label, " has invalid ordered bird model branches."))
    for (group in bird_groups) {
      validate_loess_curve_set(
        models$birds[[group]], curves, meta$loess_z_display,
        paste0(label, " bird ", group, " branch")
      )
    }
  }
  validate_stage3_gompertz_fit_table(models$gompertz_fits, models, present, curves, label)
  invisible(TRUE)
}

stage3_source_metadata <- function(inputs) {
  active <- Filter(Negate(is.null), inputs)
  list(
    paths = vapply(active, `[[`, character(1), "source"),
    md5 = vapply(active, `[[`, character(1), "md5")
  )
}

build_gompertz_loess_models <- function(gomp_params, inputs, config) {
  assert_all_gompertz_fits(gomp_params)
  curves <- config$curves
  mammals <- gomp_params[gomp_params$group == "Mammals", , drop = FALSE]
  birds <- gomp_params[gomp_params$group == "Birds", , drop = FALSE]

  mammal_models <- if (config$selected_mammals) fit_loess_set(mammals, curves, config$loess$span) else NULL
  bird_model_spec <- if (config$selected_birds) inputs$birds$bird_model else NULL
  bird_groups <- if (config$selected_birds) bird_model_spec$model_groups else character()
  bird_models <- if (config$selected_birds) stats::setNames(lapply(bird_groups, function(group) {
    fit_loess_set(
      birds[birds$bird_sigma_model_group == group, , drop = FALSE],
      curves,
      config$loess$span
    )
  }), bird_groups) else NULL

  source_metadata <- stage3_source_metadata(inputs)
  provenance_source <- if (!is.null(inputs$mammals)) inputs$mammals$provenance else inputs$birds$provenance
  present <- c(if (config$selected_mammals) "Mammals", if (config$selected_birds) "Birds")
  models <- list(
    meta = list(
      schema_version = gompertz_model_schema(),
      created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
      K0 = config$k0,
      curves = curves,
      curve_probabilities = config$curve_probabilities,
      main_curve = config$main_curve,
      selected_taxa = present,
      present_taxa = present,
      predictors = c(Mammals = "Mass_g", Birds = "GenLength"),
      bird_model_spec = bird_model_spec,
      loess_span = config$loess$span,
      loess_degree = config$loess$degree,
      loess_family = config$loess$family,
      loess_z_display = config$loess$z,
      x_transform = "log10(trait_value)",
      y_transform = "log(parameter)",
      source_files = source_metadata$paths,
      source_md5 = source_metadata$md5,
      stage2_provenance = provenance_source
    ),
    gompertz_fits = gomp_params,
    mammals = mammal_models,
    birds = bird_models
  )

  validate_gompertz_loess_models(models, config$selected_mammals, config$selected_birds, curves, config$k0,
                                 "computed Stage 3 model")
  models
}

predict_saved_stage3_models <- function(models, config) {
  validate_gompertz_loess_models(
    models,
    config$selected_mammals,
    config$selected_birds,
    config$curves,
    config$k0
  )
  fits <- models$gompertz_fits
  bird_groups <- if (config$selected_birds) {
    models$meta$bird_model_spec$model_groups
  } else {
    character()
  }
  dplyr::bind_rows(
    if (config$selected_mammals) {
      mammal <- fits[fits$group == "Mammals", , drop = FALSE]
      predict_loess_set(
        mammal,
        models$mammals,
        config$curves,
        config$loess$z,
        config$loess$prediction_points
      ) |>
        dplyr::mutate(
          group = "Mammals",
          predictor = "Mass_g",
          bird_sigma_model_group = NA_character_
        )
    } else {
      tibble::tibble()
    },
    if (config$selected_birds) {
      dplyr::bind_rows(lapply(bird_groups, function(group) {
        branch <- fits[
          fits$group == "Birds" & fits$bird_sigma_model_group == group,
          ,
          drop = FALSE
        ]
        predict_loess_set(
          branch,
          models$birds[[group]],
          config$curves,
          config$loess$z,
          config$loess$prediction_points
        ) |>
          dplyr::mutate(
            group = "Birds",
            predictor = "GenLength",
            bird_sigma_model_group = .env$group
          )
      }))
    } else {
      tibble::tibble()
    }
  )
}

validate_stage3_reuse_provenance <- function(models, inputs, config) {
  active <- Filter(Negate(is.null), inputs)
  source_md5 <- vapply(active, `[[`, character(1), "md5")
  saved_md5 <- models$meta$source_md5[names(source_md5)]
  assert(
    identical(names(saved_md5), names(source_md5)) &&
      identical(unname(saved_md5), unname(source_md5)),
    paste0(
      "The saved Stage 3 model does not match the current Stage 2 inputs. ",
      "Run Stage 3 with mode: 'fit'."
    )
  )
  assert(
    stage3_numeric_equal(models$meta$loess_span, config$loess$span),
    paste0(
      "The saved Stage 3 model uses loess_span=",
      format(models$meta$loess_span, digits = 17),
      ", but params$loess_span=", format(config$loess$span, digits = 17),
      ". Run Stage 3 with mode: 'fit' to use the requested span."
    )
  )
  invisible(TRUE)
}

# ---- Validated model/figure transaction --------------------------------------
#
# Stage 3 produces one scientific model plus a figure set whose membership
# depends on the selected taxa and dynamic bird branches. All requested files
# are rendered and validated in one staging directory, then committed together.
# Reuse mode omits the model from the transaction and regenerates figures only.
