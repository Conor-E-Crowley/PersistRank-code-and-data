# Shared Stage 7 persistence calculations.
#
# Patch areas are aggregated to population units, evaluated with the selected
# Gompertz coefficients, then combined to species persistence using a fixed
# retained-species denominator. Inputs remain in canonical stage/species order;
# absent species receive persistence zero in assemblage summaries. The Stage
# 7.1 workflow owns loading this definition-only scientific module; sourcing it
# performs no package attachment, filesystem access, or native compilation.


# Scientific persistence equations shared by the pipeline-only and benchmark
# comparison reports.
stage7_pu_persistence <- function(area_km2, density, alpha, beta, threshold_area,
                                  label = "PU persistence inputs") {
  area <- as.numeric(area_km2)
  density <- as.numeric(density)
  alpha <- as.numeric(alpha)
  beta <- as.numeric(beta)
  delta <- area - as.numeric(threshold_area)
  ok <- is.finite(area) & is.finite(delta) & delta > 0 &
    is.finite(density) & density > 0 & is.finite(alpha) & alpha > 0 &
    is.finite(beta) & beta > 0
  assert(
    all(ok),
    paste0(label, " must be finite, positive, and strictly above the population-area threshold.")
  )
  out <- exp(-exp(log(alpha) - beta * log(density) - beta * log(delta)))
  assert(
    all(is.finite(out) & out >= 0 & out <= 1),
    paste0(label, " produced invalid persistence values.")
  )
  out
}

stage7_species_persistence <- function(probability) {
  p <- as.numeric(probability)
  assert(
    all(is.finite(p) & p >= 0 & p <= 1),
    "PU persistence values must be finite and in [0,1]."
  )
  if (!length(p)) return(0)
  if (any(p >= 1)) return(1)
  out <- -expm1(sum(log1p(-p)))
  assert(is.finite(out) && out >= 0 && out <= 1, "Species persistence is invalid.")
  out
}

# Normalize the compact taxon and SDM metadata stored in shared initialization.
stage71_taxon_from_bundle <- function(x, label = "bundle taxon_class") {
  value <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(value))
  out[value == "mammalia"] <- "mammals"
  out[value == "aves"] <- "birds"
  bad <- is.na(out)
  assert(
    !any(bad),
    paste0(label, " contains unsupported value(s): ", paste(unique(value[bad]), collapse = ", "))
  )
  out
}

stage71_sdm_from_bundle <- function(x, label = "bundle sdm_method") {
  value <- tolower(trimws(as.character(x)))
  assert(
    all(value %in% c("ppm", "rangebag")),
    paste0(label, " must contain only ppm or rangebag.")
  )
  value
}

# Build the fixed species denominator from the exact shared Stage 6 snapshot,
# selecting only the requested curve's compact coefficient pair in memory.
prepare_stage71_bundle_parameters <- function(bundle, config) {
  contract <- config$contract %||% canonical_analysis_contract()
  contract <- validate_analysis_contract(contract, "Stage 7.1 abundance contract")
  validate_priority_initialization_schema(bundle, "Stage 7.1 shared initialization")
  retained <- validate_priority_initialization_metadata(
    bundle,
    taxa = config$taxa_tag,
    sdm = config$sdm,
    contract = contract,
    label = "Stage 7.1 shared initialization"
  )
  selected_species_params <- priority_species_parameters_for_curve(
    bundle$species_params,
    config$curve,
    label = "Stage 7.1 initialization species_params",
    contract = contract
  )

  params <- selected_species_params[, .(
    scientificName = trimws(as.character(species)),
    species = species_id(trimws(as.character(species))),
    taxon = stage71_taxon_from_bundle(taxon_class),
    sdm_method = stage71_sdm_from_bundle(sdm_method),
    density = as.numeric(density),
    threshold_area_km2 = as.numeric(min_population_area_km2),
    alpha = as.numeric(a_pred),
    beta = as.numeric(b_pred)
  )]
  assert(
    nrow(params) == length(retained) &&
      !anyDuplicated(params$scientificName) &&
      !anyDuplicated(params$species) &&
      setequal(params$scientificName, retained),
    "Stage 7.1 initialization species parameters do not exactly match metadata$retained_species."
  )

  metadata <- data.table::as.data.table(bundle$metadata$retained_species_table)
  need_cols(
    metadata,
    c("species", "taxon_class", "sdm_method"),
    "Stage 7.1 bundle metadata$retained_species_table"
  )
  metadata <- metadata[, .(
    scientificName = trimws(as.character(species)),
    metadata_taxon = stage71_taxon_from_bundle(taxon_class, "retained-species taxon_class"),
    metadata_sdm = stage71_sdm_from_bundle(sdm_method, "retained-species sdm_method")
  )]
  assert(
    nrow(metadata) == nrow(params) &&
      !anyDuplicated(metadata$scientificName) &&
      setequal(metadata$scientificName, params$scientificName),
    "Stage 7.1 retained-species metadata does not exactly match bundle species_params."
  )
  metadata <- metadata[match(params$scientificName, scientificName)]
  assert(
    identical(metadata$metadata_taxon, params$taxon) &&
      identical(metadata$metadata_sdm, params$sdm_method),
    "Stage 7.1 retained-species metadata disagrees with bundle species_params."
  )

  selectors <- priority_flags_from_tags(config$taxa_tag, config$sdm)
  expected_taxa <- c(
    if (selectors$do_mammals) "mammals",
    if (selectors$do_birds) "birds"
  )
  expected_sdm <- c(
    if (selectors$do_ppm) "ppm",
    if (selectors$do_rangebag) "rangebag"
  )
  assert(
    all(params$taxon %in% expected_taxa),
    "Stage 7.1 bundle species include a taxon outside the configured selection."
  )
  assert(
    all(params$sdm_method %in% expected_sdm),
    "Stage 7.1 bundle species include an SDM branch outside the configured selection."
  )

  expected_threshold <- quasi_extinction_abundance(contract) / params$density
  tolerance <- 1e-8 * pmax(1, expected_threshold)
  assert(
    all(is.finite(params$density) & params$density > 0) &&
      all(is.finite(params$threshold_area_km2) & params$threshold_area_km2 > 0) &&
      all(is.finite(params$alpha) & params$alpha > 0) &&
      all(is.finite(params$beta) & params$beta > 0) &&
      all(abs(params$threshold_area_km2 - expected_threshold) <= tolerance),
    "Stage 7.1 bundle parameters are invalid or inconsistent with the shared quasi-extinction threshold."
  )

  params[, parameter_order := .I]
  params[]
}

# Prepare all five figure-only curve parameter sets while retaining initialization-owned
# density, thresholds, species membership, taxon, and SDM metadata.
prepare_stage71_curve_parameters <- function(params, coefficient_table,
                                             configured_curve) {
  base <- data.table::copy(data.table::as.data.table(params))
  curves <- persistence_curves()
  configured_curve <- validate_persistence_curve(
    configured_curve,
    "Stage 7.1 configured curve"
  )
  coefficient_columns <- unlist(lapply(
    curves,
    function(curve) c(paste0("alpha_", curve), paste0("beta_", curve))
  ), use.names = FALSE)
  coefficients <- data.table::copy(data.table::as.data.table(coefficient_table))
  if (!"scientificName" %in% names(coefficients) && "species" %in% names(coefficients)) {
    data.table::setnames(coefficients, "species", "scientificName")
  }
  need_cols(
    coefficients,
    c("scientificName", coefficient_columns),
    "Stage 7.1 shared initialization coefficients"
  )
  coefficients[, scientificName := trimws(as.character(scientificName))]
  coefficients <- coefficients[scientificName %in% base$scientificName]
  assert(
    nrow(coefficients) == nrow(base) &&
      !anyDuplicated(coefficients$scientificName) &&
      setequal(coefficients$scientificName, base$scientificName),
    paste0(
      "The shared Stage 6 initialization must contain exactly one coefficient row for every retained species."
    )
  )
  coefficients <- coefficients[match(base$scientificName, scientificName)]

  out <- data.table::rbindlist(lapply(curves, function(curve) {
    curve_alpha <- suppressWarnings(as.numeric(coefficients[[paste0("alpha_", curve)]]))
    curve_beta <- suppressWarnings(as.numeric(coefficients[[paste0("beta_", curve)]]))
    assert(
      all(is.finite(curve_alpha) & curve_alpha > 0) &&
        all(is.finite(curve_beta) & curve_beta > 0),
      paste0(
        "The shared Stage 6 initialization contains missing, nonpositive, or nonfinite ",
        curve, " Gompertz coefficients for retained Stage 6 species."
      )
    )
    base[, .(
      scientificName, species, taxon, sdm_method, density,
      threshold_area_km2, parameter_order,
      curve = curve, alpha = curve_alpha, beta = curve_beta
    )]
  }), use.names = TRUE)
  out[, curve_order := match(curve, curves)]
  data.table::setorder(out, curve_order, parameter_order)

  configured <- out[curve == configured_curve][order(parameter_order)]
  # CSV and RDS serialization can preserve slightly different final digits
  # when Stage 4 is deterministically rebuilt. A tight relative tolerance
  # accepts that representation noise while still rejecting changed curves.
  tolerance_alpha <- 1e-8 * pmax(1, abs(base$alpha))
  tolerance_beta <- 1e-8 * pmax(1, abs(base$beta))
  assert(
    all(abs(configured$alpha - base$alpha) <= tolerance_alpha) &&
      all(abs(configured$beta - base$beta) <= tolerance_beta),
    paste0(
      "The selected ", configured_curve,
      " coefficients do not match the shared Stage 6 initialization."
    )
  )
  assert(
    nrow(out) == nrow(base) * length(curves) &&
      !anyDuplicated(out[, .(curve, species)]),
    "Stage 7.1 uncertainty parameters are incomplete or duplicated."
  )
  out[]
}

standardize_stage71_patch_lookup <- function(x, params, stage,
                                             label = "Stage 6 patch lookup") {
  dt <- data.table::as.data.table(x)
  need_cols(dt, c("species", "patch_id", "pu_id", "patch_area_km2"), label)
  if (!nrow(dt)) {
    return(data.table::data.table(
      scientificName = character(),
      species = character(),
      patch_id = integer(),
      pu_id = integer(),
      patch_area_km2 = numeric()
    ))
  }

  raw_species <- trimws(as.character(dt$species))
  if (all(raw_species %in% params$scientificName)) {
    scientific_name <- raw_species
  } else if (all(raw_species %in% params$species)) {
    scientific_name <- params$scientificName[match(raw_species, params$species)]
  } else {
    unknown <- unique(raw_species[
      !raw_species %in% params$scientificName & !raw_species %in% params$species
    ])
    stop(
      label, " for stage ", stage,
      " contains species outside the shared Stage 6 initialization: ",
      paste(utils::head(unknown, 10L), collapse = ", "),
      call. = FALSE
    )
  }

  patch_id <- suppressWarnings(as.numeric(dt$patch_id))
  pu_id <- suppressWarnings(as.numeric(dt$pu_id))
  area <- suppressWarnings(as.numeric(dt$patch_area_km2))
  assert(
    all(is.finite(patch_id) & patch_id > 0 & patch_id == floor(patch_id)) &&
      all(is.finite(pu_id) & pu_id > 0 & pu_id == floor(pu_id)),
    paste0(label, " for stage ", stage, " contains invalid patch_id or pu_id values.")
  )
  assert(
    all(is.finite(area) & area > 0),
    paste0(label, " for stage ", stage, " contains non-positive or non-finite patch areas.")
  )

  out <- data.table::data.table(
    scientificName = scientific_name,
    species = species_id(scientific_name),
    patch_id = as.integer(patch_id),
    pu_id = as.integer(pu_id),
    patch_area_km2 = area
  )
  duplicate_patch <- duplicated(out, by = c("species", "patch_id")) |
    duplicated(out, by = c("species", "patch_id"), fromLast = TRUE)
  if (any(duplicate_patch)) {
    preview <- unique(out[duplicate_patch, .(scientificName, patch_id, pu_id)])
    stop(
      label, " for stage ", stage,
      " contains duplicate species/patch rows. First duplicates: ",
      paste(
        paste0(utils::head(preview$scientificName, 5L), "|patch", utils::head(preview$patch_id, 5L)),
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  out[]
}

# Calculate one complete post-repair ecological state. PU rows are temporary;
# only the compact fixed-denominator species table is returned.
compute_stage71_stage_persistence <- function(lookup, params, stage_row,
                                              retain_pu_areas = FALSE) {
  stage <- as.integer(stage_row$stage[[1L]])
  standardized <- standardize_stage71_patch_lookup(
    lookup,
    params,
    stage,
    label = if (stage == 0L) {
      "Stage 6 shared-initialization patch_table"
    } else {
      "Stage 6 stage patch lookup"
    }
  )
  if (stage == 0L) {
    assert(
      setequal(unique(standardized$scientificName), params$scientificName),
      paste0(
        "Stage 6 shared-initialization patch_table does not exactly cover ",
        "the retained initialization species."
      )
    )
  }

  if (nrow(standardized)) {
    pu <- standardized[, .(
      pu_area_km2 = sum(patch_area_km2),
      n_patches = .N
    ), by = .(scientificName, species, pu_id)]
    pu_areas <- if (isTRUE(retain_pu_areas)) {
      pu[, .(scientificName, species, pu_id, pu_area_km2, n_patches)]
    } else {
      NULL
    }
    pu <- merge(
      pu,
      params[, .(scientificName, species, density, threshold_area_km2, alpha, beta)],
      by = c("scientificName", "species"),
      all.x = TRUE,
      sort = FALSE
    )
    assert(
      all(!is.na(pu$density)),
      paste0("Stage ", stage, " patch lookup could not be joined to all bundle parameters.")
    )
    invalid_threshold <- !is.finite(pu$pu_area_km2) |
      pu$pu_area_km2 <= pu$threshold_area_km2
    if (any(invalid_threshold)) {
      bad <- pu[invalid_threshold][seq_len(min(.N, 5L))]
      details <- paste0(
        "stage=", stage,
        " species=", bad$scientificName,
        " pu_id=", bad$pu_id,
        " area_km2=", format(bad$pu_area_km2, digits = 12, trim = TRUE),
        " threshold_km2=", format(bad$threshold_area_km2, digits = 12, trim = TRUE)
      )
      stop(
        "Represented PUs must remain strictly above their population-area threshold. First violations: ",
        paste(details, collapse = "; "),
        call. = FALSE
      )
    }
    pu[, P_pu := stage7_pu_persistence(
      area_km2 = pu_area_km2,
      density = density,
      alpha = alpha,
      beta = beta,
      threshold_area = threshold_area_km2,
      label = paste0("Stage 7.1 stage ", stage, " PU inputs")
    )]
    species_stage <- pu[, .(
      sp_persist = stage7_species_persistence(P_pu),
      n_pu = .N,
      total_pu_area_km2 = sum(pu_area_km2)
    ), by = .(scientificName, species)]
  } else {
    pu <- data.table::data.table()
    pu_areas <- if (isTRUE(retain_pu_areas)) {
      data.table::data.table(
        scientificName = character(), species = character(),
        pu_id = integer(), pu_area_km2 = numeric(), n_patches = integer()
      )
    } else {
      NULL
    }
    species_stage <- data.table::data.table(
      scientificName = character(), species = character(),
      sp_persist = numeric(), n_pu = integer(), total_pu_area_km2 = numeric()
    )
  }

  full <- merge(
    params[, .(scientificName, species, taxon, sdm_method, parameter_order)],
    species_stage,
    by = c("scientificName", "species"),
    all.x = TRUE,
    sort = FALSE
  )
  full[is.na(sp_persist), `:=`(
    sp_persist = 0,
    n_pu = 0L,
    total_pu_area_km2 = 0
  )]
  data.table::setorder(full, parameter_order)
  full[, `:=`(
    stage = stage,
    stage_order = as.integer(stage_row$stage_order[[1L]]),
    pct_cells_removed_end = as.numeric(stage_row$pct_cells_removed_end[[1L]]),
    pct_area_removed_end = as.numeric(stage_row$pct_area_removed_end[[1L]])
  )]
  data.table::setcolorder(full, c(
    "stage", "stage_order", "pct_cells_removed_end", "pct_area_removed_end",
    "scientificName", "species", "taxon", "sdm_method", "n_pu",
    "total_pu_area_km2", "sp_persist", "parameter_order"
  ))
  assert(
    nrow(full) == nrow(params) &&
      !anyDuplicated(full[, .(stage, species)]) &&
      all(is.finite(full$sp_persist) & full$sp_persist >= 0 & full$sp_persist <= 1),
    paste0("Stage 7.1 stage ", stage, " species persistence output is invalid.")
  )

  list(
    species = full[, setdiff(names(full), "parameter_order"), with = FALSE],
    lookup_rows = nrow(standardized),
    pu_count = nrow(pu),
    represented_species = sum(full$n_pu > 0L),
    pu_areas = pu_areas
  )
}

# Calculate compact mean trajectories for every Gompertz uncertainty curve from
# PU areas already aggregated for the primary Stage 7.1 calculation.
stage71_curve_means_from_pu <- function(pu_areas, curve_params, stage_row) {
  pu_areas <- data.table::as.data.table(pu_areas)
  curve_params <- data.table::as.data.table(curve_params)
  curves <- persistence_curves()
  stage <- as.integer(stage_row$stage[[1L]])
  taxa <- unique(curve_params$taxon)
  scopes <- if (length(taxa) == 1L) taxa else c("all", "mammals", "birds")

  pieces <- lapply(curves, function(curve_i) {
    par <- curve_params[curve == curve_i][order(parameter_order)]
    assert(
      nrow(par) > 0L && !anyDuplicated(par$species),
      paste0("Stage 7.1 curve parameters are invalid for ", curve_i, ".")
    )
    if (nrow(pu_areas)) {
      pu <- merge(
        pu_areas,
        par[, .(
          scientificName, species, density, threshold_area_km2, alpha, beta
        )],
        by = c("scientificName", "species"),
        all.x = TRUE,
        sort = FALSE
      )
      assert(
        all(!is.na(pu$density)),
        paste0("Stage ", stage, " PU areas could not be joined to ", curve_i, " parameters.")
      )
      pu[, P_pu := stage7_pu_persistence(
        area_km2 = pu_area_km2,
        density = density,
        alpha = alpha,
        beta = beta,
        threshold_area = threshold_area_km2,
        label = paste0("Stage 7.1 stage ", stage, " curve ", curve_i, " PU inputs")
      )]
      species_curve <- pu[, .(
        sp_persist = stage7_species_persistence(P_pu)
      ), by = .(scientificName, species)]
    } else {
      species_curve <- data.table::data.table(
        scientificName = character(), species = character(), sp_persist = numeric()
      )
    }
    full <- merge(
      par[, .(scientificName, species, taxon, parameter_order)],
      species_curve,
      by = c("scientificName", "species"),
      all.x = TRUE,
      sort = FALSE
    )
    full[is.na(sp_persist), sp_persist := 0]
    data.table::setorder(full, parameter_order)
    assert(
      nrow(full) == nrow(par) &&
        all(is.finite(full$sp_persist) & full$sp_persist >= 0 & full$sp_persist <= 1),
      paste0("Stage 7.1 stage ", stage, " curve ", curve_i, " species persistence is invalid.")
    )
    data.table::rbindlist(lapply(scopes, function(scope_i) {
      values <- if (scope_i == "all") full$sp_persist else full[taxon == scope_i, sp_persist]
      assert(length(values) > 0L, paste0("Stage 7.1 uncertainty scope is empty: ", scope_i))
      data.table::data.table(
        stage = stage,
        stage_order = as.integer(stage_row$stage_order[[1L]]),
        pct_cells_removed_end = as.numeric(stage_row$pct_cells_removed_end[[1L]]),
        scope = scope_i,
        curve = curve_i,
        n_species = length(values),
        mean_persist = mean(values)
      )
    }), use.names = TRUE)
  })
  out <- data.table::rbindlist(pieces, use.names = TRUE)
  out[, `:=`(
    scope_order__ = match(scope, c("all", "mammals", "birds")),
    curve_order__ = match(curve, curves)
  )]
  data.table::setorder(out, stage_order, scope_order__, curve_order__)
  out[, c("scope_order__", "curve_order__") := NULL]
  assert(
    nrow(out) == length(scopes) * length(curves) &&
      !anyDuplicated(out[, .(stage, scope, curve)]) &&
      all(is.finite(out$mean_persist) & out$mean_persist >= 0 & out$mean_persist <= 1),
    paste0("Stage 7.1 stage ", stage, " curve-mean summary is invalid.")
  )
  out[]
}

stage71_summarize_persistence <- function(species) {
  dt <- data.table::as.data.table(species)
  scope_summary <- function(x, scope_value) {
    x[, .(
      scope = scope_value,
      n_species = .N,
      n_species_with_pu = sum(n_pu > 0L),
      mean_persist = mean(sp_persist),
      median_persist = stats::median(sp_persist),
      q10_persist = as.numeric(stats::quantile(sp_persist, 0.10, names = FALSE, type = 7)),
      q25_persist = as.numeric(stats::quantile(sp_persist, 0.25, names = FALSE, type = 7)),
      q75_persist = as.numeric(stats::quantile(sp_persist, 0.75, names = FALSE, type = 7)),
      q90_persist = as.numeric(stats::quantile(sp_persist, 0.90, names = FALSE, type = 7)),
      min_persist = min(sp_persist),
      max_persist = max(sp_persist),
      frac_gt_05 = mean(sp_persist > 0.5),
      total_pu_area_km2 = sum(total_pu_area_km2)
    ), by = .(stage, stage_order, pct_cells_removed_end, pct_area_removed_end)]
  }

  taxa <- unique(dt$taxon)
  if (length(taxa) == 1L) {
    out <- scope_summary(dt, taxa[[1L]])
  } else {
    pieces <- list(scope_summary(dt, "all"))
    for (taxon_name in c("mammals", "birds")) {
      pieces[[length(pieces) + 1L]] <- scope_summary(dt[taxon == taxon_name], taxon_name)
    }
    out <- data.table::rbindlist(pieces, use.names = TRUE)
  }
  out[, scope_order__ := match(scope, c("all", "mammals", "birds"))]
  data.table::setorder(out, stage_order, scope_order__)
  out[, scope_order__ := NULL]
  assert(
      all(is.finite(out$mean_persist) & out$mean_persist >= 0 & out$mean_persist <= 1) &&
      all(is.finite(out$median_persist) & out$median_persist >= 0 & out$median_persist <= 1) &&
      all(is.finite(out$q10_persist) & out$q10_persist >= 0 & out$q10_persist <= 1) &&
      all(is.finite(out$q25_persist) & out$q25_persist >= 0 & out$q25_persist <= 1) &&
      all(is.finite(out$q75_persist) & out$q75_persist >= 0 & out$q75_persist <= 1) &&
      all(is.finite(out$q90_persist) & out$q90_persist >= 0 & out$q90_persist <= 1) &&
      all(
        out$q10_persist <= out$q25_persist &
          out$q25_persist <= out$median_persist &
          out$median_persist <= out$q75_persist &
          out$q75_persist <= out$q90_persist
      ) &&
      all(is.finite(out$frac_gt_05) & out$frac_gt_05 >= 0 & out$frac_gt_05 <= 1),
    "Stage 7.1 community persistence summary is invalid."
  )
  out[]
}

stage71_persistence_changes <- function(species) {
  dt <- data.table::as.data.table(species)
  first_stage <- min(dt$stage_order)
  last_stage <- max(dt$stage_order)
  auc <- dt[
    order(stage_order),
    {
      x <- as.numeric(pct_cells_removed_end)
      y <- as.numeric(sp_persist)
      span <- max(x) - min(x)
      area <- if (length(x) > 1L) {
        sum(diff(x) * (utils::head(y, -1L) + utils::tail(y, -1L)) / 2)
      } else {
        0
      }
      list(
        persistence_auc = area,
        normalized_persistence_auc = if (span > 0) area / span else y[[1L]],
        trajectory_cell_removal_span = span
      )
    },
    by = .(scientificName, species)
  ]
  initial <- dt[stage_order == first_stage, .(
    scientificName, species, taxon, sdm_method,
    initial_persistence = sp_persist
  )]
  latest <- dt[stage_order == last_stage, .(
    scientificName, species,
    latest_persistence = sp_persist,
    latest_n_pu = n_pu,
    latest_pu_area_km2 = total_pu_area_km2
  )]
  out <- merge(initial, latest, by = c("scientificName", "species"), sort = FALSE)
  out <- merge(out, auc, by = c("scientificName", "species"), sort = FALSE)
  out[, persistence_loss := initial_persistence - latest_persistence]
  assert(
    all(is.finite(out$persistence_auc) & out$persistence_auc >= 0) &&
      all(is.finite(out$normalized_persistence_auc) &
        out$normalized_persistence_auc >= 0 & out$normalized_persistence_auc <= 1),
    "Stage 7.1 species persistence AUC values are invalid."
  )
  data.table::setorder(out, -persistence_loss, scientificName)
  out[]
}

# Read one ecological state at a time and retain only species-level trajectories.
compute_stage71_pipeline_persistence <- function(stage_meta, params, lookup_loader,
                                                 progress = NULL,
                                                 curve_params = NULL) {
  meta <- data.table::as.data.table(stage_meta)
  need_cols(
    meta,
    c("stage", "patch_lookup_path", "pct_cells_removed_end", "pct_area_removed_end"),
    "Stage 7.1 stage metadata"
  )
  meta <- data.table::copy(meta)
  data.table::setorder(meta, stage)
  meta[, stage_order := seq_len(.N) - 1L]
  assert(
    identical(as.integer(meta$stage), seq.int(0L, nrow(meta) - 1L)),
    "Stage 7.1 stages must be contiguous from zero."
  )

  species_by_stage <- vector("list", nrow(meta))
  curve_summary_by_stage <- if (is.null(curve_params)) NULL else vector("list", nrow(meta))
  for (i in seq_len(nrow(meta))) {
    started <- proc.time()[["elapsed"]]
    stage <- as.integer(meta$stage[[i]])
    path <- as.character(meta$patch_lookup_path[[i]])
    lookup <- lookup_loader(stage, path)
    result <- compute_stage71_stage_persistence(
      lookup,
      params,
      meta[i],
      retain_pu_areas = !is.null(curve_params)
    )
    species_by_stage[[i]] <- result$species
    if (!is.null(curve_params)) {
      curve_summary_by_stage[[i]] <- stage71_curve_means_from_pu(
        result$pu_areas,
        curve_params,
        meta[i]
      )
    }
    if (is.function(progress)) {
      progress(list(
        stage = stage,
        lookup_rows = result$lookup_rows,
        pu_count = result$pu_count,
        represented_species = result$represented_species,
        mean_persist = mean(result$species$sp_persist),
        median_persist = stats::median(result$species$sp_persist),
        elapsed_seconds = proc.time()[["elapsed"]] - started
      ))
    }
    rm(lookup, result)
  }

  species <- data.table::rbindlist(species_by_stage, use.names = TRUE)
  expected_rows <- nrow(meta) * nrow(params)
  assert(
    nrow(species) == expected_rows &&
      !anyDuplicated(species[, .(stage, species)]),
    "Stage 7.1 persistence trajectory does not contain exactly one row per retained species and stage."
  )
  curve_summary <- if (is.null(curve_summary_by_stage)) {
    NULL
  } else {
    out <- data.table::rbindlist(curve_summary_by_stage, use.names = TRUE)
    expected_scopes <- if (data.table::uniqueN(params$taxon) == 1L) 1L else 3L
    assert(
      nrow(out) == nrow(meta) * expected_scopes * length(persistence_curves()) &&
        !anyDuplicated(out[, .(stage, scope, curve)]),
      "Stage 7.1 uncertainty curve summary is incomplete or duplicated."
    )
    out
  }
  list(
    species = species,
    summary = stage71_summarize_persistence(species),
    changes = stage71_persistence_changes(species),
    curve_summary = curve_summary
  )
}
