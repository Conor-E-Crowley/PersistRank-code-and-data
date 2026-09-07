# Stage 7.3 figure-data preparation.
#
# The Stage 7.3 workflow owns this definition-only module. It validates and
# prepares shared focal and assemblage data once without building plot copies,
# accessing artifacts, attaching packages, or writing files.

# Shared stage selection and uncertainty helpers ---------------------------

stage73_x_percent <- function(x) {
  values <- as.numeric(x)
  if (length(values) && max(values, na.rm = TRUE) <= 1.0001) {
    100 * values
  } else {
    values
  }
}

stage73_pu_bar_style <- function() {
  list(
    segment_border_colour = c(
      pipe = "#789CB5",
      rank = "#C98965"
    ),
    segment_border_linewidth = 0.26,
    outer_border_colour = c(
      pipe = "#2E78AF",
      rank = "#B7561D"
    ),
    outer_border_linewidth = 0.52
  )
}

stage73_stage_key <- function(stage_summary, label) {
  stage_summary <- data.table::as.data.table(data.table::copy(stage_summary))
  need_cols(
    stage_summary,
    c("stage", "stage_order", "x"),
    label
  )

  key <- unique(data.table::copy(stage_summary)[, .(
    stage = as.integer(stage),
    stage_order = as.integer(stage_order),
    x = as.numeric(x)
  )])
  data.table::setorder(key, stage_order)
  assert(
    nrow(key) > 0L &&
      !anyDuplicated(key$stage) &&
      !anyDuplicated(key$stage_order) &&
      all(is.finite(key$x)) &&
      all(diff(key$x) >= 0),
    paste0(label, " must contain a unique, ordered, nondecreasing stage sequence.")
  )
  key[, x_percent := stage73_x_percent(x)]
  key[]
}

select_stage73_panel_c_stages <- function(
  stage_summary,
  target_percentages = seq(50, 100, by = 5)
) {
  targets <- suppressWarnings(as.numeric(target_percentages))
  assert(
    length(targets) > 0L &&
      all(is.finite(targets) & targets >= 0 & targets <= 100),
    "Stage 7.3 retained-area targets must be percentages in [0,100]."
  )
  key <- stage73_stage_key(
    stage_summary,
    "Stage 7.3 retained-area stage input"
  )
  positions <- vapply(
    targets,
    function(target_i) which.min(abs(key$x_percent - target_i)),
    integer(1L)
  )
  selected <- key[sort(unique(positions))]
  assert(
    nrow(selected) > 0L,
    "Stage 7.3 retained-area stage selection returned no stages."
  )
  selected[]
}

stage73_focal_spec <- function() {
  list(
    figure_id = "persistence_comparison_focal_late_stage",
    minimum_removed_percent = 50,
    target_percentages = seq(50, 95, by = 5),
    persistence_area_reference_removed_percent = 50
  )
}

stage73_uncertainty_styles <- function() {
  styles <- data.table::data.table(
    interval = c("95% range", "68% range"),
    interval_label = c("95% (q2.5–q97.5)", "68% (q16–q84)"),
    lower_curve = c("q025", "q16"),
    upper_curve = c("q975", "q84"),
    alpha = c(0.10, 0.22)
  )
  assert(
    !anyNA(styles) &&
      !anyDuplicated(styles$interval) &&
      all(styles$lower_curve %in% persistence_curves()) &&
      all(styles$upper_curve %in% persistence_curves()),
    "Stage 7.3 uncertainty styles are incomplete."
  )
  styles
}

stage73_uncertainty_bands <- function(data, value_col, id_cols) {
  data <- data.table::as.data.table(data.table::copy(data))
  value_col <- validate_scalar_string(value_col, "value_col")
  id_cols <- as.character(id_cols)
  need_cols(
    data,
    c(id_cols, "curve", value_col),
    "Stage 7.3 uncertainty trajectory"
  )
  assert(
    length(id_cols) > 0L &&
      !anyDuplicated(data[, c(id_cols, "curve"), with = FALSE]) &&
      setequal(unique(as.character(data$curve)), persistence_curves()),
    "Stage 7.3 uncertainty trajectory must contain one row per identifier and curve."
  )

  formula <- stats::as.formula(
    paste(paste(id_cols, collapse = " + "), "~ curve")
  )
  wide <- data.table::dcast(
    data[, c(id_cols, "curve", value_col), with = FALSE],
    formula,
    value.var = value_col
  )
  curve_cols <- persistence_curves()
  need_cols(wide, curve_cols, "Stage 7.3 uncertainty trajectory")
  curve_values <- as.matrix(wide[, ..curve_cols])
  assert(
    all(vapply(wide[, ..curve_cols], is.numeric, logical(1L))) &&
      all(is.finite(curve_values)),
    "Stage 7.3 uncertainty trajectory contains invalid persistence values."
  )
  wide[, `:=`(
    central = q50,
    outer_ymin = pmin(q025, q975),
    outer_ymax = pmax(q025, q975),
    inner_ymin = pmin(q16, q84),
    inner_ymax = pmax(q16, q84)
  )]

  styles <- stage73_uncertainty_styles()
  outer <- wide[, ..id_cols]
  outer[, `:=`(
    ymin = wide$outer_ymin,
    ymax = wide$outer_ymax,
    interval = styles$interval[[1L]]
  )]
  inner <- wide[, ..id_cols]
  inner[, `:=`(
    ymin = wide$inner_ymin,
    ymax = wide$inner_ymax,
    interval = styles$interval[[2L]]
  )]
  bands <- data.table::rbindlist(list(outer, inner), use.names = TRUE)
  bands[, interval := factor(
    interval,
    levels = rev(styles$interval),
    labels = rev(styles$interval_label)
  )]
  list(central = wide, bands = bands)
}

# Figure data shared by every Stage 7.3 plot -------------------------------

prepare_stage73_figure_data <- function(core) {
  assert(is.list(core), "core must be the result returned by run_stage73_comparison().")
  required <- c(
    "config", "pu_long", "sp_long", "stage_sum", "curve_parameters"
  )
  assert(
    all(required %in% names(core)),
    paste0("Stage 7.3 figure preparation requires: ", paste(required, collapse = ", "), ".")
  )
  configured_curve <- validate_persistence_curve(
    core$config$curve,
    "Stage 7.3 configured curve"
  )
  curves <- persistence_curves()
  curve_parameters <- data.table::copy(core$curve_parameters)
  need_cols(
    curve_parameters,
    c("species", "density", "c_th", "curve", "alpha", "beta"),
    "Stage 7.3 figure curve parameters"
  )
  curve_parameters[, curve := as.character(curve)]
  assert(
    setequal(curve_parameters$curve, curves) &&
      !anyDuplicated(curve_parameters[, .(species, curve)]) &&
      all(
        is.finite(curve_parameters$density) & curve_parameters$density > 0 &
          is.finite(curve_parameters$c_th) & curve_parameters$c_th > 0 &
          is.finite(curve_parameters$alpha) & curve_parameters$alpha > 0 &
          is.finite(curve_parameters$beta) & curve_parameters$beta > 0
      ),
    "Stage 7.3 figure curve parameters are incomplete or invalid."
  )

  species_template <- data.table::copy(core$sp_long)[, .(
    method = as.character(method),
    stage = as.integer(stage),
    stage_order = as.integer(stage_order),
    x = as.numeric(x),
    scientificName = as.character(scientificName),
    species = as.character(species),
    className = as.character(className),
    sp_persist = as.numeric(sp_persist),
    n_pu = as.integer(n_pu),
    total_pu_area_km2 = as.numeric(total_pu_area_km2)
  )]
  species_template[, x_percent := stage73_x_percent(x)]
  assert(
    !anyDuplicated(species_template[, .(method, stage, species)]) &&
      all(
        is.finite(species_template$sp_persist) &
          species_template$sp_persist >= 0 &
          species_template$sp_persist <= 1
      ),
    "Stage 7.3 configured-curve species trajectories are invalid."
  )

  curve_species <- vector("list", length(curves))
  names(curve_species) <- curves
  configured <- data.table::copy(species_template)
  configured[, curve := configured_curve]
  curve_species[[configured_curve]] <- configured

  remaining_curves <- setdiff(curves, configured_curve)
  if (length(remaining_curves)) {
    pu_work <- core$pu_long[, .(
      method = as.character(method),
      stage = as.integer(stage),
      stage_order = as.integer(stage_order),
      x = as.numeric(x),
      species = as.character(species),
      pu_area_km2 = as.numeric(pu_area_km2)
    )]
    assert(
      all(is.finite(pu_work$pu_area_km2) & pu_work$pu_area_km2 > 0),
      "Stage 7.3 PU areas are invalid for uncertainty reconstruction."
    )

    template_columns <- setdiff(names(species_template), "sp_persist")
    full_template <- species_template[, ..template_columns]

    for (curve_i in remaining_curves) {
      params_i <- curve_parameters[curve == curve_i]
      parameter_row <- match(pu_work$species, params_i$species)
      assert(
        all(!is.na(parameter_row)),
        paste0("Missing figure parameters for curve ", curve_i, ".")
      )
      pu_work[, curve_pu_persistence := stage7_pu_persistence(
        area_km2 = pu_area_km2,
        density = params_i$density[parameter_row],
        alpha = params_i$alpha[parameter_row],
        beta = params_i$beta[parameter_row],
        threshold_area = params_i$c_th[parameter_row],
        label = paste0("Stage 7.3 figure curve ", curve_i)
      )]
      observed <- pu_work[, .(
        sp_persist = stage7_species_persistence(curve_pu_persistence)
      ), by = .(method, stage, stage_order, x, species)]
      assert(
        all(
          is.finite(observed$sp_persist) &
            observed$sp_persist >= 0 &
            observed$sp_persist <= 1
        ),
        paste0(
          "Stage 7.3 figure curve ", curve_i,
          " species persistence calculation returned invalid values."
        )
      )
      reconstructed <- merge(
        full_template,
        observed,
        by = c("method", "stage", "stage_order", "x", "species"),
        all.x = TRUE,
        sort = FALSE
      )
      reconstructed[is.na(sp_persist), sp_persist := 0]
      reconstructed[, curve := curve_i]
      data.table::setcolorder(reconstructed, c(names(species_template), "curve"))
      curve_species[[curve_i]] <- reconstructed
      pu_work[, curve_pu_persistence := NULL]
      rm(observed, reconstructed, parameter_row, params_i)
    }
    rm(pu_work, full_template)
  }

  species_curves <- data.table::rbindlist(
    curve_species[curves],
    use.names = TRUE
  )
  data.table::setorder(species_curves, method, curve, species, stage_order)
  expected_rows <- nrow(species_template) * length(curves)
  assert(
    nrow(species_curves) == expected_rows &&
      !anyDuplicated(species_curves[, .(method, curve, stage, species)]) &&
      all(
        is.finite(species_curves$sp_persist) &
          species_curves$sp_persist >= 0 &
          species_curves$sp_persist <= 1
      ),
    "Stage 7.3 figure trajectories are incomplete or invalid."
  )
  configured_check <- species_curves[curve == configured_curve, sp_persist]
  expected_configured <- species_template[
    order(method, species, stage_order),
    sp_persist
  ]
  assert(
    identical(configured_check, expected_configured),
    "Stage 7.3 configured-curve figure trajectory changed during uncertainty preparation."
  )

  stage_curves <- species_curves[, .(
    n_species = .N,
    mean_persist = mean(sp_persist)
  ), by = .(method, curve, stage, stage_order, x, x_percent)]
  styles <- stage73_uncertainty_styles()

  log_msg(
    "STAGE 7.3 FIGURE DATA | curves=", length(curves),
    " species_rows=", nrow(species_curves),
    " focal_grid=50_to_100"
  )
  list(
    species = species_curves,
    stage = stage_curves,
    uncertainty_styles = styles,
    configured_curve = configured_curve
  )
}

# Focal-species panels used by the combined full-page figure ----------------

# Point styles are a shared report configuration choice: "hollow"
# (white-filled ring), "solid" (plain filled point), or "none" (no point
# layer at all).

# Hollow, white-filled points (shape 21) and solid filled points (shape 19)
# are the two point-bearing symbologies; "none" is handled by callers, which
# skip adding a point layer entirely rather than calling this helper.
prepare_stage73_distribution_data <- function(
  core,
  figure_data = NULL,
  central_stat = "mean"
) {
  central_stat <- validate_scalar_choice(
    central_stat, report_central_statistics(), "Stage 7.3 central_stat"
  )
  assert(is.list(core), "core must be the result returned by run_stage73_comparison().")
  if (is.null(figure_data)) {
    figure_data <- prepare_stage73_figure_data(core)
  }
  assert(
    is.list(figure_data) && data.table::is.data.table(figure_data$species),
    "Stage 7.3 species-distribution figure data are unavailable."
  )
  required <- c(
    "method", "stage", "stage_order", "x_percent",
    "scientificName", "species", "sp_persist", "curve"
  )
  need_cols(
    figure_data$species,
    required,
    "Stage 7.3 species-persistence distribution trajectories"
  )
  species_source <- data.table::copy(figure_data$species)
  species <- species_source[curve == "q50", ..required]
  rm(species_source)
  species[, `:=`(
    method = as.character(method),
    stage = as.integer(stage),
    stage_order = as.integer(stage_order),
    x_percent = as.numeric(x_percent),
    scientificName = as.character(scientificName),
    species = as.character(species),
    sp_persist = as.numeric(sp_persist),
    curve = as.character(curve)
  )]
  assert(
    nrow(species) > 0L &&
      identical(sort(unique(species$method)), c("pipe", "rank")) &&
      !anyDuplicated(species[, .(method, stage, species)]) &&
      all(!is.na(species$stage) & !is.na(species$stage_order)) &&
      all(is.finite(species$x_percent) & species$x_percent >= 0 & species$x_percent <= 100) &&
      all(is.finite(species$sp_persist) &
            species$sp_persist >= 0 & species$sp_persist <= 1),
    "Stage 7.3 q50 species-persistence distributions are incomplete or invalid."
  )
  identity_key <- unique(species[, .(species, scientificName)])
  assert(
    !anyDuplicated(identity_key$species) &&
      !anyDuplicated(identity_key$scientificName) &&
      all(nzchar(identity_key$species)) && all(nzchar(identity_key$scientificName)),
    "Stage 7.3 species-distribution identifiers are inconsistent or blank."
  )
  stage_key <- unique(species[, .(method, stage, stage_order, x_percent)])
  reference_stages <- stage_key[
    method == "pipe", .(stage, stage_order, x_percent)
  ][order(stage_order)]
  stage_signatures <- stage_key[order(stage_order), .(
    signature = paste(stage, stage_order, sprintf("%.12f", x_percent), collapse = "\r")
  ), by = method]
  coverage <- species[, .(
    n_species = data.table::uniqueN(species),
    signature = paste(sort(unique(species)), collapse = "\r")
  ), by = .(method, stage)]
  assert(
    nrow(reference_stages) > 0L &&
      !anyDuplicated(stage_key[, .(method, stage)]) &&
      data.table::uniqueN(stage_signatures$signature) == 1L &&
      all(coverage$n_species == nrow(identity_key)) &&
      data.table::uniqueN(coverage$signature) == 1L,
    "Stage 7.3 species-distribution figure requires identical stages and one fixed species denominator for both methods."
  )
  counts <- list(all = as.integer(nrow(identity_key)))

  distribution <- species[, {
    quantiles <- stats::quantile(
      sp_persist,
      probs = c(0.10, 0.25, 0.50, 0.75, 0.90),
      names = FALSE,
      type = 7
    )
    list(
      n_species = data.table::uniqueN(species),
      mean_persist = mean(sp_persist),
      q10 = quantiles[[1L]],
      q25 = quantiles[[2L]],
      median_persist = quantiles[[3L]],
      q75 = quantiles[[4L]],
      q90 = quantiles[[5L]]
    )
  }, by = .(method, stage, stage_order, x_percent)]
  distribution[, central := if (identical(central_stat, "median")) {
    median_persist
  } else {
    mean_persist
  }]
  assert(
    !anyDuplicated(distribution[, .(method, stage)]) &&
      all(distribution$n_species == counts$all) &&
      all(
        is.finite(distribution$mean_persist) &
          distribution$q10 <= distribution$q25 &
          distribution$q25 <= distribution$median_persist &
          distribution$median_persist <= distribution$q75 &
          distribution$q75 <= distribution$q90 &
          distribution$q10 >= 0 & distribution$q90 <= 1
      ),
    "Stage 7.3 species-persistence distribution summaries are invalid."
  )

  data.table::setorder(distribution, method, stage_order)
  list(
    species = species,
    central = distribution,
    counts = counts,
    stages = reference_stages
  )
}
