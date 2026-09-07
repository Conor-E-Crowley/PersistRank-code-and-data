# Stage 7.3 scientific comparison statistics.
#
# Calculations retain the complete species denominator, use strict threshold
# comparisons, and preserve canonical stage and species order.
# Removal-sequence integration uses the percentage-of-cells-removed axis;
# retained-area integration uses km2 over the methods' shared support.

# The returned groups support the main-text and supplementary summaries:
#   1. late-stage mean persistence and checkpoint threshold counts
#   2. species-level signed AUC
#   3. area-normalized species persistence advantage
#   4. interpolated threshold-crossing asymmetry
#   5. PU redundancy-loss counts
#   6. demographic-curve sensitivity of late-stage contrasts and
#      species-level signed AUC across the complete removal sequence

compute_stage73_statistics <- function(
    core,
    species_name,
    figure_data = NULL
) {
  load_packages("data.table")
  stage_sum <- core$stage_sum
  sp_long <- core$sp_long
  pu_long <- core$pu_long
  rank_label <- core$rank_label
  persistence_threshold <- core$config$persistence_threshold
  species_name <- validate_scalar_string(
    species_name,
    "species_name"
  )

assert(
  length(persistence_threshold) == 1L &&
    is.numeric(persistence_threshold) &&
    is.finite(persistence_threshold) &&
    persistence_threshold > 0 &&
    persistence_threshold < 1,
  paste0(
    "core$config$persistence_threshold must be one ",
    "finite numeric value strictly between 0 and 1."
  )
)

# Only directional classification uses this tolerance; persistence values,
# thresholds, axes, and integrated differences are otherwise left unchanged.
metric_eps <- 1e-9

metric_median <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  median(x, na.rm = TRUE)
}

metric_quantile <- function(x, p) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(quantile(x, probs = p, na.rm = TRUE, names = FALSE))
}

trapz_xy <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]

  if (length(x) < 2L) return(NA_real_)

  ord <- order(x)
  x <- x[ord]
  y <- y[ord]

  dd <- data.table(x = x, y = y)
  dd <- dd[, .(y = mean(y, na.rm = TRUE)), by = x]
  setorder(dd, x)

  if (nrow(dd) < 2L) return(NA_real_)

  sum(diff(dd$x) * (head(dd$y, -1L) + tail(dd$y, -1L)) / 2)
}

fraction_axis_above_zero <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  dd <- data.table(x = as.numeric(x[ok]), y = as.numeric(y[ok]))
  if (nrow(dd) < 2L) return(NA_real_)
  dd <- dd[, .(y = mean(y)), by = x]
  setorder(dd, x)
  if (nrow(dd) < 2L || max(dd$x) <= min(dd$x)) return(NA_real_)

  x1 <- head(dd$x, -1L)
  x2 <- tail(dd$x, -1L)
  y1 <- head(dd$y, -1L)
  y2 <- tail(dd$y, -1L)
  widths <- x2 - x1
  positive_widths <- numeric(length(widths))
  both_positive <- y1 > 0 & y2 > 0
  positive_widths[both_positive] <- widths[both_positive]
  crossing <- (y1 > 0 & y2 <= 0) | (y1 <= 0 & y2 > 0)
  if (any(crossing)) {
    zero_x <- x1[crossing] - y1[crossing] * widths[crossing] /
      (y2[crossing] - y1[crossing])
    starts_positive <- y1[crossing] > 0
    positive_widths[crossing] <- ifelse(
      starts_positive,
      zero_x - x1[crossing],
      x2[crossing] - zero_x
    )
  }
  sum(pmax(0, pmin(widths, positive_widths))) / sum(widths)
}

taxon_scope <- function(class_name) {
  class_name <- trimws(as.character(class_name))
  fcase(
    tolower(class_name) == "mammalia", "Mammals",
    tolower(class_name) == "aves", "Birds",
    default = class_name
  )
}

x_to_pct <- function(x) {
  if (max(stage_sum$x, na.rm = TRUE) <= 1.0001) 100 * x else x
}

target_to_x <- function(pct) {
  if (max(stage_sum$x, na.rm = TRUE) <= 1.0001) pct / 100 else pct
}

resolve_species_id <- function(scientific_name, dt = sp_long) {
  candidate <- species_id(scientific_name)

  if (candidate %in% unique(dt$species)) {
    return(candidate)
  }

  hit <- unique(dt[
    tolower(trimws(scientificName)) == tolower(trimws(scientific_name)),
    species
  ])

  if (length(hit)) return(hit[1L])

  NA_character_
}

focal_name_results <- species_name
focal_id_results <- resolve_species_id(focal_name_results, sp_long)

# ---- Late-stage mean persistence and loss reduction -------------------------

# Every summary retains the complete species denominator from `sp_long`; a
# species is never removed merely because its persistence reaches zero.

target_pcts <- c(90, 95, 99)

stage_key <- unique(stage_sum[, .(stage, stage_order, x)])
setorder(stage_key, stage_order)
assert(nrow(stage_key) > 0L,
       "stage_sum contains no stages for checkpoint statistics.")

all_stage_wide <- dcast(
  stage_sum,
  stage + stage_order + x ~ method,
  value.var = c(
    "mean_persist",
    "median_persist",
    "frac_gt_05",
    "total_pu_area_km2"
  )
)
all_stage_wide[, `:=`(
  actual_pct_removed = x_to_pct(x),
  persistence_loss_pipe = 1 - mean_persist_pipe,
  persistence_loss_rank = 1 - mean_persist_rank,
  mean_persist_diff_pipe_minus_rank = mean_persist_pipe - mean_persist_rank,
  median_persist_diff_pipe_minus_rank = median_persist_pipe - median_persist_rank,
  persistence_loss_reduction_abs =
    (1 - mean_persist_rank) - (1 - mean_persist_pipe),
  persistence_loss_reduction_pct_vs_rank = fifelse(
    (1 - mean_persist_rank) > metric_eps,
    100 * ((1 - mean_persist_rank) - (1 - mean_persist_pipe)) /
      (1 - mean_persist_rank),
    NA_real_
  ),
  frac_gt_05_diff_pipe_minus_rank = frac_gt_05_pipe - frac_gt_05_rank,
  retained_area_diff_pipe_minus_rank =
    total_pu_area_km2_pipe - total_pu_area_km2_rank
)]
setorder(all_stage_wide, stage_order)

# Requested percentages are labels on the removal-sequence axis. The nearest
# canonical stage is selected deterministically without interpolating outcomes.
nearest_stage_checkpoints <- function(target_pcts, stage_key, label) {
  targets <- suppressWarnings(as.numeric(target_pcts))
  assert(
    length(targets) > 0L &&
      all(is.finite(targets)) &&
      all(targets >= 0) &&
      all(targets <= 100),
    paste0(label, " targets must be percentages in [0, 100].")
  )

  out <- rbindlist(lapply(targets, function(pct_i) {
    target_x <- target_to_x(pct_i)
    ii <- which.min(abs(stage_key$x - target_x))

    row <- copy(stage_key[ii])
    row[, target_pct_removed := pct_i]
    row[, actual_pct_removed := x_to_pct(x)]
    row
  }), use.names = TRUE, fill = TRUE)

  shared_stage_targets <- out[
    ,
    .(n_targets = uniqueN(target_pct_removed)),
    by = .(stage, stage_order, x)
  ][n_targets > 1L]

  if (nrow(shared_stage_targets)) {
    message(
      "[", label, "] Multiple requested removal targets map to the same ",
      "nearest stage; retaining every target label for transparent reporting."
    )
  }

  out[]
}

checkpoint_key <- nearest_stage_checkpoints(
  target_pcts = target_pcts,
  stage_key = stage_key,
  label = "late-stage persistence checkpoints"
)

checkpoint_stage_sum <- merge(
  stage_sum,
  checkpoint_key[
    ,
    .(
      stage,
      stage_order,
      x,
      target_pct_removed,
      actual_pct_removed
    )
  ],
  by = c("stage", "stage_order", "x"),
  all.x = FALSE,
  sort = FALSE
)

checkpoint_stage_sum[
  ,
  persistence_loss := 1 - mean_persist
]

checkpoint_stage_wide <- dcast(
  checkpoint_stage_sum,
  target_pct_removed + actual_pct_removed + stage + stage_order + x ~ method,
  value.var = c(
    "mean_persist",
    "persistence_loss",
    "median_persist",
    "frac_gt_05",
    "total_pu_area_km2"
  )
)

checkpoint_stage_wide[
  ,
  `:=`(
    mean_persist_diff_pipe_minus_rank =
      mean_persist_pipe - mean_persist_rank,
    persistence_loss_reduction_abs =
      persistence_loss_rank - persistence_loss_pipe,
    persistence_loss_reduction_pct_vs_rank =
      fifelse(
        persistence_loss_rank > metric_eps,
        100 * (persistence_loss_rank - persistence_loss_pipe) /
          persistence_loss_rank,
        NA_real_
      )
  )
]

setorder(checkpoint_stage_wide, target_pct_removed)

checkpoint_species <- merge(
  sp_long,
  checkpoint_key[, .(
    stage,
    stage_order,
    x,
    target_pct_removed,
    actual_pct_removed
  )],
  by = c("stage", "stage_order", "x"),
  all.x = FALSE,
  sort = FALSE
)
checkpoint_species_wide <- dcast(
  checkpoint_species,
  target_pct_removed + actual_pct_removed + stage + stage_order + x +
    species + scientificName + className ~ method,
  value.var = "sp_persist"
)
assert(
  all(is.finite(checkpoint_species_wide$pipe)) &&
    all(is.finite(checkpoint_species_wide$rank)),
  "Stage 7.3 checkpoint species persistence is incomplete."
)
checkpoint_species_wide[, taxon := taxon_scope(className)]
checkpoint_scoped <- rbindlist(list(
  copy(checkpoint_species_wide)[, scope := "All species"],
  copy(checkpoint_species_wide)[, scope := taxon]
), use.names = TRUE)
# Reporting focuses on the single configured manuscript threshold; other
# persistence levels are covered separately by the continuous
# threshold-crossing analysis below.
checkpoint_thresholds <- persistence_threshold
checkpoint_threshold_summary <- rbindlist(lapply(checkpoint_thresholds, function(thr) {
  checkpoint_scoped[, .(
    threshold = thr,
    n_species = .N,
    n_above_pipe = sum(pipe > thr),
    n_above_rank = sum(rank > thr),
    net_species_above_pipe_minus_rank = sum(pipe > thr) - sum(rank > thr),
    n_pipe_only_above = sum(pipe > thr & rank <= thr),
    n_rank_only_above = sum(rank > thr & pipe <= thr),

    n_below_pipe = sum(pipe < thr),
    n_below_rank = sum(rank < thr),
    net_species_below_rank_minus_pipe =
      sum(rank < thr) - sum(pipe < thr),

    n_pipe_only_below =
      sum(pipe < thr & rank >= thr),

    n_rank_only_below =
      sum(rank < thr & pipe >= thr),

    n_both_below =
      sum(pipe < thr & rank < thr),

    n_neither_below =
      sum(pipe >= thr & rank >= thr),

    mean_species_gap = mean(pipe - rank),
    median_species_gap = metric_median(pipe - rank)
  ), by = .(
    target_pct_removed,
    actual_pct_removed,
    stage,
    stage_order,
    x,
    scope
  )]
}), use.names = TRUE)
setorder(checkpoint_threshold_summary, target_pct_removed, scope, -threshold)

assert(
  all(checkpoint_threshold_summary[
    ,
    n_species ==
      n_pipe_only_below +
      n_rank_only_below +
      n_both_below +
      n_neither_below
  ]),
  "Checkpoint below-threshold categories do not partition the species set."
)

# ---- Species-level signed AUC ------------------------------------------------

# Integrate method differences along the canonical cell-removal sequence (`x`),
# not by stage number; both methods therefore share identical sequence support.

auc_dt <- dcast(
  sp_long[
    ,
    .(
      species,
      scientificName,
      className,
      stage,
      stage_order,
      x,
      method,
      sp_persist
    )
  ],
  species + scientificName + className + stage + stage_order + x ~ method,
  value.var = "sp_persist"
)

auc_dt <- auc_dt[is.finite(pipe) & is.finite(rank)]
auc_dt[, gap_pipe_minus_rank := pipe - rank]
setorder(auc_dt, species, stage_order)

species_auc <- auc_dt[
  ,
  .(
    n_stages = .N,
    auc_diff_pipe_minus_rank =
      trapz_xy(x, gap_pipe_minus_rank),
    mean_gap_over_sequence =
      trapz_xy(x, gap_pipe_minus_rank) /
      (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)),
    mean_abs_gap_over_sequence =
      trapz_xy(x, abs(gap_pipe_minus_rank)) /
      (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)),
    fraction_removal_axis_pipe_better =
      fraction_axis_above_zero(x, gap_pipe_minus_rank),
    max_gap_pipe_minus_rank =
      max(gap_pipe_minus_rank, na.rm = TRUE),
    min_gap_pipe_minus_rank =
      min(gap_pipe_minus_rank, na.rm = TRUE),
    x_at_max_gap =
      x[which.max(gap_pipe_minus_rank)]
  ),
  by = .(species, scientificName, className)
]

species_auc <- species_auc[is.finite(auc_diff_pipe_minus_rank)]
species_auc[, taxon := taxon_scope(className)]
species_auc[, dominant_method := fcase(
  mean_gap_over_sequence > metric_eps, "Persistence-based",
  mean_gap_over_sequence < -metric_eps, rank_label,
  default = "No material difference"
)]
setorder(species_auc, -auc_diff_pipe_minus_rank)
species_auc[, positive_auc_rank := seq_len(.N)]

n_species_auc <- nrow(species_auc)
n_positive_auc <- species_auc[
  ,
  sum(auc_diff_pipe_minus_rank > metric_eps, na.rm = TRUE)
]
n_negative_auc <- species_auc[
  ,
  sum(auc_diff_pipe_minus_rank < -metric_eps, na.rm = TRUE)
]

focal_auc_row <- species_auc[species == focal_id_results]

auc_scope_input <- rbindlist(list(
  copy(species_auc)[, scope := "All species"],
  copy(species_auc)[, scope := taxon]
), use.names = TRUE)
auc_scope_summary <- auc_scope_input[, .(
  n_species = .N,
  n_pipe_advantaged = sum(mean_gap_over_sequence > metric_eps),
  n_rank_advantaged = sum(mean_gap_over_sequence < -metric_eps),
  n_near_zero = sum(abs(mean_gap_over_sequence) <= metric_eps),
  frac_pipe_advantaged = mean(mean_gap_over_sequence > metric_eps),
  frac_rank_advantaged = mean(mean_gap_over_sequence < -metric_eps),
  frac_near_zero = mean(abs(mean_gap_over_sequence) <= metric_eps),
  mean_sequence_gap = mean(mean_gap_over_sequence),
  median_sequence_gap = metric_median(mean_gap_over_sequence),
  q10_sequence_gap = metric_quantile(mean_gap_over_sequence, 0.10),
  q90_sequence_gap = metric_quantile(mean_gap_over_sequence, 0.90),
  median_fraction_axis_pipe_better =
    metric_median(fraction_removal_axis_pipe_better),
  # Proportion of total directional-difference magnitude favoring each
  # method: positive and negative signed-AUC magnitudes summed separately
  # so opposing species responses do not cancel before being compared.
  total_positive_magnitude =
    sum(auc_diff_pipe_minus_rank[auc_diff_pipe_minus_rank > metric_eps]),
  total_negative_magnitude =
    -sum(auc_diff_pipe_minus_rank[auc_diff_pipe_minus_rank < -metric_eps]),
  frac_total_magnitude_pipe = {
    pos_mag <- sum(auc_diff_pipe_minus_rank[auc_diff_pipe_minus_rank > metric_eps])
    neg_mag <- -sum(auc_diff_pipe_minus_rank[auc_diff_pipe_minus_rank < -metric_eps])
    if ((pos_mag + neg_mag) > metric_eps) pos_mag / (pos_mag + neg_mag) else NA_real_
  }
), by = scope]
auc_scope_summary[, scope_order := match(scope, c("All species", "Mammals", "Birds"))]
setorder(auc_scope_summary, scope_order, scope)
auc_scope_summary[, scope_order := NULL]

# ---- Retained-area-normalized persistence advantage -------------------------

# Each species is compared on the overlap of the methods' retained-area curves.
# Integration uses km2 knots over that shared support and never extrapolates.

area_norm_one_species <- function(dd) {
  sp_id <- dd$species[1L]
  sp_name <- dd$scientificName[1L]
  class_name <- dd$className[1L]

  empty_result <- function() {
    data.table(
      species = sp_id,
      scientificName = sp_name,
      className = class_name,
      n_area_knots = NA_integer_,
      common_area_min_km2 = NA_real_,
      common_area_max_km2 = NA_real_,
      common_area_span_km2 = NA_real_,
      auc_diff_area_pipe_minus_rank = NA_real_,
      mean_gap_over_common_area = NA_real_,
      mean_abs_gap_over_common_area = NA_real_,
      median_gap_over_common_area = NA_real_,
      max_gap_pipe_minus_rank = NA_real_,
      min_gap_pipe_minus_rank = NA_real_,
      fraction_common_area_pipe_better = NA_real_
    )
  }

  pipe_curve <- dd[
    method == "pipe",
    .(
      area = as.numeric(total_pu_area_km2),
      persistence = as.numeric(sp_persist)
    )
  ]
  rank_curve <- dd[
    method == "rank",
    .(
      area = as.numeric(total_pu_area_km2),
      persistence = as.numeric(sp_persist)
    )
  ]
  comparison <- benchmark_exact_area_comparison(
    pipe_curve$area,
    pipe_curve$persistence,
    rank_curve$area,
    rank_curve$persistence
  )
  if (!isTRUE(comparison$valid)) {
    return(empty_result())
  }

  area_knots <- comparison$area_knots
  gap <- comparison$difference

  data.table(
    species = sp_id,
    scientificName = sp_name,
    className = class_name,
    n_area_knots = as.integer(length(area_knots)),
    common_area_min_km2 = comparison$lower,
    common_area_max_km2 = comparison$upper,
    common_area_span_km2 = comparison$span,
    auc_diff_area_pipe_minus_rank = comparison$signed_auc,
    mean_gap_over_common_area = comparison$mean_difference,
    mean_abs_gap_over_common_area =
      benchmark_absolute_trapezoid(area_knots, gap) / comparison$span,
    median_gap_over_common_area = metric_median(gap),
    max_gap_pipe_minus_rank = max(gap, na.rm = TRUE),
    min_gap_pipe_minus_rank = min(gap, na.rm = TRUE),
    fraction_common_area_pipe_better =
      fraction_axis_above_zero(area_knots, gap)
  )
}

area_norm_dt <- rbindlist(
  lapply(
    split(
      sp_long[
        ,
        .(
          method,
          species,
          scientificName,
          className,
          stage,
          stage_order,
          x,
          total_pu_area_km2,
          sp_persist
        )
      ],
      by = "species",
      keep.by = TRUE
    ),
    area_norm_one_species
  ),
  use.names = TRUE,
  fill = TRUE
)

area_norm_dt <- area_norm_dt[is.finite(mean_gap_over_common_area)]
area_norm_dt[, taxon := taxon_scope(className)]
area_norm_dt[, dominant_method := fcase(
  mean_gap_over_common_area > metric_eps, "Persistence-based",
  mean_gap_over_common_area < -metric_eps, rank_label,
  default = "No material difference"
)]
setorder(area_norm_dt, -mean_gap_over_common_area)
area_norm_dt[, positive_area_rank := seq_len(.N)]

n_area_species <- nrow(area_norm_dt)
n_area_positive <- area_norm_dt[
  ,
  sum(mean_gap_over_common_area > metric_eps, na.rm = TRUE)
]
n_area_negative <- area_norm_dt[
  ,
  sum(mean_gap_over_common_area < -metric_eps, na.rm = TRUE)
]

area_scope_input <- rbindlist(list(
  copy(area_norm_dt)[, scope := "All species"],
  copy(area_norm_dt)[, scope := taxon]
), use.names = TRUE)
area_scope_summary <- area_scope_input[, .(
  n_species = .N,
  n_pipe_advantaged = sum(mean_gap_over_common_area > metric_eps),
  n_rank_advantaged = sum(mean_gap_over_common_area < -metric_eps),
  n_near_zero = sum(abs(mean_gap_over_common_area) <= metric_eps),
  frac_pipe_advantaged = mean(mean_gap_over_common_area > metric_eps),
  frac_rank_advantaged = mean(mean_gap_over_common_area < -metric_eps),
  frac_near_zero = mean(abs(mean_gap_over_common_area) <= metric_eps),
  mean_standardized_gap = mean(mean_gap_over_common_area),
  median_standardized_gap = metric_median(mean_gap_over_common_area),
  q10_standardized_gap = metric_quantile(mean_gap_over_common_area, 0.10),
  q90_standardized_gap = metric_quantile(mean_gap_over_common_area, 0.90),
  median_fraction_area_pipe_better =
    metric_median(fraction_common_area_pipe_better),
  total_positive_magnitude =
    sum(auc_diff_area_pipe_minus_rank[auc_diff_area_pipe_minus_rank > metric_eps]),
  total_negative_magnitude =
    -sum(auc_diff_area_pipe_minus_rank[auc_diff_area_pipe_minus_rank < -metric_eps]),
  frac_total_magnitude_pipe = {
    pos_mag <- sum(auc_diff_area_pipe_minus_rank[auc_diff_area_pipe_minus_rank > metric_eps])
    neg_mag <- -sum(auc_diff_area_pipe_minus_rank[auc_diff_area_pipe_minus_rank < -metric_eps])
    if ((pos_mag + neg_mag) > metric_eps) pos_mag / (pos_mag + neg_mag) else NA_real_
  }
), by = scope]
area_scope_summary[, scope_order := match(scope, c("All species", "Mammals", "Birds"))]
setorder(area_scope_summary, scope_order, scope)
area_scope_summary[, scope_order := NULL]

# ---- Continuous/interpolated threshold crossing -----------------------------

# Crossing locations are linearly interpolated on the removal-sequence axis;
# threshold equality follows the strict classification contract above.

interp_first_crossing <- function(x, y, stage_order, threshold) {
  ok <- is.finite(x) & is.finite(y) & is.finite(stage_order)
  x <- x[ok]
  y <- y[ok]
  stage_order <- stage_order[ok]

  if (!length(x)) {
    return(data.table(
      crossed = FALSE,
      x_cross_continuous = NA_real_,
      x_cross_or_end = NA_real_
    ))
  }

  ord <- order(stage_order)
  x <- x[ord]
  y <- y[ord]

  if (y[1L] < threshold) {
    return(data.table(
      crossed = TRUE,
      x_cross_continuous = x[1L],
      x_cross_or_end = x[1L]
    ))
  }

  trans <- which(head(y, -1L) >= threshold & tail(y, -1L) < threshold)

  if (!length(trans)) {
    return(data.table(
      crossed = FALSE,
      x_cross_continuous = NA_real_,
      x_cross_or_end = max(x, na.rm = TRUE)
    ))
  }

  ii <- trans[1L]
  x1 <- x[ii]
  x2 <- x[ii + 1L]
  y1 <- y[ii]
  y2 <- y[ii + 1L]

  x_cross <- if (!is.finite(y2 - y1) || abs(y2 - y1) < .Machine$double.eps) {
    x2
  } else {
    x1 + (threshold - y1) * (x2 - x1) / (y2 - y1)
  }

  data.table(
    crossed = TRUE,
    x_cross_continuous = x_cross,
    x_cross_or_end = x_cross
  )
}

threshold_values <- c(0.9, 0.5, 0.1)

threshold_long <- rbindlist(lapply(threshold_values, function(thr) {
  sp_long[
    ,
    {
      out <- interp_first_crossing(
        x = x,
        y = sp_persist,
        stage_order = stage_order,
        threshold = thr
      )
      out[, threshold := thr]
      out
    },
    by = .(method, species, scientificName, className)
  ]
}), use.names = TRUE, fill = TRUE)

threshold_wide <- dcast(
  threshold_long,
  threshold + species + scientificName + className ~ method,
  value.var = c(
    "crossed",
    "x_cross_continuous",
    "x_cross_or_end"
  )
)

threshold_wide[
  ,
  crossing_state := fcase(
    crossed_pipe & crossed_rank, "both_crossed",
    !crossed_pipe & crossed_rank, "rank_only",
    crossed_pipe & !crossed_rank, "pipe_only",
    !crossed_pipe & !crossed_rank, "neither_crossed",
    default = NA_character_
  )
]

threshold_wide[
  ,
  delay_x_pipe_minus_rank_both_crossed :=
    fifelse(
      crossed_pipe & crossed_rank,
      x_cross_continuous_pipe - x_cross_continuous_rank,
      NA_real_
    )
]

threshold_summary <- threshold_wide[
  ,
  .(
    n_species = .N,
    n_both_crossed = sum(crossing_state == "both_crossed", na.rm = TRUE),
    n_rank_only = sum(crossing_state == "rank_only", na.rm = TRUE),
    n_pipe_only = sum(crossing_state == "pipe_only", na.rm = TRUE),
    n_neither_crossed = sum(crossing_state == "neither_crossed", na.rm = TRUE),
    mean_delay_both_crossed =
      mean(delay_x_pipe_minus_rank_both_crossed, na.rm = TRUE),
    median_delay_both_crossed =
      metric_median(delay_x_pipe_minus_rank_both_crossed)
  ),
  by = threshold
]

setorder(threshold_summary, -threshold)

# PU redundancy-loss timing

method_stage_grid <- unique(sp_long[, .(method, stage, stage_order, x)])
species_grid <- unique(sp_long[, .(species, scientificName, className)])

method_stage_grid[, tmp_join := 1L]
species_grid[, tmp_join := 1L]

full_species_stage_grid <- merge(
  method_stage_grid,
  species_grid,
  by = "tmp_join",
  allow.cartesian = TRUE,
  sort = FALSE
)[, tmp_join := NULL]

pu_redundancy_obs <- pu_long[
  ,
  .(
    n_pu_current = .N,
    n_pu_gt_09 = sum(P_pu > 0.9, na.rm = TRUE),
    n_pu_gt_05 = sum(P_pu > 0.5, na.rm = TRUE),
    expected_persisting_pu = sum(P_pu, na.rm = TRUE),
    max_P_pu = max(P_pu, na.rm = TRUE)
  ),
  by = .(
    method,
    stage,
    stage_order,
    x,
    species,
    scientificName
  )
]

pu_redundancy_full <- merge(
  full_species_stage_grid,
  pu_redundancy_obs,
  by = c(
    "method",
    "stage",
    "stage_order",
    "x",
    "species",
    "scientificName"
  ),
  all.x = TRUE,
  sort = FALSE
)

for (cc in c(
  "n_pu_current",
  "n_pu_gt_09",
  "n_pu_gt_05",
  "expected_persisting_pu",
  "max_P_pu"
)) {
  pu_redundancy_full[is.na(get(cc)), (cc) := 0]
}

first_condition_time <- function(x, y, stage_order, condition_fun) {
  ok <- is.finite(x) & is.finite(y) & is.finite(stage_order)
  x <- x[ok]
  y <- y[ok]
  stage_order <- stage_order[ok]

  if (!length(x)) {
    return(data.table(
      condition_reached = FALSE,
      x_condition = NA_real_,
      x_condition_or_end = NA_real_
    ))
  }

  ord <- order(stage_order)
  x <- x[ord]
  y <- y[ord]

  cond <- condition_fun(y)
  reached <- any(cond, na.rm = TRUE)

  if (reached) {
    ii <- which(cond)[1L]
    x_condition <- x[ii]
  } else {
    x_condition <- NA_real_
  }

  data.table(
    condition_reached = reached,
    x_condition = x_condition,
    x_condition_or_end = if (reached) x_condition else max(x, na.rm = TRUE)
  )
}

redundancy_specs <- data.table(
  target_id = c(
    "no_PU_gt_0.9",
    "no_PU_gt_0.5",
    "expected_persisting_PU_lt_1"
  ),
  metric = c(
    "n_pu_gt_09",
    "n_pu_gt_05",
    "expected_persisting_pu"
  ),
  threshold = c(0, 0, 1),
  condition_type = c("le", "le", "lt"),
  manuscript_label = c(
    "lost all PUs with P_PU > 0.9",
    "lost all PUs with P_PU > 0.5",
    "fell below one expected persisting PU"
  )
)

redundancy_long <- rbindlist(lapply(seq_len(nrow(redundancy_specs)), function(ii) {
  metric_i <- redundancy_specs$metric[ii]
  threshold_i <- redundancy_specs$threshold[ii]
  condition_type_i <- redundancy_specs$condition_type[ii]

  condition_fun <- switch(
    condition_type_i,
    le = function(z) z <= threshold_i,
    lt = function(z) z < threshold_i,
    stop("Unknown condition type.")
  )

  out <- pu_redundancy_full[
    ,
    first_condition_time(
      x = x,
      y = get(metric_i),
      stage_order = stage_order,
      condition_fun = condition_fun
    ),
    by = .(method, species, scientificName, className)
  ]

  out[, target_id := redundancy_specs$target_id[ii]]
  out[, metric := metric_i]
  out[, threshold := threshold_i]
  out[, condition_type := condition_type_i]
  out[, manuscript_label := redundancy_specs$manuscript_label[ii]]
  out
}), use.names = TRUE, fill = TRUE)

redundancy_wide <- dcast(
  redundancy_long,
  target_id + metric + threshold + condition_type + manuscript_label +
    species + scientificName + className ~ method,
  value.var = c(
    "condition_reached",
    "x_condition",
    "x_condition_or_end"
  )
)

redundancy_wide[
  ,
  loss_state := fcase(
    condition_reached_pipe & condition_reached_rank, "both_lost",
    !condition_reached_pipe & condition_reached_rank, "rank_only_lost",
    condition_reached_pipe & !condition_reached_rank, "pipe_only_lost",
    !condition_reached_pipe & !condition_reached_rank, "neither_lost",
    default = NA_character_
  )
]

redundancy_summary <- redundancy_wide[
  ,
  .(
    n_species = .N,
    n_both_lost = sum(loss_state == "both_lost", na.rm = TRUE),
    n_rank_only_lost = sum(loss_state == "rank_only_lost", na.rm = TRUE),
    n_pipe_only_lost = sum(loss_state == "pipe_only_lost", na.rm = TRUE),
    n_neither_lost = sum(loss_state == "neither_lost", na.rm = TRUE),
    mean_censored_delay =
      mean(x_condition_or_end_pipe - x_condition_or_end_rank, na.rm = TRUE),
    median_censored_delay =
      metric_median(x_condition_or_end_pipe - x_condition_or_end_rank)
  ),
  by = .(target_id, manuscript_label)
]

# ---- Demographic-curve sensitivity of the key late-stage contrast -----------

# Alternative curves are evaluated along the already completed fixed rankings.
# This is sensitivity analysis of outcomes, not re-optimization under each curve.

uncertainty_checkpoints <- data.table()
uncertainty_summary <- data.table()
uncertainty_species_auc <- data.table()
uncertainty_auc_summary <- data.table()
uncertainty_auc_stability <- data.table()
if (!is.null(figure_data)) {
  assert(
    is.list(figure_data) &&
      all(c("stage", "species") %in% names(figure_data)),
    "figure_data must be the result returned by prepare_stage73_figure_data()."
  )
  curve_stage <- copy(figure_data$stage)
  need_cols(
    curve_stage,
    c("method", "curve", "stage", "stage_order", "x", "mean_persist"),
    "Stage 7.3 uncertainty checkpoint data"
  )
  curve_checkpoint_long <- merge(
    curve_stage,
    checkpoint_key[, .(
      stage,
      stage_order,
      x,
      target_pct_removed,
      actual_pct_removed
    )],
    by = c("stage", "stage_order", "x"),
    all.x = FALSE,
    sort = FALSE
  )
  uncertainty_checkpoints <- dcast(
    curve_checkpoint_long,
    target_pct_removed + actual_pct_removed + stage + stage_order + x + curve ~ method,
    value.var = "mean_persist"
  )
  uncertainty_checkpoints[, `:=`(
    mean_gap_pipe_minus_rank = pipe - rank,
    favoured_method = fcase(
      pipe - rank > metric_eps, "Persistence-based",
      pipe - rank < -metric_eps, rank_label,
      default = "No material difference"
    )
  )]
  uncertainty_checkpoints[, curve_order := match(curve, persistence_curves())]
  setorder(uncertainty_checkpoints, target_pct_removed, curve_order)
  uncertainty_checkpoints[, curve_order := NULL]
  uncertainty_summary <- uncertainty_checkpoints[, .(
    q50_mean_gap = mean_gap_pipe_minus_rank[curve == "q50"][1L],
    minimum_curve_gap = min(mean_gap_pipe_minus_rank),
    maximum_curve_gap = max(mean_gap_pipe_minus_rank),
    n_curves_pipe_advantaged = sum(mean_gap_pipe_minus_rank > metric_eps),
    n_curves_rank_advantaged = sum(mean_gap_pipe_minus_rank < -metric_eps),
    n_curves_near_zero = sum(abs(mean_gap_pipe_minus_rank) <= metric_eps),
    direction_consistent =
      all(mean_gap_pipe_minus_rank > metric_eps) ||
      all(mean_gap_pipe_minus_rank < -metric_eps) ||
      all(abs(mean_gap_pipe_minus_rank) <= metric_eps)
  ), by = .(
    target_pct_removed,
    actual_pct_removed,
    stage,
    stage_order,
    x
  )]
  setorder(uncertainty_summary, target_pct_removed)

  curve_species <- copy(figure_data$species)
  need_cols(
    curve_species,
    c(
      "method", "curve", "species", "scientificName", "className",
      "stage", "stage_order", "x", "sp_persist"
    ),
    "Stage 7.3 uncertainty species trajectories"
  )
  curve_species[, curve := as.character(curve)]
  assert(
    setequal(unique(curve_species$curve), persistence_curves()) &&
      !anyDuplicated(curve_species[, .(method, curve, species, stage)]) &&
      all(
        is.finite(curve_species$sp_persist) &
          curve_species$sp_persist >= 0 &
          curve_species$sp_persist <= 1
      ),
    "Stage 7.3 uncertainty species trajectories are incomplete or invalid."
  )

  uncertainty_auc_wide <- dcast(
    curve_species[
      ,
      .(
        curve,
        species,
        scientificName,
        className,
        stage,
        stage_order,
        x,
        method,
        sp_persist
      )
    ],
    curve + species + scientificName + className +
      stage + stage_order + x ~ method,
    value.var = "sp_persist"
  )
  assert(
    all(c("pipe", "rank") %in% names(uncertainty_auc_wide)) &&
      all(is.finite(uncertainty_auc_wide$pipe)) &&
      all(is.finite(uncertainty_auc_wide$rank)),
    "Stage 7.3 uncertainty AUC trajectories are missing a method or persistence value."
  )
  uncertainty_auc_wide[, gap_pipe_minus_rank := pipe - rank]
  setorder(uncertainty_auc_wide, curve, species, stage_order)

  uncertainty_species_auc <- uncertainty_auc_wide[
    ,
    {
      axis_span <- max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
      signed_auc <- trapz_xy(x, gap_pipe_minus_rank)
      absolute_auc <- trapz_xy(x, abs(gap_pipe_minus_rank))
      list(
        n_stages = .N,
        removal_axis_span = axis_span,
        auc_diff_pipe_minus_rank = signed_auc,
        mean_gap_over_sequence = signed_auc / axis_span,
        mean_abs_gap_over_sequence = absolute_auc / axis_span,
        fraction_removal_axis_pipe_better =
          fraction_axis_above_zero(x, gap_pipe_minus_rank),
        max_gap_pipe_minus_rank = max(gap_pipe_minus_rank, na.rm = TRUE),
        min_gap_pipe_minus_rank = min(gap_pipe_minus_rank, na.rm = TRUE)
      )
    },
    by = .(curve, species, scientificName, className)
  ]
  uncertainty_species_auc <- uncertainty_species_auc[
    is.finite(removal_axis_span) &
      removal_axis_span > 0 &
      is.finite(auc_diff_pipe_minus_rank) &
      is.finite(mean_gap_over_sequence)
  ]
  uncertainty_species_auc[, `:=`(
    taxon = taxon_scope(className),
    auc_direction = fcase(
      mean_gap_over_sequence > metric_eps, "Positive",
      mean_gap_over_sequence < -metric_eps, "Negative",
      default = "Near zero"
    )
  )]
  uncertainty_species_auc[, curve_order := match(curve, persistence_curves())]
  setorder(
    uncertainty_species_auc,
    curve_order,
    -auc_diff_pipe_minus_rank,
    scientificName
  )
  uncertainty_species_auc[, curve_order := NULL]

  uncertainty_auc_summary <- uncertainty_species_auc[
    ,
    {
      positive_auc <- auc_diff_pipe_minus_rank[
        mean_gap_over_sequence > metric_eps
      ]
      negative_auc_magnitude <- -auc_diff_pipe_minus_rank[
        mean_gap_over_sequence < -metric_eps
      ]
      list(
        n_species = .N,
        n_positive_auc = length(positive_auc),
        n_negative_auc = length(negative_auc_magnitude),
        n_near_zero_auc = sum(abs(mean_gap_over_sequence) <= metric_eps),
        pct_positive_auc = 100 * length(positive_auc) / .N,
        pct_negative_auc = 100 * length(negative_auc_magnitude) / .N,
        total_positive_auc = sum(positive_auc),
        total_negative_auc_magnitude = sum(negative_auc_magnitude),
        mean_positive_auc = if (length(positive_auc)) mean(positive_auc) else NA_real_,
        mean_negative_auc_magnitude =
          if (length(negative_auc_magnitude)) {
            mean(negative_auc_magnitude)
          } else {
            NA_real_
          },
        net_signed_auc = sum(auc_diff_pipe_minus_rank),
        median_signed_auc = metric_median(auc_diff_pipe_minus_rank),
        mean_standardized_gap = mean(mean_gap_over_sequence),
        median_standardized_gap = metric_median(mean_gap_over_sequence)
      )
    },
    by = curve
  ]
  uncertainty_auc_summary[, curve_order := match(curve, persistence_curves())]
  setorder(uncertainty_auc_summary, curve_order)
  uncertainty_auc_summary[, curve_order := NULL]

  uncertainty_auc_stability <- uncertainty_species_auc[
    ,
    {
      n_positive <- sum(mean_gap_over_sequence > metric_eps)
      n_negative <- sum(mean_gap_over_sequence < -metric_eps)
      n_near_zero <- sum(abs(mean_gap_over_sequence) <= metric_eps)
      q50_auc <- auc_diff_pipe_minus_rank[curve == "q50"][1L]
      list(
        n_curves_positive = n_positive,
        n_curves_negative = n_negative,
        n_curves_near_zero = n_near_zero,
        q50_signed_auc = q50_auc,
        minimum_signed_auc = min(auc_diff_pipe_minus_rank),
        maximum_signed_auc = max(auc_diff_pipe_minus_rank),
        signed_auc_range =
          max(auc_diff_pipe_minus_rank) - min(auc_diff_pipe_minus_rank),
        direction_stability = fcase(
          n_positive == length(persistence_curves()), "Positive under every curve",
          n_negative == length(persistence_curves()), "Negative under every curve",
          n_near_zero == length(persistence_curves()), "Near zero under every curve",
          n_negative == 0L, "Nonnegative with near-zero curves",
          n_positive == 0L, "Nonpositive with near-zero curves",
          default = "Direction changes across curves"
        )
      )
    },
    by = .(species, scientificName, className, taxon)
  ]
  setorder(
    uncertainty_auc_stability,
    -signed_auc_range,
    scientificName
  )
}

list(
  focus = list(
    species_name = focal_name_results,
    species_id = focal_id_results,
    auc_row = focal_auc_row
  ),
  late_stage = list(
    all_stages = all_stage_wide,
    checkpoints = checkpoint_stage_wide,
    species = checkpoint_species_wide,
    thresholds = checkpoint_threshold_summary
  ),
  auc = list(
    species = species_auc,
    scope_summary = auc_scope_summary,
    n_species = n_species_auc,
    n_positive = n_positive_auc,
    n_negative = n_negative_auc
  ),
  area_normalized = list(
    species = area_norm_dt,
    scope_summary = area_scope_summary,
    n_species = n_area_species,
    n_positive = n_area_positive,
    n_negative = n_area_negative
  ),
  threshold_crossings = list(
    species = threshold_wide,
    summary = threshold_summary
  ),
  redundancy = list(
    species = redundancy_wide,
    summary = redundancy_summary
  ),
  uncertainty = list(
    checkpoints = uncertainty_checkpoints,
    summary = uncertainty_summary,
    species_auc = uncertainty_species_auc,
    auc_summary = uncertainty_auc_summary,
    auc_stability = uncertainty_auc_stability
  )
)
}
