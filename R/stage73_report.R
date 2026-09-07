# Stage 7.3 format-neutral report tables.
#
# Converts validated detailed-comparison statistics to readable data tables.
# Reads and writes no files, has no rendering dependency, and changes no
# scientific value.

build_stage73_report_tables <- function(core, statistics) {
  assert(is.list(core) && is.list(statistics),
         "core and statistics must be Stage 7.3 result objects.")
  late_stage_threshold <- data.table::copy(statistics$late_stage$thresholds)[
    scope == "All species",
    .(
      `Target removal (%)` = target_pct_removed,
      `Actual removal (%)` = actual_pct_removed,
      Stage = stage,
      `Persistence-based below threshold` = n_below_pipe,
      `Benchmark below threshold` = n_below_rank,
      `Persistence-based only` = n_pipe_only_below,
      `Benchmark only` = n_rank_only_below,
      `Both below` = n_both_below
    )
  ]
  data.table::setorder(late_stage_threshold, `Target removal (%)`)

  scope_table <- function(table, comparison_label) {
    table <- data.table::copy(table)
    table[, Comparison := comparison_label]
    table[, .(
      Comparison, Scope = scope, Species = n_species,
      `Favoring persistence-based` = n_pipe_advantaged,
      `Favoring persistence-based (%)` = 100 * frac_pipe_advantaged,
      `Favoring benchmark` = n_rank_advantaged,
      `Favoring benchmark (%)` = 100 * frac_rank_advantaged,
      `Near-zero` = n_near_zero,
      `Near-zero (%)` = 100 * frac_near_zero,
      `Persistence-based magnitude share (%)` = 100 * frac_total_magnitude_pipe
    )]
  }
  favoring_by_scope <- data.table::rbindlist(list(
    scope_table(statistics$auc$scope_summary, "Removal sequence"),
    scope_table(statistics$area_normalized$scope_summary, "Retained area")
  ), use.names = TRUE)
  top3_positive <- data.table::copy(statistics$auc$species)[
    mean_gap_over_sequence > 1e-9
  ][order(-mean_gap_over_sequence)][seq_len(min(3L, .N)), .(
    Rank = seq_len(.N), Species = scientificName, Taxon = taxon,
    `Standardized mean gap` = mean_gap_over_sequence,
    `Signed AUC` = auc_diff_pipe_minus_rank,
    `Removal axis favoring pipeline (%)` = 100 * fraction_removal_axis_pipe_better
  )]
  sensitivity_checkpoints <- data.table::copy(statistics$uncertainty$summary)
  if (nrow(sensitivity_checkpoints)) {
    sensitivity_checkpoints <- sensitivity_checkpoints[, .(
      `Target removal (%)` = target_pct_removed,
      `Actual removal (%)` = actual_pct_removed,
      `q50 mean difference` = q50_mean_gap,
      `Curves favoring pipeline` = n_curves_pipe_advantaged,
      `Curves favoring benchmark` = n_curves_rank_advantaged,
      `Direction consistent` = ifelse(direction_consistent, "Yes", "No")
    )]
  }
  sensitivity_stability <- data.table::copy(
    statistics$uncertainty$auc_stability
  )
  sensitivity_stability <- if (nrow(sensitivity_stability)) {
    sensitivity_stability[, .(Species = .N),
      by = .(`Direction across five curves` = direction_stability)
    ][order(-Species, `Direction across five curves`)]
  } else data.table::data.table()
  attr(late_stage_threshold, "benchmark_label") <- core$rank_label
  list(
    late_stage_threshold = late_stage_threshold,
    favoring_by_scope = favoring_by_scope,
    top3_positive = top3_positive,
    sensitivity_checkpoints = sensitivity_checkpoints,
    sensitivity_stability = sensitivity_stability
  )
}
