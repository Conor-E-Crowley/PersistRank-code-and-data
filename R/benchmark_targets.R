# Shared benchmark target discovery and union planning.
#
# Each consumer row identifies one completed Stage 6 state. Equal retained-cell
# counts share one union target, so the benchmark state is reconstructed once
# and its compact persistence rows are reused by every matching curve/stage.


build_benchmark_target_plan <- function(run_states) {
  run_meta <- run_states
  curves <- names(run_meta)
  assert(length(curves) > 0L && !anyDuplicated(curves) &&
           all(curves %in% persistence_curves()) && identical(names(run_meta), curves),
         "Target inputs must use a nonempty canonical curve subset.")
  consumers <- data.table::rbindlist(lapply(curves, function(curve) {
    x <- run_meta[[curve]]
    x[, .(
      optimization_curve, run_id, stage = as.integer(stage),
      stage_order = as.integer(stage_order), keep_n = as.integer(keep_n),
      target_cells_retained = as.integer(target_cells_retained),
      target_area_retained_km2 = as.numeric(target_area_retained_km2),
      pct_cells_removed = as.numeric(pct_cells_removed)
    )]
  }), use.names = TRUE)
  consumers[, curve_order := match(optimization_curve, curves)]
  data.table::setorder(consumers, curve_order, stage_order)
  assert(!anyDuplicated(consumers[, .(optimization_curve, stage)]),
         "Curve/stage target consumers must be unique.")
  assert(all(consumers[, !is.unsorted(-keep_n), by = optimization_curve]$V1),
         "Retained-cell targets must be non-increasing within every run.")

  unique_keep <- sort(unique(consumers$keep_n), decreasing = TRUE)
  union <- data.table::data.table(
    union_target_id = seq_along(unique_keep) - 1L,
    keep_n = as.integer(unique_keep)
  )
  union[, previous_keep_n := data.table::shift(keep_n, fill = keep_n[[1L]])]
  union[, rank_cells_removed := as.numeric(previous_keep_n) - as.numeric(keep_n)]
  consumers[, union_target_id := union$union_target_id[match(keep_n, union$keep_n)]]
  summary <- consumers[, .(
    n_consumers = .N,
    consumer_curves = paste(optimization_curve, collapse = ","),
    consumer_stages = paste(stage, collapse = ",")
  ), by = union_target_id]
  union <- merge(union, summary, by = "union_target_id", sort = FALSE)
  data.table::setorder(union, union_target_id)

  initial_cells <- union$keep_n[[1L]]
  independent <- sum(vapply(run_meta, function(x) {
    initial_cells - min(as.numeric(x$keep_n))
  }, numeric(1L)))
  union_traversal <- initial_cells - min(union$keep_n)
  reduction <- if (independent > 0) 100 * (1 - union_traversal / independent) else 0
  diagnostics <- data.table::data.table(
    requested_targets = nrow(consumers),
    unique_targets = nrow(union),
    duplicate_targets = nrow(consumers) - nrow(union),
    initial_cells = initial_cells,
    independent_rank_cells_traversed = independent,
    union_rank_cells_traversed = union_traversal,
    theoretical_reduction_pct = reduction
  )
  assert(
    nrow(consumers) == sum(vapply(run_meta, nrow, integer(1L))) &&
      nrow(union) == data.table::uniqueN(consumers$keep_n),
    "Benchmark target union is incomplete."
  )
  consumers[, curve_order := NULL]
  list(consumers = consumers[], union = union[], diagnostics = diagnostics)
}


benchmark_target_plans_equal <- function(observed, expected) {
  observed <- data.table::as.data.table(observed)
  expected <- data.table::as.data.table(expected)
  if (!identical(names(observed), names(expected)) || nrow(observed) != nrow(expected)) {
    return(FALSE)
  }
  all(vapply(names(expected), function(field) {
    a <- observed[[field]]
    b <- expected[[field]]
    if (is.numeric(a) || is.numeric(b)) {
      length(a) == length(b) && all(abs(as.numeric(a) - as.numeric(b)) <= 1e-10)
    } else {
      identical(as.character(a), as.character(b))
    }
  }, logical(1L)))
}
