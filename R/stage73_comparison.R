# Canonical in-memory Stage 7.3 persistence comparison.
#
# Inputs are one completed Stage 6 curve and one complete exact-target benchmark
# library. Returns PU-, species-, and stage-level comparison tables required by
# the detailed figure. Reads the shared initialization, pipeline lookups, exact
# benchmark lookups, and persistence coefficients. Writes nothing. Cost is
# moderate and entirely report-time; no rank raster, kernel, checkpoint, or
# derived cache is touched.


compute_stage73_method_trajectory <- function(stage_table, method, load_lookup,
                                              species_parameters) {
  population_units <- vector("list", nrow(stage_table))
  species_states <- vector("list", nrow(stage_table))
  parameter_columns <- species_parameters[, .(
    scientificName, species, className, density, a_pred, b_pred, c_th
  )]
  for (index in seq_len(nrow(stage_table))) {
    stage <- as.integer(stage_table$stage[[index]])
    lookup <- load_lookup(stage, stage_table$path[[index]])
    units <- if (nrow(lookup)) lookup[, .(
      pu_area_km2 = sum(patch_area_km2), n_patches = .N
    ), by = .(scientificName, species, pu_id)] else data.table::data.table()
    if (nrow(units)) {
      units <- merge(units, parameter_columns,
                     by = c("scientificName", "species"), all.x = TRUE, sort = FALSE)
      assert(all(!is.na(units$density)),
             paste0("Missing parameters for ", method, " stage ", stage, "."))
      units[, P_pu := stage7_pu_persistence(
        pu_area_km2, density, a_pred, b_pred, c_th,
        paste0(method, " stage ", stage, " population units")
      )]
      unit_output <- units[, .(
        method, stage, stage_order = stage_table$stage_order[[index]],
        x = stage_table$x[[index]], scientificName, species, pu_id,
        pu_area_km2, n_patches, P_pu
      )]
      represented <- unit_output[, .(
        sp_persist = stage7_species_persistence(P_pu), n_pu = .N,
        total_pu_area_km2 = sum(pu_area_km2)
      ), by = .(scientificName, species)]
    } else {
      unit_output <- data.table::data.table(
        method = character(), stage = integer(), stage_order = integer(), x = numeric(),
        scientificName = character(), species = character(), pu_id = integer(),
        pu_area_km2 = numeric(), n_patches = integer(), P_pu = numeric()
      )
      represented <- data.table::data.table(
        scientificName = character(), species = character(),
        sp_persist = numeric(), n_pu = integer(), total_pu_area_km2 = numeric()
      )
    }
    complete <- merge(
      species_parameters[, .(scientificName, species, className)], represented,
      by = c("scientificName", "species"), all.x = TRUE, sort = FALSE
    )
    complete[is.na(sp_persist), `:=`(
      sp_persist = 0, n_pu = 0L, total_pu_area_km2 = 0
    )]
    complete[, `:=`(
      method = method, stage = stage,
      stage_order = stage_table$stage_order[[index]], x = stage_table$x[[index]]
    )]
    data.table::setcolorder(complete, c(
      "method", "stage", "stage_order", "x", "scientificName", "species",
      "className", "sp_persist", "n_pu", "total_pu_area_km2"
    ))
    population_units[[index]] <- unit_output
    species_states[[index]] <- complete
    rm(lookup, units, unit_output, represented, complete)
  }
  list(
    pu_long = data.table::rbindlist(population_units, use.names = TRUE, fill = TRUE),
    sp_long = data.table::rbindlist(species_states, use.names = TRUE, fill = TRUE)
  )
}

summarize_stage73_species_trajectory <- function(species_trajectory) {
  data.table::as.data.table(species_trajectory)[, .(
    n_species = .N, n_species_with_pu = sum(n_pu > 0),
    mean_persist = mean(sp_persist), median_persist = stats::median(sp_persist),
    q10_persist = as.numeric(stats::quantile(sp_persist, 0.10, names = FALSE)),
    q90_persist = as.numeric(stats::quantile(sp_persist, 0.90, names = FALSE)),
    min_persist = min(sp_persist), max_persist = max(sp_persist),
    frac_gt_05 = mean(sp_persist > 0.5),
    total_pu_area_km2 = sum(total_pu_area_km2)
  ), by = .(method, stage, stage_order, x)][order(method, stage_order)]
}

compare_stage73_stage_summaries <- function(stage_summary) {
  wide <- data.table::dcast(
    stage_summary,
    stage + stage_order + x ~ method,
    value.var = c(
      "mean_persist", "median_persist", "frac_gt_05", "total_pu_area_km2"
    )
  )
  wide[, `:=`(
    mean_persist_diff_pipe_minus_rank = mean_persist_pipe - mean_persist_rank,
    median_persist_diff_pipe_minus_rank = median_persist_pipe - median_persist_rank,
    total_pu_area_diff_pipe_minus_rank =
      total_pu_area_km2_pipe - total_pu_area_km2_rank
  )]
  data.table::setorder(wide, stage_order)
  wide[]
}

stage73_parameter_tables <- function(bundle, config) {
  base <- prepare_stage71_bundle_parameters(bundle, list(
    curve = config$curve, taxa_tag = config$taxa_tag,
    sdm = config$sdm, contract = config$contract
  ))
  curves <- prepare_stage71_curve_parameters(
    base, bundle$species_params, config$curve
  )
  configured <- curves[curve == config$curve][order(parameter_order)]
  species_parameters <- configured[, .(
    scientificName, species,
    className = ifelse(taxon == "mammals", "Mammalia", "Aves"),
    density, a_pred = alpha, b_pred = beta, c_th = threshold_area_km2
  )]
  curve_parameters <- curves[, .(
    scientificName, species,
    className = ifelse(taxon == "mammals", "Mammalia", "Aves"),
    density, c_th = threshold_area_km2, curve, alpha, beta
  )]
  list(species = species_parameters, curves = curve_parameters)
}

run_stage73_comparison <- function(config) {
  validate_stage73_inputs(config)
  method_config <- benchmark_method_config(list(
    taxa = config$taxa, sdm = config$sdm,
    cells_to_remove_per_iteration = config$cells_to_remove_per_iteration,
    pruning_iterations_per_stage = config$pruning_iterations_per_stage,
    rank_method = config$rank_method, optimization_curves = config$curve
  ), config$project_paths, config$contract)
  states <- read_completed_stage6_states(method_config)
  plan <- build_benchmark_target_plan(states)
  assert(benchmark_library_complete(method_config, plan),
         "The detailed figure requires a complete exact-target lookup library.")
  state_index <- data.table::copy(states[[config$curve]])
  state_index[, x := pct_cells_removed]

  bundle <- readRDS(config$paths$initialization_bundle)
  parameters <- stage73_parameter_tables(bundle, config)
  baseline_patch_table <- bundle$patch_table
  rm(bundle)
  missing_focal <- setdiff(config$focal_species, parameters$species$scientificName)
  assert(!length(missing_focal), paste0(
    "Focal species are absent from the shared Stage 6 initialization: ",
    paste(missing_focal, collapse = ", "), "."
  ))
  pipeline_stages <- state_index[, .(
    method = "pipe", stage, stage_order, x,
    path = patch_lookup_path
  )]
  benchmark_stages <- state_index[, .(
    method = "rank", stage, stage_order, x,
    path = benchmark_lookup_file(
      config$paths$lookup_library, as.integer(keep_n)
    )
  )]
  load_pipeline <- function(stage, path) standardize_rank_lut(
    if (stage == 0L) baseline_patch_table else data.table::fread(path),
    "pipe", stage, parameters$species
  )
  load_benchmark <- function(stage, path) {
    target_stage <- as.integer(stage)
    retained_cells <- state_index[stage == target_stage, keep_n][[1L]]
    standardize_rank_lut(
      read_benchmark_lookup(
        method_config, retained_cells, validate_manifest = FALSE
      )$patch_table,
      "rank", target_stage, parameters$species
    )
  }
  pipeline <- compute_stage73_method_trajectory(
    pipeline_stages, "pipe", load_pipeline, parameters$species
  )
  benchmark <- compute_stage73_method_trajectory(
    benchmark_stages, "rank", load_benchmark, parameters$species
  )
  population_units <- data.table::rbindlist(
    list(pipeline$pu_long, benchmark$pu_long), use.names = TRUE, fill = TRUE
  )
  species <- data.table::rbindlist(
    list(pipeline$sp_long, benchmark$sp_long), use.names = TRUE
  )
  expected_rows <- 2L * nrow(state_index) * nrow(parameters$species)
  assert(nrow(species) == expected_rows &&
           !anyDuplicated(species[, .(method, stage, species)]) &&
           all(is.finite(species$sp_persist) & species$sp_persist >= 0 &
                 species$sp_persist <= 1),
         "Stage 7.3 species trajectories are incomplete or invalid.")
  stage_summary <- summarize_stage73_species_trajectory(species)
  list(
    config = config, pu_long = population_units, sp_long = species,
    stage_sum = stage_summary,
    stage_compare = compare_stage73_stage_summaries(stage_summary),
    stage_meta = state_index, species_params = parameters$species,
    curve_parameters = parameters$curves,
    rank_label = config$rank_label
  )
}
