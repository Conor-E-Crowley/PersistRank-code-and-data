# Stage 5.3 initial-state persistence orchestration.
#
# This workflow owns preflight, package loading, input reads, analysis/figure
# sequencing, and the existing result object. Focused modules own calculations
# and presentation. Sourcing is definition-only and defers every package load.

if (!exists("project_source", mode = "function")) {
  source(file.path("R", "project_loader.R"))
}
project_source(c(
  # Shared contracts, application adapters, and presentation.
  "R/project_utils.R", "R/project_artifacts.R", "R/project_rmd.R",
  "R/project_runtime_log.R", "R/project_transactions.R",
  "R/project_paths.R", "R/project_config.R", "R/demographic_contract.R",
  "R/analysis_contract.R", "R/priority_run_config.R",
  "R/species_table_config.R", "R/patch_contract.R", "R/patch_config.R",
  "R/application_storage.R", "R/storage_artifacts.R", "R/model_lifecycle.R",
  "R/application_state.R",
  "R/application_context.R", "R/application_rmd.R", "R/figure_utils.R",
  "R/stage53_config.R", "R/stage53_analysis.R", "R/stage53_figures.R"
))

run_stage53 <- function(config) {
  assert(
    is.list(config),
    "run_stage53() requires a validated stage53_config() result."
  )
  check_stage53_preflight(config)
  started <- Sys.time()
  load_packages(c("readr", "dplyr", "stringr", "ggplot2", "cowplot", "scales"))
  species_selected <- read_stage53_species(config)

  lookup <- readRDS(config$paths$patch_lookup_rds)
  persistence <- calculate_stage53_persistence(species_selected, lookup)
  represented_names <- persistence$represented_names
  absent <- persistence$absent_species
  pu <- persistence$pu
  species <- persistence$species

  coverage <- species_selected |>
    dplyr::mutate(represented = scientificName %in% represented_names) |>
    dplyr::group_by(className) |>
    dplyr::summarise(
      selected = dplyr::n(),
      represented = sum(represented),
      absent = selected - represented,
      coverage_pct = 100 * represented / selected,
      .groups = "drop"
    )
  pu_summary <- pu |>
    dplyr::group_by(className) |>
    dplyr::summarise(
      species = dplyr::n_distinct(scientificName), n_pu = dplyr::n(),
      median_area_km2 = stats::median(pu_area_km2),
      median_persistence = stats::median(P_pu), .groups = "drop"
    )
  species_summary <- species |>
    dplyr::group_by(className) |>
    dplyr::summarise(
      species = dplyr::n(), median_pus = stats::median(n_pu),
      median_area_km2 = stats::median(total_pu_area_km2),
      median_persistence = stats::median(species_persistence), .groups = "drop"
    )
  correlation_groups <- split(species, species$className)
  if (dplyr::n_distinct(species$className) > 1L) {
    correlation_groups$All <- species
  }
  correlations <- dplyr::bind_rows(lapply(names(correlation_groups), function(label) {
    x <- correlation_groups[[label]]
    rho <- if (nrow(x) >= 3L) suppressWarnings(stats::cor(
      x$total_pu_area_km2, x$species_persistence, method = "spearman"
    )) else NA_real_
    data.frame(className = label, n_species = nrow(x), spearman_rho = rho)
  }))
  ranked <- species |>
    dplyr::select(
      scientificName, className, n_pu,
      total_pu_area_km2, species_persistence
    ) |>
    dplyr::arrange(species_persistence, scientificName)
  top_n <- min(10L, nrow(ranked))
  ranked_species <- list(
    lowest = utils::head(ranked, top_n),
    highest = utils::head(dplyr::arrange(ranked, dplyr::desc(species_persistence), scientificName), top_n)
  )
  figure <- build_stage53_figure(species)
  write_stage53_figure_atomic(figure, config)
  persistence_range <- range(species$species_persistence)
  log_msg(
    "Stage 5.3 complete | taxa=", config$taxa,
    "| curve=", config$curve,
    "| selected=", nrow(species_selected),
    "| represented=", nrow(species),
    "| coverage_gaps=", nrow(absent),
    "| population_units=", nrow(pu),
    "| figure=", config$paths$figure,
    "| persistence_range=", paste(signif(persistence_range, 4), collapse = "-"),
    "| elapsed=", round(as.numeric(difftime(Sys.time(), started, units = "secs")), 2), "s"
  )
  list(
    coverage = coverage,
    absent_species = absent,
    pu = pu,
    species = species,
    pu_summary = pu_summary,
    species_summary = species_summary,
    correlations = correlations,
    ranked_species = ranked_species,
    figure = figure,
    curve = config$curve,
    taxa = config$taxa
  )
}
