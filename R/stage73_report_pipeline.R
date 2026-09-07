# In-memory Stage 7.3 report assembly.
#
# The Stage 7.3 workflow loads this after comparison, statistics, report, and
# focused figure modules. It constructs each large report object once, writes
# no artifact, creates no cache, and relies on R copy-on-write for returned
# objects.

# Build the report payload without changing the caller's configuration.
build_stage73_report_result <- function(config) {
  assert(is.list(config), "config must be a validated Stage 7.3-compatible configuration.")
  core <- run_stage73_comparison(config)

  figure_data <- prepare_stage73_figure_data(core)
  distribution_data <- prepare_stage73_distribution_data(
    core, figure_data,
    central_stat = config$assemblage_central_stat
  )
  statistics <- compute_stage73_statistics(
    core,
    config$focal_species[[1L]],
    figure_data = figure_data
  )
  report_tables <- build_stage73_report_tables(core, statistics)

  assemblage_components <- prepare_stage73_assemblage_components(
    core,
    figure_data,
    distribution_data,
    central_stat = config$assemblage_central_stat
  )
  focal_components <- prepare_stage73_focal_figure_components(
    core,
    config$focal_species,
    figure_data,
    point_style = config$focal_point_style
  )
  figure <- build_stage73_combined_figure(
    core = core,
    focus_species = config$focal_species,
    figure_data = figure_data,
    distribution_data = distribution_data,
    assemblage_components = assemblage_components,
    focal_components = focal_components,
    assemblage_point_style = config$assemblage_point_style,
    focal_point_style = config$focal_point_style,
    assemblage_central_stat = config$assemblage_central_stat
  )

  list(
    statistics = statistics,
    report_tables = report_tables,
    focal_species = config$focal_species,
    figure = figure
  )
}
