# Pure Stage 3 configuration and presentation identifiers.
#
# Loaded first by stage3_workflow.R. It validates parameters and derives paths
# without filesystem access or package loading. Input handoff and preflight are
# owned by stage3_inputs.R; fitting and expensive work remain deferred.


gompertz_model_schema <- function() "persistence_curve_models_v6_dynamic_bird"

required_gompertz_packages <- function(include_figures = TRUE) {
  packages <- c("dplyr", "readr", "tibble")
  if (isTRUE(include_figures)) {
    packages <- c(packages, "ggplot2", "cowplot", "scales")
  }
  packages
}

validate_gompertz_span <- function(x, label = "params$loess_span") {
  validate_scalar_number(
    x,
    label,
    minimum = 0,
    maximum = 1,
    minimum_open = TRUE
  )
}

stage3_figure_files <- function() {
  c(
    gompertz_wolff = "gompertz_wolff.png",
    gompertz_parameters = "S3_gompertz_parameters.png",
    abundance_persistence = "S3_abundance_persistence.png"
  )
}

stage3_selected_figure_ids <- function(taxa) {
  taxa <- validate_taxa_selector(taxa, "taxa")$value
  c(
    if (taxa %in% c("mammals", "both")) "gompertz_wolff",
    "gompertz_parameters",
    "abundance_persistence"
  )
}

stage3_branch_figure_dimensions <- function(
  include_mammals,
  bird_model_groups = character(),
  height
) {
  assert(
    is.logical(include_mammals) &&
      length(include_mammals) == 1L &&
      !is.na(include_mammals),
    "include_mammals must be TRUE or FALSE."
  )
  bird_model_groups <- as.character(bird_model_groups)
  assert(
    all(!is.na(bird_model_groups) & nzchar(trimws(bird_model_groups))) &&
      !anyDuplicated(bird_model_groups),
    "bird_model_groups must contain unique non-empty names."
  )
  height <- validate_scalar_number(
    height, "height", minimum = 0, minimum_open = TRUE
  )
  branch_count <- as.integer(include_mammals) + length(bird_model_groups)
  assert(branch_count > 0L, "At least one Stage 3 figure branch is required.")

  # Shared outer spacing plus 3.2 inches per branch keeps composite panels
  # readable without introducing a separate size setting for each taxon mix.
  c(width = 0.8 + 3.2 * branch_count, height = height)
}

# Validate Stage 3 selectors and derive model/figure paths without reading
# artifacts or loading packages. The workflow owns preflight and active work.
stage3_config <- function(params, paths = project_paths()) {
  mode <- validate_scalar_choice(params$mode, c("reuse", "fit"), "params$mode")
  selection <- validate_taxa_selector(params$taxa, "params$taxa")
  taxa <- selection$value
  write_figures <- validate_scalar_logical(
    params$write_figures %||% TRUE, "params$write_figures"
  )

  curves <- persistence_curves()
  results_dir <- paths$results
  clean_dir <- paths$clean
  figure_dir <- paths$figures
  si_dir <- paths$si_figures
  figure_files <- stage3_figure_files()
  figure_ids <- stage3_selected_figure_ids(taxa)
  paths <- list(
    data_dir = paths$data,
    mammals_points = file.path(results_dir, "persistence_points_mammals.csv"),
    birds_points = file.path(results_dir, "persistence_points_birds.csv"),
    loess_models_rds = file.path(clean_dir, "persistence_curve_models.rds"),
    figure_dir = figure_dir,
    si_dir = si_dir,
    selected_figures = stats::setNames(
      file.path(
        ifelse(figure_ids == "gompertz_wolff", figure_dir, si_dir),
        figure_files[figure_ids]
      ),
      figure_ids
    )
  )

  list(
    mode = mode,
    taxa = taxa,
    selected_mammals = selection$selected_mammals,
    selected_birds = selection$selected_birds,
    write_figures = write_figures,
    curves = curves,
    curve_probabilities = persistence_quantiles()[curves],
    main_curve = main_persistence_curve(),
    heatmap_curve = main_persistence_curve(),
    # Stage 2 is the source of truth for these two simulation settings.
    # resolve_stage3_input_settings() fills them after input validation.
    k0 = NULL,
    expected_horizon = NULL,
    loess = list(
      span = validate_gompertz_span(params$loess_span),
      z = 1.96,
      degree = 2L,
      family = "gaussian",
      prediction_points = 700L
    ),
    paths = paths
  )
}
