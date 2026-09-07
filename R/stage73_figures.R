# Stage 7.3 final figure composition.
#
# The Stage 7.3 workflow loads this definition-only module last. It composes
# already-prepared focal and assemblage components once, preserving the
# manuscript dimensions and avoiding duplicate large plot objects.

stage73_combined_figure_dimensions <- function() {
  list(
    width = 7.0,
    height = 7.45,
    assemblage_height = 2.70,
    focal_height = 4.75
  )
}

build_stage73_combined_figure <- function(
  core,
  focus_species,
  figure_data = NULL,
  distribution_data = NULL,
  assemblage_components = NULL,
  focal_components = NULL,
  assemblage_point_style = "hollow",
  focal_point_style = "hollow",
  assemblage_central_stat = "mean"
) {
  assemblage_point_style <- validate_scalar_choice(
    assemblage_point_style, report_point_styles(), "Stage 7.3 assemblage_point_style"
  )
  focal_point_style <- validate_scalar_choice(
    focal_point_style, report_point_styles(), "Stage 7.3 focal_point_style"
  )
  assemblage_central_stat <- validate_scalar_choice(
    assemblage_central_stat, report_central_statistics(),
    "Stage 7.3 assemblage_central_stat"
  )
  if (is.null(figure_data)) {
    figure_data <- prepare_stage73_figure_data(core)
  }
  if (is.null(distribution_data)) {
    distribution_data <- prepare_stage73_distribution_data(
      core, figure_data,
      central_stat = assemblage_central_stat
    )
  }
  if (is.null(assemblage_components)) {
    assemblage_components <- prepare_stage73_assemblage_components(
      core,
      figure_data,
      distribution_data,
      central_stat = assemblage_central_stat
    )
  }
  if (is.null(focal_components)) {
    focal_components <- prepare_stage73_focal_figure_components(
      core,
      focus_species,
      figure_data,
      point_style = focal_point_style
    )
  }

  assemblage <- stage73_label_assemblage_for_composite(
    stage73_assemblage_section(assemblage_components, point_style = assemblage_point_style)
  )
  focal <- stage73_focal_section(focal_components, point_style = focal_point_style)
  assert(
    data.table::is.data.table(focal$panel_key) &&
      identical(focal$panel_key$tag, paste0("(", letters[2:7], ")")),
    "Stage 7.3 focal panels must be labelled (b)-(g)."
  )
  assert(
    identical(
      focal$minimum_removed_percent,
      stage73_focal_spec()$minimum_removed_percent
    ),
    "Stage 7.3 focal panels must use the canonical late-stage range."
  )

  dimensions <- stage73_combined_figure_dimensions()
  out <- cowplot::plot_grid(
    assemblage,
    focal$plot,
    ncol = 1L,
    align = "v",
    axis = "lr",
    rel_heights = c(
      dimensions$assemblage_height,
      dimensions$focal_height
    )
  )
  attr(out, "stage73_combined_dimensions_inches") <- unlist(
    dimensions[c("width", "height")],
    use.names = TRUE
  )
  attr(out, "stage73_combined_component_heights_inches") <- c(
    assemblage = dimensions$assemblage_height,
    focal = dimensions$focal_height
  )
  attr(out, "stage73_combined_panel_tags") <- paste0("(", letters[1:7], ")")
  attr(out, "stage73_combined_selected_stages") <- focal$selected_stages
  attr(out, "stage73_combined_focal_species") <- focus_species
  out
}
