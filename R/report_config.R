# Shared lightweight presentation configuration for Stages 7.3, 8.1, and 9.
#
# These validators operate only on public scalar/vector settings. They read no
# manifests or scientific artifacts, attach no packages, and write nothing.

validate_report_persistence_threshold <- function(
  value, label = "persistence_threshold"
) {
  value <- as.numeric(value)
  assert(length(value) == 1L && is.finite(value) && value > 0 && value < 1,
         paste0(label, " must be one value in (0,1)."))
  value
}

report_point_styles <- function() c("solid", "hollow", "none")

report_central_statistics <- function() c("mean", "median")

normalize_report_focal_species <- function(x, label = "focal_species") {
  values <- trimws(as.character(unlist(x, use.names = FALSE)))
  assert(length(values) == 2L && all(nzchar(values)) && !anyDuplicated(values),
         paste0(label, " must contain exactly two unique scientific names."))
  values
}

report_presentation_settings <- function(params) list(
  persistence_threshold = validate_report_persistence_threshold(
    params$persistence_threshold
  ),
  focal_species = normalize_report_focal_species(params$focal_species),
  assemblage_point_style = validate_scalar_choice(
    params$assemblage_point_style, report_point_styles(),
    "params$assemblage_point_style"
  ),
  focal_point_style = validate_scalar_choice(
    params$focal_point_style, report_point_styles(), "params$focal_point_style"
  ),
  assemblage_central_stat = validate_scalar_choice(
    params$assemblage_central_stat, report_central_statistics(),
    "params$assemblage_central_stat"
  )
)
