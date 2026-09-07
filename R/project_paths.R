# Repository roots and the common path interface used by stage configurations.
#
# Manuscript Rmds are rendered from the repository root and do not expose
# internal data or figure directories as YAML parameters. The optional root
# argument exists for tests and programmatic use; it is not a public analysis
# setting. Scientific Rmd execution replaces the generic clean/results members
# below with the model or application hierarchy built by application_storage.R;
# no production Rmd writes a second repository-wide scientific output tree.

# Resolve one public input path relative to the repository root.
#
# This helper normalizes path syntax only. It does not require the target to
# exist, distinguish files from directories, calculate fingerprints, or write
# anything. Callers retain ownership of those stage-specific checks.
project_input_path <- function(value, default = NULL, root, label = "path") {
  if (is.null(value) || !length(value)) value <- default
  value <- validate_path_param(value, label)
  root <- normalizePath(
    validate_path_param(root, "root"), winslash = "/", mustWork = FALSE
  )
  absolute <- grepl("^/|^[A-Za-z]:[/\\\\]", value)
  if (!absolute) value <- file.path(root, value)
  normalizePath(value, winslash = "/", mustWork = FALSE)
}


project_paths <- function(root = ".") {
  assert(
    is.character(root) && length(root) == 1L && !is.na(root) && nzchar(trimws(root)),
    "root must be one non-empty directory path."
  )
  root <- normalizePath(root, winslash = "/", mustWork = FALSE)

  data <- file.path(root, "Data")
  raw <- file.path(data, "Raw")
  clean <- file.path(data, "Clean")
  results <- file.path(data, "Results")
  checkpoints <- file.path(results, "Checkpoints")
  figures <- file.path(root, "Figures")
  si_figures <- file.path(figures, "SI")
  zonation <- file.path(clean, "ZonationOutputs")
  list(
    root = root,
    data = data,
    raw = raw,
    clean = clean,
    results = results,
    checkpoints = checkpoints,
    figures = figures,
    si_figures = si_figures,
    rcpp_cache = file.path(data, "Cache", "Rcpp"),
    cache = file.path(data, "Cache"),
    iucn_cache = file.path(data, "Cache", "iucn_habitat_cache.rds"),
    models = file.path(data, "Models"),
    persistence_models = file.path(data, "Models", "Persistence"),
    applications = file.path(data, "Applications"),
    model_figures = file.path(figures, "Models"),
    application_figures = file.path(figures, "Applications"),
    patches = file.path(clean, "Patches"),
    binary_patches = file.path(clean, "Patches_binary"),
    zonation = zonation,
    zonation_feature_list =
      file.path(zonation, "zonation_binary_patch_feature_list.txt"),
    stage5_checkpoints = file.path(checkpoints, "stage5_patches"),
    priority_runs = file.path(results, "PriorityRuns")
  )
}
