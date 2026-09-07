# Transactional persistence of explicit Stage 1 fit results.
#
# Loaded last by stage1_workflow.R. It requires validated posterior tables,
# trait grids, figures, paths, and project transactions. Sourcing writes
# nothing; a requested commit stages and validates the complete output set.

# ---- Atomic fitted-output transaction -----------------------------------------

write_trait_grid_if_changed <- function(x, path, label) {
  same_grid <- function(a, b) {
    is.data.frame(a) && is.data.frame(b) &&
      identical(names(a), names(b)) &&
      nrow(a) == nrow(b) &&
      isTRUE(all.equal(
        as.data.frame(a),
        as.data.frame(b),
        tolerance = 1e-14,
        check.attributes = FALSE
      ))
  }
  existing <- if (file.exists(path)) {
    tryCatch(
      readr::read_csv(path, show_col_types = FALSE, progress = FALSE),
      error = function(e) NULL
    )
  } else {
    NULL
  }
  if (!is.null(existing) && same_grid(existing, x)) {
    log_msg("Stage 1 | trait grid unchanged | ", label, " | path=", path)
    return(FALSE)
  }
  ensure_writable_dir(dirname(path), paste0(label, " grid directory"))
  temporary <- tempfile(".stage1_grid_", tmpdir = dirname(path), fileext = ".csv")
  on.exit(unlink(temporary), add = TRUE)
  readr::write_csv(x, temporary)
  check <- readr::read_csv(temporary, show_col_types = FALSE, progress = FALSE)
  assert(same_grid(check, x), paste0("Temporary ", label, " grid validation failed."))
  project_file_set_transaction(
    temporary, path, overwrite = TRUE, label = paste0("Stage 1 ", label, " grid")
  )
  log_msg("Stage 1 | trait grid written | ", label, " | path=", path)
  TRUE
}

# Keep every posterior-density column at a stable physical width. The
# one-diet-effect layout is 7.0 × 4.2 inches; additional effect columns expand
# the canvas rather than shrinking text.
stage1_figure_specs <- function(bird_model) {
  validate_bird_model_spec(bird_model)
  diet_count <- length(bird_model$separate_intercepts)
  one_effect_relative_width <- 0.58 + 1 + 1 + 0.82
  requested_relative_width <- 0.58 + 1 + 1 + 0.82 * diet_count
  posterior_width <- 7.0 * requested_relative_width / one_effect_relative_width

  list(
    demographic_calibration = list(width = 5.2, height = 2.85),
    demographic_posteriors = list(
      width = posterior_width,
      height = 4.2
    )
  )
}

# Persist reuse-mode grids and figures in the same order as the Stage 1 Rmd.
# Scientific objects are consumed without modification; only the two changed
# grids and the two requested figures may be written.
write_stage1_reuse_outputs <- function(config, prepared, figures) {
  assert(is.list(config) && is.list(config$paths),
         "write_stage1_reuse_outputs() requires a Stage 1 config.")
  assert(identical(names(prepared), c("mammal", "bird")),
         "Stage 1 reuse trait grids must be named mammal and bird.")
  assert(identical(names(figures), c("demographic_calibration", "demographic_posteriors")),
         "Stage 1 reuse figures have unexpected names or order.")
  paths <- config$paths
  figure_specs <- stage1_figure_specs(config$bird_model)
  write_trait_grid_if_changed(
    prepared$mammal, paths$clean$mammal_mass_grid, "mammal mass"
  )
  write_trait_grid_if_changed(
    prepared$bird, paths$clean$bird_generation_length_grid,
    "bird generation length"
  )
  ggplot2::ggsave(
    paths$figures$demographic_calibration,
    figures$demographic_calibration,
    width = figure_specs$demographic_calibration$width,
    height = figure_specs$demographic_calibration$height,
    dpi = 600, units = "in", bg = "white"
  )
  ggplot2::ggsave(
    paths$figures$demographic_posteriors,
    figures$demographic_posteriors,
    width = figure_specs$demographic_posteriors$width,
    height = figure_specs$demographic_posteriors$height,
    dpi = 600, units = "in", bg = "white"
  )
  invisible(c(
    mammal_mass_grid = paths$clean$mammal_mass_grid,
    bird_generation_length_grid = paths$clean$bird_generation_length_grid,
    demographic_calibration = paths$figures$demographic_calibration,
    demographic_posteriors = paths$figures$demographic_posteriors
  ))
}

# Write every output of an explicit Stage 1 fit as one validated transaction.
# This prevents a failed figure or grid build from leaving a mixed artifact set.
write_stage1_fit_outputs <- function(
  draws_by_model,
  trait_grids,
  figures,
  paths,
  bird_model,
  expected_draws
) {
  specs <- stage1_posterior_specs(paths, bird_model)
  expected_models <- names(specs)
  assert(
    identical(names(draws_by_model), expected_models),
    "Stage 1 posterior set has unexpected model names or order."
  )
  assert(
    identical(names(trait_grids), c("mammal", "bird")),
    "Stage 1 trait grids must be named mammal and bird."
  )
  assert(
    identical(names(figures), c("demographic_calibration", "demographic_posteriors")),
    "Stage 1 figure set has unexpected names or order."
  )

  ensure_writable_dir(
    dirname(specs[[1L]]$path),
    "Stage 1 fitted-output staging directory"
  )
  staging <- tempfile("stage1_outputs_", tmpdir = dirname(specs[[1L]]$path))
  dir.create(staging)
  on.exit(unlink(staging, recursive = TRUE, force = TRUE), add = TRUE)

  staged <- character()
  targets <- character()
  for (model_id in expected_models) {
    columns <- specs[[model_id]]$cols
    validate_stage1_posterior_table(
      draws_by_model[[model_id]], columns, model_id, expected_draws
    )
    staged_path <- file.path(staging, basename(specs[[model_id]]$path))
    readr::write_csv(draws_by_model[[model_id]], staged_path)
    staged_check <- readr::read_csv(
      staged_path, show_col_types = FALSE, progress = FALSE
    )
    validate_stage1_posterior_table(
      staged_check, columns, paste0("staged ", model_id), expected_draws
    )
    staged <- c(staged, staged_path)
    targets <- c(targets, specs[[model_id]]$path)
  }

  grid_specs <- list(
    mammal = list(
      data = trait_grids$mammal,
      path = paths$clean$mammal_mass_grid,
      column = "Mass_g"
    ),
    bird = list(
      data = trait_grids$bird,
      path = paths$clean$bird_generation_length_grid,
      column = "GenLength"
    )
  )
  for (taxon in names(grid_specs)) {
    grid_spec <- grid_specs[[taxon]]
    assert(
      identical(names(grid_spec$data), grid_spec$column) &&
        nrow(grid_spec$data) >= 5L &&
        all(is.finite(grid_spec$data[[grid_spec$column]]) &
              grid_spec$data[[grid_spec$column]] > 0) &&
        all(diff(grid_spec$data[[grid_spec$column]]) > 0),
      paste0("Stage 1 ", taxon, " trait grid is invalid.")
    )
    staged_path <- file.path(staging, basename(grid_spec$path))
    readr::write_csv(grid_spec$data, staged_path)
    staged_check <- readr::read_csv(
      staged_path, show_col_types = FALSE, progress = FALSE
    )
    assert(
      isTRUE(all.equal(
        as.data.frame(staged_check),
        as.data.frame(grid_spec$data),
        tolerance = 1e-14,
        check.attributes = FALSE
      )),
      paste0("Staged Stage 1 ", taxon, " trait grid failed validation.")
    )
    staged <- c(staged, staged_path)
    targets <- c(targets, grid_spec$path)
  }

  figure_specs <- stage1_figure_specs(bird_model)
  for (figure_id in names(figure_specs)) {
    figure_spec <- figure_specs[[figure_id]]
    staged_path <- file.path(
      staging,
      basename(paths$figures[[figure_id]])
    )
    ggplot2::ggsave(
      staged_path,
      figures[[figure_id]],
      width = figure_spec$width,
      height = figure_spec$height,
      dpi = 600,
      units = "in",
      bg = "white"
    )
    assert(
      file.exists(staged_path) && file.info(staged_path)$size > 0,
      paste0("Staged Stage 1 figure is empty: ", figure_id, ".")
    )
    staged <- c(staged, staged_path)
    targets <- c(targets, paths$figures[[figure_id]])
  }

  project_file_set_transaction(
    staged,
    targets,
    overwrite = TRUE,
    label = "Stage 1 fitted output"
  )
  log_msg("Stage 1 | fitted output set committed | files=", length(targets))
  invisible(TRUE)
}
