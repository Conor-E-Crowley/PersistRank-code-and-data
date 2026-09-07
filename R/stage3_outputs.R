# Transactional persistence of Stage 3 models and figures.
#
# Loaded last by stage3_workflow.R. It requires validated model/figure objects,
# configuration, and project transactions. Sourcing writes nothing; calls stage
# and atomically replace the requested artifact set.

write_stage3_outputs <- function(
  models,
  figures,
  figure_targets,
  config
) {
  assert(
    is.list(figures) && identical(names(figures), names(figure_targets)),
    "Stage 3 figures and figure targets must have identical names."
  )
  assert(
    !anyDuplicated(unname(figure_targets)),
    "Stage 3 figure targets must be unique."
  )

  bird_groups <- if (config$selected_birds) {
    models$meta$bird_model_spec$model_groups
  } else {
    character()
  }
  figure_specs <- list(
    gompertz_wolff = c(width = 6.2, height = 2.85),
    gompertz_parameters = stage3_branch_figure_dimensions(
      config$selected_mammals,
      bird_groups,
      height = 5.4
    ),
    # One row of branch panels is sized as a complete manuscript figure.
    abundance_persistence = stage3_branch_figure_dimensions(
      config$selected_mammals,
      bird_groups,
      height = 3.4
    )
  )
  missing_specs <- setdiff(names(figures), names(figure_specs))
  assert(
    !length(missing_specs),
    paste0(
      "Stage 3 has no output dimensions for figure(s): ",
      paste(missing_specs, collapse = ", "), "."
    )
  )

  staging_dir <- tempfile(
    "stage3_outputs_",
    tmpdir = dirname(config$paths$loess_models_rds)
  )
  dir.create(staging_dir)
  on.exit(unlink(staging_dir, recursive = TRUE, force = TRUE), add = TRUE)

  staged <- character()
  targets <- character()
  if (identical(config$mode, "fit")) {
    staged_model <- file.path(
      staging_dir,
      basename(config$paths$loess_models_rds)
    )
    saveRDS(models, staged_model, version = 3L)
    validate_gompertz_loess_models(readRDS(staged_model))
    staged <- c(staged, staged_model)
    targets <- c(targets, config$paths$loess_models_rds)
  }

  for (id in names(figures)) {
    spec <- figure_specs[[id]]
    staged_figure <- file.path(
      staging_dir,
      basename(figure_targets[[id]])
    )
    ggplot2::ggsave(
      staged_figure,
      figures[[id]],
      width = spec[["width"]],
      height = spec[["height"]],
      dpi = 600,
      units = "in",
      bg = "white"
    )
    assert(
      file.exists(staged_figure) && file.info(staged_figure)$size > 0,
      paste0("Staged Stage 3 figure is empty: ", id, ".")
    )
    staged <- c(staged, staged_figure)
    targets <- c(targets, figure_targets[[id]])
  }

  project_file_set_transaction(
    staged,
    targets,
    overwrite = TRUE,
    label = "Stage 3"
  )
  invisible(targets)
}
