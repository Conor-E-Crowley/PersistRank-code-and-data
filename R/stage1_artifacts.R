# Stage 1 posterior-artifact contracts and operational preflight.
#
# Loaded by stage1_workflow.R after pure configuration. Sourcing is definition
# only; explicit reads validate posterior tables, and preflight alone checks
# input files and prepares output directories. Model fitting remains deferred.

stage1_posterior_specs <- function(paths, bird_model) {
  list(
    mammal_growth = list(
      path = paths$clean$mammal_growth_posterior,
      cols = c("alpha", "beta_logM", "residual_sd")
    ),
    mammal_environmental_variation = list(
      path = paths$clean$mammal_environmental_variation_posterior,
      cols = c("alpha", "beta_logM", "residual_sd")
    ),
    bird_growth = list(
      path = paths$clean$bird_growth_posterior,
      cols = c("alpha", "beta_logGenLength", "residual_sd")
    ),
    bird_environmental_variation = list(
      path = paths$clean$bird_environmental_variation_posterior,
      cols = c(
        "alpha", "beta_logGenLength",
        bird_model$coefficient_names, "residual_sd"
      )
    )
  )
}

validate_stage1_posterior_table <- function(
  draws,
  columns,
  label,
  expected_rows = NULL
) {
  assert(
    identical(names(draws), columns),
    paste0(
      label, " must contain exactly these columns in order: ",
      paste(columns, collapse = ", "), "."
    )
  )
  assert(nrow(draws) > 0L, paste0(label, " contains no posterior draws."))
  if (!is.null(expected_rows)) {
    assert(
      nrow(draws) == expected_rows,
      paste0(label, " contains an unexpected number of posterior draws.")
    )
  }
  assert(
    all(vapply(draws, function(x) is.numeric(x) && all(is.finite(x)), logical(1))),
    paste0(label, " must contain only finite numeric values.")
  )
  assert(
    all(draws$residual_sd > 0),
    paste0(label, " residual_sd must be positive.")
  )
  invisible(TRUE)
}

read_stage1_posteriors <- function(paths, bird_model) {
  specs <- stage1_posterior_specs(paths, bird_model)
  out <- lapply(names(specs), function(model) {
    spec <- specs[[model]]
    need_file(spec$path, paste0(model, " posterior"))
    draws <- readr::read_csv(spec$path, show_col_types = FALSE, progress = FALSE)
    validate_stage1_posterior_table(draws, spec$cols, spec$path)
    draws
  })
  names(out) <- names(specs)
  counts <- vapply(out, nrow, integer(1))
  assert(
    length(unique(counts)) == 1L,
    "The four Stage 1 posterior files must contain the same number of draws."
  )
  attr(out, "posterior_draws") <- unname(counts[[1L]])
  out
}

check_stage1_preflight <- function(config) {
  invisible(lapply(
    unlist(config$paths$raw, use.names = FALSE),
    need_file,
    label = "Stage 1 raw input"
  ))
  ensure_writable_dir(dirname(config$paths$clean$mammal_mass_grid), "Stage 1 data directory")
  ensure_writable_dir(config$paths$figure_dir, "Stage 1 figure directory")
  ensure_writable_dir(config$paths$si_dir, "Stage 1 SI figure directory")
  invisible(TRUE)
}
