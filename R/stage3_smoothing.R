# Positive LOESS smoothing for Stage 3 Gompertz parameters.
#
# Loaded by stage3_workflow.R after exact Gompertz fitting. It operates on small
# parameter tables and returns fitted/predicted curve sets. Sourcing has no side
# effects and performs no artifact I/O.

fit_loess_logparam <- function(trait_value, parameter, span = 0.5) {
  trait_value <- as.numeric(trait_value); parameter <- as.numeric(parameter)
  assert(length(trait_value) == length(parameter) && length(trait_value) >= 3L,
         "LOESS fitting requires at least three paired trait and parameter values.")
  assert(all(is.finite(trait_value) & trait_value > 0), "LOESS trait values must be positive and finite.")
  assert(all(is.finite(parameter) & parameter > 0), "LOESS parameters must be positive and finite.")
  assert(length(unique(trait_value)) == length(trait_value), "LOESS trait values must be unique within a curve.")
  data <- data.frame(logTrait = log10(trait_value), logParameter = log(parameter))
  stats::loess(
    logParameter ~ logTrait,
    data = data,
    span = span,
    degree = 2,
    family = "gaussian",
    control = stats::loess.control(surface = "interpolate")
  )
}

predict_loess_exp <- function(fit, trait_value, z = 1.96) {
  trait_value <- as.numeric(trait_value)
  assert(all(is.finite(trait_value) & trait_value > 0), "LOESS prediction trait values must be positive and finite.")
  prediction <- stats::predict(fit, newdata = data.frame(logTrait = log10(trait_value)), se = TRUE)
  out <- tibble::tibble(
    trait_value = trait_value,
    mid = exp(as.numeric(prediction$fit)),
    lo = exp(as.numeric(prediction$fit) - z * as.numeric(prediction$se.fit)),
    hi = exp(as.numeric(prediction$fit) + z * as.numeric(prediction$se.fit))
  )
  assert(all(vapply(out, function(x) all(is.finite(x)), logical(1))) && all(out$mid > 0 & out$lo > 0 & out$hi > 0),
         "LOESS produced nonfinite or nonpositive predictions within its training range.")
  out
}

fit_curve_loess <- function(data, span) {
  need_cols(data, c("trait_value", "alpha", "beta"), "Gompertz coefficients for LOESS")
  assert(all(is.finite(data$trait_value) & data$trait_value > 0), "LOESS input contains invalid trait values.")
  assert(all(is.finite(data$alpha) & data$alpha > 0 & is.finite(data$beta) & data$beta > 0),
         "LOESS input contains invalid Gompertz coefficients.")
  list(
    alpha = fit_loess_logparam(data$trait_value, data$alpha, span),
    beta = fit_loess_logparam(data$trait_value, data$beta, span),
    n = nrow(data),
    trait_range = range(data$trait_value)
  )
}

fit_loess_set <- function(data, curves, span) {
  assert(identical(sort(unique(as.character(data$curve))), sort(as.character(curves))),
         "LOESS input does not contain the complete curve set.")
  stats::setNames(lapply(curves, function(curve) {
    block <- data[data$curve == curve, , drop = FALSE]
    fit_curve_loess(block, span)
  }), curves)
}

predict_one_curve <- function(fit, trait_grid, curve, z) {
  alpha <- predict_loess_exp(fit$alpha, trait_grid, z) |>
    dplyr::transmute(.data$trait_value, curve = curve, alpha_mid = .data$mid, alpha_lo = .data$lo, alpha_hi = .data$hi)
  beta <- predict_loess_exp(fit$beta, trait_grid, z) |>
    dplyr::transmute(.data$trait_value, curve = curve, beta_mid = .data$mid, beta_lo = .data$lo, beta_hi = .data$hi)
  dplyr::left_join(alpha, beta, by = c("trait_value", "curve"))
}

loess_prediction_grid <- function(trait_range, n_grid) {
  # The interpolated LOESS surface can return NA exactly on its numerical hull.
  # Move display-grid endpoints one machine step inside the training range.
  trait_range <- as.numeric(trait_range) * c(
    1 + .Machine$double.eps,
    1 - .Machine$double.eps
  )
  log_space(trait_range[1], trait_range[2], n_grid)
}

predict_loess_set <- function(data, fits, curves, z, n_grid = 700L) {
  assert(nrow(data) > 0L, "Cannot create a LOESS prediction grid from an empty table.")
  trait_grid <- loess_prediction_grid(range(data$trait_value), n_grid)
  dplyr::bind_rows(lapply(curves, function(curve) predict_one_curve(fits[[curve]], trait_grid, curve, z)))
}

validate_loess_curve_set <- function(models, curves, z, label) {
  assert(is.list(models) && identical(names(models), curves),
         paste0(label, " must contain the ordered curve set: ", paste(curves, collapse = ", "), "."))
  for (curve in curves) {
    fit <- models[[curve]]
    assert(is.list(fit) && inherits(fit$alpha, "loess") && inherits(fit$beta, "loess") &&
             length(fit$n) == 1L && is.finite(fit$n) && fit$n > 0 && fit$n == floor(fit$n) &&
             length(fit$trait_range) == 2L && all(is.finite(fit$trait_range)) &&
             all(fit$trait_range > 0) && fit$trait_range[[1L]] < fit$trait_range[[2L]],
           paste0(label, " has no usable alpha/beta LOESS fit for ", curve, "."))
    grid <- loess_prediction_grid(fit$trait_range, 101L)
    predict_loess_exp(fit$alpha, grid, z = z)
    predict_loess_exp(fit$beta, grid, z = z)
  }
  invisible(TRUE)
}

# ---- Saved-model contract and provenance --------------------------------------
