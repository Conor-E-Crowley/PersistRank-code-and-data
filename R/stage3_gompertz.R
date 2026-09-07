# Shifted-Gompertz fitting for Stage 3 persistence curves.
#
# Loaded by stage3_workflow.R after configuration. It requires validated Stage
# 2 probability/abundance inputs. Sourcing defines functions only; optimization
# occurs only when fitting is explicitly requested.

clamp01 <- function(p) pmin(pmax(p, 1e-8), 1 - 1e-8)

predict_p <- function(alpha, beta, K, K0) {
  exp(-as.numeric(alpha) * pmax(as.numeric(K) - as.numeric(K0), 1e-12)^(-as.numeric(beta)))
}

validate_gompertz_fit_inputs <- function(K, p_eval, K0, label = "Gompertz fit") {
  K <- as.numeric(K); p_eval <- as.numeric(p_eval); K0 <- as.numeric(K0)
  assert(length(K) == length(p_eval) && length(K) >= 3L, paste0(label, " requires at least three paired points."))
  assert(length(K0) == 1L && is.finite(K0) && K0 > 0, paste0(label, " requires one positive finite Nq."))
  assert(all(is.finite(K) & K > K0), paste0(label, " contains nonfinite K or K <= Nq."))
  assert(all(diff(K) > 0), paste0(label, " requires strictly increasing K."))
  assert(all(is.finite(p_eval) & p_eval >= 0 & p_eval <= 1), paste0(label, " contains persistence outside [0, 1]."))
  assert(all(diff(p_eval) >= -sqrt(.Machine$double.eps)), paste0(label, " requires nondecreasing persistence."))
  invisible(TRUE)
}

failed_gompertz_fit <- function(n_pts, message) {
  tibble::tibble(
    alpha = NA_real_, beta = NA_real_, fit_ok = FALSE, n_pts = as.integer(n_pts),
    rmse_lin = NA_real_, rmse_prob = NA_real_, r2_prob = NA_real_, max_abs_error = NA_real_,
    fit_msg = as.character(message)
  )
}

fit_gompertz_ab <- function(K, p_eval, K0) {
  validation <- tryCatch(
    validate_gompertz_fit_inputs(K, p_eval, K0),
    error = function(e) e
  )
  if (inherits(validation, "error")) return(failed_gompertz_fit(length(K), conditionMessage(validation)))

  df <- tibble::tibble(K = as.numeric(K), p_raw = as.numeric(p_eval)) |>
    dplyr::mutate(p = clamp01(.data$p_raw))
  x <- df$K - K0
  y <- log(-log(df$p))
  lx <- log(x)
  lm0 <- stats::lm(y ~ lx)
  beta0 <- -as.numeric(stats::coef(lm0)[["lx"]])
  alpha0 <- exp(as.numeric(stats::coef(lm0)[["(Intercept)"]]))
  if (!is.finite(beta0) || beta0 <= 0) beta0 <- 1
  if (!is.finite(alpha0) || alpha0 <= 0) alpha0 <- 1

  formula <- p ~ exp(-exp(log_alpha) * pmax(K - K0, 1e-12)^(-exp(log_beta)))
  control <- stats::nls.control(maxiter = 400, tol = 1e-7, minFactor = 1 / 2048, warnOnly = TRUE)
  middle <- which.min(abs(df$p - 0.5))
  beta_starts <- unique(pmax(1e-4, pmin(80, beta0 * c(0.6, 0.85, 1, 1.25, 1.7, 2.5, 4))))
  messages <- character()
  fit <- NULL

  for (beta_start in beta_starts) {
    alpha_start <- -log(df$p[middle]) * (df$K[middle] - K0)^beta_start
    if (!is.finite(alpha_start) || alpha_start <= 0) alpha_start <- alpha0
    attempt <- tryCatch(
      suppressWarnings(stats::nls(
        formula,
        data = df,
        start = list(log_alpha = log(alpha_start), log_beta = log(beta_start)),
        algorithm = "port",
        control = control
      )),
      error = function(e) e
    )
    if (inherits(attempt, "error")) {
      messages <- c(messages, conditionMessage(attempt))
      next
    }
    if (is.null(attempt$convInfo) || !isTRUE(attempt$convInfo$isConv)) {
      messages <- c(messages, attempt$convInfo$stopMessage %||% "NLS did not report convergence")
      next
    }
    fit <- attempt
    break
  }

  if (is.null(fit)) {
    detail <- unique(messages)
    return(failed_gompertz_fit(nrow(df), paste(c("NLS failed to converge", head(detail, 2L)), collapse = ": ")))
  }

  coefficients <- stats::coef(fit)
  alpha_hat <- exp(as.numeric(coefficients[["log_alpha"]]))
  beta_hat <- exp(as.numeric(coefficients[["log_beta"]]))
  p_hat <- predict_p(alpha_hat, beta_hat, df$K, K0)
  finite_fit <- is.finite(alpha_hat) && alpha_hat > 0 && is.finite(beta_hat) && beta_hat > 0 &&
    all(is.finite(p_hat) & p_hat >= 0 & p_hat <= 1)
  if (!finite_fit) return(failed_gompertz_fit(nrow(df), "NLS returned invalid coefficients or predictions"))

  residual <- df$p_raw - p_hat
  total_ss <- sum((df$p_raw - mean(df$p_raw))^2)
  tibble::tibble(
    alpha = alpha_hat,
    beta = beta_hat,
    fit_ok = TRUE,
    n_pts = nrow(df),
    rmse_lin = sqrt(mean((log(-log(df$p)) - log(-log(clamp01(p_hat))))^2)),
    rmse_prob = sqrt(mean(residual^2)),
    r2_prob = if (is.finite(total_ss) && total_ss > 0) 1 - sum(residual^2) / total_ss else NA_real_,
    max_abs_error = max(abs(residual)),
    fit_msg = ""
  )
}

bind_stage3_points <- function(inputs) {
  dplyr::bind_rows(
    if (!is.null(inputs$mammals)) inputs$mammals$points else tibble::tibble(),
    if (!is.null(inputs$birds)) inputs$birds$points else tibble::tibble()
  )
}

fit_all_gompertz_parameters <- function(points, k0) {
  need_cols(
    points,
    c("group", "predictor", "bird_sigma_model_group", "trait_idx", "trait_value", "curve", "K", "p_eval"),
    "validated Stage 3 persistence points"
  )
  points |>
    dplyr::arrange(.data$group, .data$bird_sigma_model_group, .data$trait_idx, .data$curve, .data$K) |>
    dplyr::group_by(
      .data$group, .data$predictor, .data$bird_sigma_model_group,
      .data$trait_idx, .data$trait_value, .data$curve
    ) |>
    dplyr::group_modify(~fit_gompertz_ab(.x$K, .x$p_eval, K0 = k0)) |>
    dplyr::ungroup()
}

format_gompertz_failures <- function(fits) {
  failed <- fits |> dplyr::filter(!.data$fit_ok)
  if (!nrow(failed)) return(character())
  apply(failed, 1L, function(row) {
    paste0(
      "- taxon=", row[["group"]],
      "; trait_idx=", row[["trait_idx"]],
      "; trait_value=", row[["trait_value"]],
      if (!is.na(row[["bird_sigma_model_group"]]) && nzchar(row[["bird_sigma_model_group"]])) paste0("; bird_group=", row[["bird_sigma_model_group"]]) else "",
      "; curve=", row[["curve"]],
      "; n=", row[["n_pts"]],
      "; error=", row[["fit_msg"]]
    )
  })
}

assert_all_gompertz_fits <- function(fits) {
  failed <- fits |> dplyr::filter(!.data$fit_ok)
  if (nrow(failed)) {
    stop(paste0(
      "Shifted-Gompertz fitting failed for ", nrow(failed), " block(s):\n",
      paste(format_gompertz_failures(failed), collapse = "\n"),
      "\nNo LOESS models or output files were produced."
    ), call. = FALSE)
  }
  numeric_columns <- c("alpha", "beta", "rmse_lin", "rmse_prob", "r2_prob", "max_abs_error")
  assert(all(vapply(fits[numeric_columns], function(x) all(is.finite(x)), logical(1))),
         "Successful Gompertz fits contain nonfinite coefficients or diagnostics.")
  assert(all(fits$alpha > 0 & fits$beta > 0), "Successful Gompertz coefficients must be positive.")
  invisible(TRUE)
}

# ---- Trait smoothing on transformed scales ------------------------------------
#
# LOESS uses log10(trait) as x and log(parameter) as y. Predictions are
# back-transformed to positive alpha/beta values, with the configured z value
# used only for display intervals.
