# JAGS fitting and diagnostics for Stage 1 demographic models.
#
# Loaded by stage1_workflow.R after validated calibration inputs. Sourcing only
# defines the model string and fitting helpers; JAGS is checked and invoked only
# for an explicit fit.

jags_lm_model <- "
model {
  for (i in 1:N) {
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- alpha + inprod(beta[1:P], X[i,1:P])
  }
  alpha ~ dnorm(prior_mean, tau_beta)
  for (j in 1:P) { beta[j] ~ dnorm(prior_mean, tau_beta) }
  sigma ~ dunif(sigma_min, sigma_max)
  tau <- pow(sigma, -2)
}
"

# ---- Fitting prerequisites and shared diagnostics -----------------------------

assert_jags_available <- function(fit_models) {
  if (!isTRUE(fit_models)) return(invisible(TRUE))

  if (!requireNamespace("rjags", quietly = TRUE)) {
    stop(
      "Stage 1 mode is 'fit', but rjags is unavailable. Install rjags and ",
      "JAGS, or use mode: 'reuse' with existing posterior CSVs.",
      call. = FALSE
    )
  }
  if (!requireNamespace("coda", quietly = TRUE)) {
    stop(
      "Stage 1 mode is 'fit', but coda is unavailable. Install coda or use ",
      "mode: 'reuse'.",
      call. = FALSE
    )
  }
  if (!tryCatch(length(rjags::jags.version()) > 0L, error = function(e) FALSE)) {
    stop(
      "Stage 1 mode is 'fit', but rjags cannot find a working JAGS installation.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

validate_jags_regression_inputs <- function(y, X, coef_names, model_id) {
  X <- as.matrix(X)
  assert(length(y) > 1L, paste0(model_id, " requires at least two observations."))
  assert(nrow(X) == length(y), paste0(model_id, " response and design-matrix rows differ."))
  assert(ncol(X) == length(coef_names), paste0(model_id, " coefficient names do not match its predictors."))
  assert(
    identical(colnames(X), coef_names),
    paste0(model_id, " design-matrix column names must match its coefficient names.")
  )
  assert(all(is.finite(y)), paste0(model_id, " response contains non-finite values."))
  assert(all(is.finite(X)), paste0(model_id, " design matrix contains non-finite values."))
  assert(!anyDuplicated(coef_names), paste0(model_id, " coefficient names must be unique."))
  assert(
    qr(cbind(`(Intercept)` = 1, X))$rank == ncol(X) + 1L,
    paste0(model_id, " design matrix is rank deficient when combined with the intercept.")
  )
  invisible(X)
}

initial_residual_sd <- function(bayes) {
  width <- bayes$sigma_max - bayes$sigma_min
  lower <- bayes$sigma_min + 0.1 * width
  upper <- bayes$sigma_max - 0.1 * width
  if (!is.finite(lower) || !is.finite(upper) || lower >= upper) {
    return((bayes$sigma_min + bayes$sigma_max) / 2)
  }
  stats::runif(1L, lower, upper)
}

summarize_mcmc_diagnostics <- function(mcmc, bayes, model_id) {
  rhat <- coda::gelman.diag(
    mcmc,
    autoburnin = FALSE,
    multivariate = FALSE
  )$psrf[, "Point est."]
  ess <- coda::effectiveSize(mcmc)

  assert(all(is.finite(rhat)), paste0(model_id, " produced non-finite R-hat values."))
  assert(all(is.finite(ess)), paste0(model_id, " produced non-finite effective sample sizes."))

  max_rhat <- max(rhat)
  min_ess <- min(ess)
  assert(
    max_rhat <= bayes$rhat_max,
    paste0(model_id, " failed convergence: max R-hat=", signif(max_rhat, 5),
           " exceeds ", bayes$rhat_max, ".")
  )
  assert(
    min_ess >= bayes$ess_min,
    paste0(model_id, " failed convergence: min ESS=", signif(min_ess, 5),
           " is below ", bayes$ess_min, ".")
  )

  list(rhat = rhat, ess = ess, max_rhat = max_rhat, min_ess = min_ess)
}

posterior_interval_table <- function(draws, model_id) {
  rows <- lapply(names(draws), function(parameter) {
    q <- stats::quantile(draws[[parameter]], c(0.025, 0.5, 0.975), names = FALSE)
    tibble::tibble(
      model = model_id,
      parameter = parameter,
      lower_95 = q[[1L]],
      median = q[[2L]],
      upper_95 = q[[3L]]
    )
  })
  dplyr::bind_rows(rows)
}

# Fit y on the log10 scale using an intercept plus the supplied design matrix.
# Returned draws use transparent coefficient names and retain residual SD on
# the fitted log10 scale for optional Stage 2 posterior prediction.
fit_jags_lm <- function(y, X, coef_names, bayes, model_id, seed_offset) {
  X <- validate_jags_regression_inputs(y, X, coef_names, model_id)
  model_seed <- as.integer(bayes$seed + seed_offset)
  set.seed(model_seed)

  data_jags <- list(
    N = length(y),
    P = ncol(X),
    y = as.numeric(y),
    X = X,
    prior_mean = bayes$prior_mean,
    tau_beta = bayes$tau_beta,
    sigma_min = bayes$sigma_min,
    sigma_max = bayes$sigma_max
  )
  inits <- lapply(seq_len(bayes$n_chains), function(chain) {
    list(
      alpha = stats::rnorm(1L),
      beta = stats::rnorm(data_jags$P),
      sigma = initial_residual_sd(bayes),
      .RNG.name = "base::Wichmann-Hill",
      .RNG.seed = model_seed + chain
    )
  })

  started <- Sys.time()
  log_msg(
    "Stage 1 | model fit started | model=", model_id,
    "| n=", data_jags$N,
    "| predictors=", paste(coef_names, collapse = ",")
  )
  connection <- textConnection(jags_lm_model)
  on.exit(close(connection), add = TRUE)
  model <- rjags::jags.model(
    connection,
    data = data_jags,
    inits = inits,
    n.chains = bayes$n_chains,
    n.adapt = bayes$n_adapt,
    quiet = TRUE
  )
  beta_monitors <- paste0("beta[", seq_len(data_jags$P), "]")
  mcmc <- rjags::coda.samples(
    model,
    variable.names = c("alpha", beta_monitors, "sigma"),
    n.iter = bayes$n_iter,
    thin = bayes$thin,
    progress.bar = "none"
  )

  diagnostics <- summarize_mcmc_diagnostics(mcmc, bayes, model_id)
  matrix_draws <- as.data.frame(as.matrix(mcmc), check.names = FALSE)
  draws <- tibble::as_tibble(matrix_draws[, c("alpha", beta_monitors, "sigma"), drop = FALSE])
  names(draws) <- c("alpha", paste0("beta_", coef_names), "residual_sd")
  assert(
    nrow(draws) == bayes$n_chains * bayes$n_iter / bayes$thin,
    paste0(model_id, " produced an unexpected number of posterior draws.")
  )
  assert(
    all(vapply(draws, function(x) is.numeric(x) && all(is.finite(x)), logical(1L))),
    paste0(model_id, " produced non-numeric or non-finite posterior draws.")
  )

  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  log_msg(
    "Stage 1 | model fit complete | model=", model_id,
    "| elapsed_s=", sprintf("%.1f", elapsed),
    "| max_rhat=", sprintf("%.4f", diagnostics$max_rhat),
    "| min_ess=", sprintf("%.0f", diagnostics$min_ess)
  )

  list(
    model_id = model_id,
    draws = draws,
    mcmc = mcmc,
    diagnostics = diagnostics,
    elapsed_seconds = elapsed,
    seed = model_seed
  )
}

# ---- Four demographic allometries ---------------------------------------------

fit_loglog_allometry <- function(x, y, bayes, coef_name, model_id, seed_offset) {
  assert(all(is.finite(x) & x > 0), paste0(model_id, " predictor must be positive and finite."))
  assert(all(is.finite(y) & y > 0), paste0(model_id, " response must be positive and finite."))
  fit_jags_lm(
    y = log10(y),
    X = matrix(log10(x), ncol = 1L, dimnames = list(NULL, coef_name)),
    coef_names = coef_name,
    bayes = bayes,
    model_id = model_id,
    seed_offset = seed_offset
  )
}

fit_bird_sigma_model <- function(
  bird_sigma_data,
  bird_model,
  bayes,
  seed_offset = 4000L
) {
  design <- bird_model_design_matrix(
    bird_sigma_data$GenLength,
    bird_sigma_data$diet5_group,
    bird_model
  )
  fit_jags_lm(
    y = log10(bird_sigma_data$sigma),
    X = design,
    coef_names = colnames(design),
    bayes = bayes,
    model_id = "bird_sigma",
    seed_offset = seed_offset
  )
}

descriptive_stage1_r2 <- function(cal, bird_model) {
  bird_design <- bird_model_design_matrix(
    cal$bird_sigma_data$GenLength,
    cal$bird_sigma_data$diet5_group,
    bird_model
  )
  bird_response <- log10(cal$bird_sigma_data$sigma)
  bird_fit <- stats::lm.fit(
    x = cbind(`(Intercept)` = 1, bird_design),
    y = bird_response
  )
  bird_r2 <- 1 - sum(bird_fit$residuals^2) /
    sum((bird_response - mean(bird_response))^2)
  c(
    mammal_growth =
      summary(stats::lm(log10(rm) ~ log10(Mass_g), data = cal$mammal_rmax))$r.squared,
    mammal_environmental_variation =
      summary(stats::lm(log10(sigma) ~ log10(Mass_g), data = cal$mammal_sigma))$r.squared,
    bird_growth =
      summary(stats::lm(log10(rm) ~ log10(GenLength), data = cal$bird_rmax))$r.squared,
    bird_environmental_variation = bird_r2
  )
}

fit_all_demographic_models <- function(cal, bird_model, bayes) {
  list(
    mammal_growth = fit_loglog_allometry(
      cal$mammal_rmax$Mass_g, cal$mammal_rmax$rm, bayes,
      coef_name = "logM", model_id = "mammal_rm", seed_offset = 1000L
    ),
    mammal_environmental_variation = fit_loglog_allometry(
      cal$mammal_sigma$Mass_g, cal$mammal_sigma$sigma, bayes,
      coef_name = "logM", model_id = "mammal_sigma", seed_offset = 2000L
    ),
    bird_growth = fit_loglog_allometry(
      cal$bird_rmax$GenLength, cal$bird_rmax$rm, bayes,
      coef_name = "logGenLength", model_id = "bird_rm", seed_offset = 3000L
    ),
    bird_environmental_variation = fit_bird_sigma_model(
      cal$bird_sigma_data, bird_model, bayes, seed_offset = 4000L
    )
  )
}

log_posterior_intervals <- function(draws_by_model) {
  for (model_id in names(draws_by_model)) {
    intervals <- posterior_interval_table(draws_by_model[[model_id]], model_id)
    for (i in seq_len(nrow(intervals))) {
      row <- intervals[i, ]
      log_msg(
        "Stage 1 | posterior | model=", model_id,
        "| parameter=", row$parameter,
        "| median=", signif(row$median, 5),
        "| 95%_interval=[", signif(row$lower_95, 5), ",", signif(row$upper_95, 5), "]"
      )
    }
  }
  invisible(TRUE)
}

