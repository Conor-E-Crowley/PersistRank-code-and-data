# Simulation helpers for persistence-point generation.
#
# The public settings arrive through one validated `sim`/`grid` contract.
# Expensive C++ evaluations are memoized by K for the active trait/curve block;
# each cache miss is checkpointed immediately. Completed blocks are appended
# in canonical trait × curve order and validated before final promotion.


# ---- Demographic parameter samplers -------------------------------------------
#
# Stage 1 coefficients are on log10 scales. In posterior-predictive mode the
# stored residual SD multiplies a deterministic standard-normal deviate before
# back-transformation. Coefficient-only mode never constructs or consumes
# residual deviates.

make_mammal_sampler <- function(inputs, sim) {
  function(trait_row) {
    mass_g <- as.numeric(trait_row$Mass_g[[1L]])
    logM <- log10(mass_g)
    r_eta <-
      inputs$post_rm$alpha[inputs$idx$rm] +
        inputs$post_rm$beta_logM[inputs$idx$rm] * logM
    sigma_eta <-
      inputs$post_sigma$alpha[inputs$idx$sigma] +
        inputs$post_sigma$beta_logM[inputs$idx$sigma] * logM
    if (identical(sim$demographic_uncertainty, "posterior_predictive")) {
      r_eta <- r_eta +
        inputs$post_rm$residual_sd[inputs$idx$rm] * inputs$residual$growth
      sigma_eta <- sigma_eta +
        inputs$post_sigma$residual_sd[inputs$idx$sigma] *
          inputs$residual$environmental_variation
    }
    r <- sim$r_buffer * 10^r_eta
    s <- 10^sigma_eta
    list(r = as.numeric(r), sigma = as.numeric(s))
  }
}

make_bird_sampler <- function(inputs, sim) {
  function(trait_row) {
    GenLength <- as.numeric(trait_row$GenLength[[1L]])
    model_group <- validate_bird_model_group(
      trait_row$bird_sigma_model_group[[1L]],
      inputs$bird_model,
      "bird simulation model group"
    )
    logG <- log10(GenLength)
    r_eta <-
      inputs$post_rm$alpha[inputs$idx$rm] +
        inputs$post_rm$beta_logGenLength[inputs$idx$rm] * logG
    sigma_eta <-
      inputs$post_sigma$alpha[inputs$idx$sigma] +
        inputs$post_sigma$beta_logGenLength[inputs$idx$sigma] * logG
    for (group in inputs$bird_model$separate_intercepts) {
      coefficient <- bird_diet_coefficient_name(group)
      sigma_eta <- sigma_eta +
        inputs$post_sigma[[coefficient]][inputs$idx$sigma] *
          as.numeric(model_group == group)
    }
    if (identical(sim$demographic_uncertainty, "posterior_predictive")) {
      r_eta <- r_eta +
        inputs$post_rm$residual_sd[inputs$idx$rm] * inputs$residual$growth
      sigma_eta <- sigma_eta +
        inputs$post_sigma$residual_sd[inputs$idx$sigma] *
          inputs$residual$environmental_variation
    }
    r <- sim$r_buffer * 10^r_eta
    s <- 10^sigma_eta
    list(r = as.numeric(r), sigma = as.numeric(s))
  }
}

# ---- Abundance search and final-grid construction -----------------------------

# Expand the upper K bound by doubling until the target is bracketed, then
# bisect. The supplied evaluator reuses one CRN context throughout the search.
find_K_for_target <- function(p_target, get_q_at_K, low0, high0, rel_tol,
                              max_k = Inf, context = "K search") {
  low <- as.numeric(low0)
  high <- as.numeric(high0)
  if (high < low) {
    stop("Bisection search requires high0 >= low0.", call. = FALSE)
  }

  q_high <- get_q_at_K(high)
  assert(is.finite(q_high), paste0(context, " returned non-finite persistence at K=", high, "."))
  while (q_high < p_target) {
    if (high >= max_k) {
      stop(
        context, " could not bracket p_target=", p_target,
        " before k_search_max_k=", format(max_k, scientific = FALSE),
        "; last persistence=", signif(q_high, 6), ".",
        call. = FALSE
      )
    }
    low <- high
    high <- min(high * 2, max_k)
    q_high <- get_q_at_K(high)
    assert(is.finite(q_high), paste0(context, " returned non-finite persistence at K=", high, "."))
  }
  hi_bracket <- high

  while ((high - low) > rel_tol * high) {
    mid <- (low + high) / 2
    q_mid <- get_q_at_K(mid)
    assert(is.finite(q_mid), paste0(context, " returned non-finite persistence at K=", mid, "."))
    if (q_mid >= p_target) high <- mid else low <- mid
  }

  list(K = (low + high) / 2, hi_bracket = hi_bracket)
}

# The Gompertz fit and inversion only construct a monotone K grid. Reported
# persistence values are still evaluated by the compiled stochastic simulator.
fit_anchor_gompertz_ab <- function(K_anchor, p_anchor, sim) {
  x <- pmax(as.numeric(K_anchor) - sim$quasi_extinction_abundance, 1e-12)
  p <- as.numeric(p_anchor)
  assert(
    length(x) == length(p) &&
      length(p) >= 3L &&
      all(is.finite(x)) &&
      all(is.finite(p)) &&
      all(p > 0 & p < 1),
    "Gompertz prefit anchors must contain matching finite K and probability values in (0, 1)."
  )

  y <- log(-log(p))
  lx <- log(x)
  lm0 <- stats::lm(y ~ lx)
  b0 <- -as.numeric(stats::coef(lm0)[["lx"]])
  a0 <- exp(as.numeric(stats::coef(lm0)[["(Intercept)"]]))

  if (!is.finite(b0) || b0 <= 0) b0 <- 1
  if (!is.finite(a0) || a0 <= 0) a0 <- 1
  b0 <- max(1e-4, min(80, b0))
  a0 <- max(1e-12, a0)

  fml <- p ~ exp(-exp(loga) * x^(-exp(logb)))
  ctrl <- stats::nls.control(maxiter = 500, tol = 1e-7, minFactor = 1 / 2048, warnOnly = TRUE)

  mult <- c(0.5, 0.8, 1, 1.25, 1.6, 2, 3, 5)
  b_grid <- unique(pmax(1e-4, pmin(80, b0 * mult)))
  i_mid <- which.min(abs(p - 0.50))
  assert(
    length(i_mid) == 1L && p[i_mid] == 0.50,
    "Gompertz prefit requires p_anchor to include the 0.50 median anchor."
  )

  fit <- NULL
  for (b_try in b_grid) {
    a_try <- -log(p[i_mid]) * x[i_mid]^b_try
    if (!is.finite(a_try) || a_try <= 0) a_try <- a0

    fit <- tryCatch(
      stats::nls(
        fml,
        data = data.frame(p = p, x = x),
        start = list(loga = log(a_try), logb = log(b_try)),
        algorithm = "port",
        control = ctrl
      ),
      error = function(e) NULL
    )
    if (!is.null(fit) && isTRUE(fit$convInfo$isConv)) break
    fit <- NULL
  }

  if (is.null(fit)) stop("Gompertz pre-fit failed while building the K grid.", call. = FALSE)

  co <- stats::coef(fit)
  a <- exp(as.numeric(co[["loga"]]))
  b <- exp(as.numeric(co[["logb"]]))
  assert(is.finite(a) && a > 0 && is.finite(b) && b > 0, "Gompertz pre-fit returned invalid parameters.")
  list(a = a, b = b)
}

invert_gompertz <- function(a, b, p_grid, sim, k_search) {
  p <- pmin(pmax(as.numeric(p_grid), 1e-6), 1 - 1e-8)
  K <- sim$quasi_extinction_abundance + (a / -log(p))^(1 / b)
  K <- if (isTRUE(k_search$round_k)) round(K) else ceiling(K)
  K <- pmax(K, sim$quasi_extinction_abundance + 1)
  if (length(K) > 1L) {
    for (i in 2:length(K)) K[i] <- max(K[i], K[i - 1L] + 1)
  }
  assert(all(is.finite(K) & K <= k_search$max_k), "Inverted K grid exceeds k_search_max_k.")
  K
}

# ---- Partial-output inspection and recovery -----------------------------------

persist_resume_abort <- function(label, reason) {
  stop(
    label,
    " cannot be resumed from the existing CSV. ",
    reason,
    " Use mode: 'restart' to regenerate this taxon output.",
    call. = FALSE
  )
}

expected_persist_metadata <- function(sim, metadata) {
  posterior_rm_seed <- if (is.null(metadata$posterior_rm_seed)) {
    NA_integer_
  } else {
    as.integer(metadata$posterior_rm_seed)
  }
  posterior_sigma_seed <- if (is.null(metadata$posterior_sigma_seed)) {
    NA_integer_
  } else {
    as.integer(metadata$posterior_sigma_seed)
  }

  list(
    persistence_horizon_years = as.integer(sim$years),
    quasi_extinction_abundance = as.integer(sim$quasi_extinction_abundance),
    cap_factor = as.numeric(sim$cap_factor),
    r_buffer = as.numeric(sim$r_buffer),
    n_draws = as.integer(sim$n_draws),
    reps = as.integer(sim$reps),
    chunk_size = as.integer(sim$chunk_size),
    base_seed = as.integer(sim$base_seed),
    posterior_rm_seed = posterior_rm_seed,
    posterior_sigma_seed = posterior_sigma_seed,
    demographic_uncertainty = as.character(metadata$demographic_uncertainty),
    residual_growth_seed = suppressWarnings(as.integer(
      metadata$residual_growth_seed
    )),
    residual_environmental_variation_seed = suppressWarnings(as.integer(
      metadata$residual_environmental_variation_seed
    )),
    bird_model_signature = as.character(metadata$bird_model_signature),
    grid_signature = as.character(metadata$grid_signature),
    posterior_sampling = as.character(metadata$posterior_sampling),
    posterior_rm_md5 = as.character(metadata$posterior_rm_md5),
    posterior_sigma_md5 = as.character(metadata$posterior_sigma_md5),
    simulator_contract = as.character(metadata$simulator_contract)
  )
}

normalize_persist_trait_table <- function(trait_table, idx_col, value_col,
                                          extra_cols = character()) {
  assert(!is.null(trait_table), "trait_table is required.")
  trait_table <- as.data.frame(trait_table, stringsAsFactors = FALSE)
  required <- c(idx_col, value_col, extra_cols)
  missing_cols <- setdiff(required, names(trait_table))
  assert(
    !length(missing_cols),
    paste0("Trait table is missing column(s): ", paste(missing_cols, collapse = ", "), ".")
  )

  idx_num <- suppressWarnings(as.numeric(trait_table[[idx_col]]))
  value_num <- suppressWarnings(as.numeric(trait_table[[value_col]]))
  assert(
    all(
      is.finite(idx_num) &
        idx_num > 0 &
        idx_num == floor(idx_num) &
        is.finite(value_num) &
        value_num > 0
    ),
    paste0("Trait table must contain positive integer ", idx_col, " and positive finite ", value_col, ".")
  )

  trait_table[[idx_col]] <- as.integer(idx_num)
  trait_table[[value_col]] <- as.numeric(value_num)
  for (col in extra_cols) {
    trait_table[[col]] <- trimws(as.character(trait_table[[col]]))
    assert(
      all(!is.na(trait_table[[col]]) & nzchar(trait_table[[col]])),
      paste0("Trait table column ", col, " must be non-missing.")
    )
  }

  key_cols <- c(idx_col, extra_cols)
  trait_table$.trait_key <- do.call(
    paste,
    c(lapply(key_cols, function(col) as.character(trait_table[[col]])), sep = "\r")
  )
  assert(
    !anyDuplicated(trait_table$.trait_key),
    paste0("Trait table has duplicate key rows for: ", paste(key_cols, collapse = ", "), ".")
  )

  trait_table
}

inspect_existing_persist_points <- function(
  path,
  trait_table,
  sim,
  grid,
  metadata,
  idx_col,
  value_col,
  extra_cols = character(),
  label,
  require_complete = FALSE
) {
  trait_table <- normalize_persist_trait_table(
    trait_table = trait_table,
    idx_col = idx_col,
    value_col = value_col,
    extra_cols = extra_cols
  )
  key_cols <- c(idx_col, extra_cols)

  need_file(path, label)
  required_cols <- c(
    idx_col,
    value_col,
    extra_cols,
    "curve",
    "curve_probability",
    "p_target",
    "K",
    "p_eval",
    persistence_point_metadata_columns()
  )
  required_cols <- unique(required_cols)
  hdr <- tryCatch(
    data.table::fread(path, nrows = 0),
    error = function(e) {
      persist_resume_abort(label, paste0("The CSV could not be read: ", conditionMessage(e)))
    }
  )
  missing_cols <- setdiff(required_cols, names(hdr))
  if (length(missing_cols)) {
    persist_resume_abort(
      label,
      paste0("Missing current schema column(s): ", paste(missing_cols, collapse = ", "), ".")
    )
  }

  dt <- tryCatch(
    data.table::fread(path, select = required_cols),
    error = function(e) {
      persist_resume_abort(label, paste0("The CSV could not be read: ", conditionMessage(e)))
    }
  )

  if (!nrow(dt)) {
    if (isTRUE(require_complete)) {
      persist_resume_abort(label, "The CSV contains no rows.")
    }
    return(list(
      rows = dt,
      completed_keys = character(),
      missing_keys = character(),
      incomplete_keys = character(),
      expected_rows_per_key = length(grid$p_grid)
    ))
  }

  tryCatch(
    validate_persistence_point_metadata(
      dt[, c("curve", persistence_point_metadata_columns()), with = FALSE],
      label,
      expected_horizon = sim$years,
      expected_quasi_extinction = sim$quasi_extinction_abundance
    ),
    error = function(e) {
      persist_resume_abort(label, conditionMessage(e))
    }
  )

  expected_meta <- expected_persist_metadata(sim, metadata)
  check_integer_metadata <- function(col, expected) {
    values <- suppressWarnings(as.numeric(dt[[col]]))
    if (is.na(expected)) {
      bad <- !is.na(values)
    } else {
      bad <- is.na(values) | !is.finite(values) | values != expected | values != floor(values)
    }
    if (any(bad)) {
      persist_resume_abort(
        label,
        paste0("Column ", col, " does not match the current Stage 2 setting.")
      )
    }
  }
  check_numeric_metadata <- function(col, expected) {
    if (!all(persist_numeric_close(dt[[col]], expected))) {
      persist_resume_abort(
        label,
        paste0("Column ", col, " does not match the current Stage 2 setting.")
      )
    }
  }

  check_integer_metadata("persistence_horizon_years", expected_meta$persistence_horizon_years)
  check_integer_metadata("quasi_extinction_abundance", expected_meta$quasi_extinction_abundance)
  check_numeric_metadata("cap_factor", expected_meta$cap_factor)
  check_numeric_metadata("r_buffer", expected_meta$r_buffer)
  check_integer_metadata("n_draws", expected_meta$n_draws)
  check_integer_metadata("reps", expected_meta$reps)
  check_integer_metadata("chunk_size", expected_meta$chunk_size)
  check_integer_metadata("base_seed", expected_meta$base_seed)
  check_integer_metadata("posterior_rm_seed", expected_meta$posterior_rm_seed)
  check_integer_metadata("posterior_sigma_seed", expected_meta$posterior_sigma_seed)
  check_text_metadata <- function(col, expected) {
    values <- trimws(as.character(dt[[col]]))
    if (!all(!is.na(values) & nzchar(values) & values == expected)) {
      detail <- if (identical(col, "grid_signature")) {
        paste0(
          "Column grid_signature does not match the selected curves or current abundance-grid settings. ",
          "Set mode: 'restart' to regenerate this taxon output."
        )
      } else {
        paste0("Column ", col, " does not match the current Stage 2 setting.")
      }
      persist_resume_abort(label, detail)
    }
  }
  check_text_metadata(
    "demographic_uncertainty", expected_meta$demographic_uncertainty
  )
  check_text_metadata("bird_model_signature", expected_meta$bird_model_signature)
  for (col in c(
    "residual_growth_seed",
    "residual_environmental_variation_seed"
  )) {
    check_integer_metadata(col, expected_meta[[col]])
  }
  check_text_metadata("grid_signature", expected_meta$grid_signature)
  check_text_metadata("posterior_sampling", expected_meta$posterior_sampling)
  check_text_metadata("posterior_rm_md5", expected_meta$posterior_rm_md5)
  check_text_metadata("posterior_sigma_md5", expected_meta$posterior_sigma_md5)
  check_text_metadata("simulator_contract", expected_meta$simulator_contract)

  curve <- trimws(as.character(dt$curve))
  bad_curve <- is.na(curve) | !(curve %in% names(grid$q_levels))
  if (any(bad_curve)) {
    persist_resume_abort(
      label,
      paste0(
        "Existing rows contain curve labels outside the current contract: ",
        paste(names(grid$q_levels), collapse = ", "),
        "."
      )
    )
  }
  expected_curve_prob <- as.numeric(grid$q_levels[curve])
  if (!all(persist_numeric_close(dt$curve_probability, expected_curve_prob))) {
    persist_resume_abort(label, "curve_probability values do not match the current curve labels.")
  }
  p_target <- suppressWarnings(as.numeric(dt$p_target))
  valid_p_target <- persist_values_in_grid(p_target, grid$p_grid)
  if (any(!valid_p_target)) {
    persist_resume_abort(
      label,
      paste0(
        "p_target values do not match the current probability grid. ",
        format_invalid_persist_grid_values(p_target, grid$p_grid)
      )
    )
  }

  dt[[idx_col]] <- as.integer(suppressWarnings(as.numeric(dt[[idx_col]])))
  dt[[value_col]] <- suppressWarnings(as.numeric(dt[[value_col]]))
  for (col in extra_cols) dt[[col]] <- trimws(as.character(dt[[col]]))

  bad_trait_fields <- !is.finite(dt[[idx_col]]) |
    dt[[idx_col]] <= 0 |
    dt[[idx_col]] != floor(dt[[idx_col]]) |
    !is.finite(dt[[value_col]]) |
    dt[[value_col]] <= 0
  for (col in extra_cols) {
    bad_trait_fields <- bad_trait_fields | is.na(dt[[col]]) | !nzchar(dt[[col]])
  }
  if (any(bad_trait_fields)) {
    persist_resume_abort(label, "Existing rows contain invalid trait key or trait value fields.")
  }

  row_trait_key <- do.call(
    paste,
    c(lapply(key_cols, function(col) as.character(dt[[col]])), sep = "\r")
  )
  expected_idx <- match(row_trait_key, trait_table$.trait_key)
  if (any(is.na(expected_idx))) {
    persist_resume_abort(label, "Existing rows contain trait/diet keys outside the current trait grid.")
  }

  expected_trait_value <- trait_table[[value_col]][expected_idx]
  if (!all(persist_numeric_close(dt[[value_col]], expected_trait_value))) {
    persist_resume_abort(label, paste0(value_col, " values do not match the current trait grid."))
  }

  k_value <- suppressWarnings(as.numeric(dt$K))
  p_eval <- suppressWarnings(as.numeric(dt$p_eval))
  bad_values <- !is.finite(k_value) |
    k_value <= sim$quasi_extinction_abundance |
    !is.finite(p_eval) |
    p_eval < 0 |
    p_eval > 1
  if (any(bad_values)) {
    persist_resume_abort(label, "Existing rows contain invalid K or p_eval values.")
  }

  duplicated_rows <- duplicated(dt)
  if (any(duplicated_rows)) {
    persist_resume_abort(
      label,
      paste0("Existing rows contain ", sum(duplicated_rows), " exact duplicate row(s).")
    )
  }

  work <- data.table::copy(dt)
  work[, resume_trait_key := row_trait_key]
  work[, resume_curve := curve]
  work[, resume_key := paste(resume_trait_key, resume_curve, sep = "\r")]

  ordered <- data.table::copy(work)
  data.table::setorderv(ordered, c("resume_key", "p_target"))
  monotone <- ordered[, .(
    targets_increase = all(diff(p_target) > 0),
    K_increases = all(diff(as.numeric(K)) > 0),
    persistence_nondecreasing = all(diff(as.numeric(p_eval)) >= -sqrt(.Machine$double.eps))
  ), by = resume_key]
  if (any(!monotone$targets_increase | !monotone$K_increases | !monotone$persistence_nondecreasing)) {
    persist_resume_abort(label, "At least one trait/curve block is not monotone in p_target, K, and p_eval.")
  }

  expected_rows_per_key <- length(grid$p_grid)
  counts <- work[, .N, by = .(resume_trait_key, resume_curve, resume_key)]

  target_counts <- work[, .(rows = .N, n_targets = data.table::uniqueN(p_target)), by = resume_key]
  if (any(target_counts$n_targets != target_counts$rows)) {
    persist_resume_abort(label, "At least one trait/curve block repeats a configured p_target.")
  }

  overfull <- counts[N > expected_rows_per_key]
  if (nrow(overfull)) {
    overfull_labels <- gsub("\r", "/", overfull$resume_key, fixed = TRUE)
    persist_resume_abort(
      label,
      paste0(
        "Existing rows contain duplicate or over-complete trait/curve block(s): ",
        paste(head(overfull_labels, 8L), collapse = ", "),
        if (length(overfull_labels) > 8L) paste0(", ... and ", length(overfull_labels) - 8L, " more") else "",
        "."
      )
    )
  }

  complete <- counts[N == expected_rows_per_key]
  incomplete <- counts[N < expected_rows_per_key]

  expected_keys <- as.vector(outer(
    trait_table$.trait_key,
    names(grid$q_levels),
    paste,
    sep = "\r"
  ))
  completed_keys <- intersect(expected_keys, complete$resume_key)
  incomplete_keys <- intersect(expected_keys, incomplete$resume_key)
  missing_keys <- setdiff(expected_keys, completed_keys)

  if (isTRUE(require_complete) && length(missing_keys)) {
    persist_resume_abort(
      label,
      paste0(
        "The CSV is missing complete trait/curve block(s): ",
        paste(head(gsub("\r", "/", missing_keys), 8L), collapse = ", "),
        if (length(missing_keys) > 8L) paste0(", ... and ", length(missing_keys) - 8L, " more") else "",
        "."
      )
    )
  }

  complete_rows <- work[resume_key %in% completed_keys]
  complete_rows[, c("resume_trait_key", "resume_curve", "resume_key") := NULL]

  list(
    rows = complete_rows,
    completed_keys = completed_keys,
    missing_keys = missing_keys,
    incomplete_keys = incomplete_keys,
    expected_rows_per_key = expected_rows_per_key
  )
}

rewrite_persist_points_file <- function(rows, out_file) {
  if (!nrow(rows)) {
    if (file.exists(out_file)) file.remove(out_file)
    return(FALSE)
  }

  temp <- tempfile(".stage2_rewrite_", tmpdir = dirname(out_file), fileext = ".csv")
  on.exit(unlink(temp), add = TRUE)
  data.table::fwrite(rows, temp, append = FALSE, col.names = TRUE)
  check <- data.table::fread(temp)
  assert(nrow(check) == nrow(rows) && identical(names(check), names(rows)), "Temporary Stage 2 rewrite validation failed.")
  project_file_set_transaction(
    temp,
    out_file,
    overwrite = TRUE,
    label = "Stage 2 resumable output"
  )
  TRUE
}

validate_persist_points_complete <- function(
  path,
  trait_table,
  sim,
  grid,
  metadata,
  idx_col,
  value_col,
  extra_cols = character(),
  label
) {
  inspect_existing_persist_points(
    path = path,
    trait_table = trait_table,
    sim = sim,
    grid = grid,
    metadata = metadata,
    idx_col = idx_col,
    value_col = value_col,
    extra_cols = extra_cols,
    label = label,
    require_complete = TRUE
  )
  invisible(TRUE)
}

# Execute one taxon as ordered trait × curve blocks. The partial CSV stores
# completed blocks; the checkpoint stores only evaluated K values for the
# active block. Neither trajectories nor the compiled CRN context are
# serialized, because both are deterministically reconstructed on resume.
run_one_scenario <- function(
  trait_table,
  sampler,
  out_file,
  partial_file,
  checkpoint_file,
  taxon,
  crn_ctx,
  sim,
  grid,
  metadata = list(),
  idx_col,
  value_col,
  extra_cols = character(),
  existing_output,
  label = "persistence points",
  verbose = FALSE
) {
  existing_output <- validate_scalar_choice(
    existing_output,
    c("resume", "restart"),
    "existing_output"
  )
  ensure_writable_dir(dirname(partial_file), paste0(label, " output directory"))
  trait_table <- normalize_persist_trait_table(
    trait_table = trait_table,
    idx_col = idx_col,
    value_col = value_col,
    extra_cols = extra_cols
  )
  key_cols <- c(idx_col, extra_cols)

  completed_keys <- character()
  wrote_any <- FALSE
  checkpoint <- NULL

  # Recovery accepts only complete canonical blocks; a restart discards staged
  # transaction state without changing the requested trait or curve order.
  if (identical(existing_output, "restart")) {
    if (file.exists(partial_file)) unlink(partial_file)
    remove_stage2_checkpoint(checkpoint_file)
  } else {
    checkpoint <- read_stage2_checkpoint(checkpoint_file, sim, grid, metadata)
  }

  # A complete final file is already a successful result. Otherwise resume
  # from the partial file, or seed a new partial file with compatible complete
  # blocks found in an interrupted final artifact.
  source_file <- if (file.exists(partial_file)) {
    partial_file
  } else if (identical(existing_output, "resume") && file.exists(out_file)) {
    out_file
  } else {
    NA_character_
  }
  if (!is.na(source_file)) {
    resume <- inspect_existing_persist_points(
        path = source_file,
        trait_table = trait_table,
        sim = sim,
        grid = grid,
        metadata = metadata,
        idx_col = idx_col,
        value_col = value_col,
        extra_cols = extra_cols,
        label = label,
        require_complete = FALSE
      )

    completed_keys <- resume$completed_keys
    if (!length(resume$missing_keys) && !length(resume$incomplete_keys) &&
        identical(source_file, out_file)) {
      log_msg("RESUME | ", label, " | final output already complete", flush = TRUE)
      return(invisible(out_file))
    }
    wrote_any <- rewrite_persist_points_file(resume$rows, partial_file)
    log_msg(
      "RESUME | ", label,
      " | complete_blocks=", length(completed_keys),
      " | remaining_blocks=", length(resume$missing_keys),
      " | restored_evaluations=",
      if (is.null(checkpoint)) 0L else length(checkpoint$evaluations),
      flush = TRUE
    )
  } else {
    log_msg("RESUME | ", label, " | starting new staged output", flush = TRUE)
  }

  posterior_rm_seed <- if (is.null(metadata$posterior_rm_seed)) NA_integer_ else as.integer(metadata$posterior_rm_seed)
  posterior_sigma_seed <- if (is.null(metadata$posterior_sigma_seed)) NA_integer_ else as.integer(metadata$posterior_sigma_seed)

  write_points <- function(dt) {
    data.table::fwrite(dt, partial_file, append = wrote_any, col.names = !wrote_any)
    wrote_any <<- TRUE
  }

  # Nested trait then named-curve iteration defines the canonical append order.
  total_blocks <- nrow(trait_table) * length(grid$q_levels)
  initial_completed_blocks <- length(unique(completed_keys))
  run_started <- Sys.time()

  for (ti in seq_len(nrow(trait_table))) {
    trait_row <- trait_table[ti, , drop = FALSE]
    trait_value <- trait_row[[value_col]][[1L]]
    trait_key <- do.call(
      paste,
      c(lapply(key_cols, function(col) as.character(trait_row[[col]][[1L]])), sep = "\r")
    )
    trait_log <- paste(
      c(
        paste0(idx_col, "=", trait_row[[idx_col]][[1L]], "/", nrow(trait_table)),
        paste0(value_col, "=", formatC(trait_value, format = "f", digits = 4)),
        vapply(extra_cols, function(col) paste0(col, "=", trait_row[[col]][[1L]]), character(1))
      ),
      collapse = " "
    )
    needed_curves <- names(grid$q_levels)[
      !(paste(trait_key, names(grid$q_levels), sep = "\r") %in% completed_keys)
    ]

    if (!length(needed_curves)) {
      log_msg(
        "SKIP | ", trait_log,
        " | all curves already complete",
        flush = verbose
      )
      next
    }

    # The active trait owns one memoized K cache shared across its curves. Every
    # cache miss uses the same CRN context and is checkpointed immediately.
    pars <- sampler(trait_row)
    cache <- new.env(parent = emptyenv())
    evaluation_count <- 0L
    current_curve <- NULL
    last_checkpoint_at <- if (!is.null(checkpoint)) {
      as.character(checkpoint$saved_utc)
    } else {
      NA_character_
    }
    trait_started <- Sys.time()

    evalK <- function(K, who = "") {
      key <- sprintf("%.17g", K)
      v <- cache[[key]]
      if (!is.null(v)) return(v)

      evaluation_count <<- evaluation_count + 1L
      quantiles <- eval_persistence_quantiles_crn(
        pars$r,
        pars$sigma,
        K,
        crn_ctx,
        sim,
        q_levels = grid$q_levels
      )
      if (isTRUE(verbose)) {
        log_msg(
          "CRN", who,
          "| ", trait_log,
          "K=", formatC(K, format = "f", digits = 0),
          paste0(
            names(quantiles),
            "=",
            sprintf("%.4f", quantiles),
            collapse = " "
          ),
          flush = verbose
        )
      }

      cache[[key]] <- quantiles
      checkpoint_state <- write_stage2_checkpoint(
        checkpoint_file,
        taxon = taxon,
        trait_key = trait_key,
        curve = current_curve,
        cache = cache,
        sim = sim,
        grid = grid,
        metadata = metadata
      )
      last_checkpoint_at <<- checkpoint_state$saved_utc
      quantiles
    }

    for (qtag in needed_curves) {
      current_curve <- qtag
      restored <- restore_stage2_evaluation_cache(
        cache, checkpoint, trait_key, qtag
      )
      if (restored > 0L) {
        log_msg(
          "RESUME | ", trait_log, " | curve=", qtag,
          " | cached_evaluations=", restored,
          flush = TRUE
        )
      }
      curve_started <- Sys.time()
      evaluations_before_curve <- evaluation_count
      K_anchor <- numeric(length(grid$p_anchor))
      hi_prev <- as.numeric(grid$k_search$start_k)

      for (ai in seq_along(grid$p_anchor)) {
        pt <- grid$p_anchor[ai]
        who <- paste0("| stage=anchor q=", qtag, " p_target=", sprintf("%.3f", pt))
        get_q <- function(K) evalK(K, who = who)[[qtag]]

        low0 <- hi_prev / 2
        high0 <- hi_prev

        res <- find_K_for_target(pt, get_q, low0 = low0, high0 = high0,
                                 rel_tol = grid$k_search$rel_tol,
                                 max_k = grid$k_search$max_k,
                                 context = paste0(trait_log, " curve=", qtag))
        K_anchor[ai] <- res$K
        hi_prev <- res$hi_bracket
      }

      ab <- fit_anchor_gompertz_ab(K_anchor, grid$p_anchor, sim)
      K_guess <- invert_gompertz(ab$a, ab$b, grid$p_grid, sim, grid$k_search)
      if (any(!is.finite(K_guess) | K_guess <= sim$quasi_extinction_abundance)) {
        stop(
          "Gompertz prefit produced K values at or below the quasi-extinction threshold. ",
          "This indicates that the configured anchor prefit is outside the valid persistence domain.",
          call. = FALSE
        )
      }

      p_hat <- numeric(length(K_guess))
      for (ki in seq_along(K_guess)) {
        who <- paste0("| stage=grid q=", qtag, " p_target=", sprintf("%.3f", grid$p_grid[ki]))
        p_hat[ki] <- evalK(K_guess[ki], who = who)[[qtag]]
      }
      assert(all(diff(K_guess) > 0), paste0(trait_log, " curve=", qtag, " produced non-increasing K values."))
      assert(
        all(is.finite(p_hat) & p_hat >= 0 & p_hat <= 1) && all(diff(p_hat) >= -sqrt(.Machine$double.eps)),
        paste0(trait_log, " curve=", qtag, " produced invalid or non-monotone persistence values.")
      )

      dt <- data.table::data.table(
        curve = qtag,
        curve_probability = as.numeric(grid$q_levels[[qtag]]),
        p_target = as.numeric(grid$p_grid),
        K = as.numeric(K_guess),
        p_eval = as.numeric(p_hat),
        persistence_horizon_years = as.integer(sim$years),
        quasi_extinction_abundance = as.integer(sim$quasi_extinction_abundance),
        cap_factor = as.numeric(sim$cap_factor),
        r_buffer = as.numeric(sim$r_buffer),
        n_draws = as.integer(sim$n_draws),
        reps = as.integer(sim$reps),
        chunk_size = as.integer(sim$chunk_size),
        base_seed = as.integer(sim$base_seed),
        posterior_rm_seed = posterior_rm_seed,
        posterior_sigma_seed = posterior_sigma_seed,
        demographic_uncertainty = as.character(metadata$demographic_uncertainty),
        residual_growth_seed = suppressWarnings(as.integer(
          metadata$residual_growth_seed
        )),
        residual_environmental_variation_seed = suppressWarnings(as.integer(
          metadata$residual_environmental_variation_seed
        )),
        bird_model_signature = as.character(metadata$bird_model_signature),
        grid_signature = as.character(metadata$grid_signature),
        posterior_sampling = as.character(metadata$posterior_sampling),
        posterior_rm_md5 = as.character(metadata$posterior_rm_md5),
        posterior_sigma_md5 = as.character(metadata$posterior_sigma_md5),
        simulator_contract = as.character(metadata$simulator_contract)
      )
      dt[, (idx_col) := as.integer(trait_row[[idx_col]][[1L]])]
      dt[, (value_col) := as.numeric(trait_value)]
      for (col in extra_cols) {
        dt[, (col) := as.character(trait_row[[col]][[1L]])]
      }
      data.table::setcolorder(dt, c(
        idx_col, value_col, extra_cols, "curve", "curve_probability", "p_target", "K", "p_eval",
        "persistence_horizon_years", "quasi_extinction_abundance",
        "cap_factor", "r_buffer", "n_draws", "reps", "chunk_size",
        "base_seed", "posterior_rm_seed", "posterior_sigma_seed",
        "demographic_uncertainty", "residual_growth_seed",
        "residual_environmental_variation_seed", "bird_model_signature",
        "grid_signature", "posterior_sampling", "posterior_rm_md5",
        "posterior_sigma_md5", "simulator_contract"
      ))
      # Append one complete trait/curve block in canonical order; disk validation
      # below precedes removal of its recoverable evaluation checkpoint.
      write_points(dt)
      completed_key <- paste(trait_key, qtag, sep = "\r")

      # Validate the appended block from disk before considering it complete.
      # A crash during append is therefore recovered as an incomplete block,
      # while a successful append is confirmed before its checkpoint is removed.
      partial_state <- inspect_existing_persist_points(
        path = partial_file,
        trait_table = trait_table,
        sim = sim,
        grid = grid,
        metadata = metadata,
        idx_col = idx_col,
        value_col = value_col,
        extra_cols = extra_cols,
        label = paste0(label, " partial output"),
        require_complete = FALSE
      )
      assert(
        completed_key %in% partial_state$completed_keys,
        paste0(label, " appended block failed immediate validation.")
      )
      completed_keys <- partial_state$completed_keys
      remove_stage2_checkpoint(checkpoint_file)
      checkpoint <- NULL
      elapsed <- as.numeric(difftime(Sys.time(), run_started, units = "secs"))
      completed_n <- length(unique(completed_keys))
      completed_this_run <- completed_n - initial_completed_blocks
      eta <- if (completed_this_run > 0L) {
        elapsed / completed_this_run * (total_blocks - completed_n)
      } else {
        NA_real_
      }
      runtime_log_event(
        "stage2_curve_timing",
        taxon = taxon,
        trait_field = value_col,
        trait_value = sprintf("%.6f", trait_value),
        trait_index = as.integer(trait_row[[idx_col]][[1L]]),
        trait_position = as.integer(ti),
        trait_total = as.integer(nrow(trait_table)),
        trait_group = if (length(extra_cols)) {
          paste(vapply(extra_cols, function(column) {
            paste0(column, "=", trait_row[[column]][[1L]])
          }, character(1L)), collapse = ",")
        } else {
          "none"
        },
        curve = qtag,
        completed_blocks = as.integer(completed_n),
        total_blocks = as.integer(total_blocks),
        evaluations = as.integer(evaluation_count - evaluations_before_curve),
        elapsed_seconds = sprintf(
          "%.3f", as.numeric(difftime(Sys.time(), curve_started, units = "secs"))
        ),
        eta_seconds = if (is.finite(eta)) sprintf("%.0f", eta) else "NA",
        last_checkpoint_utc = last_checkpoint_at,
        k_min = sprintf("%.6f", min(K_guess)),
        k_max = sprintf("%.6f", max(K_guess)),
        persistence_min = sprintf("%.6f", min(p_hat)),
        persistence_max = sprintf("%.6f", max(p_hat))
      )

      if (isTRUE(verbose)) {
        x_anchor <- as.numeric(K_anchor) - sim$quasi_extinction_abundance
        p_fit <- exp(-ab$a * x_anchor^(-ab$b))
        delta <- p_fit - as.numeric(grid$p_anchor)
        log_msg(
          "GOMP | stage=pregrid q=", qtag,
          "| ", trait_log,
          "| deltas(p_fit - p_target)=",
          paste(sprintf("p=%.3f:%+.5f", grid$p_anchor, delta), collapse = " "),
          flush = verbose
        )
      }
    }

    runtime_log_event(
      "stage2_trait_timing",
      taxon = taxon,
      trait_field = value_col,
      trait_value = sprintf("%.6f", trait_value),
      trait_index = as.integer(trait_row[[idx_col]][[1L]]),
      trait_position = as.integer(ti),
      trait_total = as.integer(nrow(trait_table)),
      curves_completed = as.integer(length(needed_curves)),
      evaluations = as.integer(evaluation_count),
      elapsed_seconds = sprintf(
        "%.3f", as.numeric(difftime(Sys.time(), trait_started, units = "secs"))
      )
    )
  }

  # Promote the staged CSV only after the complete canonical block set validates.
  validate_persist_points_complete(
    partial_file,
    trait_table = trait_table,
    sim = sim,
    grid = grid,
    metadata = metadata,
    idx_col = idx_col,
    value_col = value_col,
    extra_cols = extra_cols,
    label = paste0(label, " staged output")
  )
  project_file_set_transaction(
    partial_file,
    out_file,
    overwrite = TRUE,
    label = paste0("Stage 2 ", taxon)
  )
  remove_stage2_checkpoint(checkpoint_file)
  invisible(out_file)
}
