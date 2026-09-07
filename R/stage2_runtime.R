# Deferred native runtime for Stage 2 persistence simulation.
#
# Loaded by stage2_workflow.R after pure configuration. Sourcing only defines
# the simulator/OpenMP contract and runtime calls: compilation, native symbol
# validation, CRN allocation, and self-tests occur only when explicitly invoked.

validate_persist_openmp_info <- function(info) {
  assert(
    is.list(info) && isTRUE(info$enabled),
    "Stage 2 requires OpenMP, but the compiled simulator reports that OpenMP is disabled."
  )
  assert(
    length(info$max_threads) == 1L && is.finite(info$max_threads) && info$max_threads >= 1L &&
      length(info$num_procs) == 1L && is.finite(info$num_procs) && info$num_procs >= 1L,
    "The compiled simulator returned invalid OpenMP runtime information."
  )
  if (as.integer(info$max_threads) == 1L) {
    warning(
      "OpenMP is enabled, but the runtime currently exposes only one thread. Stage 2 will run serially.",
      call. = FALSE
    )
  }
  invisible(info)
}

stage2_simulator_contract <- function() {
  "initial-threshold-v4-validated-quantiles"
}

compile_persist_cpp <- function(paths, rebuild_cpp = FALSE) {
  cache_dir <- paths$rcpp_cache
  if (isTRUE(rebuild_cpp)) {
    cache_dir <- tempfile("sourceCpp_", tmpdir = paths$rcpp_cache)
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }

  compiled <- tryCatch(
    {
      Rcpp::sourceCpp(
        file = paths$cpp_file,
        cacheDir = cache_dir,
        rebuild = isTRUE(rebuild_cpp),
        verbose = FALSE
      )
      TRUE
    },
    error = function(e) e
  )
  if (inherits(compiled, "error")) {
    stop(
      "Stage 2 requires an OpenMP-capable C++ compiler and runtime, but the simulator failed to compile. ",
      "Install or enable OpenMP for the server's R toolchain, then rebuild the simulator.\n",
      "Original compilation error: ", conditionMessage(compiled),
      call. = FALSE
    )
  }

  assert(exists("crn_context_create"), "C++ symbol missing: crn_context_create")
  assert(
    exists("simulate_persist_quantiles_cpp"),
    "C++ symbol missing: simulate_persist_quantiles_cpp"
  )
  assert(exists("persist_cpp_contract"), "C++ symbol missing: persist_cpp_contract")
  assert(exists("persist_openmp_info"), "C++ symbol missing: persist_openmp_info")
  assert(
    identical(
      as.character(persist_cpp_contract()),
      stage2_simulator_contract()
    ),
    "Compiled C++ simulator is stale. Rebuild src/simulate_persist_probs_cpp.cpp."
  )

  info <- persist_openmp_info()
  validate_persist_openmp_info(info)
  info
}

create_crn_context <- function(sim) {
  crn_context_create(
    seed = as.integer(sim$base_seed),
    n_draws = as.integer(sim$n_draws),
    reps = as.integer(sim$reps),
    years = as.integer(sim$years),
    chunk_size = as.integer(sim$chunk_size)
  )
}

eval_persistence_quantiles_crn <- function(r, sigma, K, crn_ctx, sim,
                                           q_levels = persistence_quantiles()) {
  assert(
    is.numeric(q_levels) &&
      length(q_levels) > 0L &&
      !is.null(names(q_levels)) &&
      all(!is.na(names(q_levels)) & nzchar(names(q_levels))) &&
      !anyDuplicated(names(q_levels)) &&
      all(is.finite(q_levels)) &&
      all(q_levels >= 0 & q_levels <= 1),
    "q_levels must be a non-empty, uniquely named numeric vector in [0, 1]."
  )
  assert(
    is.numeric(r) && length(r) == sim$n_draws && all(is.finite(r) & r > 0),
    "r must contain one positive finite value per posterior draw."
  )
  assert(
    is.numeric(sigma) && length(sigma) == sim$n_draws && all(is.finite(sigma) & sigma >= 0),
    "sigma must contain one non-negative finite value per posterior draw."
  )
  assert(length(K) == 1L && is.finite(K) && K > 0, "K must be one positive finite value.")

  out <- simulate_persist_quantiles_cpp(
    r = r,
    sigma = sigma,
    K = as.double(K),
    ext_thr = as.integer(sim$quasi_extinction_abundance),
    cap_factor = as.double(sim$cap_factor),
    crn_ctx = crn_ctx,
    quantile_probs = unname(q_levels)
  )

  assert(
    length(out) == length(q_levels) && all(is.finite(out)),
    "C++ persistence quantile output has the wrong length or contains non-finite values."
  )
  stats::setNames(as.numeric(out), names(q_levels))
}

run_crn_self_test <- function(sim, q_levels = persistence_quantiles()) {
  tiny_sim <- sim
  tiny_sim$n_draws <- 2L
  tiny_ctx <- crn_context_create(
    seed = as.integer(sim$base_seed),
    n_draws = 2L,
    reps = 3L,
    years = 2L,
    chunk_size = 2L
  )

  q1 <- eval_persistence_quantiles_crn(
    r = c(0.1, 0.2),
    sigma = c(0.2, 0.3),
    K = 1000,
    crn_ctx = tiny_ctx,
    sim = tiny_sim,
    q_levels = q_levels
  )
  q2 <- eval_persistence_quantiles_crn(
    r = c(0.1, 0.2),
    sigma = c(0.2, 0.3),
    K = 1000,
    crn_ctx = tiny_ctx,
    sim = tiny_sim,
    q_levels = q_levels
  )

  assert(isTRUE(all.equal(q1, q2, tolerance = 0)),
         "CRN self-test failed: repeated calls did not match exactly.")
  ordered_quantiles <- q1[order(q_levels)]
  assert(
    all(diff(ordered_quantiles) >= 0),
    "CRN self-test failed: persistence quantiles are not monotone."
  )

  q_thr <- eval_persistence_quantiles_crn(
    r = c(0.1, 0.2),
    sigma = c(0.2, 0.3),
    K = as.double(sim$quasi_extinction_abundance),
    crn_ctx = tiny_ctx,
    sim = tiny_sim,
    q_levels = q_levels
  )
  q_below <- eval_persistence_quantiles_crn(
    r = c(0.1, 0.2),
    sigma = c(0.2, 0.3),
    K = as.double(sim$quasi_extinction_abundance - 1L),
    crn_ctx = tiny_ctx,
    sim = tiny_sim,
    q_levels = q_levels
  )

  assert(all(q_thr == 0) && all(q_below == 0),
         "CRN self-test failed: K <= ext_thr must have zero persistence.")
}
