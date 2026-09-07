# Common formatting and Wolff/Gompertz figures for Stage 3.
#
# Loaded by stage3_workflow.R after scientific model definitions. Builders take
# prepared tables and return plot objects only. Sourcing has no side effects and
# performs no artifact reads or writes.

fmt_mass_label <- function(mass_g) {
  if (!is.finite(mass_g)) return("NA")
  if (mass_g < 1000) return(paste0(scales::label_number(accuracy = 1)(mass_g), " g"))
  kg <- mass_g / 1000
  acc <- if (kg >= 100) 1 else 0.1
  paste0(scales::label_number(accuracy = acc, big.mark = ",")(kg), " kg")
}

fmt_mass_axis <- function(x) {
  x <- as.numeric(x)
  out <- rep(NA_character_, length(x))
  ok <- is.finite(x)
  if (!any(ok)) return(out)
  x_ok <- x[ok]
  is_g <- x_ok < 1000
  out_ok <- character(length(x_ok))
  if (any(is_g)) out_ok[is_g] <- paste0(scales::label_number(accuracy = 1, big.mark = ",")(x_ok[is_g]), " g")
  if (any(!is_g)) {
    kg <- x_ok[!is_g] / 1000
    acc <- ifelse(kg >= 100, 1, 0.1)
    out_ok[!is_g] <- paste0(scales::label_number(accuracy = acc, big.mark = ",")(kg), " kg")
  }
  out[ok] <- out_ok
  out
}

fmt_generation_axis <- function(x) {
  x <- as.numeric(x)
  out <- rep(NA_character_, length(x))
  ok <- is.finite(x)
  if (!any(ok)) return(out)
  acc <- ifelse(x[ok] >= 10, 1, 0.1)
  out[ok] <- paste0(scales::label_number(accuracy = acc, big.mark = ",")(x[ok]), " y")
  out
}

p_wolff <- function(K, c, A) exp(-A * exp(-c * K))

r2_prob <- function(y, yhat) {
  rss <- sum((y - yhat)^2, na.rm = TRUE)
  tss <- sum((y - mean(y, na.rm = TRUE))^2, na.rm = TRUE)
  if (is.finite(tss) && tss > 0) 1 - rss / tss else NA_real_
}

validate_gompertz_figure_parameters <- function(gomp_params, label = "Gompertz parameter table") {
  need_cols(gomp_params, c("trait_value", "alpha", "beta", "fit_ok"), label)
  assert(all(gomp_params$fit_ok), paste0(label, " contains failed fits."))
  assert(
    all(is.finite(gomp_params$trait_value) & gomp_params$trait_value > 0 &
          is.finite(gomp_params$alpha) & gomp_params$alpha > 0 &
          is.finite(gomp_params$beta) & gomp_params$beta > 0),
    paste0(label, " contains invalid trait values or Gompertz coefficients.")
  )
  invisible(TRUE)
}

validate_loess_figure_predictions <- function(loess_pred, label = "LOESS prediction grid") {
  required <- c("trait_value", "alpha_mid", "alpha_lo", "alpha_hi", "beta_mid", "beta_lo", "beta_hi")
  need_cols(loess_pred, required, label)
  assert(
    all(vapply(loess_pred[required], function(x) all(is.finite(x) & x > 0), logical(1))),
    paste0(label, " contains nonpositive or nonfinite values.")
  )
  invisible(TRUE)
}

fmt_r2 <- function(x) {
  if (is.finite(x)) sprintf("%.3f", x) else "NA"
}

fit_wolff_c_start <- function(K, p, A) {
  p <- clamp01(p)
  y <- log(-log(p)) - log(A)
  fit <- stats::lm(y ~ 0 + K)
  c0 <- -as.numeric(stats::coef(fit)[["K"]])
  if (!is.finite(c0) || c0 <= 0) c0 <- 1e-8
  c0
}

fit_wolff_c_nls <- function(K, p, A = 10) {
  p <- clamp01(p)
  c0 <- fit_wolff_c_start(K, p, A = A)
  c_grid <- c0 * c(0.25, 0.5, 1, 2, 4)

  attempt <- function(c_start, iters = 250) {
    suppressWarnings(
      try(
        stats::nls(
          p ~ exp(-A * exp(-c * K)),
          start = list(c = c_start), algorithm = "port", lower = c(c = 1e-12),
          control = stats::nls.control(warnOnly = TRUE, maxiter = iters)
        ),
        silent = TRUE
      )
    )
  }

  for (cc in c_grid) {
    fit <- attempt(cc)
    if (!inherits(fit, "try-error") && !is.null(fit$convInfo) && isTRUE(fit$convInfo$isConv)) {
      return(as.numeric(stats::coef(fit)[["c"]]))
    }
  }

  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else NULL
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(123)
  for (j in 1:20) {
    cc <- c0 * exp(stats::runif(1, log(0.2), log(5)))
    fit <- attempt(cc, iters = 350)
    if (!inherits(fit, "try-error") && !is.null(fit$convInfo) && isTRUE(fit$convInfo$isConv)) {
      return(as.numeric(stats::coef(fit)[["c"]]))
    }
  }

  stop("Wolff NLS fit failed for this panel.")
}

representative_trait_rows <- function(gomp_params,
                                      group = "Mammals",
                                      curve = main_persistence_curve(),
                                      probs = c(0, 0.50, 1)) {
  need_cols(
    gomp_params,
    c("group", "trait_idx", "trait_value", "curve", "fit_ok"),
    "Gompertz parameter table"
  )
  curve <- validate_persistence_curve(curve)

  vals <- gomp_params |>
    dplyr::filter(
      .data$group == .env$group,
      .data$curve == .env$curve,
      .data$fit_ok
    ) |>
    dplyr::distinct(.data$trait_idx, .data$trait_value) |>
    dplyr::arrange(.data$trait_value)

  assert(
    nrow(vals) >= length(probs),
    paste0(
      "Need at least ", length(probs), " successful ", group,
      " Gompertz fits for representative mass panels."
    )
  )

  log_vals <- log10(vals$trait_value)
  targets <- stats::quantile(log_vals, probs = probs, names = FALSE, type = 7)
  idx <- vapply(
    targets,
    function(target) which.min(abs(log_vals - target)),
    integer(1)
  )
  selected <- vals[unique(idx), , drop = FALSE]

  if (nrow(selected) < length(probs)) {
    fill <- setdiff(seq_len(nrow(vals)), idx)
    selected <- dplyr::bind_rows(
      selected,
      vals[fill[seq_len(length(probs) - nrow(selected))], , drop = FALSE]
    )
  }

  selected |>
    dplyr::arrange(.data$trait_value)
}

left_limit_for_threshold_fraction <- function(threshold, x_max, fraction = 0.065) {
  assert(
    all(is.finite(c(threshold, x_max, fraction))) &&
      threshold > 0 &&
      x_max > threshold &&
      fraction > 0 &&
      fraction < 1,
    "left_limit_for_threshold_fraction() requires 0 < threshold < x_max and fraction in (0, 1)."
  )

  10^((log10(threshold) - fraction * log10(x_max)) / (1 - fraction))
}

wolff_display_points <- function(x, stride = 3L) {
  need_cols(x, c("K", "p_eval"), "Gompertz/Wolff display points")
  assert(
    length(stride) == 1L &&
      is.numeric(stride) &&
      is.finite(stride) &&
      stride >= 1L &&
      stride == floor(stride),
    "Gompertz/Wolff display-point stride must be one positive integer."
  )
  stride <- as.integer(stride)

  x <- x |>
    dplyr::arrange(.data$K)

  assert(nrow(x) > 0L, "Need at least one Gompertz/Wolff point to display.")
  keep <- unique(c(seq.int(1L, nrow(x), by = stride), nrow(x)))
  x[keep, , drop = FALSE]
}

build_gompertz_wolff_figure <- function(persist_points,
                                         gomp_params,
                                         k0,
                                         curve = main_persistence_curve(),
                                         wolff_A = 10) {
  k0 <- validate_scalar_number(
    k0, "k0", minimum = 0, minimum_open = TRUE
  )
  curve <- validate_persistence_curve(curve)
  validate_gompertz_figure_parameters(gomp_params)
  need_cols(
    persist_points,
    c("group", "predictor", "trait_idx", "trait_value", "curve", "K", "p_eval"),
    "persistence points for Gompertz/Wolff figure"
  )
  need_cols(
    gomp_params,
    c("group", "predictor", "trait_idx", "trait_value", "curve", "alpha", "beta", "fit_ok"),
    "Gompertz parameter table for Gompertz/Wolff figure"
  )

  reps <- representative_trait_rows(gomp_params, group = "Mammals", curve = curve)
  pal <- methods_figure_palette()

  params <- gomp_params |>
    dplyr::filter(
      .data$group == "Mammals",
      .data$curve == .env$curve,
      .data$trait_idx %in% reps$trait_idx,
      .data$fit_ok
    ) |>
    dplyr::select("group", "predictor", "trait_idx", "trait_value", "curve", "alpha", "beta")

  df <- persist_points |>
    dplyr::filter(
      .data$group == "Mammals",
      .data$curve == .env$curve,
      .data$trait_idx %in% reps$trait_idx
    ) |>
    dplyr::inner_join(
      params,
      by = c("group", "predictor", "trait_idx", "trait_value", "curve")
    )

  assert(nrow(df) > 0L, "No mammal q50 persistence points are available for the Gompertz/Wolff figure.")

  key_gomp <- paste0("Shifted Gompertz (", curve, ")")
  key_wolff <- "Wolff"
  key_pts <- "Simulation points"
  series_levels <- c(key_gomp, key_wolff, key_pts)
  series_cols <- stats::setNames(
    c(pal[["main"]], pal[["sensitivity"]], pal[["sensitivity_light"]]),
    series_levels
  )
  x_breaks_by_panel <- list(
    c(1500, 15000, 150000),
    c(1000, 3000, 9000),
    c(500, 1000, 2000, 4000)
  )

  make_panel <- function(trait_idx, x_breaks, panel_label,
                         show_y = TRUE, show_legend = FALSE) {
    dfi <- df |>
      dplyr::filter(.data$trait_idx == .env$trait_idx) |>
      dplyr::arrange(.data$K)

    assert(nrow(dfi) > 1L, "Representative Gompertz/Wolff panel has too few persistence points.")
    assert(
      dplyr::n_distinct(dfi$trait_idx) == 1L,
      "Each Gompertz/Wolff panel must contain exactly one representative mass."
    )
    mass_g <- dfi$trait_value[[1L]]
    alpha <- dfi$alpha[[1L]]
    beta <- dfi$beta[[1L]]
    c_hat <- fit_wolff_c_nls(dfi$K, dfi$p_eval, A = wolff_A)

    k_max <- max(dfi$K, na.rm = TRUE)
    k_left <- left_limit_for_threshold_fraction(k0, k_max, fraction = 0.065)
    k_grid_gomp <- exp(seq(log(k0 * 1.0005), log(k_max), length.out = 650))
    k_grid_wolff <- exp(seq(log(k_left), log(k_max), length.out = 720))

    line_gomp <- tibble::tibble(
      K = k_grid_gomp,
      p = predict_p(alpha, beta, k_grid_gomp, K0 = k0),
      Series = factor(key_gomp, levels = series_levels)
    )
    line_wolff <- tibble::tibble(
      K = k_grid_wolff,
      p = p_wolff(k_grid_wolff, c_hat, A = wolff_A),
      Series = factor(key_wolff, levels = series_levels)
    )
    r2_gomp <- r2_prob(dfi$p_eval, predict_p(alpha, beta, dfi$K, K0 = k0))
    r2_wolff <- r2_prob(dfi$p_eval, p_wolff(dfi$K, c_hat, A = wolff_A))
    legend_x <- function(p) 10^(log10(k_left) + p * (log10(k_max) - log10(k_left)))
    r2_header <- tibble::tibble(
      text_x = legend_x(0.245),
      y = 0.960,
      label = "R^2"
    )
    r2_key <- tibble::tibble(
      x = legend_x(0.13),
      xend = legend_x(0.215),
      text_x = legend_x(0.245),
      y = c(0.900, 0.840),
      label = c(fmt_r2(r2_gomp), fmt_r2(r2_wolff))
    )
    pts_here <- dfi |>
      wolff_display_points() |>
      dplyr::transmute(K = .data$K, p = .data$p_eval, Series = factor(key_pts, levels = series_levels))

    ggplot2::ggplot() +
      ggplot2::geom_vline(
        xintercept = k0,
        linewidth = 0.35,
        linetype = "dotted",
        colour = pal[["threshold"]]
      ) +
      ggplot2::geom_line(
        data = line_gomp,
        ggplot2::aes(x = .data$K, y = .data$p, color = .data$Series),
        linewidth = 0.95,
        lineend = "round",
        na.rm = TRUE
      ) +
      ggplot2::geom_line(
        data = line_wolff,
        ggplot2::aes(x = .data$K, y = .data$p, color = .data$Series, linetype = .data$Series),
        linewidth = 0.68,
        lineend = "round",
        na.rm = TRUE
      ) +
      ggplot2::geom_point(
        data = pts_here,
        ggplot2::aes(x = .data$K, y = .data$p, color = .data$Series),
        shape = 21,
        fill = scales::alpha("white", 0.78),
        stroke = 0.40,
        alpha = 0.82,
        size = 1.48
      ) +
      ggplot2::geom_segment(
        data = r2_key[1L, , drop = FALSE],
        ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$y),
        colour = series_cols[[key_gomp]],
        linewidth = 0.95,
        lineend = "round"
      ) +
      ggplot2::geom_segment(
        data = r2_key[2L, , drop = FALSE],
        ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$y),
        colour = series_cols[[key_wolff]],
        linewidth = 0.68,
        linetype = "22",
        lineend = "round"
      ) +
      ggplot2::geom_text(
        data = r2_header,
        ggplot2::aes(x = .data$text_x, y = .data$y, label = .data$label),
        hjust = 0,
        vjust = 0.5,
        size = 2.95,
        colour = pal[["text"]],
        parse = TRUE
      ) +
      ggplot2::geom_text(
        data = r2_key,
        ggplot2::aes(x = .data$text_x, y = .data$y, label = .data$label),
        hjust = 0,
        vjust = 0.5,
        size = 2.95,
        colour = pal[["text"]]
      ) +
      ggplot2::scale_linetype_manual(
        values = stats::setNames(c("solid", "22", "blank"), series_levels),
        guide = "none"
      ) +
      ggplot2::scale_x_log10(
        limits = c(k_left, k_max),
        breaks = x_breaks,
        labels = methods_compact_number_labels,
        expand = ggplot2::expansion(mult = c(0, 0))
      ) +
      ggplot2::scale_y_continuous(
        limits = c(0, 1),
        breaks = methods_probability_breaks(),
        labels = methods_probability_labels(),
        expand = ggplot2::expansion(mult = c(0, 0.01))
      ) +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::scale_color_manual(
        values = series_cols,
        breaks = series_levels,
        name = NULL
      ) +
      ggplot2::guides(
        color = ggplot2::guide_legend(
          nrow = 1,
          byrow = TRUE,
          keywidth = grid::unit(1.15, "lines"),
          keyheight = grid::unit(0.58, "lines"),
          override.aes = list(
            linetype = c("solid", "22", "blank"),
            shape = c(NA, NA, 21),
            linewidth = c(0.95, 0.68, 0),
            alpha = c(1, 1, 1),
            size = c(NA, NA, 1.48),
            fill = c(NA, NA, "white")
          )
        )
      ) +
      ggplot2::labs(
        title = fmt_mass_label(mass_g),
        tag = panel_label,
        x = NULL,
        y = if (isTRUE(show_y)) "Persistence probability" else NULL
      ) +
      theme_methods_figure(base_size = 14, legend_position = if (isTRUE(show_legend)) "bottom" else "none") +
      theme_methods_probability_guides() +
      ggplot2::theme(
        plot.margin = ggplot2::margin(7, 4, 2, 4),
        plot.title.position = "panel",
        plot.tag.position = c(0, 1.03),
        plot.tag = ggplot2::element_text(
          size = 10.4,
          face = "bold",
          hjust = 0,
          vjust = 1
        ),
        plot.title = ggplot2::element_text(
          size = 11.4,
          face = "bold",
          hjust = 0.5,
          margin = ggplot2::margin(b = 2)
        ),
        legend.box.margin = ggplot2::margin(0),
        legend.spacing.x = grid::unit(0.28, "lines"),
        legend.text = ggplot2::element_text(
          size = 9.4,
          colour = pal[["text"]],
          margin = ggplot2::margin(r = 5)
        )
      )
  }

  p_leg <- make_panel(
    reps$trait_idx[[1L]],
    x_breaks = x_breaks_by_panel[[1L]],
    panel_label = "(a)",
    show_y = TRUE,
    show_legend = TRUE
  )
  legend <- cowplot::get_legend(p_leg)

  panels <- lapply(seq_len(nrow(reps)), function(i) {
    make_panel(
      reps$trait_idx[[i]],
      x_breaks = x_breaks_by_panel[[i]],
      panel_label = paste0("(", letters[[i]], ")"),
      show_y = i == 1L,
      show_legend = FALSE
    )
  })

  row_panels <- cowplot::plot_grid(
    plotlist = panels,
    ncol = length(panels),
    align = "hv",
    axis = "tblr"
  )

  shared_x <- cowplot::ggdraw() +
    cowplot::draw_label(
      expression("Abundance, " * italic(K)),
      x = 0.5,
      y = 0.5,
      size = 11.7,
      colour = pal[["text"]]
    )

  cowplot::plot_grid(
    row_panels,
    shared_x,
    legend,
    ncol = 1,
    rel_heights = c(1, 0.08, 0.10)
  )
}
