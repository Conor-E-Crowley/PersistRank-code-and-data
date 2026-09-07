# Gompertz-parameter and LOESS presentation for Stage 3.
#
# Loaded by stage3_workflow.R after common Stage 3 figure helpers. It requires
# prepared fit and prediction tables and returns plot objects only. Sourcing is
# side-effect free and performs no artifact I/O.

gompertz_loess_taxon_available <- function(gomp_params,
                                           loess_pred,
                                           group,
                                           curve = NULL,
                                           bird_sigma_model_group = NULL) {
  need_cols(gomp_params, c("group", "curve", "fit_ok", "trait_value", "alpha", "beta"), "Gompertz parameter table")
  need_cols(loess_pred, c("group", "curve", "trait_value", "alpha_mid", "beta_mid"), "LOESS prediction grid")
  validate_gompertz_figure_parameters(gomp_params)
  validate_loess_figure_predictions(loess_pred)

  params <- gomp_params |>
    dplyr::filter(
      .data$group == .env$group,
      .data$fit_ok
    )
  pred <- loess_pred |>
    dplyr::filter(.data$group == .env$group)

  if (!is.null(curve)) {
    params <- params |> dplyr::filter(.data$curve == as.character(.env$curve))
    pred <- pred |> dplyr::filter(.data$curve == as.character(.env$curve))
  }
  if (identical(group, "Birds") && !is.null(bird_sigma_model_group)) {
    bird_sigma_model_group <- validate_scalar_string(
      bird_sigma_model_group,
      "bird_sigma_model_group"
    )
    need_cols(params, "bird_sigma_model_group", "Gompertz parameter table")
    need_cols(pred, "bird_sigma_model_group", "LOESS prediction grid")
    params <- params |> dplyr::filter(.data$bird_sigma_model_group == .env$bird_sigma_model_group)
    pred <- pred |> dplyr::filter(.data$bird_sigma_model_group == .env$bird_sigma_model_group)
  }

  isTRUE(nrow(params) > 0L && nrow(pred) > 0L)
}

gompertz_alpha_axis_labels <- function() {
  function(x) {
    labels <- rep(NA_character_, length(x))
    valid <- is.finite(x) & x > 0
    if (!any(valid)) return(scales::parse_format()(labels))

    exponent <- floor(log10(x[valid]))
    mantissa <- x[valid] / (10 ^ exponent)
    is_power_of_ten <- abs(mantissa - 1) < sqrt(.Machine$double.eps)
    mantissa_label <- scales::label_number(accuracy = 0.1, trim = TRUE)(mantissa)

    labels_valid <- paste0(mantissa_label, " %*% 10^", exponent)
    labels_valid[is_power_of_ten] <- paste0("10^", exponent[is_power_of_ten])
    labels[valid] <- labels_valid
    scales::parse_format()(labels)
  }
}

build_gompertz_loess_combined_figure <- function(gomp_params,
                                                 loess_pred,
                                                 curves = persistence_curves()) {
  curves <- as.character(curves)
  validate_gompertz_figure_parameters(gomp_params)
  validate_loess_figure_predictions(loess_pred)
  need_cols(
    gomp_params,
    c("group", "trait_value", "curve", "alpha", "beta", "fit_ok"),
    "Gompertz parameter table for combined LOESS figure"
  )
  need_cols(
    loess_pred,
    c("group", "trait_value", "curve", "alpha_mid", "alpha_lo", "alpha_hi", "beta_mid", "beta_lo", "beta_hi"),
    "LOESS prediction grid for combined LOESS figure"
  )

  has_mammals <- gompertz_loess_taxon_available(gomp_params, loess_pred, "Mammals")
  bird_groups <- unique(as.character(
    loess_pred$bird_sigma_model_group[
      loess_pred$group == "Birds" &
        !is.na(loess_pred$bird_sigma_model_group)
    ]
  ))
  has_birds <- stats::setNames(
    vapply(
      bird_groups,
      function(diet_group) {
        gompertz_loess_taxon_available(
          gomp_params,
          loess_pred,
          "Birds",
          bird_sigma_model_group = diet_group
        )
      },
      logical(1)
    ),
    bird_groups
  )
  assert(has_mammals || any(has_birds), "No mammal or bird LOESS results are available for the combined LOESS figure.")

  curve_styles <- methods_curve_styles(curves)
  legend_curves <- c("q025", "q16", "q50", "q84", "q975")
  legend_curves <- legend_curves[legend_curves %in% curves]
  # Compact quantile notation preserves the complete uncertainty meaning while
  # remaining readable when the six-panel figure is reduced to page width.
  curve_labs <- c(
    q025 = "q0.025",
    q16 = "q0.16",
    q50 = "q0.50 (median)",
    q84 = "q0.84",
    q975 = "q0.975"
  )
  pal <- methods_figure_palette()
  main_curve <- main_persistence_curve()
  log_axis_breaks <- scales::breaks_log(n = 4)
  linear_axis_breaks <- scales::breaks_pretty(n = 4)

  x_axis_label <- function(group) {
    if (identical(group, "Mammals")) "Body mass (g)" else "Generation length (years)"
  }
  x_axis_labels <- function(group) {
    if (identical(group, "Mammals")) fmt_mass_axis else fmt_generation_axis
  }
  plot_panel <- function(group, param, show_x, show_y, title = NULL,
                         show_legend = FALSE, bird_group = NULL) {
    if (identical(group, "Birds")) {
      bird_group <- validate_scalar_string(bird_group, "bird_group")
      assert(
        bird_group %in% bird_groups,
        paste0("Unknown bird model group: ", bird_group, ".")
      )
    }

    param <- match.arg(param, c("alpha", "beta"))
    mid_col <- paste0(param, "_mid")
    lo_col <- paste0(param, "_lo")
    hi_col <- paste0(param, "_hi")

    pts <- gomp_params |>
      dplyr::filter(
        .data$group == .env$group,
        .data$curve %in% .env$curves,
        .data$fit_ok
      ) |>
      dplyr::mutate(curve = factor(.data$curve, levels = curves))
    if (identical(group, "Birds")) {
      pts <- pts |> dplyr::filter(.data$bird_sigma_model_group == .env$bird_group)
    }

    pred <- loess_pred |>
      dplyr::filter(
        .data$group == .env$group,
        .data$curve %in% .env$curves
      ) |>
      dplyr::mutate(curve = factor(.data$curve, levels = curves))
    if (identical(group, "Birds")) {
      pred <- pred |> dplyr::filter(.data$bird_sigma_model_group == .env$bird_group)
    }

    pred_main <- pred |> dplyr::filter(.data$curve == main_curve)
    pred_sens <- pred |> dplyr::filter(.data$curve != main_curve)
    pts_main <- pts |> dplyr::filter(.data$curve == main_curve)
    pts_sens <- pts |> dplyr::filter(.data$curve != main_curve)

    p <- ggplot2::ggplot() +
      ggplot2::geom_ribbon(
        data = pred_main,
        ggplot2::aes(x = .data$trait_value, ymin = .data[[lo_col]], ymax = .data[[hi_col]]),
        fill = pal[["main_fill"]],
        alpha = 0.18,
        colour = NA
      ) +
      ggplot2::geom_point(
        data = pts_sens,
        ggplot2::aes(x = .data$trait_value, y = .data[[param]], color = .data$curve),
        size = 1.10,
        alpha = 0.30,
        show.legend = FALSE
      ) +
      ggplot2::geom_point(
        data = pts_main,
        ggplot2::aes(x = .data$trait_value, y = .data[[param]], color = .data$curve),
        size = 1.50,
        alpha = 0.62,
        show.legend = FALSE
      ) +
      ggplot2::geom_line(
        data = pred_sens,
        ggplot2::aes(x = .data$trait_value, y = .data[[mid_col]], color = .data$curve, linetype = .data$curve),
        linewidth = 0.86,
        alpha = 0.80,
        lineend = "round"
      ) +
      ggplot2::geom_line(
        data = pred_main,
        ggplot2::aes(x = .data$trait_value, y = .data[[mid_col]], color = .data$curve, linetype = .data$curve),
        linewidth = 1.18,
        alpha = 1,
        lineend = "round"
      ) +
      ggplot2::scale_x_log10(
        breaks = log_axis_breaks,
        labels = x_axis_labels(group),
        expand = ggplot2::expansion(mult = c(0.025, 0.060))
      ) +
      ggplot2::scale_color_manual(
        values = stats::setNames(curve_styles$color, curve_styles$curve),
        breaks = legend_curves,
        labels = curve_labs[legend_curves],
        name = NULL
      ) +
      ggplot2::scale_linetype_manual(
        values = stats::setNames(curve_styles$linetype, curve_styles$curve),
        breaks = legend_curves,
        labels = curve_labs[legend_curves],
        name = NULL
      ) +
      ggplot2::guides(
        color = ggplot2::guide_legend(
          nrow = 1,
          byrow = TRUE,
          keywidth = grid::unit(1.55, "lines"),
          keyheight = grid::unit(0.68, "lines"),
          override.aes = list(shape = NA, linewidth = 1.05, alpha = 1)
        )
      ) +
      ggplot2::labs(
        title = title,
        x = if (isTRUE(show_x)) x_axis_label(group) else NULL,
        y = if (isTRUE(show_y)) {
          if (identical(param, "alpha")) expression(alpha) else expression(beta)
        } else NULL
      ) +
      theme_methods_figure(base_size = 15.5, legend_position = if (isTRUE(show_legend)) "bottom" else "none") +
      ggplot2::theme(
        plot.margin = ggplot2::margin(6, 5, 4, 5),
        plot.title = ggplot2::element_text(
          size = 14.5,
          face = "bold",
          hjust = 0.5,
          lineheight = 0.98
        ),
        legend.text = ggplot2::element_text(size = 11.8),
        legend.spacing.x = grid::unit(0.20, "lines"),
        legend.box.margin = ggplot2::margin(0),
        axis.text.x = if (isTRUE(show_x)) ggplot2::element_text() else ggplot2::element_blank(),
        axis.ticks.x = if (isTRUE(show_x)) ggplot2::element_line(linewidth = 0.34, colour = pal[["axis"]]) else ggplot2::element_blank()
      )

    if (identical(param, "alpha")) {
      p + ggplot2::scale_y_log10(
        breaks = log_axis_breaks,
        labels = gompertz_alpha_axis_labels(),
        expand = ggplot2::expansion(mult = c(0.02, 0.06))
      )
    } else {
      p + ggplot2::scale_y_continuous(
        breaks = linear_axis_breaks,
        labels = scales::label_number(),
        expand = ggplot2::expansion(mult = c(0.02, 0.06))
      )
    }
  }

  if (has_mammals) {
    legend <- cowplot::get_legend(plot_panel("Mammals", "alpha", show_x = TRUE, show_y = TRUE, show_legend = TRUE))
  } else {
    donor_bird_group <- names(has_birds)[which(has_birds)[[1L]]]
    legend <- cowplot::get_legend(
      plot_panel("Birds", "alpha", show_x = TRUE, show_y = TRUE, show_legend = TRUE, bird_group = donor_bird_group)
    )
  }

  branches <- c(
    if (has_mammals) "Mammals",
    names(has_birds)[has_birds]
  )
  branch_title <- function(branch) {
    if (identical(branch, "Mammals")) return("Mammals")
    if (identical(branch, "Other")) return("Birds — other diets")
    paste0("Birds — ", branch)
  }
  make_branch_panel <- function(branch, parameter, column) {
    group <- if (identical(branch, "Mammals")) "Mammals" else "Birds"
    plot_panel(
      group = group,
      param = parameter,
      show_x = identical(parameter, "beta"),
      show_y = column == 1L,
      title = if (identical(parameter, "alpha")) branch_title(branch) else NULL,
      bird_group = if (identical(group, "Birds")) branch else NULL
    )
  }
  panels <- c(
    lapply(seq_along(branches), function(i) {
      make_branch_panel(branches[[i]], "alpha", i)
    }),
    lapply(seq_along(branches), function(i) {
      make_branch_panel(branches[[i]], "beta", i)
    })
  )

  grid <- panel_tag_grid(
    plotlist = panels,
    labels = paste0("(", letters[seq_along(panels)], ")"),
    ncol = length(branches),
    label_size = 14.5
  )

  cowplot::plot_grid(
    grid,
    legend,
    ncol = 1,
    rel_heights = c(1, 0.11)
  )
}
