# Allometry presentation helpers for 1_demographic_models.Rmd.
#
# Loaded by stage1_workflow.R after Stage 1 fitting definitions. Builders accept
# validated calibration data and posterior draws and return plot objects.
# Sourcing has no side effects and performs no artifact I/O.

band_from_posterior <- function(x_grid, draws, coef_col, extra_terms = NULL) {
  draw_cols <- c("alpha", coef_col, names(extra_terms))
  need_cols(draws, draw_cols, "posterior draws used for a Stage 1 figure")
  log_x <- log10(x_grid)
  pred_log10 <- draws$alpha + tcrossprod(draws[[coef_col]], log_x)
  if (length(extra_terms)) {
    for (col in names(extra_terms)) {
      pred_log10 <- pred_log10 + draws[[col]] * as.numeric(extra_terms[[col]])
    }
  }
  pred <- 10^pred_log10

  qs <- apply(pred, 2, stats::quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
  tibble::tibble(x = x_grid, lo = qs[1, ], mid = qs[2, ], hi = qs[3, ])
}

log_breaks_in_range <- function(x_limits, candidates) {
  out <- candidates[candidates >= x_limits[1] & candidates <= x_limits[2]]
  if (length(out) >= 2L) return(out)
  scales::breaks_log(n = 4)(x_limits)
}

mass_axis_labels_g <- function(x_g) {
  x_g <- as.numeric(x_g)
  out <- scales::label_number(
    accuracy = 1,
    big.mark = ",",
    drop0trailing = TRUE
  )(x_g)
  use_k <- is.finite(x_g) & x_g >= 1000
  out[use_k] <- paste0(
    scales::label_number(
      scale = 1 / 1000,
      accuracy = 0.1,
      drop0trailing = TRUE
    )(x_g[use_k]),
    "k"
  )
  out
}

validate_allometry_figure_data <- function(df, label) {
  need_cols(df, c("x", "y"), label)
  assert(
    is.numeric(df$x) && is.numeric(df$y) &&
      all(is.finite(df$x) & df$x > 0) &&
      all(is.finite(df$y) & df$y > 0),
    paste0(label, " requires positive finite x and y values.")
  )
  invisible(TRUE)
}

bird_sigma_posterior_band <- function(x_grid, draws, model_group, bird_model) {
  extra_terms <- stats::setNames(
    as.numeric(bird_model$separate_intercepts == model_group),
    bird_model$coefficient_names
  )
  band_from_posterior(
    x_grid,
    draws,
    coef_col = "beta_logGenLength",
    extra_terms = extra_terms
  ) |>
    dplyr::mutate(bird_sigma_model_group = model_group)
}

demographic_calibration_style <- function() {
  list(
    base_size = 12,
    aspect_ratio = 0.62,
    panel_margin = ggplot2::margin(4, 5, 2, 5),
    axis_title_size = 10.8,
    axis_text_size = 9.6,
    annotation_size = 3.0,
    panel_label_size = 10.2,
    column_header_size = 10.8,
    header_fraction = 0.065
  )
}

demographic_fit_annotation <- function(n, r2) {
  sprintf("n = %d\nR\u00b2 = %.2f", n, r2)
}

theme_demographic_calibration_panel <- function(legend_position = "none") {
  style <- demographic_calibration_style()
  theme_methods_figure(
    base_size = style$base_size,
    legend_position = legend_position
  ) +
    ggplot2::theme(
      aspect.ratio = style$aspect_ratio,
      plot.margin = style$panel_margin,
      axis.title = ggplot2::element_text(
        size = style$axis_title_size,
        face = "plain"
      ),
      axis.title.x = ggplot2::element_text(
        margin = ggplot2::margin(t = 3)
      ),
      axis.title.y = ggplot2::element_text(
        margin = ggplot2::margin(r = 3)
      ),
      axis.text = ggplot2::element_text(size = style$axis_text_size),
      axis.ticks.length = grid::unit(1.8, "pt")
    )
}

make_allometry_panel <- function(df, y_lab, x_lab, x_limits, x_breaks,
                                 y_limits = NULL, y_breaks = NULL,
                                 draws, r2, coef_col = "beta_logM",
                                 x_labels = scales::label_number(big.mark = ",")) {
  validate_allometry_figure_data(df, "Stage 1 allometry panel data")

  x_grid <- log_space(x_limits[1], x_limits[2], n = 300)
  band <- band_from_posterior(x_grid, draws, coef_col)

  style <- demographic_calibration_style()
  ann <- demographic_fit_annotation(nrow(df), r2)
  pal <- methods_figure_palette()

  ggplot2::ggplot(df, ggplot2::aes(x, y)) +
    ggplot2::geom_ribbon(
      data = band,
      ggplot2::aes(x = x, ymin = lo, ymax = hi),
      inherit.aes = FALSE,
      fill = pal[["main_fill"]],
      alpha = 0.24
    ) +
    ggplot2::geom_line(
      data = band,
      ggplot2::aes(x = x, y = mid),
      inherit.aes = FALSE,
      linewidth = 1.02,
      color = pal[["ink"]],
      lineend = "round"
    ) +
    ggplot2::geom_point(
      shape = 21,
      size = 1.28,
      stroke = 0.38,
      alpha = 0.64,
      colour = pal[["ink"]],
      fill = scales::alpha("white", 0.64)
    ) +
    ggplot2::annotate(
      "label",
      x = Inf, y = Inf,
      label = ann,
      hjust = 1.03,
      vjust = 1.05,
      linewidth = 0,
      label.padding = grid::unit(0.08, "lines"),
      label.r = grid::unit(0, "pt"),
      fill = scales::alpha("white", 0.88),
      color = pal[["text"]],
      size = style$annotation_size,
      lineheight = 0.94
    ) +
    ggplot2::scale_x_log10(
      limits = x_limits,
      breaks = x_breaks,
      labels = x_labels,
      expand = ggplot2::expansion(mult = c(0.02, 0.04))
    ) +
    ggplot2::scale_y_log10(
      limits = y_limits,
      breaks = y_breaks,
      labels = scales::label_number(big.mark = ","),
      expand = ggplot2::expansion(mult = c(0.04, 0.08))
    ) +
    ggplot2::labs(x = x_lab, y = y_lab) +
    theme_demographic_calibration_panel()
}

make_bird_sigma_allometry_panel <- function(df, y_lab, x_lab, x_limits, x_breaks,
                                            y_limits = NULL, y_breaks = NULL,
                                            draws, r2, bird_model) {
  validate_allometry_figure_data(df, "Stage 1 bird sigma panel data")
  validate_bird_diet5(df$diet5_group, "bird sigma allometry Diet-5Cat")

  model_groups <- bird_model$model_groups
  df$bird_sigma_model_group <- map_bird_diet_to_model_group(
    df$diet5_group, bird_model
  )
  df$bird_sigma_model_group <- factor(df$bird_sigma_model_group, levels = model_groups)
  x_grid <- log_space(x_limits[1], x_limits[2], n = 300)
  band <- dplyr::bind_rows(lapply(
    model_groups,
    bird_sigma_posterior_band,
    x_grid = x_grid,
    draws = draws,
    bird_model = bird_model
  ))
  band$bird_sigma_model_group <- factor(
    band$bird_sigma_model_group, levels = model_groups
  )

  style <- demographic_calibration_style()
  ann <- demographic_fit_annotation(nrow(df), r2)
  pal <- methods_figure_palette()

  group_cols <- if (identical(model_groups, c("Other", "VertFishScav"))) {
    c(Other = pal[["sensitivity"]], VertFishScav = pal[["main"]])
  } else {
    stats::setNames(
      c(pal[["sensitivity"]], grDevices::hcl.colors(
        max(1L, length(model_groups) - 1L), "Blues 3", rev = TRUE
      )),
      model_groups
    )
  }
  group_shapes <- stats::setNames(
    rep(c(21, 24, 22, 23, 25), length.out = length(model_groups)),
    model_groups
  )
  group_linetypes <- stats::setNames(
    rep(c("solid", "31", "42", "13", "22"), length.out = length(model_groups)),
    model_groups
  )

  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = x,
      y = y,
      color = bird_sigma_model_group,
      shape = bird_sigma_model_group
    )
  ) +
    ggplot2::geom_ribbon(
      data = band,
      ggplot2::aes(
        x = x,
        ymin = lo,
        ymax = hi,
        fill = bird_sigma_model_group
      ),
      inherit.aes = FALSE,
      alpha = 0.13,
      colour = NA
    ) +
    ggplot2::geom_line(
      data = band,
      ggplot2::aes(
        x = x,
        y = mid,
        color = bird_sigma_model_group,
        linetype = bird_sigma_model_group
      ),
      inherit.aes = FALSE,
      linewidth = 0.92,
      lineend = "round"
    ) +
    ggplot2::geom_point(
      size = 1.22,
      stroke = 0.36,
      alpha = 0.62,
      fill = scales::alpha("white", 0.35)
    ) +
    ggplot2::annotate(
      "label",
      x = Inf, y = Inf,
      label = ann,
      hjust = 1.03,
      vjust = 1.05,
      linewidth = 0,
      label.padding = grid::unit(0.08, "lines"),
      label.r = grid::unit(0, "pt"),
      fill = scales::alpha("white", 0.88),
      color = pal[["text"]],
      size = style$annotation_size,
      lineheight = 0.94
    ) +
    ggplot2::scale_color_manual(values = group_cols, breaks = model_groups, name = NULL) +
    ggplot2::scale_fill_manual(values = group_cols, breaks = model_groups, name = NULL) +
    ggplot2::scale_shape_manual(values = group_shapes, breaks = model_groups, name = NULL) +
    ggplot2::scale_linetype_manual(
      values = group_linetypes,
      breaks = model_groups,
      name = NULL
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(
          shape = unname(group_shapes),
          linetype = unname(group_linetypes),
          linewidth = 0.82,
          alpha = 1
        )
      ),
      fill = "none",
      shape = "none",
      linetype = "none"
    ) +
    ggplot2::scale_x_log10(
      limits = x_limits,
      breaks = x_breaks,
      labels = scales::label_number(big.mark = ","),
      expand = ggplot2::expansion(mult = c(0.02, 0.04))
    ) +
    ggplot2::scale_y_log10(
      limits = y_limits,
      breaks = y_breaks,
      labels = scales::label_number(big.mark = ","),
      expand = ggplot2::expansion(mult = c(0.04, 0.08))
    ) +
    ggplot2::labs(x = x_lab, y = y_lab) +
    theme_demographic_calibration_panel(legend_position = "inside") +
    ggplot2::theme(
      legend.position.inside = c(0.5, 0.022),
      legend.justification = c(0.5, 0),
      legend.direction = "horizontal",
      legend.background = ggplot2::element_rect(
        fill = scales::alpha("white", 0.84),
        colour = NA
      ),
      legend.text = ggplot2::element_text(size = 7.2, colour = "#333333"),
      legend.spacing.x = grid::unit(0.20, "lines"),
      legend.spacing.y = grid::unit(0.02, "lines"),
      legend.margin = ggplot2::margin(0.5, 1.5, 0.5, 1.5),
      legend.key.height = grid::unit(0.34, "lines"),
      legend.key.width = grid::unit(0.68, "lines")
    )
}

build_demographic_calibration_figure <- function(
    cal, trait_grids, draws_by_model, bird_model,
    stage1_r2 = descriptive_stage1_r2(cal, bird_model)) {
  assert(
    is.numeric(stage1_r2) &&
      all(c(
        "mammal_growth", "mammal_environmental_variation",
        "bird_growth", "bird_environmental_variation"
      ) %in% names(stage1_r2)) &&
      all(is.finite(stage1_r2)),
    "Stage 1 figure R-squared values must be a complete finite named vector."
  )
  dat_mammal_rm <- dplyr::transmute(cal$mammal_rmax, x = Mass_g, y = rm)
  dat_mammal_sigma <- dplyr::transmute(cal$mammal_sigma, x = Mass_g, y = sigma)
  dat_bird_rm <- dplyr::transmute(cal$bird_rmax, x = GenLength, y = rm)
  dat_bird_sigma <- dplyr::transmute(
    cal$bird_sigma_data,
    x = GenLength,
    y = sigma,
    diet5_group
  )

  mammal_xlim <- range(trait_grids$mammal$Mass_g, finite = TRUE)
  bird_xlim <- range(trait_grids$bird$GenLength, finite = TRUE)
  x_breaks_mammal <- log_breaks_in_range(mammal_xlim, c(10, 1e3, 1e5))
  x_breaks_bird <- log_breaks_in_range(bird_xlim, c(1, 2, 5, 10, 20, 50))
  sigma_limits <- c(0.05, 4)
  sigma_breaks <- c(0.1, 0.3, 1, 3)

  pA <- make_allometry_panel(
    dat_mammal_rm,
    expression(r[m]),
    "Body mass (g)",
    x_limits = mammal_xlim,
    x_breaks = x_breaks_mammal,
    y_breaks = c(0.1, 1, 10),
    draws = draws_by_model$mammal_growth,
    r2 = stage1_r2[["mammal_growth"]],
    coef_col = "beta_logM",
    x_labels = mass_axis_labels_g
  )
  pB <- make_allometry_panel(
    dat_bird_rm,
    expression(r[m]~~"(from"~~lambda[max]~~")"),
    "Generation length (yr)",
    x_limits = bird_xlim,
    x_breaks = x_breaks_bird,
    y_breaks = c(0.1, 0.3, 1),
    draws = draws_by_model$bird_growth,
    r2 = stage1_r2[["bird_growth"]],
    coef_col = "beta_logGenLength"
  )
  pC <- make_allometry_panel(
    dat_mammal_sigma,
    expression(sigma[r]),
    "Body mass (g)",
    x_limits = mammal_xlim,
    x_breaks = x_breaks_mammal,
    y_limits = sigma_limits,
    y_breaks = sigma_breaks,
    draws = draws_by_model$mammal_environmental_variation,
    r2 = stage1_r2[["mammal_environmental_variation"]],
    coef_col = "beta_logM",
    x_labels = mass_axis_labels_g
  )
  pD <- make_bird_sigma_allometry_panel(
    dat_bird_sigma,
    expression(sigma[r]),
    "Generation length (yr)",
    x_limits = bird_xlim,
    x_breaks = x_breaks_bird,
    y_limits = sigma_limits,
    y_breaks = sigma_breaks,
    draws = draws_by_model$bird_environmental_variation,
    r2 = stage1_r2[["bird_environmental_variation"]],
    bird_model = bird_model
  )

  pA <- pA + ggplot2::theme(axis.title.x = ggplot2::element_blank())
  pB <- pB + ggplot2::theme(
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank()
  )
  pD <- pD + ggplot2::theme(axis.title.y = ggplot2::element_blank())

  style <- demographic_calibration_style()
  headers <- cowplot::ggdraw() +
    cowplot::draw_label(
      "(a)",
      x = 0.008, y = 0.44,
      hjust = 0, vjust = 0.5,
      fontface = "bold",
      size = style$panel_label_size,
      colour = methods_figure_palette()[["text"]]
    ) +
    cowplot::draw_label(
      "Mammals",
      x = 0.255, y = 0.44,
      hjust = 0.5, vjust = 0.5,
      fontface = "bold",
      size = style$column_header_size,
      colour = methods_figure_palette()[["text"]]
    ) +
    cowplot::draw_label(
      "(b)",
      x = 0.508, y = 0.44,
      hjust = 0, vjust = 0.5,
      fontface = "bold",
      size = style$panel_label_size,
      colour = methods_figure_palette()[["text"]]
    ) +
    cowplot::draw_label(
      "Birds",
      x = 0.755, y = 0.44,
      hjust = 0.5, vjust = 0.5,
      fontface = "bold",
      size = style$column_header_size,
      colour = methods_figure_palette()[["text"]]
    )
  grid <- cowplot::plot_grid(
    pA, pB, pC, pD,
    ncol = 2,
    labels = c("", "", "(c)", "(d)"),
    label_fontface = "bold",
    label_size = style$panel_label_size,
    label_colour = methods_figure_palette()[["text"]],
    label_x = 0.012,
    label_y = 0.997,
    hjust = 0,
    vjust = 1,
    align = "hv",
    axis = "tblr"
  )
  out <- cowplot::plot_grid(
    headers,
    grid,
    ncol = 1,
    rel_heights = c(style$header_fraction, 1)
  )
  attr(out, "stage1_calibration_panel_labels") <- c(
    "(a) Mammals", "(b) Birds", "(c)", "(d)"
  )
  out
}

