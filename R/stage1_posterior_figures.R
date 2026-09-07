# Posterior-coefficient presentation helpers for Stage 1.
#
# Loaded by stage1_workflow.R after the allometry figure layer. It requires the
# shared bird-model and figure contracts. Builders return plot objects only;
# sourcing and figure construction do not write artifacts.

tidy_demographic_coefficient_draws <- function(draws_by_model) {
  specs <- list(
    list(
      model = "mammal_growth",
      relationship = "mammal_rm",
      coef_cols = c(intercept = "alpha", slope = "beta_logM")
    ),
    list(
      model = "mammal_environmental_variation",
      relationship = "mammal_sigma",
      coef_cols = c(intercept = "alpha", slope = "beta_logM")
    ),
    list(
      model = "bird_growth",
      relationship = "bird_rm",
      coef_cols = c(intercept = "alpha", slope = "beta_logGenLength")
    ),
    list(
      model = "bird_environmental_variation",
      relationship = "bird_sigma",
      coef_cols = c(intercept = "alpha", slope = "beta_logGenLength")
    )
  )

  relationship_levels <- vapply(specs, `[[`, character(1), "relationship")
  out <- lapply(specs, function(spec) {
    draws <- draws_by_model[[spec$model]]
    need_cols(draws, unname(spec$coef_cols), paste0(spec$model, " posterior draws"))

    dplyr::bind_rows(lapply(names(spec$coef_cols), function(coef_label) {
      tibble::tibble(
        relationship = spec$relationship,
        coefficient = coef_label,
        value = as.numeric(draws[[spec$coef_cols[[coef_label]]]])
      )
    }))
  })

  dplyr::bind_rows(out) |>
    dplyr::filter(is.finite(value)) |>
    dplyr::mutate(
      relationship = factor(relationship, levels = relationship_levels),
      coefficient = factor(
        coefficient,
        levels = c("intercept", "slope")
      )
    )
}

tidy_bird_sigma_diet_coefficient_draws <- function(draws_by_model) {
  draws <- draws_by_model$bird_environmental_variation
  coefficient_columns <- names(draws)[startsWith(names(draws), "beta_diet_")]
  if (!length(coefficient_columns)) {
    return(tibble::tibble(
      coefficient = factor(character()),
      value = numeric()
    ))
  }
  model_spec <- bird_model_spec_from_posterior(names(draws))
  assert(
    identical(coefficient_columns, model_spec$coefficient_names),
    "Bird diet coefficients are not in canonical Diet-5Cat order."
  )

  dplyr::bind_rows(lapply(coefficient_columns, function(coefficient_column) {
    tibble::tibble(
      coefficient = sub("^beta_diet_", "", coefficient_column),
      value = as.numeric(draws[[coefficient_column]])
    )
  })) |>
    dplyr::filter(is.finite(value)) |>
    dplyr::mutate(
      coefficient = factor(
        .data$coefficient,
        levels = model_spec$separate_intercepts
      )
    )
}

make_coefficient_density_panel <- function(draws, x_limits = NULL, x_breaks = NULL,
                                           density_adjust = 1.05) {
  if (!nrow(draws)) return(cowplot::ggdraw())

  median_value <- stats::median(draws$value)
  x_breaks <- if (is.null(x_breaks)) ggplot2::waiver() else x_breaks
  pal <- methods_figure_palette()

  ggplot2::ggplot(draws, ggplot2::aes(x = value)) +
    ggplot2::geom_density(
      fill = pal[["main_fill"]],
      color = pal[["ink"]],
      linewidth = 0.55,
      alpha = 0.32,
      adjust = density_adjust
    ) +
    ggplot2::geom_vline(
      xintercept = median_value,
      color = pal[["main"]],
      linewidth = 0.72,
      lineend = "round"
    ) +
    ggplot2::scale_x_continuous(
      breaks = x_breaks,
      labels = scales::label_number(accuracy = 0.01, trim = TRUE),
      expand = ggplot2::expansion(mult = c(0.08, 0.08))
    ) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.08))) +
    ggplot2::coord_cartesian(xlim = x_limits, clip = "off") +
    ggplot2::labs(x = NULL, y = NULL) +
    theme_methods_figure(base_size = 10.4, legend_position = "none") +
    ggplot2::theme(
      plot.margin = ggplot2::margin(2.5, 8, 2.5, 8),
      panel.grid = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 8.8, colour = "#333333"),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_line(linewidth = 0.34, colour = "#454545"),
      axis.ticks.y = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(linewidth = 0.42, colour = "#454545"),
      axis.line.y = ggplot2::element_blank(),
      axis.ticks.length = grid::unit(2.1, "pt")
    )
}

coefficient_figure_label <- function(label, size = 9.4, x = 0.5, hjust = 0.5, parse = FALSE) {
  if (!isTRUE(parse)) {
    return(
      cowplot::ggdraw() +
        cowplot::draw_label(label, x = x, y = 0.5, hjust = hjust, fontface = "bold", size = size, colour = "grey15")
    )
  }

  ggplot2::ggplot() +
    ggplot2::annotate(
      "text",
      x = x,
      y = 0.5,
      label = label,
      parse = TRUE,
      hjust = hjust,
      fontface = "bold",
      size = size / 2.85,
      colour = "grey15"
    ) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    ggplot2::theme_void()
}

build_stage1_coefficient_figure <- function(draws_by_model) {
  expected_models <- c(
    "mammal_growth", "mammal_environmental_variation",
    "bird_growth", "bird_environmental_variation"
  )
  assert(
    identical(names(draws_by_model), expected_models),
    "Stage 1 posterior set has unexpected model names or order."
  )

  # Configurable diet effects add columns only to the bird
  # environmental-variation row.
  draws <- tidy_demographic_coefficient_draws(draws_by_model)
  diet_draws <- tidy_bird_sigma_diet_coefficient_draws(draws_by_model)
  diet_groups <- levels(diet_draws$coefficient)
  row_order <- c("mammal_rm", "mammal_sigma", "bird_rm", "bird_sigma")
  row_labels <- c(
    mammal_rm = "atop(Mammal~r[m], body~mass)",
    mammal_sigma = "atop(Mammal~sigma[r], body~mass)",
    bird_rm = "atop(Bird~r[m], generation~length)",
    bird_sigma = "atop(Bird~sigma[r], generation~length)"
  )

  blank <- cowplot::ggdraw()
  panel_count <- 2L + length(diet_groups)
  # The wider label column accommodates unabbreviated trait names; diet-effect
  # columns remain slightly narrower so the figure reproduces cleanly at page width.
  rel_widths <- c(0.58, 1, 1, rep(0.82, length(diet_groups)))
  header_labels <- c(
    "Intercept~'(' * alpha * ')'",
    "Slope~'(' * beta * ')'",
    paste0(diet_groups, " effect")
  )
  header <- cowplot::plot_grid(
    plotlist = c(
      list(blank),
      lapply(seq_along(header_labels), function(i) {
        coefficient_figure_label(
          header_labels[[i]],
          size = 11.0,
          parse = i <= 2L
        )
      })
    ),
    ncol = panel_count + 1L,
    rel_widths = rel_widths
  )

  rows <- lapply(row_order, function(relationship) {
    row_draws <- draws[draws$relationship == relationship, , drop = FALSE]
    diet_panels <- lapply(diet_groups, function(group) {
      if (!identical(relationship, "bird_sigma")) return(blank)
      make_coefficient_density_panel(
        diet_draws[diet_draws$coefficient == group, , drop = FALSE]
      )
    })
    cowplot::plot_grid(
      plotlist = c(
        list(coefficient_figure_label(
          row_labels[[relationship]],
          size = 9.8,
          x = 0.98,
          hjust = 1,
          parse = TRUE
        )),
        lapply(levels(draws$coefficient), function(coefficient) {
          make_coefficient_density_panel(
            row_draws[row_draws$coefficient == coefficient, , drop = FALSE]
          )
        }),
        diet_panels
      ),
      ncol = panel_count + 1L,
      rel_widths = rel_widths,
      align = "h",
      axis = "tb"
    )
  })

  bottom_label <- cowplot::ggdraw() +
    cowplot::draw_label(
      "Posterior coefficient value",
      x = 0.585,
      y = 0.5,
      size = 10.5,
      fontface = "plain",
      colour = "grey15"
    )

  cowplot::plot_grid(
    plotlist = c(list(header), rows, list(bottom_label)),
    ncol = 1,
    rel_heights = c(0.20, rep(1, length(rows)), 0.22)
  )
}
