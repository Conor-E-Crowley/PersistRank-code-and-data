# Stage 7.3 assemblage figure construction.
#
# Loaded by the Stage 7.3 workflow after shared and focal helpers. This
# definition-only module owns assemblage panels, legends, section composition,
# and composite labelling; it performs no work when sourced.

stage73_distribution_panel <- function(prepared, rank_label, point_style = "hollow") {
  point_style <- validate_scalar_choice(
    point_style, report_point_styles(), "Stage 7.3 point_style"
  )
  point_params <- stage73_point_style_params(point_style)
  point_size <- stage73_point_size(point_style, hollow_size = 1.12, solid_size = 1.30)
  show_points <- !identical(point_style, "none")
  palette <- methods_figure_palette()
  method_labels <- c(pipe = "PersistRank", rank = rank_label)
  method_colours <- c(
    "PersistRank" = palette[["main"]],
    stats::setNames(palette[["benchmark"]], rank_label)
  )
  rows <- data.table::copy(prepared$central)
  assert(
    nrow(rows) > 0L &&
      identical(sort(unique(rows$method)), c("pipe", "rank")) &&
      all(is.finite(rows$central) & rows$central >= 0 & rows$central <= 1),
    "Stage 7.3 overlaid species-distribution trajectories are incomplete or invalid."
  )
  rows[, method_lab := factor(method_labels[method], levels = unname(method_labels))]
  inner_rows <- data.table::copy(rows)
  inner_rows[, `:=`(band_min = q25, band_max = q75)]
  bounds <- data.table::melt(
    rows,
    id.vars = c("method", "method_lab", "stage", "stage_order", "x_percent"),
    measure.vars = c("q10", "q90"),
    variable.name = "bound",
    value.name = "persistence"
  )
  data.table::setorder(bounds, method, bound, stage_order)
  data.table::setorder(rows, method, stage_order)
  data.table::setorder(inner_rows, method, stage_order)
  # Mark every completed stage rather than a subsampled subset, so the
  # assemblage panel's points reflect the true resolution of the underlying
  # discrete removal sequence.
  point_rows <- data.table::copy(rows)

  plot <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = inner_rows[method == "rank"],
      ggplot2::aes(x = x_percent, ymin = band_min, ymax = band_max, fill = method_lab),
      alpha = 0.16, colour = NA, show.legend = TRUE
    ) +
    ggplot2::geom_ribbon(
      data = inner_rows[method == "pipe"],
      ggplot2::aes(x = x_percent, ymin = band_min, ymax = band_max, fill = method_lab),
      alpha = 0.16, colour = NA, show.legend = TRUE
    )
  plot +
    ggplot2::geom_line(
      data = bounds[method == "rank"],
      ggplot2::aes(
        x = x_percent, y = persistence, colour = method_lab,
        group = interaction(method_lab, bound)
      ),
      linewidth = 0.48, linetype = "22", alpha = 0.64,
      show.legend = FALSE
    ) +
    ggplot2::geom_line(
      data = bounds[method == "pipe"],
      ggplot2::aes(
        x = x_percent, y = persistence, colour = method_lab,
        group = interaction(method_lab, bound)
      ),
      linewidth = 0.48, linetype = "22", alpha = 0.64,
      show.legend = FALSE
    ) +
    ggplot2::geom_line(
      data = rows[method == "rank"],
      ggplot2::aes(x = x_percent, y = central, colour = method_lab),
      linewidth = 1.12, lineend = "round", show.legend = TRUE
    ) +
    (if (show_points) {
      do.call(ggplot2::geom_point, c(
        list(
          data = point_rows[method == "rank"],
          mapping = ggplot2::aes(x = x_percent, y = central, colour = method_lab),
          size = point_size, alpha = 0.90, show.legend = FALSE
        ),
        point_params
      ))
    } else {
      NULL
    }) +
    ggplot2::geom_line(
      data = rows[method == "pipe"],
      ggplot2::aes(x = x_percent, y = central, colour = method_lab),
      linewidth = 1.12, lineend = "round", show.legend = TRUE
    ) +
    (if (show_points) {
      do.call(ggplot2::geom_point, c(
        list(
          data = point_rows[method == "pipe"],
          mapping = ggplot2::aes(x = x_percent, y = central, colour = method_lab),
          size = point_size, alpha = 0.90, show.legend = FALSE
        ),
        point_params
      ))
    } else {
      NULL
    }) +
    ggplot2::scale_colour_manual(
      values = method_colours,
      breaks = unname(method_labels),
      name = "Method"
    ) +
    ggplot2::scale_fill_manual(values = method_colours, guide = "none") +
    ggplot2::scale_x_continuous(
      limits = c(0, 100),
      breaks = seq(0, 100, by = 20),
      labels = function(x) paste0(x, "%"),
      # Keep the endpoint label comfortably inside the exported canvas.
      expand = ggplot2::expansion(mult = c(0.005, 0.04))
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1), breaks = methods_probability_breaks(),
      labels = methods_probability_labels(),
      expand = ggplot2::expansion(add = c(0.015, 0.015))
    ) +
    ggplot2::labs(
      x = "Cells removed",
      y = "Species persistence"
    ) +
    theme_methods_figure(base_size = 10.1, legend_position = "none") +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(
        colour = palette[["pale_grid"]], linewidth = 0.30
      ),
      panel.grid.major.x = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(
        face = "plain", size = 9.4, colour = palette[["text"]]
      ),
      axis.text = ggplot2::element_text(size = 8.5),
      legend.position = "none",
      plot.margin = ggplot2::margin(8, 10, 2, 8)
    )
}

# Focal-page figure: late-stage or full-sequence focal-species grid --------

prepare_stage73_assemblage_components <- function(
  core,
  figure_data = NULL,
  distribution_data = NULL,
  central_stat = "mean"
) {
  if (is.null(figure_data)) {
    figure_data <- prepare_stage73_figure_data(core)
  }
  if (is.null(distribution_data)) {
    distribution_data <- prepare_stage73_distribution_data(
      core,
      figure_data,
      central_stat = central_stat
    )
  }
  assert(
    data.table::is.data.table(distribution_data$central),
    "Stage 7.3 assemblage components are incomplete."
  )
  list(
    distribution = distribution_data,
    rank_label = validate_scalar_string(
      core$rank_label,
      "Stage 7.3 benchmark label"
    )
  )
}

stage73_assemblage_panel <- function(distribution_data, rank_label, point_style = "hollow") {
  stage73_distribution_panel(distribution_data, rank_label, point_style = point_style) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = 9.9),
      axis.text = ggplot2::element_text(size = 9.0),
      plot.margin = ggplot2::margin(3, 9, 1, 5)
    )
}

stage73_assemblage_legend <- function(rank_label, point_style = "hollow") {
  point_style <- validate_scalar_choice(
    point_style, report_point_styles(), "Stage 7.3 point_style"
  )
  point_params <- stage73_point_style_params(point_style)
  point_size <- stage73_point_size(point_style, hollow_size = 1.25, solid_size = 1.45)
  show_points <- !identical(point_style, "none")
  palette <- methods_figure_palette()
  rank_display <- sub("^Zonation\\s+", "", rank_label)
  ggplot2::ggplot() +
    ggplot2::annotate(
      "segment", x = 0.02, xend = 0.10, y = 0.70, yend = 0.70,
      linewidth = 0.90, colour = palette[["main"]], lineend = "round"
    ) +
    (if (show_points) {
      do.call(ggplot2::annotate, c(
        list("point", x = 0.06, y = 0.70, size = point_size, colour = palette[["main"]]),
        point_params
      ))
    } else {
      NULL
    }) +
    ggplot2::annotate(
      "text", x = 0.12, y = 0.70, label = "PersistRank",
      hjust = 0, size = 2.70, colour = palette[["text"]]
    ) +
    ggplot2::annotate(
      "segment", x = 0.55, xend = 0.63, y = 0.70, yend = 0.70,
      linewidth = 0.90, colour = palette[["benchmark"]], lineend = "round"
    ) +
    (if (show_points) {
      do.call(ggplot2::annotate, c(
        list("point", x = 0.59, y = 0.70, size = point_size, colour = palette[["benchmark"]]),
        point_params
      ))
    } else {
      NULL
    }) +
    ggplot2::annotate(
      "text", x = 0.65, y = 0.70, label = rank_display,
      hjust = 0, size = 2.70, colour = palette[["text"]]
    ) +
    ggplot2::annotate(
      "rect", xmin = 0.02, xmax = 0.06, ymin = 0.15, ymax = 0.35,
      fill = palette[["main"]], alpha = 0.16, colour = NA
    ) +
    ggplot2::annotate(
      "rect", xmin = 0.06, xmax = 0.10, ymin = 0.15, ymax = 0.35,
      fill = palette[["benchmark"]], alpha = 0.16, colour = NA
    ) +
    ggplot2::annotate(
      "text", x = 0.12, y = 0.25, label = "Middle 50%",
      hjust = 0, size = 2.64, colour = palette[["text"]]
    ) +
    ggplot2::annotate(
      "segment", x = 0.55, xend = 0.59, y = 0.25, yend = 0.25,
      linewidth = 0.48, linetype = "22", colour = palette[["main"]]
    ) +
    ggplot2::annotate(
      "segment", x = 0.59, xend = 0.63, y = 0.25, yend = 0.25,
      linewidth = 0.48, linetype = "22", colour = palette[["benchmark"]]
    ) +
    ggplot2::annotate(
      "text", x = 0.65, y = 0.25, label = "10th / 90th percentiles",
      hjust = 0, size = 2.64, colour = palette[["text"]]
    ) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    ggplot2::theme_void()
}

stage73_assemblage_section <- function(prepared_components, point_style = "hollow") {
  point_style <- validate_scalar_choice(
    point_style, report_point_styles(), "Stage 7.3 point_style"
  )
  assert(
    is.list(prepared_components) &&
      data.table::is.data.table(prepared_components$distribution$central),
    "Prepared Stage 7.3 assemblage figure components are invalid."
  )
  panel <- stage73_assemblage_panel(
    prepared_components$distribution,
    prepared_components$rank_label,
    point_style = point_style
  )
  legend <- stage73_assemblage_legend(prepared_components$rank_label, point_style = point_style)
  out <- cowplot::plot_grid(
    panel,
    legend,
    ncol = 1L,
    rel_heights = c(2.28, 0.37)
  )
  out
}

stage73_label_assemblage_for_composite <- function(assemblage_figure) {
  palette <- methods_figure_palette()
  assemblage_figure +
    ggplot2::labs(
      title = "Assemblage-level persistence",
      tag = "(a)"
    ) +
    ggplot2::theme(
      plot.title.position = "plot",
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 10.1,
        hjust = 0.5,
        colour = palette[["text"]],
        margin = ggplot2::margin(b = 1)
      ),
      plot.tag = ggplot2::element_text(
        face = "bold",
        size = 9.4,
        hjust = 0,
        vjust = 1,
        colour = palette[["text"]],
        margin = ggplot2::margin(r = 3)
      ),
      # Keep the composite's first panel label on the title row, but inset it
      # enough to associate it with the assemblage plot and avoid edge clipping.
      plot.tag.position = c(0.02, 0.975),
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )
}
