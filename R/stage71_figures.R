# Stage 7.1 plot construction.
#
# The Stage 7.1 workflow owns loading this definition-only module after the
# scientific and report layers. Expensive plotting dependencies are checked
# only when a builder runs; sourcing this file performs no reads or writes.
stage71_trajectory_band_colours <- function(palette = methods_figure_palette()) {
  c(
    outer = grDevices::colorRampPalette(
      c(palette[["land"]], palette[["main_fill"]])
    )(5L)[[2L]],
    inner = grDevices::colorRampPalette(
      c(palette[["land"]], palette[["main"]])
    )(5L)[[3L]]
  )
}

stage71_trajectory_theme <- function(palette = methods_figure_palette()) {
  theme_methods_figure(base_size = 11, legend_position = "bottom") +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(
        colour = palette[["pale_grid"]],
        linewidth = 0.28
      ),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(
        face = "bold",
        colour = palette[["text"]],
        size = 9.4,
        margin = ggplot2::margin(1, 2, 4, 2)
      ),
      # Extra separation prevents 100%/0% boundary labels from colliding when
      # the three-facet figure is reduced to a single manuscript column.
      panel.spacing.x = grid::unit(1.8, "lines"),
      legend.box = "horizontal",
      legend.box.just = "center",
      legend.margin = ggplot2::margin(t = 3, b = 0),
      legend.box.margin = ggplot2::margin(0, 0, 0, 0),
      legend.key.width = grid::unit(1.25, "lines"),
      legend.spacing.x = grid::unit(0.5, "lines"),
      plot.margin = ggplot2::margin(5, 6, 4, 6)
    )
}

stage71_plot_scope_levels <- function(scopes) {
  scopes <- unique(as.character(scopes))
  taxon_scopes <- intersect(c("mammals", "birds"), scopes)
  if (length(taxon_scopes) == 1L && "all" %in% scopes) {
    scopes <- setdiff(scopes, "all")
  }
  canonical <- c("all", "mammals", "birds")
  canonical[canonical %in% scopes]
}

build_stage71_persistence_summary_plot <- function(summary) {
  assert(requireNamespace("ggplot2", quietly = TRUE), "ggplot2 is required for Stage 7.1 persistence plots.")
  dt <- data.table::copy(data.table::as.data.table(summary))
  need_cols(
    dt,
    c(
      "scope", "pct_cells_removed_end", "mean_persist",
      "q10_persist", "q25_persist", "q75_persist", "q90_persist"
    ),
    "Stage 7.1 community persistence summary"
  )
  plot_scopes <- stage71_plot_scope_levels(dt$scope)
  dt <- dt[scope %in% plot_scopes]
  dt[, scope_label := factor(
    stage71_scope_label(scope),
    levels = stage71_scope_label(plot_scopes)
  )]
  pal <- methods_figure_palette()
  ribbon_data <- dt[, .(
    pct_cells_removed_end, scope_label,
    lower = q25_persist, upper = q75_persist
  )]
  quantile_data <- data.table::rbindlist(list(
    dt[, .(
      pct_cells_removed_end, scope_label,
      persistence = q10_persist,
      quantile = "10th percentile"
    )],
    dt[, .(
      pct_cells_removed_end, scope_label,
      persistence = q90_persist,
      quantile = "90th percentile"
    )]
  ))
  band_colours <- stage71_trajectory_band_colours(pal)

  ggplot2::ggplot(dt, ggplot2::aes(x = pct_cells_removed_end)) +
    ggplot2::geom_ribbon(
      data = ribbon_data,
      ggplot2::aes(ymin = lower, ymax = upper, fill = "Middle 50%"),
      colour = NA,
      alpha = 0.70
    ) +
    ggplot2::geom_line(
      data = quantile_data,
      ggplot2::aes(
        y = persistence,
        linetype = "10th / 90th percentiles",
        group = interaction(scope_label, quantile)
      ),
      colour = pal[["main"]],
      linewidth = 0.52,
      alpha = 0.78
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = mean_persist, linetype = "Mean"),
      colour = pal[["main"]], linewidth = 1.0
    ) +
    ggplot2::facet_wrap(~scope_label, nrow = 1) +
    ggplot2::scale_fill_manual(
      values = c("Middle 50%" = band_colours[["inner"]]),
      breaks = "Middle 50%",
      name = NULL
    ) +
    ggplot2::scale_linetype_manual(
      values = c(Mean = "solid", "10th / 90th percentiles" = "22"),
      breaks = c("Mean", "10th / 90th percentiles"),
      name = NULL
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, 100), breaks = seq(0, 100, 20),
      labels = function(x) paste0(x, "%"),
      expand = ggplot2::expansion(mult = c(0, 0.008))
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      breaks = methods_probability_breaks(),
      labels = methods_probability_labels(),
      expand = ggplot2::expansion(mult = c(0, 0.012))
    ) +
    ggplot2::labs(
      x = "Cells removed",
      y = "Species persistence"
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(order = 1, override.aes = list(colour = NA)),
      linetype = ggplot2::guide_legend(order = 2)
    ) +
    stage71_trajectory_theme(pal)
}

build_stage71_mean_curve_uncertainty_plot <- function(curve_summary) {
  assert(requireNamespace("ggplot2", quietly = TRUE), "ggplot2 is required for Stage 7.1 persistence plots.")
  dt <- data.table::copy(data.table::as.data.table(curve_summary))
  need_cols(
    dt,
    c(
      "stage", "stage_order", "pct_cells_removed_end", "scope",
      "curve", "n_species", "mean_persist"
    ),
    "Stage 7.1 mean uncertainty-curve summary"
  )
  assert(
    setequal(unique(dt$curve), persistence_curves()) &&
      !anyDuplicated(dt[, .(stage, scope, curve)]) &&
      all(is.finite(dt$mean_persist) & dt$mean_persist >= 0 & dt$mean_persist <= 1),
    "Stage 7.1 mean uncertainty-curve summary is invalid."
  )
  wide <- data.table::dcast(
    dt,
    stage + stage_order + pct_cells_removed_end + scope + n_species ~ curve,
    value.var = "mean_persist"
  )
  need_cols(wide, persistence_curves(), "Stage 7.1 mean uncertainty-curve summary")
  wide[, `:=`(
    outer_ymin = pmin(q025, q975),
    outer_ymax = pmax(q025, q975),
    inner_ymin = pmin(q16, q84),
    inner_ymax = pmax(q16, q84)
  )]
  plot_scopes <- stage71_plot_scope_levels(wide$scope)
  wide <- wide[scope %in% plot_scopes]
  wide[, scope_label := factor(
    stage71_scope_label(scope),
    levels = stage71_scope_label(plot_scopes)
  )]
  pal <- methods_figure_palette()
  band_colours <- stage71_trajectory_band_colours(pal)
  bands <- data.table::rbindlist(list(
    wide[, .(
      stage, stage_order, pct_cells_removed_end, scope_label,
      ymin = outer_ymin, ymax = outer_ymax,
      interval = "95% range (q2.5–q97.5)"
    )],
    wide[, .(
      stage, stage_order, pct_cells_removed_end, scope_label,
      ymin = inner_ymin, ymax = inner_ymax,
      interval = "68% range (q16–q84)"
    )]
  ))
  bands[, interval := factor(
    interval,
    levels = c("95% range (q2.5–q97.5)", "68% range (q16–q84)")
  )]

  ggplot2::ggplot(wide, ggplot2::aes(x = pct_cells_removed_end)) +
    ggplot2::geom_ribbon(
      data = bands,
      ggplot2::aes(ymin = ymin, ymax = ymax, fill = interval),
      colour = NA
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = q50, colour = "q50 mean"),
      linewidth = 1.0
    ) +
    ggplot2::facet_wrap(~scope_label, nrow = 1) +
    ggplot2::scale_fill_manual(
      values = c(
        "95% range (q2.5–q97.5)" = band_colours[["outer"]],
        "68% range (q16–q84)" = band_colours[["inner"]]
      ),
      breaks = c("68% range (q16–q84)", "95% range (q2.5–q97.5)"),
      name = NULL
    ) +
    ggplot2::scale_colour_manual(
      values = c("q50 mean" = pal[["main"]]),
      breaks = "q50 mean",
      name = NULL
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, 100), breaks = seq(0, 100, 20),
      labels = function(x) paste0(x, "%"),
      expand = ggplot2::expansion(mult = c(0, 0.008))
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      breaks = methods_probability_breaks(),
      labels = methods_probability_labels(),
      expand = ggplot2::expansion(mult = c(0, 0.012))
    ) +
    ggplot2::labs(
      x = "Cells removed",
      y = "Mean persistence"
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(order = 1, override.aes = list(colour = NA)),
      colour = ggplot2::guide_legend(
        order = 2,
        override.aes = list(colour = pal[["main"]], linewidth = 0.95)
      )
    ) +
    stage71_trajectory_theme(pal)
}

build_stage71_persistence_heatmap <- function(species) {
  assert(requireNamespace("ggplot2", quietly = TRUE), "ggplot2 is required for Stage 7.1 persistence plots.")
  dt <- data.table::copy(data.table::as.data.table(species))
  changes <- stage71_persistence_changes(dt)
  ordering <- changes[, .(
    taxon,
    scientificName,
    normalized_persistence_auc,
    latest_persistence
  )]
  data.table::setorder(
    ordering,
    taxon,
    normalized_persistence_auc,
    latest_persistence,
    scientificName
  )
  dt[, species_label := factor(scientificName, levels = ordering$scientificName)]
  dt[, taxon_label := factor(
    stage71_scope_label(taxon),
    levels = stage71_scope_label(c("mammals", "birds"))
  )]

  stage_key <- unique(dt[, .(stage_order, pct_cells_removed_end)])
  data.table::setorder(stage_key, stage_order)
  break_targets <- seq(0, 100, by = 20)
  break_rows <- unique(vapply(
    break_targets,
    function(target) which.min(abs(stage_key$pct_cells_removed_end - target)),
    integer(1L)
  ))
  breaks <- stage_key$stage_order[break_rows]
  labels <- paste0(round(stage_key$pct_cells_removed_end[break_rows]), "%")
  pal <- methods_figure_palette()

  ggplot2::ggplot(dt, ggplot2::aes(x = stage_order, y = species_label, fill = sp_persist)) +
    ggplot2::geom_tile(width = 1.01, height = 1.01, colour = NA) +
    ggplot2::facet_wrap(~taxon_label, nrow = 1, scales = "free_y") +
    ggplot2::scale_x_continuous(breaks = breaks, labels = labels, expand = c(0, 0)) +
    ggplot2::scale_fill_gradientn(
      colours = c(
        pal[["pale_grid"]], "#D7E4E9",
        pal[["main_fill"]], pal[["main"]]
      ),
      values = c(0, 0.25, 0.5, 1),
      limits = c(0, 1),
      breaks = c(0, 0.5, 1),
      labels = c("0", "0.5", "1"),
      name = "Species persistence",
      guide = ggplot2::guide_colourbar(
        title.position = "top", title.hjust = 0.5,
        barwidth = grid::unit(4.4, "cm"),
        barheight = grid::unit(0.32, "cm")
      )
    ) +
    ggplot2::labs(
      x = "Completed stage (cells removed)",
      y = NULL
    ) +
    theme_methods_figure(base_size = 9.5, legend_position = "bottom") +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(
        size = 6.2,
        face = "italic",
        colour = pal[["muted_text"]],
        margin = ggplot2::margin(r = 2)
      ),
      axis.text.x = ggplot2::element_text(size = 7.8),
      axis.title.x = ggplot2::element_text(size = 9),
      panel.spacing.x = grid::unit(0.8, "lines"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(
        face = "bold",
        colour = pal[["text"]],
        size = 9
      ),
      legend.margin = ggplot2::margin(t = 3, b = 0),
      plot.margin = ggplot2::margin(5, 5, 4, 5)
    )
}
