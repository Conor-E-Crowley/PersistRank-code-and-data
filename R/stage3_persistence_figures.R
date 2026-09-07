# Abundance/persistence heatmap presentation for Stage 3.
#
# Loaded by stage3_workflow.R after common figure helpers. It builds taxon-level
# panels from prepared LOESS predictions and combines them before the single
# manuscript figure is rendered. Sourcing and construction perform no I/O.

normalize_heatmap_taxon <- function(taxon) {
  taxon <- tolower(validate_scalar_string(taxon, "taxon"))
  if (taxon %in% c("mammal", "mammals")) {
    return(list(
      group = "Mammals",
      title = "Mammals",
      trait_label = "Body\nmass",
      labels = fmt_mass_axis,
      palette = "viridis"
    ))
  }
  if (taxon %in% c("bird", "birds")) {
    return(list(
      group = "Birds",
      title = "Birds",
      trait_label = "Generation\nlength",
      labels = fmt_generation_axis,
      palette = "plasma"
    ))
  }
  stop("taxon must be 'mammals' or 'birds'.", call. = FALSE)
}

build_abundance_persistence_heatmap <- function(loess_pred,
                                                k0,
                                                taxon = "mammals",
                                                curve = main_persistence_curve(),
                                                trait_n = 720,
                                                k_n = 1100,
                                                p_n = 1200,
                                                target_p = 0.995,
                                                bird_sigma_model_group = NULL) {
  k0 <- validate_scalar_number(
    k0, "k0", minimum = 0, minimum_open = TRUE
  )
  curve <- validate_persistence_curve(curve)
  taxon_spec <- normalize_heatmap_taxon(taxon)
  if (identical(taxon_spec$group, "Birds")) {
    bird_sigma_model_group <- validate_scalar_string(
      bird_sigma_model_group,
      "bird_sigma_model_group"
    )
    need_cols(loess_pred, "bird_sigma_model_group", "LOESS prediction grid for abundance-persistence heatmap")
  }
  need_cols(
    loess_pred,
    c("group", "trait_value", "curve", "alpha_mid", "beta_mid"),
    "LOESS prediction grid for abundance-persistence heatmap"
  )
  validate_loess_figure_predictions(loess_pred)

  pred <- loess_pred |>
    dplyr::filter(
      .data$group == taxon_spec$group,
      .data$curve == .env$curve
  )
  if (identical(taxon_spec$group, "Birds")) {
    pred <- pred |>
      dplyr::filter(.data$bird_sigma_model_group == .env$bird_sigma_model_group)
  }
  pred <- pred |>
    dplyr::arrange(.data$trait_value)

  assert(
    nrow(pred) > 0L,
    paste0(taxon_spec$title, " LOESS predictions are unavailable for curve ", curve, ".")
  )

  if (nrow(pred) > trait_n) {
    keep <- unique(round(seq(1, nrow(pred), length.out = trait_n)))
    pred <- pred[keep, , drop = FALSE]
  }

  k_target <- k0 + (pred$alpha_mid / -log(target_p))^(1 / pred$beta_mid)
  k_max <- max(k_target[is.finite(k_target)], k0 * 4, na.rm = TRUE) * 1.05
  x_max <- if (
    identical(taxon_spec$group, "Birds") &&
      identical(bird_sigma_model_group, "Other")
  ) {
    k_max * 1.20
  } else {
    k_max
  }
  k_min <- left_limit_for_threshold_fraction(k0, x_max, fraction = 0.065)
  k_breaks <- scales::breaks_log(n = 5)(c(k_min, x_max))
  log_k_grid <- seq(log10(k_min), log10(k_max), length.out = k_n)
  k_grid <- 10^log_k_grid
  p_grid <- seq(0, 1, length.out = p_n)

  heat <- lapply(seq_along(k_grid), function(i) {
    k <- k_grid[[i]]
    trait_at_p <- rep(NA_real_, length(p_grid))

    if (k > k0) {
      p_at_trait <- exp(-pred$alpha_mid * pmax(k - k0, 1e-12)^(-pred$beta_mid))
      ok <- is.finite(p_at_trait) & is.finite(pred$trait_value)
      if (length(unique(p_at_trait[ok])) >= 2L) {
        trait_at_p <- stats::approx(
          x = p_at_trait[ok],
          y = pred$trait_value[ok],
          xout = p_grid,
          ties = mean,
          rule = 1
        )$y
      }
    }

    tibble::tibble(
      logK = log_k_grid[[i]],
      persistence = p_grid,
      trait_value = trait_at_p
    )
  }) |>
    dplyr::bind_rows()

  pal <- methods_figure_palette()
  trait_range <- range(pred$trait_value, finite = TRUE)
  trait_breaks <- scales::breaks_log(n = 4)(trait_range)
  trait_breaks <- trait_breaks[trait_breaks >= trait_range[1] & trait_breaks <= trait_range[2]]
  plot_title <- if (identical(taxon_spec$group, "Birds")) {
    if (identical(bird_sigma_model_group, "Other")) {
      "Birds — other diets"
    } else {
      paste0("Birds — ", bird_sigma_model_group)
    }
  } else {
    "Mammals"
  }

  ggplot2::ggplot(
    heat,
    ggplot2::aes(
      x = .data$logK,
      y = .data$persistence,
      fill = .data$trait_value
    )
  ) +
    ggplot2::geom_raster(na.rm = TRUE, interpolate = TRUE) +
    ggplot2::geom_vline(
      xintercept = log10(k0),
      linewidth = 0.58,
      linetype = "dotted",
      colour = pal[["threshold"]]
    ) +
    ggplot2::annotate(
      "text",
      x = log10(k0) + 0.012 * diff(log10(c(k_min, x_max))),
      y = 0.975,
      label = paste0("K\u2080 = ", scales::label_comma()(k0)),
      hjust = 0,
      vjust = 1,
      size = 4.0,
      colour = pal[["threshold"]]
    ) +
    ggplot2::scale_x_continuous(
      limits = log10(c(k_min, x_max)),
      breaks = log10(k_breaks),
      labels = methods_compact_number_labels(k_breaks),
      expand = ggplot2::expansion(mult = c(0, 0.02))
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      breaks = methods_probability_breaks(),
      labels = methods_probability_labels(),
      expand = ggplot2::expansion(mult = c(0, 0.01))
    ) +
    ggplot2::scale_fill_viridis_c(
      option = taxon_spec$palette,
      trans = "log10",
      breaks = trait_breaks,
      labels = taxon_spec$labels,
      name = taxon_spec$trait_label,
      na.value = "transparent"
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        direction = "vertical",
        barheight = grid::unit(2.00, "cm"),
        barwidth = grid::unit(0.25, "cm"),
        ticks.colour = pal[["axis"]],
        frame.colour = pal[["panel_border"]]
      )
    ) +
    ggplot2::coord_cartesian(ylim = c(0, 1), clip = "off") +
    ggplot2::labs(
      title = plot_title,
      x = "Abundance K",
      y = "Persistence probability"
    ) +
    theme_methods_figure(base_size = 15.5, legend_position = "inside") +
    ggplot2::theme(
      plot.margin = ggplot2::margin(6, 7, 6, 7),
      panel.grid.major.y = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14.5),
      axis.title = ggplot2::element_text(size = 13.5),
      axis.text = ggplot2::element_text(size = 11.8),
      legend.position.inside = c(0.85, 0.29),
      legend.justification = c(0.5, 0.5),
      legend.background = ggplot2::element_rect(
        fill = scales::alpha("white", 0.92),
        colour = pal[["panel_border"]],
        linewidth = 0.25
      ),
      legend.title = ggplot2::element_text(face = "bold", size = 10.2),
      legend.text = ggplot2::element_text(size = 9.8),
      legend.margin = ggplot2::margin(1.5, 2.5, 1.5, 2.5),
      legend.box.just = "center"
    )
}

build_abundance_persistence_figure <- function(loess_pred,
                                               k0,
                                               curve,
                                               include_mammals,
                                               bird_model_groups = character()) {
  include_mammals <- validate_scalar_logical(
    include_mammals, "include_mammals"
  )
  bird_model_groups <- as.character(bird_model_groups)
  assert(
    all(!is.na(bird_model_groups) & nzchar(trimws(bird_model_groups))) &&
      !anyDuplicated(bird_model_groups),
    "bird_model_groups must contain unique non-empty names."
  )

  # The workflow supplies model groups in manuscript order. Build every panel
  # as a plot object so the complete figure is rendered only once.
  panels <- list()
  if (include_mammals) {
    panels$mammals <- build_abundance_persistence_heatmap(
      loess_pred,
      k0 = k0,
      taxon = "mammals",
      curve = curve
    )
  }
  for (model_group in bird_model_groups) {
    panels[[model_group]] <- build_abundance_persistence_heatmap(
      loess_pred,
      k0 = k0,
      taxon = "birds",
      curve = curve,
      bird_sigma_model_group = model_group
    )
  }
  assert(length(panels) > 0L, "No abundance-persistence panels were selected.")

  panel_tag_grid(
    plotlist = unname(panels),
    labels = paste0("(", letters[seq_along(panels)], ")"),
    ncol = length(panels),
    label_size = 14.5
  )
}
