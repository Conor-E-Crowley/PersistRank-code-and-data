# Stage 7.3 focal-species figure construction.
#
# Loaded by the Stage 7.3 workflow after shared figure-data preparation. This
# definition-only module owns focal panels, headers, strips, legends, grids,
# and focal-section composition; it performs no work when sourced.

stage73_point_style_params <- function(point_style) {
  point_style <- validate_scalar_choice(
    point_style, report_point_styles(), "Stage 7.3 point_style"
  )
  if (identical(point_style, "hollow")) {
    list(shape = 21, fill = "white")
  } else {
    list(shape = 19)
  }
}

# A hollow point's white interior contrasts against the line beneath it, so
# even a modestly sized ring reads clearly as a discrete marker. A solid
# point is the same colour as its line and has no such contrast to rely on,
# so it only becomes visible where its diameter clearly exceeds the line's
# rendered width; each point-bearing panel/legend therefore uses a *larger*
# size for "solid" than for "hollow" to keep both symbologies legible.
stage73_point_size <- function(point_style, hollow_size, solid_size) {
  if (identical(point_style, "hollow")) hollow_size else solid_size
}

prepare_stage73_focal_panels <- function(
  core,
  focus_species,
  figure_data = NULL,
  area_stage_ids,
  minimum_removed_percent = 50,
  persistence_area_reference_removed_percent = minimum_removed_percent,
  area_target_interval_percent = 5,
  point_style = "hollow"
) {
  suppressPackageStartupMessages(library(data.table))
  suppressPackageStartupMessages(library(ggplot2))
  point_style <- validate_scalar_choice(
    point_style, report_point_styles(), "Stage 7.3 point_style"
  )
  point_params <- stage73_point_style_params(point_style)
  point_size <- stage73_point_size(point_style, hollow_size = 0.85, solid_size = 1.05)
  show_points <- !identical(point_style, "none")
  if (is.null(figure_data)) {
    figure_data <- prepare_stage73_figure_data(core)
  }
  focus_species <- validate_scalar_string(focus_species, "focus_species")
  focal_id <- species_id(focus_species)
  rank_label <- core$rank_label
  palette <- methods_figure_palette()

  method_labels <- c(
    pipe = "PersistRank",
    rank = rank_label
  )
  method_colours <- c(
    "PersistRank" = palette[["main"]],
    stats::setNames(palette[["benchmark"]], rank_label)
  )
  method_factor <- function(x) {
    # Draw the benchmark first so persistence-based trajectories remain visible.
    factor(
      unname(method_labels[as.character(x)]),
      levels = rev(unname(method_labels))
    )
  }

  styles <- data.table::copy(figure_data$uncertainty_styles)
  interval_alphas <- stats::setNames(styles$alpha, styles$interval_label)
  outer_label <- styles$interval_label[[1L]]
  inner_label <- styles$interval_label[[2L]]
  available_stage_key <- unique(figure_data$species[, .(
    stage = as.integer(stage),
    stage_order = as.integer(stage_order),
    x = as.numeric(x),
    x_percent = as.numeric(x_percent)
  )])
  data.table::setorder(available_stage_key, stage_order)
  normalize_stage_ids <- function(values, label) {
    numeric_values <- suppressWarnings(as.numeric(values))
    assert(
      length(numeric_values) > 0L &&
        all(is.finite(numeric_values) & numeric_values >= 0 & numeric_values == floor(numeric_values)),
      paste0(label, " must contain positive-integer or zero stage IDs.")
    )
    values <- sort(unique(as.integer(numeric_values)))
    assert(
      all(values %in% available_stage_key$stage),
      paste0(label, " contains a stage absent from the Stage 7.3 trajectory.")
    )
    values
  }
  area_stage_ids <- normalize_stage_ids(
    area_stage_ids,
    "Stage 7.3 retained-area stages"
  )
  area_stage_key <- available_stage_key[match(area_stage_ids, stage)]
  panel_min_removed_percent <- suppressWarnings(as.numeric(
    minimum_removed_percent
  ))
  area_reference_removed_percent <- suppressWarnings(as.numeric(
    persistence_area_reference_removed_percent
  ))
  area_target_interval_percent <- suppressWarnings(as.numeric(
    area_target_interval_percent
  ))
  assert(
    length(panel_min_removed_percent) == 1L &&
      is.finite(panel_min_removed_percent) &&
      panel_min_removed_percent >= 0 &&
      panel_min_removed_percent <= 100,
    "Stage 7.3 focal-panel removal threshold must be a percentage in [0,100]."
  )
  assert(
    length(area_reference_removed_percent) == 1L &&
      is.finite(area_reference_removed_percent) &&
      area_reference_removed_percent >= 0 &&
      area_reference_removed_percent <= 100,
    "Stage 7.3 persistence-area reference threshold must be a percentage in [0,100]."
  )
  assert(
    length(area_target_interval_percent) == 1L &&
      (
        is.na(area_target_interval_percent) ||
          (
            is.finite(area_target_interval_percent) &&
              area_target_interval_percent > 0 &&
              area_target_interval_percent <= 100
          )
      ),
    "Stage 7.3 retained-area target interval must be NA or in (0,100]."
  )
  persistence_breaks <- methods_probability_breaks()
  compact_probability_labels <- methods_probability_labels()

  theme_persist <- function(base_size = 10.5, legend_position = "none",
                            probability_guides = TRUE) {
    out <- theme_methods_figure(
      base_size = base_size,
      legend_position = legend_position
    ) +
      theme(
        plot.title = element_text(
          face = "bold",
          hjust = 0,
          size = base_size * 0.94,
          margin = margin(b = 1)
        ),
        plot.subtitle = element_text(
          hjust = 0,
          size = base_size * 0.74
        ),
        plot.tag = element_text(face = "bold", size = base_size * 1.08),
        plot.tag.position = "topleft",
        plot.title.position = "plot",
        plot.margin = margin(6, 6, 5, 8),
        legend.box = "vertical",
        legend.title = element_text(face = "bold")
      )
    if (isTRUE(probability_guides)) {
      out <- out + theme_methods_probability_guides()
    }
    out
  }

  add_line_scales <- function(plot, x_scale) {
    plot +
      scale_color_manual(
        values = method_colours,
        breaks = unname(method_labels),
        name = "Method"
      ) +
      scale_fill_manual(
        values = method_colours,
        breaks = unname(method_labels),
        name = "Method",
        guide = "none"
      ) +
      scale_alpha_manual(
        values = interval_alphas,
        breaks = rev(styles$interval_label),
        name = "Gompertz ranges"
      ) +
      x_scale +
      scale_y_continuous(
        limits = c(0, 1),
        breaks = persistence_breaks,
        labels = compact_probability_labels,
        # Point-bearing styles ("hollow", "solid") mark every stage with a
        # point, so the top of the panel needs enough headroom that a point
        # drawn at persistence = 1 isn't clipped by the panel border. The
        # "none" style has no points on these panels, so that extra headroom
        # is unnecessary and the tighter expansion is used instead.
        expand = expansion(mult = c(0, if (show_points) 0.03 else 0.01))
      )
  }

  focus_lines <- data.table::copy(figure_data$species[species == focal_id])
  assert(nrow(focus_lines) > 0L, paste0("Species not found in figure trajectories: ", focus_species))
  focus_uncertainty <- stage73_uncertainty_bands(
    focus_lines,
    value_col = "sp_persist",
    id_cols = c(
      "method", "stage", "stage_order", "x", "x_percent",
      "scientificName", "species", "className", "n_pu", "total_pu_area_km2"
    )
  )
  focus_central <- focus_uncertainty$central
  focus_bands <- focus_uncertainty$bands
  focus_central[, method_lab := method_factor(method)]
  focus_bands[, method_lab := method_factor(method)]
  data.table::setorder(focus_central, method_lab, stage_order)
  data.table::setorder(focus_bands, interval, method_lab, stage_order)
  focus_central_display <- focus_central[x_percent >= panel_min_removed_percent]
  focus_bands_display <- focus_bands[x_percent >= panel_min_removed_percent]
  assert(
    nrow(focus_central_display) > 0L && nrow(focus_bands_display) > 0L,
    paste0(
      "Stage 7.3 focal persistence panel has no states at or above ",
      panel_min_removed_percent,
      "% cell removal."
    )
  )
  # Mark every completed stage in the displayed window rather than a
  # 5%-nearest subset, so points reflect the true resolution of the
  # underlying discrete removal sequence. Point markers are omitted entirely
  # under the "none" symbology (see show_points above).
  focus_points <- focus_central_display

  p_focus <- ggplot() +
    geom_ribbon(
      data = focus_bands_display[as.character(interval) == outer_label],
      aes(
        x = x_percent,
        ymin = ymin,
        ymax = ymax,
        fill = method_lab,
        alpha = interval,
        group = method_lab
      ),
      color = NA
    ) +
    geom_ribbon(
      data = focus_bands_display[as.character(interval) == inner_label],
      aes(
        x = x_percent,
        ymin = ymin,
        ymax = ymax,
        fill = method_lab,
        alpha = interval,
        group = method_lab
      ),
      color = NA
    ) +
    geom_line(
      data = focus_central_display[method == "rank"],
      aes(x = x_percent, y = central, color = method_lab, group = method_lab),
      linewidth = 1.05,
      lineend = "round"
    ) +
    (if (show_points) {
      do.call(geom_point, c(
        list(
          data = focus_points[method == "rank"],
          mapping = aes(x = x_percent, y = central, color = method_lab),
          size = point_size,
          alpha = 0.95
        ),
        point_params
      ))
    } else {
      NULL
    }) +
    geom_line(
      data = focus_central_display[method == "pipe"],
      aes(x = x_percent, y = central, color = method_lab, group = method_lab),
      linewidth = 1.05,
      lineend = "round"
    ) +
    (if (show_points) {
      do.call(geom_point, c(
        list(
          data = focus_points[method == "pipe"],
          mapping = aes(x = x_percent, y = central, color = method_lab),
          size = point_size,
          alpha = 0.95
        ),
        point_params
      ))
    } else {
      NULL
    }) +
    labs(
      x = "Cells removed (%)",
      y = "Species persistence",
      title = "Species persistence through removal",
      subtitle = paste0(
        focus_species,
        "; shown after at least ",
        format(panel_min_removed_percent, trim = TRUE),
        "% cell removal"
      )
    ) +
    theme_persist()
  p_focus <- add_line_scales(
    p_focus,
    scale_x_continuous(
      limits = c(panel_min_removed_percent, 100),
      breaks = pretty(c(panel_min_removed_percent, 100), n = 6),
      labels = scales::label_number(accuracy = 1),
      expand = expansion(mult = c(0.005, 0.03))
    )
  ) + theme(legend.position = "none")

  q50_parameter_rows <- which(
    core$curve_parameters$curve == "q50" &
      core$curve_parameters$species == focal_id
  )
  q50_parameters <- core$curve_parameters[q50_parameter_rows]
  assert(nrow(q50_parameters) == 1L, paste0("Missing q50 parameters for ", focus_species, "."))
  area_rows <- which(
    core$pu_long$species == focal_id &
      core$pu_long$stage %in% area_stage_ids
  )
  area_plot <- core$pu_long[
    area_rows,
    .(
      method = as.character(method),
      stage = as.integer(stage),
      stage_order = as.integer(stage_order),
      x = as.numeric(x),
      pu_id = as.integer(pu_id),
      pu_area_km2 = as.numeric(pu_area_km2)
    )
  ]
  assert(nrow(area_plot) > 0L, paste0("Species not found in PU trajectories: ", focus_species))
  area_plot[, q50_pu_persistence := stage7_pu_persistence(
    area_km2 = pu_area_km2,
    density = q50_parameters$density[[1L]],
    alpha = q50_parameters$alpha[[1L]],
    beta = q50_parameters$beta[[1L]],
    threshold_area = q50_parameters$c_th[[1L]],
    label = paste0("Stage 7.3 q50 PU persistence for ", focus_species)
  )]
  area_plot[, stack_score := {
    area_max <- max(pu_area_km2)
    area_scaled <- if (is.finite(area_max) && area_max > 0) pu_area_km2 / area_max else 0
    0.5 * area_scaled + 0.5 * q50_pu_persistence
  }, by = .(method, stage)]
  data.table::setorder(
    area_plot,
    method, stage_order, -stack_score, -pu_area_km2, -q50_pu_persistence, pu_id
  )
  area_plot[, `:=`(
    ymin = c(0, head(cumsum(pu_area_km2), -1L)),
    ymax = cumsum(pu_area_km2)
  ), by = .(method, stage)]

  stage_show <- area_stage_key[order(stage_order), stage]
  stage_label_interval <- if (is.finite(area_target_interval_percent)) {
    area_target_interval_percent
  } else {
    5
  }
  stage_show_labels <- as.character(
    stage_label_interval * round(
      area_stage_key[order(stage_order), x_percent] / stage_label_interval
    )
  )
  stage_label_positions <- seq(
    1L,
    length(stage_show_labels),
    by = 2L
  )
  area_plot[, stage_index := match(stage, stage_show)]
  assert(
    all(!is.na(area_plot$stage_index)),
    "Stage 7.3 PU-area rows do not match the representative-stage sequence."
  )
  area_plot[, x_plot := stage_index + data.table::fifelse(method == "pipe", -0.22, 0.22)]
  bar_half_width <- 0.19
  area_plot[, `:=`(
    xmin = x_plot - bar_half_width,
    xmax = x_plot + bar_half_width
  )]
  bar_outline <- area_plot[pu_area_km2 > 0, .(
    xmin = min(xmin),
    xmax = max(xmax),
    ymin = 0,
    ymax = max(ymax)
  ), by = .(method, stage, stage_order)]
  y_max <- max(area_plot$ymax, na.rm = TRUE)
  assert(is.finite(y_max) && y_max > 0, paste0("No retained PU area for ", focus_species, "."))

  pipe_fill_cols <- c(
    "#D9E8F2", "#C8DDEB", "#AFCDE1", "#82B0D1", "#5794C0", "#2E78AF"
  )
  rank_fill_cols <- c(
    "#F6DDD0", "#F2CDB9", "#ECB397", "#E38E62", "#D36E39", "#B7561D"
  )
  pipe_pal <- scales::col_numeric(pipe_fill_cols, domain = c(0, 1), na.color = pipe_fill_cols[[1L]])
  rank_pal <- scales::col_numeric(rank_fill_cols, domain = c(0, 1), na.color = rank_fill_cols[[1L]])
  pu_bar_style <- stage73_pu_bar_style()
  area_plot[, fill_col := data.table::fifelse(
    method == "pipe",
    pipe_pal(q50_pu_persistence),
    rank_pal(q50_pu_persistence)
  )]
  area_plot[, segment_border_col := unname(
    pu_bar_style$segment_border_colour[method]
  )]
  bar_outline[, outline_col := unname(
    pu_bar_style$outer_border_colour[method]
  )]

  p_area <- ggplot() +
    geom_rect(
      data = area_plot[pu_area_km2 > 0],
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax,
        fill = fill_col,
        colour = segment_border_col
      ),
      linewidth = pu_bar_style$segment_border_linewidth,
      linejoin = "mitre"
    ) +
    geom_rect(
      data = bar_outline,
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax,
        colour = outline_col
      ),
      fill = NA,
      linewidth = pu_bar_style$outer_border_linewidth,
      linejoin = "mitre"
    ) +
    scale_fill_identity() +
    scale_colour_identity() +
    scale_x_continuous(
      limits = c(0.5, length(stage_show) + 0.5),
      breaks = stage_label_positions,
      labels = stage_show_labels[stage_label_positions],
      expand = expansion(mult = c(0.03, 0.03))
    ) +
    scale_y_continuous(
      breaks = pretty(c(0, y_max), n = 3),
      labels = methods_compact_number_labels,
      expand = expansion(mult = c(0, 0.03))
    ) +
    labs(
      x = "Cells removed (%)",
      y = expression("PU area (km"^2*")"),
      title = "Retained PU area",
      subtitle = if (is.finite(area_target_interval_percent)) {
        paste0(
          focus_species,
          "; nearest completed states every ",
          format(area_target_interval_percent, trim = TRUE),
          " percentage points from ",
          format(panel_min_removed_percent, trim = TRUE),
          "%"
        )
      } else {
        paste0(
          focus_species,
          "; nearest completed states at selected removal levels"
        )
      }
    ) +
    coord_cartesian(ylim = c(0, y_max * 1.05), clip = "off") +
    theme_persist(probability_guides = FALSE)

  area_central <- data.table::copy(focus_central)
  area_bands <- data.table::copy(focus_bands)
  central_origin <- unique(area_central[, .(method, method_lab)])
  central_origin[, `:=`(
    stage = NA_integer_,
    stage_order = Inf,
    total_pu_area_km2 = 0,
    central = 0
  )]
  band_origin <- unique(area_bands[, .(method, method_lab, interval)])
  band_origin[, `:=`(
    stage = NA_integer_,
    stage_order = Inf,
    total_pu_area_km2 = 0,
    ymin = 0,
    ymax = 0
  )]
  area_central_path <- data.table::rbindlist(
    list(area_central, central_origin),
    use.names = TRUE,
    fill = TRUE
  )
  area_band_path <- data.table::rbindlist(
    list(area_bands, band_origin),
    use.names = TRUE,
    fill = TRUE
  )
  data.table::setorder(
    area_central_path,
    method_lab,
    total_pu_area_km2,
    stage_order
  )
  data.table::setorder(
    area_band_path,
    interval,
    method_lab,
    total_pu_area_km2,
    stage_order
  )
  first_window_stage <- available_stage_key[
    x_percent >= area_reference_removed_percent
  ][order(stage_order)][1L]
  assert(
    nrow(first_window_stage) == 1L,
    paste0(
      "Stage 7.3 persistence-area panel has no stage at or above ",
      area_reference_removed_percent,
      "% cell removal."
    )
  )
  window_areas <- area_central[
    stage == first_window_stage$stage,
    total_pu_area_km2
  ]
  area_xmax <- max(window_areas)
  assert(
    length(window_areas) == length(method_labels) &&
      all(is.finite(window_areas) & window_areas >= 0) &&
      is.finite(area_xmax) && area_xmax > 0,
    "Stage 7.3 persistence-area x-axis reference areas are incomplete or invalid."
  )
  # Mark every completed stage rather than a 5%-nearest subset, matching the
  # persistence-through-removal panel above. Point markers are omitted
  # entirely under the "none" symbology (see show_points above).
  area_points <- area_central
  area_x_scale <- scale_x_continuous(
    breaks = pretty(c(0, area_xmax), n = 3),
    labels = methods_compact_number_labels,
    expand = expansion(mult = 0)
  )
  area_coord <- coord_cartesian(xlim = c(0, area_xmax), clip = "on")

  p_area_persist <- ggplot() +
    geom_ribbon(
      data = area_band_path[as.character(interval) == outer_label],
      aes(
        x = total_pu_area_km2,
        ymin = ymin,
        ymax = ymax,
        fill = method_lab,
        alpha = interval,
        group = method_lab
      ),
      color = NA
    ) +
    geom_ribbon(
      data = area_band_path[as.character(interval) == inner_label],
      aes(
        x = total_pu_area_km2,
        ymin = ymin,
        ymax = ymax,
        fill = method_lab,
        alpha = interval,
        group = method_lab
      ),
      color = NA
    ) +
    geom_path(
      data = area_central_path[method == "rank"],
      aes(
        x = total_pu_area_km2,
        y = central,
        color = method_lab,
        group = method_lab
      ),
      linewidth = 1.05,
      lineend = "round"
    ) +
    (if (show_points) {
      do.call(geom_point, c(
        list(
          data = area_points[method == "rank"],
          mapping = aes(x = total_pu_area_km2, y = central, color = method_lab),
          size = point_size,
          alpha = 0.95
        ),
        point_params
      ))
    } else {
      NULL
    }) +
    geom_path(
      data = area_central_path[method == "pipe"],
      aes(
        x = total_pu_area_km2,
        y = central,
        color = method_lab,
        group = method_lab
      ),
      linewidth = 1.05,
      lineend = "round"
    ) +
    (if (show_points) {
      do.call(geom_point, c(
        list(
          data = area_points[method == "pipe"],
          mapping = aes(x = total_pu_area_km2, y = central, color = method_lab),
          size = point_size,
          alpha = 0.95
        ),
        point_params
      ))
    } else {
      NULL
    }) +
    labs(
      x = expression("Retained PU area (km"^2*")"),
      y = "Species persistence",
      title = "Persistence-area relationship",
      subtitle = focus_species
    ) +
    theme_persist()
  p_area_persist <- add_line_scales(
    p_area_persist,
    area_x_scale
  ) + theme(legend.position = "none")
  if (!is.null(area_coord)) {
    p_area_persist <- p_area_persist + area_coord
  }

  list(
    persistence = p_focus,
    retained_area = p_area,
    persistence_area = p_area_persist
  )
}

# Focal-page panel styling ---------------------------------------------

stage73_focal_column_header <- function(title) {
  palette <- methods_figure_palette()
  cowplot::ggdraw() +
    cowplot::draw_label(
      title,
      x = 0.5,
      y = 0.5,
      hjust = 0.5,
      vjust = 0.5,
      fontface = "bold",
      size = 10.1,
      colour = palette[["text"]]
    )
}

stage73_focal_row_strip <- function(scientific_name) {
  palette <- methods_figure_palette()
  cowplot::ggdraw() +
    cowplot::draw_label(
      scientific_name,
      x = 0.012,
      y = 0.58,
      hjust = 0,
      vjust = 0.5,
      fontface = "bold.italic",
      size = 10.0,
      colour = palette[["text"]]
    )
}

stage73_normalize_focal_panel <- function(plot, tag, uncertainty = FALSE) {
  palette <- methods_figure_palette()
  plot <- plot +
    labs(title = NULL, subtitle = NULL, tag = tag) +
    theme(
      legend.position = "none",
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.tag = element_text(
        face = "bold", size = 9.4, colour = palette[["text"]],
        hjust = 0, vjust = 1.35,
        margin = margin(r = 3)
      ),
      plot.tag.position = "topleft",
      plot.margin = margin(8, 5, 2, 6),
      axis.title = element_text(face = "plain", size = 9.2),
      axis.text = element_text(size = 8.3)
    )
  if (isTRUE(uncertainty)) {
    styles <- stage73_uncertainty_styles()
    plot$scales$scales <- Filter(
      function(scale_i) !"alpha" %in% scale_i$aesthetics,
      plot$scales$scales
    )
    plot <- plot + scale_alpha_manual(
      values = stats::setNames(c(0.08, 0.18), styles$interval_label),
      guide = "none"
    )
  }
  plot
}

stage73_spread_focal_area_bars <- function(plot) {
  assert(
    length(plot$layers) >= 2L,
    "Stage 7.3 retained-area panel is missing its bar layers."
  )
  stage_levels <- sort(unique(plot$layers[[1L]]$data$stage_order))
  reposition <- function(rows, half_width) {
    rows <- data.table::as.data.table(data.table::copy(rows))
    rows[, stage_index := match(get("stage_order"), stage_levels)]
    rows[, x_plot := stage_index + data.table::fifelse(method == "pipe", -0.23, 0.23)]
    rows[, `:=`(
      xmin = x_plot - half_width,
      xmax = x_plot + half_width
    )]
    rows[]
  }
  plot$layers[[1L]]$data <- reposition(plot$layers[[1L]]$data, 0.19)
  plot$layers[[2L]]$data <- reposition(plot$layers[[2L]]$data, 0.19)
  plot
}

prepare_stage73_focal_components <- function(
  core,
  focus_species,
  figure_data = NULL,
  point_style = "hollow"
) {
  point_style <- validate_scalar_choice(
    point_style, report_point_styles(), "Stage 7.3 point_style"
  )
  if (is.null(figure_data)) {
    figure_data <- prepare_stage73_figure_data(core)
  }
  focus_species <- trimws(as.character(focus_species))
  assert(
    length(focus_species) == 2L &&
      !anyNA(focus_species) &&
      all(nzchar(focus_species)) &&
      !anyDuplicated(focus_species),
    "The Stage 7.3 combined figure requires exactly two unique focal species."
  )
  focal_spec <- stage73_focal_spec()
  minimum_removed_percent <- focal_spec$minimum_removed_percent
  targets <- suppressWarnings(as.numeric(focal_spec$target_percentages))
  assert(
    length(targets) > 0L &&
      all(is.finite(targets) & targets >= 0 & targets <= 100) &&
      all(diff(targets) > 0),
    "Stage 7.3 focal-panel targets must be increasing percentages in [0,100]."
  )
  selected_stages <- select_stage73_panel_c_stages(
    core$stage_sum,
    target_percentages = targets
  )
  target_intervals <- if (length(targets) > 1L) {
    unique(diff(targets))
  } else {
    100
  }
  target_interval <- if (length(target_intervals) == 1L) {
    target_intervals
  } else {
    NA_real_
  }
  component_sets <- lapply(focus_species, function(species_name) {
    prepare_stage73_focal_panels(
      core = core,
      focus_species = species_name,
      figure_data = figure_data,
      area_stage_ids = selected_stages$stage,
      minimum_removed_percent = minimum_removed_percent,
      persistence_area_reference_removed_percent =
        focal_spec$persistence_area_reference_removed_percent,
      area_target_interval_percent = target_interval,
      point_style = point_style
    )
  })
  names(component_sets) <- focus_species
  list(
    focus_species = focus_species,
    minimum_removed_percent = as.numeric(minimum_removed_percent),
    target_percentages = targets,
    selected_stages = selected_stages,
    component_sets = component_sets,
    point_style = point_style
  )
}

# Across-species distribution panel used by the assemblage figure ----------

stage73_focal_panel_names <- function() {
  c("persistence", "retained_area", "persistence_area")
}

prepare_stage73_focal_figure_components <- function(
  core,
  focus_species,
  figure_data = NULL,
  focal_components = NULL,
  point_style = "hollow"
) {
  point_style <- validate_scalar_choice(
    point_style, report_point_styles(), "Stage 7.3 point_style"
  )
  if (is.null(figure_data)) {
    figure_data <- prepare_stage73_figure_data(core)
  }
  if (is.null(focal_components)) {
    focal_components <- prepare_stage73_focal_components(
      core,
      focus_species,
      figure_data,
      point_style = point_style
    )
  }
  assert(
    length(focal_components$focus_species) == 2L &&
      length(focal_components$component_sets) == 2L,
    "Stage 7.3 focal-page components are incomplete."
  )
  list(
    focal = focal_components,
    rank_label = validate_scalar_string(
      core$rank_label,
      "Stage 7.3 benchmark label"
    ),
    point_style = point_style
  )
}

stage73_focal_panel_key <- function(focus_species, tag_offset = 0L) {
  tag_offset <- suppressWarnings(as.integer(tag_offset))
  panels <- stage73_focal_panel_names()
  n_tags <- length(focus_species) * length(panels)
  assert(
    length(tag_offset) == 1L &&
      !is.na(tag_offset) &&
      tag_offset >= 0L &&
      tag_offset + n_tags <= length(letters),
    "Stage 7.3 focal-page tag offset is invalid."
  )
  key <- data.table::CJ(
    row = seq_along(focus_species),
    column = seq_along(panels),
    sorted = TRUE
  )
  key[, `:=`(
    scientificName = focus_species[row],
    panel = panels[column]
  )]
  key[, tag := paste0(
    "(",
    letters[tag_offset + seq_len(.N)],
    ")"
  )]
  key[]
}

stage73_focal_grid <- function(prepared, panel_key) {
  panels <- stage73_focal_panel_names()
  column_widths <- c(0.94, 1.12, 0.94)
  titles <- c(
    persistence = "Persistence through removal",
    retained_area = "Retained PU area",
    persistence_area = "Persistence–area relationship"
  )
  column_headers <- cowplot::plot_grid(
    plotlist = lapply(titles[panels], stage73_focal_column_header),
    nrow = 1L,
    rel_widths = column_widths
  )

  species_rows <- lapply(seq_along(prepared$focus_species), function(i) {
    components <- prepared$component_sets[[i]]
    plots <- lapply(panels, function(panel_i) {
      tag <- panel_key[
        row == i & panel == panel_i,
        tag
      ][[1L]]
      plot <- switch(
        panel_i,
        persistence = components$persistence,
        retained_area = stage73_spread_focal_area_bars(
          components$retained_area
        ),
        persistence_area = components$persistence_area
      )
      normalized <- stage73_normalize_focal_panel(
        plot,
        tag,
        uncertainty = panel_i %in% c("persistence", "persistence_area")
      )
      # Suppress repeated labels without changing their grob dimensions; this
      # keeps all six data panels at the approved relative size.
      if (i == 1L) {
        normalized <- normalized + ggplot2::theme(
          axis.title.x = ggplot2::element_text(colour = "transparent")
        )
      }
      normalized
    })
    panels_row <- cowplot::plot_grid(
      plotlist = plots,
      nrow = 1L,
      align = "hv",
      axis = "tblr",
      rel_widths = column_widths
    )
    cowplot::plot_grid(
      stage73_focal_row_strip(prepared$focus_species[[i]]),
      panels_row,
      ncol = 1L,
      rel_heights = c(0.060, 1)
    )
  })

  cowplot::plot_grid(
    plotlist = c(list(column_headers), species_rows),
    ncol = 1L,
    rel_heights = c(0.060, rep(0.320, length(species_rows)))
  )
}

stage73_focal_persistence_legend <- function(rank_label, point_style = "hollow") {
  point_style <- validate_scalar_choice(
    point_style, report_point_styles(), "Stage 7.3 point_style"
  )
  point_params <- stage73_point_style_params(point_style)
  point_size <- stage73_point_size(point_style, hollow_size = 0.93, solid_size = 1.10)
  # This legend keys the individual-species persistence/persistence-area
  # panels (b, d, e, g); it drops its point marker only when those panels
  # have no points at all ("none").
  show_points <- !identical(point_style, "none")
  palette <- methods_figure_palette()
  rank_display <- sub("^Zonation\\s+", "", rank_label)
  ggplot2::ggplot() +
    ggplot2::annotate(
      "text", x = 0.01, y = 0.94, label = "Species persistence (q50)",
      hjust = 0, fontface = "bold", size = 2.96, colour = palette[["text"]]
    ) +
    # Method keys (row 1)
    ggplot2::annotate(
      "segment", x = 0.02, xend = 0.08, y = 0.62, yend = 0.62,
      linewidth = 0.76, colour = palette[["main"]], lineend = "round"
    ) +
    (if (show_points) {
      do.call(ggplot2::annotate, c(
        list("point", x = 0.05, y = 0.62, size = point_size, colour = palette[["main"]]),
        point_params
      ))
    } else {
      NULL
    }) +
    ggplot2::annotate(
      "text", x = 0.095, y = 0.62, label = "PersistRank",
      hjust = 0, size = 2.54, colour = palette[["text"]]
    ) +
    ggplot2::annotate(
      "segment", x = 0.56, xend = 0.62, y = 0.62, yend = 0.62,
      linewidth = 0.76, colour = palette[["benchmark"]], lineend = "round"
    ) +
    (if (show_points) {
      do.call(ggplot2::annotate, c(
        list("point", x = 0.59, y = 0.62, size = point_size, colour = palette[["benchmark"]]),
        point_params
      ))
    } else {
      NULL
    }) +
    ggplot2::annotate(
      "text", x = 0.635, y = 0.62, label = rank_display,
      hjust = 0, size = 2.54, colour = palette[["text"]]
    ) +
    # Interval keys (row 2)
    ggplot2::annotate(
      "rect", xmin = 0.02, xmax = 0.05, ymin = 0.10, ymax = 0.32,
      fill = palette[["main"]], alpha = 0.25, colour = NA
    ) +
    ggplot2::annotate(
      "rect", xmin = 0.05, xmax = 0.08, ymin = 0.10, ymax = 0.32,
      fill = palette[["benchmark"]], alpha = 0.25, colour = NA
    ) +
    ggplot2::annotate(
      "text", x = 0.095, y = 0.21, label = "68% interval",
      hjust = 0, size = 2.54, colour = palette[["text"]]
    ) +
    ggplot2::annotate(
      "rect", xmin = 0.56, xmax = 0.59, ymin = 0.10, ymax = 0.32,
      fill = palette[["main"]], alpha = 0.10, colour = NA
    ) +
    ggplot2::annotate(
      "rect", xmin = 0.59, xmax = 0.62, ymin = 0.10, ymax = 0.32,
      fill = palette[["benchmark"]], alpha = 0.10, colour = NA
    ) +
    ggplot2::annotate(
      "text", x = 0.635, y = 0.21, label = "95% interval",
      hjust = 0, size = 2.54, colour = palette[["text"]]
    ) +
    ggplot2::coord_cartesian(xlim = c(0, 1.08), ylim = c(0.02, 1.02), clip = "off") +
    ggplot2::theme_void()
}

stage73_focal_pu_legend <- function(rank_label) {
  palette <- methods_figure_palette()
  rank_display <- sub("^Zonation\\s+", "", rank_label)
  gradient_x <- seq(0, 1, length.out = 1024L)
  pipe_cols <- c("#D9E8F2", "#C8DDEB", "#AFCDE1", "#82B0D1", "#5794C0", "#2E78AF")
  rank_cols <- c("#F6DDD0", "#F2CDB9", "#ECB397", "#E38E62", "#D36E39", "#B7561D")
  pipe_raster <- as.raster(matrix(
    scales::col_numeric(pipe_cols, domain = c(0, 1))(gradient_x), nrow = 1L
  ))
  rank_raster <- as.raster(matrix(
    scales::col_numeric(rank_cols, domain = c(0, 1))(gradient_x), nrow = 1L
  ))
  ggplot2::ggplot() +
    ggplot2::annotation_raster(pipe_raster, xmin = 0.30, xmax = 0.96, ymin = 0.61, ymax = 0.75) +
    ggplot2::annotation_raster(rank_raster, xmin = 0.30, xmax = 0.96, ymin = 0.34, ymax = 0.48) +
    ggplot2::annotate(
      "rect", xmin = 0.30, xmax = 0.96, ymin = 0.61, ymax = 0.75,
      fill = NA, colour = palette[["main"]], linewidth = 0.42
    ) +
    ggplot2::annotate(
      "rect", xmin = 0.30, xmax = 0.96, ymin = 0.34, ymax = 0.48,
      fill = NA, colour = palette[["benchmark"]], linewidth = 0.42
    ) +
    ggplot2::annotate(
      "text", x = 0.01, y = 0.96, label = "PU persistence (q50)",
      hjust = 0, fontface = "bold", size = 2.94, colour = palette[["text"]]
    ) +
    ggplot2::annotate(
      "text", x = 0.27, y = 0.68, label = "PersistRank",
      hjust = 1, size = 2.47, colour = palette[["text"]]
    ) +
    ggplot2::annotate(
      "text", x = 0.27, y = 0.41, label = rank_display,
      hjust = 1, size = 2.47, colour = palette[["text"]]
    ) +
    ggplot2::annotate(
      "text", x = c(0.30, 0.63, 0.96), y = 0.16,
      label = c("0", "0.5", "1"),
      size = 2.30, colour = palette[["muted_text"]]
    ) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0.08, 1), clip = "off") +
    ggplot2::theme_void()
}

stage73_focal_legend <- function(rank_label, point_style = "hollow") {
  cowplot::plot_grid(
    stage73_focal_persistence_legend(rank_label, point_style = point_style),
    stage73_focal_pu_legend(rank_label),
    nrow = 1L,
    rel_widths = c(0.60, 0.40)
  )
}

stage73_focal_section <- function(prepared_components, point_style = "hollow") {
  point_style <- validate_scalar_choice(
    point_style, report_point_styles(), "Stage 7.3 point_style"
  )
  assert(
    is.list(prepared_components) &&
      length(prepared_components$focal$focus_species) == 2L,
    "Prepared Stage 7.3 focal figure components are invalid."
  )
  panel_key <- stage73_focal_panel_key(
    prepared_components$focal$focus_species,
    tag_offset = 1L
  )
  focal_grid <- stage73_focal_grid(
    prepared_components$focal,
    panel_key
  )
  legend <- stage73_focal_legend(prepared_components$rank_label, point_style = point_style)
  out <- cowplot::plot_grid(
    focal_grid,
    legend,
    ncol = 1L,
    rel_heights = c(0.455, 0.090)
  )
  list(
    plot = out,
    panel_key = panel_key,
    selected_stages = prepared_components$focal$selected_stages,
    minimum_removed_percent = prepared_components$focal$minimum_removed_percent,
    target_percentages = prepared_components$focal$target_percentages
  )
}

# Assemblage figure: across-species persistence distribution ---------------
