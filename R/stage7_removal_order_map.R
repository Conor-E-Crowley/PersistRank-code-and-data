# Stage 7.1 removal surface and main-text map.
#
# The validator proves that Stage 6 events, integer removal steps, continuous
# ranks, and the bundle's initial occupied domain describe the same cells. The
# plot builder consumes only the validated rank map and does not alter its
# values, extent, or scientific interpretation.


stage71_removal_map_packages <- function() {
  c("data.table", "terra", "sf", "rnaturalearth", "ggplot2", "tidyterra", "scales")
}

check_stage71_removal_map_packages <- function() {
  required <- stage71_removal_map_packages()
  missing <- required[!vapply(required, requireNamespace, logical(1L), quietly = TRUE)]
  assert(
    !length(missing),
    paste0(
      "Stage 7.1 removal-order map requires missing package(s): ",
      paste(missing, collapse = ", "),
      "."
    )
  )
  invisible(TRUE)
}

validate_stage71_removal_surface <- function(
  rankmap,
  removal_order,
  events,
  alive0,
  tolerance = 1e-6
) {
  assert(
    inherits(rankmap, "SpatRaster") && terra::nlyr(rankmap) == 1L,
    "Stage 7.1 rank map must be a single-layer SpatRaster."
  )
  assert(
    inherits(removal_order, "SpatRaster") && terra::nlyr(removal_order) == 1L,
    "Stage 7.1 removal-order map must be a single-layer SpatRaster."
  )
  assert(
    is.numeric(tolerance) && length(tolerance) == 1L &&
      is.finite(tolerance) && tolerance > 0,
    "Stage 7.1 validation tolerance must be a positive finite scalar."
  )
  assert(
    isTRUE(terra::compareGeom(rankmap, removal_order, stopOnError = FALSE)),
    "Stage 7.1 rank map geometry differs from removal_order.tif."
  )

  events <- data.table::as.data.table(events)
  need_cols(
    events,
    c(
      "removal_step", "cells_removed", "cum_cells_removed",
      "cum_prop_cells_removed"
    ),
    "validated Stage 7.1 removal events"
  )
  assert(
    identical(as.integer(events$removal_step), seq_len(nrow(events))),
    "Stage 7.1 removal steps must be contiguous positive integers."
  )
  assert(
    length(alive0) == terra::ncell(removal_order) &&
      is.logical(alive0) &&
      !anyNA(alive0),
    "Stage 7.1 initial occupied-cell state is incompatible with the removal raster."
  )
  mask_mismatch <- as.numeric(terra::global(
    terra::ifel(is.na(rankmap) != is.na(removal_order), 1, 0),
    "sum",
    na.rm = TRUE
  )[1, 1])
  assert(
    identical(mask_mismatch, 0),
    "Stage 7.1 rank map and removal-order map have different ranked-cell domains."
  )

  alive_raster <- terra::rast(removal_order)
  terra::values(alive_raster) <- as.integer(alive0)
  domain_mismatch <- as.numeric(terra::global(
    terra::ifel(
      (alive_raster == 1L) != !is.na(removal_order),
      1,
      0
    ),
    "sum",
    na.rm = TRUE
  )[1, 1])
  assert(
    identical(domain_mismatch, 0),
    "Stage 7.1 ranked-cell domain differs from the bundle's initial occupied domain."
  )

  n_ranked <- as.numeric(terra::global(
    terra::ifel(!is.na(rankmap), 1, 0),
    "sum",
    na.rm = TRUE
  )[1, 1])
  assert(
    is.finite(n_ranked) && n_ranked > 0,
    "Stage 7.1 rank map contains no ranked cells."
  )
  assert(
    identical(n_ranked, as.numeric(sum(events$cells_removed))),
    "Stage 7.1 ranked-cell count does not match removal_events.csv."
  )

  rank_range <- terra::global(rankmap, c("min", "max"), na.rm = TRUE)
  min_rank <- as.numeric(rank_range[1, "min"])
  max_rank <- as.numeric(rank_range[1, "max"])
  assert(
    all(is.finite(c(min_rank, max_rank))) &&
      min_rank >= -tolerance && max_rank <= 1 + tolerance,
    "Stage 7.1 rank map values must be finite and in [0, 1]."
  )

  order_range <- terra::global(removal_order, c("min", "max"), na.rm = TRUE)
  min_order <- as.numeric(order_range[1, "min"])
  max_order <- as.numeric(order_range[1, "max"])
  order_fraction <- as.numeric(terra::global(
    abs(removal_order - round(removal_order)),
    "max",
    na.rm = TRUE
  )[1, 1])
  assert(
    all(is.finite(c(min_order, max_order, order_fraction))) &&
      min_order > 0 && order_fraction <= tolerance,
    "Stage 7.1 removal-order map must contain positive integer steps."
  )

  observed_counts <- data.table::as.data.table(terra::freq(removal_order))
  data.table::setnames(observed_counts, c("layer", "value", "count"))
  observed_counts <- observed_counts[, .(
    removal_step = as.integer(value),
    observed_cells = as.integer(count)
  )]
  expected_counts <- events[cells_removed > 0L, .(
    removal_step,
    expected_cells = cells_removed
  )]
  count_check <- merge(
    expected_counts,
    observed_counts,
    by = "removal_step",
    all = TRUE,
    sort = TRUE
  )
  assert(
    nrow(count_check) == nrow(expected_counts) &&
      !anyNA(count_check) &&
      identical(count_check$observed_cells, count_check$expected_cells),
    "Stage 7.1 removal-order raster frequencies do not match removal events."
  )

  zonal_ranks <- data.table::as.data.table(terra::zonal(
    rankmap,
    removal_order,
    fun = function(x, ...) c(min = min(x, ...), max = max(x, ...)),
    na.rm = TRUE
  ))
  data.table::setnames(
    zonal_ranks,
    seq_len(ncol(zonal_ranks)),
    c("removal_step", "rank_min", "rank_max")
  )
  zonal_ranks[, removal_step := as.integer(removal_step)]
  rank_check <- merge(
    events[cells_removed > 0L, .(
      removal_step,
      expected_rank = cum_prop_cells_removed
    )],
    zonal_ranks,
    by = "removal_step",
    all = TRUE,
    sort = TRUE
  )
  assert(
    !anyNA(rank_check) &&
      nrow(rank_check) == nrow(expected_counts) &&
      all(abs(rank_check$rank_min - rank_check$rank_max) <= tolerance) &&
      all(abs(rank_check$rank_min - rank_check$expected_rank) <= tolerance),
    paste0(
      "Stage 7.1 rank map is inconsistent with cumulative cell-removal ",
      "proportions in removal_events.csv."
    )
  )

  n_terminal <- as.numeric(terra::global(
    terra::ifel(abs(rankmap - 1) <= tolerance, 1, 0),
    "sum",
    na.rm = TRUE
  )[1, 1])
  assert(
    is.finite(n_terminal) && n_terminal > 0 && abs(max_rank - 1) <= tolerance,
    "Stage 7.1 rank map must contain a terminal retained layer at 1."
  )

  names(rankmap) <- "removal_rank"
  list(
    rankmap = rankmap,
    summary = list(
      n_ranked = as.integer(n_ranked),
      n_steps = nrow(events),
      min_rank = min_rank,
      max_rank = max_rank,
      n_terminal = as.integer(n_terminal)
    )
  )
}

stage71_removal_map_palette <- function() {
  c("#CAD9DC", "#ABC6CE", "#7FA9B9", "#568EA6", "#2E78AF", "#1F5269")
}

stage71_madagascar_outline <- function(template_raster) {
  outline <- rnaturalearth::ne_countries(
    country = "Madagascar",
    scale = "medium",
    returnclass = "sf"
  )
  sf::st_transform(outline, crs = terra::crs(template_raster, proj = TRUE))
}

stage71_removal_map_extent <- function(outline, pad_x = 0.020, pad_y = 0.015) {
  bounds <- sf::st_bbox(outline)
  dx <- unname(bounds[["xmax"]] - bounds[["xmin"]])
  dy <- unname(bounds[["ymax"]] - bounds[["ymin"]])
  terra::ext(
    bounds[["xmin"]] - pad_x * dx,
    bounds[["xmax"]] + pad_x * dx,
    bounds[["ymin"]] - pad_y * dy,
    bounds[["ymax"]] + pad_y * dy
  )
}

build_stage71_removal_order_map <- function(
  rankmap,
  outline = stage71_madagascar_outline(rankmap),
  maxcell = 3e6
) {
  assert(
    inherits(rankmap, "SpatRaster") && terra::nlyr(rankmap) == 1L,
    "Stage 7.1 plotting requires a single-layer rank map."
  )
  assert(
    inherits(outline, "sf") && nrow(outline) > 0L,
    "Stage 7.1 plotting requires a non-empty sf outline."
  )
  assert(
    is.numeric(maxcell) && length(maxcell) == 1L &&
      is.finite(maxcell) && maxcell > 0,
    "Stage 7.1 maxcell must be a positive finite scalar."
  )

  palette <- methods_figure_palette()
  map_extent <- stage71_removal_map_extent(outline)
  plot_raster <- terra::crop(rankmap, map_extent, snap = "out")
  names(plot_raster) <- "removal_rank"

  ggplot2::ggplot() +
    ggplot2::geom_sf(
      data = outline,
      fill = palette[["land"]],
      colour = NA,
      inherit.aes = FALSE
    ) +
    suppressMessages(tidyterra::geom_spatraster(
      data = plot_raster,
      ggplot2::aes(fill = removal_rank),
      maxcell = maxcell,
      na.rm = TRUE
    )) +
    ggplot2::scale_fill_gradientn(
      colours = stage71_removal_map_palette(),
      values = seq(0, 1, length.out = length(stage71_removal_map_palette())),
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.25),
      labels = scales::label_percent(accuracy = 1),
      oob = scales::squish,
      na.value = NA,
      name = "Removal\norder",
      guide = ggplot2::guide_colourbar(
        direction = "vertical",
        title.position = "top",
        title.hjust = 1,
        label.position = "left",
        barwidth = grid::unit(3.8, "mm"),
        barheight = grid::unit(38, "mm"),
        ticks.colour = palette[["axis"]],
        frame.colour = palette[["panel_border"]]
      )
    ) +
    ggplot2::geom_sf(
      data = outline,
      fill = NA,
      colour = palette[["axis"]],
      linewidth = 0.42,
      inherit.aes = FALSE
    ) +
    ggplot2::coord_sf(
      crs = terra::crs(plot_raster, proj = TRUE),
      xlim = c(map_extent$xmin, map_extent$xmax),
      ylim = c(map_extent$ymin, map_extent$ymax),
      expand = FALSE
    ) +
    theme_methods_map(base_size = 10.8) +
    ggplot2::theme(
      legend.position = "inside",
      legend.position.inside = c(0.982, 0.035),
      legend.justification = c(1, 0),
      legend.direction = "vertical",
      legend.title = ggplot2::element_text(
        face = "bold",
        colour = palette[["text"]],
        size = 9.6,
        hjust = 1,
        lineheight = 0.92,
        margin = ggplot2::margin(b = 4)
      ),
      legend.text = ggplot2::element_text(
        colour = palette[["text"]],
        size = 8.8,
        hjust = 1
      ),
      legend.background = ggplot2::element_rect(
        fill = scales::alpha("white", 0.88),
        colour = NA
      ),
      legend.margin = ggplot2::margin(2, 3, 2, 3),
      plot.margin = ggplot2::margin(3, 4, 3, 4)
    )
}
