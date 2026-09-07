# Stage 4 area-curve figure construction.
#
# Figure mode reads only the validated canonical species table. It extracts the
# internally fixed focal species, evaluates its q50 Gompertz curve and scaled
# marginal sensitivity over habitat area (km^2), and writes the finalized
# main-text figure through the shared rollback-capable file transaction. It
# never reads raw trait data, Stage 3 models, SDM rasters, or the IUCN API.

theme_species_area_curve <- function(base_size = 11.6) {
  theme_methods_figure(base_size = base_size, legend_position = "inside") +
    ggplot2::theme(
      plot.margin = ggplot2::margin(2, 6, 2, 5),
      panel.grid = ggplot2::element_blank(),
      legend.title = ggplot2::element_blank(),
      legend.key.width = grid::unit(1.35, "lines"),
      legend.position.inside = c(0.125, 0.975),
      legend.justification = c(0, 1),
      legend.direction = "vertical",
      legend.background = ggplot2::element_rect(
        fill = scales::alpha("white", 0.90),
        colour = NA
      ),
      legend.margin = ggplot2::margin(1, 3, 1, 3)
    )
}

build_area_curve_figure <- function(species_table,
                                    focal_species,
                                    curve = main_persistence_curve()) {
  curve <- validate_persistence_curve(curve)
  alpha_col <- paste0("alpha_", curve)
  beta_col <- paste0("beta_", curve)
  need_cols(
    species_table,
    c("scientificName", "density", "min_pop_size", alpha_col, beta_col),
    "species_table for area-curve figure"
  )

  focal_species <- as.character(focal_species)
  sp_name <- if (length(focal_species)) stringr::str_squish(focal_species[[1L]]) else ""
  assert(!is.na(sp_name) && nzchar(sp_name), "focal_species must identify one species.")

  species_names <- stringr::str_squish(as.character(species_table$scientificName))
  sp_rows <- which(species_names == sp_name)

  if (!length(sp_rows)) {
    stop(paste0("Species not found in species_table: ", sp_name), call. = FALSE)
  }
  if (length(sp_rows) > 1L) {
    detail_cols <- intersect(c("scientificName", "sdm_method", "raster_file"), names(species_table))
    details <- apply(
      as.data.frame(species_table[sp_rows, detail_cols, drop = FALSE]),
      1,
      function(row) paste(names(row), row, sep = "=", collapse = " | ")
    )
    stop(
      paste0(
        "Species matches multiple rows in species_table: ", sp_name, ".\n",
        "Use a table with one row for the focal species before building the area-curve figure.\n",
        paste(details, collapse = "\n")
      ),
      call. = FALSE
    )
  }

  sp_row <- species_table[sp_rows, , drop = FALSE]

  dens <- as.numeric(sp_row$density)
  alpha <- as.numeric(sp_row[[alpha_col]])
  beta <- as.numeric(sp_row[[beta_col]])
  A0 <- as.numeric(sp_row$min_pop_size)

  validate_species_table_curve_parameters(
    sp_row,
    curves = curve,
    label = paste0("species_table area-curve row for ", sp_name),
    rows_label = "the focal species"
  )

  if (!is.finite(dens) || dens <= 0 || !is.finite(A0) || A0 <= 0) {
    stop(
      paste0(
        "The required area-curve inputs are missing or invalid for ", sp_name, ".\n",
        "Check density and min_pop_size in species_table.csv."
      ),
      call. = FALSE
    )
  }

  c_eff <- alpha * dens^(-beta)
  delta_star <- c_eff^(1 / beta)
  delta_grid <- 10^seq(log10(delta_star / 80), log10(delta_star * 350), length.out = 1500)

  A <- A0 + delta_grid
  log_z <- log(c_eff) - beta * log(delta_grid)
  z <- exp(log_z)
  P <- exp(-z)

  log_dP <- (-z) + log(c_eff) + log(beta) - (beta + 1) * log(delta_grid)
  dP <- exp(log_dP)

  df <- tibble::tibble(
    area_km2 = c(A0, A),
    P = c(0, P),
    dP = c(0, dP)
  ) |>
    dplyr::arrange(.data$area_km2)

  s_deriv <- max(df$dP[is.finite(df$dP)], na.rm = TRUE)
  assert(is.finite(s_deriv) && s_deriv > 0, "Could not compute a finite marginal-sensitivity curve.")
  df <- df |>
    dplyr::mutate(dP_scaled = .data$dP / s_deriv)

  band_cut <- 0.35
  band_df <- df |>
    dplyr::filter(.data$dP_scaled >= band_cut)
  assert(nrow(band_df) > 0L, "Could not locate the high-sensitivity area band.")

  x_band_lo <- min(band_df$area_km2, na.rm = TRUE)
  x_band_hi <- max(band_df$area_km2, na.rm = TRUE)

  saturated_area <- df$area_km2[df$P >= 0.995 & is.finite(df$area_km2)]
  x_sat <- if (length(saturated_area)) {
    min(saturated_area)
  } else {
    max(df$area_km2, na.rm = TRUE)
  }

  x_min <- max(A0 * 0.85, x_band_lo / 1.35)
  x_max <- max(x_sat * 1.18, x_band_hi * 1.12)
  # Leave enough log-scale separation for the threshold line to remain visible
  # after the figure is reduced.
  x_min <- min(x_min, A0 * 0.70)
  x_max <- max(x_max, A0 * 1.04)

  plot_df <- df |>
    dplyr::filter(.data$area_km2 >= x_min, .data$area_km2 <= x_max)

  x_breaks <- scales::breaks_log(n = 5)(c(x_min, x_max))
  x_breaks <- x_breaks[x_breaks >= x_min & x_breaks <= x_max]
  if (!length(x_breaks)) x_breaks <- scales::breaks_log(n = 4)(c(x_min, x_max))

  pal <- methods_figure_palette()
  shade <- ggplot2::annotate(
    "rect",
    xmin = x_band_lo,
    xmax = x_band_hi,
    ymin = -Inf,
    ymax = Inf,
    fill = pal[["main_fill"]],
    alpha = 0.16
  )
  threshold_line <- ggplot2::geom_vline(
    xintercept = A0,
    linewidth = 0.52,
    linetype = "dotted",
    colour = pal[["threshold"]]
  )

  persistence_label <- "Persistence"
  sensitivity_label <- "Sensitivity (scaled)"
  series_levels <- c(persistence_label, sensitivity_label)
  line_df <- dplyr::bind_rows(
    plot_df |>
      dplyr::transmute(
        area_km2 = .data$area_km2,
        response = .data$P,
        series = persistence_label
      ),
    plot_df |>
      dplyr::transmute(
        area_km2 = .data$area_km2,
        response = .data$dP_scaled,
        series = sensitivity_label
      )
  ) |>
    dplyr::mutate(series = factor(.data$series, levels = series_levels))

  ggplot2::ggplot(
    line_df,
    ggplot2::aes(
      x = .data$area_km2,
      y = .data$response,
      colour = .data$series,
      linetype = .data$series
    )
  ) +
    shade +
    threshold_line +
    ggplot2::geom_line(
      linewidth = 1.12,
      lineend = "round",
      na.rm = TRUE
    ) +
    ggplot2::scale_colour_manual(
      values = stats::setNames(
        c(pal[["main"]], "#5E6B71"),
        series_levels
      ),
      breaks = series_levels
    ) +
    ggplot2::scale_linetype_manual(
      values = stats::setNames(c("solid", "22"), series_levels),
      breaks = series_levels
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(ncol = 1, byrow = TRUE),
      linetype = ggplot2::guide_legend(ncol = 1, byrow = TRUE)
    ) +
    ggplot2::scale_x_log10(
      limits = c(x_min, x_max),
      breaks = x_breaks,
      labels = methods_compact_number_labels,
      expand = ggplot2::expansion(mult = c(0.01, 0.03))
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      breaks = methods_probability_breaks(),
      labels = methods_probability_labels(),
      expand = ggplot2::expansion(mult = c(0, 0.01))
    ) +
    ggplot2::labs(
      x = expression("Habitat area (km"^2*")"),
      y = "Persistence probability"
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    theme_species_area_curve()
}

rebuild_area_curve_figure <- function(
    species_csv = file.path("Data", "Clean", "species_table.csv"),
    focal_species = "Fossa fossana",
    curve = main_persistence_curve(),
    figure_dir = "Figures",
    overwrite = FALSE
) {
  need_file(species_csv, "Stage 4 species table for area-curve figure")
  species_table <- readr::read_csv(
    species_csv,
    show_col_types = FALSE,
    progress = FALSE
  )
  figure <- build_area_curve_figure(
    species_table = species_table,
    focal_species = focal_species,
    curve = curve
  )
  write_area_curve_figure(figure, figure_dir = figure_dir, overwrite = overwrite)
  figure
}

install_area_curve_figure <- function(figure, figure_dir, overwrite = FALSE,
                                      rename_file = file.rename) {
  path <- file.path(figure_dir, "area_curve.png")
  assert(
    !file.exists(path) || isTRUE(overwrite),
    paste0("Refusing to overwrite Stage 4 area-curve figure: ", path)
  )
  ensure_writable_dir(figure_dir, "Stage 4 figure directory")
  staging_dir <- tempfile("stage4_area_curve_", tmpdir = figure_dir)
  assert(dir.create(staging_dir), "Could not create a temporary Stage 4 figure directory.")
  on.exit(unlink(staging_dir, recursive = TRUE, force = TRUE), add = TRUE)

  save_manuscript_figure(figure, "area_curve", figure_dir = staging_dir)
  staged <- file.path(staging_dir, "area_curve.png")
  assert(file.exists(staged) && is.finite(file.info(staged)$size) && file.info(staged)$size > 0,
         "Temporary Stage 4 area-curve figure was not written correctly.")

  project_file_set_transaction(
    staged_paths = staged,
    target_paths = path,
    overwrite = overwrite,
    rename_file = rename_file,
    label = "Stage 4 area-curve figure"
  )
  path
}

write_area_curve_figure <- function(figure, figure_dir = "Figures", overwrite = FALSE) {
  path <- install_area_curve_figure(figure, figure_dir, overwrite = overwrite)
  log_msg("Stage 4 | area-curve figure written | path=", path)
  invisible(path)
}
