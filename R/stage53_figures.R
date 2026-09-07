# Stage 5.3 presentation and atomic SI-figure persistence.
#
# stage53_workflow.R loads this module after Stage 5.3 analysis and the shared
# figure/transaction contracts. Plot construction is transient; only the writer
# changes disk state. Sourcing attaches no package and builds no figure.

stage53_taxon_colors <- function(species) {
  palette <- methods_figure_palette()
  colors <- c(AVES = palette[["main"]], MAMMALIA = palette[["sensitivity"]])
  colors[names(colors) %in% unique(as.character(species$className))]
}

validate_stage53_plot_species <- function(species) {
  need_cols(
    species,
    c("className", "sdm_method", "total_pu_area_km2", "species_persistence"),
    "Stage 5.3 species plot data"
  )
  area <- suppressWarnings(as.numeric(species$total_pu_area_km2))
  persistence <- suppressWarnings(as.numeric(species$species_persistence))
  taxon <- stringr::str_squish(as.character(species$className))
  sdm <- stringr::str_squish(as.character(species$sdm_method))
  assert(
    all(is.finite(area) & area > 0),
    "Stage 5.3 species plot areas must be positive and finite."
  )
  assert(
    all(is.finite(persistence) & persistence >= 0 & persistence <= 1),
    "Stage 5.3 species plot persistence values must be finite and in [0, 1]."
  )
  assert(
    all(taxon %in% c("MAMMALIA", "AVES")),
    "Stage 5.3 species plot data contain an unsupported taxon class."
  )
  assert(
    all(!is.na(sdm) & nzchar(sdm) & sdm %in% c("PPM", "RangeBag")),
    "Stage 5.3 species plot data must use SDM method PPM or RangeBag."
  )
  invisible(TRUE)
}

build_stage53_area_scatterplot <- function(species) {
  validate_stage53_plot_species(species)
  plot_data <- species |>
    dplyr::mutate(
      sdm_method = factor(
        stringr::str_squish(as.character(sdm_method)),
        levels = c("PPM", "RangeBag")
      )
    )
  colors <- stage53_taxon_colors(plot_data)

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = total_pu_area_km2,
      y = species_persistence,
      colour = className,
      shape = sdm_method
    )
  ) +
    ggplot2::geom_point(size = 2.0, alpha = 0.78, stroke = 0.4) +
    ggplot2::scale_colour_manual(
      values = colors,
      labels = c(AVES = "Birds", MAMMALIA = "Mammals"),
      name = "Taxon"
    ) +
    ggplot2::scale_shape_manual(
      values = c(PPM = 16, RangeBag = 17),
      breaks = c("PPM", "RangeBag"),
      name = "SDM branch"
    ) +
    ggplot2::scale_x_log10(labels = scales::label_comma()) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      breaks = methods_probability_breaks(),
      labels = methods_probability_labels()
    ) +
    ggplot2::labs(
      title = "Persistence and retained area",
      x = expression("Total retained PU area (km"^2*")"),
      y = "Initial species persistence"
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(order = 1, nrow = 1),
      shape = ggplot2::guide_legend(order = 2, nrow = 1)
    ) +
    theme_methods_figure(base_size = 12, legend_position = "bottom") +
    theme_methods_probability_guides() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 11.5),
      plot.title.position = "panel",
      axis.text = ggplot2::element_text(size = 9),
      axis.title = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 9),
      legend.title = ggplot2::element_text(size = 9)
    )
}

build_stage53_figure <- function(species) {
  validate_stage53_plot_species(species)
  colors <- stage53_taxon_colors(species)
  p_distribution <- ggplot2::ggplot(
    species,
    ggplot2::aes(x = className, y = species_persistence, colour = className)
  ) +
    ggplot2::geom_boxplot(width = 0.52, outlier.shape = NA, linewidth = 0.5) +
    ggplot2::geom_jitter(width = 0.14, height = 0, size = 1.45, alpha = 0.72) +
    ggplot2::scale_colour_manual(values = colors, guide = "none") +
    ggplot2::scale_x_discrete(
      labels = c(AVES = "Birds", MAMMALIA = "Mammals")
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      breaks = methods_probability_breaks(),
      labels = methods_probability_labels()
    ) +
    ggplot2::labs(
      title = "Persistence by taxon",
      x = NULL,
      y = "Initial species persistence"
    ) +
    theme_methods_figure(base_size = 12, legend_position = "none") +
    theme_methods_probability_guides() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 11.5),
      plot.title.position = "panel",
      axis.text = ggplot2::element_text(size = 9),
      axis.title = ggplot2::element_text(size = 10)
    )

  p_area <- build_stage53_area_scatterplot(species)
  shared_legend <- cowplot::get_legend(
    p_area +
      ggplot2::theme(
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.box.just = "center",
        legend.spacing.x = grid::unit(0.45, "cm")
      )
  )
  panels <- panel_tag_grid(
    list(p_distribution, p_area + ggplot2::theme(legend.position = "none")),
    labels = c("(a)", "(b)"),
    ncol = 2,
    label_size = 11
  )
  cowplot::plot_grid(
    panels,
    shared_legend,
    ncol = 1,
    rel_heights = c(1, 0.12)
  )
}

write_stage53_figure_atomic <- function(figure, config) {
  target <- config$paths$figure
  ensure_writable_dir(dirname(target), "Stage 5.3 SI figure directory")
  staging_root <- tempfile("stage53_figure_", tmpdir = dirname(target))
  assert(dir.create(staging_root), "Could not create Stage 5.3 figure staging directory.")
  on.exit(unlink(staging_root, recursive = TRUE, force = TRUE), add = TRUE)

  staged <- save_manuscript_figure(
    figure,
    "initial_persistence",
    figure_dir = staging_root
  )
  project_file_set_transaction(
    staged,
    target,
    overwrite = TRUE,
    label = "Stage 5.3 SI figure"
  )
  invisible(target)
}
