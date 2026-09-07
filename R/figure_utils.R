# Canonical manuscript-figure outputs.
#
# The manifest is the single authority for figure filenames, dimensions, and
# resolution. Shared palettes, probability breaks, and themes keep main-text
# and SI graphics consistent; builders remain in their stage-specific modules.


manuscript_figure_manifest <- function() {
  data.frame(
    id = c(
      "allometry",
      "stage1_coefficients",
      "gompertz_wolff",
      "gompertz_loess_combined",
      "abundance_persistence",
      "area_curve",
      "spatial_process",
      "initial_persistence",
      "stage7_pipeline_persistence",
      "stage7_curve_uncertainty",
      "stage7_species_persistence",
      "persistence_comparison_assemblage_focal",
      "persistence_comparison_assemblage_focal_caz1",
      "persistence_comparison_assemblage_focal_caz2",
      "persistence_comparison_assemblage_focal_cazmax",
      "cell_removal_order"
    ),
    filename = c(
      "demographic_calibration.png",
      "SI/S1_demographic_posteriors.png",
      "gompertz_wolff.png",
      "SI/S3_gompertz_parameters.png",
      "SI/S3_abundance_persistence.png",
      "area_curve.png",
      "spatial_process.png",
      "SI/S5_initial_persistence.png",
      "SI/S7_pipeline_persistence.png",
      "SI/S7_curve_uncertainty.png",
      "SI/S7_species_persistence.png",
      "results_comparison_assemblage_focal.png",
      "SI/S7_persistence_comparison_caz1.png",
      "SI/S7_persistence_comparison_caz2.png",
      "SI/S7_persistence_comparison_cazmax.png",
      "cell_removal_order.png"
    ),
    stage = c(
      "1", "1", "3", "3", "3", "4", "5.1", "5.3",
      "7.1", "7.1", "7.1", "7.3", "7.3", "7.3", "7.3",
      "7.1"
    ),
    content = c(
      "Mammal and bird demographic allometries.",
      "Posterior densities for Stage 1 growth, environmental-variation, and configured bird-diet coefficients.",
      "Shifted-Gompertz and Wolff comparison.",
      "Mammal and bird alpha/beta LOESS relationships for the available bird model groups.",
      "Combined median abundance-persistence relationships for mammals and the configured bird-diet groups.",
      "Focal-species persistence probability and marginal sensitivity in area units.",
      "Habitat-to-population-unit workflow.",
      "Initial species persistence by taxon and retained population-unit area.",
      "Mean persistence with middle-50% and 10th/90th percentile ranges through cell removal.",
      "Mean species persistence under five demographic coefficient sets.",
      "Species-level persistence through cell removal, ordered within taxon.",
      "Full-range assemblage distribution above the late-stage two-species comparison grid.",
      "CAZ1 comparison using the seven-panel assemblage and focal-species layout.",
      "CAZ2 comparison using the seven-panel assemblage and focal-species layout.",
      "CAZMAX comparison using the seven-panel assemblage and focal-species layout.",
      "Continuous persistence-based cell removal order across Madagascar."
    ),
    width = c(
      5.2, 7.0, 6.2, 10.4, 10.4, 4.8, 7.2, 7.0,
      7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0,
      4.8
    ),
    height = c(
      2.85, 4.2, 2.85, 5.4, 3.4, 3.0, 7.2, 3.4,
      3.3, 3.3, 9.0,
      7.45, 7.45, 7.45, 7.45,
      6.5
    ),
    dpi = c(
      600L, 600L, 600L, 600L, 600L, 600L, 300L, 600L,
      600L, 600L, 600L, 600L, 600L, 600L, 600L, 600L
    ),
    stringsAsFactors = FALSE
  )
}

validate_manuscript_figure_manifest <- function(manifest) {
  need_cols(
    manifest,
    c("id", "filename", "stage", "content", "width", "height", "dpi"),
    "manuscript figure manifest"
  )
  assert(nrow(manifest) > 0L, "manuscript figure manifest contains no rows.")
  assert(!anyDuplicated(manifest$id), "manuscript figure manifest contains duplicate ids.")
  assert(!anyDuplicated(manifest$filename), "manuscript figure manifest contains duplicate filenames.")
  assert(
    all(nzchar(as.character(manifest$id))) &&
      all(grepl("\\.png$", as.character(manifest$filename), ignore.case = TRUE)),
    "manuscript figure manifest ids must be non-empty and filenames must end in .png."
  )
  filenames <- as.character(manifest$filename)
  assert(
    all(!startsWith(filenames, "/")) &&
      all(!grepl("(^|[/\\\\])\\.\\.([/\\\\]|$)", filenames)),
    "manuscript figure manifest filenames must be safe relative paths."
  )
  assert(
    all(nzchar(as.character(manifest$stage))) &&
      all(nzchar(as.character(manifest$content))),
    "manuscript figure manifest stage and content fields must be non-empty."
  )
  assert(
    all(is.finite(manifest$width) & manifest$width > 0) &&
      all(is.finite(manifest$height) & manifest$height > 0) &&
      all(is.finite(manifest$dpi) & manifest$dpi > 0 & manifest$dpi == floor(manifest$dpi)),
    "manuscript figure manifest width, height, and dpi must be positive; dpi must be an integer."
  )
  invisible(TRUE)
}

save_manuscript_figure <- function(plot, id, figure_dir = "Figures") {
  assert(inherits(plot, c("gg", "ggplot", "grob", "gtable")), paste0("Figure '", id, "' is not a ggplot-compatible object."))
  manifest <- manuscript_figure_manifest()
  validate_manuscript_figure_manifest(manifest)
  spec <- manifest[manifest$id == id, , drop = FALSE]
  assert(nrow(spec) == 1L, paste0("Unknown manuscript figure id: ", id))

  path <- file.path(figure_dir, spec$filename)
  ensure_dir(dirname(path))
  ggplot2::ggsave(
    filename = path,
    plot = plot,
    width = spec$width,
    height = spec$height,
    dpi = spec$dpi,
    units = "in",
    bg = "white"
  )
  need_file(path, paste0("saved figure '", id, "'"))
  assert(file.info(path)$size > 0, paste0("Saved figure is empty: ", path))
  message("Wrote figure: ", path)
  invisible(path)
}

methods_figure_palette <- function() {
  c(
    ink = "#23343B",
    text = "#1F292E",
    muted_text = "#5A666D",
    axis = "#526068",
    grid = "#E8EDF0",
    pale_grid = "#F2F5F6",
    panel_border = "#D5DDE1",
    land = "#FAFBF9",
    main = "#2E78AF",
    main_dark = "#23343B",
    main_fill = "#91AFBE",
    sensitivity = "#69757B",
    sensitivity_light = "#B6C0C5",
    benchmark = "#B7561D",
    warm = "#C47B4C",
    habitat = "#8DB59B",
    distribution = "#C986B1",
    intersection = "#315F6B",
    retained = "#429477",
    discarded = "#D0D2CE",
    threshold = "#68747A"
  )
}

methods_probability_breaks <- function() {
  c(0, 0.2, 0.4, 0.6, 0.8, 1)
}

methods_probability_labels <- function() {
  c("0", "0.2", "0.4", "0.6", "0.8", "1")
}

methods_compact_number_labels <- function(x) {
  gsub(
    "K",
    "k",
    scales::label_number(
      accuracy = 0.1,
      drop0trailing = TRUE,
      scale_cut = scales::cut_short_scale()
    )(x),
    fixed = TRUE
  )
}

methods_curve_styles <- function(curves = NULL) {
  if (is.null(curves)) {
    curves <- persistence_curves()
  }
  curves <- as.character(curves)

  main_curve <- main_persistence_curve()

  pal <- methods_figure_palette()
  color_values <- stats::setNames(rep(pal[["sensitivity"]], length(curves)), curves)
  color_values[curves == main_curve] <- pal[["main"]]
  # Named styles remain stable when curves are reordered or subsetted.
  linetype_values <- c(
    q50 = "solid",
    q16 = "dotdash",
    q025 = "dotted",
    q84 = "twodash",
    q975 = "longdash"
  )[curves]

  data.frame(
    curve = curves,
    color = unname(color_values),
    linetype = unname(linetype_values),
    stringsAsFactors = FALSE
  )
}

theme_methods_figure <- function(base_size = 11, legend_position = "bottom") {
  pal <- methods_figure_palette()
  ggplot2::theme_classic(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      text = ggplot2::element_text(colour = pal[["text"]]),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin = ggplot2::margin(8, 9, 7, 9),
      axis.title = ggplot2::element_text(
        face = "plain",
        colour = pal[["text"]],
        size = base_size * 0.92
      ),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 5)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 5)),
      axis.text = ggplot2::element_text(colour = "#333333", size = base_size * 0.78),
      axis.line = ggplot2::element_line(linewidth = 0.38, colour = pal[["axis"]]),
      axis.ticks = ggplot2::element_line(linewidth = 0.34, colour = pal[["axis"]]),
      axis.ticks.length = grid::unit(2.1, "pt"),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = legend_position,
      legend.title = ggplot2::element_text(
        face = "bold",
        colour = pal[["text"]],
        size = base_size * 0.78
      ),
      legend.text = ggplot2::element_text(colour = pal[["text"]], size = base_size * 0.76),
      legend.key = ggplot2::element_blank(),
      legend.spacing.x = grid::unit(0.55, "lines"),
      legend.box.margin = ggplot2::margin(t = 2),
      plot.title.position = "plot",
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0,
        colour = pal[["text"]],
        size = base_size * 1.04,
        margin = ggplot2::margin(b = 2)
      ),
      plot.subtitle = ggplot2::element_text(
        hjust = 0,
        colour = pal[["muted_text"]],
        size = base_size * 0.78,
        lineheight = 1.05,
        margin = ggplot2::margin(b = 6)
      ),
      plot.caption = ggplot2::element_text(
        hjust = 0,
        colour = pal[["muted_text"]],
        size = base_size * 0.70,
        margin = ggplot2::margin(t = 5)
      ),
      strip.background = ggplot2::element_rect(
        fill = pal[["land"]],
        colour = pal[["panel_border"]],
        linewidth = 0.32
      ),
      strip.text = ggplot2::element_text(
        face = "bold",
        colour = pal[["text"]],
        size = base_size * 0.82,
        margin = ggplot2::margin(4, 4, 4, 4)
      )
    )
}

theme_methods_probability_guides <- function(linewidth = 0.25) {
  pal <- methods_figure_palette()
  ggplot2::theme(
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_line(
      colour = pal[["pale_grid"]],
      linewidth = linewidth
    ),
    panel.grid.minor = ggplot2::element_blank()
  )
}

theme_methods_map <- function(base_size = 10.5) {
  pal <- methods_figure_palette()
  ggplot2::theme_void(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      text = ggplot2::element_text(colour = pal[["text"]]),
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        size = base_size * 1.02,
        face = "bold",
        margin = ggplot2::margin(b = 1.5)
      ),
      plot.subtitle = ggplot2::element_text(
        hjust = 0.5,
        size = base_size * 0.76,
        colour = pal[["muted_text"]],
        margin = ggplot2::margin(b = 2.2)
      ),
      plot.title.position = "plot",
      plot.margin = ggplot2::margin(3, 3, 3, 3),
      panel.border = ggplot2::element_rect(color = pal[["panel_border"]], fill = NA, linewidth = 0.30),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA)
    )
}

panel_tag_grid <- function(plotlist,
                           labels,
                           ncol,
                           label_size = 11,
                           label_x = 0.018,
                           label_y = 0.988,
                           ...) {
  cowplot::plot_grid(
    plotlist = plotlist,
    ncol = ncol,
    labels = labels,
    label_fontface = "bold",
    label_size = label_size,
    label_colour = methods_figure_palette()[["text"]],
    label_x = label_x,
    label_y = label_y,
    hjust = 0,
    vjust = 1,
    align = "hv",
    axis = "tblr",
    ...
  )
}
