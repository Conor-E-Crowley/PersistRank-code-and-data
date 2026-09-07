# Stage 5.1 spatial-process main-text figure.
#
# This module converts the fixed focal species' retained intermediate rasters,
# patches, and population units into the finalized multipanel map and installs
# the PNG atomically. It does not modify Stage 5 lookup, connectivity, or patch
# raster artifacts.


madagascar_outline <- function(template_rast) {
  mg <- rnaturalearth::ne_countries(
    country = "Madagascar",
    scale = "medium",
    returnclass = "sf"
  )
  sf::st_transform(mg, crs = terra::crs(template_rast, proj = TRUE))
}

tight_extent_from_sf <- function(x, pad_x = 0.035, pad_y = 0.025) {
  bb <- sf::st_bbox(x)
  dx <- unname(bb["xmax"] - bb["xmin"])
  dy <- unname(bb["ymax"] - bb["ymin"])
  terra::ext(
    bb["xmin"] - pad_x * dx,
    bb["xmax"] + pad_x * dx,
    bb["ymin"] - pad_y * dy,
    bb["ymax"] + pad_y * dy
  )
}

geom_spatraster_silent <- function(...) suppressMessages(tidyterra::geom_spatraster(...))

candidate_pu_pal <- function(n) {
  if (n <= 0) return(character(0))
  grDevices::hcl(
    h = seq(15, 375, length.out = n + 1)[seq_len(n)],
    c = rep(c(46, 40, 43), length.out = n),
    l = rep(c(68, 63, 66), length.out = n),
    fixup = TRUE
  )
}

final_pu_pal <- function(n) {
  if (n <= 0) return(character(0))
  base_pal <- c("#4E79A7", "#59A14F", "#7B6FB2", "#3F7F93")
  if (n <= length(base_pal)) return(base_pal[seq_len(n)])
  c(base_pal, grDevices::hcl.colors(n - length(base_pal), palette = "Dark 3"))
}

theme_map <- function() {
  theme_methods_map(base_size = 11.2) +
    ggplot2::theme(plot.margin = ggplot2::margin(2, 2, 2, 2))
}

coord_ext <- function(r, e) {
  ggplot2::coord_sf(
    crs = terra::crs(r, proj = TRUE),
    xlim = c(e$xmin, e$xmax),
    ylim = c(e$ymin, e$ymax),
    expand = FALSE
  )
}

add_mg_fill <- function(mg_sf, fill = "white") {
  ggplot2::geom_sf(data = mg_sf, fill = fill, color = NA, inherit.aes = FALSE)
}

add_mg_outline <- function(mg_sf, col = "grey15", lwd = 0.34) {
  ggplot2::geom_sf(data = mg_sf, fill = NA, color = col, linewidth = lwd, inherit.aes = FALSE)
}

bin_plot <- function(r01, title, fill_col, mg_sf, e, maxcell = 8e5, land_fill = "white") {
  r01 <- bin01(r01)
  ggplot2::ggplot() +
    add_mg_fill(mg_sf, fill = land_fill) +
    geom_spatraster_silent(data = r01, maxcell = maxcell, na.rm = TRUE) +
    ggplot2::scale_fill_gradient(low = fill_col, high = fill_col, na.value = NA, guide = "none") +
    add_mg_outline(mg_sf) +
    coord_ext(r01, e) +
    ggplot2::labs(title = title) +
    theme_map()
}

patch_filter_plot <- function(status_r, title, mg_sf, e, retained_col, maxcell = 8e5, land_fill = "white") {
  # Panel (c) already shows the complete mapped-habitat layer. Showing only
  # retained patches here makes panel (d) the unambiguous post-filter state.
  retained_r <- terra::ifel(status_r == 2, 1, NA)
  ggplot2::ggplot() +
    add_mg_fill(mg_sf, fill = land_fill) +
    geom_spatraster_silent(data = retained_r, maxcell = maxcell, na.rm = TRUE) +
    ggplot2::scale_fill_gradient(
      low = retained_col,
      high = retained_col,
      na.value = NA,
      guide = "none"
    ) +
    add_mg_outline(mg_sf) +
    coord_ext(retained_r, e) +
    ggplot2::labs(title = title) +
    theme_map()
}

unit_plot <- function(r_id,
                      title,
                      mg_sf,
                      e,
                      palette,
                      outline_sf = NULL,
                      maxcell = 8e5,
                      land_fill = "white") {
  vals <- sort(unique(terra::values(r_id, mat = FALSE)))
  vals <- vals[!is.na(vals)]
  labs <- paste0("PU_", vals)
  r_plot <- terra::as.factor(r_id)
  if (length(vals)) levels(r_plot) <- data.frame(value = vals, label = labs)

  pal <- if (is.null(names(palette))) {
    stats::setNames(rep(palette, length.out = length(vals)), labs)
  } else {
    palette[labs]
  }
  missing_cols <- is.na(pal)
  if (any(missing_cols)) pal[missing_cols] <- candidate_pu_pal(sum(missing_cols))

  p <- ggplot2::ggplot() +
    add_mg_fill(mg_sf, fill = land_fill) +
    geom_spatraster_silent(data = r_plot, maxcell = maxcell, na.rm = TRUE) +
    ggplot2::scale_fill_manual(values = pal, guide = "none", na.translate = FALSE) +
    add_mg_outline(mg_sf)

  p <- p +
    coord_ext(r_plot, e) +
    ggplot2::labs(title = title) +
    theme_map()

  if (!is.null(outline_sf)) {
    p <- p + ggplot2::geom_sf(
      data = outline_sf,
      fill = NA,
      color = scales::alpha("grey10", 0.12),
      linewidth = 0.11,
      inherit.aes = FALSE
    )
  }

  p
}

build_single_species_layers <- function(target_row, paths, roi,
                                        clump_backend = "terra") {
  context <- load_patch_raster_context(paths$landcover_tif, roi)
  mapped <- mapped_habitat_for_species(
    target_row,
    context$habitat_masks,
    context$template,
    cell_area_km2 = context$cell_area_km2
  )
  assert(identical(mapped$status, "ok"),
         "Target species has no habitat labels represented in the land-cover crosswalk.")
  assert(mapped$mapped_habitat_cells > 0L,
         "Target species has empty mapped habitat after intersecting habitat and distribution.")
  components <- patch_components_for_species(
    target_row,
    mapped$mapped_habitat,
    species_name = target_row$scientificName,
    mapped_habitat_cells = mapped$mapped_habitat_cells,
    mapped_habitat_area_km2 = mapped$mapped_habitat_area_km2,
    cell_area_km2 = context$cell_area_km2,
    clump_backend = clump_backend,
    retain_intermediates = TRUE
  )
  assert(
    identical(components$status, "retained"),
    paste0("Target species could not produce final population units: ", components$stop_reason)
  )
  mg_sf <- madagascar_outline(mapped$habitat01)
  plot_ext <- tight_extent_from_sf(mg_sf, pad_x = 0.035, pad_y = 0.025)
  list(
    target_row = target_row,
    labels_mixed = mapped$labels,
    habitat01 = mapped$habitat01,
    presence01 = mapped$presence01,
    mapped_habitat01 = mapped$mapped_habitat,
    mapped_habitat_area_km2 = mapped$mapped_habitat_area_km2,
    patch_filter_status = components$patch_filter_status,
    pu_candidate = components$pu_candidate,
    pu_final = components$pu_final,
    final_units_sf = components$final_units_sf,
    mg_sf = mg_sf,
    plot_ext = plot_ext,
    n_raw_patches = components$n_raw_patches,
    n_kept_patches = components$n_patches_after_patch_filter,
    n_candidate_pus = components$n_candidate_pus,
    n_final_pus = components$n_final_pus,
    keep_pu = components$kept_candidate_pus,
    patch_lookup = components$patch_lookup,
    connectivity = components$connectivity,
    patch_final = components$patch_final,
    habitat_c = terra::crop(mapped$habitat01, plot_ext, snap = "out"),
    presence_c = terra::crop(mapped$presence01, plot_ext, snap = "out"),
    mapped_habitat_c = terra::crop(mapped$mapped_habitat, plot_ext, snap = "out"),
    patch_filter_c = terra::crop(components$patch_filter_status, plot_ext, snap = "out"),
    pu_candidate_c = terra::crop(components$pu_candidate, plot_ext, snap = "out"),
    pu_final_c = terra::crop(components$pu_final, plot_ext, snap = "out")
  )
}

make_single_species_process_figure <- function(layers) {
  pal <- methods_figure_palette()
  # Slightly stronger than the shared habitat green so fine cells survive
  # reduction to a narrow manuscript column.
  col_hab <- "#76A489"
  col_dist <- pal[["distribution"]]
  col_aoh <- pal[["intersection"]]
  col_keep <- pal[["retained"]]
  land_fill <- pal[["land"]]

  final_cols <- stats::setNames(final_pu_pal(layers$n_final_pus), paste0("PU_", seq_len(layers$n_final_pus)))
  candidate_cols <- stats::setNames(candidate_pu_pal(layers$n_candidate_pus), paste0("PU_", seq_len(layers$n_candidate_pus)))

  candidate_keep_labels <- paste0("PU_", layers$keep_pu)
  final_keep_labels <- paste0("PU_", seq_len(layers$n_final_pus))
  candidate_cols[candidate_keep_labels] <- unname(final_cols[final_keep_labels])

  panels <- list(
    bin_plot(layers$habitat_c, "Suitable habitat", col_hab, layers$mg_sf, layers$plot_ext, land_fill = land_fill),
    bin_plot(layers$presence_c, "Species distribution", col_dist, layers$mg_sf, layers$plot_ext, land_fill = land_fill),
    bin_plot(layers$mapped_habitat_c, "Mapped habitat", col_aoh, layers$mg_sf, layers$plot_ext, land_fill = land_fill),
    patch_filter_plot(
      layers$patch_filter_c,
      "Retained patches",
      layers$mg_sf,
      layers$plot_ext,
      col_keep,
      land_fill = land_fill
    ),
    unit_plot(
      layers$pu_candidate_c,
      "Candidate units",
      layers$mg_sf,
      layers$plot_ext,
      candidate_cols,
      land_fill = land_fill
    ),
    unit_plot(
      layers$pu_final_c,
      "Final units",
      layers$mg_sf,
      layers$plot_ext,
      final_cols,
      outline_sf = layers$final_units_sf,
      land_fill = land_fill
    )
  )

  panel_tag_grid(
    plotlist = panels,
    labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"),
    ncol = 3,
    label_size = 11.5,
    label_x = c(0.036, 0.015, 0.036, 0.036, 0.036, 0.036),
    label_y = 0.982
  )
}

write_stage51_figure_atomic <- function(plot, config, .rename_file = file.rename) {
  staging_dir <- tempfile("stage51_figure_", tmpdir = config$paths$figure_dir)
  dir.create(staging_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(staging_dir, recursive = TRUE, force = TRUE), add = TRUE)
  suppressMessages(save_manuscript_figure(plot, "spatial_process", staging_dir))
  staged <- file.path(staging_dir, basename(config$paths$figure))
  assert(file.exists(staged) && file.info(staged)$size > 0,
         "Stage 5.1 produced an empty figure file.")
  transaction <- stage5_file_set_transaction(
    staged, config$paths$figure,
    overwrite = TRUE,
    rename_file = .rename_file
  )
  transaction$finalize()
  invisible(config$paths$figure)
}
