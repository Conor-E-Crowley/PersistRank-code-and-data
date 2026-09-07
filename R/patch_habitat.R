# Shared Stage 5 land-cover and SDM raster preparation.
#
# Land cover is cropped once to the fixed ROI. Reusable habitat masks, the
# common raster template, and cell areas are then supplied to each sequential
# species calculation.

landcover_class_map <- function() {
  class_map <- list(
    forest                      = c(50L, 60L, 61L, 62L, 70L, 71L, 72L, 90L, 160L),
    savanna                     = c(120L, 121L, 122L),
    shrubland                   = c(120L, 121L, 122L, 150L, 151L, 152L, 153L, 200L, 201L, 202L),
    grassland                   = c(130L, 140L, 150L, 151L, 152L, 153L),
    wetlands_inland             = c(20L, 80L, 81L, 82L, 160L, 170L, 180L, 210L),
    rocky_areas                 = c(70L, 71L, 72L, 130L, 150L, 151L, 152L, 153L, 200L, 201L, 202L),
    desert                      = c(150L, 151L, 152L, 153L, 200L, 201L, 202L),
    arable_pastureland          = c(11L, 20L, 190L),
    plantations_degraded_forest = integer(0),
    urban_rural_gardens         = c(190L),
    artificial_aquatic          = c(20L, 160L, 170L, 180L, 190L, 210L)
  )
  class_map$artificial_terrestrial <- unique(c(
    class_map$arable_pastureland,
    class_map$urban_rural_gardens,
    class_map$plantations_degraded_forest
  ))
  class_map
}

habitat_label_map <- function() {
  c(
    "Forest"                                       = "forest",
    "Savanna"                                      = "savanna",
    "Shrubland"                                    = "shrubland",
    "Grassland"                                    = "grassland",
    "Wetlands (inland)"                            = "wetlands_inland",
    "Rocky Areas"                                  = "rocky_areas",
    "Desert"                                       = "desert",
    "Arable & Pastureland"                         = "arable_pastureland",
    "Plantations & Heavily Degraded Former Forest" = "plantations_degraded_forest",
    "Urban & Rural Gardens"                        = "urban_rural_gardens",
    "Artificial - Aquatic"                         = "artificial_aquatic",
    "Artificial - Terrestrial"                     = "artificial_terrestrial"
  )
}

parse_habitats_mixed <- function(x) {
  out <- unlist(strsplit(as.character(x), ",", fixed = TRUE))
  out <- stringr::str_squish(out)
  out <- out[nzchar(out)]
  unique(out)
}

build_habitat_masks_from_landcover <- function(landcover, class_map = landcover_class_map()) {
  base <- landcover[[1]]
  layers <- lapply(names(class_map), function(mask_name) {
    codes <- class_map[[mask_name]]
    r <- if (length(codes) == 0L || all(is.na(codes))) {
      base * NA_integer_
    } else {
      terra::ifel(base %in% codes, 1L, NA_integer_)
    }
    names(r) <- mask_name
    r
  })
  terra::rast(layers)
}

load_patch_raster_context <- function(landcover_tif, study_area) {
  landcover <- terra::rast(landcover_tif)
  assert(terra::nlyr(landcover) >= 1L, "ESA CCI land-cover raster contains no layers.")
  assert(nzchar(terra::crs(landcover, proj = TRUE)), "ESA CCI land-cover raster has no valid CRS.")
  assert(all(is.finite(terra::res(landcover)) & terra::res(landcover) > 0),
         "ESA CCI land-cover raster has invalid resolution.")
  assert(is.list(study_area) && study_area$mode %in%
           c("bounds", "vector", "full_raster"),
         "Invalid Stage 5 study-area definition.")
  crop_target <- if (identical(study_area$mode, "bounds")) {
    patch_roi(study_area$bounds)
  } else if (identical(study_area$mode, "vector")) {
    need_file(study_area$file, "study-area vector")
    vector <- sf::st_read(study_area$file, layer = study_area$layer, quiet = TRUE)
    assert(nrow(vector) > 0L && all(sf::st_geometry_type(vector) %in% c("POLYGON", "MULTIPOLYGON")),
           "Study-area vector must contain at least one polygon or multipolygon.")
    vector <- sf::st_make_valid(vector)
    vector <- sf::st_transform(sf::st_union(vector), terra::crs(landcover, proj = TRUE))
    terra::vect(vector)
  } else NULL
  landcover_roi <- tryCatch(
    if (is.null(crop_target)) landcover[[1]] else terra::crop(landcover[[1]], crop_target, snap = "out"),
    error = function(e) patch_abort("ESA CCI land-cover raster does not overlap the configured ROI: ", conditionMessage(e))
  )
  if (identical(study_area$mode, "vector")) {
    landcover_roi <- terra::mask(landcover_roi, crop_target)
  }
  assert(terra::ncell(landcover_roi) > 0L, "ESA CCI land-cover crop contains no cells.")
  observed <- terra::freq(landcover_roi, digits = 0)
  observed_codes <- if (is.null(observed) || !nrow(observed)) numeric() else as.numeric(observed$value)
  mapped_codes <- unique(unlist(landcover_class_map(), use.names = FALSE))
  assert(
    any(observed_codes %in% mapped_codes),
    "ESA CCI land-cover crop contains no class codes represented in the Stage 5 habitat crosswalk."
  )
  habitat_masks <- build_habitat_masks_from_landcover(landcover_roi, landcover_class_map())
  template <- habitat_masks[[1]]
  cell_area_km2 <- terra::cellSize(template, unit = "km")
  names(cell_area_km2) <- "cell_area_km2"
  list(
    habitat_masks = habitat_masks,
    template = template,
    cell_area_km2 = cell_area_km2,
    observed_landcover_codes = sort(unique(observed_codes))
  )
}

non_na_cell_count <- function(r) {
  counts <- terra::global(!is.na(r), "sum", na.rm = TRUE)
  as.numeric(sum(as.numeric(counts), na.rm = TRUE))
}

load_presence_aligned01_diagnostics <- function(path, template) {
  r <- terra::rast(path)
  assert(terra::nlyr(r) == 1L, paste0("Species SDM must contain exactly one layer: ", path))
  assert(nzchar(terra::crs(r, proj = TRUE)), paste0("Species SDM has no valid CRS: ", path))
  source_total <- as.numeric(terra::ncell(r))
  source_presence <- terra::global(r == 1L, "sum", na.rm = TRUE)
  source_non_na <- non_na_cell_count(r)
  projected <- FALSE
  resampled <- FALSE

  if (!terra::compareGeom(r, template, stopOnError = FALSE)) {
    if (!terra::same.crs(r, template)) {
      r <- terra::project(r, template, method = "near")
      projected <- TRUE
    }
    if (!terra::compareGeom(r, template, stopOnError = FALSE)) {
      r <- terra::resample(r, template, method = "near")
      resampled <- TRUE
    }
  }

  presence01 <- terra::ifel(r == 1L, 1L, NA_integer_)
  list(
    presence01 = presence01,
    source_sdm_layers = 1L,
    source_sdm_total_cells = source_total,
    source_sdm_non_na_cells = source_non_na,
    source_sdm_presence_cells = as.numeric(source_presence[[1L]]),
    alignment_action = if (projected) {
      if (resampled) "projected+resampled" else "projected"
    } else if (resampled) {
      "resampled"
    } else {
      "unchanged"
    },
    aligned_presence_cells = non_na_cell_count(presence01)
  )
}

habitat_union <- function(masks, mask_names) {
  x <- masks[[mask_names]]
  if (terra::nlyr(x) == 1L) return(x)

  terra::app(x, fun = function(v) {
    if (is.matrix(v)) {
      out <- rep.int(NA_integer_, nrow(v))
      out[rowSums(v == 1L, na.rm = TRUE) > 0] <- 1L
      out
    } else {
      if (any(v == 1L, na.rm = TRUE)) 1L else NA_integer_
    }
  })
}

mapped_habitat_for_species <- function(row, habitat_masks, template,
                                       cell_area_km2 = NULL,
                                       label_map = habitat_label_map()) {
  habitat_labels <- parse_habitats_mixed(row$habitats_mixed)
  mask_names <- unname(label_map[habitat_labels])
  mask_names <- unique(mask_names[!is.na(mask_names)])
  mask_names <- mask_names[mask_names %in% names(habitat_masks)]

  if (!length(mask_names)) {
    return(list(status = "skipped_no_mapped_habitat_labels", labels = habitat_labels, mask_names = mask_names))
  }

  habitat01 <- habitat_union(habitat_masks, mask_names)
  presence_diagnostics <- load_presence_aligned01_diagnostics(row$raster_path, template)
  presence01 <- presence_diagnostics$presence01
  mapped_habitat <- terra::ifel(habitat01 == 1L & presence01 == 1L, 1L, NA_integer_)

  list(
    status = "ok",
    labels = habitat_labels,
    mask_names = mask_names,
    alignment_action = presence_diagnostics$alignment_action,
    aligned_presence_cells = presence_diagnostics$aligned_presence_cells,
    mapped_habitat_cells = non_na_cell_count(mapped_habitat),
    mapped_habitat_area_km2 = raster_area_km2(mapped_habitat, cell_area_km2),
    source_sdm_layers = presence_diagnostics$source_sdm_layers,
    source_sdm_total_cells = presence_diagnostics$source_sdm_total_cells,
    source_sdm_non_na_cells = presence_diagnostics$source_sdm_non_na_cells,
    source_sdm_presence_cells = presence_diagnostics$source_sdm_presence_cells,
    habitat01 = habitat01,
    presence01 = presence01,
    mapped_habitat = mapped_habitat
  )
}

raster_area_km2 <- function(r01, cell_area_km2 = NULL) {
  cell_area <- if (is.null(cell_area_km2)) terra::cellSize(r01, unit = "km") else cell_area_km2
  if (!terra::compareGeom(cell_area, r01, stopOnError = FALSE)) {
    cell_area <- terra::crop(cell_area, r01, snap = "near")
  }
  area <- terra::global(terra::ifel(is.na(r01), NA_real_, cell_area), "sum", na.rm = TRUE)[1, 1]
  as.numeric(area)
}

bin01 <- function(r) terra::ifel(r == 1L, 1L, NA_integer_)
