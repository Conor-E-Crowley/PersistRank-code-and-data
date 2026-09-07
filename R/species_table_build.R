# Final Stage 4 species-table assembly, validation, and output.
#
# The functions in this module join the already resolved ecological/Gompertz
# rows to the queried habitat rows without changing species order, enforce the
# complete public CSV schema and scientific units, and install a validated CSV
# through the shared rollback-capable transaction. They perform no name
# matching, raster discovery, model fitting, or IUCN requests.

species_table_keep_cols <- function(curves) {
  c(
    "scientificName", "className", "orderName", "familyName", "genusName", "speciesName",
    "redlistCategory",
    "raster_stem", "name_raster", "match_raster", "sdm_method", "raster_dir", "raster_file", "raster_path",
    "name_trait", "match_trait",
    "BodyMass.Value", "GenLength", "bird_genlength_scientific", "genlength_match_source",
    "Diet", "diet5_group", "bird_sigma_model_group",
    "density", "home_range_size", "dispersal_dist", "min_patch_size", "min_pop_size",
    species_table_curve_columns(curves),
    "habitat_codes_suitable", "habitats_level1", "habitats_mixed"
  )
}

validate_stage4_species_parameters <- function(species_data, curves,
                                               use_gompertz_mammals,
                                               use_gompertz_birds,
                                               label = "Stage 4 species parameters") {
  class_col <- if ("taxon_class" %in% names(species_data)) "taxon_class" else "className"
  required <- c(
    "scientificName", class_col, "BodyMass.Value", "Diet",
    "diet5_group", "bird_sigma_model_group",
    "density", "home_range_size", "dispersal_dist", "min_patch_size", "min_pop_size",
    "match_raster", "match_trait", "raster_stem", "raster_dir", "raster_file", "raster_path",
    "GenLength", "bird_genlength_scientific", "genlength_match_source"
  )
  need_cols(species_data, required, label)
  assert(nrow(species_data) > 0L, paste0(label, " contains no retained species."))

  scientific_names <- stringr::str_squish(as.character(species_data$scientificName))
  bad_names <- is.na(scientific_names) | !nzchar(scientific_names)
  assert(!any(bad_names), paste0(label, " contains missing or blank scientificName values."))
  duplicates <- unique(scientific_names[duplicated(clean_name(scientific_names))])
  if (length(duplicates)) {
    species_abort(label, " contains duplicate retained scientific names:\n",
                  paste(utils::head(duplicates, 20L), collapse = "\n"))
  }

  class_key <- clean_name(species_data[[class_col]])
  is_mammal <- class_key == "mammalia"
  is_bird <- class_key == "aves"
  assert(all(is_mammal | is_bird), paste0(label, " contains unsupported taxon classes."))

  positive <- c("BodyMass.Value", "density", "home_range_size", "min_patch_size", "min_pop_size")
  for (column in positive) {
    value <- suppressWarnings(as.numeric(species_data[[column]]))
    assert(all(is.finite(value) & value > 0),
           paste0(label, " ", column, " must be positive and finite for every retained species."))
  }
  dispersal <- suppressWarnings(as.numeric(species_data$dispersal_dist))
  assert(all(is.finite(dispersal) & dispersal >= 0),
         paste0(label, " dispersal_dist must be finite and nonnegative for every retained species."))

  diet <- trimws(as.character(species_data$Diet))
  assert(all(diet[is_mammal] %in% names(mammal_density_coefficients$diet)),
         paste0(label, " mammal Diet must contain only herbivore, omnivore, or carnivore."))
  assert(all(diet[is_bird] %in% names(bird_density_coefficients$diet)),
         paste0(label, " bird Diet contains a category unsupported by the bird density model."))

  diet5_group <- trimws(as.character(species_data$diet5_group))
  model_group <- trimws(as.character(species_data$bird_sigma_model_group))
  mammal_diet_present <- is_mammal & !is.na(diet5_group) & nzchar(diet5_group)
  mammal_model_present <- is_mammal & !is.na(model_group) & nzchar(model_group)
  assert(
    !any(mammal_diet_present | mammal_model_present),
    paste0(label, " mammal bird diet/model group values must be absent.")
  )
  if (any(is_bird)) {
    validate_bird_diet5(diet5_group[is_bird], paste0(label, " diet5_group"))
    assert(
      identical(diet5_group[is_bird], diet[is_bird]),
      paste0(label, " diet5_group is inconsistent with Diet.")
    )
    if (isTRUE(use_gompertz_birds)) {
      assert(
        all(!is.na(model_group[is_bird]) & nzchar(model_group[is_bird])),
        paste0(label, " birds require a Stage 3 model branch.")
      )
    }
    generation_length <- suppressWarnings(as.numeric(species_data$GenLength[is_bird]))
    assert(all(is.finite(generation_length) & generation_length > 0),
           paste0(label, " bird GenLength must be positive and finite."))
    resolved_generation_name <- stringr::str_squish(as.character(
      species_data$bird_genlength_scientific[is_bird]
    ))
    assert(all(!is.na(resolved_generation_name) & nzchar(resolved_generation_name)),
           paste0(label, " birds must have a resolved Bird et al. generation-length scientific name."))
    generation_source <- as.character(species_data$genlength_match_source[is_bird])
    assert(all(generation_source %in% c("direct", "synonym", "input_synonym")),
           paste0(label, " bird generation-length match provenance is invalid."))
  }

  for (column in intersect(c("Effect_Order", "Effect_Family", "Effect_Species"), names(species_data))) {
    assert(all(is.finite(suppressWarnings(as.numeric(species_data[[column]])))),
           paste0(label, " ", column, " must be finite after zero fallback."))
  }

  density <- suppressWarnings(as.numeric(species_data$density))
  validate_density_area_thresholds(
    density = density,
    min_patch_area_km2 = suppressWarnings(as.numeric(species_data$min_patch_size)),
    min_population_area_km2 = suppressWarnings(as.numeric(species_data$min_pop_size)),
    label = paste0(label, " area thresholds")
  )

  assert(all(as.character(species_data$match_raster) %in% c("original", "synonym", "input_synonym")),
         paste0(label, " contains invalid SDM raster match provenance."))
  assert(all(as.character(species_data$match_trait) %in% c("original", "synonym", "input_synonym")),
         paste0(label, " contains invalid EltonTraits match provenance."))
  for (column in c("raster_stem", "raster_dir", "raster_file", "raster_path")) {
    value <- stringr::str_squish(as.character(species_data[[column]]))
    assert(all(!is.na(value) & nzchar(value)), paste0(label, " ", column, " contains blank values."))
  }
  missing_rasters <- !file.exists(as.character(species_data$raster_path))
  if (any(missing_rasters)) {
    species_abort(label, " references missing raster files:\n",
                  paste(utils::head(species_data$raster_path[missing_rasters], 20L), collapse = "\n"))
  }

  enabled <- (is_mammal & isTRUE(use_gompertz_mammals)) |
    (is_bird & isTRUE(use_gompertz_birds))
  disabled <- (is_mammal & !isTRUE(use_gompertz_mammals)) |
    (is_bird & !isTRUE(use_gompertz_birds))
  curve_columns <- species_table_curve_columns(curves)
  need_cols(species_data, curve_columns, label)
  if (any(enabled)) {
    validate_species_table_curve_parameters(
      species_data, curves, rows = enabled, label = label,
      rows_label = "Gompertz-enabled species"
    )
  }
  if (any(disabled)) {
    populated <- vapply(curve_columns, function(column) {
      any(disabled & !is.na(species_data[[column]]))
    }, logical(1))
    assert(
      !any(populated),
      paste0(
        label,
        " contains Gompertz parameters for taxa not selected by Stage 4: ",
        paste(curve_columns[populated], collapse = ", "), "."
      )
    )
  }
  invisible(TRUE)
}

validate_species_table <- function(species_table, curves, use_gompertz_mammals, use_gompertz_birds) {
  if (!nrow(species_table)) {
    species_abort("No species rows remained after resolving rasters and trait records.")
  }

  expected_cols <- species_table_keep_cols(curves)
  need_cols(species_table, expected_cols, "species_table.csv")
  assert(
    identical(names(species_table), expected_cols),
    "species_table.csv columns must exactly match the documented Stage 4 schema and order."
  )
  validate_stage4_species_parameters(
    species_table, curves, use_gompertz_mammals, use_gompertz_birds,
    label = "species_table.csv"
  )
  valid_methods <- c("PPM", "RangeBag")
  assert(
    all(species_table$sdm_method %in% valid_methods),
    "species_table.csv contains an unsupported sdm_method value."
  )
  required_text <- c("scientificName", "className", "name_trait", "name_raster", "habitats_mixed")
  for (col in required_text) {
    bad <- is.na(species_table[[col]]) | !nzchar(trimws(as.character(species_table[[col]])))
    if (any(bad)) {
      species_abort(col, " contains missing or blank values for ", sum(bad), " retained species.")
    }
  }
  missing_habitat <- is.na(species_table$habitats_mixed) | !nzchar(species_table$habitats_mixed)
  if (any(missing_habitat)) {
    species_abort(
      "No suitable IUCN habitat label was returned for the following retained species:\n",
      paste(utils::head(species_table$scientificName[missing_habitat], 25L), collapse = "\n"),
      if (sum(missing_habitat) > 25L) "\n..." else ""
    )
  }

  invisible(TRUE)
}


assemble_species_table <- function(species_inputs, habitats, curves,
                                   use_gompertz_mammals, use_gompertz_birds) {
  n_before <- nrow(species_inputs)
  joined <- species_inputs |>
    dplyr::left_join(habitats, by = c("genusName", "speciesName"))
  assert(nrow(joined) == n_before, "Joining IUCN habitats unexpectedly changed the species row count.")
  need_cols(joined, species_table_keep_cols(curves), "assembled Stage 4 species data")
  species_table <- joined |>
    dplyr::select(dplyr::all_of(species_table_keep_cols(curves)))
  validate_species_table(species_table, curves, use_gompertz_mammals, use_gompertz_birds)
  species_table
}

write_species_table <- function(species_table,
                                table_path,
                                curves,
                                use_gompertz_mammals,
                                use_gompertz_birds,
                                overwrite = FALSE,
                                rename_file = file.rename) {
  ensure_writable_dir(dirname(table_path), "species table output directory")
  staged_table <- tempfile("stage4_species_table_", tmpdir = dirname(table_path), fileext = ".csv")
  on.exit(unlink(staged_table, force = TRUE), add = TRUE)
  readr::write_csv(species_table, staged_table)
  written <- readr::read_csv(staged_table, show_col_types = FALSE, progress = FALSE)
  validate_species_table(written, curves, use_gompertz_mammals, use_gompertz_birds)
  assert(nrow(written) == nrow(species_table), "Temporary species table has an unexpected row count.")
  project_file_set_transaction(
    staged_paths = staged_table,
    target_paths = table_path,
    overwrite = overwrite,
    rename_file = rename_file,
    label = "Stage 4 species table and metadata"
  )
  log_msg(
    "Stage 4 | species table written | rows=", nrow(written),
    "| table=", table_path
  )
  invisible(written)
}
