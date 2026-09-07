# Stage 5 species-table input contract.
#
# This module reads only the ecological and raster fields required by patch
# construction. Gompertz coefficients are deliberately absent from the
# contract, so taxon selection here is independent of Stage 4 coefficient
# assignment.

required_patch_species_cols <- function(extra = character()) {
  unique(c(
    "scientificName",
    "className",
    "habitats_mixed",
    "raster_path",
    "sdm_method",
    "raster_dir",
    "raster_file",
    "density",
    "min_patch_size",
    "min_pop_size",
    "dispersal_dist",
    extra
  ))
}

read_patch_species_base <- function(species_csv,
                                    contract = canonical_analysis_contract()) {
  contract <- validate_analysis_contract(contract, "Stage 5 abundance contract")
  required_cols <- required_patch_species_cols()
  out <- readr::read_csv(
    species_csv,
    show_col_types = FALSE,
    progress = FALSE,
    col_select = dplyr::all_of(required_cols)
  ) |>
    dplyr::mutate(
      scientificName = stringr::str_squish(as.character(scientificName)),
      class_lc       = stringr::str_to_lower(stringr::str_squish(as.character(className))),
      habitats_mixed = stringr::str_squish(as.character(habitats_mixed)),
      raster_path    = stringr::str_squish(as.character(raster_path)),
      raster_dir     = stringr::str_squish(as.character(raster_dir)),
      raster_file    = stringr::str_squish(as.character(raster_file)),
      sdm_method     = stringr::str_squish(as.character(sdm_method)),
      disp_km        = suppressWarnings(readr::parse_number(as.character(dispersal_dist)))
    ) |>
    add_canonical_area_thresholds(
      validate = TRUE,
      label = "stage-5 species-table area thresholds",
      contract = contract
    ) |>
    dplyr::mutate(
      min_patch_km2 = min_patch_area_km2,
      min_pop_km2 = min_population_area_km2
    )
  need_cols(out, required_cols, "species_table.csv")
  out
}

filter_patch_species_taxa <- function(species_table, do_mammals = TRUE, do_birds = TRUE) {
  selected_species <- species_table |>
    dplyr::filter(
      (class_lc == "mammalia" & isTRUE(do_mammals)) |
        (class_lc == "aves" & isTRUE(do_birds))
    )

  validate_selected_patch_species(selected_species)
  selected_species
}

validate_selected_patch_species <- function(selected_species) {
  assert(nrow(selected_species) > 0L, "No species selected after applying filters.")
  assert(all(nzchar(selected_species$scientificName)), "Selected species contain blank scientificName values.")
  assert(
    all(selected_species$class_lc %in% c("mammalia", "aves")),
    "Selected species contain unsupported className values."
  )
  assert(
    all(selected_species$sdm_method %in% c("PPM", "RangeBag")),
    "Selected species contain unsupported sdm_method values."
  )
  duplicate_species <- sort(unique(selected_species$scientificName[duplicated(selected_species$scientificName)]))
  assert(
    length(duplicate_species) == 0L,
    paste0(
      "Selected species_table.csv rows contain duplicate scientificName values. ",
      "For mixed PPM+RangeBag inputs, rebuild Stage 4 with params$sdm = \"ppm_rangebag\" ",
      "so ambiguous PPM/RangeBag raster matches are diagnosed before Stage 5. First duplicate species: ",
      paste(head(duplicate_species, 20), collapse = ", ")
    )
  )
  assert(all(nzchar(selected_species$habitats_mixed)), "Selected species contain blank habitats_mixed values.")
  assert(all(nzchar(selected_species$raster_path)), "Selected species contain blank raster_path values.")
  assert(
    all(file.exists(selected_species$raster_path)),
    paste0(
      "At least one selected raster_path does not exist. First missing paths:\n",
      paste(head(selected_species$raster_path[!file.exists(selected_species$raster_path)], 20), collapse = "\n")
    )
  )
  assert(
    all(is.finite(selected_species$min_patch_km2) & selected_species$min_patch_km2 > 0),
    "Selected species contain non-positive or non-numeric min_patch_size values."
  )
  assert(
    all(is.finite(selected_species$min_pop_km2) & selected_species$min_pop_km2 > 0),
    "Selected species contain non-positive or non-numeric min_pop_size values."
  )
  assert(
    all(is.finite(selected_species$disp_km) & selected_species$disp_km >= 0),
    "Selected species contain negative or non-numeric dispersal_dist values."
  )

  patch_files <- patch_filename_from_scientific(selected_species$scientificName)
  duplicate_patch_files <- unique(patch_files[duplicated(patch_files)])
  assert(
    length(duplicate_patch_files) == 0L,
    paste0("Duplicate patch-raster filenames would be created: ", paste(duplicate_patch_files, collapse = ", "))
  )

  invisible(TRUE)
}
