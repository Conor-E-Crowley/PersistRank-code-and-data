# Stage 4 input readers and raster-manifest construction.
#
# Candidate construction and diagnostics live in the focused species-name
# modules; species_table_names.R coordinates their results. This module only
# reads canonical raw files, validates schemas and units, inventories selected
# SDM rasters, and returns the inputs needed by the build workflow.

clean_name <- function(x) stringr::str_squish(stringr::str_to_lower(as.character(x)))

as_number <- function(x) suppressWarnings(readr::parse_number(as.character(x)))

species_to_stem <- function(x) {
  paste0(
    stringr::str_replace_all(stringr::str_squish(stringr::str_to_lower(x)), "\\s+", "_"),
    "_bin"
  )
}

read_synonyms <- function(path) {
  raw <- tryCatch(
    readr::read_csv(path, show_col_types = FALSE),
    error = function(e) readr::read_delim(path, delim = "\t", show_col_types = FALSE)
  )
  need_cols(raw, c("scientificName", "genusName", "speciesName"), "synonyms.csv")
  raw |>
    dplyr::mutate(synonym_source_row = dplyr::row_number()) |>
    dplyr::transmute(
      synonym_source_row,
      scientificName = stringr::str_squish(scientificName),
      synonym = stringr::str_squish(paste(genusName, speciesName))
    ) |>
    dplyr::filter(nzchar(scientificName), nzchar(synonym)) |>
    dplyr::distinct(scientificName, synonym, .keep_all = TRUE)
}

read_input_synonyms <- function(path) {
  raw <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  need_cols(raw, c("scientificName", "synonym_scientificName"), "input_synonyms.csv")
  out <- raw |>
    dplyr::mutate(input_synonym_row = dplyr::row_number()) |>
    dplyr::transmute(
      input_synonym_row,
      scientificName = stringr::str_squish(scientificName),
      synonym = stringr::str_squish(synonym_scientificName)
    ) |>
    dplyr::filter(nzchar(scientificName), nzchar(synonym))
  assert(
    nrow(out) == nrow(raw),
    "input_synonyms.csv contains missing or blank scientific-name pairs."
  )
  assert(
    !anyDuplicated(out[c("scientificName", "synonym")]),
    "input_synonyms.csv contains duplicate scientific-name pairs."
  )
  out
}

read_traits <- function(path, taxon_class) {
  raw <- readr::read_tsv(
    path,
    col_types = readr::cols(.default = readr::col_character()),
    na = c("", "NA", "NaN", "NULL"),
    show_col_types = FALSE,
    progress = FALSE
  )
  need_cols(raw, c("Scientific", "BodyMass-Value"), paste0(taxon_class, " EltonTraits"))
  if (identical(taxon_class, "Aves")) {
    assert(
      any(startsWith(names(raw), "Diet-")) || any(names(raw) %in% c("Diet.5Cat", "Diet_5Cat")),
      "Bird EltonTraits must contain Diet-* or Diet-5Cat information."
    )
  }
  out <- raw |>
    dplyr::mutate(
      taxon_class = taxon_class,
      trait_key = clean_name(Scientific),
      BodyMass.Value = as_number(`BodyMass-Value`)
    )
  assert(
    all(!is.na(out$trait_key) & nzchar(out$trait_key)),
    paste0(taxon_class, " EltonTraits contains missing scientific names.")
  )
  assert(
    !anyDuplicated(out$trait_key),
    paste0(taxon_class, " EltonTraits contains duplicate normalized scientific names.")
  )
  assert(
    all(is.finite(out$BodyMass.Value) & out$BodyMass.Value > 0),
    paste0(taxon_class, " EltonTraits contains missing, nonfinite, or nonpositive body masses.")
  )
  out
}

list_raster_manifest <- function(raster_sources) {
  out <- lapply(seq_len(nrow(raster_sources)), function(i) {
    src <- raster_sources[i, , drop = FALSE]
    files <- list.files(src$dir_path, full.names = TRUE)
    manifest <- tibble::tibble(
      taxon_class = src$taxon_class,
      sdm_method = src$sdm_method,
      raster_dir = src$dir_path,
      raster_path = files,
      raster_file = basename(files),
      raster_stem = stringr::str_to_lower(tools::file_path_sans_ext(basename(files))),
      extension = stringr::str_to_lower(tools::file_ext(files))
    ) |>
      dplyr::filter(extension == "tif", stringr::str_ends(raster_stem, "_bin")) |>
      dplyr::select(-extension)
    assert(
      nrow(manifest) > 0L,
      paste0(
        "Selected ", src$taxon_class, " ", src$sdm_method,
        " raster directory contains no qualifying _bin.tif files: ", src$dir_path
      )
    )
    manifest
  })

  manifest <- dplyr::bind_rows(out) |>
    dplyr::arrange(taxon_class, sdm_method, raster_file)
  validate_raster_manifest(manifest)
  manifest
}

read_raster_index <- function(path, selected_methods, selected_classes) {
  need_file(path, "SDM index CSV")
  index <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  need_cols(
    index, c("scientificName", "taxon_class", "sdm_method", "raster_path"),
    "SDM index CSV"
  )
  index <- index |>
    dplyr::mutate(
      scientificName = stringr::str_squish(as.character(scientificName)),
      taxon_class = as.character(taxon_class),
      sdm_method = as.character(sdm_method),
      raster_path = resolve_sdm_index_raster_paths(raster_path, path),
      raster_dir = dirname(raster_path),
      raster_file = basename(raster_path),
      # The index is an explicit species-to-raster mapping; filenames no longer
      # need to follow the folder-mode naming convention.
      raster_stem = species_to_stem(scientificName)
    ) |>
    dplyr::filter(taxon_class %in% selected_classes, sdm_method %in% selected_methods) |>
    dplyr::select(taxon_class, sdm_method, raster_dir, raster_path, raster_file, raster_stem)
  assert(nrow(index) > 0L, "SDM index contains no rows for the selected taxa and SDM methods.")
  assert(!anyDuplicated(index[c("taxon_class", "sdm_method", "raster_stem")]),
         "SDM index contains duplicate taxon/method/species raster rows.")
  missing <- index$raster_path[!file.exists(index$raster_path)]
  assert(!length(missing), paste0("SDM index references missing raster(s):\n",
                                  paste(utils::head(missing, 20L), collapse = "\n")))
  validate_raster_manifest(index)
  index
}

validate_raster_manifest <- function(manifest) {
  need_cols(
    manifest,
    c("taxon_class", "sdm_method", "raster_stem", "raster_file"),
    "Stage 4 raster manifest"
  )
  duplicate_keys <- manifest |>
    dplyr::count(taxon_class, sdm_method, raster_stem, name = "n") |>
    dplyr::filter(n > 1L)
  if (nrow(duplicate_keys)) {
    details <- manifest |>
      dplyr::inner_join(duplicate_keys, by = c("taxon_class", "sdm_method", "raster_stem")) |>
      dplyr::transmute(detail = paste(taxon_class, sdm_method, raster_stem, raster_file, sep = " | ")) |>
      dplyr::pull(detail)
    species_abort(
      "Selected raster directories contain duplicate _bin.tif stems:\n",
      paste(utils::head(details, 25L), collapse = "\n"),
      if (length(details) > 25L) "\n..." else ""
    )
  }
  invisible(TRUE)
}

read_random_effects <- function(path) {
  raw <- readr::read_csv(path, show_col_types = FALSE)
  required <- c("Class", "Order", "Family", "Species", "Effect_Order", "Effect_Family", "Effect_Species")
  need_cols(raw, required, "random_effects.csv")
  out <- raw |>
    dplyr::mutate(
      Class = clean_name(Class),
      Order = clean_name(Order),
      Family = clean_name(Family),
      Species = clean_name(Species),
      dplyr::across(dplyr::all_of(c("Effect_Order", "Effect_Family", "Effect_Species")), as_number)
    )
  effect_specs <- list(
    list(keys = c("Class", "Order"), value = "Effect_Order"),
    list(keys = c("Class", "Family"), value = "Effect_Family"),
    list(keys = c("Class", "Species"), value = "Effect_Species")
  )
  for (spec in effect_specs) {
    assert(
      all(is.finite(out[[spec$value]])),
      paste0("random_effects.csv contains nonfinite ", spec$value, " values.")
    )
    conflicts <- out |>
      dplyr::group_by(dplyr::across(dplyr::all_of(spec$keys))) |>
      dplyr::summarise(n_values = dplyr::n_distinct(.data[[spec$value]]), .groups = "drop") |>
      dplyr::filter(n_values > 1L)
    assert(
      nrow(conflicts) == 0L,
      paste0("random_effects.csv maps at least one ", paste(spec$keys, collapse = "/"),
             " key to multiple ", spec$value, " values.")
    )
  }
  out
}

read_species_core_inputs <- function(config) {
  paths <- config$paths
  summary_raw <- readr::read_csv(paths$raw$summary, show_col_types = FALSE, progress = FALSE)
  need_cols(
    summary_raw,
    c("scientificName", "className", "orderName", "familyName", "genusName", "speciesName", "redlistCategory"),
    "simple_summary.csv"
  )
  summary_df <- summary_raw |>
    dplyr::mutate(
      row_id = dplyr::row_number(),
      scientificName = stringr::str_squish(scientificName),
      taxon_class = dplyr::case_when(
        clean_name(className) == "mammalia" ~ "Mammalia",
        clean_name(className) == "aves" ~ "Aves",
        TRUE ~ NA_character_
      )
    )
  assert(
    all(!is.na(summary_df$scientificName) & nzchar(summary_df$scientificName)),
    "simple_summary.csv contains missing or blank scientificName values."
  )
  assert(
    !anyDuplicated(clean_name(summary_df$scientificName)),
    "simple_summary.csv contains duplicate normalized scientificName values."
  )
  bad_class <- is.na(summary_df$taxon_class)
  unknown_classes <- sort(unique(dplyr::coalesce(as.character(summary_df$className[bad_class]), "<missing>")))
  assert(
    length(unknown_classes) == 0L,
    paste0("simple_summary.csv contains unsupported className values: ", paste(unknown_classes, collapse = ", "))
  )

  input_synonyms <- read_input_synonyms(paths$raw$input_synonyms)

  list(
    summary_df = summary_df,
    synonyms = read_synonyms(paths$raw$synonyms),
    raster_manifest = if (identical(config$sdm_input_mode, "index")) {
      read_raster_index(
        config$sdm_index_file,
        species_table_sdm_methods(config$sdm),
        unique(paths$rasters$taxon_class)
      )
    } else list_raster_manifest(paths$rasters),
    traits = dplyr::bind_rows(
      read_traits(paths$raw$mammal_traits, "Mammalia"),
      read_traits(paths$raw$bird_traits, "Aves")
    ),
    random_effects = read_random_effects(paths$raw$random_effects),
    bird_generation_lengths = read_bird_generation_lengths(paths$raw$bird_generation_lengths),
    input_synonyms = input_synonyms,
    models = if (config$fit_mammal_curves || config$fit_bird_curves) {
      readRDS(paths$clean$gompertz_models)
    } else NULL
  )
}
